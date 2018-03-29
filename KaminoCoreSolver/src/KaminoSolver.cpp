# include "../include/KaminoQuantity.h"


// CONSTRUCTOR / DESTRUCTOR >>>>>>>>>>


KaminoSolver::KaminoSolver(size_t nx, size_t ny, fReal gridLength, fReal frameDuration) :
	nx(nx), ny(ny), gridLen(gridLength), frameDuration(frameDuration),
	timeStep(0.0), timeElapsed(0.0)
{
	addAttr("u", 0.5, 0.0);		// u velocity
	addAttr("v", 0.0, 0.5);		// v velocity
	addAttr("p");				// p pressure
	addAttr("test");			// test scalar field

	this->gridTypes = new gridType[nx * ny];
	memset(reinterpret_cast<void*>(this->gridTypes), FLUIDGRID, nx * ny);

	initialize_velocity();
	initialize_pressure();
	precomputeLaplacian();
	initialize_test();

	initialize_boundary();
}

KaminoSolver::~KaminoSolver()
{
	for (auto& attr : this->attributeTable)
	{
		delete attr.second;
	}
	delete[] this->gridTypes;
}


// <<<<<<<<<<
// CORE FLUID SOLVER >>>>>>>>>>


void KaminoSolver::stepForward(fReal timeStep)
{
	this->timeStep = timeStep;
	advection();
	// bodyForce();
	projection();
}

void KaminoSolver::advection()
{
	for (auto quantity : this->attributeTable)
	{
		KaminoQuantity* attr = quantity.second;
		for (size_t gridX = 0; gridX < this->nx; ++gridX)
		{
			for (size_t gridY = 0; gridY < this->ny; ++gridY)
			{
				fReal gX = attr->getXCoordAtIndex(gridX);
				fReal gY = attr->getYCoordAtIndex(gridY);

				fReal uG = (*this)["u"]->sampleAt(gX, gY);
				fReal vG = (*this)["v"]->sampleAt(gX, gY);

				fReal midX = gX - 0.5 * timeStep * uG;
				fReal midY = gY - 0.5 * timeStep * vG;

				fReal uMid = (*this)["u"]->sampleAt(midX, midY);
				fReal vMid = (*this)["v"]->sampleAt(midX, midY);

				fReal pX = gX - timeStep * uMid;
				fReal pY = gY - timeStep * vMid;
				
				fReal advectedVal = attr->sampleAt(pX, pY);
				attr->writeValueTo(gridX, gridY, advectedVal);
			}
		}
	}
	this->swapAttrBuffers();
}

void KaminoSolver::projection()
{
	const fReal density = 1000.0;
	fReal rhsScaleB = -gridLen * density / timeStep;
	fReal scaleP = 1.0 / rhsScaleB;

	Eigen::VectorXd b(nx * ny);
	b.setZero();

	KaminoQuantity* u = attributeTable["u"];
	KaminoQuantity* v = attributeTable["v"];
	KaminoQuantity* p = attributeTable["p"];

	for (size_t j = 0; j < ny; ++j)
	{
		for (size_t i = 0; i < nx; ++i)
		{
			// the (unscaled) divergence at grid i, j
			fReal bij = 0.0;
			// oot : one over two
			// a grid's adjacent neighbours along x axis will always be valid
			fReal uPlus, uMinus, vPlus, vMinus;

			size_t ipoot = (i + 1) % nx;
			if (getGridTypeAt(ipoot, j) == FLUIDGRID)
			{
				uPlus = u->getValueAt(ipoot, j);
			}
			else
			{
				uPlus = 0.0;
			}

			size_t imoot = i;
			if (getGridTypeAt(imoot, j) == FLUIDGRID)
			{
				uMinus = u->getValueAt(imoot, j);
			}
			else
			{
				uMinus = 0.0;
			}
			
			// but that's not the case for y axis
			size_t jpoot = j + 1;
			if (j != ny && getGridTypeAt(i, jpoot) == FLUIDGRID)
			{
				vPlus = v->getValueAt(i, jpoot);
			}
			else
			{
				vPlus = 0.0;
			}

			size_t jmoot = j;
			if (getGridTypeAt(i, jmoot) == FLUIDGRID)
			{
				vMinus = v->getValueAt(i, jmoot);
			}
			else
			{
				vMinus = 0.0;
			}

			b(getIndex(i, j)) = (uPlus - uMinus + vPlus - vMinus);
		}
	}
	b = b * rhsScaleB;

	Eigen::VectorXd pVector(nx * ny);
	
	Eigen::ConjugateGradient<Eigen::SparseMatrix<fReal>, Eigen::Lower | Eigen::Upper> cg;
	//cg.setTolerance(pow(10, -1));
	cg.compute(Laplacian);
	pVector = cg.solve(b);

	//std::cout << "#iterations:     " << cg.iterations() << std::endl;
	//std::cout << "estimated error: " << cg.error()      << std::endl;

	// Populate updated pressure values
	for (size_t j = 0; j < ny; ++j) 
	{
		for (size_t i = 0; i < nx; ++i) 
		{
			p->writeValueTo(i, j, pVector(getIndex(i, j)));
		}
	}
	p->swapBuffer();

	// V is nx by ny + 1
	for (size_t j = 0; j < ny + 1; ++j)
	{
		for (size_t i = 0; i < nx; ++i)
		{
			const fReal usolid = 0.0;
			const fReal vsolid = 0.0;

			fReal uBeforeUpdate = u->getValueAt(i, j);
			fReal vBeforeUpdate = v->getValueAt(i, j);
			size_t iRhs = i;
			size_t iLhs = (i == 0 ? nx - 1 : i - 1);
			if (getGridTypeAt(iRhs, j) == FLUIDGRID && getGridTypeAt(iLhs, j) == FLUIDGRID)
			{
				fReal pressureSummedU = p->getValueAt(iRhs, j) - p->getValueAt(iLhs, j);
				fReal deltaU = scaleP * pressureSummedU;
				u->writeValueTo(i, j, uBeforeUpdate + deltaU);
			}
			else
			{
				u->writeValueTo(i, j, usolid);
			}
			
			if (j != ny && j != 0)
			{
				size_t jUpper = j;
				size_t jLower = j - 1;
				fReal pressureSummedV = p->getValueAt(i, jUpper) - p->getValueAt(i, jLower);
				if (getGridTypeAt(i, jLower) == FLUIDGRID && getGridTypeAt(i, jUpper) == FLUIDGRID)
				{
					fReal deltaV = scaleP * pressureSummedV;
					v->writeValueTo(i, j, vBeforeUpdate + deltaV);
				}
				else
				{
					v->writeValueTo(i, j, vsolid);
				}
			}
			else
			{
				if (j == 0)
				{
					size_t jUpper = j;
					fReal pressureSummedV = p->getValueAt(i, jUpper);
					if (getGridTypeAt(i, jUpper) == FLUIDGRID)
					{
						fReal deltaV = scaleP * pressureSummedV;
						v->writeValueTo(i, j, vBeforeUpdate + deltaV);
					}
					else
					{
						v->writeValueTo(i, j, vsolid);
					}
				}
				if (j == ny)
				{
					size_t jLower = j - 1;
					fReal pressureSummedV = p->getValueAt(i, jLower);
					if (getGridTypeAt(i, jLower) == FLUIDGRID)
					{
						fReal deltaV = scaleP * pressureSummedV;
						v->writeValueTo(i, j, vBeforeUpdate + deltaV);
					}
					else
					{
						v->writeValueTo(i, j, vsolid);
					}
				}
			}
		}
	}

	u->swapBuffer();
	v->swapBuffer();
}


// <<<<<<<<<<
// INITIALIZATION >>>>>>>>>>

/* Duplicate of getIndex() in KaminoQuantity */
size_t KaminoSolver::getIndex(size_t x, size_t y)
{
	return y * nx + x;
}

gridType KaminoSolver::getGridTypeAt(size_t x, size_t y)
{
	return gridTypes[getIndex(x, y)];
}

/* Compute Laplacian done right...probably */
void KaminoSolver::precomputeLaplacian()
{
	Laplacian = Eigen::SparseMatrix<fReal>(nx*ny, nx*ny);
	Laplacian.setZero();

	for (size_t j = 0; j < ny; ++j)
	{
		for (size_t i = 0; i < nx; ++i)
		{
			size_t numNeighbors = 0;
			size_t rowNumber = getIndex(i, j);
			size_t ip1 = (i + 1) % nx;
			size_t im1 = (i == 0 ? nx - 1 : i - 1);

			if(getGridTypeAt(ip1, j) == FLUIDGRID){
				Laplacian.coeffRef(rowNumber, getIndex(ip1, j)) = -1;
				numNeighbors++;
			}
			if(getGridTypeAt(im1, j) == FLUIDGRID){
				Laplacian.coeffRef(rowNumber, getIndex(im1, j)) = -1;
				numNeighbors++;				
			}

			if (j != ny - 1)
			{
				size_t jp1 = j + 1;
				if (getGridTypeAt(i, jp1) == FLUIDGRID)
				{
					Laplacian.coeffRef(rowNumber, getIndex(i, jp1)) = -1;
					numNeighbors++;
				}
			}
			if (j != 0)
			{
				size_t jm1 = j - 1;
				if (getGridTypeAt(i, jm1) == FLUIDGRID)
				{
					Laplacian.coeffRef(rowNumber, getIndex(i, jm1)) = -1;
					numNeighbors++;
				}
			}
			Laplacian.coeffRef(rowNumber, getIndex(i, j)) = numNeighbors;
		}
	}
}

void KaminoSolver::initialize_pressure()
{
	for(size_t i = 0; i < nx; ++i){
		for(size_t j = 0; j < ny; ++j){
			attributeTable["p"]->setValueAt(i, j, 0.0);
		}
	}
}

void KaminoSolver::initialize_velocity()
{
	fReal val = 0.0;
	size_t sizeX = attributeTable["u"]->getNx();
	size_t sizeY = attributeTable["u"]->getNy();
	for (size_t j = 0; j < sizeY; ++j) {
		for (size_t i = 0; i < sizeX; ++i) {
			val = FBM(sin(2 * M_PI * i / nx), sin(2 * M_PI * j / ny));
			attributeTable["u"]->setValueAt(i, j, val);
		}
	}
	sizeX = attributeTable["v"]->getNx();
	sizeY = attributeTable["v"]->getNy();
	for (size_t j = 0; j < sizeY; ++j) {
		for (size_t i = 0; i < sizeX; ++i) {
			val = FBM(cos(2 * M_PI * i / nx), cos(2 * M_PI * j / ny));
			attributeTable["v"]->setValueAt(i, j, 0.0);
		}
	}
}

fReal KaminoSolver::FBM(const fReal x, const fReal y) {
	fReal total = 0.0f;
	fReal resolution = 1.0;
	fReal persistance = 0.5;
	int octaves = 4;

	for (int i = 0; i < octaves; i++) {
		fReal freq = std::pow(2.0f, i);
		fReal amp = std::pow(persistance, i);
		total += amp * interpNoise2D(x * freq / resolution, y * freq / resolution);
	}
	fReal a = 1 - persistance;  // normalization

	return a * total / 2.0f;  // normalized, pseudorandom number between -1 and 1
}

fReal KaminoSolver::interpNoise2D(const fReal x, const fReal y) const {
	fReal intX = std::floor(x);
	fReal fractX = x - intX;
	fReal intY = std::floor(y);
	fReal fractY = y - intY;

	fReal v1 = rand(Eigen::Matrix<fReal, 2, 1>(intX, intY));
	fReal v2 = rand(Eigen::Matrix<fReal, 2, 1>(intX + 1, intY));
	fReal v3 = rand(Eigen::Matrix<fReal, 2, 1>(intX, intY + 1));
	fReal v4 = rand(Eigen::Matrix<fReal, 2, 1>(intX + 1, intY + 1));

	// interpolate for smooth transitions
	fReal i1 = KaminoLerp(v1, v2, fractX);
	fReal i2 = KaminoLerp(v3, v4, fractX);
	return KaminoLerp(i1, i2, fractY);
}

fReal KaminoSolver::rand(const Eigen::Matrix<fReal, 2, 1> vecA) const {
	// return pseudorandom number between -1 and 1
	Eigen::Matrix<fReal, 2, 1> vecB = Eigen::Matrix<fReal, 2, 1>(12.9898, 4.1414);
	fReal val = sin(vecA.dot(vecB) * 43758.5453);
	return val - std::floor(val);
}

void KaminoSolver::initialize_test()
{
	for(size_t i = 0; i < nx; ++i){
		for(size_t j = 0; j < ny; ++j){
			attributeTable["test"]->setValueAt(i, j, 0.0);
		}
	}

	KaminoQuantity* test = attributeTable["test"];
	size_t midX = nx / 2;
	size_t midY = ny / 2;
	size_t kernelSize = 11;
	Eigen::Matrix<fReal, 11, 11> Gaussian;
	Gaussian << 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
				1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1,
				1, 2, 3, 3, 3, 3, 3, 3, 3, 2, 1,
				1, 2, 3, 5, 5, 5, 5, 5, 3, 2, 1,
				1, 2, 3, 5, 8, 8, 8, 5, 3, 2, 1,
				1, 2, 3, 5, 8, 13, 8, 5, 3, 2, 1,
				1, 2, 3, 5, 8, 8, 8, 5, 3, 2, 1,
				1, 2, 3, 5, 5, 5, 5, 5, 3, 2, 1,
				1, 2, 3, 3, 3, 3, 3, 3, 3, 2, 1,
				1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1,
				1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1;
	for(size_t i = 0; i < kernelSize; ++i){
		for(size_t j = 0; j < kernelSize; ++j){
			test->setValueAt(i, j, 10*Gaussian(i,j));
		}
	}
}

void KaminoSolver::initialize_boundary()
{
	for (size_t gridX = 0; gridX != this->nx; ++gridX)
	{
		this->gridTypes[getIndex(gridX, 0)] = SOLIDGRID;
		this->gridTypes[getIndex(gridX, this->ny - 1)] = SOLIDGRID;
	}
}

// <<<<<<<<<<
// OUTPUT >>>>>>>>>>


void KaminoSolver::write_data_bgeo(const std::string& s, const int frame)
{
# ifndef _MSC_VER
	std::string file = s + std::to_string(frame) + ".bgeo";

	Partio::ParticlesDataMutable* parts = Partio::create();
	Partio::ParticleAttribute pH, vH, psH, test;
	pH = parts->addAttribute("position", Partio::VECTOR, 3);
	vH = parts->addAttribute("v", Partio::VECTOR, 3);
	psH = parts->addAttribute("pressure", Partio::VECTOR, 1);
	test = parts->addAttribute("test", Partio::VECTOR, 1);

	Eigen::Matrix<float, 3, 1> pos;
	Eigen::Matrix<float, 3, 1> vel;
	fReal pressure, testVal;
	fReal velX, velY;

	KaminoQuantity* u = attributeTable["u"];
	KaminoQuantity* v = attributeTable["v"];
	fReal uRight, uLeft, vUp, vDown;

	size_t upi, vpi;

	for (size_t j = 0; j < ny; ++j) {
		for (size_t i = 0; i < ny; ++i) {
			uLeft = u->getValueAt(i, j);
			i == (nx - 1) ? upi = 0 : upi = i + 1;
			vDown = v->getValueAt(i, j);
			j == (ny - 1) ? vpi = 0 : vpi = j + 1;
			uRight = u->getValueAt(upi, j);
			vUp = u->getValueAt(i, vpi);

			velX = (uLeft + uRight) / 2.0;
			velY = (vUp + vDown) / 2.0;

			pos = Eigen::Matrix<float, 3, 1>(i * gridLen, j * gridLen, 0.0);
			//mapToSphere(pos);
			mapToCylinder(pos);
			vel = Eigen::Matrix<float, 3, 1>(velX, velY, 0.0);
			pressure = attributeTable["p"]->getValueAt(i, j);
			testVal = attributeTable["test"]->getValueAt(i, j);
			
			int idx = parts->addParticle();
			float* p = parts->dataWrite<float>(pH, idx);
			float* v = parts->dataWrite<float>(vH, idx);
			float* ps = parts->dataWrite<float>(psH, idx);
			float* ts = parts->dataWrite<float>(test, idx);
			ps[0] = pressure;
			ts[0] = testVal;

			for (int k = 0; k < 3; ++k) {
				p[k] = pos(k, 0);
				v[k] = vel(k, 0);
			}
		}
	}
	Partio::write(file.c_str(), *parts);
	parts->release();
# endif
}

void KaminoSolver::mapToSphere(Eigen::Matrix<float, 3, 1>& pos) const
{
	float radius = 5.0;
	float theta = M_PI*pos[1] / (ny * gridLen);
	float phi = 2*M_PI*pos[0] / (nx * gridLen);
	pos[0] = radius * sin(theta) * cos(phi);
	pos[1] = radius * cos(theta);
	pos[2] = radius * sin(theta) * sin(phi);
}

void KaminoSolver::mapToCylinder(Eigen::Matrix<float, 3, 1>& pos) const
{
	float radius = 5.0;
	float phi = 2*M_PI*pos[0] / (nx * gridLen);
	float z = pos[1];
	pos[0] = radius * cos(phi);
	pos[1] = radius * sin(phi);
	pos[2] = z;
}

// <<<<<<<<<<
// ACCESS >>>>>>>>>>


void KaminoSolver::addAttr(std::string name, fReal xOffset, fReal yOffset)
{
	size_t attrNx = this->nx;
	size_t attrNy = this->ny;
	if (name == "v")
		attrNy += 1;
	KaminoQuantity* ptr = new KaminoQuantity(name, attrNx, attrNy, this->gridLen, xOffset, yOffset);
	this->attributeTable.emplace(std::pair<std::string, KaminoQuantity*>(name, ptr));
}

KaminoQuantity* KaminoSolver::getAttributeNamed(std::string name)
{
	return (*this)[name];
}

KaminoQuantity* KaminoSolver::operator[](std::string name)
{
	return this->attributeTable.at(name);
}

void KaminoSolver::swapAttrBuffers()
{
	for (auto quantity : this->attributeTable)
	{
		quantity.second->swapBuffer();
	}
}


// <<<<<<<<<<