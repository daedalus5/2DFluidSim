# include "../include/KaminoQuantity.h"


// CONSTRUCTOR / DESTRUCTOR >>>>>>>>>>


KaminoSolver::KaminoSolver(size_t nx, size_t ny, fReal gridLength, fReal frameDuration) :
	nx(nx), ny(ny), gridLen(gridLength), frameDuration(frameDuration),
	timeStep(0.0), timeElapsed(0.0)
{
	addAttr("u", 0.5, 0.0);		// u velocity
	addAttr("v", 0.0, 0.5);		// v velocity
	addAttr("p");				// p density

	initialize_velocity();
	initialize_pressure();
	precomputeLaplacian();
}

KaminoSolver::~KaminoSolver()
{
	for (auto& attr : this->attributeTable)
	{
		delete attr.second;
	}
}


// <<<<<<<<<<
// CORE FLUID SOLVER >>>>>>>>>>


void KaminoSolver::stepForward(fReal timeStep)
{
	this->timeStep = timeStep;
	advection();
	// bodyForce();
	// projection();
# ifdef DEBUGBUILD
	/*for (unsigned gridX = 0; gridX != nx; ++gridX)
	{
		for (unsigned gridY = 0; gridY != ny; ++gridY)
		{
			std::cout << attributeTable["u"]->getValueAt(gridX, gridY) << '\t';
		}
		std::cout << '\n';
	}*/
# endif
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

/*void KaminoSolver::projection()
{
	// fReal density = 1000;	// rest fluid density
	// fReal scale = timeStep / (density * gridLen * gridLen);
	// fReal invGridLen = 1 / gridLen;
	fReal density = 1.0;
	fReal scale = 1.0;
	fReal invGridLen = 1.0;

	// construct the vector b
	Eigen::VectorXd b(nx * ny);
	b.setZero();
	for(size_t j = 0; j < ny; ++j){
		for(size_t i = 0; i < nx; ++i){
			fReal uPlus, uMinus, vPlus, vMinus;
			i > (nx - 2) ? (uPlus = attributeTable["u"]->getValueAt(0, j)) : (uPlus = attributeTable["u"]->getValueAt(i + 1, j));
			uMinus = attributeTable["u"]->getValueAt(i, j); 
			j > (ny - 2) ? (vPlus = attributeTable["v"]->getValueAt(i, 0)) : (vPlus = attributeTable["v"]->getValueAt(i, j + 1));
			vMinus = attributeTable["v"]->getValueAt(i, j);
			b(j*nx + i) = -1 * ((uPlus - uMinus) * invGridLen + (vPlus - vMinus) * invGridLen);
		}
	}

	// pressure vector
	Eigen::VectorXd p(nx * ny);
	p.setZero();

	// solving Ax = b
	Eigen::ConjugateGradient<Eigen::SparseMatrix<fReal>, Eigen::Lower|Eigen::Upper> cg;
	cg.setTolerance(pow(10, -1));
	cg.compute(Laplacian * scale);
	p = cg.solve(b);
	p *= 0.1;

	//std::cout << "#iterations:     " << cg.iterations() << std::endl;
	//std::cout << "estimated error: " << cg.error()      << std::endl;

	// Populate updated pressure values
    for(size_t j = 0; j < ny; ++j){
    	for(size_t i = 0; i < nx; ++i){
    		attributeTable["p"]->writeValueTo(i, j, p(j*nx + i));
    	}	
    }

    // Populate updated u values
    for(size_t j = 0; j < ny; ++j){
    	for(size_t i = 0; i < nx; ++i){
    		if(i == 0){
    			attributeTable["u"]->writeValueTo(i, j, attributeTable["u"]->getValueAt(i, j) -
    			scale * gridLen * (attributeTable["p"]->getNextValueAt(i, j) - attributeTable["p"]->getNextValueAt(nx - 1, j)));
    		}
    		else{
    			attributeTable["u"]->writeValueTo(i, j, attributeTable["u"]->getValueAt(i, j) -
    			scale * gridLen * (attributeTable["p"]->getNextValueAt(i, j) - attributeTable["p"]->getNextValueAt(i - 1, j)));
    		}
    	}
    }

    // Populate updated v values
    for(size_t j = 0; j < ny; ++j){
    	for(size_t i = 0; i < nx; ++i){
    		if(j == 0){
    			attributeTable["v"]->writeValueTo(i, j, attributeTable["v"]->getValueAt(i, j) -
    			scale * gridLen * (attributeTable["p"]->getNextValueAt(i, j) - attributeTable["v"]->getNextValueAt(i, ny - 1)));
    		}
    		else{
    			attributeTable["v"]->writeValueTo(i, j, attributeTable["v"]->getValueAt(i, j) -
    			scale * gridLen * (attributeTable["p"]->getNextValueAt(i, j) - attributeTable["v"]->getNextValueAt(i, j - 1)));
    		}
    	}
    }

    this->swapAttrBuffers();
}*/

void KaminoSolver::projection()
{
	const fReal density = 1000.0;
	fReal rhsScaleB = -gridLen * density / timeStep;
	fReal scaleP = 1.0 / rhsScaleB;

	Eigen::VectorXd b(nx * ny);
	b.setZero();

	for (size_t j = 0; j < ny; ++j)
	{
		for (size_t i = 0; i < nx; ++i)
		{
			// oot : one over two
			size_t ipoot = (i + 1) % nx;
			size_t imoot = i;
			size_t jpoot = (j + 1) & ny;
			size_t jmoot = j;

			fReal uPlus = attributeTable["u"]->getValueAt(ipoot, j);
			fReal uMinus = attributeTable["u"]->getValueAt(imoot, j);
			fReal vPlus = attributeTable["v"]->getValueAt(i, jpoot);
			fReal vMinus = attributeTable["v"]->getValueAt(i, jmoot);

			b(getIndex(i, j)) = (uPlus - uMinus + vPlus - vMinus);
		}
	}
	b = b * rhsScaleB;

	Eigen::VectorXd p(nx * ny);
	
	Eigen::ConjugateGradient<Eigen::SparseMatrix<fReal>, Eigen::Lower | Eigen::Upper> cg;
	cg.setTolerance(pow(10, -1));
	cg.compute(Laplacian);
	p = cg.solve(b);

	//std::cout << "#iterations:     " << cg.iterations() << std::endl;
	//std::cout << "estimated error: " << cg.error()      << std::endl;

	// Populate updated pressure values
	for (size_t j = 0; j < ny; ++j) 
	{
		for (size_t i = 0; i < nx; ++i) 
		{
			attributeTable["p"]->writeValueTo(i, j, p(getIndex(i, j)));
		}
	}
	attributeTable["p"]->swapBuffer();

	for (size_t j = 0; j < ny; ++j)
	{
		for (size_t i = 0; i < nx; ++i)
		{
			size_t iRhs = i;
			size_t iLhs = (i == 0 ? nx - 1 : i - 1);
			size_t jUpper = j;
			size_t jLower = (j == 0 ? ny - 1 : j - 1);

			fReal uBeforeUpdate = attributeTable["u"]->getValueAt(i, j);
			fReal deltaU = scaleP * (attributeTable["p"]->getValueAt(iRhs, j) - attributeTable["p"]->getValueAt(iLhs, j));
			attributeTable["u"]->writeValueTo(i, j, uBeforeUpdate + deltaU);
			
			fReal vBeforeUpdate = attributeTable["v"]->getValueAt(i, j);
			fReal deltaV = scaleP * (attributeTable["p"]->getValueAt(i, jUpper) - attributeTable["p"]->getValueAt(i, jLower));
			attributeTable["v"]->writeValueTo(i, j, vBeforeUpdate + deltaV);
		}
	}

	this->swapAttrBuffers();
}


// <<<<<<<<<<
// INITIALIZATION >>>>>>>>>>


/*void KaminoSolver::precomputeLaplacian()
{
	// present construction assumes 2D fluid in every cell and toroidal BCs
	Laplacian = Eigen::SparseMatrix<fReal>(nx*ny, nx*ny);
	Laplacian.setZero();

	// construct Laplacian row-by-row
	size_t k = 0;

	for(size_t j = 0; j < ny; ++j){
		for(size_t i = 0; i < nx; ++i){
			Laplacian.coeffRef(k, j*nx + i) = 4;
			i > (nx - 2) ? (Laplacian.coeffRef(k, j*nx) = -1) : (Laplacian.coeffRef(k, j*nx + i + 1) = -1);
			i < 1 ? (Laplacian.coeffRef(k, j*nx + nx - 1) = -1) : (Laplacian.coeffRef(k, j*nx + i - 1) = -1);
			j > (ny - 2) ? (Laplacian.coeffRef(k, i) = -1) : (Laplacian.coeffRef(k, (j + 1)*nx + i) = -1);
			j < 1 ? (Laplacian.coeffRef(k, (ny - 1)*nx + i) = -1) : (Laplacian.coeffRef(k, (j - 1)*nx + i) = -1);
			k++;
		}
	}
}*/

size_t KaminoSolver::getIndex(size_t x, size_t y)
{
# ifdef DEBUGBUILD
	if (x >= this->nx || y >= this->ny)
	{
		std::cerr << "Index out of bound at x: " << x << " y: " << y << std::endl;
	}
# endif
	return y * nx + x;
}

void KaminoSolver::precomputeLaplacian()
{
	Laplacian = Eigen::SparseMatrix<fReal>(nx*ny, nx*ny);
	Laplacian.setZero();

	for (size_t j = 0; j < ny; ++j)
	{
		for (size_t i = 0; i < nx; ++i)
		{
			size_t rowNumber = getIndex(i, j);
			size_t ip1 = (i + 1) % nx;
			size_t im1 = (i == 0 ? nx - 1 : i - 1);
			size_t jp1 = (j + 1) % ny;
			size_t jm1 = (j == 0 ? ny - 1 : j - 1);

			Laplacian.coeffRef(rowNumber, getIndex(i, j)) = 4;
			Laplacian.coeffRef(rowNumber, getIndex(ip1, j)) = -1;
			Laplacian.coeffRef(rowNumber, getIndex(i, jp1)) = -1;
			Laplacian.coeffRef(rowNumber, getIndex(im1, j)) = -1;
			Laplacian.coeffRef(rowNumber, getIndex(i, jm1)) = -1;
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
	for (size_t j = 0; j < ny; ++j) {
		for (size_t i = 0; i < nx; ++i) {
			val = FBM(sin(2 * M_PI*i / nx), sin(2 * M_PI*j / ny));
			attributeTable["u"]->setValueAt(i, j, val);
			//val = FBM(cos(2 * M_PI*i / nx), cos(2 * M_PI*j / ny));
			attributeTable["v"]->setValueAt(i, j, 0.0);
		}
	}
}

fReal KaminoSolver::FBM(const fReal x, const fReal y) {
	fReal total = 0.0f;
	fReal resolution = 100.0;
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


// <<<<<<<<<<
// OUTPUT >>>>>>>>>>


void KaminoSolver::write_data_bgeo(const std::string& s, const int frame)
{
# ifndef _MSC_VER
	std::string file = s + std::to_string(frame) + ".bgeo";

	Partio::ParticlesDataMutable* parts = Partio::create();
	Partio::ParticleAttribute pH, vH, psH;
	pH = parts->addAttribute("position", Partio::VECTOR, 3);
	vH = parts->addAttribute("v", Partio::VECTOR, 3);
	psH = parts->addAttribute("pressure", Partio::VECTOR, 1);

	Eigen::Matrix<float, 3, 1> pos;
	Eigen::Matrix<float, 3, 1> vel;
	fReal pressure;
	fReal velX, velY;

	for (size_t j = 0; j < ny; ++j) {
		for (size_t i = 0; i < ny; ++i) {
			if(i == (nx - 1)){
				velX = (attributeTable["u"]->getValueAt(i, j) + attributeTable["u"]->getValueAt(0, j)) / 2.0;
			}
			else{
				velX = (attributeTable["u"]->getValueAt(i, j) + attributeTable["u"]->getValueAt(i + 1, j)) / 2.0;
			}
			if(j == (ny - 1)){
				velY = (attributeTable["v"]->getValueAt(i, j) + attributeTable["v"]->getValueAt(i, 0)) / 2.0;
			}
			else{
				velY = (attributeTable["v"]->getValueAt(i, j) + attributeTable["v"]->getValueAt(i, j + 1)) / 2.0;
			}
			pos = Eigen::Matrix<float, 3, 1>(i * gridLen, j * gridLen, 0.0);
			mapToSphere(pos);
			vel = Eigen::Matrix<float, 3, 1>(velX, velY, 0.0);
			pressure = attributeTable["p"]->getValueAt(i, j);
			int idx = parts->addParticle();
			float* p = parts->dataWrite<float>(pH, idx);
			float* v = parts->dataWrite<float>(vH, idx);
			float* ps = parts->dataWrite<float>(psH, idx);
			ps[0] = pressure;
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

// <<<<<<<<<<
// ACCESS >>>>>>>>>>


void KaminoSolver::addAttr(std::string name, fReal xOffset, fReal yOffset)
{
	KaminoQuantity* ptr = new KaminoQuantity(name, this->nx, this->ny, this->gridLen, xOffset, yOffset);
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