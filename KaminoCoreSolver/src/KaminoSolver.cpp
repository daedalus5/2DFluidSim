# include "../include/KaminoQuantity.h"
# include <boost/math/tools/roots.hpp>

// CONSTRUCTOR / DESTRUCTOR >>>>>>>>>>


KaminoSolver::KaminoSolver(size_t nPhi, size_t nTheta, fReal radius, fReal gridLength, fReal frameDuration) :
	nPhi(nPhi), nTheta(nTheta), radius(radius), gridLen(gridLength), frameDuration(frameDuration),
	timeStep(0.0), timeElapsed(0.0)
{
	addAttr("u", 0.5, 0.0);		// u velocity
	addAttr("v", 0.0, 0.5);		// v velocity
	addAttr("p");				// p pressure
	addAttr("test");			// test scalar field

	this->gridTypes = new gridType[nPhi * nTheta];
	memset(reinterpret_cast<void*>(this->gridTypes), FLUIDGRID, nPhi * nTheta);

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
	// geometric();
	// bodyForce();
	projection();
}

void KaminoSolver::advection()
{
	for (auto quantity : this->attributeTable)
	{
		KaminoQuantity* attr = quantity.second;
		for (size_t gridX = 0; gridX < this->nPhi; ++gridX)
		{
			for (size_t gridY = 0; gridY < this->nTheta; ++gridY)
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

template <class T>
struct cubicFunctor
{
	T G_val;
	T u_prev;
	T operator()(T u)
	{
		return G_val * G_val * pow(u, 3.0) + (G_val * u_prev + 1.0) * u - u_prev;
	}
};

template <class T>
T cbrt_noderiv(T x, cubicFunctor<T> functor)
{
	// return cube root of x using bracket_and_solve (no derivatives).
	using namespace std;                          // Help ADL of std functions.
	using namespace boost::math::tools;           // For bracket_and_solve_root.

	int exponent;
	frexp(x, &exponent);                          // Get exponent of z (ignore mantissa).
	T guess = ldexp(1.0, exponent / 3.0);         // Rough guess is to divide the exponent by three.
	T factor = 2.0;                               // How big steps to take when searching.
	
	const boost::uintmax_t maxit = 20.0;          // Limit to maximum iterations.
	boost::uintmax_t it = maxit;                  // Initally our chosen max iterations, but updated with actual.
	bool is_rising = true;                        // So if result if guess^3 is too low, then try increasing guess.
	int digits = std::numeric_limits<T>::digits;  // Maximum possible binary digits accuracy for type T.
												  // Some fraction of digits is used to control how accurate to try to make the result.
	int get_digits = digits - 3.0;                // We have to have a non-zero interval at each step, so
                                                  // maximum accuracy is digits - 1.  But we also have to
                                                  // allow for inaccuracy in f(x), otherwise the last few
                                                  // iterations just thrash around.
	eps_tolerance<T> tol(get_digits);             // Set the tolerance.
	std::pair<T, T> r = bracket_and_solve_root(functor, guess, factor, is_rising, tol, it);
	return r.first + (r.second - r.first) / 2.0;  // Midway between brackets is our result, if necessary we could
												  // return the result as an interval here.
}

void KaminoSolver::geometric()
{
	KaminoQuantity* u = attributeTable["u"];
	KaminoQuantity* v = attributeTable["v"];

	// Poles unshifted
	for (size_t phiI = 0; phiI < nPhi; ++phiI)
	{
		size_t northPole = 0;
		size_t southPole = nTheta - 1;

		u->writeValueTo(phiI, northPole, u->getValueAt(phiI, northPole));
		u->writeValueTo(phiI, southPole, u->getValueAt(phiI, southPole));
		v->writeValueTo(phiI, northPole, v->getValueAt(phiI, northPole));
		v->writeValueTo(phiI, southPole, v->getValueAt(phiI, southPole));
	}

	for (size_t thetaJ = 1; thetaJ < nTheta - 1; ++thetaJ)
	{
		for (size_t phiI = 0; phiI < nPhi; ++phiI)
		{
			fReal G = timeStep * std::cos(thetaJ * gridLen) / (radius * sin(thetaJ * gridLen));
			
			cubicFunctor<fReal> c;
			c.G_val = G;
			fReal uPrev = u->getValueAt(phiI, thetaJ);
			c.u_prev = uPrev;
			fReal uNext = cbrt_noderiv<fReal>(uPrev, c);

			fReal vPrev = v->getValueAt(phiI, thetaJ);
			fReal vNext = vPrev + G * uNext * uNext;

			u->writeValueTo(phiI, thetaJ, uNext);
			v->writeValueTo(phiI, thetaJ, vNext);
		}
	}

	u->swapBuffer();
	v->swapBuffer();
}

void KaminoSolver::projection()
{
	const fReal density = 1000.0;
	fReal rhsScaleB = -gridLen * radius * density / timeStep;
	fReal scaleP = 1.0 / rhsScaleB;

	Eigen::VectorXd b(nPhi * nTheta);
	b.setZero();

	KaminoQuantity* u = attributeTable["u"];
	KaminoQuantity* v = attributeTable["v"];
	KaminoQuantity* p = attributeTable["p"];

	for (size_t j = 0; j < nTheta; ++j)
	{
		for (size_t i = 0; i < nPhi; ++i)
		{
			// the (unscaled) divergence at grid i, j
			fReal bij = 0.0;
			// oot : one over two
			// a grid's adjacent neighbours along x axis will always be valid
			fReal uPlus, uMinus, vPlus, vMinus;

			size_t ipoot = (i + 1) % nPhi;
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
			if (j != nTheta && getGridTypeAt(i, jpoot) == FLUIDGRID)
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

	Eigen::VectorXd pVector(nPhi * nTheta);
	
	Eigen::ConjugateGradient<Eigen::SparseMatrix<fReal>, Eigen::Lower | Eigen::Upper> cg;
	//cg.setTolerance(pow(10, -1));
	cg.compute(Laplacian);
	pVector = cg.solve(b);

	//std::cout << "#iterations:     " << cg.iterations() << std::endl;
	//std::cout << "estimated error: " << cg.error()      << std::endl;

	// Populate updated pressure values
	for (size_t j = 0; j < nTheta; ++j) 
	{
		for (size_t i = 0; i < nPhi; ++i) 
		{
			p->writeValueTo(i, j, pVector(getIndex(i, j)));
		}
	}
	p->swapBuffer();

	const fReal usolid = 0.0;
	const fReal vsolid = 0.0;
	// V is nPhi by nTheta + 1
	for (size_t j = 0; j < nTheta; ++j)
	{
		for (size_t i = 0; i < nPhi; ++i)
		{
			fReal uBeforeUpdate = u->getValueAt(i, j);
			fReal vBeforeUpdate = v->getValueAt(i, j);
			size_t iRhs = i;
			size_t iLhs = (i == 0 ? nPhi - 1 : i - 1);
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
			
			if (j != 0)
			{
				fReal invSin = 1 / sin(j * gridLen);
				size_t jUpper = j;
				size_t jLower = j - 1;
				fReal pressureSummedV = p->getValueAt(i, jUpper) - p->getValueAt(i, jLower);
				if (getGridTypeAt(i, jLower) == FLUIDGRID && getGridTypeAt(i, jUpper) == FLUIDGRID)
				{
					fReal deltaV = scaleP * invSin * pressureSummedV;
					v->writeValueTo(i, j, vBeforeUpdate + deltaV);
				}
				else
				{
					v->writeValueTo(i, j, vsolid);
				}
			}
			else
			{
				size_t jUpper = j;
				fReal pressureSummedV = p->getValueAt(i, jUpper) - 0.0;
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
		}
	}
	// j = nTheta case
	for (size_t i = 0; i < nPhi; ++i)
	{
		size_t j = nTheta;
		size_t jLower = j - 1;
		fReal vBeforeUpdate = v->getValueAt(i, j);
		fReal pressureSummedV = 0.0 - p->getValueAt(i, jLower);
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

	u->swapBuffer();
	v->swapBuffer();
}


// <<<<<<<<<<
// INITIALIZATION >>>>>>>>>>

/* Duplicate of getIndex() in KaminoQuantity */
size_t KaminoSolver::getIndex(size_t x, size_t y)
{
	return y * nPhi + x;
}

gridType KaminoSolver::getGridTypeAt(size_t x, size_t y)
{
	return gridTypes[getIndex(x, y)];
}

/* Compute Laplacian done right...probably */
void KaminoSolver::precomputeLaplacian()
{
	Laplacian = Eigen::SparseMatrix<fReal>(nPhi*nTheta, nPhi*nTheta);
	Laplacian.setZero();

	for (size_t j = 0; j < nTheta; ++j)
	{
		for (size_t i = 0; i < nPhi; ++i)
		{
			size_t numPhiNeighbors = 0;
			size_t numThetaNeighbors = 0;
			size_t rowNumber = getIndex(i, j);
			size_t ip1 = (i + 1) % nPhi;
			size_t im1 = (i == 0 ? nPhi - 1 : i - 1);

			if(getGridTypeAt(ip1, j) == FLUIDGRID){
				Laplacian.coeffRef(rowNumber, getIndex(ip1, j)) = -1;
				numPhiNeighbors++;
			}
			if(getGridTypeAt(im1, j) == FLUIDGRID){
				Laplacian.coeffRef(rowNumber, getIndex(im1, j)) = -1;
				numPhiNeighbors++;				
			}

			fReal invSin;
			if(j == 0 || j == 2*M_PI){
				invSin = 1.0;
			}
			else{
				invSin = 1 / sin(j * gridLen);
			}

			if (j != nTheta - 1)
			{
				size_t jp1 = j + 1;
				if (getGridTypeAt(i, jp1) == FLUIDGRID)
				{
					Laplacian.coeffRef(rowNumber, getIndex(i, jp1)) = -1 * invSin;
					numThetaNeighbors++;
				}
			}
			if (j != 0)
			{
				size_t jm1 = j - 1;
				if (getGridTypeAt(i, jm1) == FLUIDGRID)
				{
					Laplacian.coeffRef(rowNumber, getIndex(i, jm1)) = -1 * invSin;
					numThetaNeighbors++;
				}
			}
			Laplacian.coeffRef(rowNumber, getIndex(i, j)) = numPhiNeighbors + numThetaNeighbors * invSin;
		}
	}
}

void KaminoSolver::initialize_pressure()
{
	for(size_t i = 0; i < nPhi; ++i){
		for(size_t j = 0; j < nTheta; ++j){
			attributeTable["p"]->setValueAt(i, j, 0.0);
		}
	}
}

void KaminoSolver::initialize_velocity()
{
	fReal val = 0.0;
	size_t sizePhi = attributeTable["u"]->getNPhi();
	size_t sizeTheta = attributeTable["u"]->getNTheta();
	for (size_t j = 0; j < sizeTheta; ++j) {
		for (size_t i = 0; i < sizePhi; ++i) {
			val = FBM(sin(i * gridLen), sin(j * gridLen));
			attributeTable["u"]->setValueAt(i, j, val);
		}
	}
	sizePhi = attributeTable["v"]->getNPhi();
	sizeTheta = attributeTable["v"]->getNTheta();
	for (size_t j = 0; j < sizeTheta; ++j) {
		for (size_t i = 0; i < sizePhi; ++i) {
			val = FBM(cos(i * gridLen), cos(j * gridLen));
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
	for(size_t i = 0; i < nPhi; ++i){
		for(size_t j = 0; j < nTheta; ++j){
			attributeTable["test"]->setValueAt(i, j, 0.0);
		}
	}

	KaminoQuantity* test = attributeTable["test"];
	size_t midX = nPhi / 2;
	size_t midY = nTheta / 2;
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
			test->setValueAt(i + midX, j + midY, Gaussian(i,j));
		}
	}
}

void KaminoSolver::initialize_boundary()
{
	for (size_t gridX = 0; gridX != this->nPhi; ++gridX)
	{
		this->gridTypes[getIndex(gridX, 0)] = SOLIDGRID;
		this->gridTypes[getIndex(gridX, this->nTheta - 1)] = SOLIDGRID;
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

	for (size_t j = 0; j < nTheta; ++j) {
		for (size_t i = 0; i < nPhi; ++i) {
			uLeft = u->getValueAt(i, j);
			i == (nPhi - 1) ? upi = 0 : upi = i + 1;
			vDown = v->getValueAt(i, j);
			j == (nTheta - 1) ? vpi = 0 : vpi = j + 1;
			uRight = u->getValueAt(upi, j);
			vUp = u->getValueAt(i, vpi);

			velX = (uLeft + uRight) / 2.0;
			velY = (vUp + vDown) / 2.0;

			pos = Eigen::Matrix<float, 3, 1>(i * gridLen, j * gridLen, 0.0);
			vel = Eigen::Matrix<float, 3, 1>(0.0, velY, velX);
			mapVToSphere(pos, vel);
			mapPToSphere(pos);

			pressure = attributeTable["p"]->getValueAt(i, j);
			testVal = attributeTable["test"]->getValueAt(i, j);
			
			int idx = parts->addParticle();
			float* p = parts->dataWrite<float>(pH, idx);
			float* v = parts->dataWrite<float>(vH, idx);
			float* ps = parts->dataWrite<float>(psH, idx);
			float* ts = parts->dataWrite<float>(test, idx);

			ps[0] = pressure / 5000.0;
			ts[0] = testVal / 1.0;

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

void KaminoSolver::mapPToSphere(Eigen::Matrix<float, 3, 1>& pos) const
{
	float theta = pos[1];
	float phi = pos[0];
	pos[0] = radius * sin(theta) * cos(phi);
	pos[2] = radius * sin(theta) * sin(phi);
	pos[1] = radius * cos(theta);
}

void KaminoSolver::mapVToSphere(Eigen::Matrix<float, 3, 1>& pos, Eigen::Matrix<float, 3, 1>& vel) const
{
	float theta = pos[1];
	float phi = pos[0];

	float u_theta = vel[1];
	float u_phi = vel[2];

	vel[0] = cos(theta) * cos(phi) * u_theta - sin(phi) * u_phi;
	vel[2] = cos(theta) * sin(phi) * u_theta + cos(phi) * u_phi;
	vel[1] = -sin(theta) * u_theta;
}

void KaminoSolver::mapToCylinder(Eigen::Matrix<float, 3, 1>& pos) const
{
	float radius = 5.0;
	float phi = 2*M_PI*pos[0] / (nPhi * gridLen);
	float z = pos[1];
	pos[0] = radius * cos(phi);
	pos[1] = radius * sin(phi);
	pos[2] = z;
}

// <<<<<<<<<<
// ACCESS >>>>>>>>>>


void KaminoSolver::addAttr(std::string name, fReal xOffset, fReal yOffset)
{
	size_t attrnPhi = this->nPhi;
	size_t attrnTheta = this->nTheta;
	if (name == "v")
		attrnTheta += 1;
	KaminoQuantity* ptr = new KaminoQuantity(name, attrnPhi, attrnTheta, this->gridLen, xOffset, yOffset);
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