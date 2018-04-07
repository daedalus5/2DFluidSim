# include "../include/KaminoQuantity.h"
# include "../include/CubicSolver.h"

// CONSTRUCTOR / DESTRUCTOR >>>>>>>>>>


KaminoSolver::KaminoSolver(size_t nPhi, size_t nTheta, fReal radius, fReal gridLength, fReal frameDuration) :
	nPhi(nPhi), nTheta(nTheta), radius(radius), gridLen(gridLength), frameDuration(frameDuration),
	timeStep(0.0), timeElapsed(0.0)
{
	addStaggeredAttr("u", 0.0, 0.5);		// u velocity
	addStaggeredAttr("v", 0.5, 0.0);		// v velocity
	addCenteredAttr("p", 0.5, 0.5);				// p pressure
	addCenteredAttr("test", 0.5, 0.5);			// test scalar field

	this->gridTypes = new gridType[nPhi * nTheta];
	memset(reinterpret_cast<void*>(this->gridTypes), FLUIDGRID, nPhi * nTheta);

	initialize_velocity();
	initialize_pressure();
	precomputeLaplacian();
	initialize_test();

	//initialize_boundary();
}

KaminoSolver::~KaminoSolver()
{
	for (auto& attr : this->centeredAttr)
	{
		delete attr.second;
	}
	for (auto& attr : this->staggeredAttr)
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
	advectionScalar();
	advectionSpeed();
	this->swapAttrBuffers();

	// geometric();
	// bodyForce();
	// projection();
}

// Phi: 0 - 2pi  Theta: 0 - pi
void validatePhiTheta(fReal & phi, fReal & theta)
{
	int loops = static_cast<int>(std::floor(theta / M_2PI));
	theta = theta - loops * M_2PI;
	if (theta > M_PI)
	{
		theta = M_2PI - theta;
		phi += M_PI;
	}
	loops = static_cast<int>(std::floor(phi / M_2PI));
	phi = phi - loops * M_2PI;
}

void KaminoSolver::advectAttrAt(KaminoQuantity* attr, size_t gridPhi, size_t gridTheta)
{
	KaminoQuantity* uPhi = (*this)["u"];
	KaminoQuantity* uTheta = (*this)["v"];

	fReal gTheta = attr->getThetaCoordAtIndex(gridTheta);
	fReal gPhi = attr->getPhiCoordAtIndex(gridPhi);

	fReal guPhi = uPhi->sampleAt(gPhi, gTheta);
	fReal guTheta = uTheta->sampleAt(gPhi, gTheta);

	fReal latRadius = this->radius * std::sin(gTheta);
	fReal cofPhi = timeStep / latRadius;
	fReal cofTheta = timeStep / radius;

	fReal deltaPhi = guPhi * cofPhi;
	fReal deltaTheta = guTheta * cofTheta;

	fReal midPhi = gPhi - 0.5 * deltaPhi;
	fReal midTheta = gTheta - 0.5 * deltaTheta;
	validatePhiTheta(midPhi, midTheta);

	fReal muPhi = uPhi->sampleAt(midPhi, midTheta);
	fReal muTheta = uTheta->sampleAt(midPhi, midTheta);

	deltaPhi = muPhi * cofPhi;
	deltaTheta = muTheta * cofTheta;

	fReal pPhi = gPhi - deltaPhi;
	fReal pTheta = gTheta - deltaTheta;
	validatePhiTheta(pPhi, pTheta);

	fReal advectedVal = attr->sampleAt(pPhi, pTheta);
	attr->writeValueTo(gridPhi, gridTheta, advectedVal);
}

void KaminoSolver::advectionScalar()
{
	for (auto quantity : this->centeredAttr)
	{
		KaminoQuantity* cenAttr = quantity.second;
		for (size_t gridTheta = 0; gridTheta < cenAttr->getNTheta(); ++gridTheta)
		{
			for (size_t gridPhi = 0; gridPhi < cenAttr->getNPhi(); ++gridPhi)
			{
				advectAttrAt(cenAttr, gridPhi, gridTheta);
			}
		}
	}
}

void KaminoSolver::advectionSpeed()
{
	//Advect as is for uPhi
	KaminoQuantity* uPhi = (*this)["u"];
	for (size_t gridTheta = 0; gridTheta < uPhi->getNTheta(); ++gridTheta)
	{
		for (size_t gridPhi = 0; gridPhi < uPhi->getNPhi(); ++gridPhi)
		{
			advectAttrAt(uPhi, gridPhi, gridTheta);
		}
	}

	//Trend carefully for uTheta...
	KaminoQuantity* uTheta = (*this)["v"];
	// Apart from the poles...
	for (size_t gridTheta = 1; gridTheta < uTheta->getNTheta() - 1; ++gridTheta)
	{
		for (size_t gridPhi = 0; gridPhi < uPhi->getNPhi(); ++gridPhi)
		{
			advectAttrAt(uTheta, gridPhi, gridTheta);
		}
	}
	/// TODO
	// First we derive velocity at the poles...
	// At north pole...
	// At south pole...
}

void KaminoSolver::geometric()
{
	KaminoQuantity* u = staggeredAttr["u"];
	KaminoQuantity* v = staggeredAttr["v"];

	// Poles unshifted
	for (size_t phiI = 0; phiI < nPhi; ++phiI)
	{
		size_t northPole = 0;
		size_t southPole = v->getNTheta() - 1;

		v->writeValueTo(phiI, northPole, v->getValueAt(phiI, northPole));
		v->writeValueTo(phiI, southPole, v->getValueAt(phiI, southPole));
	}
	/// TODO: Determine the upper and lower bounds
	for (size_t thetaJ = 1; thetaJ < nTheta - 1; ++thetaJ)
	{
		for (size_t phiI = 0; phiI < nPhi; ++phiI)
		{
			fReal thetaAtJ = thetaJ * gridLen;
			fReal uPrev = u->getValueAt(phiI, thetaJ);
			fReal vPrev = v->getValueAt(phiI, thetaJ);

			fReal G = 0.0;
			fReal uNext = 0.0;
			if (std::abs(thetaAtJ - M_PI / 2.0) < 1e-8)
			{
				G = 0.0;
				uNext = uPrev;
			}
			else
			{
				G = timeStep * std::cos(thetaJ * gridLen) / (radius * sin(thetaJ * gridLen));
				fReal cof = G * G;
				fReal A = 0.0;
				fReal B = (G * vPrev + 1.0) / cof;
				fReal C = -uPrev / cof;

				fReal solution[3];
				SolveP3(solution, A, B, C);

				uNext = solution[0];
			}
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

	KaminoQuantity* u = staggeredAttr["u"];
	KaminoQuantity* v = staggeredAttr["v"];
	KaminoQuantity* p = centeredAttr["p"];

	// south pole / j = 0

	for (size_t i = 0; i < nPhi; ++i)
	{
		// oot : one over two

		fReal uPlus, uMinus, vPlus, vMinus;
		// right
		size_t ipoot = (i + 1) % nPhi;
		if (getGridTypeAt(ipoot, 0) == FLUIDGRID)
		{
			uPlus = u->getValueAt(ipoot, 0);
		}
		else
		{
			uPlus = 0.0;
		}
		// left
		size_t imoot = i;
		if (getGridTypeAt(imoot, 0) == FLUIDGRID)
		{
			uMinus = u->getValueAt(imoot, 0);
		}
		else
		{
			uMinus = 0.0;
		}
		// top
		size_t jpoot = 1;
		if (getGridTypeAt(i, jpoot) == FLUIDGRID)
		{
			vPlus = v->getValueAt(i, jpoot);
		}
		else
		{
			vPlus = 0.0;
		}
		vMinus = 0.0;
		b(getIndex(i, 0)) = (uPlus - uMinus + vPlus - vMinus);
	}

	// interior of sphere grid

	for (size_t j = 0; j < nTheta; ++j)
	{
		for (size_t i = 0; i < nPhi; ++i)
		{
			// oot : one over two

			fReal uPlus, uMinus, vPlus, vMinus;
			// right
			size_t ipoot = (i + 1) % nPhi;
			if (getGridTypeAt(ipoot, j) == FLUIDGRID)
			{
				uPlus = u->getValueAt(ipoot, j);
			}
			else
			{
				uPlus = 0.0;
			}
			// left
			size_t imoot = i;
			if (getGridTypeAt(imoot, j) == FLUIDGRID)
			{
				uMinus = u->getValueAt(imoot, j);
			}
			else
			{
				uMinus = 0.0;
			}
			// top
			size_t jpoot = j + 1;
			if (getGridTypeAt(i, jpoot) == FLUIDGRID)
			{
				vPlus = v->getValueAt(i, jpoot);
			}
			else
			{
				vPlus = 0.0;
			}
			// bottom
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

	// north pole / j = nTheta - 1

	for (size_t i = 0; i < nPhi; ++i)
	{
		// oot : one over two

		fReal uPlus, uMinus, vPlus, vMinus;
		// right
		size_t ipoot = (i + 1) % nPhi;
		if (getGridTypeAt(ipoot, nTheta - 1) == FLUIDGRID)
		{
			uPlus = u->getValueAt(ipoot, nTheta - 1);
		}
		else
		{
			uPlus = 0.0;
		}
		// left
		size_t imoot = i;
		if (getGridTypeAt(imoot, nTheta - 1) == FLUIDGRID)
		{
			uMinus = u->getValueAt(imoot, nTheta - 1);
		}
		else
		{
			uMinus = 0.0;
		}
		// bottom
		size_t jmoot = nTheta - 1;
		if (getGridTypeAt(i, jmoot) == FLUIDGRID)
		{
			vMinus = v->getValueAt(i, jmoot);
		}
		else
		{
			vMinus = 0.0;
		}
		vPlus = 0.0;
		b(getIndex(i, nTheta - 1)) = (uPlus - uMinus + vPlus - vMinus);
	}

	b = b * rhsScaleB;

	Eigen::VectorXd pVector(nPhi * nTheta);
	
	Eigen::ConjugateGradient<Eigen::SparseMatrix<fReal>, Eigen::Lower | Eigen::Upper> cg;
	//cg.setTolerance(pow(10, -1));
	cg.compute(Laplacian);
	pVector = cg.solve(b);

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

	// south pole / j = 0

	fReal averageSouthV = 0.0;

	for(size_t i = 0; i < nPhi; ++i){
		fReal uBeforeUpdate = u->getValueAt(i, 0);
		fReal vBeforeUpdate = v->getValueAt(i, 0);
		fReal invSin = 1 / sin(gridLen / 2.0);
		size_t iRhs = i;
		size_t iLhs = (i == 0 ? nPhi - 1 : i - 1);
		if (getGridTypeAt(iRhs, 0) == FLUIDGRID && getGridTypeAt(iLhs, 0) == FLUIDGRID)
		{
			fReal pressureSummedU = p->getValueAt(iRhs, 0) - p->getValueAt(iLhs, 0);
			fReal deltaU = scaleP * invSin * pressureSummedU;
			u->writeValueTo(i, 0, uBeforeUpdate + deltaU);
		}
		else
		{
			u->writeValueTo(i, 0, usolid);
		}
		size_t iOpposite = (i + nPhi / 2) % nPhi;
		if (getGridTypeAt(i, 0) == FLUIDGRID && getGridTypeAt(iOpposite, 0) == FLUIDGRID)
		{
			fReal pressureSummedV = p->getValueAt(i, 0) - p->getValueAt(iOpposite, 0);
			fReal deltaV = scaleP * pressureSummedV;
			averageSouthV += uBeforeUpdate + deltaV;
			//v->writeValueTo(i, 0, vBeforeUpdate + deltaV);
		}
		else
		{
			averageSouthV += vsolid;
			//v->writeValueTo(i, 0, vsolid);
		}
	}

	// interior of sphere grid

	for(size_t j = 1; j < nTheta; ++j){
		for(size_t i = 0; i < nPhi; ++i){
			fReal uBeforeUpdate = u->getValueAt(i, j);
			fReal vBeforeUpdate = v->getValueAt(i, j);
			fReal invSin = 1 / sin(j*gridLen + gridLen / 2.0);
			size_t iRhs = i;
			size_t iLhs = (i == 0 ? nPhi - 1 : i - 1);
			if (getGridTypeAt(iRhs, j) == FLUIDGRID && getGridTypeAt(iLhs, j) == FLUIDGRID)
			{
				fReal pressureSummedU = p->getValueAt(iRhs, j) - p->getValueAt(iLhs, j);
				fReal deltaU = scaleP * invSin * pressureSummedU;
				u->writeValueTo(i, j, uBeforeUpdate + deltaU);
			}
			else
			{
				u->writeValueTo(i, j, usolid);
			}
			size_t jUpper = j;
			size_t jLower = j - 1;
			if (getGridTypeAt(i, jUpper) == FLUIDGRID && getGridTypeAt(i, jLower) == FLUIDGRID)
			{
				fReal pressureSummedV = p->getValueAt(i, jUpper) - p->getValueAt(i, jLower);
				fReal deltaV = scaleP * pressureSummedV;
				v->writeValueTo(i, j, vBeforeUpdate + deltaV);
			}
			else
			{
				v->writeValueTo(i, j, vsolid);
			}
		}
	}

	// north pole / j = nTheta

	fReal averageNorthV = 0.0;

	for(size_t i = 0; i < nPhi; ++i){
		fReal vBeforeUpdate = v->getValueAt(i, nTheta);
		fReal invSin = 1 / sin(2*M_PI - gridLen / 2.0);
		size_t iOpposite = (i + nPhi / 2) % nPhi;
		if (getGridTypeAt(i, nTheta - 1) == FLUIDGRID && getGridTypeAt(iOpposite, nTheta - 1) == FLUIDGRID)
		{
			fReal pressureSummedV = p->getValueAt(i, nTheta - 1) - p->getValueAt(iOpposite, nTheta - 1);
			fReal deltaV = scaleP * pressureSummedV;
			//v->writeValueTo(i, nTheta, vBeforeUpdate + deltaV);
			averageNorthV += vBeforeUpdate + deltaV;
		}
		else
		{
			//v->writeValueTo(i, nTheta, vsolid);
			averageNorthV += vsolid;
		}
	}

	/* average polar v velocities */

	// south pole

	averageSouthV /= nPhi;

	for(size_t i = 0; i < nPhi; ++i){
		v->writeValueTo(i, 0, averageSouthV);
	}

	// north pole

	averageNorthV /= nPhi;

	for(size_t i = 0; i < nPhi; ++i){
		v->writeValueTo(i, nTheta, averageNorthV);
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

	// south pole / j = 0

	for(size_t i = 0; i < nPhi; ++i){
		size_t numPhiNeighbors = 0;
		size_t numThetaNeighbors = 0;
		size_t rowNumber = getIndex(i, 0);
		fReal invSin = 1 / sin(gridLen / 2.0);
		// right of cell
		size_t ip1 = (i + 1) % nPhi;
		if(getGridTypeAt(ip1, 0) == FLUIDGRID){
			Laplacian.coeffRef(rowNumber, getIndex(ip1, 0)) = -1 * invSin;
			numPhiNeighbors++;
		}
		// left of cell
		size_t im1 = (i == 0 ? nPhi - 1 : i - 1);
		if(getGridTypeAt(im1, 0) == FLUIDGRID){
			Laplacian.coeffRef(rowNumber, getIndex(im1, 0)) = -1 * invSin;
			numPhiNeighbors++;
		}
		// above cell
		size_t jp1 = 1;	
		if (getGridTypeAt(i, jp1) == FLUIDGRID){
			Laplacian.coeffRef(rowNumber, getIndex(i, jp1)) = -1;
			numThetaNeighbors++;
		}
		// below cell
		// size_t iOpposite = (i + nPhi / 2) % nPhi;
		// size_t jm1 = 0;
		// if (getGridTypeAt(iOpposite, jm1) == FLUIDGRID){
		// 	Laplacian.coeffRef(rowNumber, getIndex(iOpposite, jm1)) = -1;
		// 	numThetaNeighbors++;
		// }
		Laplacian.coeffRef(rowNumber, getIndex(i, 0)) = numPhiNeighbors * invSin + numThetaNeighbors;
	}

	// interior of sphere grid

	for(size_t j = 1; j < nTheta - 1; ++j){
		for(size_t i = 0; i < nPhi; ++i){
			size_t numPhiNeighbors = 0;
			size_t numThetaNeighbors = 0;
			size_t rowNumber = getIndex(i, j);
			fReal invSin = 1 / sin(j*gridLen + gridLen / 2.0);
			// right of cell
			size_t ip1 = (i + 1) % nPhi;
			if(getGridTypeAt(ip1, j) == FLUIDGRID){
				Laplacian.coeffRef(rowNumber, getIndex(ip1, j)) = -1 * invSin;
				numPhiNeighbors++;
			}
			// left of cell
			size_t im1 = (i == 0 ? nPhi - 1 : i - 1);
			if(getGridTypeAt(im1, j) == FLUIDGRID){
				Laplacian.coeffRef(rowNumber, getIndex(im1, j)) = -1 * invSin;
				numPhiNeighbors++;
			}	
			// above cell
			size_t jp1 = j + 1;
			if(getGridTypeAt(i, jp1) == FLUIDGRID){
				Laplacian.coeffRef(rowNumber, getIndex(i, jp1)) = -1;
				numPhiNeighbors++;
			}	
			// below cell
			size_t jm1 = j - 1;
			if(getGridTypeAt(i, jm1) == FLUIDGRID){
				Laplacian.coeffRef(rowNumber, getIndex(i, jm1)) = -1;
				numPhiNeighbors++;
			}
			Laplacian.coeffRef(rowNumber, getIndex(i, 0)) = numPhiNeighbors * invSin + numThetaNeighbors;	
		}
	}

	// north pole / j = nTheta - 1

	for(size_t i = 0; i < nPhi; ++i){
		size_t numPhiNeighbors = 0;
		size_t numThetaNeighbors = 0;
		size_t rowNumber = getIndex(i, nTheta - 1);
		fReal invSin = 1 / sin(2*M_PI - gridLen / 2.0);
		// right of cell
		size_t ip1 = (i + 1) % nPhi;
		if(getGridTypeAt(ip1, nTheta - 1) == FLUIDGRID){
			Laplacian.coeffRef(rowNumber, getIndex(ip1, nTheta - 1)) = -1 * invSin;
			numPhiNeighbors++;
		}
		// left of cell
		size_t im1 = (i == 0 ? nPhi - 1 : i - 1);
		if(getGridTypeAt(im1, nTheta - 1) == FLUIDGRID){
			Laplacian.coeffRef(rowNumber, getIndex(im1, nTheta - 1)) = -1 * invSin;
			numPhiNeighbors++;
		}
		// above cell
		// size_t iOpposite = (i + nPhi / 2) % nPhi;
		// size_t jp1 = nTheta - 1;	
		// if (getGridTypeAt(iOpposite, jp1) == FLUIDGRID){
		// 	Laplacian.coeffRef(rowNumber, getIndex(iOpposite, jp1)) = -1;
		// 	numThetaNeighbors++;
		// }
		// below cell
		size_t jm1 = nTheta - 2;
		if (getGridTypeAt(i, jm1) == FLUIDGRID){
			Laplacian.coeffRef(rowNumber, getIndex(i, jm1)) = -1;
			numThetaNeighbors++;
		}
		Laplacian.coeffRef(rowNumber, getIndex(i, 0)) = numPhiNeighbors * invSin + numThetaNeighbors;
	}	
}

void KaminoSolver::initialize_pressure()
{
	for(size_t i = 0; i < nPhi; ++i){
		for(size_t j = 0; j < nTheta; ++j){
			centeredAttr["p"]->setValueAt(i, j, 0.0);
		}
	}
}

void KaminoSolver::initialize_velocity()
{
	fReal val = 0.0;
	KaminoQuantity* u = this->staggeredAttr["u"];
	KaminoQuantity* v = this->staggeredAttr["v"];

	size_t sizePhi = u->getNPhi();
	size_t sizeTheta = u->getNTheta();
	for (size_t j = 0; j < sizeTheta; ++j) {
		for (size_t i = 0; i < sizePhi; ++i) {
			val = FBM(sin(i * gridLen), sin(j * gridLen));
			u->setValueAt(i, j, val);
		}
	}
	sizePhi = v->getNPhi();
	sizeTheta = v->getNTheta();
	for (size_t j = 0; j < sizeTheta; ++j) {
		for (size_t i = 0; i < sizePhi; ++i) {
			val = FBM(cos(i * gridLen), cos(j * gridLen));
			v->setValueAt(i, j, 0.0);
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
			centeredAttr["test"]->setValueAt(i, j, 0.0);
		}
	}

	KaminoQuantity* test = centeredAttr["test"];
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

	KaminoQuantity* u = staggeredAttr["u"];
	KaminoQuantity* v = staggeredAttr["v"];
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

			pressure = centeredAttr["p"]->getValueAt(i, j);
			testVal = centeredAttr["test"]->getValueAt(i, j);
			
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


void KaminoSolver::addCenteredAttr(std::string name, fReal xOffset, fReal yOffset)
{
	size_t attrnPhi = this->nPhi;
	size_t attrnTheta = this->nTheta;
	
	KaminoQuantity* ptr = new KaminoQuantity(name, attrnPhi, attrnTheta, this->gridLen, xOffset, yOffset);
	this->centeredAttr.emplace(std::pair<std::string, KaminoQuantity*>(name, ptr));
}

void KaminoSolver::addStaggeredAttr(std::string name, fReal xOffset, fReal yOffset)
{
	size_t attrnPhi = this->nPhi;
	size_t attrnTheta = this->nTheta;
	// Is the staggered attribute u?
	if (xOffset == 0.5)
	{
		attrnTheta += 1;
	}
	KaminoQuantity* ptr = new KaminoQuantity(name, attrnPhi, attrnTheta, this->gridLen, xOffset, yOffset);
	this->staggeredAttr.emplace(std::pair<std::string, KaminoQuantity*>(name, ptr));
}

KaminoQuantity* KaminoSolver::getAttributeNamed(std::string name)
{
	return (*this)[name];
}

KaminoQuantity* KaminoSolver::operator[](std::string name)
{
	if (centeredAttr.find(name) == centeredAttr.end())
	{
		return staggeredAttr.at(name);
	}
	else
	{
		return centeredAttr.at(name);
	}
}

void KaminoSolver::swapAttrBuffers()
{
	for (auto quantity : this->centeredAttr)
	{
		quantity.second->swapBuffer();
	}
	for (auto quantity : this->staggeredAttr)
	{
		quantity.second->swapBuffer();
	}
}


// <<<<<<<<<<