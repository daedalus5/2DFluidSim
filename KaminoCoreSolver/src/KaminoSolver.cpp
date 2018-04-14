# include "../include/KaminoQuantity.h"
# include "../include/CubicSolver.h"

// CONSTRUCTOR / DESTRUCTOR >>>>>>>>>>

KaminoSolver::KaminoSolver(size_t nPhi, size_t nTheta, fReal radius, fReal gridLength, fReal frameDuration) :
	nPhi(nPhi), nTheta(nTheta), radius(radius), gridLen(gridLength), invGridLen(1.0 / gridLength), frameDuration(frameDuration),
	timeStep(0.0), timeElapsed(0.0), trc(M_PI / 2.0, M_PI / 2.0, radius)
{
	this->beffourierF = new fReal[nPhi * nTheta];
	this->fourieredFReal = new fReal[nPhi * nTheta];
	this->fourieredFImag = new fReal[nPhi * nTheta];
	this->fourierUReal = new fReal[nPhi * nTheta];
	this->fourierUImag = new fReal[nPhi * nTheta];

	this->a = new fReal[nTheta];
	this->b = new fReal[nTheta];
	this->c = new fReal[nTheta];
	this->dReal = new fReal[nTheta];
	this->dImag = new fReal[nTheta];

	// Our new staggered grid...
	addStaggeredAttr("u", -0.5, 0.5);		// u velocity
	addStaggeredAttr("v", 0.0, 0.0);		// v velocity
	addCenteredAttr("p", 0.0, 0.5);				// p pressure
	addCenteredAttr("test", 0.0, 0.5);			// test scalar field

	this->gridTypes = new gridType[nPhi * nTheta];
	memset(reinterpret_cast<void*>(this->gridTypes), FLUIDGRID, nPhi * nTheta);

	initialize_velocity();
	initialize_pressure();

	//precomputeLaplacian();
	initialize_test();
	initialize_boundary();
}

KaminoSolver::~KaminoSolver()
{
	delete[] beffourierF;
	delete[] fourieredFReal;
	delete[] fourieredFImag;
	delete[] fourierUReal;
	delete[] fourierUImag;

	delete[] a;
	delete[] b;
	delete[] c;
	delete[] dReal;
	delete[] dImag;

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

void KaminoSolver::updateTracer()
{
	fReal uPhi = (*this)["u"]->sampleAt(trc.phi, trc.theta);
	fReal uTheta = (*this)["v"]->sampleAt(trc.phi, trc.theta);
	trc.tracerStepForward(uPhi, uTheta, timeStep);
}

void KaminoSolver::stepForward(fReal timeStep)
{
	this->timeStep = timeStep;
	advectionScalar();
	advectionSpeed();
	std::cout << "Advection completed" << std::endl;
	this->swapAttrBuffers();

	geometric();
	std::cout << "Geometric completed" << std::endl;
	//bodyForce();
	projection();
	std::cout << "Projection completed" << std::endl;
	updateTracer();
}

// Phi: 0 - 2pi  Theta: 0 - pi
bool validatePhiTheta(fReal & phi, fReal & theta)
{
	/*int loops = static_cast<int>(std::floor(theta / M_2PI));
	theta = theta - loops * M_2PI;
	if (theta > M_PI)
	{
		theta = M_2PI - theta;
		phi += M_PI;
	}
	loops = static_cast<int>(std::floor(phi / M_2PI));
	phi = phi - loops * M_2PI;*/
	bool isFlipped = false;
	if (theta < 0.0)
	{
		theta = -theta;
		phi += M_PI;
		isFlipped = true;
	}
	if (theta > M_PI)
	{
		theta = M_2PI - theta;
		phi += M_PI;
		isFlipped = true;
	}
	if (phi > M_2PI)
		phi -= M_2PI;
	if (phi < 0.0)
		phi += M_2PI;
	return isFlipped;
}

void KaminoSolver::advectAttrAt(KaminoQuantity* attr, size_t gridPhi, size_t gridTheta)
{
	KaminoQuantity* uPhi = (*this)["u"];
	KaminoQuantity* uTheta = (*this)["v"];

	fReal gPhi = attr->getPhiCoordAtIndex(gridPhi);
	fReal gTheta = attr->getThetaCoordAtIndex(gridTheta);

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
	bool isFlipped = validatePhiTheta(pPhi, pTheta);

	fReal advectedVal = attr->sampleAt(pPhi, pTheta);
	if (isFlipped && attr->getName() == "v")
		advectedVal = -advectedVal;
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

enum Coord {x, y};

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

	//Tread carefully for uTheta...
	KaminoQuantity* uTheta = (*this)["v"];
	// Apart from the poles...
	for (size_t gridTheta = 1; gridTheta < uTheta->getNTheta() - 1; ++gridTheta)
	{
		for (size_t gridPhi = 0; gridPhi < uTheta->getNPhi(); ++gridPhi)
		{
			advectAttrAt(uTheta, gridPhi, gridTheta);
		}
	}
	/// TODO
	// First we derive velocity at the poles...
	size_t northernBelt = 0;
	size_t southernBelt = uPhi->getNTheta() - 1; // uTheta->getNTheta() - 2
	size_t northernPinch = 0;
	size_t southernPinch = uTheta->getNTheta() - 1;
	/*for (size_t gridPhi = 0; gridPhi < uPhi->getNPhi(); ++gridPhi)
	{
		size_t phiLower = gridPhi;
		size_t phiHigher = (gridPhi + 1) % uPhi->getNPhi();
		// duPhi / dPhi = -uTheta
		fReal uPhiDiff = uPhi->getValueAt(phiLower, northernBelt) - uPhi->getValueAt(phiHigher, northernBelt);
		fReal uThetaNorth = uPhiDiff * this->invGridLen;
		uTheta->writeValueTo(gridPhi, northernPinch, uThetaNorth);

		// duPhi / dPhi = uTheta
		uPhiDiff = uPhi->getValueAt(phiHigher, southernBelt) - uPhi->getValueAt(phiHigher, southernBelt);
		fReal uThetaSouth = uPhiDiff * this->invGridLen;
		uTheta->writeValueTo(gridPhi, southernPinch, uThetaSouth);
	}*/
	resetPoleVelocities();
	for (size_t gridPhi = 0; gridPhi < this->nPhi; ++gridPhi)
	{
		fReal gPhi = uPhi->getPhiCoordAtIndex(gridPhi);

		fReal uPhiVal = uPhi->getValueAt(gridPhi, northernBelt);
		uPhiNorthP[x] += -uPhiVal * std::sin(gPhi);
		uPhiNorthP[y] += uPhiVal * std::cos(gPhi);
		
		uPhiVal = uPhi->getValueAt(gridPhi, southernBelt);
		uPhiSouthP[x] += -uPhiVal * std::sin(gPhi);
		uPhiSouthP[y] += uPhiVal * std::cos(gPhi);
	}
	averageVelocities();
	fReal phiOfuPhiN = std::atan2(uPhiNorthP[y], uPhiNorthP[x]);
	if (phiOfuPhiN < 0.0)
	{
		phiOfuPhiN = M_2PI + phiOfuPhiN;
	}
	size_t northSplit = uTheta->getPhiIndexAtCoord(phiOfuPhiN);

	fReal phiOfuPhiS = std::atan2(uPhiSouthP[y], uPhiSouthP[x]);
	if (phiOfuPhiS < 0.0)
	{
		phiOfuPhiS = M_2PI + phiOfuPhiS;
	}
	size_t southSplit = uTheta->getPhiIndexAtCoord(phiOfuPhiS);
	
	fReal uAmplituteN = std::sqrt(uPhiNorthP[x] * uPhiNorthP[x] + uPhiNorthP[y] * uPhiNorthP[y]);
	fReal uAmplituteS = std::sqrt(uPhiSouthP[x] * uPhiSouthP[x] + uPhiSouthP[y] * uPhiSouthP[y]);

	size_t beltHalved = this->nPhi / 2;
	for (size_t i = 0; i < beltHalved; ++i)
	{
		size_t indexPhi = (i + northSplit) % this->nPhi;
		// North: neg
		uTheta->writeValueTo(indexPhi, northernPinch, -uAmplituteN);

		indexPhi = (i + southSplit) % this->nPhi;
		// South: pos
		uTheta->writeValueTo(indexPhi, southernPinch, uAmplituteS);
	}
	for (size_t i = beltHalved; i < this->nPhi; ++i)
	{
		size_t indexPhi = (i + northSplit) % this->nPhi;
		// North: pos
		uTheta->writeValueTo(indexPhi, northernPinch, uAmplituteN);

		indexPhi = (i + southSplit) % this->nPhi;
		// South: neg
		uTheta->writeValueTo(indexPhi, southernPinch, -uAmplituteS);
	}
}

void KaminoSolver::resetPoleVelocities()
{
	for (unsigned i = 0; i < 2; ++i)
	{
		uPhiNorthP[i] = 0.0;
		uPhiSouthP[i] = 0.0;
	}
}

void KaminoSolver::averageVelocities()
{
	for (unsigned i = 0; i < 2; ++i)
	{
		uPhiNorthP[i] /= this->nPhi;
		uPhiSouthP[i] /= this->nPhi;
	}
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

void KaminoSolver::bodyForce()
{
	fReal gravity = 9.8;
	KaminoQuantity* v = staggeredAttr["v"];

	for(size_t j = 0; j < nTheta + 1; ++j){
		for(size_t i = 0; i < nPhi; ++i){
			fReal vBeforeUpdate = v->getValueAt(i, j);
			fReal theta = j*gridLen;
			v->writeValueTo(i, j, vBeforeUpdate + gravity * sin(theta) * timeStep);
		}
	}

	v->swapBuffer();
}

void KaminoSolver::fillDivergence()
{
	KaminoQuantity* u = staggeredAttr["u"];
	KaminoQuantity* v = staggeredAttr["v"];

	fReal scaleDiv = density * radius / timeStep;
	/// TODO: Fill the fourierF buffer with divergence
	for (size_t j = 0; j < nTheta; ++j)
	{
		fReal thetaOftheBelt = (j + 0.5) * gridLen;
		fReal sine = std::sin(thetaOftheBelt);
		fReal invSine = 1.0 / sine;

		for (size_t i = 0; i < nPhi; ++i)
		{
			if (getGridTypeAt(i, j) != FLUIDGRID)
			{
				beffourierF[getIndex(i, j)] = 0.0;
				continue;//Leave it as 0 = 0 trivial problem
			}

			size_t rowNumber = getIndex(i, j);

			size_t imoot = i;
			size_t ipoot = (i + 1) % nPhi;
			size_t jmoot = j;
			size_t jpoot = j + 1;

			size_t grid2tRight = (i + 1) % nPhi;
			size_t grid2tLeft = i == 0 ? nPhi - 1 : i - 1;

			fReal uLeft = uSolid;
			fReal uRight = uSolid;
			fReal vUnder = vSolid;
			fReal vAbove = vSolid;

			if (getGridTypeAt(grid2tLeft, j) == FLUIDGRID)
			{
				uLeft = u->getValueAt(imoot, j);
			}
			if (getGridTypeAt(grid2tRight, j) == FLUIDGRID)
			{
				uRight = u->getValueAt(ipoot, j);
			}
			if (j != 0)
			{
				size_t gridUnder = j - 1;
				if (getGridTypeAt(i, gridUnder) == FLUIDGRID)
				{
					vUnder = v->getValueAt(i, jmoot);
				}
			}
			if (j != nTheta - 1)
			{
				size_t gridAbove = j + 1;
				if (getGridTypeAt(i, gridAbove) == FLUIDGRID)
				{
					vAbove = v->getValueAt(i, jpoot);
				}
			}
			fReal sinUpper = std::sin(thetaOftheBelt + 0.5 * gridLen);
			fReal sinLower = std::sin(thetaOftheBelt - 0.5 * gridLen);
			fReal termTheta = invSine * invGridLen * (vAbove * sinUpper - vUnder * sinLower);
			fReal termPhi = invSine * invGridLen * (uRight - uLeft);

			fReal div = termTheta + termPhi;
			//Additional divergence scaling goes here
			div *= scaleDiv;
			beffourierF[getIndex(i, j)] = div;
		}
	}
}

void KaminoSolver::transformDivergence()
{
	for (size_t thetaI = 0; thetaI < nTheta; ++thetaI)
	{
		for (int nIndex = 0; nIndex < nPhi; ++nIndex)
		{
			int n = nIndex - nPhi / 2;
			fReal accumulatedReal = 0.0;
			fReal accumulatedImag = 0.0;
			for (size_t j = 0; j < nPhi; ++j)
			{
				fReal phiJ = (M_2PI / nPhi) * j;
				fReal phase = -n * phiJ;
				fReal fThetaN = beffourierF[getIndex(j, thetaI)];
				accumulatedReal += fThetaN * std::cos(phase);
				accumulatedImag += fThetaN * std::sin(phase);
			}
			accumulatedReal = accumulatedReal / nPhi;
			accumulatedImag = accumulatedReal / nPhi;
			fourieredFReal[getIndex(nIndex, thetaI)] = accumulatedReal;
			fourieredFImag[getIndex(nIndex, thetaI)] = accumulatedImag;
		}
	}
}

void KaminoSolver::invTransformPressure()
{
	KaminoQuantity* p = (*this)["p"];
	for (size_t gTheta = 0; gTheta < nTheta; ++gTheta)
	{
		for (size_t gPhi = 0; gPhi < nPhi; ++gPhi)
		{
			fReal Phi = M_2PI / nPhi * gPhi;
			fReal accumulatedPressure = 0.0;
			for (int nIndex = 0; nIndex < nPhi; ++nIndex)
			{
				int n = nIndex - nPhi / 2;
				fReal phase = n * Phi;
				fReal pressureFourierCoefReal = fourierUReal[getIndex(nIndex, gTheta)];
				fReal pressureFourierCoefImag = fourierUImag[getIndex(nIndex, gTheta)];
				accumulatedPressure += pressureFourierCoefReal * std::cos(phase);
				accumulatedPressure -= pressureFourierCoefImag * std::sin(phase);
			}
			p->writeValueTo(gPhi, gTheta, accumulatedPressure);
		}
	}
}

void KaminoSolver::projection()
{
	KaminoQuantity* u = staggeredAttr["u"];
	KaminoQuantity* v = staggeredAttr["v"];
	KaminoQuantity* p = centeredAttr["p"];

	/// TODO: Fill these divergence
	fillDivergence();
	/// TODO: Perform forward FFT on fourierF to make them fourier coefficients
	transformDivergence();

	//fReal scaleD = density * radius * gridLen * gridLen / timeStep;
	fReal scaleD = gridLen * gridLen;
	for (int nIndex = 0; nIndex < nPhi; ++nIndex)
	{
		int n = nIndex - nPhi / 2;
		fReal nSqgridSq = n * gridLen;
		nSqgridSq = nSqgridSq * nSqgridSq;

		for (int i = 0; i < nTheta; ++i)
		{
			fReal thetaI = (i + 0.5) * gridLen;
			fReal sine = std::sin(thetaI);
			fReal sinSq = sine * sine;
			fReal sincos = std::cos(thetaI) * sine;
			fReal ip1im1Term2 = 0.5 * sincos * gridLen;
			scaleD *= sinSq;

			b[i] = -2.0 * sinSq - nSqgridSq;
			a[i] = sinSq - ip1im1Term2;
			c[i] = sinSq + ip1im1Term2;

			if (i == 0)
			{
				fReal coef = std::pow(-1.0, n);
				b[i] += coef * a[i];
				a[i] = 0.0;
			}
			if (i == nTheta - 1)
			{
				fReal coef = std::pow(-1.0, n);
				b[i] += coef * c[i];
				c[i] = 0.0;
			}
			fReal fTabled = this->fourieredFReal[getIndex(nIndex, i)];
			dReal[i] = fTabled * scaleD;
			fTabled = this->fourieredFImag[getIndex(nIndex, i)];
			dImag[i] = fTabled * scaleD;
		}
		
		//When n == 0, d = 0, whole system degenerates to Ax = 0 where A is singular
		if (n != 0)
		{
			TDMSolve(this->a, this->b, this->c, this->dReal);
			TDMSolve(this->a, this->b, this->c, this->dImag);
		}
		//d now contains Ui
		for (size_t UiIndex = 0; UiIndex < nTheta; ++UiIndex)
		{
			this->fourierUReal[getIndex(nIndex, UiIndex)] = dReal[UiIndex];
			this->fourierUImag[getIndex(nIndex, UiIndex)] = dImag[UiIndex];
		}
	}

	invTransformPressure();
	p->swapBuffer();

	// Update velocities accordingly: uPhi
	fReal factorTheta = -invGridLen * timeStep / (density * radius);
	for (size_t j = 0; j < u->getNTheta(); ++j)
	{
		for (size_t i = 0; i < u->getNPhi(); ++i)
		{
			fReal uBefore = u->getValueAt(i, j);
			fReal thetaBelt = (j + 0.5) * gridLen;
			fReal invSine = 1.0 / std::sin(thetaBelt);
			fReal factorPhi = factorTheta * invSine;

			size_t gridLeftI = (i == 0 ? u->getNPhi() - 1 : i - 1);
			size_t gridRightI = i;

			if (getGridTypeAt(gridLeftI, j) == SOLIDGRID ||
				getGridTypeAt(gridRightI, j) == SOLIDGRID)
			{
				u->writeValueTo(i, j, uSolid);
			}
			else
			{
				fReal pressurePhi = p->getValueAt(gridRightI, j) - p->getValueAt(gridLeftI, j);
				fReal deltauPhi = factorPhi * pressurePhi;
				u->writeValueTo(i, j, uBefore + deltauPhi);
			}
		}
	}
	// Update velocities accordingly: uTheta
	for (size_t j = 0; j < v->getNTheta(); ++j)
	{
		for (size_t i = 0; i < v->getNPhi(); ++i)
		{
			fReal vBefore = v->getValueAt(i, j);
			if (j == 0)
			{
				size_t thetaLeft = i;
				size_t thetaRight = (i + 1) % u->getNPhi();
				// At north pole : duTheta/dTheta = -uPhi
				fReal diff = u->getValueAt(thetaLeft, j) - u->getValueAt(thetaRight, j);
				diff *= invGridLen;
				v->writeValueTo(i, j, diff);
			}
			else if (j == v->getNTheta() - 1)
			{
				size_t thetaLeft = i;
				size_t thetaRight = (i + 1) % u->getNPhi();
				// At south pole : duTheta/dTheta = uPhi
				fReal diff = u->getValueAt(thetaRight, j - 1) - u->getValueAt(thetaLeft, j - 1);
				diff *= invGridLen;
				v->writeValueTo(i, j, diff);
			}
			else
			{
				size_t gridAboveJ = j;
				size_t gridBelowJ = j - 1;

				if (getGridTypeAt(i, gridBelowJ == SOLIDGRID) ||
					getGridTypeAt(i, gridAboveJ) == SOLIDGRID)
				{
					v->writeValueTo(i, j, vSolid);
				}
				else
				{
					fReal pressureTheta = p->getValueAt(i, gridAboveJ) - p->getValueAt(i, gridBelowJ);
					fReal deltauTheta = factorTheta * pressureTheta;
					v->writeValueTo(i, j, deltauTheta + vBefore);
				}
			}
		}
	}

	u->swapBuffer();
	v->swapBuffer();
}


// <<<<<<<<<<
// INITIALIZATION >>>>>>>>>>

/* Tri-diagonal matrix solver */
void KaminoSolver::TDMSolve(fReal* a, fReal* b, fReal* c, fReal* d)
{
	// |b0 c0 0 ||x0| |d0|
 	// |a1 b1 c1||x1|=|d1|
 	// |0  a2 b2||x2| |d2|

    int n = nTheta;
    n--; // since we index from 0
    c[0] /= b[0];
    d[0] /= b[0];

	for (int i = 1; i < n; i++) {
        c[i] /= b[i] - a[i]*c[i-1];
        d[i] = (d[i] - a[i]*d[i-1]) / (b[i] - a[i]*c[i-1]);
    }

    d[n] = (d[n] - a[n]*d[n-1]) / (b[n] - a[n]*c[n-1]);

    for (int i = n; i-- > 0;) {
        d[i] -= c[i]*d[i+1];
    }
}

/* Load diagonal element arrays */
void KaminoSolver::loadABC(size_t n)
{
	// A and B can be precomputed later to optimize
	std::vector<fReal> A;
	std::vector<fReal> B;
	for(size_t i = 0; i < nTheta; ++i){
		fReal theta = i * gridLen + gridLen / 2.0;
		fReal Aval = (gridLen / 2.0) * (cos(theta) / sin(theta));
		fReal Bval = gridLen * gridLen * n * n / (sin(theta) * sin(theta));
		A.push_back(Aval);
		B.push_back(Bval);
	}

	a[0] = 0.0;
	b[0] = pow(-1.0, n) - 2 - A[0] * pow(-1.0, n) - B[0];
	c[0] = 1 + B[0];
	for(size_t i = 1; i < nTheta - 1; ++i){
		a[i] = 1 - A[i];
		b[i] = -2 - B[i];
		c[i] = 1 + B[i];
	}
	a[nTheta - 1] = -2 - B[nTheta - 1];
	b[nTheta - 1] = pow(-1.0, n) - 2 + A[nTheta - 1] * pow(-1.0 , n) - B[nTheta - 1];
	c[nTheta - 1] = 0;
}

/* Duplicate of getIndex() in KaminoQuantity */
size_t KaminoSolver::getIndex(size_t x, size_t y)
{
	return y * nPhi + x;
}

gridType KaminoSolver::getGridTypeAt(size_t x, size_t y)
{
	return gridTypes[getIndex(x, y)];
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
			//val = FBM(sin(i * gridLen), sin(j * gridLen));
			val = sinSum(i*gridLen, j*gridLen);
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

fReal KaminoSolver::sinSum(const fReal x, const fReal y)
{
	fReal arg = x / (2*M_PI);
	fReal sum = sin(arg) + sin(2*arg) + sin(5*arg) + sin(3*y) + sin(7*y);
	return sum / 4.0;
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
		this->gridTypes[getIndex(gridX, nTheta / 2)] = SOLIDGRID;
		//this->gridTypes[getIndex(gridX, this->nTheta - 1)] = SOLIDGRID;
	}
}

// <<<<<<<<<<
// OUTPUT >>>>>>>>>>


void KaminoSolver::write_data_bgeo(const std::string& s, const int frame)
{
# ifndef _MSC_VER
	std::string file = s + std::to_string(frame) + ".bgeo";
	std::cout << "Writing to: " << file << std::endl;

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
			ts[0] = testVal / 13.0 * 255.0;

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

void KaminoSolver::write_data_tracer(const std::string& s, const int frame)
{
# ifndef _MSC_VER
	std::string file = s + std::to_string(frame) + ".bgeo";

	Partio::ParticlesDataMutable* parts = Partio::create();	
	Partio::ParticleAttribute pH, tracer;
	pH = parts->addAttribute("position", Partio::VECTOR, 3);
	tracer = parts->addAttribute("tracer", Partio::VECTOR, 3);
	Eigen::Matrix<fReal, 3, 1> tracerPos;

	int idx = parts->addParticle();
	trc.getCartesianXYZ(tracerPos[0], tracerPos[1], tracerPos[2]);
	float *p = parts->dataWrite<float>(pH, idx);
	float *tr = parts->dataWrite<float>(tracer, idx);

	for(int k = 0; k < 3; ++k){
		p[k] = 0.0;
		tr[k] = float (tracerPos(k, 0));
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
	//float radius = 5.0;
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
	// Is the staggered attribute uTheta?
	if (name == "v")
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