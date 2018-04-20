# include "../include/KaminoQuantity.h"


void KaminoSolver::fillDivergence()
{
	KaminoQuantity* u = staggeredAttr["u"];
	KaminoQuantity* v = staggeredAttr["v"];

	//fReal scaleDiv = density * radius / timeStep;
	fReal scaleDiv = 1.0;
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
			size_t grid2tLeft = (i == 0 ? nPhi - 1 : i - 1);

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
			accumulatedImag = accumulatedImag / nPhi;
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
				if (n == 0)
					pressureFourierCoefReal = 0.0;
				fReal pressureFourierCoefImag = fourierUImag[getIndex(nIndex, gTheta)];
				if (n == 0)
					pressureFourierCoefImag = 0.0;
				accumulatedPressure += pressureFourierCoefReal * std::cos(phase);
				accumulatedPressure -= pressureFourierCoefImag * std::sin(phase);
				if (abs(accumulatedPressure) > 1e2)
				{
					std::cerr << "exp" << std::endl;
				}
			}
			p->writeValueTo(gPhi, gTheta, accumulatedPressure);
		}
	}
}
