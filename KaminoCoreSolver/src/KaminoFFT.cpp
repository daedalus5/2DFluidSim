# include "../include/KaminoQuantity.h"
# include <fstream>

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
	/*std::fstream fNaiveFourier;
	fNaiveFourier.open("./naive.txt", std::ios::out);
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
			fNaiveFourier << accumulatedReal << " + " << accumulatedImag << "i" << std::endl;
		}
	}
	fNaiveFourier.close();*/

	for (size_t thetaI = 0; thetaI < nTheta; ++thetaI)
	{
		std::vector<std::complex<fReal>> output;
		std::vector<std::complex<fReal>> input;
		for (size_t j = 0; j < nPhi; ++j)
		{
			input.push_back(std::complex<float>(beffourierF[getIndex(j, thetaI)], 0.0));
		}
		fft.inv(output, input);
		for (int naiveIndex = 0; naiveIndex < nPhi; ++naiveIndex)
		{
			int fftIndex = nPhi / 2 - naiveIndex;
			if (fftIndex < 0)
				fftIndex += nPhi;
			fourieredFReal[getIndex(naiveIndex, thetaI)] = output.at(fftIndex).real();
			fourieredFImag[getIndex(naiveIndex, thetaI)] = output.at(fftIndex).imag();
		}
	}
}

void KaminoSolver::invTransformPressure()
{
	/*std::fstream fNaive;
	fNaive.open("./naiveFWD.txt", std::ios::out);
	KaminoQuantity* p = (*this)["p"];
	for (size_t gTheta = 0; gTheta < nTheta; ++gTheta)
	{
		for (size_t gPhi = 0; gPhi < nPhi; ++gPhi)
		{
			fReal Phi = M_2PI / nPhi * gPhi;
			fReal accumulatedPressure = 0.0;
			fReal accumulatedImagPres = 0.0;
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
				accumulatedImagPres += pressureFourierCoefReal * std::sin(phase) + pressureFourierCoefImag * std::cos(phase);
			}
			fNaive << accumulatedPressure << " + " << accumulatedImagPres << "i" << std::endl;
			p->writeValueTo(gPhi, gTheta, accumulatedPressure);
		}
	}
	fNaive.close();*/

	KaminoQuantity* p = (*this)["p"];
	for (size_t gTheta = 0; gTheta < nTheta; ++gTheta)
	{
		std::vector<std::complex<fReal>> output;
		std::vector<std::complex<fReal>> input;

		for (size_t j = 0; j < nPhi; ++j)
		{
			input.push_back(std::complex<float>(fourierUReal[getIndex(j, gTheta)], fourierUImag[getIndex(j, gTheta)]));
		}
		fft.fwd(output, input);
		fReal realNZeroComponent = output.at(nPhi / 2).real();
		for (size_t gPhi = 0; gPhi < nPhi; ++gPhi)
		{
			fReal pressure = 0.0;
			size_t fftIndex = 0;
			if (gPhi != 0)
				fftIndex = nPhi - gPhi;
			if (gPhi % 2 == 0)
			{
				pressure = output.at(fftIndex).real() - realNZeroComponent;
			}
			else
			{
				pressure = -output.at(fftIndex).real() - realNZeroComponent;
			}
			p->writeValueTo(gPhi, gTheta, pressure);
		}
	}
}
