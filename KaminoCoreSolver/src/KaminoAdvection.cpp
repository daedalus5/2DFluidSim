# include "../include/KaminoQuantity.h"


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

enum Coord { x, y };

void KaminoSolver::solvePolarVelocities()
{
	KaminoQuantity* uPhi = (*this)["u"];
	KaminoQuantity* uTheta = (*this)["v"];

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
	solvePolarVelocities();
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
