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

	fReal muPhi = uPhi->sampleAt(midPhi, midTheta);
	fReal muTheta = uTheta->sampleAt(midPhi, midTheta);

	deltaPhi = muPhi * cofPhi;
	deltaTheta = muTheta * cofTheta;

	fReal pPhi = gPhi - deltaPhi;
	fReal pTheta = gTheta - deltaTheta;

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
	for (size_t gridPhi = 0; gridPhi < nPhi; ++gridPhi)
	{
		fReal phi = (M_2PI / nPhi) * gridPhi;

		size_t phiLower = gridPhi;
		size_t phiHigher = (phiLower + 1) % nPhi;
		size_t phiLowerOppo = (phiLower + nPhi / 2) % nPhi;
		size_t phiHigherOppo = (phiHigher + nPhi / 2) % nPhi;

		fReal uBeltThisside = KaminoLerp(uPhi->getValueAt(phiLower, northernBelt), uPhi->getValueAt(phiHigher, northernBelt), 0.5);
		fReal uBeltOpposide = KaminoLerp(uPhi->getValueAt(phiLowerOppo, northernBelt), uPhi->getValueAt(phiHigherOppo, northernBelt), 0.5);
		fReal uPhiNorthernPolar = KaminoLerp(uBeltThisside, uBeltOpposide, 0.5);
		//uTheta is behind uPhi for pi/2 at north pole
		uTheta->writeValueTo(gridPhi, northernPinch, uPhiNorthernPolar);

		
		uBeltThisside = KaminoLerp(uPhi->getValueAt(phiLower, southernBelt), uPhi->getValueAt(phiHigher, southernBelt), 0.5);
		uBeltOpposide = KaminoLerp(uPhi->getValueAt(phiLowerOppo, southernBelt), uPhi->getValueAt(phiHigherOppo, southernBelt), 0.5);
		fReal uPhiSouthernPolar = KaminoLerp(uBeltThisside, uBeltOpposide, 0.5);
		//uTheta is before uPhi for pi/2 at south pole
		uTheta->writeValueTo(gridPhi, southernPinch, uPhiSouthernPolar);
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
