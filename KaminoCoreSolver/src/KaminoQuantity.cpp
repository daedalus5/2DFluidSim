# include "../include/KaminoQuantity.h"

KaminoQuantity::KaminoQuantity(std::string attributeName, size_t nPhi, size_t nTheta, fReal gridLen, fReal xOffset, fReal yOffset)
	: nPhi(nPhi), nTheta(nTheta), gridLen(gridLen), invGridLen(1.0 / gridLen), attrName(attributeName), xOffset(xOffset), yOffset(yOffset)
{
	thisStep = new fReal[nPhi * nTheta];
	nextStep = new fReal[nPhi * nTheta];
}

KaminoQuantity::~KaminoQuantity()
{
	delete[] thisStep;
	delete[] nextStep;
}

size_t KaminoQuantity::getNPhi()
{
	return this->nPhi;
}

size_t KaminoQuantity::getNTheta()
{
	return this->nTheta;
}

void KaminoQuantity::swapBuffer()
{
	fReal* tempPtr = this->thisStep;
	this->thisStep = this->nextStep;
	this->nextStep = tempPtr;
}

fReal KaminoQuantity::getValueAt(size_t x, size_t y)
{
	return this->accessValueAt(x, y);
}

void KaminoQuantity::setValueAt(size_t x, size_t y, fReal val)
{
	this->accessValueAt(x, y) = val;
}

fReal& KaminoQuantity::accessValueAt(size_t x, size_t y)
{
	return this->thisStep[getIndex(x, y)];
}

void KaminoQuantity::writeValueTo(size_t x, size_t y, fReal val)
{
	this->nextStep[getIndex(x, y)] = val;
}

size_t KaminoQuantity::getIndex(size_t x, size_t y)
{
# ifdef DEBUGBUILD
	if (x >= this->nPhi || y >= this->nTheta)
	{
		std::cerr << "Index out of bound at x: " << x << " y: " << y << std::endl;
	}
# endif
	return y * nPhi + x;
}

fReal KaminoQuantity::getPhiCoordAtIndex(size_t x)
{
	fReal xFloat = static_cast<fReal>(x) + xOffset;
	return xFloat * this->gridLen;
}

fReal KaminoQuantity::getThetaCoordAtIndex(size_t y)
{
	fReal yFloat = static_cast<fReal>(y) + yOffset;
	return yFloat * this->gridLen;
}

/*
Bilinear interpolated for now.
*/
fReal KaminoQuantity::sampleAt(fReal x, fReal y)
{
	int phiIndex = std::floor(x * invGridLen - this->xOffset);
	int thetaIndex = std::floor(y * invGridLen - this->yOffset);

	size_t lowerX = phiIndex < 0 ? this->nPhi - 1 : phiIndex;
	size_t upperX = lowerX + 1;
	upperX = upperX >= nPhi ? 0 : upperX;

	//This wouldn't go below 0 but incase there's floating point error...
	size_t lowerY = thetaIndex < 0 ? 0 : thetaIndex;
	size_t upperY = lowerY + 1;
	upperY = upperY >= nTheta ? nTheta - 1 : upperY;

	fReal lowerLeft = getValueAt(lowerX, lowerY);
	fReal upperLeft = getValueAt(lowerX, upperY);
	fReal lowerRight = getValueAt(upperX, lowerY);
	fReal upperRight = getValueAt(upperX, upperY);

	fReal alphaX = x - static_cast<fReal>(std::floor(x));
	fReal alphaY = y - static_cast<fReal>(std::floor(y));

	fReal lerpedLower = KaminoLerp<fReal>(lowerLeft, lowerRight, alphaX);
	fReal lerpedUpper = KaminoLerp<fReal>(upperLeft, upperRight, alphaX);
	fReal lerped = KaminoLerp<fReal>(lerpedLower, lerpedUpper, alphaY);

	return lerped;
}

fReal KaminoQuantity::getThetaOffset()
{
	return this->xOffset;
}

fReal KaminoQuantity::getPhiOffset()
{
	return this->yOffset;
}