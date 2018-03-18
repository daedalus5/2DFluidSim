# include "../include/KaminoQuantity.h"

KaminoQuantity::KaminoQuantity(std::string attributeName, size_t nx, size_t ny, fReal gridLen)
	: nx(nx), ny(ny), gridLen(gridLen), attrName(attributeName)
{
	thisStep = new fReal[nx * ny];
	nextStep = new fReal[nx * ny];
}

KaminoQuantity::~KaminoQuantity()
{
	delete[] this->thisStep;
	delete[] this->nextStep;
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
	// Handle exception
# endif
	return x * nx + y;
}

/*
Bilinear interpolated for now.
*/
fReal KaminoQuantity::sampleAt(fReal x, fReal y)
{
	size_t lowerX = getWarpedXIndex(x);
	size_t lowerY = getWarpedYIndex(y);
	size_t upperX = (lowerX + 1) % nx;
	size_t upperY = (lowerY + 1) % ny;

	fReal lowerLeft = getValueAt(lowerX, lowerY);
	fReal upperLeft = getValueAt(lowerX, upperY);
	fReal lowerRight = getValueAt(upperX, lowerY);
	fReal upperRight = getValueAt(upperX, upperY);

	fReal alphaX = x - lowerX;
	fReal alphaY = y - lowerY;

	fReal lerpedLower = KaminoLerp<fReal>(lowerLeft, lowerRight, alphaX);
	fReal lerpedUpper = KaminoLerp<fReal>(upperLeft, upperRight, alphaX);
	fReal lerped = KaminoLerp<fReal>(lerpedLower, lerpedUpper, alphaY);

	return lerped;
}

size_t KaminoQuantity::getWarpedXIndex(fReal x)
{
	int loops = static_cast<int>(std::floor(x / static_cast<fReal>(this->nx)));
	int flooredX = static_cast<int>(std::floor(x));
	int warpedX = flooredX - loops * static_cast<int>(nx);

	return static_cast<size_t>(warpedX);
}

size_t KaminoQuantity::getWarpedYIndex(fReal y)
{
	int loops = static_cast<int>(std::floor(y / static_cast<fReal>(this->ny)));
	int flooredY = static_cast<int>(std::floor(y));
	int warpedY = flooredY - loops * static_cast<int>(ny);

	return static_cast<size_t>(warpedY);
}