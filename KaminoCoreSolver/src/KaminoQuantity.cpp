# include "../include/KaminoQuantity.h"

KaminoQuantity::KaminoQuantity(std::string attributeName, size_t nx, size_t ny, fReal gridLen, fReal xOffset, fReal yOffset)
	: nx(nx), ny(ny), gridLen(gridLen), invGridLen(1.0 / gridLen), attrName(attributeName), xOffset(xOffset), yOffset(yOffset)
{
	thisStep = new fReal[nx * ny];
	nextStep = new fReal[nx * ny];
}

KaminoQuantity::~KaminoQuantity()
{
	delete[] thisStep;
	delete[] nextStep;
}

size_t KaminoQuantity::getNx()
{
	return this->nx;
}

size_t KaminoQuantity::getNy()
{
	return this->ny;
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
	if (x >= this->nx || y >= this->ny)
	{
		std::cerr << "Index out of bound at x: " << x << " y: " << y << std::endl;
	}
# endif
	return y * nx + x;
}

fReal KaminoQuantity::getXCoordAtIndex(size_t x)
{
	fReal xFloat = static_cast<fReal>(x) - xOffset;
	return xFloat * this->gridLen;
}

fReal KaminoQuantity::getYCoordAtIndex(size_t y)
{
	fReal yFloat = static_cast<fReal>(y) - yOffset;
	return yFloat * this->gridLen;
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

	fReal alphaX = x - static_cast<fReal>(std::floor(x));
	fReal alphaY = y - static_cast<fReal>(std::floor(y));

	fReal lerpedLower = KaminoLerp<fReal>(lowerLeft, lowerRight, alphaX);
	fReal lerpedUpper = KaminoLerp<fReal>(upperLeft, upperRight, alphaX);
	fReal lerped = KaminoLerp<fReal>(lerpedLower, lerpedUpper, alphaY);

	return lerped;
}

size_t KaminoQuantity::getWarpedXIndex(fReal x)
{
	x = x * invGridLen;
	x += xOffset;
	int loops = static_cast<int>(std::floor(x / static_cast<fReal>(this->nx)));
	int flooredX = static_cast<int>(std::floor(x));
	int warpedX = flooredX - loops * static_cast<int>(nx);

	return static_cast<size_t>(warpedX);
}

size_t KaminoQuantity::getWarpedYIndex(fReal y)
{
	y = y * invGridLen;
	y += yOffset;
	int loops = static_cast<int>(std::floor(y / static_cast<fReal>(this->ny)));
	int flooredY = static_cast<int>(std::floor(y));
	int warpedY = flooredY - loops * static_cast<int>(ny);

	return static_cast<size_t>(warpedY);
}