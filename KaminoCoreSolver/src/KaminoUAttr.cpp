# include "../include/KaminoSolver.h"

KaminoUAttr::KaminoUAttr(std::string attributeName, size_t nx, size_t ny, fReal gridLen)
	: KaminoAttribute(attributeName, nx, ny, gridLen)
{
	thisStep = new fReal[(nx + 1) * ny];
	nextStep = new fReal[(nx + 1) * ny];
}

KaminoUAttr::~KaminoUAttr(){}

fReal& KaminoUAttr::accessValueAt(size_t x, size_t y)
{
# ifdef DEBUGBUILD
	// Handle exception
# endif
	return this->thisStep[x * (nx + 1) + y];
}

fReal KaminoUAttr::sampleAtGC(fReal x, fReal y)
{
	fReal xOffset = 0.5;
	fReal yOffset = 0.5;
	x = x + xOffset;
	y = y + yOffset;

	size_t lowerX = std::floor(x);
	size_t lowerY = std::floor(y);
	size_t upperX = (lowerX + 1) % (nx + 1);
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