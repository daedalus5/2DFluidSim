# include "../include/KaminoSolver.h"

KaminoCenteredAttr::KaminoCenteredAttr(std::string attributeName, size_t nx, size_t ny, fReal gridLen)
	: KaminoAttribute(attributeName, nx, ny, gridLen)
{
	thisStep = new fReal[nx * ny];
	nextStep = new fReal[nx * ny];
}

KaminoCenteredAttr::~KaminoCenteredAttr(){}

fReal& KaminoCenteredAttr::accessValueAt(size_t x, size_t y)
{
# ifdef DEBUGBUILD
	// Handle exception
# endif
	return this->thisStep[x * nx + y];
}

/*
	Bilinear interpolated for now.
*/
fReal KaminoCenteredAttr::sampleAtGC(fReal x, fReal y)
{
	size_t lowerX = std::floor(x);
	size_t lowerY = std::floor(y);
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