# include "../include/KaminoSolver.h"

KaminoVAttr::KaminoVAttr(std::string attributeName, size_t nx, size_t ny, fReal gridLen)
	: KaminoAttribute(attributeName, nx, ny, gridLen)
{
	thisStep = new fReal[nx * (ny + 1)];
	nextStep = new fReal[nx * (ny + 1)];
}

KaminoVAttr::~KaminoVAttr(){}

size_t KaminoVAttr::getIndex(size_t x, size_t y)
{
# ifdef DEBUGBUILD
	// Handle exception
# endif
	return x * nx + y;
}

/*
	Interpolated with Hermite spline.
*/
fReal KaminoVAttr::sampleAtGC(fReal x, fReal y)
{
	fReal xOffset = 0.5;
	fReal yOffset = 0.5;
	x = x + xOffset;
	y = y + yOffset;

	size_t lowerX = getWarpedXIndex(x);
	size_t lowerY = getWarpedYIndex(y);
	size_t upperX = (lowerX + 1) % nx;
	size_t upperY = (lowerY + 1) % (ny + 1);

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


// x has no offset yet.
size_t KaminoVAttr::getWarpedXIndex(fReal x)
{
	int loops = std::floor(x / static_cast<fReal>(this->nx));
	int flooredX = std::floor(x);
	int warpedX = flooredX - loops * static_cast<int>(nx);

	return static_cast<size_t>(warpedX);
}

// y has no offset yet either.
size_t KaminoVAttr::getWarpedYIndex(fReal y)
{
	int loops = std::floor(y / static_cast<fReal>(this->ny + 1));
	int flooredY = std::floor(y);
	int warpedY = flooredY - loops * static_cast<int>(ny);

	return static_cast<size_t>(warpedY);
}