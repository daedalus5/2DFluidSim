# include "../include/KaminoSolver.h"

KaminoVAttr::KaminoVAttr(std::string attributeName, size_t nx, size_t ny, fReal gridLen)
	: KaminoAttribute(attributeName, nx, ny, gridLen)
{
	thisStep = new fReal[nx * (ny + 1)];
	nextStep = new fReal[nx * (ny + 1)];
}

KaminoVAttr::~KaminoVAttr(){}

fReal& KaminoVAttr::accessValueAt(size_t x, size_t y)
{
# ifdef DEBUGBUILD
	// Handle exception
# endif
	return this->thisStep[x * nx + y];
}

/*
	Interpolated with Hermite spline.
*/
fReal KaminoVAttr::sampleAtGC(fReal x, fReal y)
{

}