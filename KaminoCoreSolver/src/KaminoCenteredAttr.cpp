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