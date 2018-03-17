# include "KaminoSolver.h"

KaminoFacedAttr::KaminoFacedAttr(std::string attributeName, size_t nx, size_t ny, fReal gridLen)
	: KaminoAttribute(attributeName, nx, ny, gridLen)
{
	thisStep = new fReal[(nx + 1) * (ny + 1)];
	nextStep = new fReal[(nx + 1) * (ny + 1)];
}

fReal& KaminoFacedAttr::accessValueAt(size_t x, size_t y)
{
# ifdef DEBUGBUILD
	// Handle exception
# endif
	return this->thisStep[x * nx + y];
}