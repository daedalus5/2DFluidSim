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

}