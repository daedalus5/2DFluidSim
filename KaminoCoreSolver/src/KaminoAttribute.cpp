# include "KaminoSolver.h"
# include <iostream>
# include <cmath>

KaminoAttribute::KaminoAttribute(std::string attributeName, size_t nx, size_t ny, fReal gridLen)
	: nx(nx), ny(ny), gridLen(gridLen), attrName(attributeName)
{}

void KaminoAttribute::swapBuffer()
{
	fReal* tempPtr = this->thisStep;
	this->thisStep = this->nextStep;
	this->nextStep = tempPtr;
}

KaminoAttribute::~KaminoAttribute()
{
	delete[] this->thisStep;
	delete[] this->nextStep;

}

fReal KaminoAttribute::getValueAt(size_t x, size_t y)
{
	return this->accessValueAt(x, y);
}

void KaminoAttribute::setValueAt(size_t x, size_t y, fReal val)
{
	this->accessValueAt(x, y) = val;
}
// TODO
fReal KaminoAttribute::sampleAt(fReal x, fReal y)
{

}

int main(int argc, char** argv)
{
	std::cout << "This worked" << std::endl;
	return 0;
}