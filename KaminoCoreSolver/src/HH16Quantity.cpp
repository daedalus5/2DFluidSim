# include "../include/HH16Quantity.h"

HH16Quantity::HH16Quantity(std::string attributeName, size_t nPhi, size_t nTheta, fReal gridLen)
    : nPhi(nPhi), nTheta(nTheta), gridLen(gridLen), invGridLen(1.0 / gridLen), attrName(attributeName)
{
    thisStep = new fReal[nPhi * nTheta];
    nextStep = new fReal[nPhi * nTheta];
}

HH16Quantity::~HH16Quantity()
{
    delete[] thisStep;
    delete[] nextStep;
}

std::string HH16Quantity::getName()
{
    return this->attrName;
}

size_t HH16Quantity::getNPhi()
{
    return this->nPhi;
}

size_t HH16Quantity::getNTheta()
{
    return this->nTheta;
}

void HH16Quantity::swapBuffer()
{
    fReal* tempPtr = this->thisStep;
    this->thisStep = this->nextStep;
    this->nextStep = tempPtr;
}

fReal HH16Quantity::getValueAt(size_t x, size_t y)
{
    return this->accessValueAt(x, y);
}

void HH16Quantity::setValueAt(size_t x, size_t y, fReal val)
{
    this->accessValueAt(x, y) = val;
}

fReal& HH16Quantity::accessValueAt(size_t x, size_t y)
{
    return this->thisStep[getIndex(x, y)];
}

void HH16Quantity::writeValueTo(size_t x, size_t y, fReal val)
{
/*# ifdef DEBUGBUILD
    if (val > 1e4)
    {
        std::cerr << "Explosion detected " << std::endl;
    }
# endif*/
    this->nextStep[getIndex(x, y)] = val;
}

size_t HH16Quantity::getIndex(size_t x, size_t y)
{
# ifdef DEBUGBUILD
    if (x >= this->nPhi || y >= this->nTheta)
    {
        std::cerr << "Index out of bound at x: " << x << " y: " << y << std::endl;
    }
# endif
    return y * nPhi + x;
}

fReal HH16Quantity::getPhiCoordAtIndex(size_t x)
{
    fReal xFloat = static_cast<fReal>(x);
    return xFloat * this->gridLen;
}

size_t HH16Quantity::getPhiIndexAtCoord(fReal phi)
{
    fReal phiInt = phi * this->invGridLen;
    return static_cast<size_t>(phiInt);
}

fReal HH16Quantity::getThetaCoordAtIndex(size_t y)
{
    fReal yFloat = static_cast<fReal>(y);
    return yFloat * this->gridLen;
}

size_t HH16Quantity::getThetaIndexAtCoord(fReal theta)
{
    fReal thetaInt = theta * this->invGridLen;
    return static_cast<size_t>(thetaInt);
}

/*
Bilinear interpolated for now.
*/
fReal HH16Quantity::sampleAt(fReal x, fReal y, fReal uNorthP[2], fReal uSouthP[2])
{
    return 0;
    //TODO According to HH16

    /*
    fReal phi = x;
    fReal theta = y;
    // Phi and Theta are now shifted back to origin

    bool isFlippedPole = validatePhiTheta(phi, theta);
    fReal normedPhi = phi * invGridLen;
    fReal normedTheta = theta * invGridLen;

    int phiIndex = static_cast<int>(std::floor(normedPhi));
    int thetaIndex = static_cast<int>(std::floor(normedTheta));
    fReal alphaPhi = normedPhi - static_cast<fReal>(phiIndex);
    fReal alphaTheta = normedTheta - static_cast<fReal>(thetaIndex);

    if (thetaIndex == 0 && isFlippedPole) // If it's not flipped the theta+1 belt would just be belt 1
    {
        if (attrName == "u")
        {
            alphaTheta = 2.0 * alphaTheta;
            size_t phiLower = (phiIndex + nPhi / 2) % nPhi;
            size_t phiHigher = (phiLower + 1) % nPhi;
            fReal lowerBelt = Lerp(getValueAt(phiLower, 0), getValueAt(phiHigher, 0), alphaPhi);

            fReal lowerPhi = (phiLower - 0.5) * gridLen;
            fReal higherPhi = (phiLower + 0.5) * gridLen;
            fReal loweruPhi = -uNorthP[0] * std::cos(lowerPhi) + uNorthP[1] * std::sin(lowerPhi);
            fReal higheruPhi = -uNorthP[0] * std::cos(higherPhi) + uNorthP[1] * std::sin(higherPhi);
            fReal higherBelt = Lerp(loweruPhi, higheruPhi, alphaPhi);

            fReal lerped = Lerp(lowerBelt, higherBelt, alphaTheta);
            return lerped;
        }
        else
        {
            //Lower is to the opposite, higher is on this side
            alphaTheta = 1.0 - alphaTheta;
            size_t phiLower = phiIndex % nPhi;
            size_t phiHigher = (phiLower + 1) % nPhi;
            size_t phiLowerOppo = (phiLower + nPhi / 2) % nPhi;
            size_t phiHigherOppo = (phiHigher + nPhi / 2) % nPhi;

            fReal lowerBelt = Lerp<fReal>(getValueAt(phiLower, 0), getValueAt(phiHigher, 0), alphaPhi);
            fReal higherBelt = Lerp<fReal>(getValueAt(phiLowerOppo, 0), getValueAt(phiHigherOppo, 0), alphaPhi);
            fReal lerped = Lerp<fReal>(lowerBelt, higherBelt, alphaTheta);
            return lerped;
        }
    }
    else if (thetaIndex == nTheta - 1)
    {
        if (attrName == "u")
        {
            alphaTheta = 2.0 * alphaTheta;
            size_t phiLower = phiIndex % nPhi;
            size_t phiHigher = (phiLower + 1) % nPhi;
            fReal lowerBelt = Lerp(getValueAt(phiLower, thetaIndex), getValueAt(phiHigher, thetaIndex), alphaPhi);
            
            fReal lowerPhi = (phiLower - 0.5) * gridLen;
            fReal higherPhi = (phiLower + 0.5) * gridLen;
            fReal loweruPhi = -uSouthP[0] * std::cos(lowerPhi) + uSouthP[1] * std::sin(lowerPhi);
            fReal higheruPhi = -uSouthP[0] * std::cos(higherPhi) + uSouthP[1] * std::sin(higherPhi);
            fReal higherBelt = Lerp(loweruPhi, higheruPhi, alphaPhi);

            fReal lerped = Lerp(lowerBelt, higherBelt, alphaTheta);
            return lerped;
        }
        else
        {
            //Lower is on this side, higher is to the opposite
            size_t phiLower = phiIndex % nPhi;
            size_t phiHigher = (phiLower + 1) % nPhi;
            size_t phiLowerOppo = (phiLower + nPhi / 2) % nPhi;
            size_t phiHigherOppo = (phiHigher + nPhi / 2) % nPhi;

            fReal lowerBelt = Lerp<fReal>(getValueAt(phiLower, nTheta - 1), getValueAt(phiHigher, nTheta - 1), alphaPhi);
            fReal higherBelt = Lerp<fReal>(getValueAt(phiLowerOppo, nTheta - 1), getValueAt(phiHigherOppo, nTheta - 1), alphaTheta);

            fReal lerped = Lerp<fReal>(lowerBelt, higherBelt, alphaTheta);
            return lerped;
        }
    }
    else
    {
        size_t phiLower = phiIndex % nPhi;
        size_t phiHigher = (phiLower + 1) % nPhi;
        size_t thetaLower = thetaIndex;
        size_t thetaHigher = thetaIndex + 1;

        fReal lowerBelt = Lerp<fReal>(getValueAt(phiLower, thetaLower), getValueAt(phiHigher, thetaLower), alphaPhi);
        fReal higherBelt = Lerp<fReal>(getValueAt(phiLower, thetaHigher), getValueAt(phiHigher, thetaHigher), alphaPhi);

        fReal lerped = Lerp<fReal>(lowerBelt, higherBelt, alphaTheta);
        return lerped;
    }
    */
}