# include "../include/HH16Quantity.h"


void HH16Solver::advectAttrAt(HH16Quantity* attr, size_t gridPhi, size_t gridTheta)
{
    HH16Quantity* uPhi = (*this)["u"];
    HH16Quantity* uTheta = (*this)["v"];

    fReal gPhi = attr->getPhiCoordAtIndex(gridPhi);
    fReal gTheta = attr->getThetaCoordAtIndex(gridTheta);

    fReal guPhi = uPhi->sampleAt(gPhi, gTheta, this->uNorthP, this->uSouthP);
    fReal guTheta = uTheta->sampleAt(gPhi, gTheta, this->uNorthP, this->uSouthP);

    fReal latRadius = this->radius * std::sin(gTheta);
    fReal cofPhi = timeStep / latRadius;
    fReal cofTheta = timeStep / radius;

    fReal deltaPhi = guPhi * cofPhi;
    fReal deltaTheta = guTheta * cofTheta;

    fReal midPhi = gPhi - 0.5 * deltaPhi;
    fReal midTheta = gTheta - 0.5 * deltaTheta;

    fReal muPhi = uPhi->sampleAt(midPhi, midTheta, this->uNorthP, this->uSouthP);
    fReal muTheta = uTheta->sampleAt(midPhi, midTheta, this->uNorthP, this->uSouthP);

    fReal averuPhi = 0.5 * (muPhi + guPhi);
    fReal averuTheta = 0.5 * (muTheta + guTheta);

    deltaPhi = averuPhi * cofPhi;
    deltaTheta = averuTheta * cofTheta;

    fReal pPhi = gPhi - deltaPhi;
    fReal pTheta = gTheta - deltaTheta;

    fReal advectedVal = attr->sampleAt(pPhi, pTheta, this->uNorthP, this->uSouthP);
    
    attr->writeValueTo(gridPhi, gridTheta, advectedVal);
}

void HH16Solver::advection()
{
    for (auto quantity : this->attr)
    {
        HH16Quantity* scalarAttr = quantity.second;

        for (size_t gridTheta = 1; gridTheta < scalarAttr->getNTheta() - 1; ++gridTheta)
        {
            for (size_t gridPhi = 0; gridPhi < scalarAttr->getNPhi(); ++gridPhi)
            {
                advectAttrAt(scalarAttr, gridPhi, gridTheta);
            }
        }
    }

	solvePolarVelocitiesAdvection();
	solvePolarScalarsAdvection();
}

void HH16Solver::solvePolarVelocitiesAdvection()
{
	HH16Quantity* uPhi = (*this)["u"];
	HH16Quantity* uTheta = (*this)["v"];

	// HH16 central differencing scheme

	for (size_t gridPhi = 0; gridPhi < nPhi; ++gridPhi)
	{
		size_t gridShift = (gridPhi + nPhi / 2) % nPhi;

		// North pole
		fReal u_n = uTheta->getValueAt(gridPhi, 0);
		fReal u_down = uTheta->getValueAt(gridPhi, 1);
		fReal u_up = uTheta->getValueAt(gridShift, 1);
		fReal u_star = u_n + timeStep * (u_n / radius) * ((u_up - u_down) / (2.0 * gridLen));
		this->NPBuffer[gridPhi] = u_star;

		// South pole
		u_n = uTheta->getValueAt(gridPhi, nTheta - 1);
		u_down = uTheta->getValueAt(gridPhi, nTheta - 2);
		u_up = uTheta->getValueAt(gridShift, nTheta - 2);
		u_star = u_n + timeStep * (u_n / radius) * ((u_up - u_down) / (2.0 * gridLen));
		this->SPBuffer[gridPhi] = u_star;
	}

	applyPolarBoundaryCondition();
}

void HH16Solver::solvePolarScalarsAdvection()
{
	// HH16 central differencing scheme

	HH16Quantity* d = (*this)["density"];

	// North pole
	size_t phi_NP = std::floor(atan2(uNorthP[1], uNorthP[0]));
	size_t gridShift = (phi_NP + nPhi / 2) % nPhi;
	fReal d_n = d->getValueAt(phi_NP, 0);
	fReal d_down = d->getValueAt(phi_NP, 1);
	fReal d_up = d->getValueAt(gridShift, 1);
	fReal coeff = sqrt(uNorthP[0] * uNorthP[0] + uNorthP[1] * uNorthP[1]);
	fReal d_star = d_n + timeStep * (coeff / radius) * ((d_up - d_down) / (2.0 * gridLen));

	// South pole
	size_t phi_SP = std::floor(atan2(uSouthP[1], uSouthP[0]));
	gridShift = (phi_SP + nPhi / 2) % nPhi;
	d_n = d->getValueAt(phi_SP, nTheta - 1);
	d_down = d->getValueAt(phi_SP, nTheta - 2);
	d_up = d->getValueAt(gridShift, nTheta - 2);
	coeff = sqrt(uSouthP[0] * uSouthP[0] + uSouthP[1] * uSouthP[1]);
	d_star = d_n + timeStep * (coeff / radius) * ((d_up - d_down) / (2.0 * gridLen));

	// all pole grid points get same density value
	for (size_t gridPhi = 0; gridPhi < nPhi; ++gridPhi)
	{
		d->writeValueTo(gridPhi, 0, d_star);
		d->writeValueTo(gridPhi, nTheta - 1, d_star);
	}
}

void HH16Solver::resetPoleVelocities()
{
    for (unsigned i = 0; i < 2; ++i)
    {
        uNorthP[i] = 0.0;
        uSouthP[i] = 0.0;
    }
}
