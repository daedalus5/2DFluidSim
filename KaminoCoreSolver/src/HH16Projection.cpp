# include "../include/HH16Quantity.h"

void HH16Solver::projection()
{
    HH16Quantity* u = attr["u"];
    HH16Quantity* v = attr["v"];
    HH16Quantity* p = attr["p"];

    /// TODO: Run MKL routine to solve for pressure values

    MKLSphericalPoisson();

    p->swapBuffer();

    fReal factorTheta = -invGridLen;

    // TODO: Implement according to HH16
    /*
    // Update velocities accordingly: uPhi
    for (size_t j = 0; j < u->getNTheta(); ++j)
    {
        for (size_t i = 0; i < u->getNPhi(); ++i)
        {
            fReal uBefore = u->getValueAt(i, j);
            fReal thetaBelt = (j + 0.5) * gridLen;
            fReal invSine = 1.0 / std::sin(thetaBelt);
            fReal factorPhi = factorTheta * invSine;

            size_t gridLeftI = (i == 0 ? u->getNPhi() - 1 : i - 1);
            size_t gridRightI = i;

            if (getGridTypeAt(gridLeftI, j) == SOLIDGRID ||
                getGridTypeAt(gridRightI, j) == SOLIDGRID)
            {
                u->writeValueTo(i, j, uSolid);
            }
            else
            {
                fReal pressurePhi = 0.0;
                if (getGridTypeAt(gridLeftI, j) == FLUIDGRID)
                pressurePhi -= p->getValueAt(gridLeftI, j);
                if (getGridTypeAt(gridRightI, j) == FLUIDGRID)
                    pressurePhi += p->getValueAt(gridRightI, j);
                fReal deltauPhi = factorPhi * pressurePhi;
                u->writeValueTo(i, j, uBefore + deltauPhi);
            }
        }
    }

    u->swapBuffer();

    // Update velocities accordingly: uTheta
    for (size_t j = 1; j < v->getNTheta() - 1; ++j)
    {
        for (size_t i = 0; i < v->getNPhi(); ++i)
        {
            fReal vBefore = v->getValueAt(i, j);
            size_t gridAboveJ = j;
            size_t gridBelowJ = j - 1;
            
            if (getGridTypeAt(i, gridBelowJ) == SOLIDGRID ||
                getGridTypeAt(i, gridAboveJ) == SOLIDGRID)
            {
                v->writeValueTo(i, j, vSolid);
            }
            else
            {
                fReal pressureTheta = 0.0;
                if (getGridTypeAt(i, gridBelowJ) == FLUIDGRID)
                    pressureTheta -= p->getValueAt(i, gridBelowJ);
                if (getGridTypeAt(i, gridAboveJ) == FLUIDGRID)
                    pressureTheta += p->getValueAt(i, gridAboveJ);
                fReal deltauTheta = factorTheta * pressureTheta;
                v->writeValueTo(i, j, deltauTheta + vBefore);
            }
        }
    }
    
    */

    solvePolarVelocities();
    v->swapBuffer();
}

void HH16Solver::MKLSphericalPoisson()
{

}
