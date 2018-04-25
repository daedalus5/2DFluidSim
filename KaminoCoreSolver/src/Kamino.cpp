# include "../include/Kamino.h"

Kamino::Kamino(fReal radius, size_t nTheta, fReal particleDensity,
        float dt, float DT, int frames,
        std::string testPath, std::string particlePath) :
        radius(radius), nTheta(nTheta), nPhi(2 * nTheta), gridLen(M_PI / nTheta),
        particleDensity(particleDensity),
        dt(dt), DT(DT), frames(frames),
        testPath(testPath), particlePath(particlePath)
{
    fReal A1 = -1.0; fReal B1 = 0.5; fReal C1 = 0.5; fReal D1 = -0.9; fReal E1 = 1.0;
    fReal A2 = 1.0; fReal B2 = -0.3; fReal C2 = -0.7; fReal D2 = 0.8; fReal E2 = -0.8;
    fPhiCoeff = {A1, B1, C1, D1, E1};
    gThetaCoeff = {0.0};
    lPhiCoeff = {0.0};
    mThetaCoeff = {A2, B2, C2, D2, E2};
}

Kamino::~Kamino()
{
}

void Kamino::run()
{
    KaminoSolver solver(nPhi, nTheta, radius, gridLen, dt, fPhiCoeff, mThetaCoeff);
   
    KaminoParticles particles(particleDensity, radius, &solver);
    KaminoQuantity* u = solver.getAttributeNamed("u");
    KaminoQuantity* v = solver.getAttributeNamed("v");

    particles.write_data_bgeo(particlePath, 0);

    float T = 0.0;              // simulation time
    for(int i = 1; i <= frames; i++){
        while(T < i*DT){
            solver.stepForward(dt);
            particles.updatePositions(u, v, dt);
            T += dt;
        }
        solver.stepForward(dt + i*DT - T);
        particles.updatePositions(u, v, dt);
        T = i*DT;

        particles.write_data_bgeo(particlePath, i);
    }
}

void Kamino::test()
{
    KaminoSolver solver(nPhi, nTheta, radius, gridLen, dt, fPhiCoeff, mThetaCoeff);
    solver.write_data_bgeo(testPath, 0);

    float T = 0.0;  // simulation time
    for(int i = 1; i <= 5; i++){
        while(T < i*DT){
            solver.stepForward(dt);
            T += dt;
        }
        solver.stepForward(dt + i*DT - T);
        T = i*DT;
        solver.write_data_bgeo(testPath, i);
    }
}