# include "../include/HH16Sphere.h"
# include "../include/KaminoTimer.h"

HH16Sphere::HH16Sphere(fReal radius, size_t nTheta,
        float dt, float DT, int frames,
        std::string gridPath) :
        radius(radius), nTheta(nTheta + 1), nPhi(2 * nTheta), gridLen(M_PI / nTheta),
        dt(dt), DT(DT), frames(frames),
        gridPath(gridPath)
{
}

HH16Sphere::~HH16Sphere()
{
}

void HH16Sphere::run()
{
    HH16Solver solver(nPhi, nTheta, radius, gridLen, dt);

    HH16Quantity* u = solver.getAttributeNamed("u");
    HH16Quantity* v = solver.getAttributeNamed("v");

    solver.write_data_bgeo(gridPath, 0);

    float T = 0.0;              // simulation time
    KaminoTimer timer;
    timer.startTimer();
    for(int i = 1; i <= frames; i++){
        while(T < i*DT){
            solver.stepForward(dt);
            T += dt;
        }
        solver.stepForward(dt + i*DT - T);
        T = i*DT;

        solver.write_data_bgeo(gridPath, i);
    }
    float cpu_time = timer.stopTimer();
    std::cout << "Time spent: " << cpu_time << " seconds" << std::endl;
    std::cout << "Performance: " << frames / cpu_time << " frames per second" << std::endl;

    std::cout << "check" << std::endl;
}