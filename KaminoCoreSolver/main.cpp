#include "include/KaminoQuantity.h"

const size_t nTheta = 100;               // number of grid cells in u direction
const size_t nPhi = 2 * nTheta;         // number of grid cells in v direction
const fReal gridLen = M_PI / nTheta;    // grid spacing (square in uv plane)
const fReal radius = 5.0;               // radius of sphere
/* practical condition: dt <= 5*dx / u_max */
/* dt should be less than DT as well */

const float dt = 0.005;                 // simulation time step size
const float DT = 1.0 / 24.0;            // framerate @ 24 fps = 0.0147
const int frames = 50;                  // number of frames to output
const std::string filepath = "output/frame";
const std::string tracerPath = "tracer/trace";

int main(int argc, char** argv)
{
    KaminoSolver solver(nPhi, nTheta, radius, gridLen, dt);
	solver.write_data_bgeo(filepath, 0);
    float T = 0.0;              // simulation time
    for(int i = 1; i <= frames; i++){
        while(T < i*DT){
            solver.stepForward(dt);
            T += dt;
        }
        solver.stepForward(dt + i*DT - T);
        T = i*DT;
        solver.write_data_bgeo(filepath, i);
        solver.write_data_tracer(tracerPath, i);
    }
    return 0;
}

