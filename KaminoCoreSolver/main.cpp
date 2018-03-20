#include "include/KaminoQuantity.h"

const size_t nx = 100;          // number of grid cells in u direction
const size_t ny = 100;          // number of grid cells in v direction
const fReal gridLen = 0.1;      // grid spacing (square in uv plane)
/* practical condition: dt <= 5*dx / u_max */
/* dt should be less than DT as well */
const float dt = 0.005;         // simulation time step size
const float DT = 1.0 / 24.0;    // framerate @ 24 fps = 0.0147
const int frames = 10;          // number of frames to output
const std::string filepath = "output/frame";

int main(int argc, char** argv)
{
    KaminoSolver solver(nx, ny, gridLen, dt);
    float T = 0.0;              // simulation time
    for(int i = 1; i <= frames; i++){
        solver.write_data_bgeo(filepath, i);
        while(T < i*DT){
            solver.stepForward(dt);
            T += dt;
        }
        solver.stepForward(dt + i*DT - T);
        T = i*DT;
    }
    return 0;
}

