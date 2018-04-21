#include "include/KaminoQuantity.h"

const size_t nTheta = 512;               // number of grid cells in u direction
const size_t nPhi = 2 * nTheta;         // number of grid cells in v direction
const fReal gridLen = M_PI / nTheta;    // grid spacing (square in uv plane)
const fReal radius = 5.0;               // radius of sphere
/* practical condition: dt <= 5*dx / u_max */
/* dt should be less than DT as well */

const float dt = 0.005;                 // simulation time step size
const float DT = 1.0 / 24.0;            // framerate @ 24 fps = 0.0147
const int frames = 1000;                  // number of frames to output
const std::string filepath = "output/frame";
const std::string tracerPath = "tracer/trace";
const std::string particlePath = "particles/frame";

int main(int argc, char** argv)
{
    KaminoSolver solver(nPhi, nTheta, radius, gridLen, dt);
# ifndef _MSC_VER
    solver.write_data_bgeo(filepath, 0);
# endif
   
    KaminoParticles particles(5000, radius, &solver);
    KaminoQuantity* u = solver.getAttributeNamed("u");
    KaminoQuantity* v = solver.getAttributeNamed("v");
# ifndef _MSC_VER
    particles.write_data_bgeo(particlePath, 0);
# endif

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
# ifndef _MSC_VER
        solver.write_data_bgeo(filepath, i);
        particles.write_data_bgeo(particlePath, i);
        //solver.write_data_tracer(tracerPath, i);
# endif
    }
    return 0;
}

