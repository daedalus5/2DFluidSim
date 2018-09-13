# pragma once

# include "GlobalIncludes.h"
# include "HH16Quantity.h"

class HH16Sphere
{
private:
    size_t nTheta;              // number of grid cells in u direction // can be odd
    size_t nPhi;                // number of grid cells in v direction // MUST be even for MKL
    fReal gridLen;              // grid spacing (square in uv plane)
    fReal radius;               // radius of sphere
    
    /* practical condition: dt <= 5*dx / u_max */
    /* dt should be less than DT as well */
    float dt;                   // time step
    float DT;                   // frame rate @24 fps = 0.0147
    int frames;                 // number of frames to export

    std::string gridPath;       // folder destination for grid velocity bgeo files
    
public:
    HH16Sphere(fReal radius = 5.0, size_t nTheta = 256,
        float dt = 0.005, float DT = 1.0 / 24.0, int frames = 200,
        std::string gridPath = "HH16grid/frame");
    ~HH16Sphere(); 

    /* run the solver */
    void run();
};

