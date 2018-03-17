#include "include/KaminoSolver.h"

const size_t nx = 100;
const size_t ny = 100;
const fReal gridLen = 0.1;
const float dt = 1.0;
const int frames = 10;
const std::string filepath = "output/frame";

int main(int argc, char** argv)
{
    KaminoGrid grid(nx, ny, gridLen, dt);
    for(int i = 1; i <= frames; i++){
        grid.write_data_bgeo(filepath, i);
    }
    return 0;
}