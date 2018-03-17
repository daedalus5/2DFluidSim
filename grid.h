#pragma once

#include <iostream>
#include <vector>
#include <string>
#include "Partio.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <math.h>

class Grid
{
public:
    const int nx, ny;
    const float dx, dy;
    float** u;
    float** v;
    Eigen::Matrix<float, 2, 1>** positions;
    Eigen::Matrix<float, 2, 1>** velocities;
    //std::vector<std::vector<int>> u;
    //std::vector<std::vector<int>> v;
    //std::vector<std::vector<Eigen::Matrix<float, 2, 1>>> positions;
    //std::vector<std::vector<Eigen::Matrix<float, 2, 1>>> velocities;

    Grid(int nx, int ny, float dx, float dy);
    ~Grid();

    void write_data_bgeo(const std::string& s, const int frame);
private:
    void build_particle_grid();
    void distribute_velocity();
    Eigen::Matrix<float, 2, 1> FBM() const;
};