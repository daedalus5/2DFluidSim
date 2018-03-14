#pragma once

#include <iostream>
#include <vector>
#include <string>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <math.h>

class Grid{
public:
    std::vector<std::vector<int>> p;
    std::vector<std::vector<int>> u;
    std::vector<std::vector<int>> v;

    Grid(int nx, int ny);
    ~Grid();

private:
}