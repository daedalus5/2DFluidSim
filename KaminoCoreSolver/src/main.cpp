# include "../include/KaminoSolver.h"

int main(int argc, char** argv)
{
	KaminoGrid solver = KaminoGrid(4, 4, 0.1);
	solver.stepForward(0.001);
}