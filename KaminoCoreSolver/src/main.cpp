# include "../include/KaminoQuantity.h"

int main(int argc, char** argv)
{
	KaminoSolver solver = KaminoSolver(2, 2, 0.1);
	solver.stepForward(0.001);
}