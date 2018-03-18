# include "../include/KaminoQuantity.h"

int main(int argc, char** argv)
{
	KaminoSolver solver = KaminoSolver(4, 4, 0.1);
	solver.stepForward(0.001);
}