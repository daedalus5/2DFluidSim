# include "../include/HH16Quantity.h"

void HH16Solver::projection()
{
    HH16Quantity* u = attr["u"];
    HH16Quantity* v = attr["v"];
    HH16Quantity* p = attr["pressure"];

    // TODO: Run MKL routine to solve for pressure values

	// updates pressure buffer 
    MKLSphericalPoisson();

    p->swapBuffer();

	fReal invRad = 1.0 / radius;

    // Update velocities accordingly: uPhi
    for (size_t j = 1; j < u->getNTheta() - 1; ++j)
    {
        for (size_t i = 0; i < u->getNPhi(); ++i)
        {
            fReal uBefore = u->getValueAt(i, j);
            fReal thetaBelt = j * gridLen;
            fReal invSine = 1.0 / std::sin(thetaBelt);
            fReal factorPhi = -invGridLen * invRad * invSine;

            size_t gridLeftI = (i == 0 ? u->getNPhi() - 1 : i - 1);
			size_t gridRightI = (i == u->getNPhi() - 1 ? 0 : i + 1);

			fReal pressureGrad = p->getValueAt(gridRightI, j) - p->getValueAt(gridLeftI, j);
			fReal deltauPhi = factorPhi * pressureGrad;
			u->writeValueTo(i, j, uBefore + deltauPhi);
        }
    }
	
	// pressure is a scalar and so should be singularly defined at the poles
	// this means that there won't be a pressure gradient in the u direction 
	for (size_t i = 0; i < u->getNPhi(); ++i)
	{
		fReal uBefore_NP = u->getValueAt(i, 0);
		fReal uBefore_SP = u->getValueAt(i, u->getNTheta() - 1);
		u->writeValueTo(i, 0, uBefore_NP);
		u->writeValueTo(i, u->getNTheta() - 1, uBefore_SP);
	}

    u->swapBuffer();

    // Update velocities accordingly: uTheta
    for (size_t j = 1; j < v->getNTheta() - 1; ++j)
    {
        for (size_t i = 0; i < v->getNPhi(); ++i)
        {
            fReal vBefore = v->getValueAt(i, j);
			fReal factorTheta = -invGridLen * invRad;

			size_t gridAboveJ = j + 1;
            size_t gridBelowJ = j - 1;

			fReal pressureGrad = p->getValueAt(i, gridAboveJ) - p->getValueAt(i, gridBelowJ);
			fReal deltauTheta = factorTheta * pressureGrad;
			v->writeValueTo(i, j, deltauTheta + vBefore);
        }
    }

	for (size_t i = 0; i < v->getNPhi(); ++i)
	{
		// North Pole
		fReal vBefore = v->getValueAt(i, 0);
		fReal factorTheta = -invGridLen * invRad;

		size_t gridShift = (i + nPhi / 2) % nPhi;
		fReal pressureGrad = p->getValueAt(i, 1) - p->getValueAt(gridShift, 1);
		fReal deltauTheta = factorTheta * pressureGrad;
		v->writeValueTo(i, 0, deltauTheta + vBefore);

		// South Pole
		vBefore = v->getValueAt(i, v->getNTheta() - 1);
		pressureGrad = p->getValueAt(i, v->getNTheta() - 1) - p->getValueAt(gridShift, v->getNTheta() - 1);
		deltauTheta = factorTheta * pressureGrad;
		v->writeValueTo(i, v->getNTheta() - 1, vBefore + deltauTheta);
	}

	// uPhi vals at poles calculated using BCs
    v->swapBuffer();
}


/*******************************************************************************
* Copyright 2006-2018 Intel Corporation.
*
* This software and the related documents are Intel copyrighted  materials,  and
* your use of  them is  governed by the  express license  under which  they were
* provided to you (License).  Unless the License provides otherwise, you may not
* use, modify, copy, publish, distribute,  disclose or transmit this software or
* the related documents without Intel's prior written permission.
*
* This software and the related documents  are provided as  is,  with no express
* or implied  warranties,  other  than those  that are  expressly stated  in the
* License.
*******************************************************************************/

/*
*  Content:
*  C double precision example of solving Helmholtz problem on a whole sphere
*  using Intel(R) MKL Poisson Library
*
*******************************************************************************/
void HH16Solver::MKLSphericalPoisson()
{
	// MKL assumes grid defined as follows:
	// nTheta : (0, ..., n_theta)
	// nPhi   : (0, ..., n_phi)
	// where nPhi is the number of cells in the phi direction
	//		 nTheta is the number of cells in the theta direction
	// hence, mesh grid is size (nt + 1) x (np + 1)

	MKL_INT np = nPhi, nt = nTheta - 1;

	MKL_INT ip, it, i, stat;
	MKL_INT ipar[128];
	double ap, bp, at, bt, lp, lt, hp, ht, theta_i, ct, c1;
	double *dpar = NULL, *f = NULL, *u = NULL;
	double q;
	DFTI_DESCRIPTOR_HANDLE handle_s = 0;
	DFTI_DESCRIPTOR_HANDLE handle_c = 0;
	MKL_INT mem_error, error;

	error = 0;
	/* memory allocation */
	mem_error = 1;
	dpar = (double*)mkl_malloc((5 * np / 2 + nt + 10) * sizeof(double), 64);
	if (dpar == NULL) goto end;
	f = (double*)mkl_malloc((np + 1)*(nt + 1) * sizeof(double), 64);
	if (f == NULL) goto end;
	/* memory allocated correctly */
	mem_error = 0;

	/* Defining the rectangular domain on a sphere 0<p<2*pi, 0<t<pi
	for Helmholtz Solver on a sphere
	Poisson Library will automatically detect that this problem is on a whole sphere! */
	ap = 0.0E0;
	bp = 2 * M_PI;
	at = 0.0E0;
	bt = M_PI;

	/* Setting the coefficient q to 0.0E0 for Poisson problem
	If you like to solve Helmholtz problem, please set q to 1.0E0 */
	q = 0.0E0;

	/* Computing the mesh size hp in phi-direction */
	lp = bp - ap;
	hp = lp / np;
	/* Computing the mesh size ht in theta-direction */
	lt = bt - at;
	ht = lt / nt;

	/* 
	Filling in the right-hand side f(p,t)=(2+q)*cos(t)
	in the mesh points into the array f.
	We choose the right-hand side to correspond to the TRUE solution
	of Helmholtz equation on a sphere.
	Here we are using the mesh sizes hp and ht computed before to compute
	the coordinates (phi_i,theta_i) of the mesh points */
	for (it = 0; it <= nt; it++)
	{
		for (ip = 0; ip <= np; ip++)
		{
			theta_i = ht*it;
			ct = cos(theta_i);
			f[ip + it*(np + 1)] = ct*(2. + q);
		}
	}

	/* Initializing ipar array to make it free from garbage */
	for (i = 0; i<128; i++)
	{
		ipar[i] = 0;
	}

	/* Initializing simple data structures of Poisson Library
	for Poisson Solver on a sphere
	As we are looking for the solution on a whole sphere, this is a PERIODIC problem
	Therefore, the routines ending with "_p" are used to find the solution */
	d_init_sph_p(&ap, &bp, &at, &bt, &np, &nt, &q, ipar, dpar, &stat);
	if (stat != 0) {
		error = 1;
		goto end;
	}

	/* Initializing complex data structures of Poisson Library
	for Poisson Solver on a sphere
	NOTE: Right-hand side f may be altered after the Commit step. If you want
	to keep it, you should save it in another memory location! */
	d_commit_sph_p(f, &handle_s, &handle_c, ipar, dpar, &stat);
	if (stat != 0) {
		error = 1;
		goto end;
	}
	/* Computing the approximate solution of Poisson problem on a whole sphere */
	// note: solution is stored in f
	d_sph_p(f, &handle_s, &handle_c, ipar, dpar, &stat);
	if (stat != 0) {
		error = 1;
		goto end;
	}
	/* Cleaning the memory used by handle_s and handle_c */
	free_sph_p(&handle_s, &handle_c, ipar, &stat);
	if (stat != 0) {
		error = 1;
		goto end;
	}
	/* Now we can use handle_s and handle_c to solve another Poisson problem */
	/* after a proper initialization */

	/* Printing the results */
	printf("The number of mesh intervals in phi-direction is np=%d\n", np);
	printf("The number of mesh intervals in theta-direction is nt=%d\n\n", nt);

	for (it = 0; it <= nt; it++)
	{
		for (ip = 0; ip <= np; ip++)
		{
			printf("%10.3f", f[ip + it*(np + 1)]);
		}
	}

	// TODO: check divergence free here

end:
	/* Free Intel(R) MKL memory if any was allocated */
	mkl_free(dpar);
	mkl_free(f);
	mkl_free(u);
	MKL_Free_Buffers();
	/* Failure message to print if something went wrong */
	if (mem_error == 1)
	{
		printf("| insufficient memory \n");
	}
	if (error != 0)
	{
		printf("\nDouble precision Helmholtz example on a whole sphere has ");
		printf("FAILED to compute the solution...\n");
	}
	/* Success message to print if everything is OK */
	printf("\n Double precision Helmholtz example on a whole sphere has ");
	printf("successfully PASSED\n through all steps of computation!\n");

	printf("press any key to continue...");
	getchar();
}
