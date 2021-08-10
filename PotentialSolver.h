#ifndef _SOLVER_H
#define _SOLVER_H

#include "World.h"
#include <cmath>
#include <fftw3.h>
#include "Field.h"

class PotentialSolver
{
public:
	/*constructor, sets world*/
	PotentialSolver(World &world): 
		world(world) {}
	
	/*solves potential using Gauss-Seidel*/
	virtual bool solve() = 0;
	
	/*computes electric field = -gradient(phi)*/
	void computeEF();

protected:
	World &world;
};

class GaussSeidelSolver: public PotentialSolver
{
public:
	/*constructor, sets world*/
	GaussSeidelSolver(World &world, int max_it, double tol): 
		PotentialSolver(world), max_solver_it(max_it), tolerance(tol) {}
	
	/*solves potential using Gauss-Seidel*/
	bool solve();
	
protected:
	unsigned max_solver_it;	//maximum number of solver iterations
	double tolerance;		//solver tolerance
};

class FourierSolver: public PotentialSolver
{
public:
	/*constructor, sets world*/
	FourierSolver(World &world);

	~FourierSolver();
	
	/* solves potential using 3D FFT */
	bool solve();

private:
    /// Grid size
	int Nx, Ny, Nz;

    /// Declare FFTW components.
    //fftw_complex *mem;
    //fftw_complex *out;
    //double *in;
	
	std::vector<double> in1;
	std::vector<double> in2;
	std::vector<double> out1;
	std::vector<double> out2;
	
	fftw_plan fwrd;
    fftw_plan bwrd;
	Field phiF;			//potential

};

#endif
