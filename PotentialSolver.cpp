#include "PotentialSolver.h"
#include "Field.h"
#include <math.h>
#include <iostream>
#include <iomanip>
#include "World.h"

#include <cmath>
#include <fftw3.h>

#define SIZE 8
#define NXX SIZE
#define NYY SIZE
#define NZZ SIZE


using namespace std;
using namespace Const;

/*solves Poisson equation using Gauss-Seidel*/
bool GaussSeidelSolver::solve()
{
    //references to avoid having to write world.phi
	Field &phi = world.phi;
    Field &rho = world.rho;

	//precompute 1/(dx^2)
    double3 dh = world.getDh();
    double idx2 = 1.0/(dh[0]*dh[0]);
    double idy2 = 1.0/(dh[1]*dh[1]);
    double idz2 = 1.0/(dh[2]*dh[2]);

    double L2=0;			//norm
    bool converged= false;

    /*solve potential*/
    for (unsigned it=0;it<max_solver_it;it++)
    {
		 for (int i=1;i<world.ni-1;i++)
            for (int j=1;j<world.nj-1;j++)
                for (int k=1;k<world.nk-1;k++)
                {
					//standard internal open node
					double phi_new = (rho[i][j][k]/Const::EPS_0 +
									idx2*(phi[i-1][j][k] + phi[i+1][j][k]) +
									idy2*(phi[i][j-1][k]+phi[i][j+1][k]) +
									idz2*(phi[i][j][k-1]+phi[i][j][k+1]))/(2*idx2+2*idy2+2*idz2);

					/*SOR*/
					phi[i][j][k] = phi[i][j][k] + 1.4*(phi_new-phi[i][j][k]);
				}

		 /*check for convergence*/
		 if (it%25==0)
		 {
			double sum = 0;
			for (int i=1;i<world.ni-1;i++)
				for (int j=1;j<world.nj-1;j++)
					for (int k=1;k<world.nk-1;k++)
					{
						double R = -phi[i][j][k]*(2*idx2+2*idy2+2*idz2) +
									rho[i][j][k]/Const::EPS_0 +
									idx2*(phi[i-1][j][k] + phi[i+1][j][k]) +
									idy2*(phi[i][j-1][k]+phi[i][j+1][k]) +
									idz2*(phi[i][j][k-1]+phi[i][j][k+1]);

						sum += R*R;
					}

			L2 = sqrt(sum/(world.ni*world.nj*world.nk));
			if (L2<tolerance) {converged=true;break;}
		}
    }

    if (!converged) cerr<<"GS failed to converge, L2="<<L2<<endl;
    return converged;
}

FourierSolver::FourierSolver(World &world): PotentialSolver(world), phiF(world.ni,world.nj,world.nk) {
	Nx = world.ni-2;     // Don't include the boundary nodes for Fourier Solver
	Ny = world.nj-2;
	Nz = world.nk-2;
	//Nzh = (Nz/2+1);
    in1.resize(Nx*Ny*Nz,0.0);
    in2.resize(Nx*Ny*Nz,0.0);
    out1.resize(Nx*Ny*Nz,0.0);
    out2.resize(Nx*Ny*Nz,0.0);
	//mem = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * Nx * Ny * Nzh);
	//out = mem;
	//in = mem[0];

	//fwrd = fftw_plan_dft_r2c_3d(Nx,Ny,Nz,in,out,FFTW_MEASURE);
	//bwrd = fftw_plan_dft_c2r_3d(Nx,Ny,Nz,out,in,FFTW_MEASURE);
    fwrd = fftw_plan_r2r_3d(Nx,Ny,Nz, in1.data(), out1.data(), FFTW_RODFT00, FFTW_RODFT00, FFTW_RODFT00, FFTW_EXHAUSTIVE);
    bwrd = fftw_plan_r2r_3d(Nx,Ny,Nz, in2.data(), out2.data(), FFTW_RODFT00, FFTW_RODFT00, FFTW_RODFT00, FFTW_EXHAUSTIVE);

}

FourierSolver::~FourierSolver() {
	//fftw_free(mem);
	fftw_destroy_plan(fwrd);
	fftw_destroy_plan(bwrd);
	fftw_cleanup();
}
	

/*solves Poisson equation using FFT*/
bool FourierSolver::solve()
{
    double Pi = M_PI;
	Field &phi = world.phi;
    Field &rho = world.rho;
	//Field tmp(world.ni,world.nj,world.nk);
   
	//cout << "Fourier Solver \n";
	/*
    /// Grid size
    int Nx=NXX,Ny=NYY,Nz=NZZ,Nzh=(Nz/2+1);

    /// Declare FFTW components.
    fftw_complex *mem;
    fftw_complex *out;
    double *in;
    mem = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * Nx * Ny * Nzh);
    out = mem;
    in = mem[0];

    fftw_plan fwrd = fftw_plan_dft_r2c_3d(Nx,Ny,Nz,in,out,FFTW_MEASURE);
    fftw_plan bwrd = fftw_plan_dft_c2r_3d(Nx,Ny,Nz,out,in,FFTW_MEASURE);
    */
	
    //double L0=0.0;
    //double L1=2*Pi;
    //double xlen = (L1-L0);
	double3 X0 = world.getX0();
	double3 Xm = world.getXm();
    double xlen = Xm[0]-X0[0];
	double ylen = Xm[1]-X0[1];
	double zlen = Xm[2]-X0[2];
	
    double dx = xlen/(double)(Nx+1);
    double dy=ylen/(double)(Ny+1);
    double dz=zlen/(double)(Nz+1);

    int l=0;
	// set some values
	for (int i=0; i<Nx; i++)
		for (int j=0;j<Ny;j++)
			for(int k=0;k<Nz;k++) {				
				// consecutive ordering
				//size_t u = k*Nx*Ny + j*Nx + i;
				//in[u] = u;
				in1[l] = rho[i+1][j+1][k+1]/Const::EPS_0;
				l=l+1;
			}

    fftw_execute(fwrd);
	
    l=-1;
    for (int i = 0; i < Nx; i++){  
        for(int j = 0; j < Ny; j++){
			for(int k = 0; k < Nz; k++){

				l=l+1;
				double fact=0;

				fact=(2-2*cos((i+1)*Pi/(Nx+1)))/(dx*dx);

				fact+= (2-2*cos((j+1)*Pi/(Ny+1)))/(dy*dy);

				fact+= (2-2*cos((k+1)*Pi/(Nz+1)))/(dz*dz);

				in2[l] = out1[l]/fact;
			}
        }
    }


    fftw_execute(bwrd);
	
	//cout << "\n\n";
    //cout<<"Executed bwrd transform " << endl;
	
    l=-1;
    //double erl1 = 0.;
    //double erl2 = 0.;
    //double erl3 = 0.;
	for (int i=0; i<Nx; i++)
		for (int j=0;j<Ny;j++)
			for(int k=0;k<Nz;k++) {				
				// consecutive ordering
				l=l+1;
				//tmp[i+1][j+1][k+1] = phi[i+1][j+1][k+1];
				//phiF[i+1][j+1][k+1] = 0.125*out2[l]/((double)(Nx+1))/((double)(Ny+1))/((double)(Nz+1));
				//erl1 +=pow(fabs(phi[i+1][j+1][k+1]-phiF[i+1][j+1][k+1]),2);
				//erl2 +=pow(fabs(phi[i+1][j+1][k+1]),2);
				//erl3 +=pow(fabs(phiF[i+1][j+1][k+1]),2);
				
				/*
				if ((l > Nx*Ny*Nz/2) && (l < Nx*Ny*Nz/2 + Nz*2))
				{
					cout<< setprecision(7) << "phi[" << i << "][" << j << "]["<<k<<"] = " << tmp[i+1][j+1][k+1] << " , phiF[" << i << "][" << j << "]["<<k<<"] = " << phiF[i+1][j+1][k+1] << "\n";
					//cout<< "phi[" << i << "][" << j << "]["<<k<<"] = " << phi[i][j][k] << " , phiF[" << i << "][" << j << "]["<<k<<"] = " << phiF[i][j][k] << "\n";
				}
				*/
				
				phi[i+1][j+1][k+1] = 0.125*out2[l]/((double)(Nx+1))/((double)(Ny+1))/((double)(Nz+1));
								
			}
    //cout<< setprecision(7) << "\n phiF error=" <<sqrt(erl1) << " , " << sqrt(erl2) << " , " << sqrt(erl3) << endl ;  
      

    return true;
	
}	

/*computes electric field = -gradient(phi) using 2nd order differencing*/
void PotentialSolver::computeEF()
{
	//reference to phi to avoid needing to write world.phi
	Field &phi = world.phi;

	double3 dh = world.getDh();
	double dx = dh[0];
	double dy = dh[1];
	double dz = dh[2];

	for (int i=0;i<world.ni;i++)
		for (int j=0;j<world.nj;j++)
			for (int k=0;k<world.nk;k++)
			{
				double3 &ef = world.ef[i][j][k]; //reference to (i,j,k) ef vec3

				/*x component*/
				if (i==0)	/*forward*/
					ef[0] = -(-3*phi[i][j][k]+4*phi[i+1][j][k]-phi[i+2][j][k])/(2*dx);	
				else if (i==world.ni-1)  /*backward*/
					ef[0] = -(phi[i-2][j][k]-4*phi[i-1][j][k]+3*phi[i][j][k])/(2*dx);	
				else  /*central*/
					ef[0] = -(phi[i+1][j][k] - phi[i-1][j][k])/(2*dx);	

				/*y component*/
				if (j==0)
					ef[1] = -(-3*phi[i][j][k] + 4*phi[i][j+1][k]-phi[i][j+2][k])/(2*dy);
				else if (j==world.nj-1)
					ef[1] = -(phi[i][j-2][k] - 4*phi[i][j-1][k] + 3*phi[i][j][k])/(2*dy);
				else
					ef[1] = -(phi[i][j+1][k] - phi[i][j-1][k])/(2*dy);

				/*z component*/
				if (k==0)
					ef[2] = -(-3*phi[i][j][k] + 4*phi[i][j][k+1]-phi[i][j][k+2])/(2*dz);
				else if (k==world.nk-1)
					ef[2] = -(phi[i][j][k-2] - 4*phi[i][j][k-1]+3*phi[i][j][k])/(2*dz);
				else
					ef[2] = -(phi[i][j][k+1] - phi[i][j][k-1])/(2*dz);
			}
}
