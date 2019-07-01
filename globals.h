#include <complex>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include "omp.h"
#include "fftw3.h"
using namespace std ;

#define PI   3.141592653589793238462643383

#define pow2(x) ((x)*(x))
#define pow3(x) ((x)*(x)*(x))
#define pow4(x) ((x)*(x)*(x)*(x))
#define min(A,B) ((A)<(B) ? (A) : (B) )
#define Kdelta(i,j) ((i==j) ? 1 : 0 )

#define Dim 3


#ifndef MAIN
extern
#endif
double Range, Range2,verlet_a,verlet_b,**gn_bac,**x, **x_bac,  **f , *tmp, **grid_W, V, L[Dim], dx[Dim], gvol, Lh[Dim], 
       **rho , **w ,  *r_dudr , *tmp2, 
       *Diff, delt , chiAB, kappa, new_dt,new_C ,C, rho0 , num_averages , 
       Ubond , Uchi, Ukappa, Pvir , phiHA, phiHB,
       *rhoha, *rhohb, *rhoda, *rhodb , *rhot, *rhop, *smrhop,
       **rhoha_t, **rhohb_t, **rhoda_t, **rhodb_t , **rhop_t,
       CG_ratio, *Stress_bonds_t[Dim][Dim], Stress_bonds[Dim][Dim], Stress_nb[Dim][Dim],
       ***sts_buf , Rg, Rg3, phiP, Rp, Xi, Vp, *gammaP,
       *gradwA[Dim], *gradwB[Dim], *gradwP[Dim], 
       *uG, *grad_uG[Dim], *uP, *grad_uP[Dim], *uPG, *grad_uPG[Dim] , mem_use,
       U_chi_gg, U_chi_pg, U_chi_pp, U_kappa_gg, U_kappa_pg, U_kappa_pp ,
       *vir_func[Dim][Dim] , Ptens[Dim][Dim], Pscalar ;

#ifndef MAIN
extern
#endif
int nstot, *tp, nA, nB, nD, Nha, Nhb, Nda, Ndb, nP, 
    Nx[Dim], M, nsteps, step, print_freq , nsD, nsA, nsB, 
    **grid_inds, pmeorder, spline_weights, lagrange_weights, 
    grid_per_partic, ntypes , stress_freq , buff_size, buff_ind ,
    sample_wait, sample_freq, A_partics , pre_equil_steps ,
    nthreads ;


#ifndef MAIN
extern
#endif
char **xc, tt[80] ;


#ifndef MAIN
extern
#endif
complex<double> *ktmp2, *ktmp, I, **avg_sk , *grad_uG_hat[Dim] , *grad_uP_hat[Dim], 
  *grad_uPG_hat[Dim], *vir_func_hat[Dim][Dim] ;

#ifndef MAIN
extern
#endif
long idum;



#ifndef MAIN
extern
#endif
fftw_complex *fin, *fout ;

#ifndef MAIN
extern
#endif
fftw_plan ft_fwd, ft_bck;







double ran2(void ) ;
double gasdev2( void ) ;
int cell_stack( int* );
void cell_unstack( int , int* );
void field_gradient( double* , double* , int ) ;
void field_gradient_cdif( double* , double* , int ) ;
void convolve_fields( double*, double*, double* ) ;
void write_grid_data( const char* , double* ) ;
void write_kspace_data( const char* , complex<double>* ) ;

int unstack_stack( int ) ;
void unstack_local( int, int* ) ;
int stack( int* ) ;
int stack_local( int* ) ;
double integrate( double* );
void unstack( int , int* );
double get_k( int , double* ) ;
double get_k_alias( int , double* ) ;
void get_r( int , double* ) ;

void fftw_fwd( double* , complex<double>* );
void fftw_back( complex<double>* , double* );
int remove_dupes( int* , int );

void pbc_vdr( double*, double* , double* );
double pbc_mdr2( double*, double* , double* );
void die(const char *kill);
//Quick routine to kill the program
