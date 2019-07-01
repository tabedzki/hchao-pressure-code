#include "globals.h"
void charge_grid( void ) ;
void bonds( void ) ;


void forces() {

  int i,j, m, gind, t1, t2 ;

  charge_grid() ;


  ///////////////////////////////////////////////
  // Reset the particle forces and grid grad w //
  ///////////////////////////////////////////////
#pragma omp parallel for
  for ( i=0 ; i<nstot ; i++ )
    for ( j=0 ; j<Dim ; j++ )
      f[i][j] = 0.0 ;
 
#pragma omp parallel for
  for ( i=0 ; i<M ; i++ )
    for ( j=0 ; j<Dim ; j++ ) 
      gradwA[j][i] = 0.0 ;



  //////////////////////////////////////////////////
  // Accumulate the monomer-monomer contributions //
  //////////////////////////////////////////////////


  // Compressibility contribution //
#pragma omp parallel for
  for ( i=0 ; i<M ; i++ )
    tmp[i] = rhoha[i] ;

  fftw_fwd( tmp , ktmp ) ;

  for ( j=0 ; j<Dim ; j++ ) {
#pragma omp parallel for
    for ( i=0 ; i<M ; i++ )
      ktmp2[i] = grad_uG_hat[j][i] * ktmp[i] ;

    fftw_back( ktmp2 , tmp ) ;

#pragma omp parallel for
    for ( i=0 ; i<M ; i++ ) 
      gradwA[j][i] += tmp[i] * kappa ;
  }





  /////////////////////////////////////////////////////
  // Accumulate the nonbonded forces on each segment //
  /////////////////////////////////////////////////////
#pragma omp parallel for
  for ( i=0 ; i<nstot ; i++ ) {
    for ( m=0 ; m < grid_per_partic ; m++ ) {
      
      gind = grid_inds[ i ][ m ] ;

      for ( j=0 ; j<Dim ; j++ ) {
        if ( tp[i] == 0 )
          f[i][j] -= gradwA[ j ][ gind ] * grid_W[i][m] ;
        else if ( tp[i] == 1 )
          f[i][j] -= gradwB[ j ][ gind ] * grid_W[i][m] ;
        else if ( tp[i] == 2 )
          f[i][j] -= gradwP[ j ][ gind ] * grid_W[i][m] ;
      }
    }

    for ( j=0 ; j<Dim ; j++ )
      f[i][j] *= gvol ;
  }


  ////////////////////////////
  // Call the bonded forces //
  ////////////////////////////
  bonds() ;

}
