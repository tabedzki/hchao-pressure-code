#include "globals.h"

void calc_Unb() {

  int i ;
  // Take Gaussian potential to k-space
  fftw_fwd( uG , ktmp2 ) ;
  
  
  // Polymer chi A-B contribution //
  fftw_fwd( rho[0] , ktmp ) ;


  // Polymer kappa contribution //
#pragma omp parallel for
  for ( i=0 ; i<M ; i++ )
    tmp2[i] = rhot[i] ;

  fftw_fwd( tmp2 , ktmp ) ;

#pragma omp parallel for
  for ( i=0 ; i<M ; i++ )
    ktmp[i] *= ktmp2[i] ;

  fftw_back( ktmp , tmp ) ;

#pragma omp parallel for
  for ( i=0 ; i<M ; i++ )
    tmp[i] *= tmp2[i] ;

  U_kappa_gg = integrate( tmp ) * kappa / 2.0 ;

}

