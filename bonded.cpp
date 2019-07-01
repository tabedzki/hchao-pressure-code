#include "globals.h"


void bonds( ) {

  int i, j, k, m, ind ;
  double ratio = double(Nda+Ndb)/double(Nda + Ndb-1),mdr2, mdr, dr[Dim] ;
  
  Ubond = 0.0 ;

  // Diblock bonds //
#pragma omp parallel for \
  private(ind, i, j, m, dr, mdr2) \
  reduction(+:Ubond)
  for ( i=0 ; i<nD ; i++ ) {
    for ( m=0 ; m<Nda + Ndb - 1 ; m++ ) {

      ind = i * (Nda + Ndb) + m ;

      mdr2 = pbc_mdr2( x[ind] , x[ind+1] , dr ) ;

      Ubond += mdr2 * 1.5 ;

      for ( j=0 ; j<Dim ; j++ ) {
        f[ind][j] -= ratio*3.0 * dr[j] ;
        f[ind+1][j] += ratio*3.0 * dr[j] ;
      }

    }

  } // for ( i=0 ; i<nT[k]

ratio = double(Nha)/double(Nha-1);

  // Homopolymer A bonds //
#pragma omp parallel for \
  private(ind, j, m, dr, mdr2) \
  reduction(+:Ubond)
  for ( i=0 ; i<nA ; i++ ) {
    for ( m=0 ; m<Nha - 1 ; m++ ) {

      ind = nD * (Nda + Ndb) + i * Nha + m ;

      mdr2 = pbc_mdr2( x[ind] , x[ind+1] , dr ) ;

      Ubond += mdr2 * 1.5 ;

      for ( j=0 ; j<Dim ; j++ ) {
        f[ind][j] -= ratio*3.0 * dr[j] ;
        f[ind+1][j] += ratio*3.0 * dr[j] ;
      }

    }

  } // for ( i=0 ; i<nT[k]

  ratio  = double(Nhb)/double(Nhb-1);
  // Homopolymer B bonds //
#pragma omp parallel for \
  private(ind, i, j, m, dr, mdr2) \
  reduction(+:Ubond)
  for ( i=0 ; i<nB ; i++ ) {
    for ( m=0 ; m<Nhb - 1 ; m++ ) {

      ind = nD * (Nda + Ndb) + nA * Nha + i * Nhb + m ;

      mdr2 = pbc_mdr2( x[ind] , x[ind+1] , dr ) ;

      Ubond += mdr2 * 1.5 ;

      for ( j=0 ; j<Dim ; j++ ) {
        f[ind][j] -= ratio*3.0 * dr[j] ;
        f[ind+1][j] += ratio*3.0 * dr[j] ;
      }

    }

  } // for ( i=0 ; i<nT[k]


}
