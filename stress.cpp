#include "globals.h"

void nonbond_stress() ;
void bond_stress() ;

void calc_stress() {


  int j1, j2 ;
  nonbond_stress() ;
  bond_stress() ;


  for ( j1=0 ; j1<Dim ; j1++ )
    for ( j2=0 ; j2<Dim ; j2++ )
      Ptens[j1][j2] = -(Stress_nb[j1][j2] + Stress_bonds[j1][j2]) + double(nA*Nha)/V * Kdelta(j1,j2) ;
      //Ptens[j1][j2] = -Stress_nb[j1][j2]
        //+ double(nA)/V * Kdelta(j1,j2) ;

  Pscalar = 0.0 ;
  for ( j1=0 ; j1<Dim ; j1++ )
    Pscalar += Rg3*Ptens[j1][j1] / double( Dim ) ;

}


void nonbond_stress() {

  int i, j1, j2 ;

  fftw_fwd( rhoha , ktmp ) ;

  for ( j1=0 ; j1<Dim ; j1++ ) {
    for ( j2=0 ; j2<Dim ; j2++ ) {
#pragma omp parallel for private(i)
      for ( i=0 ; i<M ; i++ )
        ktmp2[i] = ktmp[i] * vir_func_hat[j1][j2][i] ;

      fftw_back( ktmp2 , tmp ) ;

#pragma omp parallel for private(i)
      for ( i=0 ; i<M ; i++ )
        tmp2[i] = rhoha[i] * tmp[i] ; 

      Stress_nb[j1][j2] = -integrate( tmp2 ) * kappa /2.0 / V ;
    }
  }
}


void bond_stress() {
  int i, j1, j2, m, ind ;
  double ratio=double(Nha)/double(Nha-1), mdr2, dr[Dim] ;

  ind = 0 ;

  for (j1=0 ; j1<Dim ; j1++ ) 
    for ( j2=0 ; j2<Dim ; j2++ ) {
      Stress_bonds[j1][j2] = 0.0 ;
      for ( m=0; m<nthreads ; m++ )
        Stress_bonds_t[j1][j2][m] = 0.0 ;
    }


  // Homopolymer A bonds //
#pragma omp parallel for private(ind,m,j1,j2,mdr2,dr) 
  for ( i=0 ; i<nA ; i++ ) {
    for ( m=0 ; m<Nha - 1 ; m++ ) {

      ind = nD * (Nda + Ndb) + i * Nha + m ;

      mdr2 = pbc_mdr2( x[ind] , x[ind+1] , dr ) ;

      // int tid = //omp_get_thread_num() ;
      int tid =0; //omp_get_thread_num() ;
      for ( j1=0 ; j1<Dim ; j1++ )
        for ( j2=0 ; j2<Dim ; j2++ )
          Stress_bonds_t[j1][j2][tid] += dr[j1] * dr[j2] ;
      
    }
  
  } // for ( i=0 ; i<nT[k]

  for ( j1=0 ; j1<Dim ; j1++ )   
    for ( j2=0 ; j2<Dim ; j2++ )
      for ( m=0 ; m<nthreads ; m++ ) 
        Stress_bonds[j1][j2] += ratio*3.0 * Stress_bonds_t[j1][j2][m] / V ;

}

