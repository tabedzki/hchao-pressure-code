#include "globals.h"

void update_positions() {

  int i,j ;
  double tmp_d,tmp_x,tmp_gn,delt_2 = delt*delt;
//#pragma openmp parallel for \
//  private(i, j)
  for ( i=0 ; i<nstot ; i++ ) {
    for ( j=0 ; j<Dim ; j++ ) {
 
      tmp_x = x[i][j];
      tmp_gn = gasdev2();

      tmp_d = x[i][j] - x_bac[i][j];
 
      if(tmp_d >Lh[j])
         x_bac[i][j] += L[j];
      else if(tmp_d <-Lh[j])
      	 x_bac[i][j] -= L[j];
      else{
      
      }    

   

      x[i][j] = 2.0*verlet_b*x[i][j] - verlet_a*x_bac[i][j] + verlet_b*delt_2*f[i][j] + verlet_b*delt*sqrt( 2.0/Diff[ tp[i] ] * delt )*(tmp_gn+gn_bac[i][j])/2.0;
      //x[i][j] = x_bac[i][j] + 2.0*Diff[ tp[i] ] * delt * f[i][j]+ 2.0*sqrt( Diff[ tp[i] ] * delt ) *tmp_gn;
      
      if ( x[i][j] > L[j] )
        x[i][j] -= L[j] ;
      
      else if ( x[i][j] < 0.0 )
        x[i][j] += L[j] ;
     
     x_bac[i][j] = tmp_x;
     gn_bac[i][j] = tmp_gn;
    }

  }//for ( i=0 ; ...


}

void update_positions_init() {

  int i,j ;

   double tmp_gn,delt_2 = delt*delt;


//#pragma openmp parallel for \
//  private(i, j)
  for ( i=0 ; i<nstot ; i++ ) {
    for ( j=0 ; j<Dim ; j++ ) {
      x_bac[i][j] = x[i][j] ; 
     
      tmp_gn=gasdev2();

      gn_bac[i][j] = tmp_gn;

      x[i][j] = x[i][j] + verlet_b*delt_2*f[i][j] + verlet_b*delt*sqrt( 2.0/Diff[ tp[i] ] * delt )*tmp_gn/2.0;

     // x[i][j] = x[i][j] + Diff[ tp[i] ] * delt * f[i][j]
       //       + sqrt( 2.0 * Diff[ tp[i] ] * delt ) * gasdev2() ;
      
      if ( x[i][j] > L[j] )
        x[i][j] -= L[j] ;
      
      else if ( x[i][j] < 0.0 )
        x[i][j] += L[j] ;
    }

  }//for ( i=0 ; ...


}
