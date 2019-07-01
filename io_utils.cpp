#include "globals.h"

void write_stress( ) {

  int i, j, k ;

  FILE *otp ;

  if ( step <= print_freq )
    otp = fopen( "stress.dat" , "w" ) ;
  else
    otp = fopen( "stress.dat" , "a" ) ;

  for ( i=0 ; i<buff_ind ; i++ ) {
    // Write diagonals first //
    for ( j=0 ; j<Dim ; j++ )
      fprintf( otp , "%lf " , sts_buf[i][j][j] ) ;

    for ( j=0 ; j<Dim -1 ; j++ )
      for ( k=j+1 ; k<Dim ; k++ )
        fprintf( otp , "%lf " , sts_buf[i][j][k] ) ;

    fprintf( otp , "\n" ) ;

  }

  fclose( otp ) ;

  buff_ind = 0 ;

}


void write_kspace_data( const char *nm , complex<double> *kdt ) {
  int i, j , nn[Dim] ;
  FILE *otp ;
  double kv[Dim], k2 ;
  
  otp = fopen( nm , "w" ) ;

  for ( i=1 ; i<M ; i++ ) {
    unstack( i , nn ) ;

    k2 = get_k( i , kv ) ;

    for ( j=0 ; j<Dim ; j++ ) 
      fprintf( otp , "%lf " , kv[j] ) ;

    fprintf( otp , "%1.5e %1.5e %1.5e %1.5e\n" , abs(kdt[i]), sqrt(k2), 
        real(kdt[i]) , imag(kdt[i]) ) ;

    if ( Dim == 2 && nn[0] == Nx[0]-1 )
      fprintf( otp , "\n" ) ;
  }

  fclose( otp ) ;


}

void write_grid_data( const char *nm , double *dat ) {

  int i, j, nn[Dim] ;
  FILE *otp ;
  double r[Dim] ;
  otp = fopen( nm , "w" ) ;

  for ( i=0 ; i<M ; i++ ) {
    unstack( i , nn ) ;
    
    for ( j=0 ; j<Dim ; j++ )
      fprintf( otp , "%lf " , double(nn[j]) * dx[j] ) ;
    
    fprintf( otp , "%1.16e \n" , dat[i] ) ;
    
    if ( Dim == 2 && nn[0] == Nx[0]-1 )
      fprintf( otp , "\n" ) ;
  }

  fclose( otp ) ;

}


