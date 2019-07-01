#include "globals.h"

void read_input( void ) {

  FILE *inp ;
  int i,j;
  double d1 ;

  char tt[80] ;

  inp = fopen( "dyft.input" , "r" ) ;


  fscanf( inp , "%lf %d" , &phiHA , &Nha ) ;
  fgets( tt , 80 , inp ) ;

  fscanf( inp , "%lf" , &CG_ratio ) ;
  fgets( tt , 80 , inp ) ;

  
  // Blank line //
  fgets( tt , 80 , inp ) ;

  fscanf( inp , "%lf" , &phiP ) ;
  fgets( tt , 80 , inp ) ;

  fscanf( inp , "%lf" , &Rp ) ;
  fgets( tt , 80 , inp ) ;

  fscanf( inp , "%lf" , &Xi ) ;
  fgets( tt , 80 , inp ) ;

  fscanf( inp , "%d" , &A_partics ) ;
  fgets( tt , 80 , inp ) ;
  
  
  // Blank line //
  fgets( tt , 80 , inp ) ;


  fscanf( inp , "%lf" , &C ) ;
  fgets( tt , 80 , inp ) ;

  fscanf( inp , "%lf" , &kappa ) ;
  fgets( tt , 80 , inp ) ;

  Diff = ( double* ) calloc( 3 , sizeof( double ) ) ;
  fscanf( inp , "%lf %lf" , &Diff[0] , &Diff[1] ) ;
  fgets( tt , 80 , inp ) ;


  // Blank line //
  fgets( tt , 80 , inp ) ;


  for ( i=0 ; i<Dim ; i++ ) 
    fscanf( inp , "%lf" , &L[i] ) ; 
  fgets( tt , 80 , inp ) ;

  for ( i=0 ; i<Dim ; i++ ) 
    fscanf( inp , "%d" , &Nx[i] ) ; 
  fgets( tt , 80 , inp ) ;


  fscanf( inp , "%lf", &delt ) ;
  fgets( tt , 80 , inp ) ;

  fscanf( inp , "%d" , &pmeorder ) ;
  fgets( tt , 80 , inp ) ;

  cout<<"pmeorder "<<pmeorder<<endl;
  // Blank line //
  fgets( tt , 80 , inp ) ;


  fscanf( inp , "%d" , &nsteps ) ;
  fgets( tt , 80 , inp ) ;
  fscanf( inp , "%d" , &print_freq ) ;
  fgets( tt , 80 , inp ) ;
  fscanf( inp , "%d" , &pre_equil_steps ) ;
  fgets( tt , 80 , inp ) ;
  fscanf( inp , "%d %d" , &sample_wait , &sample_freq ) ;
  fgets( tt , 80 , inp ) ;
  fscanf( inp , "%d" , &stress_freq ) ;
  printf("stress_freq: %d\n" , stress_freq ) ;
  fclose( inp ) ;


}
