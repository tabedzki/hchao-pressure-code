#include "globals.h"
void charge_grid( void ) ;
void write_grid_data( const char* , double* ) ;

void random_config( void ) ;
void write_gro( void ) ;
void read_gro( FILE* ) ;



void initialize_configuration( ) {

  int i, j ;

  FILE *inp ;
  inp = fopen( "input.gro" , "r" ) ;
  if ( inp == NULL )
    random_config() ;
  else {
    read_gro( inp ) ;
    fclose( inp ) ;
    printf("input.gro read!\n" ) ;
  }


}

void read_gro( FILE *ip ) {

  int i, m, j, k, di , ind ;
  char tt[80] ;
  fgets( tt , 80 , ip ) ;
  fscanf( ip , "%d\n" , &di ) ;

  if ( di != nstot ) 
    die("Number of sites in input.gro does not match!\n");

  ind = 0 ;

  for ( k=0 ; k<nD ; k++ ) { 
    for ( m=0 ; m<Nda + Ndb ; m++ ) {
      fscanf( ip , "%5d" , &di ) ;
      fscanf( ip , "%5s" , tt ) ;
      fscanf( ip , "%5s" , tt ) ;
      fscanf( ip , "%d" , &di ) ;
     
     
      for ( j=0 ; j<Dim ; j++ ) {
        fscanf( ip , "%lf" , &x[ind][j] ) ;
        x[ind][j] *= Rg ;
      }
     
     
      fgets( tt, 80 , ip ) ;
     
      tp[ ind ] = ( m < Nda ? 0 : 1 ) ;
      
      ind++ ;
    }
  }

  for ( k=0 ; k<nA ; k++ ) {
    for ( m=0 ; m<Nha ; m++ ) {
      fscanf( ip , "%5d" , &di ) ;
      fscanf( ip , "%5s" , tt ) ;
      fscanf( ip , "%5s" , tt ) ;
      fscanf( ip , "%d" , &di ) ;
     
     
      for ( j=0 ; j<Dim ; j++ ) {
        fscanf( ip , "%lf" , &x[ind][j] ) ;
        x[ind][j] *= Rg ;
      }
     
     
      fgets( tt, 80 , ip ) ;
     
      tp[ ind ] = 0 ;
      
      ind++ ;
    }
  }

  for ( k=0 ; k<nB ; k++ ) {
    for ( m=0 ; m<Nhb ; m++ ) {
      fscanf( ip , "%5d" , &di ) ;
      fscanf( ip , "%5s" , tt ) ;
      fscanf( ip , "%5s" , tt ) ;
      fscanf( ip , "%d" , &di ) ;
     
     
      for ( j=0 ; j<Dim ; j++ ) {
        fscanf( ip , "%lf" , &x[ind][j] ) ;
        x[ind][j] *= Rg ;
      }
     
     
      fgets( tt, 80 , ip ) ;
     
      tp[ ind ] = 1 ;
      
      ind++ ;
    }
  }

  // Assign the labels //
  for ( i=0 ; i<nstot ; i++ ) {
    if ( tp[i] == 0 ) 
      xc[i] = "H" ;
    else if ( tp[i] == 1 )
      xc[i] = "He" ;
    else if ( tp[i] == 2 )
      xc[i] = "O" ;
    else if ( tp[i] == 3 )
      xc[i] = "S" ;
    else if ( tp[i] == 4 )
      xc[i] = "N" ;
    else if ( tp[i] == 5 )
      xc[i] = "Br" ;
    else if ( tp[i] == 6 )
      xc[i] = "C" ;
    else if ( tp[i] == 7 )
      xc[i] = "Na" ;
    else if ( tp[i] == 8 )
      xc[i] = "P" ;
    else if ( tp[i] == 9 )
      xc[i] = "Ca" ;
  }
  
}


void write_gro() {

  FILE *otp ;
  int i, j, ind, k, m, resind ;
  if ( step == 0 ) 
    otp = fopen( "output.gro" , "w" ) ;
  else
    otp = fopen( "output.gro" , "a" ) ;

  fprintf( otp , "%lf = steps * delt, Writing frame %d\n" , double(step) * delt , step ) ;
  fprintf( otp , "%d\n" , nstot );

  ind = 0 ;
  resind = 0 ;
  for ( k=0 ; k<nD ; k++ ) {
    for ( m=0 ; m<Nda + Ndb ; m++ ) {
      fprintf( otp , "%5d" , resind % 100000 ) ;
      fprintf( otp , "%5s" , "BCP" ) ;
      fprintf( otp , "%5s" , xc[ind] ) ;
     
      fprintf( otp , "%5d" , ind % 100000 ) ;
      
      for ( j=0 ; j<Dim ; j++ )
        fprintf( otp , "%8.3lf " , x[ind][j]  ) ;
      
      for ( j=Dim ; j<3 ; j++ )
        fprintf( otp , "%8.3lf " , 0.0 );
      
      fprintf( otp , "\n" ) ;

      ind++;
    }
    resind++ ;
  }
  
  for ( k=0 ; k<nA ; k++ ) {
    for ( m=0 ; m<Nha ; m++ ) {
      fprintf( otp , "%5d" , resind % 100000 ) ;
      fprintf( otp , "%5s" , "HA " ) ;
      fprintf( otp , "%5s" , xc[ind] ) ;
     
      fprintf( otp , "%5d" , ind % 100000 ) ;
      
      for ( j=0 ; j<Dim ; j++ )
        fprintf( otp , "%8.3lf " , x[ind][j] ) ;
      
      for ( j=Dim ; j<3 ; j++ )
        fprintf( otp , "%8.3lf " , 0.0 );
      
      fprintf( otp , "\n" ) ;

      ind++;
    }
    resind++ ;
  }
  
  for ( k=0 ; k<nB ; k++ ) {
    for ( m=0 ; m<Nhb ; m++ ) {
      fprintf( otp , "%5d" , resind % 100000 ) ;
      fprintf( otp , "%5s" , "HB " ) ;
      fprintf( otp , "%5s" , xc[ind] ) ;
     
      fprintf( otp , "%5d" , ind % 100000 ) ;
      
      for ( j=0 ; j<Dim ; j++ )
        fprintf( otp , "%8.3lf " , x[ind][j]) ;
      
      for ( j=Dim ; j<3 ; j++ )
        fprintf( otp , "%8.3lf " , 0.0 );
      
      fprintf( otp , "\n" ) ;

      ind++;
    }
    resind++ ;
  }

  for ( k=0 ; k<nP ; k++ ) {
    fprintf( otp , "%5d" , resind % 100000 ) ;
    fprintf( otp , "%5s" , "HP " ) ;
    fprintf( otp , "%5s" , xc[ind] ) ;
    
    fprintf( otp , "%5d" , ind % 100000 ) ;
    
    for ( j=0 ; j<Dim ; j++ )
      fprintf( otp , "%8.3lf " , x[ind][j]  ) ;
    
    for ( j=Dim ; j<3 ; j++ )
      fprintf( otp , "%8.3lf " , 0.0 );
    
    fprintf( otp , "\n" ) ;

    ind++ ;
    resind++ ;
  }

  for ( j=0 ; j<Dim ; j++ )
    fprintf( otp , "%lf " , L[j]  ) ;
  for ( j=Dim ; j<3 ; j++ )
    fprintf( otp , "1.0" ) ;
  fprintf( otp , "\n" ) ;


  fclose( otp ) ;

}

void random_config( void ) {

  int i, m, j, k, ind = 0 ;

  for ( k=0 ; k<nD ; k++ ) { 
    for ( j=0 ; j<Dim ; j++ ) 
      x[ind][j] = ran2() * L[j] ;
    
    tp[ ind ] = ( Nda > 0 ? 0 : 1 ) ;

    ind += 1 ;
    
    for ( m=1 ; m<Nda + Ndb ; m++ ) {
    
      for ( j=0 ; j<Dim ; j++ ) {
        x[ind][j] = x[ ind-1 ][ j ] + gasdev2() ;

        if ( x[ind][j] > L[j] )
          x[ind][j] -= L[j] ;
        else if ( x[ind][j] < 0.0 )
          x[ind][j] += L[j] ;
      }

      tp[ ind ] = ( m < Nda ? 0 : 1 ) ;

      ind++ ;
    
    }
  }

  // Random A homopolymers //
  for ( k=0 ; k<nA ; k++ ) {
    for ( j=0 ; j<Dim ; j++ ) 
      x[ind][j] = ran2() * L[j] ;
    
    tp[ ind ] = 0 ;

    ind += 1 ;
 
    for ( m=1 ; m<Nha ; m++ ) {
    
      for ( j=0 ; j<Dim ; j++ ) {
        x[ind][j] = x[ ind-1 ][ j ] + gasdev2() ;

        if ( x[ind][j] > L[j] )
          x[ind][j] -= L[j] ;
        else if ( x[ind][j] < 0.0 )
          x[ind][j] += L[j] ;
      }

      tp[ ind ] = 0 ;

      ind++ ;
    
    }
  }

  // Random B homopolymers //
  for ( k=0 ; k<nB ; k++ ) {
    for ( j=0 ; j<Dim ; j++ ) 
      x[ind][j] = ran2() * L[j] ;
    
    tp[ ind ] = 1 ;

    ind += 1 ;
 
    for ( m=1 ; m<Nhb ; m++ ) {
    
      for ( j=0 ; j<Dim ; j++ ) {
        x[ind][j] = x[ ind-1 ][ j ] + gasdev2() ;

        if ( x[ind][j] > L[j] )
          x[ind][j] -= L[j] ;
        else if ( x[ind][j] < 0.0 )
          x[ind][j] += L[j] ;
      }

      tp[ ind ] = 1 ;

      ind++ ;
    
    }
  }

  // Random particle centers //
  for ( k=0 ; k<nP ; k++ ) {
    for ( j=0 ; j<Dim ; j++ )
      x[ind][j] = ran2() * L[j] ;

    tp[ ind ] = 2 ;

    ind += 1 ;
  }


  // Assign the labels //
  for ( i=0 ; i<nstot ; i++ ) {
    if ( tp[i] == 0 ) 
      xc[i] = "H" ;
    else if ( tp[i] == 1 )
      xc[i] = "He" ;
    else if ( tp[i] == 2 )
      xc[i] = "O" ;
    else if ( tp[i] == 3 )
      xc[i] = "S" ;
    else if ( tp[i] == 4 )
      xc[i] = "N" ;
    else if ( tp[i] == 5 )
      xc[i] = "Br" ;
    else if ( tp[i] == 6 )
      xc[i] = "C" ;
    else if ( tp[i] == 7 )
      xc[i] = "Na" ;
    else if ( tp[i] == 8 )
      xc[i] = "P" ;
    else if ( tp[i] == 9 )
      xc[i] = "Ca" ;
  }
  printf("config generated!\n") ; fflush( stdout ) ;
}

