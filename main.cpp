#define MAIN
#include "globals.h"
void update_positions( void ) ;
void update_positions_init( void ) ;
void initialize( void ) ;
void write_gro( void ) ;
void write_grid( void ) ;
void forces( void ) ;
double integrate( double* ) ;
void write_stress( void ) ;
void calc_Unb( void ) ;
void calc_stress( void ) ;
void bond_stress( void ) ;


int main( int argc , char** argv ) {

  int i,j,k ;

  new_C = 0;

  new_dt =0;

  if ( argc == 3 && !strcmp( "-nt" , argv[1] ) ) {
    nthreads = atoi( argv[2] ) ;
    //omp_set_num_threads( nthreads ) ;
    printf("\nNumber of threads set to %d!\n" , nthreads ) ;
   
  }
  else if(argc == 5 && !strcmp( "-nt" , argv[1] )){
    nthreads = atoi( argv[2] ) ;
    //omp_set_num_threads( nthreads ) ;
    printf("\nNumber of threads set to %d!\n" , nthreads ) ;
    new_C = atof(argv[4]); 
    printf("\n Adjust the chain density to  %f\n" , new_C ) ;
  }
  else if(argc == 7 && !strcmp( "-nt" , argv[1] )){
    nthreads = atoi( argv[2] ) ;
    //omp_set_num_threads( nthreads ) ;
    printf("\nNumber of threads set to %d!\n" , nthreads ) ;
    new_C = atof(argv[4]); 
    printf("\n Adjust the chain density to  %f\n" , new_C ) ;
    new_dt = atof(argv[6]); 
        printf("\n Adjust the new dt  to  %f\n" , new_dt ) ;

  }
  else {
    nthreads = 1 ;
    //omp_set_num_threads( nthreads ) ;
    printf("\nNumber of threads set to %d!\n" , nthreads ) ;
  }


  initialize() ;

  write_gro( ) ;


  FILE *tt_po,*otp, *nbo, *bono ;
  otp = fopen( "data.dat" , "w" ) ;
  nbo = fopen( "nb_stress.dat" , "w" ) ;
  bono = fopen( "bond_stress.dat" , "w" ) ;
  tt_po = fopen( "total_stress.dat" , "w" ) ;

  printf("Entering main loop!\n") ; fflush( stdout ) ;
  for ( step = 0 ; step < nsteps ; step++ )  {


    forces() ;
    if(step !=0){
      update_positions() ;
    }
    else{
	 update_positions_init() ;
    }
   
    if ( stress_freq > 0 && step % stress_freq == 0 ) {
    //  bond_stress() ;
      calc_stress();
      for ( j=0 ; j<Dim ; j++ ) 
        for ( k=0 ; k<Dim ; k++ ) 
          sts_buf[buff_ind][j][k] = Rg3*Ptens[j][k];//Stress_bonds[j][k] ;

      buff_ind++ ;
    }

    
    if ( step > sample_wait && step % sample_freq == 0 ) {
      fftw_fwd( rho[0] , ktmp ) ;
      for ( i=0 ; i<M ; i++ ) {
        double k2, kv[Dim] ;
        k2 = get_k( i , kv ) ;
        avg_sk[0][i] += ktmp[i] * conj(ktmp[i]) ;
      }
      num_averages += 1.0 ;
    }

     if( ((step > sample_wait) and (step % sample_freq == 0)) or step == 0 ){
	write_gro();
     }

    if ( step % print_freq == 0 || step == nsteps-1 ) {

      if(stress_freq <= 0)calc_stress() ;


      printf("step %d of %d  Ubond: %lf P: %lf\n" , step , nsteps , Ubond , Pscalar) ;
      fflush( stdout ) ;
 
      //if ( stress_freq >  0 ) write_gro() ;

      if ( stress_freq > 0 ){
        write_stress() ;
      }

      write_grid_data( "rhoda.dat" , rhoda ) ;
      write_grid_data( "rhodb.dat" , rhodb ) ;
      
      if ( nA > 0.0 )
        write_grid_data( "rhoha.dat" , rhoha ) ;

      if ( nB > 0.0 )
        write_grid_data( "rhohb.dat" , rhohb ) ;

      if ( nP > 0.0 ) 
        write_grid_data( "rhop.dat" , rhop ) ;

      if ( step > sample_wait ) {
        for ( i=0 ; i<M ; i++ )
          ktmp2[i] = avg_sk[0][i] / num_averages ;
 
        write_kspace_data( "avg_sk.dat" , ktmp2 ) ;
      }

      calc_Unb() ;

      fprintf( otp , "%d %lf %lf %lf\n" , step , Ubond , U_kappa_gg , Pscalar) ;
      fflush( otp ) ;

      for ( i=0 ; i<Dim ; i++ )
        fprintf( nbo , "%lf " , Rg3*Stress_nb[i][i] ) ;
      for ( i=0 ; i<Dim ; i++ )
        for ( j=i+1 ; j<Dim ; j++ )
          fprintf( nbo , "%e " , Stress_nb[i][j] ) ;

      fprintf( nbo , "\n" ) ; fflush( nbo ) ;


      for ( i=0 ; i<Dim ; i++ ) 
        fprintf( bono , "%lf " , Rg3*nA*Nha/V-Rg3*Stress_bonds[i][i]) ;
      
      for ( i=0 ; i<Dim ; i++ )
        for ( j=i+1 ; j<Dim ; j++ )
          fprintf( bono, "%e " , V*Stress_bonds[i][j]/(Nha*nA ) ) ;

      fprintf( bono , "\n" ) ; fflush( bono ) ;

     for ( i=0 ; i<Dim ; i++ )
         fprintf( tt_po , "%lf " , Rg3*nA*Nha/V-Rg3*Stress_bonds[i][i] - Rg3*Stress_nb[i][i]) ;
      fprintf( tt_po , "\n" ) ; fflush( tt_po ) ;	 
    }// if step % print_Freq == 0

  }

  fclose( otp ) ;

  return 0 ;

}
