#include "globals.h"
void allocate( void ) ;
void fft_init( void ) ;
void initialize_potential( void ) ;
void initialize_configuration( void ) ;
void charge_grid( void ) ;
void read_input( void ) ;

void initialize() {
  int i,j ;

  idum =   -long( time(0) ) ; // 9 ; //

  read_input() ;
  
  if(new_C > 0 ) C = new_C;
  
   if(new_dt > 0 ) delt = new_dt;


  phiHA = 1.0 ;

  phiP = phiHB = 0.0 ;

  if ( phiP + phiHA + phiHB > 1.0 )
    die("Invalid volume fractions!\n") ;

  lagrange_weights = 0 ;
  spline_weights = 1 ;

  mem_use = 0. ;

  M = 1 ;
  V = 1.0 ;
  grid_per_partic = 1 ;
  for ( j=0 ; j<Dim ; j++ ) {
    V *= L[j] ;
    M *= Nx[j] ;
    dx[j] = L[j] / double( Nx[j] ) ;
    Lh[j] = 0.5 * L[j] ;
    grid_per_partic *= ( pmeorder + 1 ) ;
    printf("dx: %lf\n" , dx[j] ) ;
  }
  gvol = V / double( M ) ;

  // This is used for density fields //
  ntypes = 3 ;

  Rg = sqrt( double( Nha-1 ) / 6.0 ) ;
  Rg3 = Rg * Rg * Rg ;
  
  Range = Rg *0.2;
  Range2 = Range *Range;
  cout<<"Range "<<Range<<endl;

  rho0 = C * double( Nha ) / Rg3 ;
  
  kappa = kappa * Rg3 / double( Nha * Nha ) ;
  cout << "u0: " << kappa << endl;
  cout << "rho0: " << rho0 << endl;

  chiAB = 0.0 ; // chiAB / double( Nda + Ndb ) ;


  nD = 0 ;
  nB = 0 ;
  
  
  nA = int( rho0 * V / double(Nha) ) ;

  Vp = rho0 ;
  if ( Dim == 2 )
    Vp *= PI * Rp * Rp ;
  else if ( Dim == 3 )
    Vp *= 4.0 * PI * Rp * Rp * Rp / 3.0 ;

  Diff[2] = 1.0 / Vp ;



  nP = 0 ; //int( phiP * rho0 * V / Vp * CG_ratio ) ;

  
  verlet_a = (1 - delt/2.0/Diff[0])/(1 + delt/2.0/Diff[0]); 
  verlet_b = 1.0 / (1 + delt/2.0/Diff[0]);


  nsD = nD * (Nda + Ndb) ;
  nsA = nA * Nha ;
  nsB = nB * Nhb ;

  printf("Input rho0: %lf , " , rho0 ) ;
  rho0 = ( nD * (Nda + Ndb ) + nA * Nha + nB * Nhb + nP * Vp ) / V / CG_ratio ;
  printf("actual rho0: %lf\n" , rho0 ) ;

  printf("\nnD: %d\nnA: %d\nnB: %d\nnP: %d\n\n" , nD, nA, nB, nP ) ;

  // Derived quantities //
  nstot = nA * Nha + nB * Nhb + nD * ( Nda + Ndb ) + nP ;
 
  step = 0 ;
  num_averages = 0.0 ;
  
  buff_ind = 0 ;
  if ( stress_freq > print_freq )
    stress_freq = print_freq ;

  if ( stress_freq > 0 )
    buff_size = print_freq / stress_freq + 1;
  else
    buff_size = 0 ;

  I = complex<double>( 0.0 , 1.0 ) ;

    printf("Total segments: %d\n" , nstot ) ;
    printf("grid_vol: %lf\n" , gvol ) ;
    printf("Particles per grid point: %lf\n" , double(nstot) / double(M) ) ;

  fft_init() ;
    printf("FFTW-MPI Initialized\n") ; fflush( stdout ) ;

  allocate() ;
  
  printf("Memory allocated: %lf MB\n" , mem_use / 1.0E6 ) ; 
  fflush(stdout) ;
  
  initialize_configuration() ;
  
  printf("Initial config generated\n") ; fflush(stdout) ;

  charge_grid() ;

  printf("grid charged\n"); fflush(stdout); 
  
  initialize_potential() ;

  printf("potentials initialized, written\n") ; fflush( stdout ) ; 
}



void initialize_potential( ) {

  int i, j;

  double ro[Dim], rc[Dim], dr[Dim], mdr2 , pref , mdr , k2, kv[Dim] ;
  pref = V / ( pow( 2.0 * sqrt(PI) , Dim ) ) ; // Note: the factor of V comes from the FFT

  for ( j=0 ; j<Dim ; j++ )
    ro[j] = 0.0 ;

  for ( i=0 ; i<M ; i++ ) {

    get_r( i , rc ) ;

    mdr2 = pbc_mdr2( ro, rc, dr ) ;
    mdr = sqrt( mdr2 ) ;

    uG[i] = exp( -mdr2 / 4.0 /Range2) * pref ;
    r_dudr[i] = -mdr2 * exp( -mdr2 / 2.0 /Range2) ;

    tmp[i] = rho0 / 2.0 * ( 1.0 - erf( ( mdr - Rp ) / Xi ) ) * V;
    gammaP[i] = tmp[i] ;

  }

  // Set up the particle-particle potential //
  fftw_fwd( tmp , ktmp ) ;
  for ( i=0 ; i<M ; i++ ) 
    ktmp2[i] = ktmp[i] * ktmp[i] ;
  fftw_back( ktmp2 , uP ) ;

  // Set up particle-polymer potential //
  for ( i=0 ; i<M ; i++ ) {
    k2 = get_k( i , kv ) ;
    ktmp[i] *= exp( -k2 *Range2/ 2.0 ) ;
  }
  fftw_back( ktmp , uPG ) ;


  for ( j=0 ; j<Dim ; j++ ) {
    field_gradient( uG , grad_uG[j] , j ) ;
    field_gradient( uP , grad_uP[j] , j ) ;
    field_gradient( uPG , grad_uPG[j] , j ) ;
  }

  
  int j2 ;
  for ( j=0 ; j<Dim ; j++ ) {

    for ( i=0 ; i<M ; i++ ) {
      get_r( i , rc ) ;
      mdr2 = pbc_mdr2( rc , ro , dr ) ;

      for ( j2=0 ; j2<Dim ; j2++ ) 
        vir_func[j][j2][i] = dr[j2] * -grad_uG[j][i] ;

    }
  }

  for ( j=0 ; j<Dim ; j++ )
    for ( j2=0 ; j2<Dim ; j2++ )
      fftw_fwd( vir_func[j][j2], vir_func_hat[j][j2] ) ;


  write_grid_data( "ug.dat" , uG ) ;
  write_grid_data( "up.dat" , uP ) ;
  write_grid_data( "upg.dat" , uPG ) ;


  for ( j=0 ; j<Dim ; j++ ) {
    char nm[20] ;
    fftw_fwd( grad_uG[j] , grad_uG_hat[j] ) ;
    sprintf( nm , "grad_ug_%d.dat" , j ) ;
    write_grid_data( nm , grad_uG[j] ) ;

    fftw_fwd( grad_uP[j] , grad_uP_hat[j] ) ;
    sprintf( nm , "grad_up_%d.dat" , j ) ;
    write_grid_data( nm , grad_uP[j] ) ;

    fftw_fwd( grad_uPG[j] , grad_uPG_hat[j] ) ;
    sprintf( nm , "grad_upg_%d.dat" , j ) ;
    write_grid_data( nm , grad_uPG[j] ) ;
  }
}

void allocate( ) {

  int i, j ;
  
  x = ( double** ) calloc( nstot , sizeof( double* ) ) ;
  x_bac = ( double** ) calloc( nstot , sizeof( double* ) ) ;
  gn_bac = ( double** ) calloc( nstot , sizeof( double* ) ) ;

  xc = ( char** ) calloc( nstot , sizeof( char* ) ) ;
  f = ( double** ) calloc( nstot , sizeof( double* ) ) ;
  for ( i=0 ; i<nstot ; i++ ) {
    x[i] = ( double* ) calloc( Dim , sizeof( double ) ) ;
    x_bac[i] = ( double* ) calloc( Dim , sizeof( double ) ) ;
    gn_bac[i] = ( double* ) calloc( Dim , sizeof( double ) ) ;
   f[i] = ( double* ) calloc( Dim , sizeof( double ) ) ;
    xc[i] = ( char* ) calloc( 8 , sizeof( char ) ) ;
  }

  mem_use += nstot * 2 * 2 * Dim * sizeof( double ) ; 
  mem_use += nstot * 8 * sizeof( char ) ; 
  
  tp = ( int* ) calloc( nstot , sizeof( int ) );
  grid_inds = ( int** ) calloc( nstot , sizeof( int* ) );
  grid_W = ( double** ) calloc( nstot , sizeof( double* ) ) ;
  for ( i=0 ; i<nstot ; i++ ) {
    grid_inds[i] = ( int* ) calloc( grid_per_partic , sizeof( int ) ) ;
    grid_W[i] = ( double* ) calloc( grid_per_partic , sizeof( double ) ) ;
  }

  mem_use += ( nstot + nstot * grid_per_partic ) * sizeof( int ) ;
  mem_use += nstot * grid_per_partic * sizeof( double ) ;

  sts_buf = ( double*** ) calloc( buff_size , sizeof( double** ) ) ;
  for ( i=0 ; i<buff_size ; i++ ) {
    sts_buf[i] = ( double** ) calloc( Dim , sizeof( double* ) ) ;
    for ( j=0 ; j<Dim ; j++ )
      sts_buf[i][j] = ( double* ) calloc( Dim , sizeof( double ) ) ;
  }
  mem_use += Dim * Dim * buff_size * sizeof( double ) ;


  for ( i=0 ; i<Dim ; i++ ) {
    grad_uG_hat[i] = ( complex<double>* ) calloc( M , sizeof( complex<double> ) ) ;
    grad_uP_hat[i] = ( complex<double>* ) calloc( M , sizeof( complex<double> ) ) ;
    grad_uPG_hat[i] = ( complex<double>* ) calloc( M , sizeof( complex<double> ) ) ;
  }
  
  mem_use += 3 * Dim * M * sizeof( double ) ; 
  
  ktmp = ( complex<double>* ) calloc( M , sizeof( complex<double> ) ) ;
  ktmp2 = ( complex<double>* ) calloc( M , sizeof( complex<double> ) ) ;
  tmp = ( double* ) calloc( M , sizeof( double ) ) ;
  tmp2 = ( double* ) calloc( M , sizeof( double ) ) ;
  uG = ( double* ) calloc( M , sizeof( double ) ) ;
  uP = ( double* ) calloc( M , sizeof( double ) ) ;
  uPG = ( double* ) calloc( M , sizeof( double ) ) ;
  r_dudr = ( double* ) calloc( M , sizeof( double ) ) ;
  
  mem_use += 7 * M * sizeof( double ) ; 

  for ( j=0 ; j<Dim ; j++ ) 
    for ( i=0 ; i<Dim ; i++ ) {
      vir_func[j][i] = ( double* ) calloc ( M , sizeof( double ) ) ;
      vir_func_hat[j][i] = ( complex<double>* ) calloc ( M , sizeof( complex<double> ) ) ;
    }

  mem_use += Dim * Dim * M * sizeof( double ) ;

  rhot = ( double* ) calloc( M , sizeof( double ) ) ;
  rhoha = ( double* ) calloc( M , sizeof( double ) ) ;
  rhohb = ( double* ) calloc( M , sizeof( double ) ) ;
  rhoda = ( double* ) calloc( M , sizeof( double ) ) ;
  rhodb = ( double* ) calloc( M , sizeof( double ) ) ;
  rhop = ( double* ) calloc( M , sizeof( double ) ) ;
  gammaP = ( double* ) calloc( M , sizeof( double ) ) ;
  smrhop = ( double* ) calloc( M , sizeof( double ) ) ;

  mem_use += 8 * M * sizeof( double ) ; 
  
  for ( j=0 ; j<Dim ; j++ )
    for ( i=0 ; i<Dim ; i++ )
      Stress_bonds_t[j][i] = ( double* ) calloc( nthreads , sizeof( double ) ) ;
  
  rhoha_t = ( double** ) calloc( nthreads , sizeof( double* ) ) ;
  rhohb_t = ( double** ) calloc( nthreads , sizeof( double* ) ) ;
  rhoda_t = ( double** ) calloc( nthreads , sizeof( double* ) ) ;
  rhodb_t = ( double** ) calloc( nthreads , sizeof( double* ) ) ;
  rhop_t =  ( double** ) calloc( nthreads , sizeof( double* ) ) ;

  for ( i=0 ; i<nthreads ; i++ ) {
    rhoha_t[i] = ( double* ) calloc( M , sizeof( double ) ) ;
    rhohb_t[i] = ( double* ) calloc( M , sizeof( double ) ) ;
    rhoda_t[i] = ( double* ) calloc( M , sizeof( double ) ) ;
    rhodb_t[i] = ( double* ) calloc( M , sizeof( double ) ) ;
    rhop_t[i] = ( double* ) calloc( M , sizeof( double ) ) ;
  }

  mem_use += 5 * M * nthreads * sizeof( double ) ;

  rho = ( double** ) calloc( ntypes , sizeof( double* ) ) ;
  avg_sk = ( complex<double>** ) calloc( ntypes , sizeof( complex<double>* ) ) ;
  w = ( double** ) calloc( ntypes , sizeof( double* ) ) ;
  for ( i=0 ; i<ntypes ; i++ ) {
    rho[i] = ( double * ) calloc( M , sizeof( double ) ) ;
    avg_sk[i] = ( complex<double> * ) calloc( M , sizeof( complex<double> ) ) ;
    w[i] = ( double* ) calloc( M , sizeof( double ) ) ;
  }

  mem_use += 2 * ntypes * M * sizeof( double ) ;
  mem_use += ntypes * M * sizeof( complex<double> ) ;

  for ( j=0 ; j<Dim ; j++ ) {
    grad_uG[j] = ( double* ) calloc( M , sizeof( double ) ) ;
    grad_uP[j] = ( double* ) calloc( M , sizeof( double ) ) ;
    grad_uPG[j] = ( double* ) calloc( M , sizeof( double ) ) ;
    
    gradwA[j] = ( double* ) calloc( M , sizeof( double ) ) ;
    gradwB[j] = ( double* ) calloc( M , sizeof( double ) ) ;
    gradwP[j] = ( double* ) calloc( M , sizeof( double ) ) ;
  }

  mem_use += Dim * M * 6 * sizeof( double ) ;

}
