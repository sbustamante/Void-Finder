#include <allvars.h>
/*Usage:	
 	   ./Void_Finder.out 
	   *	<eigenvalues_filename> 
 	   *	<density_filename> 
 	   *	<results_folder> 
 	   *	<print_each_void(0-False, 1-True)>
*/

int main( int argc, char *argv[] )
{
    //Properties
    float *eigen1, *eigen2, *eigen3, *FA;
    double *delta;
    //Structure for the cells of simulation
    struct void_cell *void_matrix;
    //Structure of regions
    struct void_region *regions, *regions_backup;
    //Parameters
    float p[NMAX1];
    
    int fa, n_reg_tot, n_reg, n_reg_b, l, n;
    char eig_filename[NMAX1], delta_filename[NMAX1], filename_FA[NMAX1], filename_out[NMAX1];
    float FA_i;
    

    printf( "\n\n******************************** VOID FINDER *******************************\n" );
    //Loading Configuration------------------------------------------------------------------------
    read_parameters( p, "./parameters.conf" );
    
    //Loading and calculating properties of the simulation-----------------------------------------
    //Eigenvalue 1
    sprintf( eig_filename, "%s_1", argv[1] );
    read_bin32( eig_filename, &eigen1 );
    //Eigenvalue 2
    sprintf( eig_filename, "%s_2", argv[1] );
    read_bin32( eig_filename, &eigen2 );
    //Eigenvalue 3
    sprintf( eig_filename, "%s_3", argv[1] );
    read_bin32( eig_filename, &eigen3 );
    //Density contrast field
    sprintf( delta_filename, "%s", argv[2] );
    p[NBOX] = read_bin64( delta_filename, &delta );
    
    
//Density Field!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
for( n=0; n<p[NBOX]*p[NBOX]*p[NBOX]; n++ )
    delta[n] = (double)(eigen1[n] + eigen2[n] + eigen3[n]);
    
    
    //Fractional anisotropy index
    //Allocating memory of FA
    if(!(FA=malloc((int)(p[NBOX]*p[NBOX]*p[NBOX]) * sizeof(float)))){
	fprintf(stderr, "problem with array allocation\n");
	exit(1);}
    fractional_anisotropy( eigen1, eigen2, eigen3, FA, p[NBOX] );
	    
    
    //FOF SCHEME OF CLASSIFICATION*****************************************************************
    if( p[SCHM] == 0 ){
    //Allocating memory of void matrix
    void_matrix = (struct void_cell *)calloc( (int)(p[NBOX]*p[NBOX]*p[NBOX]), \
    sizeof( struct void_cell ) );
    
    //Applying median filtering to FA field
    printf( "\n  * Applying median filtering to FA field: %1.0f iterations\n", p[ITER] );
    median_filtering_16( FA, p[NBOX], p[LNNH], p[ITER] );
    
    //Constructing void regions and hierarchical structure of voids--------------------------------
    //Sweeping all values of FA
    for( fa=0; fa<p[NPFA]; fa++ ){
	//Current FA
	FA_i = p[FAMN] + ( p[FAMX] - p[FAMN] )*fa/(p[NPFA]-1);
	printf( "\n* Hierarchy %d: fractional anisotropy %1.3f\n", fa, FA_i );
		
	//Building void matrix
	matrix_voids( void_matrix, FA, delta, eigen1, FA_i, p );
	//Void finder Algorithm
	n_reg_tot = voids_finder_FOF( void_matrix, p );	//Total number of regions, including empty ones

	//Structure for void regions
	regions = (struct void_region *)calloc( n_reg_tot, sizeof( struct void_region ) );
	for( l=0; l<n_reg_tot; l++ ){
	    regions[l].Ncells = 0;
	    MY_FREE(regions[l].cells);}
	//Building the void regions
	n_reg_b = p[NREG]; //Number of regions for the previous hierarchy
	n_reg = regions_builder( void_matrix, regions, p, n_reg_tot );
			
	//Building Hierarchical structure of voids
	if( fa != 0 )
	    hierarchical_structure( void_matrix, regions, regions_backup, n_reg, n_reg_b, p, filename_FA );
	
	//Print information of each regions (0 - False, 1 - True)
	p[PREG] = atoi(argv[4]);
	//Saving data
	sprintf( filename_FA, "%s_FA_%1.2f", argv[3], FA_i );
	//Printing all information of voids
	data_out( void_matrix, regions, n_reg_tot, p, filename_FA ); //Number of non-empty regions
	
	//Copying the current information of the void regions in a backup structure
	regions_backup = (struct void_region *)calloc( n_reg_tot, sizeof( struct void_region ) );
	memcpy( &regions_backup, &regions, sizeof regions );
    }}
	
	
    //FA WATERSHED SCHEME OF CLASSIFICATION********************************************************
    if( p[SCHM] == 1 )    
    //Sweeping all values of median filtering order
    for( n=0; n<=p[ITER]; n++ ){
	//Allocating memory of void matrix
	void_matrix = (struct void_cell *)calloc( (int)(p[NBOX]*p[NBOX]*p[NBOX]), \
	sizeof( struct void_cell ) );
      
	//Applying median filtering to FA field....................................................
	printf( "\n  * Applying median filtering to FA field: %d iterations\n", n );    
	if( n>0 )
	    median_filtering_16( FA, p[NBOX], p[LNNH], 1 );
	
	//Building void matrix
	matrix_voids( void_matrix, FA, delta, eigen1, 0.0, p );

	//Finding voids............................................................................
	n_reg_tot = voids_finder_FA( void_matrix, FA, p );	
	//Structure for void regions
	MY_FREE( regions );
	regions = (struct void_region *)calloc( n_reg_tot, sizeof( struct void_region ) );
	for( l=0; l<n_reg_tot; l++ ){
	    regions[l].Ncells = 0;
	    MY_FREE(regions[l].cells);}
	//Building the void regions
	n_reg = regions_builder( void_matrix, regions, p, n_reg_tot );
	//Printing all information of voids before boundary removal
	p[PREG] = atoi(argv[4]);
	p[BNDR] = 0;
	sprintf( filename_out, "%s_%d%d", argv[3], n, 0 );
	data_out( void_matrix, regions, n_reg, p, filename_out );

	//Applying boundary removal................................................................
	neighbour_voids( void_matrix, regions, p, n_reg, FA, delta );
	boundary_removals( void_matrix, regions, p, n_reg );
	//Reordering voids
	//Structure for void regions
	MY_FREE( regions );
	regions = (struct void_region *)calloc( n_reg, sizeof( struct void_region ) );
	for( l=0; l<n_reg; l++ ){
	    regions[l].Ncells = 0;
	    MY_FREE(regions[l].cells);}
	n_reg_b = regions_builder( void_matrix, regions, p, n_reg );	
	printf( "  * Number of void regions found before boundary removals %d, and after %d\n", 
	n_reg, n_reg_b-2 );
	//Finding neighbourhood for ordered voids
	neighbour_voids( void_matrix, regions, p, n_reg, FA, delta );
	
	//Calculating density centres of voids
	density_centres( regions, p, n_reg, delta );
	
	//Calculating FA centres of voids
	fractional_anisotropy_centres( regions, p, n_reg, FA );
	
	//Print information of each regions after boundary removal(0 - False, 1 - True)
	p[PREG] = atoi(argv[4]);
	p[BNDR] = 1;
	sprintf( filename_out, "%s_%d%d", argv[3], n, 1 );
	data_out( void_matrix, regions, n_reg, p, filename_out );

	//Releasing memory of void matrix
	MY_FREE( void_matrix );
	}

	
    //DENSITY WATERSHED SCHEME OF CLASSIFICATION***************************************************
    if( p[SCHM] == 2 )
    //Sweeping all values of median filtering order
    for( n=0; n<=p[ITER]; n++ ){
	//Allocating memory of void matrix
	void_matrix = (struct void_cell *)calloc( (int)(p[NBOX]*p[NBOX]*p[NBOX]), \
	sizeof( struct void_cell ) );
	
	//Applying median filtering to the density field...........................................
	printf( "\n  * Applying median filtering to FA field: %d iterations\n", n );
	if( n>0 )
	    median_filtering_32( delta, p[NBOX], p[LNNH], 1 );
	
	//Building void matrix
	matrix_voids( void_matrix, FA, delta, eigen1, 0.0, p );
	
	//Finding voids............................................................................
	n_reg_tot = voids_finder_DL( void_matrix, delta, p );
	//Structure for void regions
	MY_FREE( regions );
	regions = (struct void_region *)calloc( n_reg_tot, sizeof( struct void_region ) );
	for( l=0; l<n_reg_tot; l++ ){
	    regions[l].Ncells = 0;
	    MY_FREE(regions[l].cells);}
	//Building the void regions
	n_reg = regions_builder( void_matrix, regions, p, n_reg_tot );
	//Printing all information of voids before boundary removal
	p[PREG] = atoi(argv[4]);
	p[BNDR] = 0;
	sprintf( filename_out, "%s_%d%d", argv[3], n, 0 );
	data_out( void_matrix, regions, n_reg, p, filename_out );

	//Applying boundary removal................................................................
	neighbour_voids( void_matrix, regions, p, n_reg, FA, delta);
	boundary_removals( void_matrix, regions, p, n_reg );
	//Reordering voids
	//Structure for void regions
	MY_FREE( regions );
	regions = (struct void_region *)calloc( n_reg, sizeof( struct void_region ) );
	for( l=0; l<n_reg; l++ ){
	    regions[l].Ncells = 0;
	    MY_FREE(regions[l].cells);}
	n_reg_b = regions_builder( void_matrix, regions, p, n_reg );
	printf( "  * Number of void regions found before boundary removals %d, and after %d\n", 
	n_reg, n_reg_b-2 );
	//Finding neighbourhood for ordered voids
	neighbour_voids( void_matrix, regions, p, n_reg, FA, delta);

	//Print information of each regions after boundary removal(0 - False, 1 - True)
	p[PREG] = atoi(argv[4]);
	p[BNDR] = 1;
	sprintf( filename_out, "%s_%d%d", argv[3], n, 1 );
	data_out( void_matrix, regions, n_reg, p, filename_out );
	
	//Releasing memory of void matrix
	MY_FREE( void_matrix );
	}

	
    return 0;
}