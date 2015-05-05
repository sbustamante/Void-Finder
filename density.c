#include <allvars.h>
/*Usage:	
 	   ./Density_Profile.out 
	   *	Folder of voids
	   *	Density filename
*/

int main( int argc, char *argv[] )
{
    double *delta;
    //Structure for the cells of simulation
    struct void_cell *void_matrix;
    //Parameters
    float p[NMAX1];
    //Density
    FILE *delta_void;

    int i_v, i_c, void_size, n_voids, N;
    long long int n;
    char void_name[NMAX1], delta_filename[NMAX1];
    
    printf( "\n\n******************************** VOID DENSITY *******************************\n" );
    //Loading Configuration------------------------------------------------------------------------
    read_parameters( p, "./parameters.conf" );
    //Number of voids
    n_voids = void_numbers( argv[1] );
    //Loading density contrast field---------------------------------------------------------------
    sprintf( delta_filename, "%s", argv[2] );
    p[NBOX] = read_bin64( delta_filename, &delta );
    N = (int)p[NBOX];
        
    //Storing density field of each void
    for( i_v=1; i_v<n_voids; i_v++ ){
	printf( "  * In void %d\n", i_v );
      
	//Label of this void
	sprintf( void_name, "%s/void_%d.dat", argv[1], i_v );
	//Number of cells of this void
	void_size = void_cells( void_name );
      	
	//Allocating memory of void matrix
	if( i_v==1 )
	    void_matrix = (struct void_cell *)calloc( void_size, sizeof( struct void_cell ) );
	//Reading file
	read_void_cell( void_matrix, void_name );
	
	//Label of the density
	sprintf( void_name, "%s/void_%d_rho.dat", argv[1], i_v );
	delta_void = fopen( void_name, "w" );
	
	for( i_c=0; i_c<void_size; i_c++ ){
	    n = void_matrix[i_c].k + N*(void_matrix[i_c].j + N*void_matrix[i_c].i);
	    fprintf( delta_void, "%1.5e\n", delta[n] );
	}
	fclose( delta_void );
    }
    return 0;
}