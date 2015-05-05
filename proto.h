// Data Module
int conf2dump( char * );
int in2dump( char * );
int read_parameters( float *, char * );
int read_bin32( char *, float ** );
int read_bin64( char *, double ** );
int data_out( struct void_cell *, struct void_region *, int, float *, char * );
read_void_cell( struct void_cell *, char * );

// Finder Module
int matrix_voids( struct void_cell *, float *, double *, float *, float , float * );
int voids_finder_FA( struct void_cell *, float *, float * );
int voids_finder_DL( struct void_cell *, double *, float * );
int voids_finder_FOF( struct void_cell *, float * );
int regions_builder( struct void_cell *, struct void_region *, float *, int );
int hierarchical_structure( struct void_cell *, struct void_region *, struct void_region *, int, int, float *, char * );

// Tools Module
float effective_radius( int, float, int );
int fractional_anisotropy( float *, float *, float *, float *, int );
int median_filtering_16( float *, int, int , int );
int median_filtering_32( double *, int, int , int );
int neighbour_voids( struct void_cell *, struct void_region *, float *, int, float *, double * );
int boundary_removals( struct void_cell *, struct void_region *, float *, int );
int void_numbers( char * );
int void_cells( char * );
int density_centres( struct void_region *, float *, int, double * );
int fractional_anisotropy_centres( struct void_region *, float *, int, float * );