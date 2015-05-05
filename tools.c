#include <allvars.h>

/**************************************************************************************************
 NAME:	     effective_radius
 FUNCTION:   This function computes the effective radius of a voids given the number of belonging 
	     cells
 INPUTS:     int with number of cells, comoving lenght of the simulation, number of cells per side
 RETURN:     r_eff
**************************************************************************************************/
float effective_radius( int Ncells, 
			float Lbox,
			int Nbox )
{ 
    float r_eff;
    
    r_eff = pow(3*Ncells/(4*M_PI),1/3.0)*Lbox/Nbox;

    return r_eff;
}


/**************************************************************************************************
 NAME:	     fractional_anisotropy
 FUNCTION:   This function computes the FA index for all the grid of the simulation. As defined by 
	     Libeskind, et.al, 2013.
 INPUTS:     The three eigenvalues (float arrays)
 RETURN:     0
**************************************************************************************************/
int fractional_anisotropy( float *eig1, 
			   float *eig2, 
			   float *eig3, 
			   float *FA, 
			   int N )
{ 
    int i,j,k;
    long long int n;

    //Sweeping the grid    
    for( i=0;i<N;i++ )
    for( j=0;j<N;j++ )
    for( k=0;k<N;k++ ){
	//Overall index
	n = k + N*(j + N*i);
	FA[n] = (1/pow(3,0.5))*pow( ( pow(eig1[n] - eig3[n],2) + pow(eig2[n] - eig3[n],2) +\
	pow(eig1[n] - eig2[n],2) )/(eig1[n]*eig1[n] + eig2[n]*eig2[n] + eig3[n]*eig3[n]), 0.5 );}

    return 0;
}


/**************************************************************************************************
 NAME:	     median_filtering_16
 FUNCTION:   This function performs a median softening filtering process over some field sample 
	     onto a given float 3D grid.
 INPUTS:     Field to be softened, box resolution, number of neighbours, number of iterations of 
	     the filtering process.
 RETURN:     0
**************************************************************************************************/
int median_filtering_16( float *field,
			 int N,
			 int neighbours,
			    int iter )
{ 
    int i,j,k;
    int ic,jc,kc;		//Neighbour counters
    int it,jt,kt;		//Corrected neighbour counters
    long long int n, nt;
    
    float *processed_field;	//Processed field
    float ngh_values[NMAX1], ngh_tmp;
    int i_iter, n_ngh, pass;

    //Allocating memory of processed_field
    if(!(processed_field=malloc(N*N*N * sizeof(float)))){
	fprintf(stderr, "problem with array allocation\n");
	exit(1);}
    
    for( i_iter=0; i_iter<iter; i_iter++ ){
	//Sweeping the grid    
	for( i=0;i<N;i++ )
	for( j=0;j<N;j++ )
	for( k=0;k<N;k++ ){
	    //Overall index
	    n = k + N*(j + N*i);
	    
	    //Setting the neighbourhood
	    n_ngh = 0;
	    for( ic=-neighbours; ic<=neighbours; ic++ )
	    for( jc=-neighbours; jc<=neighbours; jc++ )
	    for( kc=-neighbours; kc<=neighbours; kc++ ){
	      
		it = i + ic; jt = j + jc; kt = k + kc;
		//Neighbor out of limits (Periodic boundary conditions) (X direction)
		if( i+ic>=N )		it = 0;
		if( i+ic<0 )		it = N-1;
		//Neighbor out of limits (Periodic boundary conditions) (Y direction)
		if( j+jc>=N )		jt = 0;
		if( j+jc<0 )		jt = N-1;
		//Neighbor out of limits (Periodic boundary conditions) (Z direction)
		if( k+kc>=N )		kt = 0;
		if( k+kc<0 )		kt = N-1;

		nt = kt + N*(jt + N*it);
		
		ngh_values[n_ngh] = field[nt];
		n_ngh ++;}
		
	    //Finding median value of neighbourhood
	    pass = 1;
	    while( pass == 1 ){
		pass = 0;
		for( n_ngh = 1; n_ngh<pow(2*neighbours+1,3); n_ngh++ )
		    if( ngh_values[n_ngh-1] > ngh_values[n_ngh] ){
			ngh_tmp = ngh_values[n_ngh];
			ngh_values[n_ngh] = ngh_values[n_ngh-1];
			ngh_values[n_ngh-1] = ngh_tmp;
			pass = 1;}}
	    processed_field[n] = ngh_values[(int)(pow(2*neighbours+1,3)/2.0)];}

	//Applying changes    
	for( i=0;i<N;i++ )
	for( j=0;j<N;j++ )
	for( k=0;k<N;k++ ){
	    //Overall index
	    n = k + N*(j + N*i);
	    field[n] = processed_field[n];}}
	    
    //Releasing memory of processed field
    MY_FREE(processed_field);
	
    return 0;
}


/**************************************************************************************************
 NAME:	     median_filtering_32
 FUNCTION:   This function performs a median softening filtering process over some field sample 
	     onto a given double 3D grid.
 INPUTS:     Field to be softened, box resolution, number of neighbours, number of iterations of 
	     the filtering process.
 RETURN:     0
**************************************************************************************************/
int median_filtering_32( double *field,
			 int N,
			 int neighbours,
			 int iter )
{ 
    int i,j,k;
    int ic,jc,kc;		//Neighbour counters
    int it,jt,kt;		//Corrected neighbour counters
    long long int n, nt;
    
    double *processed_field;	//Processed field
    double ngh_values[NMAX1], ngh_tmp;
    int i_iter, n_ngh, pass;

    //Allocating memory of processed_field
    if(!(processed_field=malloc(N*N*N * sizeof(double)))){
	fprintf(stderr, "problem with array allocation\n");
	exit(1);}
    
    for( i_iter=0; i_iter<iter; i_iter++ ){
	//Sweeping the grid    
	for( i=0;i<N;i++ )
	for( j=0;j<N;j++ )
	for( k=0;k<N;k++ ){
	    //Overall index
	    n = k + N*(j + N*i);
	    
	    //Setting the neighbourhood
	    n_ngh = 0;
	    for( ic=-neighbours; ic<=neighbours; ic++ )
	    for( jc=-neighbours; jc<=neighbours; jc++ )
	    for( kc=-neighbours; kc<=neighbours; kc++ ){
	      
		it = i + ic; jt = j + jc; kt = k + kc;
		//Neighbor out of limits (Periodic boundary conditions) (X direction)
		if( i+ic>=N )		it = 0;
		if( i+ic<0 )		it = N-1;
		//Neighbor out of limits (Periodic boundary conditions) (Y direction)
		if( j+jc>=N )		jt = 0;
		if( j+jc<0 )		jt = N-1;
		//Neighbor out of limits (Periodic boundary conditions) (Z direction)
		if( k+kc>=N )		kt = 0;
		if( k+kc<0 )		kt = N-1;

		nt = kt + N*(jt + N*it);
		
		ngh_values[n_ngh] = field[nt];
		n_ngh ++;}
		
	    //Finding median value of neighbourhood
	    pass = 1;
	    while( pass == 1 ){
		pass = 0;
		for( n_ngh = 1; n_ngh<pow(2*neighbours+1,3); n_ngh++ )
		    if( ngh_values[n_ngh-1] > ngh_values[n_ngh] ){
			ngh_tmp = ngh_values[n_ngh];
			ngh_values[n_ngh] = ngh_values[n_ngh-1];
			ngh_values[n_ngh-1] = ngh_tmp;
			pass = 1;}}
	    processed_field[n] = ngh_values[(int)(pow(2*neighbours+1,3)/2.0)];}

	//Applying changes    
	for( i=0;i<N;i++ )
	for( j=0;j<N;j++ )
	for( k=0;k<N;k++ ){
	    //Overall index
	    n = k + N*(j + N*i);
	    field[n] = processed_field[n];}}
	    
    //Releasing memory of processed field
    MY_FREE(processed_field);
	
    return 0;
}


/**************************************************************************************************
 NAME:	     neighbour_voids
 FUNCTION:   Finds neighbour voids for each void
 INPUTS:     Matrix with voids (1-void, 0-no void), array with empty regions, array of parameters,
	     number of regions, fractional anisotropy array, density array
 RETURN:     0
**************************************************************************************************/
int neighbour_voids( struct void_cell void_matrix[],  
		     struct void_region regions[],
		     float *p,
		     int n_reg,
		     float *FA,
		     double *delta)
{
    int i,j,k,lcell;	//Original counters
    int ic,jc,kc;	//Neighbour counters
    int it,jt,kt;	//Corrected neighbour counters
    
    int N=p[NBOX], b=p[LINK], new_neigh, reeval;
    long int l, ln, m, mn, mt;
    
    long long int n, nt;
    long int id_reg_n, id_reg_nt;
  
    //Sweeping void regions for detecting neighbour voids
    for( l=0; l<n_reg; l++ ){
	//Initial number of voids set in 0
	regions[l].Nneighbours = 0;	

	//Sweeping internal cells
	for( lcell=0; lcell<regions[l].Ncells; lcell++ ){
	    //Overall index of this cell
	    n = regions[l].cells[lcell];
	    //Cartesian indexes of this cell
	    i = void_matrix[n].i;
	    j = void_matrix[n].j;
	    k = void_matrix[n].k;
	    
	    //Setting the cell neighbourhood
	    for( ic=-b; ic<=b; ic++ )
	    for( jc=-b; jc<=b; jc++ )
	    for( kc=-b; kc<=b; kc++ )
	    if( ic!=0 || jc!=0 || kc!=0 ){
	      
		it = i + ic; jt = j + jc; kt = k + kc;
		//Neighbour out of limits (Periodic boundary conditions) (X direction)
		if( i+ic>=N )		it = 0;
		if( i+ic<0 )		it = N-1;
		//Neighbour out of limits (Periodic boundary conditions) (Y direction)
		if( j+jc>=N )		jt = 0;
		if( j+jc<0 )		jt = N-1;
		//Neighbour out of limits (Periodic boundary conditions) (Z direction)
		if( k+kc>=N )		kt = 0;
		if( k+kc<0 )		kt = N-1;
		    
		nt = kt + N*(jt + N*it);
		//If the neighbour cell is also a void
		if( void_matrix[nt].isvoid == 1 ){
		    id_reg_n = void_matrix[n].id_reg;
		    id_reg_nt = void_matrix[nt].id_reg;
		    //It is a new neighbour void
		    if( id_reg_n != id_reg_nt ){
			//First neighbour
			if( regions[l].Nneighbours == 0 ){
			    regions[l].Nneighbours = 1;
			    //Initializing regions
			    regions[l].neighbours = 
			    malloc( regions[l].Nneighbours*sizeof( long int ) );
			    regions[l].neighbours[ regions[l].Nneighbours-1 ] = id_reg_nt;
			    //Initializing neighbour cell
			    regions[l].Nneigh_cells = (int *)calloc( regions[l].Nneighbours, 
			    sizeof( int ) );
			    regions[l].Nneigh_cells[regions[l].Nneighbours-1] = 1;
			    //Initializing mean value across boundaries
			    regions[l].neigh_mean = (float *)calloc( regions[l].Nneighbours, 
			    sizeof( float ) );
			    //If FA Watershed scheme
			    if( p[SCHM] == 1 )
				regions[l].neigh_mean[regions[l].Nneighbours-1] = 
				0.5*(FA[nt] + FA[n]);
			    //If density Watershed scheme
			    if( p[SCHM] == 2 )
				regions[l].neigh_mean[regions[l].Nneighbours-1] = 
				(float)(0.5*(delta[nt]+delta[n]));
			    }
			//Other voids
			else{
			    //Flag for new neighbour voids
			    new_neigh = 0;
			    for( m=0; m<regions[l].Nneighbours; m++ )
				//An already existing void
				if( regions[l].neighbours[m] == id_reg_nt ){
				    regions[l].Nneigh_cells[m] ++;
				    //If FA Watershed scheme
				    if( p[SCHM] == 1 )
					regions[l].neigh_mean[m] += 
					0.5*(FA[nt] + FA[n]);;
				    //If density Watershed scheme
				    if( p[SCHM] == 2 )
					regions[l].neigh_mean[m] += 
					(float)(0.5*(delta[nt]+delta[n]));
				    new_neigh = 1;
				    }
			    //New neighbour void
			    if( new_neigh == 0 ){
				regions[l].Nneighbours ++;
				//Initializing regions
				regions[l].neighbours = 
				(long int *)realloc( regions[l].neighbours,
				regions[l].Nneighbours*sizeof( long int ) );
				regions[l].neighbours[ regions[l].Nneighbours-1 ] = id_reg_nt;
				//Initializing neighbour cell
				regions[l].Nneigh_cells = 
				(int *)realloc( regions[l].Nneigh_cells, 
				regions[l].Nneighbours*sizeof( int ) );
				regions[l].Nneigh_cells[regions[l].Nneighbours-1] = 1;
				//Initializing mean value across boundaries
				regions[l].neigh_mean = 
				(float *)realloc( regions[l].neigh_mean,
				regions[l].Nneighbours*sizeof( float ) );
				//If FA Watershed scheme
				if( p[SCHM] == 1 )
				    regions[l].neigh_mean[regions[l].Nneighbours-1] = 
				    0.5*(FA[nt] + FA[n]);;
				//If density Watershed scheme
				if( p[SCHM] == 2 )
				    regions[l].neigh_mean[regions[l].Nneighbours-1] = 
				    	    (float)(0.5*(delta[nt]+delta[n]));
				}
			    }
			}
		    }}}}
    return 0;
}


/**************************************************************************************************
 NAME:	     boundary_removals
 FUNCTION:   Removes boundaries below some user-defined threshold 
 INPUTS:     Matrix with voids (1-void, 0-no void), array with empty regions, array of parameters,
	     number of regions, fractional anisotropy array, density array
 RETURN:     0
**************************************************************************************************/
int boundary_removals( struct void_cell void_matrix[],  
		       struct void_region regions[],
		       float *p,
		       int n_reg)
{
    int i,j,k,lcell;	//Original counters
    int ic,jc,kc;	//Neighbour counters
    int it,jt,kt;	//Corrected neighbour counters
    
    int N=p[NBOX], b=p[LINK], new_neigh, reeval;
    long int l, ln, m, mn, mt;
    
    long long int n, nt;
    long int id_reg_n, id_reg_nt;
  
  
    //Sweeping void regions for merging with neighbours
    for( l=0; l<n_reg; l++ ){
	//Flag to reevaluate this void region
	reeval = 0;
	if( regions[l].index == l )
	for( m=0; m<regions[l].Nneighbours; m++ ){
	    //Index of neighbour void region
	    ln = regions[ regions[l].neighbours[m] ].index;
	    //Detecting threshold value across boundary
	    if( ln != l && regions[l].Nneigh_cells[m] != 0 ) //Not repeated region
	    if( (p[SCHM] == 1 && regions[l].neigh_mean[m]/regions[l].Nneigh_cells[m] > p[FAMG]) || 
		(p[SCHM] == 2 && regions[l].neigh_mean[m]/regions[l].Nneigh_cells[m] > p[DLMG]) ){
		
		//This void region has to be reevaluated
		reeval = 1;
	      
		//Discarting this neighbour for future iterations
		regions[ regions[l].neighbours[m] ].index = l;
	      
		//Merging both regions (Information of neighbours).................................
		//Reallocing regions
		regions[l].neighbours = 
		(long int *)realloc( regions[l].neighbours,
		(regions[l].Nneighbours+regions[ln].Nneighbours)*sizeof( long int ) );
		//Reallocing neighbour cell
		regions[l].Nneigh_cells = 
		(int *)realloc( regions[l].Nneigh_cells, 
		(regions[l].Nneighbours+regions[ln].Nneighbours)*sizeof( int ) );
		//Reallocing mean value across boundaries
		regions[l].neigh_mean = 
		(float *)realloc( regions[l].neigh_mean,
		(regions[l].Nneighbours+regions[ln].Nneighbours)*sizeof( float ) );
		
		//Passing information
		for( mn=0; mn<regions[ln].Nneighbours; mn++ ){
		    regions[l].neighbours[ regions[l].Nneighbours + mn ] = regions[ln].neighbours[mn];
		    regions[l].Nneigh_cells[ regions[l].Nneighbours + mn ] = regions[ln].Nneigh_cells[mn];
		    regions[l].neigh_mean[ regions[l].Nneighbours + mn ] = regions[ln].neigh_mean[mn];
		    //Removing reciprocal neighbour
		    if( regions[ regions[ln].neighbours[mn] ].index == l ){
			regions[l].Nneigh_cells[ regions[l].Nneighbours + mn ] = 0;
			regions[l].neigh_mean[ regions[l].Nneighbours + mn ] = 0.0;}}
		    
		//Detecting already existing neighbours
		for( mn=0; mn<(regions[l].Nneighbours+regions[ln].Nneighbours); mn++ )
		    for( mt=mn+1; mt<(regions[l].Nneighbours+regions[ln].Nneighbours); mt++ )
			//A repeated region
			if( regions[l].neighbours[mn] == regions[l].neighbours[mn] ){
			    regions[l].Nneigh_cells[mn] += regions[l].Nneigh_cells[mt];
			    regions[l].Nneigh_cells[mt] = 0;
			    regions[l].neigh_mean[mn] += regions[l].neigh_mean[mt];
			    regions[l].neigh_mean[mt] = 0.0;}
			    
		//Releasing memory of current neighbour region
		MY_FREE( regions[ln].neighbours );
		MY_FREE( regions[ln].Nneigh_cells );
		MY_FREE( regions[ln].neigh_mean );
		//New number of neighbours
		regions[l].Nneighbours += regions[ln].Nneighbours;
		regions[ln].Nneighbours = 0;
		
		//Merging both regions (Information of cells)......................................
		//Reallocing cells
		regions[l].cells = 
		(long long int *)realloc( regions[l].cells,
		(regions[l].Ncells+regions[ln].Ncells)*sizeof( long long int ) );
		//Passing information of cells
		for( n=0; n<regions[ln].Ncells; n++ ){
		    regions[l].cells[ regions[l].Ncells + n ] = regions[ln].cells[n];
		    void_matrix[regions[ln].cells[n]].id_reg = regions[l].index;}
		//Releasing memory of regions
		MY_FREE( regions[ln].cells );
		//A bigger void region due to merging
		regions[l].Ncells += regions[ln].Ncells;
		regions[ln].Ncells = 0;}}
	//If this region is going to be reevaluated
	if( reeval == 1 )
	    l = l-1;}
  
    return 0;
}


/**************************************************************************************************
 NAME:	     void_numbers
 FUNCTION:   This function computes the number of void into certain folder
 INPUTS:     path
 RETURN:     number of voids
**************************************************************************************************/
int void_numbers( char path[] )
{
    int file_count = 0;
    DIR * dirp;
    struct dirent * entry;

    dirp = opendir( path ); /* There should be error handling after this */      
    while ((entry = readdir(dirp)) != NULL) {
	if (entry->d_type == DT_REG) { /* If the entry is a regular file */
	    file_count++;
	}
    }
    
    closedir(dirp);
    
    return file_count-4;
}


/**************************************************************************************************
 NAME:	     void_size
 FUNCTION:   This function computes the number of cell of a single void
 INPUTS:     path
 RETURN:     number of cells
**************************************************************************************************/
int void_cells( char path[] )
{
    int i=0;
    FILE *file;
    int dumb;
    
    file = fopen( path, "r" );
    
    while( getc( file ) != EOF ){
	fscanf( file, "%d %d %d", &dumb, &dumb, &dumb );
	i++;}
	
    fclose( file );
    
    return i-4;
}


/**************************************************************************************************
 NAME:	     density_centres
 FUNCTION:   Calculates the cell where the density field is minimum
 INPUTS:     Array with void regions, array of parameters, number of regions, density array
 RETURN:     0
**************************************************************************************************/
int density_centres( struct void_region regions[],
		     float *p,
		     int n_reg,
		     double *delta )
{
    int l, lcell;
    long long int n;
    double min_rho;
        
    //Sweeping void regions for detecting neighbour voids
    for( l=0; l<n_reg; l++ ){
	min_rho=1000.0;
	//Sweeping internal cells
	for( lcell=0; lcell<regions[l].Ncells; lcell++ ){
	    //Overall index of this cell
	    n = regions[l].cells[lcell];
	    //Checking for density
	    if( delta[n]<min_rho ){
		regions[l].idrhoC = n;
		min_rho = delta[n];
	    }
	}
    }
    
    return 0;
}

/**************************************************************************************************
 NAME:	     fractional_anisotropy_centres
 FUNCTION:   Calculates the cell where the density field is minimum
 INPUTS:     Array with void regions, array of parameters, number of regions, FA array
 RETURN:     0
**************************************************************************************************/
int fractional_anisotropy_centres( struct void_region regions[],
				   float *p,
				   int n_reg,
				   float *FA )
{
    int l, lcell;
    long long int n;
    float min_FA;
        
    //Sweeping void regions for detecting neighbour voids
    for( l=0; l<n_reg; l++ ){
	min_FA = 1.1;
	//Sweeping internal cells
	for( lcell=0; lcell<regions[l].Ncells; lcell++ ){
	    //Overall index of this cell
	    n = regions[l].cells[lcell];
	    //Checking for density
	    if( FA[n]<min_FA ){
		regions[l].idrhoFA = n;
		min_FA = FA[n];
	    }
	}
    }
    
    return 0;
}