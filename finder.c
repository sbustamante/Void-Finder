#include <allvars.h>

/**************************************************************************************************
 NAME:	     matrix_voids
 FUNCTION:   This function builds the void matrix by using the FA value and the density field
 INPUTS:     Structure void_cell with empty voids, FA array, density array, eigenvalue L1 array,
	     current cut FA value, array of parameters.
 RETURN:     0
**************************************************************************************************/
int matrix_voids( struct void_cell *void_matrix, 
		  float *FA, 
		  double *delta,
		  float *eig1,
		  float FA_i, 
		  float *p )
{
    int i,j,k,N;
    long long int n;

    //Side size of simulation
    N = (int)p[NBOX];
    
    //FOF SCHEME
    if( p[SCHM] < 1 )
    //Sweeping the grid    
    for( i=0;i<N;i++ )
    for( j=0;j<N;j++ )
    for( k=0;k<N;k++ ){
	//Overall index
	n = k + N*(j + N*i);
	//Assigning position
	void_matrix[n].i = i;
	void_matrix[n].j = j;
	void_matrix[n].k = k;
	//Initializing void region
	void_matrix[n].id_reg = 0;
	//Initializing check variable
	void_matrix[n].check = 0;
	
	//Detecting number of initial void cells
	if( FA[n] <= FA_i )
	    void_matrix[n].isvoid = 1;
	else
	    void_matrix[n].isvoid = 0;}
	    
    //FA WATERSHED SCHEME
    if( p[SCHM] == 1 )
    //Sweeping the grid    
    for( i=0;i<N;i++ )
    for( j=0;j<N;j++ )
    for( k=0;k<N;k++ ){
	//Overall index
	n = k + N*(j + N*i);
	//Assigning position
	void_matrix[n].i = i;
	void_matrix[n].j = j;
	void_matrix[n].k = k;
	//Initializing void region
	void_matrix[n].id_reg = 0;
	//Initializing check variable
	void_matrix[n].check = 0;
	
	//Detecting number of initial void cells
	if( FA[n] <= p[FACT] && eig1[n] <= p[L1CT] )
	    void_matrix[n].isvoid = 1;
	else
	    void_matrix[n].isvoid = 0;}
	    
    //DENSITY WATERSHED SCHEME
    if( p[SCHM] == 2 )
    //Sweeping the grid    
    for( i=0;i<N;i++ )
    for( j=0;j<N;j++ )
    for( k=0;k<N;k++ ){
	//Overall index
	n = k + N*(j + N*i);
	//Assigning position
	void_matrix[n].i = i;
	void_matrix[n].j = j;
	void_matrix[n].k = k;
	//Initializing void region
	void_matrix[n].id_reg = 0;
	//Initializing check variable
	void_matrix[n].check = 0;
	
	//Detecting number of initial void cells
	if( delta[n] <= p[DLCT] )
	    void_matrix[n].isvoid = 1;
	else
	    void_matrix[n].isvoid = 0;}
	
    return 0;
}


/**************************************************************************************************
 NAME:	     voids_finder_FA
 FUNCTION:   Find the voids regions in the simulation by using the FA watershed scheme
 INPUTS:     Matrix with voids (1-void, 0-no void), array with FA, array with L1 eigenvalues, array
	     of parameters.
 RETURN:     0
**************************************************************************************************/
int voids_finder_FA( struct void_cell *void_matrix,
		     float *FA,
		     float *p )
{
    int i,j,k,l;	//Original counters
    int ic,jc,kc;	//Neighbour counters
    int it,jt,kt;	//Corrected neighbor counters
    int im,jm,km;	//Counters of the next cell in the path toward the local minima of FA
    int ib,jb,kb;	//Backup counters
    
    int N, b;
    long int n_reg = 0, n_reg_tot = 0, id_reg_n, id_reg_nt;
    long long int n, nt, nm, Nvoids = 0, ntmp;
    float FAm;
    
    //Regions
    struct void_region *temp_reg;
    temp_reg = NULL;
        
    //Size of the matrix
    N = p[NBOX];
    //Linking lenght of FOF scheme
    b = p[LINK];
    
    //Sweeping the grid    
    for( i=0;i<N;i++ )
    for( j=0;j<N;j++ )
    for( k=0;k<N;k++ ){
	//Making a backup of the current indexes of the cycle
	ib = i; jb = j; kb = k;
      
	//Sweeping a way throughout the minimal FA path until the local minima
	while( 1 ){
	    //Overall index
	    n = k + N*(j + N*i);

	    //Only if this cell is a void
	    if( void_matrix[n].isvoid == 1 && void_matrix[n].check == 0 ){
		//Setting this cell as checked
		void_matrix[n].check = 1;
		//Counters for numbers of voids
		Nvoids ++;
		//Assigning the current region ID if cell is alone
		if( void_matrix[n].id_reg == 0 ){
		    //Id to current zone
		    n_reg ++;
		    n_reg_tot ++;
		    //Initializing the pointer to zones
		    void_matrix[n].id_reg = n_reg;
		    //Initializing this zone
		    temp_reg = \
		    (struct void_region *)realloc( temp_reg, (n_reg+1)*sizeof( struct void_region ) );
		    //Including the current cell in this zone
		    temp_reg[n_reg].Ncells = 1;
		    temp_reg[n_reg].cells = \
		    (long long int *)calloc( temp_reg[n_reg].Ncells, sizeof( long long int ) );
		    temp_reg[n_reg].cells[ temp_reg[n_reg].Ncells-1 ] = n;}
		    
		//Detecting the way to the local minima of FA
		nm = n;
		im = i; jm = j; km = k;
		FAm = FA[n];
		for( ic=-b; ic<=b; ic++ )
		for( jc=-b; jc<=b; jc++ )
		for( kc=-b; kc<=b; kc++ ){
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
		    
		    if( FA[nt] <= FAm && void_matrix[nt].isvoid == 1 ){
			nm = nt;
			im = it; jm = jt; km = kt;
			FAm = FA[nt];
		    }}

		//Detecting a local minima of FA
		if( nm == n )
		    break;
		
		//Next cell in the minimal path
		i = im; j = jm; k = km;

		id_reg_n = void_matrix[n].id_reg;
		id_reg_nt = void_matrix[nm].id_reg;
		
		//1 CASE: If the neighbor cells does not have assigned a void region yet
		if( id_reg_nt == 0 ){
		    void_matrix[nm].id_reg = id_reg_n;
		    temp_reg[id_reg_n].Ncells ++;
		    temp_reg[id_reg_n].cells = \
		    (long long int *)realloc(temp_reg[ id_reg_n].cells, \
		    temp_reg[id_reg_n].Ncells*sizeof( long long int ) );
		    temp_reg[id_reg_n].cells[ temp_reg[id_reg_n].Ncells-1 ] = nm;}
		else{
		  
		//2 CASE: If the neighbor cells have assigned another void region but this have 
		//a major id
		if( id_reg_n < id_reg_nt && 
		    id_reg_n !=0 && id_reg_nt !=0 ){
		    //Deleting a region
		    n_reg_tot --;
		    //Mixing the two regions in the id_reg_n one
		    temp_reg[id_reg_n].cells = \
		    (long long int *)realloc(temp_reg[id_reg_n].cells, \
		    (temp_reg[id_reg_n].Ncells + temp_reg[id_reg_nt].Ncells)*\
		    sizeof( long long int ) );
		    
		    //Sweeping index of nm region
		    for( ntmp=0; ntmp<temp_reg[id_reg_nt].Ncells; ntmp++ ){
			//Assigning the id of n region to those of the nm region
			void_matrix[ temp_reg[id_reg_nt].cells[ntmp] ].id_reg = id_reg_n;
			//Storing all the cell ids of n region in the nm region
			temp_reg[id_reg_n].cells[ temp_reg[id_reg_n].Ncells + ntmp ] = \
			temp_reg[id_reg_nt].cells[ntmp];}
		    //New size of region
		    temp_reg[id_reg_n].Ncells += temp_reg[id_reg_nt].Ncells;
		    //Cleaning the nm region
		    temp_reg[id_reg_nt].Ncells = 0;
		    MY_FREE( temp_reg[id_reg_nt].cells );}
		    
		//3 CASE: If the neighbor cells have assigned another void region and this have 
		//a minor id
		if( id_reg_nt < id_reg_n && 
		    id_reg_nt !=0 && id_reg_n !=0 ){
		    //Deleting a region
		    n_reg_tot --;
		    //Mixing the two regions in the id_reg_n one
		    temp_reg[id_reg_nt].cells = \
		    (long long int *)realloc(temp_reg[id_reg_nt].cells, \
		    (temp_reg[id_reg_nt].Ncells + temp_reg[id_reg_n].Ncells)*\
		    sizeof( long long int ) );
		    
		    //Sweeping index of n region
		    for( ntmp=0; ntmp<temp_reg[id_reg_n].Ncells; ntmp++ ){
			//Assigning the id of n region to those of the nm region
			void_matrix[ temp_reg[id_reg_n].cells[ntmp] ].id_reg = id_reg_nt;
			//Storing all the cell ids of nm region in the n region
			temp_reg[id_reg_nt].cells[ temp_reg[id_reg_nt].Ncells + ntmp ] = \
			temp_reg[id_reg_n].cells[ntmp];}
		    //New size of region
		    temp_reg[id_reg_nt].Ncells += temp_reg[id_reg_n].Ncells;
		    //Cleaning the nm region
		    temp_reg[id_reg_n].Ncells = 0;
		    MY_FREE( temp_reg[id_reg_n].cells );}
	    }}
	    //If this cell is not a void or is already checked
	    else
		break;
	}
	//Recovering indexes of the cycle
	i = ib; j = jb; k = kb;}
 
     //Releasing memory
    MY_FREE(temp_reg);
 
    printf( "  * Number of void-type cells %lld\n", Nvoids );
    
    return n_reg;
}


/**************************************************************************************************
 NAME:	     voids_finder_DL
 FUNCTION:   Find the voids regions in the simulation by using the density watershed scheme
 INPUTS:     Matrix with voids (1-void, 0-no void), array with delta field, array with L1 
	     eigenvalues, array of parameters.
 RETURN:     0
**************************************************************************************************/
int voids_finder_DL( struct void_cell *void_matrix,
		     double *delta,
		     float *p )
{
    int i,j,k,l;	//Original counters
    int ic,jc,kc;	//Neighbour counters
    int it,jt,kt;	//Corrected neighbor counters
    int im,jm,km;	//Counters of the next cell in the path toward the local minima of Density
    int ib,jb,kb;	//Backup counters
    
    int N, b;
    long int n_reg = 0, n_reg_tot = 0, id_reg_n, id_reg_nt;
    long long int n, nt, nm, Nvoids = 0, ntmp;
    double deltam;
    
    //Regions
    struct void_region *temp_reg;
    temp_reg = NULL;
        
    //Size of the matrix
    N = p[NBOX];
    //Linking lenght of FOF scheme
    b = p[LINK];
    
    //Sweeping the grid    
    for( i=0;i<N;i++ )
    for( j=0;j<N;j++ )
    for( k=0;k<N;k++ ){
	//Making a backup of the current indexes of the cycle
	ib = i; jb = j; kb = k;
      
	//Sweeping a way throughout the minimal density path until the local minima
	while( 1 ){
	    //Overall index
	    n = k + N*(j + N*i);

	    //Only if this cell is a void
	    if( void_matrix[n].isvoid == 1 && void_matrix[n].check == 0 ){
		//Setting this cell as checked
		void_matrix[n].check = 1;
		//Counters for numbers of voids
		Nvoids ++;
		//Assigning the current region ID if cell is alone
		if( void_matrix[n].id_reg == 0 ){
		    //Id to current zone
		    n_reg ++;
		    n_reg_tot ++;
		    //Initializing the pointer to zones
		    void_matrix[n].id_reg = n_reg;
		    //Initializing this zone
		    temp_reg = \
		    (struct void_region *)realloc( temp_reg, (n_reg+1)*sizeof( struct void_region ) );
		    //Including the current cell in this zone
		    temp_reg[n_reg].Ncells = 1;
		    temp_reg[n_reg].cells = \
		    (long long int *)calloc( temp_reg[n_reg].Ncells, sizeof( long long int ) );
		    temp_reg[n_reg].cells[ temp_reg[n_reg].Ncells-1 ] = n;}
		    
		//Detecting the way to the local minima of delta
		nm = n;
		im = i; jm = j; km = k;
		deltam = delta[n];
		for( ic=-b; ic<=b; ic++ )
		for( jc=-b; jc<=b; jc++ )
		for( kc=-b; kc<=b; kc++ ){
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
		    
		    if( delta[nt] <= deltam && void_matrix[nt].isvoid == 1 ){
			nm = nt;
			im = it; jm = jt; km = kt;
			deltam = delta[nt];
		    }}

		//Detecting a local minima of delta
		if( nm == n )
		    break;
		
		//Next cell in the minimal path
		i = im; j = jm; k = km;

		id_reg_n = void_matrix[n].id_reg;
		id_reg_nt = void_matrix[nm].id_reg;
		
		//1 CASE: If the neighbor cells does not have assigned a void region yet
		if( id_reg_nt == 0 ){
		    void_matrix[nm].id_reg = id_reg_n;
		    temp_reg[id_reg_n].Ncells ++;
		    temp_reg[id_reg_n].cells = \
		    (long long int *)realloc(temp_reg[ id_reg_n].cells, \
		    temp_reg[id_reg_n].Ncells*sizeof( long long int ) );
		    temp_reg[id_reg_n].cells[ temp_reg[id_reg_n].Ncells-1 ] = nm;}
		else{
		  
		//2 CASE: If the neighbor cells have assigned another void region but this have 
		//a major id
		if( id_reg_n < id_reg_nt && 
		    id_reg_n !=0 && id_reg_nt !=0 ){
		    //Deleting a region
		    n_reg_tot --;
		    //Mixing the two regions in the id_reg_n one
		    temp_reg[id_reg_n].cells = \
		    (long long int *)realloc(temp_reg[id_reg_n].cells, \
		    (temp_reg[id_reg_n].Ncells + temp_reg[id_reg_nt].Ncells)*\
		    sizeof( long long int ) );
		    
		    //Sweeping index of nm region
		    for( ntmp=0; ntmp<temp_reg[id_reg_nt].Ncells; ntmp++ ){
			//Assigning the id of n region to those of the nm region
			void_matrix[ temp_reg[id_reg_nt].cells[ntmp] ].id_reg = id_reg_n;
			//Storing all the cell ids of n region in the nm region
			temp_reg[id_reg_n].cells[ temp_reg[id_reg_n].Ncells + ntmp ] = \
			temp_reg[id_reg_nt].cells[ntmp];}
		    //New size of region
		    temp_reg[id_reg_n].Ncells += temp_reg[id_reg_nt].Ncells;
		    //Cleaning the nm region
		    temp_reg[id_reg_nt].Ncells = 0;
		    MY_FREE( temp_reg[id_reg_nt].cells );}
		    
		//3 CASE: If the neighbor cells have assigned another void region and this have 
		//a minor id
		if( id_reg_nt < id_reg_n && 
		    id_reg_nt !=0 && id_reg_n !=0 ){
		    //Deleting a region
		    n_reg_tot --;
		    //Mixing the two regions in the id_reg_n one
		    temp_reg[id_reg_nt].cells = \
		    (long long int *)realloc(temp_reg[id_reg_nt].cells, \
		    (temp_reg[id_reg_nt].Ncells + temp_reg[id_reg_n].Ncells)*\
		    sizeof( long long int ) );
		    
		    //Sweeping index of n region
		    for( ntmp=0; ntmp<temp_reg[id_reg_n].Ncells; ntmp++ ){
			//Assigning the id of n region to those of the nm region
			void_matrix[ temp_reg[id_reg_n].cells[ntmp] ].id_reg = id_reg_nt;
			//Storing all the cell ids of nm region in the n region
			temp_reg[id_reg_nt].cells[ temp_reg[id_reg_nt].Ncells + ntmp ] = \
			temp_reg[id_reg_n].cells[ntmp];}
		    //New size of region
		    temp_reg[id_reg_nt].Ncells += temp_reg[id_reg_n].Ncells;
		    //Cleaning the nm region
		    temp_reg[id_reg_n].Ncells = 0;
		    MY_FREE( temp_reg[id_reg_n].cells );}
	    }}
	    //If this cell is not a void or is already checked
	    else
		break;
	}
	//Recovering indexes of the cycle
	i = ib; j = jb; k = kb;}
 
    //Releasing memory
    MY_FREE(temp_reg);
 
    printf( "  * Number of void-type cells %lld\n", Nvoids );
    return n_reg;
}


/**************************************************************************************************
 NAME:	     voids_finder_FOF
 FUNCTION:   Find the voids regions in the simulation through the FOF scheme
 INPUTS:     Matrix with voids (1-void, 0-no void), array with pointer to regions id, array of 
	     parameters.
 RETURN:     0
**************************************************************************************************/
int voids_finder_FOF( struct void_cell *void_matrix,
		      float *p )
{
    int i,j,k,l;
    int ic,jc,kc;
    int it,jt,kt;
    
    int N, b;
    long int n_reg = 0, id_reg_n, id_reg_nt;
    long long int n, nt, Nvoids = 0, ntmp;
    
    //Regions
    struct void_region *temp_reg;
    temp_reg = NULL;
        
    //Size of the matrix
    N = p[NBOX];
    //Linking lenght of FOF scheme
    b = p[LINK];
    
    //Sweeping the grid    
    for( i=0;i<N;i++ )
    for( j=0;j<N;j++ )
    for( k=0;k<N;k++ ){
    //Overall index
    n = k + N*(j + N*i);

    //Only if this cell is a void
    if( void_matrix[n].isvoid == 1 ){
	//Counters for numbers of voids
	Nvoids ++;
	//Assigning the current region ID if cell is alone
	if( void_matrix[n].id_reg == 0 ){
	    //Id to current zone
	    n_reg ++;
	    //Initializing the pointer to zones
	    void_matrix[n].id_reg = n_reg;
	    //Initializing this zone
	    temp_reg = \
	    (struct void_region *)realloc( temp_reg, (n_reg+1)*sizeof( struct void_region ) );
	    //Including the current cell in this zone
	    temp_reg[n_reg].Ncells = 1;
	    temp_reg[n_reg].cells = \
	    (long long int *)calloc( temp_reg[n_reg].Ncells, sizeof( long long int ) );
	    temp_reg[n_reg].cells[ temp_reg[n_reg].Ncells-1 ] = n;}
	    
	//Setting the neighborhood
	for( ic=-b; ic<=b; ic++ )
	for( jc=-b; jc<=b; jc++ )
	for( kc=-b; kc<=b; kc++ )
	if( ic!=0 || jc!=0 || kc!=0 ){
	  
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
	    if( void_matrix[nt].isvoid == 1 ){
		id_reg_n = void_matrix[n].id_reg;
		id_reg_nt = void_matrix[nt].id_reg;

		//1 CASE: If the neighbor cells does not have assigned a void region yet
		if( id_reg_nt == 0 ){
		    void_matrix[nt].id_reg = id_reg_n;
		    temp_reg[id_reg_n].Ncells ++;
		    temp_reg[id_reg_n].cells = \
		    (long long int *)realloc(temp_reg[ id_reg_n].cells, \
		    temp_reg[id_reg_n].Ncells*sizeof( long long int ) );
		    temp_reg[id_reg_n].cells[ temp_reg[id_reg_n].Ncells-1 ] = nt;}
		else{
		  
		//2 CASE: If the neighbor cells have assigned another void region but this have 
		//a major id
		if( id_reg_n < id_reg_nt && 
		    id_reg_n !=0 && id_reg_nt !=0 ){
		    //Mixing the two regions in the id_reg_n one
		    temp_reg[id_reg_n].cells = \
		    (long long int *)realloc(temp_reg[id_reg_n].cells, \
		    (temp_reg[id_reg_n].Ncells + temp_reg[id_reg_nt].Ncells)*\
		    sizeof( long long int ) );
		    
		    //Sweeping index of nt region
		    for( ntmp=0; ntmp<temp_reg[id_reg_nt].Ncells; ntmp++ ){
			//Assigning the id of n region to those of the nt region
			void_matrix[ temp_reg[id_reg_nt].cells[ntmp] ].id_reg = id_reg_n;
			//Storing all the cell ids of n region in the nt region
			temp_reg[id_reg_n].cells[ temp_reg[id_reg_n].Ncells + ntmp ] = \
			temp_reg[id_reg_nt].cells[ntmp];}
		    //New size of region
		    temp_reg[id_reg_n].Ncells += temp_reg[id_reg_nt].Ncells;
		    //Cleaning the nt region
		    temp_reg[id_reg_nt].Ncells = 0;
 		    MY_FREE( temp_reg[id_reg_nt].cells );}
 		    
		//3 CASE: If the neighbor cells have assigned another void region and this have 
		//a minor id
		if( id_reg_nt < id_reg_n && 
		    id_reg_nt !=0 && id_reg_n !=0 ){
		    //Mixing the two regions in the id_reg_n one
		    temp_reg[id_reg_nt].cells = \
		    (long long int *)realloc(temp_reg[id_reg_nt].cells, \
		    (temp_reg[id_reg_nt].Ncells + temp_reg[id_reg_n].Ncells)*\
		    sizeof( long long int ) );
		    
		    //Sweeping index of n region
		    for( ntmp=0; ntmp<temp_reg[id_reg_n].Ncells; ntmp++ ){
			//Assigning the id of n region to those of the nt region
			void_matrix[ temp_reg[id_reg_n].cells[ntmp] ].id_reg = id_reg_nt;
			//Storing all the cell ids of nt region in the n region
			temp_reg[id_reg_nt].cells[ temp_reg[id_reg_nt].Ncells + ntmp ] = \
			temp_reg[id_reg_n].cells[ntmp];}
		    //New size of region
		    temp_reg[id_reg_nt].Ncells += temp_reg[id_reg_n].Ncells;
		    //Cleaning the nt region
		    temp_reg[id_reg_n].Ncells = 0;
 		    MY_FREE( temp_reg[id_reg_n].cells );}
                }}}}}
              
    //Releasing memory
    MY_FREE(temp_reg);
              
    printf( "Number of void-type cells %lld\n", Nvoids );
    return n_reg;
}


/**************************************************************************************************
 NAME:	     regions_builder
 FUNCTION:   Build the voids regions in the simulation
 INPUTS:     Matrix with voids (1-void, 0-no void), array with empty regions, array of parameters,
	     number of regions.
 RETURN:     Number of non-void regions
**************************************************************************************************/
int regions_builder( struct void_cell *void_matrix,  
		     struct void_region *regions,
		     float *p,
		     int n_reg)
{
    int i=0,j=0,k=0;
    int check, n_reg_real=0;
    long int l,m;
    long long int N1, N2;
    long long int n=0, Nvoids = 0;
    long int reg_ind;
    long int *sort_regions;
    sort_regions = (long int *)calloc( n_reg, sizeof( long int ) );
    
    int N;
    //Size of the matrix
    N = p[NBOX];
    
    //Temporal structures for void regions
    struct void_region *temp_reg;
    temp_reg = (struct void_region *)calloc( n_reg, sizeof( struct void_region ) );
    for( l=0; l<n_reg; l++ ){
	temp_reg[l].Ncells = 0;
	MY_FREE(temp_reg[l].cells);}

    struct void_region *temp_reg_cut;
    temp_reg_cut = (struct void_region *)calloc( n_reg, sizeof( struct void_region ) );
    for( l=0; l<n_reg; l++ ){
	temp_reg_cut[l].Ncells = 0;
	temp_reg_cut[l].cells = NULL;}

    //Sweeping the void matrix
    for( i=0;i<N;i++ )
    for( j=0;j<N;j++ )
    for( k=0;k<N;k++ ){
    //Overall index
    n = k + N*(j + N*i);
    //Only if this cell is a void
    if( void_matrix[n].isvoid != 0 ){
	reg_ind = void_matrix[n].id_reg;
	temp_reg[reg_ind].Ncells ++;
	temp_reg[reg_ind].cells = 
	(long long int *)realloc( temp_reg[reg_ind].cells, \
	temp_reg[reg_ind].Ncells*sizeof( long long int ) );
	temp_reg[reg_ind].cells[temp_reg[reg_ind].Ncells - 1] = n;
	temp_reg[reg_ind].check = 0;
    }}
    
    //Deleting empty regions
    for( l=0; l<n_reg; l++ )
	if( temp_reg[l].Ncells > 0 ){
	    n_reg_real ++;
	    //Assigning number of cells
	    temp_reg_cut[n_reg_real].Ncells = temp_reg[l].Ncells ;
	    //Assigning id of cells
	    temp_reg_cut[n_reg_real].cells = 
	    (long long int *)calloc( temp_reg_cut[n_reg_real].Ncells,\
	    sizeof( long long int ) );
	    for( n=0; n<temp_reg_cut[n_reg_real].Ncells; n++ ){
		temp_reg_cut[n_reg_real].cells[n] = temp_reg[l].cells[n];
		void_matrix[temp_reg_cut[n_reg_real].cells[n]].id_reg = n_reg_real;}
	    MY_FREE( temp_reg[l].cells );
	}
    
    //Sorting regions according to their size (volume)
    sort_regions[0] = 0;
    N1 = 0; N2 = N*N*N;
    for( l=1; l<n_reg_real; l++ ){
	sort_regions[l] = 0;
	check = -1;
	for( m=1; m<n_reg_real; m++ ){
	    if( temp_reg_cut[m].check==0 && \
		N1<temp_reg_cut[m].Ncells && temp_reg_cut[m].Ncells<=N2 ){
		sort_regions[l] = m;
		N1 = temp_reg_cut[m].Ncells;
		temp_reg_cut[m].check = 1;
		if( check != -1 )
		    temp_reg_cut[check].check = 0;
		check = m;}}
		
	N2 = N1; N1 = 0;
	if( temp_reg_cut[sort_regions[l]].Ncells >= p[MAXV]*p[MAXV]*p[MAXV] ){
	    regions[l].Nchildren = 0;
	    regions[l].index = l;
	    regions[l].children = NULL;
	    regions[l].Ncells = temp_reg_cut[sort_regions[l]].Ncells;
	    Nvoids += regions[l].Ncells;
	    regions[l].cells = (long long int *)calloc( regions[l].Ncells,\
	    sizeof( long long int ) );
	    for( n=0; n<regions[l].Ncells; n++ ){
		regions[l].cells[n] = temp_reg_cut[sort_regions[l]].cells[n];
		void_matrix[regions[l].cells[n]].id_reg = l;}
	    MY_FREE( temp_reg_cut[sort_regions[l]].cells );
	}}
	
    //Releasing memory
    MY_FREE(temp_reg_cut);
    MY_FREE(temp_reg);
    MY_FREE(sort_regions);

    return n_reg_real+1;
}


/**************************************************************************************************
 NAME:	     hierarchical_structure
 FUNCTION:   Build the hierarchical structure of each step of the tree
 INPUTS:     Matrix with voids (1-void, 0-no void), array with current regions, array with previous
	     regions, number of current regions, number of previous regions, parameters float array.
 RETURN:     0
**************************************************************************************************/
int hierarchical_structure( struct void_cell *void_matrix,
			    struct void_region *regions,
			    struct void_region *regions_backup,
			    int n_reg,
			    int n_reg_b,
			    float *p,
			    char filename_FA[] )
{
    int i_b, i_c, i_reg;
    long long int n;
    FILE *hierarchy;
    char voidname[NMAX1];
    
    //Sweeping all previous regions in order to find current parent
    for( i_b=1; i_b<n_reg_b; i_b++ ){
	//Detecting current region of the first cell that belonged to the previous void region
	n = regions_backup[i_b].cells[0];
	i_c = void_matrix[n].id_reg;
	
	//Assigning the first children to the parent cell
	if( regions[i_c].Nchildren < 1 ){
	    regions[i_c].Nchildren = 1;
	    regions[i_c].children = (int *)calloc( regions[i_c].Nchildren, sizeof( int ) );
	    regions[i_c].children[ 0 ] = i_b;
	    }
	//Detecting major mergers of voids
	else{
	    if( regions_backup[i_b].Ncells >= p[MNVL]*regions_backup[regions[i_c].children[0]].Ncells ){
		regions[i_c].Nchildren ++;
		regions[i_c].children = (int *)realloc(regions[i_c].children, \
		regions[i_c].Nchildren*sizeof( int ) );
		regions[i_c].children[ regions[i_c].Nchildren-1 ] = i_b;
		}
	    }}
	    
    //Printing information of hierarchies of each current void region
    for( i_c=1; i_c<=n_reg; i_c++ ){
	//Creating new file
	sprintf( voidname, "%s/void_%d_child.dat", filename_FA , i_c );
	hierarchy = fopen( voidname, "w" );
	//Head
	fprintf( hierarchy, "#Current void: %d\t\tNumber of Cells = %ld\n", i_c, regions[i_c].Ncells );
	fprintf( hierarchy, "#IdVoid\tVolume\tVolume/Volume0\n");
	for( i_reg=0; i_reg<regions[i_c].Nchildren; i_reg++ ){
	    i_b = regions[i_c].children[i_reg];
	    fprintf( hierarchy, "%d\t\t%ld\t\t%1.5e\n", i_b, regions_backup[i_b].Ncells,\
	    regions_backup[i_b].Ncells/(1.0*regions_backup[ regions[i_c].children[0] ].Ncells) );}
	fclose( hierarchy );}
    
    return 0;
}