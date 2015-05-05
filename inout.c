#include <allvars.h>

/**************************************************************************************************
 NAME:	     conf2dump
 FUNCTION:   To convert a data file text in plain text 
 INPUTS:     name of configuration file
 RETURN:     0
**************************************************************************************************/
int conf2dump( char filename[] )
{
    char cmd[100];
    sprintf( cmd, "grep -v \"#\" %s | grep -v \"^$\" | gawk -F\"=\" '{print $2}' > %s.dump", 
	     filename, filename );
    system( cmd );

    return 0;
}


/**************************************************************************************************
 NAME:       in2dump
 FUNCTION:   To convert a data file text in plain text 
 INPUTS:     name of configuration file
 RETURN:     0
**************************************************************************************************/
int in2dump( char filename[] )
{
    char cmd[100];
    sprintf( cmd, "grep -v \"#\" %s > %s.dump", 
	     filename, filename );
    system( cmd );

    return 0;
}


/**************************************************************************************************
 NAME:       read_parameters
 FUNCTION:   read the file with given name and load information of array given
 INPUTS:     array where it returns reading data and file name 
	     with configuration file
 RETURN:     0 if file read ok
	     1 if file dont exist
**************************************************************************************************/
int read_parameters( float parameters[],
		     char filename[] )
{
    char cmd[100], filenamedump[100];
    int i=0;
    FILE *file;

    //Load of File
    file = fopen( filename, "r" );
    if( file==NULL ){
	printf( "  * The file '%s' don't exist!\n", filename );
	return 1;}
    fclose(file);
    
    //Converting to plain text
    conf2dump( filename );
    sprintf( filenamedump, "%s.dump", filename );
    file = fopen( filenamedump, "r" );
    
    //Reading
    while( getc( file ) != EOF ){
	fscanf( file, "%f", &parameters[i] );
	i++;}
    
    fclose( file );
    
    printf( "  * The file '%s' has been loaded!\n", filename );

    sprintf( cmd, "rm -rf %s.dump", filename );
    system( cmd );
    
    return 0;
}


/**************************************************************************************************
 NAME:       read_bin32
 FUNCTION:   read a binary file with a format of float (32 bits)
 INPUTS:     input filename, empty array to store the information
 RETURN:     side size of the simulation
**************************************************************************************************/
int read_bin32( char filename[],
		float **binary )
{
    FILE *in;
    int dumb;
    char line[30];
    int n_x, n_y, n_z;
    int n_nodes;
    long long n_total;
    float dx, dy, dz, x_0, y_0, z_0;

    //Loading file
    if(!(in=fopen(filename, "r"))){
	fprintf(stderr, "Problem opening file %s\n", filename);
	exit(1);}
    printf( "  * The binary file '%s' has been loaded!\n", filename );
    fread(&dumb,sizeof(int),1,in);
    fread(line,sizeof(char)*30,1,in);
    fread(&dumb,sizeof(int),1,in);
    fread(&dumb,sizeof(int),1,in);
    fread(&n_x,sizeof(int),1,in);    
    fread(&n_y,sizeof(int),1,in);    
    fread(&n_z,sizeof(int),1,in);    
    fread(&n_nodes,sizeof(int),1,in);    
    fread(&x_0,sizeof(float),1,in);    
    fread(&y_0,sizeof(float),1,in);    
    fread(&z_0,sizeof(float),1,in);    
    fread(&dx,sizeof(float),1,in);    
    fread(&dy,sizeof(float),1,in);    
    fread(&dz,sizeof(float),1,in);    
    fread(&dumb,sizeof(int),1,in);
    n_total = n_x * n_y * n_z;

    if(!(*binary=malloc(n_nodes * sizeof(float)))){
	fprintf(stderr, "problem with array allocation\n");
	exit(1);}
    fread(&dumb,sizeof(int),1,in);
    fread(&(*binary[0]),sizeof(float), n_total, in);
    fread(&dumb,sizeof(int),1,in);
    
    fclose(in);
    
    return n_x;
}


/**************************************************************************************************
 NAME:       read_bin64
 FUNCTION:   read a binary file with a format of double (64 bits)
 INPUTS:     input filename, empty array to store the information
 RETURN:     side size of the simulation
**************************************************************************************************/
int read_bin64( char filename[],
		double **binary )
{
    FILE *in;
    int dumb;
    char line[30];
    int n_x, n_y, n_z;
    int n_nodes;
    long long n_total;
    float dx, dy, dz, x_0, y_0, z_0;

    //Loading file
    if(!(in=fopen(filename, "r"))){
	fprintf(stderr, "Problem opening file %s\n", filename);
	exit(1);}
    printf( "  * The binary file '%s' has been loaded!\n", filename );
    fread(&dumb,sizeof(int),1,in);
    fread(line,sizeof(char)*30,1,in);
    fread(&dumb,sizeof(int),1,in);
    fread(&dumb,sizeof(int),1,in);
    fread(&n_x,sizeof(int),1,in);    
    fread(&n_y,sizeof(int),1,in);    
    fread(&n_z,sizeof(int),1,in);    
    fread(&n_nodes,sizeof(int),1,in);    
    fread(&x_0,sizeof(float),1,in);    
    fread(&y_0,sizeof(float),1,in);    
    fread(&z_0,sizeof(float),1,in);    
    fread(&dx,sizeof(float),1,in);    
    fread(&dy,sizeof(float),1,in);    
    fread(&dz,sizeof(float),1,in);    
    fread(&dumb,sizeof(int),1,in);
    n_total = n_x * n_y * n_z;

    if(!(*binary=malloc(n_nodes * sizeof(double)))){
	fprintf(stderr, "problem with array allocation\n");
	exit(1);}
    fread(&dumb,sizeof(int),1,in);
    fread(&(*binary[0]),sizeof(double), n_total, in);
    fread(&dumb,sizeof(int),1,in);
    
    fclose(in);
    
    return n_x;
}


/**************************************************************************************************
 NAME:       data_out
 FUNCTION:   write all files with voids information
 INPUTS:     'void_cell' structure with matrix data of voids, 'void_region' structure with 
	     information of void regions, number of regions used, parameters and folder name 
	     Number of pairs.
 RETURN:     Number of non-empty regions
**************************************************************************************************/
int data_out( struct void_cell void_matrix[], 
	      struct void_region regions[], 
	      int n_reg,
	      float *p,
	      char filename[] )
{
    int i;
    long long int j, n;
    char voidname[NMAX1], filenamedump[NMAX1], os_command[NMAX1];
    FILE *file, *voids, *neigh, *rho_cen, *FA_cen;
    struct stat sb;

    //Detecting previous data and deleting if it is the case
    if (stat(filename, &sb) == 0 && S_ISDIR(sb.st_mode)){
	sprintf( os_command, "rm -r %s", filename );
	system( os_command );}
    //Creating folder file to store the cells catalogues of each void region
    sprintf( os_command, "mkdir %s", filename );
    system( os_command );
    
    //File with regions
    sprintf( filenamedump, "%s/void_regions.dat", filename );
    file = fopen( filenamedump, "w" );
    
    //File with density centres
    sprintf( filenamedump, "%s/DC.dat", filename );
    rho_cen = fopen( filenamedump, "w" );
    
    //File with FA centres
    sprintf( filenamedump, "%s/FAC.dat", filename );
    FA_cen = fopen( filenamedump, "w" );
    
    //Initializing the number of regions
    p[NREG] = 0;
    //Header of file
    fprintf( file, "#IdRegion\tCellsNumber\tEff_Radius\n");
    //Finding no empty regions
    for( i=0; i<n_reg; i++ ){
	if( regions[i].Ncells != 0 ){
	    p[NREG] ++;
	    //Id regions and volumes
	    fprintf(file, "%d\t\t%ld\t\t%1.6e\n",(int)p[NREG],regions[i].Ncells, 
		    effective_radius(regions[i].Ncells, p[LBOX], p[NBOX]));
	    //Density centres
	    fprintf(rho_cen, "%d\t%e\t%e\t%e\t%d\t%d\t%d\n",
		    (int)p[NREG], 
		    void_matrix[regions[i].idrhoC].i*p[LBOX]/p[NBOX],
		    void_matrix[regions[i].idrhoC].j*p[LBOX]/p[NBOX],
		    void_matrix[regions[i].idrhoC].k*p[LBOX]/p[NBOX],
		    void_matrix[regions[i].idrhoC].i,
		    void_matrix[regions[i].idrhoC].j,
		    void_matrix[regions[i].idrhoC].k);
	    //FA centres
	    fprintf(FA_cen, "%d\t%e\t%e\t%e\t%d\t%d\t%d\n",
		    (int)p[NREG], 
		    void_matrix[regions[i].idrhoFA].i*p[LBOX]/p[NBOX],
		    void_matrix[regions[i].idrhoFA].j*p[LBOX]/p[NBOX],
		    void_matrix[regions[i].idrhoFA].k*p[LBOX]/p[NBOX],
		    void_matrix[regions[i].idrhoFA].i,
		    void_matrix[regions[i].idrhoFA].j,
		    void_matrix[regions[i].idrhoFA].k);
	    
	    //Storing information of each void region
	    if( p[PREG] == 1 ){
		//Saving catalog of cells for each region
		sprintf( voidname, "%s/void_%d.dat", filename , (int)p[NREG] );
		voids = fopen( voidname, "w" );
		//Header
		fprintf( voids, "#i\tj\tk\n");
		for( j=0; j<regions[i].Ncells; j++ ){
		    n = regions[i].cells[j];
		    fprintf(voids, "%d\t%d\t%d\n",
			    void_matrix[n].i, 
			    void_matrix[n].j, 
			    void_matrix[n].k);}
		fclose(voids);
		
	        //Saving information of neighbours
		sprintf( voidname, "%s/void_%d.ngb", filename , (int)p[NREG] );
		neigh = fopen( voidname, "w" );
		//Number of neighbours
		fprintf( neigh, "%d\n", regions[i].Nneighbours );
		//ID of neighbours
		for( j=0; j<regions[i].Nneighbours; j++ ){
		    fprintf(neigh, "%ld\n", 
			    regions[i].neighbours[j]);}
		fclose(neigh);}
	  
	}}
    fclose( file );
    fclose( rho_cen );
    fclose( FA_cen );
    
    //Storing information of each void region
    if( p[PREG] == 1 ){
	//File of regions id
	sprintf( voidname, "%s/void_index.dat", filename);
	file = fopen( voidname, "w" );
	for( n=0; n<p[NBOX]*p[NBOX]*p[NBOX]; n++ ){
	    fprintf(file, "%d\n",(int)(void_matrix[n].id_reg));}
	fclose( file );}

    //Printing number of regions
    printf( "  * Number of printed bulk void regions %d\n", (int)p[NREG] );
	
    return (int)p[NREG];
}


/**************************************************************************************************
 NAME:       read_void_cell
 FUNCTION:   read cells of a given void
 INPUTS:     empty void_cell struct for storing cells, external filename
 RETURN:     0 if file read ok
	     1 if file dont exist
**************************************************************************************************/
int read_void_cell( struct void_cell void_matrix[],
		    char filename[] )
{
    int i=0;
    FILE *file;
    int dumb;

    //Load of File
    file = fopen( filename, "r" );
    if( file==NULL ){
	printf( "  * The file '%s' don't exist!\n", filename );
	return 1;}
	
    //Reading
    fscanf( file, "%d %d %d", &dumb, &dumb, &dumb );
    while( getc( file ) != EOF ){
	fscanf( file, "%d %d %d", &void_matrix[i].i, &void_matrix[i].j, &void_matrix[i].k );
	i++;}
    
    fclose( file );
    
    return 0;
}