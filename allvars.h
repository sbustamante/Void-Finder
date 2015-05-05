/**************************************************************************************************
			      MACROS
**************************************************************************************************/
//General Macros
#define MAXREGION	1010000
#define NMAX1		1000
#define MY_FREE(ptr)	free(ptr); ptr = NULL;

//Macros for Parameters
#define LBOX		0	//Comoving length of the simulation
#define LINK		1	//Linking lenght for the FOF method
#define MAXV		2	//Minimum size to consider a void region [cells^3]
#define LNNH		3	//Size of neighbourhood for median filtering [cell]
#define ITER		4	//Number of successive iterations 
#define BNDR		5	//Applyig boundary removals
#define SCHM		6	//Classification schemes for voids (FOF-0, FA gradient-1)

#define FAMN		7	//Minimum FA index for classifying voids
#define FAMX		8	//Maximum FA index for classifying voids
#define NPFA		9	//Number of hierarchies of voids
#define MNVL		10	//Minimum volume percentage to consider a merger of voids

#define FACT		11	//Cut of the FA index in order to consider a void cell
#define L1CT		12	//Cut of the Lambda_1 eigenvalue to consider a void cell
#define FAMG		13	//FA threshold for merging voids

#define DLCT		14	//Density threshold value for density watershed scheme
#define DLMG		15	//Density threshold for merging voids

#define NREG		16	//Number of found void regions
#define PREG		17	//Print detailed information of each void
#define NBOX		18	//Number of cells per side of the simulation box


/**************************************************************************************************
			      STRUCTURES
**************************************************************************************************/
struct void_region{
    //Number of cells in void region
    long int Ncells;
    //Cells ID of this void region
    long long int *cells;
    
    //Number of children voids
    int Nchildren;
    //Children voids ID
    int *children;
    
    //Checked
    int check;
    //Index of this void
    long int index;
    
    //Number of neighbour voids
    int Nneighbours;
    //Indexes of neighbours
    long int *neighbours;
    //Number of shared cells with neighbours
    int *Nneigh_cells;
    //Mean FA(density) across boundaries for each neighbour
    float *neigh_mean;
    
    //ID of the cell corresponding with the Density centre of the void
    int idrhoC;
    //ID of the cell corresponding with the FA centre of the void
    int idrhoFA;
    };
    
struct void_cell{
    //Id of this cell
    int i, j, k;
    //Is a void?
    int isvoid;
    //Id of whole region
    int id_reg;
    //Checked region
    int check;
    };

    
/**************************************************************************************************
			      HEADERS
**************************************************************************************************/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <sys/stat.h>
#include <dirent.h>

#include <proto.h>