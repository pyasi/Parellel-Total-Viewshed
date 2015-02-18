/* multiView.h */

#ifndef __llist_h
#define __llist_h

//Grid struct for elevation and viewshed grids
typedef struct _grid {
	int rows, cols, xllcorner, yllcorner, cellsize, NODATA_value;
	float** data;
} Grid;

//Initialize grids
Grid* initGridArray();

// Reads header informtation from file into grid
Grid* readHeader(Grid* g, char* filename);

// Reads header and  data from file into grid
Grid* readAll(Grid* g, char* filename);
	
//Writes grid to file
void write(Grid* g, char* filename);

//Print header of grid
void printHeader(Grid* g);

//Print values of grid
void printValues(Grid* g);

//Print header and values of grid
void printGrid(Grid* g);

//Calculates vertical slope from (axInt,ayInt) to (bxInt, byInt)
float calculateSlope(Grid* g, float ax, float ay, float bx, float by, int crossX);
	
//Returns true if (bx, by) is visible from (ax, ay)
float isVisible(Grid* g, int axInt, int ayInt, int bxInt, int byInt);

//Returns the viewshed grid from point (x,y)
float computeViewshed(Grid* egrid, int x, int y);

float* multiView(char* elevfname, int start, int end);
 
#endif
