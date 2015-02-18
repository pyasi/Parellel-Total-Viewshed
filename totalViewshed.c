/* viewshed.c */
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include "totalViewshed.h"
#include "mpi.h"

//Global variables for MPI
static int OVERLORD_MINION_TAG = 1;
static int MINION_SCRIBE_TAG = 2;
static int MINIONSLACKER_SCRIBE_TAG;
int myId, numProcs, nameLength;
char myName[MPI_MAX_PROCESSOR_NAME];

//Returns the size of the grid
int getGridSize(char* filename){

        FILE* f;
        char s[100];
        int cols, rows;

        f = fopen(filename, "r");
        if (f == NULL) {
                printf("cannot open file..");
                exit(1);
        }

        fscanf(f, "%s", s);
        fscanf(f, "%d", &cols);

        fscanf(f, "%s", s);
        fscanf(f, "%d", &rows);

        return rows * cols;
}





//Overlord function: Overlord breaks up the DEM by sending out indicies from 0 to size of grid evenly among
//the minions. 
void overlord(char* gridName) {

	int gridSize = getGridSize(gridName);
	int sectionSize = gridSize / (numProcs - 2);
	int array[2];

	int j=1;
	int i;
	//Sends the indicies to the minion
	for(i = 0; i < numProcs - 3; i++){
		array[0] = sectionSize * i;
		array[1] = array[0] + sectionSize - 1;
		MPI_Send(array, 2, MPI_INT, j, OVERLORD_MINION_TAG, MPI_COMM_WORLD);
		j++;
	}
	array[0] = sectionSize* i;
	array[1] = gridSize - 1;

	//Sends the indicies of the slacker to the minion
	MPI_Send(array, 2, MPI_INT, j, OVERLORD_MINION_TAG, MPI_COMM_WORLD);
}

//Minion: The minons are the parallel portion of the computation.  Each minion computes the total viewshed of
//all the points passed to it from the overlord.  It puts the computations into the array and sends the array
//to the scribe.
void minion(char* gridName) {

	int array[2];
	MPI_Status status;

	//Recieves the indicies from the overloard
	MPI_Recv(array, 2, MPI_INT, 0, OVERLORD_MINION_TAG, MPI_COMM_WORLD, &status);

	int computedSectionSize = array[1] - array[0] + 2;
	float *computedSection = (float *)malloc((sizeof(float) * (computedSectionSize)));
	computedSection = multiView(gridName, array[0], array[1]);
	
	//Sends the arrays to the scribe
	MPI_Send(computedSection, computedSectionSize, MPI_FLOAT, numProcs - 1, MINION_SCRIBE_TAG, MPI_COMM_WORLD);
}

//Scribe: The scribe collects the arrays from all of the minions and puts them into a grid.  Once all of the 
//arrays have been collected, scribe writes the grid to a file
void scribe(char* gridName, char* outputName) {

	Grid* multiViewGrid;
	multiViewGrid = readHeader(multiViewGrid, gridName);

	int gridSize = getGridSize(gridName);
	int computedSectionSize = (gridSize / (numProcs - 2)) + 2;
	float *computedSection = (float *)malloc((sizeof(float)*computedSectionSize));
	MPI_Status status;

	int i,j;
	for(i = 1; i < numProcs - 2; i++) {
		
		//Recieving arrays from all of the full minions
		MPI_Recv(computedSection, computedSectionSize, MPI_FLOAT, i, MINION_SCRIBE_TAG, MPI_COMM_WORLD, &status);

		int k;
		int iterator = 1;
		int startRow = (int) computedSection[0] / multiViewGrid->cols;
		int startCol = (int) computedSection[0] % multiViewGrid->cols;
		
		//Putting the array into appropriate cells in the grid
		for(j = startRow; j < multiViewGrid->rows; j++){
			for(k = startCol; k < multiViewGrid->cols; k++){
				multiViewGrid->data[j][k] = computedSection[iterator];
				iterator++;
			}
		}
	}

	//Recieve array from slacker minion
	MPI_Recv(computedSection, computedSectionSize , MPI_FLOAT, i, MINION_SCRIBE_TAG, MPI_COMM_WORLD, &status);

	int k;
	int iterator = 1;
	int startRow = (int) computedSection[0] / multiViewGrid->cols;
	int startCol = (int) computedSection[0] % multiViewGrid->cols;

	//Putting the slacker array into the appropriate cells in the grid
	for(j = startRow; j < multiViewGrid->rows; j++){
		for(k = startCol; k < multiViewGrid->cols; k++){
			multiViewGrid->data[j][k] = computedSection[iterator];
			iterator++;
		}
	}

	//Writes grid to file
	write(multiViewGrid,outputName);
}


//Initializes the grid
Grid* initGrid(int c, int r, int xll, int yll, int cs, int v) {
	Grid* newGrid = (Grid *)malloc(sizeof(Grid));
	newGrid->rows = r;
	newGrid->cols = c;
	newGrid->xllcorner = xll;
	newGrid->yllcorner = yll;
	newGrid->cellsize = cs;
	newGrid->NODATA_value =  v;
	int i,j,k;
	newGrid->data = (float **)malloc(sizeof(float *) * r);
	for(i=0; i<r; i++) {
		newGrid->data[i] = (float *)malloc(sizeof(float) * c);
	}
	return newGrid;
}


//Reads header informtation from file into grid
Grid* readHeader(Grid* g, char* filename) {

	FILE* f;
	char s[100];
	int cols, rows, yllcorner, xllcorner, cellsize, NODATA_value;

	f = fopen(filename, "r");
	if (f == NULL) {
		printf("cannot open file..");
		exit(1);
	}

	fscanf(f, "%s", s);
	fscanf(f, "%d", &cols);

	fscanf(f, "%s", s);
	fscanf(f, "%d", &rows);

	fscanf(f, "%s", s);
	fscanf(f, "%d", &yllcorner);

	fscanf(f, "%s", s);
	fscanf(f, "%d", &xllcorner);

	fscanf(f, "%s", s);
	fscanf(f, "%d", &cellsize);

	fscanf(f, "%s", s);
	fscanf(f, "%d", &NODATA_value);

	fclose(f);	

	g = initGrid(cols, rows, yllcorner, xllcorner, cellsize, NODATA_value);
	return g;
}

//Reads header and data from file into grid
Grid* readAll(Grid* g, char* filename) {

	FILE* f;
	char s[100];
	int cols, rows, yllcorner, xllcorner, cellsize, NODATA_value;

	f = fopen(filename, "r");
	if (f == NULL) {
		printf("cannot open file..");
		exit(1);
	}

	fscanf(f, "%s", s);
	fscanf(f, "%d", &cols);

	fscanf(f, "%s", s);
	fscanf(f, "%d", &rows);

	fscanf(f, "%s", s);
	fscanf(f, "%d", &yllcorner);

	fscanf(f, "%s", s);
	fscanf(f, "%d", &xllcorner);

	fscanf(f, "%s", s);
	fscanf(f, "%d", &cellsize);

	fscanf(f, "%s", s);
	fscanf(f, "%d", &NODATA_value);

	g = initGrid(cols, rows, yllcorner, xllcorner, cellsize, NODATA_value);

	float entry;
	int r = 0;
	int c = 0;
	fscanf(f, "%f", &entry);
	while(!feof(f)){
		g->data[r][c] = entry;
		c++;
		if (c >= g->cols) {
			r++;
			c=0;
		}
		fscanf(f, "%f", &entry);
	}

	fclose(f);

	return g;
}

//Writes grid to file
void write(Grid* g, char* filename) {

	FILE* f;	
	char s[100];
	int cols, rows, yllcorner, xllcorner, cellsize, NODATA_value;

	f = fopen(filename, "w");	
	if (f == NULL) {
		printf("cannot open file..");
		exit(1);
	}

	fprintf(f, "ncols         %d\n", g->cols);
	fprintf(f, "nrows         %d\n", g->rows);
	fprintf(f, "xllcorner     %d\n", g->xllcorner);
	fprintf(f, "yllcorner     %d\n", g->yllcorner);
	fprintf(f, "cellsize      %d\n", g->cellsize);
	fprintf(f, "NODATA_value  %d\n", g->NODATA_value);

	int i,j;
	for(i=0; i<g->rows; i++) {
		for(j=0; j<g->cols; j++){
			fprintf(f, "%f ", g->data[i][j]);
		}
		fprintf(f,"\n");
	}
}

//Print header of grid
void printHeader(Grid* g) {
	printf("ncols         %d\n", g->cols);
	printf("nrows         %d\n", g->rows);
	printf("xllcorner     %d\n", g->xllcorner);
	printf("yllcorner     %d\n", g->yllcorner);
	printf("cellsize      %d\n", g->cellsize);
	printf("NODATA_value  %d\n", g->NODATA_value);

}

//Print values of grid
void printValues(Grid* g) {
	int i,j;
	for(i=0; i<g->rows; i++) {
		for(j=0; j<g->cols; j++){
			printf("%f ", g->data[i][j]);
		}
		printf("\n");
	}
}

//Print header and values of grid
void printGrid(Grid* g) {
	printHeader(g);
	printValues(g);
}


//Calculates vertical slope from (axInt,ayInt) to (bxInt, byInt)
float calculateSlope(Grid* g, float ax, float ay, float bx, float by, int crossX) {

	//Variables for calculations
	float heightA, heightB, difference, distance, roundUp, slope, length;
	heightA = g->data[(int)ax][(int)ay];
	//If intersections are with vertical gridlines
	if(crossX == 1) {
		length = by - floor(by);
		//Intersection is at a gridline intersection
		if(length == 0) {
			heightB = g->data[(int)bx][(int)by];
		}
		//Intersection is between gridline intersections
		else {
			roundUp = by + 1.0;
			slope = g->data[(int)bx][(int)roundUp] - g->data[(int)bx][(int)by];
			heightB = slope * length + g->data[(int)bx][(int)by];
		}
	}
	//If intersections are with horizontal gridlines
	else {
		length = bx - floor(bx);
		//Intersection is at gridline intersection
		if(length == 0) {
			heightB = g->data[(int)bx][(int)by];
		}
		//Intersection is between gridline intersections
		else {
			roundUp = bx + 1.0;
			slope = g->data[(int)roundUp][(int)by] - g->data[(int)bx][(int)by];
			heightB = slope * length + g->data[(int)bx][(int)by];
		}
	}
	difference = heightB - heightA; 
	distance = sqrt((by-ay) * (by-ay) + (bx-ax) * (bx-ax));
	return difference/distance;
}

//Returns true if (bx, by) is visible from (ax, ay)
float isVisible(Grid* g, int axInt, int ayInt, int bxInt, int byInt) {

	//Variables for calculations
	float ax, ay, bx, by, m, b, x, y;
	ax = (float)axInt;
	ay = (float)ayInt;
	bx = (float)bxInt;
	by = (float)byInt;

	//If (bx,by) is (ax,ay)	
	if(axInt == bxInt && ayInt == byInt) {
		return 1.0;
	}

	//Slope from (ax,ay) to (bx,by)
	float mainSlope = calculateSlope(g, axInt, ayInt, bxInt, byInt, 1);

	//If the slope is infinite
	if(bx-ax == 0) {
		y = ay; //y iterator
		if(y<by) {
			y = y + 1.0;
		}
		if(y>by){
			y = y - 1.0;
		}

		while(y!=by) {
			if(calculateSlope(g, ax, ay, bx, y, 0) > mainSlope) {
				return 0.0;
			}
			if(y<by) {
				y = y + 1.0;
			}else{
				y = y - 1.0;
			}
		}
	}
	//If the slope is zero
	else if(by-ay == 0) {
		x = ax; //y iterator
		if(x<bx) {
			x = x + 1.0;
		}
		if(x>bx){
			x = x - 1.0;
		}

		while(x!=bx) {
			if(calculateSlope(g, ax, ay, x, by, 1) > mainSlope) {
				return 0.0;
			}
			if(x<bx) {
				x = x + 1.0;
			}else{
				x = x - 1.0;
			}
		}
	} 
	//If the slope is finite and nonzero
	else {

		m = (by-ay)/(bx-ax); //slope of line
		b = ay - (m*ax); //y-intercept
		x = ax; //x iterator

		if(x<bx) {
			x = x + 1.0;
		}else{
			x = x - 1.0;
		}
		while(x != bx) {
			if(calculateSlope(g, ax, ay, x, m*x+b, 1) > mainSlope) {
				return 0.0;
			}

			if(x<bx) {
				x = x + 1.0;
			}else{
				x = x - 1.0;
			}
		}

		y = ay; //y iterator
		if(y<by) {
			y = y + 1.0;
		}else{
			y = y - 1.0;
		}

		while(y != by) {
			if(calculateSlope(g, ax, ay, (y-b)/m, y, 0) > mainSlope) {
				return 0.0;
			}

			if(y<by) {
				y = y + 1.0;
			}else{
				y = y - 1.0;
			}
		}
	}

	return 1.0;
}

//Returns the viewshed grid from point (x,y)
float computeViewshed(Grid* egrid, int x, int y) {
	int i, j;
	float visible = 0;
	for(i=0; i<egrid->rows; i++) {
		for(j=0; j<egrid->cols; j++) {
			if((float)isVisible(egrid, x, y, i, j) == 1) {
				visible = visible + 1.0;
			}
		}
	}
	return visible;
}

//MultiView returns an array of the computed totalviewshed of the indicies of the grid
//from start to end.
float* multiView(char* elevfname, int start, int end) {

	Grid* elevgrid;
	elevgrid = readAll(elevgrid, elevfname);

	float *computedSection = (float *)malloc(sizeof(float)*(end-start));
	computedSection[0] = (float)start;

	int  startRow, startCol;
	startRow = start / elevgrid->cols;
	startCol = start % elevgrid->cols;

	int row, col;        
	int count = 1;
	int index = start;

	for(row=startRow; row<elevgrid->rows; row++) {
		for(col=startCol; col<elevgrid->cols; col++) {
			computedSection[count] = computeViewshed(elevgrid, row, col);
			if(index == end + 1){
				return computedSection;
			}
			index++;
			count++;

		}
		startCol = 0;
	}
	return computedSection;
}

//Main method
int main(int argc, char *argv[]) {

	clock_t begin, end;
	double executionTime;

	if(argc != 3) {
		printf("Enter with grid file name and output file name\n");
		exit(1);
	}

	begin = clock();
	
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myId);
	MPI_Get_processor_name(myName, &nameLength);
	
	if(myId == 0) {
		overlord(argv[1]);
	} else if(myId == numProcs - 1) {
		scribe(argv[1], argv[2]);
	}
	else {
		minion(argv[1]);
	}

	MPI_Finalize();

	end = clock();
	executionTime = (double)(end - begin) / CLOCKS_PER_SEC;
	printf("%f - Proccess ID: %d\n", executionTime, myId);

}

