#ifdef __cplusplus
extern "C" {
#endif

/* read README */

/*****************************************************************
  Variables of the basic features of CA. If you want to change
  the value, do that within Initial() function, not here.
*/
extern int nplane; /* # of planes (default=0) cash2-student.c*/
extern int ndispplane;	/* # of planes that should be displayed */
extern int nrow; /* # of row (default=100) cash */
extern int ncol; /* # of column (default=100) cash */
extern int scale; /* size of the window (default=2) cash */
extern int margin; /* margin between planes on diplay.*/
extern TYPE2 boundaryvalue2; /* the value of boundary (default=0) cash2*/
extern int boundary; /* the time of boundary: FIXED, WRAP, ECHO
			(default=WRAP). Note that Margolus
			diffusion is not supported for
			ECHO. cash */
extern unsigned long ulseedG; /* random seed (default=56) cash2-student.c */
extern int display; /* if 1, display on. if 0, off */

/* One should not modify the values of the following variables
   during the simulation, but one can utilize them. */
extern int Time; /* the current time step. cash2-student.c */
extern int MaxTime;

extern long long genrandreal1_offset;
extern long long genrandreal2_offset;

/*****************************************
 * function prototype of cash2-student.c *
 *****************************************/
void InitXmgrace(void);
void InitialMain(char*);
int Display(TYPE2**,...);
int DrawSlide(TYPE2**,char *);
//int InitGradationColor(int,int,int);/* ncolor, min, max*/
//int UpdateGradationColor(TYPE2 **,double,double,int); /* plane, min_float, max_float, fval_index*/
void MDiffusion(TYPE2**);
void Plot(int,...);
void PlotArray(double []);
void PlotXY(double x, double y);
void SavePlot(char *);
//void Asynchronous(void);
//void Synchronous(int, ...);

int GetNeighbor(TYPE2**,int,int,int);
int RandomMoore8(TYPE2**,int,int);
int RandomMoore9(TYPE2**,int,int);
int RandomNeumann4(TYPE2**,int,int);
int RandomNeumann5(TYPE2**,int,int);

int countGlobal(TYPE2 **, int);

TYPE2 GetNeighborS(TYPE2**,int,int,int);
TYPE2 RandomMooreS8(TYPE2**,int,int);
TYPE2 RandomMooreS9(TYPE2**,int,int);
TYPE2 RandomNeumannS4(TYPE2**,int,int);
TYPE2 RandomNeumannS5(TYPE2**,int,int);

int CountMoore8(TYPE2**,int,int,int);
int CountMoore9(TYPE2**,int,int,int);
int CountNeumann4(TYPE2**,int,int,int);
int CountNeumann5(TYPE2**,int,int,int);

int SumMoore8(TYPE2**,int,int);
int SumMoore9(TYPE2**,int,int);
int SumNeumann4(TYPE2**,int,int);
int SumNeumann5(TYPE2**,int,int);

void GetNeighborC(TYPE2**,int,int,int,int*,int*);
void RandomMooreC8(TYPE2**,int,int,int*,int*);
void RandomMooreC9(TYPE2**,int,int,int*,int*);
void RandomNeumannC4(TYPE2**,int,int,int*,int*);
void RandomNeumannC5(TYPE2**,int,int,int*,int*);

TYPE2* GetNeighborP(TYPE2**,int,int,int);
TYPE2* RandomMooreP8(TYPE2**,int,int);
TYPE2* RandomMooreP9(TYPE2**,int,int);
TYPE2* RandomNeumannP4(TYPE2**,int,int);
TYPE2* RandomNeumannP5(TYPE2**,int,int);


void MakePlane(TYPE2***,...);
void SavePlane(char*,TYPE2**);
void ReadSavedData(char*,TYPE2**);

void SpaceTimePlot(TYPE2**,TYPE2**);

int InitialSet(TYPE2**,int,int,...);
int InitialSetS(TYPE2**,int,TYPE2,...);

/* Margolus diffusion with obstacles */
void rotate_square( TYPE2 *a, TYPE2 *b, TYPE2 *c, TYPE2 *d );
void rotate_triangle( TYPE2 *a, TYPE2 *b, TYPE2 *c );
void rotate_pair( TYPE2 *a, TYPE2 *b );
void MargolusWithObstacle( TYPE2** d_plane, TYPE2** o_plane, int obstacle, MARGOLUS*** margolus, int phase );
void ObstacleMargolus( TYPE2 **a, TYPE2 **o_plane, int obstacle );

/* functions added by Hilje */
//void MargolusWithFixedPoints( TYPE2** d_plane, int* obstacle, MARGOLUS*** margolus, int phase );
//void FixedPointsMargolus( TYPE2 **a, int* obstacle );
void DiffusionBySwap(double,TYPE2**);
int MaxIntArray(int*, int);
double Max2Doubles(double,double);
int Min2Ints(int,int);
double ExpDistr(double);
int ExpWaitingTime(int);
TYPE **NewP3(void);
TYPE **New3(void);
int DrawSlide2(char*, int, TYPE2 **,...);


/****************************************************
 * function prototype of model.c (made by students) *
 ****************************************************/
void Initial(char*);
void InitialPlane(void);
void NextState(TYPE2**,TYPE2**,int,int);
void Update(void);
//Functions added by Hilje
int ColourInf(double);
void KillInd(TYPE2**,int,int);
double MutateFloat(double,double);
void ReadParameters(char*);

#ifdef __cplusplus
 }
#endif
