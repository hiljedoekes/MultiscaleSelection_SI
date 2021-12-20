/*
  CASH2-student

  CASH is "Cellular Automata Simulated Hardware" created by R de
  Boer (http://theory.bio.uu.nl/rdb/). CASH2 was created as an
  extension of CASH by N Takeuchi. CASH2-student is an easy to
  use version of CASH2 prepared for the use in the course
  Bioinformatic Processes given by P Hogeweg at Utrecht
  University. Authors hold copy right. No warranty.
*/

/* Version
   0.1 change for-loop to duff's device in Asynchronous().
   0.0 the first.
*/

#include <stdio.h>
#include <string.h>
#include <time.h>
#include <errno.h>
#include <math.h>
#include <grace_np.h>
#include <unistd.h>
#include <float.h>
#include <limits.h>
#include <signal.h>
#include <stdarg.h>

#include "cash2003.h"
#include "cash2.h"
#include "mersenne.h"

/* system parameters */
#define BUFSIZE 250

/* PlaneNei structure will keep a track of planes and their
   neigbors: we do not know how many planes user makes and which
   plane user is referring when user want a neighbor. */
struct PlaneNei {
  TYPE2 **plane;
  TYPE2 ***neighbor[9];
};

int MaxTime=INT_MAX;
int Time=0;
int nplane=2;
int ndispplane=1;
unsigned long ulseedG=56;

int display=1;
static int movie_g=0;

int margin;
static MARGOLUS ***margolusG;
static struct updorder *updorderG;
static struct point **neighborG[9];
static struct PlaneNei *planeneiG;

extern int nrow,ncol,scale,boundary;
extern TYPE2 boundaryvalue2;

static TYPE **viewG;
static int poseG=0;
static int stopG=0;

//Added by Hilje: needed for DrawSlide2
static TYPE **viewG2;
//Added by Hilje: needed for exact reproduction from saveplanes
extern long long genrandreal1_offset;
extern long long genrandreal2_offset;
int margolus_phase = 0;			//Because we only use Margolus-diffusion on one plane, we keep track of the phase with a single int.
extern int initial_flag_UpdOrdReset;

/* function prototype of cash2-student.c */

void InitXmgrace(void);
void InitialMain(char*);
int Display(TYPE2**,...);
//int InitGradationColor(int,int,int);
//int UpdateGradationColor(TYPE2 **,double,double,int);
void MDiffusion(TYPE2**);
void Plot(int,...);
void PlotArray(double []);
void PlotXY(double x,double y);
//void Synchronous(int, ...);
//void Asynchronous(void);

int GetNeighbor(TYPE2**,int,int,int);
int RandomMoore8(TYPE2**,int,int);
int RandomMoore9(TYPE2**,int,int);
int RandomNeumann4(TYPE2**,int,int);
int RandomNeumann5(TYPE2**,int,int);

TYPE2 GetNeighborS(TYPE2**,int,int,int);
TYPE2 RandomMooreS8(TYPE2**,int,int);
TYPE2 RandomMooreS9(TYPE2**,int,int);
TYPE2 RandomNeumannS4(TYPE2**,int,int);
TYPE2 RandomNeumannS5(TYPE2**,int,int);

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

int gotMouse(void);
void MakePlane(TYPE2***,...);
void SavePlane(char*,TYPE2**);
void ReadSavedData(char*,TYPE2**);
void SpaceTimePlot(TYPE2**,TYPE2**);

int InitialSet(TYPE2**,int,int,...);
int InitialSetS(TYPE2**,int,TYPE2,...);


/* function prototype of model.c (made by students) */

void Initial(char*);
void InitialPlane(void);
void NextState(TYPE2**,TYPE2**,int,int);
void Update(void);

/* functions added by Hilje */
void DiffusionBySwap(double,TYPE2**);
int MaxIntArray(int*, int);
int Max2DIntArray(int**, int, int);
double Max2Doubles(double,double);
int Min2Ints(int,int);
double ExpDistr(double);
int ExpWaitingTime(int);
TYPE** NewP3(void);
TYPE** New3(void);
int DrawSlide2(char*, int, TYPE2**,...);


int main (int argc, char *argv[])
{
  char *savefolder;	/*name of folder in which results will be saved*/
  savefolder = "nores";	/*If no foldername is given, savefolder will keep this value*/

  /* if parameter file is given as argument: make filename available */
  if(argc > 1)
  {
    savefolder = argv[1];
    fprintf(stdout, "Reading parameters from and saving results in %s \n", savefolder);
  }
	
  /* initialize */
  InitialMain(savefolder);

  while(Time<=MaxTime){
    Update();
    Time++;
    
    if(display==1)
      if(Mouse()==1||poseG==1)gotMouse();

    if(stopG==1)
      break;
  }

  if(display==1)
    CloseDisplay();
  if(GraceIsOpen()&&display==0)
    GraceClose();
  if(movie_g==1)
    ClosePNG();
  return 0;
}

void InitXmgrace(void)
{
  int i;
  if(display==1){
    if (GraceOpen(16384) == -1) {
      fprintf(stderr,"I cannot run Xmgrace.\n");
      exit(-1);
    }
  }else{
    if (GraceOpenVA("gracebat",16384,"-noask","-nosafe",NULL) == -1) {
      fprintf(stderr,"I cannot run Grace. \n");
      exit(-1);
    }
  }
  /* g0 is a graph for population. */
  
  GracePrintf("g0 on");
  GracePrintf("with g0");

  GracePrintf("focus g0");

  for(i=0;i<15;i++){
    GracePrintf("s%d on",i);
    GracePrintf("s%d color %d",i,i+1);
  }
}

void InitialMain(char* savefolder)
{
  /* set the value of nrow,ncol,boundary,scale,nplane,ulseed */
  Initial(savefolder);

  init_genrand(ulseedG);

  // Definition of colors. Take care: only 15 colors defined now.
/*  ColorRGB(0,255,255,255);
  ColorRGB(1,0,75,255);
  ColorRGB(2,150,150,150);
*/

  ColorRGB(0,0,0,0);
  ColorRGB(1,0,75,255);
  ColorRGB(2,255,255,255);

  ColorRGB(3,0,0,255);		/*Infected with low infectivity: blue*/
  ColorRGB(4,21,0,255);
  ColorRGB(5,43,0,255);
  ColorRGB(6,64,0,255);
  ColorRGB(7,86,0,255);
  ColorRGB(8,107,0,255);
  ColorRGB(9,129,0,255);
  ColorRGB(10,150,0,255);
  ColorRGB(11,172,0,255);
  ColorRGB(12,193,0,255);		
  ColorRGB(13,214,0,255);		
  ColorRGB(14,236,0,255);
  ColorRGB(15,255,0,255);
  ColorRGB(16,255,0,220);
  ColorRGB(17,255,0,185);
  ColorRGB(18,255,0,150);
  ColorRGB(19,255,0,115);
  ColorRGB(20,255,0,80);
  ColorRGB(21,255,0,40);
  ColorRGB(22,255,0,0);		/*Infected with high infectivity: red*/


/*
  // Definition of colors. Take care: only 15 colors defined now.
  ColorRGB(0,0,0,0);
  ColorRGB(1,0,75,255);
  ColorRGB(2,255,255,255);

  ColorRGB(3,255,0,0);
  ColorRGB(4,0,175,200);
  ColorRGB(5,255,255,0);
  ColorRGB(6,255,180,50);
  ColorRGB(7,255,90,50);
  ColorRGB(8,0,0,200);
  ColorRGB(9,255,150,200);
  ColorRGB(10,200,0,100);

  ColorRGB(11,255, 165, 0);
  ColorRGB(12,114, 33, 188);
  ColorRGB(13,103, 7, 72);
  ColorRGB(14,64, 224, 208);
  ColorRGB(15,0, 139, 0);
  ColorRGB(16,255, 0, 255);
*/
  // Color gradient for local selection estimates: Light blue -> black, black-> light red.
  ColorRGB(69, 0, 0, 0); 
  ColorRGB(68, 5, 12, 13); 
  ColorRGB(67, 10, 24, 26); 
  ColorRGB(66, 15, 32, 39); 
  ColorRGB(65, 20, 48, 52); 
  ColorRGB(64, 25, 60, 65); 
  ColorRGB(63, 30, 72, 78); 
  ColorRGB(62, 35, 84, 91); 
  ColorRGB(61, 40, 96, 104); 
  ColorRGB(60, 45, 108, 117); 
  ColorRGB(59, 50, 120, 130); 
  ColorRGB(58, 55, 132, 143); 
  ColorRGB(57, 60, 144, 156); 
  ColorRGB(56, 65, 156, 169); 
  ColorRGB(55, 70, 168, 182); 
  ColorRGB(54, 75, 180, 195); 
  ColorRGB(53, 80, 192, 208); 
  ColorRGB(52, 85, 204, 221); 
  ColorRGB(51, 90, 216, 234); 
  ColorRGB(50, 95, 228, 247); 

  ColorRGB(70, 0, 0, 0);
  ColorRGB(71, 13, 7, 7); 
  ColorRGB(72, 26, 14, 14); 
  ColorRGB(73, 39, 21, 21); 
  ColorRGB(74, 52, 28, 28); 
  ColorRGB(75, 65, 35, 35); 
  ColorRGB(76, 78, 42, 42); 
  ColorRGB(77, 91, 49, 49); 
  ColorRGB(78, 104, 56, 56); 
  ColorRGB(79, 117, 63, 63); 
  ColorRGB(80, 130, 70, 70); 
  ColorRGB(81, 143, 77, 77); 
  ColorRGB(82, 156, 84, 84); 
  ColorRGB(83, 169, 91, 91); 
  ColorRGB(84, 182, 98, 98); 
  ColorRGB(85, 195, 105, 105); 
  ColorRGB(86, 208, 112, 112); 
  ColorRGB(87, 221, 119, 119); 
  ColorRGB(88, 234, 126, 126); 
  ColorRGB(89, 247, 133, 133); 


  // Definition of colors: Dark Blue -> White gradient for concentration values (QSbact).
  ColorRGB(100,00,00,85);
  ColorRGB(101,6,6,89);
  ColorRGB(102,12,12,93);
  ColorRGB(103,19,19,97);
  ColorRGB(104,25,25,102);
  ColorRGB(105,31,31,106);
  ColorRGB(106,38,38,110);
  ColorRGB(107,44,44,114);
  ColorRGB(108,51,51,119);
  ColorRGB(109,57,57,123);
  ColorRGB(110,63,63,127);
  ColorRGB(111,70,70,131);
  ColorRGB(112,76,76,136);
  ColorRGB(113,82,82,140);
  ColorRGB(114,89,89,144);
  ColorRGB(115,95,95,148);
  ColorRGB(116,102,102,153);
  ColorRGB(117,108,108,157);
  ColorRGB(118,114,114,161);
  ColorRGB(119,121,121,165);
  ColorRGB(120,127,127,170);
  ColorRGB(121,133,133,174);
  ColorRGB(122,140,140,178);
  ColorRGB(123,146,146,182);
  ColorRGB(124,153,153,187);
  ColorRGB(125,159,159,191);
  ColorRGB(126,165,165,195);
  ColorRGB(127,172,172,199);
  ColorRGB(128,178,178,204);
  ColorRGB(129,184,184,208);
  ColorRGB(130,191,191,212);
  ColorRGB(131,197,197,216);
  ColorRGB(132,204,204,221);
  ColorRGB(133,210,210,225);
  ColorRGB(134,216,216,229);
  ColorRGB(135,223,223,233);
  ColorRGB(136,229,229,238);
  ColorRGB(137,235,235,242);
  ColorRGB(138,242,242,246);
  ColorRGB(139,248,248,250);
  ColorRGB(140,255,255,255);


  /* Create planes by New2() and create neighborhood planes by
     NeighSet(), NeighSetP() and MargolusNeigh(). Also initialize
     them */
  InitialPlane();   /* this must be defined by student */

  /* Create a plane for display */
  viewG = New();

  /*Create a plane for saving pictures*/
  viewG2 = New3();

  /* OpenDisplay() must be after the definition of the color
     table, thus after InitialPlane() */
  if(display==1)
    OpenDisplay("CASH2",nrow+2*margin,ndispplane*ncol+(1+ndispplane)*margin);
}
/*
void Synchronous(int npl,...)
{
  int i,j,k,nr,nc;
  va_list ap;
  static TYPE2 ***prev=NULL; /* they will get planes of user */
/*  static TYPE2 ***next=NULL; /* they will be buffer to make next
				state, which will be copied into
				the user's planes in the end */
/*  static TYPE2 *copy=NULL;


  if(npl == 0){
    fprintf(stderr,"Synchronous() npl==0, that's wrong. See manual.\n");
    exit(-1);
  }

  /* If the first call */
/*  if(next==NULL){
    prev = (TYPE2***) malloc(sizeof(TYPE2**)*npl);

    next = (TYPE2***) malloc(sizeof(TYPE2**)*npl);
    for(i=0;i<npl;i++)
      *(next+i) = New2();

    copy = (TYPE2*) malloc(sizeof(TYPE2)*npl);
  }

  /* get arguments */
/*  va_start(ap,npl);
  for(i=0;i<npl;i++){
    *(prev+i) = va_arg(ap,TYPE2**);
  }

  nr = nrow;
  nc = ncol;

  for(i=1;i<=nr;i++)
    for(j=1;j<=nc;j++){
      /* take a copy of previous states at position [i,j]*/
/*      for(k=0;k<npl;k++)
	*(copy+k) = (*(prev+k))[i][j];

      /* This will modify the state of the [i][j] cell in
	 "*((prev+k))[i][j]". */
/*      NextState(i,j);

      for(k=0;k<npl;k++){
	/* Copy the new state of the [i][j] cells into the "next" */
/*	(*(next+k))[i][j] = (*(prev+k))[i][j];

	/* Revert the state of "prev" because this is synchronous
	   update. */
/*	(*(prev+k))[i][j] = *(copy+k);
      }
    }
  
  /* Then, update all the cells in user's plane, this makes
     CASH-student slower than CASH(or 2). */
/*  for(i=0;i<npl;i++)    
    Copy2(*(prev+i),*(next+i));
}
/*
void Asynchronous(void)
{
  int i;
  int length;
  int irow,icol;
  int testVal=0;

  UpdOrdReset(&updorderG);
  
  length = nrow*ncol;

  i=(length+7)/8;
  switch(length%8)
    {
    case 0: do { irow=(updorderG+testVal)->row;
                 icol=(updorderG+testVal)->col;
                 NextState(irow,icol);
                 testVal++;
    case 7:      irow=(updorderG+testVal)->row;
                 icol=(updorderG+testVal)->col;
                 NextState(irow,icol);
                 testVal++;
    case 6:      irow=(updorderG+testVal)->row;
                 icol=(updorderG+testVal)->col;
                 NextState(irow,icol);
                 testVal++;
    case 5:      irow=(updorderG+testVal)->row;
                 icol=(updorderG+testVal)->col;
                 NextState(irow,icol);
                 testVal++;
    case 4:      irow=(updorderG+testVal)->row;
                 icol=(updorderG+testVal)->col;
                 NextState(irow,icol);
                 testVal++;
    case 3:      irow=(updorderG+testVal)->row;
                 icol=(updorderG+testVal)->col;
                 NextState(irow,icol);
                 testVal++;
    case 2:      irow=(updorderG+testVal)->row;
                 icol=(updorderG+testVal)->col;
                 NextState(irow,icol);
                 testVal++;
    case 1:      irow=(updorderG+testVal)->row;
                 icol=(updorderG+testVal)->col;
                 NextState(irow,icol);
                 testVal++;
               }while(--i>0);
    }
}
*/
/* Usage: Display(Prey,Redator) */
int Display(TYPE2** Arg1,...)
{
  int i;
  int j,k,nr,nc;
  va_list ap;
  TYPE2** a;

  if(display==0)
    return 1;

  nr = nrow;
  nc = ncol;

  if(nplane > 1){
    va_start(ap,Arg1);
  }
  for(i=0;i<ndispplane;i++){
    if(i==0){
      a = Arg1;
    }else{
      a = va_arg(ap,TYPE2**);
    }
    for(j=1;j<=nr;j++)
      for(k=1;k<=nc;k++){
	viewG[j][k] = a[j][k].val;
      }

    PlaneDisplay(viewG,margin,margin+i*(ncol+margin),0);
  }

  return 0;
}

int DrawSlide(TYPE2 **Arg1,char *dirname)
{
  int j,k,nr,nc;
  static int init_flag=0;

  if(init_flag==0){
    OpenPNG(dirname,nrow,ncol);
    init_flag = 1;
    movie_g = 1;
  }
  
  nr = nrow;
  nc = ncol;
  
  for(j=1;j<=nr;j++)
    for(k=1;k<=nc;k++){
      viewG[j][k] = Arg1[j][k].val;
    }
  PlanePNG(viewG,0);
  
  return 0;

}

/*Newly added by Hilje: DrawSlide for more planes at the same time (DrawSlide2) + functions needed to get this (NewP3() and New3())*/
TYPE **NewP3()
{
  TYPE **a;
  a = (TYPE **)calloc(nrow+2+(2*margin), sizeof(TYPE *));
  if (a == NULL) {
    fprintf(stderr,"NewP3: error in memory allocation\n");
    exit(1);
  }
  return a;
}

TYPE **New3()
{
  TYPE **a;
  int i,j;
  a = NewP3(); 
  a[0] = (TYPE *)calloc((nrow+2+(2*margin))*(nplane*(ncol+2)+margin*(1+nplane)),sizeof(TYPE));
  if (a[0] == NULL) {
    fprintf(stderr,"New3: Error in memory allocation\n");
    exit(1);
  }
  for (i=1,j=nrow+2+(2*margin); i < j; i++)
    a[i] = a[i-1] + (nplane*(ncol + 2)+margin*(1+nplane));

  return a;
}

/* PLAN
Make a new viewG2 which has all planes in it, then put that viewG2 in PlanePNG.
viewG2 should have nc*np+nog een beetje columns.
Needed: a new New(), called New3(), which makes a viewG2 which has everything in it.
*/
//Also added by Hilje: extra parameter to indicate how many planes you want to save.
int DrawSlide2(char *dirname, int nplanes, TYPE2 **Arg1, ...)
{
  int i,j,k,j2,k2,nr,nc,nrtot,nctot,np;
  static int init_flag=0;

  va_list ap;
  TYPE2** a;

  nr = nrow;
  nc = ncol;
  np = nplanes;
  nrtot = nr + (2*margin);
  nctot = nplanes*ncol + margin*(1+nplanes);

  if(init_flag==0){
    OpenPNG(dirname,nrtot,nctot);
    init_flag = 1;
    movie_g = 1;
  }

  if(np > 1){
    va_start(ap,Arg1);
  }
  for(i=0;i<nplanes;i++){
    if(i==0){
      a = Arg1;
    }else{
      a = va_arg(ap,TYPE2**);
    }
    for(j=1;j<=nr;j++){
      for(k=1;k<=nc;k++){
        j2 = j + margin;
        k2 = k + i*(ncol) + margin*(i+1);
	viewG2[j2][k2] = a[j][k].val;
      }
    } 	
  } 
  va_end(ap);

 
  PlanePNG(viewG2,0);
  
  return 0;

}


/* InitGradationColor() & UpdateGradationColor() are for coloring
   of floating number states of the cell */
/*
static int minGradColor;
static int maxGradColor;
int InitGradationColor(int ncolor,int min,int max)
{
  int i;

  if(min<0 || max > 255 || min >= max){
    fprintf(stderr,"SetGradationColor: must be 0<min<max<255\n");
    exit(1);
  }

  /* # of color is 1 */
/*  if(ncolor<2){
    for(i=min;i<=max;i++){
      ColorRGB(i,(int)(255.*(double)(i-min)/(double)(max-min)),0,0);
    }
    minGradColor = min;
    maxGradColor = max;
  }
  /* # of color is 2 */
/*  else{
    if((max-min)%2!=0){
      fprintf(stderr,"SetGradationColor: if ncolor=2, max-min must be even\n");
      exit(1);
    }

    for(i=min;i<min+(max-min)/2;i++){
      ColorRGB(i,0,(int)(255.*(double)(min+(max-min)/2-i)/(double)((max-min)/2)),0);
    }
    ColorRGB(i,0,0,0);
    for(i=min+(max-min)/2+1;i<=max;i++){
      ColorRGB(i,(int)(255.*(double)(i-min-(max-min)/2)/(double)((max-min)/2)),0,0);
    }
    minGradColor = min;
    maxGradColor = max;
  }
  return 0;
}

int UpdateGradationColor(TYPE2 **plane,double min,double max,int fnum)
{
  int i,j;
  int nr=nrow,nc=ncol;
  double clen = (double)maxGradColor-minGradColor; /* color length */
/*  double cmin = (double)minGradColor;
  int cmini = minGradColor;
  int cmaxi = maxGradColor;
  double vlen = max - min; /* value length */
/*  double absMax,absMin,abs;

  if(max < min){
    fprintf(stderr,"UpdateGradationColor: must be max > min\n");
    exit(1);
  }

  absMax = (max > -1.*max) ? max : -1.*max;
  absMin = (min > -1.*min) ? min : -1.*min;
  abs = (absMax>absMin) ? absMax : absMin;
  abs = (abs>DBL_MIN) ? abs : DBL_MIN;

  if(vlen/abs < DBL_EPSILON){
    fprintf(stderr,"UpdateGradationColor: must be max != min\n");
    exit(1);
  }


  switch(fnum){
  case 1:
    for(i=1;i<=nr;i++)
      for(j=1;j<=nc;j++){
	/* +0.5 makes x.4 to x, x.6 to x+1*/
/*	plane[i][j].val=clen*(plane[i][j].fval-min)/vlen + cmin +0.5;
	if(plane[i][j].val<cmini)
	  plane[i][j].val=cmini;
	else if(plane[i][j].val>cmaxi)
	  plane[i][j].val=cmaxi;
      }
    break;
  case 2:
    for(i=1;i<=nr;i++)
      for(j=1;j<=nc;j++){
	plane[i][j].val=clen*(plane[i][j].fval2-min)/vlen + cmin +0.5;
	if(plane[i][j].val<cmini)
	  plane[i][j].val=cmini;
	else if(plane[i][j].val>cmaxi)
	  plane[i][j].val=cmaxi;
      }
    break;
  case 3:
    for(i=1;i<=nr;i++)
      for(j=1;j<=nc;j++){
	plane[i][j].val=clen*(plane[i][j].fval3-min)/vlen + cmin +0.5;
	if(plane[i][j].val<cmini)
	  plane[i][j].val=cmini;
	else if(plane[i][j].val>cmaxi)
	  plane[i][j].val=cmaxi;
      }
    break;
  case 4:
    for(i=1;i<=nr;i++)
      for(j=1;j<=nc;j++){
	plane[i][j].val=clen*(plane[i][j].fval4-min)/vlen + cmin +0.5;
	if(plane[i][j].val<cmini)
	  plane[i][j].val=cmini;
	else if(plane[i][j].val>cmaxi)
	  plane[i][j].val=cmaxi;
      }
    break;
  case 5:
    for(i=1;i<=nr;i++)
      for(j=1;j<=nc;j++){
	plane[i][j].val=clen*(plane[i][j].fval5-min)/vlen + cmin +0.5;
	if(plane[i][j].val<cmini)
	  plane[i][j].val=cmini;
	else if(plane[i][j].val>cmaxi)
	  plane[i][j].val=cmaxi;
      }
    break;
  default:
    break;
  }

    return 0;
}
*/

void Plot(int npl,...)
{
  static int count;
  int i,j,k;
  int nc,nr;
  int pop[15]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  va_list ap;
  TYPE2 **world;
  static int initial_flag = 0;

  if(initial_flag == 0){
    InitXmgrace();
    initial_flag =1 ;
  }
  

  nr = nrow;
  nc = ncol;
  
  va_start(ap,npl);
  for(k=0;k<npl;k++){
    world = va_arg(ap,TYPE2**);
    
    for(i=1; i<=nr; i++)
      for(j=1; j<=nc; j++){
	if(world[i][j].val>0 && world[i][j].val<=15)
	  pop[world[i][j].val-1]++;
      }
  }
  
  for(i=0;i<15;i++)
    GracePrintf("g0.s%d point %d, %d",i,Time,pop[i]);


  if(count%10==0 && display==1)
    GracePrintf("redraw");
  
  if(count%100==0 && display==1){
    GracePrintf("focus g0");
    GracePrintf("autoscale");
    GracePrintf("redraw");
  }

  count++;
}

void PlotArray(double data[15])
{
  int i;
  static int count;
  static int initial_flag = 0;

  if(initial_flag == 0){
    InitXmgrace();
    initial_flag =1 ;
  }


  for(i=0;i<15;i++)
    GracePrintf("g0.s%d point %d, %f",i,Time,data[i]);

  if(count%10==0 && display==1)
    GracePrintf("redraw");
  
  if(count%100==0 && display==1){
    GracePrintf("focus g0");
    GracePrintf("autoscale");
    GracePrintf("redraw");
  }

  count++;
}

void SavePlot(char *fname)
{
  if(GraceIsOpen()){
    GraceFlush();
    GracePrintf("saveall \"%s\"",fname);
    GraceClose();
  }
}

/* phase_list will keep track of the phase of the planes: we
   don't know how many planes user have made and which plane user
   wants to do diffusion. */
struct phase_list {
  int phase;
  TYPE2 **plane;
};

void  MDiffusion(TYPE2** a)
{
  int i;
  int planeID;
  static int init_flag=0;
  static struct phase_list *list;
 
  /* initilization */
  if(init_flag==0){
    list = (struct phase_list *)malloc(sizeof(struct phase_list)*nplane);
    for(i=0;i<nplane;i++){
      (list+i)->phase = 0;
      (list+i)->plane = NULL;
    }
    init_flag = 1;
  }

  /* search "a_name" in the list */
  planeID = -1;
  i=0;
  while(i<nplane){
    if((list+i)->plane == a){
      planeID = i;
      break;
    }
    else if((list+i)->plane==NULL){
      /* insert this plane */
      (list+i)->plane = a;
      planeID = i;
      break;
    }
    i++;
  }

  /* weak error cheack */
  if(planeID == -1){
    fprintf(stderr,"MDiffusion: a bug. please report this\n");
    exit(1);
  }

  /* do Margolus diffusion */
  MargolusDiffusion(a,margolusG,(list+planeID)->phase);

  /* incliment the phase */
  (list+planeID)->phase = ((list+planeID)->phase+1)%2;
}


/************************************************************
************ NEIGHBORHOOD RETRIEVAL FUNCTIONS ***************
************************************************************/

int GetNeighbor(TYPE2** a,int row,int col,int nei)
{
  int iplane;
  if(nei<0 || nei>8){
    fprintf(stderr,"GetNeighbor(): wrong value in the forth argument\n");
    exit(1);
  }

  /* we search through planeneiG[] to find which plane user is
     refering to */
  for(iplane=0;iplane<nplane;iplane++){
    if(planeneiG[iplane].plane == a)
      break;
  }

  return((planeneiG[iplane].neighbor[nei][row][col])->val);
}

int RandomMoore8(TYPE2** a,int row,int col)
{
  int iplane;
  int nei;

  /* we search through planeneiG[] to find which plane user is
     refering to */
  for(iplane=0;iplane<nplane;iplane++){
    if(planeneiG[iplane].plane == a)
      break;
  }

  nei = genrand_int(1,8);
  return((planeneiG[iplane].neighbor[nei][row][col])->val);
}

int RandomMoore9(TYPE2** a,int row,int col)
{
  int iplane;
  int nei;

  /* we search through planeneiG[] to find which plane user is
     refering to */
  for(iplane=0;iplane<nplane;iplane++){
    if(planeneiG[iplane].plane == a)
      break;
  }

  nei = genrand_int(0,8);
  return((planeneiG[iplane].neighbor[nei][row][col])->val);
}

int RandomNeumann4(TYPE2** a,int row,int col)
{
  int iplane;
  int nei;

  /* we search through planeneiG[] to find which plane user is
     refering to */
  for(iplane=0;iplane<nplane;iplane++){
    if(planeneiG[iplane].plane == a)
      break;
  }

  nei = genrand_int(1,4);
  return((planeneiG[iplane].neighbor[nei][row][col])->val);
}

int RandomNeumann5(TYPE2** a,int row,int col)
{
  int iplane;
  int nei;

  /* we search through planeneiG[] to find which plane user is
     refering to */
  for(iplane=0;iplane<nplane;iplane++){
    if(planeneiG[iplane].plane == a)
      break;
  }

  nei = genrand_int(0,4);
  return((planeneiG[iplane].neighbor[nei][row][col])->val);
}

int CountMoore8(TYPE2 **plane,int aval,int row,int col)
{
  int iplane;
  
  /* we search through planeneiG[] to find which plane user is
     refering to */
  for(iplane=0;iplane<nplane;iplane++){
    if(planeneiG[iplane].plane == plane)
      break;
  }

  return(((planeneiG[iplane].neighbor[1][row][col])->val==aval)+((planeneiG[iplane].neighbor[2][row][col])->val==aval)+((planeneiG[iplane].neighbor[3][row][col])->val==aval)+((planeneiG[iplane].neighbor[4][row][col])->val==aval)+((planeneiG[iplane].neighbor[5][row][col])->val==aval)+((planeneiG[iplane].neighbor[6][row][col])->val==aval)+((planeneiG[iplane].neighbor[7][row][col])->val==aval)+((planeneiG[iplane].neighbor[8][row][col])->val==aval));
}

int CountMoore9(TYPE2 **plane,int aval,int row,int col)
{
  int iplane;
  
  /* we search through planeneiG[] to find which plane user is
     refering to */
  for(iplane=0;iplane<nplane;iplane++){
    if(planeneiG[iplane].plane == plane)
      break;
  }

  return(((planeneiG[iplane].neighbor[0][row][col])->val==aval)+((planeneiG[iplane].neighbor[1][row][col])->val==aval)+((planeneiG[iplane].neighbor[2][row][col])->val==aval)+((planeneiG[iplane].neighbor[3][row][col])->val==aval)+((planeneiG[iplane].neighbor[4][row][col])->val==aval)+((planeneiG[iplane].neighbor[5][row][col])->val==aval)+((planeneiG[iplane].neighbor[6][row][col])->val==aval)+((planeneiG[iplane].neighbor[7][row][col])->val==aval)+((planeneiG[iplane].neighbor[8][row][col])->val==aval));

}

int CountNeumann4(TYPE2 **plane,int aval,int row,int col)
{
  int iplane;
  
  /* we search through planeneiG[] to find which plane user is
     refering to */
  for(iplane=0;iplane<nplane;iplane++){
    if(planeneiG[iplane].plane == plane)
      break;
  }

  return(((planeneiG[iplane].neighbor[1][row][col])->val==aval)+((planeneiG[iplane].neighbor[2][row][col])->val==aval)+((planeneiG[iplane].neighbor[3][row][col])->val==aval)+((planeneiG[iplane].neighbor[4][row][col])->val==aval));
}

int CountNeumann5(TYPE2 **plane,int aval,int row,int col)
{
  int iplane;
  
  /* we search through planeneiG[] to find which plane user is
     refering to */
  for(iplane=0;iplane<nplane;iplane++){
    if(planeneiG[iplane].plane == plane)
      break;
  }

  return(((planeneiG[iplane].neighbor[0][row][col])->val==aval)+((planeneiG[iplane].neighbor[1][row][col])->val==aval)+((planeneiG[iplane].neighbor[2][row][col])->val==aval)+((planeneiG[iplane].neighbor[3][row][col])->val==aval)+((planeneiG[iplane].neighbor[4][row][col])->val==aval));
}

int SumMoore8(TYPE2**plane,int row,int col)
{
  int iplane;
  
  /* we search through planeneiG[] to find which plane user is
     refering to */
  for(iplane=0;iplane<nplane;iplane++){
    if(planeneiG[iplane].plane == plane)
      break;
  }

  return(((planeneiG[iplane].neighbor[1][row][col])->val)+((planeneiG[iplane].neighbor[2][row][col])->val)+((planeneiG[iplane].neighbor[3][row][col])->val)+((planeneiG[iplane].neighbor[4][row][col])->val)+((planeneiG[iplane].neighbor[5][row][col])->val)+((planeneiG[iplane].neighbor[6][row][col])->val)+((planeneiG[iplane].neighbor[7][row][col])->val)+((planeneiG[iplane].neighbor[8][row][col])->val));
}

int SumMoore9(TYPE2**plane,int row,int col){
  int iplane;
  
  /* we search through planeneiG[] to find which plane user is
     refering to */
  for(iplane=0;iplane<nplane;iplane++){
    if(planeneiG[iplane].plane == plane)
      break;
  }

  return(((planeneiG[iplane].neighbor[0][row][col])->val)+((planeneiG[iplane].neighbor[1][row][col])->val)+((planeneiG[iplane].neighbor[2][row][col])->val)+((planeneiG[iplane].neighbor[3][row][col])->val)+((planeneiG[iplane].neighbor[4][row][col])->val)+((planeneiG[iplane].neighbor[5][row][col])->val)+((planeneiG[iplane].neighbor[6][row][col])->val)+((planeneiG[iplane].neighbor[7][row][col])->val)+((planeneiG[iplane].neighbor[8][row][col])->val));
}

int SumNeumann4(TYPE2**plane,int row,int col)
{
  int iplane;
  
  /* we search through planeneiG[] to find which plane user is
     refering to */
  for(iplane=0;iplane<nplane;iplane++){
    if(planeneiG[iplane].plane == plane)
      break;
  }

  return(((planeneiG[iplane].neighbor[1][row][col])->val)+((planeneiG[iplane].neighbor[2][row][col])->val)+((planeneiG[iplane].neighbor[3][row][col])->val)+((planeneiG[iplane].neighbor[4][row][col])->val));
}
int SumNeumann5(TYPE2**plane,int row,int col)
{
  int iplane;
  
  /* we search through planeneiG[] to find which plane user is
     refering to */
  for(iplane=0;iplane<nplane;iplane++){
    if(planeneiG[iplane].plane == plane)
      break;
  }

  return(((planeneiG[iplane].neighbor[0][row][col])->val)+((planeneiG[iplane].neighbor[1][row][col])->val)+((planeneiG[iplane].neighbor[2][row][col])->val)+((planeneiG[iplane].neighbor[3][row][col])->val)+((planeneiG[iplane].neighbor[4][row][col])->val));
}

TYPE2 GetNeighborS(TYPE2** a,int row,int col,int nei)
{
  int iplane;
  if(nei<0 || nei>8){
    fprintf(stderr,"GetNeighbor(): wrong value in the forth argument\n");
    exit(1);
  }

  /* we search through planeneiG[] to find which plane user is
     refering to */
  for(iplane=0;iplane<nplane;iplane++){
    if(planeneiG[iplane].plane == a)
      break;
  }

  return(*(planeneiG[iplane].neighbor[nei][row][col]));
}

TYPE2 RandomMooreS8(TYPE2** a,int row,int col)
{
  int iplane;
  int nei;

  /* we search through planeneiG[] to find which plane user is
     refering to */
  for(iplane=0;iplane<nplane;iplane++){
    if(planeneiG[iplane].plane == a)
      break;
  }

  nei = genrand_int(1,8);
  return(*(planeneiG[iplane].neighbor[nei][row][col]));
}

TYPE2 RandomMooreS9(TYPE2** a,int row,int col)
{
  int iplane;
  int nei;

  /* we search through planeneiG[] to find which plane user is
     refering to */
  for(iplane=0;iplane<nplane;iplane++){
    if(planeneiG[iplane].plane == a)
      break;
  }

  nei = genrand_int(0,8);
  return(*(planeneiG[iplane].neighbor[nei][row][col]));
}

TYPE2 RandomNeumannS4(TYPE2** a,int row,int col)
{
  int iplane;
  int nei;

  /* we search through planeneiG[] to find which plane user is
     refering to */
  for(iplane=0;iplane<nplane;iplane++){
    if(planeneiG[iplane].plane == a)
      break;
  }

  nei = genrand_int(1,4);
  return(*(planeneiG[iplane].neighbor[nei][row][col]));
}

TYPE2 RandomNeumannS5(TYPE2** a,int row,int col)
{
  int iplane;
  int nei;

  /* we search through planeneiG[] to find which plane user is
     refering to */
  for(iplane=0;iplane<nplane;iplane++){
    if(planeneiG[iplane].plane == a)
      break;
  }

  nei = genrand_int(0,4);
  return(*(planeneiG[iplane].neighbor[nei][row][col]));
}

void GetNeighborC(TYPE2** a,int row,int col,int nei,int*neirow,int*neicol)
{
  if(nei<0 || nei>8){
    fprintf(stderr,"GetNeighbor(): wrong value in the forth argument\n");
    exit(1);
  }

  if(neirow!=NULL)
    *(neirow)=neighborG[nei][row][col].row;
  if(neicol!=NULL)
    *(neicol)=neighborG[nei][row][col].col;
}

void RandomMooreC8(TYPE2** a,int row,int col,int*neirow,int*neicol)
{
  int nei;
  nei = genrand_int(1,8);

  if(neirow!=NULL)
    *(neirow)=neighborG[nei][row][col].row;
  if(neicol!=NULL)
    *(neicol)=neighborG[nei][row][col].col;
}

void RandomMooreC9(TYPE2** a,int row,int col,int*neirow,int*neicol)
{
  int nei;
  nei = genrand_int(0,8);

  if(neirow!=NULL)
    *(neirow)=neighborG[nei][row][col].row;
  if(neicol!=NULL)
    *(neicol)=neighborG[nei][row][col].col;
}

void RandomNeumannC4(TYPE2** a,int row,int col,int*neirow,int*neicol)
{
  int nei;
  nei = genrand_int(1,4);

  if(neirow!=NULL)
    *(neirow)=neighborG[nei][row][col].row;
  if(neicol!=NULL)
    *(neicol)=neighborG[nei][row][col].col;
}

void RandomNeumannC5(TYPE2** a,int row,int col,int*neirow,int*neicol)
{
  int nei;
  nei = genrand_int(0,4);

  if(neirow!=NULL)
    *(neirow)=neighborG[nei][row][col].row;
  if(neicol!=NULL)
    *(neicol)=neighborG[nei][row][col].col;
}

TYPE2* GetNeighborP(TYPE2** a,int row,int col,int nei)
{
  int iplane;
  if(nei<0 || nei>8){
    fprintf(stderr,"GetNeighbor(): wrong value in the forth argument\n");
    exit(1);
  }

  /* we search through planeneiG[] to find which plane user is
     refering to */
  for(iplane=0;iplane<nplane;iplane++){
    if(planeneiG[iplane].plane == a)
      break;
  }

  return(planeneiG[iplane].neighbor[nei][row][col]);
}

TYPE2* RandomMooreP8(TYPE2** a,int row,int col)
{
  int iplane;
  int nei;

  /* we search through planeneiG[] to find which plane user is
     refering to */
  for(iplane=0;iplane<nplane;iplane++){
    if(planeneiG[iplane].plane == a)
      break;
  }

  nei = genrand_int(1,8);
  return(planeneiG[iplane].neighbor[nei][row][col]);
}

TYPE2* RandomMooreP9(TYPE2** a,int row,int col)
{
  int iplane;
  int nei;

  /* we search through planeneiG[] to find which plane user is
     refering to */
  for(iplane=0;iplane<nplane;iplane++){
    if(planeneiG[iplane].plane == a)
      break;
  }

  nei = genrand_int(0,8);
  return(planeneiG[iplane].neighbor[nei][row][col]);
}

TYPE2* RandomNeumannP4(TYPE2** a,int row,int col)
{
  int iplane;
  int nei;

  /* we search through planeneiG[] to find which plane user is
     refering to */
  for(iplane=0;iplane<nplane;iplane++){
    if(planeneiG[iplane].plane == a)
      break;
  }

  nei = genrand_int(1,4);
  return(planeneiG[iplane].neighbor[nei][row][col]);
}

TYPE2* RandomNeumannP5(TYPE2** a,int row,int col)
{
  int iplane;
  int nei;

  /* we search through planeneiG[] to find which plane user is
     refering to */
  for(iplane=0;iplane<nplane;iplane++){
    if(planeneiG[iplane].plane == a)
      break;
  }

  nei = genrand_int(0,4);
  return(planeneiG[iplane].neighbor[nei][row][col]);
}

int gotMouse(void)
{
  char answer[BUFSIZE];

  poseG = 0;
  while(1){
    printf("Enter: Contine (c), OneStep (o), Quit (q): ?");
    fgets(answer,BUFSIZE,stdin);
    switch(answer[0]){
    case 'c':
      return 0;
    case 'o':
      poseG = 1;
      return 0;
    case 'q':
      stopG = 1;
      return 0;
    default:
      printf("Unknown command: %s\n",answer);
      break;
    }
  }


}


/* Make planes by calling New2() and make also neighborhood
   planes by calling NeighSet(), NeighSetP() and
   MargolusNeigh(). */
void MakePlane(TYPE2*** Arg1,...)
{
  int i,j;
  va_list ap;
  TYPE2*** planeP;

  planeneiG = (struct PlaneNei*) malloc(sizeof(struct PlaneNei)*nplane);

  /* if there is more plane, then there should be more arguments */
  if(nplane>1)
    va_start(ap,Arg1);

  for(i=0;i<nplane;i++){
    if(i==0){
      planeP = Arg1;
    }
    else{
      planeP = va_arg(ap,TYPE2***);
    }

    *planeP = New2();

    /* here, we keep track of planes created. */
    planeneiG[i].plane = *planeP;

    /* we also make address based neighborhood planes */
    for(j=0;j<9;j++)
      planeneiG[i].neighbor[j] = NeighSetP(*planeP,j);
  }

  /* we also make coordinate based neighborhood planes */
  for(i=0;i<9;i++)
    neighborG[i]=NeighSet(i);
  margolusG = MargolusNeigh();

}

/*Changed by Hilje.
Original: could only save all planes, else segfault. Now added extra parameter, indicating the number of planes that will be saved (see eg student cash for original). */
//StructChange
void SavePlane(char *fname, TYPE2 **Plane)
{
  int i,j,k,Ntot;
  FILE *fout;
  Ntot = nrow*ncol;

  fout = fopen(fname,"w");
  if(fout == NULL){
    fprintf(stderr,"SavePlane: can't open %s\n",fname);
    exit(1);
  }

  //fprintf(fout,"nrow=%d\nncol=%d\nplane=%d\n",nrow,ncol,nplanes);
  fprintf(fout,"Time=%d\n",Time);
  fprintf(fout,"genrandreal1_offset=%lld\n",genrandreal1_offset);
  fprintf(fout,"genrandreal2_offset=%lld\n",genrandreal2_offset);
  fprintf(fout,"margolus_phase=%d\n",margolus_phase);

  for(i=1;i<=nrow;i++){
    for(j=1;j<=ncol;j++){
      fprintf(fout,"%d %.17g %.17g \n",Plane[i][j].val, Plane[i][j].inf, Plane[i][j].repr);
    }
  }

  fclose(fout);

  initial_flag_UpdOrdReset = 0;		// Reset UpdOrdReset, to make sure that a new start from this plane is the same as continuing the simulation.
}

//StructChange
void ReadSavedData(char *fname, TYPE2** Plane)
{
  int i,j,k;
  FILE *fin;
  char line[1000];
  char *val;

  long long ii, genrand_offset_saved;
  double randdummy;

  fin = fopen(fname,"r");
  if(fin == NULL){
    fprintf(stderr,"ReadSavedData: can't open %s\n",fname);
    exit(1);
  }

/*  // nrow
  fgets(line,2000,fin);
  nr=atoi(strpbrk(line,"=")+1);
  // ncol
  fgets(line,2000,fin);
  nc=atoi(strpbrk(line,"=")+1);
  // nplane
  fgets(line,2000,fin);
  if(npl!=atoi(strpbrk(line,"=")+1)){
    fprintf(stderr,"ReadSavedData: the number of planes is different?\n");
    exit(1);
  }

  if(nrow<nr||ncol<nc){
    fprintf(stderr,"ReadSavedData: nrow or ncol are too small to read the file\n");
    exit(1);
  }
*/
  //Time
  fgets(line,1000,fin);
  Time = atoi(strpbrk(line,"=")+1);

  //genrand_offset
  fgets(line,1000,fin);
  genrand_offset_saved = atoll(strpbrk(line,"=")+1);
  for(ii=1LL;ii<=genrand_offset_saved;ii++){
    randdummy = genrand_real1();
  }
  fgets(line,1000,fin);
  genrand_offset_saved = atoll(strpbrk(line,"=")+1);
  for(ii=1LL;ii<=genrand_offset_saved;ii++){		
    randdummy = genrand_real2();
  }

  //Margolus phase
  fgets(line,1000,fin);
  margolus_phase = atoi(strpbrk(line,"=")+1);

  for(i=1;i<=nrow;i++){
    for(j=1;j<=ncol;j++){
	fgets(line,1000,fin);

     	val = strtok(line," ");
	Plane[i][j].val = atoi(val);

     	val = strtok(NULL," ");
	Plane[i][j].inf = atof(val);

     	val = strtok(NULL," ");
	Plane[i][j].repr = atof(val);
    }
  }

  fclose(fin);
}

struct space_time_rec {
  TYPE2 **plane;  
  int current; /* the row of the destination, where we should draw somehting */
};
void SpaceTimePlot(TYPE2 **dest, TYPE2 **source)
{
  int nr=nrow,nc=ncol;
  int i;
  static struct space_time_rec *list;
  static int init_flag = 0;
  int planeID;

  if(init_flag == 0){
    list = (struct space_time_rec*) malloc(sizeof(struct space_time_rec)*nplane);
    for(i=0;i<nplane;i++){
      (list+i)->plane=NULL;
      (list+i)->current = 1;
    }
    init_flag = 1;
  }

  /* search through the list */
  planeID = -1;
  for(i=0;i<nplane;i++){
    if((list+i)->plane == dest){
      planeID = i;
      break;
    }
    else if((list+i)->plane==NULL){
      /* insert this plane */
      (list+i)->plane = dest;
      planeID = i;
      break;
    }
  }

  /* weak error cheack */
  if(planeID == -1){
    fprintf(stderr,"SpaceTimePlot: a bug. please report this\n");
    exit(1);
  }

  memcpy(dest[(list+planeID)->current],source[nr/2],sizeof(TYPE2)*(nc+2));
  (list+planeID)->current = ((list+planeID)->current % nr)+1;
  
}

//StructChange
//NOTE: Might be much better to use InitialSetS(), see below!
/* A easier variant of InitailSetS() */
int InitialSet(TYPE2** plane,int ntype,int empty_int,...)
{
  va_list ap;

  TYPE2 *type_p; /* the types of cell */
  double *port_p; /* portion */
  int *numb_p; /* acutual number calculated from the portion */

  double port_total=0.;
  int nr,nc,pos,maxpos;
  int i,ii,jj;
  TYPE2 empty; 			// empty={0,0,0,0,0,0,0,0,{0.,0.},{0.,0.},0.,0}
  empty.val = 0;
  empty.inf = 0.;
  empty.repr = 0.;
  empty.exp_offspring = 0;
  empty.real_offspring = 0;
  empty.real_trans = 0.;

  /* If there is no other type than "empty" */
  if(ntype<=0){
    empty.val = empty_int;
    Fill2(plane,empty);
    return 0;
  }

  nr = nrow;
  nc = ncol;

  /* get argument */
  type_p = (TYPE2*)malloc(sizeof(TYPE2)*ntype);
  port_p = (double*)malloc(sizeof(double)*ntype);
  numb_p = (int*)malloc(sizeof(int)*ntype);

  va_start(ap,empty_int);
  empty.val = empty_int;
  for(i=0;i<ntype;i++){
    *(type_p+i) = empty;
    (type_p+i)->val = va_arg(ap,int);
    *(port_p+i) =  va_arg(ap,double);

    /* if portion<0.0, make it 0.0 */
    if(*(port_p+i)<0.0)
      *(port_p+i) = 0.0;
    port_total += *(port_p+i);
  }
  
  /* if P1+P2 > 1, then P1=P1/(P1+P2) */
  if(port_total>1.0){
    for(i=0;i<ntype;i++){
      *(port_p+i) = (*(port_p+i))/port_total;
    }
  }

  /* 
     Set the number of cells which will be S1, S2, ... 
     The fractional part will be truncated so that the sum of the
     number will not be nrow*ncol.
  */
  for(i=0;i<ntype;i++){
    *(numb_p+i) = (*(port_p+i))*nr*nc;
  }


  /* From here, we will set the state of the cell to S1, S2 ...*/

  /* "pos" is the number of the current cell [i][j]. */
  pos = nc; /* if pos==ncol, then [i][j]=[1][1].
	       if pos==ncol+1, then [i][j]=[1][2]
	       if pos==ncol+ncol+2, then [i][j]=[2][3]
	       if pos==x, then [i][j]=[x/ncol][x%ncol+1]
	       if pos==ncol+ncol*nrow-1, then [i][j]=[nrow][ncol]
	    */
  maxpos = nc+nc*nr-1;

  /* fill the plane with the given state, S0 S1 S2 ..., we will
     shafful the plane later */
  for(i=0;i<ntype;i++){
    while(*(numb_p+i)>0){
      ii = pos/nc;
      jj = pos%nc + 1;

      plane[ii][jj] = *(type_p+i);

      pos++;
      (*(numb_p+i))--;
    }
  }
  if(pos>maxpos+1){
    fprintf(stderr,"InitialSet(): this is a bug. Report please");
    exit(1);
  }

  while(pos<=maxpos){
    ii = pos/nc;
    jj = pos%nc + 1;
    
    plane[ii][jj] = empty;
    
    pos++;
  }

  PerfectMix(plane);

  free(type_p);
  free(port_p);
  free(numb_p);

  return 0;
}

/* This is exaclty the same as InitialPlaneSet() in cash2.c */
int InitialSetS(TYPE2** a,int ntype,TYPE2 empty,...)
{
  va_list ap;
  TYPE2 *type_p; /* the types of cell */
  double *port_p; /* portion */
  int *numb_p; /* acutual number calculated from the portion */
  double port_total=0.;
  int nr,nc,pos,maxpos;
  int i,ii,jj;

  /* If there is no other type thant "empty" */
  if(ntype<=0){
    Fill2(a,empty);
    return 0;
  }

  nr = nrow;
  nc = ncol;

  /* get argument */
  type_p = (TYPE2*)malloc(sizeof(TYPE2)*ntype);
  port_p = (double*)malloc(sizeof(double)*ntype);
  numb_p = (int*)malloc(sizeof(int)*ntype);

  va_start(ap,empty);

  for(i=0;i<ntype;i++){
    *(type_p+i) = va_arg(ap,TYPE2);
    *(port_p+i) =  va_arg(ap,double);

    /* if portion<0.0, make it 0.0 */
    if(*(port_p+i)<0.0)
      *(port_p+i) = 0.0;
    port_total += *(port_p+i);
  }
  
  /* if P1+P2 > 1, then P1=P1/(P1+P2) */
  if(port_total>1.0){
    for(i=0;i<ntype;i++){
      *(port_p+i) = (*(port_p+i))/port_total;
    }
  }

  /* 
     Set the number of cells which will be S1, S2, ... 
     The fractional part will be truncated so that the sum of the
     number will not be nrow*ncol.
  */
  for(i=0;i<ntype;i++){
    *(numb_p+i) = (*(port_p+i))*nr*nc;
  }


  /* From here, we will set the state of the cell to S1, S2 ...*/

  /* "pos" is the number of the current cell [i][j]. */
  pos = nc; /* if pos==ncol, then [i][j]=[1][1].
	       if pos==ncol+1, then [i][j]=[1][2]
	       if pos==ncol+ncol+2, then [i][j]=[2][3]
	       if pos==x, then [i][j]=[x/ncol][x%ncol+1]
	       if pos==ncol+ncol*nrow-1, then [i][j]=[nrow][ncol]
	    */
  maxpos = nc+nc*nr-1;

  /* fill the plane with the given state, S0 S1 S2 ..., we will
     shafful the plane later */
  for(i=0;i<ntype;i++){
    while(*(numb_p+i)>0){
      ii = pos/nc;
      jj = pos%nc + 1;

      a[ii][jj] = *(type_p+i);

      pos++;
      (*(numb_p+i))--;
    }
  }
  if(pos>maxpos+1){
    fprintf(stderr,"InitialPlaneSet(): this is a bug. Report please");
    exit(1);
  }

  while(pos<=maxpos){
    ii = pos/nc;
    jj = pos%nc + 1;
    
    a[ii][jj] = empty;
    
    pos++;
  }

  PerfectMix(a);

  free(type_p);
  free(port_p);
  free(numb_p);

  return 0;
}

/***********************************
 Margolus diffusion with obstacles by Anton
**********************************/

/* Some auxilary functions for rotating the particles
 * 
 * These functions should be inlined.... as if they were macros.
 */
inline void
rotate_square( TYPE2 *a, TYPE2 *b, TYPE2 *c, TYPE2 *d ) {
    TYPE2 aux;

    if( genrand_real1() < 0.5 ) {
        aux = *a;
        *a = *b;
        *b = *c;
        *c = *d;
        *d = aux;
    } else {
        aux = *d;
        *d = *c;
        *c = *b;
        *b = *a;
        *a = aux;
    }
}

inline void
rotate_triangle( TYPE2 *a, TYPE2 *b, TYPE2 *c ) {
    TYPE2 aux;

    if( genrand_real1() < 0.5 ) {
        aux = *a;
        *a = *b;
        *b = *c;
        *c = aux;
    } else {
        aux = *c;
        *c = *b;
        *b = *a;
        *a = aux;
    }
}

inline void
rotate_pair( TYPE2 *a, TYPE2 *b ) {
    TYPE2 aux;

    aux = *a;
    *a = *b;
    *b = aux;
}


/* Modified Margolis diffusion. Now we can diffuse with obstacles
 * in our grid as well.
 * 
 * The function expects two planes, a diffusion plane and an
 * obstacle plane. For the obstacle plane one needs to supply an
 * example obstacle. (There may be more states in the plane of
 * which only one is an obstacle.)
 *
 * Assuming both planes have equal dimensions.
 */
void
MargolusWithObstacle( TYPE2** d_plane, TYPE2** o_plane, int obstacle, 
                      MARGOLUS*** margolus, int phase ) {
  int my_nrow, my_ncol;
  int i, j, ii;
  unsigned int mask;
  int sq_row[ 4 ], sq_col[ 4 ];



  my_nrow = nrow - 1;
  my_ncol = ncol - 1;

  if( boundary == WRAP ) {
    // use margolis neighbourhood and rotate according to obstacles
    for( i = 1; i <= my_nrow; i += 2 ) {
      for( j = 1; j <= my_ncol; j += 2 ) {

	// 1. determine nbh

	/* (i,j)=HERE */
	sq_row[ 0 ] = i;
	sq_col[ 0 ] = j;
	/* (cw_r,cw_c)=CW_of_(i,j) */
	sq_row[ 1 ] = margolus[phase][i][j].m[CW].row;
	sq_col[ 1 ] = margolus[phase][i][j].m[CW].col;
	/* (ccw_r,ccw_c)=CCW_f_(i,j) */
	sq_row[ 2 ] = margolus[phase][i][j].m[CCW].row;
	sq_col[ 2 ] = margolus[phase][i][j].m[CCW].col;
	/* (opp_r,opp_c)=OPP_of_(i,j) */
	sq_row[ 3 ] = margolus[phase][i][j].m[OPP].row;
	sq_col[ 3 ] = margolus[phase][i][j].m[OPP].col;

	// 2. look for obstacles
	mask = 0;
	for( ii = 0; ii < 4; ++ii ) {
	  // FIX whole structure comparison?
	  if( o_plane[ sq_row[ ii ] ][ sq_col[ ii ] ].val == obstacle ) {
	    mask = mask | ( 1 << ii );
	  }
	}
        
	// 3. rotate depending on nr and position of obstacles
	/* The neighbourhood looks like this:
	 * 
	 * ---------
	 * | 0 | 1 |
	 * ---------
	 * | 2 | 3 |
	 * ---------
	 */
	switch( mask ) {
	case 0:
	  // rotate 0, 1, 2, 3
	  rotate_square( &d_plane[ sq_row[ 0 ] ][ sq_col[ 0 ] ],
			 &d_plane[ sq_row[ 1 ] ][ sq_col[ 1 ] ],
			 &d_plane[ sq_row[ 3 ] ][ sq_col[ 3 ] ],
			 &d_plane[ sq_row[ 2 ] ][ sq_col[ 2 ] ] );
	  break;
	case 1:
	  // rotate 0, 1, 2
	  rotate_triangle( &d_plane[ sq_row[ 0 ] ][ sq_col[ 0 ] ],
			   &d_plane[ sq_row[ 1 ] ][ sq_col[ 1 ] ],
			   &d_plane[ sq_row[ 2 ] ][ sq_col[ 2 ] ] );
	  break;
	case 2:
	  // rotate 0, 1, 3
	  rotate_triangle( &d_plane[ sq_row[ 0 ] ][ sq_col[ 0 ] ],
			   &d_plane[ sq_row[ 1 ] ][ sq_col[ 1 ] ],
			   &d_plane[ sq_row[ 3 ] ][ sq_col[ 3 ] ] );
	  break;
	case 3:
	  // rotate 0, 1
	  rotate_pair( &d_plane[ sq_row[ 0 ] ][ sq_col[ 0 ] ],
		       &d_plane[ sq_row[ 1 ] ][ sq_col[ 1 ] ] );
	  break;
	case 4:
	  // rotate 0, 2, 3
	  rotate_triangle( &d_plane[ sq_row[ 0 ] ][ sq_col[ 0 ] ],
			   &d_plane[ sq_row[ 3 ] ][ sq_col[ 3 ] ],
			   &d_plane[ sq_row[ 2 ] ][ sq_col[ 2 ] ] );
	  break;
	case 5:
	  // rotate 0, 2
	  rotate_pair( &d_plane[ sq_row[ 0 ] ][ sq_col[ 0 ] ],
		       &d_plane[ sq_row[ 2 ] ][ sq_col[ 2 ] ] );
	  break;
	case 8:
	  // rotate 1, 2, 3
	  rotate_triangle( &d_plane[ sq_row[ 1 ] ][ sq_col[ 1 ] ],
			   &d_plane[ sq_row[ 2 ] ][ sq_col[ 2 ] ],
			   &d_plane[ sq_row[ 3 ] ][ sq_col[ 3 ] ] );
	  break;
	case 10:
	  // rotate 1, 3
	  rotate_pair( &d_plane[ sq_row[ 1 ] ][ sq_col[ 1 ] ],
		       &d_plane[ sq_row[ 3 ] ][ sq_col[ 3 ] ] );
	  break;
	case 12:
	  // rotate 2, 3
	  rotate_pair( &d_plane[ sq_row[ 2 ] ][ sq_col[ 2 ] ],
		       &d_plane[ sq_row[ 3 ] ][ sq_col[ 3 ] ] );
	  break;
	}
      }
    }
  } else if( boundary == FIXED ) {
    fprintf(stderr,"MargolusWithObstacle() Use WRAP boundary\n");
    exit(-1);
  }
}

void  ObstacleMargolus( TYPE2 **a, TYPE2 **o_plane, int obstacle )
{
  int i;
  int planeID;
  static int init_flag=0;
  static struct phase_list *list;
 
  /* initilization */
  if(init_flag==0){
    list = (struct phase_list *)malloc(sizeof(struct phase_list)*nplane);
    for(i=0;i<nplane;i++){
      (list+i)->phase = 0;
      (list+i)->plane = NULL;
    }
    init_flag = 1;
  }

  /* search "a_name" in the list */
  planeID = -1;
  i=0;
  while(i<nplane){
    if((list+i)->plane == a){
      planeID = i;
      break;
    }
    else if((list+i)->plane==NULL){
      /* insert this plane */
      (list+i)->plane = a;
      planeID = i;
      break;
    }
    i++;
  }

  /* weak error cheack */
  if(planeID == -1){
    fprintf(stderr,"ObstacleMargolis: a bug. please report this\n");
    exit(1);
  }

  /* do Margolus diffusion */
  MargolusWithObstacle( a, o_plane, obstacle, margolusG, (list+planeID)->phase );

  /* incliment the phase */
  (list+planeID)->phase = ((list+planeID)->phase+1)%2;
}

//VARIANT ON MARGOLUS WITH OBSTACLES (NEXT TWO FUNCTIONS), NEWLY ADDED BY HILJE
//ALSO ADDED: MULTIPLE STATES CAN BE OBSTACLES (PASSED IN ARRAY "OBSTACLE")
/*
void
MargolusWithFixedPoints( TYPE2** d_plane, int* obstacle,	 				
                      MARGOLUS*** margolus, int phase ) {
  int my_nrow, my_ncol;
  int i, j, ii, k;
  int length = (int) (sizeof(obstacle)/sizeof(int));
  unsigned int mask;
  int sq_row[ 4 ], sq_col[ 4 ];

  my_nrow = nrow - 1;
  my_ncol = ncol - 1;

  if( boundary == WRAP ) {
    // use margolus neighbourhood and rotate according to fixed points
    for( i = 1; i <= my_nrow; i += 2 ) {
      for( j = 1; j <= my_ncol; j += 2 ) {

	// 1. determine nbh

	/* (i,j)=HERE */
/*	sq_row[ 0 ] = i;
	sq_col[ 0 ] = j;
	/* (cw_r,cw_c)=CW_of_(i,j) */
/*	sq_row[ 1 ] = margolus[phase][i][j].m[CW].row;
	sq_col[ 1 ] = margolus[phase][i][j].m[CW].col;
	/* (ccw_r,ccw_c)=CCW_f_(i,j) */
/*	sq_row[ 2 ] = margolus[phase][i][j].m[CCW].row;
	sq_col[ 2 ] = margolus[phase][i][j].m[CCW].col;
	/* (opp_r,opp_c)=OPP_of_(i,j) */
/*	sq_row[ 3 ] = margolus[phase][i][j].m[OPP].row;
	sq_col[ 3 ] = margolus[phase][i][j].m[OPP].col;

	// 2. look for obstacles
	mask = 0;
	for( ii = 0; ii < 4; ii++) {
          bool obs = false; 

	  //determine whether d_plane[sq_row[ii]][sq_col[ ii ]].val is an obstacle (obs = true) or not (obs = false)
	  for(k=0; k<length; k++)
	  {
            if(d_plane[ sq_row[ ii ] ][ sq_col[ ii ] ].val == obstacle[k] ){
              obs = true;
	      break;
	    }
	  }

	  // mask gets value depending on which squares contain obstacles.
	  if( obs ) {
	    mask += pow(2,ii);
	  }
	}
        
	// 3. rotate depending on nr and position of fixed points
	/* The neighbourhood looks like this:
	 * 
	 * ---------
	 * | 0 | 1 |
	 * ---------
	 * | 2 | 3 |
	 * ---------
	 */
/*	switch( mask ) {
	case 0:
	  // rotate 0, 1, 2, 3
	  rotate_square( &d_plane[ sq_row[ 0 ] ][ sq_col[ 0 ] ],
			 &d_plane[ sq_row[ 1 ] ][ sq_col[ 1 ] ],
			 &d_plane[ sq_row[ 3 ] ][ sq_col[ 3 ] ],
			 &d_plane[ sq_row[ 2 ] ][ sq_col[ 2 ] ] );
	  break;
	case 1:
	  // rotate 1, 2, 3
	  rotate_triangle( &d_plane[ sq_row[ 1 ] ][ sq_col[ 1 ] ],
			   &d_plane[ sq_row[ 2 ] ][ sq_col[ 2 ] ],
			   &d_plane[ sq_row[ 3 ] ][ sq_col[ 3 ] ] );
	  break;
	case 2:
	  // rotate 0, 2, 3
	  rotate_triangle( &d_plane[ sq_row[ 0 ] ][ sq_col[ 0 ] ],
			   &d_plane[ sq_row[ 3 ] ][ sq_col[ 3 ] ],
			   &d_plane[ sq_row[ 2 ] ][ sq_col[ 2 ] ] );
	  break;
	case 3:
	  // rotate 2, 3
	  rotate_pair( &d_plane[ sq_row[ 2 ] ][ sq_col[ 2 ] ],
		       &d_plane[ sq_row[ 3 ] ][ sq_col[ 3 ] ] );
	  break;
	case 4:
	  // rotate 0, 1, 3
	  rotate_triangle( &d_plane[ sq_row[ 0 ] ][ sq_col[ 0 ] ],
			   &d_plane[ sq_row[ 1 ] ][ sq_col[ 1 ] ],
			   &d_plane[ sq_row[ 3 ] ][ sq_col[ 3 ] ] );
	  break;
	case 5:
	  // rotate 1, 3
	  rotate_pair( &d_plane[ sq_row[ 1 ] ][ sq_col[ 1 ] ],
		       &d_plane[ sq_row[ 3 ] ][ sq_col[ 3 ] ] );
	  break;
	case 8:
	  // rotate 0, 1, 2
	  rotate_triangle( &d_plane[ sq_row[ 0 ] ][ sq_col[ 0 ] ],
			   &d_plane[ sq_row[ 1 ] ][ sq_col[ 1 ] ],
			   &d_plane[ sq_row[ 2 ] ][ sq_col[ 2 ] ] );
	  break;
	case 10:
	  // rotate 0, 2
	  rotate_pair( &d_plane[ sq_row[ 0 ] ][ sq_col[ 0 ] ],
		       &d_plane[ sq_row[ 2 ] ][ sq_col[ 2 ] ] );
	  break;
	case 12:
	  // rotate 0, 1
	  rotate_pair( &d_plane[ sq_row[ 0 ] ][ sq_col[ 0 ] ],
		       &d_plane[ sq_row[ 1 ] ][ sq_col[ 1 ] ] );
	  break;
	}
      }
    }
  } else if( boundary == FIXED ) {
    fprintf(stderr,"MargolusWithFixedPoints() Use WRAP boundary\n");
    exit(-1);
  }
}

void  FixedPointsMargolus( TYPE2 **a, int* obstacle )
{
  int i;
  int planeID;
  static int init_flag=0;
  static struct phase_list *list;
 
  /* initilization */
/*  if(init_flag==0){
    list = (struct phase_list *)malloc(sizeof(struct phase_list)*nplane);
    for(i=0;i<nplane;i++){
      (list+i)->phase = 0;
      (list+i)->plane = NULL;
    }
    init_flag = 1;
  }

  /* search "a_name" in the list */
/*  planeID = -1;
  i=0;
  while(i<nplane){
    if((list+i)->plane == a){
      planeID = i;
      break;
    }
    else if((list+i)->plane==NULL){
      /* insert this plane */
/*      (list+i)->plane = a;
      planeID = i;
      break;
    }
    i++;
  }

  /* weak error check */
/*  if(planeID == -1){
    fprintf(stderr,"FixedPointsMargolus: a bug. please report this\n");
    exit(1);
  }

  /* do Margolus diffusion */
/*  MargolusWithFixedPoints( a, obstacle, margolusG, (list+planeID)->phase );

  /* incliment the phase */
/*  (list+planeID)->phase = ((list+planeID)->phase+1)%2;
}
*/

int countGlobal(TYPE2 **world, int num){

  int pop = 0;
  int i,j;
  for(i=1;i<=nrow;i++){                     
    for(j=1;j<=ncol;j++){                     
     if(world[i][j].val == num)
          pop++;
     }
  }
  return pop;
}

void PlotXY(double x, double y)
{
  static int count;
  static int initial_flag = 0;

  if(initial_flag == 0){
    InitXmgrace();
    initial_flag =1 ;
  }
   
  GracePrintf("g0.s%d point %f, %f",1,x,y);

  if(count%10==0)
    GracePrintf("redraw");
  
  if(count%100==0){
    GracePrintf("focus g0");
    GracePrintf("autoscale");
    GracePrintf("redraw");
  }

  count++;
}

/*************************************************
****** NEW FUNCTIONS ADDED BY HILJE **************
*************************************************/

void DiffusionBySwap(double D, TYPE2** plane)
{
	int  k, kmod, irow, icol, nswaps, ntot;
	TYPE2* nei;
	TYPE2 temp;

	// Total number of grid points
	ntot = nrow*ncol;
	// Total number of swaps to be performed
	nswaps = (int) (D*ntot) + 1;

	// Give a random order of grid points, store (row,col)-values in updorderG
	UpdOrdReset(&updorderG);

	// For the first nswaps grid points (nswap > ntot is allowed!), perform a swap if one of them in non-empty
	for(k=0;k<nswaps;k++)
	{
		kmod = k%ntot;	// Make sure that you start again at the beginning of updorder if nswap > ntot
		irow = (updorderG[kmod]).row;
		icol = (updorderG[kmod]).col;

		// Pick a random neighbour, get pointer to this neighbour
		nei = RandomNeumannP4(plane, irow, icol);

		// Swap the two grid points if either one is non-empty
		if(nei->val != 0 || plane[irow][icol].val != 0)
		{
			temp = *nei;
			*nei = plane[irow][icol];
			plane[irow][icol] = temp;
		}
	}
}


int MaxIntArray(int* intarr, int length)
{
	int max = 0;
	int i;

	for(i=0; i<length; i++)
	{
		if(intarr[i]>max)
			max = intarr[i];
	}

	return max;
}

int Max2DIntArray(int** intarr, int xlength, int ylength)
{
	int max = 0;
	int i,j;

	for(i=0; i<xlength; i++){
		for(j=0; j<ylength; j++){
			if(intarr[i][j]>max)
				max = intarr[i][j];
		}
	}

	return max;
}

double Max2Doubles(double a1, double a2)
{
	if(a1 >= a2)
		return a1;
	else
		return a2;
}

int Min2Ints(int a1, int a2)
{
	if(a1 <= a2)
		return a1;
	else
		return a2;
}

// Return number taken from exponential distribution with mean m.
double ExpDistr(double m)
{
	double ans;
	ans = -log(1-genrand_real2())*m;
	return ans;
}

// Returns exponentially distributed waiting time (as ExpDistr, but then rounded to an integer).
int ExpWaitingTime(int m)
{
	int ans;
	double m2;
	m2 = (double) m;
	ans = (int)(ExpDistr(m2)+0.5);
	return ans;
}

