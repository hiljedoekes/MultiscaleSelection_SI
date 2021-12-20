#include <stdio.h>
#include "cash2003.h"
#include <string.h>

FILE *parfile;

/* basic parameters which are not directly related to user's specific modelling */
int seed=55;
int seedset=0;
int graphics=1;
int nrow=100;
int ncol=100;
int first=0;
int last=0;
int boundary=0; /* this is used by both "Boundaries()" and "Boundaries2()".*/
int boundaryvalue=0;
int scale=1;
int ranmargolus=1;
int psborder=1;
int psreverse=1;

int InDat(char *format,char *text,int *value)
{
  FullInDat(parfile,format,text,value,1);
  return (0);
}

int FullInDat(FILE *fp,char *format,char *text,int *value,int nvalues)
{
  int j;
  char line[80];

  rewind(fp);
  while (fscanf(fp,"%s",line) != EOF) {
    if (!strcmp(line, text)) {
      for (j=0; j<nvalues; j++)
        fscanf(fp,format,&value[j]);
      return (j);
    }
    while (getc(fp) != '\n');
  }
  return (0);
}

int ReadOptions(char filename[])
{
  parfile = fopen(filename,"r");
  if (parfile == NULL) {
    fprintf(stderr,"%s: File does not exist!\n",filename);
    exit(1);
  }
  FullInDat(parfile,"%d","nrow",&nrow,1);
  FullInDat(parfile,"%d","ncol",&ncol,1);
  FullInDat(parfile,"%d","boundary",&boundary,1);
  FullInDat(parfile,"%d","boundaryvalue",&boundaryvalue,1);
  FullInDat(parfile,"%d","graphics",&graphics,1);
  FullInDat(parfile,"%d","scale",&scale,1);
  FullInDat(parfile,"%d","seed",&seed,1);
  FullInDat(parfile,"%d","psborder",&psborder,1);
  FullInDat(parfile,"%d","psreverse",&psreverse,1);
  FullInDat(parfile,"%d","ranmargolus",&ranmargolus,1);
  if (!seedset) {
    SEED(seed);
    seedset = 1;
  }
  return (0);
}

TYPE **NewP()
{
  TYPE **a;
  a = (TYPE **)calloc(nrow+2, sizeof(TYPE *));
  if (a == NULL) {
    fprintf(stderr,"NewP: error in memory allocation\n");
    exit(1);
  }
  return a;
}

TYPE **New()
{
  TYPE **a;
  int i,j;
  a = NewP(); 
  a[0] = (TYPE *)calloc((nrow+2)*(ncol+2),sizeof(TYPE));
  if (a[0] == NULL) {
    fprintf(stderr,"Error in memory allocation\n");
    exit(1);
  }
  for (i=1,j=nrow+2; i < j; i++)
    a[i] = a[i-1] + ncol + 2;

  if (first == 0) {
    first = ncol + 3;
    last = (nrow+2)*(ncol+2) - ncol - 4;
  }
  return a;
}

int PlaneFree(TYPE **a)
{
  int i,j;
  /*$dir force_vector*/
  for (i=0,j=nrow+2; i < j; i++)
    free(a[i]);
  free(a);
  return 0;
}

TYPE **Fill(TYPE **a,TYPE c)
{
  int i,j;
  /*$dir force_vector*/
  for (i=first,j=last; i <= j; i++)
    a[0][i] = c;
  return a;
}

long Total(TYPE **a)
{
  int i,j;
  long k=0;
  int nc, nr;
  nr = nrow;
  nc = ncol;
  for (i=1; i <= nr; i++)
    /*$dir force_vector*/
    for (j=1; j <= nc; j++)
    k += a[i][j];
  return k;
}

int Equal(TYPE **a,TYPE  **b)
{
  int i,j,e=1;
  /*$dir force_vector*/
  for (i=first,j=last; i <= j; i++)
    if (a[0][i] != b[0][i]) e=0;
  return e;
}

int Max(TYPE **a)
{
  int i,j,k;
  int nc, nr;
  nr = nrow;
  nc = ncol;
  k=a[1][1];
  for (i=1; i <= nr; i++)
    /*$dir force_vector*/
    for (j=1; j <= nc; j++)
    if (a[i][j] > k) k = a[i][j];
  return k;
}

int Min(TYPE **a)
{
  int i,j,k;
  int nc, nr;
  nr = nrow;
  nc = ncol;
  k=a[1][1];
  for (i=1; i <= nr; i++)
    /*$dir force_vector*/
    for (j=1; j <= nc; j++)
    if (a[i][j] < k) k = a[i][j];
  return k;
}

TYPE **Explode(TYPE **a,TYPE  **b)
{
  int i,j;
  int nc, nr;
  nr = nrow;
  nc = ncol;
  for (i=1; i <= nr; i++)
    /*$dir force_vector*/
    for (j=1; j <= nc; j++)
      a[i][j] = b[i][j];
  return a;
}

TYPE **Copy(TYPE **a,TYPE  **b)
{
  int i,j;
  /*$dir force_vector*/
  for (i=first,j=last; i <= j; i++)
    a[0][i] = b[0][i];
  return a;
}

TYPE *CopyRow(TYPE *a,TYPE  *b)
{
  int i,j;
  /*$dir force_vector*/
  for (i=1,j=ncol; i <= j; i++) 
    a[i] = b[i];
  return a;
}

static int diffinit = 0;
static TYPE **diff;

#define DIFF(A,B) ((A) > (B) ? ((A)--,(B)++) : (A) < (B) ? ((A)++,(B)--) : 0)
#define DFBD(A,B) ((A) > (B) ? (B)++ : (A) < (B) ? (B)-- : 0)

TYPE **Diffusion(TYPE **a,TYPE **b,float d)
{
  int i, j, nc, nr;
  nr = nrow;
  nc = ncol;
  if (!diffinit) {
    diff = New();
    diffinit = 1;
  }
  VonNeumann4(diff,b);
  for (i=1; i <= nr; i++)
    /*$dir force_vector*/
    for (j=1; j <= nc; j++)
      a[i][j] = b[i][j] + d*diff[i][j] - 4*d*b[i][j] + 0.5;
  return a;
}

static int motioninit = 0;
static TYPE **ran0,**ran1;

static TYPE T;
#define SWAP(A,B) {T = A; A = B; B = T;}
#define SWBD(A,B) B = A;

TYPE **Motion(TYPE **a,TYPE **b,float d,int time)
{
  int i,ii,j,jj,nr,nc,iodd,jodd;
  int mode, bnd, bndval;
  mode = time & 1;
  nr = nrow;
  nc = ncol;
  bnd = boundary;
  bndval = boundaryvalue;
  if (!motioninit) {
    ran0 = New();
    ran1 = New();
    motioninit = 1;
  }
  if (a != b) Copy(a,b);
  Random(ran0,0.5);
  Random(ran1,d);

  for (i=1; i < nr; i++) {
    iodd = i & 1;
    ii = i + 1;
    /*$dir force_vector*/
    for (j=1; j < nc; j++) {
      jodd = j & 1;
      jj = j + 1;
      if (iodd == mode && ran0[i][j] && ran1[ii][j]) SWAP(a[ii][j],a[i][j])
      else if (jodd == mode && ran1[i][jj]) SWAP(a[i][jj],a[i][j]);
    }
  }
  if (mode == (nr & 1) && bnd < 2) {
    /*$dir force_vector*/
    for (j=1; j <= nc; j++)
    if (ran0[1][j] && ran1[1][j]) {
      if (bnd == 1) {SWBD(bndval,a[nr][j]); SWBD(bndval,a[1][j]);}
      else SWAP(a[1][j],a[nr][j]);
    }
  }
  if (mode == (nc & 1) && bnd < 2) {
    /*$dir force_vector*/
    for (i=1; i <= nr; i++)
    if (ran0[i][1] && ran1[i][1]) {
      if (bnd == 1) {SWBD(bndval,a[i][nc]); SWBD(bndval,a[i][1]);}
      else SWAP(a[i][1],a[i][nc]);
    }
  }
  return a;
}
