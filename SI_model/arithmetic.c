#include <stdio.h>
#include "cash2003.h"

extern int nrow, ncol, first, last;

TYPE **Sum(TYPE **a,TYPE  **b,TYPE  **c)
{
  int i,j;
  /*$dir force_vector*/
  for (i=first,j=last; i <= j; i++)
    a[0][i] = b[0][i] + c[0][i];
  return a;
}

TYPE **SumV(TYPE **a,TYPE  **b,TYPE  c)
{
  int i,j;
  /*$dir force_vector*/
  for (i=first,j=last; i <= j; i++)
    a[0][i] = b[0][i] + c;
  return a;
}

TYPE **Minus(TYPE **a,TYPE  **b,TYPE  **c)
{
  int i,j;
  /*$dir force_vector*/
  for (i=first,j=last; i <= j; i++)
    a[0][i] = b[0][i] - c[0][i];
  return a;
}

TYPE **MinusV(TYPE **a,TYPE  **b,TYPE  c)
{
  int i,j;
  /*$dir force_vector*/
  for (i=first,j=last; i <= j; i++)
    a[0][i] = b[0][i] - c;
  return a;
}

TYPE **Div(TYPE **a,TYPE  **b,TYPE  **c)
{
  int i,j;
  /*$dir force_vector*/
  for (i=first,j=last; i <= j; i++)
    a[0][i] = b[0][i] / c[0][i];
  return a;
}

TYPE **DivV(TYPE **a,TYPE  **b,TYPE  c)
{
  int i,j;
  /*$dir force_vector*/
  for (i=first,j=last; i <= j; i++)
    a[0][i] = b[0][i] / c;
  return a;
}

TYPE **Mod(TYPE **a,TYPE  **b,TYPE  **c)
{
  int i,j;
  /*$dir force_vector*/
  for (i=first,j=last; i <= j; i++)
    a[0][i] = b[0][i] % c[0][i];
  return a;
}

TYPE **ModV(TYPE **a,TYPE  **b,TYPE  c)
{
  int i,j;
  /*$dir force_vector*/
  for (i=first,j=last; i <= j; i++)
    a[0][i] = b[0][i] % c;
  return a;
}

TYPE **Mult(TYPE **a,TYPE  **b,TYPE  **c)
{
  int i,j;
  /*$dir force_vector*/
  for (i=first,j=last; i <= j; i++)
    a[0][i] = b[0][i] * c[0][i];
  return a;
}

TYPE **MultV(TYPE **a,TYPE  **b,TYPE  c)
{
  int i,j;
  /*$dir force_vector*/
  for (i=first,j=last; i <= j; i++)
    a[0][i] = b[0][i] * c;
  return a;
}

TYPE **MultF(TYPE **a,TYPE  **b,float c)
{
  int i,j;
  /*$dir force_vector*/
  for (i=first,j=last; i <= j; i++)
    a[0][i] = b[0][i] * c;
  return a;
}

TYPE **BinSumOld(TYPE **a,TYPE  **b[], int c)
{
  int i,j,bin;
  Copy(a,b[0]);
  for (bin=1; bin<c; bin++)
    /*$dir force_vector*/
    for (i=first,j=last; i <= j; i++)
      a[0][i] += b[bin][0][i] << bin; 
  return a;
}

TYPE **BinSum(TYPE **a,int n,TYPE **b0,TYPE **b1,TYPE **b2,TYPE **b3,TYPE **b4,TYPE **b5,TYPE **b6,TYPE **b7)
{
  int i,j;/*,bin;*/
	if (n < 1 || n > 8) {
	  printf("Error in BinSum: n=%d\n",n);
		n = 0;
  }
  if (n) Copy(a,b0);
	if (n > 1)
    /*$dir force_vector*/
    for (i=first,j=last; i <= j; i++)
      a[0][i] += b1[0][i] << 1; 
	if (n > 2)
	  /*$dir force_vector*/
    for (i=first,j=last; i <= j; i++)
      a[0][i] += b2[0][i] << 2; 
	if (n > 3)
	  /*$dir force_vector*/
    for (i=first,j=last; i <= j; i++)
      a[0][i] += b3[0][i] << 3; 
	if (n > 4)
	    /*$dir force_vector*/
    for (i=first,j=last; i <= j; i++)
      a[0][i] += b4[0][i] << 4; 
	if (n > 5)
	    /*$dir force_vector*/
    for (i=first,j=last; i <= j; i++)
      a[0][i] += b5[0][i] << 5; 
	if (n > 6)
	    /*$dir force_vector*/
    for (i=first,j=last; i <= j; i++)
      a[0][i] += b6[0][i] << 6; 
	if (n > 7)
	    /*$dir force_vector*/
    for (i=first,j=last; i <= j; i++)
      a[0][i] += b7[0][i] << 7; 
  return a;
}

TYPE **RollRight(TYPE **a,TYPE **b,int c)
{
  int i,j;
  /*$dir force_vector*/
  for (i=first,j=last; i <= j; i++)
    a[0][i] = b[0][i] >> c;
  return a;
}

TYPE **RollLeft(TYPE **a,TYPE **b,int c)
{
  int i,j;
  /*$dir force_vector*/
  for (i=first,j=last; i <= j; i++)
    a[0][i] = b[0][i] << c;
  return a;
}

TYPE **GetBits(TYPE **a,TYPE **b,int f,int l)
{
  int i,j;
  unsigned mask;
  mask = ~(~0 << (l-f+1));
  if (f > l) fprintf(stderr,"Error in GetBits: first > last");
  /*$dir force_vector*/
  for (i=first,j=last; i <= j; i++)
    a[0][i] = (b[0][i] >> f) & mask;
  return a;
}

TYPE **PutBits(TYPE **a,int p,int v)
{
  int i,j;
  unsigned mask;
  mask = v << p;
  /*$dir force_vector*/
  for (i=first,j=last; i <= j; i++)
    a[0][i] = a[0][i] | mask;
  return a;
}

TYPE **Hamming(TYPE **a,TYPE  **b,TYPE  **c)
{
  int i,j,k;
  unsigned n;
  /*$dir force_vector*/
  for (i=first,j=last; i <= j; i++) {
    n = b[0][i] ^ c[0][i];
    for (k=0; n != 0; n >>= 1)
      if (n & 01) k++;
    a[0][i] = k;
  }
  return a;
}
