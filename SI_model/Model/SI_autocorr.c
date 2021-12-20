/*
 SI-model with reproduction of individuals and evolvable infectivity.
 To run: ./SI [space] path to folder in which parameter file is and data should be saved
If no path is given, default parameters are used and data is saved in folder: Results/Default.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <errno.h>
#include <math.h>
#include <grace_np.h>
#include <unistd.h>
#include <float.h>
#include <limits.h>
#include <signal.h>
#include <cash2003.h>
#include <cash2.h>
#include <mersenne.h>
#include <cash2-s.h>
#include <complex.h>
#include <fftw3.h>

// Define different types of individuals (Susceptible / Infected)
#define Sus 2	// White
#define Inf 3	// Red

// FUNCTION DECLARATIONS // #SelectionEstimates
// Code for calculating decomposition of selection is labeled with #SelectionEstimates

fftw_complex PairwiseSum_fftwcomplex(fftw_complex*,int,int);
void OutputAutocorr(TYPE2**);

// GLOBAL VARIABLES //

int N;					// To store the total size of the grid (nrow*ncol).

// MODEL PARAMETERS
double reprS = 0.05;			// Reproduction rate susceptibles (per empty neighbour square per time step)
double reprI = 0.;			// Reproduction rate of infected individuals (per empty neighbour square per time step)
double dS = 0.05;			// Death rate of susceptibles
double dI = 0.2;			// Death rate of infected individuals
double h_inf = 20;			// "Half-max" value (of Hill-function) of infection

double mut = 0.005;			// Mutation rate
double mut_step = 0.5;			// Maximum step when mutation occurs (mutation step taken from uniform distr)

double diff = 0.1;			// Diffusion(/swapping) rate of individuals

double init_inf = 5.0;			// Initial infectivity (= infection strength) of Is.

// SIMULATION SETTINGS PARAMETERS
int init = 1;				// 0: Read data from saveplane, 1: All I have infectivity init_inf
int movie = 0;				// 1: movie is recorded (pngs are saved)
int results = 1;			// 1: Saveplanes are made

int timewindow = 2000;			// Time window over which an average of the selection differential is calculated

int startmovietime = 0;			// Sets time at which recording of movie is started
int endmovietime = 100000;		// Sets time at which recording of movie is ended
int movietime = 100;			// Sets resolution of movie: once in so many timesteps a PNG is saved	
int resulttime = 50000;			// Once in so many timesteps results are saved

extern long long genrand_offset;	// Used when restarting simulations from saveplane: offset random number generator to exactly replicate runs

char* foldername;			// Name of folder in which movie and/or results will be saved
char dummy[200];			// Dummy character-array: needed for manipulations on character*s.

static TYPE2** SIeven;			// Simulation planes
static TYPE2** SIodd;

// Autocorrelation and Fourier transform declarations
int maxdist = 100;			// Maximum distance for which autocorrelation is calculated

fftw_complex *occ;			// Matrix to store occupancy (1/0); used in Fourier transformation
fftw_plan FT_plan_occ_f, FT_plan_occ_b;	// Plans for FT of occupancy matrix 

double *autocorr;			// Array to save the pairwise correlation function
double *dist2;				// Array to keep track of no. squares at distance dist^2.


/* MAIN SIMULATION FUNCTIONS */

// Initialise parameters
void Initial(char* savefolder)
{
	//Set default values of Cash parameters (in case these are not changed in parameter file)
	MaxTime = 2000000; 	
	nrow = 512; 
  	ncol = 512;
  	scale = 2;
	margin = 5;
	display = 1; 		//If display = 1, a display is opened. For remote simulations, always set display = 0!
  	boundary = WRAP;	//Use WRAP or FIXED when using Margolus-diffusion
  	ulseedG = 56;
  	boundaryvalue2 = (TYPE2){0,0.,0.,0.,0.};		// StructChange

	char* comp = "nores";
	//Set foldername
	if(strcmp(savefolder,comp))   // If foldername is given as input
	{
		foldername = savefolder;

		//Read parameters
		char* parfile;
		strcpy(dummy, foldername);
		strcat(dummy, "/parameters.txt");
		parfile = malloc(sizeof(dummy));
		strcpy(parfile, dummy);
		ReadParameters(parfile);
	}
	else if(movie == 1 || results == 1)   // If no foldername is given
	{
		foldername = "Results/Default";
		system("mkdir -p Results/Default");
	}

	nplane = 2;
	ndispplane = 1;

	// Basic checks on parameters
	if(8*reprS > 1 || 8*reprI > 1)
	{
		fprintf(stderr, "Error: reproduction contribution of a single individual should be < 1/8.\n");
		exit(1);
	}

	if(nrow != ncol)
	{
		fprintf(stderr, "Error: To use Fourier transforms nrow should be equal to ncol.\n");
		exit(1);
	}
	N = nrow*ncol;

	if(maxdist > nrow/2)
	{
		fprintf(stderr, "Warning: maxdist too large, setting maxdist = nrow/2.\n");
		maxdist = nrow/2;
	}
}

// Initialise simulation planes
void InitialPlane()
{
	int i,j,k,di,dj,dist_sq;

	MakePlane(&SIeven,&SIodd);    // Construct simulation planes

	if(init == 0)
	// Read data from saveplane
	{
		char* saveplanefile;
		strcpy(dummy, foldername);
		strcat(dummy, "/saveplane.txt");
		saveplanefile = malloc(sizeof(dummy));
		strcpy(saveplanefile, dummy);
		ReadSavedData(saveplanefile, SIeven); 
	}

	if(init == 1)
	// Initialise plane at random, all Is get infectivity init_inf
	{
		for(i=1;i<=nrow;i++){
			for(j=1;j<=ncol;j++){
				/* Individuals */		// StructChange
				if(genrand_real1() < 0.5)
				{
					// Infected
					if(genrand_real1()<0.1)
					{
						SIeven[i][j].val = Inf;
						SIeven[i][j].inf = init_inf;
						SIeven[i][j].repr = reprI;
					}
					// Susceptibles
					else
					{
						SIeven[i][j].val = Sus;
						SIeven[i][j].inf = 0.;
						SIeven[i][j].repr = reprS;
					}
					SIeven[i][j].exp_offspring = 0.;
					SIeven[i][j].real_offspring = 0.;
				}
			}
		}
	}

	Boundaries2(SIeven);
	Boundaries2(SIodd);

	// Memory allocation + plan construction for FFT
	occ = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

	FT_plan_occ_f = fftw_plan_dft_2d(nrow,ncol,occ,occ,FFTW_FORWARD,FFTW_MEASURE);
	FT_plan_occ_b = fftw_plan_dft_2d(nrow,ncol,occ,occ,FFTW_BACKWARD,FFTW_MEASURE);

	// Memory allocation for pairwise correlation function
	autocorr = malloc(sizeof(double)*(maxdist*maxdist+1));
	dist2 = malloc(sizeof(double)*(maxdist*maxdist+1));

	// Keep track of the number of neighbouring squares at given squared distance
	for(k=0;k<=(maxdist*maxdist);k++)
		dist2[k] = 0.;
	for(i=0;i<nrow;i++){
		for(j=0;j<ncol;j++){
			di = Min2Ints(i,ncol-i);
			dj = Min2Ints(j,nrow-j);
			dist_sq = di*di + dj*dj;
			if(dist_sq <= (maxdist*maxdist))
				dist2[dist_sq] += 1.;
		}
	}

	// Make autocorr.txt file for output & distances file (to save all the distances for later use)
	if(results == 1)
	{
		char* acfile;
		strcpy(dummy, foldername);
		strcat(dummy, "/autocorr.txt");
		acfile = malloc(sizeof(dummy));
		strcpy(acfile, dummy);

		FILE *ac1;
		ac1 = fopen(acfile, "w");
		fprintf(ac1, "Time");
		for(i=0;i<=(maxdist*maxdist);i++){
			if(dist2[i] > 0.)
				fprintf(ac1,"\tdist_%.4f",sqrt(i));   // Output: Time AC_at_dist1 AC_at_dist2 AC_at_dist3 ...
		}
		fprintf(ac1,"\n");
		fclose(ac1);

		char* dfile;
		strcpy(dummy, foldername);
		strcat(dummy, "/distances.txt");
		dfile = malloc(sizeof(dummy));
		strcpy(dfile, dummy);

		FILE *d1;
		d1 = fopen(dfile, "w");
		for(i=1;i<(maxdist*maxdist);i++){
			if(dist2[i] > 0.)
				fprintf(d1,"%.17g\n",sqrt(i));    // Output dist1 \n dist2 \n dist3 \n ...
		}
		fprintf(d1,"%d",maxdist);
		fclose(d1);
	}
}

// Determine state in next timestep of single square
// Current state is read from plane Now, next state is saved in plane Next
void NextState(TYPE2** Now, TYPE2** Next, int row,int col)
{
	// Reset counters of the expected and real no. offspring for the next timestep.
	Next[row][col].exp_offspring = 0.;
	Next[row][col].real_offspring = 0.;
	Next[row][col].real_trans = 0.;

	// EMPTY SQUARES
	if(Now[row][col].val == 0)
	{
		TYPE2* nei;
		int j;
		double rand, reprprob, reprprobsum;
		reprprob = 0.;

		// Calculate prob that reproduction will take place and contribution to expected no. offspring of the neighbours
		for(j=1;j<=8;j++)
		{
			nei = GetNeighborP(Now,row,col,j);
			reprprob += (nei->repr);
			if(nei->val == Inf)
				(nei->exp_offspring) += (nei->repr);
		}

		// Actual reproduction might take place...
		reprprobsum = 0.;
		if(genrand_real1() < reprprob)
		{
			rand = genrand_real1();
			for(j=1;j<=8;j++)
			{
				nei = GetNeighborP(Now,row,col,j);
				reprprobsum += (nei->repr)/reprprob;
				if(rand < reprprobsum)
				{
					Next[row][col].val = nei->val;
					Next[row][col].inf = nei->inf;
					Next[row][col].repr = nei->repr;

					(nei->real_offspring) += 1.;
					// No need to change nei->real_trans because no mutations take place during reproduction

					break;
				}
			}
		}
		// ... or the square remains empty
		else
		{
			KillInd(Next,row,col);
		}

	}
	// SUSCEPTIBLES
	if(Now[row][col].val == Sus)
	{
		TYPE2* nei;
		int i;
		double rand, inf_tot, infprob, infprobsum;
		inf_tot = 0.;

		// Calculate prob that infection will take place
		for(i=1;i<=8;i++)
			inf_tot += GetNeighborP(Now,row,col,i)->inf;
		infprob = inf_tot / (inf_tot + h_inf);

		// Calculate contribution of potential new infection to exp no. offspring of neighbours
		if(infprob > 0.)
		{
			for(i=1;i<=8;i++)
			{
				nei = GetNeighborP(Now,row,col,i);
				if(nei->val == Inf)
					(nei->exp_offspring) += (1.-dS)*infprob*((nei->inf)/inf_tot);
			}
		}

		// Individual might die...
		if(genrand_real1() < dS)
			KillInd(Next,row,col);
		else
		{
			infprobsum = 0.;
			// ... they might get infected ...
			if(genrand_real1() < infprob)
			{
				rand = genrand_real1();
				for(i=1;i<=8;i++)
				{
					nei = GetNeighborP(Now,row,col,i);
					infprobsum += (nei->inf)/inf_tot;

					if(rand < infprobsum)
					{
						Next[row][col].val = Inf;
						Next[row][col].repr = reprI;
						Next[row][col].inf = nei->inf;	// Inherit infectivity of parent
						(nei->real_offspring) += 1.;
						// No need to change nei->real_trans because no mutations take place during infection
						break;
					}
				}
			}
			// ... or they remain susceptible.
			else
			{
				Next[row][col].val = Now[row][col].val;
				Next[row][col].inf = Now[row][col].inf;
				Next[row][col].repr = Now[row][col].repr;
			}
		}
	}
	// INFECTED INDIVIDUALS
	if(Now[row][col].val == Inf)
	{
		// Calculate contribution of "not dying" to expected no. offspring.
		Now[row][col].exp_offspring += (1.-dI);

		// Individuals might die...
		if(genrand_real1() < dI)
			KillInd(Next,row,col);
		else
		{
			// .. or they stay alive...
			Next[row][col].val = Now[row][col].val;
			Next[row][col].repr = Now[row][col].repr;
			
			// Contribution of "not dying" to real no. offspring.
			Now[row][col].real_offspring += 1.0;

			// ... while their infectivity might mutate...
			if(genrand_real1() < mut){
				Next[row][col].inf = MutateFloat(Now[row][col].inf,2.0*mut_step);
				Now[row][col].real_trans += Next[row][col].inf - Now[row][col].inf;
			}
			// ... or not.
			else
				Next[row][col].inf = Now[row][col].inf;
		}
	}
}

// Update simulation planes and output data to files (called once per timestep)
void Update()
{
	int i, j, modTime;

	modTime = Time%2;
	switch(modTime)
	{
	// EVEN TIMES
	case 0:
		// Display
		if(display==1)
			Display(SIeven);

		// Update -> new states of individuals saved in SIodd
		for(i=1;i<=nrow;i++){
			for(j=1;j<=ncol;j++){
				NextState(SIeven,SIodd,i,j);
			}
		}
		
		// Mobility
		DiffusionBySwap(diff,SIodd);    // Diffusion in the new plane
		//PerfectMix(SIodd);            // Mix the new plane
		
		// Fill occupancy matrix, and use that to calculate autocorrelation
		if(results==1)
			OutputAutocorr(SIeven);

		// Save PNGs for movie (if movie == 1)
		if( (movie==1 && Time%movietime == 0) && (Time >= startmovietime && Time <= endmovietime) )
		{
			char* directory;
			strcpy(dummy, foldername);
			strcat(dummy, "/Film");
			directory = malloc(sizeof(dummy));
			strcpy(directory, dummy);
			DrawSlide2(directory, 1, SIeven);
		}

		break;
		
	// ODD TIMES
	case 1:
		// Display
		if(display==1)
			Display(SIodd);

		// Update -> new states of individuals saved in SIodd
		for(i=1;i<=nrow;i++){
			for(j=1;j<=ncol;j++){
				NextState(SIodd,SIeven,i,j);
			}
		}

		// Mobility
		DiffusionBySwap(diff,SIeven);    // Diffusion in the new plane
		//PerfectMix(SIodd);            // Mix the new plane

		// Fill occupancy matrix, and use that to calculate autocorrelation
		if(results==1)
			OutputAutocorr(SIodd);

		// Save PNGs for movie (if movie == 1)
		if( (movie==1 && Time%movietime == 0) && (Time >= startmovietime && Time <= endmovietime) )
		{
			char* directory;
			strcpy(dummy, foldername);
			strcat(dummy, "/Film");
			directory = malloc(sizeof(dummy));
			strcpy(directory, dummy);
			DrawSlide2(directory, 1, SIodd);
		}

		break;
	}

	// !! Time is increased after execution of the Update function.
}


/* EXTRA FUNCTIONS */

// Death of individual
void KillInd(TYPE2** plane, int row, int col)
{
	plane[row][col].val = 0;
	plane[row][col].repr = 0.;
	plane[row][col].inf = 0.;
}

// Mutation of float value: new value is drawn from a uniform distribution centered around "oldval" with width "step"
// If newval < 0 it is set to 0.
double MutateFloat(double oldval, double step)
{
	double newval;
	
	newval = oldval + (genrand_real1()-0.5)*step;
	if(newval<0.0)
		newval = 0.0;

	return newval;
}

// #SelectionEstimates: Returns the sum of the input array. Uses pairwise/cascade summation algorithm to calculate the sum.
fftw_complex PairwiseSum_fftwcomplex(fftw_complex* a, int sumlength, int startindex)
{
	int k, m, m2;
	fftw_complex sum;
	if(sumlength <= 100)
	{
		sum = a[startindex];
		for(k=(startindex+1); k<(startindex+sumlength); k++)
			sum += a[k];
	}
	else
	{
		m = sumlength / 2;
		m2 = sumlength - m;
		sum = PairwiseSum_fftwcomplex(a,m,startindex) + PairwiseSum_fftwcomplex(a,m2,startindex+m);
	}
	
	return sum;
}

// Fill the occupancy matrix, calculate the autocorrelation and save in autocorr.txt
void OutputAutocorr(TYPE2** plane)
{
	int i,j,di,dj,k,nind,dist_sq;
	nind = 0;
	double correct;

	// Reset the pair correlation function
	for(k=0;k<=(maxdist*maxdist);k++)
		autocorr[k] = 0.;

	// Fill the occupancy matrix
	for(i=1;i<=nrow;i++){
		for(j=1;j<=ncol;j++){
			k = (i-1)*ncol+(j-1);
			if(plane[i][j].val > 0)
			{
				occ[k] = 1.;
				nind++;
			}
			else
				occ[k] = 0.;
		}
	}

	// Calculate pair correlation function - normalised such that autocorrelation = 1 if there is no structure
	fftw_execute(FT_plan_occ_f);
	for(k=0;k<N;k++)
		occ[k] = occ[k] * conj(occ[k]);
	fftw_execute(FT_plan_occ_b);
	correct = 1./((double)nind);
	for(k=0;k<N;k++){
		occ[k] = correct*correct*occ[k];
	}

	// Calculate the pair correlation function
	for(k=0;k<N;k++)
	{
		i = k/ncol;
		j = k%ncol;
		di = Min2Ints(i,ncol-i);
		dj = Min2Ints(j,nrow-j);
		dist_sq = di*di + dj*dj;
		if(dist_sq<=(maxdist*maxdist))
		{
			autocorr[dist_sq] += creal(occ[k]);
		}
	}
	// Normalise to account for non-equal distribution of no. potential neigbours at certain distance
	for(k=0;k<=(maxdist*maxdist);k++)
	{
		if(dist2[k] > 0)
		{
			autocorr[k] = autocorr[k] / dist2[k] ;
		}
	} 

	// Save output
	char* acfile;
	strcpy(dummy, foldername);
	strcat(dummy, "/autocorr.txt");
	acfile = malloc(sizeof(dummy));
	strcpy(acfile, dummy);

	FILE *ac1;
	ac1 = fopen(acfile, "a");
	fprintf(ac1, "%d", Time);
	for(k=0;k<=(maxdist*maxdist);k++)
		if(dist2[k] > 0)
			fprintf(ac1,"\t%.17g",autocorr[k]);
	fprintf(ac1,"\n");
	fclose(ac1);	
} 



//Read parameters from file
void ReadParameters(char* parfile)
{
	FILE *fin;

	fin = fopen(parfile, "r");
	if(fin == NULL)
	{
		fprintf(stderr, "Cannot open file %s \n", parfile);
		exit(1);
	}

	FullInDat(fin, "%lf", "reprS", &reprS, 1);
	FullInDat(fin, "%lf", "reprI", &reprI, 1);
	FullInDat(fin, "%lf", "dS", &dS, 1);
	FullInDat(fin, "%lf", "dI", &dI, 1);
	FullInDat(fin, "%lf", "h_inf", &h_inf, 1);

	FullInDat(fin, "%lf", "mut", &mut, 1);
	FullInDat(fin, "%lf", "mut_step", &mut_step, 1);

	FullInDat(fin, "%lf", "diff", &diff, 1);

	FullInDat(fin, "%lf", "init_inf", &init_inf, 1);

	FullInDat(fin, "%d", "init", &init, 1);
	FullInDat(fin, "%d", "movie", &movie, 1);
	FullInDat(fin, "%d", "results", &results, 1);
	FullInDat(fin, "%d", "display", &display, 1);

	FullInDat(fin, "%d", "startmovietime", &startmovietime, 1);
	FullInDat(fin, "%d", "endmovietime", &endmovietime, 1);
	FullInDat(fin, "%d", "movietime", &movietime, 1);
	FullInDat(fin, "%d", "resulttime", &resulttime, 1);

	FullInDat(fin, "%d", "MaxTime", &MaxTime, 1);
	FullInDat(fin, "%d", "nrow", &nrow, 1);
	FullInDat(fin, "%d", "ncol", &ncol, 1);
	FullInDat(fin, "%d", "ulseedG", &ulseedG, 1);

	fclose(fin);
}



