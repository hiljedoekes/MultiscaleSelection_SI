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

double init_inf = 4.0;			// Initial infectivity (= infection strength) of Is.

// VARIABLES TO STORE GLOBAL POPULATION CHARACTERISTICS
double mean_delta_inf = 0.;   // Change in mean infectivity (per timestep)
double mean_S_exp = 0.;       // Mean selection differential, calculated based on expected fitness
double mean_S_real = 0.;      // Mean selection differential, calculated based on realised fitness
double mean_mean_real_trans = 0.;  // Mean transmission bias, calculated based on realised fitness
int timecounter = 0;			// Not really necessary, but just to make sure that I don't make mistakes in calculations

// SIMULATION SETTINGS PARAMETERS
int init = 1;				// 0: Read data from saveplane, 1: All I have infectivity init_inf
int movie = 0;				// 1: movie is recorded (pngs are saved)
int results = 0;			// 1: Saveplanes (data-dumps) are made

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


/* MAIN SIMULATION FUNCTIONS */

// Initialise parameters
void Initial(char* savefolder)
{
	//Set default values of Cash parameters (in case these are not changed in parameter file)
	MaxTime = 2000000; 	
	nrow = 256; 
  	ncol = 256;
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
		fprintf(stderr, "Error: To use Fourier transform nrow should be equal to ncol.\n");
		exit(1);
	}
	
	N = nrow*ncol;	
}

// Initialise simulation planes
void InitialPlane()
{
	int i,j;	// To use in iterations later on

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

	// Make countglobal.txt file to store basic output in
	if(results == 1)
	{
		char* countfile;
		strcpy(dummy, foldername);
		strcat(dummy, "/countglobal.txt");
		countfile = malloc(sizeof(dummy));
		strcpy(countfile, dummy);

		FILE *count1;
		count1 = fopen(countfile, "w");
		fprintf(count1, "Time \tSus \tInfected \tmean_inf \n");   // Output: Time number_of_susceptibles number_of_infected mean_infectivity
		fclose(count1);
	}

	// Make globalPrice.txt file to store global selection estimates
	if(results == 1)
	{
		char* pricefile;
		strcpy(dummy, foldername);
		strcat(dummy, "/globalPrice.txt");
		pricefile = malloc(sizeof(dummy));
		strcpy(pricefile, dummy);

		FILE *price1;
		price1 = fopen(pricefile, "w");
		fprintf(price1, "Time \tmean_delta_inf \tmean_S_exp_offspring \tmean_S_real_offspring \tmean_mean_real_transmission \n");   // Output: Time mean_change_in_inf mean_S_based_on_expected_fitness mean_S_based_on_realised_fitness mean_transmission_bias
		fclose(price1);
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
	// Infected individuals
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

		// SavePlanes
		if(results == 1 && Time%resulttime == 0)
		{
			//Data dump in saveplaneTIME.txt
			char* planefile;
			char planefilename[100];
			memset(planefilename,0,sizeof(planefilename));
			char savenumber[20];

			sprintf(savenumber,"%d",Time);
			strcat(planefilename, "/saveplane");
			strcat(planefilename, savenumber);
			strcat(planefilename, ".txt");
			strcpy(dummy, foldername);
			strcat(dummy, planefilename);
			planefile = malloc(sizeof(dummy));
			strcpy(planefile, dummy);
			SavePlane(planefile, SIeven);
		}

		// Update -> new states of individuals saved in SIodd
		for(i=1;i<=nrow;i++){
			for(j=1;j<=ncol;j++){
				NextState(SIeven,SIodd,i,j);
			}
		}

		// Mobility		
		DiffusionBySwap(diff,SIodd);    // Diffusion in the new plane
		//PerfectMix(SIodd);            // Mix the new plane

		// Count total no. individuals and calculate means and total selection, output in countglobal.txt and globalPrice.txt at certain time steps
		if(results == 1){
			OutputGlobal(SIeven,SIodd);
		}

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
		
		// Count total no. individuals and calculate means and total selection, output in countglobal.txt and globalPrice.txt at certain time steps
		if(results == 1){
			OutputGlobal(SIodd,SIeven);
		}

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

//Output global characteristics (mean_inf, total_S etc)
void OutputGlobal(TYPE2** Now, TYPE2** Next)
{
	int x, y, nSus, nInf, nInf_next;
	double mean_inf, mean_inf_next, mean_inf_sq, mean_exp_offspring, mean_real_offspring, S_exp, mean_dx, mean_dy_exp, S_real, mean_dy_real, Cov_phi2_w, mean_dx_sq, dx, dy, delta_inf, mean_real_trans;

	// Save globalprice data
	if(Time>0 && Time%timewindow == 0)
	{
		// Calculate means
		mean_delta_inf = mean_delta_inf / timecounter;
		mean_S_exp = mean_S_exp / timecounter;
		mean_S_real = mean_S_real / timecounter;
		mean_mean_real_trans = mean_mean_real_trans / timecounter;

		// Save output in globalPrice.txt
		char* pricefile;
		strcpy(dummy, foldername);
		strcat(dummy, "/globalPrice.txt");
		pricefile = malloc(sizeof(dummy));
		strcpy(pricefile, dummy);

		FILE *price1;
		price1 = fopen(pricefile, "a");
		fprintf(price1, "%d \t%.17g \t%.17g \t%.17g \t%.17g \n", Time, mean_delta_inf, mean_S_exp, mean_S_real, mean_mean_real_trans);
		fclose(price1);

		// Reset
		mean_delta_inf = 0.;
		mean_S_exp = 0.;
		mean_S_real = 0.;
		mean_mean_real_trans = 0.;
		timecounter = 0;
	}

	nSus = 0;
	nInf = 0;
	nInf_next = 0;
	mean_inf = 0.;
	mean_inf_next = 0.;
	mean_inf_sq = 0.;
	mean_exp_offspring = 0.;
	mean_real_offspring = 0.;
	S_exp = 0.;
	mean_dx = 0.;
	mean_dy_exp = 0.;
	S_real = 0.;
	mean_dy_real = 0.;
	Cov_phi2_w = 0.;
	mean_dx_sq = 0.;
	mean_real_trans = 0.;

	for(x=1;x<=nrow;x++){
		for(y=1;y<=ncol;y++){
			if(Now[x][y].val == Sus)
				nSus++;
			if(Now[x][y].val == Inf){
				nInf++;
				mean_inf += Now[x][y].inf;
				mean_inf_sq += ( Now[x][y].inf * Now[x][y].inf );
				mean_exp_offspring += Now[x][y].exp_offspring;
				mean_real_offspring += Now[x][y].real_offspring;
			}
			if(Next[x][y].val == Inf){
				nInf_next++;
				mean_inf_next += Next[x][y].inf;
			}
		}
	}
	if(nInf == 0)
	{
		fprintf(stderr,"Infection died out.\n");
		exit(1);
	}
	else
	{
		mean_inf = mean_inf / nInf;
		mean_inf_sq = mean_inf_sq / nInf;
		mean_inf_next = mean_inf_next / nInf_next;
		delta_inf = mean_inf_next - mean_inf;
		mean_exp_offspring = mean_exp_offspring / nInf;
		mean_real_offspring = mean_real_offspring / nInf;
		// Calculate selection differential - covariance between rel fitness and inf.
		// Use naive algorithm on the residuals
		// Also calculate the transmission term (2nd term Price equation), and purifying selection
		for(x=1;x<=nrow;x++){
			for(y=1;y<=ncol;y++){
				if(Now[x][y].val == Inf){
					dx = Now[x][y].inf - mean_inf;
					dy = Now[x][y].exp_offspring/mean_exp_offspring - 1;
					mean_dx += dx;
					mean_dy_exp += dy;
					S_exp += dx*dy;
					dy = Now[x][y].real_offspring/mean_real_offspring - 1;
					mean_dy_real += dy;
					S_real += dx*dy;
					//S_exp += (SIeven[x][y].inf - mean_inf)*(SIeven[x][y].exp_offspring - mean_exp_offspring);
					//S_real += (SIeven[x][y].inf - mean_inf)*(SIeven[x][y].real_offspring - mean_real_offspring);
					dx = (Now[x][y].inf * Now[x][y].inf) - mean_inf_sq;
					mean_dx_sq += dx;
					Cov_phi2_w += dx*dy;
					
					if(Now[x][y].real_offspring > 0){
						mean_real_trans += Now[x][y].real_trans;
					}
				}
			}
		}
		mean_dx = mean_dx / nInf;
		mean_dx_sq = mean_dx_sq / nInf;
		mean_dy_exp = mean_dy_exp / nInf;
		mean_dy_real = mean_dy_real / nInf;
		mean_real_trans = mean_real_trans / (mean_real_offspring * nInf);

		S_exp = S_exp / nInf - mean_dx*mean_dy_exp;
		S_real = S_real / nInf - mean_dx*mean_dy_real;
		Cov_phi2_w = Cov_phi2_w / nInf - mean_dx_sq*mean_dy_real;

		// To calculate means
		mean_delta_inf += delta_inf;
		mean_S_exp += S_exp;
		mean_S_real += S_real;
		mean_mean_real_trans += mean_real_trans;
		timecounter++;

		// Save count-output in countglobal.txt
		if(Time%100 == 0)
		{		
			char* countfile;
			strcpy(dummy, foldername);
			strcat(dummy, "/countglobal.txt");
			countfile = malloc(sizeof(dummy));
			strcpy(countfile, dummy);

			FILE *count1;
			count1 = fopen(countfile, "a");
			fprintf(count1, "%d \t%d \t%d \t%.17g \n", Time, nSus, nInf, mean_inf);
			fclose(count1);
		}
	}
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



