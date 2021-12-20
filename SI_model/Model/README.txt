README of SI-model with calculations of spatial decomposition of selection differential (S = S_local + S_interlocal)


****************
* Dependencies *
****************

Next to some standard C-libraries, the SI model depends on two libraries:

FFTW : Fast Fourier Transform. Install from within Cash-directory. Code available in Cash-directory, latest versions and manual available from http://www.fftw.org.

CASH : Utrecht-developed library for CA simulations. Code included (see .. directory and README there)


******************
* How to compile *
******************

Make sure first FFTW and cash2 have been properly compiled (see http://www.fftw.org and README in .. directory).

Then compile model using makefile in Model directory
( make -B SI )


**************
* How to run *
**************

./SI to run with default settings

./SI <folder-name> to run with parameters read from parameters.txt in <folder-name>. 
Output is also stored in <folder-name>.
See example in Test/


*****************
* Code versions *
*****************

SI : Full functionality. Simulate SI model, calculate global population characteristics and spatial decomposition of selection
SI_basic : Stripped of multiscale selection decomposition. Simulate SI model and calculate global characteristics only.
SI_autocorr : Calculate and output autocorrelation in a given plane. Run on saved simulation plane of SI or SI_basic.


******************
* Code structure *
******************

This model uses libraries cash2 and an altered version of cash2-s. Below, we describe some of the oddities of the code structure:

Each square on the grid contains a TYPE2-struct, which carries the characteristics of a square (= potentially an individual).
Currently, the struct contains:

  int val;	// 0, S or I, describes type of individual
  double inf;	// Infectivity of individual in square
  double repr;	// Reproduction rate of individual in square

  double real_offspring; // Realised no. offspring from previous timestep (including self)
  double exp_offspring;	 // Expected no. offspring in the current timestep
  double real_trans;	 // Transmission bias, to measure the (realised) second term of the Price equation

The TYPE2-struct is defined in cash2.h.
When changes to this struct are made, a few functions in cash2-s.c should be altered. These are flagged with "StructChange" (CNTRL-F StructChange in cash2-s.c does the trick). 
After alterations to cash2- or cash2-s-files, the libraries should be rebuilt (make; make install in folder with cash2 files). 
Then, the model should be rebuilt as well (cd Model; make -B QSbact).

The main functions defining the model are found in SI.c. These are (see also manual_cash_s.pdf):

Initial		Provide initial values for the parameters
InitialPlane	Initialise simulation planes
NextState	Define update rules for a single square.
Update		Contains all calls to functions that should be made once per timestep, plus calls NextState for all squares.


main() cannot be found in SI.c, but in cash2-s.c. This function does not have to be altered often. It calls the functions in SI.c. It's basic structure is:
Initial
InitialPlane
for Time<MaxTime:
	Update

The simulation uses two simulation planes: SIodd and SIeven. Whenever a plane is updated (i.e. in a timestep), data is read from one simulation plane and the next state is saved in the other. We thus keep track of the current and previous state of the simulation. This feature is used to calculate realised fitness of individuals.

Next to the four main functions described above, SI.c contains several auxiliary functions that are called during simulation time steps. Most function names are self-explanatory. A short description of their functionality is given in SI.c.


**************
* Parameters *
**************

MODEL PARAMETERS

reprS 			Reproduction rate of susceptibles (per empty neighbour site per time step)
reprI 0.0		Reproduction rate of infected individuals (per empty neighbour site per time step)
dS 0.05			Death rate of susceptibles (per time step)
dI 0.2			Death rate of infected individuals (per time step)
h_inf 20.0		Scaling factor of infection probability (if the sum of infectivities of the eight neighbours is equal to h_inf, the probability of infection is 0.5)

mut 0.005		Mutation rate (per time step)
mut_step 0.5		Maximal mutation step size (mutation step size is chosen from uniform distribution on [-mut_step; mut_step]. If new value < 0, it is set to 0.

diff 0.05		"Diffusion rate": probability with which a square is swapped with a random neighbour

init_inf 5.0		Initial value of infectivity

SIMULATION SETTING PARAMETERS

init 1			0: Data is read from previous file stored in saveplane.txt; 1: random initialisation, all infected individuals get infectivity = init_inf
movie 0			0: No movie is recorded; 1: PNGs are saved to construct movie.
results 1		0: No regular data-dumps; 1: regular data-dumps are made.
display 0		0: No display (for remote simulation); 2: real-time display of simulation plane.

startmovietime 200000	Time point from which PNGs for movie are saved.
endmovietime 205000	Time point after which no more PNGs are saved.
movietime 5		Interval at which PNGs are saved.
resulttime 100000	Interval at which data-dumps are made (stored in saveplaneTIME.txt)
startseltime 200000	Time point from which selection decomposition is calculated and outputted in selection.txt (every time step)
endseltime 210000	Time point after which selection decomposition is no longer calculated

MaxTime 250000		Maximum number of time steps in simulation
nrow 1024		Number of rows of simulation planes
ncol 1024		Number of columns of simulation planes
ulseedG 39		Random seed




