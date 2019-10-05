#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#include "consts.h"
#include "globals.h"
#include "rocks.h"
#include "utils.h"
#include "conditions.h"


// Global variables to be share with other files.

// Holds cliff heights -- will change through time, eventually...
double **topography;

// Special SWAN matrices.
double **shelf_depth = NULL; // SWAN bathymetry
double **wave_h_sig = NULL; // SWAN wave heights.
double **wave_dir = NULL;  // SWAN wave angles.

int current_time_step = 0;  // Time step of current calculation
double current_time = 0.;  // Current model time
double stop_after = 36500; // Stop after what number of time steps
double time_step = 1; // days; reflects rate of sediment transport per time step


// Global variables to be used only within this file.
static const double kAngleFactor = 1.;

// Depth array
double **cell_depth;

// Overall Shoreface Configuration Arrays - Data file information

// Flag indicating of cell is entirely beach
char **AllBeach;
// Flag indicating if cell is entirely rock LMV
char **AllRock = NULL;
// Fractional amount of cell full of sediment LMV
double **PercentFullSand = NULL;
// Fractional amount of a cell full of rock LMV
double **PercentFullRock = NULL;
// Array to control weathering rates of rock along the beach LMV
char **type_of_rock = NULL;
// Age since cell was deposited
static int Age[X_MAX][2 * Y_MAX];

// A sink is a cell that is routinely emptied (if it is on the beach)
static const int kSinkY[] = {158, 361, 158, 361, 158, 361};
static const int kSinkX[] = {190, 190, 189, 189, 191, 191};

static FILE *SaveSandFile;
static FILE *WaveFileOutput;
static FILE *InitMetaFile;
static FILE *ReadSandFile;
static FILE *ReadWaveFile;
static FILE *ReadRealWaveFile;
static FILE *ReadControlFile;

static char savefilename[2048] = "CEM";
static char readfilename[2048] = "CEM_3285.out";
// Input Wave Distribution file: no = 0, binned file = 1,
// angle/period/height file = 2
// #define WAVE_IN (0)
static char readwavename[2048] = "In_WaveData.dat";
// use a file to initialise run? Over rides setup data
// #define INITIALIZE_FILE ('n')
static char readcontrolname[2048] = "In_CEM_init.dat";
// Create a metadata file?
// #define METADATA ('n')
static char metasavename[2048] = "Metadata.out";
// Create a wavedata file?
// #define WAVE_DATA ('n')
static char wavesavename[2048] = "Wavedata.out";

// Computational Arrays (determined for each time step)  -- will eventually be
// in a structure for BMI

static int X[MaxBeachLength]; // X Position of ith beach element
static int Y[MaxBeachLength]; // Y Position of ith beach element
static int XBehind[MaxBeachLength]; // Cell that is "behind" X[i] LMV
static int YBehind[MaxBeachLength]; // Cell that is "behind" Y[i] LMV
static int XRock[MaxBeachLength]; // X Position of ith rock element LMV
static int YRock[MaxBeachLength]; // Y Position of ith rock element LMV
static int XRockBehind[MaxBeachLength]; // Cell that is "behind" XRock[i] LMV
static int YRockBehind[MaxBeachLength]; // Cell that is "behind" YRock[i] LMV
static char InShadow[MaxBeachLength]; // Is ith beach element in shadow?
static double ShorelineAngle[MaxBeachLength]; // Angle between cell and right (z+1) neighbor
static double SurroundingAngle[MaxBeachLength]; // Angle between left and right neighbors
static char UpWind[MaxBeachLength]; // Upwind or downwind condition used to calculate sediment transport
static double VolumeIn[MaxBeachLength]; // Sediment volume into ith beach element
static double VolumeOut[MaxBeachLength]; // Sediment volume out of ith beach element
static double VolumeAcrossBorder[MaxBeachLength]; // Sediment volume across border of ith beach element in m^3/day LMV
/* amount sediment needed, not necessarily amount a cell gets */
static double ActualVolumeAcross[MaxBeachLength]; // Sediment volume that actually gets passed across border LMV
/* amount sed is limited by volumes across border upstream and downstream */
static char DirectionAcrossBorder[MaxBeachLength]; // Flag to indicate if transport across border is L or R LMV
static char FlowThroughCell[MaxBeachLength]; // Is flow through ith cell Left, Right, Convergent, or Divergent LMV
static double DistanceToBeach[MaxBeachLength]; // Distance in meters from rock to beach LMV
static double MinDistToBeach[MaxBeachLength]; // From a rock cell j, min distance (in meters) to closest beach LMV
static int ClosestBeach[MaxBeachLength]; // i position of closest rock to beach LMV
static double AmountWeathered[MaxBeachLength]; // Amount of rock weathered from rock cell j LMV

#if defined(WITH_SWAN)
static const char kSwanFlag = 'y'; // Is SWAN doing wave transformations?
#else
static const char kSwanFlag = 'n';
#endif

static double BreakDepth; // Breaking wave depth found from SWAN run

static double EvaluateAngle; // Temporary angle holder for the ConvertAngle function

// for temporary debugging only, 5-5-14
// UPDATE 11/20/14 -- now using these for upwind scheme fixing, so keep 'em
// around (probably should rename...)
static double Hsigdebug[MaxBeachLength];
static double Dirdebug[MaxBeachLength];
static double Hddebug[MaxBeachLength];
static double xdebug[MaxBeachLength];
static double ydebug[MaxBeachLength];
static double WvHeight;
static double Angle;

// Miscellaneous Global Variables

static int NextX; // Global variables used to iterate FindNextCell in global array -
static int NextY; // would've used pointer but wouldn't work
static int BehindX;
static int BehindY;
static int BehindRockX;
static int BehindRockY;
static int NextRockX; // Global variables used to iterate FindNextRockCell in global array LMV
static int NextRockY;
static int TotalBeachCells; // Number of cells describing beach at particular iteration
static int TotalRockCells; // Number of cells describing rock at an iteration LMV
static int ShadowXMax; // used to determine maximum extent of beach cells
static double WaveAngle; // wave angle for current time step
static int FindStart; // Used to tell FindBeach at what Y value to start looking
static int FindRockStart; // Used to tell FindRock at what Y value to start looking LMV
static char FellOffArray; // Flag used to determine if accidentally went off array
static char FellOffRockArray; // Flag used to determine if accidentally went off rock array LMV
static double MassInitial; // For conservation of mass calcs
static double MassCurrent;
static int NumWaveBins; // For Input Wave - number of bins
static double WaveMax[36]; // Max Angle for specific bin
static double WaveProb[36]; // Probability of Certain Bin
static double WaveAngleIn;
static double WaveHeightIn;
static double WavePeriodIn;

static char StartFromFile = 'n'; // start from saved file?

static double Period = 10.0; // seconds
static double OffShoreWvHt = 2.0; // meters

// Fractional portion of waves coming from positive (left) direction
static double Asym = 0.70;
// .5 = even dist, < .5 high angle domination. NOTE Highness actually
// determines lowness!
static double Highness = 0.35;
// Number of time steps calculations loop at same wave angle
static double Duration = 1.;

// Weathering rate of slow rock e.g. 0.5 = slow weathering, 1/2 as fast LMV
static double SlowWeatherCoeff = 0. * NORMAL_WEATHERING_RATE;
// Weathering rate of fast rock e.g. 2 = fast weathering, 2 times faster
// than normal LMV
static double FastWeatherCoeff = 1 * NORMAL_WEATHERING_RATE;
// Percent of fast weathering rock lost because it is too fine to stay
// in nearshore LMV
static double PercentFineFast = 0.0;
// Percent of slow weathering rock lost because it is too fine to stay in
// nearshore LMV
static double PercentFineSlow = 0.0;
// Amount of erosion to TotalBeachCells in m/year LMV
static double ErosionRatePerYear = 0.5;

// Angle (deg) used to align coastline with real wave climate
static double coastrotation = 40.;
// Factor applied to the input wave height; 1 = no change, 0.5 = half
// wave height, and 2 = double wave height
static double waveheightchange = 1.;
// Factor applied to the input wave period; 1 = no change, 0.5 = half
// wave period, and 2 = double wave period
static double waveperiodchange = 1.;

static int OWflag = 0; // A debugging flag for overwash routines

void AdjustShore(int i);
void AgeCells(void);
void ControlFile(void);
float ConvertAngle(float EvaluateAngle, int type);
void DetermineAngles(void);
void DetermineSedTransport(void);
void DoSink(void);
void Transport(void);
void ErodeTheBeach(int i);
void FindBeachCells(int YStart);
void FindRockCells(int YStart);
char FindIfInShadow(int xin, int yin, int ShadMax);
void FindNearestBeach(int j);
void FindNextCell(int x, int y, int z);
void FindNextRockCell(int x, int y, int z);
float FindWaveAngle(void);
void FixBeach(void);
void FixFlow(void);
void FlowInCell(void);
int getIndex(int x, int y, char interface);
float MassCount(void);
void OopsImEmpty(int x, int y);
void OopsImFull(int x, int y);
void ParseSWAN(int ShoreAngleLoc, float ShoreAngle);
void PauseRun(int x, int y, int in);
void periodic_boundary_copy(void);
float Raise(float b, float e);
float RandZeroToOne(void);
void ReadSandFromFile(void);
void ReadWaveIn(void);
void RealWaveIn(void);
void RockCalculations(void);
void SaveSandToFile(void);
void SaveLineToFile(void);
void SedTrans(int From, int To, float ShoreAngle, char MaxT, int Last);
void ShadowSweep(void);
void TransportSedimentSweep(void);
void WaveOutFile(void);
void WeatherRock(int j);
int XMaxBeach(int Max);
void ZeroVars(void);

/* Overwash function prototypes */
void CheckOverwashSweep(void);
float GetOverwashDepth(int xin, int yin, float xinfl, float yinfl, int ishore);
void CheckOverwash(int icheck);
void DoOverwash(int xfrom, int yfrom, int xto, int yto, float xintto,
                float yintto, float widthin, int ishore);

void Delay(void);
void PrintLocalConds(int x, int y, int in);

int cem_initialize(void);
int cem_update(void);
int cem_finalize(void);

int cem_initialize(void) {
  ShadowXMax = X_MAX; /* RCL: was X_MAX-5; */

  if (SEED == -999)
    srand(time(NULL));
  else
    srand(SEED);

  { // Allocate memory for arrays.
    const int n_rows = X_MAX;
    const int n_cols = 2 * Y_MAX;

    topography = (double**)malloc2d(n_rows, n_cols, sizeof(double));
    AllBeach = (char**)malloc2d(n_rows, n_cols, sizeof(char));
    AllRock = (char**)malloc2d(n_rows, n_cols, sizeof(char));
    PercentFullSand = (double**)malloc2d(n_rows, n_cols, sizeof(double));
    PercentFullRock = (double**)malloc2d(n_rows, n_cols, sizeof(double));
    cell_depth = (double**)malloc2d(n_rows, n_cols, sizeof(double));
    shelf_depth = (double**)malloc2d(n_rows, n_cols, sizeof(double));
    wave_h_sig = (double**)malloc2d(n_rows, n_cols, sizeof(double));
    wave_dir = (double**)malloc2d(n_rows, n_cols, sizeof(double));
    type_of_rock = (char**)malloc2d(n_rows, n_cols, sizeof(char));
  }

  /* Start from file or not? */
  if (PROMPT_START == 'y') {
    printf("shall we start from a file (y or n)? \n");
    scanf("%c", &StartFromFile);

    if (StartFromFile == 'y') {
      printf("Starting Filename? \n");
      scanf("%24s", readfilename);
      printf("Saving Filename? \n");
      scanf("%s", savefilename);
      printf("What time step are we starting at?");
      scanf("%d", &current_time_step);
      ReadSandFromFile();
    }

    if (StartFromFile == 'n') {
      printf("Saving Filename? \n");
      scanf("%s", savefilename);
      InitConds(cell_depth, AllBeach, AllRock, PercentFullSand,  PercentFullRock, type_of_rock, topography);
      printf("InitConds OK \n");
    }
  }

  else if (StartFromFile == 'y') {
    ReadSandFromFile();
  }

  else {
    InitConds(cell_depth, AllBeach, AllRock, PercentFullSand,  PercentFullRock, type_of_rock, topography);
  }

  /* Set Periodic Boundary Conditions */
  periodic_boundary_copy();

  /* Count Initial Mass */
  MassInitial = MassCount();

  /* Fix the Beach (?) */
  FixBeach();

  /* PWL commented out to match ndelta4.c, 12/3/14 */
  /*if (SAVE_LINE) SaveLineToFile();
     if (SAVE_FILE) SaveSandToFile(); */

  /*Read in wave data? */
  if (WAVE_IN == 1) ReadWaveIn(); /*Initialise WaveMax and WaveProb */

  if (INITIALIZE_FILE == 'y')
    ControlFile(); /*Read in initialisation data to control model */

  return 0;
}

// Update the CEM by a single time step.
int cem_update(void) {
  // duration loop variable, moved from former (now deleted) 'main.c' above
  int xx;

  {
    /*  Time Step iteration - compute same wave angle for Duration time steps */

    /*  Calculate Wave Angle */
    WaveAngle = FindWaveAngle();

    /* output wave data */
    if (WAVE_DATA == 'y') WaveOutFile();

    /*  Loop for Duration at the current wave sign and wave angle */
    for (xx = 0; xx < Duration; xx++) {
      MassCurrent = MassCount();

      /* Text to Screen? */
      if (current_time_step % SCREEN_TEXT_SPACING == 0) {
        printf("==== WaveAngle: %2.2f  MASS Percent: %1.4f  Time Step: %d\n",
               180 * (WaveAngle) / PI, MassCurrent / MassInitial,
               current_time_step);
      }

      periodic_boundary_copy(); /* Copy visible data to external model - creates
                                 boundaries */

      ZeroVars();

      /* Initialize for Find Beach Cells  (make sure strange beach does not
       * cause trouble */
      FellOffArray = 'y';
      FindStart = 1;
      FellOffRockArray = 'y';
      /*LMV*/ FindRockStart = 1;
      /*LMV*/
      /*  Look for beach - if you fall off of array, bump over a little and try
       * again */
      while (FellOffArray == 'y') {
        FindBeachCells(FindStart);
        FindStart += FindCellError;
        if (FellOffArray == 'y') {
          /*printf("NOODLE  !!!!!FoundCells: %d GetO = %c \n",
           * FindStart,FellOffArray); */
          /*PauseRun(1,1,-1); */
        }

        /* Get Out if no good beach spots exist - finish program */
        if (FindStart > Y_MAX / 2 + 1) {
          printf("Stopped Finding Beach - done %d %d", FindStart, Y_MAX / 2 - 5);
          SaveSandToFile();
          getchar();
          return 1;
        }
      }

      /* LMV */
	  if (INIT_ROCK > 0)
	  {
		  while (FellOffRockArray == 'y') {
			  FindRockCells(FindRockStart);
			  /*printf("FoundRockCells: %d GetO = %c \n",
			   * FindRockStart,FellOffRockArray); */
			  FindRockStart += FindCellError;
			  if (FellOffRockArray == 'y') {
				  /*printf("NOODLE  !!!!!FoundRockCells: %d GetO = %c \n",
				   * FindRockStart,FellOffRockArray); */
				   /*PauseRun(1,1,-1); */
			  }

			  if (FindRockStart > Y_MAX / 2 + 1) {
				  printf("Stopped Finding Rock - done %d %d", FindRockStart,
					  Y_MAX / 2 - 5);
				  if (SAVE_FILE) SaveSandToFile();
				  return 1;
			  }
		  }

		  if (debug0) {
			  printf("going to pause after finding beach, rock\n");
			  PauseRun(58, 269, -1);
		  }

		  if (INITIAL_CONDITION_TYPE != 3) {
			  RockCalculations();
		  }
	  }
      /*LMV*/ ShadowSweep();

      DetermineAngles();

      DetermineSedTransport();

      DoSink(); /* clear sand from sink cells */

      FlowInCell();
      /*LMV*/ FixFlow();
      /*LMV*/ TransportSedimentSweep();

      FixBeach();

      /* Age Empty Cells */
      if ((current_time_step % AGE_UPDATE == 0) && SAVE_AGE) AgeCells();

      /* Count Mass */
      MassCurrent = MassCount();

      /* SAVE FILE ? */
      if ((current_time_step % SAVE_SPACING == 0 &&
           current_time_step > START_SAVING_AT)) {
        if (SAVE_LINE) SaveLineToFile();
        if (SAVE_FILE) SaveSandToFile();
      }

      current_time_step++;
      current_time += time_step;
    }
  }
  return 0;
}

int cem_finalize(void) {
  free2d((void**)topography);
  free2d((void**)shelf_depth);
  free2d((void**)wave_h_sig);
  free2d((void**)wave_dir);
  free2d((void**)type_of_rock);
  free2d((void**)AllBeach);
  free2d((void**)AllRock);
  free2d((void**)PercentFullSand);
  free2d((void**)PercentFullRock);
  free2d((void**)cell_depth);

  printf("Run Complete.  Output file: %s \n", savefilename);

  return 0;
}

/* -----FUNCTIONS ARE BELOW----- */

float FindWaveAngle(void)
/* calculates wave angle for given time step */
{
  float Angle;

  /* Data Input Method */

  float RandBin;   /* Random number to pick wave angle bin */
  float RandAngle; /* Random number to pick wave angle within the bin */
  int flag = 1;
  int i = 0;

  /* Method using input binned wave distribution - */
  /* variables WaveProb[], WaveMax[], previously input from file using
   * ReadWaveIn()       */

  if (WAVE_IN == 1) {
    RandBin = RandZeroToOne();
    RandAngle = RandZeroToOne();

    /*printf("Time = %d RandBin = %f RandAng = %f\n",current_time_step, RandBin,
     * RandAngle); */

    while (flag) {
      i++;

      if (RandBin < WaveProb[i]) {
        flag = 0;
      }
    }

    Angle = -(RandAngle * (WaveMax[i] - WaveMax[i - 1]) + WaveMax[i - 1]) * PI /
            180;

    /*printf("i = %d WaveMAx[i] = %f WaveMax[i-1] = %f WaveProb[i] = %f Angle=
       %f\n",
       i,WaveMax[i],WaveMax[i-1],WaveProb[i], Angle*180/PI); */
  }

  if (WAVE_IN == 2) {
    RealWaveIn(); /*Read real data (angle/period/height) */

    Angle = ((float)(WaveAngleIn)) +
            coastrotation; /*align waveclimate to coastline in CEM */
    if (Angle >= 180.0)
      Angle -= 360.0; /* conversion to CEM angles (-180 -- 0 -- +180) */
    Angle *= -1; /*-1 is there as wave angle is from 90 in the west to -90 in
                    the east*/

    Period = WavePeriodIn * waveperiodchange;
    OffShoreWvHt = WaveHeightIn * waveheightchange;

    if ((Angle > 90) || (Angle < -90)) {
      Angle = 0; /*if waves were to come from behind coast */
      OffShoreWvHt =
          1e-10; /*if waves were to come from behind coast set to inf small */
    }

    Angle /= RAD_TO_DEG;
  }

  if (WAVE_IN == 0) {
    /* Method using wave probability step function */
    /*      variable Asym will determine fractional distribution of waves coming
     * from the           */
    /*      positive direction (positive direction coming from left)  -i.e.
     * fractional wave asymmetry */
    /*  redone by RCL. note Highness actually means lowness, for some reason */

    Angle = RandZeroToOne() * PI / 4.0;

    if (RandZeroToOne() >= Highness) Angle += PI / 4.0; /* make high angle */

    if (RandZeroToOne() >= Asym) Angle = -1.0 * Angle;

    /*printf("WaveAngle: %f degrees\n",Angle*RAD_TO_DEG); */
  }

  return Angle;
}


void FindBeachCells(int YStart)
/* Determines locations of beach cells moving from left to right direction
   */
/* This function will affect and determine the global arrays:  X[] and Y[]
   */
/* This function calls FindNextCell */
/* This will define TotalBeachCells for this time step */
{
  int y, z, xstart; /* local iterators */

  /* Starting at left end, find the x - value for first cell that is 'allbeach'
   */

  xstart = X_MAX - 1;
  y = YStart;

  while (AllBeach[xstart][y] == 'n') {
    xstart -= 1;
  }

  xstart += 1; /* Step back to where partially full beach */

  X[0] = xstart;
  Y[0] = YStart;

  if (debug1) printf("FirstX: %3d  FirstY: %3d  z: 0 \n", X[0], Y[0]);

  z = 0;

  while ((Y[z] < 2 * Y_MAX - 1) && (z < MaxBeachLength - 1)) {
    NextX = -2;
    NextY = -2;
    BehindX = -2;
    BehindY = -2;

    FindNextCell(X[z], Y[z], z);
    z++;
    X[z] = NextX;
    Y[z] = NextY;
    XBehind[z] = BehindX;
    YBehind[z] = BehindY;

    if (debug1) printf("NextX: %3d  NextY: %3d  z: %d \n", NextX, NextY, z);
    if (debug1)
      printf("z: %d, BehindX: %d; BehindY: %d; \n\n", z, BehindX, BehindY);

    if (PercentFullSand[X[z]][Y[z]] == 0) {
      /*printf("\nFINDBEACH: PercentFullSand Zero x: %d y: %d\n",X[z],Y[z]); */
      /*PauseRun(X[z],Y[z],z); */
    }

    /* If return to start point or go off left side of array, going the wrong
     * direction     */
    /* Jump off and start again closer to middle */

    if ((NextY < 1) || ((NextY == Y[0]) && (NextX == X[0])) ||
        (z > MaxBeachLength - 2)) {
      /*printf("!!!!!!!Fell Off!!!!!!!!!!!!!!!!! x = %d !!!!!!!!!!!!!", NextX);
       */
      FellOffArray = 'y';
      ZeroVars();
      return;
    }

    if (z > MaxBeachLength - 3) {
      printf("????????????  went to end of MaxBeach!! ????");
    }
  }

  TotalBeachCells = z;
  FellOffArray = 'n';

  if (debug1)
    printf("Total Beach: %d  Current Time Step: %d \n \n", TotalBeachCells,
           current_time_step);
}

void FindRockCells(int YStart)
/* LMV */
/* Determines locations of rock cells moving from left to right direction
   */
/* This function will affect and determine the global arrays:  XRock[] and
   YRock[]	*/
/* This function calls FindNextRockCell
   */
/* This will define TotalRockCells for this time step
   */
{
	if (INIT_ROCK == 0)
	{
		return;
	}
  int y, z, xstart; /* local iterators */

  /* Starting at left end, find the x - value for first cell that is 'allbeach'
   */

  xstart = X_MAX - 1;
  y = YStart;

  while (AllRock[xstart][y] == 'n') {
    xstart -= 1;
  }

  xstart += 1; /* Step back to where partially full beach */

  XRock[0] = xstart;
  YRock[0] = YStart;

  if (debug1a)
    printf("FirstRockX: %3d  FirstRockY: %3d  z: 0 \n", XRock[0], YRock[0]);

  z = 0;
  XRock[z - 1] = XRock[z];
  YRock[z - 1] = YRock[z] - 1;

  while ((YRock[z] < 2 * Y_MAX - 1) && (z < MaxBeachLength - 1)) {
    z++;
    NextRockX = -2;
    NextRockY = -2;
    BehindRockX = -2;
    BehindRockY = -2;

    FindNextRockCell(XRock[z - 1], YRock[z - 1], z - 1);
    XRock[z] = NextRockX;
    YRock[z] = NextRockY;
    XRockBehind[z] = BehindRockX;
    YRockBehind[z] = BehindRockY;

    if (debug1a)
      printf("NextRockX: %3d  NextRockY: %3d  z: %d \n", NextRockX, NextRockY,
             z);
    if (debug1a)
      printf("z: %d, BehindRockX: %d; BehindRockY: %d \n\n", z, BehindRockX,
             BehindRockY);

    /*RCL takes this out  if (PercentFullRock[XRock[z]][YRock[z]] == 0)
       {
       printf("\nFINDROCK: PercentFullRock Zero x: %d y:
       %d\n",XRock[z],YRock[z]);

       }
     */
    /* If return to start point or go off left side of array, going the wrong
     * direction     */
    /* Jump off and start again closer to middle */

    if ((NextRockY < 1) ||
        ((NextRockY == YRock[0]) && (NextRockX == XRock[0])) ||
        (z > MaxBeachLength - 2)) {
      printf("!!!!!!!Fell Off Rock!!!!!!!!!!!!! x = %d !!!!!!!!!!!!!",
             NextRockX);
      FellOffRockArray = 'y';
      ZeroVars();
      return;
    }

    if (z > MaxBeachLength - 3) {
      printf("????????????  went to end of Rock MaxBeach!! ????");
      printf("(maybe the beach made a loop)\n");
    }
  }

  TotalRockCells = z;
  FellOffRockArray = 'n';

  if (debug1a)
    printf("Total Rock: %d  Current Time Step: %d \n \n", TotalRockCells,
           current_time_step);
}

void FindNextCell(int x, int y, int z)

/* Function to find next cell that is beach moving in the general positive X
   direction */
/* changes global variables NextX and NextY, coordinates for the next beach cell
   */
/* This function will use but not affect the global arrays:  AllBeach [][], X[],
   and Y[] */
/* New thinking...5/04 LMV */
{
	if (x == X_MAX - 1 || y == 2 * Y_MAX - 1)
	{
		return;
	}
  if (z == 0 ||(X[z - 1] == X[z]) && (Y[z - 1] == Y[z] - 1))
  /* came from left */
  {
    if ((PercentFullSand[x + 1][y] + PercentFullRock[x + 1][y]) > 0.0) {
      if (debug1b) printf("1\n");
      NextX = x + 1;
      NextY = y;
      BehindX = NextX - 1;
      BehindY = NextY;
      return;
    } else if ((PercentFullSand[x + 1][y + 1] + PercentFullRock[x + 1][y + 1]) >
               0.0) {
      if (debug1b) printf("2\n");
      NextX = x + 1;
      NextY = y + 1;
      BehindX = NextX - 1;
      BehindY = NextY;
      return;
    } else if ((PercentFullSand[x][y + 1] + PercentFullRock[x][y + 1]) > 0.0) {
      if (debug1b) printf("3\n");
      NextX = x;
      NextY = y + 1;
      BehindX = NextX - 1;
      BehindY = NextY;
      return;
    } else if ((PercentFullSand[x - 1][y + 1] + PercentFullRock[x - 1][y + 1]) >
               0.0) {
      if (debug1b) printf("4\n");
      NextX = x - 1;
      NextY = y + 1;
      BehindX = NextX - 1;
      BehindY = NextY;
      return;
    } else if ((PercentFullSand[x - 1][y] + PercentFullRock[x - 1][y]) > 0.0) {
      if (debug1b) printf("5\n");
      NextX = x - 1;
      NextY = y;
      BehindX = NextX - 1;
      BehindY = NextY;
      return;
    } else if ((PercentFullSand[x - 1][y - 1] + PercentFullRock[x - 1][y - 1]) >
               0.0) {
      if (debug1b) printf("6\n");
      NextX = x - 1;
      NextY = y - 1;
      BehindX = NextX - 1;
      BehindY = NextY;
      return;
    } else
      printf("Should've found next cell (Left): %d, %d \n", x, y);
  }

  else if ((X[z - 1] == X[z] + 1) && (Y[z - 1] == Y[z] - 1))
  /* came from upper left */
  {
    if ((PercentFullSand[x + 1][y + 1] + PercentFullRock[x + 1][y + 1]) > 0.0) {
      if (debug1b) printf("7\n");
      NextX = x + 1;
      NextY = y + 1;
      BehindX = NextX - 1;
      BehindY = NextY;
      return;
    } else if ((PercentFullSand[x][y + 1] + PercentFullRock[x][y + 1]) > 0.0) {
      if (debug1b) printf("8\n");
      NextX = x;
      NextY = y + 1;
      BehindX = NextX - 1;
      BehindY = NextY;
      return;
    } else if ((PercentFullSand[x - 1][y + 1] + PercentFullRock[x - 1][y + 1]) >
               0.0) {
      if (debug1b) printf("9\n");
      NextX = x - 1;
      NextY = y + 1;
      BehindX = NextX - 1;
      BehindY = NextY;
      return;
    } else if ((PercentFullSand[x - 1][y] + PercentFullRock[x - 1][y]) > 0.0) {
      if (debug1b) printf("10\n");
      NextX = x - 1;
      NextY = y;
      BehindX = NextX;
      BehindY = NextY - 1;
      return;
    } else if ((PercentFullSand[x - 1][y - 1] + PercentFullRock[x - 1][y - 1]) >
               0.0) {
      if (debug1b) printf("11\n");
      NextX = x - 1;
      NextY = y - 1;
      BehindX = NextX;
      BehindY = NextY - 1;
      return;
    } else if ((PercentFullSand[x][y - 1] + PercentFullRock[x][y - 1]) > 0.0) {
      if (debug1b) printf("12\n");
      NextX = x;
      NextY = y - 1;
      BehindX = NextX;
      BehindY = NextY - 1;
      return;
    } else {
      printf("Should've found next cell (Upper Left): %d, %d \n", x, y);
      PauseRun(x, y, -1);
    }
  }

  else if ((X[z - 1] == X[z] + 1) && (Y[z - 1] == Y[z]))
  /* came from above */
  {
    if ((PercentFullSand[x][y + 1] + PercentFullRock[x][y + 1]) > 0.0) {
      if (debug1b) printf("13\n");
      NextX = x;
      NextY = y + 1;
      BehindX = NextX;
      BehindY = NextY - 1;
      return;
    } else if ((PercentFullSand[x - 1][y + 1] + PercentFullRock[x - 1][y + 1]) >
               0.0) {
      if (debug1b) printf("14\n");
      NextX = x - 1;
      NextY = y + 1;
      BehindX = NextX;
      BehindY = NextY - 1;
      return;
    } else if ((PercentFullSand[x - 1][y] + PercentFullRock[x - 1][y]) > 0.0) {
      if (debug1b) printf("15\n");
      NextX = x - 1;
      NextY = y;
      BehindX = NextX;
      BehindY = NextY - 1;
      return;
    } else if ((PercentFullSand[x - 1][y - 1] + PercentFullRock[x - 1][y - 1]) >
               0.0) {
      if (debug1b) printf("16\n");
      NextX = x - 1;
      NextY = y - 1;
      BehindX = NextX;
      BehindY = NextY - 1;
      return;
    } else if ((PercentFullSand[x][y - 1] + PercentFullRock[x][y - 1]) > 0.0) {
      if (debug1b) printf("17\n");
      NextX = x;
      NextY = y - 1;
      BehindX = NextX;
      BehindY = NextY - 1;
      return;
    } else if ((PercentFullSand[x + 1][y - 1] + PercentFullRock[x + 1][y - 1]) >
               0.0) {
      if (debug1b) printf("18\n");
      NextX = x + 1;
      NextY = y - 1;
      BehindX = NextX;
      BehindY = NextY - 1;
      return;
    } else if ((PercentFullSand[x + 1][y - 1] + PercentFullRock[x + 1][y - 1]) <
               1.0) { /* RCL adds ad hoc case */
      if (debug1b) printf("18.5\n");
      NextX = x + 1;
      NextY = y - 1;
      BehindX = NextX;
      BehindY = NextY - 1;
      return;
    } else {
      printf("Should've found next cell (Above): %d, %d \n", x, y);
      PauseRun(x, y, -1);
    }
  }

  else if ((X[z - 1] == X[z] + 1) && (Y[z - 1] == Y[z] + 1))
  /* came from upper right */
  {
    if ((PercentFullSand[x - 1][y + 1] + PercentFullRock[x - 1][y + 1]) > 0.0) {
      if (debug1b) printf("19\n");
      NextX = x - 1;
      NextY = y + 1;
      BehindX = NextX;
      BehindY = NextY - 1;
      return;
    } else if ((PercentFullSand[x - 1][y] + PercentFullRock[x - 1][y]) > 0.0) {
      if (debug1b) printf("20\n");
      NextX = x - 1;
      NextY = y;
      BehindX = NextX;
      BehindY = NextY - 1;
      return;
    } else if ((PercentFullSand[x - 1][y - 1] + PercentFullRock[x - 1][y - 1]) >
               0.0) {
      if (debug1b) printf("21\n");
      NextX = x - 1;
      NextY = y - 1;
      BehindX = NextX + 1;
      BehindY = NextY;
      return;
    } else if ((PercentFullSand[x][y - 1] + PercentFullRock[x][y - 1]) > 0.0) {
      if (debug1b) printf("22\n");
      NextX = x;
      NextY = y - 1;
      BehindX = NextX + 1;
      BehindY = NextY;
      return;
    } else if ((PercentFullSand[x + 1][y - 1] + PercentFullRock[x + 1][y - 1]) >
               0.0) {
      if (debug1b) printf("23\n");
      NextX = x + 1;
      NextY = y - 1;
      BehindX = NextX + 1;
      BehindY = NextY;
      return;
    } else if ((PercentFullSand[x + 1][y] + PercentFullRock[x + 1][y]) > 0.0) {
      if (debug1b) printf("24\n");
      NextX = x + 1;
      NextY = y;
      BehindX = NextX + 1;
      BehindY = NextY;
      return;
    } else {
      printf("Should've found next cell (Upper Right): %d, %d \n", x, y);
      PauseRun(x, y, -1);
    }
  }

  else if ((X[z - 1] == X[z]) && (Y[z - 1] == Y[z] + 1))
  /* came from right */
  {
    if ((PercentFullSand[x - 1][y] + PercentFullRock[x - 1][y]) > 0.0) {
      if (debug1b) printf("25\n");
      NextX = x - 1;
      NextY = y;
      BehindX = NextX + 1;
      BehindY = NextY;
      return;
    } else if ((PercentFullSand[x - 1][y - 1] + PercentFullRock[x - 1][y - 1]) >
               0.0) {
      if (debug1b) printf("26\n");
      NextX = x - 1;
      NextY = y - 1;
      BehindX = NextX + 1;
      BehindY = NextY;
      return;
    } else if ((PercentFullSand[x][y - 1] + PercentFullRock[x][y - 1]) > 0.0) {
      if (debug1b) printf("27\n");
      NextX = x;
      NextY = y - 1;
      BehindX = NextX + 1;
      BehindY = NextY;
      return;
    } else if ((PercentFullSand[x + 1][y - 1] + PercentFullRock[x + 1][y - 1]) >
               0.0) {
      if (debug1b) printf("28\n");
      NextX = x + 1;
      NextY = y - 1;
      BehindX = NextX + 1;
      BehindY = NextY;
      return;
    } else if ((PercentFullSand[x + 1][y] + PercentFullRock[x + 1][y]) > 0.0) {
      if (debug1b) printf("29\n");
      NextX = x + 1;
      NextY = y;
      BehindX = NextX + 1;
      BehindY = NextY;
      return;
    } else if ((PercentFullSand[x + 1][y + 1] + PercentFullRock[x + 1][y + 1]) >
               0.0) {
      if (debug1b) printf("30\n");
      NextX = x + 1;
      NextY = y + 1;
      BehindX = NextX + 1;
      BehindY = NextY;
      return;
    } else {
      printf("Should've found next cell (Right): %d, %d \n", x, y);
      PauseRun(x, y, -1);
    }
  }

  else if ((X[z - 1] == X[z] - 1) && (Y[z - 1] == Y[z] + 1))
  /* came from lower right */
  {
    if ((PercentFullSand[x - 1][y - 1] + PercentFullRock[x - 1][y - 1]) > 0.0) {
      if (debug1b) printf("31\n");
      NextX = x - 1;
      NextY = y - 1;
      BehindX = NextX + 1;
      BehindY = NextY;
      return;
    } else if ((PercentFullSand[x][y - 1] + PercentFullRock[x][y - 1]) > 0.0) {
      if (debug1b) printf("32\n");
      NextX = x;
      NextY = y - 1;
      BehindX = NextX + 1;
      BehindY = NextY;
      return;
    } else if ((PercentFullSand[x + 1][y - 1] + PercentFullRock[x + 1][y - 1]) >
               0.0) {
      if (debug1b) printf("33\n");
      NextX = x + 1;
      NextY = y - 1;
      BehindX = NextX + 1;
      BehindY = NextY;
      return;
    } else if ((PercentFullSand[x + 1][y] + PercentFullRock[x + 1][y]) > 0.0) {
      if (debug1b) printf("34\n");
      NextX = x + 1;
      NextY = y;
      BehindX = NextX;
      BehindY = NextY + 1;
      return;
    } else if ((PercentFullSand[x + 1][y + 1] + PercentFullRock[x + 1][y + 1]) >
               0.0) {
      if (debug1b) printf("35\n");
      NextX = x + 1;
      NextY = y + 1;
      BehindX = NextX;
      BehindY = NextY + 1;
      return;
    } else if ((PercentFullSand[x][y + 1] + PercentFullRock[x][y + 1]) > 0.0) {
      if (debug1b) printf("36\n");
      NextX = x;
      NextY = y + 1;
      BehindX = NextX;
      BehindY = NextY + 1;
      return;
    } else if (PercentFullSand[x][y + 2] + PercentFullRock[x][y + 2] > 0.0 ||
               PercentFullSand[x - 1][y + 2] + PercentFullRock[x - 1][y + 2] >
                   0.0 ||
               PercentFullSand[x + 1][y + 2] + PercentFullRock[x + 1][y + 2] >
                   0.0) {
      /* RCL adds new ad hoc case */
      if (debug1b) printf("36.5\n");
      NextX = x;
      NextY = y + 1;
      BehindX = NextX;
      BehindY = NextY + 1;
      return;
    }

    else {
      printf("Should've found next cell (Lower Right): %d, %d \n", x, y);
      PauseRun(x, y, -1);
    }
  }

  else if ((X[z - 1] == X[z] - 1) && (Y[z - 1] == Y[z]))
  /* came from below */
  {
    if ((PercentFullSand[x][y - 1] + PercentFullRock[x][y - 1]) > 0.0) {
      if (debug1b) printf("37\n");
      NextX = x;
      NextY = y - 1;
      BehindX = NextX;
      BehindY = NextY + 1;
      return;
    } else if ((PercentFullSand[x + 1][y - 1] + PercentFullRock[x + 1][y - 1]) >
               0.0) {
      if (debug1b) printf("38\n");
      NextX = x + 1;
      NextY = y - 1;
      BehindX = NextX;
      BehindY = NextY + 1;
      return;
    } else if ((PercentFullSand[x + 1][y] + PercentFullRock[x + 1][y]) > 0.0) {
      if (debug1b) printf("39\n");
      NextX = x + 1;
      NextY = y;
      BehindX = NextX;
      BehindY = NextY + 1;
      return;
    } else if ((PercentFullSand[x + 1][y + 1] + PercentFullRock[x + 1][y + 1]) >
               0.0) { /* RCL fixed typo */
      if (debug1b) printf("40\n");
      NextX = x + 1;
      NextY = y + 1;
      BehindX = NextX;
      BehindY = NextY + 1;
      return;
    } else if ((PercentFullSand[x][y + 1] + PercentFullRock[x][y + 1]) >
               0.0) { /* RCL changed behind */
      if (debug1b) printf("41\n");
      NextX = x;
      NextY = y + 1;
      BehindX = NextX - 1;
      BehindY = NextY;
      return;
    } else if ((PercentFullSand[x - 1][y + 1]) + PercentFullRock[x - 1][y + 1] >
               0.0)
    /* RCL messed with all the 42.X cases */
    {
      if (debug1b) printf("42\n");
      NextX = x - 1;
      NextY = y + 1;
      BehindX = NextX - 1;
      BehindY = NextY;
      return;
    } else if ((PercentFullSand[x - 1][y - 1]) + PercentFullRock[x - 1][y - 1] >
               0.0) {
      if (debug1b) printf("42.1\n");
      NextX = x - 1;
      NextY = y - 1;
      BehindX = NextX + 1;
      BehindY = NextY - 1;
      return;
    } else if ((PercentFullSand[x][y + 1]) + PercentFullRock[x][y + 1] < 1.0)
    /* RCL added this case */
    {
      if (debug1b) printf("42.75\n");
      NextX = x;
      NextY = y + 1;
      BehindX = NextX - 1;
      BehindY = NextY;
      return;
    } else if ((PercentFullSand[x - 1][y + 1]) + PercentFullRock[x - 1][y + 1] <
               1.0)
    /* RCL added this case */
    {
      if (debug1b) printf("42.8\n");
      NextX = x - 1;
      NextY = y + 1;
      BehindX = NextX - 1;
      BehindY = NextY;
      return;
    } else {
      printf("Should've found next cell (Below): %d, %d \n", x, y);
      PauseRun(x, y, -1);
    }
  } else if ((X[z - 1] == X[z] - 1) && (Y[z - 1] == Y[z] - 1))
  /* came from lower left */
  {
    if ((PercentFullSand[x + 1][y - 1] + PercentFullRock[x + 1][y - 1]) > 0.0) {
      if (debug1b) printf("43\n");
      NextX = x + 1;
      NextY = y - 1;
      BehindX = NextX;
      BehindY = NextY + 1;
      return;
    } else if ((PercentFullSand[x + 1][y] + PercentFullRock[x + 1][y]) > 0.0) {
      if (debug1b) printf("44\n");
      NextX = x + 1;
      NextY = y;
      BehindX = NextX;
      BehindY = NextY + 1;
      return;
    } else if ((PercentFullSand[x + 1][y + 1] + PercentFullRock[x + 1][y + 1]) >
               0.0) {
      if (debug1b) printf("45\n");
      NextX = x + 1;
      NextY = y + 1;
      BehindX = NextX - 1;
      BehindY = NextY;
      return;
    } else if ((PercentFullSand[x][y + 1] + PercentFullRock[x][y + 1]) > 0.0) {
      if (debug1b) printf("46\n");
      NextX = x;
      NextY = y + 1;
      BehindX = NextX - 1;
      BehindY = NextY;
      return;
    } else if ((PercentFullSand[x - 1][y + 1] + PercentFullRock[x - 1][y + 1]) >
               0.0) {
      if (debug1b) printf("47\n");
      NextX = x - 1;
      NextY = y + 1;
      BehindX = NextX - 1;
      BehindY = NextY;
      return;
    } else if ((PercentFullSand[x - 1][y] + PercentFullRock[x - 1][y]) > 0.0) {
      if (debug1b) printf("48\n");
      NextX = x - 1;
      NextY = y;
      BehindX = NextX - 1;
      BehindY = NextY;
      return;
    } else {
      printf("Should've found next cell (Lower Left): %d, %d \n", x, y);
      PauseRun(x, y, -1);
    }
  }

  printf("Should've found next cell : %d, %d \n", x, y);
}

void FindNextRockCell(int x, int y, int z)
/* LMV */
/* Function to find next cell that is rock moving in the general positive X
   direction           */
/* changes global variables NextRockX and NextRockY, coordinates for the next
   beach cell        */
/* This function will use but not affect the global arrays:  AllRock [][],
   XRock[], and YRock[] */
{
  if ((XRock[z - 1] == XRock[z]) && (YRock[z - 1] == YRock[z] - 1))
  /* came from left */
  {
    if (PercentFullRock[x + 1][y] > 0.0) {
      if (debug1c) printf("1 r\n");
      NextRockX = x + 1;
      NextRockY = y;
      BehindRockX = NextRockX - 1;
      BehindRockY = NextRockY;
      return;
    } else if (PercentFullRock[x + 1][y + 1] > 0.0) {
      if (debug1c) printf("2 r\n");
      NextRockX = x + 1;
      NextRockY = y + 1;
      BehindRockX = NextRockX - 1;
      BehindRockY = NextRockY;
      return;
    } else if (PercentFullRock[x][y + 1] > 0.0) {
      if (debug1c) printf("3 r\n");
      NextRockX = x;
      NextRockY = y + 1;
      BehindRockX = NextRockX - 1;
      BehindRockY = NextRockY;
      return;
    } else if (PercentFullRock[x - 1][y + 1] > 0.0) {
      if (debug1c) printf("4 r\n");
      NextRockX = x - 1;
      NextRockY = y + 1;
      BehindRockX = NextRockX - 1;
      BehindRockY = NextRockY;
      return;
    } else if (PercentFullRock[x - 1][y] > 0.0) {
      if (debug1c) printf("5 r\n");
      NextRockX = x - 1;
      NextRockY = y;
      BehindRockX = NextRockX - 1;
      BehindRockY = NextRockY;
      return;
    } else if (PercentFullRock[x - 1][y - 1] > 0.0) {
      if (debug1c) printf("6 r\n");
      NextRockX = x - 1;
      NextRockY = y - 1;
      BehindRockX = NextRockX - 1;
      BehindRockY = NextRockY;
      return;
    } else
      printf("Should've found next Rock cell (Left): %d, %d \n", x, y);
  }

  else if ((XRock[z - 1] == XRock[z] + 1) && (YRock[z - 1] == YRock[z] - 1))
  /* came from upper left */
  {
    if (PercentFullRock[x + 1][y + 1] > 0.0) {
      if (debug1c) printf("7 r\n");
      NextRockX = x + 1;
      NextRockY = y + 1;
      BehindRockX = NextRockX - 1;
      BehindRockY = NextRockY;
      return;
    } else if (PercentFullRock[x][y + 1] > 0.0) {
      if (debug1c) printf("8 r\n");
      NextRockX = x;
      NextRockY = y + 1;
      BehindRockX = NextRockX - 1;
      BehindRockY = NextRockY;
      return;
    } else if (PercentFullRock[x - 1][y + 1] > 0.0) {
      if (debug1c) printf("9 r\n");
      NextRockX = x - 1;
      NextRockY = y + 1;
      BehindRockX = NextRockX - 1;
      BehindRockY = NextRockY;
      return;
    } else if (PercentFullRock[x - 1][y] > 0.0) {
      if (debug1c) printf("10 r\n");
      NextRockX = x - 1;
      NextRockY = y;
      BehindRockX = NextRockX;
      BehindRockY = NextRockY - 1;
      return;
    } else if (PercentFullRock[x - 1][y - 1] > 0.0) {
      if (debug1c) printf("11 r\n");
      NextRockX = x - 1;
      NextRockY = y - 1;
      BehindRockX = NextRockX;
      BehindRockY = NextRockY - 1;
      return;
    } else if (PercentFullRock[x][y - 1] > 0.0) {
      if (debug1c) printf("12 r\n");
      NextRockX = x;
      NextRockY = y - 1;
      BehindRockX = NextRockX;
      BehindRockY = NextRockY - 1;
      return;
    } else
      printf("Should've found next Rock cell (Upper Left): %d, %d \n", x, y);
  }

  else if ((XRock[z - 1] == XRock[z] + 1) && (YRock[z - 1] == YRock[z]))
  /* came from above */
  {
    if (PercentFullRock[x][y + 1] > 0.0) {
      if (debug1c) printf("13 r\n");
      NextRockX = x;
      NextRockY = y + 1;
      BehindRockX = NextRockX;
      BehindRockY = NextRockY - 1;
      return;
    } else if (PercentFullRock[x - 1][y + 1] > 0.0) {
      if (debug1c) printf("14 r\n");
      NextRockX = x - 1;
      NextRockY = y + 1;
      BehindRockX = NextRockX;
      BehindRockY = NextRockY - 1;
      return;
    } else if (PercentFullRock[x - 1][y] > 0.0) {
      if (debug1c) printf("15 r\n");
      NextRockX = x - 1;
      NextRockY = y;
      BehindRockX = NextRockX;
      BehindRockY = NextRockY - 1;
      return;
    } else if (PercentFullRock[x - 1][y - 1] > 0.0) {
      if (debug1c) printf("16 r\n");
      NextRockX = x - 1;
      NextRockY = y - 1;
      BehindRockX = NextRockX;
      BehindRockY = NextRockY - 1;
      return;
    } else if (PercentFullRock[x][y - 1] > 0.0) {
      if (debug1c) printf("17 r\n");
      NextRockX = x;
      NextRockY = y - 1;
      BehindRockX = NextRockX;
      BehindRockY = NextRockY - 1;
      return;
    } else if (PercentFullRock[x + 1][y - 1] > 0.0) {
      if (debug1c) printf("18 r\n");
      NextRockX = x + 1;
      NextRockY = y - 1;
      BehindRockX = NextRockX;
      BehindRockY = NextRockY - 1;
      return;
    } else
      printf("Should've found next Rock cell (Above): %d, %d \n", x, y);
  }

  else if ((XRock[z - 1] == XRock[z] + 1) && (YRock[z - 1] == YRock[z] + 1))
  /* came from upper right */
  {
    if (PercentFullRock[x - 1][y + 1] > 0.0) {
      if (debug1c) printf("19 r\n");
      NextRockX = x - 1;
      NextRockY = y + 1;
      BehindRockX = NextRockX;
      BehindRockY = NextRockY - 1;
      return;
    } else if (PercentFullRock[x - 1][y] > 0.0) {
      if (debug1c) printf("20 r\n");
      NextRockX = x - 1;
      NextRockY = y;
      BehindRockX = NextRockX;
      BehindRockY = NextRockY - 1;
      return;
    } else if (PercentFullRock[x - 1][y - 1] > 0.0) {
      if (debug1c) printf("21 r\n");
      NextRockX = x - 1;
      NextRockY = y - 1;
      BehindRockX = NextRockX + 1;
      BehindRockY = NextRockY;
      return;
    } else if (PercentFullRock[x][y - 1] > 0.0) {
      if (debug1c) printf("22 r\n");
      NextRockX = x;
      NextRockY = y - 1;
      BehindRockX = NextRockX + 1;
      BehindRockY = NextRockY;
      return;
    } else if (PercentFullRock[x + 1][y - 1] > 0.0) {
      if (debug1c) printf("23 r\n");
      NextRockX = x + 1;
      NextRockY = y - 1;
      BehindRockX = NextRockX + 1;
      BehindRockY = NextRockY;
      return;
    } else if (PercentFullRock[x + 1][y] > 0.0) {
      if (debug1c) printf("24 r\n");
      NextRockX = x + 1;
      NextRockY = y;
      BehindRockX = NextRockX + 1;
      BehindRockY = NextRockY;
      return;
    } else
      printf("Should've found next Rock cell (Upper Right): %d, %d \n", x, y);
  }

  else if ((XRock[z - 1] == XRock[z]) && (YRock[z - 1] == YRock[z] + 1))
  /* came from right */
  {
    if (PercentFullRock[x - 1][y] > 0.0) {
      if (debug1c) printf("25 r\n");
      NextRockX = x - 1;
      NextRockY = y;
      BehindRockX = NextRockX + 1;
      BehindRockY = NextRockY;
      return;
    } else if (PercentFullRock[x - 1][y - 1] > 0.0) {
      if (debug1c) printf("26 r\n");
      NextRockX = x - 1;
      NextRockY = y - 1;
      BehindRockX = NextRockX + 1;
      BehindRockY = NextRockY;
      return;
    } else if (PercentFullRock[x][y - 1] > 0.0) {
      if (debug1c) printf("27 r\n");
      NextRockX = x;
      NextRockY = y - 1;
      BehindRockX = NextRockX + 1;
      BehindRockY = NextRockY;
      return;
    } else if (PercentFullRock[x + 1][y - 1] > 0.0) {
      if (debug1c) printf("28 r\n");
      NextRockX = x + 1;
      NextRockY = y - 1;
      BehindRockX = NextRockX + 1;
      BehindRockY = NextRockY;
      return;
    } else if (PercentFullRock[x + 1][y] > 0.0) {
      if (debug1c) printf("29 r\n");
      NextRockX = x + 1;
      NextRockY = y;
      BehindRockX = NextRockX + 1;
      BehindRockY = NextRockY;
      return;
    } else if (PercentFullRock[x + 1][y + 1] > 0.0) {
      if (debug1c) printf("30 r\n");
      NextRockX = x + 1;
      NextRockY = y + 1;
      BehindRockX = NextRockX + 1;
      BehindRockY = NextRockY;
      return;
    } else
      printf("Should've found next Rock cell (Right): %d, %d \n", x, y);
  }

  else if ((XRock[z - 1] == XRock[z] - 1) && (YRock[z - 1] == YRock[z] + 1))
  /* came from lower right */
  {
    if (PercentFullRock[x - 1][y - 1] > 0.0) {
      if (debug1c) printf("31 r\n");
      NextRockX = x - 1;
      NextRockY = y - 1;
      BehindRockX = NextRockX + 1;
      BehindRockY = NextRockY;
      return;
    } else if (PercentFullRock[x][y - 1] > 0.0) {
      if (debug1c) printf("32 r\n");
      NextRockX = x;
      NextRockY = y - 1;
      BehindRockX = NextRockX + 1;
      BehindRockY = NextRockY;
      return;
    } else if (PercentFullRock[x + 1][y - 1] > 0.0) {
      if (debug1c) printf("33 r\n");
      NextRockX = x + 1;
      NextRockY = y - 1;
      BehindRockX = NextRockX + 1;
      BehindRockY = NextRockY;
      return;
    } else if (PercentFullRock[x + 1][y] > 0.0) {
      if (debug1c) printf("34 r\n");
      NextRockX = x + 1;
      NextRockY = y;
      BehindRockX = NextRockX;
      BehindRockY = NextRockY + 1;
      return;
    } else if (PercentFullRock[x + 1][y + 1] > 0.0) {
      if (debug1c) printf("35 r\n");
      NextRockX = x + 1;
      NextRockY = y + 1;
      BehindRockX = NextRockX;
      BehindRockY = NextRockY + 1;
      return;
    } else if (PercentFullRock[x][y + 1] > 0.0) {
      if (debug1c) printf("36 r\n");
      NextRockX = x;
      NextRockY = y + 1;
      BehindRockX = NextRockX;
      BehindRockY = NextRockY + 1;
      return;
    } else
      printf("Should've found next Rock cell (Lower Right): %d, %d \n", x, y);
  }

  else if ((XRock[z - 1] == XRock[z] - 1) && (YRock[z - 1] == YRock[z]))
  /* came from below */
  {
    if (PercentFullRock[x][y - 1] > 0.0) {
      if (debug1c) printf("37 r\n");
      NextRockX = x;
      NextRockY = y - 1;
      BehindRockX = NextRockX;
      BehindRockY = NextRockY + 1;
      return;
    } else if (PercentFullRock[x + 1][y - 1] > 0.0) {
      if (debug1c) printf("38 r\n");
      NextRockX = x + 1;
      NextRockY = y - 1;
      BehindRockX = NextRockX;
      BehindRockY = NextRockY + 1;
      return;
    } else if (PercentFullRock[x + 1][y] > 0.0) {
      if (debug1c) printf("39 r\n");
      NextRockX = x + 1;
      NextRockY = y;
      BehindRockX = NextRockX;
      BehindRockY = NextRockY + 1;
      return;
    } else if (PercentFullRock[x + 1][y + 1] > 0.0) {
      if (debug1c) printf("40 r\n");
      NextRockX = x + 1;
      NextRockY = y + 1;
      BehindRockX = NextRockX;
      BehindRockY = NextRockY + 1;
      return;
    } else if (PercentFullRock[x][y + 1] > 0.0)
    /* RCL messed with 41, all the 42.X cases */
    {
      if (debug1b) printf("41 r\n");
      NextRockX = x;
      NextRockY = y + 1;
      BehindRockX = NextRockX - 1;
      BehindRockY = NextRockY;
      return;
    } else if (PercentFullRock[x - 1][y + 1] > 0.0) {
      if (debug1b) printf("42 r\n");
      NextRockX = x - 1;
      NextRockY = y + 1;
      BehindRockX = NextRockX - 1;
      BehindRockY = NextRockY;
      return;
    } else if (PercentFullRock[x - 1][y - 1] > 0.0) {
      if (debug1b) printf("42.1 r\n");
      NextRockX = x - 1;
      NextRockY = y - 1;
      BehindRockX = NextRockX + 1;
      BehindRockY = NextRockY - 1;
      return;
    } else if (PercentFullRock[x][y + 1] < 1.0)
    /* RCL added this case */
    {
      if (debug1b) printf("42.75 r\n");
      NextRockX = x;
      NextRockY = y + 1;
      BehindRockX = NextRockX - 1;
      BehindRockY = NextRockY;
      return;
    } else if (PercentFullRock[x - 1][y + 1] < 1.0)
    /* RCL added this case */
    {
      if (debug1b) printf("42.8 r\n");
      NextRockX = x - 1;
      NextRockY = y + 1;
      BehindRockX = NextRockX - 1;
      BehindRockY = NextRockY;
      return;
    } else {
      printf("Should've found next rock cell (Below): %d, %d \n", x, y);
      PauseRun(x, y, -1);
    }
  }

  else if ((XRock[z - 1] == XRock[z] - 1) && (YRock[z - 1] == YRock[z] - 1))
  /* came from lower left */
  {
    if (PercentFullRock[x + 1][y - 1] > 0.0) {
      if (debug1c) printf("43 r\n");
      NextRockX = x + 1;
      NextRockY = y - 1;
      BehindRockX = NextRockX;
      BehindRockY = NextRockY + 1;
      return;
    } else if (PercentFullRock[x + 1][y] > 0.0) {
      if (debug1c) printf("44 r\n");
      NextRockX = x + 1;
      NextRockY = y;
      BehindRockX = NextRockX;
      BehindRockY = NextRockY + 1;
      return;
    } else if (PercentFullRock[x + 1][y + 1] > 0.0) {
      if (debug1c) printf("45 r\n");
      NextRockX = x + 1;
      NextRockY = y + 1;
      BehindRockX = NextRockX - 1;
      BehindRockY = NextRockY;
      return;
    } else if (PercentFullRock[x][y + 1] > 0.0) {
      if (debug1c) printf("46 r\n");
      NextRockX = x;
      NextRockY = y + 1;
      BehindRockX = NextRockX - 1;
      BehindRockY = NextRockY;
      return;
    } else if (PercentFullRock[x - 1][y + 1] > 0.0) {
      if (debug1c) printf("47 r\n");
      NextRockX = x - 1;
      NextRockY = y + 1;
      BehindRockX = NextRockX - 1;
      BehindRockY = NextRockY;
      return;
    } else if (PercentFullRock[x - 1][y] > 0.0) {
      if (debug1c) printf("48 r\n");
      NextRockX = x - 1;
      NextRockY = y;
      BehindRockX = NextRockX - 1;
      BehindRockY = NextRockY;
      return;
    } else
      printf("Should've found next Rock cell (Lower Left): %d, %d \n", x, y);
  }

  printf("Should've found next Rock cell : %d, %d \n", x, y);
}

void RockCalculations(void)
/* LMV */
/* For each Rock cell j, moves along the beach and determines distance from rock
   to beach */
/* and amount of rock to weather */
/* This function calls FindNearestBeach and WeatherRock */
/* This function will use but not adjust the variable TotalRockCells */
/* Determines the total amount of rock weathered in one time step */
{
  int j;
  float TotalAmountWeathered = 0.0;

  for (j = 0; j <= TotalRockCells; j++) {
    FindNearestBeach(j);
    WeatherRock(j);
    TotalAmountWeathered += AmountWeathered[j];
  }
  /* printf("current_time_step = %d\n", current_time_step) */;
}

void FindNearestBeach(int j)
/* LMV */
/* Calculates the Euclidean distance between Rock and Beach for entire beach
   array      */
/* Finds the minimum beach to rock distance */
/* This function will use but not adjust the global arrays: X[], Y[], XRock[],
   YRock[]  */
/* PercentFullSand[][]                  */
/* This function will use but not adjust the variables MinDistToBeach and
   TotalBeachCells */
{
  int i;
  float DeltaX, DeltaY;
  int LookStart, LookEnd; /* for short rock to beach sweep */

  MinDistToBeach[j] = BIG_DISTANCE_TO_BEACH;

  if (current_time_step % TIME_TO_SWEEP_FULL_BEACH == 0)
  /* do a full sweep every now and then */
  {
    for (i = 0; i <= TotalBeachCells; i++) { /*full sweep */
      if (X[i] != XRock[j]) {                /*vertical distance */
        DeltaX =
            ((X[i] - XRock[j] - 1) + (PercentFullSand[X[i]][Y[i]] +
                                      PercentFullSand[XRock[j]][YRock[j]])) *
            CELL_WIDTH;
        if (debug10)
          printf("Rock cell = %d, Beach cell = %d, DeltaX (1)= %f \n", XRock[j],
                 X[i], DeltaX);
      }

      else {
        DeltaX = (PercentFullSand[X[i]][Y[i]]) * CELL_WIDTH;
        if (debug10)
          printf("Rock cell = %d, Beach cell = %d, DeltaX (2)= %f \n", XRock[j],
                 X[i], DeltaX);
      }

      if (YRock[j] <= Y[i]) { /*horizontal distance */
        /*if (XRock[j] = X[i])
           {
           DeltaY = 0.0;
           if (debug10) printf("DeltaY (0) = %f \n", DeltaY);
           }
           else
           { */
        DeltaY = (Y[i] - YRock[j]) * CELL_WIDTH;
        if (debug10) printf("DeltaY (1) = %f \n", DeltaY);
        /*} */
      } else {
        DeltaY = (YRock[j] - Y[i]) * CELL_WIDTH;
        if (debug10) printf("DeltaY (2) = %f \n", DeltaY);
      }

      DistanceToBeach[j] = sqrt(
          DeltaX * DeltaX + DeltaY * DeltaY); /*RCL removed Raise w/neg base */
      if (debug10)
        printf("Distance to beach from XRock[j] %d to X[i] %d is %f \n",
               XRock[j], X[i], DistanceToBeach[j]);

      if (DistanceToBeach[j] < MinDistToBeach[j]) {
        MinDistToBeach[j] = DistanceToBeach[j];
        ClosestBeach[j] = i;
      }
    }
    if (debug11)
      printf("Min Distance for j = %d is %f \n", j, MinDistToBeach[j]);
    if (debug11) printf("XRock[j] %d to X[i] %d \n", XRock[j], X[i]);
    if (debug11)
      printf("Closest Beach is at cell i = %d \n \n", ClosestBeach[j]);
    if (Y[i] == 150)
      printf("Min Distance for j = %d is %f \n", j, MinDistToBeach[j]);
    if (Y[i] == 150) printf("XRock[j] %d to X[i] %d \n", XRock[j], X[i]);
    if (Y[i] == 150)
      printf("Closest Beach is at cell i = %d \n \n", ClosestBeach[j]);

  }

  else
  /* LOOK_DIST allows the program to skip the full beach sweep at every iteration
     by looking at */
  /* the location of the previous closest beach and sweeping out a small arc
     plus or minus LOOK_DIST */
  {
    if (ClosestBeach[j] - LOOK_DIST >= 0)
      LookStart = ClosestBeach[j] - LOOK_DIST;
    else /* Only sweep right when on left end */
      LookStart = 0;

    if (ClosestBeach[j] + LOOK_DIST <= TotalBeachCells)
      LookEnd = ClosestBeach[j] + LOOK_DIST;
    else /* Only sweep left when at right end */
      LookEnd = TotalBeachCells;

    for (i = LookStart; i < LookEnd; i++) {
      if (X[i] != XRock[j]) {
        DeltaX =
            ((X[i] - XRock[j] - 1) + (PercentFullSand[X[i]][Y[i]] +
                                      PercentFullSand[XRock[j]][YRock[j]])) *
            CELL_WIDTH;
        if (debug10)
          printf("Rock cell = %d, Beach cell = %d, DeltaX (1)= %f \n", XRock[j],
                 X[i], DeltaX);
      }

      else {
        DeltaX = (PercentFullSand[X[i]][Y[i]]) * CELL_WIDTH;
        if (debug10)
          printf("Rock cell = %d, Beach cell = %d, DeltaX (2)= %f \n", XRock[j],
                 X[i], DeltaX);
      }

      if (YRock[j] <= Y[i]) { /*horizontal distance */
        /*if (XRock[j] = X[i])
           {
           DeltaY = 0.0;
           if (debug10) printf("DeltaY (0) = %f \n", DeltaY);
           }
           else
           { */
        DeltaY = (Y[i] - YRock[j]) * CELL_WIDTH;
        if (debug10) printf("DeltaY (1) = %f \n", DeltaY);
        /*} */
      } else {
        DeltaY = (YRock[j] - Y[i]) * CELL_WIDTH;
        if (debug10) printf("DeltaY (2) = %f \n", DeltaY);
      }

      DistanceToBeach[j] = sqrt(
          DeltaX * DeltaX + DeltaY * DeltaY); /*RCL removed Raise w/neg base */
      if (debug10)
        printf("Distance to beach from XRock[j] %d to X[i] %d is %f \n",
               XRock[j], X[i], DistanceToBeach[j]);

      if (DistanceToBeach[j] < MinDistToBeach[j]) {
        MinDistToBeach[j] = DistanceToBeach[j];
        ClosestBeach[j] = i;
      }
    }
    if (debug11)
      printf("Short Sweep - Min Distance for j = %d is %f \n", j,
             MinDistToBeach[j]);
    if (debug11) printf("XRock[j] %d to X[i] %d \n", XRock[j], X[i]);
    if (debug11) printf("Closest Beach is at cell i = %d \n ", ClosestBeach[j]);
    if (debug11) printf("X[i]: %d, Y[i]: %d \n\n", X[i], Y[i]);
  }
}

void WeatherRock(int j)
/* LMV */
/* Function changes PercentFullRock[][] and PercentFullSand[][] along the rock
   beach interface only	*/
/* At this point, weathering rate is pretty arbitrary
   */
/* Function also changes AllRock[][]
   */
/* Weathering rate changes along the shore, 'f' = fast weather, 's' = slow
   weathering			*/
/* Function uses (but does not change) type_of_rock[][] to determine rate
   */
{
  double WeatheringRatePerYear, CurrentWeatherCoeff, mslope, wscale,
      CurrentPercentFine;

  /* Assign maximum weathering rate depending on rock type  */
  CurrentWeatherCoeff =
      (type_of_rock[XRock[j]][YRock[j]] == 'f' ? FastWeatherCoeff
                                             : SlowWeatherCoeff);

  if (CurrentWeatherCoeff != FastWeatherCoeff &&
      CurrentWeatherCoeff != SlowWeatherCoeff) {
    CurrentWeatherCoeff =
        (type_of_rock[XRockBehind[j]][YRockBehind[j]] == 'f' ? FastWeatherCoeff
                                                           : SlowWeatherCoeff);
  }
  if (CurrentWeatherCoeff != FastWeatherCoeff &&
      CurrentWeatherCoeff != SlowWeatherCoeff) {
    printf("weatherrock (%d): rock & behind both have no rocktype\n", j);
  }

  /* Assign proportion of fines lost PWL */
  CurrentPercentFine =
      (type_of_rock[XRock[j]][YRock[j]] == 'f' ? PercentFineFast
                                             : PercentFineSlow);
  if (CurrentPercentFine != PercentFineFast &&
      CurrentPercentFine != PercentFineSlow)
    CurrentPercentFine =
        (type_of_rock[XRockBehind[j]][YRockBehind[j]] == 'f' ? PercentFineFast
                                                           : PercentFineSlow);
  /* Flip current percent fines to make calculations below easier PWL */
  CurrentPercentFine = 1 - CurrentPercentFine;

  /* Determine weathering: abrasion or no abrasion? */
  if ((MinDistToBeach[j]) <=
      NO_WEATHERING) { /*think in terms of vertical equivalence */

    if (ABRASION == 'n') {
      WeatheringRatePerYear =
          CurrentWeatherCoeff * (exp(-MinDistToBeach[j])); /*exponential
                                                              weathering (in
                                                              meters/yr) based
                                                              on sed cover */
      /*printf("weathering rate: %f", WeatheringRatePerYear); */
    }

    if (ABRASION == 'y') {
      if (MinDistToBeach[j] <= W_CRIT) {
        mslope = (CurrentWeatherCoeff * (N - 1)) /
                 W_CRIT; /* Slope of line between N*BareRock and BareRock PWL */
        WeatheringRatePerYear =
            CurrentWeatherCoeff + mslope * MinDistToBeach[j];
        /*if (!debug12b) printf("\nabrasion weathering rate:
         * %f\n",WeatheringRatePerYear); */
      }

      else {
        wscale =
            (NO_WEATHERING - W_CRIT) /
            (log((0.01 / N))); /* Decay coefficient for the "cover" portion of
                                  abrasion curve PWL */
        WeatheringRatePerYear = kAngleFactor * CurrentWeatherCoeff *
                                (exp((MinDistToBeach[j] - W_CRIT) / wscale));
      }
    }
  }

  else {
    WeatheringRatePerYear = 0;
  }

  AmountWeathered[j] = ((WeatheringRatePerYear * time_step) / 365) / CELL_WIDTH;
  /* was AmountWeathered[j] = ((WeatheringRatePerYear *
   * time_step)/365)/CELL_WIDTH; */
  /* weathering per time step. units: cell units, or a percent weathered (which
   * correlates with erosion) */

  if (debug12)
    printf(
        "for %d (%d,%d) rock %c, angfactor %f, weathercoeff %f, mindist %f, "
        "rate %f m/yr, frac weathered %f\n",
        j, XRock[j], YRock[j], type_of_rock[XRock[j]][YRock[j]], kAngleFactor,
        CurrentWeatherCoeff, MinDistToBeach[j], WeatheringRatePerYear,
        AmountWeathered[j]);

  /* PWL, 2/16/12: added sea cliffs, which puts in an extra amount of sediment
   * to PercentFullSand. Use variable CliffHeight (top) to adjust. */
  /* Also, need to make sure model can deal with situations when there is no
   * rock behind current rock cell (i.e., back side of an island)   */
  if (PercentFullRock[XRock[j]][YRock[j]] >= AmountWeathered[j]) {
    if (debug12) printf("%d (%d,%d) weathers (1)\n", j, XRock[j], YRock[j]);
    if (debugtopo)
      printf("\n Cliff height (1): %f, j: %d\n", topography[XRock[j]][YRock[j]],
             j);

    PercentFullRock[XRock[j]][YRock[j]] -= AmountWeathered[j];
    PercentFullSand[XRock[j]][YRock[j]] +=
        AmountWeathered[j] +
        (AmountWeathered[j] *
         (topography[XRock[j]][YRock[j]] / DEPTH_SHOREFACE) *
         CurrentPercentFine);

    if (PercentFullRock[XRock[j]][YRock[j]] < 1.0)
      AllRock[XRock[j]][YRock[j]] = 'n';
  }

  else if ((PercentFullRock[XRock[j]][YRock[j]] +
            PercentFullRock[XRockBehind[j]][YRockBehind[j]]) >=
           AmountWeathered[j]) {
    if (debug12) printf("%d (%d,%d) weathers (2)\n", j, XRock[j], YRock[j]);
    if (debugtopo)
      printf("\n Cliff height (2): %f\n", topography[XRock[j]][YRock[j]]);
    if (debugtopo)
      printf("\n Cliff height behind (2): %f\n",
             topography[XRockBehind[j]][YRockBehind[j]]);

    PercentFullRock[XRockBehind[j]][YRockBehind[j]] -=
        (AmountWeathered[j] - PercentFullRock[XRock[j]][YRock[j]]);

    PercentFullSand[XRockBehind[j]][YRockBehind[j]] +=
        (AmountWeathered[j] - PercentFullRock[XRock[j]][YRock[j]]) +
        ((AmountWeathered[j] - PercentFullRock[XRock[j]][YRock[j]]) *
         (topography[XRockBehind[j]][YRockBehind[j]] / DEPTH_SHOREFACE) *
         CurrentPercentFine);

    PercentFullSand[XRock[j]][YRock[j]] +=
        PercentFullRock[XRock[j]][YRock[j]] +
        (PercentFullRock[XRock[j]][YRock[j]] *
         (topography[XRock[j]][YRock[j]] / DEPTH_SHOREFACE) *
         CurrentPercentFine);

    PercentFullRock[XRock[j]][YRock[j]] = 0.0;
    AllBeach[XRockBehind[j]][YRockBehind[j]] = 'n';
  } else { /* current + behind don't have enough rock to weather,
              so weather all there is in both cells. */
    if (debug12) printf("%d (%d,%d) weathers (3)\n", j, XRock[j], YRock[j]);
    if (debugtopo)
      printf("\n Cliff height (3): %f\n", topography[XRock[j]][YRock[j]]);
    if (debugtopo)
      printf("\n Cliff height behind (3): %f\n",
             topography[XRockBehind[j]][YRockBehind[j]]);

    if (YRock[j] > 1 && YRock[j] < 2 * Y_MAX - 2) {
      if (debug12)
        printf("not enough rock to weather for %d (%d,%d). ", j, XRock[j],
               YRock[j]);
      if (debug12)
        printf("this is handle-able but prob. shouldn't happen often.\n");
      /*PauseRun(XRock[j],YRock[j],j); */
    }
    PercentFullSand[XRock[j]][YRock[j]] +=
        PercentFullRock[XRock[j]][YRock[j]] +
        (PercentFullRock[XRock[j]][YRock[j]] *
         (topography[XRock[j]][YRock[j]] / DEPTH_SHOREFACE) *
         CurrentPercentFine);

    PercentFullRock[XRock[j]][YRock[j]] = 0.0;

    PercentFullSand[XRockBehind[j]][YRockBehind[j]] +=
        PercentFullRock[XRockBehind[j]][YRockBehind[j]] +
        (PercentFullRock[XRockBehind[j]][YRockBehind[j]] *
         (topography[XRockBehind[j]][YRockBehind[j]] / DEPTH_SHOREFACE) *
         CurrentPercentFine);

    PercentFullRock[XRockBehind[j]][YRockBehind[j]] = 0.0;
  }
}

void ShadowSweep(void)
/*  Moves along beach and tests to see if cells are in shadow 		*/
/*  This function will use and determine the Global array:  InShadow[]	*/
/*  This function will use and adjust the variable:   ShadowXMax        */
/*  This function will use but not adjust the variable:  TotalBeachCells */
{
  int i;

  /* Find maximum extent of beach to use as a limit for shadow searching */

  /*ShadowXMax = XMaxBeach(ShadowXMax) + 3; */
  /* if(ShadowXMax > X_MAX) ShadowXMax = X_MAX; RCL */

  if (debug2)
    printf("ShadowXMax: %d   XMaxBeach: %d \n", ShadowXMax,
           XMaxBeach(ShadowXMax));

  /* Determine if beach cells are in shadow */
  for (i = 0; i <= TotalBeachCells; i++) {
    if (kSwanFlag == 'n') /* Use shadows when SWAN is not involved */
    {
      InShadow[i] = FindIfInShadow(X[i], Y[i], ShadowXMax);

      /* Code to see what is happening */

      if (InShadow[i] == 'y') {
        if (debug2)
          printf("Shadow at X: %3d  Y: %3d  i: %3d \n", X[i], Y[i], i);
        else if (debug2)
          printf("No Shadow X: %3d  Y: %3d  i: %3d \n", X[i], Y[i], i);
      }
    } else /* When using SWAN, turn off shadows... */
    {
      InShadow[i] = 'n';
    }
  }
}

int XMaxBeach(int Max)

/* Finds extent of beach in x direction                                 */
/* Starts searching at a point 3 rows higher than input Max             */
/* Function returns integer value equal to max extent of 'allbeach'     */
{
  int xtest, ytest;

  xtest = Max + 3;
  ytest = 0;

  while (xtest > 0) {
    while (ytest < 2 * Y_MAX) {
      if (AllBeach[xtest][ytest] == 'y') {
        return xtest;
      }

      ytest++;
    }

    ytest = 0;

    xtest--;
  }

  printf("***** Should've found X_MAX for shadow): %d, %d ***** \n", xtest,
         ytest);

  return X_MAX;
}

char FindIfInShadow(int xin, int yin, int ShadMax)
/*  Function to determine if particular cell xin,yin is in shadow */
/*  Returns a character 'y' or 'n' to indicate so */
/*  This function will use but not affect the global arrays: */
/*      AllBeach[][] and PercentFull[][] */
/*  This function refers to global variables:  WaveAngle, ShadowStepDistance */
{
  /* Initialize local variables */

  float xdistance, ydistance;

  int xcheck = xin;
  int ycheck = yin;
  int iteration = 1;

  while ((xcheck < ShadMax) && (ycheck > -1) && (ycheck < (2 * Y_MAX))) {
    /*  Find a cell along the projection line moving against wave direction */

    xdistance = iteration * ShadowStepDistance * cos(WaveAngle);
    xcheck = xin + rint(xdistance);

    ydistance = -iteration * ShadowStepDistance * sin(WaveAngle);
    ycheck = yin + rint(ydistance);

    if (debug2a) {
      printf("\n Iteration: %d \n", iteration);
      printf("Xdist : %f  ", xdistance);
      printf("    Xcheck: %d \n", xcheck);
      printf("Ydist : %f  ", ydistance);
      printf("    Ycheck: %d \n", ycheck);
      printf("firsttest : %f  ", (xcheck + .2));
      printf("secondtest : %f  ",
             (xin + 0.3 + fabs((ycheck - yin) / tan(WaveAngle))));
      printf("term : %f  \n", fabs((ycheck - yin) / tan(WaveAngle)));
    }

    if (xcheck >= X_MAX - 1 || ycheck >= 2 * Y_MAX - 1) return 'n';

    /* If AllBeach is along the way, and not next neighbor                  */
    /* Probably won't get to this one, though                               */

    if ((AllBeach[xcheck][ycheck] == 'y') &&
        (/* (abs(ycheck-yin) !=1 && abs(xcheck-xin) !=1) && */
         ((xcheck + 1) >
          (xin + (PercentFullSand[xin][yin] + PercentFullRock[xin][yin]) +
           fabs((ycheck - yin) / tan(WaveAngle))))) &&
        (ycheck > -1) && (ycheck < 2 * Y_MAX - 1)) {
      if (debug2) {
        printf("Allbeach used Xck: %2d  Yck: %2d   \n", xcheck, ycheck);
        printf("	      Xin: %2d  Yin: %2d   ", xin, yin);
        printf("Xdist : %d  ", abs(xcheck - xin));
        printf("Ydist : %d  ", abs(ycheck - yin));
        printf("firsttest : %f  ", (xcheck + .2));
        printf("secondtest : %f  ",
               (xin + 0.3 + fabs((ycheck - yin) / tan(WaveAngle))));
        printf("term : %f  \n\n", fabs((ycheck - yin) / tan(WaveAngle)));
      }
      if (debug2a)
        printf("finishing on a full cell with Xcheck %d, Xin %d, PFS+R %f  \n",
               xcheck, xin, (PercentFullSand[xcheck][ycheck] +
                             PercentFullSand[xcheck][ycheck]));
      return 'y';
    }
/* Compare a partially full cell's x - distance to a line projected     */
/* from the starting beach cell's x-distance                            */
/* This assumes that beach projection is in x-direction (not too bad)   */
/* RCL this could be contributing to the jaggedness...need to take into account
the way it's facing
or could just ignore partially full cells completely...
*/
#ifdef FGFG
    else if (((PercentFullSand[xcheck][ycheck] +
               PercentFullRock[xcheck][ycheck]) > 0) &&
             ((xcheck + (PercentFullSand[xcheck][ycheck] +
                         PercentFullRock[xcheck][ycheck])) >
              (xin + (PercentFullSand[xin][yin] + PercentFullRock[xin][yin]) +
               fabs((ycheck - yin) / tan(WaveAngle)))) &&
             (ycheck > -1) && (ycheck < 2 * Y_MAX - 1)
             /* && !(xcheck == xin) */
             ) {
      if (debug2a)
        printf(
            "finishing on a partially full with Xcheck %d, Xin %d, PFS+R %f  "
            "\n",
            xcheck, xin, (PercentFullSand[xcheck][ycheck] +
                          PercentFullSand[xcheck][ycheck]));
      return 'y';
    }
#endif
    iteration++;
  }

  /*  No shadow was found */
  return 'n';
}

void DetermineAngles(void)

/*  Function to determine beach angles for all beach cells from left to right */
/*  By convention, the ShorelineAngle will apply to current cell and right
   neighbor     */
/*  This function will determine global arrays: */
/*              ShorelineAngle[], UpWind[], SurroundingAngle[] */
/*  This function will use but not affect the following arrays and values: */
/*              X[], Y[], PercentFull[][], AllBeach[][], WaveAngle */
/*  ADA Revised underside, SurroundingAngle 6/03 */
/*  ADA Rerevised, 3/04 */
{
  int i, j, k; /* Local loop variables */

  /* Compute ShorelineAngle[]  */
  /*      not equal to TotalBeachCells because angle between cell and rt
   * neighbor */

  for (i = 0; i < TotalBeachCells; i++) {
    /* function revised 1-04 */
    if (Y[i] > Y[i + 1])
    /* On bottom side of Spit (assuming we must be if going right)  */
    {
      /* math revised 6-03 */

      ShorelineAngle[i] =
          PI - atan((X[i + 1] - (PercentFullSand[X[i + 1]][Y[i + 1]] +
                                 PercentFullRock[X[i + 1]][Y[i + 1]])) -
                    (X[i] - (PercentFullSand[X[i]][Y[i]]) +
                     PercentFullRock[X[i]][Y[i]]));

      if (ShorelineAngle[i] > PI) {
        ShorelineAngle[i] -= 2.0 * PI;
        if (debug3) printf("Next under -180 \n");
      }

      if (debug3)
        printf(
            "(4a) i = %d  X[i]: %d Y[i]: %d Percent %3f \n           X[+]: %d "
            "Y[+]: %d Percent{+} %3f   Angle: %f  Deg Angle: %f \n",
            i, X[i], Y[i],
            (PercentFullSand[X[i]][Y[i]] + PercentFullRock[X[i]][Y[i]]),
            X[i + 1], Y[i + 1], (PercentFullSand[X[i + 1]][Y[i + 1]] +
                                 PercentFullRock[X[i + 1]][Y[i + 1]]),
            ShorelineAngle[i], ShorelineAngle[i] * 180 / PI);

      /*PauseRun(X[i],Y[i],i); */
    }

    else if ((Y[i] < Y[i + 1]))
    /* function revised 1-04 - asusme if goin right, on regular shore */
    /*  'regular' beach */
    {
      ShorelineAngle[i] = atan(
          (X[i + 1] + (PercentFullSand[X[i + 1]][Y[i + 1]] +
                       PercentFullRock[X[i + 1]][Y[i + 1]])) -
          (X[i] + (PercentFullSand[X[i]][Y[i]] + PercentFullRock[X[i]][Y[i]])));

      if (debug3)
        printf(
            "(1) i = %d  X[i]: %d Y[i]: %d Percent %3f \n           X[+]: %d "
            "Y[+]: %d Percent[+] %3f   Angle: %f  Deg Angle: %f \n",
            i, X[i], Y[i],
            (PercentFullSand[X[i]][Y[i]] + PercentFullRock[X[i]][Y[i]]),
            X[i + 1], Y[i + 1], (PercentFullSand[X[i + 1]][Y[i + 1]] +
                                 PercentFullRock[X[i + 1]][Y[i + 1]]),
            ShorelineAngle[i], ShorelineAngle[i] * 180 / PI);
    }

    else if (Y[i] == Y[i + 1] && X[i] > X[i + 1])
    /*  Shore up and down, on right side */
    {
      ShorelineAngle[i] =
          -PI / 2.0 + atan((Y[i + 1] + (PercentFullSand[X[i + 1]][Y[i + 1]] +
                                        PercentFullRock[X[i + 1]][Y[i + 1]])) -
                           (Y[i] + (PercentFullSand[X[i]][Y[i]] +
                                    PercentFullRock[X[i]][Y[i]])));

      if (debug3)
        printf(
            "(2) i = %d  X[i]: %d Y[i]: %d Percent %3f \n           X[+]: %d "
            "Y[+]: %d Percent[+] %3f   Angle: %f  Deg Angle: %f \n",
            i, X[i], Y[i],
            (PercentFullSand[X[i]][Y[i]] + PercentFullRock[X[i]][Y[i]]),
            X[i + 1], Y[i + 1], (PercentFullSand[X[i + 1]][Y[i + 1]] +
                                 PercentFullRock[X[i + 1]][Y[i + 1]]),
            ShorelineAngle[i], ShorelineAngle[i] * 180 / PI);
    }

    else if (Y[i] == Y[i + 1] && X[i] < X[i + 1])
    /* Shore up and down, on left side */
    {
      ShorelineAngle[i] =
          PI / 2.0 + atan((Y[i + 1] + (PercentFullSand[X[i + 1]][Y[i + 1]] +
                                       PercentFullRock[X[i + 1]][Y[i + 1]])) -
                          (Y[i] + (PercentFullSand[X[i]][Y[i]] +
                                   PercentFullRock[X[i]][Y[i]])));

      if (debug3)
        printf(
            "(3) i = %d  X[i]: %d Y[i]: %d Percent %3f \n           X[+]: %d "
            "Y[+]: %d Percent[+] %3f   Angle: %f  Deg Angle: %f \n",
            i, X[i], Y[i],
            (PercentFullSand[X[i]][Y[i]] + PercentFullRock[X[i]][Y[i]]),
            X[i + 1], Y[i + 1], (PercentFullSand[X[i + 1]][Y[i + 1]] +
                                 PercentFullRock[X[i + 1]][Y[i + 1]]),
            ShorelineAngle[i], ShorelineAngle[i] * 180 / PI);
    }

    else {
      printf("Should've found ShorelineAngle): %d, %d \n", X[i], Y[i]);
      PauseRun(X[i], Y[i], i);
    }

    /* printf("i = %d  X[i]: %d Y[i]: %d Percent %3f \n           X[+]: %d Y[+]:
       %d Percent %3f   Angle: %f  Deg Angle: %f \n",i, X[i], Y[i],
       (PercentFullSand[X[i]][Y[i]] + PercentFullRock[X[i]][Y[i]]) ,
       X[i+1], Y[i+1], (PercentFullSand[X[i+1]][Y[i+1]] +
       PercentFullRock[X[i+1]][Y[i+1]]), ShorelineAngle[i],
       ShorelineAngle[i]*180/PI); */
  }

  for (k = 1; k < TotalBeachCells; k++) {
    /* compute SurroundingAngle array */
    /* 02/04 AA averaging doesn't work on bottom of spits */
    /* Use trick that x is less if on bottom of spit - angles must be different
     * signs as well */

    if ((Y[k - 1] - Y[k + 1] == 2) &&
        (copysign(ShorelineAngle[k - 1], ShorelineAngle[k]) !=
         ShorelineAngle[k - 1])) {
      SurroundingAngle[k] =
          (ShorelineAngle[k - 1] + ShorelineAngle[k]) / 2 + PI;
      if (SurroundingAngle[k] > PI) {
        SurroundingAngle[k] -= 2.0 * PI;
      }
      if (debug4) printf("Under: %d\n", k);
    } else {
      SurroundingAngle[k] = (ShorelineAngle[k - 1] + ShorelineAngle[k]) / 2;
    }
  }

  /* Determine Upwind/downwind condition */
  /* Note - Surrounding angle is based upon left and right cell neighbors, */
  /* and is centered on cell, not on right boundary */

  if (debug4) printf("\nUp/Down   Wave Angle:%f\n", WaveAngle * RAD_TO_DEG);

  for (j = 1; j < TotalBeachCells; j++) {
    if (debug4)
      printf("i: %d  Shad: %c Ang[i]: %3.1f  Sur: %3.1f  Effect: %3f  ", j,
             InShadow[j], ShorelineAngle[j] * RAD_TO_DEG,
             SurroundingAngle[j] * RAD_TO_DEG,
             (WaveAngle - SurroundingAngle[j]) * RAD_TO_DEG);

    if (fabs(WaveAngle - SurroundingAngle[j]) >= 42.0 / RAD_TO_DEG) {
      UpWind[j] = 'u';
      if (debug4) printf("U(1)  ");
    } else {
      UpWind[j] = 'd';
      if (debug4) printf("D(1)  ");
    }

    if (debug4) printf("\n");
  }
}

void DetermineSedTransport(void)

/*  Loop function to determine which neigbor/situation to use for sediment
   transport calcs      */
/*  Once situation is determined, will use function SedTrans to determine actual
   transport      */
/*  This function will call SedTrans which will determine global array:
   VolumeAcrossBorder[] LMV        */
/*  This function will use but not affect the following arrays and values: */
/*              X[], Y[], InShadow[], UpWind[], ShorelineAngle[] */
/*              PercentFullSand[][], AllBeach[][], WaveAngle */
{
  int i; /* Loop variable LMV also sent to SedTrans to determine VolAcrossBorder
            */
  float ShoreAngleUsed; /* Temporary holder for shoreline angle */
  int CalcCell; /* Cell sediment coming from to go across boundary i */
  int Next, Last; /* Indicators so test can go both left/right */
  int Correction; /* Term needed for shoreline angle and i+1 case, angle stored
                     at i      */
  char UpWindLocal; /* Local holder for upwind/downwind condition */
  char MaxTrans; /* Do we need to compute using maximum transport ? */
  int DoFlux; /* Skip sed transport calcs (added 02/04 AA) */
  float DummyAngle = 1.;

  if (debug5)
    printf("\nSEDTRANS: %d  @  %f \n\n", current_time_step, WaveAngle * RAD_TO_DEG);

  for (i = 1; i < TotalBeachCells - 1; i++) {
    if (debug5) printf("\n  i: %d  ", i);

    MaxTrans = 'n';

    /*  Is littoral transport going left or right?  */
    /* #SWAN, 11/27/14: this isn't necessarily correct when SWAN is involved --
*short-period waves can refract around and break going the opposite direction of
* their deep-water direction.  */
    /* Parse SWAN here first, and make this determination with 'breaking wave'
     * angles */
    /* OY VAY */
    if (kSwanFlag == 'y') {
      ParseSWAN(i, DummyAngle); /* Use dummy shore angle because we don't need
                                   it right now...prob a poor solution to the
                                   problem */
    }

    if (kSwanFlag == 'y') {
      if ((Angle - ShorelineAngle[i]) >
          0) /* Use SWAN nearshore angle instead */
      {
        /*  Transport going right, center on cell to left side of border
         */
        /*  Next cell in positive direction, no correction term needed */
        CalcCell = i;
        Next = 1;
        Last = -1;
        Correction = 0;
        DirectionAcrossBorder[i] = 'r'; /*LMV*/
        if (debug5) printf("RT  ");
      } else {
        /*  Transport going left, center on cell to right side of border
         */
        /*  Next cell in negative direction, correction term needed */
        CalcCell = i + 1;
        Next = -1;
        Last = 1;
        Correction = -1;
        DirectionAcrossBorder[i] = 'l'; /*LMV*/
        if (debug5) printf("LT  ");
      }
    }

    else if (kSwanFlag == 'n') /* #SWAN */
    {
      if ((WaveAngle - ShorelineAngle[i]) > 0) {
        /*  Transport going right, center on cell to left side of border
         */
        /*  Next cell in positive direction, no correction term needed */
        CalcCell = i;
        Next = 1;
        Last = -1;
        Correction = 0;
        DirectionAcrossBorder[i] = 'r'; /*LMV*/
        if (debug5) printf("RT  ");
      } else {
        /*  Transport going left, center on cell to right side of border
         */
        /*  Next cell in negative direction, correction term needed */
        CalcCell = i + 1;
        Next = -1;
        Last = 1;
        Correction = -1;
        DirectionAcrossBorder[i] = 'l'; /*LMV*/
        if (debug5) printf("LT  ");
      }
    }

    if (InShadow[CalcCell] == 'n') {
      /*  Adjustment for maximum transport when passing through 45 degrees */
      /*  This adjustment is only made for moving from downwind to upwind
       * conditions  */
      /*                                                                              */
      /*  purposefully done before shadow adjustment, only use maxtran when */
      /*      transition from dw to up not because of shadow */
      /* keeping transition from uw to dw - does not seem to be big deal (04/02
       * AA) */

      if (((UpWind[CalcCell] == 'd') && (UpWind[CalcCell + Next] == 'u') &&
           (InShadow[CalcCell + Next] == 'n')) ||
          ((UpWind[CalcCell + Last] == 'u') && (UpWind[CalcCell] == 'd') &&
           (InShadow[CalcCell + Last] == 'n'))) {
        MaxTrans = 'y';
        if (debug5) printf("MAXTRAN  ");
      }

      /*  Upwind/Downwind adjustment Make sure sediment is put into shadows */
      /*  If Next cell is in shadow, use UpWind condition */

      DoFlux = 1;
      UpWindLocal = UpWind[CalcCell];

      if (InShadow[CalcCell + Next] == 'y') {
        UpWindLocal = 'u';
        if (debug5) printf("U(2)  ");
      }

      /*  If coming out of shadow, downwind should be used            */
      /*  HOWEVER- 02/04 AA - if high angle, will result in same flux in/out
       * problem */
      /*      solution  - no flux for high angle waves */

      if ((InShadow[CalcCell + Last] == 'y') && (UpWindLocal == 'u')) {
        DoFlux = 0;
        if (debug5) printf("U(X) NOFLUX \n");
      }

      /*  Use upwind or downwind shoreline angle for calcs                    */

      if (UpWindLocal == 'u') {
        ShoreAngleUsed = ShorelineAngle[CalcCell + Last + Correction];
        if (debug5)
          printf("UP  ShoreAngle: %3.1f  ", ShoreAngleUsed * RAD_TO_DEG);
      } else if (UpWindLocal == 'd') {
        ShoreAngleUsed = ShorelineAngle[CalcCell + Correction];
        if (debug5)
          printf("DN  ShoreAngle: %3.1f  ", ShoreAngleUsed * RAD_TO_DEG);
      }

      /* !!! Do not do transport on unerneath c'cause it gets all messed up */
      if (fabs(ShoreAngleUsed) > PI / 2.0) {
        DoFlux = 0;
      }

      /* Send to SedTrans to calculate VolumeIn and VolumeOut */

      /* printf("i = %d  Cell: %d NextCell: %d Angle: %f Trans Angle: %f\n",
         i, CalcCell, CalcCell+Next, ShoreAngleUsed*180/PI, (WaveAngle -
         ShoreAngleUsed)*180/PI); */

      if (debug5)
        printf("From: %d  To: %d  TransAngle %3.1f", CalcCell, CalcCell + Next,
               (WaveAngle - ShoreAngleUsed) * RAD_TO_DEG);

      if (DoFlux) {
        /* SedTrans (i, ShoreAngleUsed, MaxTrans); */
        /*LMV*/
        SedTrans(i, CalcCell, ShoreAngleUsed, MaxTrans, Last);
      }
    }
  }
}

void SedTrans(int i, int From, float ShoreAngle, char MaxT, int Last)

/*  This central function will calcualte the sediment transported from the cell
   i-1 to          */
/*  the cell at i, using the input ShoreAngle   LMV */
/*  This function will caluclate and determine the global arrays: */
/*              VolumeAcrossBorder[]    LMV */
/*  This function will use the global values defining the wave field: */
/*      WaveAngle, Period, OffShoreWvHt */
/*  Revised 6/02 - New iterative calc for refraction and breaking, parameters
   revised           */
{
  /* Coefficients - some of these are important */

  float StartDepth = 3 * OffShoreWvHt; /* m, depth to begin refraction calcs
                                          (needs to be beyond breakers) */
  float RefractStep =
      .2; /* m, step size to iterate depth for refraction calcs */
  float KBreak = 0.5; /* coefficient for wave breaking threshold */
  float rho = 1020; /* kg/m3 - density of water and dissolved matter */

  /* Variables */

  float AngleDeep;          /* rad, Angle of waves to shore at inner shelf  */
  float Depth = StartDepth; /* m, water depth for current iteration         */
  float Angle;              /* rad, calculation angle                       */
  float CDeep;              /* m/s, phase velocity in deep water            */
  float LDeep;              /* m, offhsore wavelength                       */
  float C;                  /* m/s, current step phase velocity             */
  float kh;                 /* wavenumber times depth                       */
  float n;                  /* n                                            */
  float WaveLength;         /* m, current wavelength                        */

  /* Primary assumption is that waves refract over shore-parallel contours */
  /* New algorithm 6/02 iteratively takes wiave onshore until they break, then
   * computes Qs        */
  /* See notes 06/05/02 */

  if (debug6)
    printf("Wave Angle %2.2f Shore Angle  %2.2f    ", WaveAngle * RAD_TO_DEG,
           ShoreAngle * RAD_TO_DEG);

  AngleDeep = WaveAngle - ShoreAngle;

  if (MaxT == 'y') {
    AngleDeep = PI / 4.0;
  }
  if (debug6) printf("Deep Tranport Angle %2.2f \n\n", AngleDeep * RAD_TO_DEG);

  /*  Don't do calculations if over 90 degrees, should be in shadow  */

  if (AngleDeep > 0.995 * PI / 2.0 || AngleDeep < -0.995 * PI / 2.0) {
    return;
  } else if (kSwanFlag == 'n') {
    /* Calculate Deep Water Celerity & Length, Komar 5.11 c = gT / PI, L = CT */

    CDeep = GRAVITY * Period / (2.0 * PI);
    LDeep = CDeep * Period;
    if (debug6) printf("CDeep = %2.2f LDeep = %2.2f \n", CDeep, LDeep);

    while (TRUE) {
      /* non-iterative eqn for L, from Fenton & McKee                 */

      WaveLength =
          LDeep *
          Raise(tanh(Raise(Raise(2.0 * PI / Period, 2) * Depth / GRAVITY, .75)),
                2.0 / 3.0);
      C = WaveLength / Period;
      if (debug6)
        printf("DEPTH: %2.2f Wavelength = %2.2f C = %2.2f ", Depth, WaveLength,
               C);

      /* Determine n = 1/2(1+2kh/sinh(kh)) Komar 5.21                 */
      /* First Calculate kh = 2 pi Depth/L  from k = 2 pi/L           */

      kh = 2 * PI * Depth / WaveLength;
      n = 0.5 * (1 + 2.0 * kh / sinh(2.0 * kh));
      if (debug6) printf("kh: %2.3f  n: %2.3f ", kh, n);

      /* Calculate angle, assuming shore parallel contours and no conv/div of
       * rays    */
      /* from Komar 5.47 */

      Angle = asin(C / CDeep * sin(AngleDeep));
      if (debug6) printf("Angle: %2.2f", Angle * RAD_TO_DEG);

      /* Determine Wave height from refract calcs - Komar 5.49 */

      WvHeight = OffShoreWvHt *
                 Raise(CDeep * cos(AngleDeep) / (C * 2.0 * n * cos(Angle)), .5);
      if (debug6) printf(" WvHeight : %2.3f\n", WvHeight);

      if (WvHeight > Depth * KBreak || Depth <= RefractStep)
        break;
      else
        Depth -= RefractStep;
    }
  } else if (kSwanFlag == 'y') /* Do SWAN */
  {
    /* Instead of shoaling above, use ParseSWAN function to find the breaking
     * wave characteristics. */
    ParseSWAN(From, ShoreAngle);
  }

  /* Now Determine Transport */
  /* eq. 9.6b (10.8) Komar, including assumption of sed density = 2650 kg/m3 */
  /* additional accuracy here will not improve an already suspect eqn for sed
   * transport   */
  /* (especially with poorly constrained coefficients), */
  /* so no attempt made to make this a more perfect imperfection */
  if (kSwanFlag == 'n') {
    VolumeAcrossBorder[i] =
        fabs(1.1 * rho * Raise(GRAVITY, 3.0 / 2.0) * Raise(WvHeight, 2.5) *
             cos(Angle) * sin(Angle) * time_step); /*LMV - now global array*/
  } else /* Do Qs slightly different for SWAN input */
  {
    VolumeAcrossBorder[i] =
        fabs(0.2 * rho * Raise(GRAVITY, 3.0 / 2.0) * Raise(WvHeight, 2.5) *
             cos(Angle - ShoreAngle) * sin(Angle - ShoreAngle) * time_step);

    if (UpWind[From] == 'u') /* Use wave characteristics updrift */
    {
      VolumeAcrossBorder[i] = fabs(
          0.2 * rho * Raise(GRAVITY, 3.0 / 2.0) * Raise(Hsigdebug[From + Last], 2.5) *
          cos(Dirdebug[From + Last] - ShoreAngle) *
          sin(Dirdebug[From + Last] - ShoreAngle) * time_step);
    }
  }

  /*LMV VolumeIn/Out is now calculated below in AdjustShore */

  if (debug6a) printf("\ni: %d \n", i);
  if (debug6a) printf("VolumeAcrossBorder: %f  ", VolumeAcrossBorder[i]);
}

void DoSink(void)
/* empties of sand any cell that has been designated a sink */
{
  int s, i, foundOne = 0;
  if (!HAVE_SINKS || RandZeroToOne() > SINKINESS) return;

  if (COLUMN_SINKS) /* treat sink as any beach cell at a particular y-value */
    for (s = 0; s < NUM_SINKS; s++) {
      for (i = 0; i < MaxBeachLength; i++)
        if (Y[i] == kSinkY[s]) {
          PercentFullSand[X[i]][Y[i]] = 0.0;
          foundOne = 1;
        }
      if (!foundOne) printf("couldn't find sink %d at y == %d!\n", s, kSinkY[s]);
    }
  else /* treat sink as a specific cell; empty it if it is on the beach */
    for (s = 0; s < NUM_SINKS; s++) {
      for (i = 0; i < MaxBeachLength; i++)
        if (Y[i] == kSinkY[s] && X[i] == kSinkX[s]) {
          PercentFullSand[X[i]][Y[i]] = 0.0;
          foundOne = 1;
        }
      /*if(!foundOne) printf("couldn't find sink %d at
       * (%d,%d)!\n",s,kSinkX[s],kSinkY[s]); */
    }
  if (!foundOne) printf("didn't find any sinks @ time %d\n", current_time_step);
}

void FlowInCell(void)
/* LMV */
/* This function will determine if sediment transport into a cell is from right,
   from left,     */
/*                                              converging into the cell, or
   diverging out      */
/* Purpose is to address problems with insufficent sand, fix losing sand problem
   */
/* Moves from left to right direction */
/* This function will affect and determine the global arrays:  FlowThroughCell[]
   */
/* This function uses TotalBeachCells */
/* Uses but does not affect the global arrays:  DirectionAcrossBorder[] */
{
  int i;

  /* float AverageFull; */
  float TotalPercentInBeaches = 0.;

  for (i = 1; i <= TotalBeachCells; i++) {
    if ((DirectionAcrossBorder[i - 1] == 'r') &&
        (DirectionAcrossBorder[i] == 'r'))
      FlowThroughCell[i] = 'R'; /*Right */
    if ((DirectionAcrossBorder[i - 1] == 'r') &&
        (DirectionAcrossBorder[i] == 'l'))
      FlowThroughCell[i] = 'C'; /*Convergent */
    if ((DirectionAcrossBorder[i - 1] == 'l') &&
        (DirectionAcrossBorder[i] == 'r'))
      FlowThroughCell[i] = 'D'; /*Divergent */
    if ((DirectionAcrossBorder[i - 1] == 'l') &&
        (DirectionAcrossBorder[i] == 'l'))
      FlowThroughCell[i] = 'L'; /*Left */

    TotalPercentInBeaches +=
        (PercentFullSand[X[i]][Y[i]] + PercentFullRock[X[i]][Y[i]]);
  }

  /*
     AverageFull = TotalPercentInBeaches/TotalBeachCells;

     for (i=1; i <= TotalBeachCells; i++)
     {
     if ((PercentFullSand[X[i]][Y[i]] + PercentFullRock[X[i]][Y[i]]) >
     AverageFull)
     PutPixel(Y[i]*CellPixelSize , X[i]*CellPixelSize, 50, 200, 0);
     else
     PutPixel(Y[i]*CellPixelSize, X[i]*CellPixelSize, 200, 0, 00);
     }
   */
}

void FixFlow(void)
/* LMV */
/* Checks and distributes sediment in adjacent cells based on FlowThroughCell[]
   */
/* Sweeps left to right to check cells where DirectionAcrossBorder='r' (D to R,
   D to C, R to R, R to C) */
/* Sweeps right to left to check cells where DirectionAcrossBorder='l' (D to L,
   D to C, L to L, L to C) */
/* This function uses TotalBeachCells */
/* Uses but does not affect the global arrays: FlowThroughCell[],
   VolumeAcrossBorder[], PFS[][]         */
/* reminder--VolumeAcrossBorder[] is always on right side of cell */
/* revised 5/04 LMV */
{
  int i;
  float AmountSand; /* amount of sand availible for transport (includes amount
                       in cell behind) LMV */

  float Depth; /* Depth of current cell */
  /* float        DeltaArea;       Holds change in area for cell (m^2) */
  float Distance;   /* distance from shore to intercept of equilib. profile and
                       overall slope (m) */
  float Xintercept; /* X position of intercept of equilib. profile and overall
                       slope */
  float DepthEffective; /* depth at x intercept, Brad: where slope becomes <
                           equilib slope, */
  /*                             so sand stays piled against shoreface */

  for (i = 1; i <= TotalBeachCells - 1; i++)
    ActualVolumeAcross[i] = VolumeAcrossBorder[i];

  for (i = 1; i <= TotalBeachCells - 1;
       i++) { /*get the actual volumes across for left and right borders of D
                 cells */

    /* calculate effective depth */
    Depth = INITIAL_DEPTH + ((X[i] - INIT_BEACH) * CELL_WIDTH * SHELF_SLOPE);
    Distance = Depth / (SHOREFACE_SLOPE - SHELF_SLOPE * cos(ShorelineAngle[i]));
    Xintercept = X[i] + Distance * cos(ShorelineAngle[i]) / CELL_WIDTH;
    DepthEffective =
        INITIAL_DEPTH + ((Xintercept - INIT_BEACH) * CELL_WIDTH * SHELF_SLOPE);

    if (DepthEffective < DEPTH_SHOREFACE) DepthEffective = DEPTH_SHOREFACE;

    /*set the D's first - left and right border set */

    if (FlowThroughCell[i] == 'D') { /*Flow from a Divergent cell (to a
                                        Convergent cell or a Right cell) */

      if (AllBeach[XBehind[i]][YBehind[i]] == 'y') {
        if (PercentFullRock[X[i]][Y[i]] == 0.0) {
          AmountSand = (PercentFullSand[X[i]][Y[i]] +
                        PercentFullSand[XBehind[i]][YBehind[i]]) *
                       CELL_WIDTH * CELL_WIDTH * DepthEffective;
          if (debug13) printf("(D)X: %d, Y: %d\n", X[i], Y[i]);
        } else {
          AmountSand = (PercentFullSand[X[i]][Y[i]]) * CELL_WIDTH * CELL_WIDTH *
                       DepthEffective;
          if (debug13) printf("(RockD)X: %d, Y: %d\n", X[i], Y[i]);
        }
      }

      else { /*I don't know how much I have, just figure it out and give me some
                sand */

        /* printf("Check Amount sand X: %d, Y: %d \n", X[i], Y[i]); */
        AmountSand = (PercentFullSand[X[i]][Y[i]]) * CELL_WIDTH * CELL_WIDTH *
                     DepthEffective;
        if (debug13) printf("(Dd)X: %d, Y: %d\n", X[i], Y[i]);
      }

      if (debug13)
        printf("AmountSand[%d]: %f DepthEffective: %f\n", i, AmountSand,
               DepthEffective);

      if ((VolumeAcrossBorder[i - 1] + VolumeAcrossBorder[i]) <= AmountSand) {
        ActualVolumeAcross[i] = VolumeAcrossBorder[i];
        ActualVolumeAcross[i - 1] = VolumeAcrossBorder[i - 1];

        if (debug13)
          printf("D VolumeWantL[%d] %f, D VolumeGetL[%d] %f\n", i - 1,
                 VolumeAcrossBorder[i - 1], i - 1, ActualVolumeAcross[i - 1]);
        if (debug13)
          printf("D VolumeWantR[%d] %f, D VolumeGetR[%d] %f\n", i,
                 VolumeAcrossBorder[i], i, ActualVolumeAcross[i]);
        if (debug13)
          printf("VolumeOutL[%d] %f, VolumeOutR[%d] %f\n\n", i - 1,
                 ActualVolumeAcross[i - 1], i, ActualVolumeAcross[i]);
      }
      if ((VolumeAcrossBorder[i - 1] + VolumeAcrossBorder[i]) > AmountSand)
      /*divide up proportionally */
      {
        ActualVolumeAcross[i] =
            ((VolumeAcrossBorder[i] /
              (VolumeAcrossBorder[i - 1] + VolumeAcrossBorder[i])) *
             AmountSand);
        ActualVolumeAcross[i - 1] =
            ((VolumeAcrossBorder[i - 1] /
              (VolumeAcrossBorder[i - 1] + VolumeAcrossBorder[i])) *
             AmountSand);

        if (debug13)
          printf("   D VolumeWantL[%d] %f, D VolumeGetL[%d] %f\n", i - 1,
                 VolumeAcrossBorder[i - 1], i - 1, ActualVolumeAcross[i - 1]);
        if (debug13)
          printf("   D VolumeWantR[%d] %f, D VolumeGetR[%d] %f\n", i,
                 VolumeAcrossBorder[i], i, ActualVolumeAcross[i]);
        if (debug13)
          printf("   VolumeOutL[%d] %f, VolumeOutR[%d] %f\n\n", i - 1,
                 ActualVolumeAcross[i - 1], i, ActualVolumeAcross[i]);
      }
    }
  }

  for (i = 1; i <= TotalBeachCells - 1;
       i++) { /*get the actual volumes across right borders for R cells */

    /* calculate effective depth */
    Depth = INITIAL_DEPTH + ((X[i] - INIT_BEACH) * CELL_WIDTH * SHELF_SLOPE);
    Distance = Depth / (SHOREFACE_SLOPE - SHELF_SLOPE * cos(ShorelineAngle[i]));
    Xintercept = X[i] + Distance * cos(ShorelineAngle[i]) / CELL_WIDTH;
    DepthEffective =
        INITIAL_DEPTH + ((Xintercept - INIT_BEACH) * CELL_WIDTH * SHELF_SLOPE);

    if (DepthEffective < DEPTH_SHOREFACE) DepthEffective = DEPTH_SHOREFACE;

    if (FlowThroughCell[i] == 'R') { /*Flow from a Right cell (to a Right cell
                                        or a Convergent cell) */

      if (AllBeach[XBehind[i]][YBehind[i]] == 'y') {
        if (PercentFullRock[X[i]][Y[i]] == 0.0) {
          AmountSand = (PercentFullSand[X[i]][Y[i]] +
                        PercentFullSand[XBehind[i]][YBehind[i]]) *
                       CELL_WIDTH * CELL_WIDTH * DepthEffective;
          if (debug13) printf("(R)X: %d, Y: %d\n", X[i], Y[i]);
        } else {
          AmountSand = (PercentFullSand[X[i]][Y[i]]) * CELL_WIDTH * CELL_WIDTH *
                       DepthEffective;
          if (debug13) printf("(RockR)X: %d, Y: %d\n", X[i], Y[i]);
        }
      }

      else { /*I don't know how much I have, just figure it out and give me some
                sand */

        /* printf("Check 2 Amount sand X: %d, Y: %d \n", X[i], Y[i]); */
        AmountSand = (PercentFullSand[X[i]][Y[i]]) * CELL_WIDTH * CELL_WIDTH *
                     DepthEffective;
        if (debug13) printf("(Rd)X: %d, Y: %d\n", X[i], Y[i]);
      }

      if (debug13)
        printf("AmountSand[%d]: %f DepthEffective: %f\n", i, AmountSand,
               DepthEffective);

      if ((ActualVolumeAcross[i - 1] + AmountSand) >= VolumeAcrossBorder[i]) {
        ActualVolumeAcross[i] = VolumeAcrossBorder[i];

        if (debug13)
          printf("R VolumeWant[%d] %f, R VolumeGet[%d] %f\n", i,
                 VolumeAcrossBorder[i], i, ActualVolumeAcross[i]);
        if (debug13)
          printf("VolumeIn[%d] %f, VolumeOut[%d] %f\n\n", i,
                 ActualVolumeAcross[i - 1], i, ActualVolumeAcross[i]);
      }
      if ((ActualVolumeAcross[i - 1] + AmountSand) < VolumeAcrossBorder[i]) {
        ActualVolumeAcross[i] = ActualVolumeAcross[i - 1] + AmountSand;

        if (debug13)
          printf("   R VolumeWant[%d] %f, R VolumeGet[%d] %f\n", i,
                 VolumeAcrossBorder[i], i, ActualVolumeAcross[i]);
        if (debug13)
          printf("   VolumeIn[%d] %f, VolumeOut[%d] %f\n\n", i,
                 ActualVolumeAcross[i - 1], i, ActualVolumeAcross[i]);
      }
    }
  }

  for (i = TotalBeachCells - 1; i >= 1;
       i--) { /*get the volumes across left side of cell for L cells */

    Depth = INITIAL_DEPTH + ((X[i] - INIT_BEACH) * CELL_WIDTH * SHELF_SLOPE);
    Distance = Depth / (SHOREFACE_SLOPE - SHELF_SLOPE * cos(ShorelineAngle[i]));
    Xintercept = X[i] + Distance * cos(ShorelineAngle[i]) / CELL_WIDTH;
    DepthEffective =
        INITIAL_DEPTH + ((Xintercept - INIT_BEACH) * CELL_WIDTH * SHELF_SLOPE);

    if (DepthEffective < DEPTH_SHOREFACE) DepthEffective = DEPTH_SHOREFACE;

    if (FlowThroughCell[i] ==
        'L') { /*Flow from a Left cell (to a Convergent cell or a Left cell) */

      if (AllBeach[XBehind[i]][YBehind[i]] == 'y') {
        if (PercentFullRock[X[i]][Y[i]] == 0.0) {
          AmountSand = (PercentFullSand[X[i]][Y[i]] +
                        PercentFullSand[XBehind[i]][YBehind[i]]) *
                       CELL_WIDTH * CELL_WIDTH * DepthEffective;
          if (debug13) printf("(L)X: %d, Y: %d\n", X[i], Y[i]);
        } else {
          AmountSand = (PercentFullSand[X[i]][Y[i]]) * CELL_WIDTH * CELL_WIDTH *
                       DepthEffective;
          if (debug13) printf("(RockL)X: %d, Y: %d\n", X[i], Y[i]);
        }
      }

      else { /*I don't know how much I have, just figure it out and give me some
                sand */

        /* printf("Check 3 Amount sand X: %d, Y: %d \n", X[i], Y[i]); */
        AmountSand = (PercentFullSand[X[i]][Y[i]]) * CELL_WIDTH * CELL_WIDTH *
                     DepthEffective;
        if (debug13) printf("(Ld)X: %d, Y: %d\n", X[i], Y[i]);
      }

      if (debug13)
        printf("AmountSand[%d]: %f DepthEffective: %f\n", i, AmountSand,
               DepthEffective);

      if ((ActualVolumeAcross[i] + AmountSand) >= VolumeAcrossBorder[i - 1]) {
        ActualVolumeAcross[i - 1] = VolumeAcrossBorder[i - 1];

        if (debug13)
          printf("L VolumeWant[%d] %f, L VolumeGet[%d] %f\n", i,
                 VolumeAcrossBorder[i - 1], i, ActualVolumeAcross[i - 1]);
        if (debug13)
          printf("VolumeIn[%d] %f, VolumeOut[%d] %f\n\n", i,
                 ActualVolumeAcross[i], i, ActualVolumeAcross[i - 1]);
      }
      if ((ActualVolumeAcross[i] + AmountSand) < VolumeAcrossBorder[i - 1]) {
        ActualVolumeAcross[i - 1] = ActualVolumeAcross[i] + AmountSand;

        if (debug13)
          printf("L VolumeWant[%d] %f, L VolumeGet[%d] %f\n", i,
                 VolumeAcrossBorder[i - 1], i, ActualVolumeAcross[i - 1]);
        if (debug13)
          printf("VolumeIn[%d] %f, VolumeOut[%d] %f\n\n", i,
                 ActualVolumeAcross[i], i, ActualVolumeAcross[i - 1]);
      }
    }
  }

  if (debug13)
    for (i = 1; i <= TotalBeachCells - 1; i++) {
      printf(
          "VolumeAcrossLeft[%d]: %f, VolumeAcrossRight[%d]: %f, PFS[X%d][Y%d]: "
          "%f, PFR[X%d][Y%d]: %f,\n",
          i, ActualVolumeAcross[i - 1], i, ActualVolumeAcross[i], X[i], Y[i],
          PercentFullSand[X[i]][Y[i]], X[i], Y[i], PercentFullRock[X[i]][Y[i]]);
    }
}

void TransportSedimentSweep(void)
/*  Sweep through cells to place transported sediment */
/*  Call function AdjustShore() to move sediment. */
/*  If cell full or overempty, call OopsImFull or OopsImEmpty() */
/*  This function doesn't change any values, but the functions it calls do */
/*  Uses but doesn't change:  X[], Y[], PercentFullSand[] */
/*  sweepsign added to ensure that direction of actuating changes does not */
/*      produce unwanted artifacts (e.g. make sure symmetrical */
{
  int i, ii, sweepsign;

  if (RandZeroToOne() * 2 > 1)
    sweepsign = 1;
  else
    sweepsign = 0;

  /*if (debug7) printf("\n\n TransSedSweep  Ang %f  %d\n", WaveAngle * RAD_TO_DEG,
   * current_time_step); */

  for (i = 0; i <= TotalBeachCells - 1; i++) {
    if (sweepsign == 1)
      ii = i;
    else
      ii = TotalBeachCells - 1 - i;

    if (debug7)
      printf("i: %d  ss: %d  X: %d  Y: %d  In: %.1f  Out: %.1f Across: %.1f\n",
             ii, sweepsign, X[i], Y[i], VolumeIn[i], VolumeOut[i],
             VolumeAcrossBorder[i]);
    if (debug0 && current_time_step == 4 && (ii == 485 || ii == 486)) {
      printf("pausing in tss before adjust (PFS[58][269] is %f)\n",
             PercentFullSand[58][269]);
      PauseRun(-1, -1, ii);
    }
    AdjustShore(ii);
    if (debug0 && current_time_step == 4 && (ii == 485 || ii == 486)) {
      printf("pausing in tss after adjust (PFS[58][269] is %f)\n",
             PercentFullSand[58][269]);
      PauseRun(-1, -1, ii);
    }
    if ((PercentFullSand[X[ii]][Y[ii]] +
         PercentFullRock[X[ii]]
                        [Y[ii]]) < -0.000001) { /* RCL changes 0.0 to 1E-6
                                                   whenever testing for
                                                   oopsimempty */
      if (debug0 && current_time_step == 4 && (ii == 485 || ii == 486)) {
        printf("pausing in tss after adjust, before oops (time 4, ii 485)\n");
        PauseRun(-1, -1, ii);
      }
      OopsImEmpty(X[ii], Y[ii]);
      /*printf("Empty called in TSS\n"); */
    }

    else if ((PercentFullSand[X[ii]][Y[ii]] + PercentFullRock[X[ii]][Y[ii]]) >
             1.0) {
      OopsImFull(X[ii], Y[ii]);
      /*printf("Full called in TSS\n"); */
    }
    if (debug0 && current_time_step == 4 && (ii == 485 || ii == 486)) {
      printf(
          "pausing in tss after adjust, oops, before erode, oops (PFS[58][269] "
          "is %f)\n",
          PercentFullSand[58][269]);
      PauseRun(-1, -1, ii);
    }
    if (INITIAL_CONDITION_TYPE != 3) { /* Don't erode the beach if there are no rocks... */
      ErodeTheBeach(ii);
    }
    if (PercentFullSand[X[ii]][Y[ii]] < -0.000001) { /* RCL */
      OopsImEmpty(X[ii], Y[ii]);
      printf("Empty called in TSS after eroding\n");
    }
    if (debug0 && current_time_step == 4 && (ii == 485 || ii == 486)) {
      printf(
          "pausing in tss after adjust, oops, before erode, oops (PFS[58][269] "
          "is %f)\n",
          PercentFullSand[58][269]);
      PauseRun(-1, -1, ii);
    }
  }
}

void AdjustShore(int i)

/*  Complete mass balance for incoming and ougoing sediment */
/*  This function will change the global data array PercentFullSand[][] */
/*  Uses but does not adjust arrays: */
/*              X[], Y[], ShorelineAngle[], ActualVolumeAcross  LMV */
/*  Uses global variables: SHELF_SLOPE, CELL_WIDTH, SHOREFACE_SLOPE, INITIAL_DEPTH
   */
{
  float Depth;      /* Depth of current cell */
  float DeltaArea;  /* Holds change in area for cell (m^2) */
  float Distance;   /* distance from shore to intercept of equilib. profile and
                       overall slope (m) */
  float Xintercept; /* X position of intercept of equilib. profile and overall
                       slope */
  float DepthEffective; /* depth at x intercept, Brad: where slope becomes <
                           equilib slope, */
  /*                             so sand stays piled against shoreface */

  /*  FIND EFFECTIVE DEPTH, Deff, HERE */
  /*  Deff = FROM SOLUTION TO: */
  /*      A*Distance^n = D(X) - OverallSlope*cos(ShorelineAngle[X][Y])*Distance
   */
  /*  Where Distance is from shore, perpendicular to ShorelineAngle; */
  /* Brad's stuff not being used right now (linear slope for now) : */
  /*      Xintercept = X + Distance*cos(ShorelineAngle), and n = 2/3. */
  /*      THEN Deff = D(Xintercept); */
  /*      FIND "A" FROM DEPTH OF 10METERS AT ABOUT 1000METERS. */
  /*      FOR NOW, USE n = 1: */

  Depth = INITIAL_DEPTH + ((X[i] - INIT_BEACH) * CELL_WIDTH * SHELF_SLOPE);

  Distance = Depth / (SHOREFACE_SLOPE - SHELF_SLOPE * cos(ShorelineAngle[i]));

  Xintercept = X[i] + Distance * cos(ShorelineAngle[i]) / CELL_WIDTH;

  DepthEffective =
      INITIAL_DEPTH + ((Xintercept - INIT_BEACH) * CELL_WIDTH * SHELF_SLOPE);

  if (debug7a)
    printf(
        "\n Depth %f  Distance %f XIntercept %f  DepthEffective %f "
        "DepthShoreFace %f",
        Depth, Distance, Xintercept, DepthEffective, DEPTH_SHOREFACE);

  if (DepthEffective < DEPTH_SHOREFACE) DepthEffective = DEPTH_SHOREFACE;
  /*LMV*/ if (FlowThroughCell[i] == 'L') {
    VolumeIn[i] = ActualVolumeAcross[i];
    VolumeOut[i] = ActualVolumeAcross[i - 1];
  }
  if (FlowThroughCell[i] == 'C') {
    VolumeIn[i] = ActualVolumeAcross[i];
    VolumeOut[i] = -ActualVolumeAcross[i - 1];
  }
  if (FlowThroughCell[i] == 'D') {
    VolumeIn[i] = -ActualVolumeAcross[i - 1];
    VolumeOut[i] = ActualVolumeAcross[i];
  }
  if (FlowThroughCell[i] == 'R') {
    VolumeIn[i] = ActualVolumeAcross[i - 1];
    VolumeOut[i] = ActualVolumeAcross[i];
  }

  DeltaArea = (VolumeIn[i] - VolumeOut[i]) / DepthEffective;
  /* if(you want to do sinks this way && this cell is sink) DeltaArea = 0;
     RCL addded sinks but won't do them like this probably */
  PercentFullSand[X[i]][Y[i]] += DeltaArea / (CELL_WIDTH * CELL_WIDTH);

  if (debug14)
    printf("%c  i: %d  In: %f  Out: %f\n", FlowThroughCell[i], i, VolumeIn[i],
           VolumeOut[i]);
  if (debug14)
    printf("DeltaArea[%d]: %f, PFS[%d]: %f, time_step: %d\n", i, DeltaArea, i,
           PercentFullSand[X[i]][Y[i]], current_time_step);
  if (debug14) printf("	X[i]: %d, Y[i]: %d\n", X[i], Y[i]);
}

void ErodeTheBeach(int i)

/*  This function decreases the amount of sand in each cell of the
   TotalBeachCell array */
/*  by a uniform amount (ErosionRatePerYear) */
/*  Changes PercentFullSand[][] */
/*  LMV */
{
  float PercentEroded;

  if (InShadow[i] == 'n')
    PercentEroded = ((ErosionRatePerYear * time_step) / 365) / CELL_WIDTH;
  else /* if (InShadow[i] == 'y') */
    PercentEroded = 0.0;

  if (PercentEroded > PercentFullSand[X[i]][Y[i]]) {
    if (PercentFullRock[X[i]][Y[i]] == 0.0) {
      if (PercentEroded < (PercentFullSand[X[i]][Y[i]] +
                           PercentFullSand[XBehind[i]][YBehind[i]])) {
        if (debug15)
          printf("	Bi: %d, PFSBefore: %f\n", i,
                 PercentFullSand[X[i]][Y[i]]);
        PercentFullSand[XBehind[i]][YBehind[i]] -=
            PercentEroded - PercentFullSand[X[i]][Y[i]];
        PercentFullSand[X[i]][Y[i]] = 0.0;
        if (debug15)
          printf("	BPercentEroded: %f, PFSAfter: %f\n", PercentEroded,
                 PercentFullSand[X[i]][Y[i]]);
      }

      else if (PercentEroded > (PercentFullSand[X[i]][Y[i]] +
                                PercentFullSand[XBehind[i]][YBehind[i]]))
      /* if PFR 0 && PercentEroded > PFS + PFSBehind... */
      {
        /*if(!(PercentFullSand[X[i]][Y[i]] <= 0.000001
           && PercentFullRock[XBehind[i]][YBehind[i]] >= 0.999999)) {
           printf("Check Erosion, something strange (%d,%d)\n",X[i],Y[i]); */
        /*debug15 = 1;
           PrintLocalConds(X[i],Y[i]); */
        /*} */

        if (debug15) printf("\nErode it all!\n");
        if (debug15)
          printf(
              "cell %d (%d,%d) PFSBefore: %2.20f with Behind cell (%d,%d) "
              "PFSBefore: %2.20f\n",
              i, X[i], Y[i], PercentFullSand[X[i]][Y[i]], XBehind[i],
              YBehind[i], PercentFullSand[XBehind[i]][YBehind[i]]);
        PercentFullSand[XBehind[i]][YBehind[i]] = 0.0;
        PercentFullSand[X[i]][Y[i]] = 0.0;
        if (debug15)
          printf("PercentEroded: %f (PFS for cell and behind cell set to 0)\n",
                 PercentEroded);
        /*if(debug15) PauseRun(-1,-1,i);
           debug15 = 0; */
      }
    } else {
      PercentFullSand[X[i]][Y[i]] = 0.0;
      if (debug15)
        printf("BPercentEroded: %f, PFSAfter: %f\n", PercentEroded,
               PercentFullSand[X[i]][Y[i]]);
    }
  }

  else { /*normal case, enough sand */

    if (debug15)
      printf("i: %d, PFSBefore: %f\n", i, PercentFullSand[X[i]][Y[i]]);
    PercentFullSand[X[i]][Y[i]] -= PercentEroded;
    if (debug15)
      printf("PercentEroded: %f, PFSAfter: %f\n", PercentEroded,
             PercentFullSand[X[i]][Y[i]]);
  }
}

void OopsImEmpty(int x, int y)

/*  If a cell is under-full, this will find source for desparity and move brach
   in      */
/*  Function completly changed 5/21/02 sandrevt.c */
/*              New Approach - steal from all neighboring AllBeach cells */
/*              Backup plan - steal from all neighboring percent full > 0 */
/*  Function adjusts primary data arrays: */
/*              AllBeach[][] and PercentFullSand[][] */
/* LMV Needs to be adjusted -- can't take sand from rock cells, uses AllRock[][]
   */
{
  int emptycells = 0;
  int emptycells2 = 0;
  /* && ((x == 59 && y == 268) || (x == 58 && y == 269) || (x == 59 && y ==
  269))
  if( debug0 && x == 59 && y == 269 ) {
  debug8 = 1;
  } */
  if (debug8)
    printf("\n	OOPS I'm EMPTY!  X: %d  Y: %d PFS: %f PFR: %f", x, y,
           PercentFullSand[x][y], PercentFullRock[x][y]);

  /* find out how many AllBeaches to take from */

  if ((AllBeach[x - 1][y] == 'y') && (PercentFullRock[x - 1][y] == 0.0))
    /*LMV*/ emptycells += 1;
  if ((AllBeach[x + 1][y] == 'y') && (PercentFullRock[x + 1][y] == 0.0))
    emptycells += 1;
  if ((AllBeach[x][y - 1] == 'y') && (PercentFullRock[x][y - 1] == 0.0))
    emptycells += 1;
  if ((AllBeach[x][y + 1] == 'y') && (PercentFullRock[x][y + 1] == 0.0))
    emptycells += 1;
  if (debug8)
    printf("\n(counted %d AllBeaches with no rock to take from)\n", emptycells);

  if (emptycells > 0) {
    /* Now Move Sediment */

    if ((AllBeach[x - 1][y] == 'y') && (PercentFullRock[x - 1][y] == 0.0))
    /*LMV*/ {
      PercentFullSand[x - 1][y] += PercentFullSand[x][y] / emptycells;
      AllBeach[x - 1][y] = 'n';
      if (debug8) printf("  E MOVEDBACK\n");
      if (debug8) printf("  PFR behind: %f\n", PercentFullRock[x - 1][y]);

      if (debug8a && PercentFullSand[x - 1][y] < 0.0)
        printf("PFS went Negative in OopsI'mEmpty x-1: %d y:%d %f \n", x - 1, y,
               PercentFullSand[x - 1][y]);

      if (debug8a && PercentFullSand[x - 1][y] > 1.0)
        printf("PFS went >1 in OopsI'mEmpty x-1: %d y:%d %f \n", x - 1, y,
               PercentFullSand[x - 1][y]);
    }

    if ((AllBeach[x + 1][y] == 'y') && (PercentFullRock[x + 1][y] == 0.0)) {
      PercentFullSand[x + 1][y] += PercentFullSand[x][y] / emptycells;
      AllBeach[x + 1][y] = 'n';
      if (debug8) printf("  E MOVEDUP\n");
      if (debug8) printf("  PFR up: %f\n", PercentFullRock[x + 1][y]);

      if (debug8a && PercentFullSand[x + 1][y] < 0.0)
        printf("PFS went Negative in OopsI'mEmpty x+1: %d y:%d %f \n", x + 1, y,
               PercentFullSand[x + 1][y]);

      if (debug8a && PercentFullSand[x + 1][y] > 1.0)
        printf("PFS went >1 in OopsI'mEmpty x+1: %d y:%d %f \n", x + 1, y,
               PercentFullSand[x + 1][y]);
    }

    if ((AllBeach[x][y - 1] == 'y') && (PercentFullRock[x][y - 1] == 0.0)) {
      PercentFullSand[x][y - 1] += PercentFullSand[x][y] / emptycells;
      AllBeach[x][y - 1] = 'n';
      if (debug8) printf("  E MOVEDLEFT\n");
      if (debug8) printf("  PFR left: %f\n", PercentFullRock[x][y - 1]);
      /*if (debug8) PauseRun(x,y,-1); */

      if (debug8a && PercentFullSand[x][y - 1] < 0.0)
        printf("PFS went Negative in OopsI'mEmpty x: %d y-1:%d %f \n", x, y - 1,
               PercentFullSand[x][y - 1]);

      if (debug8a && PercentFullSand[x][y - 1] > 1.0)
        printf("PFS went >1 in OopsI'mEmpty x: %d y-1:%d %f \n", x, y - 1,
               PercentFullSand[x][y - 1]);
    }

    if ((AllBeach[x][y + 1] == 'y') && (PercentFullRock[x][y + 1] == 0.0)) {
      PercentFullSand[x][y + 1] += PercentFullSand[x][y] / emptycells;
      AllBeach[x][y + 1] = 'n';
      if (debug8) printf("  E MOVEDRIGHT\n");
      if (debug8) printf("  PFR right: %f\n", PercentFullRock[x][y + 1]);
      /*if (debug8) PauseRun(x,y,-1); */

      if (debug8a && PercentFullSand[x][y + 1] < 0.0)
        printf("PFS went Negative in OopsI'mEmpty x: %d y+1:%d %f \n", x, y + 1,
               PercentFullSand[x][y + 1]);

      if (debug8a && PercentFullSand[x][y - 1] > 1.0)
        printf("PFS went >1 in OopsI'mEmpty x: %d y-1:%d %f \n", x, y - 1,
               PercentFullSand[x][y - 1]);
    }
  }

  else {
    /* No full neighbors, so take away from partially full neighbors */

    if (PercentFullSand[x - 1][y] > 0.0) emptycells2 += 1;
    if (PercentFullSand[x + 1][y] > 0.0) emptycells2 += 1;
    if (PercentFullSand[x][y - 1] > 0.0) emptycells2 += 1;
    if (PercentFullSand[x][y + 1] > 0.0) emptycells2 += 1;
    if (debug8) printf("no full neighbors, counted %d partials\n", emptycells2);

    if (emptycells2 > 0) {
      if (PercentFullSand[x - 1][y] > 0.0) {
        PercentFullSand[x - 1][y] += PercentFullSand[x][y] / emptycells2;
        if (debug8) printf("  E NOTFULL MOVEDBACK\n");
        if (debug8) printf("  PFR behind: %f\n", PercentFullRock[x - 1][y]);

        if (debug8a && PercentFullSand[x - 1][y] < 0.0)
          printf("PFS went Negative in OopsI'mEmpty x-1: %d y:%d %f \n", x - 1,
                 y, PercentFullSand[x - 1][y]);

        if (debug8a && PercentFullSand[x - 1][y] > 1.0)
          printf("PFS went >1 in OopsI'mEmpty x-1: %d y:%d %f \n", x - 1, y,
                 PercentFullSand[x - 1][y]);
      }

      if (PercentFullSand[x + 1][y] > 0.0) {
        PercentFullSand[x + 1][y] += PercentFullSand[x][y] / emptycells2;
        if (debug8) printf("  E NOTFULL MOVEDUP\n");
        if (debug8) printf("  PFR up: %f\n", PercentFullRock[x + 1][y]);

        if (debug8a && PercentFullSand[x + 1][y] < 0.0)
          printf("PFS went Negative in OopsI'mEmpty x+1: %d y:%d %f \n", x + 1,
                 y, PercentFullSand[x + 1][y]);

        if (debug8a && PercentFullSand[x + 1][y] > 1.0)
          printf("PFS went >1 in OopsI'mEmpty x+1: %d y:%d %f \n", x + 1, y,
                 PercentFullSand[x + 1][y]);
      }

      if (PercentFullSand[x][y - 1] > 0.0) {
        PercentFullSand[x][y - 1] += PercentFullSand[x][y] / emptycells2;
        if (debug8) printf("  E NOTFULL MOVEDLEFT\n");
        if (debug8) printf("  PFR left: %f\n", PercentFullRock[x][y - 1]);
        /*if (debug8) PauseRun(x,y,-1); */

        if (debug8a && PercentFullSand[x][y - 1] < 0.0)
          printf(
              "PFS went Negative in OopsI'mEmpty...Bummer x: %d y-1:%d %f \n",
              x, y - 1, PercentFullSand[x][y - 1]);

        if (debug8a && PercentFullSand[x][y - 1] > 1.0)
          printf("PFS went >1 in OopsI'mEmpty...Bummer x: %d y-1:%d %f \n", x,
                 y - 1, PercentFullSand[x][y - 1]);
      }

      if (PercentFullSand[x][y + 1] > 0.0) {
        PercentFullSand[x][y + 1] += PercentFullSand[x][y] / emptycells2;
        if (debug8) printf(" E NOTFULL MOVEDRIGHT\n");
        if (debug8) printf("  PFR right: %f\n", PercentFullRock[x][y + 1]);
        /*if (debug8) PauseRun(x,y,-1); */

        if (debug8a && PercentFullSand[x][y + 1] < 0.0)
          printf(
              "PFS went Negative in OopsI'mEmpty...Bummer x: %d y+1:%d %f \n",
              x, y + 1, PercentFullSand[x][y + 1]);
        if (debug8a && PercentFullSand[x][y - 1] > 1.0)
          printf("PFS went >1 in OopsI'mEmpty...Bummer x: %d y-1:%d %f \n", x,
                 y - 1, PercentFullSand[x][y - 1]);
      }
    } else {
      printf("@@@ Didn't find anywhere to steal sand from!! x: %d  y: %d\n", x,
             y);
      /* PauseRun(x,y,-1); */
    }
  }

  if ((PercentFullSand[x][y]) + (PercentFullRock[x][y]) < 1.0)
    AllBeach[x][y] = 'n';
  /*LMV*/ PercentFullSand[x][y] = 0.0;

  if (debug8) printf("\n");

  if (debug0 && PercentFullSand[58][269] < 0.0)
    printf("^^^^^^(58,259) made underfull (%f) after calling oops(%d,%d)^^^^\n",
           PercentFullSand[58][269], x, y);
}

void OopsImFull(int x, int y)

/*  If a cell is overfull, push beach out in new direction */
/*  Completely revised 5/20/02 sandrevt.c to resolve 0% full problems, etc. */
/*  New approach:       put sand wherever 0% full in adjacent cells */
/*                      if not 0% full, then fill all non-allbeach */
/*  Function adjusts primary data arrays: */
/*              AllBeach[][] and PercentFullSand[][] */
/*  LMV add rock adjustments */
{
  int fillcells = 0;
  int fillcells2 = 0;

  if (debug8)
    printf("\n	OOPS I'M FULL: X: %d  Y: %d PFS: %f PFR: %f", x, y,
           PercentFullSand[x][y], PercentFullRock[x][y]);
  /*if (debug8) PrintLocalConds(x,y,-1); */

  /* find out how many cells will be filled up    */

  if ((PercentFullSand[x - 1][y] + PercentFullRock[x - 1][y]) == 0.0)
    fillcells += 1;
  if ((PercentFullSand[x + 1][y] + PercentFullRock[x + 1][y]) == 0.0)
    fillcells += 1;
  if ((PercentFullSand[x][y - 1] + PercentFullRock[x][y - 1]) == 0.0)
    fillcells += 1;
  if ((PercentFullSand[x][y + 1] + PercentFullRock[x][y + 1]) == 0.0)
    fillcells += 1;

  if (fillcells != 0) {
    /* Now Move Sediment */

    if ((PercentFullSand[x - 1][y] + PercentFullRock[x - 1][y]) == 0.0)
    /*LMV*/ {
      PercentFullSand[x - 1][y] +=
          ((PercentFullSand[x][y] + PercentFullRock[x][y]) - 1.0) / fillcells;
      if (debug80) printf("  F MOVEDBACK\n");
      if (debug8) printf("  PFR behind: %f\n", PercentFullRock[x - 1][y]);

      if (debug80 && PercentFullSand[x - 1][y] < 0.0)
        printf("1 PFS went Negative in OopsI'mFull...Bummer x-1: %d y:%d %f \n",
               x - 1, y, PercentFullSand[x - 1][y]);

      if (debug80 && PercentFullSand[x - 1][y] > 1.0)
        printf("1 PFS went >1 in OopsI'mFull...Bummer x-1: %d y:%d %f \n",
               x - 1, y, PercentFullSand[x - 1][y]);
    }

    if ((PercentFullSand[x + 1][y] + PercentFullRock[x + 1][y]) == 0.0) {
      PercentFullSand[x + 1][y] +=
          ((PercentFullSand[x][y] + PercentFullRock[x][y]) - 1.0) / fillcells;
      if (debug80) printf("  F MOVEDUP\n");
      if (debug8) printf("  PFR up: %f\n", PercentFullRock[x + 1][y]);

      if (debug80 && PercentFullSand[x + 1][y] < 0.0)
        printf("1 PFS went Negative in OopsI'mFull...Bummer x+1: %d y:%d %f \n",
               x + 1, y, PercentFullSand[x + 1][y]);

      if (debug80 && PercentFullSand[x + 1][y] > 1.0)
        printf("1 PFS went >1 in OopsI'mFull...Bummer x+1: %d y:%d %f \n",
               x + 1, y, PercentFullSand[x + 1][y]);
    }

    if ((PercentFullSand[x][y - 1] + PercentFullRock[x][y - 1]) == 0.0) {
      PercentFullSand[x][y - 1] +=
          ((PercentFullSand[x][y] + PercentFullRock[x][y]) - 1.0) / fillcells;
      if (debug80) printf("  F MOVEDLEFT\n");
      if (debug8) printf("  PFR left: %f\n", PercentFullRock[x][y - 1]);
      /*if (debug8) PauseRun(x,y,-1); */

      if (debug80 && PercentFullSand[x][y - 1] < 0.0)
        printf("1 PFS went Negative in OopsI'mFull...Bummer x: %d y-1:%d %f \n",
               x, y - 1, PercentFullSand[x][y - 1]);

      if (debug80 && PercentFullSand[x][y - 1] > 1.0)
        printf("1 PFS went >1 in OopsI'mFull...Bummer x: %d y-1:%d %f \n", x,
               y - 1, PercentFullSand[x][y - 1]);
    }

    if ((PercentFullSand[x][y + 1] + PercentFullRock[x][y + 1]) == 0.0) {
      PercentFullSand[x][y + 1] +=
          ((PercentFullSand[x][y] + PercentFullRock[x][y]) - 1.0) / fillcells;
      if (debug80) printf("  F MOVEDRIGHT\n");
      if (debug8) printf("  PFR right: %f\n", PercentFullRock[x][y + 1]);
      /*if (debug8) PauseRun(x,y,-1); */

      if (debug80 && PercentFullSand[x][y + 1] < 0.0)
        printf("1 PFS went Negative in OopsI'mFull...Bummer x: %d y+1:%d %f \n",
               x, y + 1, PercentFullSand[x][y + 1]);

      if (debug80 && PercentFullSand[x][y + 1] > 1.0)
        printf("1 PFS went >1 in OopsI'mFull...Bummer x: %d y+1:%d %f \n", x,
               y + 1, PercentFullSand[x][y + 1]);
    }
  } else {
    /* No fully empty neighbors, so distribute to partially full neighbors */

    if ((PercentFullSand[x - 1][y] + PercentFullRock[x - 1][y]) < 1.0)
      fillcells2 += 1;
    if ((PercentFullSand[x + 1][y] + PercentFullRock[x + 1][y]) < 1.0)
      fillcells2 += 1;
    if ((PercentFullSand[x][y - 1] + PercentFullRock[x][y - 1]) < 1.0)
      fillcells2 += 1;
    if ((PercentFullSand[x][y + 1] + PercentFullRock[x][y + 1]) < 1.0)
      fillcells2 += 1;

    if (fillcells2 > 0) {
      if ((PercentFullSand[x - 1][y] + PercentFullRock[x - 1][y]) < 1.0) {
        PercentFullSand[x - 1][y] +=
            ((PercentFullSand[x][y] + PercentFullRock[x][y]) - 1.0) /
            fillcells2;
        if (debug80) printf("  Fa MOVEDBACK\n");
        if (debug8) printf("  PFR behind: %f\n", PercentFullRock[x - 1][y]);

        if (debug80 && PercentFullSand[x - 1][y] < 0.0)
          printf("PFS went Negative in OopsI'mFull...Bummer x-1: %d y:%d %f \n",
                 x - 1, y, PercentFullSand[x - 1][y]);

        if (debug80 && PercentFullSand[x - 1][y] > 1.0)
          printf("PFS went >1 in OopsI'mFull...Bummer x-1: %d y:%d %f \n",
                 x - 1, y, PercentFullSand[x - 1][y]);
      }

      if ((PercentFullSand[x + 1][y] + PercentFullRock[x + 1][y]) < 1.0) {
        PercentFullSand[x + 1][y] +=
            ((PercentFullSand[x][y] + PercentFullRock[x][y]) - 1.0) /
            fillcells2;
        if (debug80) printf("  Fa MOVEDUP\n");
        if (debug8) printf("  PFR up: %f\n", PercentFullRock[x + 1][y]);

        if (debug80 && PercentFullSand[x + 1][y] < 0.0)
          printf("PFS went Negative in OopsI'mFull...Bummer x+1: %d y:%d %f \n",
                 x + 1, y, PercentFullSand[x + 1][y]);

        if (debug80 && PercentFullSand[x + 1][y] > 1.0)
          printf("PFS went >1 in OopsI'mFull...Bummer x+1: %d y:%d %f \n",
                 x + 1, y, PercentFullSand[x + 1][y]);
      }

      if ((PercentFullSand[x][y - 1] + PercentFullRock[x][y - 1]) < 1.0) {
        PercentFullSand[x][y - 1] +=
            ((PercentFullSand[x][y] + PercentFullRock[x][y]) - 1.0) /
            fillcells2;
        if (debug80) printf("  Fa MOVEDLEFT\n");
        if (debug8) printf("  PFR left: %f\n", PercentFullRock[x][y - 1]);

        if (debug80 && PercentFullSand[x][y - 1] < 0.0)
          printf("PFS went Negative in OopsI'mFull...Bummer x: %d y-1:%d %f \n",
                 x, y - 1, PercentFullSand[x][y - 1]);

        if (debug80 && PercentFullSand[x][y - 1] > 1.0)
          printf("PFS went >1 in OopsI'mFull...Bummer x: %d y-1:%d %f \n", x,
                 y - 1, PercentFullSand[x][y - 1]);
      }

      if ((PercentFullSand[x][y + 1] + PercentFullRock[x][y + 1]) < 1.0) {
        PercentFullSand[x][y + 1] +=
            ((PercentFullSand[x][y] + PercentFullRock[x][y]) - 1.0) /
            fillcells2;
        if (debug80) printf("  Fa MOVEDRIGHT\n");
        if (debug8) printf("  PFR right: %f\n", PercentFullRock[x][y + 1]);

        if (debug80 && PercentFullSand[x][y + 1] < 0.0)
          printf("PFS went Negative in OopsI'mFull...Bummer x: %d y+1:%d %f \n",
                 x, y + 1, PercentFullSand[x][y + 1]);

        if (debug80 && PercentFullSand[x][y + 1] > 1.0)
          printf("PFS went >1 in OopsI'mFull...Bummer x: %d y+1:%d %f \n", x,
                 y + 1, PercentFullSand[x][y + 1]);
      }
    } else {
      if (debug8) {
        printf("Nobody wants our sand!!! x: %d  y: %d  PFS: %f  PFR: %f\n", x,
               y, PercentFullSand[x][y], PercentFullRock[x][y]);
        printf("x+1: %f, y+1: %f, x-1: %f, y-1: %f\n",
               (PercentFullSand[x + 1][y] + PercentFullRock[x + 1][y]),
               (PercentFullSand[x][y + 1] + PercentFullRock[x][y + 1]),
               (PercentFullSand[x - 1][y] + PercentFullRock[x][y + 1]),
               (PercentFullSand[x][y - 1] + PercentFullRock[x][y - 1]));
        /* PauseRun(x,y,-1); */
        printf("(excess sand will be removed; sand = 1 - rock)\n");
      }
    }
  }

  if (PercentFullSand[x][y] + PercentFullRock[x][y] >= 1.0)
    AllBeach[x][y] = 'y';
  /*LMV*/ PercentFullSand[x][y] = 1.0 - PercentFullRock[x][y];

  if (debug80 && PercentFullSand[x][y] < 0.0)
    printf("	PFS went Negative in OopsI'mFull...Bummer x: %d y:%d %f \n", x,
           y, PercentFullSand[x][y]);

  if (debug80 && PercentFullSand[x][y] > 1.0)
    printf("	PFS went >1 in OopsI'mFull...Bummer x: %d y:%d %f \n", x, y,
           PercentFullSand[x][y]);
}

void FixBeach(void)

/* Hopefully addresses strange problems caused by filling/emptying of cells */
/* Looks at entire data set */
/* Find unattached pieces of sand and moves them back to the shore */
/* Takes care of 'floating bits' of sand */
/* Also takes care of over/under filled beach pieces */
/* Revised 5/21/02 to move sand to all adjacent neighbors sandrevt.c */
/* Changes global variable PercentFullSand[][] */
/* Uses and can change AllBeach[][] */
/* sandrevx.c - added sweepsign to reduce chances of asymmetrical artifacts */
/* LMV rock fix - do not fix bare bedrock cells */
{
  int i, x, y, sweepsign, done, counter;
  int fillcells3 = 0;
  int Corner = 0;
  /*LMV*/ int xstart;

  /*if (debug9) printf("\n\nFIXBEACH      %d     %f\n", current_time_step,
   * WaveAngle*RAD_TO_DEG); */

  if (RandZeroToOne() * 2 > 1)
    sweepsign = 1;
  else
    sweepsign = 0;

  for (x = 1; x < ShadowXMax - 1; x++) {	
    for (i = 1; i < 2 * Y_MAX - 1;
         i++) { /* RCL: changing loops to ignore border */		
      if (sweepsign == 1)
        y = i;
      else
        y = 2 * Y_MAX - i;

      /*Take care of corner problem?
         if  (((AllBeach[x][y] == 'n') && (PercentFullRock[x][y] > 0.0)) &&
         (((AllBeach[x][y-1] == 'n') && (AllBeach[x-1][y] == 'n') &&
         (AllRock[x-1][y-1] == 'y'))
         ||((AllBeach[x-1][y] == 'n') && (AllBeach[x][y+1] == 'n') &&
         (AllRock[x-1][y+1] == 'y'))))

         Corner = 1;
         else */
      Corner = 0;

      /* Take care of situations that shouldn't exist */

      if ((PercentFullSand[x][y] + PercentFullRock[x][y]) < -0.000001) {
        /* RCL changed  0.0 to 10E-6--floating-pt equality is messy */
        /* printf("Too empty\n"); */
        if (debug9 && y != 0)
          printf("\nUnder 0 Percent X: %d  Y: %d PFS: %f PFR: %f\n", x, y,
                 PercentFullSand[x][y], PercentFullRock[x][y]);
        if (debug0) {
          printf("pausing in fixbeach before before it calls oops\n");
          PauseRun(x, y, -1);
        }
        AllBeach[x][y] = 'n';
        OopsImEmpty(x, y);
      }

      if ((PercentFullSand[x][y] + PercentFullRock[x][y]) > 1.0) {
        /* printf("Too full\n"); */
        AllBeach[x][y] = 'y';
        if (debug9 && y != 0)
          printf("\nOver 100 Percent X: %d  Y: %d PFS: %f PFR: %f\n", x, y,
                 PercentFullSand[x][y], PercentFullRock[x][y]);
        OopsImFull(x, y);
      }

      if ((((PercentFullSand[x][y] + PercentFullRock[x][y]) >= 0.0) &&
           ((PercentFullSand[x][y] + PercentFullRock[x][y]) < 1.0)) &&
          (AllBeach[x][y] == 'y')) {
        AllBeach[x][y] = 'n';
        if (debug9 && y != 0) {
          printf("\nALLBeachProb X: %d  Y: %d\n", x, y);
          printf("PFS: %f, PFR: %f\n", PercentFullSand[x][y],
                 PercentFullRock[x][y]);
        }
      }

      /* Take care of 'loose' bits of sand */

      if (Corner == 0) {
        fillcells3 = 0;

        if (((PercentFullSand[x][y] + PercentFullRock[x][y]) != 0.0) &&
            ((PercentFullSand[x - 1][y] + PercentFullRock[x - 1][y]) < 1.0) &&
            ((PercentFullSand[x + 1][y] + PercentFullRock[x + 1][y]) < 1.0) &&
            ((PercentFullSand[x][y + 1] + PercentFullRock[x][y + 1]) < 1.0) &&
            ((PercentFullSand[x][y - 1] + PercentFullRock[x][y - 1]) < 1.0) &&
            (AllBeach[x][y] == 'n'))

        /* Beach in cell, but bottom, top, right, and left neighbors not all
           full */
        {
          if (debug9 && y != 0)
            printf(
                "\nFB Moved loose bit of sand,  X: %d  Y: %d  PFS: %f PFR: %f ",
                x, y, PercentFullSand[x][y], PercentFullRock[x][y]);

          /* distribute to partially full neighbors */

          if (((PercentFullSand[x - 1][y] + PercentFullRock[x - 1][y]) < 1.0) &&
              ((PercentFullSand[x - 1][y] + PercentFullRock[x - 1][y]) > 0.0))
            fillcells3 += 1;
          if (((PercentFullSand[x + 1][y] + PercentFullRock[x + 1][y]) < 1.0) &&
              ((PercentFullSand[x + 1][y] + PercentFullRock[x + 1][y]) > 0.0))
            fillcells3 += 1;
          if (((PercentFullSand[x][y + 1] + PercentFullRock[x][y + 1]) < 1.0) &&
              ((PercentFullSand[x][y + 1] + PercentFullRock[x][y + 1]) > 0.0))
            fillcells3 += 1;
          if (((PercentFullSand[x][y - 1] + PercentFullRock[x][y - 1]) < 1.0) &&
              ((PercentFullSand[x][y - 1] + PercentFullRock[x][y - 1]) > 0.0))
            fillcells3 += 1;

          if ((fillcells3 > 0) && (fillcells3 <= 4)) {
            if (((PercentFullSand[x - 1][y] + PercentFullRock[x - 1][y]) <
                 1.0) &&
                ((PercentFullSand[x - 1][y] + PercentFullRock[x - 1][y]) >
                 0.0)) {
              /* PercentFullSand[x-1][y] +=
               * (PercentFullSand[x][y]-PercentFullRock[x][y])/fillcells3; */
              PercentFullSand[x - 1][y] +=
                  PercentFullSand[x][y] / fillcells3; /* RCL */
              if (debug9) printf("  FB MOVEDBACK\n");

              if (debug9a && PercentFullSand[x - 1][y] < 0.0)
                printf(
                    "PFS went Negative in FixBeach...Bummer x-1: %d y:%d %f \n",
                    x - 1, y, PercentFullSand[x - 1][y]);

              if (debug9a && PercentFullSand[x - 1][y] > 1.0)
                printf("PFS went >1 in FixBeach...Bummer x-1: %d y:%d %f \n",
                       x - 1, y, PercentFullSand[x - 1][y]);
            }

            if (((PercentFullSand[x + 1][y] + PercentFullRock[x + 1][y]) <
                 1.0) &&
                ((PercentFullSand[x + 1][y] + PercentFullRock[x + 1][y]) >
                 0.0)) {
              /* PercentFullSand[x+1][y] +=
               * (PercentFullSand[x][y]-PercentFullRock[x][y])/fillcells3; */
              PercentFullSand[x + 1][y] +=
                  PercentFullSand[x][y] / fillcells3; /* RCL */
              if (debug9) printf("  FB MOVEDUP\n");

              if (debug9a && PercentFullSand[x + 1][y] < 0.0)
                printf(
                    "PFS went Negative in FixBeach...Bummer x+1: %d y:%d %f \n",
                    x + 1, y, PercentFullSand[x + 1][y]);

              if (debug9a && PercentFullSand[x + 1][y] > 1.0)
                printf("!!!PFS went >1 in FixBeach...Bummer x+1: %d y:%d %f \n",
                       x + 1, y, PercentFullSand[x + 1][y]);
            }

            if (((PercentFullSand[x][y - 1] + PercentFullRock[x][y - 1]) <
                 1.0) &&
                ((PercentFullSand[x][y - 1] + PercentFullRock[x][y - 1]) >
                 0.0)) {
              /* PercentFullSand[x][y-1] +=
               * (PercentFullSand[x][y]-PercentFullRock[x][y])/fillcells3; */
              PercentFullSand[x][y - 1] +=
                  PercentFullSand[x][y] / fillcells3; /* RCL */
              if (debug9) printf("  FB MOVEDLEFT\n");
              /*if (debug9) PauseRun(x,y,-1); */

              if (debug9a && PercentFullSand[x][y - 1] < 0.0)
                printf(
                    "PFS went Negative in FixBeach...Bummer x+1: %d y:%d %f \n",
                    x + 1, y, PercentFullSand[x][y - 1]);

              if (debug9a && PercentFullSand[x - 1][y] > 1.0)
                printf("PFS went >1 in FixBeach...Bummer x+1: %d y:%d %f \n",
                       x + 1, y, PercentFullSand[x + 1][y]);
            }

            if (((PercentFullSand[x][y + 1] + PercentFullRock[x][y + 1]) <
                 1.0) &&
                ((PercentFullSand[x][y + 1] + PercentFullRock[x][y + 1]) >
                 0.0)) {
              /* PercentFullSand[x][y+1] +=
               * (PercentFullSand[x][y]-PercentFullRock[x][y])/fillcells3; */
              PercentFullSand[x][y + 1] +=
                  PercentFullSand[x][y] / fillcells3; /* RCL */
              if (debug9) printf("  FB MOVEDRIGHT\n");
              /*if (debug9) PauseRun(x,y,-1); */

              if (debug9a && PercentFullSand[x][y + 1] < 0.0)
                printf(
                    "PFS went Negative in FixBeach...Bummer x: %d y+1:%d %f \n",
                    x, y + 1, PercentFullSand[x][y + 1]);

              if (debug9a && PercentFullSand[x][y + 1] > 1.0)
                printf("PFS went >1 in FixBeach...Bummer x: %d y+1:%d %f \n", x,
                       y + 1, PercentFullSand[x][y + 1]);
            }

		  }
		  else {
			  if (debug9) 
			  printf("Complete fixbeach breakdown x: %d  y: %d\n", x, y);
			  /*if (debug9) PauseRun(x,y,-1); */

			  xstart = X_MAX - 1;
			  while (AllBeach[xstart][y] == 'n') {
				  xstart -= 1;
			  }
			  xstart += 1;

			  /* printf("Moving random sand from X: %d, Y: %d, to X: %d, Y: %d\n",
				 x, y, xstart, y);
				 printf("moving this much: %f\n", PercentFullSand[x][y]);
				 printf("PFS Before move: %f\n", PercentFullSand[xstart][y]); */

			  PercentFullSand[xstart][y] += PercentFullSand[x][y];
			  PercentFullSand[x][y] = 0.0;
			  AllBeach[x][y] = 'n';

			  // printf("PFS After move: %f \n", PercentFullSand[xstart][y]);

			  if (INIT_ROCK > 0)
			  {
				  xstart = X_MAX - 1;
				  while (AllRock[xstart][y] == 'n') {
					  xstart -= 1;
				  }
				  xstart += 1;

				  printf("Moving bit of rock from X: %d, Y: %d, to X: %d, Y: %d\n", x,
					  y, xstart, y);
				  printf("moving this much: %f\n", PercentFullRock[x][y]);
				  printf("PFR of receiving cell before move: %f\n",
					  PercentFullRock[xstart][y]);

				  PercentFullRock[xstart][y] += PercentFullRock[x][y];
				  PercentFullRock[x][y] = 0.0;
				  AllRock[x][y] = 'n';
				  printf("PFR After move: %f \n", PercentFullRock[xstart][y]);
				  if (PercentFullRock[xstart][y] >= 1.0) {
					  AllRock[xstart][y] = 'y';
					  PercentFullRock[xstart + 1][y] +=
						  (PercentFullRock[xstart][y] - 1.0);
					  PercentFullRock[xstart][y] = 1.0;
					  if (PercentFullRock[xstart + 1][y] >= 1.0)
						  printf("rock spilled over\n");
					  printf("after second move PFR = %f ", PercentFullRock[xstart][y]);
					  printf("and above cell PFR = %f \n",
						  PercentFullRock[xstart + 1][y]);
				  }
			  }
			  PercentFullSand[x][y] = 0.0;
			  AllBeach[x][y] = 'n';

			  if (debug9) printf("\n");
		  }

          /* If we have overfilled any of the cells in this loop, need to
           * OopsImFull() */

          if ((PercentFullSand[x - 1][y] + PercentFullRock[x - 1][y]) > 1.0) {
            /* printf("     Below Overfilled\n");
               printf("PF[x-1][y]: %f", (PercentFullSand[x-1][y] +
               PercentFullRock[x-1][y]));
               printf("x-1:%d, y:%d\n", (x-1), y); */
            OopsImFull(x - 1, y);
            if (debug9) printf("	Below Overfilled\n");
          }
          if ((PercentFullSand[x][y - 1] + PercentFullRock[x][y - 1]) > 1.0) {
            /* printf("     Left Side Overfilled\n");
               printf("PF[x][y-1]: %f", (PercentFullSand[x][y-1] +
               PercentFullRock[x][y-1]));
               printf("x:%d, y-1:%d\n", x, (y-1)); */
            OopsImFull(x, y - 1);
            if (debug9) printf("	Left Side Overfilled\n");
          }
          if ((PercentFullSand[x][y + 1] + PercentFullRock[x][y + 1]) > 1.0) {
            /* printf("     Right Side Overfilled\n");
               printf("PF[x][y+1]: %f", (PercentFullSand[x][y+1] +
               PercentFullRock[x][y+1]));
               printf("x:%d, y+1:%d\n", x, (y+1)); */
            OopsImFull(x, y + 1);
            if (debug9) printf("	Right Side Overfilled\n");
          }
          if ((PercentFullSand[x + 1][y] + PercentFullRock[x + 1][y]) > 1.0) {
            /* printf("     Top Side Overfilled\n");
               printf("PF[x+1][y]: %f", (PercentFullSand[x+1][y] +
               PercentFullRock[x+1][y]));
               printf("x+1:%d, y:%d\n", (x+1), y); */
            OopsImFull(x + 1, y);
            if (debug9) printf("	Top Overfilled\n");
          }
        }
      } else if (debug9)
        printf("Corner cell, no fixbeach please PFR: %f x: %d, y: %d\n",
               PercentFullRock[x][y], x, y);
    }
  }
  /* now make sure there aren't any overfull or underfull cells */

  if (debug0) {
    printf("pausing in fixbeach before final loop\n");
    PauseRun(58, 269, -1);
  }

  if (debug9) printf("checking grid for overfull, underfull \n");
  done = FALSE;
  for (counter = 0; counter < 10 && !done; counter++) {
    done = TRUE;
    for (x = 1; x < X_MAX - 1; x++) {
      for (y = 1; y < 2 * Y_MAX - 1; y++) {
        if ((PercentFullSand[x][y] + PercentFullRock[x][y] > 1.0) ||
            (PercentFullSand[x][y] + PercentFullRock[x][y] < 0.0) ||
            PercentFullSand[x][y] < 0.0)
          done = FALSE;
        if (PercentFullSand[x][y] + PercentFullRock[x][y] > 1.0) {
          OopsImFull(x, y);
          if (debug9)
            printf("found an overfull in final fb loop (%d, %d) \n", x, y);
        }
        if (PercentFullSand[x][y] + PercentFullRock[x][y] < -0.000001) {
          OopsImEmpty(x, y);
          if (debug9)
            printf("found an underfull in final fb loop (%d, %d) \n", x, y);
        }
        if (PercentFullSand[x][y] <
            -0.000001) { /* RCL adds this, changes 0 to 1E-6 everywhere */
          OopsImEmpty(x, y);
          if (debug9)
            printf("found a sand < 0 in final fb loop (%d, %d) \n", x, y);
        }
      }
    }
  }
  if (debug9) printf("done checking grid for overfull, underfull \n");
}

float MassCount(void)

/* Counts the total volume occupied by beach cells      */
/* Uses same algorhythm as AdjustShore                  */
/* returns a float of the total sum                     */
/* Uses AllBeach[][] and PercentFullSand[][]            */
/* and INITIAL_DEPTH, CELL_WIDTH, SHELF_SLOPE              */
{
  int x, y;
  float Depth; /* Depth of current cell */
  float Mass = 0.0;

  for (x = 0; x < X_MAX; x++) {
    Depth = INITIAL_DEPTH + ((x - INIT_BEACH) * CELL_WIDTH * SHELF_SLOPE);

    for (y = 0; y < 2 * Y_MAX; y++) {
      Mass += PercentFullSand[x][y] * Depth;
    }
  }

  return Mass;
}

float Raise(float b, float e)
/**
 * function calulates b to the e power. pow has problems if b <= 0
 * 	note you can't use this at all when b is negative--it'll be
 * 	wrong in cases when the answer returned should have been
 * 	positive, as x^2 when x < 0.
 */
{
  if (b > 0)
    return powf(b, e);
  else {
    printf("Raise: can't handle negative base \n");
    printf("Wave data: %f %lf %lf\n", WaveAngle, Period, OffShoreWvHt);
    PauseRun(-1, -1, -1);
    return 0;
  }
}

// Apply periodic boundary conditions to CEM arrays.
void periodic_boundary_copy(void)
{
  const int buff_len = 2 * Y_MAX;
  int row;

  for (row=0; row < X_MAX; row++) {
    apply_periodic_boundary(AllBeach[row], sizeof(char), buff_len);
    apply_periodic_boundary(AllRock[row], sizeof(char),  buff_len);
    apply_periodic_boundary(type_of_rock[row], sizeof(char), buff_len);
    apply_periodic_boundary(PercentFullSand[row], sizeof(double), buff_len);
    apply_periodic_boundary(PercentFullRock[row], sizeof(double),buff_len);
    apply_periodic_boundary(Age[row], sizeof(int), buff_len);
  }
}


void ZeroVars(void)
/* Resets all arrays recalculated at each time step to 'zero' conditions */
{
  int z;

  for (z = 0; z < MaxBeachLength; z++) {
    X[z] = -1;
    Y[z] = -1;
    XRock[z] = -1;
    /*LMV*/ YRock[z] = -1;
    /*LMV*/ XBehind[z] = -1;
    YBehind[z] = -1;
    XRockBehind[z] = -1;
    YRockBehind[z] = -1;
    InShadow[z] = '?';
    ShorelineAngle[z] = -999;
    UpWind[z] = '?';

    VolumeAcrossBorder[z] = 0.0;
    /*LMV*/ ActualVolumeAcross[z] = 0.0;
  /*LMV*/}
}

void ReadSandFromFile(void)
/*  Reads saved output file, AllBeach[][] & PercentFullSand[][]	 */
/*  LMV - Also reads saved output file, AllRock[][], PercentFullRock[][],
   type_of_rock[][]  */
{
  int x, y;

  ReadSandFile = fopen(readfilename, "r");
  printf("CHECK READ \n");

  if (SAVE_FILE != 2) { /*line file input */
    for (y = Y_MAX / 2; y < 3 * Y_MAX / 2; y++)
      for (x = 0; x < X_MAX; x++) fscanf(ReadSandFile, " %c", &type_of_rock[x][y]);

    for (y = Y_MAX / 2; y < 3 * Y_MAX / 2; y++) {
      for (x = 0; x < X_MAX; x++) {
        fscanf(ReadSandFile, " %lf", &PercentFullRock[x][y]);

        if (PercentFullRock[x][y] >= 1.0)
          AllRock[x][y] = 'y';
        else
          AllRock[x][y] = 'n';
      }
    }

    for (y = Y_MAX / 2; y < 3 * Y_MAX / 2; y++) {
      for (x = 0; x < X_MAX; x++) {
        fscanf(ReadSandFile, " %lf", &PercentFullSand[x][y]);

        if ((PercentFullSand[x][y] + PercentFullRock[x][y]) >= 1.0)
          AllBeach[x][y] = 'y';
        else
          AllBeach[x][y] = 'n';
      }
    }

    if (SAVE_AGE)

      for (y = Y_MAX / 2; y < 3 * Y_MAX / 2; y++) {
        for (x = 0; x < X_MAX; x++) {
          fscanf(ReadSandFile, " %d", &Age[x][y]);
        }
      }
  }

  else { /*Array file input */

    for (x = (X_MAX - 1); x >= 0; x--) {
      for (y = Y_MAX / 2; y < 3 * Y_MAX / 2; y++) {
        fscanf(ReadSandFile, " %c", &type_of_rock[x][y]);
      }
      fscanf(ReadSandFile, "\n");
    }

    for (x = (X_MAX - 1); x >= 0; x--) {
      for (y = Y_MAX / 2; y < 3 * Y_MAX / 2; y++) {
        fscanf(ReadSandFile, " %lf", &PercentFullRock[x][y]);

        if (PercentFullRock[x][y] >= 1.0)
          AllRock[x][y] = 'y';
        else
          AllRock[x][y] = 'n';
      }
      fscanf(ReadSandFile, "\n");
    }

    for (x = (X_MAX - 1); x >= 0; x--) {
      for (y = Y_MAX / 2; y < 3 * Y_MAX / 2; y++) {
        fscanf(ReadSandFile, " %lf", &PercentFullSand[x][y]);

        if ((PercentFullSand[x][y] + PercentFullRock[x][y]) >= 1.0)
          AllBeach[x][y] = 'y';
        else
          AllBeach[x][y] = 'n';
      }
      fscanf(ReadSandFile, "\n");
    }

    if (SAVE_AGE) {
      for (x = (X_MAX - 1); x >= 0; x--) {
        for (y = Y_MAX / 2; y < 3 * Y_MAX / 2; y++) {
          fscanf(ReadSandFile, " %d", &Age[x][y]);
        }
        fscanf(ReadSandFile, "\n");
      }
    }
  }

  /*PrintLocalConds(5,5,-1); */
  fclose(ReadSandFile);
  printf("file read!");

  periodic_boundary_copy();
}

void SaveSandToFile(void)
/*  Saves current AllBeach[][] and PercentFullSand[][] data arrays to file
   */
/*  Save file name will add extension '.' and the current_time_step
   */
/*  LMV - Save AllRock[][], PercentFullRock[][], type_of_rock[][]
   */
{
  int x, y;
  char savename[40];

  printf("\n saving \n ");

  sprintf(savename, "%s_%d.out", savefilename, current_time_step);
  printf("Saving as: %s \n \n \n", savename);

  SaveSandFile = fopen(savename, "w");

  if (!SaveSandFile) {
    printf("problem opening output file\n");
    exit(1);
  }
  if (SAVE_FILE == 1) {
    for (y = Y_MAX / 2; y < 3 * Y_MAX / 2; y++)
      for (x = 0; x < X_MAX; x++) fprintf(SaveSandFile, " %c", type_of_rock[x][y]);
    /*LMV*/ for (y = Y_MAX / 2; y < 3 * Y_MAX / 2; y++)
      for (x = 0; x < X_MAX; x++)
        fprintf(SaveSandFile, " %lf", PercentFullRock[x][y]);
    /*LMV*/ for (y = Y_MAX / 2; y < 3 * Y_MAX / 2; y++)
      for (x = 0; x < X_MAX; x++)
        fprintf(SaveSandFile, " %lf", PercentFullSand[x][y]);

    if (SAVE_AGE)
      for (y = Y_MAX / 2; y < 3 * Y_MAX / 2; y++)
        for (x = 0; x < X_MAX; x++) fprintf(SaveSandFile, " %d", Age[x][y]);
  }

  if (SAVE_FILE == 2) {
    /* Array output */

    /* for (x=0; x<X_MAX; x++) */
    for (x = (X_MAX - 1); x >= 0; x--) {
      for (y = Y_MAX / 2; y < 3 * Y_MAX / 2; y++) {
        fprintf(SaveSandFile, " %c", type_of_rock[x][y]);
        /*LMV*/
        /* if (type_of_rock[x][y]=='f') fprintf(SaveSandFile, " 0"); */ /* Switch
                                                                         on if
                                                                         numbers
                                                                         required*/
        /* else fprintf(SaveSandFile, " 1"); */
      }
      fprintf(SaveSandFile, "\n");
    }

    for (x = (X_MAX - 1); x >= 0; x--) {
      for (y = Y_MAX / 2; y < 3 * Y_MAX / 2; y++) {
        fprintf(SaveSandFile, " %lf", PercentFullRock[x][y]);
      /*LMV*/}
      fprintf(SaveSandFile, "\n");
    }

    for (x = (X_MAX - 1); x >= 0; x--) {
      for (y = Y_MAX / 2; y < 3 * Y_MAX / 2; y++) {
        fprintf(SaveSandFile, " %lf", PercentFullSand[x][y]);
      }
      fprintf(SaveSandFile, "\n");
    }

    if (SAVE_AGE) {
      for (x = (X_MAX - 1); x >= 0; x--) {
        for (y = Y_MAX / 2; y < 3 * Y_MAX / 2; y++) {
          fprintf(SaveSandFile, " %d", Age[x][y]);
        }
        fprintf(SaveSandFile, "\n");
      }
    }
  }

  fclose(SaveSandFile);
  printf("file saved!\n");

  periodic_boundary_copy();
}

void SaveLineToFile(void)
/*  Saves data line of shoreline position rather than entire array */
/*  Main concern is to have only one data point at each alongshore location
   */
/*  Save file name will add extension '.' and the current_time_step
   */
{
  int y, x, xtop, i;
  float xsave;
  char savename[40];

  printf("\n saving \n ");

  sprintf(savename, "ShorePos_%d.dat", current_time_step / SAVE_SPACING);
  printf("Saving as: %s                 ", savename);

  SaveSandFile = fopen(savename, "w");
  if (!SaveSandFile) {
    printf("problem opening output file\n");
    exit(1);
  }

  for (y = Y_MAX / 2; y < 3 * Y_MAX / 2; y++) {
    x = X_MAX - 1;
    xtop = X_MAX;

    /* step back to where we encounter allbeach */
    while (AllBeach[x][y] == 'n') {
      x -= 1;
    }

    /* if on side of shape, need to average */
    if (PercentFullSand[x + 2][y] > 0) {
      xtop = x + 1;
      while (PercentFullSand[xtop][y] > 0) {
        xtop += 1;
      }

      xsave = x;

      for (i = x + 1; i < xtop; i++) xsave += PercentFullSand[i][y];

    }
    /* otherwise Regular Beach Condition */
    else {
      xsave = x + PercentFullSand[x + 1][y];
    }

    /* note this assumes average of beach locations should be 0.5 percentfull */
    fprintf(SaveSandFile, " %f", xsave - INIT_BEACH + 0.5);

    /*printf("y %d , xtop = %d xsave = %f \n",y,xtop,xsave);
       if (xtop != X_MAX)
       PauseRun(x+1,y,-1);                  */
  }

  fclose(SaveSandFile);
  printf("line file saved!\n\n");
}

void PauseRun(int x, int y, int in)
/* Pauses run until the 'q' key is pressed 	*/
/* Can Print or Plot Out Useful info		*/
{
  printf("\nPaused \n");

  if (SAVE_LINE) SaveLineToFile();
  if (SAVE_FILE) SaveSandToFile();

  if (NO_PAUSE_RUN) return;

  printf("\nend Pause\n");
}

void AgeCells(void)
/* Age Cells */
{
  int x, y;

  for (y = 0; y < 2 * Y_MAX; y++)
    for (x = 0; x < X_MAX; x++)
      if (PercentFullSand[x][y] == 0) {
        Age[x][y] = current_time_step % AGE_MAX;
      }
}

void ReadWaveIn(void)
/* Input Wave Distribution */
{
  int i;

  for (i = 0; i <= 36; i++) {
    WaveMax[i] = 0;
    WaveProb[i] = 0;
  }

  ReadWaveFile = fopen(readwavename, "r");
  printf("CHECK READ WAVE\n");

  fscanf(ReadWaveFile, " %d \n", &NumWaveBins);

  printf("Wave Bins %d \n", NumWaveBins);

  WaveMax[0] = -90;
  WaveProb[0] = 0;

  for (i = 1; i <= NumWaveBins; i++) {
    fscanf(ReadWaveFile, " %lf %lf", &WaveMax[i], &WaveProb[i]);
    printf("i= %d  Wave= %f Prob= %f \n", i, WaveMax[i], WaveProb[i]);
  }

  fclose(ReadWaveFile);
  printf("wave file read! \n");
}

void RealWaveIn(void)
/* Input Wave Distribution */
{
  int countout = 10;

  WaveAngleIn = -999;
  WavePeriodIn = -999;
  WaveHeightIn = -999;

  if (current_time_step == 0) {
    ReadRealWaveFile = fopen(readwavename, "r");
    printf("CHECK READ WAVE: %s\n", readwavename);
  }

  fscanf(ReadRealWaveFile, "%lf %lf %lf\n", &WaveAngleIn, &WavePeriodIn,
         &WaveHeightIn);

  while (WaveAngleIn == -999 || WavePeriodIn == -999 ||
         WaveHeightIn == -999) { /*Recycling wave data */
    countout--;
    printf(
        "Recycling wave data.\n If this happens every timestep then there "
        "maybe a problem reading the wave data?\n");
    fclose(ReadRealWaveFile);
    ReadRealWaveFile = fopen(readwavename, "r");
    printf("CHECK READ WAVE: %s\n", readwavename);
    fscanf(ReadRealWaveFile, "%lf %lf %lf\n", &WaveAngleIn, &WavePeriodIn,
           &WaveHeightIn);
    printf("Angle= %lf  Period= %lf Height= %lf \n", WaveAngleIn, WavePeriodIn,
           WaveHeightIn);
    if (countout == 0) exit(0);
  }

  if (debug18)
    printf("Angle= %lf  Period= %lf Height= %lf \n", WaveAngleIn, WavePeriodIn,
           WaveHeightIn);

  if (current_time_step == stop_after)
    fclose(ReadRealWaveFile); /*Keep file open till end of run */
}

void WaveOutFile(void)
/* Input Wave Distribution */
{
  if (current_time_step == 0) WaveFileOutput = fopen(wavesavename, "w");

  fprintf(WaveFileOutput, "%lf %lf %lf\n", (-WaveAngle * RAD_TO_DEG), Period,
          OffShoreWvHt);

  if (current_time_step == stop_after)
    fclose(WaveFileOutput); /*Keep file open till end of run */
}

void ControlFile(void)
/* Input Wave Distribution */
{
  OffShoreWvHt = -999;
  Period = -999;
  Asym = -999;
  Highness = -999;
  time_step = -999;
  Duration = -999;
  stop_after = -999;
  SlowWeatherCoeff = -999;
  FastWeatherCoeff = -999;
  PercentFineFast = -999;
  PercentFineSlow = -999;
  ErosionRatePerYear = -999;

  ReadControlFile = fopen(readcontrolname, "r");
  printf("READING CONTROL FILE\n");
  fscanf(ReadControlFile, "%lf %*[^\n]", &OffShoreWvHt);
  fscanf(ReadControlFile, "%lf %*[^\n]", &Period);
  fscanf(ReadControlFile, "%lf %*[^\n]", &Asym);
  fscanf(ReadControlFile, "%lf %*[^\n]", &Highness);
  fscanf(ReadControlFile, "%lf %*[^\n]", &time_step);
  fscanf(ReadControlFile, "%lf %*[^\n]", &Duration);
  fscanf(ReadControlFile, "%lf %*[^\n]", &stop_after);
  fscanf(ReadControlFile, "%lf %*[^\n]", &SlowWeatherCoeff);
  fscanf(ReadControlFile, "%lf %*[^\n]", &FastWeatherCoeff);
  fscanf(ReadControlFile, "%lf %*[^\n]", &PercentFineFast);
  fscanf(ReadControlFile, "%lf %*[^\n]", &PercentFineSlow);
  fscanf(ReadControlFile, "%lf %*[^\n]", &ErosionRatePerYear);
  fscanf(ReadControlFile, "%lf %*[^\n]", &coastrotation);
  fscanf(ReadControlFile, "%lf %*[^\n]", &waveheightchange);
  fscanf(ReadControlFile, "%lf", &waveperiodchange);
  fclose(ReadControlFile);

  if (time_step == -999) printf("***Initialisation file not read!***\n");

  if (METADATA == 'y') {
    printf("Outputting Initialization Metadata \n");

    InitMetaFile = fopen(metasavename, "w");
    fprintf(InitMetaFile, "OffShoreWvHt = %lf\n", OffShoreWvHt);
    fprintf(InitMetaFile, "Period = %lf\n", Period);
    fprintf(InitMetaFile, "Asym = %lf\n", Asym);
    fprintf(InitMetaFile, "Highness = %lf\n", Highness);
    fprintf(InitMetaFile, "Timstep = %lf\n", time_step);
    fprintf(InitMetaFile, "Duration = %lf\n", Duration);
    fprintf(InitMetaFile, "stop_after = %lf\n", stop_after);
    fprintf(InitMetaFile, "SlowWeatherCoeff = %lf\n", SlowWeatherCoeff);
    fprintf(InitMetaFile, "HighWeatherCoeff = %lf\n", FastWeatherCoeff);
    fprintf(InitMetaFile, "PercentFineFast = %lf\n", PercentFineFast);
    fprintf(InitMetaFile, "PercentFineSlow = %lf\n", PercentFineSlow);
    fprintf(InitMetaFile, "ErosionRatePerYear = %lf\n", ErosionRatePerYear);
    fprintf(InitMetaFile, "coastrotation = %lf\n", coastrotation);
    fprintf(InitMetaFile, "waveheightchange = %lf\n", waveheightchange);
    fprintf(InitMetaFile, "waveperiodchange = %lf\n", waveperiodchange);
    fclose(InitMetaFile);
  }

  printf("Control file read! \n");

  /*Write Metadata */
}

int getIndex(int x, int y, char interface)
/* don't use this during normal model run, if possible--it's too slow. as of
   10/05 it's only used during a pause */
{
  int i;
  if (interface == 's') {
    if (TotalBeachCells == 0) {
      /* printf("TotalBeachCells == 0 (probably beach ");
         printf("stuff hasn't been done yet) so no index\n"); */
      return -1;
    }
    for (i = 0; i < TotalBeachCells; i++) {
      if (X[i] == x && Y[i] == y) {
        /* printf("getIndex (s) found index to be %d \n", i); */
        return i;
      }
    }
    /* printf("getIndex (s) didn't find one\n"); */
    return -1;
  } else if (interface == 'r') {
    if (TotalRockCells == 0) {
      /* printf("TotalRockCells == 0 (probably beach ");
         printf("stuff hasn't been done yet) so no index\n"); */ return -1;
    }
    for (i = 0; i < TotalRockCells; i++) {
      if (XRock[i] == x && YRock[i] == y) {
        /* printf("getIndex (r) found index to be %d \n", i); */
        return i;
      }
    }
    /* printf("getIndex (r) didn't find one\n"); */
    return -1;
  } else {
    printf("have to call getIndex for sand 's' or rock 'r'\n");
    PauseRun(-1, -1, -1);
    return -1;
  }
}

void Delay(void) {
  int i;
  for (i = 0; i < 1000000; i++) {
    i++;
    i--;
  } /* slow it down so it doesn't print lots of times */
}

void PrintLocalConds(int x, int y, int in)

/* Prints Local Array Conditions aound x,y */
{
  int i, j, k, isee;
  float vol = CELL_WIDTH * CELL_WIDTH * DEPTH_SHOREFACE;

  printf("\n x: %d  y: %d  z: %d\n\n", x, y, in);

  /* if not given location along beach, look to see if along beach */

  if (in < 0) {
    for (i = 0; i <= TotalBeachCells; i++)
      if ((X[i] == x) && (Y[i] == y)) isee = i;
  } else
    isee = in;

  for (i = x + 2; i > x - 3; i--) {
    for (j = y - 2; j < y + 3; j++) {
      printf("  	%d,%d", i, j);
    }
    printf("\n");
  }

  printf("\n");

  printf("AllBeach\n");
  for (i = x + 2; i > x - 3; i--) {
    for (j = y - 2; j < y + 3; j++) {
      printf("	%c", AllBeach[i][j]);
    }
    printf("\n");
  }

  printf("\n");

  printf("AllRock\n");
  /*LMV*/ for (i = x + 2; i > x - 3; i--) {
    for (j = y - 2; j < y + 3; j++) {
      printf("	%c", AllRock[i][j]);
    }
    printf("\n");
  }

  printf("\n");

  printf("PercentFullSand\n");
  for (i = x + 2; i > x - 3; i--) {
    for (j = y - 2; j < y + 3; j++) {
      printf("	%2.5f", PercentFullSand[i][j]);
    }
    printf("\n");
  }

  printf("\n");

  printf("PercentFullRock\n");
  /*LMV*/ for (i = x + 2; i > x - 3; i--) {
    for (j = y - 2; j < y + 3; j++) {
      printf("	%2.5f", PercentFullRock[i][j]);
    }
    printf("\n");
  }

  printf("\n");
  printf("current_time_step: %d\n", current_time_step);
  /*LMV*/ printf("InitialMass: %f\n", MassInitial);
  /*LMV*/ printf("CurrentMass: %f\n", MassCurrent);
  /*LMV*/ printf(" %d  ", in);

  if (isee >= 0) {
    for (k = in - 3; k < in + 4; k++) {
      printf("  	%2d: %2d,%2d", k, X[k], Y[k]);
    }
    printf("\n\n\n");

    printf("Wave Angle:	%f\n\n", WaveAngle * RAD_TO_DEG);
    printf(
        "i		%d		%d		!%d		%d	"
        "	%d\n",
        in - 2, in - 1, in, in + 1, in + 2);
    printf(
        "Shadow		%c		%c		%c		"
        "%c		%c\n",
        InShadow[in - 2], InShadow[in - 1], InShadow[in], InShadow[in + 1],
        InShadow[in + 2]);
    printf(
        "Upwind		%c		%c		%c		"
        "%c		%c\n",
        UpWind[in - 2], UpWind[in - 1], UpWind[in], UpWind[in + 1],
        UpWind[in + 2]);
    printf(
        "Angle		%2.2f		%2.2f		%2.2f		%2.2f	"
        "	%2.2f\n",
        ShorelineAngle[in - 2] * RAD_TO_DEG, ShorelineAngle[in - 1] * RAD_TO_DEG,
        ShorelineAngle[in] * RAD_TO_DEG, ShorelineAngle[in + 1] * RAD_TO_DEG,
        ShorelineAngle[in + 2] * RAD_TO_DEG);
    printf(
        "Vol In 		%2.2f		%2.2f		%2.2f		"
        "%2.2f		%2.2f\n",
        VolumeIn[in - 2], VolumeIn[in - 1], VolumeIn[in], VolumeIn[in + 1],
        VolumeIn[in + 2]);
    printf(
        "Vol Out		%2.0f		%2.0f		%2.0f		"
        "%2.0f		%2.0f\n",
        VolumeOut[in - 2], VolumeOut[in - 1], VolumeOut[in], VolumeOut[in + 1],
        VolumeOut[in + 2]);
    printf(
        "Diff		%2.0f		%2.0f		%2.0f		%2.0f	"
        "	%2.0f\n",
        VolumeIn[in - 2] - VolumeOut[in - 2],
        VolumeIn[in - 1] - VolumeOut[in - 1], VolumeIn[in] - VolumeOut[in],
        VolumeIn[in + 1] - VolumeOut[in + 1],
        VolumeIn[in + 2] - VolumeOut[in + 2]);
    printf(
        "Frac Diff	%2.3f		%2.3f		%2.3f		%2.3f	"
        "	%2.3f\n",
        (VolumeIn[in - 2] - VolumeOut[in - 2]) / vol,
        (VolumeIn[in - 1] - VolumeOut[in - 1]) / vol,
        (VolumeIn[in] - VolumeOut[in]) / vol,
        (VolumeIn[in + 1] - VolumeOut[in + 1]) / vol,
        (VolumeIn[in + 2] - VolumeOut[in + 2]) / vol);
  }

  printf("\n");
}

/****
Over wash functions, copied from Kenny's code (so from Andrew Ashton's code)
*****/

void CheckOverwashSweep(void)
/* Just a loop to call overwash check founction CheckOverwash
 * Nothing done here, but can be down when CheckOVerwash is called
 */
{
  int i, ii;
  int sweepsign;

  if (RandZeroToOne() * 2 > 1) {
    sweepsign = 1;
  } else {
    sweepsign = 0;
  }

  OWflag = 0;

  for (i = 1; i < TotalBeachCells - 1; i++) {
    if (sweepsign == 1)
      ii = i;
    else
      ii = TotalBeachCells - 1 - i;

    /* To do test shoreline should be facing seaward */
    /* don't worry about shadow here, as overwash is not set to a time scale
     * with AST       */

    if ((fabs(SurroundingAngle[ii]) < (OverwashLimit / RAD_TO_DEG)) &&
        (InShadow[ii] == 'n')) {
      CheckOverwash(ii);
    }
  }
}

void CheckOverwash(int icheck)
/* Step back pixelwise in direction of SurroundingAngle[] to check 'needage' (?)
 * icheck is passed from CheckOverwashSweep() and is the index of the shoreline
 * cell being checked for overwash
 * Calls DoOverwash(). 'x' and 'y' hold real-space (floating point) values, and
 * will be mapped onto integer array
 */
{
  float slope;            /* slope of zero goes staight back */
  int ysign;              /* holder for going left or right alongshore */
  float x, y;             /* holders for 'real' location of x and y */
  float xin, yin;         /* starting 'real' locations */
  int xtest, ytest;       /* cell looking at */
  float xint, yint;       /* intercepts of overwash line in overwashable cell */
  int NextXInt, NextYInt; /* holder vairables for cell to check */
  float Ydown, DistanceDown; /* when going to next x cell, what other values */
  float Xside, DistanceSide; /* when gpoing to next y cell,other values */
  float
      checkdistance; /* distance of checking line- minimum, not actual width,
                        ends loop */
  float measwidth;   /* actual barrier width between cells */
  int AllBeachFlag; /* flag to see if overwash line has passed over at least one
                       AllBeach cell */

  /* convert angle to a slope and the direction of steps */
  /* note that for case of shoreline, positive angle will be minus y direcyion
   */
  if (SurroundingAngle[icheck] == 0.0) {
    /* unlikely, but make sure no div by zero */
    slope = 0.00001;
  } else if (fabs(SurroundingAngle[icheck]) == 90.0) {
    slope = 9999.9;
  } else {
    slope = fabs(tan(SurroundingAngle[icheck]));
  }

  if (SurroundingAngle[icheck] > 0)
    ysign = 1; /* check right */
  else
    ysign = -1; /* check left */

  if (AllBeach[X[icheck] - 1][Y[icheck]] ==
          'y' /* if the cell below is all beach */
      || ((AllBeach[X[icheck]][Y[icheck] - 1] ==
           'y') /* or if the cell to the left is all beach */
          && (AllBeach[X[icheck]][Y[icheck] + 1] ==
              'y'))) /* and the cell to the right is all beach */
  /* 'regular condition' */
  /* plus 'stuck in the middle' situation (unlikely scenario) */
  {
    xin = X[icheck] + PercentFullSand[X[icheck]][Y[icheck]];
    yin = Y[icheck] + 0.5;
  } else if (AllBeach[X[icheck]][Y[icheck] - 1] ==
             'y') /* if the cell to the left is all beach */
  /* on right side */
  {
    xin = X[icheck] + 0.5;
    yin = Y[icheck] + PercentFullSand[X[icheck]][Y[icheck]];
  } else if (AllBeach[X[icheck]][Y[icheck] + 1] == 'y')
  /* on left side */
  {
    xin = X[icheck] + 0.5;
    yin = Y[icheck] + 1.0 - PercentFullSand[X[icheck]][Y[icheck]];
  } else
  /* underneath, no overwash */
  {
    return;
  }

  x = xin;
  y = yin;
  checkdistance = 0;
  AllBeachFlag = 0;

  while ((checkdistance < CritBWidth) && (y > 0) && (y < 2 * Y_MAX) && (x > 1)) {
    NextXInt = ceil(x) - 1;
    if (ysign > 0) {
      NextYInt = floor(y) + 1;
    } else {
      NextYInt = ceil(y - 1);
    }

    /* moving to next whole 'x' position, what is y position? */
    Ydown = y + (x - NextXInt) * slope * ysign;
    DistanceDown = Raise(
        ((Ydown - y) * (Ydown - y) + (NextXInt - x) * (NextXInt - x)), .5);

    /* moving to next whole 'y' position, what is x position? */
    Xside = x - fabs(NextYInt - y) / slope;
    DistanceSide = Raise(
        ((NextYInt - y) * (NextYInt - y) + (Xside - x) * (Xside - x)), .5);

    if (DistanceDown < DistanceSide)
    /* next cell is the down cell */
    {
      x = NextXInt;
      y = Ydown;
      xtest = NextXInt - 1;
      ytest = floor(y);
    } else
    /* next cell is the side cell */
    {
      x = Xside;
      y = NextYInt;
      xtest = floor(x);
      ytest = y + (ysign - 1) / 2;
    }

    checkdistance =
        Raise(((x - xin) * (x - xin) + (y - yin) * (y - yin)), .5) * CELL_WIDTH;

    if (AllBeach[xtest][ytest] == 'y') AllBeachFlag = 1;

    if ((AllBeach[xtest][ytest] == 'n') && (AllBeachFlag) &&
        !(((X[icheck] - xtest) > 1) || (abs(ytest - Y[icheck]) > 1)))
    /* if passed through an allbeach and a neighboring partial cell, jump out,
       only bad things follow */
    {
      return;
    }

    if ((AllBeach[xtest][ytest] == 'n') && (AllBeachFlag) &&
        (xtest < X[icheck]) &&
        (((X[icheck] - xtest) > 1) || (abs(ytest - Y[icheck]) > 1)))
    /* Looking for shore cells, but don't want immediate neighbors, and go
       backwards */
    /* Also mush pass though an allbeach cell along the way */
    {
      if (AllBeach[xtest + 1][ytest] == 'y')
      /* 'regular condition' - UNDERNEATH, here */
      {
        xint = (xtest + 1 - PercentFullSand[xtest][ytest]);
        yint = yin + (xin - xint) * ysign * slope;

        if ((yint > ytest + 1.0) || (yint < ytest))
        /* This cell isn't actually an overwash cell */
        {
          measwidth = CritBWidth;
        } else {
          measwidth = CELL_WIDTH * Raise((xint - xin) * (xint - xin) +
                                            (yint - yin) * (yint - yin),
                                        0.5);
        }
      } else if (AllBeach[xtest][ytest - 1] == 'y')
      /* on right side */
      {
        yint = (ytest + PercentFullSand[xtest][ytest]);
        xint = xin - fabs(yin - yint) / slope;

        if (xint < xtest)
        /* This cell isn't actually an overwash cell */
        {
          measwidth = CritBWidth;
        } else {
          measwidth = CELL_WIDTH * Raise((xint - xin) * (xint - xin) +
                                            (yint - yin) * (yint - yin),
                                        0.5);
        }
      } else if (AllBeach[xtest][ytest + 1] == 'y')
      /* on left side */
      {
        yint = (ytest + 1 - PercentFullSand[xtest][ytest]);
        xint = xin - fabs(yin - yint) / slope;

        if (xint < xtest)
        /* This cell isn't actually an overwash cell */
        {
          measwidth = CritBWidth;
        } else {
          measwidth = CELL_WIDTH * Raise((xint - xin) * (xint - xin) +
                                            (yint - yin) * (yint - yin),
                                        0.5);
        }
      } else if (AllBeach[xtest - 1][ytest] == 'y')
      /* 'regular condition' */
      /* plus 'stuck in the middle' situation */
      {
        xint = (xtest + PercentFullSand[xtest][ytest]);
        yint = yin + (xin - xint) * ysign * slope;

        if ((yint > ytest + 1.0) || (yint < ytest))
        /* This cell isn't actually an overwash cell */
        {
          measwidth = CritBWidth;
        } else {
          measwidth = CELL_WIDTH * Raise((xint - xin) * (xint - xin) +
                                            (yint - yin) * (yint - yin),
                                        0.5);
        }
      } else if (PercentFullSand[xtest][ytest] > 0)
      /* uh oh - not good situation, no allbeach on sides */
      /* assume this is an empty cell, */
      {
        xint = x;
        yint = y;

        measwidth = CELL_WIDTH * Raise((xint - xin) * (xint - xin) +
                                          (yint - yin) * (yint - yin),
                                      0.5);

        printf(
            "-- Some Odd Over  xin: %2.2f  yin: %2.2f xt:%d yt: %d xint: %f "
            "yint: %f Meas: %3.2f Ang: %f Abs: %f\n",
            xin, yin, xtest, ytest, xint, yint, measwidth,
            SurroundingAngle[icheck] * RAD_TO_DEG,
            fabs(SurroundingAngle[icheck]) * RAD_TO_DEG);
      } else
      /* empty cell - oughta fill er up  - fill max barrier width */
      {
        xint = x;
        yint = y;
        measwidth = CritBWidth - CELL_WIDTH;

        printf(
            "-- Empty Odd Over  xin: %2.2f  yin: %2.2f xt:%d yt: %d xint: %f "
            "yint: %f Meas: %3.2f Ang: %f Abs: %f\n",
            xin, yin, xtest, ytest, xint, yint, measwidth,
            SurroundingAngle[icheck] * RAD_TO_DEG,
            fabs(SurroundingAngle[icheck]) * RAD_TO_DEG);
      }

      checkdistance = measwidth;

      if (measwidth < CritBWidth) {
        DoOverwash(X[icheck], Y[icheck], xtest, ytest, xint, yint, measwidth,
                   icheck);
        /* jump out of loop */
        OWflag = 1;
        return;
      }
    }
  }
}

void DoOverwash(int xfrom, int yfrom, int xto, int yto, float xintto,
                float yintto, float widthin, int ishore)

/*  given a cell where overwash is needed ,move sediment back
 *  for 'true' overwash based on shoreline angles
 *  will change and use PercentFull[][] and AllBeach [][]
 */
{
  float BBneed, delBB, delShore; /* local variables */

  /*float DepthBackBarrier = 6.0;      m current set depth for backbarrier (temp
   * - make into function) */

  float DepthBB; /* holds effective backbarrier depth */

  DepthBB = GetOverwashDepth(xto, yto, xintto, yintto, ishore);

  /* calculated value of most that backbarrier can move given geometry (true,
   * non-iterative solution) */

  if (DepthBB == DEPTH_SHOREFACE) {
    BBneed = MaxOver;
  } else {
    BBneed =
        (CritBWidth - widthin) / CELL_WIDTH / (1 - (DepthBB / DEPTH_SHOREFACE));
  }

  if (BBneed <= MaxOver)
  /* do all overwash */
  {
    delShore = BBneed * DepthBB / DEPTH_SHOREFACE;
    delBB = BBneed;
  } else
  /* only do overwash to max change) */
  {
    delShore = MaxOver * DepthBB / DEPTH_SHOREFACE;
    delBB = MaxOver;
  }

  PercentFullSand[xto][yto] += delBB;
  PercentFullSand[xfrom][yfrom] -= delShore;

  if (PercentFullSand[xto][yto] > 1) {
    OopsImFull(xto, yto);
  }
  if (PercentFullSand[xfrom][yfrom] < 0) {
    OopsImEmpty(xfrom, yfrom);
  }
}

float GetOverwashDepth(int xin, int yin, float xinfl, float yinfl, int ishore)
/* Rountine finds corresponding overwash depths
 *   	OWType = 0 take the depth at neightbor to the backing cell
 *   	OWType = 1 geometric rule based upon distance from back to shoreline
 */
{
  int xdepth;
  float Depth;
  float BBDistance;          /* Distance from backshore to next shore */
  float slope;               /* slope of zero goes staight back */
  int ysign;                 /* holder for going left or right alongshore */
  float x, y;                /* holders for 'real' location of x and y */
  int xtest, ytest;          /* cell looking at */
  int NextXInt, NextYInt;    /* holder variables for cell to check */
  float Ydown, DistanceDown; /* when going to next x cell, what other values */
  float Xside, DistanceSide; /* when gpoing to next y cell,other values */
  int BackFlag;              /* Flag to indicate if hit backbarrier */
  int Backi = -1;            /* i for backbarrier intersection */
  int i, j;                  /* counters */
  int FoundFlag;             /* Backbarrier intersection flag */
  float AngleSin, AngleUsed;

  if (OWType == 0)
  /* Use Cell Depths for overwash depths */
  {
    xdepth = xin;
    Depth = cell_depth[xdepth][yin];

    while ((Depth < 0) && (xdepth > 0)) {
      Depth = cell_depth[xdepth][yin];
      xdepth--;
    }

    if (Depth == DEPTH_SHOREFACE) {
      Depth = 6.0;
    }

    return Depth;

  } else if (OWType == 1)
  /* Geometric relation to determine depth through intersection of shorefaces */
  /* look in line determined by shoreline slope - reuse stepping function
     (again) */
  {
    x = xinfl;
    y = yinfl;

    if (SurroundingAngle[ishore] == 0.0) {
      /* unlikely, but make sure no div by zero */
      slope = 0.00001;
    } else if (fabs(SurroundingAngle[ishore]) == 90.0) {
      slope = 9999.9;
    } else {
      slope = fabs(tan(SurroundingAngle[ishore]));
    }

    BackFlag = 0;
    if (SurroundingAngle[ishore] > 0)
      ysign = 1;
    else
      ysign = -1;

    while ((!BackFlag) && (y > 0) && (y < 2 * Y_MAX) && (x > 1)) {
      NextXInt = ceil(x) - 1;
      if (ysign > 0)
        NextYInt = floor(y) + 1;
      else
        NextYInt = ceil(y - 1);

      /* moving to next whole 'x' position, what is y position? */
      Ydown = y + (x - NextXInt) * slope * ysign;
      DistanceDown = Raise(
          ((Ydown - y) * (Ydown - y) + (NextXInt - x) * (NextXInt - x)), .5);

      /* moving to next whole 'y' position, what is x position? */
      Xside = x - fabs(NextYInt - y) / slope;
      DistanceSide = Raise(
          ((NextYInt - y) * (NextYInt - y) + (Xside - x) * (Xside - x)), .5);

      if (DistanceDown < DistanceSide)
      /* next cell is the down cell */
      {
        x = NextXInt;
        y = Ydown;
        xtest = NextXInt - 1;
        ytest = floor(y);
      } else
      /* next cell is the side cell */
      {
        x = Xside;
        y = NextYInt;
        xtest = floor(x);
        ytest = y + (ysign - 1) / 2;
      }

      if (PercentFullSand[xtest][ytest] > 0) BackFlag = 1;
    }

    /* Try to find the i for the cell found */
    /* If you have a better idea how to do this, go ahead */

    i = 2;
    FoundFlag = 0;

    while ((i < TotalBeachCells - 1) && !(FoundFlag)) {
      if ((X[i] == xtest) && (Y[i] == ytest)) {
        FoundFlag = 1;
        Backi = i;
      }
      i++;
    }

    if (!BackFlag)
    /* The search for the backbarrier went out of bounds - not good, assume big
       = depthshoreface    */
    /* Periodic B.C.'s should make this not so important */
    {
      Depth = DEPTH_SHOREFACE;
    } else {
      BBDistance = Raise(((xinfl - xtest - PercentFullSand[xtest][ytest]) *
                          (xinfl - xtest - PercentFullSand[xtest][ytest])) +
                             ((yinfl - ytest - 0.5) * (yinfl - ytest - 0.5)),
                         .5);

      if (!FoundFlag)
      /* The backbarrier intersection isn't on the shoreline */
      /* Assume 1/2 of the length applies to this case */
      {
        Depth = BBDistance / 2 * SHOREFACE_SLOPE * CELL_WIDTH;
      } else
      /* Use the fancy geometry thing */
      {
        AngleUsed = 0;
        for (j = -1; j < 2; j++) {
          AngleUsed += SurroundingAngle[Backi + j];
        }
        AngleUsed = AngleUsed / 5;

        if (fabs(AngleUsed) > PI / 4.0) {
          AngleUsed = PI / 4.0;
        }

        AngleSin = sin(PI / 2.0 - fabs(SurroundingAngle[ishore] + AngleUsed));

        Depth = BBDistance * AngleSin / (1 + AngleSin);
      }
    }

    if (Depth < OWMinDepth) {
      Depth = OWMinDepth;
    } else if (Depth > DEPTH_SHOREFACE) {
      Depth = DEPTH_SHOREFACE;
    }
    return Depth;
  }

  printf("OWDepth all broken");
  PauseRun(xin, yin, -1);
  return DEPTH_SHOREFACE;
}

float ConvertAngle(float EvaluteAngle, int type)
/* This function takes
 1) the EvaluateAngle defined in CEM and returns an azimuth, or
 2) the azimuthal angle from SWAN output and returns an angle in the CEM
 convention.
 The int 'type' separates the two actions ('type = 1' triggers the first
 statement, and
 'type = anything except 1' triggers the second statement)

 #SWAN */
{
  float value; /* What direction to look offshore to find a breaking wave? */

  /*if (debugSWAN)*/
  /*printf("\n incoming angle: %f \n", EvaluateAngle);*/

  if (type == 1) {
    /* Moving right and down along shoreline -- most common case, unless we have
     spits or cape
     tip overgrowth */
    if (EvaluateAngle > (-90 * (PI / 180)) &&
        EvaluateAngle < (0 * (PI / 180))) {
      value = EvaluateAngle + (90 * (PI / 180));
    }

    /* Moving right and up (or just right) along shoreline -- most common case,
     unless we have
     spits or cape tip overgrowth */
    else if (EvaluateAngle >
             (0 * (PI / 180) && EvaluateAngle < (90 * (PI / 180)))) {
      value = (360 * (PI / 180)) - EvaluateAngle;
    }

    /* Moving right, stright shoreline (dx = 0) */
    else if (EvaluateAngle == 0) {
      value = 0;
    }

    /* Moving straight up or straight down along shoreline */
    else if (EvaluateAngle == (-90 * (PI / 180)) ||
             EvaluateAngle == (90 * (PI / 180))) {
      if (EvaluateAngle == (-90 * (PI / 180))) {
        value = (270 * (PI / 180));
      } else if (EvaluateAngle == (90 * (PI / 180))) {
        value = (90 * (PI / 180));
      }
    }

    /* Moving left (and up, down, or straight) along shoreline */
    else if (EvaluateAngle >
             (-270 * (PI / 180) && EvaluateAngle < (-90 * (PI / 180)))) {
      value = fabs(EvaluateAngle);
    }

    /* Oops -- did not find a search direction! */
    else {
      printf(
          "\n Can't look offshore for SWANs 'cause i didn't find a surrounding "
          "angle! (in ParseSWAN) \n");
    }

    /* Print to screen, if you like */
    /*if (debugSWAN) printf("\n search direction: %f , surrounding angle: %f
     \n",
     value, EvaluateAngle);*/

    /* Return value. */
    return value;
  }

  /* -----> Types other than 1 <----- */
  /* Convert SWAN azimuthal angle */
  else {
    /* Wave going right to left */
    if (EvaluateAngle > 0 && EvaluateAngle < 90) {
      return -EvaluateAngle * (PI / 180);
    }
    /* Wave going left to right */
    else if (EvaluateAngle < 360 && EvaluateAngle > 270) {
      return (360 - EvaluateAngle) * (PI / 180);
    }
    /* Wave coming straight onshore */
    else if (EvaluateAngle == 0 || EvaluateAngle == 360) {
      return 0;
    }
  }
  return 0.;
}

/* SWAN data parse function! PWL, 10-18-13. #SWAN */
void ParseSWAN(int ShoreAngleLoc, float ShoreAngle) {
  // float LookOffshore; /* What direction to look offshore to find a breaking wave? */
  /* float	LookSlope; */ /* What slope (x,y) is the line of sight? */
  int xcoord, ycoord;
  /* int NextInLine; */
  /* int		OopsImBroke = 0;*/ /* Flag that indicates if we've found a
                            broken wave.
                            0 = keep lookin', pal; 1 = eureka! */
  /* float		interval = 0.1; */   /* How far to search in each
                                                  iteration */
  /* float	xdist = 0.0, ydist = 0.0; */ /* Distance to search for SWAN
                                                  info (units are cells). */
  /* float	Xpos, Ypos; */               /* Tracks search distance */

  /* Tracks cells that have already been searched */
  /* int	LastYCell; */
  int LastXCell;

  /*int		CurrentXCell, CurrentYCell; */ /* Tracks integer value cell
                                                    position */
  double Hd; /* Wave height divided by water depth */
  float
      CAngle; /* Angle converted from ConvertAngle function -- extraneous for
                 dubugging */
  /* float	RAngle; */ /* Angle retrieved from SWAN -- extraneous for
                                dubugging */
  int x; /* iterator for simple search function, 11-12-13 -- can delete for
            shoreline angle-based version */

  xcoord = X[ShoreAngleLoc];
  ycoord = Y[ShoreAngleLoc];

  if (ShoreAngleLoc == 0) {
    printf("\n no alongshore location, barf \n");
  }

  /* FIRST STEP! find the direction to look offshore. Should be perpendicular to
   the shoreline
   orientation defined by the variable SurroundingAngle. The convention is
   strange, so use
   'ConvertAngle' function below to separate different cases. The end result
   will be an azimuthal
   search direction, where 0 (& 360) is straight up (seaward) */
  /* NOTE, 10-22-13: i wonder if the calcs below will be precise enough to
   satisfy the conditional
   statements? */

  /* 11-14-13: This doesn't work to calculate Qs because Qs needs a sign to
   * indicate direction...*/
  /* This conversion will be useful again when the angle-based search function
   * is used. */

  /* Convert shoreline angle to azimuth */
  /*LookOffshore = ConvertAngle(SurroundingAngle[ShoreAngleLoc], 1);*/

  /* i think this is incorrect -- PWL 11/20/14 */
  /* Don't use Surrounding Angle for Qs! This method is ignoring upwind
   * conditions */
  // LookOffshore = SurroundingAngle[ShoreAngleLoc];

  /* NEXT STEP */
  /* Start at middle of cell and search seaward. Find next cell in path. Might
   consider starting
   search from actual shoreline position (PercentFull)? */

  /* Initial conditions for loop */
  /*Xpos = (xcoord) + 0.5;*/
  /*Ypos = (ycoord) + 0.5;*/
  /* Distance to search for each loop iteration */
  /*xdist = (interval*cos(LookOffshore)); NOTE cross-shore direction --
   * confusing convention! */
  /*ydist = (interval*sin(LookOffshore)); alongshore direction */
  /* Initial conditions for LastCell */
  /*LastXCell = xcoord;*/
  /*LastYCell = ycoord;*/

  /* --> Main search loop <-- */
  /*while (!OopsImBroke)*/
  /* CAN DELETE FOR LOOP AND REPLACE WITH WHILE LOOP WHEN USING ANGLE-BASED
   * FUNCTION */
  for (x = 0; x < X_MAX; x++) {
    /* What cell am i in? */
    /*Xpos = Xpos + xdist;*/
    /*Ypos = Ypos + ydist;*/

    /* Did i enter a new cell? */
    /*if (floor(Xpos) != LastXCell || floor(Ypos) != LastYCell)
     {*/
    /* FOUND A NEW CELL */

    /* Is this a breaking wave cell? */
    /*CurrentXCell = floor(Xpos);*/
    /*CurrentYCell = floor(Ypos);*/

    /* Can remove the lines below - it's for simple search function, 11-12-13 */
    LastXCell = xcoord;
    xcoord = xcoord + x;
    /* ------ */
    ydebug[ShoreAngleLoc] = 6;
    /* Is this a breaking wave cell? Also, make sure not to divide by zero or
     * else your computer will explode */
    /* May need to lower breaking wave threshold and/or change resolution of
     * nested grid -- do some tests! #SWAN */

    if (shelf_depth[xcoord][ycoord] > 0 && shelf_depth[LastXCell][ycoord] > 0 &&
        wave_h_sig[LastXCell][ycoord] > 0) {
      Hd = wave_h_sig[xcoord][ycoord] / (shelf_depth[xcoord][ycoord]);

      /*if (debugSWAN)
                   printf("H/D threshold: %f\n", Hd);*/
      ydebug[ShoreAngleLoc] = 1;

      if (Hd < WaveBreakDepth && Hd > 0 && shelf_depth[LastXCell][ycoord] > 0) {
        /* 11-12-13: When using angle-based function, replace 'ycoord' with
         * 'LastYCell' */

        /* Wave broke! debug if necessary and then step back and take values
         * from previous cell */
        /*LookOffshore = LookOffshore*(180/PI); Convert angle to degrees because
         * SWAN is in degrees */

        /* Stepping back... */
        /* UPDATE, 11/24/14: don't step back! Replaced [LastXCell] with
         * [xcoord]*/
        WvHeight = wave_h_sig[xcoord][ycoord];
        EvaluateAngle = wave_dir[xcoord][ycoord];
        CAngle = ConvertAngle(
            EvaluateAngle,
            2); /* Convert angle from azimuth (degrees) to CEM-style (rads) */
        /*Angle = CAngle - LookOffshore; Angle is made relative to shore */
        Angle = CAngle; /* NEW METHOD 11/20/14 PWL */
        BreakDepth = shelf_depth[xcoord][ycoord];
        /*OopsImBroke = 1;*/
        /*if (debugSWAN) */
        /*{*/
        /* Save to file! */
        /*printf("Saving SWAN shoreface bathy as: %s 		",
         * SWANsavename);*/
        /* Save stuff to text file, PWL 5-5-14 #SWAN */
        Hsigdebug[ShoreAngleLoc] = WvHeight;
        Dirdebug[ShoreAngleLoc] = Angle;
        Hddebug[ShoreAngleLoc] = BreakDepth;
        xdebug[ShoreAngleLoc] = xcoord;
        ydebug[ShoreAngleLoc] = ShoreAngleLoc;
        /*}*/

        break; /* delete break and use OopsImBroke with angle-based function,
                  11-12-13 */
      }
    } else {
      ydebug[ShoreAngleLoc] = 5;
    }

    /* Not a breaking wave cell, update last cell position and keep loopin' */
    /*if (floor(Xpos) != LastXCell)
     {
     LastXCell = floor(Xpos);
     }
     if (floor(Ypos) != LastYCell)
     {
     LastYCell = floor(Ypos);
     }*/
  }
  /*}*/

  /* We're done here, muchacho. Return some values to the parent function --
   currently the
   values to be returned are global, so no further action needed... */
}
