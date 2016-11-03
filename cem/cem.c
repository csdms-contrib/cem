#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#include "consts.h"
#include "globals.h"
#include "rocks.h"
#include "utils.h"


/****************************************************************
 Global variables to be share with other files.
 ***************************************************************/
double **topography;                      // Holds cliff heights -- will change through time, eventually...

/****************************************************************
 Special SWAN parameters.
****************************************************************/
#if defined(WITH_SWAN)
static const char kSwanFlag = TRUE; // Is SWAN doing wave transformations?
#else
static const char kSwanFlag = FALSE;
#endif

double **shelf_depth = NULL;                // SWAN bathymetry
double **wave_h_sig = NULL;                 // SWAN wave heights.
double **wave_dir = NULL;                   // SWAN wave angles.

int current_time_step = 0;                  // Time step of current calculation
double current_time = 0.;                   // Current model time
double stop_after = 36500;                  // Stop after what number of time steps
double time_step = 1;                       // days; reflects rate of sediment transport per time step

static double BreakDepth;                   // Breaking wave depth found from SWAN run
static double EvaluateAngle;                // Temporary angle holder for the ConvertAngle function

/****************************************************************
 Global variables to be used only within this file.
****************************************************************/
static const double kAngleFactor = 1.;
static double cell_depth[X_MAX][2 * Y_MAX]; // Depth array
const double kCliffHeightSlow = 30;         // Cliff height above sea level for slow weathering rock PWL
const double kCliffHeightFast = 0;          // Cliff height above sea level for fast weathering rock PWL

/****************************************************************
 Overall Shoreface Configuration Arrays
/****************************************************************/
static int AllBeach[X_MAX][2 * Y_MAX];      // Flag indicating of cell is entirely beach
static int AllRock[X_MAX][2 * Y_MAX];       // Flag indicating if cell is entirely rock LMV
static double PercentFullSand[X_MAX][2 * Y_MAX];// Fractional amount of cell full of sediment LMV
static double PercentFullRock[X_MAX][2 * Y_MAX];// Fractional amount of a cell full of rock LMV
static char **type_of_rock = NULL;          // Array to control weathering rates of rock along the beach LMV
static int Age[X_MAX][2 * Y_MAX];           // Age since cell was deposited

// A sink is a cell that is routinely emptied (if it is on the beach)
static const int kSinkY[] = {158, 361, 158, 361, 158, 361};
static const int kSinkX[] = {190, 190, 189, 189, 191, 191};

/****************************************************************
 Data file information
****************************************************************/
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

/****************************************************************
 Computational Arrays (determined for each time step)  
 -- will eventually be in a structure for BMI_Model
****************************************************************/
// TODO: change limited value variables to enums
static int X[MAX_BEACH_LENGTH];                       // X Position of ith beach element
static int Y[MAX_BEACH_LENGTH];                       // Y Position of ith beach element
static int XBehind[MAX_BEACH_LENGTH];                 // Cell that is "behind" X[i] LMV
static int YBehind[MAX_BEACH_LENGTH];                 // Cell that is "behind" Y[i] LMV
static int XRock[MAX_BEACH_LENGTH];                   // X Position of ith rock element LMV
static int YRock[MAX_BEACH_LENGTH];                   // Y Position of ith rock element LMV
static int XRockBehind[MAX_BEACH_LENGTH];             // Cell that is "behind" XRock[i] LMV
static int YRockBehind[MAX_BEACH_LENGTH];             // Cell that is "behind" YRock[i] LMV
static int InShadow[MAX_BEACH_LENGTH];                // Is ith beach element in shadow?
static double ShorelineAngle[MAX_BEACH_LENGTH];       // Angle between cell and right (z+1) neighbor
static double SurroundingAngle[MAX_BEACH_LENGTH];     // Angle between left and right neighbors
// TODO: make enum
static char UpWind[MAX_BEACH_LENGTH];                 // Upwind or downwind condition used to calculate sediment transport
static double VolumeIn[MAX_BEACH_LENGTH];             // Sediment volume into ith beach element
static double VolumeOut[MAX_BEACH_LENGTH];            // Sediment volume out of ith beach element
static double VolumeAcrossBorder[MAX_BEACH_LENGTH];   // Sediment volume needed across border of ith beach element in m^3/day LMV
                                                          // not necessarily amount a cell gets
static double ActualVolumeAcross[MAX_BEACH_LENGTH];   // Sediment volume that actually gets passed across border LMV
                                                          // amount sed is limited by volumes across border upstream and downstream
static char DirectionAcrossBorder[MAX_BEACH_LENGTH];  // Flag to indicate if transport across border is L or R LMV
// TODO: make enum
static char FlowThroughCell[MAX_BEACH_LENGTH];        // Is flow through ith cell Left, Right, Convergent, or Divergent LMV
static double DistanceToBeach[MAX_BEACH_LENGTH];      // Distance in meters from rock to beach LMV
static double MinDistToBeach[MAX_BEACH_LENGTH];       // From a rock cell j, min distance (in meters) to closest beach LMV
static int ClosestBeach[MAX_BEACH_LENGTH];            // i position of closest rock to beach LMV
static double AmountWeathered[MAX_BEACH_LENGTH];      // Amount of rock weathered from rock cell j LMV

/****************************************************************
 Debugging variables -- for temporary debugging only, 5-5-14
 UPDATE 11/20/14 -- now using these for upwind scheme fixing, so keep 'em
 around (probably should rename...)
****************************************************************/
static double Hsigdebug[MAX_BEACH_LENGTH];
static double Dirdebug[MAX_BEACH_LENGTH];
static double Hddebug[MAX_BEACH_LENGTH];
static double xdebug[MAX_BEACH_LENGTH];
static double ydebug[MAX_BEACH_LENGTH];
static double WvHeight;
static double Angle;
static int OWflag = FALSE; // A debugging flag for overwash routines

/****************************************************************
 Miscellaneous Global Variables
****************************************************************/
// TODO: most of these probably shouldn't be global
static char StartFromFile = FALSE; // start from saved file?

static int NextX;                 // Global variables used to iterate FindNextCell in global array -
static int NextY;                 // would've used pointer but wouldn't work
static int BehindX;
static int BehindY;
static int BehindRockX;
static int BehindRockY;
static int NextRockX;             // Global variables used to iterate FindNextRockCell in global array LMV
static int NextRockY;
static int TotalBeachCells;       // Number of cells describing beach at particular iteration
static int TotalRockCells;        // Number of cells describing rock at an iteration LMV
static int ShadowXMax;            // used to determine maximum extent of beach cells
static double WaveAngle;          // wave angle for current time step
static int FindStart;             // Used to tell FindBeach at what Y value to start looking
static int FindRockStart;         // Used to tell FindRock at what Y value to start looking LMV
static int FellOffArray;          // Flag used to determine if accidentally went off array
static int FellOffRockArray;      // Flag used to determine if accidentally went off rock array LMV
static double MassInitial;        // For conservation of mass calcs
static double MassCurrent;

/****************************************************************
 Wave Climate 
****************************************************************/
static int NumWaveBins;           // For Input Wave - number of bins
static double WaveMax[36];        // Max Angle for specific bin
static double WaveProb[36];       // Probability of Certain Bin
static double WaveAngleIn;
static double WaveHeightIn;
static double WavePeriodIn;

static double Period = 10.0;      // seconds
static double OffShoreWvHt = 2.0; // meters
static double Asym = 0.70;        // Fractional portion of waves coming from positive (left) direction
static double Highness = 0.35;    // .5 = even dist, < .5 high angle domination. NOTE Highness actually determines lowness
static double Duration = 1.;      // Number of time steps calculations loop at same wave angle


static double coastrotation = 40.;    // Angle (deg) used to align coastline with real wave climate
static double waveheightchange = 1.;  // Factor applied to the input wave height; 1 = no change, 0.5 = half wave height, and 2 = double wave height
static double waveperiodchange = 1.;  // Factor applied to the input wave period; 1 = no change, 0.5 = half wave period, and 2 = double wave period

/****************************************************************
 Weathering parameters
****************************************************************/
static double SlowWeatherCoeff = 0. * NORMAL_WEATHERING_RATE;// Weathering rate of slow rock e.g. 0.5 = slow weathering, 1/2 as fast LMV
static double FastWeatherCoeff = 1 * NORMAL_WEATHERING_RATE;// Weathering rate of fast rock e.g. 2 = fast weathering, 2 times faster than normal
static double PercentFineFast = 0.0; // Percent of fast weathering rock lost because it is too fine to stay in nearshore LMV
static double PercentFineSlow = 0.0; // Percent of slow weathering rock lost because it is too fine to stay in nearshore LMV
static double ErosionRatePerYear = 0.5; // Amount of erosion to TotalBeachCells in m/year LMV


/****************************************************************
 Functions
****************************************************************/
void AdjustShore(int i);
void AgeCells(void);
void ControlFile(void);
double ConvertAngle(double EvaluateAngle, int type);
void DetermineAngles(void);
void DetermineSedTransport(void);
void DoSink(void);
void Transport(void);
void ErodeTheBeach(int i);
void FindBeachCells(int YStart);
void FindRockCells(int YStart);
int FindIfInShadow(int xin, int yin, int ShadMax);
void FindNearestBeach(int j);
void FindNextCell(int x, int y, int z);
void FindNextRockCell(int x, int y, int z);
double FindWaveAngle(void);
void FixBeach(void);
void FixFlow(void);
void FlowInCell(void);
int getIndex(int x, int y, char interface);
void InitConds(void);
void InitPert(void);
double MassCount(void);
void OopsImEmpty(int x, int y);
void OopsImFull(int x, int y);
void ParseSWAN(int ShoreAngleLoc, double ShoreAngle);
void PauseRun(int x, int y, int in);
void periodic_boundary_copy(void);
double Raise(double b, double e);
double RandZeroToOne(void);
void ReadSandFromFile(void);
void ReadWaveIn(void);
void RealWaveIn(void);
void RockCalculations(void);
void SaveSandToFile(void);
void SaveLineToFile(void);
void SedTrans(int From, int To, double ShoreAngle, int MaxT, int Last);
void ShadowSweep(void);
void TransportSedimentSweep(void);
void WaveOutFile(void);
void WeatherRock(int j);
int XMaxBeach(int Max);
void ZeroVars(void);

// Overwash functions
void CheckOverwashSweep(void);
double GetOverwashDepth(int xin, int yin, double xinfl, double yinfl, int ishore);
void CheckOverwash(int icheck);
void DoOverwash(int xfrom, int yfrom, int xto, int yto, double xintto, double yintto, double widthin, int ishore);

// Initialize conditions
void InitNormal(void);
int initBlock(void);
int InitWiggly(void);
void PrintLocalConds(int x, int y, int in);

// BMI function
int cem_initialize(void);
int cem_update(void);
int cem_finalize(void);

/**
 * Initializes full CEM run
 * PARAMETERS: none
 * RETURN: exit status
 */
int cem_initialize(void)
{
  ShadowXMax = X_MAX; /* RCL: was X_MAX-5; */

  if (SEED == -999)
  {
    srand(time(NULL));
  }
  else
  {
    srand(SEED);
  }

  /* Allocate memory for arrays. */
  const int n_rows = X_MAX;
  const int n_cols = 2 * Y_MAX;

  topography = (double**)malloc2d(n_rows, n_cols);
  shelf_depth = (double**)malloc2d(n_rows, n_cols);
  wave_h_sig = (double**)malloc2d(n_rows, n_cols);
  wave_dir = (double**)malloc2d(n_rows, n_cols);
  /* TODO: fix this */
  type_of_rock = (char**)malloc2d(n_rows, n_cols);

  /*  Start from file or not? */
  if (PROMPT_START)
  {
    printf("shall we start from a file (y or n)? \n");
    scanf("%c", &StartFromFile);

    if (StartFromFile)
    {
      printf("Starting Filename? \n");
      scanf("%24s", readfilename);
      printf("Saving Filename? \n");
      scanf("%s", savefilename);
      printf("What time step are we starting at?");
      scanf("%d", &current_time_step);
      ReadSandFromFile();
    }

    else
    {
      printf("Saving Filename? \n");
      scanf("%s", savefilename);
      InitConds();
      if (INITIAL_PERT)
      {
        InitPert();
      }
      printf("InitConds OK \n");
    }
  }

  else if (StartFromFile)
  {
    ReadSandFromFile();
  }

  else
  {
    InitConds();
    if (INITIAL_PERT)
    {
      InitPert();
    }
  }

  /* Set Periodic Boundary Conditions */
  periodic_boundary_copy();
  /* Count Initial Mass */
  MassInitial = MassCount();
  /* Fix the Beach (?) */
  FixBeach();

  /* PWL commented out to match ndelta4.c, 12/3/14 */
  // if (SAVE_LINE) SaveLineToFile();
  // if (SAVE_FILE) SaveSandToFile();

  /* Read in wave data? */
  if (WAVE_IN == 1)
  {
     ReadWaveIn();      /* Initialise WaveMax and WaveProb */
  }

  if (INITIALIZE_FILE)
  {
    ControlFile();      /* Read in initialisation data to control model */
  } 
  return 0;
}

/**
 * Update the CEM run by a single time step (times Duration when coupled with SWAN)
 * PARAMETERS: none
 * RETURN: exit status
 */
int cem_update(void)
{
  /* duration loop variable, moved from former (now deleted) 'main.c' above */
  int xx;
  /* Calculate Wave Angle */
  WaveAngle = FindWaveAngle();

  /* output wave data */
  if (WAVE_DATA)
  {
    WaveOutFile();
  }

  /* Loop for Duration at the current wave sign and wave angle */
  for (xx = 0; xx < Duration; xx++)
  {
    MassCurrent = MassCount();
    /* Text to Screen */
    if (current_time_step % SCREEN_TEXT_SPACING == 0)
    {
      printf("==== WaveAngle: %2.2f  MASS Percent: %1.4f  Time Step: %d\n", 180 * (WaveAngle) / PI, MassCurrent / MassInitial, current_time_step);
    }

    periodic_boundary_copy(); /* Copy visible data to external model - creates boundaries */
    ZeroVars();
    /* Initialize for Find Beach Cells  (make sure strange beach does not cause trouble */
    FellOffArray = TRUE;
    FindStart = 1;
    FellOffRockArray = TRUE; /* LMV */
    FindRockStart = 1; /* LMV */

    /* Look for beach - if you fall off of array, bump over a little and try again */
    while (FellOffArray)
    {
      FindBeachCells(FindStart);
      FindStart += FIND_CELL_ERROR;
      if (FellOffArray) 
      {
        // printf("NOODLE  !!!!!FoundCells: %d GetO = %c \n", FindStart,FellOffArray);
        // PauseRun(1,1,-1);
      }

      /* Get Out if no good beach spots exist - finish program */
      if (FindStart > Y_MAX / 2 + 1)
      {
        printf("Stopped Finding Beach - done %d %d", FindStart, Y_MAX / 2 - 5);
        if (SAVE_FILE) {
          SaveSandToFile();
        }
        getchar();
        return 1;
      }
    }

    /* LMV */
    while (FellOffRockArray)
    {
      FindRockCells(FindRockStart);
      // printf("FoundRockCells: %d GetO = %c \n",* FindRockStart,FellOffRockArray);
      FindRockStart += FIND_CELL_ERROR;
      if (FellOffRockArray)
      {
        // printf("NOODLE  !!!!!FoundRockCells: %d GetO = %c \n",* FindRockStart,FellOffRockArray);
        // PauseRun(1,1,-1);
      }

      if (FindRockStart > Y_MAX / 2 + 1)
      {
        printf("Stopped Finding Rock - done %d %d", FindRockStart, Y_MAX / 2 - 5);
        if (SAVE_FILE)
        {
          SaveSandToFile();
        }
        return 1;
      }
    }

    if (debug0)
    {
      printf("going to pause after finding beach, rock\n");
      PauseRun(58, 269, -1);
    }

    if (INITIAL_CONDITION_TYPE != 3) 
    {
      RockCalculations();
    }

    ShadowSweep(); /* LMV */
    DetermineAngles();
    DetermineSedTransport();
    DoSink();     /* clear sand from sink cells */
    FlowInCell(); /* LMV */
    FixFlow(); /* LMV */
    TransportSedimentSweep();
    FixBeach();

    /* Age Empty Cells */
    if ((current_time_step % AGE_UPDATE == 0) && SAVE_AGE)
    {
      AgeCells();
    }

    /* Count Mass */
    MassCurrent = MassCount();

    /* SAVE FILE ? */
    if ((current_time_step % SAVE_SPACING == 0 && current_time_step > START_SAVING_AT))
    {
      if (SAVE_LINE)
      {
        SaveLineToFile();
      }
      if (SAVE_FILE)
      {
        SaveSandToFile();
      }
    }

    current_time_step++;
    current_time += time_step;
  }
  return 0;
}

/**
 * Finish and exit CEM run
 * PARAMETERS: none
 * RETURN: exit status
 */
int cem_finalize(void)
{
  free2d((void**)topography);
  free2d((void**)shelf_depth);
  free2d((void**)wave_h_sig);
  free2d((void**)wave_dir);
  free2d((void**)type_of_rock);

  printf("Run Complete.  Output file: %s \n", savefilename);
  return 0;
}

/*********************************************************
 CEM FUNCTIONS ARE Below
 ********************************************************/

/**
 * Calculates wave angle for given time step
 * PARAMETERS: none
 * RETURN: Wave angle in degrees? Radians?
 */
double FindWaveAngle(void)
{
  /* Data Input Method */
  double RandBin;   /* Random number to pick wave angle bin */
  double RandAngle; /* Random number to pick wave angle within the bin */
  int flag = TRUE;
  int i = 0;

  switch (WAVE_IN)
  {
      /* Method using wave probability step function
       variable Asym will determine fractional distribution of waves coming from the positive direction (left) 
       -i.e. fractional wave asymmetry
       redone by RCL. note Highness actually means lowness, for some reason */
    case 0:
      Angle = RandZeroToOne() * PI / 4.0;
      if (RandZeroToOne() >= Highness)
      {
        Angle += PI / 4.0; /* make high angle */
      }
      if (RandZeroToOne() >= Asym)
      {
        Angle = -1.0 * Angle;
      }
      // printf("WaveAngle: %f degrees\n",Angle*RAD_TO_DEG);
    
    /* Method using input binned wave distribution variables WaveProb[], WaveMax[], 
    previously input from file using ReadWaveIn()*/
    case 1:
      RandBin = RandZeroToOne();
      RandAngle = RandZeroToOne();
      // printf("Time = %d RandBin = %f RandAng = %f\n",current_time_step, RandBin,RandAngle);

      while (flag)
      {
        i++;
        if (RandBin < WaveProb[i])
        {
          flag = FALSE;
        }
      }

      Angle = -(RandAngle * (WaveMax[i] - WaveMax[i - 1]) + WaveMax[i - 1]) * PI / 180;
      // printf("i = %d WaveMAx[i] = %f WaveMax[i-1] = %f WaveProb[i] = %f Angle= %f\n", i,WaveMax[i],WaveMax[i-1],WaveProb[i], Angle*180/PI);
 
    /*Read real data (angle/period/height) */
    case 2:
      RealWaveIn();
      Angle = ((double)(WaveAngleIn)) + coastrotation; /*align waveclimate to coastline in CEM */
      if (Angle >= 180.0)
      {
        Angle -= 360.0; /* conversion to CEM angles (-180 -- 0 -- +180) */
      }
      Angle *= -1; /*-1 is there as wave angle is from 90 in the west to -90 in the east*/

      Period = WavePeriodIn * waveperiodchange;
      OffShoreWvHt = WaveHeightIn * waveheightchange;

      if ((Angle > 90) || (Angle < -90))
      {
        Angle = 0; /*if waves were to come from behind coast */
        OffShoreWvHt = 1e-10; /*if waves were to come from behind coast set to inf small */
      }
      Angle /= RAD_TO_DEG;
  }

  return Angle;
}

// TODO: find better way to find beach cells
// TODO: more descriptive debug statements
/**
 * Determines locations of beach cells moving from left to right direction
 * PARAMETERS: Y position to start searching
 * RETURN: none, writes to X[] and Y[]
 */
void FindBeachCells(int YStart)
{
  int y, z, xstart; /* local iterators */
  /* Starting at left end, find the x - value for first cell that is 'allbeach' */
  xstart = X_MAX - 1;
  y = YStart;

  while (!AllBeach[xstart][y])
  {
    xstart -= 1;
  }

  xstart += 1; /* Step back to where partially full beach */

  X[0] = xstart;
  Y[0] = YStart;

  if (debug1)
  {
    printf("FirstX: %3d  FirstY: %3d  z: 0 \n", X[0], Y[0]);
  }

  z = 0;
  X[z - 1] = X[z];
  Y[z - 1] = Y[z] - 1;

  while ((Y[z] < 2 * Y_MAX - 1) && (z < MAX_BEACH_LENGTH - 1))
  {
    z++;
    NextX = -2;
    NextY = -2;
    BehindX = -2;
    BehindY = -2;

    FindNextCell(X[z - 1], Y[z - 1], z - 1);
    X[z] = NextX;
    Y[z] = NextY;
    XBehind[z] = BehindX;
    YBehind[z] = BehindY;

    if (debug1)
    {
      printf("NextX: %3d  NextY: %3d  z: %d \n", NextX, NextY, z);
      printf("z: %d, BehindX: %d; BehindY: %d; \n\n", z, BehindX, BehindY);
    }

    if (PercentFullSand[X[z]][Y[z]] == 0)
    {
      // printf("\nFINDBEACH: PercentFullSand Zero x: %d y: %d\n",X[z],Y[z]);
      // PauseRun(X[z],Y[z],z);
    }

    /* If return to start point or go off left side of array, going the wrong direction */
    /* Jump off and start again closer to middle */
    if ((NextY < 1) || ((NextY == Y[0]) && (NextX == X[0])) ||  (z > MAX_BEACH_LENGTH - 2))
    {
      // printf("!!!!!!!Fell Off!!!!!!!!!!!!!!!!! x = %d !!!!!!!!!!!!!", NextX);
      FellOffArray = TRUE;
      ZeroVars();
      return;
    }

    if (z > MAX_BEACH_LENGTH - 3)
    {
      printf("????????????  went to end of MaxBeach!! ????");
    }
  }

  TotalBeachCells = z;
  FellOffArray = FALSE;

  if (debug1)
  {
    printf("Total Beach: %d  Current Time Step: %d \n \n", TotalBeachCells,
           current_time_step);
  }
}

// TODO: share common logic with above
/**
 * Determines locations of rock cells moving from left to right direction
 * PARAMETERS: Y position to start searching
 * RETURN: none, writes to XRock[] and YROCK[]
 */
void FindRockCells(int YStart)
{
  int y, z, xstart; /* local iterators */
  /* Starting at left end, find the x - value for first cell that is 'allbeach' */

  xstart = X_MAX - 1;
  y = YStart;

  while (AllRock[xstart][y])
  {
    xstart -= 1;
  }

  xstart += 1; /* Step back to where partially full beach */

  XRock[0] = xstart;
  YRock[0] = YStart;

  if (debug1a)
  {
    printf("FirstRockX: %3d  FirstRockY: %3d  z: 0 \n", XRock[0], YRock[0]);
  }

  z = 0;
  XRock[z - 1] = XRock[z];
  YRock[z - 1] = YRock[z] - 1;

  while ((YRock[z] < 2 * Y_MAX - 1) && (z < MAX_BEACH_LENGTH - 1))
  {
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
    {
      printf("NextRockX: %3d  NextRockY: %3d  z: %d \n", NextRockX, NextRockY, z);
      printf("z: %d, BehindRockX: %d; BehindRockY: %d \n\n", z, BehindRockX, BehindRockY);
    }

    /* RCL takes this out */
    // if (PercentFullRock[XRock[z]][YRock[z]] == 0)
    // {
    //   printf("\nFINDROCK: PercentFullRock Zero x: %d y: %d\n",XRock[z],YRock[z]);
    // }

    /* If return to start point or go off left side of array, going the wrong direction */
    /* Jump off and start again closer to middle */

    if ((NextRockY < 1) || ((NextRockY == YRock[0]) && (NextRockX == XRock[0])) || (z > MAX_BEACH_LENGTH - 2))
    {
      printf("!!!!!!!Fell Off Rock!!!!!!!!!!!!! x = %d !!!!!!!!!!!!!", NextRockX);
      FellOffRockArray = TRUE;
      ZeroVars();
      return;
    }

    if (z > MAX_BEACH_LENGTH - 3)
    {
      printf("????????????  went to end of Rock MaxBeach!! ????");
      printf("(maybe the beach made a loop)\n");
    }
  }

  TotalRockCells = z;
  FellOffRockArray = FALSE;

  if (debug1a)
  {
    printf("Total Rock: %d  Current Time Step: %d \n \n", TotalRockCells, current_time_step);
  }
}

/**
 * Finds the next cell that is beach moving in the general positive x direction
 * PARAMETERS: current x, y, and z
 * RETURN: none, changes NextX and NextY
 */
void FindNextCell(int x, int y, int z)
{
  /* came from left */
  if ((X[z - 1] == X[z]) && (Y[z - 1] == Y[z] - 1))
  {
    if ((PercentFullSand[x + 1][y] + PercentFullRock[x + 1][y]) > 0.0)
    {
      if (debug1b)
      {
         printf("1\n");
      }
      NextX = x + 1;
      NextY = y;
      BehindX = NextX - 1;
      BehindY = NextY;
      return;
    }
    else if ((PercentFullSand[x + 1][y + 1] + PercentFullRock[x + 1][y + 1]) > 0.0)
    {
      if (debug1b)
      {
        printf("2\n");
      }
      NextX = x + 1;
      NextY = y + 1;
      BehindX = NextX - 1;
      BehindY = NextY;
      return;
    }
    else if ((PercentFullSand[x][y + 1] + PercentFullRock[x][y + 1]) > 0.0)
    {
      if (debug1b)
      {
        printf("3\n");
      }
      NextX = x;
      NextY = y + 1;
      BehindX = NextX - 1;
      BehindY = NextY;
      return;
    }
    else if ((PercentFullSand[x - 1][y + 1] + PercentFullRock[x - 1][y + 1]) > 0.0)
    {
      if (debug1b)
      {
        printf("4\n");
      }
      NextX = x - 1;
      NextY = y + 1;
      BehindX = NextX - 1;
      BehindY = NextY;
      return;
    }
    else if ((PercentFullSand[x - 1][y] + PercentFullRock[x - 1][y]) > 0.0)
    {
      if (debug1b)
      {
        printf("5\n");
      }
      NextX = x - 1;
      NextY = y;
      BehindX = NextX - 1;
      BehindY = NextY;
      return;
    }
    else if ((PercentFullSand[x - 1][y - 1] + PercentFullRock[x - 1][y - 1]) > 0.0)
    {
      if (debug1b)
      {
        printf("6\n");
      }
      NextX = x - 1;
      NextY = y - 1;
      BehindX = NextX - 1;
      BehindY = NextY;
      return;
    }
    else
    {
      printf("Should've found next cell (Left): %d, %d \n", x, y);
    }
  }
  /* came from upper left */
  else if ((X[z - 1] == X[z] + 1) && (Y[z - 1] == Y[z] - 1))
  {
    if ((PercentFullSand[x + 1][y + 1] + PercentFullRock[x + 1][y + 1]) > 0.0)
    {
      if (debug1b)
      {
        printf("7\n");
      }
      NextX = x + 1;
      NextY = y + 1;
      BehindX = NextX - 1;
      BehindY = NextY;
      return;
    }
    else if ((PercentFullSand[x][y + 1] + PercentFullRock[x][y + 1]) > 0.0)
    {
      if (debug1b)
      {
        printf("8\n");
      }
      NextX = x;
      NextY = y + 1;
      BehindX = NextX - 1;
      BehindY = NextY;
      return;
    }
    else if ((PercentFullSand[x - 1][y + 1] + PercentFullRock[x - 1][y + 1]) > 0.0)
    {
      if (debug1b)
      {
        printf("9\n");
      }
      NextX = x - 1;
      NextY = y + 1;
      BehindX = NextX - 1;
      BehindY = NextY;
      return;
    }
    else if ((PercentFullSand[x - 1][y] + PercentFullRock[x - 1][y]) > 0.0)
    {
      if (debug1b)
      {
        printf("10\n");
      }
      NextX = x - 1;
      NextY = y;
      BehindX = NextX;
      BehindY = NextY - 1;
      return;
    }
    else if ((PercentFullSand[x - 1][y - 1] + PercentFullRock[x - 1][y - 1]) > 0.0)
    {
      if (debug1b)
      {
        printf("11\n");
      }
      NextX = x - 1;
      NextY = y - 1;
      BehindX = NextX;
      BehindY = NextY - 1;
      return;
    }
    else if ((PercentFullSand[x][y - 1] + PercentFullRock[x][y - 1]) > 0.0)
    {
      if (debug1b)
      {
        printf("12\n");
      }
      NextX = x;
      NextY = y - 1;
      BehindX = NextX;
      BehindY = NextY - 1;
      return;
    }
    else
    {
      printf("Should've found next cell (Upper Left): %d, %d \n", x, y);
      PauseRun(x, y, -1);
    }
  }
  /* came from above */
  else if ((X[z - 1] == X[z] + 1) && (Y[z - 1] == Y[z]))
  {
    if ((PercentFullSand[x][y + 1] + PercentFullRock[x][y + 1]) > 0.0)
    {
      if (debug1b)
      {
        printf("13\n");
      }
      NextX = x;
      NextY = y + 1;
      BehindX = NextX;
      BehindY = NextY - 1;
      return;
    }
    else if ((PercentFullSand[x - 1][y + 1] + PercentFullRock[x - 1][y + 1]) > 0.0)
    {
      if (debug1b)
      {
        printf("14\n");
      }
      NextX = x - 1;
      NextY = y + 1;
      BehindX = NextX;
      BehindY = NextY - 1;
      return;
    }
    else if ((PercentFullSand[x - 1][y] + PercentFullRock[x - 1][y]) > 0.0)
    {
      if (debug1b)
      {
        printf("15\n");
      }
      NextX = x - 1;
      NextY = y;
      BehindX = NextX;
      BehindY = NextY - 1;
      return;
    }
    else if ((PercentFullSand[x - 1][y - 1] + PercentFullRock[x - 1][y - 1]) > 0.0)
    {
      if (debug1b)
      {
        printf("16\n");
      }
      NextX = x - 1;
      NextY = y - 1;
      BehindX = NextX;
      BehindY = NextY - 1;
      return;
    }
    else if ((PercentFullSand[x][y - 1] + PercentFullRock[x][y - 1]) > 0.0)
    {
      if (debug1b)
      {
        printf("17\n");
      }
      NextX = x;
      NextY = y - 1;
      BehindX = NextX;
      BehindY = NextY - 1;
      return;
    }
    else if ((PercentFullSand[x + 1][y - 1] + PercentFullRock[x + 1][y - 1]) > 0.0)
    {
      if (debug1b)
      {
        printf("18\n");
      }
      NextX = x + 1;
      NextY = y - 1;
      BehindX = NextX;
      BehindY = NextY - 1;
      return;
    }
    else if ((PercentFullSand[x + 1][y - 1] + PercentFullRock[x + 1][y - 1]) <1.0)
    { 
      /* RCL adds ad hoc case */
      if (debug1b)
      {
        printf("18.5\n");
      }
      NextX = x + 1;
      NextY = y - 1;
      BehindX = NextX;
      BehindY = NextY - 1;
      return;
    }
    else
    {
      printf("Should've found next cell (Above): %d, %d \n", x, y);
      PauseRun(x, y, -1);
    }
  }
  /* came from upper right */
  else if ((X[z - 1] == X[z] + 1) && (Y[z - 1] == Y[z] + 1))
  {
    if ((PercentFullSand[x - 1][y + 1] + PercentFullRock[x - 1][y + 1]) > 0.0)
    {
      if (debug1b)
      {
        printf("19\n");
      }
      NextX = x - 1;
      NextY = y + 1;
      BehindX = NextX;
      BehindY = NextY - 1;
      return;
    }
    else if ((PercentFullSand[x - 1][y] + PercentFullRock[x - 1][y]) > 0.0)
    {
      if (debug1b)
      {
        printf("20\n");
      }
      NextX = x - 1;
      NextY = y;
      BehindX = NextX;
      BehindY = NextY - 1;
      return;
    }
    else if ((PercentFullSand[x - 1][y - 1] + PercentFullRock[x - 1][y - 1]) > 0.0)
    {
      if (debug1b)
      {
        printf("21\n");
      }
      NextX = x - 1;
      NextY = y - 1;
      BehindX = NextX + 1;
      BehindY = NextY;
      return;
    }
    else if ((PercentFullSand[x][y - 1] + PercentFullRock[x][y - 1]) > 0.0)
    {
      if (debug1b)
      {
        printf("22\n");
      }
      NextX = x;
      NextY = y - 1;
      BehindX = NextX + 1;
      BehindY = NextY;
      return;
    }
    else if ((PercentFullSand[x + 1][y - 1] + PercentFullRock[x + 1][y - 1]) > 0.0)
    {
      if (debug1b)
      {
        printf("23\n");
      }
      NextX = x + 1;
      NextY = y - 1;
      BehindX = NextX + 1;
      BehindY = NextY;
      return;
    }
    else if ((PercentFullSand[x + 1][y] + PercentFullRock[x + 1][y]) > 0.0)
    {
      if (debug1b)
      {
        printf("24\n");
      }
      NextX = x + 1;
      NextY = y;
      BehindX = NextX + 1;
      BehindY = NextY;
      return;
    }
    else
    {
      printf("Should've found next cell (Upper Right): %d, %d \n", x, y);
      PauseRun(x, y, -1);
    }
  }
  /* came from right */
  else if ((X[z - 1] == X[z]) && (Y[z - 1] == Y[z] + 1))
  {
    if ((PercentFullSand[x - 1][y] + PercentFullRock[x - 1][y]) > 0.0)
    {
      if (debug1b)
      {
        printf("25\n");
      }
      NextX = x - 1;
      NextY = y;
      BehindX = NextX + 1;
      BehindY = NextY;
      return;
    }
    else if ((PercentFullSand[x - 1][y - 1] + PercentFullRock[x - 1][y - 1]) > 0.0)
    {
      if (debug1b)
      {
        printf("26\n");
      }
      NextX = x - 1;
      NextY = y - 1;
      BehindX = NextX + 1;
      BehindY = NextY;
      return;
    }
    else if ((PercentFullSand[x][y - 1] + PercentFullRock[x][y - 1]) > 0.0)
    {
      if (debug1b)
      {
        printf("27\n");
      }
      NextX = x;
      NextY = y - 1;
      BehindX = NextX + 1;
      BehindY = NextY;
      return;
    }
    else if ((PercentFullSand[x + 1][y - 1] + PercentFullRock[x + 1][y - 1]) > 0.0)
    {
      if (debug1b)
      {
        printf("28\n");
      }
      NextX = x + 1;
      NextY = y - 1;
      BehindX = NextX + 1;
      BehindY = NextY;
      return;
    }
    else if ((PercentFullSand[x + 1][y] + PercentFullRock[x + 1][y]) > 0.0)
    {
      if (debug1b)
      {
        printf("29\n");
      }
      NextX = x + 1;
      NextY = y;
      BehindX = NextX + 1;
      BehindY = NextY;
      return;
    }
    else if ((PercentFullSand[x + 1][y + 1] + PercentFullRock[x + 1][y + 1]) > 0.0)
    {
      if (debug1b)
      {
        printf("30\n");
      }
      NextX = x + 1;
      NextY = y + 1;
      BehindX = NextX + 1;
      BehindY = NextY;
      return;
    }
    else
    {
      printf("Should've found next cell (Right): %d, %d \n", x, y);
      PauseRun(x, y, -1);
    }
  }
  /* came from lower right */
  else if ((X[z - 1] == X[z] - 1) && (Y[z - 1] == Y[z] + 1))
  {
    if ((PercentFullSand[x - 1][y - 1] + PercentFullRock[x - 1][y - 1]) > 0.0)
    {
      if (debug1b)
      {
        printf("31\n");
      }
      NextX = x - 1;
      NextY = y - 1;
      BehindX = NextX + 1;
      BehindY = NextY;
      return;
    }
    else if ((PercentFullSand[x][y - 1] + PercentFullRock[x][y - 1]) > 0.0)
    {
      if (debug1b)
      {
        printf("32\n");
      }
      NextX = x;
      NextY = y - 1;
      BehindX = NextX + 1;
      BehindY = NextY;
      return;
    }
    else if ((PercentFullSand[x + 1][y - 1] + PercentFullRock[x + 1][y - 1]) > 0.0)
    {
      if (debug1b)
      {
        printf("33\n");
      }
      NextX = x + 1;
      NextY = y - 1;
      BehindX = NextX + 1;
      BehindY = NextY;
      return;
    }
    else if ((PercentFullSand[x + 1][y] + PercentFullRock[x + 1][y]) > 0.0)
    {
      if (debug1b)
      {
        printf("34\n");
      }
      NextX = x + 1;
      NextY = y;
      BehindX = NextX;
      BehindY = NextY + 1;
      return;
    }
    else if ((PercentFullSand[x + 1][y + 1] + PercentFullRock[x + 1][y + 1]) > 0.0)
    {
      if (debug1b)
      {
        printf("35\n");
      }
      NextX = x + 1;
      NextY = y + 1;
      BehindX = NextX;
      BehindY = NextY + 1;
      return;
    }
    else if ((PercentFullSand[x][y + 1] + PercentFullRock[x][y + 1]) > 0.0)
    {
      if (debug1b)
      {
        printf("36\n");
      }
      NextX = x;
      NextY = y + 1;
      BehindX = NextX;
      BehindY = NextY + 1;
      return;
    }
    /* RCL adds new ad hoc case */
    else if (PercentFullSand[x][y + 2] + PercentFullRock[x][y + 2] > 0.0 ||
               PercentFullSand[x - 1][y + 2] + PercentFullRock[x - 1][y + 2] > 0.0 ||
               PercentFullSand[x + 1][y + 2] + PercentFullRock[x + 1][y + 2] > 0.0)
      {
      if (debug1b)
      {
        printf("36.5\n");
      }
      NextX = x;
      NextY = y + 1;
      BehindX = NextX;
      BehindY = NextY + 1;
      return;
    }
    else
    {
      printf("Should've found next cell (Lower Right): %d, %d \n", x, y);
      PauseRun(x, y, -1);
    }
  }
  /* came from below */
  else if ((X[z - 1] == X[z] - 1) && (Y[z - 1] == Y[z]))
  {
    if ((PercentFullSand[x][y - 1] + PercentFullRock[x][y - 1]) > 0.0)
    {
      if (debug1b)
      {
        printf("37\n");
      }
      NextX = x;
      NextY = y - 1;
      BehindX = NextX;
      BehindY = NextY + 1;
      return;
    }
    else if ((PercentFullSand[x + 1][y - 1] + PercentFullRock[x + 1][y - 1]) > 0.0)
    {
      if (debug1b)
      {
        printf("38\n");
      }
      NextX = x + 1;
      NextY = y - 1;
      BehindX = NextX;
      BehindY = NextY + 1;
      return;
    }
    else if ((PercentFullSand[x + 1][y] + PercentFullRock[x + 1][y]) > 0.0)
    {
      if (debug1b)
      {
        printf("39\n");
      }
      NextX = x + 1;
      NextY = y;
      BehindX = NextX;
      BehindY = NextY + 1;
      return;
    }
    else if ((PercentFullSand[x + 1][y + 1] + PercentFullRock[x + 1][y + 1]) > 0.0)
    {
      if (debug1b)
      {
        printf("40\n");
      }
      NextX = x + 1;
      NextY = y + 1;
      BehindX = NextX;
      BehindY = NextY + 1;
      return;
    }
    /* RCL changed behind */
    else if ((PercentFullSand[x][y + 1] + PercentFullRock[x][y + 1]) > 0.0)
    {
      if (debug1b)
      {
        printf("41\n");
      }
      NextX = x;
      NextY = y + 1;
      BehindX = NextX - 1;
      BehindY = NextY;
      return;
    }
    /* RCL messed with all the 42.X cases */
    else if ((PercentFullSand[x - 1][y + 1]) + PercentFullRock[x - 1][y + 1] > 0.0)
    {
      if (debug1b)
      {
        printf("42\n");
      }
      NextX = x - 1;
      NextY = y + 1;
      BehindX = NextX - 1;
      BehindY = NextY;
      return;
    }
    else if ((PercentFullSand[x - 1][y - 1]) + PercentFullRock[x - 1][y - 1] > 0.0)
    {
      if (debug1b)
      {
        printf("42.1\n");
      }
      NextX = x - 1;
      NextY = y - 1;
      BehindX = NextX + 1;
      BehindY = NextY - 1;
      return;
    }
    /* RCL added this case */
    else if ((PercentFullSand[x][y + 1]) + PercentFullRock[x][y + 1] < 1.0)
    {
      if (debug1b)
      {
        printf("42.75\n");
      }
      NextX = x;
      NextY = y + 1;
      BehindX = NextX - 1;
      BehindY = NextY;
      return;
    }
    /* RCL added this case */
    else if ((PercentFullSand[x - 1][y + 1]) + PercentFullRock[x - 1][y + 1] < 1.0)
    {
      if (debug1b)
      {
        printf("42.8\n");
      }
      NextX = x - 1;
      NextY = y + 1;
      BehindX = NextX - 1;
      BehindY = NextY;
      return;
    }
    else {
      printf("Should've found next cell (Below): %d, %d \n", x, y);
      PauseRun(x, y, -1);
    }
  }
  /* came from lower left */
  else if ((X[z - 1] == X[z] - 1) && (Y[z - 1] == Y[z] - 1))
  {
    if ((PercentFullSand[x + 1][y - 1] + PercentFullRock[x + 1][y - 1]) > 0.0)
    {
      if (debug1b)
      {
        printf("43\n");
      }
      NextX = x + 1;
      NextY = y - 1;
      BehindX = NextX;
      BehindY = NextY + 1;
      return;
    }
    else if ((PercentFullSand[x + 1][y] + PercentFullRock[x + 1][y]) > 0.0)
    {
      if (debug1b)
      {
        printf("44\n");
      }
      NextX = x + 1;
      NextY = y;
      BehindX = NextX;
      BehindY = NextY + 1;
      return;
    }
    else if ((PercentFullSand[x + 1][y + 1] + PercentFullRock[x + 1][y + 1]) > 0.0)
    {
      if (debug1b)
      {
        printf("45\n");
      }
      NextX = x + 1;
      NextY = y + 1;
      BehindX = NextX - 1;
      BehindY = NextY;
      return;
    }
    else if ((PercentFullSand[x][y + 1] + PercentFullRock[x][y + 1]) > 0.0)
    {
      if (debug1b)
      {
        printf("46\n");
      }
      NextX = x;
      NextY = y + 1;
      BehindX = NextX - 1;
      BehindY = NextY;
      return;
    }
    else if ((PercentFullSand[x - 1][y + 1] + PercentFullRock[x - 1][y + 1]) > 0.0)
    {
      if (debug1b)
      {
        printf("47\n");
      }
      NextX = x - 1;
      NextY = y + 1;
      BehindX = NextX - 1;
      BehindY = NextY;
      return;
    }
    else if ((PercentFullSand[x - 1][y] + PercentFullRock[x - 1][y]) > 0.0)
    {
      if (debug1b)
      {
        printf("48\n");
      }
      NextX = x - 1;
      NextY = y;
      BehindX = NextX - 1;
      BehindY = NextY;
      return;
    }
    else
    {
      printf("Should've found next cell (Lower Left): %d, %d \n", x, y);
      PauseRun(x, y, -1);
    }
  }
  printf("Should've found next cell : %d, %d \n", x, y);
}

// TODO: pull out common code
/**
 * Finds the next cell that is rock moving in the general positive x direction
 * PARAMETERS: current x, y, and z
 * RETURN: none, changes NextRockX and NextRockY
 */
void FindNextRockCell(int x, int y, int z)
{
  /* came from left */
  if ((XRock[z - 1] == XRock[z]) && (YRock[z - 1] == YRock[z] - 1))
  {
    if (PercentFullRock[x + 1][y] > 0.0)
    {
      if (debug1c)
      {
        printf("1 r\n");
      }
      NextRockX = x + 1;
      NextRockY = y;
      BehindRockX = NextRockX - 1;
      BehindRockY = NextRockY;
      return;
    }
    else if (PercentFullRock[x + 1][y + 1] > 0.0)
    {
      if (debug1c)
      {
        printf("2 r\n");
      }
      NextRockX = x + 1;
      NextRockY = y + 1;
      BehindRockX = NextRockX - 1;
      BehindRockY = NextRockY;
      return;
    }
    else if (PercentFullRock[x][y + 1] > 0.0)
    {
      if (debug1c)
      {
        printf("3 r\n");
      }
      NextRockX = x;
      NextRockY = y + 1;
      BehindRockX = NextRockX - 1;
      BehindRockY = NextRockY;
      return;
    }
    else if (PercentFullRock[x - 1][y + 1] > 0.0)
    {
      if (debug1c)
      {
        printf("4 r\n");
      }
      NextRockX = x - 1;
      NextRockY = y + 1;
      BehindRockX = NextRockX - 1;
      BehindRockY = NextRockY;
      return;
    }
    else if (PercentFullRock[x - 1][y] > 0.0)
    {
      if (debug1c)
      {
        printf("5 r\n");
      }
      NextRockX = x - 1;
      NextRockY = y;
      BehindRockX = NextRockX - 1;
      BehindRockY = NextRockY;
      return;
    }
    else if (PercentFullRock[x - 1][y - 1] > 0.0)
    {
      if (debug1c) {
        printf("6 r\n");
      }
      NextRockX = x - 1;
      NextRockY = y - 1;
      BehindRockX = NextRockX - 1;
      BehindRockY = NextRockY;
      return;
    }
    else
    {
      printf("Should've found next Rock cell (Left): %d, %d \n", x, y);
    }
  }
  /* came from upper left */
  else if ((XRock[z - 1] == XRock[z] + 1) && (YRock[z - 1] == YRock[z] - 1))
  {
    if (PercentFullRock[x + 1][y + 1] > 0.0)
    {
      if (debug1c)
      {
        printf("7 r\n");
      }
      NextRockX = x + 1;
      NextRockY = y + 1;
      BehindRockX = NextRockX - 1;
      BehindRockY = NextRockY;
      return;
    }
    else if (PercentFullRock[x][y + 1] > 0.0)
    {
      if (debug1c)
      {
        printf("8 r\n");
      }
      NextRockX = x;
      NextRockY = y + 1;
      BehindRockX = NextRockX - 1;
      BehindRockY = NextRockY;
      return;
    }
    else if (PercentFullRock[x - 1][y + 1] > 0.0)
    {
      if (debug1c)
      {
        printf("9 r\n");
      }
      NextRockX = x - 1;
      NextRockY = y + 1;
      BehindRockX = NextRockX - 1;
      BehindRockY = NextRockY;
      return;
    }
    else if (PercentFullRock[x - 1][y] > 0.0)
    {
      if (debug1c)
      {
        printf("10 r\n");
      }
      NextRockX = x - 1;
      NextRockY = y;
      BehindRockX = NextRockX;
      BehindRockY = NextRockY - 1;
      return;
    }
    else if (PercentFullRock[x - 1][y - 1] > 0.0)
    {
      if (debug1c) {
        printf("11 r\n");
      }
      NextRockX = x - 1;
      NextRockY = y - 1;
      BehindRockX = NextRockX;
      BehindRockY = NextRockY - 1;
      return;
    }
    else if (PercentFullRock[x][y - 1] > 0.0)
    {
      if (debug1c)
      {
        printf("12 r\n");
      }
      NextRockX = x;
      NextRockY = y - 1;
      BehindRockX = NextRockX;
      BehindRockY = NextRockY - 1;
      return;
    }
    else
    {
      printf("Should've found next Rock cell (Upper Left): %d, %d \n", x, y);
    }
  }
  /* came from above */
  else if ((XRock[z - 1] == XRock[z] + 1) && (YRock[z - 1] == YRock[z]))
  {
    if (PercentFullRock[x][y + 1] > 0.0)
    {
      if (debug1c)
      {
        printf("13 r\n");
      }
      NextRockX = x;
      NextRockY = y + 1;
      BehindRockX = NextRockX;
      BehindRockY = NextRockY - 1;
      return;
    }
    else if (PercentFullRock[x - 1][y + 1] > 0.0)
    {
      if (debug1c)
      {
        printf("14 r\n");
      }
      NextRockX = x - 1;
      NextRockY = y + 1;
      BehindRockX = NextRockX;
      BehindRockY = NextRockY - 1;
      return;
    }
    else if (PercentFullRock[x - 1][y] > 0.0)
    {
      if (debug1c)
      {
        printf("15 r\n");
      }
      NextRockX = x - 1;
      NextRockY = y;
      BehindRockX = NextRockX;
      BehindRockY = NextRockY - 1;
      return;
    }
    else if (PercentFullRock[x - 1][y - 1] > 0.0)
    {
      if (debug1c)
      {
        printf("16 r\n");
      }
      NextRockX = x - 1;
      NextRockY = y - 1;
      BehindRockX = NextRockX;
      BehindRockY = NextRockY - 1;
      return;
    }
    else if (PercentFullRock[x][y - 1] > 0.0)
    {
      if (debug1c)
      {
        printf("17 r\n");
      }
      NextRockX = x;
      NextRockY = y - 1;
      BehindRockX = NextRockX;
      BehindRockY = NextRockY - 1;
      return;
    }
    else if (PercentFullRock[x + 1][y - 1] > 0.0)
    {
      if (debug1c)
      {
        printf("18 r\n");
      }
      NextRockX = x + 1;
      NextRockY = y - 1;
      BehindRockX = NextRockX;
      BehindRockY = NextRockY - 1;
      return;
    }
    else
    {
      printf("Should've found next Rock cell (Above): %d, %d \n", x, y);
    }
  }
  /* came from upper right */
  else if ((XRock[z - 1] == XRock[z] + 1) && (YRock[z - 1] == YRock[z] + 1))
  {
    if (PercentFullRock[x - 1][y + 1] > 0.0)
    {
      if (debug1c)
      {
        printf("19 r\n");
      }
      NextRockX = x - 1;
      NextRockY = y + 1;
      BehindRockX = NextRockX;
      BehindRockY = NextRockY - 1;
      return;
    }
    else if (PercentFullRock[x - 1][y] > 0.0)
    {
      if (debug1c)
      {
        printf("20 r\n");
      }
      NextRockX = x - 1;
      NextRockY = y;
      BehindRockX = NextRockX;
      BehindRockY = NextRockY - 1;
      return;
    }
    else if (PercentFullRock[x - 1][y - 1] > 0.0)
    {
      if (debug1c)
      {
        printf("21 r\n");
      }
      NextRockX = x - 1;
      NextRockY = y - 1;
      BehindRockX = NextRockX + 1;
      BehindRockY = NextRockY;
      return;
    }
    else if (PercentFullRock[x][y - 1] > 0.0)
    {
      if (debug1c)
      {
        printf("22 r\n");
      }
      NextRockX = x;
      NextRockY = y - 1;
      BehindRockX = NextRockX + 1;
      BehindRockY = NextRockY;
      return;
    }
    else if (PercentFullRock[x + 1][y - 1] > 0.0)
    {
      if (debug1c)
      {
        printf("23 r\n");
      }
      NextRockX = x + 1;
      NextRockY = y - 1;
      BehindRockX = NextRockX + 1;
      BehindRockY = NextRockY;
      return;
    }
    else if (PercentFullRock[x + 1][y] > 0.0)
    {
      if (debug1c)
      {
        printf("24 r\n");
      }
      NextRockX = x + 1;
      NextRockY = y;
      BehindRockX = NextRockX + 1;
      BehindRockY = NextRockY;
      return;
    }
    else
    {
      printf("Should've found next Rock cell (Upper Right): %d, %d \n", x, y);
    }
  }
  /* came from right */
  else if ((XRock[z - 1] == XRock[z]) && (YRock[z - 1] == YRock[z] + 1))
  {
    if (PercentFullRock[x - 1][y] > 0.0) 
    {
      if (debug1c)
      {
        printf("25 r\n");
      }
      NextRockX = x - 1;
      NextRockY = y;
      BehindRockX = NextRockX + 1;
      BehindRockY = NextRockY;
      return;
    }
    else if (PercentFullRock[x - 1][y - 1] > 0.0)
    {
      if (debug1c)
      {
        printf("26 r\n");
      }
      NextRockX = x - 1;
      NextRockY = y - 1;
      BehindRockX = NextRockX + 1;
      BehindRockY = NextRockY;
      return;
    }
    else if (PercentFullRock[x][y - 1] > 0.0)
    {
      if (debug1c)
      {
        printf("27 r\n");
      }
      NextRockX = x;
      NextRockY = y - 1;
      BehindRockX = NextRockX + 1;
      BehindRockY = NextRockY;
      return;
    }
    else if (PercentFullRock[x + 1][y - 1] > 0.0)
    {
      if (debug1c)
      {
        printf("28 r\n");
      }
      NextRockX = x + 1;
      NextRockY = y - 1;
      BehindRockX = NextRockX + 1;
      BehindRockY = NextRockY;
      return;
    }
    else if (PercentFullRock[x + 1][y] > 0.0)
    {
      if (debug1c)
      {
        printf("29 r\n");
      }
      NextRockX = x + 1;
      NextRockY = y;
      BehindRockX = NextRockX + 1;
      BehindRockY = NextRockY;
      return;
    }
    else if (PercentFullRock[x + 1][y + 1] > 0.0)
    {
      if (debug1c)
      {
        printf("30 r\n");
      }
      NextRockX = x + 1;
      NextRockY = y + 1;
      BehindRockX = NextRockX + 1;
      BehindRockY = NextRockY;
      return;
    }
    else
    {
      printf("Should've found next Rock cell (Right): %d, %d \n", x, y);
    }
  }
  /* came from lower right */
  else if ((XRock[z - 1] == XRock[z] - 1) && (YRock[z - 1] == YRock[z] + 1))
  {
    if (PercentFullRock[x - 1][y - 1] > 0.0)
    {
      if (debug1c)
      {
        printf("31 r\n");
      }
      NextRockX = x - 1;
      NextRockY = y - 1;
      BehindRockX = NextRockX + 1;
      BehindRockY = NextRockY;
      return;
    }
    else if (PercentFullRock[x][y - 1] > 0.0)
    {
      if (debug1c) {
        printf("32 r\n");
      }
      NextRockX = x;
      NextRockY = y - 1;
      BehindRockX = NextRockX + 1;
      BehindRockY = NextRockY;
      return;
    }
    else if (PercentFullRock[x + 1][y - 1] > 0.0)
    {
      if (debug1c)
      {
        printf("33 r\n");
      }
      NextRockX = x + 1;
      NextRockY = y - 1;
      BehindRockX = NextRockX + 1;
      BehindRockY = NextRockY;
      return;
    }
    else if (PercentFullRock[x + 1][y] > 0.0)
    {
      if (debug1c)
      {
        printf("34 r\n");
      }
      NextRockX = x + 1;
      NextRockY = y;
      BehindRockX = NextRockX;
      BehindRockY = NextRockY + 1;
      return;
    }
    else if (PercentFullRock[x + 1][y + 1] > 0.0)
    {
      if (debug1c)
      {
        printf("35 r\n");
      }
      NextRockX = x + 1;
      NextRockY = y + 1;
      BehindRockX = NextRockX;
      BehindRockY = NextRockY + 1;
      return;
    }
    else if (PercentFullRock[x][y + 1] > 0.0)
    {
      if (debug1c)
      {
        printf("36 r\n");
      }
      NextRockX = x;
      NextRockY = y + 1;
      BehindRockX = NextRockX;
      BehindRockY = NextRockY + 1;
      return;
    }
    else
    {
      printf("Should've found next Rock cell (Lower Right): %d, %d \n", x, y);
    }
  }
  /* came from below */
  else if ((XRock[z - 1] == XRock[z] - 1) && (YRock[z - 1] == YRock[z]))
  {
    if (PercentFullRock[x][y - 1] > 0.0)
    {
      if (debug1c)
      {
        printf("37 r\n");
      }
      NextRockX = x;
      NextRockY = y - 1;
      BehindRockX = NextRockX;
      BehindRockY = NextRockY + 1;
      return;
    }
    else if (PercentFullRock[x + 1][y - 1] > 0.0)
    {
      if (debug1c)
      {
        printf("38 r\n");
      }
      NextRockX = x + 1;
      NextRockY = y - 1;
      BehindRockX = NextRockX;
      BehindRockY = NextRockY + 1;
      return;
    }
    else if (PercentFullRock[x + 1][y] > 0.0)
    {
      if (debug1c)
      {
        printf("39 r\n");
      }
      NextRockX = x + 1;
      NextRockY = y;
      BehindRockX = NextRockX;
      BehindRockY = NextRockY + 1;
      return;
    }
    else if (PercentFullRock[x + 1][y + 1] > 0.0)
    {
      if (debug1c)
      {
        printf("40 r\n");
      }
      NextRockX = x + 1;
      NextRockY = y + 1;
      BehindRockX = NextRockX;
      BehindRockY = NextRockY + 1;
      return;
    }
    /* RCL messed with 41, all the 42.X cases */
    else if (PercentFullRock[x][y + 1] > 0.0)
    {
      if (debug1b)
      {
        printf("41 r\n");
      }
      NextRockX = x;
      NextRockY = y + 1;
      BehindRockX = NextRockX - 1;
      BehindRockY = NextRockY;
      return;
    }
    else if (PercentFullRock[x - 1][y + 1] > 0.0)
    {
      if (debug1b)
      {
        printf("42 r\n");
      }
      NextRockX = x - 1;
      NextRockY = y + 1;
      BehindRockX = NextRockX - 1;
      BehindRockY = NextRockY;
      return;
    }
    else if (PercentFullRock[x - 1][y - 1] > 0.0)
    {
      if (debug1b)
      {
        printf("42.1 r\n");
      }
      NextRockX = x - 1;
      NextRockY = y - 1;
      BehindRockX = NextRockX + 1;
      BehindRockY = NextRockY - 1;
      return;
    }
    /* RCL added this case */
    else if (PercentFullRock[x][y + 1] < 1.0)
    {
      if (debug1b)
      {
        printf("42.75 r\n");
      }
      NextRockX = x;
      NextRockY = y + 1;
      BehindRockX = NextRockX - 1;
      BehindRockY = NextRockY;
      return;
    }
    /* RCL added this case */
    else if (PercentFullRock[x - 1][y + 1] < 1.0)
    {
      if (debug1b)
      {
        printf("42.8 r\n");
      }
      NextRockX = x - 1;
      NextRockY = y + 1;
      BehindRockX = NextRockX - 1;
      BehindRockY = NextRockY;
      return;
    }
    else
    {
      printf("Should've found next rock cell (Below): %d, %d \n", x, y);
      PauseRun(x, y, -1);
    }
  }
  /* came from lower left */
  else if ((XRock[z - 1] == XRock[z] - 1) && (YRock[z - 1] == YRock[z] - 1))
  {
    if (PercentFullRock[x + 1][y - 1] > 0.0)
    {
      if (debug1c)
      {
        printf("43 r\n");
      }
      NextRockX = x + 1;
      NextRockY = y - 1;
      BehindRockX = NextRockX;
      BehindRockY = NextRockY + 1;
      return;
    }
    else if (PercentFullRock[x + 1][y] > 0.0)
    {
      if (debug1c)
      {
        printf("44 r\n");
      }
      NextRockX = x + 1;
      NextRockY = y;
      BehindRockX = NextRockX;
      BehindRockY = NextRockY + 1;
      return;
    }
    else if (PercentFullRock[x + 1][y + 1] > 0.0)
    {
      if (debug1c)
      {
        printf("45 r\n");
      }
      NextRockX = x + 1;
      NextRockY = y + 1;
      BehindRockX = NextRockX - 1;
      BehindRockY = NextRockY;
      return;
    }
    else if (PercentFullRock[x][y + 1] > 0.0)
    {
      if (debug1c)
      {
        printf("46 r\n");
      }
      NextRockX = x;
      NextRockY = y + 1;
      BehindRockX = NextRockX - 1;
      BehindRockY = NextRockY;
      return;
    }
    else if (PercentFullRock[x - 1][y + 1] > 0.0)
    {
      if (debug1c)
      {
        printf("47 r\n");
      }
      NextRockX = x - 1;
      NextRockY = y + 1;
      BehindRockX = NextRockX - 1;
      BehindRockY = NextRockY;
      return;
    }
    else if (PercentFullRock[x - 1][y] > 0.0)
    {
      if (debug1c)
      {
        printf("48 r\n");
      }
      NextRockX = x - 1;
      NextRockY = y;
      BehindRockX = NextRockX - 1;
      BehindRockY = NextRockY;
      return; 
    }
    else
    {
      printf("Should've found next Rock cell (Lower Left): %d, %d \n", x, y);
    }
  }
  printf("Should've found next Rock cell : %d, %d \n", x, y);
}

/* LMV */
/* For each Rock cell j, moves along the beach and determines distance from rock
   to beach */
/* and amount of rock to weather */
/* This function calls FindNearestBeach and WeatherRock */
/* This function will use but not adjust the variable TotalRockCells */
/* Determines the total amount of rock weathered in one time step */
void RockCalculations(void)
{
  int j;
  double TotalAmountWeathered = 0.0;

  for (j = 0; j <= TotalRockCells; j++)
  {
    FindNearestBeach(j);
    WeatherRock(j);
    TotalAmountWeathered += AmountWeathered[j];
  }
  // printf("current_time_step = %d\n", current_time_step)
}

/**
 * Calculates the Euclidean distance between Rock and Beach for entire beach array
 * Finds the minimum beach to rock distance (LMV)
 * PARAMETERS: j = edge of rock in the Y direction
 * RETURN: none, changes X[], Y[], XRock[], and YRock[]
 */
void FindNearestBeach(int j)
{
  int i;
  double DeltaX, DeltaY;
  int LookStart, LookEnd; /* for short rock to beach sweep */

  MinDistToBeach[j] = BIG_DISTANCE_TO_BEACH;

  /* do a full sweep every now and then */
  if (current_time_step % TIME_TO_SWEEP_FULL_BEACH == 0)
  {
    /*full sweep */
    for (i = 0; i <= TotalBeachCells; i++)
    {
      /*vertical distance */
      if (X[i] != XRock[j])
      {                
        DeltaX = ((X[i] - XRock[j] - 1) + (PercentFullSand[X[i]][Y[i]] + PercentFullSand[XRock[j]][YRock[j]])) *CELL_WIDTH;
        if (debug10)
        {
          printf("Rock cell = %d, Beach cell = %d, DeltaX (1)= %f \n", XRock[j], X[i], DeltaX);
        }
      }

      else
      {
        DeltaX = (PercentFullSand[X[i]][Y[i]]) * CELL_WIDTH;
        if (debug10)
        {
          printf("Rock cell = %d, Beach cell = %d, DeltaX (2)= %f \n", XRock[j], X[i], DeltaX);
        }
      }
      /*horizontal distance */
      if (YRock[j] <= Y[i])
      {
        DeltaY = (Y[i] - YRock[j]) * CELL_WIDTH;
        if (debug10)
        {
          printf("DeltaY (1) = %f \n", DeltaY);
        }
      }
      else
      {
        DeltaY = (YRock[j] - Y[i]) * CELL_WIDTH;
        if (debug10)
        {
          printf("DeltaY (2) = %f \n", DeltaY);
        }
      }

      DistanceToBeach[j] = sqrt(DeltaX * DeltaX + DeltaY * DeltaY); /*RCL removed Raise w/neg base */
      if (debug10)
      {
        printf("Distance to beach from XRock[j] %d to X[i] %d is %f \n", XRock[j], X[i], DistanceToBeach[j]);
      }
      if (DistanceToBeach[j] < MinDistToBeach[j])
      {
        MinDistToBeach[j] = DistanceToBeach[j];
        ClosestBeach[j] = i;
      }
    }
    if (debug11)
    {
      printf("Min Distance for j = %d is %f \n", j, MinDistToBeach[j]);
      printf("XRock[j] %d to X[i] %d \n", XRock[j], X[i]);
      printf("Closest Beach is at cell i = %d \n \n", ClosestBeach[j]);
    }
    if (Y[i] == 150)
    {
      printf("Min Distance for j = %d is %f \n", j, MinDistToBeach[j]);
      printf("XRock[j] %d to X[i] %d \n", XRock[j], X[i]);
      printf("Closest Beach is at cell i = %d \n \n", ClosestBeach[j]);
    }
  }

  else
  {
    /* LOOK_DIST allows the program to skip the full beach sweep at every iteration
      by looking at the location of the previous closest beach and sweeping out a small arc
      plus or minus LOOK_DIST */

    /* Only sweep right when on left end */
    LookStart = ClosestBeach[j] - LOOK_DIST >= 0 ? ClosestBeach[j] - LOOK_DIST : 0;
    /* Only sweep left when at right end */
    LookEnd = ClosestBeach[j] + LOOK_DIST <= TotalBeachCells ? ClosestBeach[j] + LOOK_DIST : TotalBeachCells;

    for (i = LookStart; i < LookEnd; i++)
    {
      if (X[i] != XRock[j])
      {
        DeltaX = ((X[i] - XRock[j] - 1) + (PercentFullSand[X[i]][Y[i]] + PercentFullSand[XRock[j]][YRock[j]])) * CELL_WIDTH;
        if (debug10)
        {
          printf("Rock cell = %d, Beach cell = %d, DeltaX (1)= %f \n", XRock[j], X[i], DeltaX);
        }
      }
      else
      {
        DeltaX = (PercentFullSand[X[i]][Y[i]]) * CELL_WIDTH;
        if (debug10)
        {
          printf("Rock cell = %d, Beach cell = %d, DeltaX (2)= %f \n", XRock[j], X[i], DeltaX);
        }
      }
      /*horizontal distance */
      if (YRock[j] <= Y[i])
      {
        DeltaY = (Y[i] - YRock[j]) * CELL_WIDTH;
        if (debug10)
        {
          printf("DeltaY (1) = %f \n", DeltaY);
        }
      }
      else{
        DeltaY = (YRock[j] - Y[i]) * CELL_WIDTH;
        if (debug10)
        {
          printf("DeltaY (2) = %f \n", DeltaY);
        }
      }

      DistanceToBeach[j] = sqrt(DeltaX * DeltaX + DeltaY * DeltaY); /*RCL removed Raise w/neg base */
      if (debug10)
      {
        printf("Distance to beach from XRock[j] %d to X[i] %d is %f \n", XRock[j], X[i], DistanceToBeach[j]);
      }
      if (DistanceToBeach[j] < MinDistToBeach[j])
      {
        MinDistToBeach[j] = DistanceToBeach[j];
        ClosestBeach[j] = i;
      }
    }
    if (debug11)
    {
      printf("Short Sweep - Min Distance for j = %d is %f \n", j, MinDistToBeach[j]);
      printf("XRock[j] %d to X[i] %d \n", XRock[j], X[i]);
      printf("Closest Beach is at cell i = %d \n ", ClosestBeach[j]);
      printf("X[i]: %d, Y[i]: %d \n\n", X[i], Y[i]);
    }
  }
}

/**
 * Simulates rock weathering at each time step using fast and slow weathering rates
 * At this point, weathering rate is pretty arbitrary (LMV)
 * PARAMETERS: j = edge of rock in the Y direction
 * RETURN: none, changes PercentFullRock[][], PercentFullSand[][], AllRock[][]
 */
void WeatherRock(int j)
{
  double WeatheringRatePerYear, CurrentWeatherCoeff, mslope, wscale, CurrentPercentFine;

  /* Assign maximum weathering rate depending on rock type  */
  CurrentWeatherCoeff = type_of_rock[XRock[j]][YRock[j]] == 'f' ? FastWeatherCoeff : SlowWeatherCoeff;

  // TODO: This block will never be true
  if (CurrentWeatherCoeff != FastWeatherCoeff && CurrentWeatherCoeff != SlowWeatherCoeff)
  {
    CurrentWeatherCoeff = type_of_rock[XRockBehind[j]][YRockBehind[j]] == 'f' ? FastWeatherCoeff : SlowWeatherCoeff;
  }
  // TODO: also false logic
  if (CurrentWeatherCoeff != FastWeatherCoeff && CurrentWeatherCoeff != SlowWeatherCoeff)
  {
    printf("weatherrock (%d): rock & behind both have no rocktype\n", j);
  }

  // TODO: same as above
  /* Assign proportion of fines lost PWL */
  CurrentPercentFine = type_of_rock[XRock[j]][YRock[j]] == 'f' ? PercentFineFast : PercentFineSlow;
  if (CurrentPercentFine != PercentFineFast && CurrentPercentFine != PercentFineSlow)
  {
    CurrentPercentFine = type_of_rock[XRockBehind[j]][YRockBehind[j]] == 'f' ? PercentFineFast : PercentFineSlow;
  }

  /* Flip current percent fines to make calculations below easier PWL */
  CurrentPercentFine = 1 - CurrentPercentFine;

  /* Determine weathering: abrasion or no abrasion? */
  /* think in terms of vertical equivalence */
  if ((MinDistToBeach[j]) <= NO_WEATHERING)
  {
    if (ABRASION)
    {
      if (MinDistToBeach[j] <= W_CRIT)
      {
        mslope = (CurrentWeatherCoeff * (N - 1)) / W_CRIT; /* Slope of line between N*BareRock and BareRock PWL */
        WeatheringRatePerYear = CurrentWeatherCoeff + mslope * MinDistToBeach[j];
        // if (!debug12b)
        // {
        //   printf("\nabrasion weathering rate:%f\n",WeatheringRatePerYear);
        // }
      }
      else
      {
        /* Decay coefficient for the "cover" portion of abrasion curve PWL */
        wscale = (NO_WEATHERING - W_CRIT) / (log((0.01 / N)));
        WeatheringRatePerYear = kAngleFactor * CurrentWeatherCoeff * (exp((MinDistToBeach[j] - W_CRIT) / wscale));
      }
    }
    else
    {
      /* exponential weathering(in meters/yr) based on sed cover */
      WeatheringRatePerYear = CurrentWeatherCoeff * (exp(-MinDistToBeach[j]));
      // printf("weathering rate: %f", WeatheringRatePerYear);
    }
  }
  else
  {
    WeatheringRatePerYear = 0;
  }

  /* weathering per time step. units: cell units, or a percent weathered (which correlates with erosion) */
  AmountWeathered[j] = ((WeatheringRatePerYear * time_step) / 365) / CELL_WIDTH;

  if (debug12)
  {
    printf("for %d (%d,%d) rock %c, angfactor %f, weathercoeff %f, mindist %f, rate %f m/yr, frac weathered %f\n",
        j, XRock[j], YRock[j], type_of_rock[XRock[j]][YRock[j]], kAngleFactor,
        CurrentWeatherCoeff, MinDistToBeach[j], WeatheringRatePerYear, AmountWeathered[j]);
  }

  /* PWL, 2/16/12: added sea cliffs, which puts in an extra amount of sediment
   * to PercentFullSand. Use variable CliffHeight (top) to adjust. */
  /* Also, need to make sure model can deal with situations when there is no
   * rock behind current rock cell (i.e., back side of an island)   */
  if (PercentFullRock[XRock[j]][YRock[j]] >= AmountWeathered[j])
  {
    if (debug12)
    {
      printf("%d (%d,%d) weathers (1)\n", j, XRock[j], YRock[j]);
    }
    if (debugtopo)
    {
      printf("\n Cliff height (1): %f, j: %d\n", topography[XRock[j]][YRock[j]], j);
    }

    PercentFullRock[XRock[j]][YRock[j]] -= AmountWeathered[j];
    PercentFullSand[XRock[j]][YRock[j]] += AmountWeathered[j] + (AmountWeathered[j] * (topography[XRock[j]][YRock[j]] / DEPTH_SHOREFACE) * CurrentPercentFine);

    if (PercentFullRock[XRock[j]][YRock[j]] < 1.0)
    {
      AllRock[XRock[j]][YRock[j]] = FALSE;
    }
  }

  else if ((PercentFullRock[XRock[j]][YRock[j]] + PercentFullRock[XRockBehind[j]][YRockBehind[j]]) >= AmountWeathered[j])
  {
    if (debug12)
    {
      printf("%d (%d,%d) weathers (2)\n", j, XRock[j], YRock[j]);
    }
    if (debugtopo)
    {
      printf("\n Cliff height (2): %f\n", topography[XRock[j]][YRock[j]]);
      printf("\n Cliff height behind (2): %f\n", topography[XRockBehind[j]][YRockBehind[j]]);
    }

    PercentFullRock[XRockBehind[j]][YRockBehind[j]] -= (AmountWeathered[j] - PercentFullRock[XRock[j]][YRock[j]]);

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
    AllBeach[XRockBehind[j]][YRockBehind[j]] = FALSE;
  } 
  /* current + behind don't have enough rock to weather, so weather all there is in both cells. */
  else
  {
    if (debug12)
    {
      printf("%d (%d,%d) weathers (3)\n", j, XRock[j], YRock[j]);
    }
    if (debugtopo)
    {
      printf("\n Cliff height (3): %f\n", topography[XRock[j]][YRock[j]]);
      printf("\n Cliff height behind (3): %f\n", topography[XRockBehind[j]][YRockBehind[j]]);
    }

    if (YRock[j] > 1 && YRock[j] < 2 * Y_MAX - 2)
    {
      if (debug12)
      {
        printf("not enough rock to weather for %d (%d,%d). ", j, XRock[j], YRock[j]);
        printf("this is handle-able but prob. shouldn't happen often.\n");
      }
      // PauseRun(XRock[j],YRock[j],j);
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

/**
 * Moves along beach and tests to see if cells are in shadow
 * PARAMETERS: none
 * RETURN: none, changes TotalBeachCells, ShadowXMax, and InShadow[]
 */
void ShadowSweep(void)
{
  int i;

  /* Find maximum extent of beach to use as a limit for shadow searching */
  // ShadowXMax = XMaxBeach(ShadowXMax) + 3;
  // if(ShadowXMax > X_MAX)
  // {
  //   ShadowXMax = X_MAX; /* RCL */
  // }

  if (debug2)
  {
    printf("ShadowXMax: %d   XMaxBeach: %d \n", ShadowXMax, XMaxBeach(ShadowXMax));
  }

  /* Determine if beach cells are in shadow */
  for (i = 0; i <= TotalBeachCells; i++)
  {
    if (kSwanFlag)
    {
      InShadow[i] = FALSE; /* When using SWAN, turn off the shadows */
    }
    else
    {
      InShadow[i] = FindIfInShadow(X[i], Y[i], ShadowXMax);

      /* Code to see what is happening */
      if (InShadow[i])
      {
        if (debug2)
        {
          printf("Shadow at X: %3d  Y: %3d  i: %3d \n", X[i], Y[i], i);
        }
        else if (debug2)
        {
          printf("No Shadow X: %3d  Y: %3d  i: %3d \n", X[i], Y[i], i);
        }
      }
    }
  }
}
/**
 * Finds extens of beach in x direction
 * PARAMETERS: Max: starts searching at point 3 rows higher than max
 * RETURN: integer value equal to max exten of 'allbeach'
 */
int XMaxBeach(int Max)
{
  int xtest, ytest;
  xtest = Max + 3;
  ytest = 0;

  while (xtest > 0)
  {
    while (ytest < 2 * Y_MAX)
    {
      if (AllBeach[xtest][ytest])
      {
        return xtest;
      }
      ytest++;
    }
    ytest = 0;
    xtest--;
  }

  printf("***** Should've found X_MAX for shadow): %d, %d ***** \n", xtest, ytest);
  return X_MAX;
}

/**
 * Determines if a particular cell is in shadow
 * PARAMETERS: xin, yin: x and y coordinates of the cell
 * ShadMax: max extent of shadowed area
 * RETURN: returns whether in shadow is TRUE or FALSE
 */
int FindIfInShadow(int xin, int yin, int ShadMax)
{
  /* Initialize local variables */
  double xdistance, ydistance;

  int xcheck = xin;
  int ycheck = yin;
  int iteration = 1;

  /*  Find a cell along the projection line moving against wave direction */
  while ((xcheck < ShadMax) && (ycheck > -1) && (ycheck < (2 * Y_MAX)))
  {
    xdistance = iteration * SHADOW_STEP_DISTANCE * cos(WaveAngle);
    xcheck = xin + (int)rint(xdistance);

    ydistance = -iteration * SHADOW_STEP_DISTANCE * sin(WaveAngle);
    ycheck = yin + (int)rint(ydistance);

    if (debug2a)
    {
      printf("\n Iteration: %d \n", iteration);
      printf("Xdist : %f  ", xdistance);
      printf("    Xcheck: %d \n", xcheck);
      printf("Ydist : %f  ", ydistance);
      printf("    Ycheck: %d \n", ycheck);
      printf("firsttest : %f  ", (xcheck + .2));
      printf("secondtest : %f  ",(xin + 0.3 + fabs((ycheck - yin) / tan(WaveAngle))));
      printf("term : %f  \n", fabs((ycheck - yin) / tan(WaveAngle)));
    }

    if (xcheck >= X_MAX - 1 || ycheck >= 2 * Y_MAX - 1)
    {
      return FALSE;
    }

    /* If AllBeach is along the way, and not next neighbor */
    /* Probably won't get to this one, though */
    if ((AllBeach[xcheck][ycheck]) &&
        (// (abs(ycheck-yin) !=1 && abs(xcheck-xin) !=1) &&
        ((xcheck + 1) > (xin + (PercentFullSand[xin][yin] + PercentFullRock[xin][yin]) +fabs ((ycheck - yin) / tan(WaveAngle))))) &&
        (ycheck > -1) && 
        (ycheck < 2 * Y_MAX - 1))
    {
      if (debug2)
      {
        printf("Allbeach used Xck: %2d  Yck: %2d   \n", xcheck, ycheck);
        printf("	      Xin: %2d  Yin: %2d   ", xin, yin);
        printf("Xdist : %d  ", abs(xcheck - xin));
        printf("Ydist : %d  ", abs(ycheck - yin));
        printf("firsttest : %f  ", (xcheck + .2));
        printf("secondtest : %f  ", (xin + 0.3 + fabs((ycheck - yin) / tan(WaveAngle))));
        printf("term : %f  \n\n", fabs((ycheck - yin) / tan(WaveAngle)));
      }
      if (debug2a)
      {
        printf("finishing on a full cell with Xcheck %d, Xin %d, PFS+R %f  \n",
               xcheck, xin, (PercentFullSand[xcheck][ycheck] + PercentFullSand[xcheck][ycheck]));
      }
      return TRUE;
    }
    /* Compare a partially full cell's x - distance to a line projected from the starting beach cell's x-distance */
    /* This assumes that beach projection is in x-direction (not too bad)   */
    /* RCL this could be contributing to the jaggedness...need to take into account
       the way it's facing or could just ignore partially full cells completely... */
    #ifdef FGFG
    else if (((PercentFullSand[xcheck][ycheck] + PercentFullRock[xcheck][ycheck]) > 0) &&
            ((xcheck + (PercentFullSand[xcheck][ycheck] + PercentFullRock[xcheck][ycheck])) >
              (xin + (PercentFullSand[xin][yin] + PercentFullRock[xin][yin]) + fabs((ycheck - yin) / tan(WaveAngle)))) &&
            (ycheck > -1) && (ycheck < 2 * Y_MAX - 1))
            // && !(xcheck == xin))
      {
      if (debug2a)
      {
        printf("finishing on a partially full with Xcheck %d, Xin %d, PFS+R %f\n",
            xcheck, xin, (PercentFullSand[xcheck][ycheck] + PercentFullSand[xcheck][ycheck]));
      }
      return TRUE;
    }
    #endif
    iteration++;
  }

  /*  No shadow was found */
  return FALSE;
}

/**
 * Determines beach angles for all beach cells from left to right
 * ShorelineAnlge will apply to current cell and right neighbor
 * PARAMETERS: none
 * RETURN: none, writes to ShorelineAngle[], UpWind[], SurroundingAngle[]
 */
void DetermineAngles(void)
{
  int i, j, k; /* Local loop variables */

  /* Compute ShorelineAngle[]  */
  /* not equal to TotalBeachCells because angle between cell and rt neighbor */

  for (i = 0; i < TotalBeachCells; i++)
  {    
    /* On bottom side of Spit (assuming we must be if going right)  */
    if (Y[i] > Y[i + 1])
    {
      ShorelineAngle[i] = PI - atan((X[i + 1] - (PercentFullSand[X[i + 1]][Y[i + 1]] + 
                    PercentFullRock[X[i + 1]][Y[i + 1]])) - (X[i] - (PercentFullSand[X[i]][Y[i]]) + PercentFullRock[X[i]][Y[i]]));

      if (ShorelineAngle[i] > PI)
      {
        ShorelineAngle[i] -= 2.0 * PI;
        if (debug3)
        {
          printf("Next under -180 \n");
        }
      }

      if (debug3)
      {
        printf("(4a) i = %d  X[i]: %d Y[i]: %d Percent %3f \n\tX[+]: %d Y[+]: %d Percent{+}%3f   Angle: %f  Deg Angle: %f\n",
            i, X[i], Y[i], (PercentFullSand[X[i]][Y[i]] + PercentFullRock[X[i]][Y[i]]),
            X[i + 1], Y[i + 1], (PercentFullSand[X[i + 1]][Y[i + 1]] + PercentFullRock[X[i + 1]][Y[i + 1]]),
            ShorelineAngle[i], ShorelineAngle[i] * 180 / PI);
      }
      // PauseRun(X[i],Y[i],i);
    }
    /* Assume if going right, on regular shore */
    else if ((Y[i] < Y[i + 1]))
    {
      ShorelineAngle[i] = atan((X[i + 1] + (PercentFullSand[X[i + 1]][Y[i + 1]] +
                     PercentFullRock[X[i + 1]][Y[i + 1]])) - (X[i] + (PercentFullSand[X[i]][Y[i]] + PercentFullRock[X[i]][Y[i]])));

      if (debug3)
      {
        printf("(1) i = %d  X[i]: %d Y[i]: %d Percent %3f \n\tX[+]: %dY[+]: %d Percent[+] %3f   Angle: %f  Deg Angle: %f\n",
            i, X[i], Y[i], (PercentFullSand[X[i]][Y[i]] + PercentFullRock[X[i]][Y[i]]),
            X[i + 1], Y[i + 1], (PercentFullSand[X[i + 1]][Y[i + 1]] + PercentFullRock[X[i + 1]][Y[i + 1]]),
            ShorelineAngle[i], ShorelineAngle[i] * 180 / PI);
      }
    }
    /*  Shore up and down, on right side */
    else if (Y[i] == Y[i + 1] && X[i] > X[i + 1])
    {
      ShorelineAngle[i] = -PI / 2.0 + atan((Y[i + 1] + (PercentFullSand[X[i + 1]][Y[i + 1]] +  PercentFullRock[X[i + 1]][Y[i + 1]])) -
                           (Y[i] + (PercentFullSand[X[i]][Y[i]] + PercentFullRock[X[i]][Y[i]])));

      if (debug3)
      {
        printf(
            "(2) i = %d  X[i]: %d Y[i]: %d Percent %3f \n           X[+]: %d Y[+]: %d Percent[+] %3f   Angle: %f  Deg Angle: %f \n",
            i, X[i], Y[i], (PercentFullSand[X[i]][Y[i]] + PercentFullRock[X[i]][Y[i]]),
            X[i + 1], Y[i + 1], (PercentFullSand[X[i + 1]][Y[i + 1]] + PercentFullRock[X[i + 1]][Y[i + 1]]),
            ShorelineAngle[i], ShorelineAngle[i] * 180 / PI);
      }
    }
    /* Shore up and down, on left side */
    else if (Y[i] == Y[i + 1] && X[i] < X[i + 1])
    {
      ShorelineAngle[i] = PI / 2.0 + atan((Y[i + 1] + (PercentFullSand[X[i + 1]][Y[i + 1]] + PercentFullRock[X[i + 1]][Y[i + 1]])) -
                          (Y[i] + (PercentFullSand[X[i]][Y[i]] + PercentFullRock[X[i]][Y[i]])));

      if (debug3)
      {
        printf(
            "(3) i = %d  X[i]: %d Y[i]: %d Percent %3f \n           X[+]: %d Y[+]: %d Percent[+] %3f   Angle: %f  Deg Angle: %f \n",
            i, X[i], Y[i], (PercentFullSand[X[i]][Y[i]] + PercentFullRock[X[i]][Y[i]]),
            X[i + 1], Y[i + 1], (PercentFullSand[X[i + 1]][Y[i + 1]] + PercentFullRock[X[i + 1]][Y[i + 1]]),
            ShorelineAngle[i], ShorelineAngle[i] * 180 / PI);
      }
    }
    else
    {
      printf("Should've found ShorelineAngle): %d, %d \n", X[i], Y[i]);
      PauseRun(X[i], Y[i], i);
    }

    // printf("i = %d  X[i]: %d Y[i]: %d Percent %3f \n           X[+]: %d Y[+]: %d Percent %3f   Angle: %f  Deg Angle: %f \n",i, X[i], Y[i],
    // (PercentFullSand[X[i]][Y[i]] + PercentFullRock[X[i]][Y[i]]), X[i+1], Y[i+1], (PercentFullSand[X[i+1]][Y[i+1]] +
    // PercentFullRock[X[i+1]][Y[i+1]]), ShorelineAngle[i],ShorelineAngle[i]*180/PI);
  }

  /* compute SurroundingAngle array */
  /* 02/04 AA averaging doesn't work on bottom of spits */
  /* Use trick that x is less if on bottom of spit - angles must be different signs as well */
  for (k = 1; k < TotalBeachCells; k++)
  {
    if ((Y[k - 1] - Y[k + 1] == 2) && (copysign(ShorelineAngle[k - 1], ShorelineAngle[k]) != ShorelineAngle[k - 1]))
    {
      SurroundingAngle[k] = (ShorelineAngle[k - 1] + ShorelineAngle[k]) / 2 + PI;
      if (SurroundingAngle[k] > PI)
      {
        SurroundingAngle[k] -= 2.0 * PI;
      }
      if (debug4)
      {
        printf("Under: %d\n", k);
      }
    }
    else
    {
      SurroundingAngle[k] = (ShorelineAngle[k - 1] + ShorelineAngle[k]) / 2;
    }
  }

  /* Determine Upwind/downwind condition */
  /* Note - Surrounding angle is based upon left and right cell neighbors, and is centered on cell, not on right boundary */

  if (debug4)
  {
    printf("\nUp/Down   Wave Angle:%f\n", WaveAngle * RAD_TO_DEG);
  }

  for (j = 1; j < TotalBeachCells; j++)
  {
    if (debug4)
    {
      printf("i: %d  Shad: %c Ang[i]: %3.1f  Sur: %3.1f  Effect: %3f  ", j, InShadow[j], ShorelineAngle[j] * RAD_TO_DEG,
             SurroundingAngle[j] * RAD_TO_DEG, (WaveAngle - SurroundingAngle[j]) * RAD_TO_DEG);
    }

    if (fabs(WaveAngle - SurroundingAngle[j]) >= 42.0 / RAD_TO_DEG)
    {
      UpWind[j] = 'u';
      if (debug4)
      {
        printf("U(1)  ");
      }
    }
    else
    {
      UpWind[j] = 'd';
      if (debug4)
      {
        printf("D(1)  ");
      }
    }
    if (debug4)
    { 
      printf("\n");
    }
  }
}

/**
 * Loop function to determine which neigbor/situation to use for sediment transport calcs
 * Once situation is determined, will use function SedTrans to determine actual transport 
 * PARAMETERS: none
 * RETURN: none
 */
void DetermineSedTransport(void)
{
  int i;                  /* Loop variable LMV also sent to SedTrans to determine VolAcrossBorder */
  double ShoreAngleUsed;  /* Temporary holder for shoreline angle */
  int CalcCell;           /* Cell sediment coming from to go across boundary i */
  int Next, Last;         /* Indicators so test can go both left/right */
  int Correction;         /* Term needed for shoreline angle and i+1 case, angle stored at i */
  char UpWindLocal;       /* Local holder for upwind/downwind condition */
  int MaxTrans;           /* Do we need to compute using maximum transport? */
  int DoFlux;             /* Skip sed transport calcs (added 02/04 AA) */
  double DummyAngle = 1.;

  if (debug5) 
  {
    printf("\nSEDTRANS: %d  @  %f \n\n", current_time_step, WaveAngle * RAD_TO_DEG);
  }

  for (i = 1; i < TotalBeachCells - 1; i++)
  {
    if (debug5)
    {
      printf("\n  i: %d  ", i);
    }
    MaxTrans = FALSE;

    /*  Is littoral transport going left or right?  */
    /* #SWAN, 11/27/14: this isn't necessarily correct when SWAN is involved --
     *short-period waves can refract around and break going the opposite direction of their deep-water direction.  */
    /* Parse SWAN here first, and make this determination with 'breaking wave' angles */
    /* OY VAY */
    if (kSwanFlag)
    {
      /* Use dummy shore angle because we don't need it right now...prob a poor solution to the problem */
      ParseSWAN(i, DummyAngle); 
      /* Use SWAN nearshore angle instead */
      if ((Angle - ShorelineAngle[i]) > 0)
      {
        /*  Transport going right, center on cell to left side of borde */
        /*  Next cell in positive direction, no correction term needed */
        CalcCell = i;
        Next = 1;
        Last = -1;
        Correction = 0;
        DirectionAcrossBorder[i] = 'r'; /*LMV*/
        if (debug5)
        {
          printf("RT  ");
        }
      }
      else
      {
        /*  Transport going left, center on cell to right side of border */
        /*  Next cell in negative direction, correction term needed */
        CalcCell = i + 1;
        Next = -1;
        Last = 1;
        Correction = -1;
        DirectionAcrossBorder[i] = 'l'; /*LMV*/
        if (debug5)
        {
          printf("LT  ");
        }
      }
    }
    else
    {
      if ((WaveAngle - ShorelineAngle[i]) > 0)
      {
        /*  Transport going right, center on cell to left side of border */
        /*  Next cell in positive direction, no correction term needed */
        CalcCell = i;
        Next = 1;
        Last = -1;
        Correction = 0;
        DirectionAcrossBorder[i] = 'r'; /*LMV*/
        if (debug5)
        {
          printf("RT  ");
        }
      }
      else
      {
        /*  Transport going left, center on cell to right side of border */
        /*  Next cell in negative direction, correction term needed */
        CalcCell = i + 1;
        Next = -1;
        Last = 1;
        Correction = -1;
        DirectionAcrossBorder[i] = 'l'; /*LMV*/
        if (debug5)
        {
          printf("LT  ");
        }
      }
    }

    /*  Adjustment for maximum transport when passing through 45 degrees */
    /*  This adjustment is only made for moving from downwind to upwind conditions  */
    /*  purposefully done before shadow adjustment, only use maxtran when transition from dw to up not because of shadow */
    /*  keeping transition from uw to dw - does not seem to be big deal (04/02 AA) */
    if (!InShadow[CalcCell])
    {
      if ((UpWind[CalcCell] == 'd' && UpWind[CalcCell + Next] == 'u' && !InShadow[CalcCell + Next]) ||
          (UpWind[CalcCell + Last] == 'u' && UpWind[CalcCell] == 'd' && !InShadow[CalcCell + Last]))
          {
        MaxTrans = TRUE;
        if (debug5)
        {
          printf("MAXTRAN  ");
        }
      }

      /*  Upwind/Downwind adjustment Make sure sediment is put into shadows */
      /*  If Next cell is in shadow, use UpWind condition */
      DoFlux = TRUE;
      UpWindLocal = UpWind[CalcCell];

      if (InShadow[CalcCell + Next])
      {
        UpWindLocal = 'u';
        if (debug5)
        {
          printf("U(2)  ");
        }
      }

      /*  If coming out of shadow, downwind should be used */
      /*  HOWEVER- 02/04 AA - if high angle, will result in same flux in/out problem */
      /*  solution  - no flux for high angle waves */
      if ((InShadow[CalcCell + Last]) && (UpWindLocal == 'u'))
      {
        DoFlux = FALSE;
        if (debug5)
        {
          printf("U(X) NOFLUX \n");
        }
      }

      /*  Use upwind or downwind shoreline angle for calcs */
      if (UpWindLocal == 'u')
      {
        ShoreAngleUsed = ShorelineAngle[CalcCell + Last + Correction];
        if (debug5)
        {
          printf("UP  ShoreAngle: %3.1f  ", ShoreAngleUsed * RAD_TO_DEG);
        }
      }
      else if (UpWindLocal == 'd')
      {
        ShoreAngleUsed = ShorelineAngle[CalcCell + Correction];
        if (debug5)
        {
          printf("DN  ShoreAngle: %3.1f  ", ShoreAngleUsed * RAD_TO_DEG);
        }
      }

      /* !!! Do not do transport on unerneath c'cause it gets all messed up */
      if (fabs(ShoreAngleUsed) > PI / 2.0)
      {
        DoFlux = FALSE;
      }

      /* Send to SedTrans to calculate VolumeIn and VolumeOut */
      // printf("i = %d  Cell: %d NextCell: %d Angle: %f Trans Angle: %f\n",
      // i, CalcCell, CalcCell+Next, ShoreAngleUsed*180/PI, (WaveAngle - ShoreAngleUsed)*180/PI);

      if (debug5)
      {
        printf("From: %d  To: %d  TransAngle %3.1f", CalcCell, CalcCell + Next, (WaveAngle - ShoreAngleUsed) * RAD_TO_DEG);
      }

      if (DoFlux)
      {
        /*LMV*/
        // SedTrans (i, ShoreAngleUsed, MaxTrans);
        SedTrans(i, CalcCell, ShoreAngleUsed, MaxTrans, Last);
      }
    }
  }
}

/**
 * This central function will calcualte the sediment transported from the cell i-1 to
 *  the cell at i, using the input ShoreAngle (LMV)
 * PARAMETERS: i: current cell index, From: ???, ShoreAngle: angle between i and i -1, 
 * MaxT: max amount of transport, Last: ??? 
 * RETURN: none, writes to VolumeAcrossBorder[]
 */
void SedTrans(int i, int From, double ShoreAngle, int MaxT, int Last)
{
  /* Coefficients - some of these are important */
  double StartDepth = 3 * OffShoreWvHt; /* m, depth to begin refraction calcs (needs to be beyond breakers) */
  double RefractStep = .2;              /* m, step size to iterate depth for refraction calcs */
  double KBreak = 0.5;                  /* coefficient for wave breaking threshold */
  double rho = 1020;                    /* kg/m3 - density of water and dissolved matter */

  /* Variables */
  double AngleDeep;          /* rad, Angle of waves to shore at inner shelf  */
  double Depth = StartDepth; /* m, water depth for current iteration         */
  // double Angle;              /* rad, calculation angle                       */
  double CDeep;              /* m/s, phase velocity in deep water            */
  double LDeep;              /* m, offhsore wavelength                       */
  double C;                  /* m/s, current step phase velocity             */
  double kh;                 /* wavenumber times depth                       */
  double n;                  /* n                                            */
  double WaveLength;         /* m, current wavelength                        */

  /* Primary assumption is that waves refract over shore-parallel contours */
  /* New algorithm 6/02 iteratively takes wave onshore until they break, then computes Qs */
  /* See notes 06/05/02 */
  if (debug6)
  {
    printf("Wave Angle %2.2f Shore Angle  %2.2f    ", WaveAngle * RAD_TO_DEG, ShoreAngle * RAD_TO_DEG);
  }

  AngleDeep = WaveAngle - ShoreAngle;
  if (MaxT)
  {
    AngleDeep = PI / 4.0;
  }
  if (debug6)
  {
    printf("Deep Tranport Angle %2.2f \n\n", AngleDeep * RAD_TO_DEG);
  }

  /*  Don't do calculations if over 90 degrees, should be in shadow  */
  if (AngleDeep > 0.995 * PI / 2.0 || AngleDeep < -0.995 * PI / 2.0)
  {
    return;
  }
  else if (kSwanFlag) /* Do SWAN */
  {
    /* Instead of shoaling above, use ParseSWAN function to find the breaking wave characteristics. */
    ParseSWAN(From, ShoreAngle);
    /* Do Qs slightly different for SWAN input */
    VolumeAcrossBorder[i] = fabs(0.2 * rho * Raise(GRAVITY, 3.0 / 2.0) * Raise(WvHeight, 2.5) * cos(Angle - ShoreAngle) * sin(Angle - ShoreAngle) * time_step);
    /* Use wave characteristics updrift */
    if (UpWind[From] == 'u')
    {
      VolumeAcrossBorder[i] = fabs(0.2 * rho * Raise(GRAVITY, 3.0 / 2.0) * Raise(Hsigdebug[From + Last], 2.5) *
          cos(Dirdebug[From + Last] - ShoreAngle) * sin(Dirdebug[From + Last] - ShoreAngle) * time_step);
    }
  }
  else
  {
    /* Calculate Deep Water Celerity & Length, Komar 5.11 c = gT / PI, L = CT */
    CDeep = GRAVITY * Period / (2.0 * PI);
    LDeep = CDeep * Period;
    if (debug6) {
      printf("CDeep = %2.2f LDeep = %2.2f \n", CDeep, LDeep);
    }

    while (TRUE)
    {
      /* non-iterative eqn for L, from Fenton & McKee */
      WaveLength = LDeep * Raise(tanh(Raise(Raise(2.0 * PI / Period, 2) * Depth / GRAVITY, .75)), 2.0 / 3.0);
      C = WaveLength / Period;
      if (debug6)
      {
        printf("DEPTH: %2.2f Wavelength = %2.2f C = %2.2f ", Depth, WaveLength, C);
      }

      /* Determine n = 1/2(1+2kh/sinh(kh)) Komar 5.21 */
      /* First Calculate kh = 2 pi Depth/L  from k = 2 pi/L */
      kh = 2 * PI * Depth / WaveLength;
      n = 0.5 * (1 + 2.0 * kh / sinh(2.0 * kh));
      if (debug6)
      {
        printf("kh: %2.3f  n: %2.3f ", kh, n);
      }

      /* Calculate angle, assuming shore parallel contours and no conv/div of rays */
      /* from Komar 5.47 */
      Angle = asin(C / CDeep * sin(AngleDeep));
      if (debug6)
      {
        printf("Angle: %2.2f", Angle * RAD_TO_DEG);
      }

      /* Determine Wave height from refract calcs - Komar 5.49 */
      WvHeight = OffShoreWvHt * Raise(CDeep * cos(AngleDeep) / (C * 2.0 * n * cos(Angle)), .5);
      if (debug6)
      {
        printf(" WvHeight : %2.3f\n", WvHeight);
      }
      // TODO: verify this logic for breaking the loop
      if (WvHeight > Depth * KBreak || Depth <= RefractStep)
      {
        break;
      }
      Depth -= RefractStep;
    }

    /* Now Determine Transport */
    /* eq. 9.6b (10.8) Komar, including assumption of sed density = 2650 kg/m3 */
    /* additional accuracy here will not improve an already suspect eqn for sed transport */
    /* (especially with poorly constrained coefficients), */
    /* so no attempt made to make this a more perfect imperfection */
    VolumeAcrossBorder[i] = fabs(1.1 * rho * Raise(GRAVITY, 3.0 / 2.0) * Raise(WvHeight, 2.5) * cos(Angle) * sin(Angle) * time_step); /*LMV - now global array*/
  }

  /*LMV VolumeIn/Out is now calculated below in AdjustShore */
  if (debug6a)
  {
    printf("\ni: %d \n", i);
  }
  if (debug6a)
  {
    printf("VolumeAcrossBorder: %f  ", VolumeAcrossBorder[i]);
  }
}

/**
 * Empties sand from any cell that has been designated a sink
 * PARAMETERS: none
 * RETURN: none
 */
void DoSink(void)
{
  if (!HAVE_SINKS || RandZeroToOne() > SINKINESS)
  {
    return;
  }

  int s, i, foundOne = FALSE;

  /* treat sink as any beach cell at a particular y-value */
  if (COLUMN_SINKS)
  {
    for (s = 0; s < NUM_SINKS; s++) 
    {
      for (i = 0; i < MAX_BEACH_LENGTH; i++)
      {
        if (Y[i] == kSinkY[s])
        {
          PercentFullSand[X[i]][Y[i]] = 0.0;
          foundOne = TRUE;
        }
      }
      if (!foundOne)
      {
        printf("couldn't find sink %d at y == %d!\n", s, kSinkY[s]);
      }
    }
  }
  /* treat sink as a specific cell; empty it if it is on the beach */
  else
  {
    for (s = 0; s < NUM_SINKS; s++)
    {
      for (i = 0; i < MAX_BEACH_LENGTH; i++)
      {
        if (Y[i] == kSinkY[s] && X[i] == kSinkX[s])
        {
          PercentFullSand[X[i]][Y[i]] = 0.0;
          foundOne = TRUE;
        }
      }
      // if(!foundOne)
      // {
      //   printf("couldn't find sink %d at (%d,%d)!\n",s,kSinkX[s],kSinkY[s]);
      // }
    }
  }
  if (!foundOne)
  {
    printf("didn't find any sinks @ time %d\n", current_time_step);
  }
}

/**
 * Determines is sediment transport into a cell is from right, left, converging into the cell
 * or diverging out. Purpose is to address problems with insufficent sand, fix losing sand problem (LMV)
 * PARAMETERS: none
 * RETURN: none, writes to FlowThroughCell[]
 */
void FlowInCell(void)
{
  int i;
  // double AverageFull;
  double TotalPercentInBeaches = 0.;

  for (i = 1; i <= TotalBeachCells; i++)
  {
    if ((DirectionAcrossBorder[i - 1] == 'r') && (DirectionAcrossBorder[i] == 'r'))
    {
      FlowThroughCell[i] = 'R'; /*Right */
    }
    if ((DirectionAcrossBorder[i - 1] == 'r') && (DirectionAcrossBorder[i] == 'l'))
    {
      FlowThroughCell[i] = 'C'; /*Convergent */
    }
    if ((DirectionAcrossBorder[i - 1] == 'l') && (DirectionAcrossBorder[i] == 'r'))
    {
      FlowThroughCell[i] = 'D'; /*Divergent */
    }
    if ((DirectionAcrossBorder[i - 1] == 'l') && (DirectionAcrossBorder[i] == 'l'))
    {
      FlowThroughCell[i] = 'L'; /*Left */
    }

    TotalPercentInBeaches += (PercentFullSand[X[i]][Y[i]] + PercentFullRock[X[i]][Y[i]]);
  }

    // AverageFull = TotalPercentInBeaches/TotalBeachCells;
    // for (i=1; i <= TotalBeachCells; i++)
    // {
    //   if ((PercentFullSand[X[i]][Y[i]] + PercentFullRock[X[i]][Y[i]]) > AverageFull)
    //   {
    //     PutPixel(Y[i]*CellPixelSize , X[i]*CellPixelSize, 50, 200, 0);
    //   }
    //   else
    //   {
    //     PutPixel(Y[i]*CellPixelSize, X[i]*CellPixelSize, 200, 0, 00);
    //   }
    // }
}

/**
 * Checks and distributes sediment in adjacent cells based on FlowThroughCell[] (LMV)
 * Sweeps left to right to check cells where DirectionAcrossBorder='r' (D to R, D to C, R to R, R to C)
 * Sweeps right to left to check cells where DirectionAcrossBorder='l' (D to L, D to C, L to L, L to C)
 * PARAMETERS: none
 * RETURN: none
 */
void FixFlow(void)
{
  int i;
  double AmountSand;    /* amount of sand availible for transport (includes amount in cell behind) LMV */

  double Depth;         /* Depth of current cell */
  // double DeltaArea;    /* Holds change in area for cell (m^2) */
  double Distance;      /* distance from shore to intercept of equilib. profile and overall slope (m) */
  double Xintercept;    /* X position of intercept of equilib. profile and overall slope */
  double DepthEffective;/* depth at x intercept, Brad: where slope becomes < equilib slope, so sand stays piled against shoreface */

  /* get the actual volumes across for left and right borders of D cells */
  for (i = 1; i <= TotalBeachCells - 1; i++)
  {
    ActualVolumeAcross[i] = VolumeAcrossBorder[i];
    /* calculate effective depth */
    Depth = INITIAL_DEPTH + ((X[i] - INIT_BEACH) * CELL_WIDTH * SHELF_SLOPE);
    Distance = Depth / (SHOREFACE_SLOPE - SHELF_SLOPE * cos(ShorelineAngle[i]));
    Xintercept = X[i] + Distance * cos(ShorelineAngle[i]) / CELL_WIDTH;
    DepthEffective = INITIAL_DEPTH + ((Xintercept - INIT_BEACH) * CELL_WIDTH * SHELF_SLOPE);
    if (DepthEffective < DEPTH_SHOREFACE)
    {
      DepthEffective = DEPTH_SHOREFACE;
    }

    /*set the D's first - left and right border set */
    /*Flow from a Divergent cell (to a Convergent cell or a Right cell) */
    if (FlowThroughCell[i] == 'D')
    {
      if (AllBeach[XBehind[i]][YBehind[i]])
      {
        if (PercentFullRock[X[i]][Y[i]] == 0.0)
        {
          AmountSand = (PercentFullSand[X[i]][Y[i]] + PercentFullSand[XBehind[i]][YBehind[i]]) * CELL_WIDTH * CELL_WIDTH * DepthEffective;
          if (debug13)
          {
            printf("(D)X: %d, Y: %d\n", X[i], Y[i]);
          }
        }
        else
        {
          AmountSand = (PercentFullSand[X[i]][Y[i]]) * CELL_WIDTH * CELL_WIDTH * DepthEffective;
          if (debug13)
          {
            printf("(RockD)X: %d, Y: %d\n", X[i], Y[i]);
          }
        }
      }
      /* I don't know how much I have, just figure it out and give me some sand */
      else{
        // printf("Check Amount sand X: %d, Y: %d \n", X[i], Y[i]);
        AmountSand = (PercentFullSand[X[i]][Y[i]]) * CELL_WIDTH * CELL_WIDTH * DepthEffective;
        if (debug13)
        {
          printf("(Dd)X: %d, Y: %d\n", X[i], Y[i]);
        }
      }

      if (debug13)
      {
        printf("AmountSand[%d]: %f DepthEffective: %f\n", i, AmountSand, DepthEffective);
      }

      if ((VolumeAcrossBorder[i - 1] + VolumeAcrossBorder[i]) <= AmountSand)
      {
        ActualVolumeAcross[i] = VolumeAcrossBorder[i];
        ActualVolumeAcross[i - 1] = VolumeAcrossBorder[i - 1];

        if (debug13)
        {
          printf("D VolumeWantL[%d] %f, D VolumeGetL[%d] %f\n", i - 1, VolumeAcrossBorder[i - 1], i - 1, ActualVolumeAcross[i - 1]);
          printf("D VolumeWantR[%d] %f, D VolumeGetR[%d] %f\n", i, VolumeAcrossBorder[i], i, ActualVolumeAcross[i]);
          printf("VolumeOutL[%d] %f, VolumeOutR[%d] %f\n\n", i - 1, ActualVolumeAcross[i - 1], i, ActualVolumeAcross[i]);
        }
      }

      /* divide up proportionally */
      else
      {
        ActualVolumeAcross[i] = ((VolumeAcrossBorder[i] / (VolumeAcrossBorder[i - 1] + VolumeAcrossBorder[i])) * AmountSand);
        ActualVolumeAcross[i - 1] = ((VolumeAcrossBorder[i - 1] / (VolumeAcrossBorder[i - 1] + VolumeAcrossBorder[i])) * AmountSand);

        if (debug13)
        {
          printf("   D VolumeWantL[%d] %f, D VolumeGetL[%d] %f\n", i - 1, VolumeAcrossBorder[i - 1], i - 1, ActualVolumeAcross[i - 1]);
          printf("   D VolumeWantR[%d] %f, D VolumeGetR[%d] %f\n", i, VolumeAcrossBorder[i], i, ActualVolumeAcross[i]);
          printf("   VolumeOutL[%d] %f, VolumeOutR[%d] %f\n\n", i - 1, ActualVolumeAcross[i - 1], i, ActualVolumeAcross[i]);
        }
      }
    }
    /* Flow from a Right cell (to a Right cell or a Convergent cell) */
    else if (FlowThroughCell[i] == 'R')
    {
      if (AllBeach[XBehind[i]][YBehind[i]])
      {
        if (PercentFullRock[X[i]][Y[i]] == 0.0)
        {
          AmountSand = (PercentFullSand[X[i]][Y[i]] + PercentFullSand[XBehind[i]][YBehind[i]]) * CELL_WIDTH * CELL_WIDTH * DepthEffective;
          if (debug13)
          {
            printf("(R)X: %d, Y: %d\n", X[i], Y[i]);
          }
        }
        else
        {
          AmountSand = (PercentFullSand[X[i]][Y[i]]) * CELL_WIDTH * CELL_WIDTH * DepthEffective;
          if (debug13)
          {
            printf("(RockR)X: %d, Y: %d\n", X[i], Y[i]);
          }
        }
      }
      /* I don't know how much I have, just figure it out and give me some sand */
      else
      {
        // printf("Check 2 Amount sand X: %d, Y: %d \n", X[i], Y[i]);
        AmountSand = (PercentFullSand[X[i]][Y[i]]) * CELL_WIDTH * CELL_WIDTH * DepthEffective;
        if (debug13)
        {
          printf("(Rd)X: %d, Y: %d\n", X[i], Y[i]);
        }
      }

      if (debug13)
      {
        printf("AmountSand[%d]: %f DepthEffective: %f\n", i, AmountSand, DepthEffective);
      }

      if ((ActualVolumeAcross[i - 1] + AmountSand) >= VolumeAcrossBorder[i])
      {
        ActualVolumeAcross[i] = VolumeAcrossBorder[i];

        if (debug13)
        {
          printf("R VolumeWant[%d] %f, R VolumeGet[%d] %f\n", i, VolumeAcrossBorder[i], i, ActualVolumeAcross[i]);
          printf("VolumeIn[%d] %f, VolumeOut[%d] %f\n\n", i, ActualVolumeAcross[i - 1], i, ActualVolumeAcross[i]);
        }
      }
      else
      {
        ActualVolumeAcross[i] = ActualVolumeAcross[i - 1] + AmountSand;
        if (debug13)
        {
          printf("   R VolumeWant[%d] %f, R VolumeGet[%d] %f\n", i, VolumeAcrossBorder[i], i, ActualVolumeAcross[i]);
          printf("   VolumeIn[%d] %f, VolumeOut[%d] %f\n\n", i, ActualVolumeAcross[i - 1], i, ActualVolumeAcross[i]);
        }
      }
    }
    /* Flow from a Left cell (to a Convergent cell or a Left cell) */
    else if (FlowThroughCell[i] == 'L')
    {
      if (AllBeach[XBehind[i]][YBehind[i]])
      {
        if (PercentFullRock[X[i]][Y[i]] == 0.0)
        {
          AmountSand = (PercentFullSand[X[i]][Y[i]] + PercentFullSand[XBehind[i]][YBehind[i]]) * CELL_WIDTH * CELL_WIDTH * DepthEffective;
          if (debug13)
          {
            printf("(L)X: %d, Y: %d\n", X[i], Y[i]);
          }
        }
        else
        {
          AmountSand = (PercentFullSand[X[i]][Y[i]]) * CELL_WIDTH * CELL_WIDTH * DepthEffective;
          if (debug13)
          {
            printf("(RockL)X: %d, Y: %d\n", X[i], Y[i]);
          }
        }
      }
      /*I don't know how much I have, just figure it out and give me some sand */
      else
      {
        // printf("Check 3 Amount sand X: %d, Y: %d \n", X[i], Y[i]);
        AmountSand = (PercentFullSand[X[i]][Y[i]]) * CELL_WIDTH * CELL_WIDTH * DepthEffective;
        if (debug13)
        {
          printf("(Ld)X: %d, Y: %d\n", X[i], Y[i]);
        }
      }

      if (debug13)
      {
        printf("AmountSand[%d]: %f DepthEffective: %f\n", i, AmountSand, DepthEffective);
      }

      if ((ActualVolumeAcross[i] + AmountSand) >= VolumeAcrossBorder[i - 1])
      {
        ActualVolumeAcross[i - 1] = VolumeAcrossBorder[i - 1];
        if (debug13)
        {
          printf("L VolumeWant[%d] %f, L VolumeGet[%d] %f\n", i, VolumeAcrossBorder[i - 1], i, ActualVolumeAcross[i - 1]);
          printf("VolumeIn[%d] %f, VolumeOut[%d] %f\n\n", i, ActualVolumeAcross[i], i, ActualVolumeAcross[i - 1]);
        }
      }
      else
      {
        ActualVolumeAcross[i - 1] = ActualVolumeAcross[i] + AmountSand;
        if (debug13)
        {
          printf("L VolumeWant[%d] %f, L VolumeGet[%d] %f\n", i, VolumeAcrossBorder[i - 1], i, ActualVolumeAcross[i - 1]);
          printf("VolumeIn[%d] %f, VolumeOut[%d] %f\n\n", i, ActualVolumeAcross[i], i, ActualVolumeAcross[i - 1]);
        }
      }
    }
    if(debug13)
    {
      printf(
          "VolumeAcrossLeft[%d]: %f, VolumeAcrossRight[%d]: %f, PFS[X%d][Y%d]: "
          "%f, PFR[X%d][Y%d]: %f,\n",
          i, ActualVolumeAcross[i - 1], i, ActualVolumeAcross[i], X[i], Y[i],
          PercentFullSand[X[i]][Y[i]], X[i], Y[i], PercentFullRock[X[i]][Y[i]]);
    }
  }
}

/**
 * Sweep through cells to place transported sediment
 * PARAMETERS: none
 * RETURN: none
 */
void TransportSedimentSweep(void)
{
  int i, ii, sweepsign;

  sweepsign = RandZeroToOne() * 2 > 1 ? 1 : 0;

  if (debug7)
  {
    printf("\n\n TransSedSweep  Ang %f  %d\n", WaveAngle * RAD_TO_DEG, current_time_step);
  }

  for (i = 0; i <= TotalBeachCells - 1; i++)
  {
    ii = sweepsign ? i : TotalBeachCells - 1 - i;

    if (debug7)
    {
      printf("i: %d  ss: %d  X: %d  Y: %d  In: %.1f  Out: %.1f Across: %.1f\n",
             ii, sweepsign, X[i], Y[i], VolumeIn[i], VolumeOut[i], VolumeAcrossBorder[i]);
    }
    if (debug0 && current_time_step == 4 && (ii == 485 || ii == 486))
    {
      printf("pausing in tss before adjust (PFS[58][269] is %f)\n", PercentFullSand[58][269]);
      PauseRun(-1, -1, ii);
    }
    AdjustShore(ii);
    if (debug0 && current_time_step == 4 && (ii == 485 || ii == 486))
    {
      printf("pausing in tss after adjust (PFS[58][269] is %f)\n", PercentFullSand[58][269]);
      PauseRun(-1, -1, ii);
    }
    /* RCL changes 0.0 to 1E-6 whenever testing for oopsimempty */
    if ((PercentFullSand[X[ii]][Y[ii]] + PercentFullRock[X[ii]] [Y[ii]]) < -0.000001)
    {
      if (debug0 && current_time_step == 4 && (ii == 485 || ii == 486))
      {
        printf("pausing in tss after adjust, before oops (time 4, ii 485)\n");
        PauseRun(-1, -1, ii);
      }
      OopsImEmpty(X[ii], Y[ii]);
      // printf("Empty called in TSS\n");
    }
    else if ((PercentFullSand[X[ii]][Y[ii]] + PercentFullRock[X[ii]][Y[ii]]) > 1.0)
    {
      OopsImFull(X[ii], Y[ii]);
      // printf("Full called in TSS\n");
    }
    if (debug0 && current_time_step == 4 && (ii == 485 || ii == 486))
    {
      printf("pausing in tss after adjust, oops, before erode, oops (PFS[58][269] is %f)\n", PercentFullSand[58][269]);
      PauseRun(-1, -1, ii);
    }
    /* Don't erode the beach if there are no rocks... */
    if (INITIAL_CONDITION_TYPE != 3)
    {
      ErodeTheBeach(ii);
    }
    /* RCL */
    if (PercentFullSand[X[ii]][Y[ii]] < -0.000001)
    {
      OopsImEmpty(X[ii], Y[ii]);
      printf("Empty called in TSS after eroding\n");
    }
    if (debug0 && current_time_step == 4 && (ii == 485 || ii == 486))
    {
      printf("pausing in tss after adjust, oops, before erode, oops (PFS[58][269] is %f)\n", PercentFullSand[58][269]);
      PauseRun(-1, -1, ii);
    }
  }
}


/**
 * Sweep through cells to place transported sediment
 * PARAMETERS: none
 * RETURN: none
 */
/*  Complete mass balance for incoming and ougoing sediment */
/*  This function will change the global data array PercentFullSand[][] */
/*  Uses but does not adjust arrays: */
/*              X[], Y[], ShorelineAngle[], ActualVolumeAcross  LMV */
/*  Uses global variables: SHELF_SLOPE, CELL_WIDTH, SHOREFACE_SLOPE, INITIAL_DEPTH
*/
void AdjustShore(int i)
{
  double Depth;      /* Depth of current cell */
  double DeltaArea;  /* Holds change in area for cell (m^2) */
  double Distance;   /* distance from shore to intercept of equilib. profile and overall slope (m) */
  double Xintercept; /* X position of intercept of equilib. profile and overall slope */
  double DepthEffective; /* depth at x intercept, Brad: where slope becomes < equilib slope, */
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
  DepthEffective = INITIAL_DEPTH + ((Xintercept - INIT_BEACH) * CELL_WIDTH * SHELF_SLOPE);

  if (debug7a) {
    printf("\n Depth %f  Distance %f XIntercept %f  DepthEffective %f DepthShoreFace %f",
        Depth, Distance, Xintercept, DepthEffective, DEPTH_SHOREFACE);
  }

  if (DepthEffective < DEPTH_SHOREFACE) {
    DepthEffective = DEPTH_SHOREFACE;
  }
  /*LMV*/ 
  switch (FlowThroughCell[i]) {
    case 'L':
      VolumeIn[i] = ActualVolumeAcross[i];
      VolumeOut[i] = ActualVolumeAcross[i - 1];
      break;
    case 'C':
      VolumeIn[i] = ActualVolumeAcross[i];
      VolumeOut[i] = -ActualVolumeAcross[i - 1];
      break;
    case 'D':
      VolumeIn[i] = -ActualVolumeAcross[i - 1];
      VolumeOut[i] = ActualVolumeAcross[i];
      break;
    case 'R':
      VolumeIn[i] = ActualVolumeAcross[i - 1];
      VolumeOut[i] = ActualVolumeAcross[i];
      break;
  }

  DeltaArea = (VolumeIn[i] - VolumeOut[i]) / DepthEffective;
  /* if(you want to do sinks this way && this cell is sink) DeltaArea = 0;
     RCL addded sinks but won't do them like this probably */
  PercentFullSand[X[i]][Y[i]] += DeltaArea / (CELL_WIDTH * CELL_WIDTH);

  if (debug14) {
    printf("%c  i: %d  In: %f  Out: %f\n", FlowThroughCell[i], i, VolumeIn[i], VolumeOut[i]);
    printf("DeltaArea[%d]: %f, PFS[%d]: %f, time_step: %d\n", i, DeltaArea, i, PercentFullSand[X[i]][Y[i]], current_time_step);
    printf("	X[i]: %d, Y[i]: %d\n", X[i], Y[i]);
  }
}

/*  This function decreases the amount of sand in each cell of the
   TotalBeachCell array */
/*  by a uniform amount (ErosionRatePerYear) */
/*  Changes PercentFullSand[][] */
/*  LMV */
void ErodeTheBeach(int i)
{
  double PercentEroded;
  PercentEroded = InShadow[i] ? 0.0 : ((ErosionRatePerYear * time_step) / 365) / CELL_WIDTH;

  if (PercentEroded > PercentFullSand[X[i]][Y[i]]) {
    if (PercentFullRock[X[i]][Y[i]] == 0.0) {
      if (PercentEroded < (PercentFullSand[X[i]][Y[i]] +  PercentFullSand[XBehind[i]][YBehind[i]])) {
        if (debug15) {
          printf("	Bi: %d, PFSBefore: %f\n", i, PercentFullSand[X[i]][Y[i]]);
        }
        PercentFullSand[XBehind[i]][YBehind[i]] -= PercentEroded - PercentFullSand[X[i]][Y[i]];
        PercentFullSand[X[i]][Y[i]] = 0.0;
        if (debug15) {
          printf("	BPercentEroded: %f, PFSAfter: %f\n", PercentEroded, PercentFullSand[X[i]][Y[i]]);
        }
      }
      /* if PFR 0 && PercentEroded > PFS + PFSBehind... */
      else {
        /*if(!(PercentFullSand[X[i]][Y[i]] <= 0.000001
           && PercentFullRock[XBehind[i]][YBehind[i]] >= 0.999999)) {
           printf("Check Erosion, something strange (%d,%d)\n",X[i],Y[i]); */
        /*debug15 = 1;
           PrintLocalConds(X[i],Y[i]); */
        /*}*/

        if (debug15) {
          printf("\nErode it all!\n");
          printf("cell %d (%d,%d) PFSBefore: %2.20f with Behind cell (%d,%d) PFSBefore: %2.20f\n",
              i, X[i], Y[i], PercentFullSand[X[i]][Y[i]], XBehind[i], YBehind[i], PercentFullSand[XBehind[i]][YBehind[i]]);
        }
        PercentFullSand[XBehind[i]][YBehind[i]] = 0.0;
        PercentFullSand[X[i]][Y[i]] = 0.0;
        if (debug15) {
          printf("PercentEroded: %f (PFS for cell and behind cell set to 0)\n", PercentEroded);
        }
        /*if(debug15) PauseRun(-1,-1,i);
           debug15 = 0; */
      }
    }
    else {
      PercentFullSand[X[i]][Y[i]] = 0.0;
      if (debug15) {
        printf("BPercentEroded: %f, PFSAfter: %f\n", PercentEroded, PercentFullSand[X[i]][Y[i]]);
      }
    }
  } 
  /*normal case, enough sand */
  else {
    if (debug15) {
      printf("i: %d, PFSBefore: %f\n", i, PercentFullSand[X[i]][Y[i]]);
    }
    PercentFullSand[X[i]][Y[i]] -= PercentEroded;
    if (debug15) {
      printf("PercentEroded: %f, PFSAfter: %f\n", PercentEroded, PercentFullSand[X[i]][Y[i]]);
    }
  }
}

/**
 * If a cell is under-full, this will find source for desparity and move brach in
 * LMV - Needs to be adjusted -- can't take sand from rock cells,
 * PARAMETERS: x, y: coordinates of the empty cell
 * RETURN: none, write to AllBeach[][] and PercentFullSand[][]
 */
void OopsImEmpty(int x, int y)
{
  int emptycells = 0;
  int emptycells2 = 0;
  
  if (debug8)
  {
    printf("\n	OOPS I'm EMPTY!  X: %d  Y: %d PFS: %f PFR: %f", x, y, PercentFullSand[x][y], PercentFullRock[x][y]);
  }
  /* find out how many AllBeaches to take from */ /* LMV */
  if ((AllBeach[x - 1][y]) && (PercentFullRock[x - 1][y] == 0.0))
  {
     emptycells += 1;
  }
  if ((AllBeach[x + 1][y]) && (PercentFullRock[x + 1][y] == 0.0))
  {
    emptycells += 1;
  }
  if ((AllBeach[x][y - 1]) && (PercentFullRock[x][y - 1] == 0.0))
  {
    emptycells += 1;
  }
  if ((AllBeach[x][y + 1]) && (PercentFullRock[x][y + 1] == 0.0))
  {
    emptycells += 1;
  }
  if (debug8)
  {
    printf("\n(counted %d AllBeaches with no rock to take from)\n", emptycells);
  }

  // TODO: save pointers to empty cells instead of rechecking each one
  if (emptycells > 0) {
    /* Now Move Sediment */ /* LMV */
   if ((AllBeach[x - 1][y]) && (PercentFullRock[x - 1][y] == 0.0))
   {
      PercentFullSand[x - 1][y] += PercentFullSand[x][y] / emptycells;
      AllBeach[x - 1][y] = FALSE;
      if (debug8)
      {
        printf("  E MOVEDBACK\n");
        printf("  PFR behind: %f\n", PercentFullRock[x - 1][y]);
      }
      if (debug8a && PercentFullSand[x - 1][y] < 0.0)
      {
        printf("PFS went Negative in OopsI'mEmpty x-1: %d y:%d %f \n", x - 1, y, PercentFullSand[x - 1][y]);
      }
      else if (debug8a && PercentFullSand[x - 1][y] > 1.0)
      {
        printf("PFS went >1 in OopsI'mEmpty x-1: %d y:%d %f \n", x - 1, y, PercentFullSand[x - 1][y]);
      }
    }

    if ((AllBeach[x + 1][y]) && (PercentFullRock[x + 1][y] == 0.0))
    {
      PercentFullSand[x + 1][y] += PercentFullSand[x][y] / emptycells;
      AllBeach[x + 1][y] = FALSE;
      if (debug8)
      {
        printf("  E MOVEDUP\n");
        printf("  PFR up: %f\n", PercentFullRock[x + 1][y]);
      }
      if (debug8a && PercentFullSand[x + 1][y] < 0.0)
      {
        printf("PFS went Negative in OopsI'mEmpty x+1: %d y:%d %f \n", x + 1, y, PercentFullSand[x + 1][y]);
      }
      else if (debug8a && PercentFullSand[x + 1][y] > 1.0)
      {
        printf("PFS went >1 in OopsI'mEmpty x+1: %d y:%d %f \n", x + 1, y, PercentFullSand[x + 1][y]);
      }
    }

    if ((AllBeach[x][y - 1]) && (PercentFullRock[x][y - 1] == 0.0))
    {
      PercentFullSand[x][y - 1] += PercentFullSand[x][y] / emptycells;
      AllBeach[x][y - 1] = FALSE;
      if (debug8)
      {
        printf("  E MOVEDLEFT\n");
        printf("  PFR left: %f\n", PercentFullRock[x][y - 1]);
        // PauseRun(x, y, -1);
      }
      if (debug8a && PercentFullSand[x][y - 1] < 0.0)
      {
        printf("PFS went Negative in OopsI'mEmpty x: %d y-1:%d %f \n", x, y - 1, PercentFullSand[x][y - 1]);
      }
      else if (debug8a && PercentFullSand[x][y - 1] > 1.0)
      {
        printf("PFS went >1 in OopsI'mEmpty x: %d y-1:%d %f \n", x, y - 1,  PercentFullSand[x][y - 1]);
      }
    }

    if ((AllBeach[x][y + 1]) && (PercentFullRock[x][y + 1] == 0.0)) {
      PercentFullSand[x][y + 1] += PercentFullSand[x][y] / emptycells;
      AllBeach[x][y + 1] = FALSE;
      if (debug8) {
        printf("  E MOVEDRIGHT\n");
        printf("  PFR right: %f\n", PercentFullRock[x][y + 1]);
        // PauseRun(x, y, -1);
      }
      if (debug8a && PercentFullSand[x][y + 1] < 0.0)
      {
        printf("PFS went Negative in OopsI'mEmpty x: %d y+1:%d %f \n", x, y + 1, PercentFullSand[x][y + 1]);
      }
      else if (debug8a && PercentFullSand[x][y - 1] > 1.0)
      {
        printf("PFS went >1 in OopsI'mEmpty x: %d y-1:%d %f \n", x, y - 1, PercentFullSand[x][y - 1]);
      }
    }
  }
  /* No full neighbors, so take away from partially full neighbors */
  else {
    if (PercentFullSand[x - 1][y] > 0.0)
    {
      emptycells2 += 1;
    }
    if (PercentFullSand[x + 1][y] > 0.0)
    {
      emptycells2 += 1;
    }
    if (PercentFullSand[x][y - 1] > 0.0)
    {
      emptycells2 += 1;
    }
    if (PercentFullSand[x][y + 1] > 0.0)
    {
      emptycells2 += 1;
    }
    if (debug8)
    {
      printf("no full neighbors, counted %d partials\n", emptycells2);
    }

    if (emptycells2 > 0)
    {
      if (PercentFullSand[x - 1][y] > 0.0)
      {
        PercentFullSand[x - 1][y] += PercentFullSand[x][y] / emptycells2;
        if (debug8)
        {
          printf("  E NOTFULL MOVEDBACK\n");}
          printf("  PFR behind: %f\n", PercentFullRock[x - 1][y]);
        }
        if (debug8a && PercentFullSand[x - 1][y] < 0.0)
        {
          printf("PFS went Negative in OopsI'mEmpty x-1: %d y:%d %f \n", x - 1, y, PercentFullSand[x - 1][y]);
        }
        else if (debug8a && PercentFullSand[x - 1][y] > 1.0)
        {
          printf("PFS went >1 in OopsI'mEmpty x-1: %d y:%d %f \n", x - 1, y, PercentFullSand[x - 1][y]);
        }
      }

      if (PercentFullSand[x + 1][y] > 0.0)
      {
        PercentFullSand[x + 1][y] += PercentFullSand[x][y] / emptycells2;
        if (debug8)
        {
          printf("  E NOTFULL MOVEDUP\n");
          printf("  PFR up: %f\n", PercentFullRock[x + 1][y]);
        }
        if (debug8a && PercentFullSand[x + 1][y] < 0.0)
        {
          printf("PFS went Negative in OopsI'mEmpty x+1: %d y:%d %f \n", x + 1, y, PercentFullSand[x + 1][y]);
        }
        else if (debug8a && PercentFullSand[x + 1][y] > 1.0)
        {
          printf("PFS went >1 in OopsI'mEmpty x+1: %d y:%d %f \n", x + 1, y, PercentFullSand[x + 1][y]);
        }
      }

      if (PercentFullSand[x][y - 1] > 0.0)
      {
        PercentFullSand[x][y - 1] += PercentFullSand[x][y] / emptycells2;
        if (debug8)
        {
          printf("  E NOTFULL MOVEDLEFT\n");
          printf("  PFR left: %f\n", PercentFullRock[x][y - 1]);
          // PauseRun(x, y, -1);
        }

        if (debug8a && PercentFullSand[x][y - 1] < 0.0)
        {
          printf("PFS went Negative in OopsI'mEmpty...Bummer x: %d y-1:%d %f \n", x, y - 1, PercentFullSand[x][y - 1]);
        }
        else if (debug8a && PercentFullSand[x][y - 1] > 1.0)
        {
          printf("PFS went >1 in OopsI'mEmpty...Bummer x: %d y-1:%d %f \n", x, y - 1, PercentFullSand[x][y - 1]);
        }
      }

      if (PercentFullSand[x][y + 1] > 0.0)
      {
        PercentFullSand[x][y + 1] += PercentFullSand[x][y] / emptycells2;
        if (debug8)
        {
          printf(" E NOTFULL MOVEDRIGHT\n");
          printf("  PFR right: %f\n", PercentFullRock[x][y + 1]);
          // PauseRun(x, y, -1);
        }
        if (debug8a && PercentFullSand[x][y + 1] < 0.0)
        {
          printf("PFS went Negative in OopsI'mEmpty...Bummer x: %d y+1:%d %f \n", x, y + 1, PercentFullSand[x][y + 1]);
        }
        else if (debug8a && PercentFullSand[x][y - 1] > 1.0)
        {
          printf("PFS went >1 in OopsI'mEmpty...Bummer x: %d y-1:%d %f \n", x, y - 1, PercentFullSand[x][y - 1]);
        }
      }
      else
      {
        printf("@@@ Didn't find anywhere to steal sand from!! x: %d  y: %d\n", x, y);
        // PauseRun(x, y, -1);
      }
  }

  if (PercentFullSand[x][y] + PercentFullRock[x][y] < 1.0)
  {
    AllBeach[x][y] = FALSE;
  }
  /*LMV*/
  PercentFullSand[x][y] = 0.0;

  if (debug8)
  {
    printf("\n");
  }
  if (debug0 && PercentFullSand[58][269] < 0.0)
  {
    printf("^^^^^^(58,259) made underfull (%f) after calling oops(%d,%d)^^^^\n", PercentFullSand[58][269], x, y);
  }
}

/**
 * If a cell is overfull, push beach out in new direction
 * PARAMETERS: x, y: coordinates of the empty cell
 * RETURN: none, write to AllBeach[][] and PercentFullSand[][]
 */
void OopsImFull(int x, int y)
{
  int fillcells = 0;
  int fillcells2 = 0;

  if (debug8)
  {
    printf("\n	OOPS I'M FULL: X: %d  Y: %d PFS: %f PFR: %f", x, y, PercentFullSand[x][y], PercentFullRock[x][y]);
  }
  // if (debug8)
  // {
  //   PrintLocalConds(x,y,-1);
  // }
  /* find out how many cells will be filled up    */
  if ((PercentFullSand[x - 1][y] + PercentFullRock[x - 1][y]) == 0.0)
  {
    fillcells += 1;
  }
  if ((PercentFullSand[x + 1][y] + PercentFullRock[x + 1][y]) == 0.0)
  {
    fillcells += 1;
  }
  if ((PercentFullSand[x][y - 1] + PercentFullRock[x][y - 1]) == 0.0)
  {
    fillcells += 1;
  }
  if ((PercentFullSand[x][y + 1] + PercentFullRock[x][y + 1]) == 0.0)
  {
    fillcells += 1;
  }

  if (fillcells != 0)
  {
    /* Now Move Sediment */ /* LMV */
    if ((PercentFullSand[x - 1][y] + PercentFullRock[x - 1][y]) == 0.0)
    {
      PercentFullSand[x - 1][y] += ((PercentFullSand[x][y] + PercentFullRock[x][y]) - 1.0) / fillcells;
      if (debug80)
      {
        printf("  F MOVEDBACK\n");
      }
      if (debug8)
      {
        printf("  PFR behind: %f\n", PercentFullRock[x - 1][y]);
      }

      if (debug80 && PercentFullSand[x - 1][y] < 0.0)
      {
        printf("1 PFS went Negative in OopsI'mFull...Bummer x-1: %d y:%d %f \n", x - 1, y, PercentFullSand[x - 1][y]);
      }
      else if (debug80 && PercentFullSand[x - 1][y] > 1.0)
      {
        printf("1 PFS went >1 in OopsI'mFull...Bummer x-1: %d y:%d %f \n", x - 1, y, PercentFullSand[x - 1][y]);
      }
    }

    if ((PercentFullSand[x + 1][y] + PercentFullRock[x + 1][y]) == 0.0)
    {
      PercentFullSand[x + 1][y] += ((PercentFullSand[x][y] + PercentFullRock[x][y]) - 1.0) / fillcells;
      if (debug80)
      {
        printf("  F MOVEDUP\n");
      }
      if (debug8)
      {
        printf("  PFR up: %f\n", PercentFullRock[x + 1][y]);
      }

      if (debug80 && PercentFullSand[x + 1][y] < 0.0)
      {
        printf("1 PFS went Negative in OopsI'mFull...Bummer x+1: %d y:%d %f \n", x + 1, y, PercentFullSand[x + 1][y]);
      }
      else if (debug80 && PercentFullSand[x + 1][y] > 1.0)
      {
        printf("1 PFS went >1 in OopsI'mFull...Bummer x+1: %d y:%d %f \n", x + 1, y, PercentFullSand[x + 1][y]);
      }
    }

    if ((PercentFullSand[x][y - 1] + PercentFullRock[x][y - 1]) == 0.0)
    {
      PercentFullSand[x][y - 1] += ((PercentFullSand[x][y] + PercentFullRock[x][y]) - 1.0) / fillcells;
      if (debug80)
      {
        printf("  F MOVEDLEFT\n");
      }
      if (debug8)
      {
        printf("  PFR left: %f\n", PercentFullRock[x][y - 1]);
        // PauseRun(x,y,-1)
      }

      if (debug80 && PercentFullSand[x][y - 1] < 0.0)
      {
        printf("1 PFS went Negative in OopsI'mFull...Bummer x: %d y-1:%d %f \n", x, y - 1, PercentFullSand[x][y - 1]);
      }
      else if (debug80 && PercentFullSand[x][y - 1] > 1.0)
      {
        printf("1 PFS went >1 in OopsI'mFull...Bummer x: %d y-1:%d %f \n", x, y - 1, PercentFullSand[x][y - 1]);
      }
    }

    if ((PercentFullSand[x][y + 1] + PercentFullRock[x][y + 1]) == 0.0)
    {
      PercentFullSand[x][y + 1] += ((PercentFullSand[x][y] + PercentFullRock[x][y]) - 1.0) / fillcells;
      if (debug80)
      {
        printf("  F MOVEDRIGHT\n");
      }
      if (debug8)
      {
        printf("  PFR right: %f\n", PercentFullRock[x][y + 1]);
        // PauseRun(x,y,-1)
      }

      if (debug80 && PercentFullSand[x][y + 1] < 0.0)
      {
        printf("1 PFS went Negative in OopsI'mFull...Bummer x: %d y+1:%d %f \n", x, y + 1, PercentFullSand[x][y + 1]);
      }
      else if (debug80 && PercentFullSand[x][y + 1] > 1.0)
      {
        printf("1 PFS went >1 in OopsI'mFull...Bummer x: %d y+1:%d %f \n", x, y + 1, PercentFullSand[x][y + 1]);
      }
    }
  }
  else
  {
    /* No fully empty neighbors, so distribute to partially full neighbors */
    if ((PercentFullSand[x - 1][y] + PercentFullRock[x - 1][y]) < 1.0)
    {
      fillcells2 += 1;
    }
    if ((PercentFullSand[x + 1][y] + PercentFullRock[x + 1][y]) < 1.0)
    {
      fillcells2 += 1;
    }
    if ((PercentFullSand[x][y - 1] + PercentFullRock[x][y - 1]) < 1.0)
    {
      fillcells2 += 1;
    }
    if ((PercentFullSand[x][y + 1] + PercentFullRock[x][y + 1]) < 1.0)
    {
      fillcells2 += 1;
    }

    if (fillcells2 > 0)
    {
      if ((PercentFullSand[x - 1][y] + PercentFullRock[x - 1][y]) < 1.0)
      {
        PercentFullSand[x - 1][y] += ((PercentFullSand[x][y] + PercentFullRock[x][y]) - 1.0) / fillcells2;
        if (debug80)
        {
          printf("  Fa MOVEDBACK\n");
        }
        if (debug8)
        {
          printf("  PFR behind: %f\n", PercentFullRock[x - 1][y]);
        }

        if (debug80 && PercentFullSand[x - 1][y] < 0.0)
        {
          printf("PFS went Negative in OopsI'mFull...Bummer x-1: %d y:%d %f \n", x - 1, y, PercentFullSand[x - 1][y]);
        }
        else if (debug80 && PercentFullSand[x - 1][y] > 1.0)
        {
          printf("PFS went >1 in OopsI'mFull...Bummer x-1: %d y:%d %f \n", x - 1, y, PercentFullSand[x - 1][y]);
        }
      }

      if ((PercentFullSand[x + 1][y] + PercentFullRock[x + 1][y]) < 1.0)
      {
        PercentFullSand[x + 1][y] += ((PercentFullSand[x][y] + PercentFullRock[x][y]) - 1.0) / fillcells2;
        if (debug80)
        {
          printf("  Fa MOVEDUP\n");
        }
        if (debug8)
        {
          printf("  PFR up: %f\n", PercentFullRock[x + 1][y]);
        }

        if (debug80 && PercentFullSand[x + 1][y] < 0.0)
        {
          printf("PFS went Negative in OopsI'mFull...Bummer x+1: %d y:%d %f \n", x + 1, y, PercentFullSand[x + 1][y]);
        }
        else if (debug80 && PercentFullSand[x + 1][y] > 1.0)
        {
          printf("PFS went >1 in OopsI'mFull...Bummer x+1: %d y:%d %f \n", x + 1, y, PercentFullSand[x + 1][y]);
        }
      }

      if ((PercentFullSand[x][y - 1] + PercentFullRock[x][y - 1]) < 1.0)
      {
        PercentFullSand[x][y - 1] += ((PercentFullSand[x][y] + PercentFullRock[x][y]) - 1.0) / fillcells2;
        if (debug80)
        {
          printf("  Fa MOVEDLEFT\n");
        }
        if (debug8)
        {
          printf("  PFR left: %f\n", PercentFullRock[x][y - 1]);
        }

        if (debug80 && PercentFullSand[x][y - 1] < 0.0) 
        {
          printf("PFS went Negative in OopsI'mFull...Bummer x: %d y-1:%d %f \n", x, y - 1, PercentFullSand[x][y - 1]);
        }
        else if (debug80 && PercentFullSand[x][y - 1] > 1.0)
        {
          printf("PFS went >1 in OopsI'mFull...Bummer x: %d y-1:%d %f \n", x, y - 1, PercentFullSand[x][y - 1]);
        }
      }

      if ((PercentFullSand[x][y + 1] + PercentFullRock[x][y + 1]) < 1.0)
      {
        PercentFullSand[x][y + 1] += ((PercentFullSand[x][y] + PercentFullRock[x][y]) - 1.0) / fillcells2;
        if (debug80)
        {
          printf("  Fa MOVEDRIGHT\n");
        }
        if (debug8)
        {
          printf("  PFR right: %f\n", PercentFullRock[x][y + 1]);
        }

        if (debug80 && PercentFullSand[x][y + 1] < 0.0)
        {
          printf("PFS went Negative in OopsI'mFull...Bummer x: %d y+1:%d %f \n", x, y + 1, PercentFullSand[x][y + 1]);
        }
        else if (debug80 && PercentFullSand[x][y + 1] > 1.0)
        {
          printf("PFS went >1 in OopsI'mFull...Bummer x: %d y+1:%d %f \n", x, y + 1, PercentFullSand[x][y + 1]);
        }
      }
    }
    else
    {
      if (debug8)
      {
        printf("Nobody wants our sand!!! x: %d  y: %d  PFS: %f  PFR: %f\n", x, y, PercentFullSand[x][y], PercentFullRock[x][y]);
        printf("x+1: %f, y+1: %f, x-1: %f, y-1: %f\n",
               (PercentFullSand[x + 1][y] + PercentFullRock[x + 1][y]),
               (PercentFullSand[x][y + 1] + PercentFullRock[x][y + 1]),
               (PercentFullSand[x - 1][y] + PercentFullRock[x][y + 1]),
               (PercentFullSand[x][y - 1] + PercentFullRock[x][y - 1]));
        // PauseRun(x,y,-1);
        printf("(excess sand will be removed; sand = 1 - rock)\n");
      }
    }
  }

  if (PercentFullSand[x][y] + PercentFullRock[x][y] >= 1.0)
  {
    AllBeach[x][y] = TRUE;
  }
  /*LMV*/ 
  PercentFullSand[x][y] = 1.0 - PercentFullRock[x][y];

  if (debug80 && PercentFullSand[x][y] < 0.0)
  {
    printf("	PFS went Negative in OopsI'mFull...Bummer x: %d y:%d %f \n", x, y, PercentFullSand[x][y]);
  }

  if (debug80 && PercentFullSand[x][y] > 1.0)
  {
    printf("	PFS went >1 in OopsI'mFull...Bummer x: %d y:%d %f \n", x, y, PercentFullSand[x][y]);
  }
}

/**
 * Hopefully addresses strange problems caused by filling/emptying of cells
 * Looks at entire data set
 * Find unattached pieces of sand and moves them back to the shore
 * Takes care of 'doubleing bits' of sand
 * Also takes care of over/under filled beach pieces
 * Revised 5/21/02 to move sand to all adjacent neighbors sandrevt.c
 * Changes global variable PercentFullSand[][]
 * Uses and can change AllBeach[][]
 * sandrevx.c - added sweepsign to reduce chances of asymmetrical artifacts
 * LMV rock fix - do not fix bare bedrock cells
 * PARAMETERS: x, y: coordinates of the empty cell
 * RETURN: none, write to AllBeach[][] and PercentFullSand[][]
 */
void FixBeach(void)
{
  int i, x, y, sweepsign, done, counter;
  int fillcells3 = 0;
  int corner = FALSE;
   int xstart; // LMV

  // if (debug9)
  // {
  //   printf("\n\nFIXBEACH      %d     %f\n", current_time_step, WaveAngle*RAD_TO_DEG);
  // }
  sweepsign = RandZeroToOne() * 2 > 1 ? 1 : 0;

  for (x = 1; x < ShadowXMax - 1; x++)
  {
    /* RCL: changing loops to ignore border */
    for (i = 1; i < 2 * Y_MAX - 1; i++)
    {
      y = sweepsign ? i : 2 * Y_MAX - i;

      /* Take care of corner problem?
         if  (((AllBeach[x][y] == 'n') && (PercentFullRock[x][y] > 0.0)) &&
         (((AllBeach[x][y-1] == 'n') && (AllBeach[x-1][y] == 'n') &&
         (AllRock[x-1][y-1] == 'y'))
         ||((AllBeach[x-1][y] == 'n') && (AllBeach[x][y+1] == 'n') &&
         (AllRock[x-1][y+1] == 'y'))))

         corner = 1;
         else */
          //corner = 0;

      /* Take care of situations that shouldn't exist */
      if ((PercentFullSand[x][y] + PercentFullRock[x][y]) < -0.000001)
      {
        /* RCL changed  0.0 to 10E-6--doubleing-pt equality is messy */
        // printf("Too empty\n");
        if (debug9 && y != 0)
        {
          printf("\nUnder 0 Percent X: %d  Y: %d PFS: %f PFR: %f\n", x, y, PercentFullSand[x][y], PercentFullRock[x][y]);
        }
        if (debug0)
        {
          printf("pausing in fixbeach before before it calls oops\n");
          PauseRun(x, y, -1);
        }
        AllBeach[x][y] = FALSE;
        OopsImEmpty(x, y);
      }
      else if ((PercentFullSand[x][y] + PercentFullRock[x][y]) > 1.0)
      {
        // printf("Too full\n");
        AllBeach[x][y] = TRUE;
        if (debug9 && y != 0)
        {
          printf("\nOver 100 Percent X: %d  Y: %d PFS: %f PFR: %f\n", x, y, PercentFullSand[x][y], PercentFullRock[x][y]);
        }
        OopsImFull(x, y);
      }

      if ((PercentFullSand[x][y] + PercentFullRock[x][y] >= 0.0) &&
         (PercentFullSand[x][y] + PercentFullRock[x][y] < 1.0) && AllBeach[x][y])
      {
        AllBeach[x][y] = FALSE;
        if (debug9 && y != 0)
        {
          printf("\nALLBeachProb X: %d  Y: %d\n", x, y);
          printf("PFS: %f, PFR: %f\n", PercentFullSand[x][y], PercentFullRock[x][y]);
        }
      }

      /* Take care of 'loose' bits of sand */
      if (!corner)
      {
        fillcells3 = 0;
        /* Beach in cell, but bottom, top, right, and left neighbors not all full */
        if (((PercentFullSand[x][y] + PercentFullRock[x][y]) != 0.0) &&
            ((PercentFullSand[x - 1][y] + PercentFullRock[x - 1][y]) < 1.0) &&
            ((PercentFullSand[x + 1][y] + PercentFullRock[x + 1][y]) < 1.0) &&
            ((PercentFullSand[x][y + 1] + PercentFullRock[x][y + 1]) < 1.0) &&
            ((PercentFullSand[x][y - 1] + PercentFullRock[x][y - 1]) < 1.0) &&
            !(AllBeach[x][y]))
        {
          if (debug9 && y != 0)
          {
            printf("\nFB Moved loose bit of sand,  X: %d  Y: %d  PFS: %f PFR: %f ", x, y, PercentFullSand[x][y], PercentFullRock[x][y]);
          }

          /* distribute to partially full neighbors */
          if (((PercentFullSand[x - 1][y] + PercentFullRock[x - 1][y]) < 1.0) && ((PercentFullSand[x - 1][y] + PercentFullRock[x - 1][y]) > 0.0))
          {
            fillcells3 += 1;
          }
          if (((PercentFullSand[x + 1][y] + PercentFullRock[x + 1][y]) < 1.0) && ((PercentFullSand[x + 1][y] + PercentFullRock[x + 1][y]) > 0.0))
          {
            fillcells3 += 1;
          }
          if (((PercentFullSand[x][y + 1] + PercentFullRock[x][y + 1]) < 1.0) && ((PercentFullSand[x][y + 1] + PercentFullRock[x][y + 1]) > 0.0))
          {
            fillcells3 += 1;
          }
          if (((PercentFullSand[x][y - 1] + PercentFullRock[x][y - 1]) < 1.0) && ((PercentFullSand[x][y - 1] + PercentFullRock[x][y - 1]) > 0.0))
          {
            fillcells3 += 1;
          }

          if ((fillcells3 > 0) && (fillcells3 <= 4))
          {
            if (((PercentFullSand[x - 1][y] + PercentFullRock[x - 1][y]) < 1.0) && ((PercentFullSand[x - 1][y] + PercentFullRock[x - 1][y]) > 0.0))
            {
              // PercentFullSand[x-1][y] += (PercentFullSand[x][y]-PercentFullRock[x][y])/fillcells3;
              PercentFullSand[x - 1][y] += PercentFullSand[x][y] / fillcells3; /* RCL */
              if (debug9)
              {
                printf("  FB MOVEDBACK\n");
              }
              if (debug9a && PercentFullSand[x - 1][y] < 0.0)
              {
                printf("PFS went Negative in FixBeach...Bummer x-1: %d y:%d %f \n", x - 1, y, PercentFullSand[x - 1][y]);
              }
              else if (debug9a && PercentFullSand[x - 1][y] > 1.0)
              {
                printf("PFS went >1 in FixBeach...Bummer x-1: %d y:%d %f \n", x - 1, y, PercentFullSand[x - 1][y]);
              }
            }

            if (((PercentFullSand[x + 1][y] + PercentFullRock[x + 1][y]) < 1.0) && ((PercentFullSand[x + 1][y] + PercentFullRock[x + 1][y]) > 0.0))
            {
              // PercentFullSand[x+1][y] += (PercentFullSand[x][y]-PercentFullRock[x][y])/fillcells3;
              PercentFullSand[x + 1][y] += PercentFullSand[x][y] / fillcells3; /* RCL */
              if (debug9)
              {
                printf("  FB MOVEDUP\n");
              }
              if (debug9a && PercentFullSand[x + 1][y] < 0.0)
              {
                printf("PFS went Negative in FixBeach...Bummer x+1: %d y:%d %f \n", x + 1, y, PercentFullSand[x + 1][y]);
              }
              else if (debug9a && PercentFullSand[x + 1][y] > 1.0)
              {
                printf("!!!PFS went >1 in FixBeach...Bummer x+1: %d y:%d %f \n", x + 1, y, PercentFullSand[x + 1][y]);
              }
            }

            if (((PercentFullSand[x][y - 1] + PercentFullRock[x][y - 1]) < 1.0) && ((PercentFullSand[x][y - 1] + PercentFullRock[x][y - 1]) > 0.0))
            {
              // PercentFullSand[x][y-1] += (PercentFullSand[x][y]-PercentFullRock[x][y])/fillcells3;
              PercentFullSand[x][y - 1] += PercentFullSand[x][y] / fillcells3; /* RCL */
              if (debug9)
              {
                printf("  FB MOVEDLEFT\n");
                // PauseRun(x,y,-1);
              }
              if (debug9a && PercentFullSand[x][y - 1] < 0.0)
              {
                printf("PFS went Negative in FixBeach...Bummer x+1: %d y:%d %f \n", x + 1, y, PercentFullSand[x][y - 1]);
              }
              else if (debug9a && PercentFullSand[x - 1][y] > 1.0)
              {
                printf("PFS went >1 in FixBeach...Bummer x+1: %d y:%d %f \n", x + 1, y, PercentFullSand[x + 1][y]);
              }
            }

            if (((PercentFullSand[x][y + 1] + PercentFullRock[x][y + 1]) < 1.0) && ((PercentFullSand[x][y + 1] + PercentFullRock[x][y + 1]) > 0.0))
            {
              // PercentFullSand[x][y+1] += (PercentFullSand[x][y]-PercentFullRock[x][y])/fillcells3;
              PercentFullSand[x][y + 1] += PercentFullSand[x][y] / fillcells3; /* RCL */
              if (debug9)
              {
                printf("  FB MOVEDRIGHT\n");                
                // PauseRun(x,y,-1);
              }
              if (debug9a && PercentFullSand[x][y + 1] < 0.0)
              {
                printf("PFS went Negative in FixBeach...Bummer x: %d y+1:%d %f \n", x, y + 1, PercentFullSand[x][y + 1]);
              }
              else if (debug9a && PercentFullSand[x][y + 1] > 1.0)
              {
                printf("PFS went >1 in FixBeach...Bummer x: %d y+1:%d %f \n", x, y + 1, PercentFullSand[x][y + 1]);
              }
            }
          }
          else {
            // if (debug9)
            printf("Complete fixbeach breakdown x: %d  y: %d\n", x, y);
            // if (debug9) PauseRun(x,y,-1);

            xstart = X_MAX - 1;
            while (!AllBeach[xstart][y])
            {
              xstart -= 1;
            }
            xstart += 1;

            /* printf("Moving random sand from X: %d, Y: %d, to X: %d, Y: %d\n",
               x, y, xstart, y);
               printf("moving this much: %f\n", PercentFullSand[x][y]);
               printf("PFS Before move: %f\n", PercentFullSand[xstart][y]); */

            PercentFullSand[xstart][y] += PercentFullSand[x][y];
            PercentFullSand[x][y] = 0.0;
            AllBeach[x][y] = FALSE;

            printf("PFS After move: %f \n", PercentFullSand[xstart][y]);

            xstart = X_MAX - 1;
            while (!AllRock[xstart][y])
            {
              xstart -= 1;
            }
            xstart += 1;

            printf("Moving bit of rock from X: %d, Y: %d, to X: %d, Y: %d\n", x, y, xstart, y);
            printf("moving this much: %f\n", PercentFullRock[x][y]);
            printf("PFR of receiving cell before move: %f\n", PercentFullRock[xstart][y]);

            PercentFullRock[xstart][y] += PercentFullRock[x][y];
            PercentFullRock[x][y] = 0.0;
            AllRock[x][y] = FALSE;
            printf("PFR After move: %f \n", PercentFullRock[xstart][y]);
            if (PercentFullRock[xstart][y] >= 1.0)
            {
              AllRock[xstart][y] = TRUE;
              PercentFullRock[xstart + 1][y] += (PercentFullRock[xstart][y] - 1.0);
              PercentFullRock[xstart][y] = 1.0;
              if (PercentFullRock[xstart + 1][y] >= 1.0)
              {
                printf("rock spilled over\n");
              }
              printf("after second move PFR = %f ", PercentFullRock[xstart][y]);
              printf("and above cell PFR = %f \n", PercentFullRock[xstart + 1][y]);
            }
          }
          PercentFullSand[x][y] = 0.0;
          AllBeach[x][y] = FALSE;

          if (debug9)
          {
            printf("\n");
          }

          /* If we have overfilled any of the cells in this loop, need to OopsImFull() */

          if ((PercentFullSand[x - 1][y] + PercentFullRock[x - 1][y]) > 1.0)
          {
            /* printf("     Below Overfilled\n");
               printf("PF[x-1][y]: %f", (PercentFullSand[x-1][y] +
               PercentFullRock[x-1][y]));
               printf("x-1:%d, y:%d\n", (x-1), y); */
            OopsImFull(x - 1, y);
            if (debug9)
            {
              printf("	Below Overfilled\n");
            }
          }
          if ((PercentFullSand[x][y - 1] + PercentFullRock[x][y - 1]) > 1.0)
          {
            /* printf("     Left Side Overfilled\n");
               printf("PF[x][y-1]: %f", (PercentFullSand[x][y-1] +
               PercentFullRock[x][y-1]));
               printf("x:%d, y-1:%d\n", x, (y-1)); */
            OopsImFull(x, y - 1);
            if (debug9) {
              printf("	Left Side Overfilled\n");
            }
          }
          if ((PercentFullSand[x][y + 1] + PercentFullRock[x][y + 1]) > 1.0)
          {
            /* printf("     Right Side Overfilled\n");
               printf("PF[x][y+1]: %f", (PercentFullSand[x][y+1] +
               PercentFullRock[x][y+1]));
               printf("x:%d, y+1:%d\n", x, (y+1)); */
            OopsImFull(x, y + 1);
            if (debug9)
            {
              printf("	Right Side Overfilled\n");
            }
          }
          if ((PercentFullSand[x + 1][y] + PercentFullRock[x + 1][y]) > 1.0)
          {
            /* printf("     Top Side Overfilled\n");
               printf("PF[x+1][y]: %f", (PercentFullSand[x+1][y] +
               PercentFullRock[x+1][y]));
               printf("x+1:%d, y:%d\n", (x+1), y); */
            OopsImFull(x + 1, y);
            if (debug9)
            {
              printf("	Top Overfilled\n");
            }
          }
        }
      }
      else if (debug9)
      {
        printf("Corner cell, no fixbeach please PFR: %f x: %d, y: %d\n", PercentFullRock[x][y], x, y);
      }
    }
  }
  /* now make sure there aren't any overfull or underfull cells */
  if (debug0)
  {
    printf("pausing in fixbeach before final loop\n");
    PauseRun(58, 269, -1);
  }

  if (debug9)
  {
    printf("checking grid for overfull, underfull \n");
  }
  done = FALSE;
  for (counter = 0; counter < 10 && !done; counter++)
  {
    done = TRUE;
    for (x = 1; x < X_MAX - 1; x++)
    {
      for (y = 1; y < 2 * Y_MAX - 1; y++)
      {
        if ((PercentFullSand[x][y] + PercentFullRock[x][y] > 1.0) ||
            (PercentFullSand[x][y] + PercentFullRock[x][y] < 0.0) ||
            PercentFullSand[x][y] < 0.0)
        {
          done = FALSE;
        }
        if (PercentFullSand[x][y] + PercentFullRock[x][y] > 1.0)
        {
          OopsImFull(x, y);
          if (debug9)
          {
            printf("found an overfull in final fb loop (%d, %d) \n", x, y);
          }
        }
        else if (PercentFullSand[x][y] + PercentFullRock[x][y] < -0.000001)
        {
          OopsImEmpty(x, y);
          if (debug9)
          {
            printf("found an underfull in final fb loop (%d, %d) \n", x, y);
          }
        }
        /* RCL adds this, changes 0 to 1E-6 everywhere */
        if (PercentFullSand[x][y] < -0.000001)
        {
          OopsImEmpty(x, y);
          if (debug9)
          {
            printf("found a sand < 0 in final fb loop (%d, %d) \n", x, y);
          }
        }
      }
    }
  }
  if (debug9)
  {
    printf("done checking grid for overfull, underfull \n");
  }
}

/**
 * Counts the total colume occupied by beach cells
 * Uses same algorithm as AdjustShore TODO: pull out said algorithm
 * PARAMETERS: none
 * RETURN: double of the total sum
 */
double MassCount(void)
{
  int x, y;
  double Depth; /* Depth of current cell */
  double Mass = 0.0;
  // TODO: only multiply by depth once per x value
  for (x = 0; x < X_MAX; x++)
  {
    Depth = INITIAL_DEPTH + ((x - INIT_BEACH) * CELL_WIDTH * SHELF_SLOPE);
    for (y = 0; y < 2 * Y_MAX; y++)
    {
      Mass += PercentFullSand[x][y] * Depth;
    }
  }
  return Mass;
}

/**
 * Raises b to the e power. Note: Can't use this when b is negative
 * PARAMETERS: b: base number, e: exponent number
 * RETURN: double of the total sum
 */
double Raise(double b, double e)
{
  if (b > 0)
  {
    return pow(b, e);
  }
  else
  {
    printf("Raise: can't handle negative base \n");
    printf("Wave data: %f %lf %lf\n", WaveAngle, Period, OffShoreWvHt);
    PauseRun(-1, -1, -1);
    return 0;
  }
}

/**
 * Generates a random number equally distributed between zero and one
 * PARAMETERS: none
 * RETURN: double of the rnadom number
 */
double RandZeroToOne(void)
{
  /* return random()/(Raise(2,31)-1); */
  double AB_rand = rand() % 1000;
  AB_rand = (AB_rand + 1) / 1000;
  return ((double)(AB_rand));
}

/**
 * Creates initial beach conditions
 * Flat beach with zone of AllBeach = TRUE separated by AllBeach = FALSE
 * Block of rock parallel to beach begins at INIT_ROCK LMV
 * LMV Assume that AllBeach = TRUE for AllRock cells (and All Full cells)
 * Bounding layer set to random fraction of fullness
 * PARAMETERS: none
 * RETURN: none
 */
void InitNormal(void)
{
  int x, y;
  int initial_beach;
  int initial_rock;
  const double amp = 10.; /* Amplitude of cos curve */

  printf("Condition Initial \n");

  for (y = 0; y < 2 * Y_MAX; y++)
  {
    for (x = 0; x < X_MAX; x++)
    {
      initial_beach = INIT_BEACH;
      initial_rock = INIT_ROCK;
      /* shoreline is a cosine curve, diffusive LMV */
      if (DIFFUSIVE_HUMP)
      {
        initial_rock = (int)(INIT_ROCK + (-amp * cos((2 * PI * y) / Y_MAX)) - (INIT_BEACH - INIT_ROCK - 1));
        initial_beach = (int)(INIT_BEACH + (-amp * cos((2 * PI * y) / Y_MAX)));
      }
      if (x < initial_rock)      /* LMV */
      { 
        PercentFullRock[x][y] = 1.0;
        PercentFullSand[x][y] = 0.0;
        AllRock[x][y] = TRUE;
        AllBeach[x][y] = FALSE;
        topography[x][y] = kCliffHeightSlow;
      }
      else if (x == initial_rock) /* LMV */
      {
        if (INITIAL_SMOOTH_ROCK)
        {
          PercentFullRock[x][y] = 0.20;
          PercentFullSand[x][y] = 0.80;
        }
        else
        {
          PercentFullRock[x][y] = RandZeroToOne();
          PercentFullSand[x][y] = 1 - PercentFullRock[x][y];
        }
        AllRock[x][y] = FALSE;
        AllBeach[x][y] = TRUE;
        topography[x][y] = kCliffHeightSlow;
      }
      else if ((x > initial_rock) && (x < initial_beach)) /* LMV */
      {
        PercentFullRock[x][y] = 0.0;
        PercentFullSand[x][y] = 1.0;
        AllRock[x][y] = FALSE;
        AllBeach[x][y] = TRUE;
        topography[x][y] = 0;

      }
      else if (x == initial_beach)
      {
        if (INITIAL_SMOOTH)
        {
          PercentFullSand[x][y] = 0.5;
        }
        else
        {
          PercentFullSand[x][y] = RandZeroToOne();
        }
        AllRock[x][y] = FALSE;
        AllBeach[x][y] = TRUE;
        topography[x][y] = 0;

      }
      else if (x > initial_beach)
      {
        PercentFullSand[x][y] = 0;
        AllRock[x][y] = FALSE;
        AllBeach[x][y] = FALSE;
        topography[x][y] = 0; /* No cliffs in the ocean... PWL */
      }
      Age[x][y] = 0;
    }
  }

  for (x = 0; x < initial_rock - 3; x++)
  {
    for (y = 0; y < 2 * Y_MAX; y++)
    {
        type_of_rock[x][y] = 'f';
        topography[x][y] = kCliffHeightFast;
    }
  }
  set_rock_blocks(type_of_rock + initial_rock - 3, topography + initial_rock - 3, X_MAX - (initial_rock - 3), 2 * Y_MAX, NUMBER_CHUNK);
}

/**
 * Creates blocks of rock and weathering types
 * PARAMETERS: none
 * RETURN: success status
 */
// TODO: standardize void versus return 0
int initBlock(void)
{
  int x, y, n, blockHeight = 4, blockWidth = 60, NumBlocks = 0; /* makes NumBlocks evenly spaced blocks */
  printf("init block\n");

  for (x = 0; x < X_MAX; x++)
  {
    for (y = Y_MAX / 2; y < 3 * Y_MAX / 2; y++)
    {
      if (x < INIT_ROCK)
      {
        PercentFullRock[x][y] = 1.0;
        PercentFullSand[x][y] = 0.0;
        AllRock[x][y] = TRUE;
        type_of_rock[x][y] = 'f';
        AllBeach[x][y] = TRUE;
      }
      else if (x < INIT_BEACH)
      {
        PercentFullRock[x][y] = 0.0;
        PercentFullSand[x][y] = 1.0;
        AllRock[x][y] = FALSE;
        AllBeach[x][y] = FALSE;
      }
      else if (x == INIT_BEACH)
      {
        PercentFullRock[x][y] = 0.0;
        PercentFullSand[x][y] = RandZeroToOne();
        AllRock[x][y] = FALSE;
        AllBeach[x][y] = FALSE;
      }
      else
      {
        PercentFullRock[x][y] = 0.0;
        PercentFullSand[x][y] = 0.0;
        AllRock[x][y] = FALSE;
        AllBeach[x][y] = FALSE;
      }
      Age[x][y] = 0;
    }
  }
  /* the regularly-spaced blocks */
  for (n = 1; n <= NumBlocks; n++)
  {
    for (x = INIT_ROCK; x >= INIT_ROCK - blockHeight; x--)
    {
      /* block starts @ upper left corner */
      for (y = Y_MAX / 2 + n * Y_MAX / (NumBlocks + 1); y < Y_MAX / 2 + n * Y_MAX / (NumBlocks + 1) + blockWidth &&  y < 2 * Y_MAX;y++)
      {
        PercentFullRock[x][y] = 1.0;
        PercentFullSand[x][y] = 0.0;
        AllRock[x][y] = TRUE;
        type_of_rock[x][y] = 's';
        AllBeach[x][y] = TRUE;
      }
    }
  }
  /* now here's space for ad hoc blocks */
  /* ad hoc block 1 */
  for (x = INIT_ROCK; x >= INIT_ROCK - blockHeight; x--)
  {
    for (y = Y_MAX / 2; y < Y_MAX / 2 + 10; y++)
    {
      PercentFullRock[x][y] = 1.0;
      PercentFullSand[x][y] = 0.0;
      AllRock[x][y] = TRUE;
      type_of_rock[x][y] = 's';
      AllBeach[x][y] = TRUE;
    }
  }
  /* ad hoc block 2 */
  for (x = INIT_ROCK; x >= INIT_ROCK - blockHeight; x--)
  {
    for (y = 3 * Y_MAX / 2 - 90; y < 3 * Y_MAX / 2; y++)
    {
      PercentFullRock[x][y] = 1.0;
      PercentFullSand[x][y] = 0.0;
      AllRock[x][y] = TRUE;
      type_of_rock[x][y] = 's';
      AllBeach[x][y] = TRUE;
    }
  }

  printf("done block init\n");
  return 0;
}

/**
 * Creates a sinusoidal beach
 * PARAMETERS: none
 * RETURN: success status
 */
int InitWiggly(void) {
  int x, y, curve;
  int NumCurves = 4;
  int RockLine[2 * Y_MAX];
  int curveFactor =  3; /* decrease fundamental wavelength/increase frequency (by a factor) --make it 1 for wavelength = Y_MAX */
  double Amp, AmpMax = INIT_ROCK / 6.0; /* INIT_ROCK/4.0; maximum amplitude of the component sin waves */
  double Min = 0.10 * X_MAX, Max = 0.95 *X_MAX; /* determines acceptable max amplitude of resulting curve */
  printf("init wiggly\n");

  for (y = 0; y < 2 * Y_MAX; y++)
  {
    RockLine[y] = INIT_ROCK;
  }

  for (curve = 1; curve <= NumCurves; curve++)
  {
    Amp = RandZeroToOne() * AmpMax;
    printf("for curve %d, amplitude %f\n", curve, Amp);
    for (y = 0; y < 2 * Y_MAX; y++)
    {
      RockLine[y] =  (int)RockLine[y] + Amp * sin(y * (curve * curveFactor) * 2.0 * PI / Y_MAX);
      if (RockLine[y] < Min || RockLine[y] > Max)
      {
        return -1;
      }
    }
  }

  for (x = 0; x < X_MAX; x++)
  {
    for (y = 0; y < 2 * Y_MAX; y++)
    {
      if (x < RockLine[y])
      {
        PercentFullRock[x][y] = 1.0;
        PercentFullSand[x][y] = 0.0;
        AllRock[x][y] = TRUE;
        type_of_rock[x][y] = 's';
        AllBeach[x][y] = TRUE;
      }
      else if (x == RockLine[y])
      {
        PercentFullRock[x][y] = 0.0;
        PercentFullSand[x][y] = 0.5;
        AllRock[x][y] = FALSE;
        AllBeach[x][y] = FALSE;
      }
      else
      {
        PercentFullRock[x][y] = 0.0;
        PercentFullSand[x][y] = 0.0;
        AllRock[x][y] = FALSE;
        AllBeach[x][y] = FALSE;
      }
      Age[x][y] = 0;
    }
  }

  printf("done wiggly init\n");
  return 0;
}

/**
 * Initializes beach conditions
 * PARAMETERS: none
 * RETURN: none
 */
void InitConds(void)
{
  current_time_step = 0;
  current_time = 0.;

  if (INITIAL_CONDITION_TYPE == 0)
  {
    InitNormal();
  }
  else if (INITIAL_CONDITION_TYPE == 1)
  {
    while (InitWiggly() == -1)
    {
      printf("init: amplitude too great, trying again...\n");
    }
  }
  else if (INITIAL_CONDITION_TYPE == 2)
  {
    initBlock();
  }
  /* All sand, all the time -- no rocks allowed */
  else if (INITIAL_CONDITION_TYPE == 3)
  {
    int x, y;
    printf("Condition Initial \n");

    for (y = 0; y < 2 * Y_MAX; y++)
    {
      for (x = 0; x < X_MAX; x++)
      {
        /* This is the only place where cell_depth is defined -- it needs to be updated through time
         for overwash functions (i think...), so perhaps there is something missing in this model version?? */
        cell_depth[x][y] = INITIAL_DEPTH + ((x - INIT_BEACH) * CELL_WIDTH * SHELF_SLOPE);

        if (x < INIT_BEACH)
        {
          PercentFullSand[x][y] = 1;
          PercentFullRock[x][y] = 0;
          AllBeach[x][y] = TRUE;
          cell_depth[x][y] = -LAND_HEIGHT;
        }
        else if (x == INIT_BEACH)
        {
          if (INITIAL_SMOOTH) {
            PercentFullSand[x][y] = 0.5;
            PercentFullRock[x][y] = 0;

          }
          else
          {
            PercentFullSand[x][y] = RandZeroToOne();
            PercentFullRock[x][y] = 0;
            printf("x: %d  Y: %d  Per: %f\n", x, y, PercentFullSand[x][y]);
          }
          AllBeach[x][y] = FALSE;
          cell_depth[x][y] = -LAND_HEIGHT;
        }
        else if (x > INIT_BEACH)
        {
          PercentFullSand[x][y] = 0;
          PercentFullRock[x][y] = 0;
          AllBeach[x][y] = FALSE;
          if (cell_depth[x][y] < DEPTH_SHOREFACE)
          {
            cell_depth[x][y] = DEPTH_SHOREFACE;
          }
        }
        else
        {
          printf("ugh! x: %d  Y: %d  Per: %f\n", x, y, PercentFullSand[x][y]);
          PauseRun(x, y, -1);
        }

        Age[x][y] = 0;
      }
    }
  }
  else
  {
    printf("have to pick an INITIAL_CONDITION_TYPE\n");
  }
  return;
}

/**
 * Add a perturbation to an otherwise flat line beach
 * (Andrew's initial bump)
 * PARAMETERS: none
 * RETURN: none
 */
void InitPert(void)
{
  int x, y;
  int PWidth = 5;
  int PHeight = 3;
  int PYstart = 25;

  /* Square perturbation */
  if (INITIAL_PERT == 1)
  {
    /* Fill AllBeach areas */
    for (x = INIT_BEACH; x <= INIT_BEACH + PHeight; x++)
    {
      for (y = PYstart; y <= PYstart + PWidth; y++)
      {
        PercentFullSand[x][y] = 1.0;
        AllBeach[x][y] = TRUE;
      }
    }
    /* PercentFull Top */
    for (y = PYstart - 1; y <= PYstart + PWidth + 1; y++)
    {
      PercentFullSand[INIT_BEACH + PHeight + 1][y] = RandZeroToOne();
    }
    /* PercentFull Sides */
    for (x = INIT_BEACH; x <= INIT_BEACH + PHeight; x++)
    {
      PercentFullSand[x][PYstart - 1] = RandZeroToOne();
      PercentFullSand[x][PYstart + PWidth + 1] = RandZeroToOne();
    }
  }
  /* Another Perturbation  - steep point */
  else if (INITIAL_PERT == 2)
  {
    x = INIT_BEACH;

    PercentFullSand[x][17] = 0.8;
    PercentFullSand[x][18] = 1.0;
    AllBeach[x][18] = TRUE;
    PercentFullSand[x][19] = 0.8;

    x = INIT_BEACH + 1;

    PercentFullSand[x][17] = 0.6;
    PercentFullSand[x][18] = 1.0;
    AllBeach[x][18] = TRUE;
    PercentFullSand[x][19] = 0.6;

    x = INIT_BEACH + 2;

    PercentFullSand[x][17] = 0.2;
    PercentFullSand[x][18] = 1.0;
    AllBeach[x][18] = TRUE;
    PercentFullSand[x][19] = 0.2;

    x = INIT_BEACH + 3;

    PercentFullSand[x][18] = 0.3;
  }
}

/**
 * Apply periodic boundary conditions to CEM arrays.
 * PARAMETERS: none
 * RETURN: none
 */
void periodic_boundary_copy(void)
{
  const int buff_len = 2 * Y_MAX;
  int row;

  for (row=0; row < X_MAX; row++)
  {
    apply_periodic_boundary(AllBeach[row], sizeof(char), buff_len);
    apply_periodic_boundary(AllRock[row], sizeof(char),  buff_len);
    apply_periodic_boundary(type_of_rock[row], sizeof(char), buff_len);
    apply_periodic_boundary(PercentFullSand[row], sizeof(double), buff_len);
    apply_periodic_boundary(PercentFullRock[row], sizeof(double),buff_len);
    apply_periodic_boundary(Age[row], sizeof(int), buff_len);
  }
}

/**
 * Resets all arrays recalculated at each time step to 'zero' conditions
 * PARAMETERS: none
 * RETURN: none
 */
void ZeroVars(void)
{
  int z;

  for (z = 0; z < MAX_BEACH_LENGTH; z++)
  {
    X[z] = -1;
    Y[z] = -1;
    XRock[z] = -1;
    YRock[z] = -1;
    XBehind[z] = -1; /*LMV*/ 
    YBehind[z] = -1; /*LMV*/ 
    XRockBehind[z] = -1;
    YRockBehind[z] = -1;
    InShadow[z] = FALSE;
    ShorelineAngle[z] = -999;
    UpWind[z] = '?';

    VolumeAcrossBorder[z] = 0.0; /*LMV*/ 
    ActualVolumeAcross[z] = 0.0; /*LMV*/
  }
}

/**
 * Reads saved output file, AllBeach[][] & PercentFullSand[][]
 * LMV - Also reads saved output file, AllRock[][], PercentFullRock[][], type_of_rock[][] 
 * PARAMETERS: none
 * RETURN: none
 */
void ReadSandFromFile(void)
{
  int x, y;

  ReadSandFile = fopen(readfilename, "r");
  printf("CHECK READ \n");
  /*line file input */
  if (SAVE_FILE != 2)
  {
    for (y = Y_MAX / 2; y < 3 * Y_MAX / 2; y++)
    {
      for (x = 0; x < X_MAX; x++)
      {
        fscanf(ReadSandFile, " %c", &type_of_rock[x][y]);
      }
    }

    for (y = Y_MAX / 2; y < 3 * Y_MAX / 2; y++)
    {
      for (x = 0; x < X_MAX; x++)
      {
        fscanf(ReadSandFile, " %lf", &PercentFullRock[x][y]);

        if (PercentFullRock[x][y] >= 1.0)
        {
          AllRock[x][y] = TRUE;
        }
        else
        {
          AllRock[x][y] = FALSE;
        }
      }
    }

    for (y = Y_MAX / 2; y < 3 * Y_MAX / 2; y++)
    {
      for (x = 0; x < X_MAX; x++)
      {
        fscanf(ReadSandFile, " %lf", &PercentFullSand[x][y]);

        if ((PercentFullSand[x][y] + PercentFullRock[x][y]) >= 1.0)
        {
          AllBeach[x][y] = TRUE;
        }
        else
        {
          AllBeach[x][y] = FALSE;
        }
      }
    }

    if (SAVE_AGE)
    {
      for (y = Y_MAX / 2; y < 3 * Y_MAX / 2; y++)
      {
        for (x = 0; x < X_MAX; x++)
        {
          fscanf(ReadSandFile, " %d", &Age[x][y]);
        }
      }
    }
  }
  /*Array file input */
  else
  {
    for (x = (X_MAX - 1); x >= 0; x--)
    {
      for (y = Y_MAX / 2; y < 3 * Y_MAX / 2; y++)
      {
        fscanf(ReadSandFile, " %c", &type_of_rock[x][y]);
      }
      fscanf(ReadSandFile, "\n");
    }

    for (x = (X_MAX - 1); x >= 0; x--)
    {
      for (y = Y_MAX / 2; y < 3 * Y_MAX / 2; y++)
      {
        fscanf(ReadSandFile, " %lf", &PercentFullRock[x][y]);

        if (PercentFullRock[x][y] >= 1.0)
        {
          AllRock[x][y] = TRUE;
        }
        else {
          AllRock[x][y] = FALSE;
        }
      }
      fscanf(ReadSandFile, "\n");
    }

    for (x = (X_MAX - 1); x >= 0; x--)
    {
      for (y = Y_MAX / 2; y < 3 * Y_MAX / 2; y++)
      {
        fscanf(ReadSandFile, " %lf", &PercentFullSand[x][y]);

        if ((PercentFullSand[x][y] + PercentFullRock[x][y]) >= 1.0)
        {
          AllBeach[x][y] = TRUE;
        }
        else
        {
          AllBeach[x][y] = FALSE;
        }
      }
      fscanf(ReadSandFile, "\n");
    }

    if (SAVE_AGE)
    {
      for (x = (X_MAX - 1); x >= 0; x--)
      {
        for (y = Y_MAX / 2; y < 3 * Y_MAX / 2; y++)
        {
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

/**
 * Saves current AllBeach[][] and PercentFullSand[][] data arrays to file
 * Save file name will add extension '.' and the current_time_step
 * LMV - Save AllRock[][], PercentFullRock[][], type_of_rock[][]
 * PARAMETERS: none
 * RETURN: none
 */
void SaveSandToFile(void)
{
  int x, y;
  char savename[40];
  printf("\n saving \n ");

  sprintf(savename, "%s_%d.out", savefilename, current_time_step);
  printf("Saving as: %s \n \n \n", savename);

  SaveSandFile = fopen(savename, "w");

  if (!SaveSandFile)
  {
    printf("problem opening output file\n");
    exit(1);
  }
  if (SAVE_FILE == 1)
  {
    for (y = Y_MAX / 2; y < 3 * Y_MAX / 2; y++)
    {
      for (x = 0; x < X_MAX; x++)
      {
        fprintf(SaveSandFile, " %c", type_of_rock[x][y]);
      }
    }
    for (y = Y_MAX / 2; y < 3 * Y_MAX / 2; y++) /* LMV */
    {
      for (x = 0; x < X_MAX; x++)
      {
        fprintf(SaveSandFile, " %lf", PercentFullRock[x][y]);
      }
    }
    for (y = Y_MAX / 2; y < 3 * Y_MAX / 2; y++) /* LMV */
    {
      for (x = 0; x < X_MAX; x++)
      {
        fprintf(SaveSandFile, " %lf", PercentFullSand[x][y]);
      }
    }

    if (SAVE_AGE)
    {
      for (y = Y_MAX / 2; y < 3 * Y_MAX / 2; y++)
      {
        for (x = 0; x < X_MAX; x++) fprintf(SaveSandFile, " %d", Age[x][y]);
      }
    }
  }

  /* Array output */
  if (SAVE_FILE == 2)
  {
    /* for (x=0; x<X_MAX; x++) */
    for (x = (X_MAX - 1); x >= 0; x--)
    {
      for (y = Y_MAX / 2; y < 3 * Y_MAX / 2; y++)
      {
        fprintf(SaveSandFile, " %c", type_of_rock[x][y]);
        /*LMV*/
        // if (type_of_rock[x][y]=='f')
        // {
        //   fprintf(SaveSandFile, " 0");  /* Switch on if numbers required*/
        // }
        // else
        // {
        //   fprintf(SaveSandFile, " 1");
        // }
      }
      fprintf(SaveSandFile, "\n");
    }

    for (x = (X_MAX - 1); x >= 0; x--)
    {
      for (y = Y_MAX / 2; y < 3 * Y_MAX / 2; y++) /* LMV */
      {
        fprintf(SaveSandFile, " %lf", PercentFullRock[x][y]);
      }
      fprintf(SaveSandFile, "\n");
    }

    for (x = (X_MAX - 1); x >= 0; x--)
    {
      for (y = Y_MAX / 2; y < 3 * Y_MAX / 2; y++)
      {
        fprintf(SaveSandFile, " %lf", PercentFullSand[x][y]);
      }
      fprintf(SaveSandFile, "\n");
    }

    if (SAVE_AGE)
    {
      for (x = (X_MAX - 1); x >= 0; x--)
      {
        for (y = Y_MAX / 2; y < 3 * Y_MAX / 2; y++)
        {
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

/**
 * Saves data line of shoreline position rather than entire array
 * Main concern is to have only one data point at each alongshore location
 * Save file name will add extension '.' and the current_time_step
 * PARAMETERS: none
 * RETURN: none
 */
void SaveLineToFile(void)
{
  int y, x, xtop, i;
  double xsave;
  char savename[40];

  printf("\n saving \n ");

  sprintf(savename, "ShorePos_%d.dat", current_time_step / SAVE_SPACING);
  printf("Saving as: %s                 ", savename);

  SaveSandFile = fopen(savename, "w");
  if (!SaveSandFile)
  {
    printf("problem opening output file\n");
    exit(1);
  }

  for (y = Y_MAX / 2; y < 3 * Y_MAX / 2; y++)
  {
    x = X_MAX - 1;
    xtop = X_MAX;

    /* step back to where we encounter allbeach */
    while (!AllBeach[x][y])
    {
      x -= 1;
    }

    /* if on side of shape, need to average */
    if (PercentFullSand[x + 2][y] > 0)
    {
      xtop = x + 1;
      while (PercentFullSand[xtop][y] > 0)
      {
        xtop += 1;
      }

      xsave = x;

      for (i = x + 1; i < xtop; i++)
      {
        xsave += PercentFullSand[i][y];
      }
    }
    /* otherwise Regular Beach Condition */
    else
    {
      xsave = x + PercentFullSand[x + 1][y];
    }

    /* note this assumes average of beach locations should be 0.5 percentfull */
    fprintf(SaveSandFile, " %f", xsave - INIT_BEACH + 0.5);

    // printf("y %d , xtop = %d xsave = %f \n",y,xtop,xsave);
    // if (xtop != X_MAX)
    // {
    //   PauseRun(x+1,y,-1);
    // }
  }

  fclose(SaveSandFile);
  printf("line file saved!\n\n");
}

/**
 * Pauses run until the 'q' key is pressed
 * Can Print or Plot Out Useful info
 * PARAMETERS: unused?
 * RETURN: none
 */
void PauseRun(int x, int y, int in)
{
  if (NO_PAUSE_RUN)
  {
    return;
  }
  printf("\nPaused \n");

  if (SAVE_LINE)
  {
    SaveLineToFile();
  }
  if (SAVE_FILE)
  {
    SaveSandToFile();
  }

  printf("\nend Pause\n");
}

/**
 * Age cells
 * PARAMETERS: none
 * RETURN: none
 */
void AgeCells(void)
{
  int x, y;
  for (y = 0; y < 2 * Y_MAX; y++)
  {
    for (x = 0; x < X_MAX; x++) {
      if (PercentFullSand[x][y] == 0)
      {
        Age[x][y] = current_time_step % AGE_MAX;
      }
    }
  }
}

/**
 * Input wave distribution
 * PARAMETERS: none
 * RETURN: none
 */
void ReadWaveIn(void)
{
  int i;
  for (i = 0; i <= 36; i++)
  {
    WaveMax[i] = 0;
    WaveProb[i] = 0;
  }

  ReadWaveFile = fopen(readwavename, "r");
  printf("CHECK READ WAVE\n");
  fscanf(ReadWaveFile, " %d \n", &NumWaveBins);
  printf("Wave Bins %d \n", NumWaveBins);

  WaveMax[0] = -90;
  WaveProb[0] = 0;

  for (i = 1; i <= NumWaveBins; i++)
  {
    fscanf(ReadWaveFile, " %lf %lf", &WaveMax[i], &WaveProb[i]);
    printf("i= %d  Wave= %f Prob= %f \n", i, WaveMax[i], WaveProb[i]);
  }

  fclose(ReadWaveFile);
  printf("wave file read! \n");
}

/**
 * Input wave distribution
 * PARAMETERS: none
 * RETURN: none
 */
void RealWaveIn(void)
{
  int countout = 10;

  WaveAngleIn = -999;
  WavePeriodIn = -999;
  WaveHeightIn = -999;

  if (current_time_step == 0)
  {
    ReadRealWaveFile = fopen(readwavename, "r");
    printf("CHECK READ WAVE: %s\n", readwavename);
  }

  fscanf(ReadRealWaveFile, "%lf %lf %lf\n", &WaveAngleIn, &WavePeriodIn, &WaveHeightIn);
  /* Recycling wave data */
  while (WaveAngleIn == -999 || WavePeriodIn == -999 || WaveHeightIn == -999)
  {
    countout--;
    printf("Recycling wave data.\n If this happens every timestep then there maybe a problem reading the wave data?\n");
    fclose(ReadRealWaveFile);
    ReadRealWaveFile = fopen(readwavename, "r");
    printf("CHECK READ WAVE: %s\n", readwavename);
    fscanf(ReadRealWaveFile, "%lf %lf %lf\n", &WaveAngleIn, &WavePeriodIn, &WaveHeightIn);
    printf("Angle= %lf  Period= %lf Height= %lf \n", WaveAngleIn, WavePeriodIn, WaveHeightIn);
    if (countout == 0)
    {
      exit(0);
    }
  }

  if (debug18)
  {
    printf("Angle= %lf  Period= %lf Height= %lf \n", WaveAngleIn, WavePeriodIn, WaveHeightIn);
  }

  if (current_time_step == stop_after)
  {
    fclose(ReadRealWaveFile); /* Keep file open till end of run */
  }
}

/**
 * Output wave distribution
 * PARAMETERS: none
 * RETURN: none
 */
void WaveOutFile(void)
{
  if (current_time_step == 0)
  {
    WaveFileOutput = fopen(wavesavename, "w");
  }

  fprintf(WaveFileOutput, "%lf %lf %lf\n", (-WaveAngle * RAD_TO_DEG), Period, OffShoreWvHt);

  if (current_time_step == stop_after)
  {
    fclose(WaveFileOutput); /*Keep file open till end of run */
  }
}

/**
 * Read control file
 * PARAMETERS: none
 * RETURN: none
 */
void ControlFile(void)
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

  if (time_step == -999)
  {
    printf("***Initialisation file not read!***\n");
  }

  if (METADATA)
  {
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
}

/**
 * Find index of x and y coordinates
 * don't use this during normal model run, if possible--it's too slow. as of
 * 10/05 it's only used during a pause
 * PARAMETERS: x, y: coordinates of cell, interface: sand or rock
 * RETURN: index
 */
int getIndex(int x, int y, char interface)
{
  int i;
  if (interface == 's')
  {
    if (TotalBeachCells == 0) 
    {
      // printf("TotalBeachCells == 0 (probably beach ");
      // printf("stuff hasn't been done yet) so no index\n");
      return -1;
    }
    for (i = 0; i < TotalBeachCells; i++)
    {
      if (X[i] == x && Y[i] == y)
      {
        // printf("getIndex (s) found index to be %d \n", i);
        return i;
      }
    }
    // printf("getIndex (s) didn't find one\n");
    return -1;
  }
  else if (interface == 'r')
  {
    if (TotalRockCells == 0)
    {
      // printf("TotalRockCells == 0 (probably beach ");
      // printf("stuff hasn't been done yet) so no index\n");
      return -1;
    }
    for (i = 0; i < TotalRockCells; i++)
    {
      if (XRock[i] == x && YRock[i] == y)
      {
        // printf("getIndex (r) found index to be %d \n", i);
        return i;
      }
    }
    // printf("getIndex (r) didn't find one\n");
    return -1;
  }
  else
  {
    printf("have to call getIndex for sand 's' or rock 'r'\n");
    PauseRun(-1, -1, -1);
    return -1;
  }
}

// TODO: clean up these prints
/**
 * Print local array conditions around x, y
 * PARAMETERS: x, y: coordinates of cell
 * RETURN: none
 */
void PrintLocalConds(int x, int y, int in)
{
  int i, j, k, isee;
  double vol = CELL_WIDTH * CELL_WIDTH * DEPTH_SHOREFACE;

  printf("\n x: %d  y: %d  z: %d\n\n", x, y, in);

  /* if not given location along beach, look to see if along beach */
  if (in < 0)
  {
    for (i = 0; i <= TotalBeachCells; i++)
    {
      if ((X[i] == x) && (Y[i] == y))
      {
        isee = i;
      }
    }
  }
  else
  {
    isee = in;
  }

  for (i = x + 2; i > x - 3; i--)
  {
    for (j = y - 2; j < y + 3; j++)
    {
      printf("  	%d,%d", i, j);
    }
    printf("\n");
  }

  printf("\n");

  printf("AllBeach\n");
  for (i = x + 2; i > x - 3; i--)
  {
    for (j = y - 2; j < y + 3; j++)
    {
      printf("	%c", AllBeach[i][j]);
    }
    printf("\n");
  }

  printf("\n");

  printf("AllRock\n"); /* LMV */
  for (i = x + 2; i > x - 3; i--)
  {
    for (j = y - 2; j < y + 3; j++)
    {
      printf("	%c", AllRock[i][j]);
    }
    printf("\n");
  }

  printf("\n");

  printf("PercentFullSand\n");
  for (i = x + 2; i > x - 3; i--)
  {
    for (j = y - 2; j < y + 3; j++)
    {
      printf("	%2.5f", PercentFullSand[i][j]);
    }
    printf("\n");
  }

  printf("\n");

  printf("PercentFullRock\n"); /* LMV */
  for (i = x + 2; i > x - 3; i--)
  {
    for (j = y - 2; j < y + 3; j++)
    {
      printf("	%2.5f", PercentFullRock[i][j]);
    }
    printf("\n");
  }

  printf("\n");
  printf("current_time_step: %d\n", current_time_step);
  printf("InitialMass: %f\n", MassInitial); /* LMV */
  printf("CurrentMass: %f\n", MassCurrent); /* LMV */
  printf(" %d  ", in); /* LMV */

  if (isee >= 0)
  {
    for (k = in - 3; k < in + 4; k++)
    {
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

/***************************************************************************
Over wash functions, copied from Kenny's code (so from Andrew Ashton's code)
***************************************************************************/

/**
 * Just a loop to call overwash check founction CheckOverwash
 * Nothing done here, but can be down when CheckOVerwash is called
 * PARAMETERS: none
 * RETURN: none
 */
void CheckOverwashSweep(void)
{
  int i, ii;
  int sweepsign;
  // TODO: pull out common code
  sweepsign = RandZeroToOne() * 2 > 1 ? 1 : 0;
  OWflag = FALSE;

  for (i = 1; i < TotalBeachCells - 1; i++)
  {
    ii = sweepsign ? i : TotalBeachCells - 1 -i;

    /* TODO: test shoreline should be facing seaward */
    /* don't worry about shadow here, as overwash is not set to a time scale
     * with AST       */

    if ((fabs(SurroundingAngle[ii]) < (OVERWASH_LIMT / RAD_TO_DEG)) && !InShadow[ii])
    {
      CheckOverwash(ii);
    }
  }
}

/* Step back pixelwise in direction of SurroundingAngle[] to check 'needage' (?)
 * cell being checked for overwash
 * Calls DoOverwash(). 'x' and 'y' hold real-space (doubleing point) values, and
 * will be mapped onto integer array
 * PARAMETERS: icehckL index of the shoreline
 * RETURN: none
 */
void CheckOverwash(int icheck)
{
  double slope;               /* slope of zero goes staight back */
  int ysign;                  /* holder for going left or right alongshore */
  double x, y;                /* holders for 'real' location of x and y */
  double xin, yin;            /* starting 'real' locations */
  int xtest, ytest;           /* cell looking at */
  double xint, yint;          /* intercepts of overwash line in overwashable cell */
  int NextXInt, NextYInt;     /* holder vairables for cell to check */
  double Ydown, DistanceDown; /* when going to next x cell, what other values */
  double Xside, DistanceSide; /* when gpoing to next y cell,other values */
  double checkdistance;       /* distance of checking line- minimum, not actual width, ends loop */
  double measwidth;           /* actual barrier width between cells */
  int AllBeachFlag;           /* flag to see if overwash line has passed over at least one AllBeach cell */

  /* convert angle to a slope and the direction of steps */
  /* note that for case of shoreline, positive angle will be minus y direction */
  if (SurroundingAngle[icheck] == 0.0)
  {
    /* unlikely, but make sure no div by zero */
    slope = 0.00001;
  }
  else if (fabs(SurroundingAngle[icheck]) == 90.0)
  {
    slope = 9999.9;
  }
  else
  {
    slope = fabs(tan(SurroundingAngle[icheck]));
  }
  ysign = SurroundingAngle[icheck] > 0 ? 1 : -1;

  /* 'regular condition' */
  /* plus 'stuck in the middle' situation (unlikely scenario) */
  /* if the cell below is all beach */
  /* or if the cell to the left is all beach and the cell to the right is all beach */
  if (AllBeach[X[icheck] - 1][Y[icheck]] || (AllBeach[X[icheck]][Y[icheck]] && AllBeach[X[icheck]][Y[icheck] + 1])) 
  {
    xin = X[icheck] + PercentFullSand[X[icheck]][Y[icheck]];
    yin = Y[icheck] + 0.5;
  }
  /* if the cell to the left is all beach on right side */
  else if (AllBeach[X[icheck]][Y[icheck] - 1]) 
  {
    xin = X[icheck] + 0.5;
    yin = Y[icheck] + PercentFullSand[X[icheck]][Y[icheck]];
  }
  /* on left side */
  else if (AllBeach[X[icheck]][Y[icheck] + 1])
  {
    xin = X[icheck] + 0.5;
    yin = Y[icheck] + 1.0 - PercentFullSand[X[icheck]][Y[icheck]];
  }
  /* underneath, no overwash */
  else
  {
    return;
  }

  x = xin;
  y = yin;
  checkdistance = 0;
  AllBeachFlag = FALSE;

  while ((checkdistance < CRIT_BEACH_WIDTH) && (y > 0) && (y < 2 * Y_MAX) && (x > 1))
  {
    NextXInt = ceil(x) - 1;
    NextYInt = ysign > 0 ? floor(y) + 1 : ceil(y -1);

    /* moving to next whole 'x' position, what is y position? */
    Ydown = y + (x - NextXInt) * slope * ysign;
    DistanceDown = Raise(((Ydown - y) * (Ydown - y) + (NextXInt - x) * (NextXInt - x)), .5);

    /* moving to next whole 'y' position, what is x position? */
    Xside = x - fabs(NextYInt - y) / slope;
    DistanceSide = Raise(((NextYInt - y) * (NextYInt - y) + (Xside - x) * (Xside - x)), .5);

    /* next cell is the down cell */
    if (DistanceDown < DistanceSide)
    {
      x = NextXInt;
      y = Ydown;
      xtest = NextXInt - 1;
      ytest = floor(y);
    }
    /* next cell is the side cell */
    else
    {
      x = Xside;
      y = NextYInt;
      xtest = floor(x);
      ytest = y + (ysign - 1) / 2;
    }

    checkdistance = Raise(((x - xin) * (x - xin) + (y - yin) * (y - yin)), .5) * CELL_WIDTH;

    if (AllBeach[xtest][ytest])
    {
      AllBeachFlag = TRUE;
    }

    /* if passed through an allbeach and a neighboring partial cell, jump out, only bad things follow */
    if (!AllBeach[xtest][ytest] && AllBeachFlag && !(((X[icheck] - xtest) > 1) || (abs(ytest - Y[icheck]) > 1)))
    {
      return;
    }

    /* Looking for shore cells, but don't want immediate neighbors, and go
       backwards */
    /* Also mush pass though an allbeach cell along the way */
    if (!AllBeach[xtest][ytest] && AllBeachFlag && (xtest < X[icheck]) && (((X[icheck] - xtest) > 1) || (abs(ytest - Y[icheck]) > 1)))
    {
      /* 'regular condition' - UNDERNEATH, here */
      if (AllBeach[xtest + 1][ytest])
      {
        xint = (xtest + 1 - PercentFullSand[xtest][ytest]);
        yint = yin + (xin - xint) * ysign * slope;

        /* This cell isn't actually an overwash cell */
        if ((yint > ytest + 1.0) || (yint < ytest))
        {
          measwidth = CRIT_BEACH_WIDTH;
        }
        else
        {
          measwidth = CELL_WIDTH * Raise((xint - xin) * (xint - xin) + (yint - yin) * (yint - yin), 0.5);
        }
      }
      /* on right side */
      else if (AllBeach[xtest][ytest - 1])
      {
        yint = (ytest + PercentFullSand[xtest][ytest]);
        xint = xin - fabs(yin - yint) / slope;

        /* This cell isn't actually an overwash cell */
        if (xint < xtest)
        {
          measwidth = CRIT_BEACH_WIDTH;
        }
        else
        {
          measwidth = CELL_WIDTH * Raise((xint - xin) * (xint - xin) + (yint - yin) * (yint - yin), 0.5);
        }
      }
      /* on left side */
      else if (AllBeach[xtest][ytest + 1])
      {
        yint = (ytest + 1 - PercentFullSand[xtest][ytest]);
        xint = xin - fabs(yin - yint) / slope;

        /* This cell isn't actually an overwash cell */
        if (xint < xtest)
        {
          measwidth = CRIT_BEACH_WIDTH;
        }
      }
      /* 'regular condition' */
      /* plus 'stuck in the middle' situation */
      else if (AllBeach[xtest - 1][ytest])
      {
        xint = (xtest + PercentFullSand[xtest][ytest]);
        yint = yin + (xin - xint) * ysign * slope;

        /* This cell isn't actually an overwash cell */
        if ((yint > ytest + 1.0) || (yint < ytest))
        {
          measwidth = CRIT_BEACH_WIDTH;
        }
        else
        {
          measwidth = CELL_WIDTH * Raise((xint - xin) * (xint - xin) + (yint - yin) * (yint - yin), 0.5);
        }
      }
      /* uh oh - not good situation, no allbeach on sides */
      /* assume this is an empty cell, */
      else if (PercentFullSand[xtest][ytest] > 0)
      {
        xint = x;
        yint = y;

        measwidth = CELL_WIDTH * Raise((xint - xin) * (xint - xin) + (yint - yin) * (yint - yin), 0.5);

        printf(
            "-- Some Odd Over  xin: %2.2f  yin: %2.2f xt:%d yt: %d xint: %f "
            "yint: %f Meas: %3.2f Ang: %f Abs: %f\n",
            xin, yin, xtest, ytest, xint, yint, measwidth,
            SurroundingAngle[icheck] * RAD_TO_DEG,
            fabs(SurroundingAngle[icheck]) * RAD_TO_DEG);
      }
      /* empty cell - oughta fill er up  - fill max barrier width */
      else
      {
        xint = x;
        yint = y;
        measwidth = CRIT_BEACH_WIDTH - CELL_WIDTH;

        printf(
            "-- Empty Odd Over  xin: %2.2f  yin: %2.2f xt:%d yt: %d xint: %f "
            "yint: %f Meas: %3.2f Ang: %f Abs: %f\n",
            xin, yin, xtest, ytest, xint, yint, measwidth,
            SurroundingAngle[icheck] * RAD_TO_DEG,
            fabs(SurroundingAngle[icheck]) * RAD_TO_DEG);
      }

      checkdistance = measwidth;

      if (measwidth < CRIT_BEACH_WIDTH)
      {
        DoOverwash(X[icheck], Y[icheck], xtest, ytest, xint, yint, measwidth, icheck);
        /* jump out of loop */
        OWflag = TRUE;
        return;
      }
    }
  }
}

/*
 *  Given a cell where overwash is needed, move sediment back
 *  for 'true' overwash based on shoreline angles
 * PARAMETERS: xfrom, yfrom: coordinates of cell sediment is moving out of,
 * xto, yto: coordinates of cell sediment is moving into, 
 * xintto, yintto: ???
 * widthin: ???
 * ishore: index of the shoreline
 * RETURN: none, write to PrecentFull[][] and AllBeach[][]
 */
void DoOverwash(int xfrom, int yfrom, int xto, int yto, double xintto, double yintto, double widthin, int ishore)
{
  double BBneed, delBB, delShore; /* local variables */

  // double DepthBackBarrier = 6.0;      /* m current set depth for backbarrier (temp - make into function) */

  double DepthBB; /* holds effective backbarrier depth */

  DepthBB = GetOverwashDepth(xto, yto, xintto, yintto, ishore);

  /* calculated value of most that backbarrier can move given geometry (true, non-iterative solution) */

  if (DepthBB == DEPTH_SHOREFACE)
  {
    BBneed = MAX_OVER;
  }
  else
  {
    BBneed = (CRIT_BEACH_WIDTH - widthin) / CELL_WIDTH / (1 - (DepthBB / DEPTH_SHOREFACE));
  }

  /* do all overwash */
  if (BBneed <= MAX_OVER)
  {
    delShore = BBneed * DepthBB / DEPTH_SHOREFACE;
    delBB = BBneed;
  }
  /* only do overwash to max change) */
  else
  {
    delShore = MAX_OVER * DepthBB / DEPTH_SHOREFACE;
    delBB = MAX_OVER;
  }

  PercentFullSand[xto][yto] += delBB;
  PercentFullSand[xfrom][yfrom] -= delShore;

  if (PercentFullSand[xto][yto] > 1)
  {
    OopsImFull(xto, yto);
  }
  if (PercentFullSand[xfrom][yfrom] < 0)
  {
    OopsImEmpty(xfrom, yfrom);
  }
}

/*
 * Finds corresponding overwash depths
 * OWType = 0 take the depth at neightbor to the backing cell
 * OWType = 1 geometric rule based upon distance from back to shoreline
 * PARAMETERS: xin, yin: coordinates of cell sediment overwashing into
 * xinfl, yinfl: ???
 * ishore: index of the shoreline
 * RETURN: depth of overwash
 */
double GetOverwashDepth(int xin, int yin, double xinfl, double yinfl, int ishore)
{
  int xdepth;
  double Depth;
  double BBDistance;            /* Distance from backshore to next shore */
  double slope;                 /* slope of zero goes staight back */
  int ysign;                    /* holder for going left or right alongshore */
  double x, y;                  /* holders for 'real' location of x and y */
  int xtest, ytest;             /* cell looking at */
  int NextXInt, NextYInt;       /* holder variables for cell to check */
  double Ydown, DistanceDown;   /* when going to next x cell, what other values */
  double Xside, DistanceSide;   /* when gpoing to next y cell,other values */
  int BackFlag;                 /* Flag to indicate if hit backbarrier */
  int Backi = -1;               /* i for backbarrier intersection */
  int i, j;                     /* counters */
  int FoundFlag;                /* Backbarrier intersection flag */
  double AngleSin, AngleUsed;

  /* Use Cell Depths for overwash depths */
  if (OVERWASH_TYPE == 0)
  {
    xdepth = xin;
    Depth = cell_depth[xdepth][yin];

    while ((Depth < 0) && (xdepth > 0))
    {
      Depth = cell_depth[xdepth][yin];
      xdepth--;
    }

    if (Depth == DEPTH_SHOREFACE)
    {
      Depth = 6.0;
    }

    return Depth;
  }
  /* Geometric relation to determine depth through intersection of shorefaces */
  /* look in line determined by shoreline slope - reuse stepping function (again) */
  else if (OVERWASH_TYPE == 1)
  {
    x = xinfl;
    y = yinfl;

    /* unlikely, but make sure no div by zero */
    if (SurroundingAngle[ishore] == 0.0)
    {
      slope = 0.00001;
    }
    else if (fabs(SurroundingAngle[ishore]) == 90.0)
    {
      slope = 9999.9;
    }
    else
    {
      slope = fabs(tan(SurroundingAngle[ishore]));
    }

    BackFlag = FALSE;
    ysign = SurroundingAngle[ishore] > 0 ? 1 : -1;

    while ((!BackFlag) && (y > 0) && (y < 2 * Y_MAX) && (x > 1))
    {
      NextXInt = ceil(x) - 1;
      NextYInt = ysign > 0 ? floor(y) + 1 : ceil(y - 1);

      /* moving to next whole 'x' position, what is y position? */
      Ydown = y + (x - NextXInt) * slope * ysign;
      DistanceDown = Raise(((Ydown - y) * (Ydown - y) + (NextXInt - x) * (NextXInt - x)), .5);

      /* moving to next whole 'y' position, what is x position? */
      Xside = x - fabs(NextYInt - y) / slope;
      DistanceSide = Raise(((NextYInt - y) * (NextYInt - y) + (Xside - x) * (Xside - x)), .5);

      /* next cell is the down cell */
      if (DistanceDown < DistanceSide)
      {
        x = NextXInt;
        y = Ydown;
        xtest = NextXInt - 1;
        ytest = floor(y);
      }
      /* next cell is the side cell */
      else
      {
        x = Xside;
        y = NextYInt;
        xtest = floor(x);
        ytest = y + (ysign - 1) / 2;
      }

      if (PercentFullSand[xtest][ytest] > 0)
      {
        BackFlag = TRUE;
      }
    }

    /* Try to find the i for the cell found */
    /* If you have a better idea how to do this, go ahead */
    i = 2;
    FoundFlag = FALSE;

    while ((i < TotalBeachCells - 1) && !(FoundFlag))
    {
      if ((X[i] == xtest) && (Y[i] == ytest))
      {
        FoundFlag = TRUE;
        Backi = i;
      }
      i++;
    }

    /* The search for the backbarrier went out of bounds - not good, assume big = depthshoreface    */
    /* Periodic B.C.'s should make this not so important */
    if (!BackFlag)
    {
      Depth = DEPTH_SHOREFACE;
    }
    else
    {
      BBDistance = Raise(((xinfl - xtest - PercentFullSand[xtest][ytest]) * (xinfl - xtest - PercentFullSand[xtest][ytest])) + ((yinfl - ytest - 0.5) * (yinfl - ytest - 0.5)), .5);

      /* The backbarrier intersection isn't on the shoreline */
      /* Assume 1/2 of the length applies to this case */
      if (!FoundFlag)
      {
        Depth = BBDistance / 2 * SHOREFACE_SLOPE * CELL_WIDTH;
      }
      /* Use the fancy geometry thing */
      else
      {
        AngleUsed = 0;
        for (j = -1; j < 2; j++)
        {
          AngleUsed += SurroundingAngle[Backi + j];
        }
        AngleUsed = AngleUsed / 5;

        if (fabs(AngleUsed) > PI / 4.0)
        {
          AngleUsed = PI / 4.0;
        }

        AngleSin = sin(PI / 2.0 - fabs(SurroundingAngle[ishore] + AngleUsed));

        Depth = BBDistance * AngleSin / (1 + AngleSin);
      }
    }

    if (Depth < OVERWASH_MIN_DEPTH)
    {
      Depth = OVERWASH_MIN_DEPTH;
    }
    else if (Depth > DEPTH_SHOREFACE)
    {
      Depth = DEPTH_SHOREFACE;
    }
    return Depth;
  }

  printf("OWDepth all broken");
  PauseRun(xin, yin, -1);
  return DEPTH_SHOREFACE;
}

/*
 * This function takes
 * 1) the EvaluateAngle defined in CEM and returns an azimuth, or
 * 2) the azimuthal angle from SWAN output and returns an angle in the CEM convention.
 * PARAMETERS: 'type' separates the two actions ('type = 1' triggers the first statement
 * type = anything except 1' triggers the second statement) #SWAN 
 * RETURNS: wave angle in the desired format
 */
 // TODO: move to utils?
double ConvertAngle(double EvaluateAngle, int type)
{
  double value; /* What direction to look offshore to find a breaking wave? */
  // if (debugSWAN)
  // {
  //    printf("\n incoming angle: %f \n", EvaluateAngle);
  // }

  if (type == 1)
  {
    /* Moving right and down along shoreline -- most common case, unless we have spits or cape tip overgrowth */
    if (EvaluateAngle > (-90 * (PI / 180)) && EvaluateAngle < (0 * (PI / 180)))
    {
      value = EvaluateAngle + (90 * (PI / 180));
    }

    /* Moving right and up (or just right) along shoreline -- most common case, unless we have spits or cape tip overgrowth */
    else if (EvaluateAngle > (0 * (PI / 180) && EvaluateAngle < (90 * (PI / 180)))) 
    {
      value = (360 * (PI / 180)) - EvaluateAngle;
    }

    /* Moving right, stright shoreline (dx = 0) */
    else if (EvaluateAngle == 0)
    {
      value = 0;
    }

    /* Moving straight up or straight down along shoreline */
    else if (EvaluateAngle == (-90 * (PI / 180)) || EvaluateAngle == (90 * (PI / 180)))
    {
      if (EvaluateAngle == (-90 * (PI / 180)))
      {
        value = (270 * (PI / 180));
      }
      else if (EvaluateAngle == (90 * (PI / 180)))
      {
        value = (90 * (PI / 180));
      }
    }

    /* Moving left (and up, down, or straight) along shoreline */
    else if (EvaluateAngle > (-270 * (PI / 180) && EvaluateAngle < (-90 * (PI / 180))))
    {
      value = fabs(EvaluateAngle);
    }

    /* Oops -- did not find a search direction! */
    else
    {
      printf("\n Can't look offshore for SWANs 'cause i didn't find a surrounding angle! (in ParseSWAN) \n");
    }

    /* Print to screen, if you like */
    // if (debugSWAN)
    // {
    //   printf("\n search direction: %f , surrounding angle: %f \n", value, EvaluateAngle);
    // }

    return value;
  }

  /* -----> Types other than 1 <----- */
  /* Convert SWAN azimuthal angle */
  else
  {
    /* Wave going right to left */
    if (EvaluateAngle > 0 && EvaluateAngle < 90)
    {
      return -EvaluateAngle * (PI / 180);
    }
    /* Wave going left to right */
    else if (EvaluateAngle < 360 && EvaluateAngle > 270)
    {
      return (360 - EvaluateAngle) * (PI / 180);
    }
    /* Wave coming straight onshore */
    else if (EvaluateAngle == 0 || EvaluateAngle == 360)
    {
      return 0;
    }
  }
  return 0.;
}

/*
 * SWAN data parse function! PWL, 10-18-13. #SWAN
 * PARAMETERS: ShoreAngleLoc: alongshore location ShoreAngle: angle of the shore
 * RETURNS: none
 */
void ParseSWAN(int ShoreAngleLoc, double ShoreAngle)
{
  // double LookOffshore;              /* What direction to look offshore to find a breaking wave? */
  // double	LookSlope;                 /* What slope (x,y) is the line of sight? */
  int xcoord, ycoord;
  // int NextInLine;
  // int OopsImBroke = 0;              /* Flag that indicates if we've found a broken wave. 0 = keep lookin', pal; 1 = eureka! */
  // double	interval = 0.1;            /* How far to search in each iteration */
  // double	xdist = 0.0, ydist = 0.0;  /* Distance to search for SWAN info (units are cells). */
  // double	Xpos, Ypos;                /* Tracks search distance */

  /* Tracks cells that have already been searched */
  // int	LastYCell; */
  int LastXCell;

  // int CurrentXCell, CurrentYCell;   /* Tracks integer value cell position */
  double Hd;                           /* Wave height divided by water depth */
  double CAngle;                       /* Angle converted from ConvertAngle function -- extraneous for dubugging */
  // double	RAngle;                    /* Angle retrieved from SWAN -- extraneous for dubugging */
  int x;                               /* iterator for simple search function, 11-12-13 -- can delete for shoreline angle-based version */

  xcoord = X[ShoreAngleLoc];
  ycoord = Y[ShoreAngleLoc];

  if (ShoreAngleLoc == 0)
  {
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

      if (Hd < WAVE_BREAK_DEPTH && Hd > 0 && shelf_depth[LastXCell][ycoord] > 0) {
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
        CAngle = ConvertAngle(EvaluateAngle, 2); /* Convert angle from azimuth (degrees) to CEM-style (rads) */
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

        break; /* delete break and use OopsImBroke with angle-based function, 11-12-13 */
      }
    }
    else
    {
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
