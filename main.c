/*
CEM_good
Base Version; includes Lisa Valvo's code, et al.
Supplied by Brad Murray, November 2013.

Add rock (PercentFull changed to PercentFullSand and PercentFullRock

Program To generate capes and sandwaves using wave angle relationships
Begun by Brad Murray 01/00
Refined by Olivier Arnoult 01/00 - 06/00
Revised and reformed by Andrew Ashton 06/00

This version includes:
Overwash routines:
Taken from Kenny Ells and Andrew Ashton's codes
Embedded by Chris Thomas (CWT). December 2013.
---Chris Thomas (CWT) Changes Dec 2013---
1.  - Disabled all of the embedded graphical functions as these are
platform/compiler specific

2.  - Overwash routines added
        -Overwash defines and variables copied in from Kenny's code, 24/4/13
        -#defines for the overwash functions
        -Overwash functions at lines 4624 onwards
        -#ifdef block is from another version of the CEM
        -Some of Jordan Slott's code.
        -Some of the overwash defs are required by the SedTrans() function


---Andy Barkwith (AB) Changes Jan 2014---
1. - Random value selection changed to something more robust
        -srandom() and random() replaced with srand() and rand() as these were
not recognised be the compiler
        -truly random numbers can be used by setting seed = -999

2. - Read in real wave data
        -added a third option (2) for WaveIn that reads Angle Period Height from
the text file defined in 'readwavename'
        -if the end of the file is reached before the end of the end of the run
the data is re-read
        -any wave >90 || <-90 is given an angle of 0 deg and a height of 0 m

3. - Read in model control data
        -allows model attributes to be changed without recompiling using
InitialiseFile and defining 'readcontrolname'

4. - Model attribute metadata
        -allows model attributes to be output into a file (only allowed if
control data is read from file)
        -allows wave climate used by the model to be output

5. - Wave height, rotation and period factors
        - Wave height, rotation and period data read from a file may be modified
by a set amount.
        - A wave rotation (deg) can be applied if the coastline has been rotated
for the model
        - Wave height and period may be multiplied by a factor (1 = no change, 2
= double and 0.5 = half)

6. - Array IO
        - Can output data as an array instead of a line by setting SaveFile = 2
        - Array is actually 3 or 4 arrays stacked on top of each other -->
Erodability -  Rock Percentage - Sand Percentage - Age (if selected)
        - If SaveFile=2, input rock/sand/sea must be read in as an array with
the same format as the output files
        - North is at the top of the IO arrays (useful to know if reading in
wave climate data)


Program Notes:
To end program, press 'd' key and 'ESC' key simultaneously
To save current iteration to file, press 's' and 'f' simultaneously
To update screen display, press 'p' key

THIS PROGRAM GONNA MAKE CAPES, SANDWAVES??
*/

/* standard includes */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

/*
CWT
Graphics and related libraries not included at present
#include <GL/gl.h>
#include <GL/glx.h>
#include <ncurses.h>
*/

/* Aspect Parameters */
#define CellWidth 100 /* size of cells (meters) */
#define Xmax 50       /* number of cells in x (cross-shore) direction */
#define Ymax 200      /* number of cells in y (longshore) direction */
#define MaxBeachLength                                                      \
  8 * Ymax /* maximum length of arrays that contain beach data at each time \
              step */
#define ShelfSlope 0.001    /* slope of continental shelf */
#define ShorefaceSlope 0.01 /* for now, linear slope of shoreface */
#define DepthShoreface \
  10 /* minimum depth of shoreface due to wave action (meters) */
#define InitBeach \
  20 /* cell where intial conditions changes from beach to ocean */
#define InitRock \
  5 /* cell where initial conditions change from beach to rock LMV*/
#define InitialDepth \
  10 /* theoretical depth in meters of continental shelf at x = InitBeach */
#define FindCellError \
  5 /* if we run off of array, how far over do we try again? */
#define ShadowStepDistance 0.2 /* step size for shadow cell checking */

/*Overwash Parameters*/
#define CritBWidth                                                          \
  350.0 /* Overwash - width barrier maintains due to overwash (m) important \
           scaling parameter! */
#define InitBWidth 4 /* Overwash - initial minimum width of barrier (Cells) */
#define OWType 1 /* Overwash - 0 = use depth array, 1 = use geometric rule */
#define OWMinDepth 1.0 /* Overwash - littlest overwash of all */

int OWflag = 0; /* A debugging flag for overwash routines */
float MaxOver =
    0.01; /* Maximum overwash step size (enforced at back barrier) */
float OverwashLimit = 60; /* Don't do over wash if the angle is > 60 degrees */
float CellDepth[Xmax][2 * Ymax]; /* Depth array */

/*
From Jordan Slott
Cuts out the TRUE & FALSE undeclared compiler errors in SedTrans()
*/

#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif

/*  Run Control Parameters */
#define HaveSinks 0 /* include sediment sinks in model run? */
#define Sinkiness                                                        \
  0.5 /* fractional value determines chance that sediment in the sink is \
         deleted */
double TimeStep =
    1; /* days - reflects rate of sediment transport per time step */
double OffShoreWvHt = 2.0; /* meters */
double Period = 10.0;      /* seconds */
double Asym = 0.70; /* fractional portion of waves coming from positive (left)
                       direction */
double Highness = 0.35; /* .5 = even dist, < .5 high angle domination. NOTE
                           Highness actually determines lowness! */
double Duration =
    1; /* Number of time steps calculations loop at same wave angle */
double StopAfter = 36500; /* Stop after what number of time steps */
#define NumberChunk                                                    \
  9 /* Number of chunks of rock in alongshore direction with different \
       weathering rates LMV */
double SlowWeatherCoeff = 0;   /* Weathering rate of slow rock e.g. 0.5 = slow
                                  weathering, 1/2 as fast LMV */
double FastWeatherCoeff = 1;   /* Weathering rate of fast rock e.g. 2 = fast
                                  weathering, 2 times faster than normal LMV */
double PercentFineFast = 0.00; /* Percent of fast weathering rock lost because
                                  it is too fine to stay in nearshore LMV */
double PercentFineSlow = 0.00; /* Percent of slow weathering rock lost because
                                  it is too fine to stay in nearshore LMV */
double ErosionRatePerYear =
    0.5; /* Amount of erosion to TotalBeachCells in m/year LMV*/
#define DoGraphics                                                           \
  =                                                                          \
      0; /* CWT Re-cast as a define, rather than the odd char DoGraphics = 1 \
            in the original; 0 = false */
double coastrotation =
    40; /* Angle (deg) used to align coastline with real wave climate*/
double waveheightchange =
    1; /* Factor applied to the input wave height; 1 = no change, 0.5 = half
          wave height, and 2 = double wave height*/
double waveperiodchange =
    1; /* Factor applied to the input wave period; 1 = no change, 0.5 = half
          wave period, and 2 = double wave period*/

int InitCType = 0; /* 0: normal (columns/blocks), 1: wiggly, 2: one block */
int seed = 1; /* random seed:  control value = 1 completely random = -999 */
int StartSavingAt = 0; /* time step to begin saving files */
int SaveSpacing = 365; /* space between saved files */
char savefilename[24] = "CEM";
char StartFromFile = 'n'; /* start from saved file? */
char readfilename[24] = "CEM_3285.out";
int WaveIn = 0; /* Input Wave Distribution file: no = 0, binned file = 1,
                   angle/period/height file =2 */
char readwavename[24] = "In_WaveData.dat";
char InitialiseFile =
    'n'; /* use a file to initialise run? Over rides setup data*/
char readcontrolname[24] = "In_CEM_init.dat";
int CurrentTimeStep = 0; /* Time step of current calculation */
char Metadata = 'n';     /*Create a metadata file?*/
char metasavename[24] = "Metadata.out";
char Wavedata = 'n'; /*Create a wavedata file?*/
char wavesavename[24] = "Wavedata.out";

int SaveFile = 2;           /* 1 = line output, 2 = array output*/
int SaveAge = 1;            /* Save/update age of cells? */
int SaveLine = 1;           /* Save line instead of whole array? */
char PromptStart = 'n';     /* ask prompts for file names, etc? */
char OffArray = 'n';        /* Initializing this variable for later use */
int ScreenTextSpacing = 30; /* Spacing of writing to screen in time steps */
int TimeToSweepFullBeach =
    50; /* Spacing of full beach sweep for each rock cell i=0:TotalBeachCells
           LMV */
int LookDist =
    10; /* For short rock to beach sweep, +- number of beach cells to look at
           LMV */
float BigDistanceToBeach =
    1000.0 * CellWidth;   /* Used to find minimum distance to beach in rock to
                             beach sweep LMV */
float NoWeathering = 5.0; /* weathering of rock only occurs below this amt of
                             sed cover (vert equiv in meters) LMV */
int InteractivePlot = 0;
int StartStop = 0;    /* Stop after every iteration 'Q' to move on */
int InterruptRun = 0; /* Allow run to be paused by pressing the 'A' key */
int NoPauseRun = 1;   /* Disbale PauseRun subroutine */
int InitialPert =
    0; /* Start with a bump? if =1-->square pert, if =2-->pointy pert */
int DiffusiveHump = 0;     /* Shoreline is sinusoidal LMV	*/
int InitialSmooth = 0;     /* Smooth starting conditions */
int InitialSmoothRock = 1; /* Smooth rock interface starting conditions LMV */
int ChunkLength =
    100; /* Alongshore length of a chunk of fast or slow weathering rock LMV */
         /* used to be Ymax/NumberChunk */

/* De-bugging Parameters */

int debug0 = 0;       /* misc. printf's--don't leave them in */
int debug1 = 0;       /* Find Next Beach Cell */
int debug1a = 0;      /* Find Next Rock Cell LMV*/
int debug1b = 0;      /* More find next cell LMV*/
int debug1c = 0;      /* More find next rock cell LMV*/
int debug2 = 0;       /* Shadow Routine */
int debug2a = 0;      /* In-Depth Shadow Results */
int debug2b = 0;      /* More shadow results */
int debug3 = 0;       /* Determine Angles */
int debug4 = 0;       /* Upwind/Downwind */
int debug5 = 0;       /* Sediment Transport Decisions*/
int debug6 = 0;       /* Sediment Trans Calculations */
int debug6a = 0;      /* Volume In/Out/AcrossBorder LMV */
int debug7 = 0;       /* Transport Sweep (move sediment) */
int debug7a = 0;      /* Slope Calcs */
int debug8 = 0;       /* Full/Empty */
int debug80 = 0;      /* Temp Full */
int debug8a = 0;      /* More Full/Empty LMV*/
int debug9 = 0;       /* FixBeach */
int debug9a = 0;      /* More FixBeach LMV */
int debug10 = 0;      /* Distance to Beach LMV*/
int debug11 = 0;      /* Min Distance to Beach & Closest Beach LMV */
int debug12 = 0;      /* Weathering Rate LMV */
int debug12a = 0;     /* Percent Fines lost LMV */
int debug13 = 0;      /* Sed Flow LMV*/
int debug14 = 0;      /* Actual Volume Across Border & PFS LMV */
int debug15 = 0;      /* Erosion Rate per Year LMV*/
int debug16 = 0;      /* Total Percent Full LMV */
int debug17 = 0;      /* graphics stuff */
int debug18 = 0;      /* reading wave data AB*/
int debugNearest = 0; /* graphs a pixel at nearest beach cell to a given rock
                         cell (for weathering) */

/* Universal Constants */
#define pi 3.1415927

float g = 9.80665;
float radtodeg = 180 / pi; /* transform rads to degrees */

/* Overall Shoreface Configuration Arrays - Data file information */
char AllBeach[Xmax][2 * Ymax]; /* Flag indicating of cell is entirely beach */
char AllRock[Xmax][2 * Ymax]; /* Flag indicating if cell is entirely rock LMV */
double PercentFullSand[Xmax][2 * Ymax]; /* Fractional amount of cell full of
                                           sediment LMV */
double PercentFullRock[Xmax][2 * Ymax]; /* Fractional amount of a cell full of
                                           rock LMV */
char TypeOfRock[Xmax][2 * Ymax]; /* Array to control weathering rates of rock
                                    along the beach LMV */
int Age[Xmax][2 * Ymax];         /* Age since cell was deposited */

int NumSinks = 6;
int ColumnSinks =
    0; /* if this is turned on, X vector is disregarded and a sink is considered
          to span a whole column */
int SinkY[] = {158, 361, 158, 361,
               158, 361}; /* a sink is a cell that is routinely emptied (if it
                             is on the beach) */
int SinkX[] = {190, 190, 189, 189, 191, 191};

FILE *SaveSandFile;
FILE *WaveFileOutput;
FILE *InitMetaFile;
FILE *ReadSandFile;
FILE *ReadWaveFile;
FILE *ReadRealWaveFile;
FILE *ReadControlFile;

/* Computational Arrays (determined for each time step) */

int X[MaxBeachLength];                  /* X Position of ith beach element */
int Y[MaxBeachLength];                  /* Y Position of ith beach element */
int XBehind[MaxBeachLength];            /* Cell that is "behind" X[i] LMV */
int YBehind[MaxBeachLength];            /* Cell that is "behind" Y[i] LMV */
int XRock[MaxBeachLength];              /* X Position of ith rock element LMV */
int YRock[MaxBeachLength];              /* Y Position of ith rock element LMV */
int XRockBehind[MaxBeachLength];        /* Cell that is "behind" XRock[i] LMV */
int YRockBehind[MaxBeachLength];        /* Cell that is "behind" YRock[i] LMV */
char InShadow[MaxBeachLength];          /* Is ith beach element in shadow? */
float ShorelineAngle[MaxBeachLength];   /* Angle between cell and right (z+1)
                                           neighbor	*/
float SurroundingAngle[MaxBeachLength]; /* Angle between left and right
                                           neighbors */
char UpWind[MaxBeachLength]; /* Upwind or downwind condition used to calculate
                                sediment transport */
float VolumeIn[MaxBeachLength];  /* Sediment volume into ith beach element*/
float VolumeOut[MaxBeachLength]; /* Sediment volume out of ith beach element */
float VolumeAcrossBorder[MaxBeachLength]; /* Sediment volume across border of
                                             ith beach element in m^3/day LMV*/
/* amount sediment needed, not necessarily amount a cell gets */
float
    ActualVolumeAcross[MaxBeachLength]; /* Sediment volume that actually gets
                                           passed across border LMV */
/* amount sed is limited by volumes across border upstream and downstream */
char DirectionAcrossBorder[MaxBeachLength]; /* Flag to indicate if transport
                                               across border is L or R LMV*/
char FlowThroughCell[MaxBeachLength]; /* Is flow through ith cell Left, Right,
                                         Convergent, or Divergent LMV*/
float DistanceToBeach
    [MaxBeachLength]; /* Distance in meters from rock to beach LMV*/
float MinDistanceToBeach[MaxBeachLength]; /* From a rock cell j, min distance
                                             (in meters) to closest beach LMV*/
int ClosestBeach[MaxBeachLength]; /* i position of closest rock to beach LMV */
float AmountWeathered[MaxBeachLength]; /* Amount of rock weathered from rock
                                          cell j LMV */

/* Miscellaneous Global Variables */

int NextX; /* Global variables used to iterate FindNextCell in global array - */
int NextY; /*	would've used pointer but wouldn't work	*/
int BehindX;
int BehindY;
int BehindRockX;
int BehindRockY;
int NextRockX; /* Global variables used to iterate FindNextRockCell in global
                  array LMV*/
int NextRockY;
int TotalBeachCells; /* Number of cells describing beach at particular iteration
                        */
int TotalRockCells;  /* Number of cells describing rock at an iteration LMV */
int ShadowXMax;      /* used to determine maximum extent of beach cells */
float WaveAngle;     /* wave angle for current time step */
int FindStart;     /* Used to tell FindBeach at what Y value to start looking */
int FindRockStart; /* Used to tell FindRock at what Y value to start looking LMV
                      */
char FellOffArray; /* Flag used to determine if accidentally went off array */
char FellOffRockArray; /* Flag used to determine if accidentally went off rock
                          array LMV */
float MassInitial;     /* For conservation of mass calcs */
float MassCurrent;     /* " */
int device;
short button;
long buttonback;
int NumWaveBins;    /* For Input Wave - number of bins */
float WaveMax[36];  /* Max Angle for specific bin */
float WaveProb[36]; /* Probability of Certain Bin */
double WaveAngleIn;
double WaveHeightIn;
double WavePeriodIn;
double ControlFileIn[20]; /* Initialisation data read from file*/
/* Function Prototypes */

void AdjustShore(int i);
void AgeCells(void);
void ControlFile(void);
void DetermineAngles(void);
void DetermineSedTransport(void);
void DoSink(void);
void Transport(void);
void ErodeTheBeach(int i); /*LMV*/
void FindBeachCells(int YStart);
void FindRockCells(int YStart); /*LMV*/
char FindIfInShadow(int xin, int yin, int ShadMax);
void FindNearestBeach(int j); /*LMV*/
void FindNextCell(int x, int y, int z);
void FindNextRockCell(int x, int y, int z); /*LMV*/
float FindWaveAngle(void);
void FixBeach(void);
void FixFlow(void);    /*LMV*/
void FlowInCell(void); /*LMV*/
int getIndex(int x, int y, char interface);
void InitConds(void);
void InitPert(void);
float MassCount(void);
void OopsImEmpty(int x, int y);
void OopsImFull(int x, int y);
void PauseRun(int x, int y, int in);
void PeriodicBoundaryCopy(void);
float Raise(float b, float e);
float RandZeroToOne(void);
void ReadSandFromFile(void);
void ReadWaveIn(void);
void RealWaveIn(void);
void RockCalculations(void);
void SaveSandToFile(void);
void SaveLineToFile(void);
void SedTrans(int i, float ShoreAngle, char MaxT); /*changed LMV*/
void ShadowSweep(void);
void TransportSedimentSweep(void);
void WaveOutFile(void);
void WeatherRock(int j); /*LMV*/
int XMaxBeach(int Max);
void ZeroVars(void);

/*
CWT
Overwash function prototypes
*/

void CheckOverwashSweep(void);
float GetOverwashDepth(int xin, int yin, float xinfl, float yinfl, int ishore);
void CheckOverwash(int icheck);
void DoOverwash(int xfrom, int yfrom, int xto, int yto, float xintto,
                float yintto, float widthin, int ishore);

/*
CWT
These assignments and declarations moved here - needed for compilation, ecven
when not using graohics
*/

int AgeMax = 1000000; /* Maximum 'age' of cells - loops back to zero */
int AgeUpdate = 10;   /* Time space for updating age of non-beach cells */
int AgeShadeSpacing =
    10000; /* For graphics - how many time steps means back to original shade */
int CellPixelSize =
    4; /* Size in pixes of plotted cell (a power of two, please) */
void PutPixel(float x, float y, float R, float G, float B);

/*
CWT
Graphics and display stuff commented out to line 330

function prototypes

void 	ButtonEnter(void);
void    GraphCells(int bugx, int bugy);
void    OpenWindow(void);
void    PutPixel(float x, float y, float R, float G, float B);
void 	PrintLocalConds(int x, int y, int in);
void    ScreenInit(void);
Bool    WaitForNotify(Display*, XEvent*, char*);

Display Crap

static WINDOW *mainwnd;
static WINDOW *screen;
WINDOW *my_win;
float xcellwidth;
float ycellwidth;
int current_getch;
int xplotoff;
int yplotoff;
Display *dpy;
Window root;
XVisualInfo *vi;
Colormap cmap;
XSetWindowAttributes swa;
Window win;
GLXContext cx;
XEvent event;

Plotting Controls

int		EveryPlotSpacing = 30;
int		XPlotExtent = Xmax;	// Height in cells of plot window
in x - direction
int		YPlotExtent = Ymax;	// Width in cells of plot window
in y - direction (Ymax or 2*Ymax)
int		KeysOn = 0;
int		YPlotStart = 0; // at what y-value do you start plotting? (0
or Ymax/2)

*/

/* Phew!! Lets's Go! - Main Program */

int main(void) {
  /* Initialize Variables and Device */

  int xx;            /* duration loop variable */
  ShadowXMax = Xmax; /* RCL: was Xmax-5; */

  // srandom(seed); //AB
  if (seed == -999)
    srand(time(NULL));
  else
    srand(seed);

  /* Start from file or not? */
  if (PromptStart == 'y') {
    printf("shall we start from a file (y or n)? \n");
    scanf("%c", &StartFromFile);

    if (StartFromFile == 'y') {
      printf("Starting Filename? \n");
      scanf("%24s", readfilename);
      printf("Saving Filename? \n");
      scanf("%s", savefilename);
      printf("What time step are we starting at?");
      scanf("%d", &CurrentTimeStep);
      ReadSandFromFile();
    }

    if (StartFromFile == 'n') {
      printf("Saving Filename? \n");
      scanf("%s", savefilename);
      InitConds();
      if (InitialPert) {
        InitPert();
      }
      printf("InitConds OK \n");
    }
  }

  else if (StartFromFile == 'y') {
    ReadSandFromFile();
  }

  else {
    InitConds();
    if (InitialPert) {
      InitPert();
    }
  }

  PeriodicBoundaryCopy();

  /* Count Initial Mass */
  MassInitial = MassCount();

  FixBeach();

  if (SaveLine) SaveLineToFile();
  if (SaveFile) SaveSandToFile();

  /*Read in wave data?*/
  if (WaveIn == 1) ReadWaveIn(); /*Initialise WaveMax and WaveProb*/

  if (InitialiseFile == 'y')
    ControlFile(); /*Read in initialisation data to control model*/

  /*
  CWT Graphics calls commented out at present
  If using graphics, open display and graph

if (DoGraphics)
{
  xcellwidth = 2.0 / (2.0 * (float) XPlotExtent) * (CellPixelSize)/2.0;
  ycellwidth = 2.0 / (2.0 * (float) YPlotExtent) * (CellPixelSize)/2.0;
  xplotoff = 0;
  yplotoff = Ymax/2;

  // Initialize Graphics Window
  OpenWindow();

  if(EveryPlotSpacing) {
      GraphCells(-1, -1);
      printf("Plotting initial condition\n");
      //getchar();
          }
}
*/

  /* ---------------------------------------    PRIMARY PROGRAM LOOP
   * ----------------- --------------- */

  while (CurrentTimeStep < StopAfter) {
    /*  Time Step iteration - compute same wave angle for Duration time steps */

    /*  Calculate Wave Angle */
    WaveAngle = FindWaveAngle();
    // output wave data
    if (Wavedata == 'y') WaveOutFile();

    /*  Loop for Duration at the current wave sign and wave angle */
    for (xx = 0; xx < Duration; xx++) {
      MassCurrent = MassCount();

      /* Text to Screen? */
      if (CurrentTimeStep % ScreenTextSpacing == 0) {
        printf("==== WaveAngle: %2.2f  MASS Percent: %1.4f  Time Step: %d\n",
               180 * (WaveAngle) / pi, MassCurrent / MassInitial,
               CurrentTimeStep);
      }

      PeriodicBoundaryCopy();  // Copy visible data to external model - creates
                               // boundaries

      ZeroVars();

      /* Initialize for Find Beach Cells  (make sure strange beach does not
       * cause trouble */
      FellOffArray = 'y';
      FindStart = 1;
      FellOffRockArray = 'y'; /*LMV*/
      FindRockStart = 1;      /*LMV*/

      /*  Look for beach - if you fall off of array, bump over a little and try
       * again */
      while (FellOffArray == 'y') {
        FindBeachCells(FindStart);
        FindStart += FindCellError;
        if (FellOffArray == 'y') {
          /*printf("NOODLE  !!!!!FoundCells: %d GetO = %c \n",
           * FindStart,FellOffArray); */
          /*PauseRun(1,1,-1);*/
        }

        /* Get Out if no good beach spots exist - finish program*/
        if (FindStart > Ymax / 2 + 1) {
          printf("Stopped Finding Beach - done %d %d", FindStart, Ymax / 2 - 5);
          SaveSandToFile();
          getchar();
          return 1;
        }
      }

      /* LMV */
      while (FellOffRockArray == 'y') {
        FindRockCells(FindRockStart);
        /*printf("FoundRockCells: %d GetO = %c \n",
         * FindRockStart,FellOffRockArray);*/
        FindRockStart += FindCellError;
        if (FellOffRockArray == 'y') {
          /*printf("NOODLE  !!!!!FoundRockCells: %d GetO = %c \n",
           * FindRockStart,FellOffRockArray); */
          /*PauseRun(1,1,-1);*/
        }

        if (FindRockStart > Ymax / 2 + 1) {
          printf("Stopped Finding Rock - done %d %d", FindRockStart,
                 Ymax / 2 - 5);
          if (SaveFile) SaveSandToFile();
          return 1;
        }
      }

      if (debug0) {
        printf("going to pause after finding beach, rock\n");
        PauseRun(58, 269, -1);
      }

      RockCalculations(); /*LMV*/

      ShadowSweep();

      DetermineAngles();

      DetermineSedTransport();

      DoSink(); /* clear sand from sink cells */

      FlowInCell(); /*LMV*/

      FixFlow(); /*LMV*/

      TransportSedimentSweep();

      FixBeach();

      /* Age Empty Cells */
      if ((CurrentTimeStep % AgeUpdate == 0) && SaveAge) AgeCells();

      /* Count Mass */
      MassCurrent = MassCount();

      /* SAVE FILE ? */
      if ((CurrentTimeStep % SaveSpacing == 0 &&
           CurrentTimeStep > StartSavingAt) ||
          (CurrentTimeStep == StopAfter)) {
        if (SaveLine) SaveLineToFile();
        if (SaveFile) SaveSandToFile();
      }

      /*
      CWT Cgraphics commented out at present
      GRAPHING
                      if (DoGraphics && EveryPlotSpacing &&
      (CurrentTimeStep%EveryPlotSpacing == 0)) {
                              GraphCells(-1,-1);
                      }
      */
      CurrentTimeStep++;
    }
  }

  /* --------------------------------------- END OF MAIN
   * -------------------------------------------------- */

  printf("Run Complete.  Output file: %s \n", savefilename);
  getchar();

  return (0);
}

float FindWaveAngle(void)
/* calculates wave angle for given time step */
{
  float Angle;

  /* Data Input Method */

  float RandBin;   /* Random number to pick wave angle bin */
  float RandAngle; /* Random number to pick wave angle within the bin */
  int flag = 1;
  int i = 0;

  /* Method using input binned wave distribution -
   */
  /* variables WaveProb[], WaveMax[], previously input from file using
   * ReadWaveIn()	*/

  if (WaveIn == 1) {
    RandBin = RandZeroToOne();
    RandAngle = RandZeroToOne();

    /*printf("Time = %d RandBin = %f RandAng = %f\n",CurrentTimeStep, RandBin,
     * RandAngle);*/

    while (flag) {
      i++;

      if (RandBin < WaveProb[i]) {
        flag = 0;
      }
    }

    Angle = -(RandAngle * (WaveMax[i] - WaveMax[i - 1]) + WaveMax[i - 1]) * pi /
            180;

    /*printf("i = %d WaveMAx[i] = %f WaveMax[i-1] = %f WaveProb[i] = %f Angle=
       %f\n",
            i,WaveMax[i],WaveMax[i-1],WaveProb[i], Angle*180/pi);*/
  }

  if (WaveIn == 2) {
    RealWaveIn(); /*Read real data (angle/period/height)*/

    Angle = ((float)(WaveAngleIn)) +
            coastrotation; /*align waveclimate to coastline in CEM*/
    if (Angle >= 180.0)
      Angle -= 360.0;  // conversion to CEM angles (-180 -- 0 -- +180)
    Angle *= -1; /*-1 is there as wave angle is from 90 in the west to -90 in
                    the east*/

    Period = WavePeriodIn * waveperiodchange;
    OffShoreWvHt = WaveHeightIn * waveheightchange;

    if ((Angle > 90) || (Angle < -90)) {
      Angle = 0; /*if waves were to come from behind coast*/
      OffShoreWvHt =
          1e-10; /*if waves were to come from behind coast set to inf small*/
    }

    Angle /= radtodeg;
  }

  if (WaveIn == 0) {
    /* Method using wave probability step function
     */
    /*  	variable Asym will determine fractional distribution of waves coming
     * from the  		*/
    /*  	positive direction (positive direction coming from left)  -i.e.
     * fractional wave asymmetry */
    /*  redone by RCL. note Highness actually means lowness, for some reason */

    Angle = RandZeroToOne() * pi / 4.0;

    if (RandZeroToOne() >= Highness) Angle += pi / 4.0; /* make high angle */

    if (RandZeroToOne() >= Asym) Angle = -1.0 * Angle;

    /*printf("WaveAngle: %f degrees\n",Angle*radtodeg);*/
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

  xstart = Xmax - 1;
  y = YStart;

  while (AllBeach[xstart][y] == 'n') {
    xstart -= 1;
  }

  xstart += 1; /* Step back to where partially full beach */

  X[0] = xstart;
  Y[0] = YStart;

  if (debug1) printf("FirstX: %3d  FirstY: %3d  z: 0 \n", X[0], Y[0]);

  z = 0;
  X[z - 1] = X[z];
  Y[z - 1] = Y[z] - 1;

  while ((Y[z] < 2 * Ymax - 1) && (z < MaxBeachLength - 1)) {
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

    if (debug1) printf("NextX: %3d  NextY: %3d  z: %d \n", NextX, NextY, z);
    if (debug1)
      printf("z: %d, BehindX: %d; BehindY: %d; \n\n", z, BehindX, BehindY);

    if (PercentFullSand[X[z]][Y[z]] == 0) {
      /*printf("\nFINDBEACH: PercentFullSand Zero x: %d y: %d\n",X[z],Y[z]);*/
      /*PauseRun(X[z],Y[z],z);*/
    }

    /* If return to start point or go off left side of array, going the wrong
     * direction 	*/
    /* Jump off and start again closer to middle
     */

    if ((NextY < 1) || ((NextY == Y[0]) && (NextX == X[0])) ||
        (z > MaxBeachLength - 2)) {
      /*printf("!!!!!!!Fell Off!!!!!!!!!!!!!!!!! x = %d !!!!!!!!!!!!!",
       * NextX);*/
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
           CurrentTimeStep);
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
  int y, z, xstart; /* local iterators */

  /* Starting at left end, find the x - value for first cell that is 'allbeach'
   */

  xstart = Xmax - 1;
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

  while ((YRock[z] < 2 * Ymax - 1) && (z < MaxBeachLength - 1)) {
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
     * direction 	*/
    /* Jump off and start again closer to middle
     */

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
           CurrentTimeStep);
}

void FindNextCell(int x, int y, int z)

/* Function to find next cell that is beach moving in the general positive X
   direction */
/* changes global variables NextX and NextY, coordinates for the next beach cell
   */
/* This function will use but not affect the global arrays:  AllBeach [][], X[],
   and Y[] */
/* New thinking...5/04 LMV
   */
{
  if ((X[z - 1] == X[z]) && (Y[z - 1] == Y[z] - 1))
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
               0.0) /* RCL fixed typo*/
    {
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
   direction 		*/
/* changes global variables NextRockX and NextRockY, coordinates for the next
   beach cell       	*/
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
/* and amount of rock to weather
   */
/* This function calls FindNearestBeach and WeatherRock
   */
/* This function will use but not adjust the variable TotalRockCells
   */
/* Determines the total amount of rock weathered in one time step
   */

{
  int j;
  float TotalAmountWeathered = 0.0;

  for (j = 0; j <= TotalRockCells; j++) {
    FindNearestBeach(j);
    WeatherRock(j);
    TotalAmountWeathered += AmountWeathered[j];
  }
  /* printf("CurrentTimeStep = %d\n", CurrentTimeStep) */;
}

void FindNearestBeach(int j)
/* LMV */
/* Calculates the Euclidean distance between Rock and Beach for entire beach
   array	*/
/* Finds the minimum beach to rock distance
   */
/* This function will use but not adjust the global arrays: X[], Y[], XRock[],
   YRock[]	*/
/* PercentFullSand[][] 			*/
/* This function will use but not adjust the variables MinDistanceToBeach and
   TotalBeachCells */

{
  int i;
  float DeltaX, DeltaY;
  int LookStart, LookEnd; /* for short rock to beach sweep */

  MinDistanceToBeach[j] = BigDistanceToBeach;

  if (CurrentTimeStep % TimeToSweepFullBeach == 0)
  /* do a full sweep every now and then */
  {
    for (i = 0; i <= TotalBeachCells; i++) /*full sweep*/
    {
      if (X[i] != XRock[j]) /*vertical distance*/
      {
        DeltaX =
            ((X[i] - XRock[j] - 1) + (PercentFullSand[X[i]][Y[i]] +
                                      PercentFullSand[XRock[j]][YRock[j]])) *
            CellWidth;
        if (debug10)
          printf("Rock cell = %d, Beach cell = %d, DeltaX (1)= %f \n", XRock[j],
                 X[i], DeltaX);
      }

      else {
        DeltaX = (PercentFullSand[X[i]][Y[i]]) * CellWidth;
        if (debug10)
          printf("Rock cell = %d, Beach cell = %d, DeltaX (2)= %f \n", XRock[j],
                 X[i], DeltaX);
      }

      if (YRock[j] <= Y[i]) /*horizontal distance*/
      {
        /*if (XRock[j] = X[i])
        {
                DeltaY = 0.0;
                if (debug10) printf("DeltaY (0) = %f \n", DeltaY);
        }
        else
        { */
        DeltaY = (Y[i] - YRock[j]) * CellWidth;
        if (debug10) printf("DeltaY (1) = %f \n", DeltaY);
        /*}*/
      } else {
        DeltaY = (YRock[j] - Y[i]) * CellWidth;
        if (debug10) printf("DeltaY (2) = %f \n", DeltaY);
      }

      DistanceToBeach[j] = sqrt(
          DeltaX * DeltaX + DeltaY * DeltaY); /*RCL removed Raise w/neg base*/
      if (debug10)
        printf("Distance to beach from XRock[j] %d to X[i] %d is %f \n",
               XRock[j], X[i], DistanceToBeach[j]);

      if (DistanceToBeach[j] < MinDistanceToBeach[j]) {
        MinDistanceToBeach[j] = DistanceToBeach[j];
        ClosestBeach[j] = i;
      }
    }
    if (debug11)
      printf("Min Distance for j = %d is %f \n", j, MinDistanceToBeach[j]);
    if (debug11) printf("XRock[j] %d to X[i] %d \n", XRock[j], X[i]);
    if (debug11)
      printf("Closest Beach is at cell i = %d \n \n", ClosestBeach[j]);
    if (Y[i] == 150)
      printf("Min Distance for j = %d is %f \n", j, MinDistanceToBeach[j]);
    if (Y[i] == 150) printf("XRock[j] %d to X[i] %d \n", XRock[j], X[i]);
    if (Y[i] == 150)
      printf("Closest Beach is at cell i = %d \n \n", ClosestBeach[j]);

  }

  else
  /* LookDist allows the program to skip the full beach sweep at every iteration
     by looking at */
  /* the location of the previous closest beach and sweeping out a small arc
     plus or minus LookDist */
  {
    if (ClosestBeach[j] - LookDist >= 0)
      LookStart = ClosestBeach[j] - LookDist;
    else /* Only sweep right when on left end */
      LookStart = 0;

    if (ClosestBeach[j] + LookDist <= TotalBeachCells)
      LookEnd = ClosestBeach[j] + LookDist;
    else /* Only sweep left when at right end */
      LookEnd = TotalBeachCells;

    for (i = LookStart; i < LookEnd; i++) {
      if (X[i] != XRock[j]) {
        DeltaX =
            ((X[i] - XRock[j] - 1) + (PercentFullSand[X[i]][Y[i]] +
                                      PercentFullSand[XRock[j]][YRock[j]])) *
            CellWidth;
        if (debug10)
          printf("Rock cell = %d, Beach cell = %d, DeltaX (1)= %f \n", XRock[j],
                 X[i], DeltaX);
      }

      else {
        DeltaX = (PercentFullSand[X[i]][Y[i]]) * CellWidth;
        if (debug10)
          printf("Rock cell = %d, Beach cell = %d, DeltaX (2)= %f \n", XRock[j],
                 X[i], DeltaX);
      }

      if (YRock[j] <= Y[i]) /*horizontal distance*/
      {
        /*if (XRock[j] = X[i])
        {
                DeltaY = 0.0;
                if (debug10) printf("DeltaY (0) = %f \n", DeltaY);
        }
        else
        { */
        DeltaY = (Y[i] - YRock[j]) * CellWidth;
        if (debug10) printf("DeltaY (1) = %f \n", DeltaY);
        /*}*/
      } else {
        DeltaY = (YRock[j] - Y[i]) * CellWidth;
        if (debug10) printf("DeltaY (2) = %f \n", DeltaY);
      }

      DistanceToBeach[j] = sqrt(
          DeltaX * DeltaX + DeltaY * DeltaY); /*RCL removed Raise w/neg base*/
      if (debug10)
        printf("Distance to beach from XRock[j] %d to X[i] %d is %f \n",
               XRock[j], X[i], DistanceToBeach[j]);

      if (DistanceToBeach[j] < MinDistanceToBeach[j]) {
        MinDistanceToBeach[j] = DistanceToBeach[j];
        ClosestBeach[j] = i;
      }
    }
    if (debug11)
      printf("Short Sweep - Min Distance for j = %d is %f \n", j,
             MinDistanceToBeach[j]);
    if (debug11) printf("XRock[j] %d to X[i] %d \n", XRock[j], X[i]);
    if (debug11) printf("Closest Beach is at cell i = %d \n ", ClosestBeach[j]);
    if (debug11) printf("X[i]: %d, Y[i]: %d \n\n", X[i], Y[i]);
  }
  /*
  CWT Commented out; not used at present
  (debugNearest)
  PutPixel(Y[ClosestBeach[j]]*CellPixelSize,X[ClosestBeach[j]]*CellPixelSize,50,200,50);
  */
}

void WeatherRock(int j)
/* LMV */
/* Fuction changes PercentFullRock[][] and PercentFullSand[][] along the rock
   beach interface only	*/
/* At this point, weathering rate is pretty arbitrary
   */
/* Fuction also changes AllRock[][]
   */
/* Weathering rate changes along the shore, 'f' = fast weather, 's' = slow
   weathering			*/
/* Function uses (but does not change) TypeOfRock[][] to determine rate
   */

{
  float WeatheringRatePerYear;

  if (TypeOfRock[XRock[j]][YRock[j]] == 's') /*slow rock weathering*/
  {
    if ((MinDistanceToBeach[j] / 100) <=
        NoWeathering) /*think in terms of vertical equivalence */
      WeatheringRatePerYear =
          SlowWeatherCoeff *
          (exp((-MinDistanceToBeach[j]) /
               100)); /*exponential weathering based on sed cover*/
    else
      WeatheringRatePerYear = 0;
  }

  if (TypeOfRock[XRock[j]][YRock[j]] == 'f') /*fast rock weathering*/
  {
    if ((MinDistanceToBeach[j] / 100) <=
        NoWeathering) /*think in terms of vertical equivalence */
      WeatheringRatePerYear =
          FastWeatherCoeff *
          (exp((-MinDistanceToBeach[j]) /
               100)); /*exponential weathering based on sed cover*/
    else
      WeatheringRatePerYear = 0;
  }

  if (debug12)
    printf("\nPercentSand = %f PercentRock = %f\n",
           PercentFullSand[XRock[j]][YRock[j]],
           PercentFullRock[XRock[j]][YRock[j]]);
  if (debug12) printf("XRock[j]: %d, YRock[j]: %d \n", XRock[j], YRock[j]);

  AmountWeathered[j] = ((WeatheringRatePerYear * TimeStep) / 365) /
                       CellWidth; /*amount weathering per time step*/
  /*note to self: I messed this up for a while, now it's a percent weathered
   * (which correlates with erosion)*/

  if (PercentFullRock[XRock[j]][YRock[j]] >= AmountWeathered[j]) {
    PercentFullRock[XRock[j]][YRock[j]] -= AmountWeathered[j];
    PercentFullSand[XRock[j]][YRock[j]] += AmountWeathered[j];
    if (PercentFullRock[XRock[j]][YRock[j]] < 1.0)
      AllRock[XRock[j]][YRock[j]] = 'n';

    if (debug12)
      printf("(>=AW) PFR:%f, PFS:%f\n", PercentFullRock[XRock[j]][YRock[j]],
             PercentFullSand[XRock[j]][YRock[j]]);
  }

  if (PercentFullRock[XRock[j]][YRock[j]] < AmountWeathered[j]) {
    if ((PercentFullRock[XRock[j]][YRock[j]] +
         PercentFullRock[XRockBehind[j]][YRockBehind[j]]) >=
        AmountWeathered[j]) {
      AllRock[XRockBehind[j]][YRockBehind[j]] = 'n';
      PercentFullRock[XRockBehind[j]][YRockBehind[j]] -=
          (AmountWeathered[j] - PercentFullRock[XRock[j]][YRock[j]]);
      PercentFullSand[XRockBehind[j]][YRockBehind[j]] +=
          (AmountWeathered[j] - PercentFullRock[XRock[j]][YRock[j]]);
      PercentFullSand[XRock[j]][YRock[j]] +=
          PercentFullRock[XRock[j]][YRock[j]];
      PercentFullRock[XRock[j]][YRock[j]] = 0.0;

      if (debug12)
        printf("(<AW) PFR Behind:%f, PFS Behind:%f\n",
               PercentFullRock[XRockBehind[j]][YRockBehind[j]],
               PercentFullSand[XRockBehind[j]][YRockBehind[j]]);
    }

    else {
      printf("WeatheringProblem j: %d, XRock[j]: %d, YRock[j]: %d \n", j,
             XRock[j], YRock[j]);
      printf(
          "PFR[XRockBehind[j]][YRockBehind[j]]: %f, "
          "PFS[XRockBehind[j]][YRockBehind[j]]: %f \n",
          PercentFullRock[XRockBehind[j]][YRockBehind[j]],
          PercentFullSand[XRockBehind[j]][YRockBehind[j]]);
      /* PauseRun(-1,-1,j); */
      PercentFullSand[XRock[j]][YRock[j]] +=
          PercentFullRock[XRock[j]][YRock[j]];
      PercentFullRock[XRock[j]][YRock[j]] = 0.0;
    }

    if (debug12)
      printf("\n(<AW) j: %d, XRock[j]: %d, YRock[j]: %d \n", j, XRock[j],
             YRock[j]);
    /* PauseRun(XRock[j],YRock[j],j); */
  }

  if (PercentFullSand[XRock[j]][YRock[j]] < 0.0) {
    if (debug12) printf("PFS: %f \n", PercentFullSand[XRock[j]][YRock[j]]);
    PercentFullSand[XRock[j]][YRock[j]] = 0.0;
    /*printf("No Sand WeatheringProblem j: %d, XRock[j]: %d, YRock[j]: %d \n",
     * j, XRock[j], YRock[j]); */
  }

  if (debug12)
    printf("Weathering of cell j = %d, Amount Weathered = %f\n", j,
           AmountWeathered[j]);
  if (debug12)
    printf("MinDistanceToBeach[j]: %f\n", MinDistanceToBeach[j] / 100);

  /*After weathering, lose a percentage of weathered rock in the closest beach
   * cell*/
  /*This represents fine grained material that does not stay in nearshore system
   */

  /*not using this for now*/

  if ((AmountWeathered[j] > 0.0) && (TypeOfRock[XRock[j]][YRock[j]] == 'f')) {
    if (debug12a)
      printf(
          "j: %d, XRock[j]: %d, YRock[j]: %d, ClosestBeach[j]: %d, X[CB]: %d, "
          "Y[CB]: %d, Type: %c\n",
          j, XRock[j], YRock[j], ClosestBeach[j], X[ClosestBeach[j]],
          Y[ClosestBeach[j]], TypeOfRock[XRock[j]][YRock[j]]);
    if (debug12a)
      printf("PFS[CB]Before: %f\n",
             PercentFullSand[X[ClosestBeach[j]]][Y[ClosestBeach[j]]]);

    PercentFullSand[X[ClosestBeach[j]]][Y[ClosestBeach[j]]] -=
        PercentFineFast * AmountWeathered[j];

    if (debug12a) printf("	AmountWeathered[%d]: %f\n", j, AmountWeathered[j]);
    if (debug12a)
      printf("PFS[CB]After: %f\n\n\n",
             PercentFullSand[X[ClosestBeach[j]]][Y[ClosestBeach[j]]]);
  }

  if ((AmountWeathered[j] > 0.0) && (TypeOfRock[XRock[j]][YRock[j]] == 's')) {
    if (debug12a)
      printf(
          "j: %d, XRock[j]: %d, YRock[j]: %d, ClosestBeach[j]: %d, X[CB]: %d, "
          "Y[CB]: %d, Type: %c\n",
          j, XRock[j], YRock[j], ClosestBeach[j], X[ClosestBeach[j]],
          Y[ClosestBeach[j]], TypeOfRock[XRock[j]][YRock[j]]);
    if (debug12a)
      printf("PFS[CB]Before: %f\n",
             PercentFullSand[X[ClosestBeach[j]]][Y[ClosestBeach[j]]]);

    PercentFullSand[X[ClosestBeach[j]]][Y[ClosestBeach[j]]] -=
        PercentFineSlow * AmountWeathered[j];

    if (debug12a) printf("	AmountWeathered[%d]: %f\n", j, AmountWeathered[j]);
    if (debug12a)
      printf("PFS[CB]After: %f\n\n\n",
             PercentFullSand[X[ClosestBeach[j]]][Y[ClosestBeach[j]]]);
  }
}

void ShadowSweep(void)
/*  Moves along beach and tests to see if cells are in shadow 		*/
/*  This function will use and determine the Global array:  InShadow[]	*/
/*  This function will use and adjust the variable:   ShadowXMax	*/
/*  This function will use but not adjust the variable:  TotalBeachCells */
{
  int i;

  /* Find maximum extent of beach to use as a limit for shadow searching */

  /*ShadowXMax = XMaxBeach(ShadowXMax) + 3; */
  /* if(ShadowXMax > Xmax) ShadowXMax = Xmax; RCL */

  if (debug2)
    printf("ShadowXMax: %d   XMaxBeach: %d \n", ShadowXMax,
           XMaxBeach(ShadowXMax));

  /* Determine if beach cells are in shadow */
  for (i = 0; i <= TotalBeachCells; i++) {
    InShadow[i] = FindIfInShadow(X[i], Y[i], ShadowXMax);

    /* Code to see what is happening */

    if (InShadow[i] == 'y') {
      if (debug2)
        printf("Shadow at X: %3d  Y: %3d  i: %3d \n", X[i], Y[i], i);
      else if (debug2)
        printf("No Shadow X: %3d  Y: %3d  i: %3d \n", X[i], Y[i], i);
    }
  }
}

int XMaxBeach(int Max)

/* Finds extent of beach in x direction					*/
/* Starts searching at a point 3 rows higher than input Max		*/
/* Function returns integer value equal to max extent of 'allbeach'	*/

{
  int xtest, ytest;

  xtest = Max + 3;
  ytest = 0;

  while (xtest > 0) {
    while (ytest < 2 * Ymax) {
      if (AllBeach[xtest][ytest] == 'y') {
        return xtest;
      }

      ytest++;
    }

    ytest = 0;

    xtest--;
  }

  printf("***** Should've found Xmax for shadow): %d, %d ***** \n", xtest,
         ytest);

  return Xmax;
}

char FindIfInShadow(int xin, int yin, int ShadMax)
/*  Function to determine if particular cell xin,yin is in shadow */
/*  Returns a character 'y' or 'n' to indicate so */
/*  This function will use but not affect the global arrays: */
/*	AllBeach[][] and PercentFull[][] */
/*  This function refers to global variables:  WaveAngle, ShadowStepDistance
   */
{
  /* Initialize local variables */

  float xdistance, ydistance;

  int xcheck = xin;
  int ycheck = yin;
  int iteration = 1;

  while ((xcheck < ShadMax) && (ycheck > -1) && (ycheck < (2 * Ymax))) {
    /*  Find a cell along the projection line moving against wave direction */

    xdistance = iteration * ShadowStepDistance * cos(WaveAngle);
    xcheck = xin + rint(xdistance);

    ydistance = -iteration * ShadowStepDistance * sin(WaveAngle);
    ycheck = yin + rint(ydistance);

    if (debug2a) {
      /*
      CWT: call to PutPixel() commented out at present
      */
      //              PutPixel(ycheck*CellPixelSize,xcheck*CellPixelSize,0,250,250);
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

    if (xcheck >= Xmax - 1 || ycheck >= 2 * Ymax - 1) return 'n';

    /* If AllBeach is along the way, and not next neighbor */
    /* Probably won't get to this one, though			 	*/

    if ((AllBeach[xcheck][ycheck] == 'y') &&
        (/* (abs(ycheck-yin) !=1 && abs(xcheck-xin) !=1) &&*/
         ((xcheck + 1) >
          (xin + (PercentFullSand[xin][yin] + PercentFullRock[xin][yin]) +
           fabs((ycheck - yin) / tan(WaveAngle))))) &&
        (ycheck > -1) && (ycheck < 2 * Ymax - 1)) {
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
/* Compare a partially full cell's x - distance to a line projected  	*/
/* from the starting beach cell's x-distance				*/
/* This assumes that beach projection is in x-direction (not too bad) 	*/
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
             (ycheck > -1) && (ycheck < 2 * Ymax - 1)
             /* && !(xcheck == xin)*/) {
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

/*  Function to determine beach angles for all beach cells from left to right
   */
/*  By convention, the ShorelineAngle will apply to current cell and right
   neighbor	*/
/*  This function will determine global arrays:
   */
/*		ShorelineAngle[], UpWind[], SurroundingAngle[]
   */
/*  This function will use but not affect the following arrays and values:
   */
/*		X[], Y[], PercentFull[][], AllBeach[][], WaveAngle
   */
/*  ADA Revised underside, SurroundingAngle 6/03
   */
/*  ADA Rerevised, 3/04
   */

{
  int i, j, k; /* Local loop variables */

  /* Compute ShorelineAngle[]  */
  /* 	not equal to TotalBeachCells because angle between cell and rt neighbor
   */

  for (i = 0; i < TotalBeachCells; i++) {
    /* function revised 1-04 */
    if (Y[i] > Y[i + 1])
    /* On bottom side of Spit (assuming we must be if going right)  */
    {
      /* math revised 6-03 */

      ShorelineAngle[i] =
          pi - atan((X[i + 1] - (PercentFullSand[X[i + 1]][Y[i + 1]] +
                                 PercentFullRock[X[i + 1]][Y[i + 1]])) -
                    (X[i] - (PercentFullSand[X[i]][Y[i]]) +
                     PercentFullRock[X[i]][Y[i]]));

      if (ShorelineAngle[i] > pi) {
        ShorelineAngle[i] -= 2.0 * pi;
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
            ShorelineAngle[i], ShorelineAngle[i] * 180 / pi);

      /*PauseRun(X[i],Y[i],i);*/
    }

    else if ((Y[i] < Y[i + 1]))
    /* function revised 1-04 - asusme if goin right, on regular shore*/
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
            ShorelineAngle[i], ShorelineAngle[i] * 180 / pi);
    }

    else if (Y[i] == Y[i + 1] && X[i] > X[i + 1])
    /*  Shore up and down, on right side */
    {
      ShorelineAngle[i] =
          -pi / 2.0 + atan((Y[i + 1] + (PercentFullSand[X[i + 1]][Y[i + 1]] +
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
            ShorelineAngle[i], ShorelineAngle[i] * 180 / pi);
    }

    else if (Y[i] == Y[i + 1] && X[i] < X[i + 1])
    /* Shore up and down, on left side */
    {
      ShorelineAngle[i] =
          pi / 2.0 + atan((Y[i + 1] + (PercentFullSand[X[i + 1]][Y[i + 1]] +
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
            ShorelineAngle[i], ShorelineAngle[i] * 180 / pi);
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
       ShorelineAngle[i]*180/pi); */
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
          (ShorelineAngle[k - 1] + ShorelineAngle[k]) / 2 + pi;
      if (SurroundingAngle[k] > pi) {
        SurroundingAngle[k] -= 2.0 * pi;
      }
      if (debug4) printf("Under: %d\n", k);
    } else {
      SurroundingAngle[k] = (ShorelineAngle[k - 1] + ShorelineAngle[k]) / 2;
    }
  }

  /* Determine Upwind/downwind condition
   */
  /* Note - Surrounding angle is based upon left and right cell neighbors,
   */
  /* and is centered on cell, not on right boundary */

  if (debug4) printf("\nUp/Down   Wave Angle:%f\n", WaveAngle * radtodeg);

  for (j = 1; j < TotalBeachCells; j++) {
    if (debug4)
      printf("i: %d  Shad: %c Ang[i]: %3.1f  Sur: %3.1f  Effect: %3f  ", j,
             InShadow[j], ShorelineAngle[j] * radtodeg,
             SurroundingAngle[j] * radtodeg,
             (WaveAngle - SurroundingAngle[j]) * radtodeg);

    if (fabs(WaveAngle - SurroundingAngle[j]) >= 42.0 / radtodeg) {
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
   transport calcs	*/
/*  Once situation is determined, will use function SedTrans to determine actual
   transport	*/
/*  This function will call SedTrans which will determine global array:
   VolumeAcrossBorder[] LMV	*/
/*  This function will use but not affect the following arrays and values:
   */
/*		X[], Y[], InShadow[], UpWind[], ShorelineAngle[]
   */
/*  		PercentFullSand[][], AllBeach[][], WaveAngle
   */

{
  int i; /* Loop variable LMV also sent to SedTrans to determine VolAcrossBorder
            */
  float ShoreAngleUsed; /* Temporary holder for shoreline angle
                           */
  int CalcCell;         /* Cell sediment coming from to go across boundary i
                           */
  int Next, Last;       /* Indicators so test can go both left/right
                           */
  int Correction; /* Term needed for shoreline angle and i+1 case, angle stored
                     at i 	*/
  char UpWindLocal; /* Local holder for upwind/downwind condition
                       */
  char MaxTrans;    /* Do we need to compute using maximum transport ?
                       */
  int DoFlux;       /* Skip sed transport calcs (added 02/04 AA)
                       */

  if (debug5)
    printf("\nSEDTRANS: %d  @  %f \n\n", CurrentTimeStep, WaveAngle * radtodeg);

  for (i = 1; i < TotalBeachCells - 1; i++) {
    if (debug5) printf("\n  i: %d  ", i);

    MaxTrans = 'n';

    /*  Is littoral transport going left or right?	*/

    if ((WaveAngle - ShorelineAngle[i]) > 0) {
      /*  Transport going right, center on cell to left side of border	 */
      /*  Next cell in positive direction, no correction term needed */
      CalcCell = i;
      Next = 1;
      Last = -1;
      Correction = 0;
      DirectionAcrossBorder[i] = 'r'; /*LMV*/
      if (debug5) printf("RT  ");
    } else {
      /*  Transport going left, center on cell to right side of border 	*/
      /*  Next cell in negative direction, correction term needed */
      CalcCell = i + 1;
      Next = -1;
      Last = 1;
      Correction = -1;
      DirectionAcrossBorder[i] = 'l'; /*LMV*/
      if (debug5) printf("LT  ");
    }

    if (InShadow[CalcCell] == 'n') {
      /*  Adjustment for maximum transport when passing through 45 degrees
       */
      /*  This adjustment is only made for moving from downwind to upwind
       * conditions	*/
      /* 										*/
      /*  purposefully done before shadow adjustment, only use maxtran when
       */
      /*	transition from dw to up not because of shadow
       */
      /* keeping transition from uw to dw - does not seem to be big deal (04/02
       * AA) */

      if (((UpWind[CalcCell] == 'd') && (UpWind[CalcCell + Next] == 'u') &&
           (InShadow[CalcCell + Next] == 'n')) ||
          ((UpWind[CalcCell + Last] == 'u') && (UpWind[CalcCell] == 'd') &&
           (InShadow[CalcCell + Last] == 'n'))) {
        MaxTrans = 'y';
        if (debug5) printf("MAXTRAN  ");
      }

      /*  Upwind/Downwind adjustment Make sure sediment is put into shadows
       */
      /*  If Next cell is in shadow, use UpWind condition
       */

      DoFlux = 1;
      UpWindLocal = UpWind[CalcCell];

      if (InShadow[CalcCell + Next] == 'y') {
        UpWindLocal = 'u';
        if (debug5) printf("U(2)  ");
      }

      /*  If coming out of shadow, downwind should be used		*/
      /*  HOWEVER- 02/04 AA - if high angle, will result in same flux in/out
       * problem */
      /*  	solution  - no flux for high angle waves */

      if ((InShadow[CalcCell + Last] == 'y') && (UpWindLocal == 'u')) {
        DoFlux = 0;
        if (debug5) printf("U(X) NOFLUX \n");
      }

      /*  Use upwind or downwind shoreline angle for calcs */

      if (UpWindLocal == 'u') {
        ShoreAngleUsed = ShorelineAngle[CalcCell + Last + Correction];
        if (debug5)
          printf("UP  ShoreAngle: %3.1f  ", ShoreAngleUsed * radtodeg);
      } else if (UpWindLocal == 'd') {
        ShoreAngleUsed = ShorelineAngle[CalcCell + Correction];
        if (debug5)
          printf("DN  ShoreAngle: %3.1f  ", ShoreAngleUsed * radtodeg);
      }

      /* !!! Do not do transport on unerneath c'cause it gets all messed up */
      if (fabs(ShoreAngleUsed) > pi / 2.0) {
        DoFlux = 0;
      }

      /* Send to SedTrans to calculate VolumeIn and VolumeOut*/

      /* printf("i = %d  Cell: %d NextCell: %d Angle: %f Trans Angle: %f\n",
      i, CalcCell, CalcCell+Next, ShoreAngleUsed*180/pi, (WaveAngle -
      ShoreAngleUsed)*180/pi); */

      if (debug5)
        printf("From: %d  To: %d  TransAngle %3.1f", CalcCell, CalcCell + Next,
               (WaveAngle - ShoreAngleUsed) * radtodeg);

      if (DoFlux) {
        SedTrans(i, ShoreAngleUsed, MaxTrans); /*LMV*/
        /* SedTrans(CalcCell, CalcCell+Next, ShoreAngleUsed, MaxTrans); Andrew's
         * way*/
      }
    }
  }
}

void SedTrans(int i, float ShoreAngle, char MaxT)

/*  This central function will calcualte the sediment transported from the cell
   i-1 to		*/
/*  the cell at i, using the input ShoreAngle	LMV
   */
/*  This function will caluclate and determine the global arrays:
   */
/*		VolumeAcrossBorder[]	LMV
   */
/*  This function will use the global values defining the wave field:
   */
/*	WaveAngle, Period, OffShoreWvHt
   */
/*  Revised 6/02 - New iterative calc for refraction and breaking, parameters
   revised		*/
{
  /* Coefficients - some of these are important*/

  float StartDepth = 3 * OffShoreWvHt; /* m, depth to begin refraction calcs
                                          (needs to be beyond breakers)	*/
  float RefractStep =
      .2; /* m, step size to iterate depth for refraction calcs */
  float KBreak =
      0.5; /* coefficient for wave breaking threshold */
  float rho =
      1020; /* kg/m3 - density of water and dissolved matter */

  /* Variables */

  float AngleDeep;          /* rad, Angle of waves to shore at inner shelf  */
  float Depth = StartDepth; /* m, water depth for current iteration
                               */
  float Angle;              /* rad, calculation angle			*/
  float CDeep;              /* m/s, phase velocity in deep water		*/
  float LDeep;              /* m, offhsore wavelength			*/
  float C;                  /* m/s, current step phase velocity 		*/
  float kh;                 /* wavenumber times depth			*/
  float n;                  /* n						*/
  float WaveLength;         /* m, current wavelength			*/
  float WvHeight;           /* m, current wave height			*/

  /* Primary assumption is that waves refract over shore-parallel contours
   */
  /* New algorithm 6/02 iteratively takes wiave onshore until they break, then
   * computes Qs	*/
  /* See notes 06/05/02
   */

  if (debug6)
    printf("Wave Angle %2.2f Shore Angle  %2.2f    ", WaveAngle * radtodeg,
           ShoreAngle * radtodeg);

  AngleDeep = WaveAngle - ShoreAngle;

  if (MaxT == 'y') {
    AngleDeep = pi / 4.0;
  }
  if (debug6) printf("Deep Tranport Angle %2.2f \n\n", AngleDeep * radtodeg);

  /*  Don't do calculations if over 90 degrees, should be in shadow  */

  if (AngleDeep > 0.995 * pi / 2.0 || AngleDeep < -0.995 * pi / 2.0) {
    return;
  }

  else {
    /* Calculate Deep Water Celerity & Length, Komar 5.11 c = gT / pi, L = CT
     */

    CDeep = g * Period / (2.0 * pi);
    LDeep = CDeep * Period;
    if (debug6) printf("CDeep = %2.2f LDeep = %2.2f \n", CDeep, LDeep);

    while (TRUE) {
      /* non-iterative eqn for L, from Fenton & McKee 		*/

      WaveLength =
          LDeep *
          Raise(tanh(Raise(Raise(2.0 * pi / Period, 2) * Depth / g, .75)),
                2.0 / 3.0);
      C = WaveLength / Period;
      if (debug6)
        printf("DEPTH: %2.2f Wavelength = %2.2f C = %2.2f ", Depth, WaveLength,
               C);

      /* Determine n = 1/2(1+2kh/sinh(kh)) Komar 5.21			*/
      /* First Calculate kh = 2 pi Depth/L  from k = 2 pi/L		*/

      kh = 2 * pi * Depth / WaveLength;
      n = 0.5 * (1 + 2.0 * kh / sinh(2.0 * kh));
      if (debug6) printf("kh: %2.3f  n: %2.3f ", kh, n);

      /* Calculate angle, assuming shore parallel contours and no conv/div of
       * rays 	*/
      /* from Komar 5.47
       */

      Angle = asin(C / CDeep * sin(AngleDeep));
      if (debug6) printf("Angle: %2.2f", Angle * radtodeg);

      /* Determine Wave height from refract calcs - Komar 5.49 */

      WvHeight = OffShoreWvHt *
                 Raise(CDeep * cos(AngleDeep) / (C * 2.0 * n * cos(Angle)), .5);
      if (debug6) printf(" WvHeight : %2.3f\n", WvHeight);

      if (WvHeight > Depth * KBreak || Depth <= RefractStep)
        break;
      else
        Depth -= RefractStep;
    }

    /* Now Determine Transport */
    /* eq. 9.6b (10.8) Komar, including assumption of sed density = 2650 kg/m3
     */
    /* additional accuracy here will not improve an already suspect eqn for sed
     * transport 	*/
    /* (especially with poorly constrained coefficients),
     */
    /* so no attempt made to make this a more perfect imperfection
     */

    VolumeAcrossBorder[i] =
        fabs(1.1 * rho * Raise(g, 3.0 / 2.0) * Raise(WvHeight, 2.5) *
             cos(Angle) * sin(Angle) * TimeStep); /*LMV - now global array*/

    /*LMV VolumeIn/Out is now calculated below in AdjustShore */

    if (debug6a) printf("\ni: %d \n", i);
    if (debug6a) printf("VolumeAcrossBorder: %f  ", VolumeAcrossBorder[i]);
  }
}

void DoSink(void)
/* empties of sand any cell that has been designated a sink */
{
  int s, i, foundOne = 0;
  if (!HaveSinks || RandZeroToOne() > Sinkiness) return;

  if (ColumnSinks) /* treat sink as any beach cell at a particular y-value */
    for (s = 0; s < NumSinks; s++) {
      for (i = 0; i < MaxBeachLength; i++)
        if (Y[i] == SinkY[s]) {
          PercentFullSand[X[i]][Y[i]] = 0.0;
          foundOne = 1;
        }
      if (!foundOne) printf("couldn't find sink %d at y == %d!\n", s, SinkY[s]);
    }
  else /* treat sink as a specific cell; empty it if it is on the beach */
    for (s = 0; s < NumSinks; s++) {
      for (i = 0; i < MaxBeachLength; i++)
        if (Y[i] == SinkY[s] && X[i] == SinkX[s]) {
          PercentFullSand[X[i]][Y[i]] = 0.0;
          foundOne = 1;
        }
      /*if(!foundOne) printf("couldn't find sink %d at
       * (%d,%d)!\n",s,SinkX[s],SinkY[s]);*/
    }
  if (!foundOne) printf("didn't find any sinks @ time %d\n", CurrentTimeStep);
}

void FlowInCell(void)
/* LMV */
/* This function will determine if sediment transport into a cell is from right,
   from left,	*/
/*						converging into the cell, or diverging out
   */
/* Purpose is to address problems with insufficent sand, fix losing sand problem
   */
/* Moves from left to right direction
   */
/* This function will affect and determine the global arrays:  FlowThroughCell[]
   */
/* This function uses TotalBeachCells
   */
/* Uses but does not affect the global arrays:  DirectionAcrossBorder[]
   */
{
  int i;

  /* float AverageFull; */
  float TotalPercentInBeaches;

  for (i = 1; i <= TotalBeachCells; i++) {
    if ((DirectionAcrossBorder[i - 1] == 'r') &&
        (DirectionAcrossBorder[i] == 'r'))
      FlowThroughCell[i] = 'R'; /*Right*/
    if ((DirectionAcrossBorder[i - 1] == 'r') &&
        (DirectionAcrossBorder[i] == 'l'))
      FlowThroughCell[i] = 'C'; /*Convergent*/
    if ((DirectionAcrossBorder[i - 1] == 'l') &&
        (DirectionAcrossBorder[i] == 'r'))
      FlowThroughCell[i] = 'D'; /*Divergent*/
    if ((DirectionAcrossBorder[i - 1] == 'l') &&
        (DirectionAcrossBorder[i] == 'l'))
      FlowThroughCell[i] = 'L'; /*Left*/

    TotalPercentInBeaches +=
        (PercentFullSand[X[i]][Y[i]] + PercentFullRock[X[i]][Y[i]]);
  }

  /*
  AverageFull = TotalPercentInBeaches/TotalBeachCells;

  for (i=1; i <= TotalBeachCells; i++)
  {
  if ((PercentFullSand[X[i]][Y[i]] + PercentFullRock[X[i]][Y[i]]) > AverageFull)
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
/* Sweeps left to right	to check cells where DirectionAcrossBorder='r' (D
   to R, D to C, R to R, R to C) */
/* Sweeps right to left to check cells where DirectionAcrossBorder='l' (D to L,
   D to C, L to L, L to C) */
/* This function uses TotalBeachCells
   */
/* Uses but does not affect the global arrays: FlowThroughCell[],
   VolumeAcrossBorder[], PFS[][]		*/
/* reminder--VolumeAcrossBorder[] is always on right side of cell
   */
/* revised 5/04 LMV
   */

{
  int i;
  float AmountSand; /* amount of sand availible for transport (includes amount
                       in cell behind) LMV */

  float Depth; /* Depth of current cell */
  /* float	DeltaArea;	 Holds change in area for cell (m^2)*/
  float Distance;   /* distance from shore to intercept of equilib. profile and
                       overall slope (m)*/
  float Xintercept; /* X position of intercept of equilib. profile and overall
                       slope*/
  float DepthEffective; /* depth at x intercept, Brad: where slope becomes <
                           equilib slope,*/
  /*  			       so sand stays piled against shoreface*/

  for (i = 1; i <= TotalBeachCells - 1; i++)
    ActualVolumeAcross[i] = VolumeAcrossBorder[i];

  for (i = 1; i <= TotalBeachCells - 1;
       i++) /*get the actual volumes across for left and right borders of D
               cells*/
  {
    /* calculate effective depth */
    Depth = InitialDepth + ((X[i] - InitBeach) * CellWidth * ShelfSlope);
    Distance = Depth / (ShorefaceSlope - ShelfSlope * cos(ShorelineAngle[i]));
    Xintercept = X[i] + Distance * cos(ShorelineAngle[i]) / CellWidth;
    DepthEffective =
        InitialDepth + ((Xintercept - InitBeach) * CellWidth * ShelfSlope);

    if (DepthEffective < DepthShoreface) DepthEffective = DepthShoreface;

    /*set the D's first - left and right border set*/

    if (FlowThroughCell[i] == 'D') { /*Flow from a Divergent cell (to a
                                        Convergent cell or a Right cell)*/

      if (AllBeach[XBehind[i]][YBehind[i]] == 'y') {
        if (PercentFullRock[X[i]][Y[i]] == 0.0) {
          AmountSand = (PercentFullSand[X[i]][Y[i]] +
                        PercentFullSand[XBehind[i]][YBehind[i]]) *
                       CellWidth * CellWidth * DepthEffective;
          if (debug13) printf("(D)X: %d, Y: %d\n", X[i], Y[i]);
        } else {
          AmountSand = (PercentFullSand[X[i]][Y[i]]) * CellWidth * CellWidth *
                       DepthEffective;
          if (debug13) printf("(RockD)X: %d, Y: %d\n", X[i], Y[i]);
        }
      }

      else /*I don't know how much I have, just figure it out and give me some
              sand*/
      {
        /* printf("Check Amount sand X: %d, Y: %d \n", X[i], Y[i]); */
        AmountSand = (PercentFullSand[X[i]][Y[i]]) * CellWidth * CellWidth *
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
      /*divide up proportionally*/
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
       i++) /*get the actual volumes across right borders for R cells*/
  {
    /* calculate effective depth */
    Depth = InitialDepth + ((X[i] - InitBeach) * CellWidth * ShelfSlope);
    Distance = Depth / (ShorefaceSlope - ShelfSlope * cos(ShorelineAngle[i]));
    Xintercept = X[i] + Distance * cos(ShorelineAngle[i]) / CellWidth;
    DepthEffective =
        InitialDepth + ((Xintercept - InitBeach) * CellWidth * ShelfSlope);

    if (DepthEffective < DepthShoreface) DepthEffective = DepthShoreface;

    if (FlowThroughCell[i] ==
        'R') { /*Flow from a Right cell (to a Right cell or a Convergent cell)*/

      if (AllBeach[XBehind[i]][YBehind[i]] == 'y') {
        if (PercentFullRock[X[i]][Y[i]] == 0.0) {
          AmountSand = (PercentFullSand[X[i]][Y[i]] +
                        PercentFullSand[XBehind[i]][YBehind[i]]) *
                       CellWidth * CellWidth * DepthEffective;
          if (debug13) printf("(R)X: %d, Y: %d\n", X[i], Y[i]);
        } else {
          AmountSand = (PercentFullSand[X[i]][Y[i]]) * CellWidth * CellWidth *
                       DepthEffective;
          if (debug13) printf("(RockR)X: %d, Y: %d\n", X[i], Y[i]);
        }
      }

      else /*I don't know how much I have, just figure it out and give me some
              sand*/
      {
        /* printf("Check 2 Amount sand X: %d, Y: %d \n", X[i], Y[i]); */
        AmountSand = (PercentFullSand[X[i]][Y[i]]) * CellWidth * CellWidth *
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
       i--) /*get the volumes across left side of cell for L cells*/
  {
    Depth = InitialDepth + ((X[i] - InitBeach) * CellWidth * ShelfSlope);
    Distance = Depth / (ShorefaceSlope - ShelfSlope * cos(ShorelineAngle[i]));
    Xintercept = X[i] + Distance * cos(ShorelineAngle[i]) / CellWidth;
    DepthEffective =
        InitialDepth + ((Xintercept - InitBeach) * CellWidth * ShelfSlope);

    if (DepthEffective < DepthShoreface) DepthEffective = DepthShoreface;

    if (FlowThroughCell[i] ==
        'L') { /*Flow from a Left cell (to a Convergent cell or a Left cell)*/

      if (AllBeach[XBehind[i]][YBehind[i]] == 'y') {
        if (PercentFullRock[X[i]][Y[i]] == 0.0) {
          AmountSand = (PercentFullSand[X[i]][Y[i]] +
                        PercentFullSand[XBehind[i]][YBehind[i]]) *
                       CellWidth * CellWidth * DepthEffective;
          if (debug13) printf("(L)X: %d, Y: %d\n", X[i], Y[i]);
        } else {
          AmountSand = (PercentFullSand[X[i]][Y[i]]) * CellWidth * CellWidth *
                       DepthEffective;
          if (debug13) printf("(RockL)X: %d, Y: %d\n", X[i], Y[i]);
        }
      }

      else /*I don't know how much I have, just figure it out and give me some
              sand*/
      {
        /* printf("Check 3 Amount sand X: %d, Y: %d \n", X[i], Y[i]); */
        AmountSand = (PercentFullSand[X[i]][Y[i]]) * CellWidth * CellWidth *
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
/*  This function doesn't change any values, but the functions it calls do
   */
/*  Uses but doesn't change:  X[], Y[], PercentFullSand[] */
/*  sweepsign added to ensure that direction of actuating changes does not */
/*  	produce unwanted artifacts (e.g. make sure symmetrical */
{
  int i, ii, sweepsign;

  if (RandZeroToOne() * 2 > 1)
    sweepsign = 1;
  else
    sweepsign = 0;

  /*if (debug7) printf("\n\n TransSedSweep  Ang %f  %d\n", WaveAngle * radtodeg,
   * CurrentTimeStep);*/

  for (i = 0; i <= TotalBeachCells - 1; i++) {
    if (sweepsign == 1)
      ii = i;
    else
      ii = TotalBeachCells - 1 - i;

    if (debug7)
      printf("i: %d  ss: %d  X: %d  Y: %d  In: %.1f  Out: %.1f Across: %.1f\n",
             ii, sweepsign, X[i], Y[i], VolumeIn[i], VolumeOut[i],
             VolumeAcrossBorder[i]);
    if (debug0 && CurrentTimeStep == 4 && (ii == 485 || ii == 486)) {
      printf("pausing in tss before adjust (PFS[58][269] is %f)\n",
             PercentFullSand[58][269]);
      PauseRun(-1, -1, ii);
    }
    AdjustShore(ii);
    if (debug0 && CurrentTimeStep == 4 && (ii == 485 || ii == 486)) {
      printf("pausing in tss after adjust (PFS[58][269] is %f)\n",
             PercentFullSand[58][269]);
      PauseRun(-1, -1, ii);
    }
    if ((PercentFullSand[X[ii]][Y[ii]] +
         PercentFullRock[X[ii]]
                        [Y[ii]]) < -0.000001) { /* RCL changes 0.0 to 1E-6
                                                   whenever testing for
                                                   oopsimempty */
      if (debug0 && CurrentTimeStep == 4 && (ii == 485 || ii == 486)) {
        printf("pausing in tss after adjust, before oops (time 4, ii 485)\n");
        PauseRun(-1, -1, ii);
      }
      OopsImEmpty(X[ii], Y[ii]);
      /*printf("Empty called in TSS\n");*/
    }

    else if ((PercentFullSand[X[ii]][Y[ii]] + PercentFullRock[X[ii]][Y[ii]]) >
             1.0) {
      OopsImFull(X[ii], Y[ii]);
      /*printf("Full called in TSS\n");*/
    }
    if (debug0 && CurrentTimeStep == 4 && (ii == 485 || ii == 486)) {
      printf(
          "pausing in tss after adjust, oops, before erode, oops (PFS[58][269] "
          "is %f)\n",
          PercentFullSand[58][269]);
      PauseRun(-1, -1, ii);
    }

    ErodeTheBeach(ii);
    if (PercentFullSand[X[ii]][Y[ii]] < -0.000001) /* RCL */
    {
      OopsImEmpty(X[ii], Y[ii]);
      printf("Empty called in TSS after eroding\n");
    }
    if (debug0 && CurrentTimeStep == 4 && (ii == 485 || ii == 486)) {
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
/*  This function will change the global data array PercentFullSand[][]
   */
/*  Uses but does not adjust arrays: */
/*		X[], Y[], ShorelineAngle[], ActualVolumeAcross	LMV
   */
/*  Uses global variables: ShelfSlope, CellWidth, ShorefaceSlope, InitialDepth
   */
{
  float Depth;      /* Depth of current cell */
  float DeltaArea;  /* Holds change in area for cell (m^2)*/
  float Distance;   /* distance from shore to intercept of equilib. profile and
                       overall slope (m)*/
  float Xintercept; /* X position of intercept of equilib. profile and overall
                       slope*/
  float DepthEffective; /* depth at x intercept, Brad: where slope becomes <
                           equilib slope,*/
  /*  			       so sand stays piled against shoreface*/

  /*  FIND EFFECTIVE DEPTH, Deff, HERE */
  /*  Deff = FROM SOLUTION TO: */
  /*	A*Distance^n = D(X) - OverallSlope*cos(ShorelineAngle[X][Y])*Distance
   */
  /*  Where Distance is from shore, perpendicular to ShorelineAngle;
   */
  /* Brad's stuff not being used right now (linear slope for now) :
   */
  /*	Xintercept = X + Distance*cos(ShorelineAngle), and n = 2/3.
   */
  /*	THEN Deff = D(Xintercept); */
  /*	FIND "A" FROM DEPTH OF 10METERS AT ABOUT 1000METERS. */
  /*	FOR NOW, USE n = 1: */

  Depth = InitialDepth + ((X[i] - InitBeach) * CellWidth * ShelfSlope);

  Distance = Depth / (ShorefaceSlope - ShelfSlope * cos(ShorelineAngle[i]));

  Xintercept = X[i] + Distance * cos(ShorelineAngle[i]) / CellWidth;

  DepthEffective =
      InitialDepth + ((Xintercept - InitBeach) * CellWidth * ShelfSlope);

  if (debug7a)
    printf(
        "\n Depth %f  Distance %f XIntercept %f  DepthEffective %f "
        "DepthShoreFace %d",
        Depth, Distance, Xintercept, DepthEffective, DepthShoreface);

  if (DepthEffective < DepthShoreface) DepthEffective = DepthShoreface;
  /*LMV*/
  if (FlowThroughCell[i] == 'L') {
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
  PercentFullSand[X[i]][Y[i]] += DeltaArea / (CellWidth * CellWidth);

  if (debug14)
    printf("%c  i: %d  In: %f  Out: %f\n", FlowThroughCell[i], i, VolumeIn[i],
           VolumeOut[i]);
  if (debug14)
    printf("DeltaArea[%d]: %f, PFS[%d]: %f, TimeStep: %d\n", i, DeltaArea, i,
           PercentFullSand[X[i]][Y[i]], CurrentTimeStep);
  if (debug14) printf("	X[i]: %d, Y[i]: %d\n", X[i], Y[i]);
}

void ErodeTheBeach(int i)

/*  This function decreases the amount of sand in each cell of the
   TotalBeachCell array */
/*  by a uniform amount (ErosionRatePerYear)
   */
/*  Changes PercentFullSand[][]
   */
/*  LMV */

{
  float PercentEroded;

  if (InShadow[i] == 'n')
    PercentEroded = ((ErosionRatePerYear * TimeStep) / 365) / CellWidth;
  else if (InShadow[i] == 'y')
    PercentEroded = 0.0;

  if (PercentEroded > PercentFullSand[X[i]][Y[i]]) {
    if (PercentFullRock[X[i]][Y[i]] == 0.0) {
      if (PercentEroded < (PercentFullSand[X[i]][Y[i]] +
                           PercentFullSand[XBehind[i]][YBehind[i]])) {
        if (debug15)
          printf("	Bi: %d, PFSBefore: %f\n", i, PercentFullSand[X[i]][Y[i]]);
        PercentFullSand[XBehind[i]][YBehind[i]] -=
            PercentEroded - PercentFullSand[X[i]][Y[i]];
        PercentFullSand[X[i]][Y[i]] = 0.0;
        if (debug15)
          printf("	BPercentEroded: %f, PFSAfter: %f\n", PercentEroded,
                 PercentFullSand[X[i]][Y[i]]);
      }

      else if (PercentEroded > (PercentFullSand[X[i]][Y[i]] +
                                PercentFullSand[XBehind[i]][YBehind[i]]))
      /* if PFR 0 && PercentEroded > PFS + PFSBehind...*/
      {
        /*if(!(PercentFullSand[X[i]][Y[i]] <= 0.000001
        && PercentFullRock[XBehind[i]][YBehind[i]] >= 0.999999)) {
                printf("Check Erosion, something strange
        (%d,%d)\n",X[i],Y[i]);*/
        /*debug15 = 1;
        PrintLocalConds(X[i],Y[i]);*/
        /*}*/

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
        debug15 = 0;*/
      }
    } else {
      PercentFullSand[X[i]][Y[i]] = 0.0;
      if (debug15)
        printf("BPercentEroded: %f, PFSAfter: %f\n", PercentEroded,
               PercentFullSand[X[i]][Y[i]]);
    }
  }

  else /*normal case, enough sand */
  {
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
   in	*/
/*  Function completly changed 5/21/02 sandrevt.c
   */
/*  		New Approach - steal from all neighboring AllBeach cells
   */
/*		Backup plan - steal from all neighboring percent full > 0
   */
/*  Function adjusts primary data arrays:
   */
/*		AllBeach[][] and PercentFullSand[][]
   */

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

  if ((AllBeach[x - 1][y] == 'y') && (PercentFullRock[x - 1][y] == 0.0)) /*LMV*/
    emptycells += 1;
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

    if ((AllBeach[x - 1][y] == 'y') &&
        (PercentFullRock[x - 1][y] == 0.0)) /*LMV*/
    {
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
      /*if (debug8) PauseRun(x,y,-1);*/

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
      /*if (debug8) PauseRun(x,y,-1);*/

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
        /*if (debug8) PauseRun(x,y,-1);*/

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
        /*if (debug8) PauseRun(x,y,-1);*/

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
    AllBeach[x][y] = 'n'; /*LMV*/
  PercentFullSand[x][y] = 0.0;

  if (debug8) printf("\n");

  if (debug0 && PercentFullSand[58][269] < 0.0)
    printf("^^^^^^(58,259) made underfull (%f) after calling oops(%d,%d)^^^^\n",
           PercentFullSand[58][269], x, y);
}

void OopsImFull(int x, int y)

/*  If a cell is overfull, push beach out in new direction
   */
/*  Completely revised 5/20/02 sandrevt.c to resolve 0% full problems, etc.
   */
/*  New approach: 	put sand wherever 0% full in adjacent cells
   */
/*			if not 0% full, then fill all non-allbeach
   */
/*  Function adjusts primary data arrays:
   */
/*		AllBeach[][] and PercentFullSand[][]
   */
/*  LMV add rock adjustments
   */

{
  int fillcells = 0;
  int fillcells2 = 0;

  if (debug8)
    printf("\n	OOPS I'M FULL: X: %d  Y: %d PFS: %f PFR: %f", x, y,
           PercentFullSand[x][y], PercentFullRock[x][y]);
  /*if (debug8) PrintLocalConds(x,y,-1);*/

  /* find out how many cells will be filled up	*/

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

    if ((PercentFullSand[x - 1][y] + PercentFullRock[x - 1][y]) == 0.0) /*LMV*/
    {
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
      /*if (debug8) PauseRun(x,y,-1);*/

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
      /*if (debug8) PauseRun(x,y,-1);*/

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
    AllBeach[x][y] = 'y'; /*LMV*/
  PercentFullSand[x][y] = 1.0 - PercentFullRock[x][y];

  if (debug80 && PercentFullSand[x][y] < 0.0)
    printf("	PFS went Negative in OopsI'mFull...Bummer x: %d y:%d %f \n", x,
           y, PercentFullSand[x][y]);

  if (debug80 && PercentFullSand[x][y] > 1.0)
    printf("	PFS went >1 in OopsI'mFull...Bummer x: %d y:%d %f \n", x, y,
           PercentFullSand[x][y]);
}

void FixBeach(void)

/* Hopefully addresses strange problems caused by filling/emptying of cells
   */
/* Looks at entire data set */
/* Find unattached pieces of sand and moves them back to the shore */
/* Takes care of 'floating bits' of sand */
/* Also takes care of over/under filled beach pieces */
/* Revised 5/21/02 to move sand to all adjacent neighbors sandrevt.c */
/* Changes global variable PercentFullSand[][] */
/* Uses and can change AllBeach[][] */
/* sandrevx.c - added sweepsign to reduce chances of asymmetrical artifacts
   */
/* LMV rock fix	- do not fix bare bedrock cells */

{
  int i, x, y, sweepsign, done, counter;
  int fillcells3 = 0;
  int Corner = 0; /*LMV*/
  int xstart;

  /*if (debug9) printf("\n\nFIXBEACH      %d     %f\n", CurrentTimeStep,
   * WaveAngle*radtodeg);*/

  if (RandZeroToOne() * 2 > 1)
    sweepsign = 1;
  else
    sweepsign = 0;

  for (x = 1; x < ShadowXMax - 1; x++) {
    for (i = 1; i < 2 * Ymax - 1;
         i++) /* RCL: changing loops to ignore border */
    {
      if (sweepsign == 1)
        y = i;
      else
        y = 2 * Ymax - i;

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
              /*if (debug9) PauseRun(x,y,-1);*/

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
              /*if (debug9) PauseRun(x,y,-1);*/

              if (debug9a && PercentFullSand[x][y + 1] < 0.0)
                printf(
                    "PFS went Negative in FixBeach...Bummer x: %d y+1:%d %f \n",
                    x, y + 1, PercentFullSand[x][y + 1]);

              if (debug9a && PercentFullSand[x][y + 1] > 1.0)
                printf("PFS went >1 in FixBeach...Bummer x: %d y+1:%d %f \n", x,
                       y + 1, PercentFullSand[x][y + 1]);
            }

          } else {
            /* if (debug9) */
            printf("Complete fixbeach breakdown x: %d  y: %d\n", x, y);
            /*if (debug9) PauseRun(x,y,-1);*/

            xstart = Xmax - 1;
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

            printf("PFS After move: %f \n", PercentFullSand[xstart][y]);

            xstart = Xmax - 1;
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

          /* If we have overfilled any of the cells in this loop, need to
           * OopsImFull() */

          if ((PercentFullSand[x - 1][y] + PercentFullRock[x - 1][y]) > 1.0) {
            /* printf("	Below Overfilled\n");
            printf("PF[x-1][y]: %f", (PercentFullSand[x-1][y] +
            PercentFullRock[x-1][y]));
            printf("x-1:%d, y:%d\n", (x-1), y); */
            OopsImFull(x - 1, y);
            if (debug9) printf("	Below Overfilled\n");
          }
          if ((PercentFullSand[x][y - 1] + PercentFullRock[x][y - 1]) > 1.0) {
            /* printf("	Left Side Overfilled\n");
            printf("PF[x][y-1]: %f", (PercentFullSand[x][y-1] +
            PercentFullRock[x][y-1]));
            printf("x:%d, y-1:%d\n", x, (y-1)); */
            OopsImFull(x, y - 1);
            if (debug9) printf("	Left Side Overfilled\n");
          }
          if ((PercentFullSand[x][y + 1] + PercentFullRock[x][y + 1]) > 1.0) {
            /* printf("	Right Side Overfilled\n");
            printf("PF[x][y+1]: %f", (PercentFullSand[x][y+1] +
            PercentFullRock[x][y+1]));
            printf("x:%d, y+1:%d\n", x, (y+1)); */
            OopsImFull(x, y + 1);
            if (debug9) printf("	Right Side Overfilled\n");
          }
          if ((PercentFullSand[x + 1][y] + PercentFullRock[x + 1][y]) > 1.0) {
            /* printf("	Top Side Overfilled\n");
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
    for (x = 1; x < Xmax - 1; x++) {
      for (y = 1; y < 2 * Ymax - 1; y++) {
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
            -0.000001) /* RCL adds this, changes 0 to 1E-6 everywhere*/
        {
          OopsImEmpty(x, y);
          if (debug9)
            printf("found a sand < 0 in final fb loop (%d, %d) \n", x, y);
        }
      }
    }
  }
  // if(!done) printf("final fb loop didn't fix all problems\n");
  if (debug9) printf("done checking grid for overfull, underfull \n");
}

float MassCount(void)

/* Counts the total volume occupied by beach cells 	*/
/* Uses same algorhythm as AdjustShore			*/
/* returns a float of the total sum 			*/
/* Uses AllBeach[][] and PercentFullSand[][]		*/
/* and InitialDepth, CellWidth, ShelfSlope		*/
{
  int x, y;
  float Depth; /* Depth of current cell */
  float Mass = 0.0;

  for (x = 0; x < Xmax; x++) {
    Depth = InitialDepth + ((x - InitBeach) * CellWidth * ShelfSlope);

    for (y = 0; y < 2 * Ymax; y++) {
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

float RandZeroToOne(void)

/* function will return a random number equally distributed between zero and one
   */
/* currently this function has no seed */

{
  // return random()/(Raise(2,31)-1);
  double AB_rand = rand() % 1000;
  AB_rand = (AB_rand + 1) / 1000;
  return ((float)(AB_rand));
}

void InitNormal(void)
/* Creates initial beach conditions */
/* Flat beach with zone of AllBeach = 'y' separated by AllBeach = 'n' */
/* Block of rock parallel to beach begins at InitRock LMV */
/* LMV Assume that AllBeach = 'Y' for AllRock cells (and All Full cells)
   */
/* Bounding layer set to random fraction of fullness */

{
  int x, y, n; /* y1 */
  int InitialBeach;
  int InitialRock;
  int Amp; /*Amplitude of cos curve*/
  int blocks =
      1; /* have hard stuff just as blocks near the beach, or columns? */
  printf("Condition Initial \n");

  Amp = 10;

  for (y = 0; y < 2 * Ymax; y++)
    for (x = 0; x < Xmax; x++) {
      InitialBeach = InitBeach;
      InitialRock = InitRock;

      if (DiffusiveHump)
      /* shoreline is a cosine curve, diffusive LMV */
      {
        InitialRock = (InitRock + (-Amp * cos((2 * pi * y) / Ymax)) -
                       (InitBeach - InitRock - 1));
        InitialBeach = (InitBeach + (-Amp * cos((2 * pi * y) / Ymax)));
      }

      if (x < InitialRock) /*LMV*/
      {
        PercentFullRock[x][y] = 1.0;
        PercentFullSand[x][y] = 0.0;
        AllRock[x][y] = 'y';
        AllBeach[x][y] = 'n';
      }

      else if (x == InitialRock) /*LMV*/
      {
        if (InitialSmoothRock) {
          PercentFullRock[x][y] = 0.20;
          PercentFullSand[x][y] = 0.80;
        } else {
          PercentFullRock[x][y] = RandZeroToOne();
          PercentFullSand[x][y] = 1 - PercentFullRock[x][y];
        }
        AllRock[x][y] = 'n';
        AllBeach[x][y] = 'y';
      }

      else if ((x > InitialRock) && (x < InitialBeach)) /*LMV*/
      {
        PercentFullRock[x][y] = 0.0;
        PercentFullSand[x][y] = 1.0;
        AllRock[x][y] = 'n';
        AllBeach[x][y] = 'y';
      } else if (x == InitialBeach) {
        if (InitialSmooth) {
          PercentFullSand[x][y] = 0.5;
        } else {
          PercentFullSand[x][y] = RandZeroToOne();
        }
        AllRock[x][y] = 'n';
        AllBeach[x][y] = 'n';
      }

      else if (x > InitialBeach) {
        PercentFullSand[x][y] = 0;
        AllRock[x][y] = 'n';
        AllBeach[x][y] = 'n';
      }

      Age[x][y] = 0;
    }

  for (x = 0; x <= Xmax; x++) /*LMV Assign fast and slow weathering portions*/
  {
    for (n = 0; n <= 2 * NumberChunk; n++) {
      for (y = n * ChunkLength; y < ((n + 2) * ChunkLength); y++) {
        if (n % 3 == 0 && (!blocks || x > InitialRock - 3)) /* if even */
          TypeOfRock[x][y] = 's';
        else /*if odd*/
          TypeOfRock[x][y] = 'f';
      }
    }
    /*for (y=0;y<=Ymax;y++)
    printf("x %d, y %d, Type %c\n", x, y, TypeOfRock[x][y]);*/
  }
}

int initBlock(void) {
  int x, y, n, blockHeight = 4, blockWidth = 60,
               NumBlocks = 0; /* makes NumBlocks evenly spaced blocks */
  printf("init block\n");

  for (x = 0; x < Xmax; x++) {
    for (y = Ymax / 2; y < 3 * Ymax / 2; y++) {
      if (x < InitRock) {
        PercentFullRock[x][y] = 1.0;
        PercentFullSand[x][y] = 0.0;
        AllRock[x][y] = 'y';
        TypeOfRock[x][y] = 'f';
        AllBeach[x][y] = 'y';
      } else if (x < InitBeach) {
        PercentFullRock[x][y] = 0.0;
        PercentFullSand[x][y] = 1.0;
        AllRock[x][y] = 'n';
        AllBeach[x][y] = 'y';
      } else if (x == InitBeach) {
        PercentFullRock[x][y] = 0.0;
        PercentFullSand[x][y] = RandZeroToOne();
        AllRock[x][y] = 'n';
        AllBeach[x][y] = 'n';
      } else {
        PercentFullRock[x][y] = 0.0;
        PercentFullSand[x][y] = 0.0;
        AllRock[x][y] = 'n';
        AllBeach[x][y] = 'n';
      }
      Age[x][y] = 0;
    }
  }
  for (n = 1; n <= NumBlocks; n++) { /* the regularly-spaced blocks */
    for (x = InitRock; x >= InitRock - blockHeight; x--) {
      for (y = Ymax / 2 + n * Ymax / (NumBlocks + 1);
           y < Ymax / 2 + n * Ymax / (NumBlocks + 1) + blockWidth &&
           y < 2 * Ymax;
           y++) { /* block starts @ upper left corner */
        PercentFullRock[x][y] = 1.0;
        PercentFullSand[x][y] = 0.0;
        AllRock[x][y] = 'y';
        TypeOfRock[x][y] = 's';
        AllBeach[x][y] = 'y';
      }
    }
  }
  /* now here's space for ad hoc blocks */
  /* ad hoc block 1 */
  for (x = InitRock; x >= InitRock - blockHeight; x--) {
    for (y = Ymax / 2; y < Ymax / 2 + 10; y++) {
      PercentFullRock[x][y] = 1.0;
      PercentFullSand[x][y] = 0.0;
      AllRock[x][y] = 'y';
      TypeOfRock[x][y] = 's';
      AllBeach[x][y] = 'y';
    }
  }
  /* ad hoc block 2 */
  for (x = InitRock; x >= InitRock - blockHeight; x--) {
    for (y = 3 * Ymax / 2 - 90; y < 3 * Ymax / 2; y++) {
      PercentFullRock[x][y] = 1.0;
      PercentFullSand[x][y] = 0.0;
      AllRock[x][y] = 'y';
      TypeOfRock[x][y] = 's';
      AllBeach[x][y] = 'y';
    }
  }

  printf("done block init\n");
  return 0;
}

int InitWiggly(void) {
  int x, y, curve;
  int NumCurves = 4;
  int RockLine[2 * Ymax];
  int curveFactor =
      3; /* decrease fundamental wavelength/increase frequency (by a factor)
            --make it 1 for wavelength = Ymax */
  float Amp, AmpMax = InitRock / 6.0;
      /* InitRock/4.0; */ /* maximum amplitude of the component sin waves */
  int Min = 0.10 * Xmax,
      Max = 0.95 *
            Xmax; /* determines acceptable max amplitude of resulting curve */
  printf("init wiggly\n");

  for (y = 0; y < 2 * Ymax; y++) RockLine[y] = InitRock;

  for (curve = 1; curve <= NumCurves; curve++) {
    Amp = RandZeroToOne() * AmpMax;
    printf("for curve %d, amplitude %f\n", curve, Amp);
    for (y = 0; y < 2 * Ymax; y++) {
      RockLine[y] =
          RockLine[y] + Amp * sin(y * (curve * curveFactor) * 2.0 * pi / Ymax);
      if (RockLine[y] < Min || RockLine[y] > Max) return -1;
    }
  }

  for (x = 0; x < Xmax; x++)
    for (y = 0; y < 2 * Ymax; y++) {
      if (x < RockLine[y]) {
        PercentFullRock[x][y] = 1.0;
        PercentFullSand[x][y] = 0.0;
        AllRock[x][y] = 'y';
        TypeOfRock[x][y] = 's';
        AllBeach[x][y] = 'y';
      } else if (x == RockLine[y]) {
        PercentFullRock[x][y] = 0.0;
        PercentFullSand[x][y] = 0.5;
        AllRock[x][y] = 'n';
        AllBeach[x][y] = 'n';
      } else {
        PercentFullRock[x][y] = 0.0;
        PercentFullSand[x][y] = 0.0;
        AllRock[x][y] = 'n';
        AllBeach[x][y] = 'n';
      }
      Age[x][y] = 0;
    }

  printf("done wiggly init\n");
  return 0;
}

void InitConds(void) {
  CurrentTimeStep = 0;
  if (InitCType == 0)
    InitNormal();
  else if (InitCType == 1)
    while (InitWiggly() == -1)
      printf("init: amplitude too great, trying again...\n");
  else if (InitCType == 2)
    initBlock();
  else
    printf("have to pick an InitCType\n");
  return;
}

void InitPert(void)
/* Andrew's initial bump */
{
  int x, y;
  int PWidth = 5;
  int PHeight = 3;
  int PYstart = 25;

  if (InitialPert == 1)
  /* Square perturbation */
  {
    /* Fill AllBeach areas */

    for (x = InitBeach; x <= InitBeach + PHeight; x++) {
      for (y = PYstart; y <= PYstart + PWidth; y++) {
        PercentFullSand[x][y] = 1.0;
        AllBeach[x][y] = 'y';
      }
    }

    /* PercentFull Top */

    for (y = PYstart - 1; y <= PYstart + PWidth + 1; y++) {
      PercentFullSand[InitBeach + PHeight + 1][y] = RandZeroToOne();
    }

    /* PercentFull Sides */

    for (x = InitBeach; x <= InitBeach + PHeight; x++) {
      PercentFullSand[x][PYstart - 1] = RandZeroToOne();
      PercentFullSand[x][PYstart + PWidth + 1] = RandZeroToOne();
    }
  }

  else if (InitialPert == 2)
  /* Another Perturbation  - steep point */
  {
    x = InitBeach;

    PercentFullSand[x][17] = 0.8;
    PercentFullSand[x][18] = 1.0;
    AllBeach[x][18] = 'y';
    PercentFullSand[x][19] = 0.8;

    x = InitBeach + 1;

    PercentFullSand[x][17] = 0.6;
    PercentFullSand[x][18] = 1.0;
    AllBeach[x][18] = 'y';
    PercentFullSand[x][19] = 0.6;

    x = InitBeach + 2;

    PercentFullSand[x][17] = 0.2;
    PercentFullSand[x][18] = 1.0;
    AllBeach[x][18] = 'y';
    PercentFullSand[x][19] = 0.2;

    x = InitBeach + 3;

    PercentFullSand[x][18] = 0.3;
  }
}

void PeriodicBoundaryCopy(void)
/* Simulates periodic boundary conditions by copying middle section to front and
   end of arrays */
{
  int x, y;

  for (y = Ymax; y < 3 * Ymax / 2; y++)
    for (x = 0; x < Xmax; x++) {
      AllBeach[x][y - Ymax] = AllBeach[x][y];
      AllRock[x][y - Ymax] = AllRock[x][y];       /*LMV*/
      TypeOfRock[x][y - Ymax] = TypeOfRock[x][y]; /*LMV*/
      PercentFullSand[x][y - Ymax] = PercentFullSand[x][y];
      PercentFullRock[x][y - Ymax] = PercentFullRock[x][y]; /*LMV*/
      Age[x][y - Ymax] = Age[x][y];
    }
  for (y = Ymax / 2; y < Ymax; y++)
    for (x = 0; x < Xmax; x++) {
      AllBeach[x][y + Ymax] = AllBeach[x][y];
      AllRock[x][y + Ymax] = AllRock[x][y];       /*LMV*/
      TypeOfRock[x][y + Ymax] = TypeOfRock[x][y]; /*LMV*/
      PercentFullSand[x][y + Ymax] = PercentFullSand[x][y];
      PercentFullRock[x][y + Ymax] = PercentFullRock[x][y]; /*LMV*/
      Age[x][y + Ymax] = Age[x][y];
    }
}

void ZeroVars(void)
/* Resets all arrays recalculated at each time step to 'zero' conditions */
{
  int z;

  for (z = 0; z < MaxBeachLength; z++) {
    X[z] = -1;
    Y[z] = -1;
    XRock[z] = -1; /*LMV*/
    YRock[z] = -1; /*LMV*/
    XBehind[z] = -1;
    YBehind[z] = -1;
    XRockBehind[z] = -1;
    YRockBehind[z] = -1;
    InShadow[z] = '?';
    ShorelineAngle[z] = -999;
    UpWind[z] = '?';

    VolumeAcrossBorder[z] = 0.0; /*LMV*/
    ActualVolumeAcross[z] = 0.0; /*LMV*/
  }
}

void ReadSandFromFile(void)
/*  Reads saved output file, AllBeach[][] & PercentFullSand[][]	 */
/*  LMV - Also reads saved output file, AllRock[][], PercentFullRock[][],
   TypeOfRock[][]  */
{
  int x, y;

  ReadSandFile = fopen(readfilename, "r");
  printf("CHECK READ \n");

  if (SaveFile != 2) /*line file input*/
  {
    for (y = Ymax / 2; y < 3 * Ymax / 2; y++)
      for (x = 0; x < Xmax; x++) fscanf(ReadSandFile, " %c", &TypeOfRock[x][y]);

    for (y = Ymax / 2; y < 3 * Ymax / 2; y++) {
      for (x = 0; x < Xmax; x++) {
        fscanf(ReadSandFile, " %lf", &PercentFullRock[x][y]);

        if (PercentFullRock[x][y] >= 1.0)
          AllRock[x][y] = 'y';
        else
          AllRock[x][y] = 'n';
      }
    }

    for (y = Ymax / 2; y < 3 * Ymax / 2; y++) {
      for (x = 0; x < Xmax; x++) {
        fscanf(ReadSandFile, " %lf", &PercentFullSand[x][y]);

        if ((PercentFullSand[x][y] + PercentFullRock[x][y]) >= 1.0)
          AllBeach[x][y] = 'y';
        else
          AllBeach[x][y] = 'n';
      }
    }

    if (SaveAge)

      for (y = Ymax / 2; y < 3 * Ymax / 2; y++) {
        for (x = 0; x < Xmax; x++) {
          fscanf(ReadSandFile, " %d", &Age[x][y]);
        }
      }
  }

  else /*Array file input*/
  {
    for (x = (Xmax - 1); x >= 0; x--) {
      for (y = Ymax / 2; y < 3 * Ymax / 2; y++) {
        fscanf(ReadSandFile, " %c", &TypeOfRock[x][y]);
      }
      fscanf(ReadSandFile, "\n");
    }

    for (x = (Xmax - 1); x >= 0; x--) {
      for (y = Ymax / 2; y < 3 * Ymax / 2; y++) {
        fscanf(ReadSandFile, " %lf", &PercentFullRock[x][y]);

        if (PercentFullRock[x][y] >= 1.0)
          AllRock[x][y] = 'y';
        else
          AllRock[x][y] = 'n';
      }
      fscanf(ReadSandFile, "\n");
    }

    for (x = (Xmax - 1); x >= 0; x--) {
      for (y = Ymax / 2; y < 3 * Ymax / 2; y++) {
        fscanf(ReadSandFile, " %lf", &PercentFullSand[x][y]);

        if ((PercentFullSand[x][y] + PercentFullRock[x][y]) >= 1.0)
          AllBeach[x][y] = 'y';
        else
          AllBeach[x][y] = 'n';
      }
      fscanf(ReadSandFile, "\n");
    }

    if (SaveAge) {
      for (x = (Xmax - 1); x >= 0; x--) {
        for (y = Ymax / 2; y < 3 * Ymax / 2; y++) {
          fscanf(ReadSandFile, " %d", &Age[x][y]);
        }
        fscanf(ReadSandFile, "\n");
      }
    }
  }

  /*PrintLocalConds(5,5,-1);*/
  fclose(ReadSandFile);
  printf("file read!");

  PeriodicBoundaryCopy();
}

void SaveSandToFile(void)
/*  Saves current AllBeach[][] and PercentFullSand[][] data arrays to file
   */
/*  Save file name will add extension '.' and the CurrentTimeStep
   */
/*  LMV - Save AllRock[][], PercentFullRock[][], TypeOfRock[][]
   */
{
  int x, y;
  char savename[40];

  printf("\n saving \n ");

  sprintf(savename, "%s_%d.out", savefilename, CurrentTimeStep);
  printf("Saving as: %s \n \n \n", savename);

  SaveSandFile = fopen(savename, "w");

  if (!SaveSandFile) {
    printf("problem opening output file\n");
    exit(1);
  }
  if (SaveFile == 1) {
    for (y = Ymax / 2; y < 3 * Ymax / 2; y++)
      for (x = 0; x < Xmax; x++)
        fprintf(SaveSandFile, " %c", TypeOfRock[x][y]); /*LMV*/

    for (y = Ymax / 2; y < 3 * Ymax / 2; y++)
      for (x = 0; x < Xmax; x++)
        fprintf(SaveSandFile, " %lf", PercentFullRock[x][y]); /*LMV*/

    for (y = Ymax / 2; y < 3 * Ymax / 2; y++)
      for (x = 0; x < Xmax; x++)
        fprintf(SaveSandFile, " %lf", PercentFullSand[x][y]);

    if (SaveAge)
      for (y = Ymax / 2; y < 3 * Ymax / 2; y++)
        for (x = 0; x < Xmax; x++) fprintf(SaveSandFile, " %d", Age[x][y]);
  }

  // Array output
  if (SaveFile == 2) {
    // for (x=0; x<Xmax; x++)
    for (x = (Xmax - 1); x >= 0; x--) {
      for (y = Ymax / 2; y < 3 * Ymax / 2; y++) {
        fprintf(SaveSandFile, " %c", TypeOfRock[x][y]); /*LMV*/
        // if (TypeOfRock[x][y]=='f') fprintf(SaveSandFile, " 0");	 /*
        // Switch on if numbers required*/
        // else fprintf(SaveSandFile, " 1");
      }
      fprintf(SaveSandFile, "\n");
    }

    for (x = (Xmax - 1); x >= 0; x--) {
      for (y = Ymax / 2; y < 3 * Ymax / 2; y++) {
        fprintf(SaveSandFile, " %lf", PercentFullRock[x][y]); /*LMV*/
      }
      fprintf(SaveSandFile, "\n");
    }

    for (x = (Xmax - 1); x >= 0; x--) {
      for (y = Ymax / 2; y < 3 * Ymax / 2; y++) {
        fprintf(SaveSandFile, " %lf", PercentFullSand[x][y]);
      }
      fprintf(SaveSandFile, "\n");
    }

    if (SaveAge) {
      for (x = (Xmax - 1); x >= 0; x--) {
        for (y = Ymax / 2; y < 3 * Ymax / 2; y++) {
          fprintf(SaveSandFile, " %d", Age[x][y]);
        }
        fprintf(SaveSandFile, "\n");
      }
    }
  }

  fclose(SaveSandFile);
  printf("file saved!\n");

  PeriodicBoundaryCopy();
}

void SaveLineToFile(void)
/*  Saves data line of shoreline position rather than entire array */
/*  Main concern is to have only one data point at each alongshore location
   */
/*  Save file name will add extension '.' and the CurrentTimeStep
   */
{
  int y, x, xtop, i;
  float xsave;
  char savename[40];

  printf("\n saving \n ");

  sprintf(savename, "ShorePos_%d.dat", CurrentTimeStep / SaveSpacing);
  printf("Saving as: %s                 ", savename);

  SaveSandFile = fopen(savename, "w");
  if (!SaveSandFile) {
    printf("problem opening output file\n");
    exit(1);
  }

  for (y = Ymax / 2; y < 3 * Ymax / 2; y++) {
    x = Xmax - 1;
    xtop = Xmax;

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
    fprintf(SaveSandFile, " %f", xsave - InitBeach + 0.5);

    /*printf("y %d , xtop = %d xsave = %f \n",y,xtop,xsave);
            if (xtop != Xmax)
                    PauseRun(x+1,y,-1);			*/
  }

  fclose(SaveSandFile);
  printf("line file saved!\n\n");
}

void PauseRun(int x, int y, int in)
/* Pauses run until the 'q' key is pressed 	*/
/* Can Print or Plot Out Useful info		*/
{
  int xsee = 1, ysee = -1, isee = -1, i;

  printf("\nPaused \n");

  if (SaveLine) SaveLineToFile();
  if (SaveFile) SaveSandToFile();

  if (NoPauseRun) return;

  printf("\nend Pause\n");
}

void AgeCells(void)
/* Age Cells */
{
  int x, y;

  for (y = 0; y < 2 * Ymax; y++)
    for (x = 0; x < Xmax; x++)
      if (PercentFullSand[x][y] == 0) {
        Age[x][y] = CurrentTimeStep % AgeMax;
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
    fscanf(ReadWaveFile, " %f %f", &WaveMax[i], &WaveProb[i]);
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

  if (CurrentTimeStep == 0) {
    ReadRealWaveFile = fopen(readwavename, "r");
    printf("CHECK READ WAVE: %s\n", readwavename);
  }

  fscanf(ReadRealWaveFile, "%lf %lf %lf\n", &WaveAngleIn, &WavePeriodIn,
         &WaveHeightIn);

  while (WaveAngleIn == -999 || WavePeriodIn == -999 ||
         WaveHeightIn == -999) /*Recycling wave data*/
  {
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

  if (CurrentTimeStep == StopAfter)
    fclose(ReadRealWaveFile); /*Keep file open till end of run*/
}

void WaveOutFile(void)
/* Input Wave Distribution */
{
  if (CurrentTimeStep == 0) WaveFileOutput = fopen(wavesavename, "w");

  fprintf(WaveFileOutput, "%lf %lf %lf\n", (-WaveAngle * radtodeg), Period,
          OffShoreWvHt);

  if (CurrentTimeStep == StopAfter)
    fclose(WaveFileOutput); /*Keep file open till end of run*/
}

void ControlFile(void)
/* Input Wave Distribution */
{
  OffShoreWvHt = -999;
  Period = -999;
  Asym = -999;
  Highness = -999;
  TimeStep = -999;
  Duration = -999;
  StopAfter = -999;
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
  fscanf(ReadControlFile, "%lf %*[^\n]", &TimeStep);
  fscanf(ReadControlFile, "%lf %*[^\n]", &Duration);
  fscanf(ReadControlFile, "%lf %*[^\n]", &StopAfter);
  fscanf(ReadControlFile, "%lf %*[^\n]", &SlowWeatherCoeff);
  fscanf(ReadControlFile, "%lf %*[^\n]", &FastWeatherCoeff);
  fscanf(ReadControlFile, "%lf %*[^\n]", &PercentFineFast);
  fscanf(ReadControlFile, "%lf %*[^\n]", &PercentFineSlow);
  fscanf(ReadControlFile, "%lf %*[^\n]", &ErosionRatePerYear);
  fscanf(ReadControlFile, "%lf %*[^\n]", &coastrotation);
  fscanf(ReadControlFile, "%lf %*[^\n]", &waveheightchange);
  fscanf(ReadControlFile, "%lf", &waveperiodchange);
  fclose(ReadControlFile);

  if (TimeStep == -999) printf("***Initialisation file not read!***\n");

  if (Metadata == 'y') {
    printf("Outputting Initialisation Metadata \n");

    InitMetaFile = fopen(metasavename, "w");
    fprintf(InitMetaFile, "OffShoreWvHt = %lf\n", OffShoreWvHt);
    fprintf(InitMetaFile, "Period = %lf\n", Period);
    fprintf(InitMetaFile, "Asym = %lf\n", Asym);
    fprintf(InitMetaFile, "Highness = %lf\n", Highness);
    fprintf(InitMetaFile, "Timstep = %lf\n", TimeStep);
    fprintf(InitMetaFile, "Duration = %lf\n", Duration);
    fprintf(InitMetaFile, "StopAfter = %lf\n", StopAfter);
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

  /*Write Metadata*/
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
			printf("stuff hasn't been done yet) so no index\n"); */
      return -1;
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
  float vol = CellWidth * CellWidth * DepthShoreface;

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

  printf("AllRock\n"); /*LMV*/
  for (i = x + 2; i > x - 3; i--) {
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

  printf("PercentFullRock\n"); /*LMV*/
  for (i = x + 2; i > x - 3; i--) {
    for (j = y - 2; j < y + 3; j++) {
      printf("	%2.5f", PercentFullRock[i][j]);
    }
    printf("\n");
  }

  printf("\n");
  printf("CurrentTimeStep: %d\n", CurrentTimeStep); /*LMV*/
  printf("InitialMass: %f\n", MassInitial);         /*LMV*/
  printf("CurrentMass: %f\n", MassCurrent);         /*LMV*/
  printf(" %d  ", in);

  if (isee >= 0) {
    for (k = in - 3; k < in + 4; k++) {
      printf("  	%2d: %2d,%2d", k, X[k], Y[k]);
    }
    printf("\n\n\n");

    printf("Wave Angle:	%f\n\n", WaveAngle * radtodeg);
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
        ShorelineAngle[in - 2] * radtodeg, ShorelineAngle[in - 1] * radtodeg,
        ShorelineAngle[in] * radtodeg, ShorelineAngle[in + 1] * radtodeg,
        ShorelineAngle[in + 2] * radtodeg);
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

/*-----------------Graphics-------------------------
CWT Commented out at present

void ButtonEnter(void)
// Not sure what this does
{
    char newdigit = 'z';
    int flag = 0;
    int i = 1;
    char digits[7] = "";

    printf("Press <Space> to Start\n");
    printf("Enter Digit and <Space> (<Z> to finish)\n");

    i += 1;

    sprintf(digits, "%s%c", digits, newdigit);
    printf("\nCurrent %s\n", digits);

    newdigit = 'z';
}


void ScreenInit(void)
// this is for the keyboard thingies to work
{
    mainwnd = initscr();
    noecho();
    cbreak();
    nodelay(mainwnd, TRUE);
    halfdelay(.1);
    refresh(); // 1)
    wrefresh(mainwnd);
    screen = newwin(13, 27, 1, 1);
    box(screen, ACS_VLINE, ACS_HLINE);
}


Bool WaitForNotify(Display *d, XEvent *e, char *arg) {
    return (e->type == MapNotify) && (e->xmap.window == (Window)arg);
}


void OpenWindow(void) {
    static  int  attributeListSgl[]  =  {GLX_RGBA, GLX_RED_SIZE, 1,
GLX_GREEN_SIZE, 1, GLX_BLUE_SIZE, 1, None };
    static   int   attributeListDbl[]   =   { GLX_RGBA, 1, GLX_RED_SIZE, 1,
GLX_GREEN_SIZE, 1, GLX_BLUE_SIZE, 1, None };

        //     int FALSE = 0;
        //     int TRUE =1;
    int winwidth;
    int winheight;
    int swap_flag = FALSE;

    dpy = XOpenDisplay(0);

    vi = glXChooseVisual(dpy, DefaultScreen(dpy), attributeListSgl);

    if (vi == NULL) {
        vi = glXChooseVisual(dpy, DefaultScreen(dpy), attributeListDbl);
        swap_flag = TRUE;
    }
    cx = glXCreateContext(dpy, vi, 0, GL_TRUE);
    cmap = XCreateColormap(dpy, RootWindow(dpy, vi->screen),  vi->visual,
AllocNone);

    swa.colormap = cmap;
    swa.border_pixel = 0;
    swa.event_mask = StructureNotifyMask;

    winwidth = CellPixelSize*Xmax;
    winheight = CellPixelSize*Ymax;

    win = XCreateWindow(dpy, RootWindow(dpy, vi->screen), 0, 0, YPlotExtent *
CellPixelSize, XPlotExtent * CellPixelSize,
                                                0, vi->depth, InputOutput,
vi->visual,
                                                CWBorderPixel|CWColormap|CWEventMask,
&swa);
    XMapWindow(dpy, win);
    XIfEvent(dpy, &event, WaitForNotify, (char*)win);

    glXMakeCurrent(dpy, win, cx);
    glClearColor(0.0, 0.0, 0.0, 1);
    glClear(GL_COLOR_BUFFER_BIT);
    glFlush();

    printf("WIN %d \n ", dpy);
}

void GraphCells(int bugx, int bugy)
// Plots entire Array
{
    int x, y;
    float Red, Green, Blue, backRed, backGreen, backBlue;
    float AgeFactorRed, AgeFactorGreen, AgeFactorBlue, DepthBlue;
    float DepthFactorX;

        int xbug, ybug;
        xbug = bugx;
        ybug = bugy;

    DepthFactorX = XPlotExtent;

    for (y=yplotoff; y <= YPlotExtent+yplotoff; y++) {

        for (x= xplotoff; x <= XPlotExtent+xplotoff; x++) {

                        if (x == xbug && y == ybug) {
                                Red = 1;
                                Blue = 0;
                                Green = 0;
                        }
                        else {

                                backRed = 0;
                                backBlue = (75 + 130 * (x/DepthFactorX));
                                backGreen = (165 - 125 * (x/DepthFactorX));

                                DepthBlue =  floor( (.8 -
(float)x/(XPlotExtent*3)) * 255 )/255.0;
                                PutPixel( CellPixelSize*(y-yplotoff),
CellPixelSize*x,0,0, DepthBlue);

                                AgeFactorRed =
(float)((Age[x][y])%AgeShadeSpacing)/AgeShadeSpacing;
                                AgeFactorGreen =
(float)((Age[x][y]+AgeShadeSpacing/3)%AgeShadeSpacing)/AgeShadeSpacing;
                                AgeFactorBlue =
(float)((Age[x][y]+2*AgeShadeSpacing/3)%AgeShadeSpacing)/AgeShadeSpacing;

                                if ((PercentFullSand[x][y] > 0) &&
(AllBeach[x][y] == 'n')) {
                                        Red =((((235 - 100
*(AgeFactorRed))-backRed)* PercentFullSand[x][y] )+backRed)/255.0;
                                        Green =((((235 - 95 *
(AgeFactorGreen))-backGreen) * PercentFullSand[x][y])+backGreen)/255.0;
                                        Blue =((((210 - 150 *
AgeFactorBlue)-backBlue)* PercentFullSand[x][y])+backBlue)/255.0;

                                }
                                else if (AllBeach[x][y] == 'y') {
                                        Red =((((235 - 100
*(AgeFactorRed))-backRed) * PercentFullSand[x][y]) + backRed)/255.0;
                                        Green =((((235 - 95 *
(AgeFactorGreen))-backGreen) * PercentFullSand[x][y]) + backGreen)/255.0;
                                        Blue =(((210 - 150 *
(AgeFactorBlue)-backBlue) * PercentFullSand[x][y]) + backBlue)/255.0;
                                }

                                // LMV
                                else if ((PercentFullRock[x][y] > 0.0) &&
(AllRock[x][y] == 'n'))
                                {
                                        Red = ((((225 - 100 *
(AgeFactorRed))-backRed) * PercentFullSand[x][y]) + backRed)/255.0;
                                        Green = ((((225 - 95 *
(AgeFactorGreen))-backGreen) * PercentFullSand[x][y]) + backGreen)/255.0;
                                        Blue = ((((200 - 150 *
(AgeFactorBlue))-backBlue) * PercentFullSand[x][y]) + backBlue)/255.0;

                                }

                                else if ((AllRock[x][y] == 'y') &&
(TypeOfRock[x][y] == 's'))
                                {

                                        Red = (159 *
PercentFullRock[x][y])/255.0;
                                        Green = (89 *
PercentFullRock[x][y])/255.0;
                                        Blue = (39 *
PercentFullRock[x][y])/255.0;
                                }

                                else if ((AllRock[x][y] == 'y') &&
(TypeOfRock[x][y] == 'f'))
                                {
                                        Red = (205 *
PercentFullRock[x][y])/255.0;
                                        Green = (133 *
PercentFullRock[x][y])/255.0;
                                        Blue = (63 *
PercentFullRock[x][y])/255.0;
                                }

                                else {
                                        Red = backRed/255.0;
                                        Blue = backBlue/255.0;
                                        Green = backGreen/255.0;
                                }
                        }

                        PutPixel(x-xplotoff, y-yplotoff, Red, Green, Blue);

        }
    }

    x = 0;
    y = Ymax; // was StreamSpot, which was defined as Ymax in the "deltas"
parameters, I guess for Andrew's delta stuff...

    while (AllBeach[x][y] == 'y') {
        PutPixel(x-xplotoff, y-yplotoff, 1,0,0);
        x += 1;
    }

    glFlush();

}

void PutPixel(float x, float y, float R, float G, float B)
// translate x and y integer components to the openGL grid
{

    float xstart, ystart;

    xstart = (x - Xmax/2.0)/(Xmax/2.0);
    ystart = (y - Ymax/2.0)/(Ymax/2.0);

    glXMakeCurrent(dpy, win, cx);
    glColor3f(R, G, B);
    glBegin(GL_POLYGON);
    glVertex3f(ystart, xstart, 0.0);
    glVertex3f(ystart, xstart+xcellwidth, 0.0);
    glVertex3f(ystart + ycellwidth, xstart + xcellwidth, 0.0);
    glVertex3f(ystart + ycellwidth, xstart, 0.0);
    glEnd();

}

 CWT End of commenting out of graphics routines
 */

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

    /* To do test shoreline should be facing seaward
     */
    /* don't worry about shadow here, as overwash is not set to a time scale
     * with AST 	*/

    if ((fabs(SurroundingAngle[ii]) < (OverwashLimit / radtodeg)) &&
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
    ysign = 1;  // check right
  else
    ysign = -1;  // check left

  if (AllBeach[X[icheck] - 1][Y[icheck]] ==
          'y'  // if the cell below is all beach
      || ((AllBeach[X[icheck]][Y[icheck] - 1] ==
           'y')  // or if the cell to the left is all beach
          && (AllBeach[X[icheck]][Y[icheck] + 1] ==
              'y')))  // and the cell to the right is all beach
  /* 'regular condition' */
  /* plus 'stuck in the middle' situation (unlikely scenario)*/
  {
    xin = X[icheck] + PercentFullSand[X[icheck]][Y[icheck]];
    yin = Y[icheck] + 0.5;
  } else if (AllBeach[X[icheck]][Y[icheck] - 1] ==
             'y')  // if the cell to the left is all beach
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

  while ((checkdistance < CritBWidth) && (y > 0) && (y < 2 * Ymax) && (x > 1)) {
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
        Raise(((x - xin) * (x - xin) + (y - yin) * (y - yin)), .5) * CellWidth;

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
          measwidth = CellWidth * Raise((xint - xin) * (xint - xin) +
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
          measwidth = CellWidth * Raise((xint - xin) * (xint - xin) +
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
          measwidth = CellWidth * Raise((xint - xin) * (xint - xin) +
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
          measwidth = CellWidth * Raise((xint - xin) * (xint - xin) +
                                            (yint - yin) * (yint - yin),
                                        0.5);
        }
      } else if (PercentFullSand[xtest][ytest] > 0)
      /* uh oh - not good situation, no allbeach on sides */
      /* assume this is an empty cell, */
      {
        xint = x;
        yint = y;

        measwidth = CellWidth * Raise((xint - xin) * (xint - xin) +
                                          (yint - yin) * (yint - yin),
                                      0.5);

        printf(
            "-- Some Odd Over  xin: %2.2f  yin: %2.2f xt:%d yt: %d xint: %f "
            "yint: %f Meas: %3.2f Ang: %f Abs: %f\n",
            xin, yin, xtest, ytest, xint, yint, measwidth,
            SurroundingAngle[icheck] * radtodeg,
            fabs(SurroundingAngle[icheck]) * radtodeg);
      } else
      /* empty cell - oughta fill er up  - fill max barrier width*/
      {
        xint = x;
        yint = y;
        measwidth = CritBWidth - CellWidth;

        printf(
            "-- Empty Odd Over  xin: %2.2f  yin: %2.2f xt:%d yt: %d xint: %f "
            "yint: %f Meas: %3.2f Ang: %f Abs: %f\n",
            xin, yin, xtest, ytest, xint, yint, measwidth,
            SurroundingAngle[icheck] * radtodeg,
            fabs(SurroundingAngle[icheck]) * radtodeg);
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

  /*float DepthBackBarrier = 6.0;  	 m current set depth for backbarrier
   * (temp - make into function)*/

  float DepthBB; /* holds effective backbarrier depth */

  DepthBB = GetOverwashDepth(xto, yto, xintto, yintto, ishore);

  /* calculated value of most that backbarrier can move given geometry (true,
   * non-iterative solution) */

  if (DepthBB == DepthShoreface) {
    BBneed = MaxOver;
  } else {
    BBneed =
        (CritBWidth - widthin) / CellWidth / (1 - (DepthBB / DepthShoreface));
  }

  if (BBneed <= MaxOver)
  /* do all overwash */
  {
    delShore = BBneed * DepthBB / DepthShoreface;
    delBB = BBneed;
  } else
  /* only do overwash to max change) */
  {
    delShore = MaxOver * DepthBB / DepthShoreface;
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
    Depth = CellDepth[xdepth][yin];

    while ((Depth < 0) && (xdepth > 0)) {
      Depth = CellDepth[xdepth][yin];
      xdepth--;
    }

    if (Depth == DepthShoreface) {
      Depth = 6.0;
    }

    return Depth;

  } else if (OWType == 1)
  /* Geometric relation to determine depth through intersection of shorefaces
     */
  /* look in line determined by shoreline slope - reuse stepping function
     (again)	*/
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

    while ((!BackFlag) && (y > 0) && (y < 2 * Ymax) && (x > 1)) {
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
       = depthshoreface 	*/
    /* Periodic B.C.'s should make this not so important
       */
    {
      Depth = DepthShoreface;
    } else {
      BBDistance = Raise(((xinfl - xtest - PercentFullSand[xtest][ytest]) *
                          (xinfl - xtest - PercentFullSand[xtest][ytest])) +
                             ((yinfl - ytest - 0.5) * (yinfl - ytest - 0.5)),
                         .5);

      if (!FoundFlag)
      /* The backbarrier intersection isn't on the shoreline */
      /* Assume 1/2 of the length applies to this case */
      {
        Depth = BBDistance / 2 * ShorefaceSlope * CellWidth;
      } else
      /* Use the fancy geometry thing */
      {
        AngleUsed = 0;
        for (j = -1; j < 2; j++) {
          AngleUsed += SurroundingAngle[Backi + j];
        }
        AngleUsed = AngleUsed / 5;

        if (fabs(AngleUsed) > pi / 4.0) {
          AngleUsed = pi / 4.0;
        }

        AngleSin = sin(pi / 2.0 - fabs(SurroundingAngle[ishore] + AngleUsed));

        Depth = BBDistance * AngleSin / (1 + AngleSin);
      }
    }

    if (Depth < OWMinDepth) {
      Depth = OWMinDepth;
    } else if (Depth > DepthShoreface) {
      Depth = DepthShoreface;
    }
    return Depth;
  }

  printf("OWDepth all broken");
  PauseRun(xin, yin, -1);
  return DepthShoreface;
}
