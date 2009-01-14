/*  CAPERIFFIC                                                                */
/*  Program To generate capes?? and sandwaves?? using wave angle relationships*/
/*  Begun by Brad Murray 01/00                                                */
/*  Refined by Olivier Arnoult 01/00 - 06/00                                  */
/*  Revised and reformed by Andrew Ashton 06/00 -                             */
/*                                                                            */
/*  Program Notes -                                                           */
/*      To end program, press 'd' key and 'ESC' key simultaneously            */
/*      To save current iteration to file, press 's' and 'f' simultaneously   */
/*      To update screen display, press 'p' key                               */
/*                                                                            */

#include <stdlib.h>      /*THIS PROGRAM GONNA MAKE CAPES, SANDWAVES??*/
#include <stdio.h>
#include <math.h>
#include <time.h>
#ifdef WITH_OPENGL
#include <GL/glx.h> 
#include <GL/gl.h> 
#endif
#include <unistd.h>
#include <ncurses.h>

#undef DEBUG_ON

#ifdef DEBUG_ON
# if defined(__GNUC__)
#  define DEBUG_PRINT( exp, format... ) { if (exp) fprintf(stderr,format); }
# else
static void DEBUG_PRINT( int exp, const char* format, ... )
{
   if ( exp )
   {
      va_list args;
      va_start( args , format );
      fprintf( stderr, format, args );
      va_end( args );
   }
}
# endif
#else
# define DEBUG_PRINT( exp, format... ) { }
#endif

/*  Run Control Parameters */

#define TimeStep     0.2  /* days - reflects rate of sediment transport per
                             time step */
#define OffShoreWvHt 2    /* meters */
#define Period       7    /* seconds */
#define Asym         0.7  /* Fractional portion of waves coming from positive
                             (left) direction */
#define Highness     0.1  /* All New! .5 = even dist, > .5 high angle
                             domination */
#define Duration     1    /* Number of time steps calculations loop at same wave
                            angle */
#define StopAfter    2600 /* Stop after what number of time steps */

/* init,finalize, SaveSandToFile...*/
//char            savefilename[24] = "fileout";
/* init,ReadSandFromFile */
//char            readfilename[24] = "";
/* init, FindWaveAngle */
#define SAVE_FILENAME "fileout"
#define READ_FILENAME ""

typedef struct
{
   char savefilename[24];
   char readfilename[24];
}
Deltas_io;

Deltas_io _io = { SAVE_FILENAME, READ_FILENAME };

#define WAVE_IN (0) /* Input Wave Distribution file? */
//int     WaveIn = 0;             /* Input Wave Distribution file? */


/* Delta Info */


/* Unused!!! */#define NMouths         2       /* number of river mouths */
/* Unused!!! */float   QRiver[NMouths];

#define SedRate 0.02
#define StreamSpot Ymax

/* Aspect Parameters */

#define CellWidth       100.0  /* size of cells (meters) */
#define CritBWidth    350.0    /* width barrier maintains due to overwash (m) important scaling param! */
#define Xmax            200      /* number of cells in x (cross-shore) direction */
#define Ymax            500     /* number of cells in y (longshore) direction */
#define MaxBeachLength  8*Ymax  /* maximum length of arrays that contain beach data at each time step */
#define ShelfSlope      0.001   /* slope of continental shelf */
#define ShorefaceSlope  0.01    /* for now, linear slope of shoreface */
                                /* future : shoreface exponent m^1/3, from depth of ~10m at ~1000 meters */
#define DepthShoreface  10.0    /* minimum depth of shoreface due to wave action (meters) */   
#define InitBeach       30      /* cell where intial conditions changes from beach to ocean */ 
#define InitialDepth    9.0     /* theoretical depth in meters of continental shelf at x = InitBeach */
#define LandHeight      1.0     /* elevation of land above MHW  */
#define InitCType       0       /* type of initial conds 0 = sandy, 1 = barrier */
#define InitBWidth      4       /* initial minimum width of barrier (Cells) */
#define OWType        1         /* 0 = use depth array, 1 = use geometric rule */
#define OWMinDepth	0.1	/*  littlest overwash of all */
#define FindCellError	5	/* if we run off of array, how far over do we try again? */

/* Plotting Controls */
#define CELL_PIXEL_SIZE (4)
#define XPlotExtent     (Xmax)
#define YPlotExtent     (Ymax)
//int CELL_PIXEL_SIZE = 4;
//int XPlotExtent = Xmax;
//int YPlotExtent = Ymax;

/* DeBuggin Parameters */

#define DO_GRAPHICS       (1)
#define KeysOn            (0)
#define SaveAge           (1)    /* Save/update age of cells? */
#define PromptStart       ('n')  /* ask prompts for file names, etc? */
#define ScreenTextSpacing (1000) /* Spacing of writing to screen in time steps */ 
#define EveryPlotSpacing  (100)
#define StartStop         (0)    /* Stop after every iteration 'Q' to move on */
#define NoPauseRun        (1)    /* Disbale PauseRun subroutine */
#define InitialPert       (0)    /* Start with a bump? */
#define InitialSmooth     (0)    /* Smooth starting conditions */
#define WaveAngleSign     (1)    /* used to change sign of wave angles */

//int     DoGraphics = 1;
//int     KeysOn = 0;
//int 	SaveAge = 1;		/* Save/update age of cells? */
//char	PromptStart = 'n';	/* ask prompts for file names, etc? */
/* unused !!! */char 	OffArray = 'n';		/* Initializing this variable for later use */
//int	ScreenTextSpacing = 1000;/* Spacing of writing to screen in time steps */
//int	EveryPlotSpacing = 100;
//int	StartStop = 0;		/* Stop after every iteration 'Q' to move on */
/* unused */int	InterruptRun = 0;	/* Allow run to be paused by pressing the 'A' key */
//int	NoPauseRun = 1;		/* Disbale PauseRun subroutine */
//int	InitialPert = 0;	/* Start with a bump? */
//int	InitialSmooth = 0;	/* Smooth starting conditions */
//int	WaveAngleSign = 1;	/* used to change sign of wave angles */

#define DEBUG_0   (0)  /* Main program steps */
#define DEBUG_1   (0)  /* Find Next Cell */
#define DEBUG_2   (0)  /* Shadow Routine */
#define DEBUG_3   (0)  /* Determine Angles */
#define DEBUG_4   (0)  /* Upwind/Downwind */
#define DEBUG_5   (0)  /* Sediment Transport Decisions*/
#define DEBUG_6   (0)  /* Sediment Trans Calculations */
#define DEBUG_7   (0)  /* Transport Sweep (move sediment) */
#define DEBUG_7A  (0)  /* Slope Calcs */
#define DEBUG_8   (0)  /* Full/Empty */
#define DEBUG_9   (0)  /* FixBeach */
#define DEBUG_10A (0)  /* Overwash Tests*/
#define DEBUG_10B (0)  /* doing overwash (w/screen) */
int	OWflag = 0;     /* debugger */

/* Universal Constants */
//#define	pi		3.1415927
//#define	exp		2.7182818 /* e */
#define GRAV            (9.80665)
#define radtodeg        (180.0/M_PI) /* transform rads to degrees */
//float g =		9.80665;
//float radtodeg = 	180.0/pi; /* transform rads to degrees */

typedef struct
{
   char AllBeach[Xmax][2*Ymax]; /* Flag indicating of cell is entirely beach */
   float PercentFull[Xmax][2*Ymax]; /* Fractional amount of shore cell full of
                                       sediment */
   int Age[Xmax][2*Ymax]; /* Age since cell was deposited */
   float CellDepth[Xmax][2*Ymax]; /* Depth array (m) (ADA 6/3) */
   int X[MaxBeachLength]; /* X Position of ith beach element */
   int Y[MaxBeachLength]; /* Y Position of ith beach element */
   char	InShadow[MaxBeachLength];	 /* Is ith beach element in shadow? */
   float ShorelineAngle[MaxBeachLength]; /* Angle between cell and right (z+1)
                                            neighbor */
   float SurroundingAngle[MaxBeachLength];/* Cell-orientated angle based upon
                                             left and right neighbor */
   char UpWind[MaxBeachLength]; /* Upwind or downwind condition used to
                                   calculate sediment transport */
   float VolumeIn[MaxBeachLength];  /* Sediment volume into ith beach
                                       element */	
   float VolumeOut[MaxBeachLength]; /* Sediment volume out of ith beach
                                       element */

   int CurrentTimeStep; /* Time step of current calculation */ 

   int NextX; /* used to iterate FindNextCell in global array - */
   int NextY;

   int TotalBeachCells; /* Number of cells describing beach at particular iteration */
   int ShadowXMax; /* used to determine maximum extent of beach cells */
}
Deltas_state;

Deltas_state _s;

/* Overall Shoreface Configuration Arrays - Data file information */

//char	AllBeach[Xmax][2*Ymax];		/* Flag indicating of cell is entirely beach */
//float	PercentFull[Xmax][2*Ymax];	/* Fractional amount of shore cell full of sediment */
//int	Age[Xmax][2*Ymax];		/* Age since cell was deposited */
//float	CellDepth[Xmax][2*Ymax];	/* Depth array (m) (ADA 6/3) */


/* Computational Arrays (determined for each time step) */
 
//int 	X[MaxBeachLength];		/* X Position of ith beach element */
//int	Y[MaxBeachLength];		/* Y Position of ith beach element */
//char	InShadow[MaxBeachLength];	/* Is ith beach element in shadow? */
//float	ShorelineAngle[MaxBeachLength];	/* Angle between cell and right (z+1) neighbor	*/
//float	SurroundingAngle[MaxBeachLength];/* Cell-orientated angle based upon left and right neighbor */
//char	UpWind[MaxBeachLength];		/* Upwind or downwind condition used to calculate sediment transport */
//float	VolumeIn[MaxBeachLength];	/* Sediment volume into ith beach element */	
//float 	VolumeOut[MaxBeachLength];	/* Sediment volume out of ith beach element */


/* Miscellaneous Global Variables */

//int	CurrentTimeStep = 0;  	/* Time step of current calculation */ 
///*FindBeachCells,FindNextCell*/int	NextX;			/* Global variables used to iterate FindNextCell in global array - */
///*FindBeachCells,FindNextCell*/int	NextY;			/*	would've used pointer but wouldn't work	*/
//int 	TotalBeachCells;	/* Number of cells describing beach at particular iteration */
//int 	ShadowXMax; 		/* used to determine maximum extent of beach cells */
float 	WaveAngle;		/* wave angle for current time step */	
/* run */int	FindStart;		/* Used to tell FindBeach at what Y value to start looking */
/*run,FindBeachCells*/char	FellOffArray;		/* Flag used to determine if accidentally went off array */
/* init,run */float   MassInitial;		/* For conservation of mass calcs */
/* run */float	MassCurrent;		/* " */
/* unused!!! */int	device;			
/* unused!!! */short	button;
/* unused!!! */long	buttonback;
int	NumWaveBins;		/* For Input Wave - number of bins	*/
float	WaveMax[36];		/* Max Angle for specific bin */
float	WaveProb[36];		/* Probability of Certain Bin */

//float xcellwidth;
//float ycellwidth;
///*ScreenInit*/static WINDOW *mainwnd;
///*ScreenInit*/static WINDOW *screen;
/* unused!!! */WINDOW *my_win;
/* unused!!! */int current_getch;
//int xplotoff;
//int yplotoff;

typedef struct
{
   float xcellwidth;
   float ycellwidth;
   int   xplotoff;
   int   yplotoff;
}
Deltas_graphics;

Deltas_graphics _g;

/*Unused!!! */int KEY_P = 93;


/* Function Prototypes */

void    AdjustShore(int i);
void	AgeCells(void);
void 	ButtonEnter(void);
void	CheckOverwash(int icheck);
void	CheckOverwashSweep(void);
void    DeliverSediment(void);
void	DetermineAngles(void);
void 	DetermineSedTransport(void);
void	DoOverwash(int xfrom,int yfrom, int xto, int yto,  float xintto, float yintto, float distance, int ishore);
void 	FindBeachCells(int YStart);
char 	FindIfInShadow(int icheck, int ShadMax);
void 	FindNextCell(int x, int y, int z);
float 	FindWaveAngle(void);
void	FixBeach(void);
float 	GetOverwashDepth(int xin, int yin, float xintto, float yintto, int ishore);
void    GraphCells(void);
void	InitConds(void);
void	InitPert(void);
float	MassCount(void);
void	OopsImEmpty(int x, int y);
void 	OopsImFull(int x, int y);
void    OpenWindow(void);
void	PauseRun(int x, int y, int in);
void	PeriodicBoundaryCopy(void);
void    PutPixel(float x, float y, float R, float G, float B);
void 	PrintLocalConds(int x, int y, int in);
float 	Raise(float b, float e);
float	RandZeroToOne(void);
void	ReadSandFromFile(void);
void	ReadWaveIn(void);
void	SaveSandToFile(void);
void	SaveLineToFile(void);
void    ScreenInit(void);
void 	SedTrans(int From, int To, float ShoreAngle, char MaxT);
void 	ShadowSweep(void);
void 	TransportSedimentSweep(void);
int 	XMaxBeach(int Max);
void	ZeroVars(void);

/*void	FillUpGap(int X, int Y, int LorR);*/

int deltas_init    ( void );
int deltas_run     ( void );
int deltas_finalize( void );

int
main( void )
{
   deltas_init();
   deltas_run();
   deltas_finalize();

   return EXIT_SUCCESS;
}

#define SEED             (44)  /* random seed  control value = 1 */
#define START_FROM_FILE  ('n') /* start from saved file? */

///* init */char    StartFromFile = 'n';    /* start from saved file? */

int
deltas_init()
{ /* Initialize Variables and Device */
    int seed = 44;
    char StartFromFile = 'n'; /* start from saved file? */

    _s.ShadowXMax = Xmax-5;

    srandom(seed);

    /* Start from file or not? */
    if (PromptStart == 'y')
    {
	printf("shall we start from a file (y or n)? \n");
	scanf("%c", &StartFromFile);

	if (StartFromFile == 'y')
	{
	    printf("Starting Filename? \n");
	    scanf("%24s", _io.readfilename);
	    printf("Saving Filename? \n");
	    scanf("%s", _io.savefilename);
	    printf("What time step are we starting at?");
	    scanf("%d", &_s.CurrentTimeStep);
	    ReadSandFromFile();
	}
	
	if (StartFromFile == 'n')
	{
	    printf("Saving Filename? \n");
	    scanf("%s", _io.savefilename);
	    InitConds();
	    if (InitialPert) 
	    {
		InitPert();
	    }
	    printf("InitConds OK \n");
	}
    }
    else if (StartFromFile == 'y')
    {
	ReadSandFromFile();
    }
    else
    {
	InitConds();
	if (InitialPert) 
	{
	    InitPert();
	}
    }	

#ifdef WITH_OPENGL
    if (KeysOn)
	ScreenInit();
#endif

    /* Count Initial Mass */

    PeriodicBoundaryCopy();
    FixBeach();
    MassInitial = MassCount();

    /*if (SaveLine) 
      SaveLineToFile();
      if (SaveFile) 
      SaveSandToFile();*/

    /* Open Display*/
#ifdef WITH_OPENGL
    if ( DO_GRAPHICS )
    {
	_g.xcellwidth = 2.0 / (2.0 * (float) XPlotExtent) * (CELL_PIXEL_SIZE)/2.0;
	_g.ycellwidth = 2.0 / (2.0 * (float) YPlotExtent) * (CELL_PIXEL_SIZE)/2.0; 
	_g.xplotoff = 0;
	_g.yplotoff = Ymax/2;

	OpenWindow();
     	
	if(EveryPlotSpacing)
	    GraphCells();
    }
#endif


    if( WAVE_IN ) 
	ReadWaveIn();

    return TRUE;
}

#define START_SAVING_AT   (0)    /* time step to begin saving files */
#define SAVE_SPACING      (2500) /* space between saved files */
#define SAVE_LINE_SPACING (1000) /* space between saved line files */
#define SAVE_FILE         (1)    /* save full file? */
#define SAVE_LINE         (0)    /* Save line */

///* run */int     SaveSpacing = 2500 ;    /* space between saved files */
///* run */int     SaveLineSpacing = 1000;/* space between saved line files */
///* run */int     SaveFile = 1;           /* save full file? */
///* run */int     SaveLine = 0;           /* Save line */

#define AGE_UPDATE (10) /* Time space for updating age of non-beach cells */
//int	AgeUpdate = 10;		/* Time space for updating age of non-beach cells */
int
deltas_run( void )
{ /* PRIMARY PROGRAM LOOP */
    int	xx;			/* duration loop variable */
    int StartSavingAt   = START_SAVING_AT;
    int SaveSpacing     = SAVE_SPACING;
    int SaveLineSpacing = SAVE_LINE_SPACING;
    int SaveFile        = SAVE_FILE;
    int SaveLine        = SAVE_LINE;
    int AgeUpdate       = AGE_UPDATE;

    while ( _s.CurrentTimeStep < StopAfter )
    {
	/*  Time Step iteration - compute same wave angle for Duration time steps */

	/*  Calculate Wave Angle */

	WaveAngle = FindWaveAngle();

	/*  Loop for Duration at the current wave sign and wave angle */

	for (xx = 0; xx < Duration; xx++)
	{
		
	    /* Text to Screen? */

	    if (_s.CurrentTimeStep%ScreenTextSpacing == 0)
	    {
		printf("==== WaveAngle: %2.2f  MASS Percent: %1.4f  Time Step: %d\n", 180*(WaveAngle)/M_PI, 
		       MassCurrent/MassInitial, _s.CurrentTimeStep);
	    }

	    PeriodicBoundaryCopy();

	    ZeroVars();
			
	    /* Initialize for Find Beach Cells  (make sure strange beach does not cause trouble */

	    FellOffArray = 'y';
	    FindStart = 1;
		
	    /*  Look for beach - if you fall off of array, bump over a little and try again */

	    while (FellOffArray == 'y')
	    {	
		FindBeachCells(FindStart);
		/*printf("FoundCells: %d GetO = %c \n", FindStart,FellOffArray);*/
		FindStart += FindCellError;
		if (FellOffArray == 'y')
		{
		    /*printf("NOODLE  !!!!!FoundCells: %d GetO = %c \n", FindStart,FellOffArray); */
		    /*PauseRun(1,1,-1);*/
		}
		
		/* Get Out if no good beach spots exist - finish program*/
			
		if (FindStart > Ymax/2+1)
		{
		    printf("Stopped Finding Beach - done %d %d",FindStart,Ymax/2-5); 
		    SaveSandToFile();
		    return 1;
		}
	    }
	
	    /* printf("Foundbeach!: %d \n", _s.CurrentTimeStep); */

	    ShadowSweep();
	    //DEBUG_PRINT( DEBUG_0, "Shadowswept: %d \n", _s.CurrentTimeStep);
            DEBUG_PRINT( DEBUG_0, "Shadowswept: %d \n", _s.CurrentTimeStep );

	    DetermineAngles();
	    //DEBUG_PRINT( DEBUG_0, "AngleDet: %d \n", _s.CurrentTimeStep);
	    DEBUG_PRINT( DEBUG_0, "AngleDet: %d \n", _s.CurrentTimeStep );
	    DetermineSedTransport();
	    DEBUG_PRINT( DEBUG_0, "Sed Trans: %d \n", _s.CurrentTimeStep); 
	    TransportSedimentSweep();
	    DEBUG_PRINT( DEBUG_0, "Transswept: %d \n", _s.CurrentTimeStep);

	    DeliverSediment();

	    FixBeach();
	    DEBUG_PRINT( DEBUG_0, "Fixed Beach: %d \n", _s.CurrentTimeStep);


		/* OVERWASH */
		/* because shoreline config may have been changed, need to refind shoreline and recalc angles*/

			ZeroVars();
			/* Initialize for Find Beach Cells  (make sure strange beach does not cause trouble */

			FellOffArray = 'y';
			FindStart = 1;
			
			/*  Look for beach - if you fall off of array, bump over a little and try again */

			while (FellOffArray == 'y')
			{	
				FindBeachCells(FindStart);
				/*printf("FoundCells: %d GetO = %c \n", FindStart,FellOffArray);*/
				FindStart += FindCellError;
				if (FellOffArray == 'y')
				{
					/*printf("NOODLE  !!!!!FoundCells: %d GetO = %c \n", FindStart,FellOffArray); */
					/*PauseRun(1,1,-1);*/
				}
			
				/* Get Out if no good beach spots exist - finish program*/
				
				if (FindStart > Ymax/2+1)
				{
					printf("Stopped Finding Beach - done %d %d",FindStart,Ymax/2-5); 
					SaveSandToFile();
					return 1;
				}
			}
	
				/* printf("Foundbeach!: %d \n", _s.CurrentTimeStep); */

			ShadowSweep();
				DEBUG_PRINT( DEBUG_0, "Shadowswept: %d \n", _s.CurrentTimeStep);
			DetermineAngles();
				DEBUG_PRINT( DEBUG_0, "AngleDet: %d \n", _s.CurrentTimeStep);
			CheckOverwashSweep();
			FixBeach();
	

	    if ((StartStop)   ) 
	    {
		printf("---- You Paused it, bud ---- \npp");
		PauseRun(1,1,-1);			
	    }

	    DEBUG_PRINT( DEBUG_3, "End of Time Step: %d \n", _s.CurrentTimeStep);
		
	    /* Age Empty Cells */

	    if ((_s.CurrentTimeStep%AgeUpdate == 0) && SaveAge)
		AgeCells();

	    /* Count Mass */

	    MassCurrent = MassCount();

	    /* GRAPHING */		
		
#ifdef WITH_OPENGL
	    if ( DO_GRAPHICS && EveryPlotSpacing && (_s.CurrentTimeStep%EveryPlotSpacing == 0))
		GraphCells();
#endif
	    

	    /* current_getch = getch();
	    printf("%d",current_getch);
if(DoGraphics && KeysOn && (current_getch == KEY_P))
GraphCells();*/

		
	    _s.CurrentTimeStep ++;

	    /* SAVE FILE ? */

	    if (((_s.CurrentTimeStep%SaveSpacing == 0 && _s.CurrentTimeStep >= StartSavingAt) 
		 || (_s.CurrentTimeStep == StopAfter)) && SaveFile)
	    {	
		SaveSandToFile();
	    }

	    if (((_s.CurrentTimeStep%SaveLineSpacing == 0 && _s.CurrentTimeStep >= StartSavingAt) 
		 || (_s.CurrentTimeStep == StopAfter)) && SaveLine)
	    {
		SaveLineToFile();
	    }


			
	    /*if (_s.CurrentTimeStep > 14300)
	      SaveSandToFile();*/
	}				
	
	
    }

    return TRUE;
}

int
deltas_finalize( void )
{
    printf("Run Complete.  Output file: %s\n" , _io.savefilename);
    return TRUE;
}

float FindWaveAngle(void)

/* calculates wave angle for given time step */

{

    float 	Angle;

    /* Data Input Method */

    float 	RandBin;	/* Random number to pick wave angle bin */
    float 	RandAngle;	/* Random number to pick wave angle within the bin */
    int	flag = 1;
    int 	i = 0;


    /*  Variables for Asymmetry Method */

    float	AsymRandom;   /* variable used to determine wave direction for current time step */
    float 	AngleRandom;  /* variable used to determine wave power for current time step */
    int	Sign;	      /* defines direction of incoming waves for current time step, */
    /* positive from left, negative from right */

	
    /* Method using input binned wave distribution - 					*/
    /* variables WaveProb[], WaveMax[], previously input from file using ReadWaveIn()	*/ 
	
    if ( WAVE_IN )
    {
	
	RandBin = RandZeroToOne();
	RandAngle = RandZeroToOne();
	
	/*printf("Time = %d RandBin = %f RandAng = %f\n",_s.CurrentTimeStep, RandBin, RandAngle);*/

	while (flag)
	{
	    i++;		

	    if (RandBin < WaveProb[i])
	    {
		flag = 0;
	    }
	}

	Angle = - ( RandAngle * (WaveMax[i] - WaveMax[i-1]) + WaveMax[i-1])*M_PI/180;

	/*printf("i = %d WaveMAx[i] = %f WaveMax[i-1] = %f WaveProb[i] = %f Angle= %f\n",
	  i,WaveMax[i],WaveMax[i-1],WaveProb[i], Angle*180/pi);*/

    }
    else{

	/* Method using wave probability step function 							*/

	/*  Determine sign */	
	/*  	variable Asym will determine fractional distribution of waves coming from the  		*/
	/*  	positive direction (positive direction coming from left)  -i.e. fractional wave asymmetry */

	AsymRandom = RandZeroToOne();

	if ( AsymRandom <= Asym )
	    Sign = 1;
	else
	    Sign = -1;


	/*  Determine wave angle */
	/*  New formulation - even distribution */
	
	
	AngleRandom = RandZeroToOne();
		
	if (AngleRandom > Highness)
	{
	    Angle = Sign * ((AngleRandom-Highness)/(1-Highness)) * M_PI / 4.0;
	}
	else
	{
	    Angle = Sign * ((AngleRandom/Highness)*M_PI/4.0 + M_PI/4.0);
	}

	/*printf("Highness sub: %f AngleRandom: %f \n", Highness, AngleRandom);
	  printf("WaveAngle sub: %f Angle: %f \n", 180*(Angle)/pi, AngleRandom);*/
	
    }
    Angle = WaveAngleSign*Angle;
    /*Angle =1.59/180.0*pi;*/
    return Angle;	

}


void FindBeachCells(int YStart)

/* Determines locations of beach cells moving from left to right direction 	*/
/* This function will affect and determine the global arrays:  _s.X[] and _s.Y[]	*/
/* This function calls FindNextCell   						*/
/* This will define _s.TotalBeachCells for this time step				*/

{
    int 	y, z, xstart;	/* local iterators */	


    /* Starting at left end, find the x - value for first cell that is 'allbeach' */

    xstart = Xmax -1; y = YStart;
	
    while (_s.AllBeach[xstart][y] == 'n')
    {	
	xstart -= 1;
    }
	
    xstart += 1;			/* Step back to where partially full beach */

    _s.X[0] = xstart; 	_s.Y[0] = YStart;

    DEBUG_PRINT( DEBUG_1, "FirsX: %3d  FrstY: %3d  z: 0 \n", _s.X[0], _s.Y[0]);	

    z = 0;

    while ((_s.Y[z] < 2*Ymax -1) && (z < MaxBeachLength-1)) 
    {
	z++;
	_s.NextX = -2;
	_s.NextY = -2;
			
	FindNextCell(_s.X[z-1], _s.Y[z-1], z-1);
	_s.X[z] = _s.NextX;
	_s.Y[z] = _s.NextY;
			
	DEBUG_PRINT( DEBUG_1, "_s.NextX: %3d  _s.NextY: %3d  z: %d \n", _s.NextX, _s.NextY, z);

	if (_s.PercentFull[_s.X[z]][_s.Y[z]] == 0) 
	{
	    printf("\nFINDBEACH: PercentFull Zero x: %d y: %d\n",_s.X[z],_s.Y[z]);
	    /*PauseRun(_s.X[z],_s.Y[z],z);*/
	}

	/* If return to start point or go off left side of array, going the wrong direction 	*/
	/* Jump off and start again closer to middle						*/

	if ((_s.NextY < 1) || ((_s.NextY == _s.Y[0])&&(_s.NextX==_s.X[0])) || (z > MaxBeachLength -2))
	{
	    /*printf("!!!!!!!Fell Off!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! x = %d !!!!!!!!!!!!!", _s.NextX);*/
	    FellOffArray = 'y';
	    ZeroVars();
	    return;
	}

			

	if (z > MaxBeachLength - 3)
	{
	    printf("????????????  went to end of MaxBeach!! ????");
	}

    }
	
    _s.TotalBeachCells = z; 
    FellOffArray = 'n';		

    DEBUG_PRINT( DEBUG_1, "Total Beach: %d  \n \n", _s.TotalBeachCells); 

}


void FindNextCell(int x, int y, int z)

/* Function to find next cell that is beach moving in the general positive X direction */
/* changes global variables _s.NextX and _s.NextY, coordinates for the next beach cell       */
/* This function will use but not affect the global arrays:  AllBeach [][], _s.X[], and _s.Y[] */

{

    if ( _s.AllBeach[x-1][y] == 'n')
	/* No beach directly beneath cell */
    {  	
	if ( _s.AllBeach[x][y-1] == 'y' && _s.AllBeach[x][y+1] == 'n') 
	    /* If on right side of protuberance */
	{
	    if ( _s.AllBeach[x-1][y-1] == 'y' )
	    {  	/* Move one inshore */	
		_s.NextX = x-1; _s.NextY = y; return;
	    }	
	    else if (_s.AllBeach[x-1][y-1] == 'n' )	/* This is where shadow procedure was */
	    {	/* Back and to the left */
		_s.NextX = x-1; _s.NextY = y-1; return;
	    }
	    printf("Should've found next cell (1): %d, %d \n", x, y);
	    PauseRun(x, y, z); 
	}
		

	else if ( _s.AllBeach[x][y-1] == 'n' && _s.AllBeach[x][y+1] == 'y')	
	    /* If on left side of protuberance */
	{
	    if ( _s.AllBeach[x+1][y+1] == 'n' && _s.AllBeach[x+1][y] == 'n')
		/*  Up and right - move around spit end */
	    {
		_s.NextX = x+1; _s.NextY = y+1; return;			
	    }
		
	    else if (  _s.AllBeach[x+1][y] == 'y')
		/*  On underside of regular or diagonally thin spit */
	    {
		if ( _s.AllBeach[x+1][y-1] == 'n' && _s.AllBeach[x-1][y-1] == 'n' && _s.X[z-1]>x)
		    /* Reaching end of spit - not going in circles */
		{
		    _s.NextX = x-1; _s.NextY = y; return;
		}
		else if (_s.AllBeach[x+1][y-1] == 'n')
		    /* This is reaching end of spit */
		{
		    _s.NextX = x+1; _s.NextY = y-1; return;
		}
		/* Moving along back side of spit */
		{
		    _s.NextX = x; _s.NextY = y-1; return;
		}
	    }
		
	    else if ( _s.AllBeach[x+1][y+1] == 'y')
		/* we know ( _s.AllBeach[x+1][y] == 'n') */
		/* Moving straight up */
		/* NEW - we still don't want to go in */
	    {
		_s.NextX = x+1; _s.NextY = y; return;
	    }
					
	    printf("Should've found next cell (2): %d, %d \n", x, y);
	    PauseRun(x, y ,z);
	}


	if  (_s.AllBeach[x][y-1] == 'n' && _s.AllBeach[x][y+1] == 'n')
	    /* Hanging out - nothing on sides or top - maybe on corner? */
	{
	    if (_s.AllBeach[x-1][y+1] == 'y' && _s.AllBeach[x+1][y] == 'n')
		/* On left corner of protuberence, move right*/
	    {
		_s.NextX = x; _s.NextY = y+1; return;
	    }

	    else if (_s.AllBeach[x+1][y] == 'y' && _s.AllBeach[x+1][y-1] == 'n')
		/* Under protuberance, move around to left and up  */ 
	    {
		_s.NextX = x+1; _s.NextY = y-1; return;
	    }
			
	    else if (_s.AllBeach[x+1][y] == 'y' && _s.AllBeach[x+1][y-1] == 'y')
		/* Under protuberance, move to left */ 
	    {
		_s.NextX = x; _s.NextY = y-1; return;
	    }
	    printf("Should've found next cell (3): %d, %d \n", x, y);
	    PauseRun(x, y, z); 
	}

	else if ( _s.AllBeach[x][y-1] == 'y' && _s.AllBeach[x][y+1] == 'y' )
	    /* thin entrance between spits.  Don't even think about going in there */
	    /* (Similar case to over head and underneath - don't go in */
	    /* check to see which way we were coming in - from below or from side	*/
	{	if (_s.X[z-1] > x)
	    /* coming from above */
	{
	    if (_s.AllBeach[x+1][y+1] == 'n')
		/* Move right and up*/
	    {
		_s.NextX = x+1; _s.NextY = y+1; return;
	    }
	    else if (_s.AllBeach[x+1][y] == 'n')
		/* Straight up*/
	    {
		_s.NextX = x+1; _s.NextY = y; return;
	    }
	    else if (_s.AllBeach[x+1][y-1] == 'n')
		/* Up and left*/
		/* shouldn't need this, this where coming from */
	    {
		_s.NextX = x+1; _s.NextY = y-1; return;
	    }
	}
	else if (_s.X[z-1] < x)
	    /* coming from below */
	{
	    if (_s.AllBeach[x-1][y-1] == 'n')
		/* move down and left*/
	    {
		_s.NextX = x-1; _s.NextY = y-1; return;
	    }
	    else if (_s.AllBeach[x-1][y] == 'n')
		/*move straight down*/
	    {
		_s.NextX = x-1; _s.NextY = y; return;
	    }
	    else if (_s.AllBeach[x-1][y+1] == 'n')
		/*move straight down*/
		/* shouldn't need this, this would be where coming from*/
	    {
		_s.NextX = x-1; _s.NextY = y+1; return;
	    }
	}
	printf("Should've found next cell (3.5): %d, %d \n", x, y);		
	PauseRun(x, y, z); 
	}
		
	printf("Should've found next cell (4): %d, %d \n", x, y);
	PauseRun(x, y, z); 
    }

	
    else if ( _s.AllBeach[x-1][y] == 'y' && _s.AllBeach[x+1][y] == 'n') 
	/* There is beach beneath cell, nothing over the head */
    {
	if ( _s.AllBeach[x][y+1] == 'n')
	    /*  Adjacent Cell to right is vacant */
	{
	    if ( _s.AllBeach[x-1][y+1] == 'y' )
		/* move straight right */
	    {
		_s.NextX = x; _s.NextY = y+1; return;
	    }
	    else if ( _s.AllBeach[x-1][y+1] == 'n' ) 
		/* Move down and to right */
	    {			 
		_s.NextX = x-1; _s.NextY = y+1; return;	
	    }

	    printf("Should've found next cell (5): %d, %d \n", x, y);
	    PauseRun(x, y, z); 
	}
		
	else if ( _s.AllBeach[x][y+1] == 'y')	
	    /*Brad's note : DON'T REALLY NEED TO REPEAT THIS (WORKS SAME IN BOTH CASES) */
	    /* Right neighbor occupied */
	{				
	    if ( _s.AllBeach[x+1][y+1] == 'n' )
		/* Move up and to right */
	    {
		_s.NextX = x+1; _s.NextY = y+1; return;
	    }
	    else if ( _s.AllBeach[x+1][y+1] == 'y')
		/* Move straight up */
	    {
		_s.NextX = x+1; _s.NextY = y; return;
	    }
			
	    printf("Should've found next cell (6): %d, %d \n", x, y);
	    PauseRun(x, y, z); 	
	}

	printf("Should've found next cell (7): %d, %d \n", x, y);
	PauseRun(x, y, z); 
    }


    else if ( (_s.AllBeach[x-1][y] == 'y') && (_s.AllBeach[x+1][y] == 'y')) 
	/* There is beach behind cell, and over the head don't want to go in (will be shadowed anyway */
	/* Need to use last cell to find out if going into left or right enclosure */
    {
	/* Fill up that nasty piece of work */

	/*	FillUpGap(x,y, y-_s.Y[z-1]);*/


	if (_s.Y[z-1] < y)
	    /* Moving towards right, bump up and over the problem */
	{
	    if (_s.AllBeach[x+1][y-1] == 'n')
		/* Move up and to the left */
	    {
		_s.NextX = x+1;  _s.NextY = y-1; return;
	    }
	    else if (_s.AllBeach[x][y-1] == 'n')
		/* Move directly left */	
	    {
		_s.NextX = x;  _s.NextY = y-1; return;
	    }
	    else if (_s.AllBeach[x-1][y-1] == 'n')
		/* Move left and down */	
	    {
		_s.NextX = x-1;  _s.NextY = y-1; return;
	    }
	    printf("Should've found next cell (8): %d, %d \n", x, y);
	    PauseRun(x, y, z); 
	}

	else if (_s.Y[z-1] > y)
	    /* Moving towards left, go back right */
	{
	    if (_s.AllBeach[x-1][y+1] == 'n')
		/* Move down and to the right */
	    {
		_s.NextX = x-1;  _s.NextY = y+1; return;
	    }
	    else if (_s.AllBeach[x][y+1] == 'n')
		/* Move directly right */	
	    {
		_s.NextX = x;  _s.NextY = y+1; return;
	    }
	    else if (_s.AllBeach[x+1][y+1] == 'n')
		/* Move right and up */	
	    {
		_s.NextX = x+1;  _s.NextY = y+1; return;
	    }
	    printf("Should've found next cell (8): %d, %d \n", x, y);
	    PauseRun(x, y, z); 
	}
	printf("Should've found next cell (9): %d, %d \n", x, y);
	PauseRun(x, y, z); 	
    }
    printf("Should've found next cell (10): %d, %d \n", x, y);
    PauseRun(x, y, z); 	
}


/*void FillUpGap(int X, int Y, int LorR)

*  NOT CURRENTLY BEING USED			 *
*  When thin entrance is discovered, fill it up *
	
{
if (_s.AllBeach[X][Y+2*LorR] == 'n')
{
    _s.AllBeach[X][Y+2*LorR] = 'y';
    _s.PercentFull[X][Y+2*LorR] = 1;
    printf("!!!!!!!!!!!!!!!!!\n		FILLEDERUP: %d, %d, %d \n", X, Y, LorR);
}	
}*/


void ShadowSweep(void)
	
/*  Moves along beach and tests to see if cells are in shadow 		*/
/*  This function will use and determine the Global array:  _s.InShadow[]	*/
/*  This function will use and adjust the variable:   _s.ShadowXMax	*/
/*  This function will use but not adjust the variable:  _s.TotalBeachCells */

{

    int	i;

    /* Find maximum extent of beach to use as a limit for shadow searching */
	
    _s.ShadowXMax = XMaxBeach(_s.ShadowXMax) + 3;
	
    DEBUG_PRINT( DEBUG_3, "_s.ShadowXMax: %d   XMaxBeach: %d \n", _s.ShadowXMax, XMaxBeach(_s.ShadowXMax));

    /* Determine if beach cells are in shadow */

    for (i=0;  i <= _s.TotalBeachCells; i++)
    {
	_s.InShadow[i] = FindIfInShadow(i, _s.ShadowXMax);	
    }

}


int XMaxBeach(int Max)

/* Finds extent of beach in x direction					*/
/* Starts searching at a point 3 rows higher than input Max		*/
/* Function returns integer value equal to max extent of 'allbeach'	*/

{
    int xtest, ytest;
	
    xtest = Max + 2;
    ytest = 0;

    while (xtest > 0)
    {
	while (ytest < 2 *Ymax)
	{
	    if (_s.AllBeach[xtest][ytest] == 'y')
	    {
		return xtest;
	    }
			
	    ytest++;
	}

	ytest = 0;

	xtest--;
    }

    printf("***** Should've found Xmax for shadow): %d, %d ***** \n", xtest, ytest);

    return Xmax;

}
	

char FindIfInShadow(int icheck, int ShadMax)

/*  Function to determine if particular cell xin,yin is in shadow 		*/
/*  Returns a character 'y' if yes 'n' if no					*/
/*  New 2/04 - use pixelwise march  - make it faster, more accurate - aa	*/
/*  New 3/04 - correctly take acocunt for sideways and underneath shadows - aa	*/
/*  This function will use but not affect the global arrays: 			*/
/*	_s.AllBeach[][] and _s.PercentFull[][]					*/
/*  This function refers to global variable:  WaveAngle				*/
	
{

    float	slope;			/* search line slope - slope of zero goes staight forward */
    int	ysign;			/* holder for going left or right alongshore */
    float	x,y;			/* holders for 'real' location of x and y */
    float	xin, yin;		/* starting 'real' locations */
    int	xinint, yinint;		/* integer locations for starting location */
    int	xtestint,ytestint;	/* cell looking at */
    float 	xtest,ytest;		/* 'real' location of testing */
    float	xout,yout;		/* used in AllBeach check - exit coordinates */
    int	NextXInt, NextYInt;	/* holder vairables for cell to check */
    float 	Yup, DistanceUp;	/* when going to next x cell, what other values */
    float 	Xside, DistanceSide;	/* when gpoing to next y cell,other values */
    int 	DEBUG_2a = 0;		/* local debuggers */
    int	debug2b = 0;
	
	
    /* convert angle to a slope and the direction of steps */
    /* note that for case of shoreline, positive angle will be minus y direction */
    /*if (icheck == 106) {DEBUG_2a = 1;debug2b=1;}*/

    if (WaveAngle == 0.0)
    {
	/* unlikely, but make sure no div by zero */
	slope = 0.00001;
    }
    else if (fabs(WaveAngle) == 90.0)
    {
	slope = 9999.9;
    }
    else
    {
	slope = fabs(tan(WaveAngle));
    }

    if (WaveAngle > 0)
	ysign = -1;
    else
	ysign = 1;
		
    DEBUG_PRINT( DEBUG_3, "\nI: %d----------x: %d  Y: %d  Wang:  %f Slope: %f sign: %d \n",
			icheck, _s.X[icheck],_s.Y[icheck],WaveAngle*radtodeg,slope, ysign); 
	
    /* 03/04 AA: depending on local orientations, starting point will differ */
    /* so go through scenarios */

    xinint = _s.X[icheck];
    yinint = _s.Y[icheck];

    if (_s.AllBeach[xinint-1][yinint] == 'y' || ((_s.AllBeach[xinint][yinint-1] == 'y') && 
					      (_s.AllBeach[xinint][yinint+1] == 'y')) )
	/* 'regular condition' */
	/* plus 'stuck in the middle' situation (unlikely scenario)*/
    {
	xin = xinint + _s.PercentFull[xinint][yinint];
	yin = yinint + 0.5;
	DEBUG_PRINT( DEBUG_3, "-- Regular xin: %f  yin: %f\n",xin,yin);
    }
    else if (_s.AllBeach[xinint][yinint-1] == 'y')
	/* on right side */
    {
	xin = xinint + 0.5;
	yin = yinint + _s.PercentFull[xinint][yinint];
	DEBUG_PRINT( DEBUG_3, "-- Right xin: %f  yin: %f\n",xin,yin);
    }
    else if (_s.AllBeach[xinint][yinint+1] == 'y')
	/* on left side */
    {
	xin = xinint + 0.5;
	yin = yinint + 1.0 - _s.PercentFull[xinint][yinint];
	DEBUG_PRINT( DEBUG_3, "-- Left xin: %f  yin: %f\n",xin,yin);
    }
    else if (_s.AllBeach[xinint+1][yinint] == 'y')
	/* gotta be on the bottom now */
    {
	xin = xinint + 1 - _s.PercentFull[xinint][yinint];
	yin = yinint + 0.5;
	DEBUG_PRINT( DEBUG_3, "-- Under xin: %f  yin: %f\n",xin,yin);
    }
    else
	/* debug ain't just an insect */
    {
	printf("Shadowstart Broke !!!! ");
	PauseRun(xinint,yinint,icheck);
    }

    x = xin;
    y = yin;
	
    while ((floor(x) < ShadMax) && (y > Ymax/2) && (y < 3*Ymax/2))
    {
	NextXInt = floor(x) + 1;
	if (ysign > 0)
	    NextYInt = floor(y) + 1;
	else
	    NextYInt = ceil(y-1);			 

	/* moving to next whole 'x' position, what is y position? */
	Yup = y + (NextXInt-x)*slope * ysign;
	DistanceUp = ((Yup - y)*(Yup - y) + (NextXInt - x)*(NextXInt - x));
	
	/* moving to next whole 'y' position, what is x position? */
	Xside = x + fabs(NextYInt - y) / slope;
	DistanceSide = ((NextYInt - y)*(NextYInt - y) + (Xside - x)*(Xside - x));
	
	DEBUG_PRINT( DEBUG_3, "x: %f  y: %f  X:%d  Y: %d  Yd: %f  DistD: %f Xs: %f DistS: %f\n",
			    x,y,NextXInt,NextYInt, Yup,DistanceUp,Xside,DistanceSide); 

	if (DistanceUp < DistanceSide)
	    /* next cell is the up cell */
	{
	    x = NextXInt;
	    y = Yup;
	    xtestint = NextXInt;	
	    ytestint = floor(y);
	    DEBUG_PRINT( DEBUG_3, " up ");
	}
	else
	    /* next cell is the side cell */
	{
	    x = Xside;
	    y = NextYInt;
	    xtestint = floor(x);
	    ytestint = y + (ysign-1)/2;
	    DEBUG_PRINT( DEBUG_3, " side ");
	}
		
	DEBUG_PRINT( DEBUG_3, "	x: %f  y: %f  xtesti: %d ytesti: %d \n\n",x,y,xtestint,ytestint); 
			
			
	/* Now Test */
	/* If AllBeach is along the way, will we pass through 'diamond'?	*/
	/* Trick - if crossing through the diamond, will change quadrants 	*/
	/* Probably won't get to this one, though			 	*/		

	if 	(_s.AllBeach[xtestint][ytestint] == 'y')
	{
	    /* use same approach to find exit (could make this modular)	 */
	    /* don't change 'x' or 'y' and this will be ok			*/
	    NextXInt = floor(x) + 1;
	    if (ysign > 0)
		NextYInt = floor(y) + 1;
	    else
		NextYInt = ceil(y-1);	

	    Yup = y + (NextXInt-x)*slope * ysign;
	    DistanceUp = ((Yup - y)*(Yup - y) + (NextXInt - x)*(NextXInt - x));
	    Xside = x + fabs(NextYInt - y) / slope;
	    DistanceSide = ((NextYInt - y)*(NextYInt - y) + (Xside - x)*(Xside - x));

	    if (DistanceUp < DistanceSide)
		/* next cell is the up cell */
	    {
		xout = NextXInt;
		yout = Yup;
	    }
	    else
		/* next cell is the side cell */
	    {
		xout = Xside;
		yout = NextYInt;
	    }

	    /*DEBUG_PRINT( DEBUG_3, "In Allbeach xin: %2.2f yin: %2.2f xout: %2.2f yout: %2.2f\n",
	      x,y,xout,yout);
	      DEBUG_PRINT( DEBUG_3, "In Allbeach xin: %2.2f yin: %2.2f xout: %2.2f yout: %2.2f\n",
	      (xout-xtestint-0.5),(x-xtestint-0.5),(yout-ytestint-0.5),(y-ytestint-0.5));*/

	    if(( (xout-xtestint-0.5) * (x-xtestint-0.5) < 0 ) || ((yout-ytestint-0.5) * (y-ytestint-0.5) < 0)) 
	    {
		DEBUG_PRINT( DEBUG_3, "  Shaddowded ");
		return 'y';
	    }
	}	

	

	/* Compare a partially full cell's x - distance to a line projected  	*/
	/* from the starting beach cell's x-distance				*/
	/* This assumes that beach projection is in x-direction (not too bad) 	*/

	else if ( _s.PercentFull[xtestint][ytestint] > 0 ) 
	{
	    if (_s.AllBeach[xtestint-1][ytestint] == 'y' || ((_s.AllBeach[xtestint][ytestint-1] == 'y') && 
							  (_s.AllBeach[xtestint][ytestint+1] == 'y')) )
		/* 'regular' condition */
		/* plus 'stuck in the middle' situation (unlikely scenario) */
	    {
		xtest = xtestint + _s.PercentFull[xtestint][ytestint];
		ytest = ytestint + 0.5;

		if (xtest > (xin  + fabs(ytest-yin)/slope) ) 
		{
		    if (debug2b) printf("Top: sl: %f xt: %2.2f xin: %2.2f yt: %2.2f yin: %2.2f comp: %2.2f > Thing: %2.2f\n",
					slope, xtest, xin, ytest, yin, xtest, (xin  + fabs(ytest-yin)/slope));
		    return 'y';
		}
	    }
	    else if (_s.AllBeach[xtestint][ytestint-1] == 'y')
		/* on right side */
	    {
		xtest = xtestint + 0.5;
		ytest = ytestint + _s.PercentFull[xtestint][ytestint];

		if (ytest > (yin + (xtest-xin) * slope))
		{
		    if (debug2b) printf("Right:  xt: %f  yt: %f  comp: %f > Thing: %f\n",
					xtest, ytest, ytest,(yin + (xtest-xin) * slope));
		    return 'y';
		}
	    }
	    else if (_s.AllBeach[xtestint][ytestint+1] == 'y')
		/* on left side */
	    {
		xtest = xtestint + 0.5;
		ytest = ytestint + 1.0 - _s.PercentFull[xtestint][ytestint];

		if (ytest < (yin + (xtest-xin) * slope))
		{
		    if (debug2b) printf("Left:  xt: %f  yt: %f  comp: %f < Thing: %f\n",
					xtest, ytest, ytest,(yin + (xtest-xin) * slope));
		    return 'y';
		}
	    }
	    else if (_s.AllBeach[xtestint+1][ytestint] == 'y')
		/* gotta be on the bottom now */
	    {
		xtest = xtestint + 1 - _s.PercentFull[xtestint][ytestint];
		ytest = ytestint + 0.5;
				
		if (xtest < (xin  + fabs(ytest-yin)/slope) ) 
		{
		    if (debug2b) printf("Bottom:  xt: %f  yt: %f  comp: %f < Thing: %f\n",
					xtest, ytest, xtest,(xin  + fabs(ytest-yin)/slope));
						
		    return 'y';		
		}
	    }
	    else
		/* debug ain't just an insect */
	    {
		printf("'Shaddows' not responding xin: %f yin: %f xt: %f  yt: %f  \n",
		       xin,yin,xtest, ytest);
		/*PauseRun(xtestint,ytestint,icheck);*/
	    }

	}
    }	
    return 'n';
}


void  DetermineAngles(void)
	
/*  Function to determine beach angles for all beach cells from left to right		*/
/*  By convention, the ShorelineAngle will apply to current cell and right neighbor	*/
/*  This function will determine global arrays:						*/
/*		_s.ShorelineAngle[], _s.UpWind[], _s.SurroundingAngle[]						*/
/*  This function will use but not affect the following arrays and values:		*/
/*		_s.X[], _s.Y[], _s.PercentFull[][], _s.AllBeach[][], WaveAngle			*/
/*  ADA Revised underside, SurroundingAngle 6/03, 2/04 fixed 				*/
/*  ADA Revised angle calc 5/04								*/

{

    int i,j,k;  			/* Local loop variables */
    int x2int, y2int;		/* shoreline location vars */
    float x2,x1,y2,y1;		/* shoreline location variables */
    int debug3a = 0;		/* local debuggr */

	
    /* Shoreline Angle Calcs - ADA 05/04 - use correct 'point' to do calcs (like in shaddow */
    /* Set first point									*/
    /* first angle should be regular one - periodic BC's should also take care 		*/
	
    x2 = _s.X[0] + _s.PercentFull[_s.X[0]][_s.Y[0]];
    y2 = _s.Y[0] + 0.5;

    /* Compute _s.ShorelineAngle[]  */
    /* 	not equal to _s.TotalBeachCells because angle between cell and rt neighbor */

    for (i=0 ; i < _s.TotalBeachCells ; i++)
    {
		
	x1 = x2;
	y1 = y2;

	x2int = _s.X[i+1];
	y2int = _s.Y[i+1];

	if (_s.AllBeach[x2int-1][y2int] == 'y' || ((_s.AllBeach[x2int][y2int-1] == 'y') && 
						(_s.AllBeach[x2int][y2int+1] == 'y')) && (_s.AllBeach[x2int+1][y2int] == 'n'))
	    /* 'regular condition' - if between  */
	    /* plus 'stuck in the middle' situation (unlikely scenario)*/
	{
	    x2 = x2int + _s.PercentFull[x2int][y2int];
	    y2 = y2int + 0.5;
	    if (debug3a) printf("-- Regular xin: %f  yin: %f\n",x2,y2);
	}
	else if ((_s.AllBeach[x2int+1][y2int] == 'y') && (_s.AllBeach[x2int-1][y2int] == 'y'))
	    /* in a sideways nook (or is that a cranny?) */
	{
	    x2 = x2int + 0.5;

	    if (_s.AllBeach[x2int][y2int-1] == 'y')
		/* right-facing nook */
	    {
		y2 = y2int + _s.PercentFull[x2int][y2int];
	    }
	    else
		/* left-facing nook */
	    {
		y2 = y2int + 1.0 - _s.PercentFull[x2int][y2int];
	    }
	    if (debug3a) printf("-- Nook  xin: %f  yin: %f\n",x2,y2);
	}
	else if (_s.AllBeach[x2int][y2int-1] == 'y')
	    /* on right side */
	{
	    x2 = x2int + 0.5;
	    y2 = y2int + _s.PercentFull[x2int][y2int];
	    if (debug3a) printf("-- Right xin: %f  yin: %f\n",x2,y2);
	}
	else if (_s.AllBeach[x2int][y2int+1] == 'y')
	    /* on left side */
	{
	    x2 = x2int + 0.5;
	    y2 = y2int + 1.0 - _s.PercentFull[x2int][y2int];
	    if (debug3a) printf("-- Left xin: %f  yin: %f\n",x2,y2);
	}
	else if (_s.AllBeach[x2int+1][y2int] == 'y')
	    /* gotta be on the bottom now */
	{
	    x2 = x2int + 1 - _s.PercentFull[x2int][y2int];
	    y2 = y2int + 0.5;
	    if (debug3a) printf("-- Under xin: %f  yin: %f\n",x2,y2);
	}
	else
	    /* debug ain't just an insect */
	{
	    printf("Shadowstart Broke !!!! ");
	    PauseRun(x2int,y2int,i+1);
	}


	/* compute angles */

	if (y2 > y1)
	{
	    _s.ShorelineAngle[i] = atan((x2 - x1) / (y2 - y1));
	    DEBUG_PRINT( DEBUG_3, "(R) i = %d  _s.X[i]: %d _s.Y[i]: %d Percent %3f x: %f y: %f Angle:%f  Deg Angle: %f \n",
			       i, _s.X[i], _s.Y[i], _s.PercentFull[_s.X[i]][_s.Y[i]],x2,y2,_s.ShorelineAngle[i], _s.ShorelineAngle[i]*180/M_PI);
	}
	else if (y2 == y1)
	{
	    _s.ShorelineAngle[i] = M_PI/2.0 * (x1 - x2) / fabs(x2 - x1);
	    DEBUG_PRINT( DEBUG_3, "(G) i = %d  _s.X[i]: %d _s.Y[i]: %d Percent %3f x: %f y: %f Angle:%f  Deg Angle: %f \n",
			       i, _s.X[i], _s.Y[i], _s.PercentFull[_s.X[i]][_s.Y[i]],x2,y2,_s.ShorelineAngle[i], _s.ShorelineAngle[i]*180/M_PI);
	}
	else 
	    /* y2 < y1 */
	{
	    _s.ShorelineAngle[i] = atan((x2 - x1) / (y2 - y1)) - M_PI;
		
	    if (_s.ShorelineAngle[i] < - M_PI)
	    {
		_s.ShorelineAngle[i] += 2.0 * M_PI;
	    }
	    DEBUG_PRINT( DEBUG_3, "(U) i = %d  _s.X[i]: %d _s.Y[i]: %d Percent %3f x: %f y: %f Angle:%f  Deg Angle: %f \n",
			       i, _s.X[i], _s.Y[i], _s.PercentFull[_s.X[i]][_s.Y[i]],x2,y2,_s.ShorelineAngle[i], _s.ShorelineAngle[i]*180/M_PI);
	}

    }

    for (k=1 ; k < _s.TotalBeachCells ; k++)
    {
	/* compute SurroundingAngle array */
	/* 02/04 AA averaging doesn't work on bottom of spits */
	/* Use trick that x is less if on bottom of spit - angles might be different signs as well */
		
	if ((_s.Y[k-1] - _s.Y[k+1] == 2) && 
	    (copysign(_s.ShorelineAngle[k-1],_s.ShorelineAngle[k]) != _s.ShorelineAngle[k-1]))
	{		
	    _s.SurroundingAngle[k] = (_s.ShorelineAngle[k-1] + _s.ShorelineAngle[k]) / 2 + M_PI;
	    if (_s.SurroundingAngle[k] > M_PI)
	    {
		_s.SurroundingAngle[k] -= 2.0 * M_PI;
	    }
	    DEBUG_PRINT( DEBUG_4, "Under: %d\n",k);
	}
	else
	{
	    _s.SurroundingAngle[k] = (_s.ShorelineAngle[k-1] + _s.ShorelineAngle[k]) / 2;
	}
    }
	
    /* Determine Upwind/downwind condition						*/
    /* Note - Surrounding angle is based upon left and right cell neighbors, 	*/
    /* and is centered on cell, not on right boundary				*/
	

    DEBUG_PRINT( DEBUG_4, "\nUp/Down   Wave Angle:%f\n", WaveAngle * radtodeg);

    for (j=1 ; j < _s.TotalBeachCells  ; j++)
    {
	DEBUG_PRINT( DEBUG_4, "i: %d  Shad: %c Ang[i]: %3.1f  Sur: %3.1f  Effect: %3f  ",
			   j,_s.InShadow[j], _s.ShorelineAngle[j]*radtodeg, 
			   _s.SurroundingAngle[j]*radtodeg, (WaveAngle - _s.SurroundingAngle[j])*radtodeg);

	if ( fabs(WaveAngle - _s.SurroundingAngle[j]) >= 42.0/radtodeg )
	{	
	    _s.UpWind[j] = 'u';
	    DEBUG_PRINT( DEBUG_4, "U(1)  ");
	}
	else 
	{
	    _s.UpWind[j] = 'd';
	    DEBUG_PRINT( DEBUG_4, "D(1)  ");
	}

	DEBUG_PRINT( DEBUG_4, "\n");

    }

}

#define SED_TRANS_LIMIT (90) /* beyond what absolute slope don't do sed trans (degrees)*/
//float   SedTansLimit =  90;	/* beyond what absolute slope don't do sed trans (degrees)*/

void DetermineSedTransport(void)

/*  Loop function to determine which neigbor/situation to use for sediment transport calcs	*/
/*  Once situation is determined, will use function SedTrans to determine actual transport	*/
/*  This function will call SedTrans which will determine global arrays:			*/
/*		_s.VolumeIn[], _s.VolumeOut[]								*/
/*  This function will use but not affect the following arrays and values:			*/
/*		_s.X[], _s.Y[], _s.InShadow[], _s.UpWind[], _s.ShorelineAngle[]				*/
/*  		_s.PercentFull[][], _s.AllBeach[][], WaveAngle					*/

{

    int i;			/* Loop variable */
    float ShoreAngleUsed;	/* Temporary holder for shoreline angle 				*/
    int CalcCell;		/* Cell sediment coming from to go across boundary i 			*/
    int Next,Last;		/* Indicators so test can go both left/right			 	*/
    int Correction;		/* Term needed for shoreline angle and i+1 case, angle stored at i 	*/
    char UpWindLocal;	/* Local holder for upwind/downwind condition				*/
    char MaxTrans;		/* Do we need to compute using maximum transport ? 			*/
    int DoFlux;		/* Skip sed transport calcs (added 02/04 AA)				*/
    float   SedTansLimit =  SED_TRANS_LIMIT;


    DEBUG_PRINT( DEBUG_5, "\nSEDTRANS: %d  @  %f \n\n", _s.CurrentTimeStep, WaveAngle * radtodeg);

    for (i=1 ; i < _s.TotalBeachCells-1 ; i++)
    {
	DEBUG_PRINT( DEBUG_5, "\n  i: %d  ",i);		


	MaxTrans = 'n';

	/*  Is littoral transport going left or right?	*/

	if ((WaveAngle-_s.ShorelineAngle[i]) > 0)
	{
	    /*  Transport going right, center on cell to left side of border	 */
	    /*  Next cell in positive direction, no correction term needed 		*/
	    CalcCell = i;
	    Next = 1;
	    Last = -1;
	    Correction = 0;

	    DEBUG_PRINT( DEBUG_5, "RT  %d ",CalcCell);
	}
	else
	{
	    /*  Transport going left, center on cell to right side of border 	*/
	    /*  Next cell in negative direction, correction term needed 		*/
	    CalcCell = i+1;
	    Next = -1;
	    Last = 1;
	    Correction = -1;

	    DEBUG_PRINT( DEBUG_5, "LT  %d ",CalcCell);
	}
			

	if ( _s.InShadow[CalcCell] == 'n') 
	{
			
	    /*  Adjustment for maximum transport when passing through 45 degrees		*/
	    /*  This adjustment is only made for moving from downwind to upwind conditions	*/
	    /* 										*/
	    /*  purposefully done before shadow adjustment, only use maxtran when 		*/
	    /*	transition from dw to up not because of shadow				*/
	    /* keeping transition from uw to dw - does not seem to be big deal (04/02 AA) */
			
	    if ( ((_s.UpWind[CalcCell] == 'd') && (_s.UpWind[CalcCell+Next] == 'u') &&
		  (_s.InShadow[CalcCell + Next] == 'n')) ||
		 ((_s.UpWind[CalcCell+Last] == 'u') && (_s.UpWind[CalcCell] == 'd')
		  && (_s.InShadow[CalcCell+Last] == 'n')) )
	    {
		MaxTrans = 'y';
		DEBUG_PRINT( DEBUG_5, "MAXTRAN  ");
	    } 
			

	    /*  Upwind/Downwind adjustment Make sure sediment is put into shadows		*/
	    /*  If Next cell is in shadow, use UpWind condition				*/
			
	    DoFlux = 1;
	    UpWindLocal = _s.UpWind[CalcCell];

	    if (_s.InShadow[CalcCell+Next] == 'y')  
	    {
		UpWindLocal = 'u';
		DEBUG_PRINT( DEBUG_5, "U(2)  ");
	    }

	    /*  If coming out of shadow, downwind should be used		*/
	    /*  HOWEVER- 02/04 AA - if high angle, will result in same flux in/out problem */
	    /*  	solution  - no flux for high angle waves */
		
	    if ((_s.InShadow[CalcCell+Last] == 'y') &&(UpWindLocal == 'u')) 
	    {
		DoFlux = 0;
		DEBUG_PRINT( DEBUG_5, "U(X) NOFLUX \n");
				
	    }

	    /*  Use upwind or downwind shoreline angle for calcs			*/

	    if (UpWindLocal == 'u')
	    {
		ShoreAngleUsed = _s.ShorelineAngle[CalcCell+Last+Correction];
		DEBUG_PRINT( DEBUG_5, "UP  ShoreAngle: %3.1f  ", ShoreAngleUsed * radtodeg);
	    }	
	    else if (UpWindLocal == 'd')
	    {
		ShoreAngleUsed = _s.ShorelineAngle[CalcCell+Correction];
		DEBUG_PRINT( DEBUG_5, "DN  ShoreAngle: %3.1f  ", ShoreAngleUsed *radtodeg);
	    }
			

	    /* !!! Do not do transport on unerneath c'cause it gets all messed up */
	    if (fabs(ShoreAngleUsed) > SedTansLimit/radtodeg)
	    {
		DoFlux = 0;
	    }
	
	    /* Send to SedTrans to calculate VolumeIn and VolumeOut*/


	    /* printf("i = %d  Cell: %d NextCell: %d Angle: %f Trans Angle: %f\n",
	       i, CalcCell, CalcCell+Next, ShoreAngleUsed*180/pi, (WaveAngle - ShoreAngleUsed)*180/pi); */

	    DEBUG_PRINT( DEBUG_5, "From: %d  To: %d  TransAngle %3.1f", CalcCell, CalcCell+Next, 
			       (WaveAngle - ShoreAngleUsed) * radtodeg);

	    if (DoFlux)
	    {
		SedTrans(CalcCell, CalcCell+Next, ShoreAngleUsed, MaxTrans);
	    }
	}

    }

}


void SedTrans(int From, int To, float ShoreAngle, char MaxT)

/*  This central function will calcualte the sediment transported from the cell at From to	*/
/*  the cell at To, using the input ShoreAngle							*/
/*  This function will caluclate and determine the global arrays:				*/
/*		_s.VolumeIn[] and _s.VolumeOut[]							*/
/*  This function does not use any other arrays							*/
/*  This function will use the global values defining the wave field:				*/
/*	WaveAngle, Period, OffShoreWvHt								*/
/*  Revised 6/02 - New iterative calc for refraction and breaking, parameters revised		*/
{

    /* Coefficients - some of these are important*/

    float	StartDepth = 3*OffShoreWvHt;/* m, depth to begin refraction calcs (needs to be beyond breakers)	*/
    float 	RefractStep = .2;	/* m, step size to iterate depth for refraction calcs			*/
    float 	KBreak = 0.5;		/* coefficient for wave breaking threshold 				*/
    float 	rho = 1020;		/* kg/m3 - density of water and dissolved matter				*/
	
    /* Variables */

    int 	Broken = 0;		/* is wave broken yet?				*/
    float	AngleDeep;		/* rad, Angle of waves to shore at inner shelf  */
    float 	Depth = StartDepth;	/* m, water depth for current iteration		*/
    float	Angle;			/* rad, calculation angle			*/
    float	CDeep;			/* m/s, phase velocity in deep water		*/
    float	LDeep;			/* m, offhsore wavelength			*/
    float	C;			/* m/s, current step phase velocity 		*/
    float	kh;			/* wavenumber times depth			*/
    float	n;			/* n						*/
    float	WaveLength;		/* m, current wavelength			*/
    float	WvHeight;		/* m, current wave height			*/
    float	VolumeAcrossBorder;	/* m3/day 					*/


    /* Primary assumption is that waves refract over shore-parallel contours			*/
    /* New algorithm 6/02 iteratively takes wiave onshore until they break, then computes Qs	*/
    /* See notes 06/05/02										*/
	
    DEBUG_PRINT( DEBUG_6, "Wave Angle %2.2f Shore Angle  %2.2f    ",WaveAngle*radtodeg, ShoreAngle*radtodeg);

    AngleDeep = WaveAngle - ShoreAngle;

    if (MaxT == 'y')
    {
	AngleDeep = 42.0 / radtodeg;
    }
    DEBUG_PRINT( DEBUG_6, "Deep Tranport Angle %2.2f \n\n",AngleDeep*radtodeg);	

    /*  Don't do calculations if over 90 degrees, should be in shadow  */
	
    if (AngleDeep > 0.995*M_PI/2.0 || AngleDeep < -0.995*M_PI/2.0)
    {
	return;
    }

    else
    {
	/* Calculate Deep Water Celerity & Length, Komar 5.11 c = gT / pi, L = CT	*/
		
	CDeep = GRAV * Period / (2.0 * M_PI);
	LDeep = CDeep * Period;
	DEBUG_PRINT( DEBUG_6, "CDeep = %2.2f LDeep = %2.2f \n",CDeep, LDeep);

	while(!Broken)
	{
	    /* non-iterative eqn for L, from Fenton & McKee 		*/	

	    WaveLength = LDeep * Raise(tanh(Raise(Raise(2.0*M_PI/Period,2)*Depth/GRAV,.75)),2.0/3.0);
	    C = WaveLength/Period;
	    DEBUG_PRINT( DEBUG_6, "DEPTH: %2.2f Wavelength = %2.2f C = %2.2f ", Depth, WaveLength,C);
			
	    /* Determine n = 1/2(1+2kh/tanh(kh)) Komar 5.21			*/
	    /* First Calculate kh = 2 pi Depth/L  from k = 2 pi/L		*/

	    kh =  M_PI * Depth / WaveLength;
	    n =0.5 * ( 1 + 2.0 * kh / sinh(2.0*kh));
	    DEBUG_PRINT( DEBUG_6, "kh: %2.3f  n: %2.3f ", kh, n);

	    /* Calculate angle, assuming shore parallel contours and no conv/div of rays 	*/
	    /* from Komar 5.47								*/

	    Angle = asin(C/CDeep * sin(AngleDeep));
	    DEBUG_PRINT( DEBUG_6, "Angle: %2.2f",Angle*radtodeg);

	    /* Determine Wave height from refract calcs - Komar 5.49			*/
	
	    WvHeight = OffShoreWvHt * Raise(CDeep*cos(AngleDeep)/(C*2.0*n*cos(Angle)),.5);
	    DEBUG_PRINT( DEBUG_6, " WvHeight : %2.3f\n",WvHeight);

	    if (WvHeight > Depth*KBreak)
		Broken = 1;
	    else if (Depth == RefractStep)
	    {
		Broken = 1;
		Depth -= RefractStep;
	    }			
	    else 
		Depth -= RefractStep;
	}

		
	/* Now Determine Transport */
	/* eq. 9.6b (10.8) Komar, including assumption of sed density = 2650 kg/m3  		*/
	/* additional accuracy here will not improve an already suspect eqn for sed transport 	*/
	/* (especially with poorly constrained coefficients), 					*/
	/* so no attempt made to make this a more perfect imperfection				*/
		
	VolumeAcrossBorder = 	fabs(1.1*rho*Raise(GRAV,3.0/2.0)*Raise(WvHeight,2.5)*
				     cos(Angle)*sin(Angle)*TimeStep);
				
	_s.VolumeOut[From] = _s.VolumeOut[From] + VolumeAcrossBorder;
		
	_s.VolumeIn[To] = _s.VolumeIn[To] + VolumeAcrossBorder;
	
	DEBUG_PRINT( DEBUG_6, "VolumeAcrossBorder: %f  ",VolumeAcrossBorder);
	DEBUG_PRINT( DEBUG_6, "VolumeIn : %f ",_s.VolumeIn[To]);
	DEBUG_PRINT( DEBUG_6, "VolumeOut : %f \n\n",_s.VolumeOut[From]);
		
    }
}


void TransportSedimentSweep(void)

/*  Sweep through cells to place transported sediment				*/
/*  Call function AdjustShore() to move sediment.  				*/
/*  If cell full or overempty, call OopsImFull or OopsImEmpty()			*/
/*  This function doesn't change any values, but the functions it calls do	*/
/*  Uses but doesn't change:  _s.X[], _s.Y[], _s.PercentFull[]				*/
/*  sweepsign added to ensure that direction of actuating changes does not  	*/
/*  	produce unwanted artifacts (e.g. make sure symmetrical			*/

{

    int i,ii;
    int sweepsign;

    if (RandZeroToOne()*2 > 1)
    {
	sweepsign = 1;
	DEBUG_PRINT( DEBUG_7A, "L  ");
    }
    else
    {
	sweepsign = 0;
	DEBUG_PRINT( DEBUG_7A, "R  ");
    }

    DEBUG_PRINT( DEBUG_7A, "\n\n TransSedSweep  Ang %f  %d\n", WaveAngle * radtodeg, _s.CurrentTimeStep);
	
    for (i=0; i < _s.TotalBeachCells-1 ; i++)
    {
	
	if (sweepsign == 1)
	    ii = i;
	else
	    ii = _s.TotalBeachCells-1-i;

	DEBUG_PRINT( DEBUG_7A, "i: %d  ss: %d  X: %d  Y: %d  In: %.1f  Out: %.1f\n", ii, sweepsign,
			   _s.X[i], _s.Y[i], _s.VolumeIn[i], _s.VolumeOut[i]);

	AdjustShore(ii);
				
	if (_s.PercentFull[_s.X[ii]][_s.Y[ii]] < 0)
	{
	    OopsImEmpty(_s.X[ii],_s.Y[ii]);
	}
	else if (_s.PercentFull[_s.X[ii]][_s.Y[ii]]> 1) 
	{
	    OopsImFull(_s.X[ii],_s.Y[ii]);
	}
    }

}

void AdjustShore(int i)

/*  Complete mass balance for incoming and ougoing sediment			*/
/*  This function will change the global data array _s.PercentFull[][]		*/
/*  Uses but does not adjust arrays:  						*/
/*		_s.VolumeIn[], _s.VolumeOut[], _s.X[], _s.Y[], _s.ShorelineAngle[]		*/
/*  Uses global variables: ShelfSlope, CellWidth, ShorefaceSlope, InitialDepth	*/
/*  NEW - AA 05/04 fully utilize shoreface depths				*/

{

    float	Depth;		/* Depth of convergence*/
    float	DeltaArea;	/* Holds change in area for cell (m^2)*/
    float	Distance;	/* distance from shore to intercept of equilib. profile and overall slope (m)*/
    float	PercentIn;
    float 	PercentOut;
    float	PercentSum;
    int	Xintint, Yintint;	/* integer representing location shoreface cell */
    float	Xintfloat,Yintfloat;	/* floaters for shoreface cell */
	
    /* variables for loop */
    float	slope;			/* slope of zero goes staight back */
    int	ysign;			/* holder for going left or right alongshore */
    float	x,y;			/* holders for 'real' location of x and y */
    int	xtest,ytest;		/* cell looking at */
    int	NextXInt, NextYInt;	/* holder vairables for cell to check */
    float 	Ydown, DistanceDown;	/* when going to next x cell, what other values */
    float 	Xside, DistanceSide;	/* when gpoing to next y cell,other values */
    int 	ShorefaceFlag;		/* flag to see if started intersecting shoreface cells */

    if (_s.VolumeIn[i] <= _s.VolumeOut[i])
	/* eroding, just have to use shoreface depth */
    {
	Depth = DepthShoreface;
    }
    else
	/* accreting, good god */
    {
	/* where should we intersect shoreface depth ? */
	
	/* uncomplicated way - assume starting in middle of cell */
	Distance = DepthShoreface/CellWidth/ShorefaceSlope;
	Xintfloat = _s.X[i] + 0.5 + Distance * cos(_s.SurroundingAngle[i]);
	Xintint = floor(Xintfloat);
	Yintfloat = _s.Y[i] + 0.5 - Distance * sin(_s.SurroundingAngle[i]);
	Yintint = floor(Yintfloat);

	DEBUG_PRINT( DEBUG_7A, "xs: %d  ys: %d  Xint: %f Xint:%d Yint: %f Yint: %d  Dint: %f SAng: %f Sin = %f\n",
			   _s.X[i],_s.Y[i],Xintfloat,Xintint,Yintfloat,Yintint,_s.CellDepth[Xintint][Yintint],_s.SurroundingAngle[i]*radtodeg,sin(_s.SurroundingAngle[i]));


	if ((Yintint < 0) || (Yintint > 2*Ymax))
	{
	    Depth = DepthShoreface;
	    if ((Yintint > Ymax/2) && (Yintint < 3/2*Ymax))
	    {
		printf("Periodic Boundary conditions and Depth Out of Bounds");
		PauseRun(_s.X[i],_s.Y[i],i);
	    }
	}
	else if ((Xintint < 0) || (Xintint > Xmax))
	{
	    Depth = DepthShoreface;
	    printf("-- Warning - depth location off of x array: X %d Y %d",Xintint,Yintint);
	    PauseRun(_s.X[i],_s.Y[i],i);
	}
	else if (_s.CellDepth[Xintint][Yintint] <= 0)
	    /* looking back on land */
	{
	    Depth = DepthShoreface;
	    DEBUG_PRINT( DEBUG_7A, "=== Shoreface is Shore, eh? Accreti:  xs: %d  ys: %d  Xint:%d  Yint: %d  Dint: %f \n",
				_s.X[i],_s.Y[i],Xintint,Yintint,_s.CellDepth[Xintint][Yintint]);
	}
	else if (_s.CellDepth[Xintint][Yintint] < DepthShoreface)
	{
	    printf("Shallow but underwater Depth %f",_s.CellDepth[Xintint][Yintint]);
	    PauseRun(Xintint,Yintint,01);
	}
	else
	{
	    Depth = _s.CellDepth[Xintint][Yintint];
			

	    /* That was the easy part - now we need to 'fix' all cells towards shoreface */
	    /* probably due to accretion from previous moving forward */
	    /* reuse some of the overwash checking code here */


	    if (_s.SurroundingAngle[i] == 0)
	    {
		/* unlikely, but make sure no div by zero */
		slope = 0.00001;
	    }
	    else if (fabs(_s.SurroundingAngle[i]) == 90.0)
	    {
		slope = 9999.9;
	    }
	    else
	    {
		slope = fabs(tan(_s.SurroundingAngle[i]));
	    }

	    if (_s.SurroundingAngle[i] > 0)
		ysign = 1;
	    else
		ysign = -1;
			
	    x = Xintfloat;
	    y = Yintfloat;
	    xtest = Xintint;
	    ytest = Yintint;
	    ShorefaceFlag = 0;	
	
	    while (( _s.CellDepth[xtest][ytest] > DepthShoreface) && !(ShorefaceFlag))
	    {
		NextXInt = ceil(x) -1;	
		if (ysign > 0)	
		    NextYInt = floor(y) + 1;
		else	
		    NextYInt = ceil(y-1);			 

		/* moving to next whole 'x' position, what is y position? */
		Ydown = y + (x - NextXInt)*slope * ysign;
		DistanceDown = Raise(((Ydown - y)*(Ydown - y) + (NextXInt - x)*(NextXInt - x)),.5);
		
		/* moving to next whole 'y' position, what is x position? */
		Xside = x - fabs(NextYInt - y) / slope;
		DistanceSide = Raise(((NextYInt - y)*(NextYInt - y) + (Xside - x)*(Xside - x)),.5);
	
				

		if (DistanceDown < DistanceSide)
		    /* next cell is the down cell */
		{
		    x = NextXInt;
		    y = Ydown;
		    xtest = NextXInt-1;	
		    ytest = floor(y);			
		}
		else
		    /* next cell is the side cell */
		{
		    x = Xside;
		    y = NextYInt;
		    xtest = floor(x);
		    ytest = y + (ysign-1)/2;			
		}
						
		if (_s.CellDepth[xtest][ytest] > DepthShoreface)	
		    /* Deep hole - fill 'er in - mass came from previous maths */ 
		{

		    DEBUG_PRINT( DEBUG_7A, "=== Deep Hole, eh? Accreti:  xs: %d  ys: %d  Xint:%d  Yint: %d  Dint: %f Xfill: %d Yfill: %d Dt: %f\n",
					_s.X[i],_s.Y[i],Xintint,Yintint,_s.CellDepth[Xintint][Yintint],xtest,ytest,
					_s.CellDepth[xtest][ytest]); 
		    _s.CellDepth[xtest][ytest] = DepthShoreface;
					
		    /*PauseRun(xtest,ytest,i);*/


		}
		else 
		    /* stop checking - ostensibly we have hit the shoreface or shore */
		{
		    ShorefaceFlag = 1;
				
		    if (_s.PercentFull[xtest][ytest] > 0)
			/* not good - somehow crossing the shore */
		    {
			/*printf("Shoreface is the Beach !!!??");*/
			/*PauseRun(xtest,ytest,-1);*/
		    }
		}
	    }
	}
    }		
		
    Depth += LandHeight;
	

    if (Depth < DepthShoreface)
    {
	printf("too deep");
	PauseRun(x,y,-1);
    }

    DeltaArea = (_s.VolumeIn[i] - _s.VolumeOut[i])/Depth;

    _s.PercentFull[_s.X[i]][_s.Y[i]] += DeltaArea/(CellWidth*CellWidth);
	
    PercentIn = _s.VolumeIn[i]/(CellWidth*CellWidth*Depth);
    PercentOut = _s.VolumeOut[i]/(CellWidth*CellWidth*Depth);
    PercentSum = DeltaArea/(CellWidth*CellWidth);

    DEBUG_PRINT( DEBUG_7A, "  In: %2.4f  Out: %2.4f  Sum: %2.4f\n", PercentIn, PercentOut, PercentSum);

}


void OopsImEmpty(int x, int y)

/*  If a cell is under-full, this will find source for desparity and move brach in	*/
/*  Function completly changed 5/21/02 sandrevt.c					*/
/*  		New Approach - steal from all neighboring AllBeach cells		*/
/*		Backup plan - steal from all neighboring percent full > 0		*/
/*  Function adjusts primary data arrays:						*/
/*		_s.AllBeach[][] and _s.PercentFull[][]					*/
	

{

    int emptycells = 0;
    int emptycells2 = 0;

    DEBUG_PRINT( DEBUG_8, "\n		OOPS I'm EMPTY!  X: %d  Y: %d Per: %f ", x, y, _s.PercentFull[x][y]);

    /* find out how many AllBeaches to take from */

    if (_s.AllBeach[x-1][y] == 'y')
	emptycells += 1;
    if (_s.AllBeach[x+1][y] == 'y')
	emptycells += 1;
    if (_s.AllBeach[x][y-1] == 'y')
	emptycells += 1;
    if (_s.AllBeach[x][y+1] == 'y')
	emptycells += 1;

    if (emptycells > 0)
    {
	/* Now Move Sediment */

	if (_s.AllBeach[x-1][y] == 'y')
	{
	    _s.PercentFull[x-1][y] += (_s.PercentFull[x][y])/emptycells;
	    _s.AllBeach[x-1][y] = 'n';
	    DEBUG_PRINT( DEBUG_8, "  MOVEDBACK");
	}		
	if (_s.AllBeach[x+1][y] == 'y')
	{
	    _s.PercentFull[x+1][y] += (_s.PercentFull[x][y])/emptycells;
	    _s.AllBeach[x+1][y] = 'n';
	    DEBUG_PRINT( DEBUG_8, "  MOVEDUP");
	}
	if (_s.AllBeach[x][y-1] == 'y')
	{
	    _s.PercentFull[x][y-1] += (_s.PercentFull[x][y])/emptycells;
	    _s.AllBeach[x][y-1] = 'n';
	    DEBUG_PRINT( DEBUG_8, "  MOVEDLEFT");
	    /*if (DEBUG_8) PauseRun(x,y,-1);*/
	}
	if (_s.AllBeach[x][y+1] == 'y')
	{
	    _s.PercentFull[x][y+1] += (_s.PercentFull[x][y])/emptycells;
	    _s.AllBeach[x][y+1] = 'n';
	    DEBUG_PRINT( DEBUG_8, "  MOVEDRIGHT");
	    /*if (DEBUG_8) PauseRun(x,y,-1);*/
	}
    }
    else
    {
	/* No full neighbors, so take away from partially full neighbors */
		
	if (_s.PercentFull[x-1][y] > 0)
	    emptycells2 += 1;
	if (_s.PercentFull[x+1][y] > 0)
	    emptycells2 += 1;
	if (_s.PercentFull[x][y-1] > 0)
	    emptycells2 += 1;
	if (_s.PercentFull[x][y+1] > 0)
	    emptycells2 += 1;

	if (emptycells2 > 0)
	{

	    if (_s.PercentFull[x-1][y] > 0)
	    {
		_s.PercentFull[x-1][y] += (_s.PercentFull[x][y])/emptycells2;
		DEBUG_PRINT( DEBUG_8, "  NOTFULL MOVEDBACK");
	    }		
	    if (_s.PercentFull[x+1][y] > 0)
	    {
		_s.PercentFull[x+1][y] += (_s.PercentFull[x][y])/emptycells2;
		DEBUG_PRINT( DEBUG_8, "  NOTFULL MOVEDUP");
	    }
	    if (_s.PercentFull[x][y-1] > 0)
	    {
		_s.PercentFull[x][y-1] += (_s.PercentFull[x][y])/emptycells2;
		DEBUG_PRINT( DEBUG_8, "  NOTFULL MOVEDLEFT");
		/*if (DEBUG_8) PauseRun(x,y,-1);*/
	    }
	    if (_s.PercentFull[x][y+1] > 0)
	    {
		_s.PercentFull[x][y+1] += (_s.PercentFull[x][y])/emptycells2;
		DEBUG_PRINT( DEBUG_8, "  NOTFULL MOVEDRIGHT");
		/*if (DEBUG_8) PauseRun(x,y,-1);*/
	    }	
	}
	else
	{
	    printf("@@@ Didn't find anywhere to steal sand from!! x: %d  y: %d\n",x,y);
	    PauseRun(x,y,-1);
	}

    }

    _s.AllBeach[x][y] = 'n';
    _s.PercentFull[x][y] = 0.0;
    _s.CellDepth[x][y] = DepthShoreface;

    DEBUG_PRINT( DEBUG_8, "\n");

}


void OopsImFull(int x, int y)

/*  If a cell is overfull, push beach out in new direction				*/
/*  Completely revised 5/20/02 sandrevt.c to resolve 0% full problems, etc.		*/
/*  New approach: 	put sand wherever 0% full in adjacent cells			*/
/*			if not 0% full, then fill all non-allbeach			*/
/*  Function adjusts primary data arrays:						*/
/*		_s.AllBeach[][] and _s.PercentFull[][]					*/
	

{

    int fillcells = 0;
    int fillcells2 = 0;

    DEBUG_PRINT( DEBUG_8, "\n		OOOPPPS I'M FULLL: X: %d  Y: %d Per: %f  ==", x, y, _s.PercentFull[x][y]);
    /*if (DEBUG_8) PrintLocalConds(x,y,-1);*/

    /* find out how many cells will be filled up	*/

    if (_s.PercentFull[x-1][y] == 0.0)
	fillcells += 1;
    if (_s.PercentFull[x+1][y] == 0.0)
	fillcells += 1;
    if (_s.PercentFull[x][y-1] == 0.0)
	fillcells += 1;
    if (_s.PercentFull[x][y+1] == 0.0)
	fillcells += 1;

    if (fillcells != 0)
    {
	/* Now Move Sediment */

	if (_s.PercentFull[x-1][y] == 0.0)
	{
	    _s.PercentFull[x-1][y] += (_s.PercentFull[x][y]-1)/fillcells;
	    _s.CellDepth[x-1][y] = - LandHeight;
	    DEBUG_PRINT( DEBUG_8, "  MOVEDBACK");
	}		
	if (_s.PercentFull[x+1][y] == 0.0)
	{
	    _s.PercentFull[x+1][y] += (_s.PercentFull[x][y]-1)/fillcells;
	    _s.CellDepth[x+1][y] = - LandHeight;
	    DEBUG_PRINT( DEBUG_8, "  MOVEDUP");
	}
	if (_s.PercentFull[x][y-1] == 0.0)
	{
	    _s.PercentFull[x][y-1] += (_s.PercentFull[x][y]-1)/fillcells;
	    _s.CellDepth[x][y-1] = - LandHeight;
	    DEBUG_PRINT( DEBUG_8, "  MOVEDLEFT");
	    /*if (DEBUG_8) PauseRun(x,y,-1);*/
	}
	if (_s.PercentFull[x][y+1] == 0.0)
	{
	    _s.PercentFull[x][y+1] += (_s.PercentFull[x][y]-1)/fillcells;
	    _s.CellDepth[x][y+1] = - LandHeight;
	    DEBUG_PRINT( DEBUG_8, "  MOVEDRIGHT");
	    /*if (DEBUG_8) PauseRun(x,y,-1);*/
	}
    }
    else
    {
	/* No fully empty neighbors, so distribute to partially full neighbors */
		
	if (_s.PercentFull[x-1][y] < 1)
	    fillcells2 += 1;
	if (_s.PercentFull[x+1][y] < 1)
	    fillcells2 += 1;
	if (_s.PercentFull[x][y-1] < 1)
	    fillcells2 += 1;
	if (_s.PercentFull[x][y+1] < 1)
	    fillcells2 += 1;

	if (fillcells2 > 0)
	{

	    if (_s.PercentFull[x-1][y] < 1)
	    {
		_s.PercentFull[x-1][y] += (_s.PercentFull[x][y]-1)/fillcells2;
		DEBUG_PRINT( DEBUG_8, "  MOVEDBACK");
	    }		
	    if (_s.PercentFull[x+1][y] < 1)
	    {
		_s.PercentFull[x+1][y] += (_s.PercentFull[x][y]-1)/fillcells2;
		DEBUG_PRINT( DEBUG_8, "  MOVEDUP");
	    }
	    if (_s.PercentFull[x][y-1] < 1)
	    {
		_s.PercentFull[x][y-1] += (_s.PercentFull[x][y]-1)/fillcells2;
		DEBUG_PRINT( DEBUG_8, "  MOVEDLEFT");
	    }
	    if (_s.PercentFull[x][y+1] < 1)
	    {
		_s.PercentFull[x][y+1] += (_s.PercentFull[x][y]-1)/fillcells2;
		DEBUG_PRINT( DEBUG_8, "  MOVEDRIGHT");
	    }	
	}
	else
	{
	    DEBUG_PRINT( DEBUG_8, "Nobody wants our sand!!! x: %d  y: %d Per: %f\n",x,y,_s.PercentFull[x][y]);
	    /*PauseRun(x,y,-1);*/
	}

    }

    _s.AllBeach[x][y] = 'y';
    _s.PercentFull[x][y] = 1.0;
    _s.CellDepth[x][y] = - LandHeight;

    DEBUG_PRINT( DEBUG_8, "\n");
	
	
}


void FixBeach(void)

/* Hopefully addresses strange problems caused by filling/emptying of cells	*/
/* Looks at entire data set							*/
/* Find unattached pieces of sand and moves them back to the shore 		*/
/* Takes care of 'floating bits' of sand					*/
/* Also takes care of over/under filled beach pieces				*/
/* Revised 5/21/02 to move sand to all adjacent neighbors sandrevt.c 		*/
/* Changes global variable _s.PercentFull[][]					*/
/* Uses but does not change _s.AllBeach[][]					*/
/* sandrevx.c - added sweepsign to reduce chances of asymmetrical artifacts	*/


{

    int i,x,y,sweepsign;
    int FixXMax;
    int fillcells3 = 0;

    /*DEBUG_PRINT( DEBUG_9, "\n\nFIXBEACH      %d     %f\n", _s.CurrentTimeStep, WaveAngle*radtodeg);*/

    if (RandZeroToOne()*2 > 1)
    {
	sweepsign = 1;
	DEBUG_PRINT( DEBUG_9, "fixL  ");
    }
    else
    {
	sweepsign = 0;
	DEBUG_PRINT( DEBUG_9, "fixR  ");
    }


    FixXMax = _s.ShadowXMax + ceil(DepthShoreface/CellWidth/ShorefaceSlope) +3;
    if (FixXMax > Xmax)
	FixXMax = Xmax; 

    for (x = FixXMax; x >= 0 ; x--)
    {
	for (i = 0; i <= 2*Ymax ; i++)
	{

	    if (sweepsign == 1)
		y = i;
	    else
		y = 2*Ymax-i;

	    /* ye olde depth fix */
	    if ((_s.PercentFull[x][y] <= 0) && (_s.CellDepth[x][y] > DepthShoreface) && 
		(_s.CellDepth[x-1][y] == DepthShoreface))
	    {
		if ((_s.CellDepth[x+1][y] == DepthShoreface) && (_s.CellDepth[x][y-1] == DepthShoreface)
		    && (_s.CellDepth[x][y+1] == DepthShoreface))
		{
		    /* Fill Hole */
		    _s.CellDepth[x][y] = DepthShoreface;
		}
	    }
	    if (_s.PercentFull[x][y]> 100)
	    {
		printf("too full");
		_s.PercentFull[x][y] = 0;
		PauseRun(x,y,-1);
	    }
	


	    /* Take care of situations that shouldn't exist */

	    if (_s.PercentFull[x][y] < 0)
	    {
		_s.AllBeach[x][y] = 'n';
		DEBUG_PRINT( DEBUG_9 && y != 0, "\nUnder 0 Percent X: %d  Y: %d Percent: %f\n", x,y,_s.PercentFull[x][y]);
		OopsImEmpty(x,y);
		printf("Underzerofill");
		/*PauseRun(x,y,-1);*/
	    }

	    if (_s.PercentFull[x][y] > 1)
	    {
		_s.AllBeach[x][y] = 'y';
		_s.CellDepth[x][y] = - LandHeight;
		DEBUG_PRINT( DEBUG_9 && y != 0, "\nOver 100 Percent X: %d  Y: %d Per: %f\n"
					     ,x,y, _s.PercentFull[x][y]);
		OopsImFull(x,y);
	    }

	    if (((_s.PercentFull[x][y] >=0) && (_s.PercentFull[x][y] <1)) && (_s.AllBeach[x][y] == 'y'))
	    {
		_s.AllBeach[x][y] = 'n';
		_s.CellDepth[x][y] = - LandHeight;
		DEBUG_PRINT( DEBUG_9 && y != 0, "\nALLBeachProb X: %d  Y: %d\n", x,y);
	    }


	    /* Take care of 'loose' bits of sand */

	    fillcells3 = 0;
			
	    if ( (_s.PercentFull[x][y] != 0) && (_s.PercentFull[x-1][y] < 1) && (_s.PercentFull[x+1][y] < 1) &&
		 (_s.PercentFull[x][y+1] < 1) && (_s.PercentFull[x][y-1] < 1) && (_s.AllBeach[x][y] =='n'))
		/* Beach in cell, but bottom, top, right, and left neighbors not all full */
	    {
		DEBUG_PRINT( DEBUG_9 && y != 0, "\nFB Moved loose bit of sand,  X: %d  Y: %d  Per: %f  ",
					     x, y, _s.PercentFull[x][y]);

		/* distribute to partially full neighbors */
		
		if ((_s.PercentFull[x-1][y] < 1) && (_s.PercentFull[x-1][y] > 0))
		    fillcells3 += 1;
		if ((_s.PercentFull[x+1][y] < 1) && (_s.PercentFull[x+1][y] > 0))
		    fillcells3 += 1;
		if ((_s.PercentFull[x][y-1] < 1) && (_s.PercentFull[x][y-1] > 0))
		    fillcells3 += 1;
		if ((_s.PercentFull[x][y+1] < 1) && (_s.PercentFull[x][y+1] > 0))
		    fillcells3 += 1;				

		if ((fillcells3 > 0))
		{			

		    if ((_s.PercentFull[x-1][y] < 1) && (_s.PercentFull[x-1][y] > 0))
		    {
			_s.PercentFull[x-1][y] += (_s.PercentFull[x][y])/fillcells3;
			DEBUG_PRINT( DEBUG_9, "  MOVEDBACK");
		    }		
		    if ((_s.PercentFull[x+1][y] < 1) && (_s.PercentFull[x+1][y] > 0))
		    {
			_s.PercentFull[x+1][y] += (_s.PercentFull[x][y])/fillcells3;
			DEBUG_PRINT( DEBUG_9, "  MOVEDUP");
		    }
		    if ((_s.PercentFull[x][y-1] < 1) && (_s.PercentFull[x][y-1] > 0))
		    {
			_s.PercentFull[x][y-1] += (_s.PercentFull[x][y])/fillcells3;
			DEBUG_PRINT( DEBUG_9, "  MOVEDLEFT");
			/*if (DEBUG_9) PauseRun(x,y,-1);*/
		    }
		    if ((_s.PercentFull[x][y+1] < 1) && (_s.PercentFull[x][y+1] > 0))
		    {
			_s.PercentFull[x][y+1] += (_s.PercentFull[x][y])/fillcells3;
			DEBUG_PRINT( DEBUG_9, "  MOVEDRIGHT");
			/*if (DEBUG_9) PauseRun(x,y,-1);*/
		    }	
		}
		else
		{
		    printf("Loner fixbeach breakdown - mass disintegrated x: %d  y: %d\n",x,y);
		    if (DEBUG_9) 
			PauseRun(x,y,-1);
		}	

		_s.PercentFull[x][y] = 0;
		_s.AllBeach[x][y] = 'n';
		_s.CellDepth[x][y] = DepthShoreface;

		DEBUG_PRINT( DEBUG_9, "\n");


		/* If we have overfilled any of the cells in this loop, need to OopsImFull() */

		if (_s.PercentFull[x-1][y] > 1)
		{
		    OopsImFull(x-1,y);
		    DEBUG_PRINT( DEBUG_9, "	Below Overfilled\n");
		}
		if (_s.PercentFull[x][y-1] > 1)
		{
		    OopsImFull(x,y-1);
		    DEBUG_PRINT( DEBUG_9, "	Left Side Overfilled\n");	
		}
		if (_s.PercentFull[x][y+1] > 1)
		{
		    OopsImFull(x,y+1);
		    DEBUG_PRINT( DEBUG_9, "	Right Side Overfilled\n");
		}
		if (_s.PercentFull[x+1][y+1] > 1)
		{
		    OopsImFull(x+1,y+1);
		    DEBUG_PRINT( DEBUG_9, "	Top Overfilled\n");
		}
			
	    }
		
	    /*if ((_s.AllBeach[x][y] =='y') && (_s.PercentFull[x-1][y] < 1) && (_s.PercentFull[x+1][y] < 1)
	      && (_s.PercentFull[x][y-1] < 1) && (_s.PercentFull[x][y+1] < 1)
	      && (_s.AllBeach[x-1][y-1] == 'n') && (_s.AllBeach[x-1][y+1] == 'n') &&
	      (_s.AllBeach[x+1][y+1] == 'n') && (_s.AllBeach[x+1][y-1] == 'n') )

	      {
	      printf("%% Booger !! x: %d  y: %d", x,y);
	      PauseRun(x,y,-1);
	      }*/

	}
    }
	

}


float MassCount(void)

/* Counts the total volume occupied by beach cells 	*/
/* Uses same algorhythm as AdjustShore			*/
/* returns a float of the total sum 			*/
/* Uses _s.AllBeach[][] and _s.PercentFull[][]		*/
/* and InitialDepth, CellWidth, ShelfSlope		*/

{

    int 	x,y;
    float 	Mass = 0;
    /*float  	MassHere;
      float 	refdepth;
			
      refdepth = InitialDepth;*/

    for (x=0; x < Xmax ; x++)
    {
	for(y=Ymax/2; y < 3 * Ymax /2; y++)
	{
	    /*if ((_s.PercentFull[x][y] > 0) && (_s.PercentFull[x][y] < 1.0))
	      MassHere = _s.PercentFull[x][y] * (refdepth - _s.CellDepth[x][y]) + 
	      (1 - _s.PercentFull[x][y])*(refdepth - DepthShoreface);
	      else 
	      MassHere = refdepth - _s.CellDepth[x][y];*/

	    Mass += _s.PercentFull[x][y];
	}
    }

    return Mass;

}

	
float Raise(float b, float e)		

/* function calulates b to the e power */
/* pow has problems if b <= 0 */
	
{
    if (b>0)
	return powf(b,e);
    else 
	return -powf(fabs(b),e);
}


float RandZeroToOne(void)

/* function will return a random number equally distributed between zero and one */
/* currently this function has no seed */

{
    return random()/(Raise(2,31)-1);
}


void InitConds(void)

/* Creates initial beach conditions 						*/
/* Flat beach with zone of AllBeach = 'y' separated by AllBeach = 'n' 		*/
/* Bounding layer set to random fraction of fullness 				*/
	
{
    int 	x,y;
    printf("Condition Initial \n");

    if (InitCType == 0)
	/* 'Regular Initial cons - beach backed by sandy land */
    {
	
	for (y = 0; y <= 2*Ymax; y++)
	    for (x = 0; x <= Xmax; x++)
	    {

		_s.CellDepth[x][y] = InitialDepth + ((x-InitBeach) * CellWidth * ShelfSlope);
			
		if (x < InitBeach)
		{  	 
		    _s.PercentFull[x][y] = 1;
		    _s.AllBeach[x][y] = 'y';
		    _s.CellDepth[x][y] = - LandHeight;
		}
		else if (x == InitBeach)
		{
		    if (InitialSmooth)
		    {
			_s.PercentFull[x][y] = .5;
		    }
		    else
		    {
			_s.PercentFull[x][y] = RandZeroToOne();
			printf("x: %d  Y: %d  Per: %f\n",x,y,_s.PercentFull[x][y]);
		    }
		    _s.AllBeach[x][y] = 'n';
		    _s.CellDepth[x][y] = - LandHeight;
		}
		else if (x > InitBeach)
		{
		    _s.PercentFull[x][y] = 0;
		    _s.AllBeach[x][y] = 'n';
		    if (_s.CellDepth[x][y] < DepthShoreface)
		    {
			_s.CellDepth[x][y] = DepthShoreface;
		    }
		}
		else
		{
		printf("WTF! x: %d  Y: %d  Per: %f\n",x,y,_s.PercentFull[x][y]);
		PauseRun(x,y,-1);
		}

		_s.Age[x][y] = 0;		
	    }
    }

    else if (InitCType == 1)
	/* 'Simple Barrier' type initial condition - island backed by lagoon at slope of shelf */
    {
	
	for (y = 0; y <= 2*Ymax; y++)
	{
	    for (x = 0; x < Xmax; x++)
	    {

		_s.CellDepth[x][y] = InitialDepth + ((x-InitBeach) * CellWidth * ShelfSlope);
			
		if (_s.CellDepth[x][y] <= 0)
		    /* This must be land due to continental shelf intersection */
		{
		    _s.PercentFull[x][y] = 1.0;
		    _s.AllBeach[x][y] = 'y';
		    _s.CellDepth[x][y] = - LandHeight;
		}
		else if (x > InitBeach)
		    /* Shoreward of beach - enforce ShorefaceDepth if necessary */
		{
		    _s.PercentFull[x][y] = 0;
		    _s.AllBeach[x][y] = 'n';
		    if (_s.CellDepth[x][y] < DepthShoreface)
		    {
			_s.CellDepth[x][y] = DepthShoreface;
		    }
		}
		else if (x == InitBeach)
		    /* Beach */
		{
		    if (InitialSmooth)
		    {
			_s.PercentFull[x][y] = .5;
		    }
		    else
		    {
			_s.PercentFull[x][y] = RandZeroToOne();
			/*printf("x: %d  Y: %d  Per: %f\n",x,y,_s.PercentFull[x][y]);*/
		    }
		    _s.AllBeach[x][y] = 'n';
		    _s.CellDepth[x][y] = - LandHeight;
		}
		else if ((x < InitBeach) && (x > InitBeach - InitBWidth - 1))
		    /* Island */
		{  	 
		    _s.PercentFull[x][y] = 1.0;
		    _s.AllBeach[x][y] = 'y';
		    _s.CellDepth[x][y] = - LandHeight;
		}
		else if (x == InitBeach - InitBWidth -1)
		    /* Back of Barrier */
		{  	 
		    if (InitialSmooth)
		    {
			_s.PercentFull[x][y] = .5;
		    }
		    else
		    {
			_s.PercentFull[x][y] = RandZeroToOne();
			printf("x: %d  Y: %d  Per: %f\n",x,y,_s.PercentFull[x][y]);
		    }
		    _s.AllBeach[x][y] = 'n';
		    _s.CellDepth[x][y] = - LandHeight;
		}
		else if (x < InitBeach - InitBWidth -1)
		    /* Lagoon at depth of shelf slope  */
		{
		    _s.PercentFull[x][y] = 0;
		    _s.AllBeach[x][y] = 'n';
		}
		if (_s.PercentFull[x][y] > 1)
		{
		    printf("x: %d  Y: %d  Per: %f\n",x,y,_s.PercentFull[x][y]);
		    PauseRun(x,y,-1);
		}
		_s.Age[x][y] = 0;		
	    }
	}
    }

}
	

void InitPert(void)
	
/* Andrew's initial bump */

	
{
	
    int x,y;	
    int PWidth = 40;
    int PHeight = 40;
    int PYstart = 100;

    if (InitialPert == 1)
	/* Square perturbation */	
    {

	/* Fill AllBeach areas */

	for (x = InitBeach ; x <= InitBeach + PHeight ; x++)
	{	
	    for (y = PYstart ; y <= PYstart + PWidth ; y++)
	    {
		_s.PercentFull[x][y] = 1.0;
		_s.AllBeach[x][y] = 'y';
	    } 	
	}


	/* PercentFull Top */

	for (y = PYstart -1; y <= PYstart + PWidth +1; y++)
	{
	    _s.PercentFull[InitBeach + PHeight + 1][y] = RandZeroToOne();
	}

	/* PercentFull Sides */

	for (x = InitBeach ; x <= InitBeach + PHeight ; x++)
	{
	    _s.PercentFull[x][PYstart-1] = RandZeroToOne();
	    _s.PercentFull[x][PYstart+PWidth + 1] = RandZeroToOne();
	}
    }

    else if (InitialPert == 2)
	/* Another Perturbation  - steep point */
    {
	
	x = InitBeach;

	_s.PercentFull[x][17] = 0.8;
	_s.PercentFull[x][18] = 1.0;
	_s.AllBeach[x][18] = 'y';
	_s.PercentFull[x][19] = 0.8;

	x = InitBeach + 1;

	_s.PercentFull[x][17] = 0.6;
	_s.PercentFull[x][18] = 1.0;
	_s.AllBeach[x][18] = 'y';
	_s.PercentFull[x][19] = 0.6;

	x = InitBeach + 2;

	_s.PercentFull[x][17] = 0.2;
	_s.PercentFull[x][18] = 1.0;
	_s.AllBeach[x][18] = 'y';
	_s.PercentFull[x][19] = 0.2;

	x = InitBeach + 3;

	_s.PercentFull[x][18] = 0.3;

    }

}


void PeriodicBoundaryCopy(void)	

/* Simulates periodic boundary conditions by copying middle section to front and end of arrays */


{
    int	x,y;

    for (y = Ymax; y < 3*Ymax/2; y++)
	for (x = 0; x < Xmax; x++)
	{
	    _s.AllBeach[x][y-Ymax] = _s.AllBeach[x][y];
	    _s.PercentFull[x][y-Ymax] = _s.PercentFull[x][y];
	    _s.Age[x][y-Ymax] = _s.Age[x][y];
	    _s.CellDepth[x][y-Ymax] = _s.CellDepth[x][y];
	}
    for (y = Ymax/2; y <= Ymax; y++)
	for (x = 0; x < Xmax; x++)
	{
	    _s.AllBeach[x][y+Ymax] = _s.AllBeach[x][y];
	    _s.PercentFull[x][y+Ymax] = _s.PercentFull[x][y];
	    _s.Age[x][y+Ymax] = _s.Age[x][y];
	    _s.CellDepth[x][y+Ymax] = _s.CellDepth[x][y];
	}
	
}










void ZeroVars(void)

/* Resets all arrays recalculated at each time step to 'zero' conditions */

{

    int z;

    for (z=0; z < MaxBeachLength; z++)
    {
	_s.X[z] = -1;	
	_s.Y[z] = -1;		
	_s.InShadow[z] = '?';	
	_s.ShorelineAngle[z] = -999;
	_s.SurroundingAngle[z] = -998;
	_s.UpWind[z] = '?';	
	_s.VolumeIn[z] = 0;	
	_s.VolumeOut[z] = 0;	
    }
}

			
void ReadSandFromFile(void)	
	
/*  Reads saved output file, _s.AllBeach[][] & _s.PercentFull[][]	 */

{
    int x,y;
    FILE *ReadSandFile;
	
    ReadSandFile = fopen(_io.readfilename,"r");printf("CHECK READ \n");

	
    for (y = Ymax/2; y < 3*Ymax/2; y++)
    {	

	for (x=0; x<Xmax; x++)
	{
	    fscanf(ReadSandFile, " %f", &_s.PercentFull[x][y]);

	    if (_s.PercentFull[x][y] >= 1.0)
		_s.AllBeach[x][y] = 'y';
	    else 
		_s.AllBeach[x][y] = 'n';
	}
    }

    for (y = Ymax/2; y < 3*Ymax/2; y++)
	for (x=0; x<Xmax; x++)
	    fscanf(ReadSandFile, " %f", &_s.CellDepth[x][y]);

				
    if (SaveAge)
	
	for (y = Ymax/2; y < 3*Ymax/2; y++)
	{
	    for (x=0; x<Xmax; x++)
	    {
		fscanf(ReadSandFile, " %d", &_s.Age[x][y]);
	    }
	}

    /*PrintLocalConds(5,5,-1);*/
    fclose(ReadSandFile);
    printf("file read!");

    PeriodicBoundaryCopy();

}


void SaveSandToFile(void)

/*  Saves current _s.AllBeach[][] and _s.PercentFull[][] data arrays to file 		*/
/*  Save file name will add extension '.' and the _s.CurrentTimeStep		*/

{
    int	x,y;
    char savename[40];
    FILE *SaveSandFile;

    printf("\n saving \n ");

    sprintf(savename, "%s.%d", _io.savefilename, _s.CurrentTimeStep);
    printf( "Saving as: %s 		", savename );	


    SaveSandFile = fopen(savename, "w");
    if (!SaveSandFile)
    {
	printf("problem opening output file\n");
	exit(1);
    }
	
	
    for (y= Ymax/2; y< 3*Ymax/2; y++)
	for (x=0; x<Xmax; x++)
	    fprintf(SaveSandFile, " %f", _s.PercentFull[x][y]);

    for (y= Ymax/2; y< 3*Ymax/2; y++)
	for (x=0; x<Xmax; x++)
	    fprintf(SaveSandFile, " %f", _s.CellDepth[x][y]);

    if (SaveAge)
	for (y=Ymax/2; y< 3*Ymax/2; y++)
	    for (x=0; x<Xmax; x++)
		fprintf(SaveSandFile, " %d", _s.Age[x][y]);
		

    fclose(SaveSandFile);
    printf("--- regular file saved! ----\n\n");

}


#define SAVE_LINE_NAME "lineout"
///* SaveLineToFile */char            savelinename[24] = "lineout";

void SaveLineToFile(void)

/*  Saves data line of shoreline position rather than entire array 		*/
/*  Main concern is to have only one data point at each alongshore location	*/
/*  Save file name will add extension '.' and the _s.CurrentTimeStep		*/

{

    int	y,x,xtop,i;
    float   xsave;
    char savename[40];
    char savelinename[24] = SAVE_LINE_NAME;
    FILE *SaveSandFile;

    printf("\n saving \n ");

    sprintf(savename, "%s%d", savelinename, _s.CurrentTimeStep);
    printf( "Saving as: %s                 ", savename );	

    SaveSandFile = fopen(savename, "w");
    if (!SaveSandFile)
    {
	printf("problem opening output file\n");
	exit(1);
    }
	

    for (y=Ymax/2; y < 3*Ymax/2; y++)
    {

	x = Xmax-1;	
	xtop = Xmax;	

	/* step back to where we encounter allbeach */
	while(_s.AllBeach[x][y] == 'n')
	{
	    x -= 1;
	}

	/* if on side of shape, need to average */
	if (_s.PercentFull[x+2][y] > 0)
	{
	    xtop = x+1;
	    while(_s.PercentFull[xtop][y] > 0)
	    {
		xtop +=1;
	    } 

	    xsave = x;

	    for (i=x+1; i<xtop ; i++)
		xsave += _s.PercentFull[i][y];

	}
	/* otherwise Regular Beach Condition */
	else
	{
	    xsave = x + _s.PercentFull[x+1][y];
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
	






void PrintLocalConds(int x, int y, int in)
	
/* Prints Local Array Conditions aound x,y */

{ 

    int i,j,k,isee;
    float vol = CellWidth*CellWidth*DepthShoreface;

    printf("\n x: %d  y: %d  z: %d\n\n", x,y,in);

    /* if not given location along beach, look to see if along beach */

    if (in<0)
    {
	for (i=0; i <= _s.TotalBeachCells; i++)
	    if ((_s.X[i]==x) && (_s.Y[i]==y))
		isee = i;
    }
    else
	isee = in;


    for (i= x+2 ; i > x-3 ; i--)
    {
	for (j = y-2 ; j < y+3 ; j++)
	{
	    printf("  	%d,%d", i , j);
	}
	printf("\n");
    }
    printf("\n");

    for (i= x+2 ; i > x-3 ; i--)
    {
	for (j = y-2 ; j < y+3 ; j++)
	{
	    printf("  	%f", _s.CellDepth[i][j]);
	    if (_s.CellDepth[i][j] == DepthShoreface)
		printf("y");
	    else
		printf("n");
	}
	printf("\n");
    }
    printf("\n");

    for (i= x+2 ; i > x-3 ; i--)
    {
	for (j = y-2 ; j < y+3 ; j++)
	{
	    printf("	%c", _s.AllBeach[i][j]);
	}
	printf("\n");
    }

    printf("\n");

    for (i= x+2 ; i > x-3 ; i--)
    {
	for (j = y-2 ; j < y+3 ; j++)
	{
	    printf("	%2.5f",_s.PercentFull[i][j]);
	}
	printf("\n");
    }		

    printf("\n");

    printf(" %d  ", in );

    if (isee>=0)
    {
	for (k = in-3; k <in+4; k++)
	{
	    printf("  	%2d: %2d,%2d", k , _s.X[k] , _s.Y[k]);
	}
	printf("\n\n\n");

	printf("Wave Angle:	%f\n\n",WaveAngle*radtodeg);
	printf("i		%d		%d		!%d		%d		%d\n", 
	       in-2, in-1,in,in+1,in+2);
	printf("Shadow		%c		%c		%c		%c		%c\n", 
	       _s.InShadow[in-2], _s.InShadow[in-1],_s.InShadow[in], _s.InShadow[in+1], _s.InShadow[in+2]);
	printf("Upwind		%c		%c		%c		%c		%c\n", 
	       _s.UpWind[in-2],  _s.UpWind[in-1],_s.UpWind[in],  _s.UpWind[in+1],  _s.UpWind[in+2]);
	printf("Angle		%2.2f		%2.2f		%2.2f		%2.2f		%2.2f\n",
	       _s.ShorelineAngle[in-2]*radtodeg,_s.ShorelineAngle[in-1]*radtodeg, _s.ShorelineAngle[in]*radtodeg,
	       _s.ShorelineAngle[in+1]*radtodeg,_s.ShorelineAngle[in+2]*radtodeg);
	printf("SurrAngle	%2.2f		%2.2f		%2.2f		%2.2f		%2.2f\n",
	       _s.SurroundingAngle[in-2]*radtodeg, _s.SurroundingAngle[in-1]*radtodeg,  _s.SurroundingAngle[in]*radtodeg,
	       _s.SurroundingAngle[in+1]*radtodeg, _s.SurroundingAngle[in+2]*radtodeg);
	printf("Vol In 		%2.2f		%2.2f		%2.2f		%2.2f		%2.2f\n",
	       _s.VolumeIn[in-2], _s.VolumeIn[in-1],_s.VolumeIn[in],_s.VolumeIn[in+1],_s.VolumeIn[in+2]);
	printf("Vol Out		%2.2f		%2.2f		%2.2f		%2.2f		%2.2f\n",
	       _s.VolumeOut[in-2], _s.VolumeOut[in-1], _s.VolumeOut[in],_s.VolumeOut[in+1],_s.VolumeOut[in+2]);
	printf("Diff		%2.2f		%2.2f		%2.2f		%2.2f		%2.2f\n",
	       _s.VolumeIn[in-2]-_s.VolumeOut[in-2], _s.VolumeIn[in-1]-_s.VolumeOut[in-1], _s.VolumeIn[in]-_s.VolumeOut[in],
	       _s.VolumeIn[in+1]-_s.VolumeOut[in+1],_s.VolumeIn[in+2]-_s.VolumeOut[in+2]);
	printf("Frac Diff	%2.3f		%2.3f		%2.3f		%2.3f		%2.3f\n",
	       (_s.VolumeIn[in-2]-_s.VolumeOut[in-2])/vol, (_s.VolumeIn[in-1]-_s.VolumeOut[in-1])/vol,
	       (_s.VolumeIn[in]-_s.VolumeOut[in])/vol, (_s.VolumeIn[in+1]-_s.VolumeOut[in+1])/vol,
	       (_s.VolumeIn[in+2]-_s.VolumeOut[in+2])/vol);

    }

    printf("\n");





}


void PauseRun(int x, int y, int in)

/* Pauses run intil the 'q' key is pressed 	*/
/* Can Print or Plot Out Useful info		*/
	

{
		
    int xsee=1,ysee=-1,isee=-1,i;

    printf("\nPaused x: %d  y: %d Time: %d\n",x,y,_s.CurrentTimeStep);

    /*if (SaveLine) SaveLineToFile();
      else SaveSandToFile();*/

    if(NoPauseRun)
	return;

    sleep(5);
    printf("\nend Pause\n");
	
}
 

void ButtonEnter(void)

{

    char newdigit = 'z';
    int flag = 0;	
    int i = 1;
    char digits[7] = "";

    /*printf("Flag1 %d\n",flag);*/	

    printf("Press <Space> to Start\n");

	

    printf("Enter Digit and <Space> (<Z> to finish)\n");

	
	 
    i += 1;
    sprintf(digits,"%s%c",digits,newdigit);
    printf("\nCurrent %s\n",digits);
    newdigit = 'z';
			


	

}


#define AGE_MAX (1000000) /* Maximum 'age' of cells - loops back to zero */
//int 	AgeMax = 1000000;	/* Maximum 'age' of cells - loops back to zero */

void AgeCells(void)

/* Age Cells */

{
	
    int x,y;
    int	AgeMax = AGE_MAX;

    for (y = 0; y < 2*Ymax; y++)
	for (x=0; x<Xmax; x++)
	    if (_s.PercentFull[x][y] == 0)
	    {
		_s.Age[x][y] = _s.CurrentTimeStep%AgeMax;
	    }
				
}

#define READ_WAVE_NAME "WIS_509_150.dat"
///* ReadWaveIn */char            readwavename[24] = "WIS_509_150.dat";

void ReadWaveIn(void)

/* Input Wave Distribution */

{
    int i;
    char readwavename[24] = READ_WAVE_NAME;
    FILE *ReadWaveFile;

    for (i=0 ; i<= 36; i++)
    {
	WaveMax[i] =0;
	WaveProb[i] = 0;
    }

    ReadWaveFile = fopen(readwavename,"r");printf("CHECK READ WAVE\n");

    fscanf(ReadWaveFile, " %d \n", &NumWaveBins);
	
    printf("Wave Bins %d \n",NumWaveBins);

    WaveMax[0] = -90;
    WaveProb[0] = 0;

    for (i=1 ; i<= NumWaveBins ; i++)
    {
	fscanf(ReadWaveFile, " %f %f", &WaveMax[i] , &WaveProb[i]);
	printf("i= %d  Wave= %f Prob= %f \n",i, WaveMax[i], WaveProb[i]);
    }	

    fclose(ReadWaveFile);
    printf("wave file read! \n");
				
}


#ifdef WITH_OPENGL
Bool WaitForNotify(Display *d, XEvent *e, char *arg)
{
   return (e->type == MapNotify) && (e->xmap.window == (Window)arg);
}

void OpenWindow(void)

{

    static  int  attributeListSgl[]  =  {GLX_RGBA, GLX_RED_SIZE, 1, GLX_GREEN_SIZE, 1, GLX_BLUE_SIZE, 1, None };
    static   int   attributeListDbl[]   =   { GLX_RGBA, 1, GLX_RED_SIZE, 1, GLX_GREEN_SIZE, 1, GLX_BLUE_SIZE, 1, None };
                                                                                                                      
                                                                                                                      
//    static Bool WaitForNotify(Display *d, XEvent *e, char *arg) {
//	return (e->type == MapNotify) && (e->xmap.window == (Window)arg); }

    /*int FALSE =0;
      int TRUE =1;*/
    int winwidth;
    int winheight;
	   
    Display *dpy;
    XVisualInfo *vi;
    Colormap cmap;
    XSetWindowAttributes swa;
    Window win;
    GLXContext cx;
    XEvent event;
    int swap_flag = FALSE;                                                                                                         
                                                                                                                      
    dpy = XOpenDisplay(0);                                                                                            
                                                                                                                      
    vi = glXChooseVisual(dpy, DefaultScreen(dpy), attributeListSgl); 	

if (vi == NULL) {
	vi = glXChooseVisual(dpy, DefaultScreen(dpy), attributeListDbl);
	swap_flag = TRUE;
    }
    cx = glXCreateContext(dpy, vi, 0, GL_TRUE);
                                                                                                                      
                                                                                                                      
    cmap = XCreateColormap(dpy, RootWindow(dpy, vi->screen),  vi->visual, AllocNone);
                                                                                                                      
                                                                                                                      
    swa.colormap = cmap;
    swa.border_pixel = 0;
    swa.event_mask = StructureNotifyMask;

    winwidth = CELL_PIXEL_SIZE*Xmax;
    winheight = CELL_PIXEL_SIZE*Ymax;
    

    win = XCreateWindow(dpy, RootWindow(dpy, vi->screen),0, 0, YPlotExtent * CELL_PIXEL_SIZE, XPlotExtent * CELL_PIXEL_SIZE,
			0, vi->depth, InputOutput, vi->visual,
			CWBorderPixel|CWColormap|CWEventMask, &swa);
    XMapWindow(dpy, win);
    XIfEvent(dpy, &event, WaitForNotify, (char*)win);
                                                                                                    
    glXMakeCurrent(dpy, win, cx);                                                    
    glClearColor(0.0,0.0,0.0,1);
    glClear(GL_COLOR_BUFFER_BIT);
    glFlush();
}


void PutPixel(float x, float y, float R, float G, float B)

{

    float xstart, ystart;

    /* translate x and y integer components to the openGL grid */

    xstart = (x - Xmax/2.0)/(Xmax/2.0);
    ystart = (y - Ymax/2.0)/(Ymax/2.0);

    glColor3f (R, G, B);
    glBegin(GL_POLYGON);
    glVertex3f (ystart, xstart, 0.0);
    glVertex3f (ystart, xstart+_g.xcellwidth, 0.0);
    glVertex3f (ystart + _g.ycellwidth, xstart + _g.xcellwidth, 0.0);
    glVertex3f (ystart + _g.ycellwidth, xstart, 0.0);
    glEnd();

}



#define AGE_SHADE_SPACING (10000) /* For graphics - how many time steps means back to original shade */
//int 	AgeShadeSpacing = 10000; /* For graphics - how many time steps means back to original shade */
void GraphCells(void)

/* Plots entire Array */ 

{ 
    int x,y;
    float Red,Green,Blue,backRed,backGreen,backBlue;
    float AgeFactorRed,AgeFactorGreen,AgeFactorBlue,DepthBlue;
    float DepthFactorX;
    int AgeShadeSpacing = AGE_SHADE_SPACING; /* For graphics - how many time steps means back to original shade */
	
    DepthFactorX = XPlotExtent;
    for (y=_g.yplotoff; y <= YPlotExtent+_g.yplotoff; y++)

	for (x= _g.xplotoff; x <= XPlotExtent+_g.xplotoff; x++)
	{

	    backRed = 0;
	    backBlue = (75 + 130 * (x/DepthFactorX));
	    backGreen = (165 - 125 * (x/DepthFactorX));

	    DepthBlue =  floor( (.8 - (float)x/(XPlotExtent*3)) * 255)/255.0; 
	    /*	PutPixel( CELL_PIXEL_SIZE*(y-yplotoff), CELL_PIXEL_SIZE*x,0,0, DepthBlue);*/
 
	    AgeFactorRed = (float)((_s.Age[x][y])%AgeShadeSpacing)/AgeShadeSpacing;
	    AgeFactorGreen = (float)((_s.Age[x][y]+AgeShadeSpacing/3)%AgeShadeSpacing)/AgeShadeSpacing;
	    AgeFactorBlue = (float)((_s.Age[x][y]+2*AgeShadeSpacing/3)%AgeShadeSpacing)/AgeShadeSpacing;

	    if ((_s.PercentFull[x][y] > 0) && (_s.AllBeach[x][y] == 'n'))
	    {
		Red =((((235 - 100 *(AgeFactorRed))-backRed)* _s.PercentFull[x][y] )+backRed)/255.0;
		Green =((((235 - 95 * (AgeFactorGreen))-backGreen) * _s.PercentFull[x][y])+backGreen)/255.0;
		Blue =((((210 - 150 * AgeFactorBlue)-backBlue)* _s.PercentFull[x][y])+backBlue)/255.0;

	    }
	    else if (_s.AllBeach[x][y] == 'y')
	    {
		Red =((((235 - 100 *(AgeFactorRed))-backRed)* _s.PercentFull[x][y] )+backRed)/255.0;
		Green =((((235 - 95 * (AgeFactorGreen))-backGreen) * _s.PercentFull[x][y])+backGreen)/255.0;
		Blue =((((210 - 150 * AgeFactorBlue)-backBlue)* _s.PercentFull[x][y])+backBlue)/255.0;
	    }
	    else
	    {
		Red = backRed/255.0;
		Blue = backBlue/255.0;
		Green = backGreen/255.0;
	    }

	    PutPixel(x-_g.xplotoff,y-_g.yplotoff,Red,Green,Blue);

	}

    x = 0;

    y = StreamSpot;

    while (_s.AllBeach[x][y] == 'y')
    {
	PutPixel(x-_g.xplotoff, y-_g.yplotoff, 1,0,0);
	x += 1;
    }

    glFlush();

}


void ScreenInit(void) 
/* this is for the keyboard thingies to work */

{
    WINDOW *mainwnd;
    WINDOW *screen;

    mainwnd = initscr();
/*noecho();*/
    cbreak();
    nodelay(mainwnd, TRUE);
    halfdelay(.1);
    refresh(); // 1)
    wrefresh(mainwnd);
    screen = newwin(13, 27, 1, 1);
    box(screen, ACS_VLINE, ACS_HLINE);
}
#endif

void	DeliverSediment(void)

/* Simple 'first approximation of sediment delivery	*/
/* At certain alongshore location, add certain amount of sed to the coast */

	{

	int x,y;
	
	x = 0;
	y = StreamSpot;

	while (_s.AllBeach[x][y] == 'y')
	{
		x += 1;
	}

	_s.PercentFull[x][y] += SedRate;

	}

#define OVERWASH_LIMIT (75) /* beyond what angle don't do overwash */
//float	OverwashLimit = 75;	/* beyond what angle don't do overwash */

void CheckOverwashSweep(void)

	/* Just a loop to call overwash check founction CheckOverwash				 	*/
	/* Nothing done here, but can be down when CheckOVerwash is called				*/


	{
   float OverwashLimit = OVERWASH_LIMIT;

	int i,ii;		/* local loop variable */	
	int sweepsign;

	
	if (RandZeroToOne()*2 > 1)
	{
		sweepsign = 1;
		DEBUG_PRINT( DEBUG_10A, "L  ");
	}
	else
	{
		sweepsign = 0;
		DEBUG_PRINT( DEBUG_10A, "R  ");
	}

	OWflag = 0;
	for (i=1; i < _s.TotalBeachCells-1 ; i++)
	{
		if (sweepsign == 1)
			ii = i;
		else
			ii = _s.TotalBeachCells-1-i;

		/* To do test shoreline should be facing seaward 					*/
		/* don't worry about shadow here, as overwash is not set to a time scale with AST 	*/

		if ((fabs(_s.SurroundingAngle[ii]) < (OverwashLimit/radtodeg))  && (_s.InShadow[ii] == 'n'))
		{
			CheckOverwash(ii);
		}	

	}

	
	/*if (OWflag) PauseRun(1,1,-1);*/

}


	

void CheckOverwash(int icheck)

	/* New 1/04 ADA - Step back pixelwise in direction of Surrounding Angle to check needage 	*/
	/* If too short, calls DoOverwash, which will move some sediment				*/
	/* Uses _s.AllBeach[][] and _s.PercentFull[][] (can be changed when DoOVerwash is called		*/
	/* Need to change sweepsign because filling cells should affect neighbors 			*/
	/* 'x' and 'y' hold real-space values, will be mapped onto ineger array				*/

	{

	float	slope;			/* slope of zero goes staight back */
	int	ysign;			/* holder for going left or right alongshore */
	float	x,y;			/* holders for 'real' location of x and y */
	float	xin, yin;		/* starting 'real' locations */
	int	xtest,ytest;		/* cell looking at */
	float	xint,yint;		/* intercepts of overwash line in overwashable cell */
	int	NextXInt, NextYInt;	/* holder vairables for cell to check */
	float 	Ydown, DistanceDown;	/* when going to next x cell, what other values */
	float 	Xside, DistanceSide;	/* when gpoing to next y cell,other values */
	float	checkdistance;		/* distance of checking line- minimum, not actual width, ends loop */
	float	measwidth;		/* actual barrier width between cells */
	int 	AllBeachFlag;		/* flag to see if overwash line has passed over at least one AllBeach cell */

	/* convert angle to a slope and the direction of steps */
	/* note that for case of shoreline, positive angle will be minus y direcyion */

	/*if(icheck == 122)
		DEBUG_10A = 1;
	else
		DEBUG_10A = 0;*/

	if (_s.SurroundingAngle[icheck] == 0.0)
	{
		/* unlikely, but make sure no div by zero */
		slope = 0.00001;
	}
	else if (fabs(_s.SurroundingAngle[icheck]) == 90.0)
	{
		slope = 9999.9;
	}
	else
	{
		slope = fabs(tan(_s.SurroundingAngle[icheck]));
	}

	if (_s.SurroundingAngle[icheck] > 0)
		ysign = 1;
	else
		ysign = -1;
		
		DEBUG_PRINT( DEBUG_10A, "\nI: %d------------- Surr: %f  %f Slope: %f sign: %d \n",
		 icheck, _s.SurroundingAngle[icheck],_s.SurroundingAngle[icheck]*radtodeg,slope, ysign); 
	

	if (_s.AllBeach[_s.X[icheck]-1][_s.Y[icheck]] == 'y' || ((_s.AllBeach[_s.X[icheck]][_s.Y[icheck]-1] == 'y') && 
		(_s.AllBeach[_s.X[icheck]][_s.Y[icheck]+1] == 'y')) )
	/* 'regular condition' */
	/* plus 'stuck in the middle' situation (unlikely scenario)*/
	{
		xin = _s.X[icheck] + _s.PercentFull[_s.X[icheck]][_s.Y[icheck]];
		yin = _s.Y[icheck] + 0.5;
	}
	else if (_s.AllBeach[_s.X[icheck]][_s.Y[icheck]-1] == 'y')
	/* on right side */
	{
		xin = _s.X[icheck] + 0.5;
		yin = _s.Y[icheck] + _s.PercentFull[_s.X[icheck]][_s.Y[icheck]];
		DEBUG_PRINT( DEBUG_10A, "-- Right xin: %f  yin: %f\n",xin,yin);
	}
	else if (_s.AllBeach[_s.X[icheck]][_s.Y[icheck]+1] == 'y')
	/* on left side */
	{
		xin = _s.X[icheck] + 0.5;
		yin = _s.Y[icheck] + 1.0 - _s.PercentFull[_s.X[icheck]][_s.Y[icheck]];
		DEBUG_PRINT( DEBUG_10A, "-- Left xin: %f  yin: %f\n",xin,yin);
	}
	else	
	/* underneath, no overwash */
	{
		return;
	}

	x = xin;
	y = yin;
	checkdistance = 0 ; 
	AllBeachFlag = 0;

	while ((checkdistance < CritBWidth) && (y > 0) && (y < 2*Ymax) && (x > 1))
	{
		NextXInt = ceil(x) -1;
		if (ysign > 0)
			NextYInt = floor(y) + 1;
		else
			NextYInt = ceil(y-1);			 

		/* moving to next whole 'x' position, what is y position? */
		Ydown = y + (x - NextXInt)*slope * ysign;
		DistanceDown = Raise(((Ydown - y)*(Ydown - y) + (NextXInt - x)*(NextXInt - x)),.5);
	
		/* moving to next whole 'y' position, what is x position? */
		Xside = x - fabs(NextYInt - y) / slope;
		DistanceSide = Raise(((NextYInt - y)*(NextYInt - y) + (Xside - x)*(Xside - x)),.5);
	
			DEBUG_PRINT( DEBUG_10A, "x: %f  y: %f  X:%d  Y: %d  Yd: %f  DistD: %f Xs: %f DistS: %f\n",
			x,y,NextXInt,NextYInt, Ydown,DistanceDown,Xside,DistanceSide); 

		if (DistanceDown < DistanceSide)
		/* next cell is the down cell */
		{
			x = NextXInt;
			y = Ydown;
			xtest = NextXInt-1;	
			ytest = floor(y);
				/*DEBUG_PRINT( DEBUG_10A, " down ");*/
		}
		else
		/* next cell is the side cell */
		{
			x = Xside;
			y = NextYInt;
			xtest = floor(x);
			ytest = y + (ysign-1)/2;
				/*DEBUG_PRINT( DEBUG_10A, " side ");*/
		}
		
		/*if ((DEBUG_10A) && (DoGraphics == 'y'))PutPixel(ytest*CELL_PIXEL_SIZE,xtest*CELL_PIXEL_SIZE,0,0,200);*/

		checkdistance = Raise(((x - xin)*(x - xin) +  (y - yin)*(y - yin)),.5) * CellWidth;	
		if (_s.AllBeach[xtest][ytest] == 'y')
			AllBeachFlag = 1;

		DEBUG_PRINT( DEBUG_10A, "	x: %f  y: %f  xtest: %d ytest: %d check: %f\n\n",x,y,xtest,ytest,checkdistance);

		if ((_s.AllBeach[xtest][ytest] == 'n') && (AllBeachFlag) && !(((_s.X[icheck]-xtest) > 1) || (abs(ytest - _s.Y[icheck]) > 1)))
		/* if passed through an allbeach and a neighboring partial cell, jump out, only bad things follow */
		{
			return;
		}		
		
		if((_s.AllBeach[xtest][ytest] == 'n') && (AllBeachFlag) && (xtest < _s.X[icheck]) &&
			(((_s.X[icheck]-xtest) > 1) || (abs(ytest - _s.Y[icheck]) > 1)))
		/* Looking for shore cells, but don't want immediate neighbors, and go backwards */
		/* Also mush pass though an allbeach cell along the way */
		{
		
			if (_s.AllBeach[xtest+1][ytest] == 'y')
			/* 'regular condition' - UNDERNEATH, here */
			{
				xint = (xtest + 1 - _s.PercentFull[xtest][ytest]);
				yint = yin + (xin - xint)*ysign * slope;

				if ((yint > ytest+1.0) || (yint < ytest))
				/* This cell isn't actually an overwash cell */
				{
					measwidth = CritBWidth;
					DEBUG_PRINT( DEBUG_10A, "-- Regunder Cancelled  xin: %2.2f  yin: %2.2f xt:%d yt: %d xint: %f yint: %f sl: %2.2fMMeas: %3.2f\n",
					xin,yin,xtest,ytest,xint,yint,slope,measwidth);
				}
				else
				{
					measwidth = CellWidth * Raise((xint - xin)*(xint - xin)+ (yint - yin)*(yint - yin),0.5);

					DEBUG_PRINT( DEBUG_10A, "-- Regunder Over  xin: %2.2f  yin: %2.2f xt:%d yt: %d xint: %f yint: %f sl: %2.2fMeas: %3.2f\n",
					xin,yin,xtest,ytest,xint,yint,slope,measwidth);
				}
			}
			else if (_s.AllBeach[xtest][ytest-1] == 'y')
			/* on right side */
			{
				yint = (ytest + _s.PercentFull[xtest][ytest]);
				xint = xin - fabs(yin - yint)/ slope;

				if (xint < xtest)
				/* This cell isn't actually an overwash cell */
				{
					measwidth = CritBWidth;

					DEBUG_PRINT( DEBUG_10A, "-- Right Cancelled  xin: %2.2f  yin: %2.2f xt:%d yt: %d xint: %f yint: %f sl: %2.2fMMeas: %3.2f\n",
					xin,yin,xtest,ytest,xint,yint,slope,measwidth);
				}
				else
				{
					measwidth = CellWidth * Raise((xint - xin)*(xint - xin)+ (yint - yin)*(yint - yin),0.5);

					DEBUG_PRINT( DEBUG_10A, "-- Right Over  xin: %2.2f  yin: %2.2f xt:%d yt: %d xint: %f yint: %f sl: %2.2fMMeas: %3.2f\n",
					xin,yin,xtest,ytest,xint,yint,slope,measwidth);
				}
			}
			else if (_s.AllBeach[xtest][ytest+1] == 'y')
			/* on left side */
			{
				yint = (ytest + 1 - _s.PercentFull[xtest][ytest]);
				xint = xin - fabs(yin - yint)/ slope;

				if (xint < xtest)
				/* This cell isn't actually an overwash cell */
				{
					measwidth = CritBWidth;
					
					DEBUG_PRINT( DEBUG_10A, "-- Left cancelled  xin: %2.2f  yin: %2.2f xt:%d yt: %d xint: %f yint: %f sl: %2.2fMMeas: %3.2f\n",
					xin,yin,xtest,ytest,xint,yint,slope,measwidth);
				}
				else
				{
					measwidth = CellWidth * Raise((xint - xin)*(xint - xin)+ (yint - yin)*(yint - yin),0.5);

					DEBUG_PRINT( DEBUG_10A, "-- Left Over  xin: %2.2f  yin: %2.2f xt:%d yt: %d xint: %f yint: %f sl: %2.2fMMeas: %3.2f\n",
					xin,yin,xtest,ytest,xint,yint,slope,measwidth);
				}
			}
			else if (_s.AllBeach[xtest-1][ytest] == 'y')
			/* 'regular condition' */
			/* plus 'stuck in the middle' situation */
			{
				xint = (xtest  + _s.PercentFull[xtest][ytest]);
				yint = yin + (xin - xint)*ysign * slope;

				if ((yint > ytest+1.0) || (yint < ytest))
				/* This cell isn't actually an overwash cell */
				{
					measwidth = CritBWidth;
					DEBUG_PRINT( DEBUG_10A, "-- RegularODD Cancelled  xin: %2.2f  yin: %2.2f xt:%d yt: %d xint: %f yint: %f Meas: %3.2f\n",
					xin,yin,xtest,ytest,xint,yint,measwidth);
				}
				else
				{
					measwidth = CellWidth * Raise((xint - xin)*(xint - xin)+ (yint - yin)*(yint - yin),0.5);

					DEBUG_PRINT( DEBUG_10A, "-- RegularODD Over  xin: %2.2f  yin: %2.2f xt:%d yt: %d xint: %f yint: %f Meas: %3.2f\n",
					xin,yin,xtest,ytest,xint,yint,measwidth);
					/*PauseRun(xtest,ytest,icheck);*/
				}
			}
			else if (_s.PercentFull[xtest][ytest] > 0)
			/* uh oh - not good situation, no allbeach on sides */
			/* assume this is an empty cell, */
			{
				xint = x;
				yint = y;
			
				measwidth = CellWidth * Raise((xint - xin)*(xint - xin)+ (yint - yin)*(yint - yin),0.5);

				printf("-- Some Odd Over  xin: %2.2f  yin: %2.2f xt:%d yt: %d xint: %f yint: %f Meas: %3.2f Ang: %f Abs: %f\n",
				xin,yin,xtest,ytest,xint,yint,measwidth, _s.SurroundingAngle[icheck]*radtodeg,fabs(_s.SurroundingAngle[icheck])*radtodeg);
				/*PauseRun(xtest,ytest,icheck);*/
			}
			else
			/* empty cell - oughta fill er up  - fill max barrier width*/
			{
				xint = x;
				yint = y;
				measwidth = CritBWidth - CellWidth;
				
				printf("-- Empty Odd Over  xin: %2.2f  yin: %2.2f xt:%d yt: %d xint: %f yint: %f Meas: %3.2f Ang: %f Abs: %f\n",
				xin,yin,xtest,ytest,xint,yint,measwidth, _s.SurroundingAngle[icheck]*radtodeg,fabs(_s.SurroundingAngle[icheck])*radtodeg);
				/*PauseRun(xtest,ytest,icheck); */
			}

			checkdistance = measwidth;
			
			if (measwidth < CritBWidth)
			{
				DoOverwash(_s.X[icheck],_s.Y[icheck],xtest,ytest,xint,yint,measwidth,icheck);
				/* jump out of loop */
				OWflag = 1;
				return;
			}

			
		}
	
		}
/*	while(!getbutton(GKEY)){}*/
		
}

			

void DoOverwash(int xfrom,int yfrom, int xto, int yto, float xintto, float yintto, float widthin, int ishore)

	/*  given a cell where overwash is needed ,move sediment back  *** ADA 09/03, rev 01/04 */
	/*  for 'true' overwash based on shoreline angles					*/
	/*  will change and use _s.PercentFull[][] and AllBeach [][]				*/

{
	float BBneed, delBB, delShore; 	/* local variables */
	float MaxOver = 0.2; 		/*Maximum overwash step size (enforced at backbarrier) */
	/*float DepthBackBarrier = 6.0;  	 m current set depth for backbarrier (temp - make into function)*/
	float DepthBB;			/* holds effective backbarrier depth */
	short	vertex[2];
	
	DepthBB = GetOverwashDepth(xto,yto,xintto,yintto,ishore);

	/* calculated value of most that backbarrier ca nmove given geometry (true, non-iterative solution) */
	

	if (DepthBB == DepthShoreface)
	{
		BBneed = MaxOver;
	}
	else
	{
		BBneed = (CritBWidth - widthin) / CellWidth /  (1 - (DepthBB / DepthShoreface));
	}


	if (BBneed <= MaxOver)
	/* do all overwash */
	{
		delShore = BBneed * DepthBB / DepthShoreface;
		delBB = BBneed;
	}
	else
	/* only do overwash to max change) */
	{
		delShore = MaxOver * DepthBB / DepthShoreface ;
		delBB = MaxOver;
		
	}


	DEBUG_PRINT( DEBUG_10B, "** Overwash From X: %d  Y: %d  To: X: %d Y: %d Width: %f \n"
		, xfrom, yfrom,xto,yto,widthin );
	DEBUG_PRINT( DEBUG_10B, "DepthBB: %f  BBNeed: %f DelShore: %f  DelBB: %f\n",
			DepthBB, BBneed,delShore,delBB );
	/*if (DepthBB == DepthShoreface) PauseRun(xto,yto,-1);*/

#ifdef WITH_OPENGL
#ifdef DEBUG_ON
	if (DEBUG_10B && ( DO_GRAPHICS == 'y'))	
	{
		/*bgnpolygon();
			RGBcolor(250,0,0);
			vertex[0] = (yfrom+0.2)*CELL_PIXEL_SIZE;
			vertex[1] = (xfrom+.5)*CELL_PIXEL_SIZE;
			v2s(vertex);
			vertex[0] = (yfrom+0.8)*CELL_PIXEL_SIZE;
			v2s(vertex);
			vertex[0] = (yto+0.8)*CELL_PIXEL_SIZE;
			vertex[1] = (xto+0.5)*CELL_PIXEL_SIZE;
			v2s(vertex);
			vertex[0] = (yto+0.3)*CELL_PIXEL_SIZE;;
			v2s(vertex);
			vertex[0] = (yfrom+0.2)*CELL_PIXEL_SIZE;
			vertex[1] = (xfrom+.5)*CELL_PIXEL_SIZE;
			v2s(vertex);
			endpolygon();*/
	}
#endif
#endif

	_s.PercentFull[xto][yto] += delBB;
	_s.PercentFull[xfrom][yfrom] -= delShore;

	if (_s.PercentFull[xto][yto] > 1)
	{
		OopsImFull(xto,yto);
	}
	if (_s.PercentFull[xfrom][yfrom] < 0)
	{
		OopsImEmpty(xfrom,yfrom);
	}

#ifdef DEBUG_ON
	if (DEBUG_10B) PauseRun(xto,yto,-1);
#endif

}

	 
float GetOverwashDepth(int xin, int yin, float xinfl, float yinfl, int ishore)

	/* Rountine finds corresponding overwash depths					*/
	/*   	OWType = 0 take the depth at neightbor to the backing cell		*/
	/*   	OWType = 1 geometric rule based upon distance from back to shoreline	*/
	/* AA 5/04									*/
{
	int	xdepth;
	float 	Depth;
	float	BBDistance;  /* Distance from backshore to next shore */
	float	slope;			/* slope of zero goes staight back */
	int	ysign;			/* holder for going left or right alongshore */
	float	x,y;			/* holders for 'real' location of x and y */
	int	xtest,ytest;		/* cell looking at */
	int	NextXInt, NextYInt;	/* holder vairables for cell to check */
	float 	Ydown, DistanceDown;	/* when going to next x cell, what other values */
	float 	Xside, DistanceSide;	/* when gpoing to next y cell,other values */
	int 	BackFlag;		/* Flag to indicate if hit backbarrier */
	int	Backi = -1;		/* i for backbarrier intersection */
	int	i,j;			/* counters */
	int	FoundFlag;		/* Backbarrier intersection flag */
	float	AngleSin, AngleUsed;		


	if (OWType == 0)
	/* Use Cell Depths for overwash depths */
	{
		xdepth = xin;
		Depth = _s.CellDepth[xdepth][yin];
	
		while ((Depth < 0) && (xdepth > 0))
		{
			Depth = _s.CellDepth[xdepth][yin];
			/*printf("-- Overwash depth problem - Here = %f Next = %f",_s.CellDepth[xdepth][yto],_s.CellDepth[xdepth-1][yto] );
			PauseRun(xdepth,yto,-1);*/
			xdepth --;
		}
	
		if (Depth == DepthShoreface)
		{
			Depth = 6.0;
		}

		return Depth;
	
	}
	else if (OWType == 1)
	/* Geometric relation to determine depth through intersection of shorefaces 	*/
	/* look in line determined by shoreline slope - reuse stepping function (again)	*/
	{
		x = xinfl;
		y = yinfl;
				
		if (_s.SurroundingAngle[ishore] == 0.0)
		{
			/* unlikely, but make sure no div by zero */
			slope = 0.00001;
		}
		else if (fabs(_s.SurroundingAngle[ishore]) == 90.0)
		{
			slope = 9999.9;
		}
		else
		{
			slope = fabs(tan(_s.SurroundingAngle[ishore]));
		}

		BackFlag = 0;
		if (_s.SurroundingAngle[ishore] > 0)
			ysign = 1;
		else
			ysign = -1;
		
		while ((!BackFlag) && (y > 0) && (y < 2*Ymax) && (x > 1))
		{
			NextXInt = ceil(x) -1;
			if (ysign > 0)
				NextYInt = floor(y) + 1;
			else
				NextYInt = ceil(y-1);			 
	
			/* moving to next whole 'x' position, what is y position? */
			Ydown = y + (x - NextXInt)*slope * ysign;
			DistanceDown = Raise(((Ydown - y)*(Ydown - y) + (NextXInt - x)*(NextXInt - x)),.5);
		
			/* moving to next whole 'y' position, what is x position? */
			Xside = x - fabs(NextYInt - y) / slope;
			DistanceSide = Raise(((NextYInt - y)*(NextYInt - y) + (Xside - x)*(Xside - x)),.5);
		
				DEBUG_PRINT( DEBUG_10B, "x: %f  y: %f  X:%d  Y: %d  Yd: %f  DistD: %f Xs: %f DistS: %f\n",
				x,y,NextXInt,NextYInt, Ydown,DistanceDown,Xside,DistanceSide); 		

			if (DistanceDown < DistanceSide)
			/* next cell is the down cell */
			{
				x = NextXInt;
				y = Ydown;
				xtest = NextXInt-1;	
				ytest = floor(y);
			}
			else
			/* next cell is the side cell */
			{
				x = Xside;
				y = NextYInt;
				xtest = floor(x);
				ytest = y + (ysign-1)/2;
			}
		
			if (_s.PercentFull[xtest][ytest] > 0)
				BackFlag = 1;
		}
	
		/* Try to find the i for the cell found */
		/* If you have a better idea how to do this, go ahead */

		i = 2;
		FoundFlag = 0;
		
		while ((i < _s.TotalBeachCells-1) && !(FoundFlag))
		{
			if ((_s.X[i] == xtest) && (_s.Y[i] == ytest))
			{
				FoundFlag = 1;
				Backi = i;
			}
			i ++;
		}

	
		if (!BackFlag)
		/* The search for the backbarrier went out of bounds - not good, assume big = depthshoreface 	*/
		/* Periodic B.C.'s should make this not so important 						*/
		{
			Depth = DepthShoreface;
			DEBUG_PRINT( DEBUG_10B, "\nbackbarrier out of bounds: xin: %d yin: %d xbi: %d ybi: %d xinf: %f yinf: %f Per: %f Dist:  Depth: %f\n",
			xin, yin, xtest, ytest, xinfl, yinfl,_s.PercentFull[xtest][ytest], Depth);
			/*PauseRun(xin,yin,-1);*/
		}	
		else
		{
			BBDistance =  Raise(((xinfl - xtest-_s.PercentFull[xtest][ytest])*
			(xinfl - xtest-_s.PercentFull[xtest][ytest])) + 
			((yinfl - ytest - 0.5)*(yinfl - ytest - 0.5)),.5);

 			if(!FoundFlag)
			/* The backbarrier intersection isn't on the shoreline */
			/* Assume 1/2 of the length applies to this case */
			{
				Depth = BBDistance/2 * ShorefaceSlope * CellWidth;
				DEBUG_PRINT( DEBUG_10B, "\nNot Found backi: %d bx: %d by: %d Depth:%f",
					Backi,xtest,ytest,Depth);
			}
			else
			/* Use the fancy geometry thing */
			{
				AngleUsed = 0;
				for (j = -1; j <2 ; j ++)
				{
					AngleUsed += _s.SurroundingAngle[Backi+j];
				}
				AngleUsed = AngleUsed/5;

				if (fabs(AngleUsed) > M_PI/4.0)
				{
					AngleUsed = M_PI/4.0;
					DEBUG_PRINT( DEBUG_10B, "Big Angle");
					/*PauseRun(_s.X[Backi],_s.Y[Backi],Backi);*/
				}

				AngleSin = sin(M_PI/2.0 - fabs(_s.SurroundingAngle[ishore] + AngleUsed));

				Depth = BBDistance * AngleSin / (1 + AngleSin);
		
		DEBUG_PRINT( DEBUG_10B, "\nBack Angle backi: %d bx: %d by: %d BackA: %f AngU: %f Asin: %f L/2: %f Depth:%f",
			Backi,_s.X[Backi],_s.Y[Backi],_s.SurroundingAngle[ishore]*radtodeg,AngleUsed*radtodeg,AngleSin,
					BBDistance/2.0,Depth);
	
			}
		}

		if (Depth < OWMinDepth)
		{
			Depth = OWMinDepth;
		}
		else if (Depth > DepthShoreface)
		{
			Depth = DepthShoreface;
		}
			DEBUG_PRINT( DEBUG_10B, "\nOverwash Depth2: xin: %d yin: %d xbi: %d ybi: %d xinf: %f yinf: %f Per: %f Dist: %f  Depth: %f\n",
			xin, yin, xtest, ytest, xinfl, yinfl,_s.PercentFull[xtest][ytest],BBDistance, Depth);	
		return Depth;
	} 

	printf("OWDepth all broken");
	PauseRun(xin,yin,-1);
	return DepthShoreface;
}

