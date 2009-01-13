/*  CAPERIFFIC                                                                  */
/*  Program To generate capes?? and sandwaves?? using wave angle relationships  */
/*  Begun by Brad Murray 01/00                                                  */
/*  Refined by Olivier Arnoult 01/00 - 06/00                                    */
/*  Revised and reformed by Andrew Ashton 06/00 -                               */
/*                                                                              */      
/*  Program Notes -                                                             */
/*      To end program, press 'd' key and 'ESC' key simultaneously              */
/*      To save current iteration to file, press 's' and 'f' simultaneously     */
/*      To update screen display, press 'p' key                                 */
/*                                                                              */

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


/*  Run Control Parameters */

#define TimeStep        0.2      /* days - reflects rate of sediment transport per time step */
#define OffShoreWvHt    2     /* meters */
#define Period          7     /* seconds */
#define Asym            0.7    /*ractional portion of waves coming from positive (left) direction */
#define Highness        0.1    /*ll New! .5 = even dist, > .5 high angle domination */
#define Duration        1       /* Number of time steps calculations loop at same wave angle */
#define StopAfter       2600 /* Stop after what number of time steps */

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

int	debug0 = 0;	/* Main program steps */
int	debug1 = 0 ;	/* Find Next Cell */
int 	debug2 = 0;	/* Shadow Routine */
int 	debug3 = 0;	/* Determine Angles */
int 	debug4 = 0;	/* Upwind/Downwind */
int 	debug5 = 0;	/* Sediment Transport Decisions*/
int 	debug6 = 0;	/* Sediment Trans Calculations */
int 	debug7 = 0;	/* Transport Sweep (move sediment) */
int	debug7a = 0;	/* Slope Calcs */
int 	debug8 = 0;	/* Full/Empty */
int 	debug9 = 0;	/* FixBeach */
int	debug10a = 0;	/* Overwash Tests*/
int	debug10b = 0;	/* doing overwash (w/screen) */
int	OWflag = 0;		/* debugger */

/* Universal Constants */
//#define	pi		3.1415927
//#define	exp		2.7182818 /* e */
#define GRAV            (9.80665)
#define radtodeg        (180.0/M_PI) /* transform rads to degrees */
//float g =		9.80665;
//float radtodeg = 	180.0/pi; /* transform rads to degrees */


/* Overall Shoreface Configuration Arrays - Data file information */

char	AllBeach[Xmax][2*Ymax];		/* Flag indicating of cell is entirely beach */
float	PercentFull[Xmax][2*Ymax];	/* Fractional amount of shore cell full of sediment */
int	Age[Xmax][2*Ymax];		/* Age since cell was deposited */
float	CellDepth[Xmax][2*Ymax];	/* Depth array (m) (ADA 6/3) */


/* Computational Arrays (determined for each time step) */
 
int 	X[MaxBeachLength];		/* X Position of ith beach element */
int	Y[MaxBeachLength];		/* Y Position of ith beach element */
char	InShadow[MaxBeachLength];	/* Is ith beach element in shadow? */
float	ShorelineAngle[MaxBeachLength];	/* Angle between cell and right (z+1) neighbor	*/
float	SurroundingAngle[MaxBeachLength];/* Cell-orientated angle based upon left and right neighbor */
char	UpWind[MaxBeachLength];		/* Upwind or downwind condition used to calculate sediment transport */
float	VolumeIn[MaxBeachLength];	/* Sediment volume into ith beach element */	
float 	VolumeOut[MaxBeachLength];	/* Sediment volume out of ith beach element */


/* Miscellaneous Global Variables */

int	CurrentTimeStep = 0;  	/* Time step of current calculation */ 
/*FindBeachCells,FindNextCell*/int	NextX;			/* Global variables used to iterate FindNextCell in global array - */
/*FindBeachCells,FindNextCell*/int	NextY;			/*	would've used pointer but wouldn't work	*/
int 	TotalBeachCells;	/* Number of cells describing beach at particular iteration */
int 	ShadowXMax; 		/* used to determine maximum extent of beach cells */
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

// Graphics variables???
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

    ShadowXMax = Xmax-5;

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
	    scanf("%d", &CurrentTimeStep);
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

    while ( CurrentTimeStep < StopAfter )
    {
	/*  Time Step iteration - compute same wave angle for Duration time steps */

	/*  Calculate Wave Angle */

	WaveAngle = FindWaveAngle();

	/*  Loop for Duration at the current wave sign and wave angle */

	for (xx = 0; xx < Duration; xx++)
	{
		
	    /* Text to Screen? */

	    if (CurrentTimeStep%ScreenTextSpacing == 0)
	    {
		printf("==== WaveAngle: %2.2f  MASS Percent: %1.4f  Time Step: %d\n", 180*(WaveAngle)/M_PI, 
		       MassCurrent/MassInitial, CurrentTimeStep);
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
	
	    /* printf("Foundbeach!: %d \n", CurrentTimeStep); */

	    ShadowSweep();
	    if (debug0) printf("Shadowswept: %d \n", CurrentTimeStep);
	    DetermineAngles();
	    if (debug0) printf("AngleDet: %d \n", CurrentTimeStep);
	    DetermineSedTransport();
	    if (debug0) printf("Sed Trans: %d \n", CurrentTimeStep); 
	    TransportSedimentSweep();
	    if (debug0) printf("Transswept: %d \n", CurrentTimeStep);

	    DeliverSediment();

	    FixBeach();
	    if (debug0) printf("Fixed Beach: %d \n", CurrentTimeStep);


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
	
				/* printf("Foundbeach!: %d \n", CurrentTimeStep); */

			ShadowSweep();
				if (debug0) printf("Shadowswept: %d \n", CurrentTimeStep);
			DetermineAngles();
				if (debug0) printf("AngleDet: %d \n", CurrentTimeStep);
			CheckOverwashSweep();
			FixBeach();
	

	    if ((StartStop)   ) 
	    {
		printf("---- You Paused it, bud ---- \npp");
		PauseRun(1,1,-1);			
	    }

	    if(debug0) printf("End of Time Step: %d \n", CurrentTimeStep);
		
	    /* Age Empty Cells */

	    if ((CurrentTimeStep%AgeUpdate == 0) && SaveAge)
		AgeCells();

	    /* Count Mass */

	    MassCurrent = MassCount();

	    /* GRAPHING */		
		
#ifdef WITH_OPENGL
	    if ( DO_GRAPHICS && EveryPlotSpacing && (CurrentTimeStep%EveryPlotSpacing == 0))
		GraphCells();
#endif
	    

	    /* current_getch = getch();
	    printf("%d",current_getch);
if(DoGraphics && KeysOn && (current_getch == KEY_P))
GraphCells();*/

		
	    CurrentTimeStep ++;

	    /* SAVE FILE ? */

	    if (((CurrentTimeStep%SaveSpacing == 0 && CurrentTimeStep >= StartSavingAt) 
		 || (CurrentTimeStep == StopAfter)) && SaveFile)
	    {	
		SaveSandToFile();
	    }

	    if (((CurrentTimeStep%SaveLineSpacing == 0 && CurrentTimeStep >= StartSavingAt) 
		 || (CurrentTimeStep == StopAfter)) && SaveLine)
	    {
		SaveLineToFile();
	    }


			
	    /*if (CurrentTimeStep > 14300)
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
	
	/*printf("Time = %d RandBin = %f RandAng = %f\n",CurrentTimeStep, RandBin, RandAngle);*/

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
/* This function will affect and determine the global arrays:  X[] and Y[]	*/
/* This function calls FindNextCell   						*/
/* This will define TotalBeachCells for this time step				*/

{
    int 	y, z, xstart;	/* local iterators */	


    /* Starting at left end, find the x - value for first cell that is 'allbeach' */

    xstart = Xmax -1; y = YStart;
	
    while (AllBeach[xstart][y] == 'n')
    {	
	xstart -= 1;
    }
	
    xstart += 1;			/* Step back to where partially full beach */

    X[0] = xstart; 	Y[0] = YStart;

    if (debug1) printf("FirsX: %3d  FrstY: %3d  z: 0 \n", X[0], Y[0]);	

    z = 0;

    while ((Y[z] < 2*Ymax -1) && (z < MaxBeachLength-1)) 
    {
	z++;
	NextX = -2;
	NextY = -2;
			
	FindNextCell(X[z-1], Y[z-1], z-1);
	X[z] = NextX;
	Y[z] = NextY;
			
	if (debug1) printf("NextX: %3d  NextY: %3d  z: %d \n", NextX, NextY, z);

	if (PercentFull[X[z]][Y[z]] == 0) 
	{
	    printf("\nFINDBEACH: PercentFull Zero x: %d y: %d\n",X[z],Y[z]);
	    /*PauseRun(X[z],Y[z],z);*/
	}

	/* If return to start point or go off left side of array, going the wrong direction 	*/
	/* Jump off and start again closer to middle						*/

	if ((NextY < 1) || ((NextY == Y[0])&&(NextX==X[0])) || (z > MaxBeachLength -2))
	{
	    /*printf("!!!!!!!Fell Off!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! x = %d !!!!!!!!!!!!!", NextX);*/
	    FellOffArray = 'y';
	    ZeroVars();
	    return;
	}

			

	if (z > MaxBeachLength - 3)
	{
	    printf("????????????  went to end of MaxBeach!! ????");
	}

    }
	
    TotalBeachCells = z; 
    FellOffArray = 'n';		

    if (debug1) printf("Total Beach: %d  \n \n", TotalBeachCells); 

}


void FindNextCell(int x, int y, int z)

/* Function to find next cell that is beach moving in the general positive X direction */
/* changes global variables NextX and NextY, coordinates for the next beach cell       */
/* This function will use but not affect the global arrays:  AllBeach [][], X[], and Y[] */

{

    if ( AllBeach[x-1][y] == 'n')
	/* No beach directly beneath cell */
    {  	
	if ( AllBeach[x][y-1] == 'y' && AllBeach[x][y+1] == 'n') 
	    /* If on right side of protuberance */
	{
	    if ( AllBeach[x-1][y-1] == 'y' )
	    {  	/* Move one inshore */	
		NextX = x-1; NextY = y; return;
	    }	
	    else if (AllBeach[x-1][y-1] == 'n' )	/* This is where shadow procedure was */
	    {	/* Back and to the left */
		NextX = x-1; NextY = y-1; return;
	    }
	    printf("Should've found next cell (1): %d, %d \n", x, y);
	    PauseRun(x, y, z); 
	}
		

	else if ( AllBeach[x][y-1] == 'n' && AllBeach[x][y+1] == 'y')	
	    /* If on left side of protuberance */
	{
	    if ( AllBeach[x+1][y+1] == 'n' && AllBeach[x+1][y] == 'n')
		/*  Up and right - move around spit end */
	    {
		NextX = x+1; NextY = y+1; return;			
	    }
		
	    else if (  AllBeach[x+1][y] == 'y')
		/*  On underside of regular or diagonally thin spit */
	    {
		if ( AllBeach[x+1][y-1] == 'n' && AllBeach[x-1][y-1] == 'n' && X[z-1]>x)
		    /* Reaching end of spit - not going in circles */
		{
		    NextX = x-1; NextY = y; return;
		}
		else if (AllBeach[x+1][y-1] == 'n')
		    /* This is reaching end of spit */
		{
		    NextX = x+1; NextY = y-1; return;
		}
		/* Moving along back side of spit */
		{
		    NextX = x; NextY = y-1; return;
		}
	    }
		
	    else if ( AllBeach[x+1][y+1] == 'y')
		/* we know ( AllBeach[x+1][y] == 'n') */
		/* Moving straight up */
		/* NEW - we still don't want to go in */
	    {
		NextX = x+1; NextY = y; return;
	    }
					
	    printf("Should've found next cell (2): %d, %d \n", x, y);
	    PauseRun(x, y ,z);
	}


	if  (AllBeach[x][y-1] == 'n' && AllBeach[x][y+1] == 'n')
	    /* Hanging out - nothing on sides or top - maybe on corner? */
	{
	    if (AllBeach[x-1][y+1] == 'y' && AllBeach[x+1][y] == 'n')
		/* On left corner of protuberence, move right*/
	    {
		NextX = x; NextY = y+1; return;
	    }

	    else if (AllBeach[x+1][y] == 'y' && AllBeach[x+1][y-1] == 'n')
		/* Under protuberance, move around to left and up  */ 
	    {
		NextX = x+1; NextY = y-1; return;
	    }
			
	    else if (AllBeach[x+1][y] == 'y' && AllBeach[x+1][y-1] == 'y')
		/* Under protuberance, move to left */ 
	    {
		NextX = x; NextY = y-1; return;
	    }
	    printf("Should've found next cell (3): %d, %d \n", x, y);
	    PauseRun(x, y, z); 
	}

	else if ( AllBeach[x][y-1] == 'y' && AllBeach[x][y+1] == 'y' )
	    /* thin entrance between spits.  Don't even think about going in there */
	    /* (Similar case to over head and underneath - don't go in */
	    /* check to see which way we were coming in - from below or from side	*/
	{	if (X[z-1] > x)
	    /* coming from above */
	{
	    if (AllBeach[x+1][y+1] == 'n')
		/* Move right and up*/
	    {
		NextX = x+1; NextY = y+1; return;
	    }
	    else if (AllBeach[x+1][y] == 'n')
		/* Straight up*/
	    {
		NextX = x+1; NextY = y; return;
	    }
	    else if (AllBeach[x+1][y-1] == 'n')
		/* Up and left*/
		/* shouldn't need this, this where coming from */
	    {
		NextX = x+1; NextY = y-1; return;
	    }
	}
	else if (X[z-1] < x)
	    /* coming from below */
	{
	    if (AllBeach[x-1][y-1] == 'n')
		/* move down and left*/
	    {
		NextX = x-1; NextY = y-1; return;
	    }
	    else if (AllBeach[x-1][y] == 'n')
		/*move straight down*/
	    {
		NextX = x-1; NextY = y; return;
	    }
	    else if (AllBeach[x-1][y+1] == 'n')
		/*move straight down*/
		/* shouldn't need this, this would be where coming from*/
	    {
		NextX = x-1; NextY = y+1; return;
	    }
	}
	printf("Should've found next cell (3.5): %d, %d \n", x, y);		
	PauseRun(x, y, z); 
	}
		
	printf("Should've found next cell (4): %d, %d \n", x, y);
	PauseRun(x, y, z); 
    }

	
    else if ( AllBeach[x-1][y] == 'y' && AllBeach[x+1][y] == 'n') 
	/* There is beach beneath cell, nothing over the head */
    {
	if ( AllBeach[x][y+1] == 'n')
	    /*  Adjacent Cell to right is vacant */
	{
	    if ( AllBeach[x-1][y+1] == 'y' )
		/* move straight right */
	    {
		NextX = x; NextY = y+1; return;
	    }
	    else if ( AllBeach[x-1][y+1] == 'n' ) 
		/* Move down and to right */
	    {			 
		NextX = x-1; NextY = y+1; return;	
	    }

	    printf("Should've found next cell (5): %d, %d \n", x, y);
	    PauseRun(x, y, z); 
	}
		
	else if ( AllBeach[x][y+1] == 'y')	
	    /*Brad's note : DON'T REALLY NEED TO REPEAT THIS (WORKS SAME IN BOTH CASES) */
	    /* Right neighbor occupied */
	{				
	    if ( AllBeach[x+1][y+1] == 'n' )
		/* Move up and to right */
	    {
		NextX = x+1; NextY = y+1; return;
	    }
	    else if ( AllBeach[x+1][y+1] == 'y')
		/* Move straight up */
	    {
		NextX = x+1; NextY = y; return;
	    }
			
	    printf("Should've found next cell (6): %d, %d \n", x, y);
	    PauseRun(x, y, z); 	
	}

	printf("Should've found next cell (7): %d, %d \n", x, y);
	PauseRun(x, y, z); 
    }


    else if ( (AllBeach[x-1][y] == 'y') && (AllBeach[x+1][y] == 'y')) 
	/* There is beach behind cell, and over the head don't want to go in (will be shadowed anyway */
	/* Need to use last cell to find out if going into left or right enclosure */
    {
	/* Fill up that nasty piece of work */

	/*	FillUpGap(x,y, y-Y[z-1]);*/


	if (Y[z-1] < y)
	    /* Moving towards right, bump up and over the problem */
	{
	    if (AllBeach[x+1][y-1] == 'n')
		/* Move up and to the left */
	    {
		NextX = x+1;  NextY = y-1; return;
	    }
	    else if (AllBeach[x][y-1] == 'n')
		/* Move directly left */	
	    {
		NextX = x;  NextY = y-1; return;
	    }
	    else if (AllBeach[x-1][y-1] == 'n')
		/* Move left and down */	
	    {
		NextX = x-1;  NextY = y-1; return;
	    }
	    printf("Should've found next cell (8): %d, %d \n", x, y);
	    PauseRun(x, y, z); 
	}

	else if (Y[z-1] > y)
	    /* Moving towards left, go back right */
	{
	    if (AllBeach[x-1][y+1] == 'n')
		/* Move down and to the right */
	    {
		NextX = x-1;  NextY = y+1; return;
	    }
	    else if (AllBeach[x][y+1] == 'n')
		/* Move directly right */	
	    {
		NextX = x;  NextY = y+1; return;
	    }
	    else if (AllBeach[x+1][y+1] == 'n')
		/* Move right and up */	
	    {
		NextX = x+1;  NextY = y+1; return;
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
if (AllBeach[X][Y+2*LorR] == 'n')
{
    AllBeach[X][Y+2*LorR] = 'y';
    PercentFull[X][Y+2*LorR] = 1;
    printf("!!!!!!!!!!!!!!!!!\n		FILLEDERUP: %d, %d, %d \n", X, Y, LorR);
}	
}*/


void ShadowSweep(void)
	
/*  Moves along beach and tests to see if cells are in shadow 		*/
/*  This function will use and determine the Global array:  InShadow[]	*/
/*  This function will use and adjust the variable:   ShadowXMax	*/
/*  This function will use but not adjust the variable:  TotalBeachCells */

{

    int	i;

    /* Find maximum extent of beach to use as a limit for shadow searching */
	
    ShadowXMax = XMaxBeach(ShadowXMax) + 3;
	
    if (debug2) printf("ShadowXMax: %d   XMaxBeach: %d \n", ShadowXMax, XMaxBeach(ShadowXMax));

    /* Determine if beach cells are in shadow */

    for (i=0;  i <= TotalBeachCells; i++)
    {
	InShadow[i] = FindIfInShadow(i, ShadowXMax);	
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
	    if (AllBeach[xtest][ytest] == 'y')
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
/*	AllBeach[][] and PercentFull[][]					*/
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
    int 	debug2a = 0;		/* local debuggers */
    int	debug2b = 0;
	
	
    /* convert angle to a slope and the direction of steps */
    /* note that for case of shoreline, positive angle will be minus y direction */
    /*if (icheck == 106) {debug2a = 1;debug2b=1;}*/

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
		
    if (debug2a) printf("\nI: %d----------x: %d  Y: %d  Wang:  %f Slope: %f sign: %d \n",
			icheck, X[icheck],Y[icheck],WaveAngle*radtodeg,slope, ysign); 
	
    /* 03/04 AA: depending on local orientations, starting point will differ */
    /* so go through scenarios */

    xinint = X[icheck];
    yinint = Y[icheck];

    if (AllBeach[xinint-1][yinint] == 'y' || ((AllBeach[xinint][yinint-1] == 'y') && 
					      (AllBeach[xinint][yinint+1] == 'y')) )
	/* 'regular condition' */
	/* plus 'stuck in the middle' situation (unlikely scenario)*/
    {
	xin = xinint + PercentFull[xinint][yinint];
	yin = yinint + 0.5;
	if (debug2a) printf("-- Regular xin: %f  yin: %f\n",xin,yin);
    }
    else if (AllBeach[xinint][yinint-1] == 'y')
	/* on right side */
    {
	xin = xinint + 0.5;
	yin = yinint + PercentFull[xinint][yinint];
	if (debug2a) printf("-- Right xin: %f  yin: %f\n",xin,yin);
    }
    else if (AllBeach[xinint][yinint+1] == 'y')
	/* on left side */
    {
	xin = xinint + 0.5;
	yin = yinint + 1.0 - PercentFull[xinint][yinint];
	if (debug2a) printf("-- Left xin: %f  yin: %f\n",xin,yin);
    }
    else if (AllBeach[xinint+1][yinint] == 'y')
	/* gotta be on the bottom now */
    {
	xin = xinint + 1 - PercentFull[xinint][yinint];
	yin = yinint + 0.5;
	if (debug2a) printf("-- Under xin: %f  yin: %f\n",xin,yin);
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
	
	if (debug2a) printf("x: %f  y: %f  X:%d  Y: %d  Yd: %f  DistD: %f Xs: %f DistS: %f\n",
			    x,y,NextXInt,NextYInt, Yup,DistanceUp,Xside,DistanceSide); 

	if (DistanceUp < DistanceSide)
	    /* next cell is the up cell */
	{
	    x = NextXInt;
	    y = Yup;
	    xtestint = NextXInt;	
	    ytestint = floor(y);
	    if (debug2a) printf(" up ");
	}
	else
	    /* next cell is the side cell */
	{
	    x = Xside;
	    y = NextYInt;
	    xtestint = floor(x);
	    ytestint = y + (ysign-1)/2;
	    if (debug2a) printf(" side ");
	}
		
	if (debug2a) printf("	x: %f  y: %f  xtesti: %d ytesti: %d \n\n",x,y,xtestint,ytestint); 
			
			
	/* Now Test */
	/* If AllBeach is along the way, will we pass through 'diamond'?	*/
	/* Trick - if crossing through the diamond, will change quadrants 	*/
	/* Probably won't get to this one, though			 	*/		

	if 	(AllBeach[xtestint][ytestint] == 'y')
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

	    /*if (debug2a) printf("In Allbeach xin: %2.2f yin: %2.2f xout: %2.2f yout: %2.2f\n",
	      x,y,xout,yout);
	      if (debug2a) printf("In Allbeach xin: %2.2f yin: %2.2f xout: %2.2f yout: %2.2f\n",
	      (xout-xtestint-0.5),(x-xtestint-0.5),(yout-ytestint-0.5),(y-ytestint-0.5));*/

	    if(( (xout-xtestint-0.5) * (x-xtestint-0.5) < 0 ) || ((yout-ytestint-0.5) * (y-ytestint-0.5) < 0)) 
	    {
		if (debug2a) printf("  Shaddowded ");
		return 'y';
	    }
	}	

	

	/* Compare a partially full cell's x - distance to a line projected  	*/
	/* from the starting beach cell's x-distance				*/
	/* This assumes that beach projection is in x-direction (not too bad) 	*/

	else if ( PercentFull[xtestint][ytestint] > 0 ) 
	{
	    if (AllBeach[xtestint-1][ytestint] == 'y' || ((AllBeach[xtestint][ytestint-1] == 'y') && 
							  (AllBeach[xtestint][ytestint+1] == 'y')) )
		/* 'regular' condition */
		/* plus 'stuck in the middle' situation (unlikely scenario) */
	    {
		xtest = xtestint + PercentFull[xtestint][ytestint];
		ytest = ytestint + 0.5;

		if (xtest > (xin  + fabs(ytest-yin)/slope) ) 
		{
		    if (debug2b) printf("Top: sl: %f xt: %2.2f xin: %2.2f yt: %2.2f yin: %2.2f comp: %2.2f > Thing: %2.2f\n",
					slope, xtest, xin, ytest, yin, xtest, (xin  + fabs(ytest-yin)/slope));
		    return 'y';
		}
	    }
	    else if (AllBeach[xtestint][ytestint-1] == 'y')
		/* on right side */
	    {
		xtest = xtestint + 0.5;
		ytest = ytestint + PercentFull[xtestint][ytestint];

		if (ytest > (yin + (xtest-xin) * slope))
		{
		    if (debug2b) printf("Right:  xt: %f  yt: %f  comp: %f > Thing: %f\n",
					xtest, ytest, ytest,(yin + (xtest-xin) * slope));
		    return 'y';
		}
	    }
	    else if (AllBeach[xtestint][ytestint+1] == 'y')
		/* on left side */
	    {
		xtest = xtestint + 0.5;
		ytest = ytestint + 1.0 - PercentFull[xtestint][ytestint];

		if (ytest < (yin + (xtest-xin) * slope))
		{
		    if (debug2b) printf("Left:  xt: %f  yt: %f  comp: %f < Thing: %f\n",
					xtest, ytest, ytest,(yin + (xtest-xin) * slope));
		    return 'y';
		}
	    }
	    else if (AllBeach[xtestint+1][ytestint] == 'y')
		/* gotta be on the bottom now */
	    {
		xtest = xtestint + 1 - PercentFull[xtestint][ytestint];
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
/*		ShorelineAngle[], UpWind[], SurroundingAngle[]						*/
/*  This function will use but not affect the following arrays and values:		*/
/*		X[], Y[], PercentFull[][], AllBeach[][], WaveAngle			*/
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
	
    x2 = X[0] + PercentFull[X[0]][Y[0]];
    y2 = Y[0] + 0.5;

    /* Compute ShorelineAngle[]  */
    /* 	not equal to TotalBeachCells because angle between cell and rt neighbor */

    for (i=0 ; i < TotalBeachCells ; i++)
    {
		
	x1 = x2;
	y1 = y2;

	x2int = X[i+1];
	y2int = Y[i+1];

	if (AllBeach[x2int-1][y2int] == 'y' || ((AllBeach[x2int][y2int-1] == 'y') && 
						(AllBeach[x2int][y2int+1] == 'y')) && (AllBeach[x2int+1][y2int] == 'n'))
	    /* 'regular condition' - if between  */
	    /* plus 'stuck in the middle' situation (unlikely scenario)*/
	{
	    x2 = x2int + PercentFull[x2int][y2int];
	    y2 = y2int + 0.5;
	    if (debug3a) printf("-- Regular xin: %f  yin: %f\n",x2,y2);
	}
	else if ((AllBeach[x2int+1][y2int] == 'y') && (AllBeach[x2int-1][y2int] == 'y'))
	    /* in a sideways nook (or is that a cranny?) */
	{
	    x2 = x2int + 0.5;

	    if (AllBeach[x2int][y2int-1] == 'y')
		/* right-facing nook */
	    {
		y2 = y2int + PercentFull[x2int][y2int];
	    }
	    else
		/* left-facing nook */
	    {
		y2 = y2int + 1.0 - PercentFull[x2int][y2int];
	    }
	    if (debug3a) printf("-- Nook  xin: %f  yin: %f\n",x2,y2);
	}
	else if (AllBeach[x2int][y2int-1] == 'y')
	    /* on right side */
	{
	    x2 = x2int + 0.5;
	    y2 = y2int + PercentFull[x2int][y2int];
	    if (debug3a) printf("-- Right xin: %f  yin: %f\n",x2,y2);
	}
	else if (AllBeach[x2int][y2int+1] == 'y')
	    /* on left side */
	{
	    x2 = x2int + 0.5;
	    y2 = y2int + 1.0 - PercentFull[x2int][y2int];
	    if (debug3a) printf("-- Left xin: %f  yin: %f\n",x2,y2);
	}
	else if (AllBeach[x2int+1][y2int] == 'y')
	    /* gotta be on the bottom now */
	{
	    x2 = x2int + 1 - PercentFull[x2int][y2int];
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
	    ShorelineAngle[i] = atan((x2 - x1) / (y2 - y1));
	    if (debug3) printf("(R) i = %d  X[i]: %d Y[i]: %d Percent %3f x: %f y: %f Angle:%f  Deg Angle: %f \n",
			       i, X[i], Y[i], PercentFull[X[i]][Y[i]],x2,y2,ShorelineAngle[i], ShorelineAngle[i]*180/M_PI);
	}
	else if (y2 == y1)
	{
	    ShorelineAngle[i] = M_PI/2.0 * (x1 - x2) / fabs(x2 - x1);
	    if (debug3) printf("(G) i = %d  X[i]: %d Y[i]: %d Percent %3f x: %f y: %f Angle:%f  Deg Angle: %f \n",
			       i, X[i], Y[i], PercentFull[X[i]][Y[i]],x2,y2,ShorelineAngle[i], ShorelineAngle[i]*180/M_PI);
	}
	else 
	    /* y2 < y1 */
	{
	    ShorelineAngle[i] = atan((x2 - x1) / (y2 - y1)) - M_PI;
		
	    if (ShorelineAngle[i] < - M_PI)
	    {
		ShorelineAngle[i] += 2.0 * M_PI;
	    }
	    if (debug3) printf("(U) i = %d  X[i]: %d Y[i]: %d Percent %3f x: %f y: %f Angle:%f  Deg Angle: %f \n",
			       i, X[i], Y[i], PercentFull[X[i]][Y[i]],x2,y2,ShorelineAngle[i], ShorelineAngle[i]*180/M_PI);
	}

    }

    for (k=1 ; k < TotalBeachCells ; k++)
    {
	/* compute SurroundingAngle array */
	/* 02/04 AA averaging doesn't work on bottom of spits */
	/* Use trick that x is less if on bottom of spit - angles might be different signs as well */
		
	if ((Y[k-1] - Y[k+1] == 2) && 
	    (copysign(ShorelineAngle[k-1],ShorelineAngle[k]) != ShorelineAngle[k-1]))
	{		
	    SurroundingAngle[k] = (ShorelineAngle[k-1] + ShorelineAngle[k]) / 2 + M_PI;
	    if (SurroundingAngle[k] > M_PI)
	    {
		SurroundingAngle[k] -= 2.0 * M_PI;
	    }
	    if (debug4) printf("Under: %d\n",k);
	}
	else
	{
	    SurroundingAngle[k] = (ShorelineAngle[k-1] + ShorelineAngle[k]) / 2;
	}
    }
	
    /* Determine Upwind/downwind condition						*/
    /* Note - Surrounding angle is based upon left and right cell neighbors, 	*/
    /* and is centered on cell, not on right boundary				*/
	

    if (debug4) printf("\nUp/Down   Wave Angle:%f\n", WaveAngle * radtodeg);

    for (j=1 ; j < TotalBeachCells  ; j++)
    {
	if (debug4) printf("i: %d  Shad: %c Ang[i]: %3.1f  Sur: %3.1f  Effect: %3f  ",
			   j,InShadow[j], ShorelineAngle[j]*radtodeg, 
			   SurroundingAngle[j]*radtodeg, (WaveAngle - SurroundingAngle[j])*radtodeg);

	if ( fabs(WaveAngle - SurroundingAngle[j]) >= 42.0/radtodeg )
	{	
	    UpWind[j] = 'u';
	    if (debug4) printf("U(1)  ");
	}
	else 
	{
	    UpWind[j] = 'd';
	    if (debug4) printf("D(1)  ");
	}

	if (debug4) printf("\n");

    }

}

#define SED_TRANS_LIMIT (90) /* beyond what absolute slope don't do sed trans (degrees)*/
//float   SedTansLimit =  90;	/* beyond what absolute slope don't do sed trans (degrees)*/

void DetermineSedTransport(void)

/*  Loop function to determine which neigbor/situation to use for sediment transport calcs	*/
/*  Once situation is determined, will use function SedTrans to determine actual transport	*/
/*  This function will call SedTrans which will determine global arrays:			*/
/*		VolumeIn[], VolumeOut[]								*/
/*  This function will use but not affect the following arrays and values:			*/
/*		X[], Y[], InShadow[], UpWind[], ShorelineAngle[]				*/
/*  		PercentFull[][], AllBeach[][], WaveAngle					*/

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


    if (debug5) printf("\nSEDTRANS: %d  @  %f \n\n", CurrentTimeStep, WaveAngle * radtodeg);

    for (i=1 ; i < TotalBeachCells-1 ; i++)
    {
	if (debug5) printf("\n  i: %d  ",i);		


	MaxTrans = 'n';

	/*  Is littoral transport going left or right?	*/

	if ((WaveAngle-ShorelineAngle[i]) > 0)
	{
	    /*  Transport going right, center on cell to left side of border	 */
	    /*  Next cell in positive direction, no correction term needed 		*/
	    CalcCell = i;
	    Next = 1;
	    Last = -1;
	    Correction = 0;

	    if (debug5) printf("RT  %d ",CalcCell);
	}
	else
	{
	    /*  Transport going left, center on cell to right side of border 	*/
	    /*  Next cell in negative direction, correction term needed 		*/
	    CalcCell = i+1;
	    Next = -1;
	    Last = 1;
	    Correction = -1;

	    if (debug5) printf("LT  %d ",CalcCell);
	}
			

	if ( InShadow[CalcCell] == 'n') 
	{
			
	    /*  Adjustment for maximum transport when passing through 45 degrees		*/
	    /*  This adjustment is only made for moving from downwind to upwind conditions	*/
	    /* 										*/
	    /*  purposefully done before shadow adjustment, only use maxtran when 		*/
	    /*	transition from dw to up not because of shadow				*/
	    /* keeping transition from uw to dw - does not seem to be big deal (04/02 AA) */
			
	    if ( ((UpWind[CalcCell] == 'd') && (UpWind[CalcCell+Next] == 'u') &&
		  (InShadow[CalcCell + Next] == 'n')) ||
		 ((UpWind[CalcCell+Last] == 'u') && (UpWind[CalcCell] == 'd')
		  && (InShadow[CalcCell+Last] == 'n')) )
	    {
		MaxTrans = 'y';
		if (debug5) printf("MAXTRAN  ");
	    } 
			

	    /*  Upwind/Downwind adjustment Make sure sediment is put into shadows		*/
	    /*  If Next cell is in shadow, use UpWind condition				*/
			
	    DoFlux = 1;
	    UpWindLocal = UpWind[CalcCell];

	    if (InShadow[CalcCell+Next] == 'y')  
	    {
		UpWindLocal = 'u';
		if (debug5) printf("U(2)  ");
	    }

	    /*  If coming out of shadow, downwind should be used		*/
	    /*  HOWEVER- 02/04 AA - if high angle, will result in same flux in/out problem */
	    /*  	solution  - no flux for high angle waves */
		
	    if ((InShadow[CalcCell+Last] == 'y') &&(UpWindLocal == 'u')) 
	    {
		DoFlux = 0;
		if (debug5) printf("U(X) NOFLUX \n");
				
	    }

	    /*  Use upwind or downwind shoreline angle for calcs			*/

	    if (UpWindLocal == 'u')
	    {
		ShoreAngleUsed = ShorelineAngle[CalcCell+Last+Correction];
		if (debug5) printf("UP  ShoreAngle: %3.1f  ", ShoreAngleUsed * radtodeg);
	    }	
	    else if (UpWindLocal == 'd')
	    {
		ShoreAngleUsed = ShorelineAngle[CalcCell+Correction];
		if (debug5) printf("DN  ShoreAngle: %3.1f  ", ShoreAngleUsed *radtodeg);
	    }
			

	    /* !!! Do not do transport on unerneath c'cause it gets all messed up */
	    if (fabs(ShoreAngleUsed) > SedTansLimit/radtodeg)
	    {
		DoFlux = 0;
	    }
	
	    /* Send to SedTrans to calculate VolumeIn and VolumeOut*/


	    /* printf("i = %d  Cell: %d NextCell: %d Angle: %f Trans Angle: %f\n",
	       i, CalcCell, CalcCell+Next, ShoreAngleUsed*180/pi, (WaveAngle - ShoreAngleUsed)*180/pi); */

	    if (debug5) printf("From: %d  To: %d  TransAngle %3.1f", CalcCell, CalcCell+Next, 
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
/*		VolumeIn[] and VolumeOut[]							*/
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
	
    if (debug6) printf("Wave Angle %2.2f Shore Angle  %2.2f    ",WaveAngle*radtodeg, ShoreAngle*radtodeg);

    AngleDeep = WaveAngle - ShoreAngle;

    if (MaxT == 'y')
    {
	AngleDeep = 42.0 / radtodeg;
    }
    if (debug6) printf("Deep Tranport Angle %2.2f \n\n",AngleDeep*radtodeg);	

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
	if (debug6) printf("CDeep = %2.2f LDeep = %2.2f \n",CDeep, LDeep);

	while(!Broken)
	{
	    /* non-iterative eqn for L, from Fenton & McKee 		*/	

	    WaveLength = LDeep * Raise(tanh(Raise(Raise(2.0*M_PI/Period,2)*Depth/GRAV,.75)),2.0/3.0);
	    C = WaveLength/Period;
	    if (debug6) printf("DEPTH: %2.2f Wavelength = %2.2f C = %2.2f ", Depth, WaveLength,C);
			
	    /* Determine n = 1/2(1+2kh/tanh(kh)) Komar 5.21			*/
	    /* First Calculate kh = 2 pi Depth/L  from k = 2 pi/L		*/

	    kh =  M_PI * Depth / WaveLength;
	    n =0.5 * ( 1 + 2.0 * kh / sinh(2.0*kh));
	    if (debug6) printf("kh: %2.3f  n: %2.3f ", kh, n);

	    /* Calculate angle, assuming shore parallel contours and no conv/div of rays 	*/
	    /* from Komar 5.47								*/

	    Angle = asin(C/CDeep * sin(AngleDeep));
	    if (debug6) printf("Angle: %2.2f",Angle*radtodeg);

	    /* Determine Wave height from refract calcs - Komar 5.49			*/
	
	    WvHeight = OffShoreWvHt * Raise(CDeep*cos(AngleDeep)/(C*2.0*n*cos(Angle)),.5);
	    if (debug6) printf(" WvHeight : %2.3f\n",WvHeight);

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
				
	VolumeOut[From] = VolumeOut[From] + VolumeAcrossBorder;
		
	VolumeIn[To] = VolumeIn[To] + VolumeAcrossBorder;
	
	if (debug6) printf("VolumeAcrossBorder: %f  ",VolumeAcrossBorder);
	if (debug6) printf("VolumeIn : %f ",VolumeIn[To]);
	if (debug6) printf("VolumeOut : %f \n\n",VolumeOut[From]);
		
    }
}


void TransportSedimentSweep(void)

/*  Sweep through cells to place transported sediment				*/
/*  Call function AdjustShore() to move sediment.  				*/
/*  If cell full or overempty, call OopsImFull or OopsImEmpty()			*/
/*  This function doesn't change any values, but the functions it calls do	*/
/*  Uses but doesn't change:  X[], Y[], PercentFull[]				*/
/*  sweepsign added to ensure that direction of actuating changes does not  	*/
/*  	produce unwanted artifacts (e.g. make sure symmetrical			*/

{

    int i,ii;
    int sweepsign;

    if (RandZeroToOne()*2 > 1)
    {
	sweepsign = 1;
	if (debug7) printf("L  ");
    }
    else
    {
	sweepsign = 0;
	if (debug7) printf("R  ");
    }

    if (debug7) printf("\n\n TransSedSweep  Ang %f  %d\n", WaveAngle * radtodeg, CurrentTimeStep);
	
    for (i=0; i < TotalBeachCells-1 ; i++)
    {
	
	if (sweepsign == 1)
	    ii = i;
	else
	    ii = TotalBeachCells-1-i;

	if (debug7) printf("i: %d  ss: %d  X: %d  Y: %d  In: %.1f  Out: %.1f\n", ii, sweepsign,
			   X[i], Y[i], VolumeIn[i], VolumeOut[i]);

	AdjustShore(ii);
				
	if (PercentFull[X[ii]][Y[ii]] < 0)
	{
	    OopsImEmpty(X[ii],Y[ii]);
	}
	else if (PercentFull[X[ii]][Y[ii]]> 1) 
	{
	    OopsImFull(X[ii],Y[ii]);
	}
    }

}

void AdjustShore(int i)

/*  Complete mass balance for incoming and ougoing sediment			*/
/*  This function will change the global data array PercentFull[][]		*/
/*  Uses but does not adjust arrays:  						*/
/*		VolumeIn[], VolumeOut[], X[], Y[], ShorelineAngle[]		*/
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

    if (VolumeIn[i] <= VolumeOut[i])
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
	Xintfloat = X[i] + 0.5 + Distance * cos(SurroundingAngle[i]);
	Xintint = floor(Xintfloat);
	Yintfloat = Y[i] + 0.5 - Distance * sin(SurroundingAngle[i]);
	Yintint = floor(Yintfloat);

	if (debug7a)printf("xs: %d  ys: %d  Xint: %f Xint:%d Yint: %f Yint: %d  Dint: %f SAng: %f Sin = %f\n",
			   X[i],Y[i],Xintfloat,Xintint,Yintfloat,Yintint,CellDepth[Xintint][Yintint],SurroundingAngle[i]*radtodeg,sin(SurroundingAngle[i]));


	if ((Yintint < 0) || (Yintint > 2*Ymax))
	{
	    Depth = DepthShoreface;
	    if ((Yintint > Ymax/2) && (Yintint < 3/2*Ymax))
	    {
		printf("Periodic Boundary conditions and Depth Out of Bounds");
		PauseRun(X[i],Y[i],i);
	    }
	}
	else if ((Xintint < 0) || (Xintint > Xmax))
	{
	    Depth = DepthShoreface;
	    printf("-- Warning - depth location off of x array: X %d Y %d",Xintint,Yintint);
	    PauseRun(X[i],Y[i],i);
	}
	else if (CellDepth[Xintint][Yintint] <= 0)
	    /* looking back on land */
	{
	    Depth = DepthShoreface;
	    if (debug7a) printf("=== Shoreface is Shore, eh? Accreti:  xs: %d  ys: %d  Xint:%d  Yint: %d  Dint: %f \n",
				X[i],Y[i],Xintint,Yintint,CellDepth[Xintint][Yintint]);
	}
	else if (CellDepth[Xintint][Yintint] < DepthShoreface)
	{
	    printf("Shallow but underwater Depth %f",CellDepth[Xintint][Yintint]);
	    PauseRun(Xintint,Yintint,01);
	}
	else
	{
	    Depth = CellDepth[Xintint][Yintint];
			

	    /* That was the easy part - now we need to 'fix' all cells towards shoreface */
	    /* probably due to accretion from previous moving forward */
	    /* reuse some of the overwash checking code here */


	    if (SurroundingAngle[i] == 0)
	    {
		/* unlikely, but make sure no div by zero */
		slope = 0.00001;
	    }
	    else if (fabs(SurroundingAngle[i]) == 90.0)
	    {
		slope = 9999.9;
	    }
	    else
	    {
		slope = fabs(tan(SurroundingAngle[i]));
	    }

	    if (SurroundingAngle[i] > 0)
		ysign = 1;
	    else
		ysign = -1;
			
	    x = Xintfloat;
	    y = Yintfloat;
	    xtest = Xintint;
	    ytest = Yintint;
	    ShorefaceFlag = 0;	
	
	    while (( CellDepth[xtest][ytest] > DepthShoreface) && !(ShorefaceFlag))
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
						
		if (CellDepth[xtest][ytest] > DepthShoreface)	
		    /* Deep hole - fill 'er in - mass came from previous maths */ 
		{

		    if (debug7a) printf("=== Deep Hole, eh? Accreti:  xs: %d  ys: %d  Xint:%d  Yint: %d  Dint: %f Xfill: %d Yfill: %d Dt: %f\n",
					X[i],Y[i],Xintint,Yintint,CellDepth[Xintint][Yintint],xtest,ytest,
					CellDepth[xtest][ytest]); 
		    CellDepth[xtest][ytest] = DepthShoreface;
					
		    /*PauseRun(xtest,ytest,i);*/


		}
		else 
		    /* stop checking - ostensibly we have hit the shoreface or shore */
		{
		    ShorefaceFlag = 1;
				
		    if (PercentFull[xtest][ytest] > 0)
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

    DeltaArea = (VolumeIn[i] - VolumeOut[i])/Depth;

    PercentFull[X[i]][Y[i]] += DeltaArea/(CellWidth*CellWidth);
	
    PercentIn = VolumeIn[i]/(CellWidth*CellWidth*Depth);
    PercentOut = VolumeOut[i]/(CellWidth*CellWidth*Depth);
    PercentSum = DeltaArea/(CellWidth*CellWidth);

    if (debug7) printf("  In: %2.4f  Out: %2.4f  Sum: %2.4f\n", PercentIn, PercentOut, PercentSum);

}


void OopsImEmpty(int x, int y)

/*  If a cell is under-full, this will find source for desparity and move brach in	*/
/*  Function completly changed 5/21/02 sandrevt.c					*/
/*  		New Approach - steal from all neighboring AllBeach cells		*/
/*		Backup plan - steal from all neighboring percent full > 0		*/
/*  Function adjusts primary data arrays:						*/
/*		AllBeach[][] and PercentFull[][]					*/
	

{

    int emptycells = 0;
    int emptycells2 = 0;

    if (debug8) printf("\n		OOPS I'm EMPTY!  X: %d  Y: %d Per: %f ", x, y, PercentFull[x][y]);

    /* find out how many AllBeaches to take from */

    if (AllBeach[x-1][y] == 'y')
	emptycells += 1;
    if (AllBeach[x+1][y] == 'y')
	emptycells += 1;
    if (AllBeach[x][y-1] == 'y')
	emptycells += 1;
    if (AllBeach[x][y+1] == 'y')
	emptycells += 1;

    if (emptycells > 0)
    {
	/* Now Move Sediment */

	if (AllBeach[x-1][y] == 'y')
	{
	    PercentFull[x-1][y] += (PercentFull[x][y])/emptycells;
	    AllBeach[x-1][y] = 'n';
	    if (debug8) printf ("  MOVEDBACK");
	}		
	if (AllBeach[x+1][y] == 'y')
	{
	    PercentFull[x+1][y] += (PercentFull[x][y])/emptycells;
	    AllBeach[x+1][y] = 'n';
	    if (debug8) printf ("  MOVEDUP");
	}
	if (AllBeach[x][y-1] == 'y')
	{
	    PercentFull[x][y-1] += (PercentFull[x][y])/emptycells;
	    AllBeach[x][y-1] = 'n';
	    if (debug8) printf ("  MOVEDLEFT");
	    /*if (debug8) PauseRun(x,y,-1);*/
	}
	if (AllBeach[x][y+1] == 'y')
	{
	    PercentFull[x][y+1] += (PercentFull[x][y])/emptycells;
	    AllBeach[x][y+1] = 'n';
	    if (debug8) printf ("  MOVEDRIGHT");
	    /*if (debug8) PauseRun(x,y,-1);*/
	}
    }
    else
    {
	/* No full neighbors, so take away from partially full neighbors */
		
	if (PercentFull[x-1][y] > 0)
	    emptycells2 += 1;
	if (PercentFull[x+1][y] > 0)
	    emptycells2 += 1;
	if (PercentFull[x][y-1] > 0)
	    emptycells2 += 1;
	if (PercentFull[x][y+1] > 0)
	    emptycells2 += 1;

	if (emptycells2 > 0)
	{

	    if (PercentFull[x-1][y] > 0)
	    {
		PercentFull[x-1][y] += (PercentFull[x][y])/emptycells2;
		if (debug8) printf ("  NOTFULL MOVEDBACK");
	    }		
	    if (PercentFull[x+1][y] > 0)
	    {
		PercentFull[x+1][y] += (PercentFull[x][y])/emptycells2;
		if (debug8) printf ("  NOTFULL MOVEDUP");
	    }
	    if (PercentFull[x][y-1] > 0)
	    {
		PercentFull[x][y-1] += (PercentFull[x][y])/emptycells2;
		if (debug8) printf ("  NOTFULL MOVEDLEFT");
		/*if (debug8) PauseRun(x,y,-1);*/
	    }
	    if (PercentFull[x][y+1] > 0)
	    {
		PercentFull[x][y+1] += (PercentFull[x][y])/emptycells2;
		if (debug8) printf ("  NOTFULL MOVEDRIGHT");
		/*if (debug8) PauseRun(x,y,-1);*/
	    }	
	}
	else
	{
	    printf("@@@ Didn't find anywhere to steal sand from!! x: %d  y: %d\n",x,y);
	    PauseRun(x,y,-1);
	}

    }

    AllBeach[x][y] = 'n';
    PercentFull[x][y] = 0.0;
    CellDepth[x][y] = DepthShoreface;

    if (debug8) printf("\n");

}


void OopsImFull(int x, int y)

/*  If a cell is overfull, push beach out in new direction				*/
/*  Completely revised 5/20/02 sandrevt.c to resolve 0% full problems, etc.		*/
/*  New approach: 	put sand wherever 0% full in adjacent cells			*/
/*			if not 0% full, then fill all non-allbeach			*/
/*  Function adjusts primary data arrays:						*/
/*		AllBeach[][] and PercentFull[][]					*/
	

{

    int fillcells = 0;
    int fillcells2 = 0;

    if (debug8) printf("\n		OOOPPPS I'M FULLL: X: %d  Y: %d Per: %f  ==", x, y, PercentFull[x][y]);
    /*if (debug8) PrintLocalConds(x,y,-1);*/

    /* find out how many cells will be filled up	*/

    if (PercentFull[x-1][y] == 0.0)
	fillcells += 1;
    if (PercentFull[x+1][y] == 0.0)
	fillcells += 1;
    if (PercentFull[x][y-1] == 0.0)
	fillcells += 1;
    if (PercentFull[x][y+1] == 0.0)
	fillcells += 1;

    if (fillcells != 0)
    {
	/* Now Move Sediment */

	if (PercentFull[x-1][y] == 0.0)
	{
	    PercentFull[x-1][y] += (PercentFull[x][y]-1)/fillcells;
	    CellDepth[x-1][y] = - LandHeight;
	    if (debug8) printf ("  MOVEDBACK");
	}		
	if (PercentFull[x+1][y] == 0.0)
	{
	    PercentFull[x+1][y] += (PercentFull[x][y]-1)/fillcells;
	    CellDepth[x+1][y] = - LandHeight;
	    if (debug8) printf ("  MOVEDUP");
	}
	if (PercentFull[x][y-1] == 0.0)
	{
	    PercentFull[x][y-1] += (PercentFull[x][y]-1)/fillcells;
	    CellDepth[x][y-1] = - LandHeight;
	    if (debug8) printf ("  MOVEDLEFT");
	    /*if (debug8) PauseRun(x,y,-1);*/
	}
	if (PercentFull[x][y+1] == 0.0)
	{
	    PercentFull[x][y+1] += (PercentFull[x][y]-1)/fillcells;
	    CellDepth[x][y+1] = - LandHeight;
	    if (debug8) printf ("  MOVEDRIGHT");
	    /*if (debug8) PauseRun(x,y,-1);*/
	}
    }
    else
    {
	/* No fully empty neighbors, so distribute to partially full neighbors */
		
	if (PercentFull[x-1][y] < 1)
	    fillcells2 += 1;
	if (PercentFull[x+1][y] < 1)
	    fillcells2 += 1;
	if (PercentFull[x][y-1] < 1)
	    fillcells2 += 1;
	if (PercentFull[x][y+1] < 1)
	    fillcells2 += 1;

	if (fillcells2 > 0)
	{

	    if (PercentFull[x-1][y] < 1)
	    {
		PercentFull[x-1][y] += (PercentFull[x][y]-1)/fillcells2;
		if (debug8) printf ("  MOVEDBACK");
	    }		
	    if (PercentFull[x+1][y] < 1)
	    {
		PercentFull[x+1][y] += (PercentFull[x][y]-1)/fillcells2;
		if (debug8) printf ("  MOVEDUP");
	    }
	    if (PercentFull[x][y-1] < 1)
	    {
		PercentFull[x][y-1] += (PercentFull[x][y]-1)/fillcells2;
		if (debug8) printf ("  MOVEDLEFT");
	    }
	    if (PercentFull[x][y+1] < 1)
	    {
		PercentFull[x][y+1] += (PercentFull[x][y]-1)/fillcells2;
		if (debug8) printf ("  MOVEDRIGHT");
	    }	
	}
	else
	{
	    if (debug8) printf("Nobody wants our sand!!! x: %d  y: %d Per: %f\n",x,y,PercentFull[x][y]);
	    /*PauseRun(x,y,-1);*/
	}

    }

    AllBeach[x][y] = 'y';
    PercentFull[x][y] = 1.0;
    CellDepth[x][y] = - LandHeight;

    if (debug8) printf("\n");
	
	
}


void FixBeach(void)

/* Hopefully addresses strange problems caused by filling/emptying of cells	*/
/* Looks at entire data set							*/
/* Find unattached pieces of sand and moves them back to the shore 		*/
/* Takes care of 'floating bits' of sand					*/
/* Also takes care of over/under filled beach pieces				*/
/* Revised 5/21/02 to move sand to all adjacent neighbors sandrevt.c 		*/
/* Changes global variable PercentFull[][]					*/
/* Uses but does not change AllBeach[][]					*/
/* sandrevx.c - added sweepsign to reduce chances of asymmetrical artifacts	*/


{

    int i,x,y,sweepsign;
    int FixXMax;
    int fillcells3 = 0;

    /*if (debug9) printf("\n\nFIXBEACH      %d     %f\n", CurrentTimeStep, WaveAngle*radtodeg);*/

    if (RandZeroToOne()*2 > 1)
    {
	sweepsign = 1;
	if (debug9) printf("fixL  ");
    }
    else
    {
	sweepsign = 0;
	if (debug9) printf("fixR  ");
    }


    FixXMax = ShadowXMax + ceil(DepthShoreface/CellWidth/ShorefaceSlope) +3;
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
	    if ((PercentFull[x][y] <= 0) && (CellDepth[x][y] > DepthShoreface) && 
		(CellDepth[x-1][y] == DepthShoreface))
	    {
		if ((CellDepth[x+1][y] == DepthShoreface) && (CellDepth[x][y-1] == DepthShoreface)
		    && (CellDepth[x][y+1] == DepthShoreface))
		{
		    /* Fill Hole */
		    CellDepth[x][y] = DepthShoreface;
		}
	    }
	    if (PercentFull[x][y]> 100)
	    {
		printf("too full");
		PercentFull[x][y] = 0;
		PauseRun(x,y,-1);
	    }
	


	    /* Take care of situations that shouldn't exist */

	    if (PercentFull[x][y] < 0)
	    {
		AllBeach[x][y] = 'n';
		if (debug9 && y != 0) printf("\nUnder 0 Percent X: %d  Y: %d Percent: %f\n", x,y,PercentFull[x][y]);
		OopsImEmpty(x,y);
		printf("Underzerofill");
		/*PauseRun(x,y,-1);*/
	    }

	    if (PercentFull[x][y] > 1)
	    {
		AllBeach[x][y] = 'y';
		CellDepth[x][y] = - LandHeight;
		if (debug9 && y != 0) printf("\nOver 100 Percent X: %d  Y: %d Per: %f\n"
					     ,x,y, PercentFull[x][y]);
		OopsImFull(x,y);
	    }

	    if (((PercentFull[x][y] >=0) && (PercentFull[x][y] <1)) && (AllBeach[x][y] == 'y'))
	    {
		AllBeach[x][y] = 'n';
		CellDepth[x][y] = - LandHeight;
		if (debug9 && y != 0) printf("\nALLBeachProb X: %d  Y: %d\n", x,y);
	    }


	    /* Take care of 'loose' bits of sand */

	    fillcells3 = 0;
			
	    if ( (PercentFull[x][y] != 0) && (PercentFull[x-1][y] < 1) && (PercentFull[x+1][y] < 1) &&
		 (PercentFull[x][y+1] < 1) && (PercentFull[x][y-1] < 1) && (AllBeach[x][y] =='n'))
		/* Beach in cell, but bottom, top, right, and left neighbors not all full */
	    {
		if (debug9 && y != 0) printf("\nFB Moved loose bit of sand,  X: %d  Y: %d  Per: %f  ",
					     x, y, PercentFull[x][y]);

		/* distribute to partially full neighbors */
		
		if ((PercentFull[x-1][y] < 1) && (PercentFull[x-1][y] > 0))
		    fillcells3 += 1;
		if ((PercentFull[x+1][y] < 1) && (PercentFull[x+1][y] > 0))
		    fillcells3 += 1;
		if ((PercentFull[x][y-1] < 1) && (PercentFull[x][y-1] > 0))
		    fillcells3 += 1;
		if ((PercentFull[x][y+1] < 1) && (PercentFull[x][y+1] > 0))
		    fillcells3 += 1;				

		if ((fillcells3 > 0))
		{			

		    if ((PercentFull[x-1][y] < 1) && (PercentFull[x-1][y] > 0))
		    {
			PercentFull[x-1][y] += (PercentFull[x][y])/fillcells3;
			if (debug9) printf ("  MOVEDBACK");
		    }		
		    if ((PercentFull[x+1][y] < 1) && (PercentFull[x+1][y] > 0))
		    {
			PercentFull[x+1][y] += (PercentFull[x][y])/fillcells3;
			if (debug9) printf ("  MOVEDUP");
		    }
		    if ((PercentFull[x][y-1] < 1) && (PercentFull[x][y-1] > 0))
		    {
			PercentFull[x][y-1] += (PercentFull[x][y])/fillcells3;
			if (debug9) printf ("  MOVEDLEFT");
			/*if (debug9) PauseRun(x,y,-1);*/
		    }
		    if ((PercentFull[x][y+1] < 1) && (PercentFull[x][y+1] > 0))
		    {
			PercentFull[x][y+1] += (PercentFull[x][y])/fillcells3;
			if (debug9) printf ("  MOVEDRIGHT");
			/*if (debug9) PauseRun(x,y,-1);*/
		    }	
		}
		else
		{
		    printf("Loner fixbeach breakdown - mass disintegrated x: %d  y: %d\n",x,y);
		    if (debug9) 
			PauseRun(x,y,-1);
		}	

		PercentFull[x][y] = 0;
		AllBeach[x][y] = 'n';
		CellDepth[x][y] = DepthShoreface;

		if (debug9) printf("\n");


		/* If we have overfilled any of the cells in this loop, need to OopsImFull() */

		if (PercentFull[x-1][y] > 1)
		{
		    OopsImFull(x-1,y);
		    if (debug9) printf("	Below Overfilled\n");
		}
		if (PercentFull[x][y-1] > 1)
		{
		    OopsImFull(x,y-1);
		    if (debug9) printf("	Left Side Overfilled\n");	
		}
		if (PercentFull[x][y+1] > 1)
		{
		    OopsImFull(x,y+1);
		    if (debug9) printf("	Right Side Overfilled\n");
		}
		if (PercentFull[x+1][y+1] > 1)
		{
		    OopsImFull(x+1,y+1);
		    if (debug9) printf("	Top Overfilled\n");
		}
			
	    }
		
	    /*if ((AllBeach[x][y] =='y') && (PercentFull[x-1][y] < 1) && (PercentFull[x+1][y] < 1)
	      && (PercentFull[x][y-1] < 1) && (PercentFull[x][y+1] < 1)
	      && (AllBeach[x-1][y-1] == 'n') && (AllBeach[x-1][y+1] == 'n') &&
	      (AllBeach[x+1][y+1] == 'n') && (AllBeach[x+1][y-1] == 'n') )

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
/* Uses AllBeach[][] and PercentFull[][]		*/
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
	    /*if ((PercentFull[x][y] > 0) && (PercentFull[x][y] < 1.0))
	      MassHere = PercentFull[x][y] * (refdepth - CellDepth[x][y]) + 
	      (1 - PercentFull[x][y])*(refdepth - DepthShoreface);
	      else 
	      MassHere = refdepth - CellDepth[x][y];*/

	    Mass += PercentFull[x][y];
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

		CellDepth[x][y] = InitialDepth + ((x-InitBeach) * CellWidth * ShelfSlope);
			
		if (x < InitBeach)
		{  	 
		    PercentFull[x][y] = 1;
		    AllBeach[x][y] = 'y';
		    CellDepth[x][y] = - LandHeight;
		}
		else if (x == InitBeach)
		{
		    if (InitialSmooth)
		    {
			PercentFull[x][y] = .5;
		    }
		    else
		    {
			PercentFull[x][y] = RandZeroToOne();
			printf("x: %d  Y: %d  Per: %f\n",x,y,PercentFull[x][y]);
		    }
		    AllBeach[x][y] = 'n';
		    CellDepth[x][y] = - LandHeight;
		}
		else if (x > InitBeach)
		{
		    PercentFull[x][y] = 0;
		    AllBeach[x][y] = 'n';
		    if (CellDepth[x][y] < DepthShoreface)
		    {
			CellDepth[x][y] = DepthShoreface;
		    }
		}
		else
		{
		printf("WTF! x: %d  Y: %d  Per: %f\n",x,y,PercentFull[x][y]);
		PauseRun(x,y,-1);
		}

		Age[x][y] = 0;		
	    }
    }

    else if (InitCType == 1)
	/* 'Simple Barrier' type initial condition - island backed by lagoon at slope of shelf */
    {
	
	for (y = 0; y <= 2*Ymax; y++)
	{
	    for (x = 0; x < Xmax; x++)
	    {

		CellDepth[x][y] = InitialDepth + ((x-InitBeach) * CellWidth * ShelfSlope);
			
		if (CellDepth[x][y] <= 0)
		    /* This must be land due to continental shelf intersection */
		{
		    PercentFull[x][y] = 1.0;
		    AllBeach[x][y] = 'y';
		    CellDepth[x][y] = - LandHeight;
		}
		else if (x > InitBeach)
		    /* Shoreward of beach - enforce ShorefaceDepth if necessary */
		{
		    PercentFull[x][y] = 0;
		    AllBeach[x][y] = 'n';
		    if (CellDepth[x][y] < DepthShoreface)
		    {
			CellDepth[x][y] = DepthShoreface;
		    }
		}
		else if (x == InitBeach)
		    /* Beach */
		{
		    if (InitialSmooth)
		    {
			PercentFull[x][y] = .5;
		    }
		    else
		    {
			PercentFull[x][y] = RandZeroToOne();
			/*printf("x: %d  Y: %d  Per: %f\n",x,y,PercentFull[x][y]);*/
		    }
		    AllBeach[x][y] = 'n';
		    CellDepth[x][y] = - LandHeight;
		}
		else if ((x < InitBeach) && (x > InitBeach - InitBWidth - 1))
		    /* Island */
		{  	 
		    PercentFull[x][y] = 1.0;
		    AllBeach[x][y] = 'y';
		    CellDepth[x][y] = - LandHeight;
		}
		else if (x == InitBeach - InitBWidth -1)
		    /* Back of Barrier */
		{  	 
		    if (InitialSmooth)
		    {
			PercentFull[x][y] = .5;
		    }
		    else
		    {
			PercentFull[x][y] = RandZeroToOne();
			printf("x: %d  Y: %d  Per: %f\n",x,y,PercentFull[x][y]);
		    }
		    AllBeach[x][y] = 'n';
		    CellDepth[x][y] = - LandHeight;
		}
		else if (x < InitBeach - InitBWidth -1)
		    /* Lagoon at depth of shelf slope  */
		{
		    PercentFull[x][y] = 0;
		    AllBeach[x][y] = 'n';
		}
		if (PercentFull[x][y] > 1)
		{
		    printf("x: %d  Y: %d  Per: %f\n",x,y,PercentFull[x][y]);
		    PauseRun(x,y,-1);
		}
		Age[x][y] = 0;		
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
		PercentFull[x][y] = 1.0;
		AllBeach[x][y] = 'y';
	    } 	
	}


	/* PercentFull Top */

	for (y = PYstart -1; y <= PYstart + PWidth +1; y++)
	{
	    PercentFull[InitBeach + PHeight + 1][y] = RandZeroToOne();
	}

	/* PercentFull Sides */

	for (x = InitBeach ; x <= InitBeach + PHeight ; x++)
	{
	    PercentFull[x][PYstart-1] = RandZeroToOne();
	    PercentFull[x][PYstart+PWidth + 1] = RandZeroToOne();
	}
    }

    else if (InitialPert == 2)
	/* Another Perturbation  - steep point */
    {
	
	x = InitBeach;

	PercentFull[x][17] = 0.8;
	PercentFull[x][18] = 1.0;
	AllBeach[x][18] = 'y';
	PercentFull[x][19] = 0.8;

	x = InitBeach + 1;

	PercentFull[x][17] = 0.6;
	PercentFull[x][18] = 1.0;
	AllBeach[x][18] = 'y';
	PercentFull[x][19] = 0.6;

	x = InitBeach + 2;

	PercentFull[x][17] = 0.2;
	PercentFull[x][18] = 1.0;
	AllBeach[x][18] = 'y';
	PercentFull[x][19] = 0.2;

	x = InitBeach + 3;

	PercentFull[x][18] = 0.3;

    }

}


void PeriodicBoundaryCopy(void)	

/* Simulates periodic boundary conditions by copying middle section to front and end of arrays */


{
    int	x,y;

    for (y = Ymax; y < 3*Ymax/2; y++)
	for (x = 0; x < Xmax; x++)
	{
	    AllBeach[x][y-Ymax] = AllBeach[x][y];
	    PercentFull[x][y-Ymax] = PercentFull[x][y];
	    Age[x][y-Ymax] = Age[x][y];
	    CellDepth[x][y-Ymax] = CellDepth[x][y];
	}
    for (y = Ymax/2; y <= Ymax; y++)
	for (x = 0; x < Xmax; x++)
	{
	    AllBeach[x][y+Ymax] = AllBeach[x][y];
	    PercentFull[x][y+Ymax] = PercentFull[x][y];
	    Age[x][y+Ymax] = Age[x][y];
	    CellDepth[x][y+Ymax] = CellDepth[x][y];
	}
	
}










void ZeroVars(void)

/* Resets all arrays recalculated at each time step to 'zero' conditions */

{

    int z;

    for (z=0; z < MaxBeachLength; z++)
    {
	X[z] = -1;	
	Y[z] = -1;		
	InShadow[z] = '?';	
	ShorelineAngle[z] = -999;
	SurroundingAngle[z] = -998;
	UpWind[z] = '?';	
	VolumeIn[z] = 0;	
	VolumeOut[z] = 0;	
    }
}

			
void ReadSandFromFile(void)	
	
/*  Reads saved output file, AllBeach[][] & PercentFull[][]	 */

{
    int x,y;
    FILE *ReadSandFile;
	
    ReadSandFile = fopen(_io.readfilename,"r");printf("CHECK READ \n");

	
    for (y = Ymax/2; y < 3*Ymax/2; y++)
    {	

	for (x=0; x<Xmax; x++)
	{
	    fscanf(ReadSandFile, " %f", &PercentFull[x][y]);

	    if (PercentFull[x][y] >= 1.0)
		AllBeach[x][y] = 'y';
	    else 
		AllBeach[x][y] = 'n';
	}
    }

    for (y = Ymax/2; y < 3*Ymax/2; y++)
	for (x=0; x<Xmax; x++)
	    fscanf(ReadSandFile, " %f", &CellDepth[x][y]);

				
    if (SaveAge)
	
	for (y = Ymax/2; y < 3*Ymax/2; y++)
	{
	    for (x=0; x<Xmax; x++)
	    {
		fscanf(ReadSandFile, " %d", &Age[x][y]);
	    }
	}

    /*PrintLocalConds(5,5,-1);*/
    fclose(ReadSandFile);
    printf("file read!");

    PeriodicBoundaryCopy();

}


void SaveSandToFile(void)

/*  Saves current AllBeach[][] and PercentFull[][] data arrays to file 		*/
/*  Save file name will add extension '.' and the CurrentTimeStep		*/

{
    int	x,y;
    char savename[40];
    FILE *SaveSandFile;

    printf("\n saving \n ");

    sprintf(savename, "%s.%d", _io.savefilename, CurrentTimeStep);
    printf( "Saving as: %s 		", savename );	


    SaveSandFile = fopen(savename, "w");
    if (!SaveSandFile)
    {
	printf("problem opening output file\n");
	exit(1);
    }
	
	
    for (y= Ymax/2; y< 3*Ymax/2; y++)
	for (x=0; x<Xmax; x++)
	    fprintf(SaveSandFile, " %f", PercentFull[x][y]);

    for (y= Ymax/2; y< 3*Ymax/2; y++)
	for (x=0; x<Xmax; x++)
	    fprintf(SaveSandFile, " %f", CellDepth[x][y]);

    if (SaveAge)
	for (y=Ymax/2; y< 3*Ymax/2; y++)
	    for (x=0; x<Xmax; x++)
		fprintf(SaveSandFile, " %d", Age[x][y]);
		

    fclose(SaveSandFile);
    printf("--- regular file saved! ----\n\n");

}


#define SAVE_LINE_NAME "lineout"
///* SaveLineToFile */char            savelinename[24] = "lineout";

void SaveLineToFile(void)

/*  Saves data line of shoreline position rather than entire array 		*/
/*  Main concern is to have only one data point at each alongshore location	*/
/*  Save file name will add extension '.' and the CurrentTimeStep		*/

{

    int	y,x,xtop,i;
    float   xsave;
    char savename[40];
    char savelinename[24] = SAVE_LINE_NAME;
    FILE *SaveSandFile;

    printf("\n saving \n ");

    sprintf(savename, "%s%d", savelinename, CurrentTimeStep);
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
	while(AllBeach[x][y] == 'n')
	{
	    x -= 1;
	}

	/* if on side of shape, need to average */
	if (PercentFull[x+2][y] > 0)
	{
	    xtop = x+1;
	    while(PercentFull[xtop][y] > 0)
	    {
		xtop +=1;
	    } 

	    xsave = x;

	    for (i=x+1; i<xtop ; i++)
		xsave += PercentFull[i][y];

	}
	/* otherwise Regular Beach Condition */
	else
	{
	    xsave = x + PercentFull[x+1][y];
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
	for (i=0; i <= TotalBeachCells; i++)
	    if ((X[i]==x) && (Y[i]==y))
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
	    printf("  	%f", CellDepth[i][j]);
	    if (CellDepth[i][j] == DepthShoreface)
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
	    printf("	%c", AllBeach[i][j]);
	}
	printf("\n");
    }

    printf("\n");

    for (i= x+2 ; i > x-3 ; i--)
    {
	for (j = y-2 ; j < y+3 ; j++)
	{
	    printf("	%2.5f",PercentFull[i][j]);
	}
	printf("\n");
    }		

    printf("\n");

    printf(" %d  ", in );

    if (isee>=0)
    {
	for (k = in-3; k <in+4; k++)
	{
	    printf("  	%2d: %2d,%2d", k , X[k] , Y[k]);
	}
	printf("\n\n\n");

	printf("Wave Angle:	%f\n\n",WaveAngle*radtodeg);
	printf("i		%d		%d		!%d		%d		%d\n", 
	       in-2, in-1,in,in+1,in+2);
	printf("Shadow		%c		%c		%c		%c		%c\n", 
	       InShadow[in-2], InShadow[in-1],InShadow[in], InShadow[in+1], InShadow[in+2]);
	printf("Upwind		%c		%c		%c		%c		%c\n", 
	       UpWind[in-2],  UpWind[in-1],UpWind[in],  UpWind[in+1],  UpWind[in+2]);
	printf("Angle		%2.2f		%2.2f		%2.2f		%2.2f		%2.2f\n",
	       ShorelineAngle[in-2]*radtodeg,ShorelineAngle[in-1]*radtodeg, ShorelineAngle[in]*radtodeg,
	       ShorelineAngle[in+1]*radtodeg,ShorelineAngle[in+2]*radtodeg);
	printf("SurrAngle	%2.2f		%2.2f		%2.2f		%2.2f		%2.2f\n",
	       SurroundingAngle[in-2]*radtodeg, SurroundingAngle[in-1]*radtodeg,  SurroundingAngle[in]*radtodeg,
	       SurroundingAngle[in+1]*radtodeg, SurroundingAngle[in+2]*radtodeg);
	printf("Vol In 		%2.2f		%2.2f		%2.2f		%2.2f		%2.2f\n",
	       VolumeIn[in-2], VolumeIn[in-1],VolumeIn[in],VolumeIn[in+1],VolumeIn[in+2]);
	printf("Vol Out		%2.2f		%2.2f		%2.2f		%2.2f		%2.2f\n",
	       VolumeOut[in-2], VolumeOut[in-1], VolumeOut[in],VolumeOut[in+1],VolumeOut[in+2]);
	printf("Diff		%2.2f		%2.2f		%2.2f		%2.2f		%2.2f\n",
	       VolumeIn[in-2]-VolumeOut[in-2], VolumeIn[in-1]-VolumeOut[in-1], VolumeIn[in]-VolumeOut[in],
	       VolumeIn[in+1]-VolumeOut[in+1],VolumeIn[in+2]-VolumeOut[in+2]);
	printf("Frac Diff	%2.3f		%2.3f		%2.3f		%2.3f		%2.3f\n",
	       (VolumeIn[in-2]-VolumeOut[in-2])/vol, (VolumeIn[in-1]-VolumeOut[in-1])/vol,
	       (VolumeIn[in]-VolumeOut[in])/vol, (VolumeIn[in+1]-VolumeOut[in+1])/vol,
	       (VolumeIn[in+2]-VolumeOut[in+2])/vol);

    }

    printf("\n");





}


void PauseRun(int x, int y, int in)

/* Pauses run intil the 'q' key is pressed 	*/
/* Can Print or Plot Out Useful info		*/
	

{
		
    int xsee=1,ysee=-1,isee=-1,i;

    printf("\nPaused x: %d  y: %d Time: %d\n",x,y,CurrentTimeStep);

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
	    if (PercentFull[x][y] == 0)
	    {
		Age[x][y] = CurrentTimeStep%AgeMax;
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
 
	    AgeFactorRed = (float)((Age[x][y])%AgeShadeSpacing)/AgeShadeSpacing;
	    AgeFactorGreen = (float)((Age[x][y]+AgeShadeSpacing/3)%AgeShadeSpacing)/AgeShadeSpacing;
	    AgeFactorBlue = (float)((Age[x][y]+2*AgeShadeSpacing/3)%AgeShadeSpacing)/AgeShadeSpacing;

	    if ((PercentFull[x][y] > 0) && (AllBeach[x][y] == 'n'))
	    {
		Red =((((235 - 100 *(AgeFactorRed))-backRed)* PercentFull[x][y] )+backRed)/255.0;
		Green =((((235 - 95 * (AgeFactorGreen))-backGreen) * PercentFull[x][y])+backGreen)/255.0;
		Blue =((((210 - 150 * AgeFactorBlue)-backBlue)* PercentFull[x][y])+backBlue)/255.0;

	    }
	    else if (AllBeach[x][y] == 'y')
	    {
		Red =((((235 - 100 *(AgeFactorRed))-backRed)* PercentFull[x][y] )+backRed)/255.0;
		Green =((((235 - 95 * (AgeFactorGreen))-backGreen) * PercentFull[x][y])+backGreen)/255.0;
		Blue =((((210 - 150 * AgeFactorBlue)-backBlue)* PercentFull[x][y])+backBlue)/255.0;
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

    while (AllBeach[x][y] == 'y')
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

	while (AllBeach[x][y] == 'y')
	{
		x += 1;
	}

	PercentFull[x][y] += SedRate;

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
		if (debug10a) printf("L  ");
	}
	else
	{
		sweepsign = 0;
		if (debug10a) printf("R  ");
	}

	OWflag = 0;
	for (i=1; i < TotalBeachCells-1 ; i++)
	{
		if (sweepsign == 1)
			ii = i;
		else
			ii = TotalBeachCells-1-i;

		/* To do test shoreline should be facing seaward 					*/
		/* don't worry about shadow here, as overwash is not set to a time scale with AST 	*/

		if ((fabs(SurroundingAngle[ii]) < (OverwashLimit/radtodeg))  && (InShadow[ii] == 'n'))
		{
			CheckOverwash(ii);
		}	

	}

	
	/*if (OWflag) PauseRun(1,1,-1);*/

}


	

void CheckOverwash(int icheck)

	/* New 1/04 ADA - Step back pixelwise in direction of Surrounding Angle to check needage 	*/
	/* If too short, calls DoOverwash, which will move some sediment				*/
	/* Uses AllBeach[][] and PercentFull[][] (can be changed when DoOVerwash is called		*/
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
		debug10a = 1;
	else
		debug10a = 0;*/

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

	if (SurroundingAngle[icheck] > 0)
		ysign = 1;
	else
		ysign = -1;
		
		if (debug10a) printf("\nI: %d------------- Surr: %f  %f Slope: %f sign: %d \n",
		 icheck, SurroundingAngle[icheck],SurroundingAngle[icheck]*radtodeg,slope, ysign); 
	

	if (AllBeach[X[icheck]-1][Y[icheck]] == 'y' || ((AllBeach[X[icheck]][Y[icheck]-1] == 'y') && 
		(AllBeach[X[icheck]][Y[icheck]+1] == 'y')) )
	/* 'regular condition' */
	/* plus 'stuck in the middle' situation (unlikely scenario)*/
	{
		xin = X[icheck] + PercentFull[X[icheck]][Y[icheck]];
		yin = Y[icheck] + 0.5;
	}
	else if (AllBeach[X[icheck]][Y[icheck]-1] == 'y')
	/* on right side */
	{
		xin = X[icheck] + 0.5;
		yin = Y[icheck] + PercentFull[X[icheck]][Y[icheck]];
		if (debug10a) printf("-- Right xin: %f  yin: %f\n",xin,yin);
	}
	else if (AllBeach[X[icheck]][Y[icheck]+1] == 'y')
	/* on left side */
	{
		xin = X[icheck] + 0.5;
		yin = Y[icheck] + 1.0 - PercentFull[X[icheck]][Y[icheck]];
		if (debug10a) printf("-- Left xin: %f  yin: %f\n",xin,yin);
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
	
			if (debug10a) printf("x: %f  y: %f  X:%d  Y: %d  Yd: %f  DistD: %f Xs: %f DistS: %f\n",
			x,y,NextXInt,NextYInt, Ydown,DistanceDown,Xside,DistanceSide); 

		if (DistanceDown < DistanceSide)
		/* next cell is the down cell */
		{
			x = NextXInt;
			y = Ydown;
			xtest = NextXInt-1;	
			ytest = floor(y);
				/*if (debug10a) printf(" down ");*/
		}
		else
		/* next cell is the side cell */
		{
			x = Xside;
			y = NextYInt;
			xtest = floor(x);
			ytest = y + (ysign-1)/2;
				/*if (debug10a) printf(" side ");*/
		}
		
		/*if ((debug10a) && (DoGraphics == 'y'))PutPixel(ytest*CELL_PIXEL_SIZE,xtest*CELL_PIXEL_SIZE,0,0,200);*/

		checkdistance = Raise(((x - xin)*(x - xin) +  (y - yin)*(y - yin)),.5) * CellWidth;	
		if (AllBeach[xtest][ytest] == 'y')
			AllBeachFlag = 1;

		if (debug10a) printf("	x: %f  y: %f  xtest: %d ytest: %d check: %f\n\n",x,y,xtest,ytest,checkdistance);

		if ((AllBeach[xtest][ytest] == 'n') && (AllBeachFlag) && !(((X[icheck]-xtest) > 1) || (abs(ytest - Y[icheck]) > 1)))
		/* if passed through an allbeach and a neighboring partial cell, jump out, only bad things follow */
		{
			return;
		}		
		
		if((AllBeach[xtest][ytest] == 'n') && (AllBeachFlag) && (xtest < X[icheck]) &&
			(((X[icheck]-xtest) > 1) || (abs(ytest - Y[icheck]) > 1)))
		/* Looking for shore cells, but don't want immediate neighbors, and go backwards */
		/* Also mush pass though an allbeach cell along the way */
		{
		
			if (AllBeach[xtest+1][ytest] == 'y')
			/* 'regular condition' - UNDERNEATH, here */
			{
				xint = (xtest + 1 - PercentFull[xtest][ytest]);
				yint = yin + (xin - xint)*ysign * slope;

				if ((yint > ytest+1.0) || (yint < ytest))
				/* This cell isn't actually an overwash cell */
				{
					measwidth = CritBWidth;
					if (debug10a) printf("-- Regunder Cancelled  xin: %2.2f  yin: %2.2f xt:%d yt: %d xint: %f yint: %f sl: %2.2fMMeas: %3.2f\n",
					xin,yin,xtest,ytest,xint,yint,slope,measwidth);
				}
				else
				{
					measwidth = CellWidth * Raise((xint - xin)*(xint - xin)+ (yint - yin)*(yint - yin),0.5);

					if (debug10a) printf("-- Regunder Over  xin: %2.2f  yin: %2.2f xt:%d yt: %d xint: %f yint: %f sl: %2.2fMeas: %3.2f\n",
					xin,yin,xtest,ytest,xint,yint,slope,measwidth);
				}
			}
			else if (AllBeach[xtest][ytest-1] == 'y')
			/* on right side */
			{
				yint = (ytest + PercentFull[xtest][ytest]);
				xint = xin - fabs(yin - yint)/ slope;

				if (xint < xtest)
				/* This cell isn't actually an overwash cell */
				{
					measwidth = CritBWidth;

					if (debug10a) printf("-- Right Cancelled  xin: %2.2f  yin: %2.2f xt:%d yt: %d xint: %f yint: %f sl: %2.2fMMeas: %3.2f\n",
					xin,yin,xtest,ytest,xint,yint,slope,measwidth);
				}
				else
				{
					measwidth = CellWidth * Raise((xint - xin)*(xint - xin)+ (yint - yin)*(yint - yin),0.5);

					if (debug10a) printf("-- Right Over  xin: %2.2f  yin: %2.2f xt:%d yt: %d xint: %f yint: %f sl: %2.2fMMeas: %3.2f\n",
					xin,yin,xtest,ytest,xint,yint,slope,measwidth);
				}
			}
			else if (AllBeach[xtest][ytest+1] == 'y')
			/* on left side */
			{
				yint = (ytest + 1 - PercentFull[xtest][ytest]);
				xint = xin - fabs(yin - yint)/ slope;

				if (xint < xtest)
				/* This cell isn't actually an overwash cell */
				{
					measwidth = CritBWidth;
					
					if (debug10a) printf("-- Left cancelled  xin: %2.2f  yin: %2.2f xt:%d yt: %d xint: %f yint: %f sl: %2.2fMMeas: %3.2f\n",
					xin,yin,xtest,ytest,xint,yint,slope,measwidth);
				}
				else
				{
					measwidth = CellWidth * Raise((xint - xin)*(xint - xin)+ (yint - yin)*(yint - yin),0.5);

					if (debug10a) printf("-- Left Over  xin: %2.2f  yin: %2.2f xt:%d yt: %d xint: %f yint: %f sl: %2.2fMMeas: %3.2f\n",
					xin,yin,xtest,ytest,xint,yint,slope,measwidth);
				}
			}
			else if (AllBeach[xtest-1][ytest] == 'y')
			/* 'regular condition' */
			/* plus 'stuck in the middle' situation */
			{
				xint = (xtest  + PercentFull[xtest][ytest]);
				yint = yin + (xin - xint)*ysign * slope;

				if ((yint > ytest+1.0) || (yint < ytest))
				/* This cell isn't actually an overwash cell */
				{
					measwidth = CritBWidth;
					if (debug10a) printf("-- RegularODD Cancelled  xin: %2.2f  yin: %2.2f xt:%d yt: %d xint: %f yint: %f Meas: %3.2f\n",
					xin,yin,xtest,ytest,xint,yint,measwidth);
				}
				else
				{
					measwidth = CellWidth * Raise((xint - xin)*(xint - xin)+ (yint - yin)*(yint - yin),0.5);

					if (debug10a) printf("-- RegularODD Over  xin: %2.2f  yin: %2.2f xt:%d yt: %d xint: %f yint: %f Meas: %3.2f\n",
					xin,yin,xtest,ytest,xint,yint,measwidth);
					/*PauseRun(xtest,ytest,icheck);*/
				}
			}
			else if (PercentFull[xtest][ytest] > 0)
			/* uh oh - not good situation, no allbeach on sides */
			/* assume this is an empty cell, */
			{
				xint = x;
				yint = y;
			
				measwidth = CellWidth * Raise((xint - xin)*(xint - xin)+ (yint - yin)*(yint - yin),0.5);

				printf("-- Some Odd Over  xin: %2.2f  yin: %2.2f xt:%d yt: %d xint: %f yint: %f Meas: %3.2f Ang: %f Abs: %f\n",
				xin,yin,xtest,ytest,xint,yint,measwidth, SurroundingAngle[icheck]*radtodeg,fabs(SurroundingAngle[icheck])*radtodeg);
				/*PauseRun(xtest,ytest,icheck);*/
			}
			else
			/* empty cell - oughta fill er up  - fill max barrier width*/
			{
				xint = x;
				yint = y;
				measwidth = CritBWidth - CellWidth;
				
				printf("-- Empty Odd Over  xin: %2.2f  yin: %2.2f xt:%d yt: %d xint: %f yint: %f Meas: %3.2f Ang: %f Abs: %f\n",
				xin,yin,xtest,ytest,xint,yint,measwidth, SurroundingAngle[icheck]*radtodeg,fabs(SurroundingAngle[icheck])*radtodeg);
				/*PauseRun(xtest,ytest,icheck); */
			}

			checkdistance = measwidth;
			
			if (measwidth < CritBWidth)
			{
				DoOverwash(X[icheck],Y[icheck],xtest,ytest,xint,yint,measwidth,icheck);
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
	/*  will change and use PercentFull[][] and AllBeach [][]				*/

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


	if (debug10b) printf("** Overwash From X: %d  Y: %d  To: X: %d Y: %d Width: %f \n"
		, xfrom, yfrom,xto,yto,widthin );
	if (debug10b) printf("DepthBB: %f  BBNeed: %f DelShore: %f  DelBB: %f\n",
			DepthBB, BBneed,delShore,delBB );
	/*if (DepthBB == DepthShoreface) PauseRun(xto,yto,-1);*/

#ifdef WITH_OPENGL
	if (debug10b && ( DO_GRAPHICS == 'y'))	
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

	PercentFull[xto][yto] += delBB;
	PercentFull[xfrom][yfrom] -= delShore;

	if (PercentFull[xto][yto] > 1)
	{
		OopsImFull(xto,yto);
	}
	if (PercentFull[xfrom][yfrom] < 0)
	{
		OopsImEmpty(xfrom,yfrom);
	}

	if (debug10b) PauseRun(xto,yto,-1);

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
		Depth = CellDepth[xdepth][yin];
	
		while ((Depth < 0) && (xdepth > 0))
		{
			Depth = CellDepth[xdepth][yin];
			/*printf("-- Overwash depth problem - Here = %f Next = %f",CellDepth[xdepth][yto],CellDepth[xdepth-1][yto] );
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
				
		if (SurroundingAngle[ishore] == 0.0)
		{
			/* unlikely, but make sure no div by zero */
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

		BackFlag = 0;
		if (SurroundingAngle[ishore] > 0)
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
		
				if (debug10b) printf("x: %f  y: %f  X:%d  Y: %d  Yd: %f  DistD: %f Xs: %f DistS: %f\n",
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
		
			if (PercentFull[xtest][ytest] > 0)
				BackFlag = 1;
		}
	
		/* Try to find the i for the cell found */
		/* If you have a better idea how to do this, go ahead */

		i = 2;
		FoundFlag = 0;
		
		while ((i < TotalBeachCells-1) && !(FoundFlag))
		{
			if ((X[i] == xtest) && (Y[i] == ytest))
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
			if (debug10b) printf("\nbackbarrier out of bounds: xin: %d yin: %d xbi: %d ybi: %d xinf: %f yinf: %f Per: %f Dist:  Depth: %f\n",
			xin, yin, xtest, ytest, xinfl, yinfl,PercentFull[xtest][ytest], Depth);
			/*PauseRun(xin,yin,-1);*/
		}	
		else
		{
			BBDistance =  Raise(((xinfl - xtest-PercentFull[xtest][ytest])*
			(xinfl - xtest-PercentFull[xtest][ytest])) + 
			((yinfl - ytest - 0.5)*(yinfl - ytest - 0.5)),.5);

 			if(!FoundFlag)
			/* The backbarrier intersection isn't on the shoreline */
			/* Assume 1/2 of the length applies to this case */
			{
				Depth = BBDistance/2 * ShorefaceSlope * CellWidth;
				if (debug10b) printf("\nNot Found backi: %d bx: %d by: %d Depth:%f",
					Backi,xtest,ytest,Depth);
			}
			else
			/* Use the fancy geometry thing */
			{
				AngleUsed = 0;
				for (j = -1; j <2 ; j ++)
				{
					AngleUsed += SurroundingAngle[Backi+j];
				}
				AngleUsed = AngleUsed/5;

				if (fabs(AngleUsed) > M_PI/4.0)
				{
					AngleUsed = M_PI/4.0;
					if (debug10b) printf("Big Angle");
					/*PauseRun(X[Backi],Y[Backi],Backi);*/
				}

				AngleSin = sin(M_PI/2.0 - fabs(SurroundingAngle[ishore] + AngleUsed));

				Depth = BBDistance * AngleSin / (1 + AngleSin);
		
		if (debug10b) printf("\nBack Angle backi: %d bx: %d by: %d BackA: %f AngU: %f Asin: %f L/2: %f Depth:%f",
			Backi,X[Backi],Y[Backi],SurroundingAngle[ishore]*radtodeg,AngleUsed*radtodeg,AngleSin,
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
			if (debug10b) printf("\nOverwash Depth2: xin: %d yin: %d xbi: %d ybi: %d xinf: %f yinf: %f Per: %f Dist: %f  Depth: %f\n",
			xin, yin, xtest, ytest, xinfl, yinfl,PercentFull[xtest][ytest],BBDistance, Depth);	
		return Depth;
	} 

	printf("OWDepth all broken");
	PauseRun(xin,yin,-1);
	return DepthShoreface;
}

