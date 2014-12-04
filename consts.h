#ifndef CEM_CONSTS_INCLUDED
#define CEM_CONSTS_INCLUDED

#if defined(__cplusplus)
extern "C" {
#endif

	/* Aspect Parameters */
	#define CellWidth			100 	/* size of cells (meters) */
	#define Xmax				50	/* number of cells in x (cross-shore) direction */
	#define Ymax				200	/* number of cells in y (longshore) direction */
	#define MaxBeachLength			8*Ymax	/* maximum length of arrays that contain beach data at each time step */
	#define ShelfSlope			0.001	/* slope of continental shelf */
	#define ShorefaceSlope			0.01	/* for now, linear slope of shoreface */
	#define DepthShoreface			10	/* minimum depth of shoreface due to wave action (meters) */
	#define InitBeach			20	/* cell where intial conditions changes from beach to ocean */
	#define InitRock			5	/* cell where initial conditions change from beach to rock LMV*/
	#define InitialDepth			10	/* theoretical depth in meters of continental shelf at x = InitBeach */
	#define FindCellError			5	/* if we run off of array, how far over do we try again? */
	#define ShadowStepDistance  		0.2 	/* step size for shadow cell checking */
	
	/*Overwash Parameters*/
	#define CritBWidth			350.0   /* Overwash - width barrier maintains due to overwash (m) important scaling parameter! */
	#define InitBWidth			4       /* Overwash - initial minimum width of barrier (Cells) */
	#define OWType				1       /* Overwash - 0 = use depth array, 1 = use geometric rule */
	#define OWMinDepth			1.0		/* Overwash - littlest overwash of all */
	
	#define MaxOver				0.01    /* Maximum overwash step size (enforced at back barrier) */
	#define OverwashLimit			60      /* Don't do over wash if the angle is > 60 degrees */

	#ifndef TRUE
	#define TRUE	1
	#endif
	#ifndef FALSE
	#define FALSE 	0
	#endif
	
	/*  Run Control Parameters */
	#define HaveSinks			0	/* include sediment sinks in model run? */
	#define Sinkiness			0.5 	/* fractional value determines chance that sediment in the sink is deleted */
	#define STOP_AFTER			36500	/* Stop after what number of time steps? This is replaced by variable pushed here from BMI, 12/3/14 */
	#define NumberChunk			9	/* Number of chunks of rock in alongshore direction with different weathering rates LMV */
	#define NormalWeatheringRate 		0.2	/* Baseline rock retreat rate, m/yr   PWL */

	#define VaryCliffHeight			('y')	/* Vary cliff height for different rock weathering rates? PWL */
	#define Abrasion			('y')   /* Turn abrasion on? If yes, then rock weathering is maximized at a given amount of sediment cover (wcrit).
												Otherwise, use existing exponential weathering.  PWL */
	#define Wcrit				20	/* Sediment cover that maximizes rock weathering. PWL */
	#define N				4	/* How much to maximize weathering rate above bare-rock rate at Wcrit? PWL */
	#define Emin				0.001   /* Determines weathering rate when sediment thickness = NoWeathering. 
												Needed to calculate decay constant in WeatherRock function PWL */ 
	#define	DoGraphics			0       /* CWT Re-cast as a define, rather than the odd char DoGraphics = 1 in the original; 0 = false */
	
	#define InitCType			0	/* 0: normal (columns/blocks), 1: wiggly, 2: one block */
	#define	seed				1	/* random seed:  control value = 1 completely random = -999 */
	#define	StartSavingAt			0	/* time step to begin saving files */
	#define	SaveSpacing			365	/* space between saved files */
	#define	savefilename			"CEM"
	#define	readfilename			"CEM_3285.out"
	#define	WaveIn				0	/* Input Wave Distribution file: no = 0, binned file = 1, angle/period/height file =2 */
	#define	readwavename			"In_WaveData.dat"
	#define InitialiseFile			('n')	/* use a file to initialise run? Over rides setup data*/
	#define readcontrolname			"In_CEM_init.dat"
	#define Metadata			('n')	/*Create a metadata file?*/
	#define metasavename			"Metadata.out"
	#define Wavedata			('n')	/*Create a wavedata file?*/
	#define wavesavename			"Wavedata.out"
	
	#define	SaveFile			2       /* 1 = line output, 2 = array output*/
	#define SaveAge				1	/* Save/update age of cells? */
	#define	SaveLine			1	/* Save line instead of whole array? */
	#define	PromptStart			('n')	/* ask prompts for file names, etc? */
	#define OffArray			('n')	/* Initializing this variable for later use */
	#define	ScreenTextSpacing		30	/* Spacing of writing to screen in time steps */
	#define TimeToSweepFullBeach 		50  	/* Spacing of full beach sweep for each rock cell i=0:TotalBeachCells LMV */
	#define LookDist			10	/* For short rock to beach sweep, +- number of beach cells to look at LMV */
	#define BigDistanceToBeach  		1000.0*CellWidth  /* Used to find minimum distance to beach in rock to beach sweep LMV */
	#define NoWeathering			5.0	/* weathering of rock only occurs below this amt of sed cover (vert equiv in meters) LMV */
	#define	InteractivePlot			0
	#define	StartStop			0	/* Stop after every iteration 'Q' to move on */
	#define	InterruptRun			0	/* Allow run to be paused by pressing the 'A' key */
	#define	NoPauseRun			1	/* Disbale PauseRun subroutine */
	#define	InitialPert			0	/* Start with a bump? if =1-->square pert, if =2-->pointy pert */
	#define	DiffusiveHump			0	/* Shoreline is sinusoidal LMV	*/
	#define	InitialSmooth			0	/* Smooth starting conditions */
	#define	InitialSmoothRock		1	/* Smooth rock interface starting conditions LMV */
	#define ChunkLength			100	/* Alongshore length of a chunk of fast or slow weathering rock LMV */
								/* (used to be Ymax/NumberChunk) */
	
	/* De-bugging Parameters */
	
	#define debug0				0		/* misc. printf's--don't leave them in */
	#define	debug1				0		/* Find Next Beach Cell */
	#define debug1a				0		/* Find Next Rock Cell LMV*/
	#define	debug1b				0		/* More find next cell LMV*/
	#define	debug1c				0		/* More find next rock cell LMV*/
	#define debug2				0		/* Shadow Routine */
	#define	debug2a				0		/* In-Depth Shadow Results */
	#define	debug2b				0		/* More shadow results */
	#define debug3				0		/* Determine Angles */
	#define debug4				0		/* Upwind/Downwind */
	#define debug5				0		/* Sediment Transport Decisions*/
	#define debug6				0		/* Sediment Trans Calculations */
	#define	debug6a				0		/* Volume In/Out/AcrossBorder LMV */
	#define debug7				0		/* Transport Sweep (move sediment) */
	#define	debug7a				0		/* Slope Calcs */
	#define debug8				0		/* Full/Empty */
	#define debug80				0		/* Temp Full */
	#define debug8a				0		/* More Full/Empty LMV*/
	#define debug9				0		/* FixBeach */
	#define debug9a				0		/* More FixBeach LMV */
	#define debug10				0		/* Distance to Beach LMV*/
	#define debug11				0		/* Min Distance to Beach & Closest Beach LMV */
	#define debug12				0		/* Weathering Rate LMV */
	#define debug12a			0		/* Percent Fines lost LMV */
	#define debug13				0		/* Sed Flow LMV*/
	#define	debug14				0		/* Actual Volume Across Border & PFS LMV */
	#define debug15				0		/* Erosion Rate per Year LMV*/
	#define debug16				0		/* Total Percent Full LMV */
	#define debug17				0		/* graphics stuff */
	#define debug18				0		/* reading wave data AB*/
	#define debugNearest			0		/* graphs a pixel at nearest beach cell to a given rock cell (for weathering) */
	#define	debugtopo			0		/* debug topography...not really useful at this point -- wait until land is DEM! */
	
	/* Universal Constants */
	#define	pi        (3.1415927)
	#define g					(9.80665)
	#define	radtodeg  (180. / pi)  /* transform rads to degrees */
	
	#define NumSinks			6
	#define ColumnSinks			0		/* if this is turned on, X vector is disregarded and a sink is considered
											to span a whole column */
	
	#define	AgeMax				1000000		/* Maximum 'age' of cells - loops back to zero */
	#define	AgeUpdate			10		/* Time space for updating age of non-beach cells */
	#define	AgeShadeSpacing			10000		/* For graphics - how many time steps means back to original shade */
	#define	CellPixelSize			4		/* Size in pixes of plotted cell (a power of two, please) */

#if defined(__cplusplus)
}
#endif

#endif
