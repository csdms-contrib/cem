#ifndef CEM_CONSTS_INCLUDED
#define CEM_CONSTS_INCLUDED

#if defined(__cplusplus)
extern "C" {
#endif

// IMPORTANT -- specify wave transformation routine
#undef WITH_SWAN

// Aspect Parameters
// size of cells (meters)
#define CELL_WIDTH (100.)
// number of cells in x (cross-shore) direction
#define X_MAX (50)
// number of cells in y (longshore) direction
#define Y_MAX (200)

// maximum length of arrays that contain beach data at each time step
#define MaxBeachLength (8 * Y_MAX)

// slope of continental shelf
#define SHELF_SLOPE (0.001)
// for now, linear slope of shoreface
#define SHOREFACE_SLOPE (0.01)
// minimum depth of shoreface due to wave action (meters)
#define DEPTH_SHOREFACE (10.)
// cell where intial conditions changes from beach to ocean
#define INIT_BEACH (20)
// cell where initial conditions change from beach to rock LMV
#define INIT_ROCK (5)
// theoretical depth in meters of continental shelf at x = INIT_BEACH
#define INITIAL_DEPTH (10.)
#define LAND_HEIGHT (1.)


// if we run off of array, how far over do we try again?
#define FindCellError (5)
// step size for shadow cell checking
#define ShadowStepDistance (0.2)

// Overwash Parameters
#define CritBWidth                                                          \
  350.0 /* Overwash - width barrier maintains due to overwash (m) important \
           scaling parameter! */
#define InitBWidth 4 /* Overwash - initial minimum width of barrier (Cells) */
#define OWType 1 /* Overwash - 0 = use depth array, 1 = use geometric rule */
#define OWMinDepth 1.0 /* Overwash - littlest overwash of all */

// SWAN
// Define wave breaking threshold (wave height/depth) to pull out necessary
// metrics from SWAN
#define WaveBreakDepth (0.2)

// Maximum overwash step size (enforced at back barrier)
#define MaxOver (0.01)
// Don't do over wash if the angle is > 60 degrees
#define OverwashLimit (60)

#ifndef TRUE
# define TRUE 1
#endif
#ifndef FALSE
# define FALSE 0
#endif

//  Run Control Parameters

// include sediment sinks in model run?
#define HAVE_SINKS (0)
// fractional value determines chance that sediment in the sink is deleted
#define SINKINESS (0.5)
// Number of chunks of rock in alongshore direction with different
// weathering rates LMV
#define NUMBER_CHUNK (9)
// Baseline rock retreat rate, m/yr   PWL
#define NORMAL_WEATHERING_RATE (0.2)

// Turn abrasion on? If yes, then rock weathering is maximized at a
// given amount of sediment cover (wcrit). Otherwise, use existing
// exponential weathering.  PWL
#define ABRASION ('y')
// Sediment cover that maximizes rock weathering. PWL
#define W_CRIT (20)
// How much to maximize weathering rate above bare-rock rate at Wcrit? PWL
#define N (4.)

// 0: normal (columns/blocks),
// 1: wiggly,
// 2: one block,
// 3: all sand, no  rock
#define INITIAL_CONDITION_TYPE (0)
// random seed:  control value = 1 completely random = -999
#define SEED (1973)
// time step to begin saving files
#define START_SAVING_AT (0)
// space between saved files
#define SAVE_SPACING (365)

// char savefilename[2048] = "CEM";
// char readfilename[2048] = "CEM_3285.out";
// Input Wave Distribution file: no = 0, binned file = 1,
// angle/period/height file = 2
#define WAVE_IN (0)
// char readwavename[2048] = "In_WaveData.dat";
// use a file to initialise run? Over rides setup data
#define INITIALIZE_FILE ('n')
// char readcontrolname[2048] = "In_CEM_init.dat";
// Create a metadata file?
#define METADATA ('n')
// char metasavename[2048] = "Metadata.out";
// Create a wavedata file?
#define WAVE_DATA ('n')
// char wavesavename[2048] = "Wavedata.out";

// 1 = line output,
// 2 = array output
#define SAVE_FILE (2)
// Save/update age of cells?
#define SAVE_AGE (1)
// Save line instead of whole array?
#define SAVE_LINE (1)
// ask prompts for file names, etc?
#define PROMPT_START ('n')
// Spacing of writing to screen in time steps
#define SCREEN_TEXT_SPACING (30)
// Spacing of full beach sweep for each rock cell i=0:TotalBeachCells LMV
#define TIME_TO_SWEEP_FULL_BEACH (50)
// For short rock to beach sweep, +- number of beach cells to look at LMV
#define LOOK_DIST (10)
/* Used to find minimum distance to beach in rock to beach sweep LMV */
#define BIG_DISTANCE_TO_BEACH (1000.0 * CELL_WIDTH)
// weathering of rock only occurs below this amt of sed cover
// (vert equiv in meters) LMV
#define NO_WEATHERING (5.0)
// Disbale PauseRun subroutine
#define NO_PAUSE_RUN (1)
// Start with a bump? if =1-->square pert, if =2-->pointy pert
#define INITIAL_PERT (0)
// Shoreline is sinusoidal LMV
#define DIFFUSIVE_HUMP (0)
// Smooth starting conditions
#define INITIAL_SMOOTH (0)
// Smooth rock interface starting conditions LMV
#define INITIAL_SMOOTH_ROCK (1)
// Alongshore length of a chunk of fast or slow weathering rock LMV
// (used to be Y_MAX/NumberChunk)
// #define CHUNK_LENGTH (100)

/* De-bugging Parameters */

#define debug0 0   /* misc. printf's--don't leave them in */
#define debug1 0   /* Find Next Beach Cell */
#define debug1a 0  /* Find Next Rock Cell LMV */
#define debug1b 0  /* More find next cell LMV */
#define debug1c 0  /* More find next rock cell LMV */
#define debug2 0   /* Shadow Routine */
#define debug2a 0  /* In-Depth Shadow Results */
#define debug2b 0  /* More shadow results */
#define debug3 0   /* Determine Angles */
#define debug4 0   /* Upwind/Downwind */
#define debug5 0   /* Sediment Transport Decisions */
#define debug6 0   /* Sediment Trans Calculations */
#define debug6a 0  /* Volume In/Out/AcrossBorder LMV */
#define debug7 0   /* Transport Sweep (move sediment) */
#define debug7a 0  /* Slope Calcs */
#define debug8 0   /* Full/Empty */
#define debug80 0  /* Temp Full */
#define debug8a 0  /* More Full/Empty LMV */
#define debug9 0   /* FixBeach */
#define debug9a 0  /* More FixBeach LMV */
#define debug10 0  /* Distance to Beach LMV */
#define debug11 0  /* Min Distance to Beach & Closest Beach LMV */
#define debug12 0  /* Weathering Rate LMV */
#define debug12a 0 /* Percent Fines lost LMV */
#define debug13 0  /* Sed Flow LMV */
#define debug14 0  /* Actual Volume Across Border & PFS LMV */
#define debug15 0  /* Erosion Rate per Year LMV */
#define debug16 0  /* Total Percent Full LMV */
#define debug17 0  /* graphics stuff */
#define debug18 0  /* reading wave data AB */
#define debugNearest                                                  \
  0 /* graphs a pixel at nearest beach cell to a given rock cell (for \
       weathering) */
#define debugtopo                                                            \
  0 /* debug topography...not really useful at this point -- wait until land \
       is DEM! */

/* Universal Constants */
#define PI (3.1415927)
#define GRAVITY (9.80665)
#define RAD_TO_DEG (180. / PI) /* transform rads to degrees */

#define NUM_SINKS (6)
// if this is turned on, X vector is disregarded and a sink is
// considered to span a whole column
#define COLUMN_SINKS (0)

// Maximum 'age' of cells - loops back to zero
#define AGE_MAX (1000000)
// Time space for updating age of non-beach cells
#define AGE_UPDATE (10)

#if defined(__cplusplus)
}
#endif

#endif
