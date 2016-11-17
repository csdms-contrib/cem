#ifndef CEM_CONSTS_INCLUDED
#define CEM_CONSTS_INCLUDED

#if defined(__cplusplus)
extern "C" {
#endif

// IMPORTANT -- specify wave transformation routine
#undef WITH_SWAN

// Aspect Parameters
#define CELL_WIDTH (100.)           // size of cells (meters)
#define X_MAX (50)                  // number of cells in x (cross-shore) direction
#define Y_MAX (200)                 // number of cells in y (longshore) direction


// Beach definition parameters
#define MAX_BEACH_LENGTH (8 * Y_MAX)// maximum length of arrays that contain beach data at each time step
#define SHELF_SLOPE (0.001)         // slope of continental shelf
#define SHOREFACE_SLOPE (0.01)      // for now, linear slope of shoreface
#define DEPTH_SHOREFACE (10.)       // minimum depth of shoreface due to wave action (meters)
#define INIT_BEACH (20)             // cell where intial conditions changes from beach to ocean
#define INIT_ROCK (0)               // cell where initial conditions change from beach to rock LMV
#define INITIAL_DEPTH (10.)         // theoretical depth in meters of continental shelf at x = INIT_BEACH
#define LAND_HEIGHT (1.)


// Overwash Parameters
#define CRIT_BEACH_WIDTH (350.0)    // Overwash - width barrier maintains due to overwash (m)
#define INIT_BEACH_WIDTH (4)        // Overwash - initial minimum width of barrier (Cells)
#define OVERWASH_TYPE (1)           // Overwash - 0 = use depth array, 1 = use geometric rule
#define OVERWASH_MIN_DEPTH (1.0)    // Overwash - littlest overwash of all
#define MAX_OVER (0.01)             // Maximum overwash step size (enforced at back barrier)
#define OVERWASH_LIMT (60)          // Don't do over wash if the angle is > 60 degrees


//  Run Control Parameters
#define WAVE_BREAK_DEPTH (0.2)        //wave breaking threshold (wave height/depth) (SWAN)

#define FIND_CELL_ERROR (5)           // if we run off of array, how far over do we try again?
#define SHADOW_STEP_DISTANCE (0.2)    // step size for shadow cell checking
#define AGE_MAX (1000000)             // Maximum 'age' of cells - loops back to zero
#define AGE_UPDATE (10)               // Time space for updating age of non-beach cells
// 0: normal (columns/blocks),
// 1: wiggly,
// 2: one block,
// 3: all sand, no  rock
#define INITIAL_CONDITION_TYPE (0)
#define INITIAL_PERT (0)              // Start with a bump? if =1-->square pert, if =2-->pointy pert
#define DIFFUSIVE_HUMP (0)            // Shoreline is sinusoidal LMV
#define INITIAL_SMOOTH (0)            // Smooth starting conditions
#define INITIAL_SMOOTH_ROCK (1)       // Smooth rock interface starting conditions LMV
#define SEED (1973)                    // random seed:  control value = 1 completely random = -999
// Sinks
#define HAVE_SINKS (0)                // include sediment sinks in model run?
#define NUM_SINKS (6)                 // Number of sinks to include
#define COLUMN_SINKS (0)              // disregard x vector and consider sink to span entire column
#define SINKINESS (0.5)               // fractional value determines chance that sediment in the sink is deleted
// weathering
#define NUMBER_CHUNK (9)              // Number of chunks of rock in alongshore direction with different weather rates LMV
#define NORMAL_WEATHERING_RATE (0.2)  // Baseline rock retreat rate, m/yr   PWL
#define TIME_TO_SWEEP_FULL_BEACH (50) // Spacing of full beach sweep for each rock cell i=0:TotalBeachCells LMV
#define LOOK_DIST (10)                // For short rock to beach sweep, +- number of beach cells to look at LMV
#define BIG_DISTANCE_TO_BEACH (1000.0 * CELL_WIDTH) // Used to find minimum distance to beach in rock to beach sweep LMV
#define NO_WEATHERING (5.0)           // weathering of rock only occurs below this amt of sed cover (vert equiv in meters) LMV
// Alongshore length of a chunk of fast or slow weathering rock LMV
// (used to be Y_MAX/NumberChunk)
// #define CHUNK_LENGTH (100)

// Turn abrasion on? If yes, then rock weathering is maximized at a
// given amount of sediment cover (wcrit). Otherwise, use existing
// exponential weathering.  PWL
#define ABRASION (1)
#define W_CRIT (20)                   // Sediment cover that maximizes rock weathering. PWL
#define N (4.)                        // How much to maximize weathering rate above bare-rock rate at Wcrit? PWL


// I/O Parameters
#define START_SAVING_AT (0)           // time step to begin saving files
#define SAVE_SPACING (365)            // space between saved files
// char savefilename[2048] = "CEM";
// char readfilename[2048] = "CEM_3285.out";
#define WAVE_IN (0)                   // Input Wave Distribution file: no = 0, binned file = 1, angle/period/height file = 2
// char readwavename[2048] = "In_WaveData.dat";
#define INITIALIZE_FILE (0)           // use a file to initialise run? Over rides setup data
// char readcontrolname[2048] = "In_CEM_init.dat";
#define METADATA (0)                // Create a metadata file?
// char metasavename[2048] = "Metadata.out";
#define WAVE_DATA (0)               // Create a wavedata file?
// char wavesavename[2048] = "Wavedata.out";
#define SAVE_FILE (2)               // 1 = line output, 2 = array output
#define SAVE_AGE (1)                // Save/update age of cells?
#define SAVE_LINE (1)               // Save line instead of whole array?
#define PROMPT_START (0)            // ask prompts for file names, etc?
#define SCREEN_TEXT_SPACING (30)    // Spacing of writing to screen in time steps
#define NO_PAUSE_RUN (1)            // Disable PauseRun subroutine


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

// Universal Constants
#define PI (3.1415927)
#define GRAVITY (9.80665)
#define RAD_TO_DEG (180. / PI) // transform rads to degrees

// Booleans
#ifndef TRUE
# define TRUE 1
#endif
#ifndef FALSE
# define FALSE 0
#endif


#if defined(__cplusplus)
}
#endif

#endif
