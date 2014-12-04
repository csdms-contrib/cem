#ifndef CEM_GLOBALS_INCLUDED
#define CEM_GLOBALS_INCLUDED

#if defined(__cplusplus)
extern "C" {
#endif

#include "consts.h"

	float						CellDepth[Xmax][2*Ymax]; /* Depth array */
	double CliffHeightSlow = 30;		/* Cliff height above sea level for slow weathering rock PWL */
	double CliffHeightFast = 0;			/* Cliff height above sea level for fast weathering rock PWL */

	/* Overall Shoreface Configuration Arrays - Data file information */
	char   AllBeach[Xmax][2*Ymax];			/* Flag indicating of cell is entirely beach */
	char   AllRock[Xmax][2*Ymax];			/* Flag indicating if cell is entirely rock LMV */
	double PercentFullSand[Xmax][2*Ymax];	/* Fractional amount of cell full of sediment LMV */
	double PercentFullRock[Xmax][2*Ymax];	/* Fractional amount of a cell full of rock LMV */
	char   TypeOfRock[Xmax][2*Ymax];		/* Array to control weathering rates of rock along the beach LMV */
	int    Age[Xmax][2*Ymax];				/* Age since cell was deposited */
	double Topography[Xmax][2*Ymax];		/* Holds cliff heights -- will change through time, eventually...*/	
	
	int 	SinkY[] = {158, 361, 158, 361, 158, 361};		/* a sink is a cell that is routinely emptied (if it is on the beach) */
	int 	SinkX[] = {190, 190, 189, 189, 191, 191};
	
	FILE	*SaveSandFile;
	FILE    *WaveFileOutput;
	FILE    *InitMetaFile;
	FILE	*ReadSandFile;
	FILE	*ReadWaveFile;
	FILE    *ReadRealWaveFile;
	FILE    *ReadControlFile;
	
	
	/* Computational Arrays (determined for each time step)  -- will eventually be in a structure for BMI */
	
	int	    X[MaxBeachLength];		/* X Position of ith beach element */
	int	    Y[MaxBeachLength];		/* Y Position of ith beach element */
	int 	XBehind[MaxBeachLength];	/* Cell that is "behind" X[i] LMV */
	int	    YBehind[MaxBeachLength];	/* Cell that is "behind" Y[i] LMV */
	int	    XRock[MaxBeachLength];		/* X Position of ith rock element LMV */
	int	    YRock[MaxBeachLength];		/* Y Position of ith rock element LMV */
	int 	XRockBehind[MaxBeachLength];	/* Cell that is "behind" XRock[i] LMV */
	int 	YRockBehind[MaxBeachLength];	/* Cell that is "behind" YRock[i] LMV */
	char	InShadow[MaxBeachLength];	/* Is ith beach element in shadow? */
	float	ShorelineAngle[MaxBeachLength];	/* Angle between cell and right (z+1) neighbor	*/
	float   SurroundingAngle[MaxBeachLength];   	/* Angle between left and right neighbors */
	char	UpWind[MaxBeachLength];		/* Upwind or downwind condition used to calculate sediment transport */
	float	VolumeIn[MaxBeachLength];	/* Sediment volume into ith beach element*/
	float 	VolumeOut[MaxBeachLength];	/* Sediment volume out of ith beach element */
	float 	VolumeAcrossBorder[MaxBeachLength];	/* Sediment volume across border of ith beach element in m^3/day LMV*/
	/* amount sediment needed, not necessarily amount a cell gets */
	float	ActualVolumeAcross[MaxBeachLength];	/* Sediment volume that actually gets passed across border LMV */
	/* amount sed is limited by volumes across border upstream and downstream */
	char	DirectionAcrossBorder[MaxBeachLength];	/* Flag to indicate if transport across border is L or R LMV*/
	char	FlowThroughCell[MaxBeachLength];	/* Is flow through ith cell Left, Right, Convergent, or Divergent LMV*/
	float	DistanceToBeach[MaxBeachLength]; 	/* Distance in meters from rock to beach LMV*/
	float	MinDistToBeach[MaxBeachLength];	/* From a rock cell j, min distance (in meters) to closest beach LMV*/
	int     ClosestBeach[MaxBeachLength];	/* i position of closest rock to beach LMV */
	float	AmountWeathered[MaxBeachLength];	/* Amount of rock weathered from rock cell j LMV */
	
	/* Miscellaneous Global Variables -- also will be included in the BMI structure */
	
	int	    NextX;			/* Global variables used to iterate FindNextCell in global array - */
	int	    NextY;			/*	would've used pointer but wouldn't work	*/
	int 	BehindX;
	int	    BehindY;
	int 	BehindRockX;
	int	    BehindRockY;
	int 	NextRockX;		/* Global variables used to iterate FindNextRockCell in global array LMV*/
	int	    NextRockY;
	int 	TotalBeachCells;	/* Number of cells describing beach at particular iteration */
	int	    TotalRockCells; 	/* Number of cells describing rock at an iteration LMV */
	int 	ShadowXMax; 		/* used to determine maximum extent of beach cells */
	float 	WaveAngle;		/* wave angle for current time step */
	int	    FindStart;		/* Used to tell FindBeach at what Y value to start looking */
	int 	FindRockStart;		/* Used to tell FindRock at what Y value to start looking LMV */
	char	FellOffArray;		/* Flag used to determine if accidentally went off array */
	char    FellOffRockArray;	/* Flag used to determine if accidentally went off rock array LMV */
	float   MassInitial;		/* For conservation of mass calcs */
	float	MassCurrent;		/* " */
	int	    device;
	short	button;
	long	buttonback;
	int	    NumWaveBins;		/* For Input Wave - number of bins */
	float	WaveMax[36];		/* Max Angle for specific bin */
	float	WaveProb[36];		/* Probability of Certain Bin */
	double   WaveAngleIn;
	double  WaveHeightIn;
	double  WavePeriodIn;
	double ControlFileIn[20];   /* Initialisation data read from file*/

	char StartFromFile = 'n';	/* start from saved file? */

int CurrentTimeStep = 0; /* Time step of current calculation */
double StopAfter = 36500;       /* Stop after what number of time steps */

double Period = 10.0;           /* seconds */
double OffShoreWvHt = 2.0;	/* meters */

double Asym = 0.70;	/* fractional portion of waves coming from positive (left) direction */
double Highness = 0.35;	/* .5 = even dist, < .5 high angle domination. NOTE Highness actually determines lowness! */
double Duration = 1.; /* Number of time steps calculations loop at same wave angle */
double TimeStep = 1;	/* days - reflects rate of sediment transport per time step */

double SlowWeatherCoeff = 0. * NormalWeatheringRate;		/* Weathering rate of slow rock e.g. 0.5 = slow weathering, 1/2 as fast LMV */
double FastWeatherCoeff = 1 * NormalWeatheringRate;		/* Weathering rate of fast rock e.g. 2 = fast weathering, 2 times faster than normal LMV */
double PercentFineFast = 0.0;	/* Percent of fast weathering rock lost because it is too fine to stay in nearshore LMV */
double PercentFineSlow = 0.0;	/* Percent of slow weathering rock lost because it is too fine to stay in nearshore LMV */
double ErosionRatePerYear = 0.5;	/* Amount of erosion to TotalBeachCells in m/year LMV*/

double coastrotation = 40.;	/* Angle (deg) used to align coastline with real wave climate*/
double waveheightchange = 1.;	/* Factor applied to the input wave height; 1 = no change, 0.5 = half wave height, and 2 = double wave height*/
double waveperiodchange = 1.;	/* Factor applied to the input wave period; 1 = no change, 0.5 = half wave period, and 2 = double wave period*/

int OWflag = 0; /* A debugging flag for overwash routines */

#if defined(__cplusplus)
}
#endif

#endif
