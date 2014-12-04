#ifndef CEM_GLOBALS_INCLUDED
#define CEM_GLOBALS_INCLUDED

#if defined(__cplusplus)
extern "C"
{
#endif

#include "consts.h"

extern float CellDepth[Xmax][2 * Ymax];      /* Depth array */
extern double CliffHeightSlow;  /* Cliff height above sea level for slow weathering rock PWL */
extern double CliffHeightFast;   /* Cliff height above sea level for fast weathering rock PWL */

/* Overall Shoreface Configuration Arrays - Data file information */
extern char AllBeach[Xmax][2 * Ymax];        /* Flag indicating of cell is entirely beach */
extern char AllRock[Xmax][2 * Ymax]; /* Flag indicating if cell is entirely rock LMV */
extern double PercentFullSand[Xmax][2 * Ymax];       /* Fractional amount of cell full of sediment LMV */
extern double PercentFullRock[Xmax][2 * Ymax];       /* Fractional amount of a cell full of rock LMV */
extern char TypeOfRock[Xmax][2 * Ymax];      /* Array to control weathering rates of rock along the beach LMV */
extern int Age[Xmax][2 * Ymax];      /* Age since cell was deposited */
extern double Topography[Xmax][2 * Ymax];    /* Holds cliff heights -- will change through time, eventually... */

extern int SinkY[6]; /* a sink is a cell that is routinely emptied (if it is on the beach) */
extern int SinkX[6];

extern FILE *SaveSandFile;
extern FILE *WaveFileOutput;
extern FILE *InitMetaFile;
extern FILE *ReadSandFile;
extern FILE *ReadWaveFile;
extern FILE *ReadRealWaveFile;
extern FILE *ReadControlFile;

/* Computational Arrays (determined for each time step)  -- will eventually be in a structure for BMI */

extern int X[MaxBeachLength];        /* X Position of ith beach element */
extern int Y[MaxBeachLength];        /* Y Position of ith beach element */
extern int XBehind[MaxBeachLength];  /* Cell that is "behind" X[i] LMV */
extern int YBehind[MaxBeachLength];  /* Cell that is "behind" Y[i] LMV */
extern int XRock[MaxBeachLength];    /* X Position of ith rock element LMV */
extern int YRock[MaxBeachLength];    /* Y Position of ith rock element LMV */
extern int XRockBehind[MaxBeachLength];      /* Cell that is "behind" XRock[i] LMV */
extern int YRockBehind[MaxBeachLength];      /* Cell that is "behind" YRock[i] LMV */
extern char InShadow[MaxBeachLength];        /* Is ith beach element in shadow? */
extern float ShorelineAngle[MaxBeachLength]; /* Angle between cell and right (z+1) neighbor  */
extern float SurroundingAngle[MaxBeachLength];       /* Angle between left and right neighbors */
extern char UpWind[MaxBeachLength];  /* Upwind or downwind condition used to calculate sediment transport */
extern float VolumeIn[MaxBeachLength];       /* Sediment volume into ith beach element */
extern float VolumeOut[MaxBeachLength];      /* Sediment volume out of ith beach element */
extern float VolumeAcrossBorder[MaxBeachLength];     /* Sediment volume across border of ith beach element in m^3/day LMV */
/* amount sediment needed, not necessarily amount a cell gets */
extern float ActualVolumeAcross[MaxBeachLength];     /* Sediment volume that actually gets passed across border LMV */
/* amount sed is limited by volumes across border upstream and downstream */
extern char DirectionAcrossBorder[MaxBeachLength];   /* Flag to indicate if transport across border is L or R LMV */
extern char FlowThroughCell[MaxBeachLength]; /* Is flow through ith cell Left, Right, Convergent, or Divergent LMV */
extern float DistanceToBeach[MaxBeachLength];        /* Distance in meters from rock to beach LMV */
extern float MinDistToBeach[MaxBeachLength]; /* From a rock cell j, min distance (in meters) to closest beach LMV */
extern int ClosestBeach[MaxBeachLength];     /* i position of closest rock to beach LMV */
extern float AmountWeathered[MaxBeachLength];        /* Amount of rock weathered from rock cell j LMV */

/* Miscellaneous Global Variables -- also will be included in the BMI structure */

extern int NextX;                    /* Global variables used to iterate FindNextCell in global array - */
extern int NextY;                    /*      would've used pointer but wouldn't work */
extern int BehindX;
extern int BehindY;
extern int BehindRockX;
extern int BehindRockY;
extern int NextRockX;                /* Global variables used to iterate FindNextRockCell in global array LMV */
extern int NextRockY;
extern int TotalBeachCells;          /* Number of cells describing beach at particular iteration */
extern int TotalRockCells;           /* Number of cells describing rock at an iteration LMV */
extern int ShadowXMax;               /* used to determine maximum extent of beach cells */
extern float WaveAngle;              /* wave angle for current time step */
extern int FindStart;                /* Used to tell FindBeach at what Y value to start looking */
extern int FindRockStart;            /* Used to tell FindRock at what Y value to start looking LMV */
extern char FellOffArray;            /* Flag used to determine if accidentally went off array */
extern char FellOffRockArray;        /* Flag used to determine if accidentally went off rock array LMV */
extern float MassInitial;            /* For conservation of mass calcs */
extern float MassCurrent;            /* " */
extern int device;
extern short button;
extern long buttonback;
extern int NumWaveBins;              /* For Input Wave - number of bins */
extern float WaveMax[36];            /* Max Angle for specific bin */
extern float WaveProb[36];           /* Probability of Certain Bin */
extern double WaveAngleIn;
extern double WaveHeightIn;
extern double WavePeriodIn;
extern double ControlFileIn[20];     /* Initialisation data read from file */

extern char StartFromFile;     /* start from saved file? */

extern int CurrentTimeStep;      /* Time step of current calculation */
extern double StopAfter;     /* Stop after what number of time steps */

extern double Period;         /* seconds */
extern double OffShoreWvHt;    /* meters */

extern double Asym;           /* fractional portion of waves coming from positive (left) direction */
extern double Highness;       /* .5 = even dist, < .5 high angle domination. NOTE Highness actually determines lowness! */
extern double Duration;         /* Number of time steps calculations loop at same wave angle */
extern double TimeStep;          /* days - reflects rate of sediment transport per time step */

extern double SlowWeatherCoeff;  /* Weathering rate of slow rock e.g. 0.5 = slow weathering, 1/2 as fast LMV */
extern double FastWeatherCoeff;   /* Weathering rate of fast rock e.g. 2 = fast weathering, 2 times faster than normal LMV */
extern double PercentFineFast; /* Percent of fast weathering rock lost because it is too fine to stay in nearshore LMV */
extern double PercentFineSlow; /* Percent of slow weathering rock lost because it is too fine to stay in nearshore LMV */
extern double ErosionRatePerYear;      /* Amount of erosion to TotalBeachCells in m/year LMV */

extern double coastrotation;   /* Angle (deg) used to align coastline with real wave climate */
extern double waveheightchange; /* Factor applied to the input wave height; 1 = no change, 0.5 = half wave height, and 2 = double wave height */
extern double waveperiodchange; /* Factor applied to the input wave period; 1 = no change, 0.5 = half wave period, and 2 = double wave period */

extern int OWflag;               /* A debugging flag for overwash routines */

#if defined(__cplusplus)
}
#endif

#endif
