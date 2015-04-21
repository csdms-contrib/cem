#if !defined (CEM_MODEL_H)
#define CEM_MODEL_H

#if defined(__cplusplus)
extern "C" {
#endif


#define DEFAULT_CELL_WIDTH (100.0) /**< size of cells (meters) */
/*#define MaxBeachLength  (8*Ymax) */ /**< maximum length of arrays that contain beach data at each time step */
#define TimeStep     (0.2)  /**< days - reflects rate of sediment transport per
                             time step */

typedef struct
{
  int use_sed_flux;  /**< Use SedFlux rather than SedRate */
  double SedFlux;  /**< Sediment flux in kg/s. */
  double SedRate;  /**< Sedimentation rate as percent per time step. */
  double angle_highness;  /**< Fraction of high-angle waves. */
  double angle_asymmetry;  /**< Fraction of waves coming from the left. */
  double wave_height;  /**< Height of incoming waves in meters. */
  double wave_period;  /**< Period of incoming waves in seconds. */
  double shoreface_slope;  /**< Gradient of the shoreface. */
  double shelf_slope;  /**< Gradient of the shelf. */
  double shoreface_depth;  /**< Water depth of the shoreface in meters. */

  int nx;  /**< Number of cells in x (cross-shore) direction */
  int ny;  /**< Number of cells in y (long-shore) direction */
  int max_beach_len;  /**< Max number of cells that can make up the coastline */

  int n_rivers;
  double *river_flux;
  int *river_x_ind;
  int *river_y_ind;
  double *river_x;
  double *river_y;

  int stream_spot;

  double cell_width;  /**< Size of cell in meters */

   /** Input/output file names. */
  char *savefilename;  /**< Name of save file. */
  char *readfilename;  /**< Namve of file to read input from. */

   /** Overall Shoreface Configuration Arrays - Data file information
       This grids will be of size (nx, ny).
    */

  char **AllBeach;  /**< Flag indicating of cell is entirely beach */
  double **PercentFull;  /**< Fractional amount of shore cell full of
                                       sediment */
  int **Age;  /**< Age since cell was deposited */
  double **CellDepth;  /**< Depth array (m) (ADA 6/3) */
  double **InitDepth;  /**< Save initial depths (m) (EWHH 2010/8/11) */

   /** Computational Arrays (determined for each time step) */
  int *X;  /**< X Position of ith beach element */
  int *Y;  /**< Y Position of ith beach element */
  char *InShadow;                /**< Is ith beach element in shadow? */
  double *ShorelineAngle;  /**< Angle between cell and right (z+1)
                                            neighbor */
  double *SurroundingAngle; /**< Cell-orientated angle based upon
                                             left and right neighbor */
  char *UpWind;  /**< Upwind or downwind condition used to
                                   calculate sediment transport */
  double *VolumeIn;   /**< Sediment volume into ith beach element */
  double *VolumeOut;  /**< Sediment volume out of ith beach element */

   /** Miscellaneous State Variables */
  int CurrentTimeStep;  /**< Time step of current calculation */

  int NextX;  /**< used to iterate FindNextCell in global array - */
  int NextY;

  int TotalBeachCells;  /**< Number of cells describing beach at particular iteration */
  int ShadowXMax;  /**< used to determine maximum extent of beach cells */

  int external_waves;
  double WaveAngle;  /**< wave angle for current time step */

  int FindStart;  /**< Used to tell FindBeach at what Y value to start looking */

  char FellOffArray;  /**< Flag used to determine if accidentally went off array */

  double MassInitial;  /**< For conservation of mass calcs */
  double MassCurrent;

  int NumWaveBins;     /**< For Input Wave - number of bins */
  double WaveMax[36];   /**< Max Angle for specific bin */
  double WaveProb[36];  /**< Probability of Certain Bin */

   /** Graphics variables. */
  double xcellwidth;
  double ycellwidth;
  int xplotoff;
  int yplotoff;

  /* char state[256]; */
  char *state;
}
CemModel;

int cem_initialize (const char *config_file, CemModel **handle);

int deltas_get_nx (CemModel * model);
int deltas_get_ny (CemModel * model);
double deltas_get_dx (CemModel * model);
double deltas_get_dy (CemModel * model);
double deltas_get_current_time (CemModel * model);
double deltas_get_end_time (CemModel * model);
double deltas_get_time_step (CemModel * model);

CemModel * deltas_set_shoreface_slope (CemModel * model, double shoreface_slope);
CemModel * deltas_set_shelf_slope (CemModel * model, double shelf_slope);
CemModel * deltas_set_shoreface_depth (CemModel * model, double shoreface_depth);
CemModel * deltas_set_save_file (CemModel * model, const char *name);

void deltas_avulsion (CemModel * model, double *qs, double river_flux);

#if defined(__cplusplus)
}
#endif

#endif
