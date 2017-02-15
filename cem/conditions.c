#include <stdlib.h>
#include <string.h>

#include "conditions.h"
#include "consts.h"
#include "utils.h"

static int *rock_line;
static int *beach_line;

enum ConditionType
{
    NORMAL = 0,
    DIFFUSIVE = 1,
    WIGGLY = 2,
    SQUARE_PERTURBATION = 3,
    POINTED_PERTURBATION = 4
};

enum RockPattern
{
    NONE = 0,
    STRIPE = 1,
    BLOCK = 2
};

void InitBeach(double **percent_full_sand, double **percent_full_rock, char **all_beach, char **all_rock, double **topography);
void InitDiffusive(void);
int InitWiggly(void);
void InitPert(int pert_type, double **percent_full_sand, char **all_beach);
void InitRockBlocks(char **type_of_rock);
void InitRockStripes(char **type_of_rock, double **topography);
void ApplySpecialConditions(double **percent_full_sand, double **percent_full_rock, char **all_beach, char **all_rock, double **topography);
void CopyData(int start_row, int stop_row, double **percent_full_sand, double **percent_full_rock, char **all_beach, char **all_rock, double **topography);

void InitConds (double **cell_depth, char **all_beach, char **all_rock, double **percent_full_sand, double **percent_full_rock, char **type_of_rock, double **topography)
{

    int start = Y_MAX / 2;
    int end = Y_MAX + start;

    rock_line = malloc(sizeof(int) * (end-start));
    beach_line = malloc(sizeof(int) * (end-start));

    /* Base beach conditions */
    InitBeach(percent_full_sand, percent_full_rock, all_beach, all_rock, topography);
    /* Add special characteristics */
    switch (INITIAL_CONDITION_TYPE)
    {
        case DIFFUSIVE :
            InitDiffusive();
            break;
        case WIGGLY :
            while (InitWiggly() == -1)
                printf("init: amplitude too great, trying again...\n");
            break;
        case SQUARE_PERTURBATION :
            InitPert(1, percent_full_sand, all_beach);
            break;
        case POINTED_PERTURBATION :
            InitPert(2, percent_full_sand, all_beach);
            break;
    }

    ApplySpecialConditions(percent_full_sand, percent_full_rock, all_beach, all_rock, topography);

    /* apply rock type patterns */
    switch (INITIAL_ROCK_TYPE)
    {
        case STRIPE :
            InitRockStripes(type_of_rock, topography);
            break;
        case BLOCK :
            InitRockBlocks(type_of_rock);
            break;
    }

    free(rock_line);
    free(beach_line);
}

/**
 * Fills in a flat beach with in
 */
void InitBeach(double **percent_full_sand, double **percent_full_rock, char **all_beach, char **all_rock, double **topography)
{
    int col;
    int start = Y_MAX / 2;
    int end = Y_MAX + start;

    /* Fill in conditions on transition rows */
    for(col = start; col < end; col ++)
    {
        rock_line[col - start] = INIT_ROCK;
        beach_line[col - start] = INIT_BEACH;

        if( INIT_ROCK > 0){
            /* Fill in first row of rock */
            percent_full_rock[0][col] = 1.0;
            percent_full_sand[0][col] = 0.0;
            all_rock[0][col] = 'y';
            all_beach[0][col] = 'y';
            topography[0][col] = kCliffHeightSlow;
            /* Fill in rock/sand interface row */
            float percent_rock = INITIAL_SMOOTH_ROCK ? 0.2 : RandZeroToOne();
            percent_full_rock[INIT_ROCK][col] = percent_rock;
            percent_full_sand[INIT_ROCK][col] = 1 - percent_rock;
            all_rock[INIT_ROCK][col] = 'n';
            all_beach[INIT_ROCK][col] = 'y';
            topography[INIT_ROCK][col] = kCliffHeightSlow;
        }
        
        int beach_start = INIT_ROCK+1;
        /* Fill in first row of sand */
        percent_full_rock[beach_start][col] = 0.0;
        percent_full_sand[beach_start][col] = 1.0;
        all_rock[beach_start][col] = 'n';
        all_beach[beach_start][col] = 'y';
        topography[beach_start][col] = 0;
        /* Fill in sand/ocean interface row */
        percent_full_rock[beach_start][col] = 0.0;
        percent_full_sand[INIT_BEACH][col] = INITIAL_SMOOTH ? 0.5 : RandZeroToOne();
        all_rock[INIT_BEACH][col] = 'n';
        all_beach[INIT_BEACH][col] = 'n';
        topography[INIT_BEACH][col] = 0;
        
        int ocean_start = INIT_BEACH + 1;
        /* Fill in first row of ocean */
        percent_full_rock[beach_start][col] = 0.0;
        percent_full_sand[ocean_start][col] = 0.0;
        all_rock[ocean_start][col] = 'n';
        all_beach[ocean_start][col] = 'n';
        topography[ocean_start][col] = 0;
    }

    /* Copy rock into all rows in between 0 and INIT_ROCK*/
    CopyData(0, INIT_ROCK, percent_full_sand, percent_full_rock, all_beach, all_rock, topography);
    /* Copy sand into all rows in between INIT_ROCK and INIT_BEACH*/
    CopyData(INIT_ROCK + 1, INIT_BEACH, percent_full_sand, percent_full_rock, all_beach, all_rock, topography);
    /* Copy ocean into all remaining rows */
    CopyData(INIT_BEACH + 1, X_MAX, percent_full_sand, percent_full_rock, all_beach, all_rock, topography);
}

/**
 * Initialize diffusive hump (cosine curve)
 */
void InitDiffusive(void)
{
    const double amp = 10.; // Amplitude of cos curve */
    int col;
    int start = Y_MAX / 2;
    int end = Y_MAX + start;

    for (col = start; col < end; col++)
    {
      rock_line[col-start] = (INIT_ROCK + (-amp * cos((2 * PI * col) / Y_MAX)) - (INIT_BEACH - INIT_ROCK - 1));
      beach_line[col-start] = (INIT_BEACH + (-amp * cos((2 * PI * col) / Y_MAX)));
    }
}

// TODO
/**
 * Initialize sinusoidal coastline
 */
int InitWiggly(void)
{
    int start = Y_MAX / 2;
    int end = Y_MAX + start;

    int y, curve;
    int NumCurves = 4;
    int curveFactor = 3; /* decrease fundamental wavelength/increase frequency (by a factor)--make it 1 for wavelength = Y_MAX */
    float Amp,
        AmpMax = INIT_ROCK / 6.0; /* INIT_ROCK/4.0; maximum amplitude of the component sin waves */
    int Min = 0.10 * X_MAX,
        Max = 0.95 * X_MAX; /* determines acceptable max amplitude of resulting curve */

    printf("init wiggly\n");

    for (curve = 1; curve <= NumCurves; curve++)
    {
        Amp = RandZeroToOne() * AmpMax;
        printf("for curve %d, amplitude %f\n", curve, Amp);
        for (y = start; y < end; y++)
        {
            rock_line[y] = rock_line[y] + Amp * sin(y * (curve * curveFactor) * 2.0 * PI / Y_MAX);
            beach_line[y] = rock_line[y];
            if (rock_line[y] < Min || rock_line[y] > Max) return -1;
        }
    }

    printf("done wiggly init\n");
    return 0;
}

/* Andrew's initial bump */
void InitPert(int pert_type, double **percent_full_sand, char **all_beach)
{
    int x, y;
    int PWidth = 5;
    int PHeight = 3;
    int PYstart = 25;

    int yPeak = 18;

    switch(pert_type)
    {
        /* Square perturbation */
        case 1 :
            for (y = PYstart; y <= PYstart + PWidth; y++)
            {
                beach_line[y] = INIT_BEACH + PHeight;
            }

            /* PercentFull Sides */
            for (x = INIT_BEACH; x <= INIT_BEACH + PHeight; x++)
            {
                percent_full_sand[x][PYstart - 1] =  INITIAL_SMOOTH ? 0.5 : RandZeroToOne();
                percent_full_sand[x][PYstart + PWidth + 1] =  INITIAL_SMOOTH ? 0.5 : RandZeroToOne();
            }
            break;

        /* Another Perturbation  - steep point */
        case 2 :
            x = INIT_BEACH;

            percent_full_sand[x][yPeak-1] = 0.8;
            percent_full_sand[x][yPeak] = 1.0;
            all_beach[x][yPeak] = 'y';
            percent_full_sand[x][yPeak+1] = 0.8;

            x++;

            percent_full_sand[x][yPeak - 1] = 0.6;
            percent_full_sand[x][yPeak] = 1.0;
            all_beach[x][yPeak] = 'y';
            percent_full_sand[x][yPeak + 1] = 0.6;

            x++;;

            percent_full_sand[x][yPeak - 1] = 0.2;
            percent_full_sand[x][yPeak] = 1.0;
            all_beach[x][yPeak] = 'y';
            percent_full_sand[x][yPeak + 1] = 0.2;

            x++;;

            percent_full_sand[x][yPeak] = 0.3;
            break;
    }
}

/**
 * Shift rock and beach lines if necessary
 */
void ApplySpecialConditions(double **percent_full_sand, double **percent_full_rock, char **all_beach, char **all_rock, double **topography)
{
    int col;
    int start = Y_MAX / 2;
    int end = Y_MAX + start;

    for(col = start; col < end; col++)
    {
        int row;
        int rockline = rock_line[col-start];

        if (rockline != INIT_ROCK)
        {
            /* adjust rock line */
            percent_full_rock[rockline][col] = INITIAL_SMOOTH_ROCK ? 0.2 : RandZeroToOne();
            all_rock[rockline][col] = 'n';

            for(row = rockline + 1; row <= INIT_ROCK; row ++)
            {
                percent_full_rock[row][col] = 0.0;
                all_rock[row][col] = 'n';
                topography[row][col] = 0;
            }

            for(row = INIT_ROCK; row < rockline; row ++)
            {
                percent_full_rock[row][col] = 1.0;
                all_rock[row][col] = 'y';
                topography[row][col] = kCliffHeightSlow;
            }
        }

        int sandline = beach_line[col-start];
        if(sandline != INIT_BEACH)
        {
            /* adjust sand line */
            percent_full_sand[sandline][col] = INITIAL_SMOOTH ? 0.5 : RandZeroToOne();
            all_beach[sandline][col] = 'n';

            for(row = sandline + 1; row <= INIT_BEACH; row ++)
            {
                percent_full_sand[row][col] = 0.0;
                all_beach[row][col] = 'n';
            }

            for(row = INIT_BEACH; row < sandline; row ++)
            {
                percent_full_sand[row][col] = 1.0;
                all_beach[row][col] = 'y';
            }
        }
    }
}

/**
 * Initialize blocky rock pattern
 */
void InitRockBlocks(char **type_of_rock)
{
    int start = Y_MAX / 2;
    int end = Y_MAX + start;
    int x, y, n, blockHeight = 4, blockWidth = 60, NumBlocks = 0; /* makes NumBlocks evenly spaced blocks */

    printf("init block\n");

    int col;
    for (col = start; col < end; col ++)
    {
        int row;
        for (row = 0; row < rock_line[col-start]; row ++)
        {
            type_of_rock[row][col] = 'f';
        }
    }

    /* the regularly-spaced blocks */
    for (n = 1; n <= NumBlocks; n++)
    {
        for (x = INIT_ROCK; x >= INIT_ROCK - blockHeight; x--)
        {
            for (y = Y_MAX / 2 + n * Y_MAX / (NumBlocks + 1);
                y < Y_MAX / 2 + n * Y_MAX / (NumBlocks + 1) + blockWidth &&
                y < 2 * Y_MAX;
                y++) { /* block starts @ upper left corner */
                type_of_rock[x][y] = 's';
            }
        }
    }
    /* now here's space for ad hoc blocks */
    /* ad hoc block 1 */
    for (x = INIT_ROCK; x >= INIT_ROCK - blockHeight; x--)
    {
        for (y = Y_MAX / 2; y < Y_MAX / 2 + 10; y++)
        {
            type_of_rock[x][y] = 's';
        }
    }
    /* ad hoc block 2 */
    for (x = INIT_ROCK; x >= INIT_ROCK - blockHeight; x--)
    {
        for (y = 3 * Y_MAX / 2 - 90; y < 3 * Y_MAX / 2; y++)
        {
            type_of_rock[x][y] = 's';
        }
    }

    printf("done block init\n");
}

// TODO
/**
 * Initialize stiped rock pattern
 */
void InitRockStripes(char **type_of_rock, double **topography)
{
    if(INIT_ROCK > 3)
    {
        int y;
        for (y = 0; y < 2 *Y_MAX; y++)
        {
            int x;
            for (x = 0; x < INIT_ROCK - 3; x++)
            {
                type_of_rock[x][y] = 'f';
                topography[x][y] = kCliffHeightFast;
            }
        }
        set_rock_blocks(type_of_rock + INIT_ROCK - 3, topography + INIT_ROCK - 3, X_MAX - (INIT_ROCK - 3), 2 * Y_MAX, NUMBER_CHUNK);
    }
}

/**
 * Copy one row pf beach data to a block of rows
 */
void CopyData(int start_row, int stop_row, double **percent_full_sand, double **percent_full_rock, char **all_beach, char **all_rock, double **topography)
{
    int start = Y_MAX / 2;
    int row;
    for(row = start_row + 1; row < stop_row; row++)
    {
        memcpy(percent_full_rock[row] + start, percent_full_rock[start_row] + start, Y_MAX*sizeof(double));
        memcpy(percent_full_sand[row] + start, percent_full_sand[start_row] + start, Y_MAX*sizeof(double));
        memcpy(all_rock[row] + start, all_rock[start_row] + start, Y_MAX*sizeof(char));
        memcpy(all_beach[row] + start, all_beach[start_row] + start, Y_MAX*sizeof(char));
        memcpy(topography[row] + start, topography[start_row] + start, Y_MAX*sizeof(double));
    }
}
