#include <stdlib.h>
#include <string.h>

#include "conditions.h"
#include "consts.h"
#include "utils.h"

extern double** cell_depth;
extern char** AllBeach;
extern char** AllRock;
extern double** PercentFullSand;
extern double** PercentFullRock;
extern char** type_of_rock;
extern double** topography;

extern double kCliffHeightSlow;
extern double kCliffHeightFast;

static int* rock_line;
static int* beach_line;

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

void InitBeach(void);
void InitDiffusive(void);
int InitWiggly(void);
void InitPert(int pert_type);
void InitRockBlocks(void);
void InitRockStripes(void);
void ApplySpecialConditions(void);
void CopyRow(int start_row, int stop_row);

void InitConds (void) 
{
    int start = Y_MAX / 2;
    int end = Y_MAX + start;

    rock_line = malloc(sizeof(int) * (end-start));
    beach_line = malloc(sizeof(int) * (end-start));

    /* Base beach conditions */
    InitBeach();
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
            InitPert(1);
            break;
        case POINTED_PERTURBATION :
            InitPert(2);
            break;
    }

    ApplySpecialConditions();

    /* apply rock type patterns */
    switch (INITIAL_ROCK_TYPE)
    {
        case STRIPE :
            InitRockStripes();
            break;
        case BLOCK :
            InitRockBlocks();
            break;
    }

    free(rock_line);
    free(beach_line);
}

/**
 * Fills in a flat beach with in
 */
void InitBeach(void)
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
            PercentFullRock[0][col] = 1.0;
            PercentFullSand[0][col] = 0.0;
            AllRock[0][col] = 'y';
            AllBeach[0][col] = 'y';
            topography[0][col] = kCliffHeightSlow;
            /* Fill in rock/sand interface row */
            float percent_rock = INITIAL_SMOOTH_ROCK ? 0.2 : RandZeroToOne();
            PercentFullRock[INIT_ROCK][col] = percent_rock;
            PercentFullSand[INIT_ROCK][col] = 1 - percent_rock;
            AllRock[INIT_ROCK][col] = 'n';
            AllBeach[INIT_ROCK][col] = 'y';
            topography[INIT_ROCK][col] = kCliffHeightSlow;
        }
        
        int beach_start = INIT_ROCK+1;
        /* Fill in first row of sand */
        PercentFullRock[beach_start][col] = 0.0;
        PercentFullSand[beach_start][col] = 1.0;
        AllRock[beach_start][col] = 'n';
        AllBeach[beach_start][col] = 'y';
        topography[beach_start][col] = 0;
        /* Fill in sand/ocean interface row */
        PercentFullRock[beach_start][col] = 0.0;
        PercentFullSand[INIT_BEACH][col] = INITIAL_SMOOTH ? 0.5 : RandZeroToOne();
        AllRock[INIT_BEACH][col] = 'n';
        AllBeach[INIT_BEACH][col] = 'n';
        topography[INIT_BEACH][col] = 0;
        
        int ocean_start = INIT_BEACH + 1;
        /* Fill in first row of ocean */
        PercentFullRock[beach_start][col] = 0.0;
        PercentFullSand[ocean_start][col] = 0.0;
        AllRock[ocean_start][col] = 'n';
        AllBeach[ocean_start][col] = 'n';
        topography[ocean_start][col] = 0;
    }

    /* Copy rock into all rows in between 0 and INIT_ROCK*/
    CopyRow(0, INIT_ROCK);
    /* Copy sand into all rows in between INIT_ROCK and INIT_BEACH*/
    CopyRow(INIT_ROCK + 1, INIT_BEACH);
    /* Copy ocean into all remaining rows */
    CopyRow(INIT_BEACH + 1, X_MAX);
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
            if (rock_line[y] < Min || rock_line[y] > Max) return -1;
        }
    }

    printf("done wiggly init\n");
    return 0;
}

/* Andrew's initial bump */
void InitPert(int pert_type)
{
    int x, y;
    int PWidth = 5;
    int PHeight = 3;
    int PYstart = 25;

    switch(pert_type)
    {
        /* Square perturbation */
        case 1 :
            /* Fill AllBeach areas */
            for (x = INIT_BEACH; x <= INIT_BEACH + PHeight; x++)
            {
                for (y = PYstart; y <= PYstart + PWidth; y++)
                {
                    PercentFullSand[x][y] = 1.0;
                    AllBeach[x][y] = 'y';
                }
            }

            /* PercentFull Top */
            for (y = PYstart - 1; y <= PYstart + PWidth + 1; y++)
            {
                PercentFullSand[INIT_BEACH + PHeight + 1][y] = INITIAL_SMOOTH ? 0.5 : RandZeroToOne();
            }

            /* PercentFull Sides */
            for (x = INIT_BEACH; x <= INIT_BEACH + PHeight; x++)
            {
                PercentFullSand[x][PYstart - 1] =  INITIAL_SMOOTH ? 0.5 : RandZeroToOne();
                PercentFullSand[x][PYstart + PWidth + 1] =  INITIAL_SMOOTH ? 0.5 : RandZeroToOne();
            }
            break;
        /* Another Perturbation  - steep point */
        case 2 :
            x = INIT_BEACH;

            PercentFullSand[x][17] = 0.8;
            PercentFullSand[x][18] = 1.0;
            AllBeach[x][18] = 'y';
            PercentFullSand[x][19] = 0.8;

            x = INIT_BEACH + 1;

            PercentFullSand[x][17] = 0.6;
            PercentFullSand[x][18] = 1.0;
            AllBeach[x][18] = 'y';
            PercentFullSand[x][19] = 0.6;

            x = INIT_BEACH + 2;

            PercentFullSand[x][17] = 0.2;
            PercentFullSand[x][18] = 1.0;
            AllBeach[x][18] = 'y';
            PercentFullSand[x][19] = 0.2;

            x = INIT_BEACH + 3;

            PercentFullSand[x][18] = 0.3;
            break;
    }
}

/**
 * Shift rock and beach lines if necessary
 */
void ApplySpecialConditions()
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
            PercentFullRock[rockline][col] = INITIAL_SMOOTH_ROCK ? 0.2 : RandZeroToOne();
            AllRock[rockline][col] = 'n';

            for(row = rockline + 1; row <= INIT_ROCK; row ++)
            {
                PercentFullRock[row][col] = 0.0;
                AllRock[row][col] = 'n';
                topography[row][col] = 0;
            }

            for(row = INIT_ROCK; row < rockline; row ++)
            {
                PercentFullRock[row][col] = 1.0;
                AllRock[row][col] = 'y';
                topography[row][col] = kCliffHeightSlow;
            }
        }

        int sandline = beach_line[col-start];
        if(sandline != INIT_BEACH)
        {
            /* adjust sand line */
            PercentFullSand[sandline][col] = INITIAL_SMOOTH ? 0.5 : RandZeroToOne();
            AllBeach[sandline][col] = 'n';

            for(row = sandline + 1; row <= INIT_BEACH; row ++)
            {
                PercentFullSand[row][col] = 0.0;
                AllBeach[row][col] = 'n';
            }

            for(row = INIT_BEACH; row < sandline; row ++)
            {
                PercentFullSand[row][col] = 1.0;
                AllBeach[row][col] = 'y';
            }
        }
    }
}

/**
 * Initialize blocky rock pattern
 */
void InitRockBlocks(void)
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

/**
 * Initialize stiped rock pattern
 */
void InitRockStripes()
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
 * Copy one row to a block of rows
 */
void CopyRow(int start_row, int stop_row)
{
    int start = Y_MAX / 2;
    int row;
    for(row = start_row + 1; row < stop_row; row++)
    {
        memcpy(PercentFullRock[row] + start, PercentFullRock[start_row] + start, Y_MAX*sizeof(double));
        memcpy(PercentFullSand[row] + start, PercentFullSand[start_row] + start, Y_MAX*sizeof(double));
        memcpy(AllRock[row] + start, AllRock[start_row] + start, Y_MAX*sizeof(char));
        memcpy(AllBeach[row] + start, AllBeach[start_row] + start, Y_MAX*sizeof(char));
        memcpy(topography[row] + start, topography[start_row] + start, Y_MAX*sizeof(double));
    }
}
