#ifndef BMI_WAVES_H_INCLUDED
#define BMI_WAVES_H_INCLUDED

#if defined(__cplusplus)
extern "C" {
#endif

#include "bmi.h"
#include "waves_model.h"

BMI_Model * register_bmi_waves(BMI_Model *model);
WavesModel* waves_new(void);

#if defined(__cplusplus)
}
#endif

#endif
