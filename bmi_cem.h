#ifndef BMI_CEM_H_INCLUDED
#define BMI_CEM_H_INCLUDED

#if defined(__cplusplus)
extern "C" {
#endif

#include "bmi.h"
#include "cem_model.h"

BMI_Model * register_bmi_cem(BMI_Model *model);
CemModel* new_cem_model(void);

#if defined(__cplusplus)
}
#endif

#endif
