#ifndef BMI_CEM_H_INCLUDED
#define BMI_CEM_H_INCLUDED

#if defined(__cplusplus)
extern "C" {
#endif

#include "bmi.h"
#include "cem_model.h"

Bmi * register_bmi_cem(Bmi *model);
CemModel* new_cem_model(void);

#if defined(__cplusplus)
}
#endif

#endif
