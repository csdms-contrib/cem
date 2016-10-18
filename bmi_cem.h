#ifndef BMI_CEM_H_INCLUDED
#define BMI_CEM_H_INCLUDED

#if defined(__cplusplus)
extern "C" {
#endif

#include "bmi.h"

BMI_Model * register_bmi_cem(BMI_Model *model);


typedef struct {
  double time;
} CemModel;


#if defined(__cplusplus)
}
#endif

#endif
