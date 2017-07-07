#include "bmi/bmi.h"
#include "bmi/bmi_cem.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char *argv[]) {
  BMI_Model *cem = (BMI_Model *)malloc(sizeof(BMI_Model));
  int err = BMI_SUCCESS;
  int status;

  if (argc > 1) {
    if (strcmp(argv[1], "--version") == 0) {
      fprintf(stdout, "The Coastal Evolution Model version 0.1.1\n");
      exit(0);
    } else if (strcmp(argv[1], "--help") == 0) {
      fprintf(stdout, "Usage: cem [--help] [--version]\n");
      exit(0);
    }
  }

  register_bmi_cem(cem);

  status = cem->initialize(NULL, &(cem->self));
  if (status == BMI_FAILURE) {
    fprintf(stdout, "FAIL.\n");
    fprintf(stderr, "Unable to initialize\n");
    return EXIT_FAILURE;
  }

  {
    char name[BMI_MAX_COMPONENT_NAME];
    cem->get_component_name(cem->self, name);
    fprintf(stdout, "%s\n", name);
  }

  {
    int i;
    const int n_steps = 1;
    double time;

    for (i = 0; i < n_steps; i++) {
      fprintf(stdout, "Running until t = %d... ", i + 1);
      if (cem->update(cem->self) == 0 &&
          cem->get_current_time(cem->self, &time) == 0) {
        if (fabs(time - (i + 1)) < 1e-6)
          fprintf(stdout, "PASS\n");
        else {
          return EXIT_FAILURE;
        }
      } else
        return EXIT_FAILURE;
    }

    fprintf(stdout, "Running until t = %f... ", 1000.);
    if (cem->update_until(cem->self, 1000) == 0 &&
        cem->get_current_time(cem->self, &time) == 0) {
      if (fabs(time - 20) < 1e-6)
        fprintf(stdout, "PASS\n");
      else {
        fprintf(stdout, "%f\n", time);
        return EXIT_FAILURE;
      }
    } else
      return EXIT_FAILURE;
  }

  fprintf(stdout, "Finalizing... ");
  err = cem->finalize(cem->self);
  if (err) return EXIT_FAILURE;
  fprintf(stdout, "PASS\n");

  return EXIT_SUCCESS;
}
