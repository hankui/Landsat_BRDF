
#ifndef _MODEL_H
#define _MODEL_H

#include "tiff.h"

int16 calc_refl(int16 *pars, float vzn, float szn, float raa);

double calc_refl_noround(int16 *pars, float vzn, float szn, float raa);

void CalculateKernels(double *resultsArray, double tv, double ti, double phi);

#define DE2RA 0.0174532925199432956

#endif




