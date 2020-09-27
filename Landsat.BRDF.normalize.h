/******************
Hankui Jan 15, 2015
Hankui adapted to Landsat ARD on Nov 24 2017
******************/

#ifndef _L457_BRDF_NORM_H
#define _L457_BRDF_NORM_H

#define NBANDS_L 6
#define FILL_VALUE_L -32768

void   UsageRule(char * arg_0);
int NBAR_calculate_global(int Landsat_input [], float view_zenith, float solar_zenith,float relative_azimuth, float solar_zenith_output,int Landsat_output []);


// int NBAR_calculate_global2(int*Landsat_input, float view_zenith, float solor_zenith, float relative_azimuth, float solor_zenith_output,int* Landsat_output);

// int NBAR_calculate_global(int b, float input_SRF, float view_zenith, float solor_zenith, float relative_azimuth, float solor_zenith_output);
// void l457_srf_tile_update2(l457_srf_tile_info_t *srtile, l457_srf_tile_info_t *tile);
// int l57_srf_tile_normalize(l457_srf_tile_info_t *tile);
// int l57_srf_tile_normalize_albers(l457_tile_info_t *srfComp, char* timeperoid, bool sz_change);

#endif
