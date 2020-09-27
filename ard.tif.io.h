/******************
Hankui on Nov 24 2017
******************/

#ifndef _ARD_TIF_IO_H
#define _ARD_TIF_IO_H


#include "geotiff.h"
#include "xtiffio.h"


void   ReadBands(char* *PathSR_1d, int b, int **Bands_2d);
void   SetUpTIFFDirectory(TIFF *tif, int WIDTH, int HEIGHT, int LX, int UY, int PixelSize);
void   SetUpGeoKeys(GTIF *gtif);
int write_tiff(int b, int **Bands_2d, char* *PathBR_1d, int nrow, int ncol, int LX, int UY, int PixelSize);

#endif
