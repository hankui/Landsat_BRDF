/******************
Hankui on Nov 24 2017
******************/

#include <stdlib.h>

#include "ard.tif.io.h"
#include "geovalues.h"

extern int NPixels; 

//----Routine--------------------Read tif signed 16 int-----------------------------------
void ReadBands(char* *PathSR_1d, int b, int **Bands_2d)
{
	signed int i;                // i  - сквозной порядковый номер пикселя
	unsigned int tl, tw, ninner; // tl - tif tile length, tw - tif tile width
	uint32 imageWidth, imageLength;
	uint32 tileWidth, tileLength;
	uint32 x, y;
	tdata_t buf;

	printf ("%s %s %s\n", "Reading ", PathSR_1d[b], "..." );
	TIFF* tif = TIFFOpen(PathSR_1d[b], "r");

	if (tif)
	{
		TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &imageWidth); TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &imageLength);
		TIFFGetField(tif, TIFFTAG_TILEWIDTH,  &tileWidth);  TIFFGetField(tif, TIFFTAG_TILELENGTH,  &tileLength);

		buf = _TIFFmalloc(TIFFTileSize(tif));
		for (y = 0; y < imageLength; y += tileLength)
		{
			for (x = 0; x < imageWidth; x += tileWidth)
			{
				TIFFReadTile(tif, buf, x, y, 0, 0);
				for (tl=0, ninner=0; tl<tileLength; tl++)
				{
					for (tw=0; tw<tileWidth; tw++,ninner++) 
					{
						if ((tl+y)>=imageLength || (tw+x)>=imageWidth) continue;
						i = (tl+y)*imageWidth + (tw+x);
						Bands_2d[b][i] = ((int16*)buf)[ninner];
					}
				}
			}
		}
		_TIFFfree(buf);
		TIFFClose(tif);
	}
	else {
		printf("cannot open the tiff file: %s", PathSR_1d[b]);
	}
}
//----Routine--------------------Read tif signed 16 int-----------------------------------

//----Routine----------------------Set geotiff fields-------------------------------------
void SetUpTIFFDirectory(TIFF *tif, int WIDTH, int HEIGHT, int LX, int UY, int PixelSize)
{
	double tiepoints[6]={0,0,0,(double)LX,(double)UY,0.0};
	//double tiepoints[6]={0,0,0,-2565585,3314805,0.0};
	double pixscale[3]={(double)PixelSize,(double)PixelSize,0};

	TIFFSetField(tif,TIFFTAG_IMAGEWIDTH,    WIDTH);
	TIFFSetField(tif,TIFFTAG_IMAGELENGTH,   HEIGHT);
	// TIFFSetField(tif,TIFFTAG_COMPRESSION,   COMPRESSION_DEFLATE);// NOTE THAT ENVI cannot open deflate version of tiff 
	TIFFSetField(tif,TIFFTAG_PHOTOMETRIC,   PHOTOMETRIC_MINISBLACK);
	TIFFSetField(tif,TIFFTAG_PLANARCONFIG,  PLANARCONFIG_CONTIG);
	TIFFSetField(tif,TIFFTAG_BITSPERSAMPLE, 16);
	TIFFSetField(tif,TIFFTAG_SAMPLEFORMAT,  SAMPLEFORMAT_INT); // нужно дабы сохранить в signed 16 int
	TIFFSetField(tif,TIFFTAG_ROWSPERSTRIP,  HEIGHT);
	TIFFSetField(tif,TIFFTAG_GEOTIEPOINTS,  6,tiepoints);
	TIFFSetField(tif,TIFFTAG_GEOPIXELSCALE, 3,pixscale);
}
//----Routine----------------------Set geotiff fields-------------------------------------


//----Routine-----------------Установить ещё какуюто хуйню--------------------------------
void SetUpGeoKeys(GTIF *gtif)
{
	//https://github.com/ufz/geotiff/blob/master/bin/makegeo.c
	// http://geotiff.maptools.org/spec/geotiff2.7.html#2.7
	// // * 6.3.1 GeoTIFF General Codes * // // 
	// 6.3.1.1 Model Type Codes
	GTIFKeySet(gtif, GTModelTypeGeoKey, TYPE_SHORT, 1, ModelTypeProjected);
	// 6.3.1.2 Raster Type Codes
	GTIFKeySet(gtif, GTRasterTypeGeoKey, TYPE_SHORT, 1, RasterPixelIsArea);
	// GTIFKeySet(gtif, GTCitationGeoKey, TYPE_ASCII, 0, "Just An Example");

	// // * 6.3.2 Geographic CS Codes * // // 
	//6.3.2.1 Geographic CS Type Codes
	GTIFKeySet(gtif, GeographicTypeGeoKey,   TYPE_SHORT, 1, 32767          );
	GTIFKeySet(gtif, GeogAngularUnitsGeoKey, TYPE_SHORT, 1, Angular_Degree);
	GTIFKeySet(gtif, GeogLinearUnitsGeoKey,  TYPE_SHORT, 1, Linear_Meter);

	//6.3.2.2 Geodetic Datum Codes//
	GTIFKeySet(gtif, GeogCitationGeoKey, TYPE_ASCII, 0, "WGS 84");
	GTIFKeySet(gtif, GeogGeodeticDatumGeoKey, TYPE_SHORT, 1, 6326); // Datum_WGS84
	//6.3.2.3 Ellipsoid Codes //
	GTIFKeySet(gtif, GeogEllipsoidGeoKey, TYPE_SHORT,     1, 7030 ); // Ellipse_WGS_84

	GTIFKeySet(gtif, GeogSemiMajorAxisGeoKey, TYPE_DOUBLE, 1, (double)6378140);
	GTIFKeySet(gtif, GeogInvFlatteningGeoKey, TYPE_DOUBLE, 1, (double)298.2569999999986);

	// // * 6.3.3 Projected CS Codes * // // 
	//6.3.3.1 Projected CS Type Codes
	// GTIFKeySet(gtif, PCSCitationGeoKey, TYPE_ASCII, 0, "Sinusoidal");

	//https://github.com/opengdp/opengdp/blob/master/MRTSwath/InitGeoTiff.c
	GTIFKeySet(gtif, GTCitationGeoKey, TYPE_ASCII, 0, "Albers");
	GTIFKeySet(gtif, ProjCoordTransGeoKey, TYPE_SHORT, 1, CT_AlbersEqualArea);

	GTIFKeySet(gtif, GTModelTypeGeoKey, TYPE_SHORT, 1, ModelTypeProjected);
	GTIFKeySet(gtif, GTRasterTypeGeoKey, TYPE_SHORT, 1, RasterPixelIsArea);
	GTIFKeySet(gtif, ProjStdParallel1GeoKey, TYPE_DOUBLE, 1, (double)29.5);
	GTIFKeySet(gtif, ProjStdParallel2GeoKey, TYPE_DOUBLE, 1, (double)45.5);
	GTIFKeySet(gtif, ProjNatOriginLongGeoKey, TYPE_DOUBLE, 1, (double)-96);
	GTIFKeySet(gtif, ProjNatOriginLatGeoKey, TYPE_DOUBLE, 1, (double)23);
	GTIFKeySet(gtif, ProjFalseEastingGeoKey, TYPE_DOUBLE, 1, (double)0);
	GTIFKeySet(gtif, ProjFalseNorthingGeoKey, TYPE_DOUBLE, 1, (double)0);
	GTIFKeySet(gtif, ProjLinearUnitsGeoKey, TYPE_SHORT, 1, Linear_Meter);

}
//----Routine-----------------Установить ещё какуюто хуйню--------------------------------


//----Routine--------------------------Write geotiff--------------------------------------
int write_tiff(int b, int **Bands_2d, char* *PathBR_1d, int nrow, int ncol, int LX, int UY, int PixelSize)
{
	int i;
	signed short int *raster; raster = (signed short int *)malloc(nrow*ncol * sizeof(signed short int));

	for (i=0;i<NPixels;i++)
	{
		if (Bands_2d[b][i] >  20000) {raster[i] = 20000;}
		else if (Bands_2d[b][i] <  -32768) {raster[i] = -32768;}
		else {raster[i] = Bands_2d[b][i];}
	}

	TIFF *tif=(TIFF*)0;  /* TIFF-level descriptor */
	GTIF *gtif=(GTIF*)0; /* GeoKey-level descriptor */

	tif=XTIFFOpen(PathBR_1d[b],"w");
	if (!tif) goto failure;

	gtif = GTIFNew(tif);
	if (!gtif)
	{
		printf("failed in GTIFNew\n");
		goto failure;
	}

	SetUpTIFFDirectory(tif, nrow, ncol, LX, UY, PixelSize); // http://geotiff.maptools.org/spec/geotiff2.7.html#2.7
	SetUpGeoKeys(gtif);
	
	// write
	for (i=0;i<nrow;i++)
		if (!TIFFWriteScanline(tif, raster+i*ncol, i, 0))
		{
			printf("something is wrong during writing\n");
			break;
		}

	GTIFWriteKeys(gtif);
	
	GTIFFree(gtif);
	XTIFFClose(tif);
	free(raster);
	return 0;

	failure:
	printf("failure in makegeo\n");
	if (tif) TIFFClose(tif);
	// if (gtif) GTIFFree(gtif);
	free(raster);
	
	return -1;
}
//----Routine--------------------------Write geotiff--------------------------------------
//========================================Routines========================================
