/************************************************************************************************************************************************************************
BRDF normalization for ARD Landsat 4/5/7/8 TM/ETM+/OLI reflective bands (blue, green, red, NIR, SWIR-1 & SWIR-2) reflectance (top-of-atmosphere or surface reflectance)
This setting normalizes the reflectance to NBAR with the same solar zenith as the input. 

Usage: 
.exe <hh.vv> <ard_input_solar_zenith_path> <ard_input_solar_azimuth_path> <ard_input_view_zenith_path> <ard_input_view_azimuth_path> <ard_input_band1path> <ard_input_band2path> ... <ard_input_band6path> <output_band1path> <output_band2path> ... <output_band6path>

Parameters:
	(1) hh.vv:
		Landsat ARD tile id, see https://landsat.usgs.gov/ard
	(2) <ard_input_solar_zenith_path> <ard_input_solar_azimuth_path> <ard_input_view_zenith_path> <ard_input_view_azimuth_path>:
		four Landsat geometry files:
			ARD data provide the geometry in tiff files 
			Landsat Collection 1 data provoide a tool and an associated angle coeffient files to create per-pixel geometry (https://landsat.usgs.gov/solar-illumination-and-sensor-viewing-angle-coefficient-file) 
	(3) <ard_input_band1path> <ard_input_band2path> ... <ard_input_band6path>: 
		six band top-of-atmosphere or surface reflectance tiff files in the order of blue, green, red, NIR, SWIR-1 & SWIR-2

	(4) <output_band1path> <output_band2path> ... <output_band6path>: 
		six band top-of-atmosphere or surface reflectance tiff files in the order of blue, green, red, NIR, SWIR-1 & SWIR-2

Running example (check the associated run.sh file): 
	./Landsat.BRDF.normalize h04v02 \
		/gpfs/data2/temp.test/brdf.test/LT05_CU_004002_19860819_20170711_C01_V01_SOZ4.tif \
		/gpfs/data2/temp.test/brdf.test/LT05_CU_004002_19860819_20170711_C01_V01_SOA4.tif \
		/gpfs/data2/temp.test/brdf.test/LT05_CU_004002_19860819_20170711_C01_V01_SEZ4.tif \
		/gpfs/data2/temp.test/brdf.test/LT05_CU_004002_19860819_20170711_C01_V01_SEA4.tif \
		/gpfs/data2/temp.test/brdf.test/LT05_CU_004002_19860819_20170711_C01_V01_SRB1.tif \
		/gpfs/data2/temp.test/brdf.test/LT05_CU_004002_19860819_20170711_C01_V01_SRB2.tif \
		/gpfs/data2/temp.test/brdf.test/LT05_CU_004002_19860819_20170711_C01_V01_SRB3.tif \
		/gpfs/data2/temp.test/brdf.test/LT05_CU_004002_19860819_20170711_C01_V01_SRB4.tif \
		/gpfs/data2/temp.test/brdf.test/LT05_CU_004002_19860819_20170711_C01_V01_SRB5.tif \
		/gpfs/data2/temp.test/brdf.test/LT05_CU_004002_19860819_20170711_C01_V01_SRB7.tif \
		${outputdir}/LT05_CU_004002_19860819_20170711_C01_V01_BRB1_InC.tif \
		${outputdir}/LT05_CU_004002_19860819_20170711_C01_V01_BRB2_InC.tif \
		${outputdir}/LT05_CU_004002_19860819_20170711_C01_V01_BRB3_InC.tif \
		${outputdir}/LT05_CU_004002_19860819_20170711_C01_V01_BRB4_InC.tif \
		${outputdir}/LT05_CU_004002_19860819_20170711_C01_V01_BRB5_InC.tif \
		${outputdir}/LT05_CU_004002_19860819_20170711_C01_V01_BRB7_InC.tif \
		2>warninglog.txt 
Check the following for algorithms: 
	(1) Roy, D. P., Zhang, H. K., Ju, J., Gomez-Dans, J. L., Lewis, P. E., Schaaf, C. B., Sun Q., Li J., Huang H., & Kovalskyy, V. (2016). 
	A general method to normalize Landsat reflectance data to nadir BRDF adjusted reflectance. 
	Remote Sensing of Environment, 176, 255-271.

	(2) Zhang, H. K., Roy, D. P., & Kovalskyy, V. (2016). 
	Optimal solar geometry definition for global long-term Landsat time-series bidirectional reflectance normalization. 
	IEEE Transactions on Geoscience and Remote Sensing, 54(3), 1410-1418.


Hankui Jan 15, 2015
Revised by Hankui on Aug 7, 2017
Revised by Alex on Nov 21 2017 for ARD data
************************************************************************************************************************************************************************/

#include <stdlib.h>
// #include <tiff.h>
// #include <tiffio.h>
#include <math.h>
#include <string.h>
#include <omp.h>

#include "model.h"

#include "Landsat.BRDF.normalize.h"
#include "ard.tif.io.h"


// #include "geovalues.h"


int NPixels; 


char * __basename (const char *filename)
{
  char *p = strrchr (filename, '/');
  return p ? p + 1 : (char *) filename;
}


/******************************************************************************************/
int main(int argc, char* argv[])
{
	const int NBands = 6; // six Landsat bands need to be normalized
	int NGeometries = 4; // number of geometry tiff layers (4: solar zenith & azimuth; view zenith & azimuth)
	int nrows, ncols; 
	
	int PixelSize;
	
	int b, p, a, count, status, HDFtileH, HDFtileV, LX, RX, UY, LY, Albers_ULX, Albers_ULY; 

    float solar_zenith, view_zenith, relative_azimuth, solar_zenith_output; //SOZ4, SEZ4, SOA4 - SEA4, solar_zenith respectively
    int   Landsat_input  [NBands]; 
    int   Landsat_output [NBands]; 
    
    int MODELLED_SOLAR = 0;
    
    if (argc<18) 
	{
		status = 1;
		fprintf(stderr, "status=%d\tIncorrect number of command-line parameters; \n", status);
		fprintf(stderr, "\tthe following parameters are supplied:\n");
		for (a=1; a<argc; a++) {fprintf(stderr, "\t%s\n", argv[a]);}
		UsageRule(argv[0]);
		exit(status);
	}

	//*******************************************************************************************
	// initialization 
    nrows = ncols = 5000;
    NPixels = 25000000;
    PixelSize = 30;
    Albers_ULX  = -2565585;
    Albers_ULY  =  3314805;

    char** PathSR_1d; PathSR_1d = (char**)malloc(NBands      * sizeof(char*));
    char** PathBR_1d; PathBR_1d = (char**)malloc(NBands      * sizeof(char*));
    char** PathGE_1d; PathGE_1d = (char**)malloc(NGeometries * sizeof(char*));

    sscanf(argv[1],"h%dv%d", &HDFtileH, &HDFtileV); //printf ("h: %d v: %d \n", HDFtileH, HDFtileV);
    for (a=2;  a<6;  a++) { PathGE_1d[a-2]  = argv[a]; }
    for (a=6;  a<12; a++) { PathSR_1d[a-6]  = argv[a]; }
    for (a=12; a<18; a++) { PathBR_1d[a-12] = argv[a]; }
    
    if (argc>18 )  { MODELLED_SOLAR = atoi(argv[a]);}
    
	//*******************************************************************************************
	// basic metadata 
    LX = Albers_ULX + (HDFtileH * ncols) * PixelSize;
    RX = LX + ncols * PixelSize;
    UY = Albers_ULY - (HDFtileV * nrows) * PixelSize;
    LY = UY - nrows * PixelSize;//printf ("LX: %d RX: %d UY: %d LY: %d \n", LX, RX, UY, LY);
	
    char basename [1024];
    char yearname [1024];
    char monthname[1024];
    char dayname  [1024];
    strcpy(basename ,  __basename(PathSR_1d[0])+0 );
    strcpy(yearname ,  __basename(PathSR_1d[0])+15); yearname [4] = '\0';
    strcpy(monthname,  __basename(PathSR_1d[0])+19); monthname[2] = '\0';
    strcpy(dayname  ,  __basename(PathSR_1d[0])+21); dayname  [2] = '\0';
    
    printf("%s %s %s %s\n", basename,yearname, monthname, dayname);
    
    int month, year, day;
    
    year =atoi( yearname);
    month=atoi(monthname);
    day  =atoi(  dayname);
    printf("MODELLED_SOLAR=%d %d %d %d\n", MODELLED_SOLAR,year, month, day);
    
	//*******************************************************************************************
	// allocate memory read tiff 
    int **Bands_2d; Bands_2d = (int **)malloc(sizeof(int *)*NBands);
    for (b=0; b<NBands; b++) { Bands_2d[b] = (int *)malloc(NPixels * sizeof(int)); }

    int **Geometry_2d; Geometry_2d = (int **)malloc(sizeof(int *)*NGeometries);
    for (b=0; b<NGeometries; b++) { Geometry_2d[b] = (int *)malloc(NPixels * sizeof(int)); }

    for (b = 0; b<NGeometries; b++) { ReadBands(PathGE_1d, b, Geometry_2d); }
    //Geometry_2d[0] SOZ4     float solar_zenith = 50; //SOZ4
    //Geometry_2d[1] SOA4     float view_zenith = 7;//SEZ4
    //Geometry_2d[2] SEZ4     float relative_azimuth = 100;//SOA4 - SEA4
    //Geometry_2d[3] SEA4     float solar_zenith_output = solar_zenith;

    for (b = 0; b<NBands;      b++) { ReadBands(PathSR_1d, b, Bands_2d); }

	
	//*******************************************************************************************
	// BRDF normalization
	printf("\n doing BRDF normalization....\n");
    int i,j;
    double ax, ay;
    // for (p = 0; p<NPixels; p++) // for cycle for each pixel 
    for (p=0, i = 0; i<nrows; i++) // for cycle for each pixel 
    for (j = 0; j<ncols; j++, p++) // for cycle for each pixel 
	{
        ax = LX+j*PixelSize;
        ay = LY-i*PixelSize;
        
        if (p%10000000==0) printf("line %d\n", i);
		for (b = 0; b<NBands; b++)
		{
			Landsat_input[b] = Bands_2d[b][p]; // geometry and surface reflectance //int Landsat_input [6] = {200, 300, 400, 2000, 1000, 1500}; 
		}
        
		solar_zenith = (float)Geometry_2d[0][p] * 0.01; //SOZ4
        
        if (MODELLED_SOLAR == 0) //default 
            solar_zenith_output = solar_zenith;
        else // modelled solar zenith in Zhang et al. (2016)
            solar_zenith_output = get_modelled_solar_zenith(ax, ay, year, month, day);
         
		view_zenith  = (float)Geometry_2d[2][p] * 0.01; //SEZ4
		relative_azimuth = ((float)Geometry_2d[1][p] - (float)Geometry_2d[3][p]) * 0.01; // SOA4 - SEA4
		NBAR_calculate_global(Landsat_input, view_zenith, solar_zenith, relative_azimuth, solar_zenith_output, Landsat_output);// give new value to nbar
		for (b = 0; b<NBands; b++)
		{
			if(Bands_2d[b][p] != 20000 && Bands_2d[b][p] != -9999) {Bands_2d[b][p] = Landsat_output[b];}
		}
		Bands_2d[NBands-1][p] = (int)(solar_zenith_output*100+0.5);
	}// end loop 

	//printf("%d %d %d %d %d %d \n", Landsat_output[0], Landsat_output[1], Landsat_output[2], Landsat_output[3], Landsat_output[4], Landsat_output[5]); 
	//printf("%d %d %d %d %d %d \n", Bands_2d[0][24999999], Bands_2d[1][24999999], Bands_2d[2][24999999], Bands_2d[3][24999999], Bands_2d[4][24999999], Bands_2d[5][24999999]); 
    //should be  335 563 590 1908 1789 1232

	
	//*******************************************************************************************
 	// save nbar tiff file
	printf("\nsaving nbar tiff file....\n");
    for (b = 0; b<NBands; b++) { write_tiff(b, Bands_2d, PathBR_1d, ncols, nrows, LX, UY, PixelSize); }

	
	//*******************************************************************************************
	// release memory
	free (PathSR_1d);
	free (PathBR_1d);
	free (PathGE_1d);
	
    for (b=0; b<NGeometries; b++) { free(Geometry_2d[b]); }
	free (Geometry_2d); 

    for (b=0; b<NBands; b++) { free (Bands_2d[b]); }
	free (Bands_2d); 
	
	return 0; 
}



//========================================Routines========================================
//----Routine-----------------------define input parameters-----------------------------------
void UsageRule(char * arg_0)
{
	fprintf(stderr, "\nUSAGE: %s <hh.vv> <ard_input_solar_zenith_path> <ard_input_solar_azimuth_path> <ard_input_view_zenith_path> <ard_input_view_azimuth_path> <ard_input_band1path> <ard_input_band2path> ... <ard_input_band6path> <output_band1path> <output_band2path> ... <output_band6path>\n", arg_0);
	// fprintf(stderr, "A NEW CALLING EXAMPLE :\n"
	// "%s\t CONUS.5year.2006to2010.h10v10.v1.6.lcluc.v1.0.hdf 25.15 \\\n"
	// "8bit    \t\t /band1bils/h25v15/result.bil \\ \n"
	// "8bit    \t\t /band2bils/h25v15/result.bil \\ \n"
	// "8bit    \t\t /band3bils/h25v15/result.bil \\ \n"
	// "8bit  ... \n"
	// "16bit_s \t\t /band4bils/h25v15/result.bil \\ \n"
	// "16bit_s \t\t /band5bils/h25v15/result.bil \\ \n",  arg_0);
}


//----Routine------------------------Calculate NBAR---------------------------------------
/***************************************************************************
* calculate nadir normalized reflectance (scaled by 10000)
* Landsat_input:		input six bands (in the order of 1, 2, 3, 4, 5 and 7) Landsat surface reflectance to be normalized (scaled by 10000)
* Landsat_output:		output normalized six bands (in the order of 1, 2, 3, 4, 5 and 7) Landsat surface reflectance (scaled by 10000)
* view_zenith:          viewing zenith angle of the Landsat surface refletance (in degree 0~8)
* solar_zenith:         solar zenith angle of the Landsat surface refletance (in degree 0~90)
* relative_azimuth:     relative azimuth angle of obtained by (soloar_azimuth-viewing_azimuth) or (viewing_azimuth-soloar_azimuth) (in degree -180~180)
* solar_zenith_output:  output solar zenith angle the "Landsat_input" to be normalize to
***************************************************************************/
int NBAR_calculate_global(int Landsat_input [], float view_zenith, float solar_zenith,
		float relative_azimuth, float solar_zenith_output,int Landsat_output [])
{
	// Coeffients in  Roy, D. P., Zhang, H. K., Ju, J., Gomez-Dans, J. L., Lewis, P. E., Schaaf, C. B., Sun Q., Li J., Huang H., & Kovalskyy, V. (2016). 
	// A general method to normalize Landsat reflectance data to nadir BRDF adjusted reflectance. 
	// Remote Sensing of Environment, 176, 255-271.
	
	const int pars_12m_global[NBANDS_L][3] = {
			{ 774, 372, 79},
			{1306, 580,178},
			{1690, 574,227},
			{3093,1535,330},
			{3430,1154,453},
			{2658, 639,387},
	};

	//input kernel values
	double nbarkerval[3];
	CalculateKernels(nbarkerval, view_zenith*DE2RA, solar_zenith*DE2RA, relative_azimuth*DE2RA);

	//output kernel values
	double nbarkerval_i[3];
	CalculateKernels(nbarkerval_i, 0*DE2RA, solar_zenith_output*DE2RA, 180.0*DE2RA);

	// process for each band
	double nbar,modis_srf;
	int16 pars[3];
	int i = 0;
	int b;
	for (b=0;b<NBANDS_L;b++)
	{
		for (i=0;i<3;i++)
		{
			pars[i] = pars_12m_global[b][i];
		}

		// it does not matter what is the relative_azimuth value when the view zenith is zero tested by Hankui
		nbar = nbarkerval_i[0] * pars[0] +  nbarkerval_i[1] * pars[1] +  nbarkerval_i[2] * pars[2];
		modis_srf = nbarkerval[0] * pars[0] +  nbarkerval[1] * pars[1] +  nbarkerval[2] * pars[2];

		if (nbar <=0 || modis_srf<=0 || Landsat_input[b]==FILL_VALUE_L)
			Landsat_output[b] = FILL_VALUE_L;
		else
			//ratio based normalization
			Landsat_output[b] = (int)(nbar/modis_srf*Landsat_input[b]);
	}
	return 0;
}
//----Routine------------------------Calculate NBAR---------------------------------------
















