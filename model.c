#include <math.h>
#include <stdio.h>
#include <stdlib.h>
//#include "hdf.h"
#include "model.h"



#define NUMBER_BS_KERNELALBEDOS 91

#define MODSPABR 1.0
#define MODSPAHB 2.0
//#define MAX_PAR (1000)
#define MAX_PAR (10000) //MODIFIED by hankui
#define PAR_FILL_VALUE (32767)
#ifndef PI
#define PI (3.1415926)
#endif


#define MAX(x,y) ((x) >= (y)) ? (x) : (y)
#define MIN(x,y) ((x) <= (y)) ? (x) : (y)


double bsaLiSparseRker[91] =
	{ -1.287889, -1.287944, -1.288058, -1.288243, -1.288516, -1.288874,
		-1.289333, -1.289877, -1.290501, -1.291205,
		-1.291986, -1.292855, -1.293812, -1.294861, -1.295982, -1.297172,
		-1.298447, -1.299805, -1.301228, -1.302730,
		-1.304319, -1.306023, -1.307829, -1.309703, -1.311615, -1.313592,
		-1.315674, -1.317847, -1.320090, -1.322396,
		-1.324763, -1.327212, -1.329779, -1.332436, -1.335150, -1.337899,
		-1.340710, -1.343621, -1.346625, -1.349686,
		-1.352780, -1.355925, -1.359123, -1.362429, -1.365825, -1.369290,
		-1.372774, -1.376264, -1.379825, -1.383465,
		-1.387144, -1.390834, -1.394543, -1.398269, -1.402071, -1.405931,
		-1.409818, -1.413691, -1.417533, -1.421378,
		-1.425284, -1.429204, -1.433052, -1.436887, -1.440693, -1.444474,
		-1.448251, -1.452013, -1.455691, -1.459185,
		-1.462709, -1.466247, -1.469651, -1.472869, -1.476028, -1.479219,
		-1.482235, -1.485117, -1.487871, -1.490654,
		-1.493254, -1.495838, -1.498365, -1.500982, -1.503707, -1.506758,
		-1.510596, -1.516135, -1.526204, -1.555139, -37286.245539
};

double bsaRossThickker[91] = {
		-0.021079, -0.021026, -0.020866, -0.020598, -0.020223, -0.019740,
		-0.019148, -0.018447, -0.017636, -0.016713,
		-0.015678, -0.014529, -0.013265, -0.011883, -0.010383, -0.008762,
		-0.007018, -0.005149, -0.003152, -0.001024,
		0.001236, 0.003633, 0.006170, 0.008850, 0.011676, 0.014653, 0.017785,
		0.021076, 0.024531, 0.028154,
		0.031952, 0.035929, 0.040091, 0.044444, 0.048996, 0.053752, 0.058720,
		0.063908, 0.069324, 0.074976,
		0.080873, 0.087026, 0.093443, 0.100136, 0.107117, 0.114396, 0.121987,
		0.129904, 0.138161, 0.146772,
		0.155755, 0.165127, 0.174906, 0.185112, 0.195766, 0.206890, 0.218509,
		0.230649, 0.243337, 0.256603,
		0.270480, 0.285003, 0.300209, 0.316139, 0.332839, 0.350356, 0.368744,
		0.388061, 0.408372, 0.429747,
		0.452264, 0.476012, 0.501088, 0.527602, 0.555677, 0.585457, 0.617102,
		0.650800, 0.686769, 0.725269,
		0.766607, 0.811157, 0.859381, 0.911864, 0.969367, 1.032915, 1.103966,
		1.184742, 1.279047, 1.394945, 1.568252
};

void GetPhaang(cos1, cos2, sin1, sin2, cos3, cosres, res, sinres)
	double cos1, cos2, sin1, sin2, cos3;
	double *cosres, *res, *sinres;
{
	*cosres = cos1 * cos2 + sin1 * sin2 * cos3;
	*res = acos(MAX(-1., MIN(1., *cosres)));
	*sinres = sin(*res);

	return;
}

void GetDistance(tan1, tan2, cos3, res)
	double tan1, tan2, cos3;
	double *res;
{
	double temp;

	temp = tan1 * tan1 + tan2 * tan2 - 2. * tan1 * tan2 * cos3;
	*res = sqrt(MAX(0., temp));

	return;
}

void GetpAngles(brratio, tan1, sinp, cosp, tanp)
	double brratio, tan1;
	double *sinp, *cosp, *tanp;
{
	double angp;

	*tanp = brratio * tan1;
	if (*tanp < 0)
		*tanp = 0.;
	angp = atan(*tanp);
	*sinp = sin(angp);
	*cosp = cos(angp);
	return;
}

void GetOverlap(hbratio, distance, cos1, cos2, tan1, tan2, sin3, overlap, temp)
	double hbratio, distance, cos1, cos2, tan1, tan2, sin3;
	double *overlap, *temp;
{
	double cost, sint, tvar;

	*temp = 1. / cos1 + 1. / cos2;
	cost =
	    hbratio * sqrt(distance * distance + tan1 * tan1 * tan2 * tan2 * sin3 * sin3) / *temp;
	cost = MAX(-1., MIN(1., cost));
	tvar = acos(cost);
	sint = sin(tvar);
	*overlap = 1. / PI * (tvar - sint * cost) * (*temp);
	*overlap = MAX(0., *overlap);

	return;
}

void LiKernel(hbratio, brratio, tantv, tanti, sinphi, cosphi, result, SparseFlag, RecipFlag)
	double hbratio, brratio, tantv, tanti, sinphi, cosphi;
	double *result;
	int SparseFlag, RecipFlag;
{
	double sintvp, costvp, tantvp, sintip, costip, tantip;
	double phaangp, cosphaangp, sinphaangp, distancep, overlap, temp;

	GetpAngles(brratio, tantv, &sintvp, &costvp, &tantvp);
	GetpAngles(brratio, tanti, &sintip, &costip, &tantip);
	GetPhaang(costvp, costip, sintvp, sintip, cosphi, &cosphaangp, &phaangp, &sinphaangp);
	GetDistance(tantvp, tantip, cosphi, &distancep);
	GetOverlap(hbratio, distancep, costvp, costip, tantvp, tantip, sinphi, &overlap, &temp);

	if (SparseFlag) {
		if (RecipFlag) {
			*result = overlap - temp + 1. / 2. * (1. + cosphaangp) / costvp / costip;
		} else {
			*result = overlap - temp + 1. / 2. * (1. + cosphaangp) / costvp;
		}
	} else {
		if (RecipFlag) {
			*result = (1 + cosphaangp) / (costvp * costip * (temp - overlap)) - 2.;
		} else {
			*result = (1 + cosphaangp) / (costvp * (temp - overlap)) - 2.;
		}
	}
	
	return;
}

void CalculateKernels(resultsArray, tv, ti, phi)	
	double *resultsArray;
	double tv, ti, phi;
/*tv view zenith angle, 
	ti solar zenith angle, 
	phi relative azimuth angle */
{
	int /*currker,*/ SparseFlag, RecipFlag;
	double cosphi, sinphi;
	double costv, costi, sintv, sinti, tantv, tanti;
	double phaang, sinphaang, cosphaang;
	double rosselement;// , distance;
	//int iphi;


	resultsArray[0] = 1.;

	cosphi = cos(phi);

	costv = cos(tv);
	costi = cos(ti);
	sintv = sin(tv);
	sinti = sin(ti);
	GetPhaang(costv, costi, sintv, sinti, cosphi, &cosphaang, &phaang, &sinphaang);
	rosselement = (PI / 2. - phaang) * cosphaang + sinphaang;
	resultsArray[1] = rosselement / (costi + costv) - PI / 4.;

/*finish rossthick kernal */

	sinphi = sin(phi);
	tantv = tan(tv);
	tanti = tan(ti);

	// SparseFlag = TRUE;
	// RecipFlag = TRUE;
	SparseFlag = 1;
	RecipFlag = 1;
	LiKernel(MODSPAHB, MODSPABR, tantv, tanti, sinphi, cosphi,
		 &resultsArray[2], SparseFlag, RecipFlag);

	return;
}

int16 calc_refl(int16 *pars, float vzn, float szn, float raa)
{
	if(pars[0] > MAX_PAR || pars[1] > MAX_PAR || pars[2] > MAX_PAR)
		return PAR_FILL_VALUE;

	double nbarkerval[3];
	double nbar;
	//int k;

	CalculateKernels(nbarkerval, vzn*DE2RA, szn*DE2RA, raa*DE2RA);
	
	nbar = nbarkerval[0] * pars[0] +  nbarkerval[1] * pars[1] +  nbarkerval[2] * pars[2];
	
	if(nbar < 0 && nbar > -30)
		nbar = 0.0;
	else if (nbar < 0 || nbar > MAX_PAR)
		return PAR_FILL_VALUE;

	int16 ret = rintf(nbar);

	return ret;
}

/***************************************************************************
 * Added by Hankui to implement calculation without rounding on Jan 18, 2015
***************************************************************************/
double calc_refl_noround(int16 *pars, float vzn, float szn, float raa)
{
	if(pars[0] > MAX_PAR || pars[1] > MAX_PAR || pars[2] > MAX_PAR)
		return PAR_FILL_VALUE;

	double nbarkerval[3];
	double nbar;
	//int k;

	CalculateKernels(nbarkerval, vzn*DE2RA, szn*DE2RA, raa*DE2RA);

	nbar = nbarkerval[0] * pars[0] +  nbarkerval[1] * pars[1] +  nbarkerval[2] * pars[2];

	if(nbar < 0 && nbar > -30)
		nbar = 0.0;
	else if (nbar < 0 || nbar > MAX_PAR)
		return PAR_FILL_VALUE;

	return nbar;

}
