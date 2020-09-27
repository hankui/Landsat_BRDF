#!/bin/bash

# Hankui wrote in on Nov 21, 2017
# bash run.example.sh

# outputdir=/gpfs/data2/temp.test/brdf.test/c.hank;
inputARD=/gpfs/data2/workspace/zhangh/brdf.test/
outputdir=/gpfs/data2/workspace/zhangh/brdf.test/c.hank.git;
outputdir=/gpfs/data2/workspace/zhangh/brdf.test/c.hank.git2;
mkdir -p ${outputdir};

date_start=`date|awk -F"[ :]+" '{print $3*3600*24 + $4*60*60 + $5*60 + $6}'`;

./Landsat.BRDF.normalize.v1.0 h04v02 \
	${inputARD}/LT05_CU_004002_19860819_20170711_C01_V01_SOZ4.tif \
	${inputARD}/LT05_CU_004002_19860819_20170711_C01_V01_SOA4.tif \
	${inputARD}/LT05_CU_004002_19860819_20170711_C01_V01_SEZ4.tif \
	${inputARD}/LT05_CU_004002_19860819_20170711_C01_V01_SEA4.tif \
	${inputARD}/LT05_CU_004002_19860819_20170711_C01_V01_SRB1.tif \
	${inputARD}/LT05_CU_004002_19860819_20170711_C01_V01_SRB2.tif \
	${inputARD}/LT05_CU_004002_19860819_20170711_C01_V01_SRB3.tif \
	${inputARD}/LT05_CU_004002_19860819_20170711_C01_V01_SRB4.tif \
	${inputARD}/LT05_CU_004002_19860819_20170711_C01_V01_SRB5.tif \
	${inputARD}/LT05_CU_004002_19860819_20170711_C01_V01_SRB7.tif \
	${outputdir}/LT05_CU_004002_19860819_20170711_C01_V01_BRB1_InC.tif \
	${outputdir}/LT05_CU_004002_19860819_20170711_C01_V01_BRB2_InC.tif \
	${outputdir}/LT05_CU_004002_19860819_20170711_C01_V01_BRB3_InC.tif \
	${outputdir}/LT05_CU_004002_19860819_20170711_C01_V01_BRB4_InC.tif \
	${outputdir}/LT05_CU_004002_19860819_20170711_C01_V01_BRB5_InC.tif \
	${outputdir}/LT05_CU_004002_19860819_20170711_C01_V01_BRB7_InC.tif \
    1 \
    2>warninglog.txt 

gdal_translate ${outputdir}/LT05_CU_004002_19860819_20170711_C01_V01_BRB3_InC.tif ${outputdir}/gdal.LT05_CU_004002_19860819_20170711_C01_V01_BRB3_InC.tif 
gdal_translate ${inputARD}/LT05_CU_004002_19860819_20170711_C01_V01_SRB3.tif ${outputdir}/gdal.LT05_CU_004002_19860819_20170711_C01_V01_SRB3.tif 
gdal_translate ${outputdir}/LT05_CU_004002_19860819_20170711_C01_V01_BRB7_InC.tif ${outputdir}/gdal.LT05_CU_004002_19860819_20170711_C01_V01_BRB7_InC.tif
gdal_translate ${inputARD}/LT05_CU_004002_19860819_20170711_C01_V01_SOZ4.tif ${outputdir}/gdal.LT05_CU_004002_19860819_20170711_C01_V01_SOZ4.tif
	
date_end=`date|awk -F"[ :]+" '{print $3*3600*24 + $4*60*60 + $5*60 + $6}'`;
time_diff=`echo "scale=2;($date_end-$date_start+0.01)*1.0/60.0"|bc`;
date
echo "$time_diff minutes used";
