# Landsat_BRDF
This code works only for view angle correction of Landsat BRDF impact (Roy et al. 2016).

The derived NBAR solar zenith of this code is defined by default same as the input solar zenith. However, users can choose a modelled solar zenith 
that has only <2° difference to the observed 2011 Landsat 5 and Landsat 7 solar zenith angles (Zhang et al. 2016). 

In either case, this code is *NOT* good to correct Landsat solar angle BRDF, i.e., 
    (1) *NOT* good to correct for Landsat orbit drift impact 
    (2) *NOT* good to correct for Landsat BRDF induced by sun-angle seasonal variation


Roy, D.P., Zhang, H. K., Ju, J., Gomez-Dans, J. L., Lewis, P.E., Schaaf C.B., Sun, Q., Li, J., Huang, H., Kovalskyy, V., 2016,
A general method to normalize Landsat reflectance data to nadir BRDF adjusted reflectance, Remote Sensing of Environment, 176, 255-271.   
https://www.sciencedirect.com/science/article/pii/S0034425716300220

Zhang, H. K., Roy, D. P., & Kovalskyy, V. (2016). 
Optimal solar geometry definition for global long-term Landsat time-series bidirectional reflectance normalization. 
IEEE Transactions on Geoscience and Remote Sensing, 54(3), 1410-1418.
https://ieeexplore.ieee.org/document/7295567



