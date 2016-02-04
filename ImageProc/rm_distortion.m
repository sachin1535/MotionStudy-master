function x_caltech = rm_distortion(x_p, K, fc, prin_p, skew, dist_c)

x_caltech = caltech_normalize(x_p, fc, prin_p, dist_c, skew);
x_caltech = [1, 0, 0; 0,1, 0]*K*[x_caltech;1];


