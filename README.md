# DbyD_DDFCsparsePA_2D
A local method to reduce the Gibbs oscillation on 2D image

These are the Matlab files used in Shi, R., Jung, J., Schweser, F.: Two-dimensional local Fourier image reconstruction via domain 
decomposition Fourier continuation method. PLoS ONE 14(1):e0197963. (2019)

DbyD_abs_sample_region.m: This file shows the reconstruction of sample regions by dimension-by-dimension Fourier continuation 
sparse PA method I proposed in the article. 

DbyD_abs_abs_whole.m: This file shows the reconstruction of whole image by splitting the whole image into several sub-images and 
applying dimension-by-dimension Fourier continuation sparse PA method on each sub-image, finally stitching the reconstructions together.


