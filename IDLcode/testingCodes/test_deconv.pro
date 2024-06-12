PRO TEST_DECONV
; uses output from sharp_blur_image to see if information from
; less-smoothed image can be recovered (gaussian psf's only for
; this test, so far).

FXREAD,'fake_contour_fwhm5.fits',fwhm5
FXREAD,'fake_contour_fwhm10.fits',fwhm10

; Perform deconvolution on the more-smoothed contoured image, see if
; resulting contour looks like the less-smoothed contoured image

; generate the smoothing function that hypothetically would have
; blurred the sharper image to the less sharp one:
FWHM=SQRT(10.^2-5^2)
PSF=GAUSSIAN_FUNCTION([FWHM,FWHM]/(2.*SQRT(2.*ALOG(2.))),WIDTH=101)
; force integral under PSF function to equal unity (e.g., normalize
; function)
PSF=PSF/TOTAL(PSF)
; note: fwhm --> sigma in construction of gaussian


NITER=50
FOR I=1,NITER DO MAX_LIKELIHOOD,FWHM10,PSF,DECONV,RE_CONV,FT_PSF=PSF_FT

WRITEFITS,'fake_deconv10-5.fits',DECONV
WRITEFITS,'fake_reconv10-5.fits',RE_CONV
stop,'finished'
END
