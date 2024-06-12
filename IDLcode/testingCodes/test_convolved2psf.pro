PRO TEST_CONVOLVED2PSF
; A little test to see if when one makes an image by semi-randomly
; placing down point sources, each with some gaussian profile of fixed
; FWHM, then convolves that resulting image with a gaussian kernal of
; some chosen width, if the resulting smeared image is equivalent to
; going back to those point sources and just placing down gaussians
; that have the equivalent FWHM of the final convolved image.  The
; answer is YES!  This program confirms that understanding of how
; convolution works. 
;
; Plop down several little 2-D gaussian's, just keep adding
; them up to generate an extended area.  Randomize where they land so
; that several situations can be tested.

; to insure that there is adequate clumping of points to make a
; "galaxy", force 30% of the points to land within 40 pixels of the
; middle of the 100x100 frame.  Force another 40% to land somewhere along a diagonal
; stretching from the upper left to lower right corner within a
; thickness of +/-35 degrees, to insure
; elongation.  Sprinkle the remaining 30% randomly.


; for plotting purposes, generate a user-defined circular plotting symbol
TH=FINDGEN(41)/40.*2.*!PI  ; make an array of angles defining a circle
XSYM=SIN(TH) & YSYM=COS(TH)  ; radius is implicitly unity
USERSYM,XSYM,YSYM,/FILL
  
WINDOW,4,XSIZE=800,YSIZE=800
PLOT,FINDGEN(100),FINDGEN(100),XRA=[0,99],YRA=[0,99],/XST,/YST,/ISO,/NODATA
; Get inner 30% of points
CENTERDIAM=40.  ; number of pixels of diameter for central concentration
RADCENTER=(RANDOMU(SEED,300)-0.5)*CENTERDIAM
ANGCENTER=(RANDOMU(SEED,300)-0.5)*!PI

; get diagonal 40% of points
DIAGTHICK=35.                   ; number of degrees in thickness for diagonal strip
DIAGLENGTH=80.                  ; maximum length from end to end along diagonal
DIAGPOSANG=135.   ; desired position angle of diagonal feature
RADDIAG=(RANDOMU(SEED,400)-0.5)*DIAGLENGTH
ANGDIAG=((RANDOMU(SEED,400)-0.5)*DIAGTHICK+DIAGPOSANG)*!PI/180.

; get the 30% random sprinkle
RANDLENGTH=100.    ; desired extent in x (or y) direction for random sprinkling
RADRAND=(RANDOMU(SEED,300)-0.5)*RANDLENGTH
ANGRAND=(RANDOMU(SEED,300)-0.5)*!PI

; Consolidate
RADIUS=[RADCENTER,RADDIAG,RADRAND]
ANGLE=[ANGCENTER,ANGDIAG,ANGRAND]
AMP=RANDOMU(SEED,1000)*100.      ; amplitude of gaussian

INITPSF=1.0  ; fwhm
INITSIG=INITPSF/(2.*SQRT(2.*ALOG(2)))
IMAGE=FLTARR(100,100)

; Turn into x and y coordinates
X=(RADIUS*COS(ANGLE)+50.)>1.
Y=(RADIUS*SIN(ANGLE)+50.)>1.

; pixelate by forcing rounded numbers
X=ROUND(X)<100.
Y=ROUND(Y)<100.

PLOTS,X,Y,PSYM=8

; Now turn each coordinate into a 2-D gaussian and form an image
; let the FWHM of the gaussian be 1 pixel
FOR I=0,N_ELEMENTS(X)-1 DO BEGIN
; Determine the maximal distance between this point and the furthest
; edge of frame and set the width of the gaussian to be that distance
; to insure the entire frame is filled
   DLLC=SQRT( (X(I)-0)^2 + (Y(I)-0.)^2 )
   DURC=SQRT( (X(I)-100.)^2 + (Y(I)-100.)^2 )
   DULC=SQRT( (X(I)-0.)^2 + (Y(I)-100.)^2 )
   DLRC=SQRT( (X(I)-100.)^2 + (Y(I)-0.)^2 )
   DISTANCE=MAX([DLLC,DURC,DULC,DLRC])
; Force DISTANCE to be a round number as well as an odd number
   DISTANCE=2.*ROUND(DISTANCE)+1.
   PSFIMG=GAUSSIAN_FUNCTION([INITSIG,INITSIG],WIDTH=DISTANCE,MAXIMUM=AMP(I))
; To aid in the overlay onto the image frame, construct a coordinate
; system to go with the gaussian array that is centered on the X,Y point
   XGAUSS=FIX(FINDGEN(DISTANCE)-FIX(DISTANCE/2.)+X(I))
   YGAUSS=FIX(FINDGEN(DISTANCE)-FIX(DISTANCE/2.)+Y(I))
   XGOOD=WHERE(XGAUSS GE 0 AND XGAUSS LE 99)
   YGOOD=WHERE(YGAUSS GE 0 AND YGAUSS LE 99)
   XGAUSS=XGAUSS(XGOOD)
   YGAUSS=YGAUSS(YGOOD)
   XG1=MIN(XGAUSS) & XG2=MAX(XGAUSS)
   YG1=MIN(YGAUSS) & YG2=MAX(YGAUSS)
; Now overlay this 2D gaussian into the image frame, trimming off the
; parts that fall outside the frame
   XS1=MIN(XGOOD) & XS2=MAX(XGOOD)
   YS1=MIN(YGOOD) & YS2=MAX(YGOOD)
   IMAGE(XG1:XG2,YG1:YG2)=IMAGE(XG1:XG2,YG1:YG2)+PSFIMG(XS1:XS2,YS1:YS2)
ENDFOR
WRITEFITS,'fake_galaxy_orig.fits',IMAGE


; Generate a bunch of contours, so that the contour holding 50% of all
; the light can be found
;
; Divide the levels up evenly between the 30%  of the mean value and
; up to 80% of the highest pixel value
ORIGLEV=(FINDGEN(1001)/1000.)*(MAX(IMAGE)*0.8-MEAN(IMAGE)*0.3)+MEAN(IMAGE)*0.3
CONTOUR,IMAGE,FINDGEN(100),FINDGEN(100),LEVELS=ORIGLEV,PATH_INFO=ORIGINFO,$
        PATH_XY=ORIGXY,/CLOSED,/PATH_DATA_COORDS
; Go through each contour level ...
FORIG=ORIGLEV*0.
FOR I=0,N_ELEMENTS(ORIGLEV)-1 DO BEGIN
; Grab all contours associated with this level
   ILEV=WHERE(ORIGINFO.LEVEL EQ I)
; Assume that the biggest contour in this group is the one of interest
; (a more accurate thing to do is to get the contour that encloses the
; point of interest, like the center of the image)
   IBIG=WHERE(ORIGINFO(ILEV).N EQ MAX(ORIGINFO(ILEV).N))
   IBIG=ILEV(IBIG(0))
; Get the starting and ending point for this single contour
   ISTART=ORIGINFO(IBIG).OFFSET
   IF IBIG EQ N_ELEMENTS(ORIGLEV)-1 THEN IEND=ISTART ELSE IEND=ORIGINFO(IBIG+1).OFFSET
; Extract the x/y coordinates for this single contour
   XC=REFORM(ORIGXY(0,ISTART:IEND))
   YC=REFORM(ORIGXY(1,ISTART:IEND))
; Get all the pixels within this contour and add up the total light
   IF N_ELEMENTS(XC) EQ 1 THEN BEGIN
      FORIG(I)=IMAGE(XC,YC)
   ENDIF ELSE BEGIN
      INSIDE=POLYFILLV([XC,XC(0)],[YC,YC(0)],100,100)
      FORIG(I)=TOTAL(IMAGE(INSIDE))
   ENDELSE
ENDFOR

; Now "plot up" the growth curve and figure out the contour that most
; closely matches 50% encircled energy
;
; Since we have such finely-gridded levels, let's just take
; whatever is closed to 50% rather than do some kind of fancy
; interpolation
IPICK=WHERE( ABS(FORIG-TOTAL(IMAGE)/2.) EQ MIN(ABS(FORIG-TOTAL(IMAGE)/2.)))
IPICK=IPICK(0)

; Now get the x,y coordinates for this contour
ILEV=WHERE(ORIGINFO.LEVEL EQ IPICK)
IBIG=WHERE(ORIGINFO(ILEV).N EQ MAX(ORIGINFO(ILEV).N))
IBIG=ILEV(IBIG(0))
ISTART=ORIGINFO(IBIG).OFFSET
IF IBIG EQ N_ELEMENTS(ORIGLEV)-1 THEN IEND=ISTART ELSE IEND=ORIGINFO(IBIG+1).OFFSET
ORIGX=REFORM(ORIGXY(0,ISTART:IEND))
ORIGY=REFORM(ORIGXY(1,ISTART:IEND))

PLOTS,ORIGX,ORIGY,LINEST=0,THICK=2

; Now convolve the image with a gaussian kernal that would make the
; image's effective psf be 3 instead of 1. Determine the kernal
; needed to accomplish this task:
FINPSF=3.0
FINSIG=FINPSF/(2.*SQRT(2.*ALOG(2.)))
KFWHM=SQRT(FINPSF^2-INITPSF^2)
KSIGMA=KFWHM/(2.*SQRT(2.*ALOG(2)))
NEWIMG=GAUSS_SMOOTH(IMAGE,KSIGMA,/EDGE_TRUNCATE)
WRITEFITS,'fake_galaxy_convolved.fits',NEWIMG

CONVLEV=(FINDGEN(1001)/1000.)*(MAX(NEWIMG)*0.8-MEAN(NEWIMG)*0.3)+MEAN(NEWIMG)*0.3
CONTOUR,NEWIMG,FINDGEN(100),FINDGEN(100),LEVELS=CONVLEV,PATH_INFO=CONVINFO,$
        PATH_XY=CONVXY,/CLOSED,/PATH_DATA_COORDS
; Go through each contour level ...
FCONV=CONVLEV*0.
FOR I=0,N_ELEMENTS(CONVLEV)-1 DO BEGIN
; Grab all contours associated with this level
   ILEV=WHERE(CONVINFO.LEVEL EQ I)
; Assume that the biggest contour in this group is the one of interest
; (a more accurate thing to do is to get the contour that encloses the
; point of interest, like the center of the image)
   IBIG=WHERE(CONVINFO(ILEV).N EQ MAX(CONVINFO(ILEV).N))
   IBIG=ILEV(IBIG(0))
; Get the starting and ending point for this single contour
   ISTART=CONVINFO(IBIG).OFFSET
   IF IBIG EQ N_ELEMENTS(CONVLEV)-1 THEN IEND=ISTART ELSE IEND=CONVINFO(IBIG+1).OFFSET
; Extract the x/y coordinates for this single contour
   XC=REFORM(CONVXY(0,ISTART:IEND))
   YC=REFORM(CONVXY(1,ISTART:IEND))
; Get all the pixels within this contour and add up the total light
   IF N_ELEMENTS(XC) EQ 1 THEN BEGIN
      FCONV(I)=NEWIMG(XC,YC)
   ENDIF ELSE BEGIN
      INSIDE=POLYFILLV([XC,XC(0)],[YC,YC(0)],100,100)
      FCONV(I)=TOTAL(NEWIMG(INSIDE))
   ENDELSE
ENDFOR

; Now "plot up" the growth curve and figure out the contour that most
; closely matches 50% encircled energy
;
; Since we have such finely-gridded levels, let's just take
; whatever is closed to 50% rather than do some kind of fancy
; interpolation
IPICK=WHERE( ABS(FCONV-TOTAL(NEWIMG)/2.) EQ MIN(ABS(FCONV-TOTAL(NEWIMG)/2.)))
IPICK=IPICK(0)

; Now get the x,y coordinates for this contour
ILEV=WHERE(CONVINFO.LEVEL EQ IPICK)
IBIG=WHERE(CONVINFO(ILEV).N EQ MAX(CONVINFO(ILEV).N))
IBIG=ILEV(IBIG(0))
ISTART=CONVINFO(IBIG).OFFSET
IF IBIG EQ N_ELEMENTS(CONVLEV)-1 THEN IEND=ISTART ELSE IEND=CONVINFO(IBIG+1).OFFSET
CONVX=REFORM(CONVXY(0,ISTART:IEND))
CONVY=REFORM(CONVXY(1,ISTART:IEND))

PLOTS,CONVX,CONVY,LINEST=0,THICK=4

; Just for grins, see if constructing the image using point sources
; that have a fwhm of 7.0 will duplate the image generated by
; convolving to an effective psf of 7.0.

BLURIMG=FLTARR(100,100)
FOR I=0,N_ELEMENTS(X)-1 DO BEGIN
; Determine the maximal distance between this point and the furthest
; edge of frame and set the width of the gaussian to be that distance
; to insure the entire frame is filled
   DLLC=SQRT( (X(I)-0)^2 + (Y(I)-0.)^2 )
   DURC=SQRT( (X(I)-100.)^2 + (Y(I)-100.)^2 )
   DULC=SQRT( (X(I)-0.)^2 + (Y(I)-100.)^2 )
   DLRC=SQRT( (X(I)-100.)^2 + (Y(I)-0.)^2 )
   DISTANCE=MAX([DLLC,DURC,DULC,DLRC])
; Force DISTANCE to be a round number as well as an odd number
   DISTANCE=2.*ROUND(DISTANCE)+1.
   PSFIMG=GAUSSIAN_FUNCTION([FINSIG,FINSIG],WIDTH=DISTANCE,MAXIMUM=AMP(I))
; I just found out that this function does not preserve area ... in
; other words, the total integrated flux for a PSF having FWHM should
; be exactly the same as one having FHWM that is twice as big, but
; this function does not seem to preserve that relationship.  Force it
; to, by normalizing it:
   OLDIMG=GAUSSIAN_FUNCTION([INITSIG,INITSIG],WIDTH=DISTANCE,MAXIMUM=AMP(I))
   PSFIMG=PSFIMG*TOTAL(OLDIMG)/TOTAL(PSFIMG)
; To aid in the overlay onto the image frame, construct a coordinate
; system to go with the gaussian array that is centered on the X,Y point
   XGAUSS=FIX(FINDGEN(DISTANCE)-FIX(DISTANCE/2.)+X(I))
   YGAUSS=FIX(FINDGEN(DISTANCE)-FIX(DISTANCE/2.)+Y(I))
   XGOOD=WHERE(XGAUSS GE 0 AND XGAUSS LE 99)
   YGOOD=WHERE(YGAUSS GE 0 AND YGAUSS LE 99)
   XGAUSS=XGAUSS(XGOOD)
   YGAUSS=YGAUSS(YGOOD)
   XG1=MIN(XGAUSS) & XG2=MAX(XGAUSS)
   YG1=MIN(YGAUSS) & YG2=MAX(YGAUSS)
; Now overlay this 2D gaussian into the image frame, trimming off the
; parts that fall outside the frame
   XS1=MIN(XGOOD) & XS2=MAX(XGOOD)
   YS1=MIN(YGOOD) & YS2=MAX(YGOOD)
   BLURIMG(XG1:XG2,YG1:YG2)=BLURIMG(XG1:XG2,YG1:YG2)+PSFIMG(XS1:XS2,YS1:YS2)
ENDFOR
WRITEFITS,'fake_galaxy_blur.fits',BLURIMG

BLURLEV=(FINDGEN(1001)/1000.)*(MAX(BLURIMG)*0.8-MEAN(BLURIMG)*0.3)+MEAN(BLURIMG)*0.3
CONTOUR,BLURIMG,FINDGEN(100),FINDGEN(100),LEVELS=BLURLEV,PATH_INFO=BLURINFO,$
        PATH_XY=BLURXY,/CLOSED,/PATH_DATA_COORDS
; Go through each contour level ...
FBLUR=BLURLEV*0.
FOR I=0,N_ELEMENTS(BLURLEV)-1 DO BEGIN
; Grab all contours associated with this level
   ILEV=WHERE(BLURINFO.LEVEL EQ I)
; Assume that the biggest contour in this group is the one of interest
; (a more accurate thing to do is to get the contour that encloses the
; point of interest, like the center of the image)
   IBIG=WHERE(BLURINFO(ILEV).N EQ MAX(BLURINFO(ILEV).N))
   IBIG=ILEV(IBIG(0))
; Get the starting and ending point for this single contour
   ISTART=BLURINFO(IBIG).OFFSET
   IF IBIG EQ N_ELEMENTS(BLURLEV)-1 THEN IEND=ISTART ELSE IEND=BLURINFO(IBIG+1).OFFSET
; Extract the x/y coordinates for this single contour
   XC=REFORM(BLURXY(0,ISTART:IEND))
   YC=REFORM(BLURXY(1,ISTART:IEND))
; Get all the pixels within this contour and add up the total light
   IF N_ELEMENTS(XC) EQ 1 THEN BEGIN
      FBLUR(I)=BLURIMG(XC,YC)
   ENDIF ELSE BEGIN
      INSIDE=POLYFILLV([XC,XC(0)],[YC,YC(0)],100,100)
      FBLUR(I)=TOTAL(BLURIMG(INSIDE))
   ENDELSE
ENDFOR

; Now "plot up" the growth curve and figure out the contour that most
; closely matches 50% encircled energy
;
; Since we have such finely-gridded levels, let's just take
; whatever is closed to 50% rather than do some kind of fancy
; interpolation
IPICK=WHERE( ABS(FBLUR-TOTAL(BLURIMG)/2.) EQ MIN(ABS(FBLUR-TOTAL(BLURIMG)/2.)))
IPICK=IPICK(0)

; Now get the x,y coordinates for this contour
ILEV=WHERE(BLURINFO.LEVEL EQ IPICK)
IBIG=WHERE(BLURINFO(ILEV).N EQ MAX(BLURINFO(ILEV).N))
IBIG=ILEV(IBIG(0))
ISTART=BLURINFO(IBIG).OFFSET
IF IBIG EQ N_ELEMENTS(BLURLEV)-1 THEN IEND=ISTART ELSE IEND=BLURINFO(IBIG+1).OFFSET
BLURX=REFORM(BLURXY(0,ISTART:IEND))
BLURY=REFORM(BLURXY(1,ISTART:IEND))

PLOTS,BLURX,BLURY,LINEST=0,THICK=1,color=100

DISPLAY,'fake_galaxy_orig.fits',0,/ADJ,/NOC,/NOZ,FLD=1,MIND=MIND,MAXD=MAXD
DISPLAY,'fake_galaxy_convolved.fits',0,/ADJ,/NOC,/NOZ,FLD=2,MIND=MIND,MAXD=MAXD
DISPLAY,'fake_galaxy_blur.fits',0,/ADJ,/NOC,/NOZ,FLD=3,MIND=MIND,MAXD=MAXD

; Conclusion:  the "smeared" galaxy constructed by a series of point
; sources and the convolved image are nearly identical !  Cool,
; verifies what I thought!

STOP,'Z'

END
  
