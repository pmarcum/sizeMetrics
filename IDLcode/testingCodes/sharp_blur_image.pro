PRO SHARP_BLUR_IMAGE,CENTERDIAM=CENTERDIAM,CENTERFRAC=CENTERFRAC,$
                     DIAGTHICK=DIAGTHICK,DIAGFRAC=DIAGFRAC,RANDLENGTH=RANDLENGTH
; Script to first test idea that a smaller concentric polygon can be
; created by scaling the radius of each polygon vertice, relative to
; an origin (for example, to the centroid of the shape).  If that
; approach does not work, then try subtracting off a fixed margin from
; the radius, where the amount of margin subtracted off depends on the
; angle of intersection between the line drawn from reference point to
; the vertex and the orientation of that line to the "surface" of the
; polygon's shape at that point (e.g., perpendicular would be
; full margin, orthogonal would mean that nothing gets subtracted
; off).  Probably need to recenter the re-gridded polygon in this
; case.
;
; The next thing that is done, after the above test is performed, is
; to produce a simualted image that is an extended object generated by
; summing lots of randomly-placed (and randomly-strengthed) gaussian
; functions.  Generate the isophote at some level.  Then apply
; gaussian smearing to the image, see if the adjusted isophote from
; unaltered image aligns where it should in the smeared image, given
; your knowledge of how much smearing occurred.

; Test to see if one can scale the radii from some reference point to each
; polygon vertice to make a smaller concentric polygon, or if this
; process is indeed a subtractive one instead.


;*************************************************************
; Now generate fake data and test the algorithm of the trimmer
;*************************************************************
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

FRAC=0.5  ; 50% encircled energy.  Change to 0.8 to study 80% encircled energy

; for plotting purposes, generate a user-defined circular plotting symbol
TH=FINDGEN(41)/40.*2.*!PI  ; make an array of angles defining a circle
XSYM=SIN(TH) & YSYM=COS(TH)  ; radius is implicitly unity
USERSYM,XSYM,YSYM,/FILL

NPOINTS=1000.
WINDOW,4,XSIZE=800,YSIZE=800
PLOT,FINDGEN(100),FINDGEN(100),XRA=[0,99],YRA=[0,99],/XST,/YST,/ISO,/NODATA
; Get inner 30% of points
IF N_ELEMENTS(CENTERDIAM) EQ 0 THEN CENTERDIAM=40.
IF N_ELEMENTS(CENTERFRAC) EQ 0 THEN CENTERFRAC=0.3
CENTERDIAM=40.  ; number of pixels of diameter for central concentration
RADCENTER=(RANDOMU(SEED,CENTERFRAC*NPOINTS)-0.5)*CENTERDIAM
ANGCENTER=(RANDOMU(SEED,CENTERFRAC*NPOINTS)-0.5)*!PI

; get diagonal 40% of points
IF N_ELEMENTS(DIAGTHICK) EQ 0 THEN DIAGTHICK=35.   ; number of degrees in thickness for diagonal strip
DIAGLENGTH=80.                  ; maximum length from end to end along diagonal
DIAGPOSANG=135.                 ; desired position angle of diagonal feature
IF N_ELEMENTS(DIAGFRAC) EQ 0 THEN DIAGFRAC=0.4
RADDIAG=(RANDOMU(SEED,NPOINTS*DIAGFRAC)-0.5)*DIAGLENGTH
ANGDIAG=((RANDOMU(SEED,NPOINTS*DIAGFRAC)-0.5)*DIAGTHICK+DIAGPOSANG)*!PI/180.

; get the 30% random sprinkle
IF N_ELEMENTS(RANDLENGTH) EQ 0 THEN RANDLENGTH=100.
; desired extent in x (or y) direction for random sprinkling
RANDFRAC=1.0-DIAGFRAC-CENTERFRAC
RADRAND=(RANDOMU(SEED,NPOINTS*RANDFRAC)-0.5)*RANDLENGTH
ANGRAND=(RANDOMU(SEED,NPOINTS*RANDFRAC)-0.5)*!PI

; Consolidate
RADIUS=[RADCENTER,RADDIAG,RADRAND]
ANGLE=[ANGCENTER,ANGDIAG,ANGRAND]
AMP=RANDOMU(SEED,1000)*100.      ; amplitude of gaussian

IMAGE=FLTARR(100,100)

; Determine what the total energy is for each of these point sources,
; given the amplitude.  Be lazy and just do this computation
; numerically:
ENERGY=AMP*0.
INITPSF=1.0  ; fwhm
INITSIG=INITPSF/(2.*SQRT(2.*ALOG(2)))
FOR I=0,N_ELEMENTS(AMP)-1 DO BEGIN
; Make the psf be positioned in the middle of the frame, to optimally
; extract the total energy without having to worry about effects of
; offsets of centroid:
   DISTANCE=100./2.
   PSFIMG=GAUSSIAN_FUNCTION([INITSIG,INITSIG],WIDTH=DISTANCE,MAXIMUM=AMP(I))
   ENERGY(I)=TOTAL(PSFIMG)
ENDFOR

; Turn into x and y coordinates
X=(RADIUS*COS(ANGLE)+49.)>1.
Y=(RADIUS*SIN(ANGLE)+49.)>1.

; pixelate by forcing rounded numbers
X=ROUND(X)<99.
Y=ROUND(Y)<99.

IMAGE(X,Y)=AMP
WRITEFITS,'fake_galaxy_fwhm0.fits',IMAGE

PLOTS,X,Y,PSYM=8

;
; Without any smearing, figure out the shape of the polygon that would
; enclose 50% of the total energy of the image. Perform this task by
; identifying the pixels above some threshold, determining if their
; sum is equivalent to 50% of the total energy, and seeing where the
; cross-over point is (when one more point is added that is just a
; little further away from the radius, the 50% is exceeded).  Use
; convex hull to bound the points within 50%. Sort in ascending radius
; from the "object" of interest (positioned in the center of the
; frame), so that the enclosure insures that the innermost points are
; included first.
RADIUS=SQRT( (X-49.)^2 + (Y-49.)^2 )
ISORT=SORT(RADIUS)
ETOT=0.
XIN=[-1.]
YIN=[-1.]
I=0
WHILE ETOT LT TOTAL(ENERGY)*FRAC DO BEGIN
   ETOT=ETOT+ENERGY(ISORT(I))
   XIN=[XIN,X(ISORT(I))]
   YIN=[YIN,Y(ISORT(I))]
   I=I+1
ENDWHILE
N=N_ELEMENTS(XIN)
XIN=XIN(1:N-1)
YIN=YIN(1:N-1)

; The question is now, how much over the half-way mark did that last
; point push the total?  If removing it would make the total closer to
; the half-energy, then just remove it
I=I-1
IF ABS(ETOT-TOTAL(ENERGY)*FRAC) LT ABS(ETOT-ENERGY(ISORT(I))-TOTAL(ENERGY)*FRAC) THEN BEGIN
   N=N_ELEMENTS(XIN)
   XIN(N-1)=-1
   IN=WHERE(XIN NE -1)
   XIN=XIN(IN) & YIN=YIN(IN)
ENDIF

; Now make a "contour" around these points by using concave hull
; method:
TRIANGULATE,XIN,YIN,TRIANGLES,HULL
XHULL=XIN(HULL)
YHULL=YIN(HULL)

PLOTS,[XHULL,XHULL(0)],[YHULL,YHULL(0)],LINEST=0,THICK=3,COLOR=100

; Now turn each coordinate into a 2-D gaussian and form an image
; let the FWHM of the gaussian be FWHM(i)

FWHM=[0.,1.,3.,5.,8.,10.,12.,15.,20.,25.]
CVALUE=FWHM*0.
XCONT=[-1.]
YCONT=[-1.]
FWHMCONT=[-1.]

FOR F=0,N_ELEMENTS(FWHM)-1 DO BEGIN ; for each desired image PSF
   SIG=FWHM(F)/(2.*SQRT(2.*ALOG(2.)))
   IF SIG EQ 0 THEN GOTO,SKIPSMOOTH
   
   IMAGE=FLTARR(100,100)   
   FOR I=0,N_ELEMENTS(X)-1 DO BEGIN  ; for each point producing light in the image
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
      PSFIMG=GAUSSIAN_FUNCTION([SIG,SIG],WIDTH=DISTANCE,MAXIMUM=AMP(I))
; Force the PSF to have a total energy that is equivalent to
; ENERGY(i), computed using a FWHM=1, so that there is an overall
; conservation of energy here
      PSFIMG=PSFIMG*ENERGY(I)/TOTAL(PSFIMG)
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
   WRITEFITS,'fake_galaxy_fwhm'+STRTRIM(STRING(FIX(FWHM(F))),2)+'.fits',IMAGE

SKIPSMOOTH:
   FIND_ENCLOSED_CONTOUR,IMAGE,CX,CY,CVAL,FRAC=FRAC
   XCONT=[XCONT,CX]
   YCONT=[YCONT,CY]
   FWHMCONT=[FWHMCONT,CX*0.+FWHM(F)]
; While we are at it, save the contour level value that fit the
; enclosed fraction criteria, as well
   CVALUE(F)=CVAL
ENDFOR
N=N_ELEMENTS(XCONT)
XCONT=XCONT(1:N-1)
YCONT=YCONT(1:N-1)
FWHMCONT=FWHMCONT(1:N-1)

; Plot up the percent-enclosed contours from the various smoothing
; levels to visually see how they vary
FOR I=0,N_ELEMENTS(FWHM)-1 DO BEGIN
   IN=WHERE(FWHMCONT EQ FWHM(I))
   OPLOT,XCONT(IN),YCONT(IN),LINEST=0
   PRINT,'FWHM ',FWHM(I)
   ANS='' & READ,'Hit a key to continue ',ANS
; Construct an image that consists only of the contour that was picked
; out
   INSIDE=POLYFILLV([XCONT(IN),XCONT(IN(0))],[YCONT(IN),YCONT(IN(0))],100,100)
   IMAGE=IMAGE*0. & IMAGE(INSIDE)=1
   WRITEFITS,'fake_contour_fwhm'+STRTRIM(STRING(FIX(FWHM(I))),2)+'.fits',IMAGE
ENDFOR

WINDOW,5
PLOT,FWHM,CVALUE,PSYM=4

;DO_PSF_RESIZE,BLURX,BLURY,EDGE,XTRIM,YTRIM

STOP,'LAST STOP!!!'

END
  
