PRO TEST_XRAYPETROSIAN,DEVAUC=DEVAUC,TRUNCATE=TRUNCATE

; After a lot of discussion, literature-searching and derivations,
; Mehmet and I have developed an algorithm to compute x-ray halo sizes
; that is based on the optical Petrosian radius.  Instead of using
; concentric elliptical apertures, our method uses nested isophotal
; contours. Therefore, we can't reference "radii" but rather
; "isophotal levels" which then gets converted into an effective
; radius later, by converting the isophotal area into a circular area
; and cmputing the radius of that circle that has the same area as the
; contour.

; This program generates some fake data to test the final algorithm.

; Make a fake perfectly elliptical-shaped galaxy that is
; well-resolved to represent the closest galaxy. Then make the same
; galaxy, but less resolved, to represent the most distant galaxy
; (but having the same light profile, etc.).  See how well the
; end-to-end algorithm computes the sizes of these 2 galaxies.  The
; sizes should end up to be the same value.

; Generate an elliptical galaxy having semi-minor/semi-major = 0.6,
; with a position angle of 30 degrees (relative to +y axis, eg lay to
; left of +y axis), a total luminosity of 1000, and an effective radius
; of 100.

IMSIZE=2999                   ; size of square image, along a side
;;RATIO=1./0.6                  ; major-to-minor ratio
RATIO=1.0   ; circle
BETA=0.52
POS_ANG=30.                   ; in units of degrees
XC=ROUND(IMSIZE/2.)           ; center of the image
YC=XC
DISTNEAR=10.                  ; Distance to galaxy in Mpc
DNEARTXT='10'

RMAX=ROUND(IMSIZE*SQRT(2.0))
R=FINDGEN(RMAX)

R50NEAR=100.
RENEAR=100.
IE=1.0    ; when r=re, intensity will equal IE, which equals in this case "1"
IDEVAUC=IE*EXP(-7.669*( (R/RENEAR)^0.25 - 1.0) )
RCUTOFF=R50NEAR*1.2
IN=WHERE(R GT RCUTOFF)
IDEVAUC(IN)=0.0

; I have derived the following relationship, and tested it in a
; separate program called growthcurve.pro:
;      r50 = rc sqrt( (1/2)^(1/(B+1) -1) where B=-3*beta+1/2
; Invert to compute the RC, if R50 is the above
;      rc = r50/sqrt( (0.5^(1/(b+1))-1)
RCNEAR=R50NEAR/SQRT( 0.5^(1.0/(-3.0*BETA+0.5+1.0)) -1.0 ) 
; Force the intensity to equal 1 when the radius equals the virial
; radius, meaning that I0 must equal:
;  I0=1.0/(2.0)^(-3.0*BETA+0.5)
; But all we have right now is the core radius, so need to determined
; the virial radius.  To both determine where to truncate the light
; and the radius out to which to integrate to impose flux-conservation, need a relationship
; between Rc and Rvirial.  Let's use the one provided in
; equation 117 at
; https://ned.ipac.caltech.edu/level5/Sept03/Peacock/Peacock6_5.html
;     DeltaC=200c^3/3/( ln(1+c) - c/(1+c) ) 
; where c is ratio of virial to core radius, and DeltaC is set to equation 118:
;     DeltaC =~ 3000(1+zc)^3   where zc=collapse redshift, which is
;     zero for massive systems, so DeltaC is about 3000
; I numerically deteremined the value of c that would satisfy the
; first equation above, if DeltaC were 3000.  And by "numerically", I
; mean that I literally starting picking values for "c" out of the hat
; until I could get the left hand side of that equation to match
; 3000. That value turns out to be about 3.085.  R_virial/R_core ~3.085
; From yet another reference (Figure 3 from
; https://ned.ipac.caltech.edu/level5/Sept03/Mulcheay/Muchaey3_3.html), we
; see that the X-ray halo size is about equal to 0.7 of the virial
; radius, as a median value (just eyeballing ... perhaps even a bit
; smaller, but we'll go with the 0.7). So (whew!):
CORE2VIRIAL=3.085
RVIR=RCNEAR*CORE2VIRIAL
I0=1.0/(1.0 + (RVIR/RCNEAR)^2)^(-3.0*BETA+0.5) 
IBETA=I0*(1.0 + (R/RCNEAR)^2 )^(-3.0*BETA+0.5)
;;RCUTOFF=RVIR*0.75
RCUTOFF=R50NEAR*1.2
;;;;; NOTE:  the whole cutting the thing at the virial radius
;;;;; didn't seem to work so well.  Virial radius seemed WAYYY
;;;;; too small compared to the R50.  So a question is, do people
;;;;; generally detect xray emmission out to at least the half-light
;;;;; radius? 
IN=WHERE(R GT RCUTOFF)
IBETA(IN)=0.0

IF KEYWORD_SET(DEVAUC) THEN F=IDEVAUC ELSE F=IBETA

MAKEGAL,IMSIZE,RATIO,POS_ANG,NEARGAL,XC=XC,YC=YC,FLUX=F
IZERO=WHERE(NEARGAL EQ 0.,NZERO)
IF NZERO EQ 0 THEN BEGIN
; find all pixels within the effective radius (devauc) or virial radius (beta)
   INNER1=WHERE(NEARGAL GE 1.0,NNEAR)
ENDIF ELSE INNER1=WHERE(NEARGAL GT 0.0,NNEAR) 
NEARINNER=TOTAL(NEARGAL(INNER1))
NEARPROF=NEARGAL(XC,*)
WRITEFITS,'galaxy-nonoise-'+DNEARTXT+'mpc.fits',NEARGAL
;DISPLAY,'galaxy-nonoise-'+DNEARTXT+'mpc.fits',0,/NOC,/NOZ,FLD=1,/REV,/ADJ,TITLE='no noise '+DNEARTXT+' mpc'

;-----------------------------------------------------------------------------------
;-----------------------------------------------------------------------------------
; Now generate the same galaxy, but at 30 mpc away (3 times further
; away)
DISTMID=30.
DMIDTXT='30'

REMID=RENEAR*(DISTNEAR/DISTMID)
IE=1.0    ; when r=re, intensity will equal IE, which equals in this case "1"
IDEVAUC=IE*EXP(-7.669*( (R/REMID)^0.25 - 1.0) )
; Truncate the galaxy at some distance to replicate the x-ray images
RCUTOFF=REMID*1.2
IN=WHERE(R GT RCUTOFF)
IDEVAUC(IN)=0.0

RCMID=RCNEAR*(DISTNEAR/DISTMID)
R50MID=R50NEAR*(DISTNEAR/DISTMID)
IBETA=I0*(1.0 + (R/RCMID)^2 )^(-3.0*BETA+0.5)
RVIR=RCMID*CORE2VIRIAL
;;RCUTOFF=RVIR*0.75
RCUTOFF=R50MID*1.2
IN=WHERE(R GT RCUTOFF)
IBETA(IN)=0.0

IF KEYWORD_SET(DEVAUC) THEN F=IDEVAUC ELSE F=IBETA

MAKEGAL,IMSIZE,RATIO,POS_ANG,MIDGAL,XC=XC,YC=YC,FLUX=F
IZERO=WHERE(MIDGAL EQ 0.,NZERO)
IF NZERO EQ 0 THEN BEGIN
   INNER2=WHERE(MIDGAL GE 1.0,NMID) ; pixels within effective radius
ENDIF ELSE INNER2=WHERE(MIDGAL GT 0.0,NMID)
MIDINNER=TOTAL(MIDGAL(INNER2))
; Note:  I spent nearly all day on this code.  Realized that if you
; just take the image as chucked out by MAKEGAL, the surface
; brightness is inherently the same per pixel for each of the 2 galaxy
; distances (e.g., the central pixel has same surface brightness, for
; example).  Also proved that the area within 1 effective radius for
; the near galaxy is a factor of (DISTFAR/DISTNEAR)^2 times more that
; the area within 1 effective radius for the more distant galaxy.
;
; But here's where I am confused:  the total flux within 1
; effective radius for the near galaxy should be larger by 100 as
; compared to that of the more distant galaxy, but its only about 25
; times larger.  I can't for the life of me understand why. The
; flux enclosed should scale as (distance 1/ distance 2)^2, right???? 

; Next morning thoughts:  I think that the problem is that makegal
; samples the de vauc profile rather than allowing for flux
; conservation.  In other words, when the near galaxy light profile is
; sqished into the lower resolution more distant galaxy version, what
; happens to all the flux in between the values corresponding to radii
; that happen to "land" on the fewer pixels of the more distant
; galaxy? I think they are ignored, so flux is lost.  Not convinced
; that the surface brightness profile in the center should even be the
; same as that of the near galaxy's centeral pixel.

; Force the average surface brightness w/i effective radius of near galaxy to
; equal the surface brightness w/i effective radius of far
; galaxy. YEAH, WORKS!! The integrated fluxes compare as expected,
; with difference due to distance^2 effect. 
SBNEAR=NEARINNER/FLOAT(NNEAR)
SBMID=MIDINNER/FLOAT(NMID)
MIDGAL=(MIDGAL/MIDINNER)*NEARINNER*(DISTNEAR/DISTMID)^2

PRINT,TOTAL(MIDGAL(INNER2)),TOTAL(NEARGAL(INNER1))
PRINT,TOTAL(MIDGAL(INNER2))/FLOAT(NMID), TOTAL(NEARGAL(INNER1))/FLOAT(NNEAR)
MIDPROF=MIDGAL(XC,*)
WRITEFITS,'galaxy-nonoise-'+DMIDTXT+'mpc.fits',MIDGAL
;DISPLAY,'galaxy-nonoise-'+DMIDTXT+'mpc.fits',0,/NOC,/NOZ,FLD=2,/REV,/ADJ,TITLE='no noise '+DMIDTXT+'mpc'

;-----------------------------------------------------------------------------------
;-----------------------------------------------------------------------------------
; Now generate the same galaxy, but at 80 mpc away (8 times further away)
DISTFAR=80.
DFARTXT='80'

REFAR=RENEAR*(DISTNEAR/DISTFAR)
IE=1.0    ; when r=re, intensity will equal IE, which equals in this case "1"
IDEVAUC=IE*EXP(-7.669*( (R/REFAR)^0.25 - 1.0) )
RCUTOFF=REFAR*1.2
IN=WHERE(R GT REFAR)
IDEVAUC(IN)=0.0

R50FAR=R50NEAR*(DISTNEAR/DISTFAR)
RCFAR=RCNEAR*(DISTNEAR/DISTFAR)
IBETA=I0*(1.0 + (R/RCFAR)^2 )^(-3.0*BETA+0.5)
RVIR=RCFAR*CORE2VIRIAL
;;RCUTOFF=RVIR*0.75
RCUTOFF=R50FAR*1.2
IN=WHERE(R GT RCUTOFF)
IBETA(IN)=0.0

IF KEYWORD_SET(DEVAUC) THEN F=IDEVAUC ELSE F=IBETA

MAKEGAL,IMSIZE,RATIO,POS_ANG,FARGAL,XC=XC,YC=YC,FLUX=F

IZERO=WHERE(FARGAL EQ 0.,NZERO)
IF NZERO EQ 0 THEN BEGIN
   INNER3=WHERE(FARGAL GE 1.0,NFAR) ; pixels within effective radius
ENDIF ELSE INNER3=WHERE(FARGAL GT 0.0,NFAR)
FARINNER=TOTAL(FARGAL(INNER3))
SBNEAR=NEARINNER/FLOAT(NNEAR)
SBFAR=FARINNER/FLOAT(NFAR)
FARGAL=(FARGAL/FARINNER)*NEARINNER*(DISTNEAR/DISTFAR)^2

PRINT,TOTAL(FARGAL(INNER3)),TOTAL(NEARGAL(INNER1))
PRINT,TOTAL(FARGAL(INNER3))/FLOAT(NFAR), TOTAL(NEARGAL(INNER1))/FLOAT(NNEAR)
FARPROF=FARGAL(XC,*)
WRITEFITS,'galaxy-nonoise-'+DFARTXT+'mpc.fits',FARGAL
;DISPLAY,'galaxy-nonoise-'+DFARTXT+'mpc.fits',0,/NOC,/NOZ,FLD=2,/REV,/ADJ,TITLE='no noise '+DFARTXT+'mpc'

;*************************************************************************************
;*************************************************************************************

; Now inject some random noise into the images
NOISE=RANDOMU(SEED,IMSIZE,IMSIZE)
; Force the median value of the noise to be such that the photometry
; within the effective radius of the NEAR galaxy would be 1%.
; Recall that error_phot = Area_aperture x background_noise, so
; background_noise = error_phot/Area_aperture
NOISE=NOISE/MEDIAN(NOISE)*NEARINNER*0.01/FLOAT(NNEAR)
WRITEFITS,'noise.fits',NOISE
;DISPLAY,'noise.fits',0,/NOC,/NOZ,FLD=3,/REV,/ADJ,TITLE='noise'

NEARGALNOISE=NEARGAL
; Now compute poisson noise for the pixels that have detected light
IN=WHERE(NEARGAL GT 0,N)
FOR I=0,N-1 DO BEGIN
   POIS=RANDOMU(SEED,POISSON=NEARGAL(IN(I)))
   NEARGALNOISE(IN(I))=NEARGAL(IN(I))+POIS
ENDFOR
IF KEYWORD_SET(DEVAUC) THEN NEARGALNOISE=NEARGALNOISE+NOISE
; Trim up the image so that it is smaller and less blank sky. Might
; also need to add some keywords into the fits header
NEARGALNOISE=NEARGALNOISE(XC-250:XC+250,YC-250:YC+250)

WRITEFITS,'galaxy-5noise-'+DNEARTXT+'mpc.fits',NEARGALNOISE
;DISPLAY,'galaxy-5noise-'+DNEARTXT+'mpc.fits',0,/NOC,/NOZ,FLD=4,$
;        /REV,/ADJ,TITLE='with noise '+DNEARTXT+'mpc',MIND=-0.001,MAXD=0.2

; Add the same noise array to the MID-distant galaxy, and make no
; adjustment to the noise amplitude
MIDGALNOISE=MIDGAL
IN=WHERE(MIDGAL GT 0.,N)
FOR I=0,N-1 DO BEGIN
   POIS=RANDOMU(SEED,POISSON=MIDGAL(IN(I)))
   MIDGALNOISE(IN(I))=MIDGAL(IN(I))+POIS
ENDFOR
IF KEYWORD_SET(DEVAUC) THEN MIDGALNOISE=MIDGALNOISE+NOISE
MIDGALNOISE=MIDGALNOISE(XC-250:XC+250,YC-250:YC+250)
WRITEFITS,'galaxy-5noise-'+DMIDTXT+'mpc.fits',MIDGALNOISE
DISPLAY,'galaxy-5noise-'+DMIDTXT+'mpc.fits',0,/NOC,/NOZ,FLD=5,$
        /REV,/ADJ,TITLE='with noise '+DMIDTXT+' mpc',MIND=-0.001,MAXD=0.2

; Add the same noise array to the distant galaxy, and make no
; adjustment to the noise amplitude
FARGALNOISE=FARGAL
IN=WHERE(FARGAL GT 0.,N)
FOR I=0,N-1 DO BEGIN
   POIS=RANDOMU(SEED,POISSON=FARGAL(IN(I)))
   FARGALNOISE(IN(I))=FARGAL(IN(I))+POIS
ENDFOR
IF KEYWORD_SET(DEVAUC) THEN FARGALNOISE=FARGALNOISE+NOISE
FARGALNOISE=FARGALNOISE(XC-250:XC+250,YC-250:YC+250)
WRITEFITS,'galaxy-5noise-'+DFARTXT+'mpc.fits',FARGALNOISE
DISPLAY,'galaxy-5noise-'+DFARTXT+'mpc.fits',0,/NOC,/NOZ,FLD=6,$
        /REV,/ADJ,TITLE='with noise '+DFARTXT+' mpc',MIND=-0.001,MAXD=0.2

NEARGAL=NEARGAL(XC-250:XC+250,YC-250:YC+250)
MIDGAL=MIDGAL(XC-250:XC+250,YC-250:YC+250)
FARGAL=FARGAL(XC-250:XC+250,YC-250:YC+250)

WRITEFITS,'galaxy-nonoise-'+DNEARTXT+'mpc.fits',NEARGAL
WRITEFITS,'galaxy-nonoise-'+DMIDTXT+'mpc.fits',MIDGAL
WRITEFITS,'galaxy-nonoise-'+DFARTXT+'mpc.fits',FARGAL

stop,'1'

END
