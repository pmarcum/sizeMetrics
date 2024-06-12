PRO GET_PETRO,IMGFILE,PERCENT=PERCENT,PTARGET=PTARGET
; Determines the 50% encircled energy surface brightness level,
; it's corresponsing areal circularized radius, and a petrosian
; radius by determining the surface brightness contour at which the
; enclosed energy divided by the surface brightness of the contour
; equals a specific target number (called the Petrosian index).

rcmodel=10.    ; the core radius used to create the image being processed, in arcsec
pixscale=0.492 ; the pixel scale (number arcsec on side) of image, to convert rcmodel to pix units  
  
IF N_ELEMENTS(PTARGET) EQ 0 THEN PTARGET=13.122
IF N_ELEMENTS(PERCENT) EQ 0 THEN PERCENT=80.
; ptarget=13.122 for percent=80.; 3.62 for percent=50 for beta=0.6
IF N_ELEMENTS(BLANK) EQ 0 THEN BLANK=-9999
IF N_ELEMENTS(IMGFILE) EQ 0 THEN IMGFILE='marxsims/disteffects1_final.fits'

FXREAD,IMGFILE,IMG,HDR

IGOOD=WHERE(IMG GT 0)
MED=MEDIAN(IMG(IGOOD))
ST=STDEV(IMG(IGOOD))
MINVAL=MIN(IMG(IGOOD))
MAXVAL=MAX(IMG(IGOOD))

;  ------------------------    max value
; -------------------------   Q9 (median of pix above Q8 - 96.875th percentile)
;
; --------------------------  Q8 (median of pix above Q7 - 93.75th percentile)
; --------------------------  Q7 (median of pix above Q6 - 87.5th percentile)
;
;
; --------------------------  Q6 (median of pix above median value - 75th percentile)
;
; --------------------------   median value  (50th percentile)
;
;
;
; --------------------------   Q4 (median of the pixs below median value - 25th percentile) 
;
;
; --------------------------   Q3 (median if the pixs below Q2 - 12.5th percentile)
;
;
; --------------------------   Q2 (median of pixs below Q3 - 6.25th percentile)
;
; --------------------------   Q1 (median of pixs below Q2 - 3.125th percentile)
;
; --------------------------   Q0 (median of pix below Q1 -- 1.5625 percentile
;
; --------------------------   min value

IN=WHERE(IMG(IGOOD) LT MED)
Q4=MEDIAN(IMG(IGOOD(IN)))  ; the median value of the pixels below the overall median

IN=WHERE(IMG(IGOOD) LT Q4)
Q3=MEDIAN(IMG(IGOOD(IN)))  ; the median value of the pixels below Q4

IN=WHERE(IMG(IGOOD) LT Q3)
Q2=MEDIAN(IMG(IGOOD(IN)))   ; the median values of pixels below Q3

IN=WHERE(IMG(IGOOD) LT Q2)
Q1=MEDIAN(IMG(IGOOD(IN)))

IN=WHERE(IMG(IGOOD) LT Q1)
Q0=MEDIAN(IMG(IGOOD(IN)))

IN=WHERE(IMG(IGOOD) GT MED)
Q6=MEDIAN(IMG(IGOOD(IN)))

IN=WHERE(IMG(IGOOD) GT Q6)
Q7=MEDIAN(IMG(IGOOD(IN)))

IN=WHERE(IMG(IGOOD) GT Q7)
Q8=MEDIAN(IMG(IGOOD(IN)))

IN=WHERE(IMG(IGOOD) GT Q8)
Q9=MEDIAN(IMG(IGOOD(IN)))

TMPLEV=[Q0,Q1,Q2,Q3,Q4,MED,Q6,Q7,Q8,Q9]
TMPPERC=[1.5625,3.125,6.25,12.5,25.0,50.0,75.0,87.5,93.75,96.875]

; Want to get nice even percentile numbers like 10,20,30,40,50,60th
; percentile, so interpolate the PERCENTILE vs LEVELS function to get
; these nicer numbers
PERCENTILE=[1.7,2.,2.5,3.0,3.5,4.0,4.5,5.0,6.0,$
            7.0,8.0,9.0,10.,15.,20.,25.,30.,35.,$
            40.,45.,50.,55.,60.,65.,70.,75.,80.,$
            85.,90.,95.]
LEVELS=INTERPOL(TMPLEV,TMPPERC,PERCENTILE)

PERCENTILE=[TMPPERC(0),PERCENTILE,TMPPERC(N_ELEMENTS(TMPLEV)-1)]
LEVELS=[TMPLEV(0),LEVELS,TMPLEV(N_ELEMENTS(TMPLEV)-1)]

CONTOUR,IMG,LEVELS=LEVELS,PATH_INFO=PINFO,PATH_XY=PXY,/CLOSED,/PATH_DATA_COORDS

; Now go through each contour, construct growth curve and a petrosian curve

AREA=LEVELS*0.
FLUX=LEVELS*0.
IMGSZ=SIZE(IMG)
XCENTER=LEVELS*0.
YCENTER=LEVELS*0.
FOR I=0,N_ELEMENTS(LEVELS)-1 DO BEGIN
; grab all contours associated with this ith level
   ILEV=WHERE(PINFO.LEVEL EQ I,NCONTOURS)
   AREATEMP=0
   FOR J=0,NCONTOURS-1 DO BEGIN
; get starting and stopping point for the contours
      ISTART=PINFO(ILEV(J)).OFFSET
      XC=REFORM(PXY(0,ISTART:ISTART+PINFO(ILEV(J)).N-1))
      YC=REFORM(PXY(1,ISTART:ISTART+PINFO(ILEV(J)).N-1))
; Now get area enclosed by contour, add to cumulative for this contour level
      IF N_ELEMENTS(XC) EQ 1 THEN BEGIN
         IF IMG(XC,YC) NE BLANK THEN BEGIN
            AREA(I)=AREA(I)+1
            FLUX(I)=FLUX(I)+IMG(XC,YC)
            IF AREATEMP EQ 0 THEN BEGIN
               AREATEMP=1
               XCENTER(I)=XC
               YCENTER(I)=YC
            ENDIF
         ENDIF 
      ENDIF ELSE BEGIN
         INSIDE=POLYFILLV([XC,XC(0)],[YC,YC(0)],IMGSZ(1),IMGSZ(2))
         FROM2TO1D,XC,YC,IMGSZ(1),BORDER
; comine the inside pixels with the border of the contour.  There may
; be some overlap between the 2 groups, so need to exclude any
; redundancies and just get unique pixel coordinates
         INSIDE=[INSIDE,BORDER]
         INSIDE=INSIDE(SORT(INSIDE))
         INSIDE=INSIDE(UNIQ(INSIDE))
         IN=WHERE(IMG(INSIDE) NE BLANK,NOK)
         IF NOK NE 0 THEN BEGIN
            AREA(I)=AREA(I)+NOK
            FLUX(I)=FLUX(I)+TOTAL(IMG(INSIDE(IN)))
            IF NOK GT AREATEMP THEN BEGIN
               AREATEMP=NOK
               FROM1TO2D,INSIDE,IMGSZ(1),XX,YY
               XCENTER(I)=MEAN(XX)
               YCENTER(I)=MEAN(YY)
            ENDIF
         ENDIF
      ENDELSE
   ENDFOR
ENDFOR

; Determine the centroid of galaxy by averaging the centroids of the
; largest contour associated with each level:
GET_STATS,XCENTER,XCENTROID,/AVERAGE,SIGCUT=3.0
GET_STATS,YCENTER,YCENTROID,/AVERAGE,SIGCUT=3.0
; now that we have a centroid, we can go back and figure out the
; average radius for the largest contour of each level.
RADIUS=LEVELS*0.
FOR I=0,N_ELEMENTS(LEVELS)-1 DO BEGIN
; grab all contours associated with this ith level
   ILEV=WHERE(PINFO.LEVEL EQ I,NCONTOURS)
   AREATEMP=0
   FOR J=0,NCONTOURS-1 DO BEGIN
; get starting and stopping point for the contours
      ISTART=PINFO(ILEV(J)).OFFSET
      XC=REFORM(PXY(0,ISTART:ISTART+PINFO(ILEV(J)).N-1))
      YC=REFORM(PXY(1,ISTART:ISTART+PINFO(ILEV(J)).N-1))
; Now get area enclosed by contour, add to cumulative for this contour level
      IF N_ELEMENTS(XC) EQ 1 THEN BEGIN
         IF IMG(XC,YC) NE BLANK THEN BEGIN
            IF AREATEMP EQ 0 THEN BEGIN
               AREATEMP=1
               RADIUS(I)=SQRT( (XCENTROID-XC)^2 + (YCENTROID-YC)^2 )
            ENDIF
         ENDIF 
      ENDIF ELSE BEGIN
         INSIDE=POLYFILLV([XC,XC(0)],[YC,YC(0)],IMGSZ(1),IMGSZ(2))
         FROM2TO1D,XC,YC,IMGSZ(1),BORDER
; comine the inside pixels with the border of the contour.  There may
; be some overlap between the 2 groups, so need to exclude any
; redundancies and just get unique pixel coordinates
         INSIDE=[INSIDE,BORDER]
         INSIDE=INSIDE(SORT(INSIDE))
         INSIDE=INSIDE(UNIQ(INSIDE))
         IN=WHERE(IMG(INSIDE) NE BLANK,NOK)
         IF NOK NE 0 THEN BEGIN
            IF NOK GT AREATEMP THEN BEGIN
               AREATEMP=NOK
               FROM1TO2D,INSIDE,IMGSZ(1),XX,YY
               RADIUS(I)=MEAN(SQRT( (XCENTROID-XC)^2 + (YCENTROID-YC)^2 ))
            ENDIF
         ENDIF
      ENDELSE
   ENDFOR
ENDFOR

; below is a test to make sure that the surface brightness profile
; indeed follows the function (beta model with beta=0.52,
; rocre=RCMODEL/pixscale pixels) that was
; used in CIAO/MARX to make the image from which the surface
; brightness profile was derived. I found that there is a near-perfect
; match. 
WINDOW,1,TITLE='RADIUS VS SURFACE BRIGHTNESS'
PLOT,RADIUS,LEVELS
rcore=rcmodel/pixscale  ; pixel units
oplot,[3.*rcore,3.*rcore],[min(levels),max(levels)]
sb3rc=interpol(levels,radius,3.*rcore)
plots,3.*rcore,sb3rc,psym=4,symsi=2
beta=0.60
b=-3.0*beta+0.5
sb0=sb3rc/(10.0^b)
bfn=sb0*(1.0+(radius/rcore)^2)^b
oplot,radius,bfn,thick=2

; Note: major accomplishment here, proved that the theoretical beta
; model can be recovered, and fits almost exactly on image's profile! 

; Now that we've gathered up the total areas involved per
; contour level, the total flux enclosed, etc.  we can now develop the
; growth curve, compute the 50% encircled energy surface brightness
; value, develop the petrosian curve and compute the petrosian surface
; brightness value:

; GROWTH CURVE  (integrated flux within each contour as fn of surface brightness)
WINDOW,1,TITLE='GROWTH CURVE TO DETERMINE EFFECTIVE RADIUS'
PLOT,LEVELS,FLUX,XTITLE='Isophotal surface brightness',$
   YTITLE='Isophotal integrated flux',CHARSIZE=1.7

; Now do an interpolation to determine at what surface brightness the
; PERCENT% enclosed energy is:
TOTFLUX=MAX(FLUX)

SBPERCENT=INTERPOL(LEVELS,FLUX,TOTFLUX*PERCENT/100.)

;;totflux=-!pi*sb0*(rcmodel/pixscale)^2/(b+1.)
sb50=interpol(levels,flux,TOTflux*0.5)
sb80=interpol(levels,flux,TOTflux*0.8)
plots,[sb50,sb50],[min(flux),TOTflux],thick=0.7
plots,[sb80,sb80],[min(flux),TOTflux],thick=2.0

; do the same interpolation but in radius units rather than contour sb
; units:
rad50=interpol(radius,flux,TOTflux*0.5)
rad80=interpol(radius,flux,TOTflux*0.8)
; Now link the radius with the surface brightness at that radius, and
; compare to sb50 and sb80
sbrad50=interpol(levels,radius,rad50)
sbrad80=interpol(levels,radius,rad80)

; Now make a contour at this level only, to get the final properties
; of the associated contour
CONTOUR,IMG,LEVELS=SBPERCENT,PATH_INFO=PINFO,PATH_XY=PXY,/CLOSED,/PATH_DATA_COORDS

; Determine the total area, etc of the contours associated with the SBPERCENT contour
PIC=FLTARR(IMGSZ(1),IMGSZ(2))
AREAPERCENT=0.
FLUXPERCENT=0.
FOR J=0,N_ELEMENTS(PINFO.LEVEL)-1 DO BEGIN  ; could be multiple contours at SBPERCENT
   ISTART=PINFO(J).OFFSET
   XC=REFORM(PXY(0,ISTART:ISTART+PINFO(J).N-1))
   YC=REFORM(PXY(1,ISTART:ISTART+PINFO(J).N-1))
   IF N_ELEMENTS(XY) EQ 1 THEN BEGIN
      IF IMG(XC,YC) NE BLANK THEN BEGIN
         AREAPERCENT=AREAPERCENT+1.
         FLUXPERCENT=FLUXPERCENT+IMG(XC,YC)
         XX=XC
         YY=YC
      ENDIF
   ENDIF ELSE BEGIN
      INSIDE=POLYFILLV([XC,XC(0)],[YC,YC(0)],IMGSZ(1),IMGSZ(2))
      FROM2TO1D,XC,YC,IMGSZ(1),BORDER
      INSIDE=[INSIDE,BORDER]
      INSIDE=INSIDE(SORT(INSIDE))
      INSIDE=INSIDE(UNIQ(INSIDE))
      IN=WHERE(IMG(INSIDE) NE BLANK,NOK)
      IF NOK GT 0 THEN BEGIN
         AREAPERCENT=AREAPERCENT+N_ELEMENTS(INSIDE(IN))
         FLUXPERCENT=FLUXPERCENT+TOTAL(IMG(INSIDE(IN)))
         FROM1TO2D,INSIDE(IN),IMGSZ(1),XX,YY
      ENDIF
   ENDELSE
ENDFOR

; PETROSIAN CURVE
WINDOW,2,TITLE='PETROSIAN CURVE TO DETERMINE PETROSIAN RADIUS'
PINDEX=(FLUX/AREA)/LEVELS       ; construct petrosian index
; recall that petrosian index "eta" is simply a ratio of the average
; flux within an aperture to the surface brightness at the radius of
; aperture 
PLOT,LEVELS,PINDEX,XTITLE='Isophotal surface brightness',YTITLE='Isophotal petrosian index',$
     CHARSIZE=1.7
plots,[sb50,sb50],[min(pindex),max(pindex)],THICK=0.8
plots,[sb80,sb80],[min(pindex),max(pindex)],THICK=2.0

petro50=interpol(pindex,levels,sb50)
petro80=interpol(pindex,levels,sb80)
plots,sb50,petro50,psym=4,syms=2.5
plots,sb80,petro80,psym=4,syms=2.5

; Now do an interpolation to determine at what surface brightness the
; petrosian index hits a specified value
SBPINDEX=INTERPOL(LEVELS,PINDEX,PTARGET)
plots,[sbpindex,sbpindex],[min(pindex),max(pindex)],thick=2,linest=3

; On this plot, overplot the "theoretical" petrosian index curve.
; Start with the Beta model, numberically integrate to get the
; integrated energy within some radius, then divide by the surface
; brightness of that radius.  Need to do a best-fit to the beta model
; represented by this image, though (e.g., we know the Beta value but
; not the amplitude. 
;
; The amplitude (central surface brightness) can be computed as
; u(rc)/2^B, where rc=core radius, B=-3Beta+1/2.  So interpolate to
; determine what the surface brightness is at the core radius, RC:
;;SBRC=INTERPOL(FLUXPERCENT


; Now construct a contour at this level only, to get final properties
; of the associated contour
CONTOUR,IMG,LEVELS=SBPINDEX,PATH_INFO=PINFO,PATH_XY=PXY,/CLOSED,/PATH_DATA_COORDS

AREAPINDEX=0.
FLUXPINDEX=0.                   ; Determine the total area, etc
areatemp=0.
radpindex=0.
FOR J=0,N_ELEMENTS(PINFO.LEVEL)-1 DO BEGIN ; could be more than 1 contour at this level
   ISTART=PINFO(J).OFFSET
   XC=REFORM(PXY(0,ISTART:ISTART+PINFO(J).N-1))
   YC=REFORM(PXY(1,ISTART:ISTART+PINFO(J).N-1))
   IF N_ELEMENTS(XY) EQ 1 THEN BEGIN
      IF IMG(XC,YC) NE BLANK THEN BEGIN
         AREAPINDEX=AREAPINDEX+1.
         FLUXPINDEX=FLUXPINDEX+TOTAL(IMG(XC,YC))
      ENDIF
   ENDIF ELSE BEGIN
      INSIDE=POLYFILLV([XC,XC(0)],[YC,YC(0)],IMGSZ(1),IMGSZ(2))
      FROM2TO1D,XC,YC,IMGSZ(1),BORDER
      INSIDE=[INSIDE,BORDER]
      INSIDE=INSIDE(SORT(INSIDE))
      INSIDE=INSIDE(UNIQ(INSIDE))
      IN=WHERE(IMG(INSIDE) NE BLANK,NOK)
      IF NOK GT 0 THEN BEGIN
         AREAPINDEX=AREAPINDEX+N_ELEMENTS(INSIDE(IN))
         FLUXPINDEX=FLUXPINDEX+TOTAL(IMG(INSIDE(IN)))
         if nok gt areatemp then begin
            areatemp=nok
            radpindex=mean(sqrt( (xcentroid-xc)^2 + (ycentroid-yc)^2 ))
         endif
      ENDIF 
   ENDELSE
ENDFOR

; compute the theoretical radius and surface brightness at which the
; desired petrosian index should have been met:
theorytot=-1.0*!pi*sb0*(rcore^2)/(b+1.0)
theoryr50=rcore*sqrt( (1.0-0.5)^(1.0/(b+1.0)) - 1.0 )
theoryr80=rcore*sqrt( (1.0-0.8)^(1.0/(b+1.0)) - 1.0 )
theoryr50sb=sb0*(1.0+(theoryr50/rcore)^2)^b
theoryr80sb=sb0*(1.0+(theoryr80/rcore)^2)^b
; compare theorytot to totflux, theoryr50 to rad50, theory80 to rad80
; compare theoryr50sb to sb50, theoryr80sb to sb80
   

; Now that the contours have been made for the targeted petrosian
; index and percent-enclosed eneryy, we can now generate the final
; radii by circularizing the associated areas:

PRADIUS=SQRT(AREAPINDEX/!PI)
ERADIUS=SQRT(AREAPERCENT/!PI)

stop,'check the PRADIUS and ERADIUS, see if they make sense'

END
