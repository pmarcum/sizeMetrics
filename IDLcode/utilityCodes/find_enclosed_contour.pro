PRO FIND_ENCLOSED_CONTOUR,IMAGE,XCONT,YCONT,CVALUE,FRAC=FRAC,pref=pref,suff=suff
; given an image, will determine the contour that produces the P
; percent enclosed energy (where total energy is assumed for these
; purposes to be the total summed pixel values of the image).

; Generate a bunch of contours, so that the contour holding 50% of all
; the light can be found
;
; Divide the levels up evenly between the 30%  of the mean value and
; up to 80% of the highest pixel value 

dir='/home/pmarcum/marx-master/marxsims/'
if keyword_set(pref) and keyword_set(suff) then begin
   file=dir+pref+'effects'+suff+'.fits'
   fxread,file,image
endif

IN=WHERE(IMAGE NE 0)
MEANPIX=MEAN(IMAGE(IN))
IMAGE=IMAGE/MEANPIX

IF NOT KEYWORD_SET(FRAC) THEN FRAC=0.5 ; 50% light enclosed

MAXVAL=MAX(IMAGE)*0.9
MINVAL=MIN(IMAGE)*1.1

LEV=(FINDGEN(101)/100.)*(MAXVAL-MINVAL)+MINVAL

SZ=SIZE(IMAGE)
X0=FIX(SZ(1)/2.)
Y0=FIX(SZ(2)/2.)

;;NEW=MIN_CURVE_SURF(IMAGE)
;;IMAGE=NEW & NEW=0.0

window,1
contour,image,findgen(sz(1)),findgen(sz(2)),levels=lev,/closed

window,2
in=where(image ne 0)
from1to2d,in,sz(1),xx,yy
contour,(image(in)),xx,yy,levels=lev,/closed,/irregular

stop,'1111'

CONTOUR,IMAGE,FINDGEN(SZ(1)),FINDGEN(SZ(2)),LEVELS=LEV,PATH_INFO=INFO,$
         PATH_XY=XY,/CLOSED,/PATH_DATA_COORDS

; Go through each contour level ...
FENCLOSED=LEV*0.
RADIUS=LEV*0.

WINDOW,1,XSIZE=800,YSIZE=800
PLOT,FINDGEN(SZ(1)),FINDGEN(SZ(2)),/ISO,XRA=[X0-100,X0+100],YRA=[Y0-100,Y0+100],/YST,/NODATA

FOR I=0,N_ELEMENTS(LEV)-1 DO BEGIN
; Grab all contours associated with this level
   ILEV=WHERE(INFO.LEVEL EQ I,NILEV)
   AREA=FLTARR(NILEV)
; Find the contour that encloses the central point
   FOR J=0,NILEV-1 DO BEGIN
      ISTART=INFO(ILEV(J)).OFFSET
      IEND=INFO(ILEV(J)).OFFSET+INFO(ILEV(J)).N-1     
      XC=REFORM(XY(0,ISTART:IEND))
      YC=REFORM(XY(1,ISTART:IEND))
      IF N_ELEMENTS(XC) EQ 1 THEN BEGIN
         RADTEST=SQRT( (XC-X0)^2 + (YC-Y0)^2 )
         IF RADTEST LE 3 THEN AREA(J)=1.
      ENDIF ELSE BEGIN
         INSIDE=POLYFILLV([XC,XC(0)],[YC,YC(0)],SZ(1),SZ(2))
         FROM1TO2D,INSIDE,SZ(1),XX,YY
         IN=WHERE(XX EQ X0 AND YY EQ Y0,NTEST)
         IF NTEST NE 0 THEN AREA(J)=N_ELEMENTS(INSIDE)
      ENDELSE
   ENDFOR

   if i gt 80 then stop,'1'
   
   IN=WHERE(AREA EQ MAX(AREA))
   IN=IN(0)
   IF AREA(IN) NE 0 THEN BEGIN
      IBIG=ILEV(IN)
; Get the starting and ending point for this single contour
      ISTART=INFO(IBIG).OFFSET
      IEND=INFO(IBIG).OFFSET+INFO(IBIG).N-1
; Extract the x/y coordinates for this single contour
      XC=REFORM(XY(0,ISTART:IEND))
      YC=REFORM(XY(1,ISTART:IEND))
; Get all the pixels within this contour and add up the total light
      IF N_ELEMENTS(XC) EQ 1 THEN BEGIN
        FENCLOSED(I)=IMAGE(XC,YC)
        PLOTS,XC,YC
      ENDIF ELSE BEGIN
         INSIDE=POLYFILLV([XC,XC(0)],[YC,YC(0)],SZ(1),SZ(2))
         FENCLOSED(I)=TOTAL(IMAGE(INSIDE))
         PLOTS,XC,YC
      ENDELSE
;;   PLOTS,X0,Y0,PSYM=4,SYMSI=2
      RADIUS(I)=MEAN(SQRT( (XC-X0)^2 + (YC-Y0)^2 ))
   ENDIF
ENDFOR

; Now "plot up" the growth curve and figure out the contour that most
; closely matches 50% encircled energy
;
; Since we have such finely-gridded levels, let's just take
; whatever is closed to 50% rather than do some kind of fancy
; interpolation
IPICK=WHERE( ABS(FENCLOSED-TOTAL(IMAGE)*FRAC) EQ MIN(ABS(FENCLOSED-TOTAL(IMAGE)*FRAC)))
IPICK=IPICK(0)

; Now get the x,y coordinates for this contour
ILEV=WHERE(INFO.LEVEL EQ IPICK)
IBIG=WHERE(INFO(ILEV).N EQ MAX(INFO(ILEV).N))
IBIG=ILEV(IBIG(0))
ISTART=INFO(IBIG).OFFSET
IF IBIG EQ N_ELEMENTS(LEV)-1 THEN IEND=ISTART ELSE IEND=INFO(IBIG+1).OFFSET
XCONT=REFORM(XY(0,ISTART:IEND))
YCONT=REFORM(XY(1,ISTART:IEND))

; While we are at it, save the contour level value that fit the
; enclosed fraction criteria, as well
CVALUE=INFO(IBIG).VALUE

STOP,'1'
WINDOW,1
PLOT,RADIUS,LEV

STOP,'1'
END
