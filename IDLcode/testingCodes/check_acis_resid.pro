PRO CHECK_ACIS_RESID
; Flips through a list of residual images produced by taking the
; deconvolution of a CHANDRA ACIS image, then convolving it back with
; the PSF, then subtracting from the original image to make sure that
; the residual is mostly zero, which would indicate that the
; deconvolution went pretty well.  This program plops up each image on
; the screen for visual inspection and produces some basic stats just
; to make sure there are no obscenely negative or positive values
; within the image.

DIR='/Gdrive/Pam/ACIS_Mar2017/'
RDIR=DIR+'Residuals/'
IMGDIR=DIR+'Images/'
DDIR=DIR+'Deconvolved/'

READCOL,DIR+'deconv.list',BASENAME,FORMAT='A'
BASENAME=STRTRIM(BASENAME,2)

FOR I=0,N_ELEMENTS(BASENAME)-1 DO BEGIN
   rname=rdir+basename(i)+'resid.fits'
   oname=imgdir+basename(I)+'regimg3.fits'
   dname=ddir+basename(i)+'deconv.fits'
   FXREAD,rname,IMG
   fxread,oname,orig
   fxread,dname,deconv

   in=where(orig ne 0)
   mind=min(orig(in))*1.2 & maxd=max(orig(in))*0.8
   display,oname,0,/adj,/noc,/noz,fld=1,/rev,mind=mind,maxd=maxd,title=basename(i)
;   IMG=ABS(IMG)
;   FXREAD,IMGDIR+BASENAME(I)+'regimg3.fits',ORIG
;   PERCENT=IMG*0.0-9999.
;   IN=WHERE(ORIG NE 0.)
;   PERCENT(IN)=IMG(IN)/ORIG(IN)*100.
;   WRITEFITS,'tmp.fits',PERCENT
;   DISPLAY,'tmp.fits',0,/ADJ,/NOC,/NOZ,FLD=1,$
;           TITLE=BASENAME(I)+' (RESIDUAL)',MIND=MIND,MAXD=MAXD1
;   MEANPIX=MEAN(IMG)
;   MEDIANPIX=MEDIAN(IMG)
;   STPIX=STDDEV(IMG)
;   MINMAXPIX=MINMAX(IMG)
;stop,'1'
;   FXREAD,IMGDIR+BASENAME(I)+'regimg3.fits',ORIG
;   IN=WHERE(ORIG NE 0)
;   ORIG=ORIG/MEDIAN(ORIG(IN))
;   WRITEFITS,'tmp.fits',ORIG
;   DISPLAY,'tmp.fits',0,/ADJ,/NOC,/NOZ,FLD=2,$
;           TITLE=BASENAME(I)+' (ORIG)',MIND=MIND,MAXD=MAXD2

;   FXREAD,DDIR+BASENAME(I)+'deconv.fits',DECON
;   IN=WHERE(DECON NE 0)
;   DECON=DECON/MEDIAN(DECON(IN))
;   WRITEFITS,'tmp.fits',DECON
;   DISPLAY,'tmp.fits',0,/ADJ,/NOC,/NOZ,FLD=3,$
;           TITLE=BASENAME(I)+' (DECON)',MIND=MIND3,MAXD=MAXD3
;   PRINT,MEANPIX,STPIX,MEDIANPIX,MINMAXPIX(0),MINMAXPIX(1),$
;         FORMAT='("MEAN: ",E9.2,"+/-",E9.2,2X,"MEDIAN: "'+$
;                  ',E9.2,2X,"MIN: ",E9.2,1X,"MAX: ",E9.2)'
   A='' & READ,'hit any key to continue',A

ENDFOR

END
