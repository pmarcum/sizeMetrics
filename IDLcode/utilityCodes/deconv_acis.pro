PRO DECONV_ACIS
; Constructs deconvoluted images from Chanda/ACIS using psfs obtained
; through the Chandra database. 

; Read the list of images to be processed:

NITER=50  ; USER-defined max number of iterations in the maximum likelihood analysis

IMGSUFF='regimg3'
PSFSUFF='psf3'
OPSFSUFF='psf'
DIR='/Gdrive/Pam/ACIS_Mar2017/'
PSFDIR=DIR+'PSF_soft_CSC/'
OPSFDIR=DIR+'oldPSF/'
IMGDIR=DIR+'Images/'
DDIR=DIR+'Deconvolved/'
RDIR=DIR+'Residuals/'
ODIR=DIR+'Oldpsf_deconv/'
CDIR=DIR+'ConvergenceTests/'

READCOL,DIR+'deconv.list',BASENAME,FORMAT='A'
BASENAME=STRTRIM(BASENAME,2)

; Now go through each image, one by one, apply the PSF in the
; deconvolution and produce a result, then place result into google
; drive directory

FOR I=0,N_ELEMENTS(BASENAME)-1 DO BEGIN
; Read in the image
   FXREAD,IMGDIR+BASENAME(I)+IMGSUFF+'.fits',IMG,IMGHDR
; Just to get rid of the small numbers and hopefully allow
; computations to work more smoothly, divide the image by the mean
; nonzero pixel value
   NORMVAL=MEAN(IMG(WHERE(IMG GT 0.)))
   IMG=IMG/NORMVAL
; Read in the PSF
   FXREAD,PSFDIR+BASENAME(I)+PSFSUFF+'.fits',PSF
; Insure that theintegral under PSF function is equal unity (e.g., normalize
; function)
   PSF=PSF/TOTAL(PSF)
; Perform deconvolution
; Reset the deconvolved image to a single pixel, so that it will start
; out fresh
  DECONV=FLTARR(1)
  RECONV=FLTARR(1)
  DELTA1=FLTARR(NITER)
  DELTA2=FLTARR(NITER)
  FOR J=0,NITER-1 DO BEGIN
     MAX_LIKELIHOOD,IMG,PSF,DECONV,RECONV
     IN=WHERE(DECONV NE 0)
     IF J EQ 1 THEN PREV=DECONV
     IF J GT 1 THEN BEGIN
        DELTA1(J)=MEAN((DECONV(IN)-PREV(IN))/DECONV(IN))*100.
        PREV=DECONV
     ENDIF
     IN=WHERE(IMG NE 0)
     DELTA2(J)=MEAN( (IMG(IN)-RECONV(IN))/IMG(IN))*100.
  ENDFOR

  TH=FINDGEN(41)/40.*2.*!PI
  XP=COS(TH) & YP=SIN(TH)
  USERSYM,XP,YP,/FILL
  WINDOW,1
  YMIN=MIN(DELTA1)-ABS(MIN(DELTA1))*0.05
  YMAX=MAX(DELTA1)+ABS(MAX(DELTA1))*0.05
  PLOT,FINDGEN(NITER),DELTA1,SYMSIZE=1.0,PSYM=8,$
       TITLE='CONVERGENCE TEST: DIFFERENCE BETWEEN ith and (i-1)th DECONVOLVED IMAGE',$
       XTIT='ITERATION NUMBER',YTIT='MEAN PERCENT DIFFERENCE',$
       XRA=[0,50],/XST,CHARSIZE=1.7,YRA=[YMIN,YMAX],/YST
  IM=TVRD()
  WRITE_JPEG,CDIR+BASENAME(I)+'convergence1.jpg',IM
  WINDOW,2
  YMIN=MIN(DELTA2)-ABS(MIN(DELTA2))*0.05
  YMAX=MAX(DELTA2)+ABS(MAX(DELTA2))*0.05
  PLOT,FINDGEN(NITER),DELTA2,SYMSIZE=1.0,PSYM=8,$
       TITLE='CONVERGENCE TEST: DIFFERENCE BETWEEN ORIGINAL and RECONVOLVED IMAGE',$
       XTIT='ITERATION NUMBER',YTIT='MEAN PERCENT DIFFERENCE',XRA=[0,50],/XST,CHARSIZE=1.7,$
       YRA=[YMIN,YMAX],/YST

  IM=TVRD()
  WRITE_JPEG,CDIR+BASENAME(I)+'convergence2.jpg',IM
  
  WRITEFITS,DDIR+BASENAME(I)+'deconv.fits',DECONV*NORMVAL,IMGHDR
  PRINT,'Just created '+BASENAME(I)+'deconv.fits'

; See how well the deconvoluted image, after convolving back with the
; PSF, matches the original image.
  RESID=IMG-RECONV
  IZERO=WHERE(IMG EQ 0)
  IN=WHERE(IMG NE 0)
  RESID(IN)=(RESID(IN)/IMG(IN))*100.
  RESID(IZERO)=-9999.
  WRITEFITS,RDIR+BASENAME(I)+'resid.fits',RESID,IMGHDR

; For grins, read in the old PSF and see how well a job it does
; compared to the smoother PSF
;  FXREAD,OPSFDIR+BASENAME(I)+OPSFSUFF+'.fits',OPSF
;  OPSF=OPSF/TOTAL(OPSF)
;  FOR J=1,NITER DO MAX_LIKELIHOOD,IMG,OPSF,DECONV
;  WRITEFITS,ODIR+BASENAME(I)+'deconv.fits',DECONV*NORMVAL,IMGHDR
  
ENDFOR

END
