PRO BETAPETRO
; generates the petrosian index for a range of values of Beta, R/Rc,
; and encircled energy radii.

RCMODEL=10.    ; arcseconds, the core radius used in model to generate image
PIXSCALE=0.492 ; number arcseconds on pixel side              
  
BETA=[0.52,0.55,0.60]
SELECT=2  ; index of the beta to display in plots.  0 is first beta value, 1 is second, etc. 
PERCENT=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
RCORE=RCMODEL/PIXSCALE

RRCARR=FINDGEN(10001)/0.01+600.
RRCARR=[FINDGEN(1501)/3.+0.01,RRCARR]
WTMP=RRCARR(1)-RRCARR(0)
RADII=[WTMP/2.+RRCARR(0)]
WIDTH=[WTMP]
FOR I=2,N_ELEMENTS(RRCARR)-1 DO BEGIN
   WTMP=RRCARR(I)-RRCARR(I-1)
   WIDTH=[WIDTH,WTMP]
   RADII=[RADII,WTMP/2.+RRCARR(I-1)]
ENDFOR
RRCARR=RADII & RADII=0

; Using beta model, determine the surface brightness at each of the
; above normalized radii:
SB0=1.0                         ; central surface brightness
B=-3.0*BETA+0.5
SB=DINDGEN(N_ELEMENTS(BETA),N_ELEMENTS(RRCARR))
FOR I=0,N_ELEMENTS(BETA)-1 DO SB(I,*)=SB0*(1.0+RRCARR^2)^B(I)

; Now construct growth curves by computing the integrated flux at each radius
GROWTH=SB*0.
GROWTH(*,0)=SB(*,0)*!PI*RRCARR(*,0)^2

FOR I=0,N_ELEMENTS(BETA)-1 DO $
   FOR J=1,N_ELEMENTS(RRCARR)-1 DO $
      GROWTH(I,J)=TOTAL(SB(I,0:J)*2.*!PI*RRCARR(0:J)*WIDTH(0:J)*RCORE^2)
   
;;;      GROWTH(I,J)=INT_TABULATED(RRCARR(0:J)*RCORE,SB(I,0:J)*2.*!PI*RRCARR(0:J)*RCORE)

; Now construct petrosian indice curve
PINDEX=SB*0.
FOR I=0,N_ELEMENTS(BETA)-1 DO $
   PINDEX(I,*)=(GROWTH(I,*)/(!PI*(RRCARR(*)*RCORE)^2))/SB(I,*)

ETOT=BETA*0.
FOR I=0,N_ELEMENTS(BETA)-1 DO ETOT(I)=MAX(GROWTH(I,*))
ETHEORY=-1.*!PI*SB0*RCORE^2/(B+1)

;find percent encircled energy radii
RVAL=FLTARR(N_ELEMENTS(BETA),N_ELEMENTS(PERCENT))
FOR I=0,N_ELEMENTS(BETA)-1 DO $
   FOR J=0,N_ELEMENTS(PERCENT)-1 DO RVAL(I,J)=INTERPOL(RRCARR(*)*RCORE,GROWTH(I,*),ETHEORY(I)*PERCENT(J))

PVAL=FLTARR(N_ELEMENTS(BETA),N_ELEMENTS(PERCENT))  ; find petrosian index for different percent enclosed light radii
FOR I=0,N_ELEMENTS(BETA)-1 DO $
   FOR J=0,N_ELEMENTS(PERCENT)-1 DO PVAL(I,J)=INTERPOL(PINDEX(I,*),RRCARR(*)*RCORE,RVAL(I,J))

WINDOW,1,TITLE='Beta function'
PLOT,RRCARR,SB(SELECT,*),XRA=[0,500],/XST,XTIT='R/Rc',YTIT='Beta model (r)',CHARSIZE=1.8
OPLOT,[RVAL(SELECT,4),RVAL(SELECT,4)]/RCORE,[MIN(SB(SELECT,*)),MAX(SB(SELECT,*))]
OPLOT,[RVAL(SELECT,7),RVAL(SELECT,7)]/RCORE,[MIN(SB(SELECT,*)),MAX(SB(SELECT,*))]

WINDOW,2,TITLE='Growth curve'
PLOT,RRCARR,GROWTH(SELECT,*),XRA=[0,500],/XST,XTIT='R/Rc',YTIT='Growth curve (r)',CHARSIZE=1.8
OPLOT,[RVAL(SELECT,4),RVAL(SELECT,4)]/RCORE,[MIN(GROWTH(SELECT,*)),MAX(GROWTH(SELECT,*))]
OPLOT,[RVAL(SELECT,7),RVAL(SELECT,7)]/RCORE,[MIN(GROWTH(SELECT,*)),MAX(GROWTH(SELECT,*))]

WINDOW,3,TITLE='Petrosian index curve'
PLOT,RRCARR,PINDEX(SELECT,*),XRA=[0,500],/XST,XTIT='R/Rc',$
     YTIT='Petrosian index = (Avg SB within R)/SB(R)',CHARSIZE=1.8
OPLOT,[RVAL(SELECT,4),RVAL(SELECT,4)]/RCORE,[MIN(PINDEX(SELECT,*)),MAX(PINDEX(SELECT,*))]
OPLOT,[RVAL(SELECT,7),RVAL(SELECT,7)]/RCORE,[MIN(PINDEX(SELECT,*)),MAX(PINDEX(SELECT,*))]
OPLOT,[MIN(RRCARR),MAX(RRCARR)],[PVAL(SELECT,4),PVAL(SELECT,4)]
OPLOT,[MIN(RRCARR),MAX(RRCARR)],[PVAL(SELECT,7),PVAL(SELECT,7)]

stop,'inspect the plots, investigate the relationship bet/ p50 and p80 with beta to confirm expectations'

END
