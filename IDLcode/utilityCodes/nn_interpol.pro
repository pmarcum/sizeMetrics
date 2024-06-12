FUNCTION NN_INTERPOL,IMAGE, NOISE, SNR, CTSIMAGE=CTSIMAGE, $
   CENTER=CENTER, MAXCELLSIZE=MAXCELLSIZE
; Performs a "natural neighbor" interpolation by doing the following:
; (1) forms weighted voronoi tesselation on points X, Y to generate the
;     honeycombed-shaped "Vorioni" cells,
; (2) then forms a different set of voronoi cells by including point
;     X0,Y0 to the list of X Y values.
; (3) The vooroi cell around  X0,Y0 is overlayed on the first Veronoi;
;     cell map,
; (4) the percent overlap between the X0,Y0 cell and the surrounding
;     original cells is computed.  The percent areas of overlap
;     between the neightboring cells relative to the total area of the X0 Y0 cell
;     are computed 
; (5) the X0 Y0 pixel value is computed by the following weighted sum:
;     VALUE = f1 V1 + f2 V2 + f3 V3 + ... where fi is the fractional
;     area of overlap with the original cell i and Vi is the value of
;     that cell. 
; The value X0 Y0 should not be included in X Y
; The pixel value for X0 Y0 is returned.

; Put down the Voronoi cells for the underlying, known points using WVT:

BLANK=-9999.
  
IF NOT KEYWORD_SET(MAXCELLSIZE) THEN MAXAREA=36. ELSE MAXAREA=MAXCELLSIZE^2


WVT_IMAGE,IMAGE,NOISE,SNR,WVTRESULT,XNODE,YNODE,CTSIMAGE=CNTS, $
         CENTER=CENTER, MAX_AREA=MAXAREA,BINNUMBER=BINIMAGE

INONIMAGE=WHERE(IMAGE NE 0)
NORMVAL=MEAN(IMAGE(INONIMAGE))

; Now that we have the cells for the underlying image, do the
; interpolation for all the pixels that are not the XNODE, YNODE
; values:

RESULT=IMAGE*0.+BLANK

SZ=SIZE(IMAGE)

NX=SZ(1)
NY=SZ(2)

FOR IX=0,NX-1 DO BEGIN
   FOR IY=0,NY-1 DO BEGIN
      XARR=[IX,XNODE]
      YARR=[IY,YNODE]
      TRIANGULATE,XARR,YARR,TR,CONN=C
      VORONOI,XARR,YARR,0,C,XP,YP,[0,0,NX-1,NY-1]
      
      IBAD=WHERE(XP LT 0 OR YP LT 0 OR XP GT NX-1 OR YP GT NY-1,NBAD)
      IF NBAD GT 0 THEN GOTO,NEXTY
      
; Now figure out the overlap of this polygon (the X0, Y0 cell) with
; the cells laid down with just the X Y list of points:
      ICELL=POLYFILLV(XP,YP,NX,NY)
      TAREA=FLOAT(N_ELEMENTS(ICELL))
      VSUM=TOTAL(WVTRESULT(ICELL)/NORMVAL*(1.0/TAREA))
      RESULT(IX,IY)=VSUM*NORMVAL
      NEXTY:
   ENDFOR
ENDFOR

RETURN,RESULT

END
