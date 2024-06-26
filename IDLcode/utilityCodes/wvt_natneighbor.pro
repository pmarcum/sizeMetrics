PRO WVT_NATNEIGHBOR,INFILE,CNTSFILE,EXPFILE,OUTFILE,SNR=SNR,$
   CENTER=CENTER,MAXCELLSIZE=MAXCELLSIZE

; Applies the Weighted Voronoi Tesselations (Diehl & Statler 2006,
; MNRAS 368, 497) to the fits image called INFILE that was produced by
; MARX x-ray simulation (or could be a real Chandra ACIS image file),
; then produces a "continuous" image that can be contoured by
; stretching a surface across the points corresponding to the centers
; of each tesselation cell.  We could just leave the image as-is after
; applying the WVT, but the problem is that each cell, which can be
; several pixels, is all the same value and has a honeycomb shape,
; which seriously impacts the ability to draw realistic contours using
; the resulting image.  Effectively, what we do here is to remove
; those artifical plateaus and force there to be gradients in the
; pixels in betwen the centers of the cells to create a more realistic
; light distribution that is more continous than the original x-ray
; image. This interpolation is performed by "natural neighbors"
; interpolation, which is done in the following manner: 
;
; (1) forms weighted voronoi tesselation on image to generate the
;     honeycombed-shaped "Vorioni" cells,
; (2) determines the XVor and YVor coordinates of cell centers, 
; (3) systematically goes through every pixel in the image, and forms
;     a new set of Voronoi cells by including the pixel Xi, Yi to the
;     list XVor and YVor.
; (4) The Voronoi cell around  Xi,Yi is overlayed on the
;     originally-generated Veronoi;
; (5) the percent overlap between the Xi,Yi cell and the surrounding
;     original cells is computed.  The percent areas of overlap
;     between the neightboring cells relative to the total area of the Xi Yi cell
;     are computed 
; (5) the Xi Yi pixel value is computed by the following weighted sum:
;     VALUE = f1 V1 + f2 V2 + f3 V3 + ... where fi is the fractional
;     area of overlap with the original cell i and Vi is the value of
;     that cell.
; (6) Keep going through the image one pixel at a time, Xi Yi until
;     all pixels have had values computed for them. 


BLANK=-9999.             ; Value of pixel where interpolation was not possible

FXREAD,CNTSFILE,CNTS            ; counts per pixel
FXREAD,EXPFILE,EFFEXP        ; effective area and effective exposure time map (cm^2 s photon/count)

NOISE=SQRT(CNTS)/EFFEXP
; The noise is computed by sqrt(counts)/expmap.  Note that the flux
; image, in units of photons/s/cm^2/pixel is generated by dividing the
; counts image (soft_thresh.img, output from the CIAO command
; "fluximage") by the effective exposure map (soft_thres.expmap),
; which also accounts for the effective area if fluximage is run with
; the default for the type of exposure map that is generated.  See
; Diehl & Statler 2006, MNRAS 368, 497 eqn 9.

FXREAD,INFILE,PHOTFLUX,HDR      ; photons/s/cm^2/pixel

; target signal-to-noise ratio per tesselation cell, can be changed to
; whatever:
IF NOT KEYWORD_SET(SNR) THEN SNR=4.0
; A target SNR of 4-10 is standard for x-ray, see pg 8 of WVT
; manual. A temperature map would require SNR~30-100, and integral
; field spectroscopy, SNR~50 per bin would be reasonable to extract
; stellar kinematics info from galaxy spectra.

IF NOT KEYWORD_SET(MAXCELLSIZE) THEN MAXCELLSIZE=6.0
; Maximum size of a cell, as measured on a side if cell were
; square. This value is used to compute the maximum square area of a
; cell. 

MAXAREA=MAXCELLSIZE^2

WVT_IMAGE,PHOTFLUX,NOISE,SNR,WVTRESULT,XNODE,YNODE,CTSIMAGE=CNTS, $
         CENTER=CENTER, MAX_AREA=MAXAREA,BINNUMBER=BINIMAGE

INONIMAGE=WHERE(WVTRESULT NE 0)
NORMVAL=MEAN(WVTRESULT(INONIMAGE))

; Now that we have the cells for the underlying image, do the
; interpolation for all the pixels that are not the XNODE, YNODE
; values:

RESULT=PHOTFLUX*0.+BLANK

SZ=SIZE(PHOTFLUX)

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

WRITEFITS,OUTFILE,RESULT,HDR

END
