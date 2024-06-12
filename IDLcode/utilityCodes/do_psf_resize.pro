PRO DO_PSF_RESIZE,XPOLY,YPOLY,EDGE,XTRIM,YTRIM
; Code taken from the test that figured out how to snip an edge off of
; a polygon of a given size EDGE.  This code accepts a contour
; defined by XBLUR, YBLUR and snips away EDGE margin, returns XTRIM,YTRIM

NPOLY=N_ELEMENTS(XPOLY)

NEWX=XPOLY*0. & NEWY=NEWX
; 
; Now construct an image in which pixels that are enclosed by the
; polygon are masked by the value of "1"
XDIM=MAX(XPOLY)+10.*EDGE
YDIM=MAX(YPOLY)+10*EDGE
IMAGE=INTARR(XDIM,YDIM)
INDEX=POLYFILLV([XPOLY,XPOLY(0)],[YPOLY,YPOLY(0)],XDIM,YDIM)
IMAGE(INDEX)=1
WRITEFITS,'test.fits',IMAGE

FOR I=0,NPOLY-1 DO BEGIN
; Now go through each point, fit a smooth curve locally using the
; neighbor to either side, and then use that curve to determine the
; local tangent line, from which the perpendicular can be computed.
; grab the indices of the 2 nearest neighbor points to either side of
; point I.  Have to consider the special case when I is close to the
; beginning or end, to insure that the nearest neighbors are drawn
; from the other end of the loop. 
   IF I EQ 0 THEN IND=[NPOLY-2,NPOLY-1,0,1,2] ELSE $
      IF I EQ 1 THEN IND=[NPOLY-1,0,1,2,3] ELSE $
         IF I EQ NPOLY-2 THEN IND=[I-2,I-1,I,I+1,0] ELSE $
            IF I EQ NPOLY-1 THEN IND=[I-2,I-1,I,0,1] ELSE $
               IND=[I-2,I-1,I,I+1,I+2]

; temporarily redefine the coordinate system for these nearest
; neighbor points, in which Point I is the origin
   XARR=XPOLY(IND)-XPOLY(I)
   YARR=YPOLY(IND)-YPOLY(I)

; Determine the relative angle of the position vectors, wrt the
; central point, of the 2 neighboring points
   ANG1=ATAN(XARR(1),YARR(1))
   ANG3=ATAN(XARR(3),YARR(3))

   ANGDIFF=ABS(ANG1-ANG3)
; If the angle difference is greater than 180 degrees, then the larger
; angle was chosen instead of the smaller one between the 2 vectors,
; so subtract from 360 degrees
   IF ANGDIFF GT !PI THEN ANGDIFF=ABS((2.*!PI)-ANGDIFF)

   
   MINANGLE=!PI/4.
   ANGTOLERANCE=0.6
   BISECTOR=-9999
   
   IF ABS(ANG1-ANG3)LT MINANGLE THEN BEGIN
; The angle of separation is small, such that they form a skinny
; triangle (with the central point at the pointy edge), fitting a spline
; is probably not the best approach, as the radius of curvature of a spline
; may not be sharp enough to turn a sharp corner of the triangle.  In
; this  case, the best bet is taking the bisector ande declaring the
; bisector to be the perpendicular at that location.   If the angle is
; larger, but the point orientation does not fall under the above
; categories  (which guarantee an angular separation of at least 90
; degrees), then proceed with the spline fit.  Let's say that
; if the angular separation is less than 45 degrees, go with the
; bisector, and if greater than 45 degrees, proceed with a spline
; fit. Will have to work through all the possible cases, though. 
      BISECTOR=MEAN([ANG1,ANG3])
; If the BISECTOR angle is greater than 90 degrees or smaller than -90
; degrees, the obtuse angle was chosen, and the acute angle needs to
; be selected.  So just add (or subtract) from 180 degrees
      IF BISECTOR GT !PI/2. THEN BISECTOR=!PI-BISECTOR
      IF BISECTOR LT -!PI/2. THEN BISECTOR=!PI+BISECTOR
      PERP_SLOPE=TAN(BISECTOR)
   ENDIF ELSE BEGIN
; In this case, the angle separating the 2 points is larger than 45
; degrees, and large enough for a spline to be made ... so proceed
; Now do spline fit and interpolation to points that lay
; just a wee bit to either side of the point of interest
      SPLINE_P,XARR,YARR,XSPLINE,YSPLINE ; parametric cubic spline interpolation      
; locate the 2 points to either side of the central point
      RADSPLINE=SQRT( XSPLINE^2+YSPLINE^2 ) 
      ANGSPLINE=ATAN(YSPLINE,XSPLINE) 

      IGOOD=WHERE(RADSPLINE NE 0.0) ; picking out all points except the origin
      INEAR1=WHERE(RADSPLINE(IGOOD) EQ MIN(RADSPLINE(IGOOD))) & INEAR1=IGOOD(INEAR1(0))
; Remove the neighbor just picked, so it doesn't get picked again
      IGOOD(INEAR1)=-999 & ITMP=WHERE(IGOOD NE -999) & IGOOD=IGOOD(ITMP)  

; The other nearest neighbor needs to be at a substantially different
; angle, otherwise you are not picking up the nearest neighbor on the
; other side but rather another close neighbor on the same side. We
; know that the other side is likely to be much greater than MINANGLE,
; but to account for possible curvature in the fit, we'll allow
; the angle to be down to ANGSCALE percent of MINANGLE

; Now determine the difference in position angle of all the other
; points relative to this first near neighbor that has been
; identified:
      DIFFANG=ABS(ANGSPLINE-ANGSPLINE(INEAR1))
      IFIX=WHERE(DIFFANG GT !PI,NFIX)
      IF NFIX GT 0 THEN DIFFANG(IFIX)=ABS(2.0*!PI-DIFFANG(IFIX))
      
      IBIGANG=WHERE(DIFFANG(IGOOD) GE MINANGLE*ANGTOLERANCE) 
      INEAR2=WHERE(RADSPLINE(IGOOD(IBIGANG)) EQ MIN(RADSPLINE(IGOOD(IBIGANG))))
      INEAR2=IGOOD(IBIGANG(INEAR2(0)))

      XTANG=[XSPLINE(INEAR1),0.0,XSPLINE(INEAR2)]
      YTANG=[YSPLINE(INEAR1),0.0,YSPLINE(INEAR2)]
; Now fit a line to the 2 interpolated points that are super-close and
; bracketing the central point of interest.  The line will essentially
; be the tangent to the shape at this location
      FIT=POLY_FIT(XTANG,YTANG,1)   
      TAN_SLOPE=FIT(1)
  
; The above is the equation of the line for the tangent at this point,
; but we need the equation of the line for the orthogonal direction at
; this point, the line that is perpendicular to this line.  THis problem
; is, yes, another high school problem!
;          y = const + (-1./slope)*x
; where slope is the tangent slope. the polygon point also lay on this
; perpendicular, so we can solve for the new constant of this line:
;           const = y - (-1/slope)*x
;                 = 0 + 1/slope*0
;                 = 0
; which again makes sense given our careful choice of coord system.
; Thus, for the PERPENDICULAR LINE:
;           y = -1/slope*x
;           y = -1./[FIT(1)+2*FIT(2)*x+3*FIT(3)*X^2+4*FIT(4)*X^3] * x

      IF TAN_SLOPE NE 0.0 THEN $
         PERP_SLOPE=-1.0/TAN_SLOPE $
            ELSE STOP,'the slope is exactly zero -- put in more code to work around this situation!'
         
   ENDELSE
   
; Now what we do is slide the XPOLY(i), YPOLY(i) point down along this
; perpendicular to a distance that corresponds to the trim margin, and
; call that new location the new corresponding point in the shrunk
; polygon.  In other words, the old coordinate moves to this new coordinate
; as its contribution to creating the new, smaller, concentric polygon
;
; To accomplish this task, we essentially need to determine the point of
; intersection between a circle whose radius is the amount to be
; trimmed ("margin") and this perpendicular line.
;      y=perp_slope*x  (the equation of the perpendicular line)
;      (x^2 + y^2) = (margin to be trimmed)^2
; (the latter is a circle of radius MARGIN centered on the polygon point).
; Note that all the x and y's in the above are referenced to the polygon
; point as the origin of the coordinate system. 
; Solving these 2 equations, we would normally have gotten a quadratic,
; but our choice of origin simplified the situation.  The solution is :
;
;      X = M/sqrt(1+(-1/slope)^2)  ; where slope is slope of tangent line
;      Y = -1/slope * X
;  Note that X could be either positive or negatve and still solve
; the equations (a line will have 2 points of intersection in the circle, if
; it intersects it at all).
;
; So which will be the correct choice?  The one that is closest to
; the original reference point (e.g., the one closest to the inside
; of the polygon). So at this point, will need to shift the origin
; back to the reference point:
;      X = M/sqrt(1+perp_slope^2) + XPOLY(I)
;      Y = perp_slope * X + YPOLY(i)
;  Do the above to both of the possible coordinates.  The point that
;  falls in the area previously defined as the interoior (e.g., is a
;  masked pixel) will be the one to choose.
;
; In complicated shapes, there could be masked pixels on either side
; of the perpendicular (e.g., some parts of the polygon could have
; significant curvature and put a "penisula" out there for an
; otherwise "external" perpendcular endpoint to land in).  In such a
; case, there would be some unmasked points that lay between that
; "Fake" endpoint and the point of interest.  So let the real interior
; perpendicular end be the one that lay on the half of the line that
; has the highest percentage of masked points.  Make a series of
; points that go from the point of interest in either direction along
; the perpendicular to either candidate endpoint:
   
   
; Want to define an array of points that extend from the central point
; on the polygon in both directions down the perpendicular, to the
; EDGE endpoint on either side.

; Define an array of fractions of EDGE:
   FEDGE=FINDGEN(11)/10.        ; 0,0.1,0.2,0.3,....1.0
   XPERP1=(FEDGE*EDGE)/SQRT(1.0+PERP_SLOPE^2)
   YPERP1=PERP_SLOPE*XPERP1

   XPERP2=-1.0*(FEDGE*EDGE)/SQRT(1.0+PERP_SLOPE^2)
   YPERP2=PERP_SLOPE*XPERP2

   XPERP1=XPERP1+XPOLY(I) & YPERP1=YPERP1+YPOLY(I)
   XPERP2=XPERP2+XPOLY(I) & YPERP2=YPERP2+YPOLY(I)
   
   
   XCAND1=EDGE/SQRT(1.0+PERP_SLOPE^2)
   YCAND1=PERP_SLOPE*XCAND1

   XCAND2=-1.0*EDGE/SQRT(1.0+PERP_SLOPE^2)
   YCAND2=PERP_SLOPE*XCAND2

   XCAND1=XCAND1+XPOLY(I)
   YCAND1=YCAND1+YPOLY(I)
   XCAND2=XCAND2+XPOLY(I)
   YCAND2=YCAND2+YPOLY(I)
  

   ITMP=WHERE(IMAGE(ROUND(XPERP1),ROUND(YPERP1)) EQ 1,NMASKED1)
   ITMP=WHERE(IMAGE(ROUND(XPERP1),ROUND(YPERP1)) EQ 0,NNOT1)
   ITMP=WHERE(IMAGE(ROUND(XPERP2),ROUND(YPERP2)) EQ 1,NMASKED2)
   ITMP=WHERE(IMAGE(ROUND(XPERP2),ROUND(YPERP2)) EQ 0,NNOT2)

   NPERC1=FLOAT(NMASKED1)/FLOAT(NNOT1+NMASKED1)
   NPERC2=FLOAT(NMASKED2)/FLOAT(NNOT2+NMASKED2)

; OK, which of these 2 sets of candidate points (they are points at a
; distance of EDGE away from the polygon point, to either side, laying
; along the perpendicular line to the curve of the shape at that point)
; is closest to the inside of the polygon?

; If any of the new coordinates fall off the edge of the image,
; that's a good sign that the new coordinate was not the right
   NEWX(I)=-999. & NEWY(I)=-999.
; one to pick
   IF XCAND1 LT 0. OR XCAND1 GT XDIM-1 THEN BEGIN
      NEWX(I)=XCAND2 & NEWY(I)=YCAND2
   ENDIF ELSE IF XCAND2 LT 0. OR XCAND2 GT XDIM-1 THEN BEGIN
      NEWX(I)=XCAND1 & NEWY(I)=YCAND1
   ENDIF ELSE IF YCAND1 LT 0. OR YCAND1 GT YDIM-1 THEN BEGIN
      NEWX(I)=XCAND2 & NEWY(I)=YCAND2
   ENDIF ELSE IF YCAND2 LT 0. OR YCAND2 GT YDIM-1 THEN BEGIN
      NEWX(I)=XCAND1 & NEWY(I)=YCAND1
   ENDIF
   IF NEWX(I) NE -999. AND NEWY(I) NE -999. THEN BEGIN
; Do a final check to make sure that the coordinates that were picked
; did not also fall off edge
      IF NEWX(I) LT 0. OR NEWX(I) GT XDIM-1 THEN STOP,'something bad happened with new X'
      IF NEWY(I) LT 0. OR NEWY(I) GT YDIM-1 THEN STOP,'something bad happened with new Y'
   ENDIF ELSE BEGIN
; Probably, both sets of candidates fell wthin the original frame of
; the image.  So then pick who landed on the side with the most masked
; pixels:
      IF NPERC1 GT NPERC2 THEN BEGIN
         NEWX(I)=XCAND1 & NEWY(I)=YCAND1
      ENDIF ELSE BEGIN
         NEWX(I)=XCAND2 & NEWY(I)=YCAND2
      ENDELSE
; Now do some tests.  The chosen coordinates should land on masked
; pixels
      IF IMAGE(ROUND(NEWX(I)),ROUND(NEWY(I))) NE 1 THEN STOP,'Weird shit going down ....'
   ENDELSE
ENDFOR


INDEX=POLYFILLV([NEWX,NEWX(0)],[NEWY,NEWY(0)],XDIM,YDIM)
IMAGE(INDEX)=1+IMAGE(INDEX)
WRITEFITS,'test.fits',IMAGE

; Need to prune the polygon, if there are places in which the new
; polygon border twisted because the points on one side moved past the
; location of the points on the other side of a structure during the
; new coordinate transform. Theoretically this situation should never
; happen if the structure is real and if the psf has been adequately
; characterized, but there could be artifacts that cause such a
; situation.

; To look for places in which the polygon twists/folds over itself,
; need to compare the line segment formed by 2 adjacent points while
; going along the boundary with the line segement formed by all the
; other possible adjacent pairs and determine if the line segments
; intersect.

; How to determine if line segments intersect?  One approach is to
; compute the coordinate at which the intersection would occur if the
; two line segments were instead infinitely long lines, and see if
; that intersection position lay within the bounding box of either
; line segment (doesn't really matter which bounding box is
; selected).  There are other approaches (mathematical in nature), but
; this idea seems as good as any.

; Fit every single adjacent pair of points to a line and determine
; y-intercept and slope.
SLOPE=NEWX*0.
YINT=SLOPE
XLL=SLOPE                       ; lower left corner of bounding box in X
XUR=SLOPE                       ; upper right corner of bounding box in X
YLL=SLOPE
YUR=SLOPE

FOR I=0,NPOLY-2 DO BEGIN
   IF NEWX(I+1)-NEWX(I) NE 0. THEN $
      SLOPE(I)=( NEWY(I+1)-NEWY(I) )/( NEWX(I+1) - NEWX(I) )  ELSE $
         SLOPE(I)=-9999            ; flag as a value of infinity -- a vertical line
   IF SLOPE(I) NE -9999 THEN $
      YINT(I)=NEWY(I)-SLOPE(I)*NEWX(I) ELSE $
         YINT(I)=-9999.         ; vertical line has no y intercept
   XLL(I)=MIN([NEWX(I),NEWX(I+1)])
   XUR(I)=MAX([NEWX(I),NEWX(I+1)])
   YLL(I)=MIN([NEWY(I),NEWY(I+1)])
   YUR(I)=MAX([NEWY(I),NEWY(I+1)])
ENDFOR
IF NEWX(NPOLY-1)-NEWX(0) NE 0. THEN $
   SLOPE(NPOLY-1)=( NEWY(0)-NEWY(NPOLY-1) )/( NEWX(0) - NEWX(NPOLY-1) ) ELSE $ 
      SLOPE(NPOLY-1)=-9999
IF SLOPE(NPOLY-1) NE -9999 THEN $
   YINT(NPOLY-1)=NEWY(NPOLY-1)-SLOPE(NPOLY-1)*NEWX(NPOLY-1) ELSE $
      YINT(NPOLY-1)=-9999.
XLL(NPOLY-1)=MIN([NEWX(NPOLY-1),NEWX(0)])
XUR(NPOLY-1)=MAX([NEWX(NPOLY-1),NEWX(0)])
YLL(NPOLY-1)=MIN([NEWY(NPOLY-1),NEWY(0)])
YUR(NPOLY-1)=MAX([NEWY(NPOLY-1),NEWY(0)])

; OK, we now know the line-fits to every adjacent pair. Now take one
; pair, compare to all others and determine if there are any twisted
; pieces of the new polygon.

XTRIM=NEWX
YTRIM=NEWY
DONE=XTRIM*0.

FOR I=0,NPOLY-1 DO BEGIN
; There is no need to do comparisons to segment pairs that are far away and
; have no prayer of intersecting. Weed those out by only selecting the
; pairs whose bounding boxes at least have some overlap with the
; bounding box of this pair:

   IF DONE(I) EQ -999 THEN GOTO,NEXT
; don't waste time on a segment that has already been
; identified as being part of a cross-over piece

   IND=INDGEN(NPOLY)

; indice all the other pairs except the one of interest here, and the
; pair that it makes with the neighbor to either side, as by
; definition those line segments intersect (at the endpoints)

; first find any segment pair for which a corner of the bounding box
; falls within the bounding box of another pair.  If the only common
; overlap in bounding boxes is right at the ends of the segment, then
; disregard as a crossed pair (touching bounding boxes are not
; overlapping bounding boxes)
   ILLC=WHERE(XLL(I) GT XLL(IND) AND XLL(I) LT XUR(IND) AND $
              YLL(I) GT YLL(IND) AND YLL(I) LT YUR(IND) AND $
              SLOPE(I) NE SLOPE(IND) AND DONE(IND) EQ 0,NLLC)
   IURC=WHERE(XUR(I) GT XLL(IND) AND XUR(I) LT XUR(IND) AND $
              YUR(I) GT YLL(IND) AND YUR(I) LT YUR(IND) AND $
              SLOPE(I) NE SLOPE(IND) AND DONE(IND) EQ 0,NURC)
   ILRC=WHERE(XLL(I)+(XUR(I)-XLL(I)) GT XLL(IND) AND XLL(I)+(XUR(I)-XLL(I)) LT XUR(IND) AND $
              YLL(I) GT YLL(IND) AND YLL(I) LT YUR(IND) AND $
              SLOPE(I) NE SLOPE(IND) AND DONE(IND) EQ 0,NLRC)
   IULC=WHERE(XLL(I) GT XLL(IND) AND XLL(I) LT XUR(IND) AND $
              YLL(I)+(YUR(I)-YLL(I)) GT YLL(IND) AND YLL(I)+(YUR(I)-YLL(I)) LT YUR(IND) AND $
              SLOPE(I) NE SLOPE(IND) AND DONE(IND) EQ 0,NULC)

   IF NLLC GT 0 THEN ILLC=IND(ILLC)
   IF NURC GT 0 THEN IURC=IND(IURC)
   IF NLRC GT 0 THEN ILRC=IND(ILRC)
   IF NULC GT 0 THEN IULC=IND(IULC)
   
; Now determine if there are any bounding boxes that forms a cross
; with another bounding box, such that the segments associated with
; those bounding boxes indeed cross, but the corners of either box lay
; outside the other's box  (sort of like the situation below):
;
;        ____
;       |    2
;       |    |
;       |    |
;       |    |
;   3...|....|.....
;   .   |    |    .
;   .   |    |    .
;   ....|....|....4
;       |    |
;       |    |
;       |    |
;       1____
;
; pair 1,2 cross with pair 3,4 but the corners of their bounding boxes
; are outside the other's bounding box


   ICROSS=WHERE(XLL(I) GT XLL(IND) AND XLL(I) LT XUR(IND) AND $
                XUR(I) GT XLL(IND) AND XUR(I) LT XUR(IND) AND $
                SLOPE(I) NE SLOPE(IND) AND DONE(IND) EQ 0,NCROSS)
; If such a situation exists, then the other pair should have corners
; that lay within the range of y values of the bounding box for point
; I if the 2 segments truly overlap
   IF NCROSS GT 0 THEN BEGIN
      II=WHERE(YLL(IND(ICROSS)) GT YLL(I) AND YLL(IND(ICROSS)) LT YUR(I) AND $
               YUR(IND(ICROSS)) GT YLL(I) AND YUR(IND(ICROSS)) LT YUR(I),NCROSS)
      IF NCROSS GT 0 THEN ICROSS=ICROSS(II) ELSE ICROSS=[-1]
   ENDIF

   IOVERLAP=[ILLC,IURC,ILRC,IULC,ICROSS]
   II=WHERE(IOVERLAP NE -1,NOVERLAP)
   IF NOVERLAP NE 0 THEN IOVERLAP=IOVERLAP(II)
   
   FOR J=0,NOVERLAP-1 DO BEGIN
      IF SLOPE(I) NE -9999 AND SLOPE(IOVERLAP(J)) NE -9999 THEN BEGIN
;If neither line is vertical,determine the coordinate at which the
;segments would cross if they were extended to be infinitely-long lines
         XCROSS=(YINT(I)-YINT(IOVERLAP(J)))/(SLOPE(IOVERLAP(J))-SLOPE(I))
         YCROSS=YINT(I)+SLOPE(I)*XCROSS
      ENDIF ELSE IF SLOPE(I) EQ -9999 AND SLOPE(IOVERLAP(J)) NE -9999 THEN BEGIN
; Vertical line intersecting another line that is not vertical
         XCROSS=NEWX(I)
         YCROSS=YINT(IOVERLAP(J))+SLOPE(IOVERLAP(J))*XCROSS
      ENDIF ELSE IF SLOPE(IOVERLAP(J)) EQ -9999 AND SLOPE(I) NE -9999 THEN BEGIN
; Vertical line intersecting another line that is not vertical
         XCROSS=NEWX(IOVERLAP(J))
         YCROSS=YINT(I)+SLOPE(I)*XCROSS
      ENDIF
      
      IF XCROSS GT XLL(I) AND XCROSS LT XUR(I) AND $
         YCROSS GT YLL(I) AND YCROSS LT YUR(I) AND $
         XCROSS GT XLL(IOVERLAP(J)) AND XCROSS LT XUR(IOVERLAP(J)) AND $
         YCROSS GT YLL(IOVERLAP(J)) AND YCROSS LT YUR(IOVERLAP(J)) THEN BEGIN

;;;;  The intersection must fall within BOTH bounding boxes to truly
;;;;  be a segment interection (and not just the intersection of the
;;;;  extension of the segments into infinite lines)

; need to trim out all points that lay between the 2 indices here
; 
; When a pair of segments cross, the loop could be seen as either the
; true loop that needs to be removed, but the loop could also be the
; main polygon.  The only way to remove the chaff from the wheat is to
; make the assumption that the loop will be the smaller of the 2
; loops. To determine loop size, go along the loop and sum up the
; distances between the points.  The loop with the largest
; circumference will be considered to be the "loop" that should
; remain.  The other loop is to be removed.  In a future version of
; this code, could have a keyword that declares where the centroid of
; the polygon is (e.g., location of galaxy), and force program to
; choose the loop that encircles that point.


;********************************************************************
;                               L O O P   1
;
;     
;
;
;                    >>>>> LOOP 1 <<<<<
;       ---->       |     *
;        * * *      |   *  * /\
;      *       *    \/ *   *  |
;     *           *    *  *   |
;   *               *  * *     ORDER OF POLYGON IS IN COUNTERCLOCKWISE
;/\ *  loop2           *       DIRECTION
;|  *                    * 
;|  *                    *   |
;    *                   *   |    
;     *                  *  \/
;       *               *
;            *  *  *  *
;            <------
;
;*********************************************************************
      
         ILOOP1=NEWX*0.
         IF I LT IOVERLAP(J) THEN BEGIN
; Compute distance from the point of intersection to the first inner
; loop point, point I+1
            LOOP1=SQRT( (XCROSS-NEWX(I+1))^2 + (YCROSS-NEWY(I+1))^2 )
; Now add up the rest of the distances between points within the loop as one goes
; around the loop. Treat the pairs as the current coordinate of
; interest and the coordinate that is next in line along the polygon
; (rather than the current coord of interest and the point behind
; it).
            FOR K=I+1,IOVERLAP(J)-1 DO $
               LOOP1=LOOP1+SQRT( (NEWX(K)-NEWX(K+1))^2 + (NEWY(K)-NEWY(K+1))^2 )
; As last step, we handle the joining of the point IOVERLAP(J)
; and the intersection.
            LOOP1=LOOP1+SQRT( (XCROSS-NEWX(IOVERLAP(J)))^2 + (YCROSS-NEWY(IOVERLAP(J)))^2 )
; Indicate all the points that are considered inside-loop points that
; should potentially be removed
            ILOOP1(I+1:IOVERLAP(J))=1
         ENDIF ELSE IF I GT IOVERLAP(J) THEN BEGIN
; go around polygon from IOVERLAP(J)+1 (assuming first inner point of
; loop to I (last loop point on other side)
            LOOP1=SQRT( (XCROSS-NEWX(IOVERLAP(J)+1))^2 + (YCROSS-NEWY(IOVERLAP(J)+1))^2 )
            FOR K=IOVERLAP(J)+1,I-1 DO $
               LOOP1=LOOP1+SQRT( (NEWX(K)-NEWX(K+1))^2 + (NEWY(K)-NEWY(K+1))^2 )
            LOOP1=LOOP1+SQRT( (XCROSS-NEWX(I))^2 + (YCROSS-NEWY(I))^2 )
            ILOOP1(IOVERLAP(J)+1:I)=1
         ENDIF ELSE STOP,'SOMETHING SCREWED UP AT LOOP1'

;********************************************************************
;                               L O O P   2
;
;     
;
;
;                       loop2
;       ---->       |     *
;        * * *      |   *  * /\
;      *       *    \/ *   *  |
;     *           *    *  *   |
;   *               *  * *     ORDER OF POLYGON IS IN COUNTERCLOCKWISE
;/\ *  >>> LOOP 2<<<<   *       DIRECTION
;|  *                    * 
;|  *                    *   |
;    *                   *   |    
;     *                  *  \/
;       *               *
;            *  *  *  *
;            <------
;
;*********************************************************************
      
; Now go around whatever part of the polygon that didn't get
; gone around before, and assume that was in fact the actual loop to
; be removed ("loop2" as designated in diagram above)

         ILOOP2=NEWX*0.
         IF I LT IOVERLAP(J) THEN BEGIN
; link up the intersection point to the POINT I (which we are
; considering to be an inside point to the loop)
            LOOP2=SQRT( (XCROSS-NEWX(I))^2 + (YCROSS-NEWY(I))^2 )
; Now go from 0 to point I, add up distances along the loop
            FOR K=0,I-1 DO $
               LOOP2=LOOP2+SQRT( (NEWX(K)-NEWX(K+1))^2 + (NEWY(K)-NEWY(K+1))^2 )
; Go from the IOVERLAP(J)+1 point, the first loop point on the other
; side, and go around back to the 0 point
            FOR K=IOVERLAP(J)+1,NPOLY-2 DO $
               LOOP2=LOOP2+SQRT( (NEWX(K)-NEWX(K+1))^2 + (NEWY(K)-NEWY(K+1))^2 )
; Now link up the end points of the polygon
            LOOP2=LOOP2+SQRT( (NEWX(0)-NEWX(NPOLY-1))^2 + (NEWY(0)-NEWY(NPOLY-1))^2 )
; Now link up the intersection point to the other loop end
; (IOVERLAP(J)+1). Assuming that IOVERLAP(J) is not at the end of the
; polygon
            IF IOVERLAP(J) NE NPOLY-1 THEN BEGIN
               LOOP2=LOOP2+SQRT( (XCROSS-NEWX(IOVERLAP(J)+1))^2 + (YCROSS-NEWY(IOVERLAP(J)+1))^2 )
               ILOOP2(IOVERLAP(J)+1:NPOLY-1)=1
            ENDIF ELSE $
               LOOP2=LOOP2+SQRT( (XCROSS-NEWX(0))^2 + (YCROSS-NEWY(0))^2 )
            ILOOP2(0:I)=1
         ENDIF ELSE IF I GT IOVERLAP(J) THEN BEGIN
            IF I NE NPOLY-1 THEN BEGIN
               LOOP2=SQRT( (XCROSS-NEWX(I+1))^2 + (YCROSS-NEWY(I+1))^2 )
               ILOOP2(I+1:NPOLY-1)=1
            ENDIF ELSE $
               ILOOP2=SQRT( (XCROSS-NEWX(0))^2 + (YCROSS-NEWY(0))^2 )
            FOR K=I,NPOLY-2 DO $
               LOOP2=LOOP2+SQRT( (NEWX(K)-NEWX(K+1))^2 + (NEWY(K)-NEWY(K+1))^2 )
            LOOP2=LOOP2+SQRT( (NEWX(NPOLY-1)-NEWX(0))^2 + (NEWY(NPOLY-1)-NEWY(0))^2 )
            FOR K=0,IOVERLAP(J)-1 DO $
               LOOP2=LOOP2+SQRT( (NEWX(K)-NEWX(K+1))^2 + (NEWY(K)-NEWY(K+1))^2 )
            LOOP2=LOOP2+SQRT( (XCROSS-NEWX(IOVERLAP(J)))^2 + (YCROSS-NEWY(IOVERLAP(J)))^2 )
            ILOOP2(0:IOVERLAP(J))=1
         ENDIF ELSE STOP,'SOMETHING SCREWED UP AT LOOP2'

; All we now have to ask is, which loop circumference is largest?  Then
; assume that the smaller of the two is the one that should be
; removed:

         IF LOOP1 LT LOOP2 THEN ITRIM=ILOOP1 ELSE ITRIM=ILOOP2
         IREPLACE=WHERE(ITRIM NE 0)
         XTRIM(IREPLACE)=-999
         YTRIM(IREPLACE)=-999
         DONE(IREPLACE)=1
         XTRIM(IREPLACE(0))=XCROSS
         YTRIM(IREPLACE(0))=YCROSS
      ENDIF
   ENDFOR
NEXT:
ENDFOR
   

IKEEP=WHERE(XTRIM NE -999)

INDEX=POLYFILLV([XTRIM(IKEEP),XTRIM(0)],[YTRIM(IKEEP),YTRIM(0)],XDIM,YDIM)
IMAGE(INDEX)=1+IMAGE(INDEX)
WRITEFITS,'test.fits',IMAGE


END
  
