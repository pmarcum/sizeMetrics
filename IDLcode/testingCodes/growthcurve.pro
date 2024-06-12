PRO GROWTHCURVE,BETA=BETA,I0=I0,E=E,RC=RC
; Predicts the growth curve (photometric) for a "perfect" de vauceleur
; or beta-model surface brightness profile, given an axial ratio of
; the ellipse of e (semi-minor/semi-major)

IF NOT KEYWORD_SET(BETA) THEN BETA=0.52
IF NOT KEYWORD_SET(E) THEN E=0.3
IF NOT KEYWORD_SET(RC) THEN RC=15.
IF NOT KEYWORD_SET(I0) THEN BEGIN
; Force I0 to be what it needs to be such that the profile equals 1 at
; the virial radius, which is approximately 3.085 times the core
; radius (see
; ned.ipac.caltech.edu/level5/Sept03/Peacock/Peacock6_5.html, eqns
; 117-118, to understand how I came to this conclusion

   CORE2VIRIAL=3.085
   RVIR=RC*CORE2VIRIAL
   I0=1.0/(1.0 + (RVIR/RC)^2)^(-3.0*BETA+0.5)
ENDIF

R=DINDGEN(1e9)
H=((1.0-E)/(1.0+E))^2
B=-3.0*BETA+0.5
CONST=!PI*(1.0+E)*(1.0+(3.0*H)/(10.0+SQRT(4.0-3.0*H)))/(2.0*(B+1.0))
FINT=I0*CONST*RC^2*( (1.0+(R/RC)^2)^(B+1)-1.0 )

;;;  NOTE: I am commenting out the below, because for a large number
;;   of points, it takes forever to execute.  And besides, it was only
;;   added to confirm, by a completely independent way, that the above
;;   equation is valid.   I have confirmed that when plotting R versus
;;   FINT and then R versus FSUM, the 2 plots lay exactly on top of each
;;   other!!   So the analytical expression above for FINT is the
;;   correct expression for the total flux integrated within an
;;   ellipse for a Beta model having the ellipticity, etc. as
;;   specified.  whew!!! 
;;INTENSITY=I0*(1.0+(R/RC)^2)^B
;;FSUM=R*0.
;;DAREA=CONST*(2.0*(B+1.0))*R  ; the approximatation for "perimeter" of an ellipse
;;FOR I=1,N_ELEMENTS(R)-1 DO FSUM(I)=INT_TABULATED(R(0:I),INTENSITY(0:I)*DAREA(0:I))

; note that int_tabulated works in the following way:  integral of
; function(x) dx, where dx is a simple array of radius, x-values,
; etc.  Even though in this case we are really doing function (x)
; d(area), where d(area) is also a function of x, to make
; int_tabulated work, we have to "trick" it into thinking that the
; x-dependent piece of the darea is really just part of the
; function(x).  So we group all the x-dependency together, call it
; f(x), and then let that consolidated function be integrated against
; x ("r" here). For example, if one wanted to integrate to get the
; area of a circle, you would not do the following:
; int_tabulated(2*!pi*r, r*0.+1), which is sort of the
; way you'd do it mathematically (e.g., integral{dArea} =
; integral{2 !pi r dr}, but rather you'd lump everything but
; the dr together and treat that as the function being integrated
; over, eg int_tabluated{r, 2*!pi*r} 

; Compute the half-light radius:

R50=RC*SQRT(0.5^(1.0/(B+1.0)) - 1.0 )

STOP,'1'

END
