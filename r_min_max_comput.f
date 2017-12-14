c	version 1: april, 1997

c       MODULE r_min_max_comput
c	including functions and subroutines
c	RMM
c	RSTART
C	RTBISEC
c	SDPR2

c  THE CALLING PROGRAM MUST HAVE DEFINED:
c	cutoff
c	the "common /s_d_param/" variables

C &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

	SUBROUTINE RMM (
c
C  INPUT
     *		cutoff, 
C  OUTPUT
     *          rmin, rmax)

c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c input: cutoff (scalar): threshold for the integrand function 
c			  (size dis. times pi times square radius)
c output: rmin, rmax (scalars): integration limits for Mie integrals
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

	implicit none

	real cutoff, rmin, rmax

c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c local variables
c maxrmax (scalar): max (in a physically sense) value for rmax (in micron)
c minrmin (scalar): min (in a physically sense) value for rmin (in micron)
c epsmx (scalar): relative accuracy for rmax
c epsmn (scalar): relative accuracy for rmin
c rm(2) (vector), r: temporary names for rmax and rmin
c rst (logical): to test the position of the initial value for rmax and rmin 
c		 with respect to maxrmax and minrmin
c eps(2) (vector): local names for epsmx and epsmn 
c i: index for loop to compute rmax (i=1) and rmin (i=2)
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

	real maxrmax, minrmin, epsmx, epsmn
	real mrm(2), rm(2), r, eps(2)
	integer i
	logical rst

	parameter (maxrmax=100.,  minrmin=5.e-5,  
     *             epsmx=1.e-4,  epsmn=1.e-4)
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c functions
c rtbisec: search for the root of the equation [ sdpr2 - cutoff = 0] with
c          the accuracy eps(i)
c rstart: search for the initial point to compute rmax and rmin
c sdpr2: compute the product size-distribution times pi times square radius
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

	real rtbisec, rstart, sdpr2

c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c		 EXECUTABLE STATEMENTS
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

	mrm(1) = maxrmax
	mrm(2) = minrmin
	eps(1) = epsmx
	eps(2) = epsmn

	do 100 i=1,2
		r = rstart(i)
		rst = (r.ge.mrm(i)) .eqv. (i.eq.1)
		if (sdpr2(mrm(i)) .ge. cutoff .or. rst) then
			rm(i) = mrm (i)
			goto 100
		else
			if (sdpr2(r) .le. cutoff) goto 50
			r = rtbisec(sdpr2,r,mrm(i),eps(i),
     *				    3-2*i,cutoff)
		endif
			
 50		rm(i) = r

 100 	continue

	rmax = rm(1)
	rmin = rm(2)

	return
	end


c &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

	FUNCTION RSTART (
c
C  INPUT
     *		i)

c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c input: i (scalar): index to select between rmax (i=1) and rmin (i=2)
c		     case
c output: rstart (scalar): initial radius in searching rmin and rmax
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

	implicit none

	real rstart
	integer i

c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c common variables
c block: s_d_param
c neq (scalar): size distr "type" (neq=1 sumof lognorm.'s, neq=2 gamma-mod.)
c modanum (scalar): num ber of lognorm.'s modalities (max 10)
c lognorpar (matrix): lognorm.'s modal radii, sigma and mixing ratios 
c gammodpar (vector): the four gamma-modified parameters   
c &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

	integer neq, modanum
	real lognorpar(3,10)
	real gammodpar(4)
	common /s_d_param/ lognorpar,gammodpar,neq,modanum

c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c local variables
c jmod (scalar): index for loop over modalities
c modr (scalar): modality's modal radius
c sigma (scalar): modality's standard deviation
c alpha, b, gamma (scalars): gamma-modified parameters   
c sngchk (scalar): plus (or minus) one, to drive the max (or min) search
c rr (scalar): guess for rstart
c &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

	integer jmod
	real sngchk, rr, modr, sigma, alpha, b, gamma

c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c		 EXECUTABLE STATEMENTS
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

	if (neq.eq.1) then
		if (i.eq.1) then
			sngchk=1.
			rstart=-1.e6
		else
			sngchk=-1.
			rstart=1.e6
		endif

		do jmod=1,modanum
		   modr=lognorpar(1,jmod)
		   sigma=lognorpar(2,jmod)
		   rr=modr*10.**(sigma*sigma*alog(10.))
		   if(sngchk*rr.gt.sngchk*rstart)  rstart=rr 
		enddo
	else
		alpha=gammodpar(2)
		b=gammodpar(3)
		gamma=gammodpar(4)
		rstart=( (alpha+2)/(b*gamma))**(1./gamma)
	endif

	return

	end




c &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

	FUNCTION RTBISEC(
c
C  INPUT
     *		func,x1,x2,xacc,isel,thr)

c modified version of Num. Rec's "rtbis" routine.
c by M. Cervino, April 1997
c a root betw. x1 and x2 isearched, for the equation func-thr=0
c Main new features:
c	** The function "func" is evaluated from a threshold "thr"
c	** The accuracy "xacc" is relative to the final abscissa "rtbisec"
c	** There is a NEW selection flag "isel"
c			isel=0 "rtbisec" is returned as in "rtbis" (the 
c					accuracy is at least "xacc")
c			isel<0 "rtbisec" is returned lessened by dx (the
c                                       accuracy is at least 2 * xacc)
c                       isel>0 "rtbisec" is returned increased by dx (the
c                                       accuracy is at least 2 * xacc)   
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c input: 
c func (function): the one-variable function tested
c x1, x2 (scalars): interval limits for the root
c xacc (scalar): RELATIVE accuracy for the root
c thr (scalar): a constant parameter of the equation func-thr=0
c isel (scalar): a flag to select the root (see above)
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

	implicit none

      REAL rtbisec,x1,x2,xacc,func,thr
      INTEGER isel

c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c local variables
c jmax (scalar): max number os allowed bisections
c j (scalar): index for the loop over bisections
c dx (scalar): interval width (halved each iteration)
c f, fmid (scalars): intermediate values of func-thr
c xmid (scalar): intermediate value of the searched root
c &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

      INTEGER JMAX, J
      PARAMETER (JMAX=60)
      REAL dx,f,fmid,xmid

c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c functions
c func: dummy name of the tested function
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

      EXTERNAL func

c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c		 EXECUTABLE STATEMENTS
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

      if(isel.gt.0) isel = 1 
      if(isel.lt.0) isel = -1 
      fmid=func(x2)-thr
      f=func(x1)-thr
      if(f*fmid.ge.0.) pause 'root must be bracketed in rtbisec'
      if(f.lt.0.)then
        rtbisec=x1
        dx=x2-x1
      else
        rtbisec=x2
        dx=x1-x2
      endif
      do 11 j=1,JMAX
        dx=dx*.5
        xmid=rtbisec+dx
        fmid=func(xmid)-thr
        if(fmid.le.0.)rtbisec=xmid
	if(fmid.eq.0.) return
	if(rtbisec.eq.0.) rtbisec=1.e-18 
        if(abs(dx)/abs(rtbisec).lt.xacc) go to 20
11    continue
      pause 'too many (more than 60!) bisections in rtbisec'
 20   rtbisec=rtbisec+float(isel)*abs(dx)
      return
      END

c &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

	FUNCTION sdpr2(
c
C  INPUT
     *		r)

c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c input: r (scalar): radius at which the size-dis. times pi times square 
c radius (sdpr2) is returned.
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

	implicit none

	real r, sdpr2

c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c common variables
c block: s_d_param
c neq (scalar): size distr "type" (neq=1 sumof lognorm.'s, neq=2 gamma-mod.)
c modanum (scalar): num ber of lognorm.'s modalities (max 10)
c lognorpar (matrix): lognorm.'s modal radii, sigma and mixing ratios 
c gammodpar (vector): the four gamma-modified parameters   
c &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

	integer neq, modanum
	real lognorpar(3,10)
	real gammodpar(4)
	common /s_d_param/ lognorpar,gammodpar,neq,modanum

c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c local variables
c jmod (scalar): index for loop over modalities
c modr (scalar): modality's modal radius
c sigma (scalar): modality's standard deviation
c mixrat (scalar): modality's mixing ratio
c a, alpha, b, gamma (scalars): gamma-modified parameters   
c pi, root2p, ln10 (scalars): obvious constants (see below)
c expon (scalar): exponent for the lognormal sizedis.
c &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

	integer jmod
	real modr, sigma, mixrat, a, alpha, b, gamma
	real pi, root2p, ln10, expon

c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c		 EXECUTABLE STATEMENTS
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

	pi	= acos(-1.e0)
	root2p	= sqrt(pi+pi)
	ln10	= alog(10.)

	if (neq.eq.1) then

		sdpr2=0.
		do jmod=1,modanum
		   modr=lognorpar(1,jmod)
		   sigma=lognorpar(2,jmod)
		   mixrat=lognorpar(3,jmod)
		   expon=-.5*(alog10(r/modr))**2/(sigma*sigma)
		   sdpr2=sdpr2+ mixrat/sigma*exp(expon)
		enddo
		sdpr2=sdpr2/(root2p*r*ln10)
	else
		a=gammodpar(1)
		alpha=gammodpar(2)
		b=gammodpar(3)
		gamma=gammodpar(4)
		sdpr2=a*r**alpha*exp(-b*r**gamma)
	endif

	sdpr2=sdpr2*pi*r*r

	return
	end

