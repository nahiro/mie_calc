      SUBROUTINE legendre_coeff(NWL,NTHETA,THET,WAVE,X,ELLE,ASYM)
C
c-----------------------------------------------------------------
c ngl       :  number of points for the Gauss-Legendre quadrature
c n_thetph  :  number of points in the input phase function
c	       (GAEROS_DATASET.F)
c elle      :  degree of the Legendre expansion
c-----------------------------------------------------------------
	integer ngl,flag_ctrl,elle,indc,nwl,ntheta
        parameter(ngl=1000,max_el=300,NTHET_MAX=1000,
     &            N_WL_OUT=1000)
        real    x(0:elle,nwl),ASYM(NWL),wld(100),
     &          thet(ntheta),weig(ngl),phasep_gl(ngl),ab(ngl),
     &          gamma(nthet_MAX),phase(nthet_MAX),wave(nwl),plg,
     &          control,diff,diff_mean,max_diff,integral,effe,
     &          lambda,pi,a
        common /FDIFASE/ phas_fun(1000,N_WL_OUT)
c       write(*,*)nwl,ntheta,elle
C    test con la funzione isotropa
c      do 200 iwl=1,nwl
c       do 100 i=1,ntheta
c        write(33,*)i,iwl,phas_fun(i,iwl)
c 100  continue
c 200  continue
C
        pi=acos(-1.0)
        do jj=1,nwl
            do jjj=0,elle
             x(jjj,jj)=0.0
            enddo
        enddo
	flag_ctrl=1
	pi=acos(-1.0d0)

	open(unit=77,file='legendre_coeff.dat',
     &	     form='formatted')
C
	if(flag_ctrl.ne.0)
C
     &    OPEN(UNIT=22,FILE='sens.LOG')

c----------------------------------------------------------------------
c  		input phase function interpolation
c----------------------------------------------------------------------

	do ind=1,ntheta
		gamma(ind)= thet(ind) * pi /180.0d0
	enddo
C
	call gauleg(0.0,pi,ab,weig,ngl)

c linear interpolation of the phase function in the ngl points for quadrature
C==============================
C begin loop over wavelengths
C==============================
        nintegr=0
        do 330 iwl=1,nwl
C=============================
	ind_start=2
        ph_int=0.0

	  do ind=1,ngl
 	    do jjj=ind_start,NTHETA
		if(gamma(jjj).ge.ab(ind)) then	
			ind_start=jjj
			goto 3333
		endif
       	    enddo
 3333	   continue


	   call inter2(1,ab(ind),gamma(jjj-1),gamma(jjj),
     &        phasep_gl(ind),phas_fun(jjj-1,iwl),
     &        phas_fun(jjj,iwl))

           ph_int=ph_int+phasep_gl(ind)*sin(ab(ind))*weig(ind)
                       
          enddo
           ph_integral=ph_int*2.0*pi
           err_ph=abs(ph_integral-(4.0*pi))/(4.0*pi)*100.0
           write(*,100)err_ph
           if (abs(ph_integral-(4.0*pi)) .gt. 7.0e-3) then
             nintegr=nintegr+1
           endif

          write(*,111)wave(iwl), ph_int*2.0*pi
 111    format(2x, 'wl=',f5.3,2x,'Int_ph=',e12.5)
C
	if(flag_ctrl.ne.0) THEN
           WRITE(22,*)  ' NUMBER OF LEG. MOMENTS       = ', elle
C           WRITE(22,*)  ' WAVELENGTH                  = ', lambda
           WRITE(22,*)  '  WAVELENGTH                  = ', wave(iwl)
	endif
c compute the elle Legendre coefficients and delta-m renormalize them
	do i=0,elle
		do iii=1,ngl
	           x(i,iwl)= x(i,iwl)+ phasep_gl(iii)*
     &              plg(i,0,cos(ab(iii)))*sin(ab(iii))*weig(iii)
 
    		enddo
	enddo
C========== modifica fatta per Phil ==============
         do i=0,elle
C            x(i,iwl)= x(i,iwl) * 4.0 * PI / 2.0
             x(i,iwl)= x(i,iwl) / 2.0
         enddo
C
C        do i=0,elle
C           x(i,iwl)=(2.0*i+1.0)/2.0 * x(i,iwl)
C        enddo
C=================================================

	effe = x(elle,iwl)
c	effe = 0.0

	if(flag_ctrl.ne.0)
     &     write(22,*)  ' EFFE                        = ', effe

	 do i=0,elle
		write(77,*) i, x(i,iwl)	
C		x(i,iwl) = (x(i,iwl)-effe)/(1.0 - effe)
	 enddo
C===================================
C  end loop over wavelengths
C===================================
  330    CONTINUE
C===================================
         icont=0
         do iil=1,nwl
          firstlc=(x(0,iil))-1.0
          diffcoeff=( x(1,iil)-(3.*asym(iil)*x(0,iil)))*100.
          if (firstlc .gt. 1.0e-3 .and. diffcoeff .gt. 1.0e-2) then
            icont=icont+1
            wld(icont)=wave(iil)
          endif
         enddo
         if (icont.ne.0) then 
          write(*,*)'***** WARNING on LEGENDRE COEFFICIENTS *****'
	  write(*,*)'* Problem at wavelength=',(wld(iL),iL=1,icont)
          write(*,*)'* Solution: increase the number of angles in the 
     & phase function tabulation'
         endif
          if (nintegr.ne.0) then 
          write(*,*)'****** WARNING on THE PHASE FUNCTION ****'
          WRITE(*,*)'* The integral NOT equal 4 pi) '
          WRITE(*,*)'* Solution: increse the number of angles in the
     & phase function tabulation'
          endif

	close(77)

	if(flag_ctrl.ne.0) then
	  call ctrl(nwl,ngl,NTHETA,elle,x,gamma,ab,weig)
	endif
C
	if(flag_ctrl.ne.0)
     &		close(22)

	return
 100    FORMAT('the integral over P(gamma) differs from 4PI for',
     & E12.5,'%')
	end

c--------------------------------------------------------------------

      FUNCTION PLG(L,M,X)

      IMPLICIT NONE
      INTEGER  I,L,LL,M
      REAL   FACT,PLL,PMM,PMMP1,SOMX2,X,PLG

C -----------  0<M<L  -1<X<1  -------------------------
      IF(M.LT.0.OR.M.GT.L.OR.ABS(X).GT.1.0) THEN
        WRITE(6,*) 'Bad arguments'
        STOP
      ENDIF
C ------------ Compute P(M,M)  ------------------------      
      PMM = 1.0
      IF(M.GT.0) THEN
        SOMX2 = SQRT((1.00-X)*(1.00+X))
        FACT = 1.0
        DO 11 I=1,M
           PMM = -PMM*FACT*SOMX2
           FACT = FACT+2.0
  11    CONTINUE
      ENDIF
      IF(L.EQ.M) THEN 
        PLG = PMM
      ELSE
C ------------ Compute P(M,M+1)  ---------------------
        PMMP1 = X * FLOAT(2*M+1) * PMM
        IF(L.EQ.M+1) THEN
           PLG = PMMP1
        ELSE
C ------------Compute P(L,M)  -------------------------
           DO 12 LL = M+2, L
              PLL = (X * FLOAT(2*LL-1) * PMMP1 - FLOAT(LL+M-1) * PMM)
     &              / FLOAT(LL-M)
              PMM   = PMMP1
              PMMP1 = PLL
  12       CONTINUE
           PLG = PLL
        ENDIF     
      ENDIF
      RETURN
      END  

c-------------------------------------------------------------------------

        SUBROUTINE INTER2(INTYPE,X,X1,X2,F,F1,F2)

C    "INTER2" INTERPOLATES TO DETERMINE THE VALUE OF F AT X,
C     GIVEN F1 AT X1 AND F2 AT X2.
C       INTYPE=1 FOR LINEAR INTERPOLATION
C       INTYPE=2 FOR LOGARITHMIC INTERPOLATION
C  ............................................................
        
        REAL X,X1,X2,F,F1,F2,a1,a2,a
         
        ITYPE=INTYPE
        IF(F1 .LE. 0.0  .OR. F2 .LE. 0.0)ITYPE=1
        IF(ITYPE .EQ. 2)GO TO 100
C
        F=F1+(X-X1)*(F2-F1)/(X2-X1)
        RETURN
C
 100    A1=LOG(F1)
        A2=LOG(F2)
        A=A1+(X-X1)*(A2-A1)/(X2-X1)
        F=EXP(A)
        RETURN
        END

c-------------------------------------------------------------------------


	subroutine ctrl(nwl,ngl,N_THETPH,elle,x,gamma,
     &                  absci,weights)
        PARAMETER ( N_WL_OUT=1000 )
	integer elle,max_number,nden
        common /FDIFASE/phase_in(1000,N_WL_OUT)
	real x(0:elle,nwl),phase_out(1000,N_WL_OUT),gamma(N_THETPH),
     &       absci(ngl),weights(ngl),diff,max_diff,diff_mean,plg,
     &       max_angle,integral,phase,pi
        pi=acos(-1.0)

        do ileng=1,nwl
C
	diff_mean=0.0
        diff=0.0
        max_diff=0.0
        max_angle=0.0
	integral = 0.0
	nden = 0
	do ind=1,N_THETPH
	   phase=0.0
 	   do iii=0,elle
	      phase=phase+x(iii,ileng)*PLG(iii,0,cos(gamma(ind)))
	   enddo
	   diff=abs(phase-phase_in(ind,ileng))/phase_in(ind,ileng)*100.0

	   if(gamma(ind).eq.0.0) then
		write(22,*) ' DIFF AT SC. ANGLE 0         = ', diff
		goto 3322
	   else 
	        if(diff.gt.max_diff) then
		   max_diff = diff
	           max_angle = gamma(ind)*180.0/pi
c                   WRITE(*,*)GAMMA(IND),MAX_ANGLE 
	        endif
c		diff_mean = diff_mean + dabs(diff)
		diff_mean = diff_mean + dabs(DBLE(diff))
		nden = nden + 1
	   endif
 3322	   continue
	   phase_out(ind,ileng) = phase
	enddo

c linear interpolation of the phase function in the ngl points for quadrature

	ind_start=2
	integral = 0.0

	do ind=1,ngl
	   do jjj=ind_start,N_THETPH
		if(gamma(jjj).ge.absci(ind)) then	
			ind_start=jjj
			goto 3333
		endif
	   enddo
 3333	   continue	
	   call inter2(1,absci(ind),gamma(jjj-1),gamma(jjj),
     &		phase,phase_out(jjj-1,ileng),phase_out(jjj,ileng))
	   integral = integral + phase*sin(absci(ind))*weights(ind)
	enddo
         
	integral = integral * 2 * pi
	write(22,*) ' INTEGRAL (RECONSTR. PF)     = ', integral

	diff_mean = diff_mean/ nden
	write(22,1122) diff_mean
	write(22,1133) max_diff
	write(22,1144) max_angle
        enddo
 1122   format('  MEAN DIFF PHASE FUNC.       = ',f10.3,'% ')
 1133   format('  MAX  DIFF PHASE FUNC.       = ',f10.3,'% ')
 1144   format('  AT ANGLE                    = ',f10.3,'% ')

	return
	end

