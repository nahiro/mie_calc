       SUBROUTINE VOL_PERC(NUMAEROS,NRH,C)
C-------------------------------------------------------------
C    This routine computes the mixing ratio by volume 
C    percentage starting from the mixing ratios by number density
C---------------------------------------------------------------
        INCLUDE 'AEROS_DB.PARS'
        INCLUDE 'AEROS_DB.VARS'
       DIMENSION V(4),C(4)
          PI=acos(-1.0)
          E=EXP(1.0)
          VOL=0.0
      DO 300 I=1,VEC_COMP(NUMAEROS)
       ICOMP=COMP_MAT(I,NUMAEROS)
C       S=EXP(4.5*( SIGMA(ICOMP) / LOG10(E) )**2)
C       WRITE(*,*)S
       V(I)=4.0/3.0*PI*(HUM_RADII(NRH,ICOMP))**3*(EXP(4.5*
     +   ( SIGMA(ICOMP) / LOG10(E) )**2))
       VOL=VOL+V(I)*MIX_RATIO(I,NUMAEROS)
 300  CONTINUE
C      WRITE(*,*)'VOL=',VOL
      DO 200 K=1,VEC_COMP(NUMAEROS)
       C(K)=V(K)*MIX_RATIO(K,NUMAEROS)/VOL
C      WRITE(*,*)C
 200  CONTINUE

       RETURN 
       END


       SUBROUTINE RH_RIND(NRH,NWLENGTH,NCOMP,N_RI,RI_R,RI_I)
C-----------------------------------------------------------------
C      This routine computes though Hanel formula components 
C      refractive index at any humidity level
C-----------------------------------------------------------------
        INCLUDE 'AEROS_DB.PARS'
        INCLUDE 'AEROS_DB.VARS'
C
             RI_R=REFIND(NWLENGTH,11,1)+
     +            ((REFIND(NWLENGTH,N_RI,1)-REFIND(NWLENGTH,11,1))*
     +            (HUM_RADII(1,NCOMP)/HUM_RADII(NRH,NCOMP))**3)
             RI_I=REFIND(NWLENGTH,11,2)+
     +            ((REFIND(NWLENGTH,N_RI,2)-REFIND(NWLENGTH,11,2))*
     +            (HUM_RADII(1,NCOMP)/HUM_RADII(NRH,NCOMP))**3)
        RETURN
        END
C
       SUBROUTINE LINTER(WNEW,NJ,NPART,I_RI,RI_IND)
C--------------------------------------------------------------
C     this routine does a linear interpolation used to
C     compute the refractive indices at any wavelength
C     in the range 0.3-0.860
C--------------------------------------------------------------       
        INCLUDE 'AEROS_DB.PARS'
        INCLUDE 'AEROS_DB.VARS'
      REAL WNEW
C
        RI_IND=((REFIND(NJ,I_RI,NPART)-REFIND(NJ-1,I_RI,NPART))/
     +           (WLR(NJ)-WLR(NJ-1))*(WNEW-WLR(NJ)))+
     +           REFIND(NJ,I_RI,NPART)
C
       RETURN
       END
C
       SUBROUTINE WRITEFILE(JPOLY,KPROB,IPROB,NTHETA,WAVE,EXT_COEFF,
     *      SING_S_ALB,ASYM_FAC,THET,XI,ELLE,nscalerad,NUMAEROS)
C
        real pesi(1000), absci(1000),wh(1000)
        PARAMETER (N_WL_OUT=1000)
        INTEGER ELLE,KPROB,JPOLY,IPROB,NTHETA,nscalerad
        REAL WAVE(KPROB),THET(NTHETA),EXT_COEFF(KPROB),
     *       SING_S_ALB(KPROB), ASYM_FAC(KPROB),
     *       XI(0:ELLE,KPROB)
        COMMON /FDIFASE/ FUNZ_FASE(1000,N_WL_OUT)
C
        PI=acos(-1.0)
       IF (IPROB.EQ.1) THEN
        WRITE(2,972)
        WRITE(2,973)
        WRITE(2,974)
        WRITE(2,976)
        WRITE(2,975)
        WRITE(2,977)
       ENDIF
        DEN=1.0
        DO 302 NWL=1,KPROB
c          IF (WAVE(NWL).EQ.0.550) THEN
c             DEN=EXT_COEFF(NWL)
C              DEN=1.0
C              NUMDEN=9088
c          ENDIF
 302    CONTINUE
        NDIV=INT(KPROB/4)
        MDIV=KPROB-(4*NDIV)
         IF (MDIV.EQ.0) THEN
          NDIV=NDIV-1
          MDIV=4
         ENDIF
C        WRITE(*,*)NDIV
C        WRITE(*,*)MDIV
C
C togli a simulazione finita

        if (ntheta.LT.181) then
C
        IF (IPROB.EQ.1) THEN
         DO 303 J=1,NDIV
          WRITE(2,960)(WAVE(NWL),NWL=1+(4*(J-1)),4*J)
 960      FORMAT (5X,'&',4(1X,F7.4,','))
 303     CONTINUE
          WRITE(2,971)(WAVE(NWL),NWL=(4*NDIV)+1,(4*NDIV)+MDIV )
C971      FORMAT (5X,'&',<MDIV>(1X,F7.4,','),TL1,'/')
 971      FORMAT (5X,'&',10000(:,1X,F7.4,','),TL1,'/')
        ENDIF
C
         DO 301 IN=1,KPROB
          WRITE(6,989)WAVE(IN),EXT_COEFF(IN)/DEN,
     + SING_S_ALB(IN),ASYM_FAC(IN)
 301     CONTINUE
          EXT=0.0
         WRITE(2,976)

c       warning: fog models: extinction coefficient to N0 particles/cm3
c       see FORMAT and readme.txt
         if (numaeros.eq.44.or.numaeros.eq.45) then
                write(2,1978)
         else if (numaeros.eq.46) then
                write(2,2978)
         else if (numaeros.eq.47) then
                write(2,3978)
         else
                WRITE(2,978)
         endif

         WRITE(2,979)IPROB
         DO 304 J=1,NDIV
          WRITE(2,961)(EXT_COEFF(NWL)/DEN,NWL=1+(4*(J-1)),4*J)
 961     FORMAT (5X,'&',4(1X,1PE12.5,','))
 304     CONTINUE
         WRITE(2,962)(EXT_COEFF(NWL)/DEN, NWL=(4*NDIV)+1,
     +                      (4*NDIV)+MDIV )
C962     FORMAT (5X,'&',<MDIV>(1X,1PE12.5,','),TL1,'/')
 962     FORMAT (5X,'&',10000(:,1X,1PE12.5,','),TL1,'/')
C
         WRITE(2,976)
         WRITE(2,980)
         WRITE(2,981)IPROB
         DO 305 J=1,NDIV
          WRITE(2,963)(SING_S_ALB(NWL), NWL=1+(4*(J-1)),4*J)
 963     FORMAT (5X,'&',4(1X,0PF8.6,','))
 305     CONTINUE
         WRITE(2,964)(SING_S_ALB(NWL), NWL=(4*NDIV)+1,(4*NDIV)+MDIV)
C964     FORMAT (5X,'&',<MDIV>(1X,0PF8.6,','),TL1,'/')
 964     FORMAT (5X,'&',10000(:,1X,0PF8.6,','),TL1,'/')
C
         WRITE(2,976)
         WRITE(2,982)
         WRITE(2,983)IPROB
         DO 306 J=1,NDIV
          WRITE(2,965)(ASYM_FAC(NWL), NWL=1+(4*(J-1)),4*J)
 965     FORMAT (5X,'&',4(1X,0PF8.6,','))
 306     CONTINUE
         WRITE(2,966)(ASYM_FAC(NWL),NWL=(4*NDIV)+1,(4*NDIV)+MDIV )
C966     FORMAT (5X,'&',<MDIV>(1X,0PF8.6,','),TL1,'/')
 966     FORMAT (5X,'&',10000(:,1X,0PF8.6,','),TL1,'/')
C
          NDIVF=INT(NTHETA/5)
          MDIVF=NTHETA-(5*NDIVF)
           IF (MDIVF.EQ.0) THEN
              NDIVF=NDIVF-1
              MDIVF=5
           ENDIF
C
          IF (IPROB.EQ.1) THEN
            WRITE(2,976)
            WRITE(2,984)
            WRITE(2,985)
          DO 307 J=1,NDIVF
            WRITE(2,967)(THET(NGAMMA), NGAMMA=1+(5*(J-1)),5*J)
 967        FORMAT (5X,'&',5(1X,F6.2,','))
 307      CONTINUE
           IF (MDIVF.NE.0) THEN
            WRITE(2,968)(THET(NGAMMA), NGAMMA=(5*NDIVF)+1,
     +                                  (5*NDIVF)+MDIVF )
C968        FORMAT (5X,'&',<MDIVF>(1X,F6.2,','),TL1,'/')
 968        FORMAT (5X,'&',10000(:,1X,F6.2,','),TL1,'/')
           ELSE
            WRITE(2,988)
           ENDIF
         ENDIF
C
         WRITE(2,976)
           NDIVF=INT(NTHETA/4)
           MDIVF=NTHETA-(4*NDIVF)
            IF (MDIVF.EQ.0) THEN
             NDIVF=NDIVF-1
             MDIVF=4
            ENDIF
         DO 309 NWL=1,KPROB
         WRITE(2,986)WAVE(NWL)
         WRITE(2,987)NWL,IPROB
         DO 308 J=1,NDIVF
       WRITE(2,969)(FUNZ_FASE(NGAMMA,NWL),NGAMMA=1+(4*(J-1)),4*J)
 969     FORMAT (5X,'&',4(1X,1PE12.5,','))
 308     CONTINUE
          IF (MDIVF.NE.0) THEN
        WRITE(2,970)(FUNZ_FASE(NGAMMA,NWL), NGAMMA=(4*NDIVF)+1,
     +                  (4*NDIVF)+MDIVF )
C970        FORMAT (5X,'&',<MDIVF>(1X,1PE12.5,','),TL1,'/')
 970        FORMAT (5X,'&',10000(:,1X,1PE12.5,','),TL1,'/')
          ELSE
            WRITE (2,988)
          ENDIF
 309     CONTINUE

        ENDIF
C togli a simulazione finita
C
C   scrittura coefficienti X(i)   
C 
      IF (JPOLY.EQ.1) THEN
         WRITE(2,976)
           NDIVF=INT((ELLE+1)/4)
           MDIVF=(ELLE+1)-(4*NDIVF)
           IF (MDIVF.EQ.0) THEN
            NDIVF=NDIVF-1
            MDIVF=4
           ENDIF
       DO 1001 Nlung=1,KPROB
         WRITE(2,990)WAVE(NLUNG)
         WRITE(2,991)NLUNG,IPROB
         
         DO 1004 J=1,NDIVF
           WRITE(2,1002)(XI(NELLE,Nlung),NELLE=(4*(J-1)),(4*J)-1 )
 1002      FORMAT (5X,'&',4(1X,1PE12.5,','))
       
 1004    CONTINUE
          IF (MDIVF.NE.0) THEN
            WRITE(2,1003)(XI(NELLE,Nlung), NELLE=(4*NDIVF),
     +                  (4*NDIVF)+MDIVF-1 )
C1003        FORMAT (5X,'&',<MDIVF>(1X,1PE12.5,','),TL1,'/')
 1003        FORMAT (5X,'&',10000(:,1X,1PE12.5,','),TL1,'/')
          ELSE
            WRITE (2,988)
          ENDIF
 1001   CONTINUE    
        ENDIF 
              
C
C
  
  972 FORMAT (6X,'BLOCK DATA AEROSOL_OPT_PROP')
  973 FORMAT (6X,'INCLUDE',1x,' ''AEROSOL_OPT_PROP.PARS''')
  974 FORMAT (6X,'INCLUDE',1x,' ''AEROSOL_OPT_PROP.VARS''')
  975 FORMAT ('C',5X,'Wavelength (microns)')
  976 FORMAT ('C',5X,'---------------------------------------')
  977 FORMAT (6x,'DATA ( aer_wl(i), i=1, n_aer_wl )              /')

  978 FORMAT ('C',5X,'Extintion Coefficient [km-1/cm-3]')
 1978 FORMAT ('C',5X,'Extintion Coefficient [km-1/(20*cm-3)]')
 2978 FORMAT ('C',5X,'Extintion Coefficient [km-1/(100*cm-3)]')
 3978 FORMAT ('C',5X,'Extintion Coefficient [km-1/(200*cm-3)]')
  979 FORMAT (6x,'DATA ( extval(i,',I2,'), i=1, n_aer_wl )       /')
  980 FORMAT ('C',5X,'Single Scattering Albedo')
  981 FORMAT (6x,'DATA ( omval(i,',I2,'), i=1, n_aer_wl )        /')
  982 FORMAT ('C',5X,'Asymmetry Factor')
  983 FORMAT (6x,'DATA ( asym_fac(i,',I2,'), i=1,n_aer_wl )      /')
  984 FORMAT ('C',5X,'Scattering Angle')
  985 FORMAT (6x,'DATA ( thetph(j), j=1, n_thetph )              /')
  986 FORMAT ('C',5X,'Phase Function at wavelength=', F7.4)
  987 FORMAT (6x,'DATA ( phsfnc(j,',I2,',',I2,'), j=1,n_thetph ) /')
  988 FORMAT (5X,'&',3X,'/')
  990 FORMAT ('C',5X,'Legendre polynomial terms at wavelength=',F7.4)
  991 FORMAT (6x,'DATA ( x(i,',I2,',',I2,'), i=1, elle )         /')

  989 FORMAT (/,2X,'WAVELENGTH =',F8.4,/,2X,
     *          'MIXTURE EXTINCTION COEFFICIENT =',1PE12.5,/,2X,
     *          'MIXTURE SINGLE SCATTER. ALBEDO =',0PF8.6,/,2X,
     *          'ASYMMETRY PARAMETER   = ',0PF8.6,/,2X)
        RETURN
        END

C<FF>
      SUBROUTINE GAULEG(X1,X2,X,W,N)
C
C      Routine copyed by
C  *Numerical Recepies in Fortran*
C
C  Given a lower and a upper limits of integration X1 and X2,
C  and given N, this routine returns arrays X and W of lenght N,
C  containing the abscissas and weights of the Gauss-Legendre
C  N-point quadrature formula.
C
      IMPLICIT REAL (A-H,O-Z)
      REAL X1,X2,X(N),W(N)
      PARAMETER (EPS=3.D-14)
      DOUBLE PRECISION P1,P2,P3,PP,XL,XM,Z,Z1
      M=(N+1)/2
      XM=0.5D0*(X2+X1)
      XL=0.5D0*(X2-X1)
      DO 12 I=1,M
       Z=COS(3.141592654D0*(I-0.25D0)/(N+.5D0))
 1      CONTINUE
        P1=1.D0
        P2=0.D0
        DO 11 J=1,N
         P3=P2
         P2=P1
         P1=((2.D0*J-1.D0)*Z*P2-(J-1.D0)*P3)/J
 11     CONTINUE
C
C P1 is the desired Legendre polynomial.
C
        PP=N*(Z*P1-P2)/(Z*Z-1.D0)
        Z1=Z
        Z=Z1-P1/PP
        IF(ABS(Z-Z1).GT.EPS) GO TO 1
        X(I)=XM-XL*Z
        X(N+1-I)=XM+XL*Z
        W(I)=2.D0*XL/((1.D0-Z*Z)*PP*PP)
        W(N+1-I)=W(I)
 12   CONTINUE
C
      RETURN
      END
C

C
      SUBROUTINE RRA42 (X,RENT1,RENT2,ENT3,ENT4,EXKROS,SCKROS,COSAV,PN)
C
      COMPLEX W(3),A(13000),SA(13000),SB(13000),S1(1000),S2(1000),PM,WO,
     *   WM1,S2CONJ,S1S2CJ,PMX,SA1,SB1
      COMPLEX FF,CTAN,Z,SPLUS,SMINUS,ACNJ,SBL
      DIMENSION RENT1(*),RENT2(*),ENT3(*),ENT4(*),PI(3)
      COMMON PM,NTHETA,SINS(1000),COS2(1000),COSA(1000),APM,A
C
C     MAXDIM SHOULD EQUAL DIMENSION OF A
C
      DATA MAXDIM / 13000 /
C
      J = 1
      J1 = 2
      J2 = 3
C
      SN = SIN(X)
      CN = COS(X)
      PMX = PM*X
      NN = 1.01*(X+4.05*X**(1./3.)+2.0) 
      NP = NN+2
C
C     TEST NUMBER TERMS REQUIRED AGAINST THE DIMENSION OF A 
C
      IF (NP.GT.MAXDIM) THEN
         WRITE (6,900) NP,X,PM
         STOP
      ENDIF
C
C     INITIALIZE PARAMETERS FOR RECURSION
C
      ROX = 1./X
      COSAV = 0.
      EXKROS = 0.
      SCKROS = 0.
      COEPH = 2./(X*X)
C
      WO = CMPLX(SN,CN)
      WM1 = CMPLX(CN,-SN)
C
      W(1) = (1./X)*WO-WM1
      W(2) = (3./X)*W(1)-WO
C
      WREAL = REAL(W(1))
C
      PMIMAG = ABS(AIMAG(PM)) 
      PMREAL = REAL(PM)
      XMIM = PMIMAG*X
      F2MRE = 13.78*PMREAL**2-10.8*PMREAL+3.9
      ZF = 1.0/X
      Z = ZF/PM
      FF = CTAN(PM*X)
      IF (XMIM.LE.F2MRE) THEN 
C
C        DO UPWARD RECURSION
C
         CALL UPREC (Z,FF,NP,A)
      ELSE
C
C        DO DOWNWARD RECURSION
C
         CALL DWNREC (Z,ASTART,NP,A)
      ENDIF
C
C        START UPWARD RECURSION FOR ATTENUATION CROSS-SECTIONS
C
      WREALO = REAL(WO)
      SA(1) = ((A(1)/PM+ROX)*WREAL-WREALO)/((A(1)/PM+ROX)*W(1)-WO)
      SB(1) = ((A(1)*PM+ROX)*WREAL-WREALO)/((A(1)*PM+ROX)*W(1)-WO)
C
      DO 10 N = 1, NN
         N1 = N+1
         PN = FLOAT(N)
         PN2 = FLOAT(N+2)
         PX1 = FLOAT(N1)/X
C
         W(J2) = ((2.*PN2-1.)/X)*W(J1)-W(J)
         WREAL = REAL(W(J1))
         WREAL1 = REAL(W(J))
         SA1 = A(N1)/PM+PX1
         SB1 = A(N1)*PM+PX1
         SA(N1) = (SA1*WREAL-WREAL1)/(SA1*W(J1)-W(J))
         SB(N1) = (SB1*WREAL-WREAL1)/(SB1*W(J1)-W(J))
C
         RSAB = CABS(SA(N))**2+CABS(SB(N))**2
C
         W(1) = W(2)
         W(2) = W(3)
         CONST = 2.*PN+1.
         DSC = CONST*RSAB
         DEX = CONST*REAL(SA(N)+SB(N))
         EXKROS = EXKROS+DEX
         SCKROS = SCKROS+DSC
         IF (N.NE.1) COSAV = (PN-1./PN)*REAL(ACNJ*SA(N)+CONJG(SBL)*SB(N)
     *      )+((CONST-2.0)/(PN*PN-PN))*REAL(ACNJ*SBL)+COSAV 
         ACNJ = CONJG(SA(N))
         SBL = SB(N)
C
   10 CONTINUE
C
      KEEP = IFIX(PN)
      EXKROS = EXKROS*COEPH
      COSAV = 2.0*COSAV/SCKROS
      SCKROS = SCKROS*COEPH
C
      DO 30 K = 1, NTHETA
         CTHETA = COSA(K)
         SPLUS = 1.5*(SA(1)+SB(1))*(1.+CTHETA)
         SMINUS = 1.5*(SA(1)-SB(1))*(1.-CTHETA)
C
         PI(1) = 1. 
         PI(2) = 3.*CTHETA
         DO 20 N = 2, KEEP
            FN = FLOAT(N)
            FN1 = FLOAT(N+1)
            S = CTHETA*PI(J1) 
            T = S-PI(J)
            PI(J2) = S+T*(FN1/FN)
            TAU = FN*T-PI(J)
            FACT = (FN+FN1)/(FN*FN1)
            SPLUS = FACT*(SA(N)+SB(N))*(PI(J1)+TAU)+SPLUS
            SMINUS = FACT*(SA(N)-SB(N))*(PI(J1)-TAU)+SMINUS 
            PI(1) = PI(2)
            PI(2) = PI(3)
   20    CONTINUE
         S1(K) = 0.5*(SPLUS+SMINUS)
         S2(K) = 0.5*(SPLUS-SMINUS)
   30 CONTINUE
C
      DO 40 I = 1, NTHETA
         S2CONJ = CONJG(S2(I))
         S1S2CJ = S1(I)*S2CONJ
         RENT1(I) = REAL(S1(I)*CONJG(S1(I)))
         RENT2(I) = REAL(S2(I)*S2CONJ)
         ENT3(I) = REAL(S1S2CJ)
         ENT4(I) = -AIMAG(S1S2CJ)
   40 CONTINUE
C
      RETURN
C
C
C
C
C
  900 FORMAT (' NMAX =',I20,' EXCEEDING DIMENSION OF A',/,
     *   ' X,REF INDEX=',1P3G15.5)
C
      END 
C
C
      SUBROUTINE IPLAW (KEEP,H,F,AR)
C
      DIMENSION H(1000),F(1000),AR(1000)
C
      KEEB = KEEP-1 
      DO 30 I = 1, KEEB
         IF (ABS(H(I)).EQ.0.0.OR.ABS(H(I+1)).EQ.0.0) GO TO 10
         HF = H(I)*F(I)
         HF1 = H(I+1)*F(I+1)
         IF (ABS(HF).EQ.0.0.OR.ABS(HF1).EQ.0.0) GO TO 10
         HII1 = H(I)/H(I+1)
         HFI1 = HF/HF1
         IF (HII1.LE.0.0.OR.HFI1.LE.0.0) GO TO 10 
         ALPH1 = ALOG(HFI1)/ALOG(HII1)
         AR(I+1) = (HF1-HF)/ALPH1
         GO TO 20
   10    AR(I+1) = 0.5*(F(I)+F(I+1))*(H(I+1)-H(I))
   20    AR(I+1) = -AR(I+1)
   30 CONTINUE
C
      RETURN
      END 
      SUBROUTINE UPREC (Z,FF,N,AN)
C
C     PERFORM UPWARD RECURSION FOR ATTENUATION CROSS-SECTIONS
C
      COMPLEX Z,AN(3400),FF
      AN(1) = -Z+FF/(Z*FF-1.0)
      DO 10 L = 2, N
         AN(L) = -(L*Z)+1.0/((L*Z)-AN(L-1))
   10 CONTINUE
      RETURN
      END 
      SUBROUTINE DWNREC (Z,ASTART,N,A)
C
C     CALCULATE A(N) IN FUNCTION ALENTZ THEN PERFORM DOWNWARD RECURSION
C
      COMPLEX ASTART,Z,A(13000),ALENTZ
      ASTART = ALENTZ(Z,N)
      A(N) = ASTART 
      DO 10 I = N, 2, -1
         A(I-1) = (I*Z)-1.0/((I*Z)+A(I))
   10 CONTINUE
      RETURN
      END 
      COMPLEX FUNCTION CTAN (Z)
C
C     COMPLEX TANGENT OF Z
C
      COMPLEX Z
C      REAL*8 E2Y, E4Y,TWOX
C TWOX = 2.d0*REAL(Z)
      TWOX = 2.0*REAL(Z)
      IF (AIMAG(Z).GT.-337.9) THEN
         E2Y = EXP(2.0*AIMAG(Z))
      ELSE
         E2Y = 0.0
      ENDIF
      E4Y = E2Y**2
C      DEN = 1./((2.d0*E2Y)*DCOS(TWOX)+E4Y+1.)
C      CTAN = DCMPLX((2.*E2Y)*DSIN(TWOX)*DEN,(E4Y-1.)*DEN)
         DEN = 1./((2.*E2Y)*COS(TWOX)+E4Y+1.)
         CTAN = CMPLX((2.*E2Y)*SIN(TWOX)*DEN,(E4Y-1.)*DEN)
      RETURN
      END
C
      FUNCTION ALENTZ (Z,NTERM)
C
C     USES LENTZ'S PRODUCTION REPRESENTATION (APPL. OPT., 15, P.668,
C     1976) OF THE CONTINUED FRACTION FORM FOR A SUB N TO CALCULATE
C     A STARTING VALUE OF A(N) FOR USE IN A DOWNWARD RECURSION
C
      COMPLEX FF,TNTN,DTD
      COMPLEX ALENTZ,Z,ZINV,AK,DEN,TNUM,TT
      DATA MAXIT / 10000 /
      DATA EPS1 / 1.0E-02 / ,EPS2 / 1.0E-08 /
      ZINV = Z
      FF = (NTERM+1)*ZINV
      MM = -1
      KK = 2*NTERM+3
      AK = (MM*KK)*ZINV
      DEN = AK
      TNUM = AK+1.0/FF
      ITER = 1
   10 ITER = ITER+1 
C
C     IF TOO MANY ITERATIONS, EXIT WITH COMMENT
C     IF ILL CONDITIONED CASE, SKIP CONVERGE CHECK AND DO A DOUBLE STEP
C
      IF (ITER.GT.MAXIT) THEN 
         WRITE (6,900) NTERM,Z,AK,TNUM,DEN,TT,FF
         STOP '8001'
      ENDIF
      IF (CABS(TNUM/AK).GT.EPS1.AND.CABS(DEN/AK).GT.EPS1) THEN
         TT = TNUM/DEN
         FF = TT*FF 
C
C        CONVERGENCE CHECK
C
         IF (ABS(REAL(TT)-1.0).LT.EPS2.AND.ABS(AIMAG(TT)).LT.EPS2) GO TO
     *       20
         MM = -MM
         KK = KK+2
         AK = (MM*KK)*ZINV
         TNUM = AK+1.0/TNUM
         DEN = AK+1.0/DEN
      ELSE
         MM = -MM
         KK = KK+2
         AK = (MM*KK)*ZINV
         TNTN = AK*TNUM+1.0
         DTD = AK*DEN+1.0
         FF = (TNTN/DTD)*FF
         MM = -MM
         KK = KK+2
         AK = (MM*KK)*ZINV
         TNUM = AK+TNUM/TNTN
         DEN = AK+DEN/DTD
         ITER = ITER+1
      ENDIF
      GO TO 10
C
   20 CONTINUE
      ALENTZ = FF
      RETURN
C
C
C
C
C
C
C
C
C
C
  900 FORMAT (//,'  CONTINUED FRACTION FOR A(N) FAILED TO CONVERGE',/,
     *   'NTERM, Z, AK, TNUM, DEN, TT, FF =',//,2X,I6,6(1P2E10.4))
C
      END 
C
C
      FUNCTION PROD (R,ENR,NRS,DM,MOM,MODE)
C
C     PROD FINDS THE MOMENTS OF THE SIZE DISTRIBUTION ENR
C     PROD = INTEGRAL OF ENR*R**MOM  FROM R(1) TO R(NRS)
C     WITH RESPECT TO THE RADI R
C     INTEGRATES OVER R, ASSUMING A POWER LAW BETWEEN ADJACENT POINTS 
C     MOM IS THE ORDER OF THE MOMENT
C     IF MODE = 0, DM IS NOT USED
C     IF MODE = 1, DM(I) IS THE CONTRIBUTION TO PROD CENTERED AT R(I) 
C     IF MODE = 2, DM(I) IS THE CUMULATIVE INTEGRAL TO A POINT
C     HALFWAY BETWEEN R(I) AND R(I+1)
C
      DIMENSION R(NRS),ENR(NRS),DM(NRS) 
      PROD = 0.
      print*,'mode,nrs',mode,nrs
      IF (MODE.GE.1) DM(1) = 0.
      NRS1 = NRS-1
      print*,'nrs1,r(1)',nrs1,r(1)
      IF (R(1).EQ.0.) Y2 = 0. 
      IF (R(1).GT.0.) Y2 = ENR(1)*R(1)**MOM
C
      print*,'y2',y2
      DO 10 I = 1, NRS1
      print*,'i',i
         Y1 = Y2
         Y2 = ENR(I+1)*R(I+1)**MOM
         Z1 = Y1*R(I)
         Z2 = Y2*R(I+1)
         IF (ENR(I).EQ.0.OR.ENR(I+1).EQ.0.OR.R(I).EQ.0.OR.ABS(Z1-Z2).LE.
     *      ABS(Z1)*1.E-6) THEN
C
C           USE LINEAR FIT BECAUSE POWER LAW CANNOT BE USED 
C
            DELR = R(I+1)-R(I)
            IF (MODE.GE.1) THEN
               DM(I) = DM(I)+(Y2+3*Y1)*DELR/8.
               DM(I+1) = (3*Y2+Y1)*DELR/8.
               IF (MODE.EQ.2) DM(I+1) = DM(I+1)+DM(I)
               IF (MODE.EQ.1) PROD = PROD+DM(I)
            ELSE
               PROD = PROD+(Y1+Y2)*DELR/2.
            ENDIF
         ELSE
         print*,'z2,z1,r(i,i+1)',z2,z1,r(i),r(i+1)
            ALPH1 = ALOG(Z2/Z1)/ALOG(R(I+1)/R(I)) 
         print*,'mode,alph1',mode,alph1
            IF (MODE.GE.1) THEN
               Z = (Z1/R(I)**ALPH1)*((R(I+1)+R(I))/2.)**ALPH1
          print*,'z',z
               DM(I) = DM(I)+(Z-Z1)/ALPH1
               DM(I+1) = (Z2-Z)/ALPH1
         print*,'dm',dm(i),dm(i+1)
               IF (MODE.EQ.2) DM(I+1) = DM(I+1)+DM(I)
               IF (MODE.EQ.1) PROD = PROD+DM(I)
            ELSE
               PROD = PROD+(Z2-Z1)/ALPH1
            ENDIF
         ENDIF
   10 CONTINUE
C
      IF (MODE.EQ.1) PROD = PROD+DM(NRS)
      IF (MODE.EQ.2) PROD = DM(NRS)
      RETURN
C
C
C
C
C
C
C
      END 
C
C
      SUBROUTINE NEWRAD (RI,ENR,NRS,MOM)
      DIMENSION RI(800),ENR(800),AREA(800),RIP(800)
C
C     GENERATE A NEW SET OF RADII UNIFORMLY SPACED IN N(R)**MOM
C     NOTE "AREA" REPRESENTS THE APPROPRIATE MOMENT OF N(R) 
C
      DATA PII / 3.14159265 / 
      NRS1 = NRS-1
      print*,'nrs,nrs1',nrs,nrs1
      ARETOT = PII*PROD(RI,ENR,NRS,AREA,MOM,2)
      DO 10 I = 1, NRS
         AREA(I) = AREA(I)/AREA(NRS)
   10 CONTINUE
      NDEL = AMIN1(SQRT(FLOAT(NRS)),15.) 
      NRW2 = NRS1-NDEL
      ARD = AREA(NRS)/FLOAT(NRS-2*NDEL+1)
      AS = ARD/2**(NDEL-1)
      KL = 1
      DO 30 I = 2, NRS
         IF (I.LE.NDEL) AS = AS+AS
         IF (I.GT.NDEL.AND.I.LE.NRW2) AS = AS+ARD 
         IF (I.GT.NRW2) AS = AS+ARD/2**(I-NRW2)
         DO 20 KM = KL, NRS
            IF (AREA(KM).GE.AS) THEN
               KL = KM
               IF (KM.EQ.1) THEN
                  RIP(I) = RI(1)+0.5*(RI(2)-RI(1))*AS/AREA(1)
               ELSE 
                  R1 = 0.5*(RI(KM-1)+RI(KM))
                  IF (KM.LT.NRS) R2 = 0.5*(RI(KM+1)+RI(KM)) 
                  IF (KM.EQ.NRS) R2 = RI(KM)
                  AL = ALOG(R1/R2)/ALOG(AREA(KM-1)/AREA(KM))
                  RIP(I) = R2*(AS/AREA(KM))**AL
               ENDIF
               GO TO 30
            ENDIF
   20    CONTINUE
         RIP(I) = RI(NRS)
         KL = NRS1
   30 CONTINUE
C
      DO 40 I = 2, NRS1
         RI(I) = (RIP(I)+RIP(I+1))/2.
   40 CONTINUE
      RETURN
C
C
C
C
C
C
C
      END 
      SUBROUTINE RAYLIM (X,RENT1,RENT2,ENT3,ENT4,EXKRO,SCKRO,COSAV)
C
C     CALCULATES THE MATRIX ELEMENTS, MIE SCATTERING AND ABSORPTION
C     EFFICIENCIES IN THE RAYLEIGH LIMIT
C
      DIMENSION RENT1(*),RENT2(*),ENT3(*),ENT4(*) 
      COMPLEX PM,FAC
      COMMON PM,NTHETA,SINS(1000),COS2(1000),COSA(1000)
C
      FAC = (PM**2-1.)/(PM**2+2.)
      SCKRO = (8./3.)*X**4*CABS(FAC)**2 
      ABKRO = -4*X*AIMAG(FAC) 
      EXKRO = SCKRO+ABKRO
      COSAV = 0.
C
C     COMPUTE THE MATRIX ELEMENTS
C
      DO 10 I = 1, NTHETA
         RENT1(I) = CABS(FAC)**2*X**6
         RENT2(I) = COSA(I)**2*CABS(FAC)**2*X**6
         ENT3(I) = COSA(I)*CABS(FAC)**2*X**6
         ENT4(I) = 0.
   10 CONTINUE
  
      RETURN
      END 
      SUBROUTINE RAYOPT (X,NP,THET,RENT1,RENT2,ENT3,ENT4)
C
C     MATRIX ELEMENTS FOR A SPHERE USING GEOMETRICAL OPTICS 
C
      REAL GDIFF(1000),GAIN1(1000),GAIN2(1000),J1,THETPR(181),THETA(10)
      DIMENSION THET(*),RENT1(*),RENT2(*),ENT3(*),ENT4(*)
      COMPLEX PM
      INTEGER P,PNAUT
      COMMON PM,NTHETA
      DATA PI,RADS / 3.14159265,0.017453292 /
C
      ZN = REAL(PM) 
      ZK = -AIMAG(PM)
      TRUNC = 1E-5
C
C     DIFFRACTION GAIN
C
      DO 10 I = 1, NTHETA
         IF (THET(I).EQ.0.) THEN
            GDIFF(I) = 0.5*X**2
         ELSEIF (THET(I).GT.90) THEN
            GDIFF(I) = 0.
         ELSE
            THETR = THET(I)*RADS
            SINT = SIN(THETR) 
            GDIFF(I) = 2.*(J1(X*SINT)/SINT)**2
         ENDIF
   10 CONTINUE
C
      DO 20 I = 1, NTHETA
         GAIN1(I) = GDIFF(I)
         GAIN2(I) = GDIFF(I)
   20 CONTINUE
C
C     GAIN FOR GEOMETRICAL OPTICS; FIRST EXTERNAL REFLECTION
C
      P = 0
      DO 30 I = 1, NTHETA
         THETA(1) = THET(I)*RADS
         TAU = THETA(1)/2.
         TAUPR = ACOS(COS(TAU)/ZN)
         GAIN1(I) = GAIN1(I)+GAINP(P,1,X,TAU,TAUPR,THETA(1))
         GAIN2(I) = GAIN2(I)+GAINP(P,2,X,TAU,TAUPR,THETA(1))
   30 CONTINUE
C
C     GAIN FOR P <= ZN
C
      PNAUT = INT(ZN)
      IF (PNAUT.GT.NP) PNAUT = NP
      DO 80 P = 1, PNAUT
         DO 40 NT = 1, 180
            TAU = (NT-1)*0.5*RADS
            TAUPR = ACOS(COS(TAU)/ZN)
            THETPR(NT) = 2.*TAU-2.*P*TAUPR
   40    CONTINUE
         THETPR(181) = (1-P)*PI
         DO 70 I = 1, NTHETA
            THETR = THET(I)*RADS
            TAU = -99.
            IF (MOD(P,2).EQ.1) THETA(1) = (1-P)*PI-THETR
            IF (MOD(P,2).EQ.0) THETA(1) = -P*PI+THETR
            IF (MOD(P,2).EQ.1.AND.THETA(1).LT.THETPR(1)) GO TO 80
            DO 50 NT = 0, 179 
               K = 181-NT
               IF (THETA(1)-THETPR(K).GE.-TRUNC) THEN
                  S = (THETA(1)-THETPR(K))/(THETPR(K-1)-THETPR(K))
                  TAU = ((K-1)-S)*0.5*RADS
                  IF (ZN.EQ.FLOAT(P).AND.ABS(TAU-PI/2).LE.TRUNC) THEN 
                     TAU = TAU-2.*RADS
                     TPR = ACOS(COS(TAU)/ZN)
                     T = 2.*TAU-2.*P*TPR
                     THETR = ABS(ABS(T-THETA(1))-THETR)
                  ENDIF
                  GO TO 60
               ENDIF
   50       CONTINUE
            IF (ABS(THETA(1)-THETPR(1)).LE.TRUNC) TAU = 0.
   60       IF (TAU.NE.-99) THEN
               TAUPR = ACOS(COS(TAU)/ZN)
               GAIN1(I) = GAIN1(I)+GAINP(P,1,X,TAU,TAUPR,THETR)
               GAIN2(I) = GAIN2(I)+GAINP(P,2,X,TAU,TAUPR,THETR)
            ENDIF
   70    CONTINUE
   80 CONTINUE
C
C     GAINS FOR ALL OTHER P'S (SLOW!!)
C
      IF (PNAUT.EQ.NP) GO TO 170
      DO 160 P = PNAUT+1, NP
         RAINBO = ASIN(SQRT((ZN*ZN-1.)/(P*P-1.))) 
         THETBO = 2.*RAINBO-2.*P*ACOS(COS(RAINBO)/ZN)
         DO 90 NT = 1, 180
            TAU = (NT-1)*0.5*RADS
            TAUPR = ACOS(COS(TAU)/ZN)
            THETPR(NT) = 2.*TAU-2.*P*TAUPR
   90    CONTINUE
         THETPR(181) = (1-P)*PI
         DO 150 I = 1, NTHETA 
            THETR = THET(I)*RADS
            IF (MOD(P,2).EQ.1) THEN
               DO 100 N = 1, 5, 2
                  THETA(N) = (1-P)*PI-THETR+PI*(N-1)
                  THETA(N+1) = (1-P)*PI+THETR+PI*(N-1)
  100          CONTINUE
            ELSE
               DO 110 N = 1, 5, 2
                  THETA(N) = -P*PI+THETR+PI*(N-1) 
                  THETA(N+1) = (2-P)*PI-THETR+PI*(N-1)
  110          CONTINUE
            ENDIF
            DO 140 N = 1, 6
               IF (THETA(N)-THETBO.GT.TRUNC) GO TO 150
               IF (N.EQ.1.AND.THETPR(181).LE.THETPR(1)) GO TO 140
               IF (ABS(THETA(N)-THETA(N+1)).LE.TRUNC) GO TO 140
               DO 130 NT = 1, 180
                  IF (AMAX1(THETA(N),THETPR(NT),THETPR(NT+1)).EQ.THETA(N
     *               )) GO TO 130
                  IF (AMIN1(THETA(N),THETPR(NT),THETPR(NT+1)).EQ.THETA(N
     *               ).AND.THETA(N).NE.THETPR(NT).AND.THETA(N).NE.THETPR
     *               (NT+1)) GO TO 130
                  S = (THETA(N)-THETPR(NT))/(THETPR(NT+1)-THETPR(NT)) 
                  TAU = ((NT-1)+S)*0.5*RADS
  120             IF (((ABS(THETR-PI).LE.TRUNC.OR.ABS(THETR).LE.TRUNC)
     *               .AND.ABS(TAU-PI/2).GT.TRUNC).OR.TAU.EQ.RAINBO) THEN
                     TAU = TAU-2.*RADS
                     IF (TAU.LT.0) TAU = TAU+4.*RADS
                     TPR = ACOS(COS(TAU)/ZN)
                     T = 2.*TAU-2.*P*TPR
                     THETR = ABS(ABS(T-THETA(N))-THETR)
                     GO TO 120
                  ENDIF
                  TAUPR = ACOS(COS(TAU)/ZN)
                  GAIN1(I) = GAIN1(I)+GAINP(P,1,X,TAU,TAUPR,THETR)
                  GAIN2(I) = GAIN2(I)+GAINP(P,2,X,TAU,TAUPR,THETR)
  130          CONTINUE
  140       CONTINUE
  150    CONTINUE
  160 CONTINUE
  170 CONTINUE
C
C     CONVERT THE GAINS TO MATRIX ELEMENTS
C
      DO 180 I = 1, NTHETA
         RENT1(I) = GAIN1(I)/2.*X**2
         RENT2(I) = GAIN2(I)/2.*X**2
         ENT3(I) = SQRT(GAIN1(I)*GAIN2(I))/2.*X**2
         ENT4(I) = 0.
  180 CONTINUE
      RETURN
      END 
      FUNCTION GAINP (P,NPOL,X,TAU,TAUPR,THETR)
C
C     FUNCTION TO COMPUTE THE GAIN, GAINP
C
      COMPLEX PMJ,PM,R
      INTEGER P
      COMMON PM
C
      ZN = REAL(PM) 
      ZK = -AIMAG(PM)
      PMJ = CONJG(PM)
      SINT = SIN(TAU)
C
C     FRESNEL REFLECTION COEFFICIENT
C
      IF (NPOL.EQ.1) R = (SINT-PMJ*SIN(TAUPR))/(SINT+PMJ*SIN(TAUPR))
      IF (NPOL.EQ.2) R = (PMJ*SINT-SIN(TAUPR))/(PMJ*SINT+SIN(TAUPR))
      RSQR = CABS(R)**2
      IF (P.EQ.0) THEN
         ESQR = RSQR
      ELSE
         IF (RSQR.EQ.0.) THEN 
            ESQR = 0.
         ELSE
            ESQR = (1.-RSQR)**2*(RSQR)**(P-1)
         ENDIF
      ENDIF
C
      AB = EXP(-4.*X*ZK*P*SIN(TAUPR))
      IF (ABS(TAU-ACOS(0.)).LE.1E-5) THEN
         GAINP = 0.5*ESQR/(1-P/ZN)**2*AB
      ELSEIF (TAU.EQ.0.) THEN 
         IF (P.EQ.0) THEN
            GAINP = 0.5
         ELSE
            GAINP = 0.0
         ENDIF
      ELSE
         DTDTAU = 2.-2.*P*TAN(TAU)/TAN(TAUPR)
         D = SINT*COS(TAU)/SIN(THETR)/ABS(DTDTAU) 
         GAINP = 2.*ESQR*D*AB 
         IF (GAINP.GT.100.) GAINP = 100 
      ENDIF
      RETURN
      END 
      FUNCTION J1 (Z)
C
C     RETURNS THE BESSEL FUNCTION OF FIRST ORDER
C
      REAL J1
      DATA PI / 3.14159265 /
      IF (Z.LE.15.0) THEN
         DEL = 0.05*PI
         B = 1.
         BESS = 0.
         DO 10 I = 1, 20
            A = B
            THETA = DEL*I
            ARG = Z*SIN(THETA)-THETA
            B = COS(ARG)
            DI = (A+B)*DEL/2. 
            BESS = BESS+DI
   10    CONTINUE
         J1 = BESS/PI
      ELSE
         J1 = SQRT(2./PI/Z)*COS(Z-0.75*PI)
      ENDIF
      RETURN
      END 
      REAL FUNCTION QEXAPX(PMCONJ,X)
C     ******************************************************************
C
C     ROUTINE TO CALCULATE THE MIE EXTINCTION EFFICIENCY FOR LARGE
C     SIZE PARAMETERS USING MODIFIED ANAMOLOUS DIFFRACTION THEORY
C
C     REFERENCE: ACKERMAN AND STEPHENS  JAS, VOL 44 P. 1574 
C
C     ******************************************************************
      COMPLEX PMCONJ
C
C     DEFINE SOME CONSTANTS
C
      ZN = REAL(PMCONJ)
      ZK = AIMAG(PMCONJ)
      BETA = ATAN(ZK/(ZN-1.)) 
      RHO = 2.*X*(ZN-1.)
C
      IF (RHO*TAN(BETA).GE.600.) TERM1 = 0.
      IF (RHO*TAN(BETA).LT.600.) TERM1 = 4.*EXP(-RHO*TAN(BETA))
      TERM2 = COS(BETA)/RHO
      U = 4.*X*ZK
      V = 4.*X*ZK*SQRT(1.-1./ZN**2.)
C
C     CALCULATE QEXAPX
C
      QEXAPX = 2.-TERM1*TERM2*SIN(RHO-BETA)-TERM1*TERM2**2.*COS(RHO-2.
     * *BETA)+4.*TERM2**2.*COS(2.*BETA)+2.*X**(-2./3.)
      RETURN
      END 
      REAL FUNCTION QABAPX (CREFIN,SIZPAR,GAMMA)
C
C     ******************************************************************
C
C     MIE ABSORPTION EFFICIENCY FROM ASYMPTOTIC FORMULAS
C     (GEOMETRIC OPTICS PLUS EDGE CORRECTIONS)
C
C     ROUTINE OBTAINED FROM WARREN WISCOMBE, GSFC/NASA, GREENBELT, MD.
C
C     REFERENCE:  NUSSENZWEIG, H. AND W. WISCOMBE, 1980, EFFICIENCY
C                 FACTORS IN MIE SCATTERING, PHYS. REV. LETT. 45,
C                 1490-1494
C
C     ACCURACY:   BETTER THAN 1 PERCENT FOR SIZE PARAMETER .GT. 100
C                 EXCEPT IN RIPPLE SPIKES
C
C     ******************************************************************
C
      COMPLEX CN,CREFIN
      REAL LOLIM,RINT(8),XX,GAM,SIZPAR,GAMMA
      COMMON /APPXMI/ CN,XX,GAM
      EXTERNAL QAGOMI,QAAEMI,QABEMI
      SAVE PID2
      DATA PID2 / 0.0 /
C
      IF (PID2.EQ.0.0) PID2 = ASIN(1.0) 
      CN = CREFIN
      XX = SIZPAR
      GAM = GAMMA
      IF (AIMAG(CN)*XX.LT.1.E-7) THEN
         QABAPX = 0.0
         RETURN
      ENDIF
C
C     GEOMETRIC OPTICS CONTRIBUTION
C
      LOLIM = 0.0
      UPLIM = PID2
      ACCUR = 1.E-6 
      CALL RKRONQ (LOLIM,UPLIM,QAGOMI,ACCUR,RINT,KKK)
      IF (KKK.EQ.0) CALL RINTER (KKK,XX,CN,'QABAPX','QAGOMI',LOLIM,UPLIM
     *   ,ACCUR,RINT)
      GEOMOP = RINT(KKK)
C
C     EDGE CONTRIBUTIONS
C
      TMAX = 1.0/GAM
      XMAX = (REAL(CN)-1.0)/(0.5*GAM**2)
      XUP = AMIN1(XMAX,6.0)
C
      QABAPX = GEOMOP+(0.5*GAM**2)*EDGINT(QAAEMI,QABEMI,XUP,2.,TMAX)
C
      RETURN
      END 
      REAL FUNCTION QPRAPX (CREFIN,SIZPAR,GAMMA)
C
C     ******************************************************************
C
C     MIE RADIATION PRESSURE EFFICIENCY FROM ASYMPTOTIC FORMULAS
C     (GEOMETRIC OPTICS PLUS EDGE CORRECTIONS)
C
C     ROUTINE OBTAINED FROM WARREN WISCOMBE, GSFC/NASA, GREENBELT, MD.
C
C     REFERENCE:  NUSSENZWEIG, H. AND W. WISCOMBE, 1980: EFFICIENCY
C                 FACTORS IN MIE SCATTERING, PHYS. REV. LETT. 45,
C                 1490-1494
C
C     ACCURACY:   BETTER THAN 1 PERCENT FOR SIZE PARAMETER .GT. 60
C                 EXCEPT IN RIPPLE SPIKES
C
C     ******************************************************************
C
      COMPLEX CN,CREFIN
      REAL LOLIM,RINT(8),XX,GAM,SIZPAR,GAMMA
      COMMON /APPXMI/ CN,XX,GAM
      EXTERNAL QPGOMI,QPAEMI,QPBEMI
      SAVE PID2
      DATA PID2 / 0.0 /
C
      IF (PID2.EQ.0.0) PID2 = ASIN(1.0) 
      CN = CREFIN
      XX = SIZPAR
      GAM = GAMMA
C
C     GEOMETRIC OPTICS CONTRIBUTION
C
      LOLIM = 0.0
      UPLIM = PID2
      ACCUR = 1.E-6 
      CALL RKRONQ (LOLIM,UPLIM,QPGOMI,ACCUR,RINT,KKK)
      IF (KKK.EQ.0) CALL RINTER (KKK,XX,CN,'QPRAPX','QPGOMI',LOLIM,UPLIM
     *   ,ACCUR,RINT)
      GEOMOP = RINT(KKK)
C
C     EDGE CONTRIBUTIONS
C
      TMAX = 1.0/GAM
      XMAX = (REAL(CN)-1.0)/(0.5*GAM**2)
      XUP = AMIN1(XMAX,6.0)
C
      WGASY = GEOMOP+(0.5*GAM**2)*EDGINT(QPAEMI,QPBEMI,XUP,2.,TMAX)
      QPRAPX = 1.0-WGASY
C
      RETURN
      END 
      REAL FUNCTION EDGINT (EDGUP,EDGDN,XUP,T0,TMAX)
C
C     COMPUTES EDGE CORRECTION INTEGRALS TO GEOMETRIC OPTICS
C
      COMPLEX CN
      COMMON /APPXMI/ CN,XX,GAM
      REAL LOLIM,RINT(8)
      EXTERNAL EDGUP,EDGDN
C
C     DELT = MAXIMUM INTEGRATION INCREMENT IN ITERATION TO
C     CALCULATE BELOW-EDGE INTEGRAL
C
      DATA DELT / 4.0 /
C
C     ABOVE-EDGE INTEGRAL (USUALLY LARGER THAN BELOW-EDGE)
C
      LOLIM = 0.0
      UPLIM = XUP
      ACCUR = 1.E-6 
      CALL RKRONQ (LOLIM,UPLIM,EDGUP,ACCUR,RINT,KKK)
      IF (KKK.EQ.0) CALL RINTER (KKK,XX,CN,'EDGINT','EDGUP',LOLIM,UPLIM,
     *   ACCUR,RINT)
      EDGINT = RINT(KKK)
C
C     BELOW-EDGE INTEGRAL
C
      IF (TMAX.LE.T0) THEN
C
         LOLIM = 0.0
         UPLIM = TMAX
         ACCUR = 1.E-6
         CALL RKRONQ (LOLIM,UPLIM,EDGDN,ACCUR,RINT,KKK)
         IF (KKK.EQ.0) CALL RINTER (KKK,XX,CN,'EDGINT','EDGDN',LOLIM, 
     *      UPLIM,ACCUR,RINT) 
         EDGINT = EDGINT+RINT(KKK)
C
      ELSE
C
C     INTEGRATE FROM ZERO TO T0, THEN KEEP ADDING CONTRIBUTIONS TO
C     INTEGRAL UNTIL EITHER (A) THEY BECOME NEGLIGIBLE, OR (B) THE UPPER
C     LIMIT -TMAX- IS REACHED.
C
         LOLIM = 0.0
         UPLIM = T0 
         ACCUR = 1.E-6
         CALL RKRONQ (LOLIM,UPLIM,EDGDN,ACCUR,RINT,KKK)
         IF (KKK.EQ.0) CALL RINTER (KKK,XX,CN,'EDGINT','EDGDN',LOLIM, 
     *      UPLIM,ACCUR,RINT) 
         EDGINT = EDGINT+RINT(KKK)
C
         T = T0
   10    TUP = AMIN1(T+DELT,TMAX)
         LOLIM = T
         UPLIM = TUP
         ACCUR = 1.E-3
         CALL RKRONQ (LOLIM,UPLIM,EDGDN,ACCUR,RINT,KKK)
         IF (KKK.EQ.0) CALL RINTER (KKK,XX,CN,'EDGINT','EDGDN',LOLIM, 
     *      UPLIM,ACCUR,RINT) 
         ADD = RINT(KKK)
         EDGINT = EDGINT+ADD
         T = T+DELT 
         IF (T.LT.TMAX.AND.ABS(ADD/EDGINT).GT.1.E-5) GO TO 10
      ENDIF
C
      RETURN
      END 
      REAL FUNCTION QAGOMI (THETA)
C
C     INTEGRAND IN GEOM. OPTICS PART OF ABSORPTION EFFICIENCY
C
      COMPLEX CN
      COMMON /APPXMI/ CN,XX,GAM
      COMPLEX CI,CSITHP,COSTHP,CR1,CR2,CW
      COMPLEX C
      SQ(C) = REAL(C)**2+AIMAG(C)**2
      DATA CI / (0.,1.) /
C
C
      SINTH = SIN(THETA)
      COSTH = COS(THETA)
      CSITHP = SINTH/CN
      COSTHP = CSQRT(1.0-CSITHP**2)
      CR1 = (COSTH-CN*COSTHP)/(COSTH+CN*COSTHP)
      CR2 = (COSTH*CN-COSTHP)/(COSTH*CN+COSTHP)
C
      IF (AIMAG(CN)*XX.GT.5.0) THEN
C
C        SPECIAL CASE: HI ABSORPTION
C
         QAGOMI = SINTH*COSTH*(2.-SQ(CR1)-SQ(CR2))
C
      ELSE
C
C        GENERAL CASE
C
         COSTSQ = COSTH**2
         CW = CN*COSTHP
         V = AIMAG(CW)
         VT = REAL(CN)*AIMAG(COSTHP)-AIMAG(CN)*REAL(COSTHP) 
         ZETA = 0.5*ALOG(SQ(COSTHP+CI*CSITHP))
         E = EXP(-4.0*XX*(V-ZETA*SINTH))
         R1SQ = SQ(CR1)
         R2SQ = SQ(CR2)
         TERM1 = (1.-R1SQ)-16.*R1SQ*E*(V**2*COSTSQ)/((COSTSQ-SQ(CW))**2+
     *      4.*(V**2*COSTSQ)) 
         TERM2 = (1.-R2SQ)-16.*R2SQ*E*(VT**2*COSTSQ)/((SQ(CN)*COSTSQ-SQ(
     *      COSTHP))**2+4.*(VT**2*COSTSQ))
         QAGOMI = SINTH*COSTH*(1.-E)*(TERM1/(1.-R1SQ*E)+TERM2/(1.-R2SQ*E
     *      ))
      ENDIF
C
      RETURN
      END 
      REAL FUNCTION QAAEMI (X)
C
C     ABOVE-EDGE INTEGRAND FOR ABSORPTION EFFICIENCY
C
      COMPLEX CN
      COMMON /APPXMI/ CN,XX,GAM
      COMPLEX CFITF1
      COMPLEX CI,CSITHP,COSTHP,CW,CRM11,CRM22,CRE11,CRE22,Z,ZC
      COMPLEX C
      SQ(C) = REAL(C)**2+AIMAG(C)**2
      DATA CI / (0.,1.) /
C
C
      SINTH = 1.+0.5*GAM**2*X 
      CSITHP = SINTH/CN
      COSTHP = CSQRT(1.0-CSITHP**2)
      Z = GAM*CFITF1(SQRT(X)) 
      ZC = CONJG(Z) 
      CRM22 = (ZC-CN*COSTHP)/(Z+CN*COSTHP)
      CRE22 = (ZC*CN-COSTHP)/(Z*CN+COSTHP)
C
      IF (AIMAG(CN)*XX.GT.5.0) THEN
C
C        SPECIAL CASE: HI ABSORPTION
C
         QAAEMI = 2.-SQ(CRM22)-SQ(CRE22)
C
      ELSE
C
C        GENERAL CASE
C
         CRM11 = (Z-CN*COSTHP)/(Z+CN*COSTHP)
         CRE11 = (Z*CN-COSTHP)/(Z*CN+COSTHP)
         ZETA = 0.5*ALOG(SQ(COSTHP+CI*CSITHP))
         E = EXP(-4.0*XX*(AIMAG(CN*COSTHP)-ZETA*SINTH))
         QAAEMI = (1.-E)*((1.-SQ(CRM22))/(1.-E*SQ(CRM11))+(1.-SQ(CRE22))
     *      /(1.-E*SQ(CRE11)))
      ENDIF
C
      RETURN
      END 
      REAL FUNCTION QABEMI (T)
C
C     BELOW-EDGE INTEGRAND FOR ABSORPTION EFFICIENCY
C
      COMPLEX CN
      COMMON /APPXMI/ CN,XX,GAM
      COMPLEX CFITF2
      COMPLEX CI,CSITHP,COSTHP,CW,CNSQ,CRM11,CRM22,CRE11,CRE22,CRM,CRE,Z
     *   ,ZC
      COMPLEX C
      SQ(C) = REAL(C)**2+AIMAG(C)**2
      DATA CI / (0.,1.) /
C
C
      X = T**2
      SINTH = 1.-0.5*GAM**2*X 
      CSITHP = SINTH/CN
      COSTHP = CSQRT(1.0-CSITHP**2)
      CW = CN*COSTHP
      Z = GAM*CFITF2(T)
      CRM22 = (CONJG(Z)-CW)/(Z+CW)
      CRE22 = (CN*CONJG(Z)-COSTHP)/((CN*Z)+COSTHP)
      CRM = ((GAM*T)-CW)/((GAM*T)+CW)
      CRE = ((CN*(GAM*T))-COSTHP)/((CN*(GAM*T))+COSTHP)
C
      IF (AIMAG(CN)*XX.GT.5.0) THEN
C
C        SPECIAL CASE: HI ABSORPTION
C
         QABEMI = 2.*T*(SQ(CRM)-SQ(CRM22)+SQ(CRE)-SQ(CRE22))
C
      ELSE
C
C        GENERAL CASE
C
         CRM11 = (Z-CW)/(Z+CW)
         CRE11 = ((CN*Z)-COSTHP)/((CN*Z)+COSTHP)
         RMSQ = SQ(CRM)
         RESQ = SQ(CRE)
         ZETA = 0.5*ALOG(SQ(COSTHP+CI*CSITHP))
         E = EXP(-4.0*XX*(AIMAG(CW)-ZETA*SINTH))
         QABEMI = 2.*T*(1.-E)*((1.-SQ(CRM22))/(1.-E*SQ(CRM11))+(1.-
     *      SQ(CRE22))/(1.-E*SQ(CRE11))-(1.-RMSQ)/(1.-RMSQ*E)-(1.-RESQ)/
     *      (1.-RESQ*E))
      ENDIF
C
      RETURN
      END 
      REAL FUNCTION QPGOMI (THETA)
C
C     INTEGRAND IN THE GEOM. OPTICS EXPRESSION FOR W*G, WHERE W AND G 
C     ARE THE GEOMETRICAL OPTICS CONTRIBUTIONS TO QSCA AND THE ASYMMETRY
C     FACTOR RESPECTIVELY
C
      COMPLEX CN
      COMMON /APPXMI/ CN,XX,GAM
      COMPLEX CI,CSITHP,COSTHP,CR1,CR2,CFAC,CW,CNSQ,EM2ITH
      COMPLEX C
      SQ(C) = REAL(C)**2+AIMAG(C)**2
      DATA CI / (0.,1.) /
C
C
      SINTH = SIN(THETA)
      COSTH = COS(THETA)
      CSITHP = SINTH/CN
      COSTHP = CSQRT(1.0-CSITHP**2)
      CW = CN*COSTHP
      CR1 = (COSTH-CW)/(COSTH+CW)
      CR2 = (COSTH*CN-COSTHP)/(COSTH*CN+COSTHP)
C
      IF (AIMAG(CN)*XX.GT.5.0) THEN
C
C        SPECIAL CASE: HI ABSORPTION
C
         QPGOMI = SINTH*COSTH*(2.*COSTH**2-1.0)*(-SQ(CR1)-SQ(CR2))
C
      ELSE
C
C        GENERAL CASE
C
         EM2ITH = CMPLX(2.*COSTH**2-1.,-2.*SINTH*COSTH)
         R1SQ = SQ(CR1)
         R2SQ = SQ(CR2)
         ZETA = 0.5*ALOG(SQ(COSTHP+CI*CSITHP))
         E = EXP(-4.0*XX*(AIMAG(CW)-ZETA*SINTH))
         CFAC = E*(2.*COSTHP**2-1.+(0.,2.)*CSITHP*COSTHP)
         QPGOMI = SINTH*COSTH*REAL(EM2ITH*(-R1SQ-R2SQ+CFAC*(SQ(1.-CR1**2
     *      )/(1.+CFAC*R1SQ)+SQ(1.-CR2**2)/(1.+CFAC*R2SQ))))
      ENDIF
C
      RETURN
      END 
      REAL FUNCTION QPAEMI (X)
C
C     ABOVE-EDGE INTEGRAND FOR W*G
C
      COMPLEX CN
      COMMON /APPXMI/ CN,XX,GAM
      COMPLEX CFITF1
      COMPLEX CI,CSITHP,COSTHP,CW,CWP,CWM,CRM11,CRM22,CRE11,CRE22,C11M,
     *   C22M,C11E,C22E,CE,CFF,CP,CNSQ,CIMSQ,CC0,Z,ZC
      COMPLEX C
      SQ(C) = REAL(C)**2+AIMAG(C)**2
      DATA CI / (0.,1.) /
C
C
      CNSQ = CN**2
      CIMSQ = (0.,1.)*(CNSQ-1.0)
      CC0 = 2.*CNSQ-1.0
      SINTH = 1.+0.5*GAM**2*X 
      CSITHP = SINTH/CN
      COSTHP = CSQRT(1.0-CSITHP**2)
      CW = CN*COSTHP
      CWP = CW+CIMSQ
      Z = GAM*CFITF1(SQRT(X)) 
      ZC = CONJG(Z) 
      CP = CC0+CIMSQ*CW
C
      IF (AIMAG(CN)*XX.GT.5.0) THEN
C
C        SPECIAL CASE: HI ABSORPTION
C
         QPAEMI = REAL(CONJG((CW-ZC)/(CW+Z))*(CWP-CNSQ*ZC)/(CWP+CNSQ*Z))
     *      +REAL(CONJG((CW-CNSQ*ZC)/(CW+CNSQ*Z))*(CWP-CP*ZC)/(CWP+CP*Z)
     *      )-2.0
      ELSE
C
C        GENERAL CASE
C
         CWM = CW-CIMSQ
         ZETA = 0.5*ALOG(SQ(COSTHP+CI*CSITHP))
         E = EXP(-4.0*XX*(AIMAG(CW)-ZETA*SINTH))
         CE = 2.*COSTHP**2-1.0+(0.,2.)*CSITHP*COSTHP
         CRM11 = -(Z-CW)/(Z+CW)
         CRM22 = (ZC-CW)/(Z+CW)
         C11M = (CNSQ*Z-CWM)/(CNSQ*Z+CWP)
         C22M = (CNSQ*ZC-CWP)/(CNSQ*Z+CWP)
         CRE11 = -(CNSQ*Z-CW)/(CNSQ*Z+CW)
         CRE22 = (CNSQ*ZC-CW)/(CNSQ*Z+CW)
         C11E = ((CC0-CIMSQ*CW)*Z-CWM)/(CP*Z+CWP) 
         C22E = (CP*ZC-CWP)/(CP*Z+CWP)
         CFF = (1.+CI*ZC)/(1.-CI*Z)
         QPAEMI = REAL(C22M*CONJG(CRM22))+REAL(C22E*CONJG(CRE22))-E*
     *      (REAL((CFF+C22M)*(CE+C11M)*CONJG((1.+CRM11)*(1.+CRM22))/(1.+
     *      E*C11M*CONJG(CRM11)))+REAL((CFF+C22E)*(CE+C11E)*CONJG((1.+
     *      CRE11)*(1.+CRE22))/(1.+E*C11E*CONJG(CRE11))))-2.0
      ENDIF
C
      RETURN
      END 
      REAL FUNCTION QPBEMI (T)
C
C     BELOW-EDGE INTEGRAND FOR W*G
C
      COMPLEX CN
      COMMON /APPXMI/ CN,XX,GAM
      REAL PART(2)
      COMPLEX CFITF2
      COMPLEX CI,CSITHP,COSTHP,CW,CWP,CWM,CRM11,CRM22,CRE11,CRE22,C11M,
     *   C22M,C11E,C22E,CE,CM,CP,CFF,CNSQ,CIMSQ,CC0,CQ,CQC,Z,ZC
      COMPLEX C
      SQ(C) = REAL(C)**2+AIMAG(C)**2
      DATA CI / (0.,1.) /
C
C
      X = T**2
      CNSQ = CN**2
      CIMSQ = (0.,1.)*(CNSQ-1.0)
      CC0 = 2.*CNSQ-1.0
      SINTH = 1.-0.5*GAM**2*X 
      CSITHP = SINTH/CN
      COSTHP = CSQRT(1.0-CSITHP**2)
      CW = CN*COSTHP
      CWP = CW+CIMSQ
      CP = CC0+CIMSQ*CW
      Z = GAM*CFITF2(T)
      CQ = GAM*CMPLX(X*T,0.25)
C
      IF (AIMAG(CN)*XX.GT.5.0) THEN
C
C        SPECIAL CASE: HI ABSORPTION
C
         ZC = CONJG(Z)
         CQC = CONJG(CQ)
C
         QPBEMI = REAL(CONJG((CW-ZC)/(CW+Z))*(CWP-CNSQ*ZC)/(CWP+CNSQ*Z))
     *      +REAL(CONJG((CW-CNSQ*ZC)/(CW+CNSQ*Z))*(CWP-CP*ZC)/(CWP+CP*Z)
     *      )-REAL(CONJG((X*CW-CQC)/(X*CW+CQ))*(X*CWP-CNSQ*CQC)/(X*CWP+
     *      CNSQ*CQ))-REAL(CONJG((X*CW-CNSQ*CQC)/(X*CW+CNSQ*CQ))*(X*CWP-
     *      CP*CQC)/(X*CWP+CP*CQ))
C
         QPBEMI = 2.*T*QPBEMI 
C
      ELSE
C
C        GENERAL CASE
C
         CWM = CW-CIMSQ
         CM = CC0-CIMSQ*CW
         ZETA = 0.5*ALOG(SQ(COSTHP+CI*CSITHP))
         E = EXP(-4.0*XX*(AIMAG(CW)-ZETA*SINTH))
         CE = 2.*COSTHP**2-1.0+(0.,2.)*CSITHP*COSTHP
C
         DO 10 II = 1, 2
            ZC = CONJG(Z)
            CRM11 = -(Z-CW)/(Z+CW)
            CRM22 = (ZC-CW)/(Z+CW)
            C11M = (CNSQ*Z-CWM)/(CNSQ*Z+CWP)
            C22M = (CNSQ*ZC-CWP)/(CNSQ*Z+CWP)
            CRE11 = -(CNSQ*Z-CW)/(CNSQ*Z+CW)
            CRE22 = (CNSQ*ZC-CW)/(CNSQ*Z+CW)
            C11E = (CM*Z-CWM)/(CP*Z+CWP)
            C22E = (CP*ZC-CWP)/(CP*Z+CWP)
C
            IF (II.EQ.1) CFF = (1.0+CI*ZC)/(1.0-CI*Z)
            IF (II.EQ.2) CFF = (X+CI*ZC)/(X-CI*Z) 
C
            PART(II) = REAL(C22M*CONJG(CRM22))+REAL(C22E*CONJG(CRE22))-E
     *         *(REAL((CFF+C22M)*(CE+C11M)*CONJG((1.+CRM11)*(1.+CRM22))/
     *         (1.+E*C11M*CONJG(CRM11)))+REAL((CFF+C22E)*(CE+C11E)*CONJG
     *         ((1.+CRE11)*(1.+CRE22))/(1.+E*C11E*CONJG(CRE11))))
C
            IF (II.EQ.1) THEN 
               Z = CQ
               CW = X*CW
               CWM = X*CWM
               CWP = X*CWP
            ENDIF
   10    CONTINUE
C
         QPBEMI = 2.*T*(PART(1)-PART(2))
C
      ENDIF
C
      RETURN
      END 
      COMPLEX FUNCTION CFITF1 (T)
C
C     FITS THE FUNCTION: 
C
C     F1( T ) =  - EXP( I*PI/6 ) * AI-PRIME(Z1) / AI(Z1)
C
C     WHERE AI IS THE AIRY FUNCTION,  AI-PRIME ITS DERIVATIVE,
C     I = SQRT(-1),  AND Z1 = T**2 * EXP( 2*I*PI/3 )
C
      PARAMETER (R3=1./3.,F4D3=4./3.,F5D24=5./24.)
      REAL NUMR,NUMI,PR(9),QR(8),PI(9),QI(8)
      DATA PR / 1.7315719625372,2.9897177730570,2.1675211928955,
     *   1.2069718439574,0.54713883599714,0.18584610613686, 
     *   0.045603149940072,0.0070247460663800,0.00023588144084310 /
      DATA QR / 0.65541334445671,1.0,0.75451174058638,0.36651461456878,
     *   0.15513583081292,0.043966063702531,0.0082020278256018,
     *   0.00041273922145329 /
      DATA PI / 2.0010707952241,3.3622086655917,2.5394407186989,
     *   1.4270851380827,0.6954258833838,0.25268331608842,
     *   0.073924586046895,0.014168632079269,0.0017133588160996 /
      DATA QI / 0.85209525954058,0.95591657885629,1.0,0.37848651011038,
     *   0.22670902999311,0.057589255900429,0.016077200177600,
     *   0.0021324145727573 / 
C
C
      IF (T.LE.3.0) THEN
C
C        RATIONAL FUNCTION APPROXIMATION
C
         ARG = R3*(2.*T-3.)
         TNM1 = 1.0 
         TN = ARG
         NUMR = PR(1)+PR(2)*TN
         NUMI = PI(1)+PI(2)*TN
         DENR = QR(1)+QR(2)*TN
         DENI = QI(1)+QI(2)*TN
         DO 10 K = 3, 8
            TNP1 = (2.*ARG)*TN-TNM1
            TNM1 = TN
            TN = TNP1
            NUMR = NUMR+PR(K)*TN
            NUMI = NUMI+PI(K)*TN
            DENR = DENR+QR(K)*TN
            DENI = DENI+QI(K)*TN
   10    CONTINUE
C
         TN = (2.*ARG)*TN-TNM1
         NUMR = NUMR+PR(9)*TN 
         NUMI = NUMI+PI(9)*TN 
         CFITF1 = CMPLX(EXP(-F4D3*T**3)*NUMR/DENR,NUMI/DENI)
C
      ELSE
C
C        ASYMPTOTIC FORM
C
         TCUB = T**3
         CFITF1 = T*CMPLX(EXP(-F4D3*TCUB)*(1.0-F5D24/TCUB),1.0-0.25/TCUB
     *      )
      ENDIF
C
      RETURN
      END 
      COMPLEX FUNCTION CFITF2 (T)
C
C     FITS THE FUNCTION: 
C
C     F2( T ) =  - EXP( I*PI/6 ) * AI-PRIME(Z2) / AI(Z2)
C
C     WHERE I = SQRT(-1),  AI IS THE AIRY FUNCTION,  AI-PRIME ITS
C     DERIVATIVE, AND Z2 = T**2 * EXP( - I*PI/3 ) 
C
      REAL NUMR,NUMI,PR(7),QR(6),PI(5),QI(7)
      DATA PR / 2.0599303588769,3.3534188921726,1.8877130218248,
     *   0.72699531367472,0.18767024937609,0.029189541046810,
     *   0.0021851157151614 / 
      DATA QR / 0.67306884509460,1.0,0.49837696874507,0.15350461962531,
     *   0.028703419917545,0.0024609936972022 /
      DATA PI / 0.016498766116069,0.021397695518234,0.010819984202854,
     *   0.0027705151330218,0.00046655822980077 / 
      DATA QI / 0.60230756331207,1.0,0.59913150168617,0.25428428666886,
     *   0.075618290222945,0.014257900593919,0.0014081320250867 /
C
C
      IF (T.LE.3.5) THEN
C
C        RATIONAL FUNCTION APPROXIMATION
C
         ARG = (2.*T-3.5)/3.5 
         TNM1 = 1.0 
         TN = ARG
         NUMR = PR(1)+PR(2)*TN
         NUMI = PI(1)+PI(2)*TN
         DENR = QR(1)+QR(2)*TN
         DENI = QI(1)+QI(2)*TN
         DO 10 K = 3, 5
            TNP1 = (2.*ARG)*TN-TNM1
            TNM1 = TN
            TN = TNP1
            NUMR = NUMR+PR(K)*TN
            NUMI = NUMI+PI(K)*TN
            DENR = DENR+QR(K)*TN
            DENI = DENI+QI(K)*TN
   10    CONTINUE
C
         TNP1 = (2.*ARG)*TN-TNM1
         TNM1 = TN
         TN = TNP1
         NUMR = NUMR+PR(6)*TN 
         DENR = DENR+QR(6)*TN 
         DENI = DENI+QI(6)*TN 
         TN = (2.*ARG)*TN-TNM1
         NUMR = NUMR+PR(7)*TN 
         DENI = DENI+QI(7)*TN 
         CFITF2 = CMPLX(NUMR/DENR,NUMI/DENI)
C
      ELSE
C
C        ASYMPTOTIC SERIES
C
         T6 = T**6
         CFITF2 = CMPLX(T*(1.0+0.15625/T6),0.25/T**2*(1.0-0.9375/T6)) 
      ENDIF
C
      RETURN
      END 
      SUBROUTINE RINTER (KKK,XX,CN,CALROU,INTROU,LOLIM,UPLIM,ACCUR,RINT)
C
C     SET  KKK = 8, REPORT ON CONVERGENCE FAILURE IN KRONROD SCHEME
C
      INTEGER KKK
      REAL XX,LOLIM,UPLIM,ACCUR,RINT(*) 
      COMPLEX CN
      CHARACTER*(*) CALROU,INTROU
C
C
      KKK = 8
      WRITE (*,900) CALROU,INTROU,XX,CN,ACCUR,LOLIM,UPLIM,(RINT(I),I=1,8
     *   )
      WRITE (*,901) 
      RETURN
C
C
C
C
C
C
C
C
  900 FORMAT (/,' >>>>>>>> WARNING:  KRONROD CONVERGENCE FAILURE.  ',
     *   ' CALLED BY SUBPROGRAM  ',A,/,11X,'FOR INTEGRAND  ',A,
     *   ',   SIZE PARAM =',1P,E11.3,'    REFRAC INDEX =',2E11.3,/,11X,
     *   'ACCURACY DEMANDED =',1P,E10.2,'    INTEG. LIMITS:  ',E11.3,
     *   '  TO ',E11.3,/,11X,'APPROX. INTEGRAL VALUES :',E17.8,/,(37X,
     *   E17.8))
  901 FORMAT (1X,78('='))
C
      END 
      SUBROUTINE RKRONQ (A,B,RINTGD,EPSIL,RINT,KCONVG)
C
C     >>>>>>>>   REAL INTEGRAND   <<<<<<<
C
C     NUMERICALLY INTEGRATES OVER AN INTERVAL USING A SEQUENCE OF
C     INTERLEAVING 1, 3, 7, 15, 31, 63, 127 AND 255-POINT EXTENDED
C     GAUSS-TYPE QUADRATURE FORMULAE.  SINCE EACH SUCCESSIVE
C     FORMULA EMPLOYS ALL POINTS USED BY ITS PREDECESSOR, NO
C     INTEGRAND VALUES ARE WASTED WHEN THE ORDER OF THE INTEGRATION
C     FORMULA IS INCREASED.
C
C     SPECIFICATIONS OF ARGUMENTS
C
      INTEGER KCONVG
      REAL A,B,EPSIL,RINTGD,RINT(*)
C
C     SPECIFICATIONS OF LOCAL VARIABLES 
C
      PARAMETER (MAXNO=8)
      PARAMETER (MAXORD=2**MAXNO-1,MXORD2=(MAXORD+1)/2)
      LOGICAL PASS1 
      INTEGER NORD(MAXNO)
      REAL ABK(MXORD2),WTK(MXORD2,MAXNO),X(MAXORD),W(MAXORD),
     *   RFUNCT(MAXORD)
      COMMON /KRABWT/ ABK,WTK 
      DATA NORD / 1,3,7,15,31,63,127,255 /
      DATA PASS1 / .TRUE. /
C
C
      IF (PASS1) THEN
C
C        ENSURE -ABK,WTK- DATA IS LOADED AT RUN TIME
C
         CALL KRNDAT
         PASS1 = .FALSE.
      ENDIF
C
      IF (B.EQ.A) THEN
         KCONVG = 1 
         RINT(1) = 0.0
      ENDIF
C
      CP = 0.5*(B+A)
      CM = 0.5*(B-A)
C
C     SET UP QUADRATURE ABSCISSAE
C
      X(1) = CP+CM*ABK(1)
      DO 10 I = 2, MXORD2
         X(2*I-2) = CP+CM*ABK(I)
         X(2*I-1) = CP-CM*ABK(I)
   10 CONTINUE
C
C     FIRST QUADRATURE RULE IN SEQUENCE (MID-POINT RULE)
C
      RFUNCT(1) = RINTGD(X(1))
      RINT(1) = CM*WTK(1,1)*RFUNCT(1)
C
C     LOOP OVER QUADRATURE RULES
C
      DO 50 K = 2, MAXNO
C
C        SET UP QUADRATURE WEIGHTS
C
         W(1) = CM*WTK(1,K)
         DO 20 I = 2, (NORD(K)+1)/2
            W(2*I-2) = CM*WTK(I,K)
            W(2*I-1) = W(2*I-2)
   20    CONTINUE
C
         RINT(K) = 0.0
C
C        CONTRIBUTION FROM FUNCTION VALUES ALREADY COMPUTED 
C
         DO 30 I = 1, NORD(K-1)
            RINT(K) = RINT(K)+W(I)*RFUNCT(I)
   30    CONTINUE
C
C        CONTRIBUTION FROM NEW FUNCTION VALUES
C
         DO 40 I = NORD(K-1)+1, NORD(K) 
            RFUNCT(I) = RINTGD(X(I))
            RINT(K) = RINT(K)+W(I)*RFUNCT(I)
   40    CONTINUE
C
C        CHECK FOR CONVERGENCE
C
         IF (ABS(RINT(K)-RINT(K-1)).LE.ABS(EPSIL*RINT(K))) GO TO 60
   50 CONTINUE
C
C     REQUESTED CONVERGENCE NOT ACHIEVED
C
      K = 0
C
   60 KCONVG = K
C
      RETURN
      END 

C
C
      SUBROUTINE KRNDAT
C
C     ABSCISSAE (ABK) AND WEIGHTS (WTK) FOR KRONROD QUADRATURE
C
      PARAMETER (MAXNO=8)
      PARAMETER (MAXORD=2**MAXNO-1,MXORD2=(MAXORD+1)/2)
      REAL ABK(MXORD2),WTK(MXORD2,MAXNO)
      COMMON /KRABWT/ ABK,WTK 
C
      DATA (ABK(I),I=1,58) / 0.0,0.77459666924148E+00,0.96049126870802E+
     *   00,0.43424374934680E+00,0.99383196321276E+00,0.88845923287226E+
     *   00,0.62110294673723E+00,0.22338668642897E+00,0.99909812496767E+
     *   00,0.98153114955374E+00,0.92965485742974E+00,0.83672593816887E+
     *   00,0.70249620649153E+00,0.53131974364437E+00,0.33113539325798E+
     *   00,0.11248894313319E+00,0.99987288812036E+00,0.99720625937222E+
     *   00,0.98868475754743E+00,0.97218287474858E+00,0.94634285837340E+
     *   00,0.91037115695701E+00,0.86390793819369E+00,0.80694053195022E+
     *   00,0.73975604435270E+00,0.66290966002478E+00,0.57719571005205E+
     *   00,0.48361802694584E+00,0.38335932419873E+00,0.27774982202182E+
     *   00,0.16823525155221E+00,0.56344313046593E-01,0.99998243035489E+
     *   00,0.99959879967191E+00,0.99831663531841E+00,0.99572410469841E+
     *   00,0.99149572117810E+00,0.98537149959852E+00,0.97714151463970E+
     *   00,0.96663785155842E+00,0.95373000642576E+00,0.93832039777959E+
     *   00,0.92034002547001E+00,0.89974489977694E+00,0.87651341448471E+
     *   00,0.85064449476835E+00,0.82215625436498E+00,0.79108493379985E+
     *   00,0.75748396638051E+00,0.72142308537010E+00,0.68298743109108E+
     *   00,0.64227664250976E+00,0.59940393024224E+00,0.55449513263193E+
     *   00,0.50768775753372E+00,0.45913001198983E+00,0.40897982122989E+
     *   00,0.35740383783153E+00 /
      DATA (ABK(I),I=59,115) / 0.30457644155671E+00,0.25067873030348E+00
     *   ,0.19589750271110E+00,0.14042423315256E+00,0.84454040083711E-01
     *   ,0.28184648949746E-01,0.99999759637975E+00,0.99994399620705E+00
     *   ,0.99976049092443E+00,0.99938033802502E+00,0.99874561446810E+00
     *   ,0.99780535449596E+00,0.99651414591489E+00,0.99483150280062E+00
     *   ,0.99272134428279E+00,0.99015137040077E+00,0.98709252795403E+00
     *   ,0.98351865757863E+00,0.97940628167086E+00,0.97473445975240E+00
     *   ,0.96948465950246E+00,0.96364062156981E+00,0.95718821610986E+00
     *   ,0.95011529752129E+00,0.94241156519108E+00,0.93406843615773E+00
     *   ,0.92507893290708E+00,0.91543758715576E+00,0.90514035881326E+00
     *   ,0.89418456833556E+00,0.88256884024734E+00,0.87029305554811E+00
     *   ,0.85735831088623E+00,0.84376688267271E+00,0.82952219463740E+00
     *   ,0.81462878765514E+00,0.79909229096084E+00,0.78291939411828E+00
     *   ,0.76611781930376E+00,0.74869629361694E+00,0.73066452124218E+00
     *   ,0.71203315536225E+00,0.69281376977911E+00,0.67301883023042E+00
     *   ,0.65266166541002E+00,0.63175643771119E+00,0.61031811371519E+00
     *   ,0.58836243444766E+00,0.56590588542365E+00,0.54296566649831E+00
     *   ,0.51955966153746E+00,0.49570640791876E+00,0.47142506587166E+00
     *   ,0.44673538766203E+00,0.42165768662616E+00,0.39621280605762E+00
     *   ,0.37042208795008E+00 /
      DATA (ABK(I),I=116,128) / 0.34430734159944E+00,0.31789081206848E+
     *   00,0.29119514851825E+00,0.26424337241093E+00,0.23705884558983E+
     *   00,0.20966523824318E+00,0.18208649675925E+00,0.15434681148138E+
     *   00,0.12647058437230E+00,0.98482396598119E-01,0.70406976042855E-
     *   01,0.42269164765364E-01,0.14093886410783E-01 /
C
      DATA WTK(1,1) / 2.0 /
      DATA (WTK(I,2),I=1,2) / 0.88888888888889E+00,0.55555555555556E+00
     *   /
      DATA (WTK(I,3),I=1,4) / 0.45091653865847E+00,0.26848808986833E+00,
     *   0.10465622602647E+00,0.40139741477596E+00 /
      DATA (WTK(I,4),I=1,8) / 0.22551049979821E+00,0.13441525524378E+00,
     *   0.51603282997080E-01,0.20062852937699E+00,0.17001719629940E-01,
     *   0.92927195315125E-01,0.17151190913639E+00,0.21915685840159E+00
     *   /
      DATA (WTK(I,5),I=1,16) / 0.11275525672077E+00,0.67207754295991E-01
     *   ,0.25807598096177E-01,0.10031427861179E+00,0.84345657393211E-02
     *   ,0.46462893261758E-01,0.85755920049990E-01,0.10957842105593E+00
     *   ,0.25447807915619E-02,0.16446049854388E-01,0.35957103307129E-01
     *   ,0.56979509494123E-01,0.76879620499004E-01,0.93627109981265E-01
     *   ,0.10566989358023E+00,0.11195687302095E+00 /
      DATA (WTK(I,6),I=1,32) / 0.56377628360385E-01,0.33603877148208E-01
     *   ,0.12903800100351E-01,0.50157139305899E-01,0.42176304415588E-02
     *   ,0.23231446639910E-01,0.42877960025008E-01,0.54789210527963E-01
     *   ,0.12651565562301E-02,0.82230079572359E-02,0.17978551568128E-01
     *   ,0.28489754745834E-01,0.38439810249456E-01,0.46813554990628E-01
     *   ,0.52834946790117E-01,0.55978436510476E-01,0.36322148184553E-03
     *   ,0.25790497946857E-02,0.61155068221173E-02,0.10498246909621E-01
     *   ,0.15406750466559E-01,0.20594233915913E-01,0.25869679327215E-01
     *   ,0.31073551111688E-01,0.36064432780783E-01,0.40715510116944E-01
     *   ,0.44914531653632E-01,0.48564330406673E-01,0.51583253952048E-01
     *   ,0.53905499335266E-01,0.55481404356559E-01,0.56277699831254E-01
     *    /
      DATA (WTK(I,7),I=1,56) / 0.28188814180192E-01,0.16801938574104E-01
     *   ,0.64519000501757E-02,0.25078569652950E-01,0.21088152457266E-02
     *   ,0.11615723319955E-01,0.21438980012504E-01,0.27394605263981E-01
     *   ,0.63260731936263E-03,0.41115039786547E-02,0.89892757840641E-02
     *   ,0.14244877372917E-01,0.19219905124728E-01,0.23406777495314E-01
     *   ,0.26417473395058E-01,0.27989218255238E-01,0.18073956444539E-03
     *   ,0.12895240826104E-02,0.30577534101755E-02,0.52491234548089E-02
     *   ,0.77033752332797E-02,0.10297116957956E-01,0.12934839663607E-01
     *   ,0.15536775555844E-01,0.18032216390391E-01,0.20357755058472E-01
     *   ,0.22457265826816E-01,0.24282165203337E-01,0.25791626976024E-01
     *   ,0.26952749667633E-01,0.27740702178280E-01,0.28138849915627E-01
     *   ,0.50536095207863E-04,0.37774664632698E-03,0.93836984854238E-03
     *   ,0.16811428654215E-02,0.25687649437940E-02,0.35728927835173E-02
     *   ,0.46710503721143E-02,0.58434498758356E-02,0.70724899954336E-02
     *   ,0.83428387539682E-02,0.96411777297025E-02,0.10955733387838E-01
     *   ,0.12275830560083E-01,0.13591571009765E-01,0.14893641664815E-01
     *   ,0.16173218729578E-01,0.17421930159464E-01,0.18631848256139E-01
     *   ,0.19795495048098E-01,0.20905851445812E-01,0.21956366305318E-01
     *   ,0.22940964229388E-01,0.23854052106039E-01,0.24690524744488E-01
     *    /
      DATA (WTK(I,7),I=57,64) / 0.25445769965465E-01,0.26115673376706E-
     *   01,0.26696622927450E-01,0.27185513229625E-01,0.27579749566482E-
     *   01,0.27877251476614E-01,0.28076455793817E-01,0.28176319033017E-
     *   01 /
      DATA (WTK(I,8),I=1,58) / 0.14094407090096E-01,0.84009692870519E-02
     *   ,0.32259500250879E-02,0.12539284826475E-01,0.10544076228633E-02
     *   ,0.58078616599776E-02,0.10719490006252E-01,0.13697302631991E-01
     *   ,0.31630366082226E-03,0.20557519893274E-02,0.44946378920321E-02
     *   ,0.71224386864584E-02,0.96099525623639E-02,0.11703388747657E-01
     *   ,0.13208736697529E-01,0.13994609127619E-01,0.90372734658751E-04
     *   ,0.64476204130573E-03,0.15288767050878E-02,0.26245617274044E-02
     *   ,0.38516876166399E-02,0.51485584789782E-02,0.64674198318037E-02
     *   ,0.77683877779220E-02,0.90161081951957E-02,0.10178877529236E-01
     *   ,0.11228632913408E-01,0.12141082601668E-01,0.12895813488012E-01
     *   ,0.13476374833816E-01,0.13870351089140E-01,0.14069424957813E-01
     *   ,0.25157870384281E-04,0.18887326450651E-03,0.46918492424785E-03
     *   ,0.84057143271072E-03,0.12843824718970E-02,0.17864463917587E-02
     *   ,0.23355251860572E-02,0.29217249379178E-02,0.35362449977168E-02
     *   ,0.41714193769841E-02,0.48205888648513E-02,0.54778666939189E-02
     *   ,0.61379152800414E-02,0.67957855048828E-02,0.74468208324076E-02
     *   ,0.80866093647888E-02,0.87109650797321E-02,0.93159241280694E-02
     *   ,0.98977475240488E-02,0.10452925722906E-01,0.10978183152659E-01
     *   ,0.11470482114694E-01,0.11927026053019E-01,0.12345262372244E-01
     *   ,0.12722884982732E-01,0.13057836688353E-01 /
      DATA (WTK(I,8),I=59,115) / 0.13348311463725E-01,0.13592756614812E-
     *   01,0.13789874783241E-01,0.13938625738307E-01,0.14038227896909E-
     *   01,0.14088159516508E-01,0.69379364324108E-05,0.53275293669781E-
     *   04,0.13575491094923E-03,0.24921240048300E-03,0.38974528447328E-
     *   03,0.55429531493037E-03,0.74028280424450E-03,0.94536151685853E-
     *   03,0.11674841174300E-02,0.14049079956552E-02,0.16561127281544E-
     *   02,0.19197129710139E-02,0.21944069253638E-02,0.24789582266576E-
     *   02,0.27721957645934E-02,0.30730184347026E-02,0.33803979910869E-
     *   02,0.36933779170256E-02,0.40110687240750E-02,0.43326409680930E-
     *   02,0.46573172997569E-02,0.49843645647655E-02,0.53130866051871E-
     *   02,0.56428181013845E-02,0.59729195655082E-02,0.63027734490858E-
     *   02,0.66317812429019E-02,0.69593614093904E-02,0.72849479805538E-
     *   02,0.76079896657190E-02,0.79279493342949E-02,0.82443037630329E-
     *   02,0.85565435613077E-02,0.88641732094825E-02,0.91667111635608E-
     *   02,0.94636899938301E-02,0.97546565363174E-02,0.10039172044057E-
     *   01,0.10316812330948E-01,0.10587167904885E-01,0.10849844089337E-
     *   01,0.11104461134007E-01,0.11350654315981E-01,0.11588074033044E-
     *   01,0.11816385890830E-01,0.12035270785280E-01,0.12244424981612E-
     *   01,0.12443560190714E-01,0.12632403643542E-01,0.12810698163877E-
     *   01,0.12978202239537E-01 /
      DATA (WTK(I,8),I=116,128) / 0.13134690091960E-01,0.13279951743930E
     *   -01,0.13413793085110E-01,0.13536035934956E-01,0.13646518102571E
     *   -01,0.13745093443002E-01,0.13831631909506E-01,0.13906019601325E
     *   -01,0.13968158806517E-01,0.14017968039457E-01,0.14055382072650E
     *   -01,0.14080351962554E-01,0.14092845069160E-01 /
C
      RETURN
      END 


