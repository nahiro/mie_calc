      SUBROUTINE SIZEDISDN(KN,IA,ID,NSD,NMD,CM,SM,RMM,
     &RMIN,RMAX,RRR,AR,AC)
C*********************************************************************
C**     Determining bi-modal LogNormal size distribution:           **
C**         d(...)/dlnR in KN - grid radius points                  **
C*********************************************************************
C  INPUT:
C************
C  KN   I(NSD) - number of radius points
C        <0 logarithmic intervals
C        >0 linear intervals
C  IA   I  - defines quadrature 
C     =0 - rectangualar approximation
C     =1 - trapezoidal approximation
C  ID  I - dimension of d(...)/dlnR or d(...)/dR
C         = 0 - number
C         = 1 - radius
C         = 2 - area
C         = 3 - volume
C  NSD I    - number of sizedistributions (up to 3)
C  NMD I    - number of modes (up to 2)
C  CM  R(NMD,NSD) - concentrations
C  SM  R(NMD,NSD) - halfwidths 
C  RM  R(NMD,NSD) - mean radii (mkm)
C  RMIN  R(NSD) - minimum radius (mkm)
C  RMAX  R(NSD) - maximum radius (mkm)
C*****************************************************
C  OUTPUT:
C  RRR    R(NSD,KPAR) - Radii for Size Dsitribution
C  AR     R(KPAR) - d()/dlnR  or d()/dR (in M)
C  AC     R    - total concentration (M3/M3)
C*****************************************************
      include 'optchar.par'
      PARAMETER (KPAR=KNpar)
      INTEGER KN(NSD),KNN(NSD)
      REAL CM(NMD,NSD),SM(NMD,NSD),RMM(NMD,NSD),RMIN(NSD),RMAX(NSD)
      REAL AR(NSD,KPAR)
cl      REAL CM0(NMD,NSD),SM0(NMD,NSD),RMM0(NMD,NSD)
      REAL CM0(NMD,NSD),RMM0(NMD,NSD)
      REAL CM1(NMD,NSD),RMM1(NMD,NSD)
      REAL CM2(NMD,NSD),RMM2(NMD,NSD)
      REAL CM3(NMD,NSD),RMM3(NMD,NSD)
	REAL CMR(NMD,NSD,0:3),RMMR(NMD,NSD,0:3)
     &    ,AAR(NSD,0:3),ACR(NSD,0:3),ARR(NSD,KPAR,0:3)
      REAL AA(NSD),AC(NSD),RRR(NSD,KPAR)
C*****************************************************
C  recalculation Size distribution d(...)/dlnR 
C  parameters in different dimensions:
C*****************************************************
 
      PI= ACOS( -1.0 )
      DO ISD=1,NSD
      KNN(ISD)=KN(ISD)
      IF(KN(ISD).LT.0) KNN(ISD)=-KN(ISD)
      DO JJ=1,KNN(ISD)
      AR(ISD,JJ)=0.
	DO I1=0,3
      ARR(ISD,JJ,I1)=0.
	ENDDO
      ENDDO
      ENDDO
      DO IM=1,NMD
       DO I=1,NSD
        SM(IM,I)=EXP(SM(IM,I))
       ENDDO
      ENDDO
        IF (ID.EQ.0) THEN
        DO IM=1,NMD
        DO I=1,NSD
        CM0(IM,I)=CM(IM,I)
        RMM0(IM,I)=RMM(IM,I)
        ENDDO
        ENDDO
        GO TO 5
        ENDIF
      IF (ID.EQ.1) THEN
      DO IM=1,NMD
      DO I=1,NSD
      RMM0(IM,I)=EXP(LOG(RMM(IM,I))-LOG(SM(IM,I))*LOG(SM(IM,I)))
      CM0(IM,I)=CM(IM,I)/RMM0(IM,I)*EXP(0.5*LOG(SM(IM,I))*LOG(SM(IM,I)))
      ENDDO
      ENDDO
      GO TO 5
      ENDIF
        IF (ID.EQ.2) THEN
        DO IM=1,NMD
        DO I=1,NSD
        RMM0(IM,I)=EXP(LOG(RMM(IM,I))-2.0*LOG(SM(IM,I))*LOG(SM(IM,I)))
        CM0(IM,I)=CM(IM,I)/(PI*RMM0(IM,I)*RMM0(IM,I))/
     &   EXP(2.0*LOG(SM(IM,I))*LOG(SM(IM,I)))
        ENDDO
        ENDDO
        GO TO 5
        ENDIF
      IF (ID.EQ.3) THEN
      DO IM=1,NMD
      DO I=1,NSD
      RMM0(IM,I)=EXP(LOG(RMM(IM,I))-3.0*LOG(SM(IM,I))*LOG(SM(IM,I)))
      CM0(IM,I)=CM(IM,I)/(4.0/3.0)/
     &(PI*RMM0(IM,I)*RMM0(IM,I)*RMM0(IM,I))
     &/EXP(9.0/2.0*LOG(SM(IM,I))*LOG(SM(IM,I)))
      ENDDO
      ENDDO
      GO TO 5
      ENDIF
    5 CONTINUE
      DO IM=1,NMD
      DO I=1,NSD
	  RMMR(IM,I,0)=RMM0(IM,I)
        CMR(IM,I,0) =CM0(IM,I)
	ENDDO
	ENDDO

      I=0
      DO 4 IM=1,NMD
      DO 4 I=1,NSD
      CM1(IM,I)=CM0(IM,I)*RMM0(IM,I)
     &*EXP(0.5*LOG(SM(IM,I))*LOG(SM(IM,I)))
      RMM1(IM,I)=EXP(LOG(RMM0(IM,I))+LOG(SM(IM,I))*LOG(SM(IM,I)))
      CM2(IM,I)=CM0(IM,I)*PI*RMM0(IM,I)*RMM0(IM,I)*EXP(2.0*LOG(SM(IM,I))
     &*LOG(SM(IM,I)))
      RMM2(IM,I)=EXP(LOG(RMM0(IM,I))+2.0*LOG(SM(IM,I))*LOG(SM(IM,I)))
      CM3(IM,I)=CM0(IM,I)*4.0/3.0*PI*RMM0(IM,I)*RMM0(IM,I)*RMM0(IM,I)
     &*EXP(9.0/2.0*LOG(SM(IM,I))*LOG(SM(IM,I)))
      RMM3(IM,I)=EXP(LOG(RMM0(IM,I))+3.0*LOG(SM(IM,I))*LOG(SM(IM,I))) 
        RMMR(IM,I,1)=RMM1(IM,I)      
        CMR(IM,I,1)=CM1(IM,I)
	  RMMR(IM,I,2)=RMM2(IM,I) 
        CMR(IM,I,2)=CM2(IM,I)     
	  RMMR(IM,I,3)=RMM3(IM,I)  
        CMR(IM,I,3)=CM3(IM,I)
    4 CONTINUE
      OPEN (7,FILE='LNPAR.dat',status='unknown')
      DO I=1,NSD 
      WRITE(7,*)'  halfwidth:'
      WRITE(7,17) (LOG(SM(IM,I)),IM=1,NMD)
      WRITE(7,*) ' number concentration, mean radius:'
      WRITE(7,17) (CM0(IM,I),RMM0(IM,I),IM=1,NMD)
      WRITE(7,*) ' radius concentration, mean radius:'
      WRITE(7,17) (CM1(IM,I),RMM1(IM,I),IM=1,NMD)
      WRITE(7,*) ' area concentration, mean radius:'
      WRITE(7,17) (CM2(IM,I),RMM2(IM,I),IM=1,NMD)
      WRITE(7,*) ' volume concentration, mean radius:'
      WRITE(7,17) (CM3(IM,I),RMM3(IM,I),IM=1,NMD)
      ENDDO
17    FORMAT (4E15.5) 
C******************************************
       DO I1=0,3
      DO II=1,NSD
      AAR(II,I1)=0.0
      ACR(II,I1)=0.0
      ENDDO
cl      DO IM=1,NMD  ! delete after asking Oleg
cl      AA(IM,)=0.0
cl      AC(IM)=0.0
cl      ENDDO
      DO I=1,NSD
      DO IM=1,NMD
      AI=0.
      CALL SDNORM(CMR(IM,I,I1),SM(IM,I),RMMR(IM,I,I1),
     &                                  RMIN(I),RMAX(I),AI)
      AAR(I,I1)=AAR(I,I1)+AI
      ACR(I,I1)=ACR(I,I1)+CMR(IM,I,I1)
      ENDDO
      ENDDO
	 ENDDO ! I1
C***** KN<0 logarithmic intervals *******
      OPEN (9,FILE='SizeDis.dat',status='unknown')
      DO 1 II=1,NSD
      DO I=1,KNN(II)
      AR(II,I)=0.0
      ENDDO
      WRITE(9,*) II,' number of size distributions'
      IF(IA.EQ.1) THEN
      IF(KN(II).LT.0) RH=(LOG(RMAX(II))-LOG(RMIN(II)))/(KNN(II)-1)
      IF(KN(II).GT.0) RH=((RMAX(II))-(RMIN(II)))/(KNN(II)-1)
      ENDIF
      IF(IA.EQ.0) THEN
      IF(KN(II).LT.0) RH=(LOG(RMAX(II))-LOG(RMIN(II)))/(KNN(II)-1)
      IF(KN(II).GT.0) RH=((RMAX(II))-(RMIN(II)))/(KNN(II)-1)
      IF(KN(II).LT.0) RH1=(LOG(RMAX(II))-LOG(RMIN(II)))/(KNN(II)-1)
      IF(KN(II).GT.0) RH1=((RMAX(II))-(RMIN(II)))/(KNN(II)-1)
      ENDIF
      DO 1 IM=1,NMD
      DO 1 I=1,KNN(II)
      IF(IA.EQ.1) THEN
      IF(KN(II).LT.0) RR=EXP(LOG(RMIN(II))+(I-1)*RH)
      IF(KN(II).GT.0) RR=RMIN(II)+(I-1)*RH
      ENDIF
      IF(IA.EQ.0) THEN
CD      IF(I.EQ.1.OR.I.EQ.KNN(II)) THEN
CD      RH=RH1/2
CD      ELSE
CD      RH=RH1
CD      ENDIF
CD      IF(KN(II).LT.0) RR=EXP(LOG(RMIN(II))+I*RH-0.5*RH)
CD      IF(KN(II).GT.0) RR=RMIN(II)+I*RH-0.5*RH
      IF(KN(II).LT.0) RR=EXP(LOG(RMIN(II))+(I-1)*RH)
      IF(KN(II).GT.0) RR=RMIN(II)+(I-1)*RH
      ENDIF
      RRR(II,I)=RR      
	 DO I1=0,3
       ARR(II,I,I1)=ARR(II,I,I1)
     &           +SLOG(CMR(IM,II,I1),SM(IM,II),RMMR(IM,II,I1),RR)
       ENDDO ! I1
       AR(II,I)=AR(II,I)
     &           +SLOG(CM(IM,II),SM(IM,II),RMM(IM,II),RR)
1     CONTINUE 
cl      OPEN (9,FILE='SizeDis.dat',status='unknown')
c **  Write into 'SizeDis.dat'
      DO ISD=1,NSD
       DO I=1,KNN(ISD)
       WRITE(9,*) RRR(ISD,I),AR(ISD,I)
       ENDDO
      ENDDO ! ISD
      CLOSE (9)
c **  Write into 'LNPAR.dat'
      DO I1=0,3
      WRITE(7,*)
      IF (I1.EQ.0) WRITE(7,*) ' ******* number distribution:'
      IF (I1.EQ.1) WRITE(7,*) ' ******* radius distribution:'
	IF (I1.EQ.2) WRITE(7,*) ' ******* area   distribution:'
	IF (I1.EQ.3) WRITE(7,*) ' ******* volume distribution:'

      DO ISD=1,NSD
      IF (ISD.EQ.1) WRITE(7,*) 'first componet:'
      IF (ISD.EQ.2) WRITE(7,*) 'second componet:'
       DO I=1,KNN(ISD)
       WRITE(7,*) RRR(ISD,I),ARR(ISD,I,I1)
       ENDDO
       WRITE(7,*) (ACR(ISD,I1)-AAR(ISD,I1))/ACR(ISD,I1), 
     &'relative error in concentation due to limit: Rmin-Rmax'
      ENDDO ! ISD
      ENDDO ! I1                
      CLOSE (7)
      DO ISD=1,NSD
	 AC(ISD)=ACR(ISD,3)
       DO I=1,KNN(ISD)
       AR(ISD,I)=ARR(ISD,I,3)
       ENDDO
      ENDDO ! ISD

      DO IM=1,NMD
       DO I=1,NSD
        SM(IM,I)=LOG(SM(IM,I))
       ENDDO
      ENDDO

      RETURN
      END

      SUBROUTINE SDNORM(C,S,RM,RMIN,RMAX,AC)
C******************************************************************
C**  Normalization of lognormal function d(...)/dlnR             ** 
C**  for descrete interval of sizes: RMIN-RMAX                   **
C******************************************************************
C INPUT:
C**********
C  C    R  - concentration
C  S    R  - halfwidth
C  RM   R  - mean radius
C  RMIN R  - minimum radius
C  RMAX R  - maximum radius
C******************************************************************* 
C OUTPUT:
C**********
C  AC   R  - normalization constant
C*******************************************************************
      REAL SLOG
      KNSIM=301
      PI = ACOS( -1.0 )
      AHR=(LOG(RMAX)-LOG(RMIN))/(KNSIM-1)
      AA=0.
      DO I=1,KNSIM
      RR=EXP(LOG(RMIN)+(I-1)*AHR)
      ADD=SLOG(C,S,RM,RR)
      AA=AA+ADD*AHR
      ENDDO
      AC=AA
      RETURN
      END
 
cl      REAL FUNCTION SLOG(C,S,RM,RR)
clC******************************************************************
clC**     Lognormal function d(RR)/dlnR                            ** 
clC******************************************************************     
clC  C   R  - concentration
clC  S   R  - halfwidth
clC  RM  R  - mean radius
clC  RR  R  - value
clC****************************************************************** 
cl      PI = ACOS( -1.0 )
cl      IF(S.LE.1) WRITE(*,*) S,'TROUBLES in SLOG:S.LE.1'
cl      SLOG=C/SQRT(2.0*PI)/LOG(S)*EXP((-0.5)*((DLOG(RR/RM)*
cl     & DLOG(RR/RM))/(DLOG(S)*DLOG(S))))
cl      RETURN
cl      END

      REAL FUNCTION SLOG(C,S,RM,RR)
C******************************************************************
C**     Lognormal function d(RR)/dlnR                            ** 
C******************************************************************     
C  C   R  - concentration
C  S   R  - halfwidth
C  RM  R  - mean radius
C  RR  R  - value
C****************************************************************** 
      PI = ACOS( -1.0 )
      IF(S.LE.1) THEN
      WRITE(*,*) S,'TROUBLES in SLOG:S.LE.1'
	STOP 'STOP: in SLOG:S.LE.1'
	ENDIF
	A1=LOG(RR/RM)
	A2=LOG(S)
      SLOG=C/SQRT(2.0*PI)/A2*EXP(-0.5*((A1*A1)/(A2*A2)))
      RETURN
      END
 