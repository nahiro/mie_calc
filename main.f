      PROGRAM DBASE
      PARAMETER (N_WL_OUT=1000, NANG=1000)
      INTEGER KPROB,JPOLY,ELLE,nscale,numaeros
      COMPLEX PM,PMCONJ,A(13000)
      REAL*8 R_IND_I, R_IND_R
      REAL      P(110,4),R(110),F(NANG),H(NANG),XX(20),
     *          WAVE_NEW(N_WL_OUT)
      REAL      ab(900),size_d(900),weig(900)
      REAL      THET(NANG),ENTAVG(NANG),RENT1(NANG),RENT2(NANG),
     *   ENT3(NANG),ENT4(NANG),PHAS1(NANG),PHAS2(NANG),PHAS3(NANG),
     *   PHAS4(NANG),PHAS5(NANG),POL(NANG),APHAS(NANG),RI(800),
     *   ENR(800),ELL(NANG),ANGOLO(NANG),FUNZ_FASE(NANG,10),
     *   XI(0:300,N_WL_OUT)
C ********* MODIFIED by Chiara Levoni 02-04-96 ****************
      REAL      EXT(N_WL_OUT),SCA(N_WL_OUT),ASY(N_WL_OUT),
     *          SS_ALB(N_WL_OUT),ASYM(N_WL_OUT), SCALE(9),
     *          EXT_COEFF(N_WL_OUT,4),SING_S_ALB(N_WL_OUT,4),
     *          ASYM_FAC(N_WL_OUT,4),PHAS_F(NANG,N_WL_OUT,4),
     *          WAVE(N_WL_OUT), SCAT_COEFF(N_WL_OUT,4),
     *          WH(998),PESI(1000),absci(998),COORD(998),C(4),
     *          AMBDA(N_WL_OUT),PARTRE(N_WL_OUT,4),PARTIM(N_WL_OUT,4)
C **************************************************************
      COMMON PM,NTHETA,SINS(1000),COS2(1000),COSA(1000),APM,A
      COMMON /FDIFASE/ PHAS_FUN(1000,N_WL_OUT)
      COMMON /A1R0/ NRS,PII,RMIN,RMAX,SQRPI2,NSTEP
      DATA PII,PI2,SQRPI2,RAD,FACTDB / 3.14159265,6.28318531,5.771724899
     *   ,0.017453292,4.3429448 /
      DATA SCALE /0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0 /
      DATA MAXDIM,NPRNT / 13000,0 /
C
      INCLUDE 'AEROS_DB.PARS'
      INCLUDE 'AEROS_DB.VARS'
C     *******************************************************************
C     /
C     /           This program...........
C     /
C     /
C     /
      character*128 release_str
      character*70 header

C     Name or Version of Aerosol DB
      release_str=' user-ok (11-12-02)'

      OPEN (UNIT=1,FILE='entry.dat',FORM='FORMATTED',STATUS='old')
      OPEN (UNIT=2,FILE='out.out',FORM='FORMATTED')
      OPEN (UNIT=6,FILE='mie6.out',FORM='FORMATTED')
C     OPEN (UNIT=8,FILE='mie8.out',FORM='FORMATTED',STATUS='SCRATCH')
      OPEN (UNIT=8,FILE='mie8.out',FORM='FORMATTED')
      OPEN (UNIT=9,FILE='reff.out',FORM='FORMATTED')
C
C Write aerosol DB version name
      write(2,'(a)')'C AEROSOL DB :'//release_str
C
      CALL FINDPAR (1, 'Type of class')
      READ (1,*)NUMAEROS
      IF (NUMAEROS.EQ.1) THEN
        WRITE (6,*)' AEROSOL CLASS: Clean-Continental'
        WRITE (2,'(a,i3)')'C CLASS / DB INDEX : Clean-Continental  '
     +   ,numaeros

      ELSEIF (NUMAEROS.EQ.2) THEN
        WRITE (6,*)' AEROSOL CLASS: Average-Continental'
        WRITE (2,'(a,i3)')'C CLASS / DB INDEX : Average-Continental / '
     +   ,numaeros

      ELSEIF (NUMAEROS.EQ.3) THEN
        WRITE (6,*)' AEROSOL CLASS: Urban'
        WRITE (2,'(a,i3)')'C CLASS / DB INDEX : Urban / '
     +   ,numaeros


      ELSEIF (NUMAEROS.EQ.4) THEN
        WRITE (6,*)' AEROSOL CLASS: Clean-Maritime'
        WRITE (2,'(a,i3)')'C CLASS / DB INDEX : Clean-Maritime / '
     +   ,numaeros


      ELSEIF (NUMAEROS.EQ.5) THEN
        WRITE (6,*)' AEROSOL CLASS: Maritime-Polluted'
        WRITE (2,'(a,i3)')'C CLASS / DB INDEX : Maritime-Polluted / '
     +   ,numaeros


      ELSEIF (NUMAEROS.EQ.6) THEN
        WRITE (6,*)' AEROSOL CLASS: Desert Background (wintertime)'
        WRITE (2,'(a,i3)')'C CLASS / DB INDEX : Desert Background
     +   (wintertime) / ',numaeros


      ELSEIF (NUMAEROS.EQ.7) THEN
        WRITE (6,*)' AEROSOL CLASS: Desert Wind-Carry (summertime)'
        WRITE (2,'(a,i3)')'C CLASS / DB INDEX : Desert Wind-Carry
     +   (summertime) / ',numaeros


      ELSEIF (NUMAEROS.EQ.8) THEN
        WRITE (6,*)' AEROSOL CLASS: Background Stratospheric 1'
        WRITE (2,'(a,i3)')'C CLASS / DB INDEX : Background
     +   Stratospheric 1 ',numaeros


      ELSEIF (NUMAEROS.EQ.9) THEN
        WRITE (6,*)' AEROSOL CLASS: Volcanic 1'
        WRITE (2,'(a,i3)')'C CLASS / DB INDEX : Volcanic 1 / '
     +   ,numaeros

      ELSEIF (NUMAEROS.EQ.10) THEN
        WRITE (6,*)' AEROSOL CLASS: Maritime'
        WRITE (2,'(a,i3)')'C CLASS / DB INDEX : Maritime / '
     +   ,numaeros

      ELSEIF (NUMAEROS.EQ.11) THEN
        WRITE (6,*)' AEROSOL CLASS: Continental'
        WRITE (2,'(a,i3)')'C CLASS / DB INDEX : Continental / '
     +   ,numaeros

      ELSEIF (NUMAEROS.EQ.12) THEN
        WRITE (6,*)' AEROSOL CLASS: Urban-Industrial'
        WRITE (2,'(a,i3)')'C CLASS / DB INDEX : Urban-Industrial / '
     +   ,numaeros

      ELSEIF (NUMAEROS.EQ.13) THEN
        WRITE (6,*)' AEROSOL CLASS: Background Stratospheric 2'
        WRITE (2,'(a,i3)')'C CLASS / DB INDEX : Background
     +   Stratospheric 2 / ',numaeros

      ELSEIF (NUMAEROS.EQ.14) THEN
        WRITE (6,*)' AEROSOL CLASS: Volcanic 2'
        WRITE (2,'(a,i3)')'C CLASS / DB INDEX : Volcanic 2 / '
     +   ,numaeros


      ELSEIF (NUMAEROS.EQ.15) THEN
        WRITE (6,*)' AEROSOL CLASS: Rural (LOWTRAN)'
        WRITE (2,'(a,i3)')'C CLASS / DB INDEX : Rural (LOWTRAN) / '
     +   ,numaeros


      ELSEIF (NUMAEROS.EQ.16) THEN
        WRITE (6,*)' AEROSOL CLASS: Urban (LOWTRAN)'
        WRITE (2,'(a,i3)')'C CLASS / DB INDEX : Urban (LOWTRAN) / '
     +   ,numaeros


      ELSEIF (NUMAEROS.EQ.17) THEN
        WRITE (6,*)' AEROSOL CLASS: Maritime (LOWTRAN)'
        WRITE (2,'(a,i3)')'C CLASS / DB INDEX : Maritime (LOWTRAN) / '
     +   ,numaeros


      ELSEIF (NUMAEROS.EQ.18) THEN
        WRITE (6,*)' AEROSOL CLASS: Tropospheric (LOWTRAN)'
        WRITE (2,'(a,i3)')'C CLASS / DB INDEX : Tropospheric (LOWTRAN)
     +    / ',numaeros


      ELSEIF (NUMAEROS.EQ.19) THEN
        WRITE (6,*)' AEROSOL CLASS: Biomass Burning'
        WRITE (2,'(a,i3)')'C CLASS / DB INDEX : Biomass Burning / '
     +   ,numaeros


      ELSEIF (NUMAEROS.EQ.43) THEN
        WRITE (6,*)' AEROSOL CLASS: User defined'
        WRITE (2,'(a,i3)')'C CLASS / DB INDEX : User defined / '
     +  ,numaeros


      ELSEIF (NUMAEROS.EQ.44) THEN
        WRITE (6,*)' AEROSOL CLASS: Advection Fog 1 (heavy)'
        WRITE (2,'(a,i3)')'C CLASS / DB INDEX : Advection Fog 1
     +  (heavy) / ',numaeros


      ELSEIF (NUMAEROS.EQ.45) THEN
        WRITE (6,*)' AEROSOL CLASS: Advection Fog 2 (moderate)'
        WRITE (2,'(a,i3)')'C CLASS / DB INDEX : Advection Fog 2
     +  (moderate) / ',numaeros


      ELSEIF (NUMAEROS.EQ.46) THEN
        WRITE (6,*)' AEROSOL CLASS: Radiation Fog 3 (heavy)'
        WRITE (2,'(a,i3)')'C CLASS / DB INDEX : Radiation Fog 3 (heavy)
     +  / ',numaeros


      ELSEIF (NUMAEROS.EQ.47) THEN
        WRITE (6,*)' AEROSOL CLASS: Radiation Fog 4 (moderate)'
        WRITE (2,'(a,i3)')'C CLASS / DB INDEX : Radiation Fog 4
     +  (moderate) / ',numaeros

      ELSE

        WRITE (6,*)' DB COMPONENT INDEX (see entry.dat): ', numaeros
        WRITE (2,'(a,i3)')'C DB AEROSOL COMPONENT INDEX (see entry.dat)
     +    :',numaeros

      ENDIF

      integ=0
      ISEE=0
      iunits=0
      idb=0

C====================================================
C modified to add  radii scaling factor
C for Phil Watts (01-04-1999)
C===================================================

      CALL FINDPAR (1, 'Scaling radii factor')
      READ(1,*)NSCALE
c       write(*,*)scale(nscale)
c       write(*,*)(hum_radii(i,25),i=1,7)
      do iyhum=1,7
       do iscale=1,26
         hum_radii(iyhum,iscale)=HUM_RADII(iyhum,ISCALE)*
     + SCALE(NSCALE)
       enddo
      enddo
C
C==================================================
C
C  change the mixing ratios
C

      CALL FINDPAR (1, 'Type of mixing')
      READ(1,*)MIXTYPE
      WRITE(6,*)'TIPO DI MIXING=',MIXTYPE
      
      IF (NUMAEROS.NE.43) THEN
        CALL FINDPAR (1, 'User mixing ratio')
        READ(1,*)MIXR
        IF (MIXR.NE.0) THEN
          NNPROB=VEC_COMP ( NUMAEROS )
          WRITE(*,*)'The old mixing ratios are:'
            WRITE(*,*)(MIX_RATIO ( I,NUMAEROS ),i=1,NNPROB)
            WRITE(*,*)'Type the new mixing ratios:'
      WRITE(*,*)'in particle number density - normalized at 1part/cm^3'
            read(*,*)(MIX_RATIO ( I,NUMAEROS ),i=1,NNPROB)
           ENDIF
           
           IF (MIXTYPE.EQ.2.AND.NUMAEROS.LE.19) THEN
             NNPROB=VEC_COMP ( NUMAEROS )
             write(6,99)NNPROB
 99          FORMAT (' THIS CLASS IS MADE UP OF',1X,I2,1X,'COMPONENTS')
           ELSE
             NNPROB=1
           ENDIF
           
           CALL FINDPAR (1, 'Relative humidity level')
           READ (1,*)NRH
         WRITE(6,*)' HUMIDITY LEVEL=',NRH
         
       ELSE                                    !LEGGE LA USER DEFINED!!!
         CALL FINDPAR (1, 'User parameters')
         read(1,'(a)')header
c         print*,'C nome:  ',header
         write(2,'(a)')'C nome:  '//header(1:lnblnk(header))
         READ(1,*)NCOMP, NLAMBDA
         DO 201 IICOM=1,NCOMP
           READ(1,*)ANUMCOM, RMCOM, SIGMACOM
           MIX_RATIO(IICOM,1)=ANUMCOM
           HUM_RADII(1,IICOM)=RMCOM
           SIGMA(IICOM)=SIGMACOM
           DO 202 ILAMBDA=1,NLAMBDA
             READ(1,*) AMBDA(ILAMBDA), PARTRE(ILAMBDA,IICOM),
     &         PARTIM(ILAMBDA,IICOM)
 202       CONTINUE
 201     CONTINUE
         
        DO 204 IICOM=1,NCOMP
         DO 203 ILAMBDA=1,NLAMBDA
c         write(*,*) MIX_RATIO(IICOM,1),HUM_RADII(1,IICOM),SIGMA(IICOM)
c        write(*,*)AMBDA(ILAMBDA), PARTRE(ILAMBDA,IICOM),
c     +  PARTIM(ILAMBDA,IICOM)
 203    CONTINUE
 204    CONTINUE
           IF (MIXTYPE.EQ.2) THEN
             NNPROB=NCOMP
           ELSE
             NNPROB=1
           ENDIF
      ENDIF


C
      CALL FINDPAR (1, 'Problem number')
      READ (1,*)IPROB

C      WRITE(*,*)IPROB
C
C
      CALL FINDPAR (1,'Phase function tabulation')
      READ (1,*)NTYPEPH
C
      IF (NTYPEPH.EQ.1) THEN
         CALL FINDPAR (1,'Step size')
         READ( 1,*)DTHET1, CHANG1, DTHET2, CHANG2, DTHET3
C         WRITE(*,*)DTHET1, CHANG1, DTHET2, CHANG2, DTHET3
         KFG = (CHANG1+0.005)/DTHET1
         KFH = (180.0-CHANG2+0.005)/DTHET3
         KFF = (CHANG2-CHANG1+0.005)/DTHET2
         KFF = KFF+KFG
         KFH = KFH+KFF
         KFG1 = KFG+1
         KFG2 = KFG+2
         KFF1 = KFF+1
         KFF2 = KFF+2
         NTHETA = KFH+1
         THET(1) = 0.0
         DO 10 I = 2, KFG1
            THET(I) = THET(I-1)+DTHET1
  10    CONTINUE
         DO 20 I = KFG2, KFF1
            THET(I) = THET(I-1)+DTHET2
  20    CONTINUE
         DO 30 I = KFF2, NTHETA
            THET(I) = THET(I-1)+DTHET3
  30    CONTINUE

         DO 40 I = 1, NTHETA
            ANGLE = THET(I)*RAD
            COSA(I) = COS(ANGLE)
  40    CONTINUE
         THET(NTHETA) = 180.0
         COSA(NTHETA) = -1.0
C===================================================================
C  28-01-98 Modified to compute an other array of theta more dense
C           near theta=0 in order to compute more precisely the
C           Legendre coefficient of the phase function.
C====================================================================
         WRITE(6,96)NTHETA
  96     FORMAT (' PHASE FUNCTION TABULATION IS IN',1X,I4,1X,'ANGLES 
     + IN 3 INTERVALS' )
        ELSE IF (NTYPEPH.EQ.2) THEN
           CALL FINDPAR (1,'Gauss abscissae')
           READ (1,*)NTHETA
C           WRITE(*,*)NTHETA
             X1=0.0
             X2=180.0
              CALL GAULEG(X1,X2,COORD,WH,NTHETA-2)
             DO 275 J=2,NTHETA-1
                THET(J)=COORD(J-1)
  275        CONTINUE
              THET(1)=0.0
              THET(NTHETA)=180.0
              KFG1=INT(NTHETA/3.0)
              KFF1=INT(2*NTHETA/3.0)
         DO 41 I = 1, NTHETA
            ANGLE = THET(I)*RAD
            COSA(I) = COS(ANGLE)
  41    CONTINUE
         THET(NTHETA) = 180.0
         COSA(NTHETA) = -1.0
C
        WRITE(6,91)NTHETA
 91     FORMAT ('PHASE FUNCTION TABULATION IS IN',1X,I4,1X,
     + 'GAUSS ABSCISSAE')       
C
         ELSE
           CALL FINDPAR (1,'Number abscissae')
           READ (1,*)NTHETA
C           WRITE(*,*)NTHETA
          CALL FINDPAR (1,'Abscissae')
           IF (ntheta.LE.8) THEN
             READ (1,*)(THET(I), I=1,NTHETA)
           ELSE
             ND=INT(NTHETA/8)
             MD=NTHETA-(8*ND)
            DO 122 K=1,ND
               READ(1,*)(THET(I), I=1+(8*(K-1)),8*K)
  122       CONTINUE
               READ(1,*)(THET(I), I=(8*ND)+1,(8*ND)+MD)
           ENDIF
              KFG1=INT(NTHETA/3.0)
              KFF1=INT(2*NTHETA/3.0)
        DO 42 I = 1, NTHETA
            ANGLE = THET(I)*RAD
            COSA(I) = COS(ANGLE)
  42    CONTINUE
           WRITE(*,*)(THET(I), I=1,NTHETA)
           WRITE (6,*)'PHASE FUNCTION TABULATION IS READ FROM INPUT'
         ENDIF
C
         ZAPRXS=2.0E+5
         NSTEP=799
         IRADSP=1
         MOM=0

         IF (NUMAEROS.LE.19) THEN
          NNEQ= NEQVEC ( NUMAEROS )
            IF (NNEQ.eq.1) THEN
              WRITE(6,*)' LOG NORMAL SIZE DISTRIBUTION'
            ELSE
              WRITE(6,*)' MODIFIED GAMMA SIZE DISTRIBUTION'
            ENDIF
C                  ELSEIF (NUMAEROS.LE.42.OR.MIXTYPE.EQ.2) THEN
C                   NEQ=3

C         ELSEIF (MIXTYPE.EQ.1.AND.NCOMP.EQ.3) THEN
C           NEQ=4
C         ELSE
C           NEQ=3
         ENDIF
C
         NRS=NSTEP+1
         MODR=INT(NSTEP/100.0+1.)
         MNRS=MOD(NRS,MODR)

C
C  begin loop over componets
C

             size_int1=0.0
             size_int2=0.0
             size_int=0.0
        DO 200 JYZ=1,NNPROB

         IF (IRADSP.EQ.0) THEN
C            READ (15,905) (RI(I),ENR(I),I=1,NRS)
         ELSE
C         write(6,651) nrs
C 651         format(/,'up to size',i4)
c           print*,'up to size'
            CALL SIZE (RI,ENR,NUMAEROS,NRH,JYZ,mixtype)
c          print*,'past size'

C==========================================================
C    Modified for Phil Watts in order to compute the
C    effective radius relative to the chosen distribution

             call gauleg(rmin,rmax,ab,weig,900)
             indstart=2

          do i=1,900

               do j=indstart,nrs
                  if (ri(j).ge.ab(i)) then
                        ind_start=j
                        goto 3333
                  endif
               enddo
 3333         continue

           call inter2(1,ab(i),ri(j-1),ri(j),
     &        size_d(i),enr(j-1),enr(j))

         SIZE_INT=SIZE_INT+mix_ratio(jyz,numaeros)*size_d(i)*weig(i)
      size_int1=size_int1 +mix_ratio(jyz,numaeros)*size_d(i)*
     &              weig(i)*(ab(i)**3)
      size_int2=size_int2 +mix_ratio(jyz,numaeros)*size_d(i)*
     &              weig(i)*(ab(i)**2)


          enddo

C          write(*,*)mix_ratio(jyz,numaeros)
C          WRITE(*,*)size_int

C========================================================
            IF (IRADSP.EQ.2) THEN
               DO 50 KTIMES = 1, 2
c            print*,'up to newrad',ktimes
                  CALL NEWRAD (RI,ENR,NRS,MOM)
c            print*,'up to siz2',ktimes
             CALL SIZE (RI,ENR,NUMAEROS,NRH,JYZ,mixtype)
 50          CONTINUE
               WRITE (6,920) MOM
           ENDIF
         ENDIF

      IF (NUMAEROS.NE.43) THEN
            CALL FINDPAR (1, 'Type of wavelenghts')
            READ (1,*) NWLTYPE

         IF (NWLTYPE.EQ.2) THEN
             KPROB=9
            IF (JYZ.EQ.1) THEN
                WRITE(6,9)
  9             FORMAT('THE COMPUTATION WILL BE PERFORMED AT
     + THE 9 STANDARD WAVELENGTHS')
            ENDIF
         ELSE

            CALL FINDPAR (1, 'Number wavelengths')
            READ (1,*)KPROB
              IF (JYZ.EQ.1) THEN
              WRITE(6,92)KPROB
  92      FORMAT (' THE COMPUTATION WILL BE PERFOEMED FOR THESE
     +',1X,I2,1X,'WAVELENGTH VALUES')
           ENDIF
C
           CALL FINDPAR (1, 'Wavelength values')
           IF (KPROB.LE.10) THEN
             READ(1,*)(WAVE_NEW(I), I=1,KPROB)
           ELSE
             NDIME=INT(KPROB/10)
             MDIME=KPROB-(10*NDIME)
              DO 121 K=1,NDIME
               READ(1,*)(WAVE_NEW(I), I=1+(10*(K-1)),10*K)
  121         CONTINUE
               READ(1,*)(WAVE_NEW(I), I=(10*NDIME)+1,(10*NDIME)+MDIME)
           ENDIF
              IF (JYZ.EQ.1) THEN
                WRITE(6,*)(WAVE_NEW(I), I=1,KPROB)
              ENDIF
        ENDIF
C
      ELSE
          KPROB=NLAMBDA
      ENDIF
         IF (MIXTYPE.EQ.2)   THEN
             ICOMP=COMP_MAT(JYZ,NUMAEROS)
             I_RI=COMP_TO_RI(ICOMP)
         ENDIF

C
C    Loop over wavelength
C
      DO 300 NIII=1,KPROB

      IF (NUMAEROS.NE.43) THEN
C
       IF (NWLTYPE.EQ.2) THEN
C
C computation of complex refractive index at any humidity level
C
C
         IF (MIXTYPE.EQ.2) THEN
            CALL RH_RIND(NRH,NIII,ICOMP,I_RI,RI_R,RI_I)
             PM=CMPLX(RI_R,RI_I)
             write(*,*)PM
             WAVEL=WLR(NIII)
             WRITE(6,960)WAVEL,PM
 960      FORMAT(/,'WAVELENGTH =',F8.4,' UM, INDEX OF REFRACTION = (',
     *      e12.5,',',e12.5,')  ')
         ELSE
           RI_REAL=0.0
           RI_IMAG=0.0
              CALL VOL_PERC(NUMAEROS,NRH,C)
            DO 961 J=1,VEC_COMP(NUMAEROS)
              ICOMP=COMP_MAT(J,NUMAEROS)
              I_RI=COMP_TO_RI(ICOMP)
              CALL RH_RIND(NRH,NIII,ICOMP,I_RI,RI_R,RI_I)
C              write(*,*)c(j)
              RI_REAL = RI_REAL + (RI_R * C(J))
              RI_IMAG = RI_IMAG + (RI_I * C(J))
 961        CONTINUE
              WAVEL=WLR(NIII)
              PM=CMPLX(RI_REAL,RI_IMAG)
              WRITE(6,960)WAVEL,PM
         ENDIF

       ELSE
         IF (MIXTYPE.EQ.2) THEN
          DO 400 JJ=1,N_WLR
           IF (wlr(jj).gt.WAVE_NEW(NIII)) THEN
C
             W_NEW=WAVE_NEW(NIII)
C
C  computation of complex refractive indices at user wavelength
C                 for components and water
C
             CALL LINTER(W_NEW,JJ,1,I_RI,RI_IND_R)
             CALL LINTER(W_NEW,JJ,2,I_RI,RI_IND_I)
             CALL LINTER(W_NEW,JJ,1,11,RI_H2_R)
             CALL LINTER(W_NEW,JJ,2,11,RI_H2_I)
C
C   Hanel formula at user wavelngth
C

C             write(*,*)RI_IND_R,RI_IND_I,RI_H2_R,RI_H2_I

              RI_R = RI_H2_R + ((RI_IND_R - RI_H2_R)*
     +            (HUM_RADII(1,ICOMP)/HUM_RADII(NRH,ICOMP))**3)
              RI_I = RI_H2_I + ((RI_IND_I - RI_H2_I)*
     +            (HUM_RADII(1,ICOMP)/HUM_RADII(NRH,ICOMP))**3)
c              write(*,*)HUM_RADII(1,ICOMP),HUM_RADII(NRH,ICOMP)

C              write(*,*)ri_r, ri_i
C
              PM=CMPLX(RI_R,RI_I)
              WAVEL=WAVE_NEW(NIII)
              WRITE(6,960)WAVEL,PM
              GO TO 325
          ENDIF
 400          CONTINUE
       ELSE
          RI_R=0.0
          RI_I=0.0
          DO 405 JJ=1,N_WLR
            IF (wlr(jj).gt.WAVE_NEW(NIII)) THEN
C
              W_NEW=WAVE_NEW(NIII)
C
C  computation of complex refractive indices at user wavelength
C                 for components and water
C
              CALL VOL_PERC(NUMAEROS,NRH,C)
              DO 404 J=1,VEC_COMP(NUMAEROS)
               ICOMP=COMP_MAT(J,NUMAEROS)
               I_RI=COMP_TO_RI(ICOMP)

               CALL LINTER(W_NEW,JJ,1,I_RI,RI_IND_R)
               CALL LINTER(W_NEW,JJ,2,I_RI,RI_IND_I)
               CALL LINTER(W_NEW,JJ,1,11,RI_H2_R)
               CALL LINTER(W_NEW,JJ,2,11,RI_H2_I)
C
C   Hanel formula at user wavelength
C
C              write(*,*)c(j)
              RI_R = RI_R +(( RI_H2_R + ((RI_IND_R-RI_H2_R)*
     +          (HUM_RADII(1,ICOMP)/HUM_RADII(NRH,ICOMP))**3))* C(J))
              RI_I = RI_I +(( RI_H2_I + ((RI_IND_I - RI_H2_I)*
     +          (HUM_RADII(1,ICOMP)/HUM_RADII(NRH,ICOMP))**3))*C(J))
 404          CONTINUE
C
               PM=CMPLX(RI_R,RI_I)
               WAVEL=WAVE_NEW(NIII)
               WRITE(6,960)WAVEL,PM
               GO TO 325
               ENDIF
 405       CONTINUE
          ENDIF
 325     ENDIF

       ELSE
            WAVEL=AMBDA(NIII)
            PM=CMPLX(PARTRE(NIII,JYZ),PARTIM(NIII,JYZ))
       ENDIF
C********************* MIE2 ****************************************
c            write(*,*)MIX_RATIO(jyz,1)
c            write(*,*)HUM_RADII(1,jyz)
c            write(*,*)SIGMA(jyz)
            IF (ZAPRXS.GT.0.AND.REAL(PM).LE.1.) WRITE(*,943)PM
            IF (AIMAG(PM).GT.0.) PM = CONJG(PM)
            APM = CABS(PM)
            PLOGIM = -99.
            IF (AIMAG(PM).NE.0) PLOGIM = ALOG10(-AIMAG(PM))
            ALPH = WAVEL/PI2
            WOP = 1.0E-03*ALPH**2
            CEXT = 0.0
            SCOE = 0.0
            ASYM1 = 0.
            KR = 0
            DO 60 I = 1, NTHETA
               PHAS1(I) = 0.0
               PHAS2(I) = 0.0
               PHAS3(I) = 0.0
               PHAS4(I) = 0.0
   60       CONTINUE
C
C           LOOP FOR BASIC MIE CALCULATIONS AND INTEGRATION OVER SIZE
C
            DO 90 K = 1, NRS
               X = RI(K)/ALPH
               IF (X.LE.0.0) GO TO 90
               XDIM = 1.01*(X+4.05*X**(1./3.)+2.0)
C
C              CALCULATE PHASE FUNCTION AND CROSS SECTIONS
C
               IF (PLOGIM.LE.-4.*ALOG10(X)-14.) THEN
                  CALL RAYLIM (X,RENT1,RENT2,ENT3,ENT4,EXKRO,SCKRO,COSAV
     *               )
                  RECN = 0.
               ELSEIF (REAL(PM).LE.1) THEN
                  IF (XDIM.GT.MAXDIM) THEN
                     WRITE (*,944)
                     STOP
                  ELSE
                     CALL RRA42 (X,RENT1,RENT2,ENT3,ENT4,EXKRO,SCKRO,
     *                  COSAV,RECN)
                  ENDIF
               ELSEIF (XDIM.LE.MAXDIM.AND.(X.LT.ZAPRXS.OR.ZAPRXS.LE.0))
     *            THEN
                  CALL RRA42 (X,RENT1,RENT2,ENT3,ENT4,EXKRO,SCKRO,COSAV,
     *               RECN)
               ELSE
                  IF (XDIM.GT.MAXDIM.AND.ZAPRXS.LE.0.AND.NPRNT.EQ.0)
     *                THEN
                     NPRNT = 1
                     WRITE (*,945)
                  ENDIF
                  PMCONJ = CONJG(PM)
                  GAMMA = (2./X)**(1./3.)
                  EXKRO = QEXAPX(PMCONJ,X)
                  SCKRO = EXKRO-QABAPX(PMCONJ,X,GAMMA)
                  COSAV = (EXKRO-QPRAPX(PMCONJ,X,GAMMA))/SCKRO
                  NP = 3
                  CALL RAYOPT (X,NP,THET,RENT1,RENT2,ENT3,ENT4)
                  RECN = 0.
               ENDIF
               IF (ISEE.GT.0) THEN
C
C                 PRINT RESULTS FOR EACH AEROSOL RADIUS
C
                  WRITE (6,926) IPROB
                  IF (IUNITS.EQ.0) THEN
                     WRITE (6,946) WAVEL,PM,PML
                     WRITE (6,928) RI(K),X
                  ENDIF
                  IF (IUNITS.EQ.1) THEN
                     WRITE (6,947) WAVEL,PM,PML
                     WRITE (6,929) RI(K),X
                  ENDIF
C                  WRITE (6,927)
C                  DO 70 I = 1, NTHETA
C                     ENTAVG(I) = (RENT1(I)+RENT2(I))/2.0
C                     WRITE (6,909) THET(I),RENT1(I),RENT2(I),ENT3(I),
C     *                  ENT4(I),ENTAVG(I)
C   70             CONTINUE
                  ABKRO = EXKRO-SCKRO
                  SSA = SCKRO/EXKRO
                  WRITE (6,930) EXKRO,SCKRO,ABKRO,SSA,COSAV,RECN
              ENDIF
C              CUMULATIVE INTEGRATION OVER SIZE DISTRIBUTION
C
               IF (K.GT.1.AND.K.LT.NRS) DELR = (RI(K+1)-RI(K-1))/2.0
               IF (K.EQ.NRS) DELR = (RI(K)-RI(K-1))/2.0
               IF (K.EQ.1) DELR = (RI(K+1)-RI(K))/2.0
               EDW = ENR(K)*DELR*WOP
               DO 80 I = 1, NTHETA
                  PHAS1(I) = PHAS1(I)+EDW*RENT1(I)
                  PHAS2(I) = PHAS2(I)+EDW*RENT2(I)
                  PHAS3(I) = PHAS3(I)+EDW*ENT3(I)
                  PHAS4(I) = PHAS4(I)+EDW*ENT4(I)
   80          CONTINUE
               EDRR = ENR(K)*DELR*RI(K)*RI(K)
               CEXT = CEXT+EXKRO*EDRR
               SCOE = SCOE+SCKRO*EDRR
               ASYM1 = ASYM1+COSAV*SCKRO*EDRR
               IF (MOD(K,MODR).EQ.MNRS) THEN
                  KR = KR+1 
                  R(KR) = RI(K)
                  P(KR,1) = PHAS1(1)+PHAS2(1)
                  P(KR,2) = PHAS1(KFG1)+PHAS2(KFG1)
                  P(KR,3) = PHAS1(KFF1)+PHAS2(KFF1)
                  P(KR,4) = PHAS1(NTHETA)+PHAS2(NTHETA)
c                  write(*,*)kr, P(KR,1),p(kr,2),p(kr,3),p(kr,4)
               ENDIF
   90       CONTINUE
C
C           END LOOP FOR MIE CALCULATIONS
C
            KRMAX = KR
            ASYM1 = ASYM1/SCOE
            CEXT = PII*CEXT*1.0E-03
            SCOE = PII*SCOE*1.0E-03
            IF (IDB.EQ.1) THEN
               CEXT = CEXT*FACTDB
               SCOE = SCOE*FACTDB
               DO 95 I = 1, NTHETA
                  PHAS1(I) = PHAS1(I)*FACTDB
                  PHAS2(I) = PHAS2(I)*FACTDB
                  PHAS3(I) = PHAS3(I)*FACTDB
                  PHAS4(I) = PHAS4(I)*FACTDB
   95          CONTINUE
            ENDIF
            ABOE = CEXT-SCOE
C                  WRITE (6,930) EXKRO,SCKRO,ABKRO,SSA,COSAV,RECN
            SSA = SCOE/CEXT
C
            DO 100 KR = 1, KRMAX
               DO 100 I = 1, 4
                  P(KR,I) = P(KR,I)/P(KRMAX,I)
  100       CONTINUE
C
C           OUTPUT THE CUMULATIVE INTEGRAL OF THE PHASE FUNCTION AND
C           THE PHASE FUNCTION
C
C            WRITE (6,926) IPROB
C            IF (IUNITS.EQ.0) WRITE (6,946) WAVEL,PM,PML
C            IF (IUNITS.EQ.1) WRITE (6,947) WAVEL,PM,PML
C            WRITE (6,931) THET(1),THET(KFG1),THET(KFF1),THET(NTHETA)
C            WRITE (6,932) INT(THET(1)),INT(THET(KFG1)),INT(THET(KFF1)),
C     *         INT(THET(NTHETA)),INT(THET(1)),INT(THET(KFG1)),
C     *         INT(THET(KFF1)),INT(THET(NTHETA)),INT(THET(1)),
C     *         INT(THET(KFG1)),INT(THET(KFF1)),INT(THET(NTHETA))
C            WRITE (6,912) (R(KR),(P(KR,I),I=1,4),KR=1,KRMAX)
C
C            WRITE (6,926) IPROB
C            IF (IUNITS.EQ.0) WRITE (6,946) WAVEL,PM,PML
C            IF (IUNITS.EQ.1) WRITE (6,947) WAVEL,PM,PML
C            IF (IDB.EQ.0) WRITE (6,933)
C            IF (IDB.EQ.1) WRITE (6,934)
C            WRITE (6,927)
            DO 110 I = 1, NTHETA
               PHAS5(I) = (PHAS1(I)+PHAS2(I))/2.0
C               WRITE (6,909) THET(I),PHAS1(I),PHAS2(I),PHAS3(I),PHAS4(I)
C     *            ,PHAS5(I)
  110       CONTINUE
C
C           COMPUTE THE ASYMMETRY PARAMETER BY INTEGRATING OVER THET
C
            APHAS(1) = 0.0
            ASYM2 = 0.0
            IF (INTEG.GT.0) THEN
               KEEP = NTHETA
               DO 120 I = 1, NTHETA
                  F(I) = PI2*PHAS5(I)
  120          CONTINUE
               CALL IPLAW (KEEP,COSA,F,APHAS)
               F(1) = PHAS5(1)
               DO 130 I = 2, NTHETA
                  APHAS(I) = APHAS(I)+APHAS(I-1)
                  F(I) = PHAS5(I)*COSA(I)
  130          CONTINUE
               CALL IPLAW (KEEP,COSA,F,POL)
               H(1) = 0.0
               DO 140 I = 2, NTHETA
                  H(I) = POL(I)+H(I-1)
  140          CONTINUE
               TPHAS = APHAS(NTHETA)
               ASYM2 = H(NTHETA)*PI2/TPHAS
            ELSE
               DO 150 I = 2, NTHETA
                  DCOS = COSA(I-1)-COSA(I)
                  PDC = PII*DCOS
                  ASYM2 = ASYM2+PDC*(PHAS5(I-1)*COSA(I-1)+PHAS5(I)*COSA(
     *               I))
                  APHAS(I) = APHAS(I-1)+PDC*(PHAS5(I-1)+PHAS5(I))
  150          CONTINUE
               TPHAS = APHAS(NTHETA)
               ASYM2 = ASYM2/TPHAS
            ENDIF
              
C
C           COMPUTE NORMALIZED PHASE FUNCTION, DEGREE OF POLARIZATION,
C           AND ELLIPTICITY. OUTPUT THE RESULTS
C
            DO 160 I = 1, NTHETA
               PHAS1(I) = PHAS1(I)/SCOE
               PHAS2(I) = PHAS2(I)/SCOE
               PHAS3(I) = PHAS3(I)/SCOE
               PHAS4(I) = PHAS4(I)/SCOE
               PHAS5(I) = PHAS5(I)/SCOE
               Q = (PHAS1(I)-PHAS2(I))/2.0
               POL(I) = 100.0*Q/PHAS5(I)
               W = SQRT(Q**2+PHAS3(I)**2+PHAS4(I)**2)
               IF (ABS(W).GT.1.0E-15) THEN
                  BETA = 0.5*ASIN(PHAS4(I)/W)
                  ELL(I) = TAN(BETA)
               ELSE
                  ELL(I) = 0.0
               ENDIF
  160       CONTINUE
C
C            WRITE (6,926) IPROB
C            IF (IUNITS.EQ.0) WRITE (6,946) WAVEL,PM,PML
C            IF (IUNITS.EQ.1) WRITE (6,947) WAVEL,PM,PML
C            WRITE (6,935)
C            WRITE (6,927)
            DO 170 I = 1, NTHETA
C               WRITE (6,909) THET(I),PHAS1(I),PHAS2(I),PHAS3(I),PHAS4(I)
C     *            ,PHAS5(I)
               ANGOLO(I)=THET(I)
C******************** MODIFIED BY CHIARA ON THE 24-08-96 ********
               IF (NNPROB.NE.1) THEN
                 PHAS_F(I,NIII,JYZ)=PHAS5(I)
               ELSE
                 PHAS_FUN(I,NIII)=PHAS5(I)*4.0*PII
               ENDIF
  170       CONTINUE
C            WRITE (6,936) CEXT,SCOE,ABOE
C************ MODIFIED BY CHIARA 2-04-96 *******************
            IF (NNPROB.NE.1) THEN
                EXT_COEFF(NIII,JYZ)=CEXT
                SCAT_COEFF(NIII,JYZ)=SCOE
                SING_S_ALB(NIII,JYZ)=SSA
                ASYM_FAC(NIII,JYZ)=ASYM1
            ELSE
             EXT(NIII)=CEXT
             SCA(NIII)=SCOE
             SS_ALB(NIII)=SSA
             ASYM(NIII)=ASYM1
            ENDIF
C            WRITE (2,956) CEXT,SSA,ASYM1
C**********************************************************
C            WRITE (6,937) SSA,ASYM1
C            WRITE (6,938) TPHAS,ASYM2
C            IF (IDB.EQ.0) WRITE (6,954)
C            IF (IDB.EQ.1) WRITE (6,955)
C            CALL SECOND (TIMN)
C            WRITE (6,925) TIMN-TIMS,TIMN
C
C            WRITE (6,926) IPROB
C            IF (IUNITS.EQ.0) WRITE (6,946) WAVEL,PM,PML
C            IF (IUNITS.EQ.1) WRITE (6,947) WAVEL,PM,PML
C            WRITE (6,940)
C            WRITE (6,941)
C            DO 180 I = 1, NTHETA
C               WRITE (6,910) THET(I),POL(I),ELL(I)
C  180       CONTINUE
C
c            IF (IUNITS.EQ.0) WRITE (8,948) WAVEL,PM,PML
c            IF (IUNITS.EQ.1) WRITE (8,949) WAVEL,PM,PML
             WRITE (6,939) CEXT,SCOE,ABOE,ASYM1
            DO 190 II = 1, NTHETA
               WRITE (8,911) THET(II),PHAS1(II),PHAS2(II),PHAS3(II),
     *            PHAS4(II)
  190       CONTINUE
C
C**** MODIFIED by Chiara Levoni 02-04-96 *****************
        WAVE(NIII)=WAVEL
C********************************************************
C  integral of the phase function
c
c          x1=thet(1)
c          x2=PII
c          call GAULEG(X1,X2,absci,WH,NTHETA-2)
c          do 175 j=2,ntheta-1
c           pesi(j)=wh(j-1)
c  175     continue
c           pesi(1)=0.0
c           pesi(ntheta)=0.0
c          gapmie2=0.0
c           do 173 k=1,ntheta
c             gapmie2=gapmie2+phas5(k)*pesi(k)*
c     +          (sin(thet(k)*PII/180.0))
c  173     continue
c          gapmie2=gapmie2*2.0d0*PII
c          write(*,*)'INTEGRALE DI MIE2 (ultima colonna)',gapmie2
C******************* MIE2(the end) ************************************
 300    CONTINUE
C  End of loop over wavelength             
                         
 200    CONTINUE        
C  End of loop over problem
       WRITE(9,*)'CLASSE',NUMAEROS,'SCALE=',NSCALE
       WRITE(9,*)'RAGGIO EFFETTIVO=', SIZE_INT1/SIZE_INT2

C END LOOP OVER SCALING FACTOR
C
C  if the aerosol is made of more than 1 component 
C  the weighted average of optical properties is performed
C
       IF (NNPROB.NE.1) THEN
         DO 301 Icnt=1,kprob
           
c       DO 303 Jcnt=1,VEC_COMP(NUMAEROS)
C
           
c MODFICA FATTA il 12-2-2003 per far funzionare
c l'user defined con piu' componenti
c
C IPOTESI: che NPROB sia il numero di componenti!!!
c sia in verione con le classi sia con l'user defined.

           DO 303 Jcnt=1,NNPROB
             
             IF (NUMAEROS.NE.43) THEN
               EXT(Icnt)=EXT(Icnt)+EXT_COEFF(Icnt,Jcnt)*
     &           MIX_RATIO(Jcnt,NUMAEROS)
               SCA(Icnt)=SCA(Icnt)+SCAT_COEFF(Icnt,Jcnt)*
     &           MIX_RATIO(Jcnt,NUMAEROS)
               ASY(Icnt)=ASY(Icnt)+(ASYM_FAC(Icnt,Jcnt)*
     &           MIX_RATIO(Jcnt,NUMAEROS)*SCAT_COEFF(Icnt,Jcnt))
             ELSE
               EXT(Icnt)=EXT(Icnt)+EXT_COEFF(Icnt,Jcnt)*
     &           MIX_RATIO(Jcnt,1)
               SCA(Icnt)=SCA(Icnt)+SCAT_COEFF(Icnt,Jcnt)*
     &           MIX_RATIO(Jcnt,1)
               ASY(Icnt)=ASY(Icnt)+(ASYM_FAC(Icnt,Jcnt)*
     &           MIX_RATIO(Jcnt,1)*SCAT_COEFF(Icnt,Jcnt))
             ENDIF
C
             
 303       CONTINUE
           SS_ALB(Icnt)=SCA(Icnt)/EXT(Icnt)
           ASYM(Icnt)=ASY(Icnt)/SCA(Icnt)
C        WRITE(2,959)WAVE(Icnt),EXT,SCA,SS_ALB,ASYM
 301     CONTINUE
         DO 402 NW=1,kprob
           DO 401 IANG=1,NTHETA

c riga vecchia!!!
c            DO 403 jcnt=1,VEC_COMP(NUMAEROS)
c
c riga nuova (12-2-2003)
             DO 403 jcnt=1,NNPROB
               

               IF (NUMAEROS.NE.43) THEN
             PHAS_FUN(IANG,NW)=PHAS_FUN(IANG,NW)+PHAS_F(IANG,NW,JCNT)*
     &             MIX_RATIO(Jcnt,NUMAEROS)*SCAT_COEFF(NW,Jcnt)
               ELSE
             PHAS_FUN(IANG,NW)=PHAS_FUN(IANG,NW)+PHAS_F(IANG,NW,JCNT)*
     &             MIX_RATIO(Jcnt,1)*SCAT_COEFF(NW,Jcnt)
               ENDIF 
 403         CONTINUE
             PHAS_FUN(IANG,NW)=PHAS_FUN(IANG,NW)/SCA(NW)*4.0*PII
 401       CONTINUE
 402     CONTINUE
       ENDIF
C phase function expantion in Legendre polynomial
c       write(*,*)ntheta
       CALL FINDPAR(1, 'Legendre polynomial')
       READ(1,*) JPOLY
       CALL FINDPAR(1, 'Terms')
       READ(1,*)ELLE
       IF (ELLE.GT.300) THEN
           ELLE=68
           WRITE(*,*)'the number of Legendre coeff. you have chosen'
           WRITE(*,*)'is greater than 300, so it has been set to 68'
       ENDIF
       IF (JPOLY.EQ.1) THEN
C        do iwl=1,KPROB
C           do i=1,NTHETA
C            phas_fun(i,iwl)=1.0
C           ENDDO
C        ENDDO
             CALL LEGENDRE_COEFF(KPROB,NTHETA,THET,
     *            WAVE,XI,ELLE,ASYM)
C

C 
             CALL WRITEFILE(JPOLY,KPROB,IPROB,NTHETA,WAVE,EXT,
     *            SS_ALB,ASYM,THET,XI,ELLE,nscale,NUMAEROS)



c        do iwl=1,KPROB
c           do i=1,NTHETA                           
c            write(*,*)I,IWL,phas_fun(i,iwl)
c           ENDDO
c        ENDDO

       ELSE   
 
C  calcolo raggi effettivi per watts. Dopo scommenta writefile      
         CALL WRITEFILE(JPOLY,KPROB,IPROB,NTHETA,WAVE,EXT,SS_ALB,
     *            ASYM,THET,XI,ELLE,nscale,NUMAEROS)
       ENDIF

C Append input file entry.dat
       rewind(1)
       call append_f(2,1,"entry.dat")
       rewind(9)
       call append_f(2,9,"reff.out")

       STOP
C
C**********************************************************
C

  901 FORMAT (3I5)
  902 FORMAT (4I5,20A4)
  903 FORMAT (5F10.2)
  904 FORMAT (3E10.3,4I5,A20)
  905 FORMAT (2E15.5)
C  906 FORMAT (3E15.8,A20)
  906 FORMAT (3E15.8)
  907 FORMAT (5(1X,0PF10.5,1PE13.5,2X))
  908 FORMAT (1P3E10.3,4I5,A20)
  909 FORMAT (8X,F6.2,5X,1P5E17.6)
  910 FORMAT (23X,F6.2,14X,F11.6,15X,1PE13.6)
  911 FORMAT (1P,6(E10.3,2x))
  912 FORMAT (3(F10.4,4F7.4))
  913 FORMAT (20A4)
  914 FORMAT ('1 TOTAL TIME =',F12.4,//,'  INPUT DATA',/)
  915 FORMAT ('  NNPROB, INTEG, ISEE =',3I5,/)
  916 FORMAT ('  IPROB, KPROB, IUNITS, IDB =',4I5,//,2X,20A4,/)
  917 FORMAT ('  DTHET1, CHANG1, DTHET2, CHANG2, DTHET3 =',5F8.2,/)
  918 FORMAT ('  RMIN, RMAX, ZAPRXS =',1P3E10.3,/,
     *   '  NSTEP, IRADSP, MOM, NNEQ =',4I5,2X,A20,/)
  919 FORMAT (5F10.2,2X,'NO. OF ANGLES =',I5)
  920 FORMAT (/,'  RADII UNIFORMLY SPACED IN N(R)*R**',I2)
  921 FORMAT (/,30X,'  INPUT AEROSOL SIZE DISTRIBUTION FOR',I4,
     *   ' RADII (RADII EXPRESSED IN UM)',/)
  922 FORMAT (/,30X,'  INPUT AEROSOL SIZE DISTRIBUTION FOR',I4,
     *   ' RADII (RADII EXPRESSED IN MM)',/)
  923 FORMAT (5(7X,'R',9X,'N(R)',5X))
  924 FORMAT (/,'  TOTAL TIME =',F12.4,/)
  925 FORMAT (//,'  TIME DIFF =',F12.4,'  TOTAL TIME =',F12.4)
  926 FORMAT ('1',16X,'AIR FORCE GEOPHYSICS LABORATORY ',6X,
     *   'PROGRAM MIE2',6X,'PROBLEM',I7,/)
  927 FORMAT (/,3X,'SCATTERING ANGLE',10X,'I1',15X,'I2',15X,'I3',15X,
     *   'I4',14X,'IAVG'//)
  928 FORMAT (/,2X,'PHASE MATRIX FOR AEROSOL WITH RADIUS',1PE12.5,
     *   ' UM AND SIZE PARAMETER',1PE12.5,2X,'(UNITS ARE STR-1)')
  929 FORMAT (/,2X,'PHASE MATRIX FOR AEROSOL WITH RADIUS',1PE12.5,
     *   ' MM AND SIZE PARAMETER',1PE12.5,2X,'(UNITS ARE STR-1)')
  930 FORMAT (/,2X,'EXTINCTION EFFICIENCY = ',1PE12.6,/,2X,
     *   'SCATTERING EFFICIENCY = ',1PE12.6,/,2X,
     *   'ABSORPTION EFFICIENCY = ',1PE12.6,/,2X,
     *   'SINGLE SCAT. ALBEDO   = ',0PF8.6,/,2X,
     *   'ASYMMETRY PARAMETER   = ',0PF8.6,/,2X,
     *   'RECURSION CYCLES      = ',0PF6.0)
  931 FORMAT (/,29X,' CUMULATIVE INTEGRAL OF PHASE FUNCTION VS SIZE',/,
     *   30X,'FOR THETA =',3(F7.2,','),F7.2,' DEGREES',/)
  932 FORMAT (3(7X,'R',2X,4(' P(',I3,')')))
  933 FORMAT (/,10X,'UNNORMALIZED PHASE MATRIX FOR AEROSOL SIZE',
     *   ' DISTRIBUTION  (UNITS ARE KM-1*STR-1)')
  934 FORMAT (/,10X,'UNNORMALIZED PHASE MATRIX FOR AEROSOL SIZE',
     *   ' DISTRIBUTION  (UNITS ARE DB*KM-1*STR-1)')
  935 FORMAT (/,16X,'NORMALIZED PHASE MATRIX FOR AEROSOL SIZE',
     *   ' DISTRIBUTION  (UNITS ARE STR-1)')
  936 FORMAT (/,2X,'MACROSCOPIC EXTINCTION COEFFICIENT =',1PE12.5,/,2X
     *   ,'MACROSCOPIC SCATTERING COEFFICIENT =',1PE12.5,/,2X,
     *   'MACROSCOPIC ABSORPTION COEFFICIENT =',1PE12.5)
  937 FORMAT (2X,'SINGLE SCATTERING ALBEDO',11X,'=',5X,F7.5,/,2X,
     *   'ASYMMETRY PARAMETER',16X,'=',5X,F7.5)
  938 FORMAT (/,'  VALUES BY INTEGRATING PHASE FUNCTION:',/,2X,
     *   'MACROSCOPIC SCATTERING COEFFICIENT =',1PE12.5,/,2X,
     *   'ASYMMETRY PARAMETER',16X,'=',5X,0PF7.5)
  939 FORMAT (' EXT,SCAT,ABS,ASYM =',1P,3(G12.5,1X),0PF7.5)
  940 FORMAT (/,29X,'POLARIZATION DATA FOR AEROSOL SIZE DISTRIBUTION')
  941 FORMAT (/,18X,'SCATTERING ANGLE',4X,'DEGREE OF POLARIZATION',4X,
     *   'ELLIPTICITY (TAN(BETA))',//)
  942 FORMAT ('  ZAPRXS SET TO 10, WHICH IS THE MINIMUM SIZE PARAMETER',
     *   /,'  WHERE THE MIE APPROXIMATIONS ARE VALID',/)
  943 FORMAT ('  THE MIE APPROXIMATIONS ARE NOT USED FOR PM =',2F10.4,/,
     *   '  SINCE THEY ARE NOT VALID FOR REAL PARTS LESS THAN ONE',/)
  944 FORMAT ('  FMIE2 ABORTED:  REAL INDEX LESS THAN ONE AND NUMBER ',
     *   'OF TERMS',/,'  IN MIE EXPANSION IS GREATER THAN MAXDIM')
  945 FORMAT ('  WARNING: ZAPRXS = 0, BUT THE APPROXIMATIONS ARE USED ',
     *   'SINCE',/,'  THE NUMBER OF TERMS IN THE MIE EXPANSION',
     *   ' EXCEEDS MAXDIM')
  946 FORMAT (18X,'WAVELENGTH =',F8.4,' UM, INDEX OF REFRACTION = (',
     *   F9.6,',',F10.6,')  ',A20,/)
  947 FORMAT (18X,'WAVELENGTH =',F8.4,' MM, INDEX OF REFRACTION = (',
     *   F9.6,',',F10.6,')  ',A20,/)
  948 FORMAT (' WAVELENGTH =',F8.3,' UM, INDEX OF REF. = ',2F9.5,2X,A20)
  949 FORMAT (' WAVELENGTH =',F8.3,' MM, INDEX OF REF. = ',2F9.5,2X,A20)
  950 FORMAT (2I5,'   WAVELENGTHS IN UM AND COEFFICIENTS IN KM-1')
  951 FORMAT (2I5,'   WAVELENGTHS IN UM AND COEFFICIENTS IN DB*KM-1')
  952 FORMAT (2I5,'   WAVELENGTHS IN MM AND COEFFICIENTS IN KM-1')
  953 FORMAT (2I5,'   WAVELENGTHS IN MM AND COEFFICIENTS IN DB*KM-1')
  954 FORMAT (/,2X,'UNITS OF COEFFICIENTS ARE KM-1')
  955 FORMAT (/,2X,'UNITS OF COEFFICIENTS ARE DB*KM-1')
C***** MODIFIED by Chiara Levoni 02-04-96 ************************
  956 FORMAT (/,2X,'MACROSCOPIC EXTINCTION COEFFICIENT =',1PE12.5,/,2X,
     *         'SINGLE SCATTERING ALBEDO',11X,'=',5X,e12.5,/,2X,
     *         'ASYMMETRY PARAMETER',16X,'=',5X,e12.5)
  957 FORMAT (/,2X,'WAVELENGTH =',F8.4,/,2X,
     *          'MIXTURE EXTINCTION COEFFICIENT =',1PE12.5,/,2X,
     *          'MIXTURE SCATTERING COEFFICIENT =',1PE12.5,/,2X,
     *          'MIXTURE SINGLE SCATTER. ALBEDO =',0PF8.6,/,2X,
     *          'ASYMMETRY PARAMETER   = ',0PF8.6,/,2X)
C 958 FORMAT ( '--- WL --- EXT --- SCA -- S.S.ALB
C     * -- ASYM ---')
  959 FORMAT (/2X,F8.4,2X,1PE12.5,2X,1PE12.5,2X,0PF8.6,2X,0PF8.6,2X)
C******************************************************************
C

       END

