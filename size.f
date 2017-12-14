
      SUBROUTINE SIZE (RI,ENR,NUMAEROS,NRH,NELEM,MIXTYPE)
C
C     ROUTINE TO READ IN THE AEROSOL SIZE DISTRIBUTION PARAMETERS BASED
C     ON THE VALUE OF NEQ. THE ROUTINE THEN COMPUTES THE NUMBER OF
C     PARTICLES AT RI(I) BEGINNING AT ENTRY POINT SIZ2
C     NEQ = 1: MODIFIED GAMMA DISTRIBUTION
C     NEQ = 2: TRUNCATED POWER LAW
C     NEQ = 3: SUM OF TWO LOG NORMAL DISTRIBUTIONS
C     NEQ = 4: SUM OF THREE LOG NORMAL DISTRIBUTIONS
C     NEQ = 5: MARSHALL-PALMER RAINDROP DISTRIBUTION
C     NEQ = 6: GENERALIZED EXPONENTIAL RAINDROP DISTRIBUTION
C**************************************************************
C  This routine has been modified on the 23-04-97.
C  The common block "s_d_param" has been added in
C  order to share the size distribution parameters
C  with the new subroutine RMM that computes the
C  minimum and maximum radius for the integration
C  over the size distribution.
C**************************************************************

c size-distrubution cutoff value to select rmin and rmax
        integer modanum,neq
        real lognorpar(3,10)
        real gammodpar(4),top(10),bi(10), abs(900),
     +     size_d(900),weig(900)

	 parameter (cutoff=1.e-6)
         parameter (dummy=-99.)

      INCLUDE 'AEROS_DB.PARS'
      INCLUDE 'AEROS_DB.VARS'
      DIMENSION RI(*),ENR(*)

      COMMON /A1R0/ NRS,PII,RMIN,RMAX,SQRPI2,NSTEP
      COMMON /s_d_param/ lognorpar,gammodpar,neq,modanum
C
C
      do k=1,nrs
        enr(k)=0.0
      enddo
C
      NEQ=NEQVEC(NUMAEROS)
C
      IF (NUMAEROS.EQ.43) THEN
          NRH=1
      ENDIF

C       IF (NUMAEROS.LE.19.and.mixtype.eq.2) THEN
C          NNPROB=VEC_COMP ( NUMAEROS )
C       ELSE
C          NNPROB=1
C       ENDIF

C       write(*,*)(MIX_RATIO ( I,NUMAEROS ),i=1,NNPROB)

C modified gamma size distribtuon
C
        IF (NEQ.EQ.2) THEN
C
           IF (NUMAEROS.EQ.13) THEN          !75% H2SO4 2 (Background Stratospheric 2)
              modanum=1                      !because it's  gamma_mod(i,1) - see aeros_bb.f
              do i=1,4
                gammodpar(i)=gamma_mod(i,1)
              enddo
c              write (6,910) A,alpha,b,gamma
              WRITE(6,910)(gammodpar(i), i=1,4)
	      call rmm (cutoff, rmin, rmax )
              write(*,*)'rmin found',rmin,'rmax found',rmax


           ELSEIF (NUMAEROS.EQ.14) THEN      !Volcanic ash 2
              modanum=1
              do i=1,4
                gammodpar(i)=gamma_mod(i,2)
              enddo
c              write (6,910) A,alpha,b,gamma
              WRITE(6,910)(gammodpar(i), i=1,4)
              call rmm (cutoff, rmin, rmax )
              write(*,*)'rmin found',rmin,'rmax found',rmax


           ELSEIF (NUMAEROS.EQ.44) THEN      !Advection Fog 1
              modanum=1
              do i=1,4
                gammodpar(i)=gamma_mod(i,3)
              enddo
c              write (6,910) A,alpha,b,gamma
              WRITE(6,910)(gammodpar(i), i=1,4)
              call rmm (cutoff, rmin, rmax )
              write(*,*)'rmin found',rmin,'rmax found',rmax

           ELSEIF (NUMAEROS.EQ.45) THEN      !Advection Fog 2
              modanum=1
              do i=1,4
                gammodpar(i)=gamma_mod(i,4)
              enddo
c              write (6,910) A,alpha,b,gamma
              WRITE(6,910)(gammodpar(i), i=1,4)
              call rmm (cutoff, rmin, rmax )
              write(*,*)'rmin found',rmin,'rmax found',rmax

           ELSEIF (NUMAEROS.EQ.46) THEN      !Radiation Fog 3
              modanum=1
              do i=1,4
                gammodpar(i)=gamma_mod(i,5)
              enddo
c              write (6,910) A,alpha,b,gamma
              WRITE(6,910)(gammodpar(i), i=1,4)
              call rmm (cutoff, rmin, rmax )
              write(*,*)'rmin found',rmin,'rmax found',rmax

           ELSEIF (NUMAEROS.EQ.47) THEN      !Radiation Fog 4
              modanum=1
              do i=1,4
                gammodpar(i)=gamma_mod(i,6)
              enddo
c              write (6,910) A,alpha,b,gamma
              WRITE(6,910)(gammodpar(i), i=1,4)
              call rmm (cutoff, rmin, rmax )
              write(*,*)'rmin found',rmin,'rmax found',rmax

           ENDIF
C Truncated power law size distribution
C
C        ELSEIF (NEQ.EQ.3) THEN
C            CONST=0.0
C            RPMAX=0.0
C            D=0.0
C            ALPHA=0.0
C
C Log-normal size distr (NEQ=1)
        ELSE
C mixing esterno
C
            IF ( MIXTYPE.EQ.2 ) THEN
              modanum=1
              lognorpar(1,1)=HUM_RADII(NRH,COMP_MAT(NELEM,NUMAEROS))
              lognorpar(2,1)=SIGMA(COMP_MAT(NELEM,NUMAEROS))
              lognorpar(3,1)=1.0
              call rmm ( cutoff, rmin, rmax )
              write(*,*)'rmin found',rmin,'rmax found',rmax

            ELSE
              modanum=vec_comp(numaeros)
              do ip=1,modanum
               lognorpar(1,ip)=HUM_RADII(NRH,COMP_MAT(ip,NUMAEROS))
               lognorpar(2,ip)=SIGMA(COMP_MAT(ip,NUMAEROS))
               lognorpar(3,ip)=MIX_RATIO(ip,numaeros)
              enddo
              call rmm ( cutoff, rmin, rmax )
              write(*,*)'rmin found',rmin,'rmax found',rmax
            ENDIF


C            ANUM2=MIX_RATIO(2,NUMAEROS)
C           ENDIF
C            R1=HUM_RADII(NRH,COMP_MAT(1,NUMAEROS))
C            SIG1=SIGMA(COMP_MAT(1,NUMAEROS))
C         IF((NUMAEROS.EQ.8).OR.(NUMAEROS.EQ.9).OR.(NUMAEROS.EQ.18)) THEN
C              R2=R1
C              SIG2=SIG1
C         ELSE
C              R2=HUM_RADII(NRH,COMP_MAT(2,NUMAEROS))
C              SIG2=SIGMA(COMP_MAT(2,NUMAEROS))
C         ENDIF

C

C in una seconda fase scommenta !
C                 IF (R1.LE.0.0.OR.R2.LE.0.0) THEN
C                     WRITE (6,904)
C                     GO TO 70
C                 ENDIF
C            WRITE(6,911)R1,SIG1
C 911     FORMAT (/,'VALUES OF SIZE DISTRIBUTION PARAMETERS:',/,
C     + 'MODAL RADIUS=',F8.4,/,'STANDARD DEVIATION=',F8.4)

C check wether the aerosol is dry or humidity indipendet
C                  scommentare opportunamente
C            IF (NRH.NE.1) THEN
C             IF ((R1.EQ.99.).OR.(R2.EQ.99.).OR.(R3.EQ.99.)) THEN
C                     WRITE(*,97)NUMAEROS
C 97    format (' ERROR: For aerosol number',1x,i2,1x,'humidity level
C     + should be 1',/,'       resubmit a correct input file')
C        STOP
C       ELSEIF ((NUMAEROS.GE.6.AND.NUMAEROS.LE.9).OR.
C     +             (NUMAEROS.GE.13.AND.NUMAEROS.LE.14)) THEN
C       WRITE(*,98)NUMAEROS
C 98    FORMAT (' WARNING: Aerosol number',1X,I2,1X,'is humidity
C     + indipendent')
C       ENDIF
C       ENDIF
C end of the check
C       WRITE(*,*)neq, anum1, anum2, anum3,r1, r2, r3, sig1,
C     + sig2, sig3
C


C      ELSEIF (NEQ.EQ.5) THEN
C        RAIN=0.0
C      ELSE
C         ZNUM = 0.0
C         VLAM = 0.0

      ENDIF

        RI(1) = RMIN
        DELR = 2.*(RMAX-RMIN)/(NSTEP*(NSTEP+1.))
        do  I = 2, NRS
           RI(I) = RI(I-1)+FLOAT(I-1)*DELR
        enddo


C
C     CONTROL POINT TO CALCULATE SIZE DISTRIBUTION AT RI(I)
C
      ENTRY SIZ2
C
       IF (NEQ.EQ.2) THEN
          DO 20 I = 1, NRS
             ENR(I) = gammodpar(1)*RI(I)**gammodpar(2)*
     &               EXP(-gammodpar(3)*RI(I)**gammodpar(4))
          write(6,*) 'i,enr(i)',i,enr(i)
   20    CONTINUE
C      ELSEIF (NEQ.EQ.2) THEN
C         DO 30 I = 1, NRS
C            IF (RI(I).GT.RPMAX) THEN
C               ENR(I) = D*RI(I)**ALPHA
C            ELSE
C               ENR(I) = CONST
C            ENDIF
C   30    CONTINUE
       ELSE
         DO 40 I = 1, NRS
            IF (RI(I).LE.0.0) THEN
               ENR(I) = 0.
               GO TO 40
            ENDIF
C            write(*,*)anum1,sig1

             do ip=1,modanum
              bi(ip)=lognorpar(3,ip)/(RI(i)*lognorpar(2,ip)*SQRPI2)
              top(ip) = -(ALOG10(RI(I)/lognorpar(1,ip)))**2
              enr(i) = enr(i)+bi(ip)*EXP(top(ip)/
     *                (2.0*lognorpar(2,ip)**2))
             enddo
   40    CONTINUE

       ENDIF
C
      RETURN
C
C   70 CALL EXIT
C
C
C
C
C
C
C
C
C
   10 FORMAT (E12.5,E12.5)
  900 FORMAT ('  MODIFIED GAMMA DISTRIBUTION ')
  901 FORMAT ('  TRUNCATED POWER LAW DISTRIBUTION')
  902 FORMAT ('  LOG NORMAL DISTRIBUTION')
  903 FORMAT (5E15.8)
  904 FORMAT ('  RADIUS IS ZERO OR NEGATIVE ')
  905 FORMAT (1P9E10.3)
  906 FORMAT (/)
  907 FORMAT ('  MARSHALL-PALMER RAINDROP DISTRIBUTION')
  908 FORMAT ('  GENERALIZED RAINDROP DISTRIBUTION')
  909 FORMAT (3E15.8)
  910 FORMAT (/,'MODIFIED GAMMA PARAMETERS:',1X,F10.5,1X,F10.5,1X,
     * F10.5,1X,F10.5,/ )
C
      END
