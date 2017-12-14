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
      DOUBLE PRECISION z1,z2,y1,y2
cc    real*16          z1,z2,y1,y2
      PROD = 0.
      IF (MODE.GE.1) DM(1) = 0.
      NRS1 = NRS-1
      IF (R(1).EQ.0.) Y2 = 0.
      IF (R(1).GT.0.) Y2 = ENR(1)*R(1)**MOM
C
      DO 10 I = 1, NRS1
         Y1 = Y2
         Y2 = ENR(I+1)*R(I+1)**MOM
cc       print*,' y1,y2 ',y1,y2
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
C           zrat = z2/z1
            ZRAT=(Y2/Y1)*(R(I+1)/R(I))
cc          print*,' z2 z1 r ',zrat,r(i+1),r(i)
            ALPH1 = ALOG(zrat)/ALOG(R(I+1)/R(I))
            IF (MODE.GE.1) THEN
c              Z = (Z1/R(I)**ALPH1)*((R(I+1)+R(I))/2.)**ALPH1
                 if(alph1.gt.-15.)then
cc                   print*,' z1 r alph1 ',z1,r(i),alph1
                     Za = (Z1/R(I)**ALPH1)
                     zb = ((R(I+1)+R(I))/2.)**ALPH1
cc                   print*,' za,zb',za,zb
c
cc                   zb1 = ((R(I+1)+R(I))/2.)
cc                    zp1 = (zb1/r(i))**alph1
cc                   zp = zp1 * z1
                     z = za * zb
cc                   print*,' z,zp ',z,zp,zp1
               else
                     zb1 = ((R(I+1)+R(I))/2.)
                     zp1 = (zb1/r(i))**alph1
cc                   print*,' zp1 zb1 r ',zp1,zb1,r(i)
                     zp = zp1 * z1
                     z = zp
               endif
               DM(I) = DM(I)+(Z-Z1)/ALPH1
               DM(I+1) = (Z2-Z)/ALPH1
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
      SUBROUTINE NEWRAD (RI,ENR,NRS,MOM)
      DIMENSION RI(800),ENR(800),AREA(800),RIP(800)
C
C     GENERATE A NEW SET OF RADII UNIFORMLY SPACED IN N(R)**MOM
C     NOTE "AREA" REPRESENTS THE APPROPRIATE MOMENT OF N(R)
C
      DATA PII / 3.14159265 /
      NRS1 = NRS-1
      ARETOT = PII*PROD(RI,ENR,NRS,AREA,MOM,2)
      DO 10 I = 1, NRS
         AREA(I) = AREA(I)/AREA(NRS)
   10 CONTINUE
      NDEL = MIN1(SQRT(FLOAT(NRS)),15.)
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
