      MODULE alloc
        ALLOCATABLE UF11(:,:,:,:), UF12(:,:,:,:)
     &             ,UF22(:,:,:,:), UF33(:,:,:,:)
     &             ,UF34(:,:,:,:), UF44(:,:,:,:)
     &             ,UFEA(:,:,:,:)
      END MODULE alloc

      MODULE alloc1
        ALLOCATABLE U11(:,:,:,:,:), U12(:,:,:,:,:)
     &             ,U22(:,:,:,:,:), U33(:,:,:,:,:)
     &             ,U34(:,:,:,:,:), U44(:,:,:,:,:)
     &             ,UEA(:,:,:,:,:)
      END MODULE alloc1

      PROGRAM OPTCHAR1
c ** ANGLE & PO Sline interpolation; n, k, R Linear interpolation
c ** 11/29/05 Size Distribution option: SD table or
c ** LogNormal parameters (only KN=-KN (logarithmic intervals) &
c ** IA=1 - trapezoidal approximation case)
c ** 10/09/03 don't calculate f34 and f44 if key_f344=1
c ** 08/27/03 smooth f33 and f44 for 40< angle <60 
c ** 05/05/03 this version can be used to retrieve an axis ratio or
c **          axis ratio distribution
c **************************************************************** c
c **   02/28/03                                                 ** c
c **   Program calculates optical characteristics for given     ** c
c **   size distribution, refractive index, axis ratio          ** c
c **   distribution and wavelength                              ** c
c **************************************************************** c
c **                                                            ** c
c ** INPUT:                                                     ** c
c **                                                            ** c 
c **   key  = 1 - create fixed kernels (for fixed axis          ** c 
c **              ratio distr.) and save them                   ** c
c **              into 'Rke...fix...' files and calculate       ** c
c **              opt.characteristics                           ** c
c **          2 - read fixed kernels from 'Rke...fix...'  files ** c
c **          3 - create fixed kernels but don't save them      ** c
c **          4 - don't create fixed kernels, calculate         ** c
c **              opt.characteristics from original kernels     ** c
c **   key_RD =1 - volume mixture of spheroids                  ** c
c **           2 - surface area  mixture of spheroids           ** c
c **   key_f11= 0 - calculate scattering matrix                 ** c
c **            1 - calculate only phase function               ** c
c **   key_f344=0 - calculate scattering matrix                 ** c
c **            1 - don't calculate f34 and f44                 ** c
c **   key_org=0 - read original kernels simultaniously make    ** c 
c **               angle interpolation and calculate opt.char.  ** c
c **           1 -  -"-, save new kernels in                    ** c
c **                /distname_N/ directory, STOP                ** c
c **   key_fx works when key=1                                  ** c
c **        =0 -  create fixed kernels (for fixed axis          ** c 
c **              ratio distr.) and save them                   ** c
c **              into 'Rke...fix...' files and calculate       ** c
c **              opt.characteristics                           ** c
c **         1 - save fixed kernels with original kernel format ** c
c **             in order to be used as input kernels;          ** c
c **             'Rke...fix' kernels have to be renamed and moved* c
c **             into directory 'dir_name'(see 'matrix_fixget.f')* c
c **             The files can not be used if key=2.            ** c
c **                                                            ** c
c **   WL   - wavelength                                        ** c
c **   RN   - real part of the refractive index                 ** c
c **   RK   - imaginary part of the refractive index            ** c
c **   rgmin,rgmax &                                            ** c 
c **   wlmin,wlmax - min,max radii and wlmin,wlmax wavelengths  ** c
c **                 that are used to recalculate grid radii for** c
c **                 fixed kernels. New input file              ** c
c **                'grid1.dat.new' will be created if key=1    ** c
c **                 or key_org=1. Use key_grid1 to choose      ** c
c **                 'grid1.dat' or 'grid1.dat.new' will be read** c
c **                 for further calculations                   ** c  
c **   key_SD=0 - read Size Distribution table dV/dlnR          ** c
c **         =1 - calculate Size Distribution for grid radii    ** c
c **              using Log Normal function                     ** c
c **   ID    - dimension of d(...)/dlnR or d(...)/dR            ** c
c **       = 0 - number                                         ** c
c **       = 1 - radius                                         ** c
c **       = 2 - area                                           ** c
c **       = 3 - volume                                         ** c
c **   NMD   - number of modes (up to 2)                        ** c
c **   KN   - number of grid radii                              ** c
c **   grid(KN) - grid radii                                    ** c
c **   SD(KN)   - size distribution for grid radii              ** c
c **   (CM(i),SM(i),RMM(i),i=1,NMD) - size distribution         ** c
c **   function (LogNormal) parameters:                         ** c
c **                         CM - concentration                 ** c
c **                         SM - standard deviation            ** c
c **                         RMM - median radius                ** c
c **   distname_O - original kernel directory name              ** c
c **   distname_F - .fix kernel directory name                  ** c
c **   distname_N - new original kernel directory               ** c
c **                                      name (key_org=1)      ** c
c **   KR  - number of axis ratios                              ** c
c **   R(KR)  - grid axis ratios                                ** c
c **   RD(KR) - axis ratio distribution for grid axis ratios    ** c
c **   KM   - number of scattering angles                       ** c
c **   ANGLE(KM) - scattering angles                            ** c
c **                                                            ** c
c ** OUTPUT:                                                    ** c
c **                                                            ** c
c **   ext     - extinction                                     ** c
c **   albedo  - albedo                                         ** c
c **   f... - scattering matrix elements                        ** c
c **************************************************************** c
      use alloc
      use alloc1
      include 'optchar.par'
      integer KN, KM, KR, NDP,
     &       key, key_RD, key_f11, key_f344, key_org, KSIMP
      real WL, RN, RK
      dimension grid(KNpar), SD(KNpar)
     &            , RD(KRpar), ANGLE(KMpar)
      real*4 R(KRpar)
      dimension f11(KMpar) 
     &         ,f12(KMpar)
     &         ,f22(KMpar)
     &         ,f33(KMpar)
     &         ,f34(KMpar)
     &         ,f44(KMpar)
      real ext, albedo
	dimension CM(KMD),SM(KMD),RMM(KMD),
     &          RRR(KNpar),AR(KNpar)
      CHARACTER(60) distname_O,distname_F,distname_N
c ----- Timer ------
      real*4 T_TTL
      real*4 tarray(2)
      real*4 external dtime
      real*4 external etime
      
	NDP=0
c
c ** READ INPUT
c       
      open(10,file='dls_input.dat',status='old')
      read(10,*) key,key_RD,key_f11,key_f344,
     &           key_org,key_fx
      write(*,*) key,key_RD,key_f11,key_f344,
     &           key_org,key_fx
      read(10,*) WL,RN,RK,rgmin,rgmax,wlmin,wlmax
	write(*,*) WL,RN,RK,rgmin,rgmax,wlmin,wlmax
	read(10,*) key_SD,ID,NMD
	if(key_SD.eq.0) then
      read(10,*) 
      read(10,*) KN
	if(KN.gt.KNpar) then
      write(*,*) 'in input.dat KN=',KN,' .gt. ',
     & 'KNpar=',KNpar,' in optchar.par'
	STOP 'STOP: KN should be < or = KNpar'
	endif ! KN&KNpar
      write(*,*) 'KN=',KN
      do i=1,KN
      read(10,*) grid(i), SD(i)
      enddo ! i
	else ! key_SD=1
	read(10,*) (CM(i),SM(i),RMM(i),i=1,NMD)
	read(10,*) KN
      write(*,*) 'KN=',KN
      do i=1,KN
      read(10,*) grid(i)
      enddo ! i
	endif ! key_SD
      write(*,*) 'after grid'
      read(10,*) distname_O
	write(*,*) 'distname_O=',TRIM(distname_O)
	read(10,*) distname_F
	write(*,*) 'distname_F=',TRIM(distname_F)
	read(10,*) distname_N
	write(*,*) 'distname_N=',TRIM(distname_N)
      read(10,*) KR
	if(KR.gt.KRpar) then
      write(*,*) 'in input.dat KR=',KR,' .gt. ',
     & 'KRpar=',KRpar,' in optchar.par'
	STOP 'STOP: KR should be < or = KRpar'
	endif ! KR&KRpar
      do i=1,KR
      read(10,*) R(i), RD(i)
      enddo ! i
      read(10,*) KM
	if(KM.gt.KMpar) then
      write(*,*) 'in input.dat KM=',KM,' .gt. ',
     & 'Kmpar=',KMpar,' in optchar.par'
	STOP 'STOP: KM should be < or = KMpar'
	endif ! KM&KMpar
      do j=1,KM
      read(10,*) ANGLE(j)
      enddo ! j
      close(10)
      if(key.eq.4.and.key_org.eq.1) then
      write(*,*) ' if key=4, key_org=',key_org,
     & ' should be 0.'
	stop 'STOP in OPTCHAR1 (matrix_optchr.f)'
	endif 
	if(key_org.eq.1.and.key_fx.eq.1) then
      write(*,*) 'STOP: key_org=1 & key_fx=1'
      write(*,*) 'If you want to use key_fx=1',
     &           ' change key_org to 0'
      write(*,*) 'If key_org=0 is your choice,',
     &           ' check dir_name1 in matrix_fixget',
     &           ' and run the code again.'
	stop 
	endif

c **   key_grid1 read grid radii and scat.angles which were used** c
c **             for kernel look up tables or fixed kernels     ** c
c **          =0 - 'grid1.dat'                                  ** c
c **           1 - 'grid1.dat.fix'                              ** c

      if(key.eq.2) then
	 key_grid1=1
	else 
	 key_grid1=0
	endif

c
c ** CALCULATE SD if key_SD=1
c      
      if(key_SD.eq.1) then
      call SIZEDISDN(-KN,1,ID,1,NMD,CM(1:NMD),SM(1:NMD),RMM(1:NMD),
     &                      grid(1),grid(KN),RRR,AR,AC)
      write(*,*) 'CM=',CM(1:NMD),' AC=',AC

      do i=1,KN,10
	write(*,13) i,grid(i),RRR(i),SD(i),AR(i) 
      enddo ! i
      do i=1,KN
      SD(i)=AR(i)
      enddo ! i
	endif ! key_SD
c
c ** CHECK INPUT
c
      write(*,*) 'CHECK INPUT in OPTCHAR1:'
      write(*,*) 'key, key_RD, key_RD, key_f11',
     &   ' key_f344, key_org, key_fx, key_grid1'
      write(*,'(i3,7i6)') key, key_RD, key_RD,
     &   key_f11, key_f344, key_org, key_fx, key_grid1
      write(*,*) 'WL,RN,RK,rgmin,rgmax,wlmin,wlmax:'
      write(*,'(f6.4,4e12.5,2f7.4)')  WL,RN,RK,rgmin,rgmax,wlmin,wlmax
      write(*,*) 'size distribution:'
	if(key_SD.eq.0) then
       do i=1,KN,10
       write(*,13) i,grid(i), SD(i)
       enddo ! i
	else 
       do i=1,KN,10
       write(*,13) i,grid(i),RRR(i),SD(i)
      enddo ! i
	endif 
      write(*,*) 'axis ratio distribution:'
      do i=1,KR
      write(*,13) i,R(i), RD(i)
      enddo ! i
      write(*,*) 'SCATTERING ANGLES: KM=',KM
      write(*,14) ANGLE(:KM)
      if(key_f11.eq.0) then
      write(*,*) 'Scattering matrix will be calculated'
      else      
      write(*,*) 'Phase function will be calculated'
      endif ! key_f11
       if(key.eq.1.and.key_fx.eq.1) 
     &  write(*,*) 'Fixed kernels will be used as original ones !'
c
c ** rgrid min & max calculation for fixed or NEWORG kernels
c
      if(key.eq.1.or.key_org.eq.1) then
      pi2=2*ACOS(-1.)
	pomin=pi2*rgmin/wlmax
	pomax=pi2*rgmax/wlmin  
	endif ! key or key_org
	         
c
c ** GET OPTICAL CARACTERISTICS
c
c *** ALLOCATE ARRAYS to be used in subroutine USMATRIX 
c *** (in matrix_intrpl.f)
c
      IF(key.lt.4) THEN 
        ALLOCATE(UF11(KMpar,KN1par,KIMpar,KREpar),stat=ierr)
        if(ierr/=0) stop 'Can not allocate UF11 array'
        ALLOCATE(UFEA(2,KN1par,KIMpar,KREpar),    stat=ierr)
        if(ierr/=0) stop 'Can not allocate UFEA array'
        UF11=0.
        UFEA=0.
        if(key_f11.eq.0) then
          ALLOCATE(UF12(KMpar,KN1par,KIMpar,KREpar),stat=ierr)
          if(ierr/=0) stop 'Can not allocate UF12 array'
          ALLOCATE(UF22(KMpar,KN1par,KIMpar,KREpar),stat=ierr)
          if(ierr/=0) stop 'Can not allocate UF22 array'
          ALLOCATE(UF33(KMpar,KN1par,KIMpar,KREpar),stat=ierr)
          if(ierr/=0) stop 'Can not allocate UF33 array'
          UF12=0.
          UF22=0.
          UF33=0.
	    if(key_f344.ne.1) then
          ALLOCATE(UF34(KMpar,KN1par,KIMpar,KREpar),stat=ierr)
          if(ierr/=0) stop 'Can not allocate UF34 array'
          ALLOCATE(UF44(KMpar,KN1par,KIMpar,KREpar),stat=ierr)
          if(ierr/=0) stop 'Can not allocate UF44 array'
          UF34=0.
          UF44=0.
	    endif ! key_f344
        endif ! key_f11
      ENDIF ! key<4

      IF(key.gt.2) THEN 
c
c *** ALLOCATE ARRAYS (alloc1.mod)
c
	  ALLOCATE(U11(KR1par,KMpar,KN1par,KIMpar,KREpar),stat=ierr)
        if(ierr/=0) stop 'Can not allocate U11 array'
        ALLOCATE(UEA(KR1par,2,KN1par,KIMpar,KREpar),    stat=ierr)
        if(ierr/=0) stop 'Can not allocate UEA array'
        U11=0. 
        UEA=0.
        if(key_f11.eq.0) then
          ALLOCATE(U12(KR1par,KMpar,KN1par,KIMpar,KREpar),stat=ierr)
          if(ierr/=0) stop 'Can not allocate U12 array'
          ALLOCATE(U22(KR1par,KMpar,KN1par,KIMpar,KREpar),stat=ierr)
          if(ierr/=0) stop 'Can not allocate U22 array'
		ALLOCATE(U33(KR1par,KMpar,KN1par,KIMpar,KREpar),stat=ierr)
          if(ierr/=0) stop 'Can not allocate U33 array'
          U12=0.
          U22=0.
          U33=0.
		if(key_f344.ne.1) then
          ALLOCATE(U34(KR1par,KMpar,KN1par,KIMpar,KREpar),stat=ierr)
          if(ierr/=0) stop 'Can not allocate U34 array'
          ALLOCATE(U44(KR1par,KMpar,KN1par,KIMpar,KREpar),stat=ierr)
          if(ierr/=0) stop 'Can not allocate U44 array'
          U34=0.
          U44=0.
	    endif ! key_f344
        endif ! key_f11
      ENDIF ! key.gt.2    

c      DO i=1,2  ! loop for test

      CALL OPTCHAR(key,key_RD,key_f11,key_f344,
     &                     key_org,key_fx,key_grid1,
     &                     WL,RN,RK,KN,grid,SD,
     &                     KR,R,RD,KM,ANGLE,ext,albedo,
     &                     f11,f12,f22,f33,f34,f44,pomin,pomax
     &                    ,distname_O,distname_F,distname_N,NDP)
c      write(*,*) 'LOOP TEST:   after  OPTCHAR  i=',i
c      WL=0.64 ! delete after test
c      RN=1.34 ! delete after test
c      RK=0.01 ! delete after test
c      ENDDO ! i

13    FORMAT(i4,4E14.6)
14    FORMAT(7F8.2)
c
c *** DEALLOCATE ARRAYS (alloc.mod) in subroutine USMATRIX 
c *** (in matrix_intrpl.f)
c
      IF(key.ne.4) THEN
        DEALLOCATE(UF11,stat=ierr)
        if(ierr/=0) stop 'Can not deallocate UF11 array'
        DEALLOCATE(UFEA,stat=ierr)
        if(ierr/=0) stop 'Can not deallocate UFEA array'
        if(key_f11.eq.0) then
          DEALLOCATE(UF12,stat=ierr)
          if(ierr/=0) stop 'Can not deallocate UF12 array'
          DEALLOCATE(UF22,stat=ierr)
          if(ierr/=0) stop 'Can not deallocate UF22 array'
          DEALLOCATE(UF33,stat=ierr)
          if(ierr/=0) stop 'Can not deallocate UF33 array'
	     if(key_f344.ne.1) then
          DEALLOCATE(UF34,stat=ierr)
          if(ierr/=0) stop 'Can not deallocate UF34 array'
          DEALLOCATE(UF44,stat=ierr)
          if(ierr/=0) stop 'Can not deallocate UF44 array'
	     endif ! key_f344
        endif ! key_f11
      ENDIF ! key.ne.4
      IF(key.gt.2) THEN 
c
c *** DEALLOCATE ARRAYS (alloc1.mod)
c
      DEALLOCATE(U11,stat=ierr)
      if(ierr/=0) stop 'Can not deallocate U11 array'
      DEALLOCATE(UEA,stat=ierr)
      if(ierr/=0) stop 'Can not deallocate UEA array'
      if(key_f11.eq.0) then
        DEALLOCATE(U12,stat=ierr)
        if(ierr/=0) stop 'Can not deallocate U12 array'
        DEALLOCATE(U22,stat=ierr)
        if(ierr/=0) stop 'Can not deallocate U22 array'
        DEALLOCATE(U33,stat=ierr)
        if(ierr/=0) stop 'Can not deallocate U33 array'
	   if(key_f344.ne.1) then
        DEALLOCATE(U34,stat=ierr)
        if(ierr/=0) stop 'Can not deallocate U34 array'
        DEALLOCATE(U44,stat=ierr)
        if(ierr/=0) stop 'Can not deallocate U44 array'
	   endif ! key_f344
      endif ! key_f11 
      ENDIF ! key>2    
                                  T_TTL=etime(tarray)
c      write(*,61) T_TTL/60.
      T_TTL=tarray(1)+tarray(2)  
      write(*,61) T_TTL/60.
      write(*,*) '-------------------------------------------------' 
61    format(' Total execution time ..... ',f8.3,' min.')
      STOP
      END PROGRAM OPTCHAR1

c ****************************************************************

      SUBROUTINE OPTCHAR(key,key_RD,key_f11,key_f344,
     &                       key_org,key_fx,key_grid1,
     &                       WL,RN,RK,KN,grid,SD,
     &                       KR,R,RD,KM,ANGLE,ext,albedo,
     &                       f11,f12,f22,f33,f34,f44,pomin,pomax
     &                      ,distname_O,distname_F,distname_N, NDP)

c **************************************************************** c
c **   08/26/02                                                 ** c
c **   Subroutine calculates optical characteristics for given  ** c
c **   size distribution, refractive index, axis ratio          ** c
c **   distribution and wavelength.                             ** c
c **                                                            ** c
c **   In case RN or RK or R() is out of correspoding ranges    ** c 
c **   subroutine changes RN or RK or R() for edge value and    ** c
c **   gives optical characteristics for new values             ** c 
c **************************************************************** c
c **                                                            ** c
c ** INPUT:                                                     ** c
c **                                                            ** c 
c **   key  = 1 - create fixed kernels (for fixed axis          ** c 
c **              ratio distr.) and save them                   ** c
c **              into 'Rke...fix...' files and calculate       ** c
c **              opt.characteristics                           ** c
c **          2 - read fixed kernels from 'Rke...fix...'  files ** c
c **          3 - create fixed kernels but don't save them      ** c
c **          4 - don't create fixed kernels, calculate         ** c
c **              opt.characteristics from original kernels     ** c
c **   key_RD =1 - volume mixture of spheroids                  ** c
c **           2 - surface area  mixture of spheroids           ** c
c **   key_f11=0 - calculate scattering matrix                  ** c
c **           1 - calculate only phase function                ** c
c **   key_f344=0 - calculate scattering matrix                 ** c
c **            1 - don't calculate f34 and f44                 ** c
c **   key_org=0 - read original kernels simultaniously make    ** c 
c **               angle interpolation and calculate opt.char.  ** c
c **           1 -  -"-, save new kernels in                    ** c
c **                /NEWORG/ directory, STOP                    ** c
c **        =0 -  create fixed kernels (for fixed axis          ** c 
c **              ratio distr.) and save them                   ** c
c **              into 'Rke...fix...' files and calculate       ** c
c **              opt.characteristics                           ** c
c **         1 - save fixed kernels with original kernel format ** c
c **             in order to be used as input kernels;          ** c
c **             'Rke...fix' kernels have to be renamed and moved* c
c **             into directory 'dir_name'(see 'matrix_fixget.f')* c
c **             The files can not be used if key=2.            ** c
c **                                                            ** c
c **   key_grid1 read grid radii and scat.angles which were used** c
c **             for kernel look up tables or fixed kernels     ** c
c **          =0 - 'grid1.dat'                                  ** c
c **           1 - 'grid1.dat.new'                              ** c
c **   WL   - wavelength                                        ** c
c **   RN   - real part of the refractive index                 ** c
c **   RK   - imaginary part of the refractive index            ** c
c **   KN   - number of grid radii                              ** c
c **   grid(KN) - grid radii                                    ** c
c **   SD(KN)   - size distribution for grid radii              ** c
c **   distname_O - original kernel directory name              ** c
c **   distname_F - .fix kernel directory name                  ** c
c **   distname_N - new original kernel directory name          ** c
c **   KR  - number of axis ratios                              ** c
c **   R(KR)  - grid axis ratios                                ** c
c **   RD(KR) - axis ratio distribution for grid axis ratios    ** c
c **   KM   - number of scattering angles                       ** c
c **   ANGLE(KM) - scattering angles                            ** c
c **                                                            ** c
c ** OUTPUT:                                                    ** c
c **                                                            ** c
c **   ext     - extinction                                     ** c
c **   albedo  - absorption                                     ** c
c **   f... - scattering matrix elements                        ** c
c **************************************************************** c
      include 'optchar.par'
      integer KN, KM, KR, NDP,
     &       key, key_RD, key_f11, key_f344, key_org
      real WL, RN, RK, dlnr, xnorm
      dimension grid(KNpar), SD(KNpar)
     &            , RD(KRpar), ANGLE(KMpar)
      real*4 R(KRpar)
      COMMON /US1/ US11(KMpar,KNpar) 
      COMMON /US2/ US12(KMpar,KNpar)
      COMMON /US3/ US22(KMpar,KNpar)
      COMMON /US4/ US33(KMpar,KNpar)
      COMMON /US5/ US34(KMpar,KNpar)
      COMMON /US6/ US44(KMpar,KNpar)
      COMMON /US0/ USEA(2,KNpar)
      dimension f11(KMpar) 
     &         ,f12(KMpar)
     &         ,f22(KMpar)
     &         ,f33(KMpar)
     &         ,f34(KMpar)
     &         ,f44(KMpar)
      real ext, abs, sca, albedo
      dimension X(2), Y(2) 
      real LINEAR
      CHARACTER(60) distname_O,distname_F,distname_N
c
c ** GET MATRICES
c
c      write(*,*) 'before USMATRIX'
c      do itest=1,2  ! delete after test
      CALL USMATRIX(key,key_RD,key_f11,key_f344,
     &                  key_org,key_fx,key_grid1,
     &                  WL,RN,RK,KN,grid,
     &                  KR,R,RD,KM,ANGLE,dlnr,pomin,pomax
     &                 ,distname_O,distname_F,distname_N,NDP)
cl      write(*,*) 'after USMATRIX'

cl         write(*,*) 'in OPTCHAR11: US11=',US11(1,1),
cl     &    ' US11=',US11(KM,KN)
c
c ** CALCULATE APPROXIMATED OPTICAL CHARACTERISTICS      
c
cl      f11(:)=0.
cl      f12(:)=0.
cl      f22(:)=0.
cl      f33(:)=0.
cl      f34(:)=0.
cl      f44(:)=0.
cl      ext=0.
cl      abs=0.
c      coeff=1.
 
c	write(*,*) 'in OPTCHAR - SD1: dlnr=',dlnr
c	write(*,*) SD(1:KN)

      coeff=1e+3                  ! for AERONET dV/dlnR
      SD(1:KN)=SD(1:KN)*coeff*dlnr


c	write(*,*) 'in OPTCHAR - SD2:'
c	write(*,*) SD(1:KN)
     
      ext=DOT_PRODUCT(USEA(1,1:KN),SD(1:KN))
      abs=DOT_PRODUCT(USEA(2,1:KN),SD(1:KN))
      do j=1,KM
      f11(j)=DOT_PRODUCT(US11(j,1:KN),SD(1:KN))
      enddo ! j
      sca=ext-abs
      albedo=sca/ext
      f11(1:KM)=f11(1:KM)/sca
        
      if(key_f11.eq.0) then
      do j=1,KM
      f12(j)=DOT_PRODUCT(US12(j,1:KN),SD(1:KN))
      f22(j)=DOT_PRODUCT(US22(j,1:KN),SD(1:KN))
      f33(j)=DOT_PRODUCT(US33(j,1:KN),SD(1:KN))
      enddo ! j

      f12(1:KM)=-f12(1:KM)/sca/f11(1:KM)
      f22(1:KM)= f22(1:KM)/sca/f11(1:KM)
      f33(1:KM)= f33(1:KM)/sca/f11(1:KM)

cl      f12(1:KM)=-f12(1:KM)/sca
cl      f22(1:KM)= f22(1:KM)/sca
cl      f33(1:KM)= f33(1:KM)/sca      

	if(key_f344.ne.1) then
      do j=1,KM
      f34(j)=DOT_PRODUCT(US34(j,1:KN),SD(1:KN))
      f44(j)=DOT_PRODUCT(US44(j,1:KN),SD(1:KN))
	enddo ! j

      f34(1:KM)= f34(1:KM)/sca/f11(1:KM)
      f44(1:KM)= f44(1:KM)/sca/f11(1:KM)

cl      f34(1:KM)= f34(1:KM)/sca
cl      f44(1:KM)= f44(1:KM)/sca
	 endif ! key_f344
      endif ! key_f11

      SD(1:KN)=SD(1:KN)/coeff/dlnr
c
c ** SMOOTH f33 and f44 
c
      if(key_f11.eq.0) then
      j=1
      do while(ANGLE(j).lt.(40.).or.j.gt.KM)
	j1=j
	j=j+1
	enddo
      j=1
      do while(ANGLE(j).lt.(50.).or.j.gt.KM)
	j2=j
	j=j+1
	enddo

      if(j2.gt.j1) then
	X(1)=ANGLE(j1-1)
	X(2)=ANGLE(j2+1)
	Y(1)=f33(j1-1)
	Y(2)=f33(j2+1)
	Y(1)=LOG(Y(1))
	Y(2)=LOG(Y(2))
      do j=j1,j2
      f33(j)=EXP(LINEAR(X, Y, 2, ANGLE(j)))
	enddo ! j
      endif ! j2&j1
      endif ! key_f11
c
      if(key_f344.ne.1) then
      j=1
      do while(ANGLE(j).lt.(50.).or.j.gt.KM)
	j1=j
	j=j+1
	enddo
      j=1
      do while(ANGLE(j).lt.(60.).or.j.gt.KM)
	j2=j
	j=j+1
	enddo

      if(j2.gt.j1) then
	X(1)=ANGLE(j1-1)
	X(2)=ANGLE(j2+1)
	Y(1)=f44(j1-1)
	Y(2)=f44(j2+1)
	Y(1)=LOG(Y(1))
	Y(2)=LOG(Y(2))
      do j=j1,j2
      f44(j)=EXP(LINEAR(X, Y, 2, ANGLE(j)))
	enddo ! j
      endif ! j2&j1
	endif ! key_f344

c ** Check f11 norma
c
      if((KM.eq.180.or.KM.eq.181).and.
     &       (ANGLE(1).eq.(0.).and.ANGLE(KM).eq.(180.))) then
       call SINT  ( ANGLE, f11, KM, xnorm )
       call SINT  ( ANGLE, f11, KM, xnorm )
       call ASYPAR( ANGLE, f11, KM, g )
       endif ! KM 
c ** WRITE APPROXIMATED OPTICAL CHARACHTERISTICS
c
      open(10,file='sca_mtrx.dat',status='unknown')

      write(10,*)
c      write(10,*)'<<<     WELCOME TO NONSPHERICAL',
c     &      ' AEROSOL WORLD     >>>'
c      write(10,*)
      write(10,*)'Size distribution:'
      write(10,*) '      r(mkm)             Sd(r)'
      do i=1,KN
      write(10,'(2e25.16)') grid(i),SD(i)
      enddo ! i
      write(10,*)
      write(10,*)'Axis ratio distribution:'
      write(10,*) '         R                Rd(R)'
      do i=1,KR
      write(10,'(2e25.16)') R(i),RD(i)
      enddo ! i
      write(10,*) 
      write(10,10) WL, RN, RK
      write(10,*)
      if(key_RD.eq.1) then
      write(10,*) '        APPROXIMATED OPTICAL CHARACTERISTICS'
      write(10,*) '                 (volume mixture)'
      endif ! key_RD
      if(key_RD.eq.2) then 
      write(10,*) '        APPROXIMATED OPTICAL CHARACTERISTICS'
      write(10,*) '               (surface area mixture)'
      endif ! key_RD
      write(10,*) 
      write(10,11) ext, abs, sca, albedo
	write(10,*) 
      if(key_f11.eq.0) then
	if(key_f344.ne.1) then
      write(10,'(2x,''ANGLE'',7x,''F11'',8x,''-F12/F11'',7x,
     &''F22/F11'',7x,''F33/F11'',7x,''F34/F11'',7x,''F44/F11'')') 
      do j=1,KM       
      write(10,12) ANGLE(j),f11(j),f12(j),f22(j)
     &                     ,f33(j),f34(j),f44(j)
      enddo ! j
	else
      write(10,*) 
     &' ANGLE     F11         -F12/F11      F22/F11      F33/F11' 
      do j=1,KM       
      write(10,12) ANGLE(j),f11(j),f12(j),f22(j),f33(j)
      enddo ! j
	endif ! key_f344
      else
      write(10,*) 
     &' ANGLE       F11'
      do j=1,KM       
      write(10,12) ANGLE(j),f11(j)
      enddo ! j
      endif ! key_f11
cl	write(10,13) xnorm,1.
c
      if((KM.eq.180.or.KM.eq.181).and.
     &       (ANGLE(1).eq.(0.).and.ANGLE(KM).eq.(180.)))
     &  write(10,'(''asymmetry parameter='',e25.16)') g
      close(10)
13    format('check: f11 norma=',f7.3,'   one=',f7.3) 
10    FORMAT('wavelength=',e11.4,'  n=',e15.8,'  k=',e15.8) 
11    FORMAT('ext=',e25.16,'  abs=',e25.16,
     &       '  sca=',e25.16,'  albedo=',e25.16) 
12    FORMAT(F7.2,6E25.16)
c      enddo ! itest  - delete after test
      RETURN 
      END SUBROUTINE OPTCHAR

c **************************************************************** c
      REAL FUNCTION LINEAR( X, Y, M, X1 )
c ***************************************************************** c
c **       Linear interpolation function.                        ** c
c ***************************************************************** c
      IMPLICIT  REAL   (A-H, O-Z)
      REAL      X( * ), Y( * )
C
      ind=0
      IF ( X1.LT.X(1) )  THEN
         LINEAR = Y(1) - ( X1-X(1) )*( Y(2)-Y(1) )/( X(2)-X(1) )
      ELSE IF ( X1.GT.X(M) )  THEN
         LINEAR = Y(M) + ( X1-X(M) )*( Y(M)-Y(M-1) )/( X(M)-X(M-1) )
      ELSE
      DO N=2,M
       IF ( X1.GE.X(N-1).AND.X1.LE.X(N).AND.X(N-1).NE.X(N) ) THEN
       LINEAR = Y(N-1) + ( X1-X(N-1) )*( Y(N)-Y(N-1) )/
     $                                        ( X(N)-X(N-1) )
c         IF ( X1.GE.X(N-1) .AND. X1.LE.X(N) ) LINEAR = Y(N-1) +
c     $         ( X1-X(N-1) )*( Y(N)-Y(N-1) )/( X(N)-X(N-1) )
       ind=ind+1
       ENDIF
C*****************************
        IF(X1.EQ.X(N-1)) LINEAR =Y(N-1)
        IF(X1.EQ.X(N))   LINEAR =Y(N)
C*****************************

	if(ind.eq.1) goto 10
      ENDDO ! N
      ENDIF

10    RETURN
      END
c ***************************************************************** c
      subroutine SINT(X1, Y1, KM, xnorm)
c ** Simpson itegration
c ** xh - angle step (degree)
      dimension X1( KM ), Y1( KM ), Y(2)
      real LINEAR, xnorm
      KSIMP=721
      xh=180./float(KSIMP-1)
      KSIMP1=(180.-X1(1))/xh+1
c	write(*,*) 'KSIMP1=',KSIMP1
	pi=ACOS(-1.)
      Y(1)=Y1(1)
	Y(2)=Y1(KM)
      Y1(1:KM)=LOG(Y1(1:KM))
	      
	xnorm=0.
	XA=X1(1)
      do IK=1,KSIMP1
       II=IK/2
       IF(IK.EQ.1.OR.IK.EQ.KSIMP) THEN
	  if(IK.eq.1) then
	   F11=1./3.*Y(1)
	  else
         F11=1./3.*Y(2)
	  endif
	 GOTO 1
	 ENDIF
       IF(II*2.LT.IK) THEN
	 F11=2./3.*EXP(LINEAR(X1, Y1, KM, XA))
	 GOTO 1
	 ENDIF
       IF(II*2.EQ.IK) THEN
	 F11=4./3.*EXP(LINEAR(X1, Y1, KM, XA))
	 GOTO 1
	 ENDIF
1      CONTINUE
       xnorm=xnorm+F11*SIN(XA*pi/180.)
	 XA=XA+xh
	enddo ! IK
	xnorm=0.5*xnorm*xh*pi/180.
cl      write(*,'(''xnorm='',e13.5)') xnorm
	Y1(1:KM)=EXP(Y1(1:KM))/xnorm

	return
	end subroutine SINT
c
c ***************************************************************** c
      subroutine ASYPAR(X1, Y1, KM, g)
c ** Asymmetry parameter
c ** xh - angle step (degree)
      dimension X1( KM ), Y1( KM ), Y(2)
      real LINEAR, g, XX
      KSIMP=721
      xh=180./float(KSIMP-1)
	pi=ACOS(-1.)
      xpi=pi/180.
      Y(1)=Y1(1)
	Y(2)=Y1(KM)
      Y1(1:KM)=LOG(Y1(1:KM))	      
	g=0.
	XA=0.
      do IK=1,KSIMP
       II=IK/2
       IF(IK.EQ.1.OR.IK.EQ.KSIMP) THEN
	  if(IK.eq.1) then
	   F11=1./3.*Y(1)
	  else
         F11=1./3.*Y(2)
	  endif
	 GOTO 1
	 ENDIF
       IF(II*2.LT.IK) THEN
	 F11=2./3.*EXP(LINEAR(X1, Y1, KM, XA))
	 GOTO 1
	 ENDIF
       IF(II*2.EQ.IK) THEN
	 F11=4./3.*EXP(LINEAR(X1, Y1, KM, XA))
	 GOTO 1
	 ENDIF
1      CONTINUE
       XX=XA*xpi
       g=g+F11*SIN(XX)*COS(XX)
	 XA=XA+xh
	enddo ! IK
	g=0.5*g*(xh*xpi)

	Y1(1:KM)=EXP(Y1(1:KM))

	return
	end subroutine ASYPAR
c
