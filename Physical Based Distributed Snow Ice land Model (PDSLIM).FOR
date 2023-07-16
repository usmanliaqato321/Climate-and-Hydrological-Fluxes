$debug 
C
C   This PDSLIM_Karakoram.FOR code is a development of PDSLIMCOMPLETE2021_ETPP.FOR compiled 20210711 
C   which is a development of the PDSLIM code developed for the Mandrone Glacier 
C   and adapted to Naltar catchment in the Karakoram mountain range.
C   KARAKORAM 2.15
C   THE MASS BALANCE unit 22 MASSBAL.RST IS COMPUTED
C
C   KARAKORAM 2.14
C
C   READ THE DATE OF START OF RECORDING THE ANIMATION 
C
C   KARAKORAM 2.13 
C   Zw=6.5 m
C   Zt=6.0 m 
C   KARAKORAM_2.10
C 2.10   we save in Gv5d(NX,NY,10) either tice, Tsoil, PACT1
C
C	OPEN(70,FILE='RLOWXY.dat',STATUS='unknown',FORM='binary')
C     OPEN(85,FILE='HLATEN.dat',STATUS='unknown',FORM='binary')
C
C    
C   KARAKORAM_2.9
C 2.9   we add saving the field of surface temperature
C
C   KARAKORAM_2.8
C
C 2.8   we go back to etppm  with zt=6 m
C   KARAKORAM_2.7
C
C 2.7   we go back to etppm but with zt=2 m
C
C   KARAKORAM_2.6
C
C 2.6   write SWE and xymelt in binary format on Unit 27 and 30 respectively with wribin1
C       cancel ISTAMP -2 tesidtm tesispck 
C
C   KARAKORAM_2.5
C
C   write HCONVXY and SWAVXY in binary format wribin1
C
C   KARAKORAM_2.4
C
C   write HCONVXY and SWAVXY in ascii format wriascF7
C
C   PDSLIM KARAKORAM_2.3
C
C      dw=amin1(zw-1.,d)
C      log1=alog((zW-dw)/zom)
C      dt=amin1(zt-0.5,d)
C      log2=amax1(alog((zt-dt)/zov),0.)
C 
C	Gv5d(nx,ny,6)=SWAVBA on Unit 60
C	Gv5d(nx,ny,7)=RLOWAB
C	Gv5d(nx,ny,8)=HCONVE on Unit 80
C	Gv5d(nx,ny,9)=HLATEN
C  
C     PRINTING HCONVE IN PERC1.OUT      
C   PDSLIM KARAKORAM_2.2
C 
C   check
C
C   PDSLIM KARAKORAM_2.1
C   CHANGE OF SUBROUTINE ETPPM2
C   and DWind
C   PDSLIM_KARAKORAM_2.0
C   20211226 RR changed subroutine Dw(.) to take into account wind and temperature taken at height Zw and Zt respectively
C  
C
C
C   The program starts with a BLOCK DATA definition: all the variables,
C    included or not in a COMMON statement can be inizialized in this block
C    
C
C	Version updated by Gitti and Grossi on 20060124 changing in ETP the input Z0
C     20060201 by Ranzi and Gitti  
C	Zw =6 m
C     Zt =6 m for Mandrone
C     Zt =2 m for Naltar
C	xyicem  amax(0.,xyicem)
c     melt spatial average computed in statxy
c	Version with Albedo ice(Albedi) from ADAMELLO GLACIER 2003 08 23 ASTER IMAGE (about 0.37 see also Regine Hock, 1998)
C     26/06/07 
CRRAG	20070831 V01.02
C     landuse flag includes 14 vegetazione boschiva ed arbustiva in evoluzione
c              15 vegetazione rada 
c              treated as 13 (Alder) brughiera e cespuglieti   
CRRMUL WE HAVE TO CHECK WHY IT STOPS ON 29 SEPTEMBER 2012 PROBABLY BECAUSE OF BISEXTILE YEARS
CRRMUL
      BLOCK DATA
      PARAMETER(NNR=353,NNC=262)
	REAL RatiKh,CKa,RatiKe,VLh
      REAL CJOMME,VVCFCT(12)
	REAL KATM
      COMMON /RECORD/ XX(10,6,8784),NTIMST
      COMMON /PARAME/ XS(20),X(20),INDX(20),N,NPAR,STEPS(20),TESTPS(20),
     +TESSTP(20),STEP(20)
      COMMON /HEAT/CPw,Wden,ROhu,Cpa,RatiKh,CKa,Z0,VLh,RatiKe,Rgas,ZW,ZT
      COMMON /SOIL/ ASOIL,SLAMDA,Tmsoil,Cv(NNR,NNC),DSOIL1(NNR,NNC)
      COMMON /SNOW/A0,A1,TSPRIN,ASPRIN,TRNCF,AM0,CAM,GANESN,ALFAS,
     +BETAS,STRN,CJOMME,HSNLA1,FHI,Swi,FHIeKs,BEERSN,CMELTM,SPCFCT
      COMMON /ATMOSP/ DOX,SATM0,TRAATM,PERCLO,VVCFCT,KATM
      COMMON /DATE/ DAYN(13),STPPER
      COMMON /SUNS/ SCI,SCE,ANRAD,RADAN,PI,PI12,PI2,PASSO,SOLST,CALJOU
C
c      COMMON /meteomo/SP10M,RH2M,TPR,RITOT,ALB,T2M,BLWB
      COMMON /meteogr/day,T2Mgr,RH2Mgr,TPRgr,RADIgr,Pmed,SPeed
      COMMON /date/YYYY,MM,DD,HH
      DATA DAYN/0.,31.,28.,31.,30.,31.,30.,31.,31.,30.,31.,30.,31./
      DATA CPw,Wden,ROhu,Cpa,RatiKh,CKa,VLh,RatiKe,Rgas,Zw,Zt
     +/4217.8,1002.28,1.15,1005.,1.,0.41,2510000.,1.,2.87287,6.5,6./
c CPw = specific heat of water at constant pressure [J/(Kgm °C)]
c Wden = water density [kgm/m³]
c ROhu = density of humid air [Kgm/m³]
c Cpa = specific heat of air at constant pressure [J/(Kgm °C)]
c RatiKh = ratio of Kh/Km,turbulent eddy diffusivity for heat and momentum,respectively
c Cka = von Karman's constant
c VLh = latent heat of vaporisation in J/kg
c VLh =2.501-0.002361*Ta   ???
c RatiKe = ratio of Ke/Km,turbulent eddy diffusivity for mass? and momentum,respectively
c Rgas = constant of perfect gas (air)[m³ mb/(Kgm K)]
c Zw = height of temperature and humidity sensor (m)
c Zt = height of temperature and humidity sensor (m)  (2 m for Naltar)
      DATA CJOMME,HSNLA1,FHI,FHIeKs,BEERSN,GANESN/334944.,0.1,0.517,
     &0.002788,10.6,150./
      END
C
*****************************************
C      program DISMA
C	created  15/3/99 
C	updated  13/9/99 
*****************************************
C
      CALL SIMULA
	STOP
	END
C
CC
C     SIMULA  
C             
C     DEM :Digital Elevation Model: Matrix of elevations NR,NC             
CC

	SUBROUTINE SIMULA
c Calculates Potential Evapotranspiration by Penman-Monteith Equation.
c The numbers after formulas are a reference to 'Handbook of Hydrology', Maidment, 1993.
C Give the choice between Bolam (B) and Ground (G) meteorological data.
c      character*2 choice
c      write(*,*) 'Bolam (B) or Ground (G) meteorological data?'
c      read(*,'(a2)') choice
c      if (choice.eq.b) goto 100
c      if (choice.eq.g) goto 200

Cc Bolam
c      parameter(nhour=30,nc=17,nr=16)
c      REAL T2M(nr,nc),SP10M(nr,nc),RH2M(nr,nc),TPR(nr,nc),RITOT(nr,nc),
c     +ALB(nr,nc),BLWB(NR,NC)

Cc Ground
      PARAMETER (NOBS=5140,NSTRMX=900,NNR=353,NNC=262,NBAS=10)
C
      character*12 filemet(6),namefdem,namefdoc,iname,filebas(nbas)
CGG      CHARACTER*80 HEADER,FORMAL,HEADE1
c      character*12 filetp(NHOUR),filetp1(NHOUR),filetp2(NHOUR),
c     + fileta(nhour),filepr(nhour),fileu(NHOUR)
CRRMUL NAMES OF SCA FILES
      character*8 FILESCA(220)

      INTEGER T,TIN,TIME(4,8784),YEAR(8784),TMPTIM(4),STPDTH,STPTHE
	INTEGER NSENSP(2,10),NSENST(2,10),IMASK(10)
      INTEGER NESPO(nnr,nnc),NHEAD(6),NUMSE(6),idsub(10),NCOLMET(10,6)
      INTEGER (kind=2) IND(NNR,NNC),luse(nnr,nnc)
      INTEGER YYYY,MM,DD,HH,hsensp(10),hsenst(10),risp
      LOGICAL MASK(nnr,nnc),GROUND       
C
      real aspect(nnr,nnc),aspsmoo(nnr,nnc),day,T2Mgr,RH2Mgr,TPRgr,
     +RADIgr,Pmed,SPeed
CGG      real h0c,lai,laiout(nnr,nnc),DEMm(nnr,nnc),deltaH(nnr,nnc)
      real lai,laiout(nnr,nnc),hdebri(nnr,nnc)
CGG      real ETPrad(nnr,nnc),ETPaer(nnr,nnc),ETP(nnr,nnc),
CGG     +tamatr(nnr,nnc),prmatr(nnr,nnc),dtesp(nnr,nnc),viewft(nnr,nnc)
      real ETPrad(nnr,nnc),ETPaer(nnr,nnc),ETP(nnr,nnc),
     +dtesp(nnr,nnc),viewft(nnr,nnc)
      REAL SLTL(nnr,nnc),SLTR(nnr,nnc),SLL(8,nnr,nnc),SLR(8,nnr,nnc)
      REAL MELTXY,MEMUT,MESMUT,MESTT,MESSTT
      REAL MELT(NBAS),MESTXY(NBAS),LONGM(NBAS),LONGST(NBAS),LATM(NBAS)
	REAL LATST(NBAS),meltice (NBAS),pemed(nbas),PAREA(NBAS)
	REAL EVRMEDS(NBAS),EVRMEDI(NBAS),SNOWMEDS(nbas), NSNOWXY(NBAS)
	REAL RAINMEDS(NBAS),RAINMEDI(NBAS), MELTSNOW(NBAS),NSNOWICE(NBAS)
CGG      REAL DAYS(13),TMPXX(6),METEOXX(6,22),ALBTM(NNR,NNC),
CGG	+ALBSMOO(nnr,nnc)
      REAL DAYS(13),METEOXX(6,22),ALBTM(NNR,NNC),ALBSMOO(nnr,nnc)
      REAL DEM(nnr,nnc),XNORM(3,nnr,nnc),TMAXIM(nnr,nnc),
     +TEMCU(nnr,nnc),T2112(nnr,nnc)
C   Snowpack State variables
      REAL SWST(2,nnr,nnc),PDKEF(2,nnr,nnc),PACT(2,nnr,nnc),
     +WEQUI(nnr,nnc)
      REAL HSNPCK(nnr,nnc),DENSN(2,nnr,nnc),XYMELT(nnr,nnc),HEAFLO(2)
	real xyicem(nnr,nnc),pe(nnr,nnc),PXY(NNR,NNC),EVRXY(NNR,NNC),
     +SNOWXY(NNR,NNC)
Cgg      REAL HEABXY(nnr,nnc),VVCFCT(12),SAT(nnr,nnc),CMAX(nnr,nnc),
Cgg	+tk(nnr,nnc)
      REAL HEABXY(nnr,nnc),VVCFCT(12)
	REAL KATM,TSOIL(NNR,NNC)
	real Gv5d(NNR,NNC,10)
CRRMUL dates of the screenshot images
      integer idatess(214)
      COMMON /RECORD/ XX(10,6,8784),NTIMST
      COMMON /PARAME/ XS(20),X(20),INDX(20),N,NPAR,STEPS(20),TESTPS(20),
     +TESSTP(20),STEP(20)
      COMMON /SOIL/ ASOIL,SLAMDA,Tmsoil,Cv(NNR,NNC),DSOIL1(NNR,NNC)
      COMMON /SNOW/A0,A1,TSPRIN,ASPRIN,TRNCF,AM0,CAM,GANESN,ALFAS,
     +BETAS,STRN,CJOMME,HSNLA1,FHI,Swi,FHIeKs,BEERSN,CMELTM,SPCFCT
      COMMON /ATMOSP/ DOX,SATM0,TRAATM,PERCLO,VVCFCT,KATM
      COMMON /DATE/ DAYN(13),STPPER
      COMMON /PERIOD/ NMONIN,NDAYIN,NMONFI,NDAYFI,ISTAMP,NCYEAR,NFYEAR,
     +NMONAN,NDAYAN
      COMMON /SUNS/ SCI,SCE,ANRAD,RADAN,PI,PI12,PI2,PASSO,SOLST,CALJOU
      COMMON /HEAT/CPw,Wden,ROhu,Cpa,RatiKh,CKa,Z0,VLh,RatiKe,Rgas,ZW,ZT
      COMMON /STEF/ STEFAA
      COMMON /DEM/ DEM,XNORM,HSENSO,NXSENS,NYSENS,NR,NC,DEMDX,DEMDY,
     +AREAA,NGRIDS
      COMMON /SNPK/ SWST,PDKEF,PACT,HSNPCK,DENSN
      COMMON /MASKD/ MASK,NXPRNT,NYPRNT,IND,NUMSUB,IDSUB,IMASK
      COMMON/RADIA/SWAVXY(nnr,nnc),RLOWXY(nnr,nnc),HCONVXY(nnr,nnc),
     +HLATXY(nnr,nnc),HPXY(nnr,nnc),SNOGRXY(nnr,nnc)
      COMMON /AGGRES/ RAMUXY(NBAS),RASTXY(NBAS),HPMUXY(NBAS),
     +HPSTXY(NBAS),LONGM,LONGST,CONVM(NBAS),CONVST(NBAS),LATM,
     +LATST,MESTXY,SWEM(NBAS),SWEST(NBAS),SVM(7)
      COMMON /MOMENT/SCSW(3),SCRL(3),SCHC(3),SCHL(3)
Cgg    COMMON /PORTV/IDEF,JDEF,CSTR,FHIT,cinv,TK,SKS,CMAX,ww
      COMMON /VEG/FOR,ALBTM,ALBSMOO,LUSE,h0c(nnr,nnc)
      COMMON/SLOPEV/AGRID,VX,VY,VD,BDEF,VXR,VYR,DXY,VDR,FDEF,DDEF,
     +SLTL,SLTR,SLL,SLR
c      COMMON /meteomo/SP10M,RH2M,TPR,RITOT,ALB,T2M,BLWB
      COMMON /meteogr/day,T2Mgr,RH2Mgr,TPRgr,RADIgr,Pmed,SPeed
      COMMON /date/YYYY,MM,DD,HH
	COMMON/rhid/NCMAT,NRMAT,DELTAMAT,AMINXMAT,AMAXXMAT,AMINYMAT,
     +AMAXYMAT
      COMMON /THIESS/ TH(NSTRMX),XP(NOBS,NSTRMX),IPU(NNR,NNC),
     +COORDT(2,NSTRMX),COORDP(2,NSTRMX)
C
	data filebas/'medie_01.out','medie_02.out','medie_03.out',
     +'medie_04.out','medie_05.out','medie_06.out','medie_07.out',
     +'medie_08.out','medie_09.out','medie_10.out'/
      DATA DAyS,tcoresp/0.,31.,59.,90.,120.,151.,181.,212.,243.,
     +273.,304.,334.,365.,0.5/
c
      DATA gv5d/924860*0./
c       CS    =Specific heat of ice   (Eagleson) (2093.4)        [J/(Kg C)]
      DATA CS/2093.4/

C
CCC
CCC
C       INPUT METEOROLOGICAL DATA
CCC
      WRITE(*,*) '             INPUT DATA :'
C
c      OPEN(3,FILE='METEODAT.INP',STATUS='OLD',FORM='FORMATTED')
Cgg   File 39 contains filenames of the meteo variables and the association of
Cgg   each hydrologic unit to a meteo station
      OPEN(39,FILE='INFMETEO.INP',STATUS='OLD',FORM='FORMATTED')
 	do iskip=1,5
       read(39,*)
	enddo 
        READ(39,*) NSIMMX,NMISUR,DT,NDT,NMISP,NDTP,jumlam,nserb
C     nsimmx  NUMERO INTERVALLI SIMULAZIONE
C     nmisur  NUMERO INTERVALLI PIOGGIA misurata
C     DT      PASSO DI MISURA [ORE]
C     NDT     N. PASSI DI CALCOLO/PASSO DI MISURA []
C     NMISP   N. INTERVALLI DI DURATA DT DI MISURA PLUVIOMETRI []
C     NDTP    N. INTERVALLI MISURA PLUVIOMETRI/PASSO DI MISURA []
c     JUMLAM  n. intervalli di simulazione LAM da saltare in lettura
c     NSERB   n. SERBATOI
      DTT=DT*3600./NDT
C     DTT     PASSO DI CALCOLO [sec]
        READ(39,*) RISP
c RISP 0-RASTER,1-PLUVIOMETRI */
CGGG Inizializzo la variabile GROUND
        ISTART=0  
Cggg        IF(NMISP.GT.0) THEN
          GROUND=.TRUE.
Cggg        ELSE
Cggg          GROUND=.FALSE.
Cggg        ENDIF
cc
c     READING BOLAM OUTPUT FILES               (Z)
cc
      WRITE(*,*) 'NOME FILE PIOGGE LAM?_'
      READ(39,'(A12)') INAME
	WRITE(*,'(1x,A12)') INAME
      OPEN(52,FILE=INAME,IOSTAT=IOS,ERR=999,STATUS='OLD')
      WRITE(*,*) 'NOME FILE NEVE LAM?_'
      READ(39,'(A12)') INAME
	WRITE(*,'(1x,A12)') INAME
      OPEN(53,FILE=INAME,IOSTAT=IOS,ERR=999,STATUS='OLD')
      WRITE(*,*) 'NOME FILE T2M LAM?_'
      READ(39,'(A12)') INAME
	WRITE(*,'(1x,A12)') INAME
      OPEN(54,FILE=INAME,IOSTAT=IOS,ERR=999,STATUS='OLD')
CGGG  salta l'intestazione del file di neve
			do nnskip=1,18
				READ(53,*)
			enddo
CGGG  salta l'intestazione del file di temperatura
			do nnskip=1,18
				READ(54,*)
			enddo
c
	DO IREAD=1,6
       read(39,'(A12)') FILEMET(IREAD)
      write(*,'(A12)') FILEMET(IREAD)
	ENDDO
      read(39,*) (NHEAD(JH),JH=1,6)
      read(39,*) (NUMSE(JS),JS=1,6)
      read(39,*) NUMSUB
      read(39,*) 
	DO ISUB=1,NUMSUB
 	read(39,*) IDSUB(isub),(NCOLMET(ISUB,JS),JS=1,6),hsensp(isub),
     +hsenst(isub)
	 	WRITE(*,*) IDSUB(isub),(NCOLMET(ISUB,JS),JS=1,6),hsensp(isub),
     +hsenst(isub)
	ENDDO
      read(39,*) NUMSRIF
C      WRITE(*,*) NUMSRIF
      DO IMETEO=1,6
	  MUNIT=30+IMETEO
        OPEN(MUNIT,FILE=FILEMET(IMETEO))
	  IF(IMETEO.GT.2) THEN
	  NSKIP=NHEAD(IMETEO)
        DO ISKIP=1,NSKIP
          read(MUNIT,*)
	  enddo
	  ENDIF
	ENDDO 
C       INPUT PARAMETERS
CCC
c	write(*,*) 'simula'
      CALL INITIA(SMPLI,STPDTH)
      CALL SUNSET(XLAT,XLONG,STD,DEL,DELL,PASEC)      
C
       OPEN(1,FILE='DISMASIM.OUT',STATUS='unknown',FORM='FORMATTED')
Cgg	file disma2.out contains columns 'melt, discharge, port'
      OPEN(20,FILE='DISMA2.OUT',STATUS='unknown',FORM='FORMATTED')
      OPEN(22,FILE='MASSBAL.RST',STATUS='unknown',FORM='binary')

c      OPEN(25,FILE='ESPO.DAT',STATUS='unknown',FORM='FORMATTED')
Cgg	file .dat contains columns 'discharge, port, nsnow' (daily values)
      OPEN(26,FILE='DAILY.DAT',STATUS='unknown',FORM='FORMATTED')
Cgg	file perc1.out=BALANCE.rst and perc2.out=SWE.rst can be used to write a matrix at one tme step
Cgg	originally they were introduced to check percolation through thre first
Cgg    and the second layer perc3.OUT=LWC.rst.
      OPEN(27,FILE='SWE.DAT',STATUS='unknown',FORM='binary')
      OPEN(28,FILE='HSNOWPck.rst',STATUS='unknown',FORM='FORMATTED')
      open(71,file='DATESCA.OUT',STATUS='UNKNOWN',FORM='FORMATTED')
      OPEN(30,FILE='XYMELT.DAT',STATUS='unknown',FORM='binary')
C
      IF(ISTAMP.EQ.-1)THEN
CGG      writing snowdepth matrix at a fixed time (search for tx1)
CGG      OPEN(18,FILE='NEVE.rst',STATUS='unknown',FORM='FORMATTED')
       OPEN(64,FILE='SNOGRXY.rst',STATUS='unknown',FORM='FORMATTED')
c       OPEN(19,FILE='RAD.MAT',STATUS='unknown',FORM='FORMATTED')
CGG      OPEN(61,FILE='RLOWXY.rst',STATUS='unknown',FORM='FORMATTED')
CGG      OPEN(63,FILE='HLATXY.rst',STATUS='unknown',FORM='FORMATTED')

      ELSEIF(ISTAMP.EQ.-4)THEN
Cgg   writing statistics time series
C       OPEN(18,FILE='ORD2.DAT',STATUS='unknown',FORM='FORMATTED')
C       OPEN(19,FILE='ORDSUP.DAT',STATUS='unknown',FORM='FORMATTED')
       do iss=1,numsub
	ifs=21+iss
	 OPEN(ifs,FILE=filebas(iss),STATUS='unknown',FORM='FORMATTED')
	enddo
      ELSEIF(ISTAMP.EQ.-3)THEN
       OPEN(8,FILE='TESIDTM.DAT',STATUS='unknown',FORM='FORMATTED')
       OPEN(9,FILE='TESISPCK.DAT',STATUS='unknown',FORM='FORMATTED')
       OPEN(21,FILE='TESIBAL.DAT',STATUS='unknown',FORM='FORMATTED')      
c
      ELSEIF(ISTAMP.EQ.-2)THEN
C       OPEN(4,FILE='TESISUN.DAT',STATUS='unknown',FORM='FORMATTED')
       OPEN(7,FILE='TESIHEAT.DAT',STATUS='unknown',FORM='FORMATTED')
       OPEN(8,FILE='TESIDTM.DAT',STATUS='unknown',FORM='FORMATTED')
       OPEN(9,FILE='TESISPCK.DAT',STATUS='unknown',FORM='FORMATTED')
       OPEN(10,FILE='TESIRALB.DAT',STATUS='unknown',FORM='FORMATTED')      
CRR KARAKORAM 2.4 OPEN(60,FILE='SWAVXY.dat',STATUS='unknown',FORM='formatted')
CRR KARAKORAM 2.5 
	OPEN(60,FILE='SWAVXY.dat',STATUS='unknown',FORM='binary')
	OPEN(77,FILE='RLOWXY.dat',STATUS='unknown',FORM='binary')
	OPEN(80,FILE='HCONVXY.dat',STATUS='unknown',FORM='binary')
      OPEN(85,FILE='HLATEN.dat',STATUS='unknown',FORM='binary')
	OPEN(90,FILE='SURFTEMP.dat',STATUS='unknown',FORM='binary')
CGG      OPEN(61,FILE='RLOWXY.rst',STATUS='unknown',FORM='FORMATTED')      
CGG      OPEN(63,FILE='HLATXY.rst',STATUS='unknown',FORM='FORMATTED')

      ELSEIF(ISTAMP.EQ.2)THEN
       OPEN(10,FILE='TESIRALB.DAT',STATUS='unknown',FORM='FORMATTED')
      ENDIF
CCC
C       Headers of OUTPUT FILES are written
CCC
c	CALL WRITET(ISTAMP)

	IF(ISTAMP.EQ.-1) THEN
c	  WRITE(18,*)' Stato del manto nevoso in forma distribuita ' 
c        WRITE(19,*)' Componenti radiative in forma distribuita'
	ELSEIF(ISTAMP.EQ.-4)THEN
	do iss=1,numsub
	   ifs=21+iss
c    ICEm	    = Ice melt 
c    Smel     = Snow melt
c    P-Ei&t   = Rain-Pot Ev. over ice and ground
c    P-E+SNOW	= Rain+Snow-Pot. Ev.
	WRITE(ifs,'(2A35,3A52)')'Anno Me Gi Ora  Shortwave   +/-StD ',
     &'Longwave   +/-StD Sensible +/-StD  ',
     &'Latent     +/-StD RAINnosnow Evnosnow ICEmelt SnowCA',
     &' Rainonsnow Snowonsnow Evonsnow Meltonsnow Snowmelt ',
     &'NetRain+m Netpre SWE   SCA+ICE                      '
      WRITE(ifs,'(2A35,3A52)')'                  (W/m²)           ',
     &' (W/m²)            (W/m²)          ',
     &' (W/m²)             (mm)     ----------------(cells)',
     &' (mm) ----------------------------------------------',
     &'     -----------------------------------------  (mm)'
	enddo
c	WRITE(IFS,996) NYY,(TIME(J,T),J=1,4),RAMUXY(ISS),RASTXY(ISS),
c     +LONGM(ISS),LONGST(ISS),CONVM(ISS),CONVST(ISS),LATM(ISS),
c     +LATST(ISS),RAINMEDI(ISS),EVRMEDI(ISS),MELTICE(ISS),
c     +IFIX(NSNOWXY(ISS)),RAINMEDS(ISS),SNOWMEDS(ISS),EVRMEDS(ISS),
c     +MELTSNOW(ISS),MELT(ISS),PEMED(ISS),PAREA(ISS),SWEM(ISS)
c	  WRITE(18,*)'AnMe Gi Ora Shortwave   Longwave   Sensible heat ',
c     &'   Latent Heat  Nsnow'
c	  WRITE(18,*)'            CENTR - NORM [W2/m4]  '
c	  WRITE(19,*)'AnMe Gi Ora Shortwave   Longwave   Sensible heat ',
c     &'   Latent Heat   Nsnow'
c	  WRITE(19,*)'           TERZO - QUARTO  [W3/m6] - [W4/m8]  '
	ELSEIF(ISTAMP.EQ.2) THEN
	  WRITE(10,*)'      POINT Simulation of RADIATION & ALBEDO ' 
	  WRITE(10,*)'AnMeGi Ora  Altit  Azim PROJEC Radmis Rclsky ',
     &' AlbedoAlbedfDifradCloudi'
        WRITE(10,*)'               [deg]           [W/m²] [W/m²] ',
     &'[W/m²]   []    []    []    []  '
	ELSEIF(ISTAMP.EQ.-2) THEN
        WRITE(1,*)'   CELL(Station)  METEOROLOGICAL INPUT & MELT'
        WRITE(1,'(A37,A55)')'AnMeGi Ora   Pioggia Temper.  Radiaz.',
     &'Velven. Humid.   HSup   Wequi   Melt  Tsurf Portata'
        WRITE(1,'(A37,A55)')'             [mm]      [°C]    [W/m²]',
     &' [m/s]   [%]   [W/m²]    [m]    [mm]  [°C]   [m³/s]'
        WRITE(7,*)'   CELL(STATION)  HEAT EXCHANGE COMPONENTS   '
        WRITE(7,*)'AnMeGi Ora Swavbal Lwavbal Hconvec Hlatent  Hpdt',
     &'  Snogro Melt  Portata PotEv'
        WRITE(7,*)'           ----------------[W/m²]---------------',
     &'-------- [mm]  [mm]    [mm] '
        WRITE(8,*)'DISTRIBUTED HEAT EXCHANGE COMPONENTS [W/m²]-Melt[mm]'
        WRITE(8,*)'AnMeGi Ora X Y Ta RADNET Swavba Lwavba Hconve',
     &' Hlaten Hpdt Snogr(HG12) Heabal Melt'
        WRITE(8,*)'              [°C]   -------------[W/m²]-------',
     &'------------------------ [mm]'
        WRITE(9,*)'                       DISTRIBUTED SNOWPACK MODEL'
        WRITE(9,*)'AnMeGiOra   X  Y Hsnow Wequi Tspck Temp1 Temp2',
     &'Dens1Dens2 Sat1 Sat2 Watin Psnow Melt'
        WRITE(9,*)'             [m]   [mm]        [°C]       ',
     &' [kg/m³]  [1/1000]  [mm]  [mm] [mm]'
        WRITE(10,*)'      POINT Simulation of RADIATION & ALBEDO ' 
        WRITE(10,*)'AnMeGi Ora  Altit  Azim PROJEC Radmis Rclsky ',
     &'Rcldif AlbedoAlbedfDifradCloudi'
        WRITE(10,*)'               [deg]           [W/m²] [W/m²] ',
     &'[W/m²]   []    []    []    []  '      
	ELSEIF(ISTAMP.EQ.-3)THEN
        WRITE(8,*)' CELL(STATION)   EXCHANGE COMPONENTS [W/m²]-Melt[mm]'
        WRITE(8,*)'AnMeGi Ora  X  Y   Ta  RADNET Swavba Lwavba Hconve',
     &' Hlaten Hpdt Snogr(HG12) Heabal  Melt'
        WRITE(8,*)'              [ ° C] -------------[W/m²]-------',
     &'------------------------ [mm]'      
        WRITE(9,*)'                   CELL(STATION)   SNOWPACK MODEL'
        WRITE(9,*)'AnMeGiOra   X    Y Hsnow Wequi Tspck Temp1 Temp2',
     &'Dens1Dens2 Sat1 Sat2 Watin Psnow Melt'
        WRITE(9,*)'             [m]   [m]        [ ° C]        ',
     &' [Kg/m3]  [1/1000]  [mm]  [mm] [mm]'
        WRITE(21,*)'AnMeGiOra WequiM1 WEQUIV RAIN    SNOW    EVAP    ',
     &'MELT    DIFSTO'
        WRITE(21,*)'             [mm] '
	ENDIF
CC
C    Finished writing headers of OUTPUT files
CC     
C
C  Correction for bissextile years
C
CRRMUL
      IF ((NCYEAR.EQ.1996).OR.(NCYEAR.EQ.2000).OR.(NCYEAR.EQ.2004).OR.
     +(NCYEAR.EQ.2008).OR.(NCYEAR.EQ.2012).OR.(NCYEAR.EQ.2016).OR.
     +(NCYEAR.EQ.2020)) THEN
	   write(*,*) ' Anno bisestile'
	   DAYN(3)=29.
	   DO IM=3,13
	      DAYS(IM)=DAYS(IM)+1.
	   ENDDO
      ENDIF
CCC
C       SETUP OF INITIAL STATE OF SHORT WAVE AND HEAT BALANCE
C     DTM normal vector is computed once only
C     Initial state for albedo,temperature and snowmelt computations
CCC
C	 Reading aspect file
	namefdoc=' adaaspt.rdc'
	namefdem=' adaaspt.rst'
	call rhidri32 (namefdoc)
	call ridin3r (namefdem,aspect)
C	 Reading aspect file (smoothed 750 m or 1250 m)
	namefdoc=' aspsmoo.rdc'
	namefdem=' aspsmoo.rst'
	call rhidri32 (namefdoc)
	call ridin3r (namefdem,aspsmoo)

      DO 110 NX=1,NR
	DO 105 NY=1,NC
	 CALL VNORM(NX,NY)
C
C     Computing aspect     
C         
Cgg	 IF (ISTAMP.EQ.3)  THEN
Cgg	    ETANG=0.
Cgg	    NESPO(NX,NY)=0
Cgg	    E1=XNORM(2,NX,NY)
Cgg	    E2=XNORM(1,NX,NY)
Cgg	    CALL ESPO(E1,E2,ETANG,NESPO(NX,NY))


        alfaesp=aspsmoo(nx,ny)/180.*4.*atan(1.)
Cgg	  dtesp is negative for north aspect and positive for south aspect
        dtesp(nx,ny)=tcoresp*sin(alfaesp-2.*atan(1.))
c     	write(27,*) dtesp(nx,ny)
c	 ENDIF
C
	 SWAVXY(NX,NY)=0.
	 HEABXY(NX,NY)=0.
	 TEMCU(NX,NY)=0.
Cgg	 XYMELT(NX,NY)=-999.
	 XYMELT(NX,NY)=0.
	 TMAXIM(NX,NY)=-10. 
	 T2112(NX,NY)=10.+DAYS(NMONIN)+NDAYIN
105     CONTINUE
110   CONTINUE
      
CCC
C     DATA INPUT
CCC
c      READ(3,'(A80)') HEADE1
c      READ(3,'(A80)') HEADER
c      READ(3,'(A80)') FORMAL
ccc
ccc
ccc
C READ INPUT DATA
c	Reading lai matrix
	namefdoc='  laiinp.rdc'
	namefdem='  laiinp.rst'
	call rhidri32 (namefdoc)
	call ridri32 (namefdem,laiout)
c	Reading view factor matrix
	namefdoc='vtfactor.rdc'
	namefdem='vtfactor.rst'
	call rhidri32 (namefdoc)
	call ridri32 (namefdem,viewft)
ccc
C      I=1
      T=1
cgg      NTIMST=5856
      NTIMST=8784
C
      DAYINI=DAYS(NMONIN)+FLOAT(NDAYIN)
C      
      IDYEAR=NFYEAR-NCYEAR
      DPLUS=0.
      NYCUR=0
      DO IDY=1,IDYEAR       
	  NYCUR=NCYEAR+1
	  IF((NYCUR.EQ.1996).OR.(NYCUR.EQ.2000).OR.(NYCUR.EQ.2004).OR.
     +(NYCUR.eq.2008).OR.(nycur.eq.2012).or.(NYCUR.EQ.2016).OR.
     +(NYCUR.EQ.2020)) DPLUS=DPLUS+1.
      ENDDO
C      
      DAYFIN=FLOAT(IDYEAR)*365.+DPLUS+DAYS(NMONFI)+FLOAT(NDAYFI)
C      time loop
      DO 500 TIN=1,NTIMST
C
Cggg       IF(T.GT.NMISP) GROUND=.FALSE.
      IF (RISP.EQ.1) THEN
C     LE PIOGGE SONO LETTE DAI PLUVIOMETRI NEL FILE munit=31
CGGG	 opzione osservazioni a terra
         IF(GROUND) THEN
CGGG   passo anche la variabile ISTART per fare in modo che
CGGG   le coordinate dei pluviometri a terra siano lette una sola volta
     	 IF(ISTART.EQ.0) THEN
	 ISTART=-1
	 CALL RAINR
     +(NRMAT,NCMAT,32,NSTR,T,AMINXMAT,AMAXYMAT,DELTAMAT,ISTART,MASK)
	ISTART=0
	 CALL RAINR
     +(NRMAT,NCMAT,31,NSTR,T,AMINXMAT,AMAXYMAT,DELTAMAT,ISTART,MASK)
	 do isub=1,numsub
	   XCOORDP=COORDP(1,NCOLMET(ISUB,1))
	   YCOORDP=COORDP(2,NCOLMET(ISUB,1))
	   XCOORDT=COORDT(1,NCOLMET(ISUB,2))
	   YCOORDT=COORDT(2,NCOLMET(ISUB,2))
	   NSENSP(2,ISUB)=NINT((XCOORDP-AMINXMAT)/DELTAMAT+0.5)
	   NSENSP(1,ISUB)=NINT((AMAXYMAT-YCOORDP)/DELTAMAT+0.5)
	   NSENST(2,ISUB)=NINT((XCOORDT-AMINXMAT)/DELTAMAT+0.5)
	   NSENST(1,ISUB)=NINT((AMAXYMAT-YCOORDT)/DELTAMAT+0.5)
	 enddo
	 ENDIF
c
C	READ(3,FORMAL,END=550,ERR=999) NYY,(TMPTIM(J),J=1,4),
C     +(TMPXX(J),J=1,6)
C	READ(3,*) NJJ,NYY,(TMPTIM(J),J=1,3),(TMPXX(J),J=1,6)
Cgg Meteorological variables to be read are supposed to be 6
Cgg	 1. precipitation      [mm/h]
Cgg	 2.	temperature		     [°C]
Cgg	 3.	humidity              [%]
Cgg	 4.	pressure           [mbar]
Cgg	 5.	wind speed - 10 m  [km/h]
Cgg	 6.	solar radiation	   [W/m2]
C
	DO IMET=1,6
	  MUNIT=30+IMET
Cgg NUMSE is the number of columns to be read for each variable (one file per
Cgg variable)
	  NSERIE=NUMSE(IMET)

	  READ(MUNIT,*)NYY,(TMPTIM(J),J=1,3),trash,
     +(METEOXX(IMET,J),J=1,NSERIE)

c	if(imet.eq.1) write(*,*)NYY,(TMPTIM(J),J=1,3),trash,
c     +(METEOXX(IMET,J),J=1,NSERIE)
	IF (IMET.EQ.1) THEN
	   DO KS=1,NSERIE
           XP(T,KS)=METEOXX(IMET,KS)
	   ENDDO
	ENDIF
C

	ENDDO
Cgg tmptim(4) are the minutes (set to zero)
	TMPTIM(4)=0
CGGc	calcolo pesi di thissen
Cggg!!!!!
Cggg      CALL RAINR
Cggg	+(NRMAT,NCMAT,31,NSTR,T,AMINXMAT,AMAXYMAT,DELTAMAT,ISTART,MASK)
cc	 
cc ELSE
C         GROUND=.FALSE.  SI LEGGE DAL LAM
CGGG IL prossimo if e' necessario nel caso che si parta subito da BOLAM 
CRRR SI DEVONO LEGGERE UNA TANTUM LE COORDINATE E CALCOLARE  I PESI DI THISSEN
cc			IF((ISTART.EQ.0).or.(T.eq.nmisp+1)) THEN 
cc	            istart=0
cc				CALL RAINR
cc    +(IDIM,JDIM,52,NSTR,IMISUR,XMIN,YMAX,DX,ISTART)
cc		        CALL RAINR
cc    +(IDIM,JDIM,52,NSTR,IMISUR,XMIN,YMAX,DX,ISTART)
cc		if(jumlam.gt.0) then
cc				do ijum=1,jumlam  
cc					READ(52,*)
cc					READ(53,*)
cc					READ(54,*)
cc			end do  
cc		end IF 

			ENDIF     
Cgggg		fine if sulla variabile logica GROUND

cc			READ(52,*,err=999,END=2000)IANNO,IMESE,IGG,IHH,ITZ,IFC,
ccc    +(XP(IMISUR,KS),KS=1,NSTR)
CGGG  leggo la snowfall e la sottraggo subito dalla p totale
cc			READ(53,*,END=2000)IANNO,IMESE,IGG,IHH,ITZ,IFC,
ccc    +(sn(KS),KS=1,NSTR)
ccc			READ(54,*,END=2000)IANNO,IMESE,IGG,IHH,ITZ,IFC,
ccc     +(T2M(IMISUR,KS),KS=1,NSTR)
c	write(*,*) (xp(imisur,ks),sn(ks),t2m(imisur,ks),ks=1,2,3)
ccc			do ks=1,nstr
ccc				xp(imisur,ks)=xp(imisur,ks)-sn(ks)
ccc			enddo
	ENDIF
CGGGG	 fine di RISP=1
Cggg	ENDIF
CGGGG ???????
C
      IDYY=NYY-NCYEAR
      IF (IDYY.GT.0) THEN
	  DPLUS=0.
	  NYCUR=0
	  DO IYY=1,IDYY
	     NYCUR=NCYEAR+1 
	     IF ((NYCUR.EQ.96).OR.(NYCUR.EQ.00).OR.(NYCUR.EQ.04).
     +OR.(NYCUR.EQ.08).OR.(NYCUR.EQ.12).OR.(NYCUR.EQ.16).
     +OR.(NYCUR.EQ.20)) THEN
		DPLUS=DPLUS+1.
		IF ((NYCUR.EQ.NYY).AND.(DAYN(3).EQ.28)) THEN
		     DAYN(3)=29.
		     DO IM=3,13
				DAYS(IM)=DAYS(IM)+1. 
		     ENDDO
		ENDIF
	    
		   ELSEIF ((NYCUR.EQ.NYY).AND.(DAYN(3).EQ.29)) THEN
	     	 DAYN(3)=28
				DO IM=3,13
					DAYS(IM)=DAYS(IM)-1.
				ENDDO
	       ENDIF
	   ENDDO
      ENDIF
C
       DAYCUR=FLOAT(IDYY)*365.+DPLUS+DAYS(TMPTIM(1))+FLOAT(TMPTIM(2))
      IF(DAYCUR.GE.DAYINI) THEN
	  IF(DAYCUR.LE.DAYFIN) THEN
	  YEAR(T)=NYY
C  CHANGING METEO VARIABLE DEFINITION FOLLOWING INFMETEO.INP
C
C	  DO  490 IASS=1,4
C	   TIME(IASS,T)=TMPTIM(IASS)
C	   XX(I,IASS,T)=TMPXX(IASS)
C490       CONTINUE 
C	  XX(I,5,T)=TMPXX(5)
C	  XX(I,6,T)=TMPXX(6)
C
	  DO  490 IASS=1,4
	   TIME(IASS,T)=TMPTIM(IASS)
490       CONTINUE 
C       XX(I,1,T)=TMPXX(4)
C	  XX(I,2,T)=TMPXX(1)
C	  XX(I,3,T)=TMPXX(6)
C	  XX(I,4,T)=TMPXX(5)
C	  XX(I,5,T)=TMPXX(2)
C	  XX(I,6,T)=TMPXX(3)
	DO ISUB=1,NUMSUB
       XX(ISUB,1,T)=METEOXX(1,NCOLMET(ISUB,1))
       XX(ISUB,2,T)=METEOXX(2,NCOLMET(ISUB,2))
       XX(ISUB,3,T)=METEOXX(3,NCOLMET(ISUB,3))
       XX(ISUB,4,T)=METEOXX(4,NCOLMET(ISUB,4))
       XX(ISUB,5,T)=METEOXX(5,NCOLMET(ISUB,5))
       XX(ISUB,6,T)=METEOXX(6,NCOLMET(ISUB,6))
      ENDDO
	 T=T+1
	 ELSE

	   GOTO 550
	 ENDIF
	ENDIF
C1390  FORMAT(5X,3I2,1X,2I2,5F8.2,8X,F8.3)
500   CONTINUE
550   IF(T.LT.NTIMST) NTIMST=T-1
c	CLOSE(3)
      DO IMETEO=1,6
	  MUNIT=30+IMETEO
	  CLOSE(MUNIT)
	ENDDO
c
C  FINITA LA LETTURA DEI DATI METEO PER IL PASSO TEMPORALE CORRENTE
C
c hc = reference crop mean height [m]
C
C
      DIFRAD=0.
CGGG	meglio partire da CLOUDM=1, se la simulazione parte da mezzanotte
CGGG	!!!!!
CGGG      CloudM=0.
      CloudM=1.
      CLOUDN=1.
      SECDT=SMPLI*60.
      SECDTH=SECDT/FLOAT(STPDTH)
      J=1
C
C       Reset snowmelt variables     
C
C     Initial water storage at sensor cell [m]
       Wequin=Wequ(NXPRNT,NYPRNT,HSNPCK,DENSN)
       WRITE(*,*)'WWW=',WequIN,HSNPCK(NXPRNT,NYPRNT),
     + DENSN(1,NXPRNT,NYPRNT),DENSN(2,NXPRNT,NYPRNT)
       WEQUM1=Wequin
C     Cumulated melt at sensor
       CUMELT=0. 
C     Cumulated rain/snow -evaporation at sensor
       CUMWIN=0. 
C     Snow surface temperature
       Tsnow=0.
C     Snowpack mean temperature
       Tspack=0.
C
CGG      read(2,*)
CGG	READ (2,*) IDEF,JDEF,CSTR,FHIT,cinv,SATI,ww   
      read(2,*)itsobs
      read(2,*)iglacio
CRRMUL
CRRMUL  number of images (nimages) and idatess=mmdd  of the output 'screenshots' of 2D fields 
CRRMUL
      read(2,*) 
	read (2,*) nimages
	write (*,*) nimages
	do iim=1,nimages
	 read (2,'(I8)') idatess(iim)
CRRMUL      write (*,'(I9)') idatess(iim)
	write (71,'(I8)') idatess(iim)
	end do
CRRMUL im= counter of the screenshot images
	iim=1
	CLOSE (2)
	close (71)
      open(71,file='DATESCA.OUT',STATUS='UNKNOWN',FORM='FORMATTED')

	if(itsobs.gt.0) then
        OPEN(29,FILE='tsuolo.prn',STATUS='unknown',FORM='FORMATTED')
	  read(29,*)
	endif
	if(iglacio.gt.0) then
CgggRRMUL	reading debris depth formerly [cm] -> [m], now (m)

	  namefdoc='  debris.rdc'
	  namefdem='  debris.rst'
	  call rhidri32 (namefdoc)
	  call ridri32 (namefdem,hdebri)
	endif
Crrr
CCRRMUL	DO  IR=1,nrmat
CCRRMUL        DO IC=1,ncmat
CC	   IF((IR.EQ.85).AND.(IC.EQ.42)) WRITE(*,*)'HD1=',HDEBRI(IR,IC)
CCMUL now we read hdebris in m	   hdebri(ir,ic)=hdebri(ir,ic)/100.
CC	   IF((IR.EQ.85).AND.(IC.EQ.42)) THEN
CC	   WRITE(*,*)'HD2=',HDEBRI(IR,IC)
CC	PAUSE
CC	   ENDIF
CCRRMUL	   ENDDO
CCRRMUL      enddo
CGG  Initializing soil water content
CGG	 DO NX=1,NR
CGG	    DO NY=1,NC
CGG	       SAT(NX,NY)=SATI*CMAX(NX,NY)
CGG	    ENDDO
CGG	 ENDDO
CGG      SKSK=(FHIEKS*FHI*(1-SWI))**3
CGG       WRITE (*,*)'SKSK=', SKSK
CGG      SKS=SKSK*5.47*(10**6)*1000.
CGG       WRITE (*,*)'SKS=', SKS
CGG       WRITE(*,*)'SATI,CMAX,SAT=',SATI,CMAX(IDEF,JDEF),SAT(IDEF,JDEF)
CCC
C
c      CALL SLOPE(MASK)
CCC
C
C   ****  MAIN TEMPORAL LOOP ***
C								 
CC
	ISTEP=0
      DAYP=0.
      DAYD=0.
      NSNOW=0
	  CI=0.
	oldhconv=0.
      WRITE (26,*) '       discharge[mm/d] port[mm/d] nsnow'
C      time loop
      write(*,*)'ntimst=',ntimst
      DO 1000 T=1,NTIMST           
C
      NYY=YEAR(T)
      MONTH=TIME(1,T)
      IDAY=TIME(2,T)
      DAYCUR=DAYS(MONTH)+FLOAT(IDAY)
      IHOUR=TIME(3,T)
c  ????????????????? local time
      IMINU=TIME(4,T)      
      TIMENO=1104.+24.*(DAYS(MONTH)+IDAY)
      IF((IHOUR.EQ.1).AND.(IMINU.EQ.0)) THEN    
C
C     New day initialization
C
	  WRITE(*,*) NYY,MONTH,IDAY  
	  DO 610 NY=1,NC
	  DO 605 NX=1,NR
	   TEMCU(NX,NY)=TEMCU(NX,NY)+10.+TMAXIM(NX,NY)
	   TMAXIM(NX,NY)=-10. 
605       CONTINUE
610      CONTINUE
C
	   DATE=FLOAT(IDAY)
	   JDAY=INT(DAYS(MONTH))+IDAY
C                                                                          
	   CALL SOL(DATE,DAYS(MONTH),ET,ALFA,R)                          
C                  *** CALCOLO COSTANTE SOLARE DEL GIORNO ***
	   ETCMV=-ET
	   ALFAG=ALFA*RADAN                                                
	   DOX=SCE/(R*R)
	   ETIM=60.*ET                                                     
C   TCORR: Tempo vero - Tempo civile medio [ore]
	   TCORR=ET+DELL                                                   
C                                                                          
	   CALL DAILY(ALFA,XLAT,TR,TS,BETAM,RISANG,PSIR)                   
C                                                                          
	   TR=TR+TCORR                                                     
	   TS=TS+TCORR                                                     
	   TSR=TS-TR                                                       
C                  *** ESTREMI ROUTINE ORARIA 
C
	   IATR=NINT(TR-2.*PASSO)                                          
	   IATS=NINT(TS)                                                   
C     ATRS: Durata, in ore, della giornata (d)
	   ATRS=FLOAT(IATS-IATR)                                           
	   ILONG=INT(ATRS/PASSO)                                           
cccx           IF(ISTAMP.GE.1) THEN
C
C   End of new day inizialization
	 ENDIF  
C                      *** ROUTINE ORARIA ***
	 TCM=IHOUR
	 AST=TCM+TCORR
	 IF((AST.LT.TR).OR.(AST.GT.TS)) THEN
C
C        IT IS NIGHT
C
	   THETR=0.
	   THETD=0.
	  ELSE                                                                           
	      CALL ORA(XLAT,ALFA,AST,PSIZ,THETR)                           
C                                                                          
	      BETA=PI2-THETR                                               
	      COSB=COS(BETA)                                               
	      SINB=SIN(BETA)                   
	      PSIZD=PSIZ*RADAN                            
	      THETD=THETR*RADAN                                            
	      IF(THETD.GE.90.)THETD=90.                                    
	      IF(THETD.LE.0.)THEN
		 THETD=0.                                      
	      ENDIF
C          
cccx              IF(ISTAMP.GE.1) WRITE(4,1410) TCM,AST,PSIZD,THETD
ccx1410   FORMAT(3X,F6.2,6X,F6.2,4X,F7.3,4X,F7.3)
CCC                                                                          
C          ***  FINE ROUTINE ORARIA  ***                                   
C                                                                          
	 ENDIF         


C**********************************************
CC      
C     Set up Meteorological input
C     at sensor station
CC
C     Height of snowpack at sensor station 
      Hsnst=HSNPCK(NXSENS,NYSENS)
      Solrad=XX(NUMSRIF,6,T)
c
c**************************************************
C
C     Cloudiness computations now varying! 
C**************************************************
      IF(Solrad.LE.0) THEN 
c	SOLDIF=0.
	soldxy=0.
	solskydif=0.
	SOLBEA=0. 
	CLOUDI=CLOUDM
      ELSE
CGG
	vfsens=viewft(nxsens,nysens)
	asens=albsmoo(nxsens,nysens)/100.
CGG	CALL RADSUN(Solrad,THETR,SolDIF,SolBEA,DIFRAD,RCLSKY)
	CALL RADSUN(Solrad,THETR,SOLSKYDIF,SolBEA,DIFRAD,RCLSKY,vfsens,
     +asens,hsenso,SATM)
	IF((IHOUR.Gt.TR).AND.(IHOUR.Lt.TS)) THEN
	  IF((IHOUR.EQ.INT(TR+1)).AND.(IMINU.EQ.0)) Then
C
C    Reset for Cloudiness computations (Used in night hours)
C
	    CloudM=0.
	    CLOUDN=1.
	  ENDIF
	  CloudM=(CloudM*(CLOUDN-1)+Difrad)/CLOUDN
	  CLOUDN=CLOUDN+1.
	ENDIF 
	CLOUDI=DIFRAD
      ENDIF
Cgg  subroutine dtweight computes first gammat in clear sky conditions,
Cgg  following Ca'Zorzi, Dalla Fontana
Cgg  gammat in overcast conditions is later computed as a weighed average, where
Cgg  the weight is difrad (cloudiness)  
	call dtweight (ast,difrad,gammat)
      PSNOWP=0.         
      RAIN=0.
c
      RADNET=0.
      POTEV=0.
 	ETPW=0.
      Swavba=0.      
      Rlowab=0.      
      Hconve=0.      
      Hlaten=0.      
      Hpdt=0.        
C      Melt=0.
      ETPrs=0.
      ETPas=0.
      ETPsum=0.
Cgg	write(27,'(f8.3)') difrad
CCC
C    *** DEM Loop for SPATIAL VARIABLES COMPUTATIONS
CCC   
CCRRMUL
CCRRMUL	if(ihour.eq.1) then
CCRRMUL	  write(*,*)'iday,ast,tr,ts,atrs=',iday,ast,tr,ts,atrs
CCRRMUL	else
CCRRMUL	  write(*,*)'ast=',ast
CCRRMUL	endif

c      WRITE (*,*) 'DEM LOOP'
c Centro Mandron
c       DO 900 NX=249,249
c       DO 800 NY=295,295
c Dosson di Genova    
c       DO 900 NX=297,297
c       DO 800 NY=325,325
c Fronte Mandrone
c       DO 900 NX=164,164
c       DO 800 NY=357,357
c Passo Brizio
c       DO 900 NX=227,227
c       DO 800 NY=249,249
c Passo Adamè
c       DO 900 NX=314,314
c       DO 800 NY=287,287
c Ai caduti dell'Adamello
c       DO 900 NX=232,232
c       DO 800 NY=356,356
c Traccia Lobbia Alta
c       DO 900 NX=178,178
c       DO 800 NY=349,349
c Ghiacciaio Venezia
c       DO 900 NX=202,202
c       DO 800 NY=273,273
c SIM CELLA coniferous
C       DO 900 NX=160,160
C       DO 800 NY=435,435
c SIM CELLA grassland
c       DO 900 NX=200,200
c       DO 800 NY=500,500

       DO 900 NX=1,NR
       DO 800 NY=1,NC
CGG   Evapotranspiration of wet surfaces is set to zero initially
          ETP(NX,NY)=ETPW

      IF(MASK(NX,NY)) THEN
C        SIAMO DENTRO IL BACINO     
	  CALL DTM(NX,NY,THETR,PSIZ,SOLBEA,SOLBXY,HGRID,PROJEC,SATM)
Cgg	WRITE(*,*) HGRID,LAIOUT(NX,NY),LUSE(NX,NY)
CGG	correzione per SIBERIA
Cgg	albtm(nx,ny)=18
Cgg      HGRID=252.
Cgg	LAIOUT(NX,NY)=1.
Cgg	LUSE(NX,NY)=5
CGG
C
Cgg  subroutine assmeteo assigns to each subbasin the meteo time series 
Cgg  as defined in the input filre infmeteo.inp
Cgg  - ncheck=ind(nx,ny) is the subbasin ID read by the subroutine maskde in
Cgg    the mask.rst input file
Cgg  - numsub is the number of subbbasins included
Cgg  - idsub is the subbasin ID read in the infmeteo.inp file 
Cgg  - P,Tasns,Vw, rhumi,pmed are the precipitation, temperature, wind speed,
Cgg    humidity, pressure observed in the associated station
Cgg  - hsensp(numsub) are the altitudes of the precipitation sensors
Cgg  - hsenst(numsub) are the altitudes of the thermal sensors
Cgg  - hspmis is the altitude of the associated precipitation sensor
Cgg  - hstmis is the altitude of the associated thermal sensor
Cgg  - numsrif is the serial number of the reference subbasin (for missing data) 
C
        ncheck=ind(nx,ny)
	  call assmeteo(t,ncheck,NUMSUB,idsub,P,TASNS,VW,RHUMI,PMED,
     +hsensp,hsenst,hspmis,hstmis,numsrif,nsensp,nsenst,isp,jsp,ist,jst)
cggg   correggo precipitazione Passo Moro per avere totale evento tt4 registrato a C. Morelli
C?!@#[:-(        if((month.eq.6).and.(iday.eq.29)) P=P*4.
C	   Snow- measurement correction  
	if(hspmis.ne.hstmis) then
	  Tasp=TLAPSE(Tasns,hspmis,hstmis,gammat)
CrrESP
        if(isp.ge.1.and.isp.le.nr.and.jsp.ge.1.and.jsp.le.nc.and.
     +ist.ge.1.and.ist.le.nr.and.jst.ge.1.and.jst.le.nc) 
     +tasp=tasp+dtesp(isp,jsp)-dtesp(ist,jst)
	else
	  tasp=tasns
	endif
        IF((P.GT.0).AND.(Tasp.LT.(STRN+1.)))THEN
CGG	  not the total any more, but only the solid fraction of precipitation is
CGG     corrected; the correction factor of the rain fraction of precipitation is
CGG     assumed to be equal to unity
CGG
CGG 21/05/09 correzione gradiente altitudinale preecipitazione!! 
CGG          il gradiente è 0.312*Dh a scala annuale
CGG          P=P+0.000324*(hgrid-hspmis)
CGG          si preferisce convertire questo incremento in un fattore moltiplicativo
CGG 
CGG          P=P*TPCFCT
          P=P*PRLAPSE (HGRID, hspmis)
	    SNOW=RAINSN(Tasp,P,Rain)                         
	    P=SNOW*SPCFCT+RAIN
CGG      	P=P*SPCFCT
        ENDIF
c      Vw=XX(1,4,T)*VVCFCT(MONTH)
      Vw=Vw*VVCFCT(MONTH)
C	  write(*,*)'vw=',vw
C      DIS=XX(1,6,T)
	DIS=-999.
C
CCCCCC
C       WE ARE INSIDE THE BASIN
CCCCCC
C       RESET VARIABLES FOR THERMAL LOOP
	  STPTHE=STPDTH
	  SECTHE=SECDTH
c
Cgg a linear vertical distribution of air pressure is assumed (gammap in sub.for)
	  PRESS=PLAPSE(Pmed,HGRID,Hsenso)
c	  Ta=TLAPSE(Tasns,HGRID,Hsenso)
Cgg distribution of air tempaerature is assumed to result by the superimposition
Cgg of different effects:
Cgg   - linear vertical distribution with gammat varying following Ca'Zorzi, 
Cgg     Dalla Fontana and weighed by DIFRAD
Cgg   - bias induced by the aspect, described with an harmonic function of a
Cgg     smoothed aspect
	  Ta=TLAPSE(Tasns,HGRID,hstmis,gammat)
CGG	Ta=Tasns
CESP
	if(ist.ge.1.and.ist.le.nr.and.jst.ge.1.and.jst.le.nc) then
     	  ta=ta+dtesp(nx,ny)-dtesp(ist,jst)
	else
        ta=ta+dtesp(nx,ny)
	endif
c
	  IF(Ta.GT.-10.) TmAXIM(NX,NY)=AMAX1(TMAXIM(NX,NY),Ta) 
c	write(27,*)ta,P,gammat
	  IF(P.GT.0) THEN                                  
	    SNOW=RAINSN(Ta,P,Rain)                         
	    PSNOWP=SNOW/1000.
	    IF(Snow.GT.0) THEN
	      TEMCU(NX,NY)=0.                        
	      T2112(NX,NY)=10.+DAYS(MONTH)+IDAY
	    ENDIF
	  ELSE                                             
	     RAIN=0.                                       
	     SNOW=0.                                       
	  ENDIF                                            

c
c    relationships for TRNCF (fraction of shortwave radiation transmitted 
c    by the canopy) and FOR (fractional vegetation cover) are taken from
c    Kustas & Jackson, WRR, 1999
c
        lai=laiout(NX,NY)
	if (luse(nx,ny).eq.1) lai=0.
Cgg	if((luse(nx,ny).eq.1).and.(lai.gt.0.01))then
Cgg	  write(*,*) nx,ny,luse(nx,ny),lai
Cgg	endif
  	  TRNCF=exp(-0.67*lai)
	  FOR=1.-TRNCF
      if (rain.gt.0.) then
          rs=0.
       else 
          if(lai.gt.0) then
		  rs=200./lai
	    else
	      rs=0.
		endif   
	    rss=25.
	    rsn=0.
      endif
C  Computing distributed diffused radiation
	 SOLDXY=SOLSKYDIF*(1.-VIEWFT(NX,NY))+
     +ALBSMOO(NX,NY)/100.*SOLRAD*VIEWFT(NX,NY)
CGGG	  RADNET=SOLDIF+SOLBXY
       RADNET=SOLDXY+SOLBXY
	  STEFAA=STEFAN(Ta+273.14)
	  Ea=Rhumi/100.*EAS(Ta)
	IF(ITSOBS.GT.0) READ(29,*)NYY,(TMPTIM(J),J=1,3),TSOBS
c
	  IF((HSNPCK(NX,NY).LE.0).AND.(SNOW.LE.0)) THEN
C  No snow and no snowfall goto line 650
C            WRITE(*,*)'no snow',nx,ny,ihour,iday,month  
	    WATIN=0.
	    WATOU=0.
	    VAPINT=0.
Cgg  ICE!! !!!!!!!(att. per far simulare solo il Belvedere in maskde
Cgg  viene posto mask=.false. se id=19)
ccccc	    niter=0
		  IF(luse(nx,ny).eq.1)THEN
Cglacier
 	      HDEBRIS=HDEBRI(NX,NY)
Crr090613!!!?? Check Tice 
	      Tice=pact(1,nx,ny)
Crrice	   THERE IS ICE/debris
crrice    See Regine Hock, ZGS, 1998 page 96
	       if(hdebris.lt.0.02) then
crr090625	   TICE=amin1(0.,Ta)
		   Dwind=Dw(Vw,Ta,Tice,Hsnst,0.007)
Crr		THERE IS ICE only ALBEDI=albedo ice ESTIMATED AFTER ADAMELLO GLACIER 2003 08 23 ASTER IMAGE
Cag		Visto il risultato non del tutto positivo dal confronto con le misure di ablazione(Alberti,2005),  
Cag		provo ad aumentare (compatibilmente con l'immagine ADAMELLO GLACIER 2003 08 23 ASTER IMAGE) l'albedo Ice da 0.31 a 0.37 (Hock,1998)
C20060308			ALBEDI=0.37
			ALBEDI=albtm(nx,ny)/100.
C   Albedo in presence of Diffused Radiation Petzold-like formula
			VAR=0.499+0.0097*(10*DIFRAD)**3.
crrice Aged-ice Clear-sky albedo supposed constant
crrice Aged ice albedo in presence of diffused radiation and/or Cloudiness
			ALBEDO=AMIN1(ALBEDI*(1+VAR/100),1.)
			VAR=0.499+0.0097*(10.)**3.
			ALBEDF=AMIN1(ALBEDI*(1+VAR/100),1.)
	       else
Cdebris	there is debris albedo 
			ALBEDO=0.20
			ALBEDF=ALBEDO
			Dwind=Dw(Vw,Ta,Tice,Hsnst,0.01)
		   endif

	   	   Hpdt=HP(rain,Ta,Tice)/SECDT                           
Crr090626
		HeawhcDT=0.
	    HeaconDT=0.
	    HeameltDT=0.

	  DO 701 ISHEAT=1,STPTHE
          IF(RADNET.GT.0) THEN                             
C     Short-wave radiation balance at the ground                      
		    SWavba=GminS0(RADNET,ALBEDO,FOR,ALBEDO)
		   ELSE                                     
	        SWavba=0.                              
	    ENDIF                                            

C     Long-wave radiation balance at the ground  
	       if((hdebris.lt.0.001).and.(tice.gt.0.)) tice=0.
Cg		if (tice.gt.30.) tice=30.
     	   RLowab=RLio(FOR,Cloudi,Ta,Tice,Ea)
Cg	pause
Cg	 endif
	energy=swavba+rlowab
      if (hdebris.gt.0.02) then
	   hfict=0.
	   rsfict=rss
	  else
	   hfict=-88.
	   rsfict=rsn
	endif

cRR20211226
C2,7	
c	call etppm(Ta,Tice,RHumi,press,Vw,energy,h0c(nx,ny),ETP1,
c     +ETP2,ETPW,ra,deltag,gamma,hfict)
c
c2.13
	call etppm2(Ta,Tice,RHumi,press,Vw,energy,h0c(nx,ny),ETP1,
     +ETP2,ETPW,ra,deltag,gamma,hfict,hsnst,h0c(nxsens,nysens))
c
c
 	if(ra.gt.0.) then
        etpice=etpw*(deltag+gamma)/(deltag+gamma*(1+rsfict/ra))
	 else
	  etpice=etpw
	endif
    
          Eice=Rhumi/100.*EAS(Tice)
c          Evrx=Ev(Dwind,Eice,Ta,Tice)                        
          Evrx=Ev(Dwind,Ea,Ta,Tice)                        
CGG	PM-ETP
c20060206	evr=etpice
	evr=evrx
CGGice
 			VAPINT=VAPINT+Evr
CGG		if((iday.eq.10).and.(ihour.eq.9)) then
Crr090613 latent heat of sublimation instead of evaporation for ice
			Hlaten=Hl(Evr)*(CJOMME+Vlh)/VLh
CGG	IF (HLATEN.LT.-100.) HLATEN=-100.
CGG	IF (HLATEN.GT.100.) HLATEN=100.
CGG	 write(*,*)'ta,tice,hlaten=',ta,tice,hlaten
CGG	pause
CGG	 endif
	if (ra.gt.0.) then
        Hconve=Hc(1/ra,Ta,tice)
	else
	  hconve=0.
	endif
C
C     Ice Surface Heat balance
Cggrr090626	  Heabal=SWavba+RLowab+Hconve+Hlaten+Hpdt  
	  HG12=0.
	  Heawhc=SWavba+RLowab+Hconve+Hlaten+Hpdt  
	  HeawhcDT=HeawhcDT+Heawhc
	  if(hdebris.gt.0.001) then
cc	   HEACON=2.*SQRT(800.*2700.*2./PI/SECDTH)*tice
ccc   vecchia versione (21/06/01):  
cggrr 800 calore specific<o gneiss J/kg/K (Ministero, Colombo)
cggrr 1350 kg/m3
cggrr 0.35 W/K/m
cggrr nuova versione 24.2.2004  
cggrr 1821 bulk density gneiss+aria Kg/m3 da volumometro di Illy
cggrr 919 calore specifico gneiss Barontini da letteratura J/kg/K 
cggrr 0.58 conduttivita' termica detrito belvedere (da misure 2003) W/K/m  v. anche hpsuolo.xls
cggrr tice temperatura superficiale (doppia di quella media dello strato) 
         ticemed=tice/2.
	   Conk=2.*SQRT(919.*1821.*0.58/PI/SECDTH)
crr090626
         Cvzddt=919.*1821.*hdebris/secdth
         ticemed=(Heawhc+ticemed*(Cvzddt-Conk/2.))/(Cvzddt+Conk/2.)
	   ticemed=amin1(ticemed,0.)
	   HEACON=Conk*(ticemed+tice/2.)/2.
	   
cc    	tice=(Cvzddt*Tice+Heawhc-heacon)/Cvzddt
ccc   nuova versione (21/06/01):  
    	   tice=Ticemed*2.
Cgg    	tice=((Cvzddt-2./hdebris)*Tice+Heawhc)/Cvzddt
Cgg	   Tice=Heawhc/(2./hdebris)
Cgg	   if(hdebris.gt.0.05) then
Cgg	      if(niter.le.10) goto 123
Cgg	    else
Cgg	      if(niter.le.10) goto 123
Cgg	   endif
      else
CC         PURE ICE
C
crr090614 controllare metterei un ciclo sulla temperatura ghiaccio
cc    firn density =830 kg/m3 (Paterson, 1994; page 9)
CC    ThKICE (eq.2 page 205 Paterson, 1994) = 9,828 *exp(-0.0057 TK)
CC    A surface ice/firn layer 0.1 m deep is assumed to exchange conductive heat with a 0°C isothermal glacier
c       CS    =Specific heat of ice   (Eagleson) (2093.4)        [J/(Kg C)]

	   ticemed=tice/2.
	   TK=TICE+273.14
	   ThKICE=9.828 *exp(-0.0057*TK)
	   Conk=2.*SQRT(CS*830.*ThKICE/PI/SECDTH)
crr090626	   HEACON=Conk*ticemed
     	   Cvzddt=CS*830.*0.1/secdth
	   ticemed=(Heawhc+ticemed*(Cvzddt-Conk/2.))/(Cvzddt+Conk/2.)
Crr090626 heat content available for ice fusion
	   Heamelt=amax1(Cvzddt*Ticemed,0.)
	   ticemed=amin1(ticemed,0.)
	   HEACON=Conk*(ticemed+tice/2.)/2.
    	   tice=Ticemed*2.
crr090626         tice=amin1(Tice+Heabal/Cvzddt*2.,0.)
	endif
	HEACONDT=HEACONDT+HEACON
	HeameltDT=HeameltDT+Heamelt
701    continue
      Heamelt=HeameltDT/STPTHE
      HEACON=HEACONDT/STPTHE
      Heawhc=HeawhcDT/STPTHE

      PACT(1,NX,NY)=TICE
CGG icemelt is computed in mm
CGGRR
      if(hdebris.gt.0.001) then
c20060201 no snow, debris over ice
      	xyicem(nx,ny)=amax1(HeaCON*SECDT/cjomme,0.0)
	else
c20060201
crr090626	no snow , just ice melt only if ice is at 0°C    
          Heabal=Heawhc-Heacon
          if(ticemed.ge.0.0) then 
			xyicem(nx,ny)=amax1(Heamelt*SECDT/cjomme,0.0)
	    else
			xyicem(nx,ny)=0.
          endif 
	endif
CGG watin la produzione totale di defusso (pioggia + icemelt), depurata delle perdite di evaporazione viene 
CGG viene invece calcolata in m 
c20060201
	watin=rain/1000.-evr/Wden*SECDT
	meltxy=watin*1000.+xyicem(nx,ny)
C Karakoram 2.10
	Gv5d(nx,ny,1)=wequi(nx,ny)*1000.
	Gv5d(Nx,Ny,2)=hsnpck(nx,ny)
	Gv5d(nx,Ny,3)=meltxy
	Gv5d(nx,ny,4)=lwc(nx,ny,hsnpck,swst)
	Gv5d(nx,ny,6)=SWAVBA
	Gv5d(nx,ny,7)=RLOWAB
	Gv5d(nx,ny,8)=HCONVE
	Gv5d(nx,ny,9)=HLATEN
	Gv5d(NX,NY,10)=TICE

Cggg ICE!!!!!!!!!!!!
	  ELSE
Cggg  (ICE version!)there is ground goto 640

		albedo=albtm(NX,NY)/100.
		IF(T.eq.1) Tsoil(NX,NY)=Ta
		IF ((ITSOBS.GT.0).AND.(TSOBS.GT.-999.)) tsoil(nx,ny)=tsobs

          Dwind=Dw(Vw,Ta,Tsoil(nx,ny),Hsnst,0.01)
c		write(*,*)'dw,vw,ta,tsoil,hsnst=',dwind,vw,ta,tsoil,hsnst                                       
CGG	Rainfall heat contribution is now computed after the interception,
CGG    taking into account only what is left, after evaporation from interception
CGG	    Hpdt=HP(Rain,Ta,Tsoil)/SECDT                           
C	    POTEV=POTEV+Potevr(DAYCUR)

		IF(RADNET.GT.0) THEN                             
C
C     Short-wave radiation balance at the ground                      
			SWavba=GminS0(RADNET,ALBEDO,FOR,ALBEDO)
		ELSE                                     
			SWavba=0.                              
		ENDIF                                            
C     Long-wave radiation balance at the ground  
 	  RLowab=RLio(FOR,Cloudi,Ta,Tsoil(nx,ny),Ea)      
C CALLS THE SUBROUTINE THAT CALCULATES ETP FOLLOWING PENMAN-MONTEITH EQUATION.
	energy=swavba+rlowab
C2.13
	call etppm2(Ta,Tsoil(nx,ny),RHumi,press,Vw,energy,h0c(nx,ny),ETP1,
     +ETP2,ETPW,ra,deltag,gamma,0.,hsnst,h0c(nxsens,nysens))

CRR20211226      
c2.7	call etppm(Ta,Tsoil(nx,ny),RHumi,press,Vw,energy,h0c(nx,ny),ETP1,
c2.7     +ETP2,ETPW,ra,deltag,gamma,0.)
c	write(*,*)'etpw,etp1,energy=',etpw,etp1,energy,ta,tsoil(nx,ny)
c
cccc	ETP in kg/m2/s (=mm/m2/s)
c	in INTERCET values are expressed in mm/h
c
c	write(*,*)'rs,ra,deltag,gamma=',rs,ra,deltag,gamma
		if(lai.gt.0.) then
	
			CALL  INTERCET(T,rain*FOR,ETPW*3600.*FOR,
     +CI,QI,EI,ET,ETA,TR,rs,ra,deltag,gamma)
c	write(*,*)'etpw,ei,et,eta=',etpw*3600.,ei,et,eta
		else
			ei=0.
			et=0.
			eta=0.
			tr=rain*FOR
		endif
CGGG heat contribution due to the rainfall (through the trees)
	    HPINP=TR+RAIN*(1-FOR)
Cggg	    Hpdt=HP(HPINP,Ta,Tsurf)/SECDT                           
	    Hpdt=HP(HPINP,Ta,Tsoil(nx,ny))/SECDT                           
c 
		if(ra.gt.0.) then
      etpsoil=(etpw-eta/3600.)*(deltag+gamma)/(deltag+gamma*(1+rss/ra))
		else
	etpsoil=etpw-eta/3600.
		endif
Cgg       Esoil=Rhumi/100.*EAS(Tsurf)
       Esoil=Rhumi/100.*EAS(Tsoil(nx,ny))
c                              
Cgg       EvrX=Ev(Dwind,Esoil,Ta,Tsurf)                        
       EvrX=Ev(Dwind,Esoil,Ta,Tsoil(nx,ny))                        
c	  write(27,'(F12.10,F8.6,3(F9.5))')evrx,dwind,ea,ta,tsoil(NX,NY)                         
c	WRITE(*,*)'SUOLO, ETP1,ETP2,ETPW,EVRX=',ETP1,ETP2,ETPW,EVRX
c	  evr=amin1(et/3600.,evrx)
	evr=et/3600.+Etpsoil
	  VAPINT=VAPINT+Evr
	  Hlaten=Hl(Evr)
C
C     Soil Surface Heat balance
Cgg    All except the convective heat which is now computed after the 
Cgg    solution of T soil           
Cgg	  Heabal=SWavba+RLowab+Hconve+Hlaten+Hpdt  
	  Julday=DAYS(MONTH)+IDAY
	  HG12=Hg(dsoil1(NX,NY),julday,rain)
c
	  HeaWHC=SWavba+RLowab+Hlaten+Hpdt+HG12
c	EXTINC=EXP(-BEERSN*dsoil1(nx,ny))        
	  EXTINC=0.       
C  2 Soil Layers
	  heawhc=Heawhc-Swavba*EXTINC
c
      tsold=tsoil(nx,ny)
	  Tsoil(NX,NY)=tsc(heaWHC,Ta,Cv(NX,NY),dsoil1(NX,NY),
     +Tsold,secdt,1./ra)

	  if(tsoil(nx,ny).lt.-30) tsoil(nx,ny)=-30.
		if (ra.gt.0.) then
			 Hconve=Hc(1/ra,Ta,TSoil(nx,ny))
		else
			hconve=0.
		endif
Cgg
    	xyicem(nx,ny)=0.

C Water Input [m/h]
      WATIN=TR/1000.+TRNCF*RAIN/1000.-Evr/Wden*SECDT
Cg	 else
Cg	 hconve=0.
Cg	endif 
Cggg
	    MELTXY=watin*1000.
C Karakoram 2.10
	Gv5d(nx,ny,1)=wequi(nx,ny)*1000.
	Gv5d(Nx,Ny,2)=hsnpck(nx,ny)
	Gv5d(nx,Ny,3)=meltxy
	Gv5d(nx,ny,4)=lwc(nx,ny,hsnpck,swst)
	Gv5d(nx,ny,6)=SWAVBA
	Gv5d(nx,ny,7)=RLOWAB
	Gv5d(nx,ny,8)=HCONVE
	Gv5d(nx,ny,9)=HLATEN
	Gv5d(NX,NY,10)=TSOIL(NX,NY)

640	ENDIF
CGGG - fine if sul landuse=1 del ghiacciaio   da qui abbiamo o terreno o ice             
CGGG	calcolo l pioggia efficace (oppure pioggia - evt +icemelt) in mm
	pe(nx,ny)=watin*1000.
CCC
c      call soilbal(LAI,W,Wmax,smax,smin1,smin2,ratioGR,tpr,TReff,
c     + ETeff,ESUeff,ESGeff,EIeff,ETdef,ESUdef,ESGdef,ESdef,EVAPO,Q)
CCC
	    WEQUIV=0.
	    WEQUI(NX,NY)=0.
    	  xymelt(nx,ny)=meltxy
	  SWAVXY(NX,NY)=SWAVBA
	  RLOWXY(NX,NY)=RLOWAB
	  HCONVXY(NX,NY)=HCONVE
	  HLATXY(NX,NY)=HLATEN
	  HEABXY(NX,NY)=HEAwhc
       HPXY(NX,NY)=HPDT
       SNOGRXY(NX,NY)=HG12
	Gv5d(nx,ny,1)=wequi(nx,ny)*1000.
	Gv5d(Nx,Ny,2)=hsnpck(nx,ny)
	Gv5d(nx,Ny,3)=xymelt(nx,ny)
	Gv5d(nx,ny,4)=lwc(nx,ny,hsnpck,swst)
C KARAKORAM 2.10
	  IF(ISTAMP.EQ.-2) THEN
	ivar=0
			IF(LUSE(NX,NY).EQ.1) THEN
cggrr ghiacciaio
		    if(hdebris.lt.0.001) IVAR=2
       WRITE(8,790) NYY,(TIME(J,T),J=1,3),NX,NY,Ta,RADNET,SWavba,RLowab,
     +Hconve,Hlaten,Hpdt,hg12,Heawhc,TICE,ivar,XYICEM(NX,NY)
Cggrr Heacon invece di Hg12, XYICEM(X,Y) invece di PE ????
c2.6     +Hconve,Hlaten,Hpdt,Heacon,Heawhc,TICE,ivar,XYICEM(NX,NY)
			ELSE
	WRITE(8,790) NYY,(TIME(J,T),J=1,3),NX,NY,Ta,RADNET,SWavba,RLowab,
     +Hconve,Hlaten,Hpdt,hg12,Heawhc,tsoil(nx,ny),ivar,PE(NX,NY)
			ENDIF

		  WRITE(9,888) NYY,(TIME(J,T),J=1,3),NX,NY,HSNPCK(NX,NY)
888         FORMAT(I4,3I3,2I4,F5.2)
	  ENDIF

C  Goto 798
C		   ELSE THERE IS SNOW ON THE GROUND
C
650	  ELSE
	  tsoil(nx,ny)=0.
	  ALBEDO=ALBDIF(DIFRAD,T2112(NX,NY),TEMCU(NX,NY))
C
C   Albedo in presence of Diffused Radiation (as supposed in forest
	  ALBEDF=ALBDIF(1.,T2112(NX,NY),TEMCU(NX,NY))                            
C
C  Diffused radiation measured at station and BEAM radiation computed by DTM
	  Dwind=Dw(Vw,Ta,PACT(1,NX,NY),Hsnst,Z0)                                       
	  Hpdt=HP(rain,Ta,PACT(1,NX,NY))/SECDT                           
C	  
	  IF(RADNET.GT.0) THEN                             
C
C     Short-wave radiation balance                       
	    SWavba=GminS0(RADNET,ALBEDO,FOR,ALBEDF)
	  ELSE                                     
	    SWavba=0.                              
	  ENDIF                                            
	  WATOU=0.
C       INTEGRATED EVAPORATION
	  VAPINT=0.
c
	  IF(HSNPCK(NX,NY).LE.0.05) THEN
		STPTHE=STPDTH*4
		SECTHE=SECDTH/4.
	  ENDIF
	  DO 700 ISHEAT=1,STPTHE
C
C     Long-wave radiation balance   
	  Tspack=Tspk(NX,NY,HSNPCK,DENSN,PACT)
	  if(tspack.lt.-30) tspack=-30.
	  if(tspack.gt.0) tspack=0.
	  RLowab=RLio(FOR,Cloudi,Ta,Tspack,Ea)      
C                                                  
C                                                  
C     Latent heat                              
C
C CALLS THE SUBROUTINE THAT CALCULATES ETP FOLLOWING PENMAN-MONTEITH EQUATION.
	  energy=swavba+rlowab
CRR20211226
c2,7
c      call etppm(Ta,pact(1,nx,ny),RHumi,press,Vw,energy,h0c(nx,ny),
c     +ETP1,ETP2,ETPW,ra,deltag,gamma,hsnpck(nx,ny))
c2.13
      call etppm2(Ta,pact(1,nx,ny),RHumi,press,Vw,energy,h0c(nx,ny),ETP1
     +,ETP2,ETPW,ra,deltag,gamma,hsnpck(nx,ny),hsnst,h0c(nxsens,nysens))

c
cccc	ETP in kg/m2/s (=mm/m2/s)
c	in INTERCET values are expressed in mm/h
c
	  if(hsnpck(nx,ny).ge.0.9*h0c(nx,ny))then
	      lai=0.
	      for=0.
	      trncf=1.
	  endif
	  if(lai.gt.0.) then
	     CALL  INTERCET(T,RAIN*FOR,ETPW*3600.*FOR,
     +CI,QI,EI,ET,ETA,TR,rs,ra,deltag,gamma)
	  else
	     ei=0.
	     et=0.
	     eta=0.
           tr=rain*FOR
	  endif
c
	 if(ra.gt.0.)then
	 etpsnow=(etpw-eta/3600.)*(deltag+gamma)/(deltag+gamma*(1+rsn/ra))
	  else
	 etpsnow=etpw-eta/3600. 
	endif
	  Esnow=Rhumi/100.*EAS(pact(1,nx,ny))
c        EvrX=Ev(Dwind,Ea,Ta,PACT(1,NX,NY))                        
        EvrX=Ev(Dwind,Esnow,Ta,PACT(1,NX,NY))                        
c	  write(27,'(F12.10,F8.6,3(F9.5))')evrx,dwind,ea,ta,pact(1,NX,NY)                         
c      evr=eta/3600.*FOR+EVRX*trncf
c	write(*,*)'snow:ev,etp,etpw,eta=',evrx,etpsnow,etpw,eta/3600.
c20060206        evr=et/3600.+Etpsnow
        evr=evrx 
c	  evr=amin1(et/3600.,evrx)
	  VAPINT=VAPINT+Evr
	  Hlaten=Hl(Evr)                           
C     Convective heat             
Cgg	  Hconve=Hc(dwind,Ta,PACT(1,NX,NY))
	  if (ra.gt.0.) then
	    Hconve=Hc(1/ra,Ta,PACT(1,NX,NY))
 	  ELSE
	    hconve=0.
	  ENDIF
	    
Cg	if(hconve.lt.-700) then
Cg	  write(*,*)'ra,ta,pact1=',ra,ta,pact1
Cg	stop
Cg	endif
C
C     Snow Surface Heat balance                
	  Heabal=SWavba+RLowab+Hconve+Hlaten+Hpdt  
CCC
C       SNOWPACK ROUTING
CCC
	  HSNST1=AMIN1(HSNPCK(NX,NY),HSNLA1)
	  EXTINC=EXP(-BEERSN*HSNST1)        
c	  if(isheat.eq.1)
c	 +write(27,'(7f11.3)')Heabal,Hpdt,pact(1,nx,ny),rain,tr
	  IF(HSNPCK(NX,NY).GT.HSNLA1) THEN 
C  2 Layers
	  HEAFLO(1)=Heabal-Swavba*EXTINC
	  HEAFLO(2)=Swavba*(EXTINC-EXP(-BEERSN*(HSNPCK(NX,NY))))
	  ELSE
C   1 Layer only
	  HEAFLO(1)=Heabal-Swavba*EXP(-BEERSN*HSNPCK(NX,NY))
	  HEAFLO(2)=0.
	  ENDIF
C Water Input [m/measurement step e.g. h]
	  WATIN=TR/1000.+TRNCF*RAIN/1000.-Evr/Wden*SECDT
C
	  CALL SNPACK(NX,NY,SwST,PDKEF,HSNPCK,PACT,DENSN,WATIN,WATOUT,
     &HEAFLO,PSNOWP,Ta,SECDT,SECTHE,STPTHE,Snogro)
	  WATOU=WATOU+WATOUT
C  End of Heat equation temporal loop
700    CONTINUE


c
CCC
C         Melt equivalent [mm H20/(m2 interval)]
CCC
	xyicem(nx,ny)=0.
      MeltXY=WATOU*1000.                                               
      XYMELT(NX,NY)=meltxy
C20060201  over snow we are not interested in net rainfall 
	pe(nx,ny)=0.
	Gv5d(nx,ny,1)=wequi(nx,ny)*1000.
	Gv5d(Nx,Ny,2)=hsnpck(nx,ny)
	Gv5d(nx,Ny,3)=meltxy
	Gv5d(nx,ny,4)=lwc(nx,ny,hsnpck,swst)
	  Wequiv=Wequ(NX,NY,HSNPCK,DENSN)
	  IF((ihour.EQ.11).AND.(WEQUIV.GT.0)) NSNOW=NSNOW+1
	  WEQUI(NX,NY)=Wequiv
	  SWAVXY(NX,NY)=SWAVBA
	  RLOWXY(NX,NY)=RLOWAB
	  HCONVXY(NX,NY)=HCONVE
	  HLATXY(NX,NY)=HLATEN
        HPXY(NX,NY)=HPDT
        SNOGRXY(NX,NY)=SNOGRO
	  HEABXY(NX,NY)=HEABAL
	  IF(ISTAMP.EQ.-2) THEN
	  IEQUIP=IFIX(Wequiv*1000.)
	  IWST1P=IFIX(SWST(1,NX,NY)*1000.)
	  IWST2P=IFIX(SWST(2,NX,NY)*1000.)
	  WATINW=WATIN*1000.
	ivar=1
CRR KARAKORAM 2.3
C Karakoram 2.10
	Gv5d(nx,ny,6)=SWAVBA
	Gv5d(nx,ny,7)=RLOWAB
	Gv5d(nx,ny,8)=HCONVE
	Gv5d(nx,ny,9)=HLATEN
	Gv5d(NX,NY,10)=PACT(1,NX,NY)


c2.6      WRITE(8,790) NYY,(TIME(J,T),J=1,3),NX,NY,Ta,RADNET,SWavba,RLowab,
c2.6     +Hconve,Hlaten,Hpdt,Snogro,Heabal,PACT(1,NX,NY),ivar,MELTXY
C     +Hconve,Hlaten,Hpdt,Snogro,Heabal,MeltXY,ivar
      
c2.6	  WRITE(9,795) NYY,(TIME(J,T),J=1,3),NX,NY,HSNPCK(NX,NY),
c2.6     &IEQUIP,Tspack,(PACT(IL,NX,NY),IL=1,2),
c2.6     &(DENSN(IL,NX,NY),IL=1,2),IWST1P,IWST2P,WATINW,SNOW,MELTXY
795         FORMAT(I4,3I3,2I4,F5.2,I5,5F6.1,2I4,F7.2,F6.1,F6.2)
    	ENDIF
      ENDIF
C                  ENDIF for HSNPCK(NX,NY)=0. AND SNOW=0.
C     Melt integrated over the basin area
790         FORMAT(I4,3I3,2I4,F6.1,F8.1,F7.1,3F8.1,2F7.1,F9.1,F9.3,i4,
     +f8.2)
  
       ETPRAD(nx,ny)=ETP1
       ETPAER(nx,ny)=ETP2
       ETP(nx,ny)=ETPw
	 PXY(NX,NY)=WATIN*1000.+SNOW
	SNOWXY(NX,NY)=SNOW
	EVRXY(NX,NY)=EVR*SECDT
	  IF((HSNPCK(NX,NY).LE.0).AND.(SNOW.LE.0)) THEN
c karakoram 2.14
c	Gv5d(Nx,Ny,5)=Gv5d(nx,ny,5)+PXY(nx,ny)-xymelt(nx,ny)-f,NY)
      if(luse(nx,ny).eq.1) then
	 Gv5d(Nx,Ny,5)=Gv5d(nx,ny,5)-xyicem(nx,ny)-amax1(EVRXY(NX,NY),0.)
	endif
	  else
c	Gv5d(Nx,Ny,5)=Gv5d(nx,ny,5)+PXY(nx,ny)-xymelt(nx,ny)
        endif
          ETPrs=ETPrs+ETP1
          ETPas=ETPas+ETP2
          ETPsum=ETPsum+ETPw
CC798     Melt=Melt+MELTXY
798     MeltXY=MELTXY
C
C                 Output for CELL at OUTPUT station
C
       IF((NX.EQ.NXPRNT).AND.(NY.EQ.NYPRNT))THEN 
	  ALBEDS=ALBEDO
	  ALBDFS=ALBEDF
	  HEABAS=HEABAL
	  SMELTS=MELTXY
CCC
C       Water balance at sensor [m]
CCC
c	  CUMWIN=CUMWIN+p/1000.-VAPINT/Wden*SECTHE
	  CUMWIN=CUMWIN+TR/1000.+P/1000*TRNCF-VAPINT/Wden*SECTHE
	  CUMELT=CUMELT+MELTXY/1000.
	  IF((ISTAMP.EQ.-2).OR.(ISTAMP.EQ.2)) THEN
	  WRITE(10,796) NYY,(TIME(J,T),J=1,3),THETD,PSIZD,PROJEC,SOLRAD,
     &RCLSKY,SOLBEA,ALBEDO,ALBEDF,DIFRAD,CLOUDI
796        FORMAT(I2,3I3,F5.1,1X,F5.1,1X,f5.2,3(1X,F6.1),4F6.3)        
	  SWavbs=Swavba
	  RLowas=RLowab
	  Hconvs=Hconve
	  Hlates=Hlaten
	  Snogrs=Snogro
	  Hpdts=Hpdt
	  ELSEIF(ISTAMP.EQ.-3) THEN
	  WEQUIP=Wequiv*1000.
	  WEQUMP=Wequm1*1000.
	  IEQUIP=IFIX(Wequip)
	  IWST1P=IFIX(SWST(1,NX,NY)*1000.)
	  IWST2P=IFIX(SWST(2,NX,NY)*1000.)
	  WATINW=WATIN*1000.
	  EVAPOR=VAPINT/WDEN*SECTHE*1000.
	  DIFSTO=WEQUMP-WEQUIP+P-EVAPOR
      WRITE(8,790) NYY,(TIME(J,T),J=1,3),NX,NY,Ta,RADNET,SWavba,RLowab,
     +Hconve,Hlaten,Hpdt,Snogro,Heabal,MeltXY
      WRITE(9,795) NYY,(TIME(J,T),J=1,3),NX,NY,HSNPCK(NX,NY),
     &IEQUIP,Tspack,(PACT(IL,NX,NY),IL=1,2),
     &(DENSN(IL,NX,NY),IL=1,2),IWST1P,IWST2P,WATINW,SNOW,MELTXY
      WRITE(21,799) NYY,(TIME(J,T),J=1,3),WEQUMP,WEQUIP,RAIN,SNOW,
     &EVAPOR,MELTXY,DIFSTO
799         FORMAT(I2,3I3,7F8.3)
	  WEQUM1=WEQUIV
	  ENDIF 
       ENDIF
C  ENDIF FOR CELL(NXPRNT,NYPRNT)
       ENDIF
C  ENDIF fOR HGRID(NX,NY)=0.
c	if(iday.eq.15) then
c          IF(IHOUR.EQ.9) THEN
c            WRITE(27,'(2F8.3)') meltxy   
c	      stop
c          ENDIF
c	endif
c	IF(IHOUR.EQ.13) STOP
800     CONTINUE
900    CONTINUE
C       WRITE (*,*)'END OD DEM LOOP'
CC
C                   ******** END OF (DTM) SPATIAL LOOP **********
CC
C     Melt averaged over the basin area
Cgg	WRITE(28,*)RAIN,EVR*1000.,XYICEM(85,42)
	
c       WRITE(1,*)'ETPRS=',ETPRS
c       WRITE(1,*)'ETPAS=',ETPAS
c       WRITE(1,*)'ETPSUM=',ETPSUM
c	  Melt=Melt/FLOAT(NGRIDS)
	  POTEV=POTEV/FLOAT(NGRIDS)
C	 CALL PORTATA(MASK,SAT,XYMELT,PORT,HSNPCK,FHI,ihour,iday)
c	  XX(1,4,T)=MELT
CRRMUL	  XX(1,5,T)=PORT
CRRMUL	  DAYP=DAYP+PORT       
CRRMUL	  DAYD=DAYD+XX(1,6,t)*3.6/areaa
	  IF(ihour.EQ.24) THEN
CRRMUL		DAYP=DAYP/24.
CRRMUL		DAYD=DAYD/24.
		WRITE(26,545) month,iday,DAYD,DAYP,NSNOW
CRRMUL		DAYP=0.
CRRMUL		DAYD=0.
CRRMUL		NSNOW=0
	  ENDIF
545    format (2x,2(i2,2x),2(f10.6,2x),i8)
	  if (istamp.eq.-2) then 	
c karakoram 2.3          call wriBIN2(6,Gv5d,nr,nc,60)
c          call wriBIN2(8,Gv5d,nr,nc,80)
c karakoram 2.4          call wriASCF7(6,Gv5d,nr,nc,60)
c karakoram 2.4          call wriASCF7(8,Gv5d,nr,nc,80)
c karakoram 2.5          call wribin1(6,Gv5d,nr,nc,60)
c 80 Convective sensible heat
c 85 Latent heat
c 90 Surface temperature 
	     call wribin1(1,Gv5d,nr,nc,27)
          call wribin1(3,Gv5d,nr,nc,30)
		call wribin1(6,Gv5d,nr,nc,60)
		call wribin1(7,Gv5d,nr,nc,77)
		call wribin1(8,Gv5d,nr,nc,80)
		call wribin1(9,Gv5d,nr,nc,85)
		call wribin1(10,Gv5d,nr,nc,90)

	  endif 
cRRMUL
        idate=month*100+iday+ncyear*10000
CRRMUL	write(*,*) idate
        
	  if((idate.EQ.idatess(iim)).AND.(ihour.eq.9)) then
	  read(71,'(A8)') fileSCA(iim)
cRRMUL	write(*,'(2A8)') ' FileSC=',fileSCA(iim)
cRRMUL        open(70,file=fileSCA(iim),status='unknown',form='formatted')
        open(70,file=fileSCA(iim),status='unknown',form='binary')

c write a binary matrix SWE (mm) on perc1.out
    	 CALL WRIv5d (1,Gv5d,nr,nc,27)
CRRMUL scrive 1 matrice binaria (WEQ) su SWE,RST ex perc2.out and an ascii single file of SWE
CRRMUL CALL WRIv5d (1,Gv5d,nr,nc,28)
       CALL WRIvASC (2,Gv5d,nr,nc,28) 
CRRMUL scrive 1 matrice ascii (SCA=1) o ICE=2 
CRRMUL       CALL WRIvASC2 (1,Gv5d,luse,nr,nc,28) 
CRRMUL writes binary matrix  (SCA=1) or ICE=2 

           CALL WRIbin3 (1,Gv5d,luse,nr,nc,70) 
c scrive matrice binaria lwc su perc3.out
    	 CALL WRIv5d (3,Gv5d,nr,nc,30)
	write(*,*) ' Ciao Usman! I am storing image for',idatess(iim)
       iim=iim+1
	close(70)
	endif


CCC   ISTAMP=-1 Optional OUTPUT of MATRIX of WATER EQUIVALENT and MELT
	  IF(ISTAMP.EQ.-1) THEN
	    TX1=TIME(1,T)
	    TX2=TIME(2,T)
	    TX3=TIME(3,T)
	    TV11=6
	    TV12=29
	    TV21=6
	    TV22=29
	    TV31=0
	    TV32=0
c
Cggg	     IF(((TX1.EQ.TV11).AND.(TX2.EQ.TV12)).OR.
Cggg     +((TX1.EQ.TV21).AND.(TX2.EQ.TV22)).OR.
Cggg     +((TX1.EQ.TV31).AND.(TX2.EQ.TV32))) THEN
Cggg	      IF (TX3.EQ.12) THEN
Cggg	  scrivo gli if di controllo per far scrivere le mtarici di
Cggg     scioglimento per un intero evento
	 if(tx1.eq.tv11) then
Cgg	   if ((tx2.ge.tv12).and.(tx2.le.tv22)) then
CGG	      IF (((tx2.eq.tv12).and.(tx3.ge.18)).or.
CGG	+((tx2.eq.tv22).and.(tx3.le.18)))then
	      IF ((tx2.eq.tv12).and.(tx3.EQ.12))then
	ISTEP=1
CGG		istep=istep+1
Cgg	      WRITE(18,*) NYY,(TIME(J,T),J=1,3)
c	      CALL OUTPUT(NR,NC,XYMELT,istep,18)
CGG	      CALL OUTPUT(NR,NC,SWAVXY,istep,60)
CGG	      CALL OUTP1(NR,NC,RLOWXY,istep,61)
CGG	      CALL OUTP2(NR,NC,HCONVXY,istep,62)
CGG	      CALL OUTP3(NR,NC,HLATXY,istep,63)
	      CALL OUTP4(NR,NC,tsoil,istep,64)
Cgg	      CALL OUTPUT(NR,NC,hsnpck,istep)
c              WRITE(19,*) NYY,(TIME(J,T),J=1,3)
c              CALL OUTPUS(NR,NC,HSNPCK,HEABXY)
Cgg	    endif
	    ENDIF
	  ENDIF
C

C  OUTPUT OF METEOROLOGICAL INPUT DATA AND HEAT EXCHANGE COMPONENTS AT STATION
C
C     ISTAMP=-4 Optional OUTPUT of MATRIX of SWAVXY & STATISTICS
CCC               about RADIATION, HEAT BALANCE & MELT IN SPACE AND TIME
c       IN THIS VERSION THE HEAT BALANCE DISTRIBUTED MATRIX HEABXY IS 
C       SUBSTITUTED
C       BY THE WATER EQUIVALENT DISTRIBUTED MATRIX WEQUI
C       BUT THE STATISTICS OF SPACE AVERAGE (HEMUXY), SPACE STANDARD DEVIATION
C       (HESTXY) KEEP THE PREVIOUS NAME
C
      ELSEIF(ISTAMP.EQ.-4) THEN
CC
C spatial mean
CC
ccc         CALL STATXY(SWAVXY,HEABXY,XYMELT,MELT,RAMUXY,RASTXY,HEMUXY,
      CALL 
     +STATXY(WEQUI,XYMELT,MELT,XYICEM,PE,MELTICE,PEMED,PXY,PAREA,EVRXY,
     +SNOWXY, EVRMEDS, EVRMEDI, SNOWMEDS, RAINMEDS,RAINMEDI,MELTSNOW,
     +NXY,NSNOWXY,LUSE,NSNOWICE)

	  IF(NXY.LE.0) THEN
	   SWAVSG=0.
	   ELSEIF(HSNPCK(NXPRNT,NYPRNT).LE.0) THEN
	     SWAVSG=0.
	   ELSE 
	     SWAVSG=GminS0(Solrad,ALBEDS,FOR,ALBDFS)
	  ENDIF
CC
C
C		BILANCIO A SCALA DI BACINO
cccc
	  DO ISS=1,NUMSUB
	    IFS=21+ISS
       WRITE(IFS,996) NYY,(TIME(J,T),J=1,4),RAMUXY(ISS),RASTXY(ISS),
     +LONGM(ISS),LONGST(ISS),CONVM(ISS),CONVST(ISS),LATM(ISS),
     +LATST(ISS),RAINMEDI(ISS),EVRMEDI(ISS),MELTICE(ISS),
     +IFIX(NSNOWXY(ISS)),RAINMEDS(ISS),SNOWMEDS(ISS),EVRMEDS(ISS),
     +MELTSNOW(ISS),MELT(ISS),PEMED(ISS),PAREA(ISS),SWEM(ISS),
     +IFIX(NSNOWICE(ISS))
	  ENDDO
C       Input
C       WEQUI   Distributed water equivalent
C       XYMELT  Distributed snowmelt matrix
C       XYICEM  Distributed icemelt
C       PE	  Distributed net rainfall + snow/icemelt
C       PXY     Distributed net precipitation Rain-Evrpot+Snow 
C

C
C       Output
C
C       MELT    Spatial average Snowmelt
C	  MELTICE Spatial average Icemelt
C       PEMED   Spatial average net rainfall + snow/icemelt
c       PAREA   Spatial average of net precipitation
c       EVRXY   Distributed Evaporation
c       SNOWXY  Snowfall Precipitation
c       EVRMEDS Spatial Mean Evaporation over Snow 
c       EVRMEDI Spatial Mean Evaporation over snowfree area  
c       SNOWMEDS Spatial Average of snowfall over snow
c       RAINMEDS Spatial Average of rainfall over snow
c       RAINMEDI Spatial Average of rainfall over snowfree area
c       MELTSNOW Spatial Average of Snowmelt
c       NXY     # grids cells in the basin  
c       NSNOWXY # snowcover cells 
c 
c       Area Common   
c     
C       SWEM    Spatial mean Water EQUIVALENT over total basin area 
C       SWAVXY  Distributed Net shortwave radiation matrix
C       RAMUXY  Spatial mean Net shortwave radiation over area covered by snow
C       LONGM   Spatial mean Net longwave radiation over area covered by snow
C       CONVM   Spatial mean CONVECTIVE FLUX over area covered by snow
C       LATM    Spatial mean LATENT FLUX over area covered by snow
C       HEMUXY  Spatial mean HEAT BALANCE over area covered by snow
C       SWEST   Spatial standard deviation Water EQUIVALENT 
C       RASTXY  Spatial standard deviation Net shortwave radiation
C       LONGST  Spatial standard deviation Net longwave radiation
C       CONVST  Spatial standard deviation CONVECTIVE FLUX
C       LATST   Spatial standard deviation LATENT FLUX
C       HESTXY  Spatial standard deviation HEAT BALANCE 
C       MESTXY  Spatial standard deviation Melt over BASIN area
C       NXY     # grids in the basin 
C

996      FORMAT(2I4,2I3,I2,11(1X,F8.2),I9,8(1X,f8.2),I9)                
CRRMUL996      FORMAT(2I4,2I3,I2,11(1X,F8.1),1x,F6.1,1x,f8.1)                
995      FORMAT(2I4,2I3,I2,8(1X,E9.3),1x,I6)                
CCC	  XX(1,1,T)=RAMUXY
CRR20211226	  XX(1,2,T)=HEMUXY
	  XX(1,3,T)=SWAVSG
CC
C       remind that   XX(1,4,T)=MELT
CC
c         XX(1,5,T)=HEABAS   !!!!!ATTENZIONE,COMMENTO TEMPORANEO
C         XX(1,6,T)=smelts   !!!!!ATTENZIONE,COMMENTO TEMPORANEO 
	  ELSEIF(ISTAMP.EQ.-2) THEN
       WRITE(1,310) NYY,(TIME(J,T),J=1,3),P,Tasns,Solrad,Vw,Rhumi,
     +Heabal,wequi(nxsens,nysens),Melt(1),Gv5d(nxsens,nysens,10),DIS
       WRITE(7,797) NYY,(TIME(J,T),J=1,3),SWavbs,RLowas,Hconvs,
     &Hlates,Hpdts,Snogrs,SMelts,DIS,POTEV
310         FORMAT(I4,3I3,6F8.2,F8.3,F8.2,F6.1,F8.3)
797         FORMAT(I4,3I3,5F8.3,2F5.2,F8.3,F5.2)   
CCC
      ENDIF
CCC
1000  CONTINUE
c	write summer mass balance
      DO NX=1,NR
       DO NY=1,NC
	   if(luse(nx,ny).eq.1) then
			Gv5d(nx,ny,5)=Gv5d(nx,ny,5)+Gv5d(nx,ny,1)
	   else
		Gv5d(nx,ny,5)=-9999.0
	  
	 end do
	end do

      call wribin1(5,Gv5d,nr,nc,22) 


CCC
C                   ********+END OF TEMPORAL LOOP **********
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      
C  Scrive sul file TESIOUT.DAT  
C        Melt    Discharge     Potev      
C       
      IF(ISTAMP.EQ.3) THEN
	  DO NX=1,NR
	     WRITE (25,'(200(1x,I2))') (NESPO(NX,NY),NY=1,NC)
	  ENDDO
	  ENDIF
       WRITE(20,*) ' melt  discharge(PE)   port [mm] NSNOW'
       DO 1100 T=1,NTIMST
1100      WRITE(20,1110) XX(1,4,T),XX(1,6,T)*3.6/areaa,
     +XX(1,5,T),NSNOWXY
1110   FORMAT(1X,F5.2,2(1X,F10.3),1X,I10)
C
       IF(ISTAMP.EQ.-3) THEN
CCC
C    WATER SNOWPACK BALANCE at STATION CELL [meters of Water equivelent]
CCC
	  WEQUFI=WEQUIN+CUMWIN-CUMELT     
	  WRITE(9,1120) WEQUIN,CUMWIN,CUMELT,WEQUFI,WEQUI(NXPRNT,NYPRNT)
1120    FORMAT(//,' THEORICAL WATER BALANCE OF SNOWPACK [m]',//,
     &' WINITIAL + WATIN-EVAP - SNOWMELT= WFINAL     SIMULATED',
     &/,4F11.4,6X,F11.4)
       ENDIF
      IF(ISTAMP.EQ.-4) THEN
CCC
C       Temporal Mean and STD DEV about RADIATION, HEAT BALANCE & SNOWMELT
C       HERE ALSO HEAT BALANCE STATISTICS ARE SUBSTITUTED BY WATER EQUIVALENT
C       ONES :  HEMUT (TEMPORAL MEAN OF WEQUIV SPATIAL MEAN)
C               HESTT (TEMPORAL MEAN OF WEQUIV SPATIAL STANDARD DEVIATION)      
CCC
	  CALL STATT(NSTSWA,NSTRA,SWAMUT,SWASTT,RAMUT,RASTT,
     &HEMUT,HESTT,MEMUT,MESTT,HESMUT,HESSTT,MESMUT,MESSTT)
C 
	  WRITE(22,997) NSTSWA,SWAMUT,RAMUT,HESMUT,HEMUT,MESMUT,MEMUT,
     &SWASTT,RASTT,HESSTT,HESTT,MESSTT,MESTT

997     FORMAT(/,6X,'N.ore SWAMUT,RAMUT(snow on sensor) HESMUT ,  HEMUT'
     &'    MESMUT,  MEMUT',/,
     &'               [W/m²]                           [W/m²]     [mm]',
     &/,'Mean  ',2X,I4,2F6.1,13X,F8.1,F8.3,4X,2F7.2,
     &/,'St.Dev',6X,2F6.1,13X,F8.1,F8.3,4X,2F7.2)

      ENDIF
1001  IF(ISTAMP.EQ.-4)THEN
       CLOSE(22)
      ELSEIF(ISTAMP.EQ.-3)THEN
       CLOSE(8)
       CLOSE(9)

	   CLOSE(21)
      ELSEIF(ISTAMP.EQ.-2)THEN

       CLOSE(1)
       CLOSE(7)
       CLOSE(8)
       CLOSE(9)
       CLOSE(10)
      ELSEIF(ISTAMP.EQ.2)THEN
       CLOSE(10)
	  ENDIF

      CLOSE(20)

      CLOSE(26)
      RETURN
999   WRITE(*,*) ' Errore di lettura !' 
      RETURN
      END

C
C
C
	  SUBROUTINE WRITET(ISTAMP)
C
C  Scrive le intestazioni dei files di output
C        
	  RETURN
	  END
CCCCC
C     INITIA  This routine assignes parameters values read from the 
C             parameter file disma.inp to parameter variables    
C             without optimizing them                               
CCCCC
      SUBROUTINE INITIA(SMPLI,STPDTH)
CCCC
      PARAMETER (NNR=353,NNC=262)
	  CHARACTER*12 namefdem,namefdoc
CGG      CHARACTER*80 HEADER,HEADE1
      INTEGER (kind=2) ind(NNR,NNC),luse(nnr,nnc)
      INTEGER STPDTH,IUTMN(NNR,NNC),
     +IUTME(NNR,NNC),IDSUB(10),IMASK(10)
      REAL DEM(nnr,nnc),XNORM(3,nnr,nnc),VTMP(NNR,NNC)
Cgg    REAL SWST(2,nnr,nnc),PDKEF(2,nnr,nnc),PACT(2,nnr,nnc),tk(nnr,nnc),cmax(nnr,nnc)
      REAL SWST(2,nnr,nnc),PDKEF(2,nnr,nnc),PACT(2,nnr,nnc)
      REAL HSNPCK(nnr,nnc),DENSN(2,nnr,nnc),VVCFCT(12),ALBTM(nnr,nnc),
     +ALBSMOO(nnr,nnc)
	REAL KATM 
      LOGICAL MASK(nnr,nnc)       
CCCC
      COMMON /PERIOD/ NMONIN,NDAYIN,NMONFI,NDAYFI,ISTAMP,NCYEAR,NFYEAR,
     +NMONAN,NDAYAN
      COMMON /PARAME/ XS(20),X(20),INDX(20),N,NPAR,STEPS(20),TESTPS(20),
     +TESSTP(20),STEP(20)
      COMMON /HEAT/CPw,Wden,ROhu,Cpa,RatiKh,CKa,Z0,VLh,RatiKe,Rgas,ZW,ZT
      COMMON /SOIL/ ASOIL,SLAMDA,Tmsoil,Cv(NNR,NNC),DSOIL1(NNR,NNC)
      COMMON /SNOW/A0,A1,TSPRIN,ASPRIN,TRNCF,AM0,CAM,GANESN,ALFAS,
     +BETAS,STRN,CJOMME,HSNLA1,FHI,Swi,FHIeKs,BEERSN,CMELTM,SPCFCT
      COMMON /ATMOSP/ DOX,SATM0,TRAATM,PERCLO,VVCFCT,KATM
      COMMON /DEM/ DEM,XNORM,HSENSO,NXSENS,NYSENS,NR,NC,DEMDX,DEMDY,
     +AREAA,NGRIDS
      COMMON /SNPK/ SWST,PDKEF,PACT,HSNPCK,DENSN
      COMMON /MASKD/ MASK,NXPRNT,NYPRNT,IND,NUMSUB,IDSUB,IMASK
Cgg    COMMON /PORTV/IDEF,JDEF,CSTR,FHIT,cinv,TK,SKS,CMAX,ww
      COMMON /VEG/FOR,ALBTM,ALBSMOO,LUSE,h0c(nnr,nnc)
	  COMMON/rhid/NCMAT,NRMAT,DELTAMAT,AMINXMAT,AMAXXMAT,AMINYMAT,
     +AMAXYMAT
CCCC
      DATA CS/2093.4/
CCCC
C      WRITE(*,*) '                INPUT PARAMETERS FILE:'
C      WRITE(*,*)
C      WRITE(*,*) '            if disma.inp  type  RETURN,'
C      WRITE(*,*) '            if NEW FILE   type a 12 character name:' 
C      READ(*,'(A12)') INP
c
C        Open input file disma.inp or an input file named INP

C      IF(INP.EQ.BLANK) THEN
C       INP='DISMA.INP'
C      ENDIF  
      WRITE(*,*)  ' Apertura del file input'
      WATEQU=0.
      HMEAN=0.
      OPEN(2,FILE='DISMA.INP',STATUS='OLD',ERR=999)
      READ(2,'(A80)') HEADER 
	  WRITE(*,'(A80)') HEADER 
C dates of start of simulation and of start and end of animation
      READ(2,*) NMONIN,NDAYIN,NCYEAR,NMONAI,NDAYAI,NMONAF,NDAYAF 
      READ(2,*) NMONFI,NDAYFI,NFYEAR 
      WRITE(*,*) NMONIN,NDAYIN,NCYEAR 
      WRITE(*,*) NMONFI,NDAYFI,NFYEAR 
C         Max. number of parameters, Printing flag -4;-3;-2;-1;0;1 
      READ(2,*) NPAR,ISTAMP
C         Snow parameters header 
      READ(2,'(a80)') HEADE1
      READ(2,'(a80)') HEADER
CCC
C   Z0,A0,A1,TSPRIN,ASPRIN,TRNCF,AM0,CAM,STRN,ALFAS,BETAS,Swi
C   XS () PARAMETERS;  INDX() 0=Parameter Fixed 1=Parameter To be OPTIMIZED
C   STEPS() Dimension of Initial SIMPLEX; TESTPS() Simplex convergence Range
CCC
      READ(2,*) (XS(I),I=1,12)
      READ(2,*) (INDX(I),I=1,12)
      READ(2,*) (STEPS(I),I=1,12)
      READ(2,*) (TESTPS(I),I=1,12)
c      WRITE(*,*) HEADE1
c      WRITE(*,*) HEADER
c      WRITE(*,*) (XS(I),I=1,12)
C                               Soil parameters
C                               ASOIL,SLAMDA,Tmsoil 
      READ(2,'(a80)') HEADE1
      READ(2,'(a80)') HEADER
      READ(2,*) (XS(I),I=13,15)
      READ(2,*) (INDX(I),I=13,15)
      READ(2,*) (STEPS(I),I=13,15)
      READ(2,*) (TESTPS(I),I=13,15)
c      WRITE(*,*) HEADE1
c      WRITE(*,*) HEADER
c      WRITE(*,*) (XS(I),I=13,15)
C                               Atmosphere parameters
C                               SATM0,PERCLO,TRAATM,KATM
      READ(2,'(a80)') HEADE1
	READ(2,'(a80)') HEADER
	READ(2,*) (XS(I),I=16,19)
	READ(2,*) (INDX(I),I=16,19)
      READ(2,*) (STEPS(I),I=16,19)
	READ(2,*) (TESTPS(I),I=16,19)
c       WRITE(*,*) HEADE1
c      WRITE(*,*) HEADER
c      WRITE(*,*) (XS(I),I=16,19)
c           Other parameters
	  READ(2,'(a80)') HEADE1
	  READ(2,'(a80)') HEADER
	  READ(2,*) XS(20)
      READ(2,*) INDX(20)
      READ(2,*) STEPS(20)
      READ(2,*) TESTPS(20)
c      WRITE(*,*) HEADE1
c      WRITE(*,*) HEADER
c      WRITE(*,*) XS(20)
C                                SMPLI  STPDTH SPCKCF
      

	  
	  READ(2,'(a80)') HEADE1
      READ(2,'(a80)') HEADER
      
	  READ(2,*) SMPLI,STPDTH,SPCKCF
c      WRITE(*,*) HEADE1
c      WRITE(*,*) HEADER
c      WRITE(*,*) SMPLI,STPDTH,SPCKCF
CCCC
C
C       Digital Elevation Model
C
C       HSENSO  ELEVATION AT SENSOR STATION
C       NXSENS, NYSENS  LOCATION OF SENSOR STATION
C       NXPRNT, NYPRNT  LOCATION OF OUTPUT CELL
CCCC
      READ(2,'(a80)') HEADE1
      READ(2,'(a80)') HEADER
      READ(2,*) NR,NC,NGRIDS,DEMDX,DEMDY,AREAA
      READ(2,*) HSENSO,NXSENS,NYSENS,NXPRNT,NYPRNT
      WRITE(*,*) HEADE1
      WRITE(*,*) HEADER
      WRITE(*,*) NR,NC,NGRIDS,DEMDX,DEMDY,AREAA
      WRITE(*,*) HSENSO,NXSENS,NYSENS
C
	  namefdoc='  deminp.rdc'
	  namefdem='  deminp.rst'
	  call rhidri32 (namefdoc)
	  call ridin3r (namefdem,dem)
	  idmask=0
	  DO 100 IR=1,nrmat
      	DO IC=1,ncmat
	  IUTMN(IR,IC)=AMAXYMAT-DELTAMAT/2-DELTAMAT*(IR-1)
C	WRITE(*,*)'UTM-N=',IUTMN(IR,IC)
	  IUTME(IR,IC)=AMINXMAT+DELTAMAT/2+DELTAMAT*(IC-1)
C	WRITE(*,*)'UTM-E=',IUTME(IR,IC)
	       if ((ir.gt.1).and.(ir.lt.nrmat)) then
           	  if ((ic.gt.1).and.(ic.lt.ncmat)) then
      	         if (DEM(IR,IC).gt.0) idmask=idmask+1
	          endif
	       endif
	    ENDDO
100   CONTINUE
      write(*,*)'idmask=',idmask  
C
      CALL MASKDE(DEM,NR,NC)       
C
c      OPEN(52,FILE='HSNOW.INP',STATUS='OLD',ERR=999)
c      OPEN(53,FILE='PACT1.INP',STATUS='OLD',ERR=999)
c      OPEN(54,FILE='SWST1.INP',STATUS='OLD',ERR=999)
c      OPEN(55,FILE='SNOWDEN1.INP',STATUS='OLD',ERR=999)
c      OPEN(56,FILE='PACT2.INP',STATUS='OLD',ERR=999)
c      OPEN(57,FILE='SWST2.INP',STATUS='OLD',ERR=999)
c      OPEN(58,FILE='SNOWDEN2.INP',STATUS='OLD',ERR=999)
C
C       Snowpack initial state
C
C      READ(2,'(a80)') HEADE1
C      WRITE(*,*) HEADE1
C
      DO 120 NLAY=1,2
c
	  IF(NLAY.EQ.1) THEN
C
	namefdoc='  hsnow.rdc'
	namefdem='  hsnow.rst'
	call rhidri32 (namefdoc)
	call ridri32 (namefdem,vtmp)
	DO  nx=1,nrmat
      	DO ny=1,ncmat
CMUL
    	    hsnpck(nx,ny)=vtmp(nx,ny)
	    ENDDO
      enddo
	namefdoc='  pact1.rdc'
	namefdem='  pact1.rst'
	call rhidri32 (namefdoc)
	call ridin3r (namefdem,vtmp)
	DO  nx=1,nrmat
      	DO ny=1,ncmat
    	    pact(1,nx,ny)=vtmp(nx,ny)
	    ENDDO
      enddo

	namefdoc='  swst1.rdc'
	namefdem='  swst1.rst'
	call rhidri32 (namefdoc)
	call ridin3r (namefdem,vtmp)
	DO  nx=1,nrmat
      	DO ny=1,ncmat
    	    swst(1,nx,ny)=vtmp(nx,ny)
	    ENDDO
      enddo
	namefdoc='snowden1.rdc'
	namefdem='snowden1.rst'
	call rhidri32 (namefdoc)
	call ridin3r (namefdem,vtmp)
	DO  nx=1,nrmat
      	DO ny=1,ncmat
    	    densn(1,nx,ny)=vtmp(nx,ny)
	    ENDDO
      enddo
c	   NPUNIT=53
c	   NSUNIT=54
c	   NDUNIT=55
	  ELSE
c
	namefdoc='  pact2.rdc'
	namefdem='  pact2.rst'
	call rhidri32 (namefdoc)
	call ridin3r (namefdem,vtmp)
	DO  nx=1,nrmat
      	DO ny=1,ncmat
    	    pact(2,nx,ny)=vtmp(nx,ny)
	    ENDDO
      enddo
	namefdoc='  swst2.rdc'
	namefdem='  swst2.rst'
	call rhidri32 (namefdoc)
	call ridin3r (namefdem,vtmp)
	DO  nx=1,nrmat
      	DO ny=1,ncmat
    	    swst(2,nx,ny)=vtmp(nx,ny)
	    ENDDO
      enddo
	namefdoc='snowden2.rdc'
	namefdem='snowden2.rst'
	call rhidri32 (namefdoc)
	call ridin3r (namefdem,vtmp)
	DO  nx=1,nrmat
      	DO ny=1,ncmat
    	    densn(2,nx,ny)=vtmp(nx,ny)
	    ENDDO
      enddo

	  ENDIF

	  WRITE(*,*) ' INITIAL HEAT STORAGE'

134     DO 125 NY=1,NC
	  DO 130 NXPDK=1,NR
	  IF(NLAY.EQ.1) THEN
C  Calibration of initial state by multiplying HSNPCK by SPCKCF
	    HSNPCK(NXPDK,NY)=AMAX1(HSNPCK(NXPDK,NY)*SPCKCF,0.)
C
	    HSNOW=AMIN1(HSNPCK(NXPDK,NY),HSNLA1)
	    HMEAN=HMEAN+DEM(NXPDk+1,NY+1)
	    IF(HSNPCK(NXPDK,NY).GT.0) 
     +WATEQU=WATEQU+HSNOW*DENSN(NLAY,NXPDK,NY)
	  ELSE
	    HSNOW=AMAX1(HSNPCK(NXPDK,NY)-HSNLA1,0.)   
	    WATEQU=WATEQU+HSNOW*DENSN(NLAY,NXPDK,NY)
	  ENDIF
	    PDKEF(NLAY,NXPDK,NY)=CS*DENSN(NLAY,NXPDK,NY)*
     +HSNOW*PACT(NLAY,NXPDK,NY)
C          WRITE(*,*) 'Contenuto di freddo',PDKEF(NLAY,NXPDK,NY),
c     +nlay,nxpdk,ny
130       CONTINUE
125     CONTINUE 
120   CONTINUE
	HMEAN=HMEAN/FLOAT(NGRIDS)
	VOLEQU=WATEQU/WDEN*AREAA/FLOAT(NGRIDS)
	WATEQU=VOLEQU/AREAA
	   WRITE(*,*)' Altezza media bacino= ',HMEAN,' m.s.l.m.'
    	WRITE(*,*)' Volume Eq. in acqua neve = ',VOLEQU,'x 10^6 [m^3]'
     	WRITE(*,*)' Equivalente in acqua neve= ',WATEQU,'       [m]'
C
C       Measurements CORRECTION Factors
c       SPCFCT=Snow Precipitation Correction Factor VVCFCT(12) Wind Velocity
	  READ(2,'(a80)') HEADER
	  READ(2,*) SPCFCT
	  READ(2,*) (VVCFCT(IV),IV=1,12)
c	  OPEN(30,FILE='LUSE.INP',STATUS='OLD',ERR=999)
c	  OPEN(31,FILE='ALB.INP',STATUS='OLD',ERR=999)
c	  OPEN(40,FILE='SOILCMAX.INP',STATUS='OLD',ERR=999)
c	  OPEN(41,FILE='SOILTK.INP',STATUS='OLD',ERR=999)
c	 DO NX=1,NR
c	  do ny=1,nc
c	    READ(30,*) LUSE(NX,NY)
c	  enddo
c	 ENDDO    
	namefdoc=' luseinp.rdc'
	namefdem=' luseinp.rst'
	call rhidri32 (namefdoc)
	call ridin32 (namefdem,luse)
c
c    landuse classification:
c	0  = Shadow/no DTM values
c	1  = Snow/Glacier
c	2  = Settlements
c	3  = Water
c	4  = Grassland
c	5  = Agriculture
c	6  = Coniferous
c	7  = Rock
c	8  = Gravel
c	9  = Larches
c	10 = Deciduous
c	11 = Mixed Forest
c	12 = Alpine Meadows
c	13 = Alder
c	14 = vegetazione boschiva ed arbustiva in evoluzione
c    15 vegetazione rada
c
   	do nx=1,nr
	  do ny=1,nc
      landuse=luse(nx,ny)
	IF((LANDUSE.EQ.7).OR.(LANDUSE.LE.3)) THEN
cggrrli vedi hpsuolo.xls CV J/m3/K roccia e suolo
	   CV(NX,NY)=2160000.
	   DSOIL1(NX,NY)=0.43
	 ELSE
	   CV(NX,NY)=2389500.
	   DSOIL1(NX,NY)=0.24
	ENDIF
	if(landuse.eq.4) then
	 h0c(nx,ny)=0.12
	elseif (landuse.eq.5) then
	 h0c(nx,ny)=0.5
	elseif ((landuse.eq.6).or.(landuse.eq.11)) then
	 h0c(nx,ny)=7.5
	elseif ((landuse.eq.9).or.(landuse.eq.10)) then
	 h0c(nx,ny)=5.0
CRRAG	20070831 V01.02
	elseif ((landuse.eq.12).or.(landuse.eq.13).or.(landuse.eq.14).or.
     +(landuse.eq.15)) then
	 h0c(nx,ny)=0.5
	else
	 h0c(nx,ny)=0.
	endif
	enddo
	enddo
Cggg	reading resampled Landsat derived albedo [multiplied by 100]
	namefdoc='  albinp.rdc'
	namefdem='  albinp.rst'
	call rhidri32 (namefdoc)
	call ridin3r (namefdem,albtm)
Cggg	reading smoothed (750 m or 1250 m) albedo [multiplied by 100]
	namefdoc=' albsmoo.rdc'
	namefdem=' albsmoo.rst'
	call rhidri32 (namefdoc)
	call ridri32 (namefdem,albsmoo)
Cggg   reading soil maximum capacity
CGGG	namefdoc='soilcmax.rdc'
CGGG	namefdem='soilcmax.rst'
CGGG	call rhidri32 (namefdoc)
CGGG	call ridri32 (namefdem,cmax)
Cggg   reading soil permeability
CGGG	namefdoc='  soiltk.rdc'
CGGG	namefdem='  soiltk.rst'
CGGG	call rhidri32 (namefdoc)
CGGG	call ridri32 (namefdem,tk)

	 write(*,*)'END OF SOIL INPUT READING'   
C
C     Parameters are set up just once for simulation
C
      CALL SETUP
C
      RETURN
999   WRITE(*,*) ' Error opening FILE input'
      STOP      
      END
C
CCC
C
	SUBROUTINE MASKDE(DEM,NR,NC)
	PARAMETER (NNR=353,NNC=262,nbas=10)
	LOGICAL MASK(nnr,nnc)
	INTEGER (kind=2) IND(nnr,nnc)
	INTEGER IDSUB(nbas),IMASK(nbas)
	REAL DEM(nnr,nnc)
      COMMON /MASKD/ MASK,NXPRNT,NYPRNT,IND,NUMSUB,IDSUB,IMASK
	DATA  IMASK/10*0/
C	
      IMASK0=0
C	IMASK1=0
C	IMASK19=0
C	IMASK20=0
C	IMASK21=0
C	IMASK22=0
C	IMASK23=0
C	IMASK24=0
C	IMASK25=0
C	IMASK26=0

      OPEN (50,FILE='MASK.rst', form='binary')
c	read (50)	 IDRISI32 file (*.rst)
        DO I=1,NR
		do j=1,nc
            READ(50) IND(I,j)
	  enddo
        ENDDO
	DO 100 I=1,NR
	DO 200 J=1,NC
		MASK(I,J)=.FALSE.
          IF ((IND(I,J).le.0).or.(dem(i,j).le.0)) THEN
		IMASK0=IMASK0+1
	  ELSE
	     DO ISS=1,NUMSUB
 	        NUM=IDSUB(ISS)
			if(ind(i,j).eq.NUM) then
		     MASK(I,J)=.TRUE.
		     IMASK(iss)=IMASK(iss)+1
	     GOTO 200
	        ENDIF
	     ENDDO
	 ENDIF

200     CONTINUE
100     CONTINUE
	close(50)
	WRITE(*,*) ' NUMERO CELLE fuori bacino= ', IMASK0
       DO ISS=1,NUMSUB
	    NUM=IDSUB(ISS)
	 WRITE(*,*) ' NUMERO CELLE DEL BACINO ',NUM,' = ', IMASK(iss)
	ENDDO
C	WRITE(*,*) ' NUMERO CELLE DEL BACINO 1= ', IMASK1
C	WRITE(*,*) ' NUMERO CELLE DEL BACINO 19= ', IMASK19
C	WRITE(*,*) ' NUMERO CELLE DEL BACINO 20= ', IMASK20
C	WRITE(*,*) ' NUMERO CELLE DEL BACINO 21= ', IMASK21
C	WRITE(*,*) ' NUMERO CELLE DEL BACINO 22= ', IMASK22
C	WRITE(*,*) ' NUMERO CELLE DEL BACINO 23= ', IMASK23
C	WRITE(*,*) ' NUMERO CELLE DEL BACINO 24= ', IMASK24
C	WRITE(*,*) ' NUMERO CELLE DEL BACINO 25= ', IMASK25
C	WRITE(*,*) ' NUMERO CELLE DEL BACINO 26= ', IMASK26
	RETURN
	END
C
C       SUBROUTINE DTM(NX,NY,THETR,PSIZ,SOLBEA,SOLBXY,HGRID,PROJEC)
C
C       Input
C
C       NX      = X coordinate position of actual square    
C       NY      = Y coordinate position of actual square
C       THETR   = SUN ELEVATION                         [rad]
C       PSIZ    = SUN AZIMUT (respect to NORTH)         [rad]      
C       SolBEA  =Direct beam radiation [W/m²] on sensor station grid
C
C       Output
C       SOLBXY  = Direct beam radiation on DTM grid     [W/m²]
C       HGRID   = Grid square mean elevation            [m]
C       PROJEC  = COSINE PROJECTION
C
C       Common input
C       DEM(x,y)= ARRAY of elevations                   [m]
C       Note:     External [DEM(1,i), DEM(NC,i), DEM(i,1) & DEM(i,NR) ]
C                 elevations ane needed only to compute topographic gradient 
C                 but are not part of actual squares
C                 Example: Upper left corner of basin = DEM(2,2)
C                          but NX=1, NY=1
C                      or  DEM(4,5) -> NX=3, NY=4 ....
C    
C       DEMDX      = X Grid spacing  (North-South)        [m]
C       DEMDY      = Y Grid spacing  (West-East)      [m]
C       NC         = n. of grid cells in Y coordinate (number of columns) [-]
C       NR         = n. of grid cells in X coordinate (number of rows) [-]
C
		   
	SUBROUTINE DTM(NX,NY,THETR,PSIZ,SOLBEA,SOLBXY,HGRID,PROJEC,SATM)
	PARAMETER (NNR=353,nnc=262)
	REAL DEM(nnr,nnc),XNORM(3,nnr,nnc)
	REAL VVCFCT(12),KATM
	INTEGER XCROSS,YCROSS
	COMMON /DEM/ DEM,XNORM,HSENSO,NXSENS,NYSENS,NR,NC,DEMDX,DEMDY,
     +AREAA,NGRIDS
       COMMON /SUNS/ SCI,SCE,ANRAD,RADAN,PI,PI12,PI2,PASSO,SOLST,CALJOU
      COMMON /ATMOSP/ DOX,SATM0,TRAATM,PERCLO,VVCFCT,KATM
C       OLD STATEMENTS, EQUIVALENT TO MASK(I,J)
CC          HGRID=0.
C          IF(DEM(NX,NY).LE.0) RETURN   
C          IF(DEM(NX+2,NY).LE.0) RETURN    
C          IF(DEM(NX,NY+2).LE.0) RETURN   
CC          IF(DEM(NX+2,NY+2).LE.0) RETURN   
C          WRITE (*,*) 'DTM'
	  HGRID=DEM(NX+1,NY+1)
C          WRITE(*,*) HGRID
CG!!!	  IF(THETR.LE.0) THEN
	  IF(THETR.LE.0.1) THEN
C    -2 Flag: SUN IS NOT RISEN ABOVE THE HORIZON 
	    PROJEC=-2. 
	    SOLBXY=0.
	    RETURN
	  ELSEIF(SOLBEA.LE.0) THEN
C    -3 Flag: SENSOR MEASURES 0 VALUE OF RADIATION (night)
	    PROJEC=-3. 
	    SOLBXY=0.
	    RETURN
	  ENDIF
	TINTEX=0.
	TINTEY=0.
	PSIZEA=PSIZ-PI2
C	Da qui in poi inverto tutte le x con le y!!!!!

C	DTINTX=ABS(DEMDX/COS(PSIZEA))!!!
C	DTINTY=ABS(DEMDY/SIN(PSIZEA))!!!
	DTINTX=ABS(DEMDX/SIN(PSIZEA))
	DTINTY=ABS(DEMDY/COS(PSIZEA))
C
C XCROSS YCROSS Maximum possible number of X-lines and Y-lines crossings
C    NXINT  N. of crossed lines counter   
C
	  NYINT=NY+1

C
	IF(PSIZ.GT.PI) THEN
	  YCROSS=NY
	  NYSGN=-1
	  IF(PSIZEA.GT.PI)THEN
	    XCROSS=NX
	    NXSGN=-1 
	  ELSE
	    NXSGN=1
	    XCROSS=NR+1-NX 
	  ENDIF
	ELSE
	  YCROSS=NC+1-NY
	  NYSGN=1
	  IF(PSIZEA.LT.0) THEN
	    XCROSS=NX
	    NXSGN=-1 
	  ELSE
	    NXSGN=1
	    XCROSS=NR+1-NX 
	  ENDIF
	ENDIF 
C
C       X0,Y0 Distance from grid origin of X-lines and Y-lines crossings 
C       TINTEX,TINTEX Distance of grid crossing from current grid center
C
	X0=DEMDX*(FLOAT(NX))
	Y0=DEMDY*(FLOAT(NY))
	DO 100 NYCROS=1,YCROSS
	    TINTEY=TINTEY+DTINTY
C    X-Y Coordinates of X-lines intersections
	  XINTER=X0+TINTEY*SIN(PSIZEA)
	  NYINT=NYINT+NYSGN
	  NXINT=INT(XINTER/DEMDX)+1
C    Terrain elevation at intersection
	  HINTER=HINTY(NXINT,XINTER,NYINT)
	  HSBEAM=HGRID+TINTEY*TAN(THETR)
	  IF(HSBEAM.LT.HINTER) THEN
C   SHADOW
	    SOLBXY=0.
	    PROJEC=-1.
	    RETURN
	  ENDIF         
100     CONTINUE
C
	  NXINT=NX+1
C
	DO 200 NXCROS=1,XCROSS
	    TINTEX=TINTEX+DTINTX
C    X-Y Coordinates of Y-lines intersections
	  YINTER=Y0+TINTEX*COS(PSIZEA)
	  NYINT=INT(YINTER/DEMDY)+1
	  NXINT=NXINT+NXSGN
	  HINTER=HINTX(NXINT,YINTER,NYINT)
	  HSBEAM=HGRID+TINTEX*TAN(THETR)
	  IF(HSBEAM.LT.HINTER) THEN
C   SHADOW
	    SOLBXY=0.
	    PROJEC=-1.5
	    RETURN
	  ENDIF
200     CONTINUE
C                               Sun-grid vector
	    COSTHE=COS(THETR)
C	    XSUNVE=COSTHE*SIN(PSIZEA)!!!!
C	    YSUNVE=COSTHE*COS(PSIZEA)!!!!
	    XSUNVE=COSTHE*COS(PSIZEA)
	    YSUNVE=COSTHE*SIN(PSIZEA)
	    ZSUNVE=SIN(THETR)
C
C    COSINE (Scalar product) between Normal vector and Sun vector
C
	  PROJEC=XNORM(2,NX,NY)*XSUNVE+XNORM(1,NX,NY)*YSUNVE+
     +XNORM(3,NX,NY)*ZSUNVE
C   Approximation for direct beam radiation !!!
	  IF(PROJEC.LT.0.) THEN
	     SOLBXY=0.
	     RETURN
	  ELSE
C Beam radiation is attenuated as a function of the altitude
	     COEFF=EXP(-1./SIN(THETR)*(SATM0*EXP(-HGRID/KATM)-SATM))
		 SOLBXY=(SOLBEA*COEFF/ZSUNVE)*PROJEC
	     RETURN
	  ENDIF 
	  END
C
C       Computes elevation of sun-beam projection on terrain at intersection
C       NXINT=nr. of DTM ROW of segment to be intersected 
C       NYINT=nr. of DTM COLUMN of upper part of segment to be intersected 
C       YINTER Y-coordinate of point at intersection
C 
	FUNCTION HINTX(NXINT,YINTER,NYINT)
	PARAMETER (NNR=353,NNC=262)
	REAL DEM(nnr,nnc),XNORM(3,nnr,nnc)
	COMMON /DEM/ DEM,XNORM,HSENSO,NXSENS,NYSENS,NR,NC,DEMDX,DEMDY,
     +AREAA,NGRIDS
	IF((NXINT.LT.1).OR.(NXINT.GT.NR+2).OR.
     &(NYINT.LT.1).OR.(NYINT.GT.NC+1)) THEN
	  HINTX=0.
	  RETURN
	ENDIF 
	DH=YINTER/DEMDY-FLOAT(NYINT-1)
	NYINP1=NYINT+1
	HINTX=DEM(NXINT,NYINT)+(DEM(NXINT,NYINP1)-DEM(NXINT,NYINT))*DH
	RETURN
	END
C
C       NXINT=nr. of DTM ROW of left part of segment to be intersected 
C       NYINT=nr. of DTM COLUMN of segment to be intersected 
C       XINTER X-coordinate of point at intersection
C 
	FUNCTION HINTY(NXINT,XINTER,NYINT)
	PARAMETER (NNR=353,NNC=262)
	REAL DEM(nnr,nnc),XNORM(3,nnr,nnc)
	COMMON /DEM/ DEM,XNORM,HSENSO,NXSENS,NYSENS,NR,NC,DEMDX,DEMDY,
     +AREAA,NGRIDS        
	IF((NXINT.LT.1).OR.(NXINT.GT.NR+1).OR.
     &(NYINT.LT.1).OR.(NYINT.GT.NC+2)) THEN
	  HINTY=0.
	  RETURN
 	ENDIF 
	DH=XINTER/DEMDX-FLOAT(NXINT-1)
	NXINP1=NXINT+1
	HINTY=DEM(NXINT,NYINT)+(DEM(NXINP1,NYINT)-DEM(NXINT,NYINT))*DH
	RETURN
	END
C
C    VNORM Computes Normal vector respect to terrain surface
c
	  SUBROUTINE VNORM(NX,NY)
	PARAMETER (NNR=353,NNC=262)
	REAL DEM(nnr,nnc),XNORM(3,nnr,nnc)
	COMMON /DEM/ DEM,XNORM,HSENSO,NXSENS,NYSENS,NR,NC,DEMDX,DEMDY,
     +AREAA,NGRIDS          
	  DEM11=DEM(NX,NY)
	  IF(DEM11.LE.0) RETURN   
	  DEM13=DEM(NX+2,NY)
	  IF(DEM13.LE.0) RETURN    
	  DEM31=DEM(NX,NY+2)
	  IF(DEM31.LE.0) RETURN   
	  DEM33=DEM(NX+2,NY+2)
	  IF(DEM33.LE.0) RETURN   
	  DEM12=DEM(NX+1,NY)
	  DEM32=DEM(NX+1,NY+2)
	 XNORMA=-Gy(DEM11,DEM12,DEM13,DEM31,DEM32,DEM33,DEMDY)
	  DEM21=DEM(NX,NY+1)
	  DEM23=DEM(NX+2,NY+1)
	 YNORMA=-Gx(DEM11,DEM13,DEM21,DEM23,DEM31,DEM33,DEMDX)
	 ZNORMA=1.
	  VECNOR=SQRT(xnorma*xnorma+ynorma*ynorma+znorma*znorma)
	  XNORM(1,NX,NY)=xnorma/vecnor
	  XNORM(2,NX,NY)=ynorma/vecnor
	  XNORM(3,NX,NY)=znorma/vecnor
	  RETURN
  	END
C
C       Gx and Gy Sobel algorithm for x,y gradient
C
	FUNCTION Gy(DEM11,DEM12,DEM13,DEM31,DEM32,DEM33,Dy)
	  Gy=(DEM31+2*(DEM32-DEM12)+DEM33-DEM11-DEM13)/(8*DY)
	RETURN
	END
C
	FUNCTION Gx(DEM11,DEM13,DEM21,DEM23,DEM31,DEM33,DX)
	  Gx=(DEM13+2*(DEM23-DEM21)+DEM33-DEM11-DEM31)/(8*DX)
	RETURN
	END
C
C
C
	    SUBROUTINE ESPO(E1,E2,ETANG,NESPO)
c
C   Calcola l'esposizione per ogni cella
C            
      IF (E2.EQ.0) THEN
	  IF (E1.GT.0) THEN
	      NESPO=1
	   ELSEIF (E1.LT.0) THEN
	      NESPO=5
	   ELSE
	      NESPO=99
	   ENDIF
       ELSEIF (E1.EQ.0) THEN
	  IF (E2.GT.0) THEN
	      NESPO=3
	   ELSEIF (E2.LT.0) THEN
	      NESPO=7
	   ELSE
	      NESPO=99
	  ENDIF
       ELSE
	  ETANG=E2/E1
	  IF ((ETANG.GT.-0.4142).AND.(ETANG.LT.0.4142)) THEN
	       IF (E1.GT.0) THEN
		   NESPO=1
		ELSE
		   NESPO=5
	       ENDIF
	   ELSEIF ((ETANG.GT.2.4142).OR.(ETANG.LT.-2.4142)) THEN
	       IF (E2.GT.0) THEN
		   NESPO=3
		ELSE
		   NESPO=7
	       ENDIF
	   ELSEIF (E2.GT.0) THEN
	       IF (E1.GT.0) THEN
		   NESPO=2
		ELSE
		   NESPO=4
	       ENDIF
	   ELSEIF (E1.GT.0) THEN
	       NESPO=8
	   ELSE
	       NESPO=6
	   ENDIF   
      ENDIF
      RETURN
      END
C

	subroutine assmeteo(t,ncheck,NUMSUB,idsub,P,TASNS,VW,RHUMI,PMED,
     +hsensp,hsenst,hspmis,hstmis,numsrif,nsensp,nsenst,isp,jsp,ist,jst)
      integer T,idsub(10),hsenst(10),hsensp(10)
	integer nsenst(2,10),nsensp(2,10)
      COMMON /RECORD/ XX(10,6,8784),NTIMST
c
	ifound=0
      do isub=1,numsub
	   nid=idsub(isub)
          if(ncheck.eq.nid) then
             ifound=1
	       P=XX(isub,1,T)
		   isp=nsensp(1,isub)
		   jsp=nsensp(2,isub)
 		   hspmis=hsensp (isub)
		   if (P.eq.-99.9) then
		      P=XX(numsrif,1,T)
		      isp=nsensp(1,numsrif)
		      jsp=nsensp(2,numsrif)
 		      hspmis=hsensp(numsrif)
		   endif               
             Tasns=XX(isub,2,T) 
		   ist=nsenst(1,isub)
		   jst=nsenst(2,isub)
      	   hstmis=hsenst(isub)
		   if (Tasns.eq.-99.9) then
		   Tasns=XX(numsrif,2,T)
 		   ist=nsenst(1,numsrif)
 		   jst=nsenst(2,numsrif)
  		   hstmis=hsenst(numsrif)
		   endif
             Vw=XX(isub,5,T)
		   if (Vw.eq.(-99.9)) then
		      Vw=XX(isub,5,T-1)
	 			  if (Vw.eq.(-99.9)) Vw=1.
             endif				               
C
            
	        Rhumi=XX(isub,3,T)
		    if (Rhumi.eq.-99.9) Rhumi=XX(numsrif,3,T)             
            
	       PMED=XX(isub,4,T)
	if (pmed.eq.-99.9) pmed=1013.
c               pressure in mbar
	    endif
	enddo
	if(ifound.eq.0) then
		  write(*,*)'maschera vuota ncheck=',ncheck
	      stop
	endif

	return
	end
ccc
cc
cc	 
	subroutine dtweight (ast,difrad,gammat)
Cgg Valori per MAP-SOP 99
Cgg	DATA GAMMACC /-0.0051/
Cgg	DATA tm,hv,d1,d2 /-0.00574,-0.00153,9.,15./
Cgg Valori per giugno 1997
	DATA GAMMACC /-0.0052/
	DATA tm,hv,d1,d2 /-0.0065,-0.00195,9.,15./
	DATA fi1,fi2 /5.,14./
	if (ast.lt.fi1) then
	  gammacl=tm+hv*(cos((ast+24-fi2)*4.*atan(1.)/d2)-
     +0.5*sin((ast+24-fi2)*4.*atan(1.)/d2)*
     +(1+cos((ast+24-fi2)*4.*atan(1.)/d2)))
	elseif (ast.gt.fi2) then
	  gammacl=tm+hv*(cos((ast-fi2)*4.*atan(1.)/d2)-
     +0.5*sin((ast-fi2)*4.*atan(1.)/d2)*
     +(1+cos((ast-fi2)*4.*atan(1.)/d2)))
	else
	  gammacl=tm-hv*cos((ast-fi1)*4.*atan(1.)/d1)
	endif
	GAMMAT=DIFRAD*GAMMACC+(1.-DIFRAD)*GAMMACL
	return
	end

C
C  Subroutine RHIDRISI reads an IDRISI Document file
C
	subroutine rhidrisi (namefdoc)
	character*12 namefdoc
	COMMON/rhid/NCMAT,NRMAT,DELTAMAT,AMINXMAT,AMAXXMAT,AMINYMAT,
     +AMAXYMAT

      OPEN(51,FILE=namefdoc,STATUS='OLD')
	WRITE(*,'(1X,12A)')namefdoc
      READ(51,*)
	READ(51,*)
	READ(51,*)
	READ(51,'(14X,I4)') NCMAT
	READ(51,'(14X,I4)') NRMAT
	WRITE(*,*)'NC,NR=',NCMAT,NRMAT
	READ(51,*)
	READ(51,*)
	READ(51,'(14X,F10.0)') DELTAMAT
	READ(51,'(14X,F10.0)') AMINXMAT
	READ(51,'(14X,F10.0)') AMAXXMAT
	READ(51,'(14X,F10.0)') AMINYMAT
	READ(51,'(14X,F10.0)') AMAXYMAT
	WRITE(*,*)DELTAMAT
	WRITE(*,*)AMINXMAT,AMAXXMAT,AMINYMAT,AMAXYMAT
	CLOSE(51)

	RETURN 
	END


C
C  Subroutine RIDRISI reads an IDRISI Image file (real numbers)
C
	subroutine ridrisi (namefdem,varf)
	parameter(NNR=353,NNC=262)
	character*12 namefdem
	real varf(nnr,nnc)
 	COMMON/rhid/NCMAT,NRMAT,DELTAMAT,AMINXMAT,AMAXXMAT,AMINYMAT,
     +AMAXYMAT

      OPEN(51,FILE=namefdem,STATUS='OLD')
	write(*,'(1x,12A)') namefdem
	DO  IR=1,nrmat
      	DO IC=1,ncmat
		   READ(51,*) varf(IR,IC)
	    ENDDO
      enddo
	CLOSE (51)

	RETURN 
	END

C
C  Subroutine RIDRint reads an IDRISI Image file (integer numbers)
C
	subroutine ridrint (namefdem,nvarf)
	parameter(NNR=353,NNC=262)
	character*12 namefdem
	integer nvarf(nnr,nnc)
 	COMMON/rhid/NCMAT,NRMAT,DELTAMAT,AMINXMAT,AMAXXMAT,AMINYMAT,
     +AMAXYMAT

      OPEN(51,FILE=namefdem,STATUS='OLD')
	write(*,'(1x,12A)') namefdem
	DO  IR=1,nrmat
      	DO IC=1,ncmat
		   READ(51,*) nvarf(IR,IC)
	    ENDDO
      enddo
	CLOSE (51)

	RETURN 
	END

C
C  Subroutine RHIDRI32 reads an IDRISI32 Document file (*rdc)
C
	subroutine rhidri32 (namefdoc)
	character*12 namefdoc
	COMMON/rhid/NCMAT,NRMAT,DELTAMAT,AMINXMAT,AMAXXMAT,AMINYMAT,
     +AMAXYMAT

      OPEN(51,FILE=namefdoc,STATUS='OLD')
	WRITE(*,'(1X,12A)')namefdoc
      READ(51,*)
c
      READ(51,*)
	READ(51,*)
	READ(51,*)
	READ(51,'(14X,I4)') NCMAT
	READ(51,'(14X,I4)') NRMAT
	WRITE(*,*)'NC,NR=',NCMAT,NRMAT
	READ(51,*)
	READ(51,*)
	READ(51,'(14X,F10.0)') DELTAMAT
	READ(51,'(14X,F10.0)') AMINXMAT
	READ(51,'(14X,F10.0)') AMAXXMAT
	READ(51,'(14X,F10.0)') AMINYMAT
	READ(51,'(14X,F10.0)') AMAXYMAT
	WRITE(*,*)DELTAMAT
	WRITE(*,*)AMINXMAT,AMAXXMAT,AMINYMAT,AMAXYMAT
	CLOSE(51)

	RETURN 
	END


C
C  Subroutine RIDRI32 reads an IDRISI32 Image file *.rst binary (real numbers)
C
	subroutine ridri32 (namefdem,varf)
	parameter(NNR=353,NNC=262)
	character*12 namefdem
	real varf(nnr,nnc)
 	COMMON/rhid/NCMAT,NRMAT,DELTAMAT,AMINXMAT,AMAXXMAT,AMINYMAT,
     +AMAXYMAT

      OPEN(51,FILE=namefdem,FORM='binary')
	write(*,'(1x,12A)') namefdem
	DO  IR=1,nrmat
      	DO IC=1,ncmat
		   READ(51) varf(IR,IC)
	    ENDDO
      enddo
	CLOSE (51)

	RETURN 
	END

C
C  Subroutine RIDin32 reads an IDRISI32 Image file (*.rst)(integer numbers, converted to real)
C
	subroutine ridin3r (namefdem,varf)
	parameter(NNR=353,NNC=262)
	character*12 namefdem
	integer (kind=2) nvarf(nnr,nnc)
	real varf(nnr,nnc)
 	COMMON/rhid/NCMAT,NRMAT,DELTAMAT,AMINXMAT,AMAXXMAT,AMINYMAT,
     +AMAXYMAT

      OPEN(51,FILE=namefdem,STATUS='OLD',FORM='BINARY')
	write(*,'(1x,12A)') namefdem
	DO  IR=1,nrmat
      	DO IC=1,ncmat
		   READ(51) nvarf(IR,IC)
	varf(ir,ic)=nvarf(ir,ic)
	    ENDDO
      enddo
	CLOSE (51)

	RETURN 
	END

C  Subroutine RIDin32 reads an IDRISI32 Image file (*.rst)(integer numbers)
C
	subroutine ridin32 (namefdem,nvarf)
	parameter(NNR=353,NNC=262)
	character*12 namefdem
	integer (kind=2) nvarf(nnr,nnc)
 	COMMON/rhid/NCMAT,NRMAT,DELTAMAT,AMINXMAT,AMAXXMAT,AMINYMAT,
     +AMAXYMAT

      OPEN(51,FILE=namefdem,STATUS='OLD',FORM='BINARY')
	write(*,'(1x,12A)') namefdem
	DO  IR=1,nrmat
      	DO IC=1,ncmat
		   READ(51) nvarf(IR,IC)
	    ENDDO
      enddo
	CLOSE (51)

	RETURN 
	END
C       NELLA SOUBROUTINE RAINR VENGONO LETTE LE COORDINATE DEI
C       PLUVIOMETRI E NELLA SUB. THISSEN VENGONO CALCOLATI I
C       COEFFICIENTI DI THISSEN
C       VENGONO INOLTRE LETTI I VALORI DI AMCC RELATIVI A CIASCUN
C       PLUVIOGRAFO E CON I PESI DI THISSEN VIENE ASSEGNATO IL
C       VALORE CORRISPONDENTE AD OGNI CELLA
CC
      SUBROUTINE RAINR(IDIM,JDIM,IDCB,NSTR,IMISUR,XMIN,YMAX,DX,
     +ISTART,MASK)
CGGG Passo anche IMISUR per poter controllare la validita' del dato
CGGG di pioggia
      PARAMETER (NOBS=5140,NSTRMX=900,NNR=353,NNC=262)
C      COMMON /COM/ Z(NNR,NNC),ZZ(NNR,NNC),MAT(NNR,NNC),Q(NNR,NNC),
C     .QQ(NNR,NNC),Y(NNR,NNC),YY(NNR,NNC),IPU(NNR,NNC),EF(NNR,NNC),
C     .PRC(NNR,NNC),PE(NNR,NNC),AREA(NNR,NNC),
C     .VOLIN(NNR,NNC),DEFIC(NNR,NNC)
C      COMMON /THIESS/ TH(900),AMCC(900)
      COMMON /THIESS/ TH(NSTRMX),XP(NOBS,NSTRMX),IPU(NNR,NNC),
     +COORDT(2,NSTRMX),COORDP(2,NSTRMX)
	LOGICAL MASK(NNR,NNC) 
C      COMMON /CONST/ QPDT36,AREAU
C      COMMON /FORMLE/ FL
CGGG   Aggiungo il COMMON delle piogge per il calcolo dinamico 
CGGG   dei poligoni di Thissen
C      COMMON /PIO/ xp(nobs,NSTRMX),XPMED(nobs),T2M(NOBS,NSTRMX)
C       CHARACTER*80 FL
C       REAL COORD(2,900)

CC
C  INPUT:
C  NSTR=NUMERO STRUMENTI COSIDERATI NEL BACINO
C  NMISUR=NUMERO DI INTERVALLI CONSECUTIVI DI MISURAZIONE DELLE PIOGGE
C  (O PORTATE)
CC
CGGG	se le coordinate dei pluviometri a terra sono gia' state lette,
CGGG   per il calcolo dinamico dei poligoni di Thissen si rimanda
CGGG   subito alla chiamata della subroutine THISSEN
CRRR   altrimenti si leggono solo le coordinate la prima volta

      IF (ISTART.LE.0) then
c  LETTURA FILE BOLAM (salta 13 righe,
c   legge nst1  colonne e nst2 righe
c   salta riga
c   legge coord1 (longitudine) 
c   legge coord2 (latitudine)
c   salta riga
	    do nnskip=1,13
              READ(IDCB,*)
	    enddo
          READ(IDCB,'(20x,i3,1x,i3)') nst1,nst2
c	write(*,*)'nst1,nst2=',nst1,nst2
	    NSTR=NST1*NST2
	    WRITE(*,*)'NST. MISURA=',NSTR
	    READ(IDCB,*)
	    IF(ISTART.EQ.0) THEN
           READ(IDCB,*) (COORDP(1,I),I=1,NSTR)
	WRITE(*,*) 'C1=',COORDP(1,1),COORDP(1,2),COORDP(1,3)
		 READ(IDCB,*) (COORDP(2,I),I=1,NSTR)
	WRITE(*,*) 'C2=',COORDP(2,1),COORDP(2,2),COORDP(2,3)
	    ELSE
           READ(IDCB,*) (COORDT(1,I),I=1,NSTR)
	WRITE(*,*) 'C1=',COORDT(1,1),COORDT(1,2),COORDT(1,3)
		 READ(IDCB,*) (COORDT(2,I),I=1,NSTR)
	WRITE(*,*) 'C2=',COORDT(2,1),COORDT(2,2),COORDT(2,3)
	    ENDIF
      READ(IDCB,*)
	  ISTART=1           
      ELSE  
CC
C       CALCOLO DEI PESI DI THISSEN
CC
      CALL THISSEN (NSTR,IMISUR,IDIM,JDIM,XMIN,YMAX,DX,MASK)
c      WRITE(*,*)
c      WRITE(*,*) ' STAZIONE  RIGA COLONNA  PESO THISSEN AMCC '
c        DO 10 I=1,NSTR
c          WRITE(*,20) I,COORD(1,I),COORD(2,I),TH(I),AMCC(I)
c20         FORMAT(I5,2F10.2,F7.4,F6.1)
c10      CONTINUE       
25     FORMAT(A26,2I4,A2,F9.2,A4)
      endif
       RETURN
       END

        SUBROUTINE THISSEN(NSTR,IMISUR,IDIM,JDIM,XMIN,YMAX,DX,MASK)
CGGG Passo anche IMISUR per poter controllare la validita' del dato
CGGG di pioggia
      PARAMETER (NOBS=5140,NSTRMX=900,NNR=353,NNC=262)
      COMMON /THIESS/ TH(NSTRMX),XP(NOBS,NSTRMX),IPU(NNR,NNC),
     +COORDT(2,NSTRMX),COORDP(2,NSTRMX)
	LOGICAL MASK(NNR,NNC)
CGGG	 aggiungo il COMMON delle piogge per il calcolo dinamico
CGGG     dei poligoni di Thissen
C      COMMON /PIO/ xp(nobs,NSTRMX),XPMED(nobs),T2M(NOBS,NSTRMX)
C      DIMENSION COORD(2,900)

C      calcolo di IPU solo per aree con pendenza=ZZ nota
C      IPU PUNTATORI DELLE PIOGGE
C      COORD OGNI STRUMENTO HA DUE COORDIINATE E CI SONO 900 STRUMENTI
C
C      INPUT:NSTR = NUMERO DI STRUMENTI A DISPOSIZIONE
C            IDIM,JDIM= DIMENSIONE DELLE MATRICI
C            COORD=COORD(1,K) E' LA COORDINATA "I" DEL K-ESIMO STRUMENTO
C                  COORD(2,K) E' LA COORDINATA "J" DEL K-ESIMO STRUMENTO
C
C      OUTPUT:TH =AREA PERCENTUALE DEL POLIGONO DI THISSEN ASSOCIATA ALLO
C                 STRUMENTO I-ESIMO (MATRICE DEI PESI)

       WRITE(*,*) ' SUBROUTINE THISSEN'
C!?
      DO 20 I=1,NSTR
         TH(I)=0
20    CONTINUE

CC
C     DMIN DISTANZA AL QUADRATO FISSATA TRA PUNTO E STRUMENTO
C     D2 DISTANZA AL QUADRATO TRA LA CELLA IN CONSIDERAZIONE
C     ( COORDINATE (I,J) ) E LO STRUMENTO K-ESIMO
C     CON IPU SI DEFINISCE IL NUMERO KKK DELLA STAZIONE PIU'
C     VICINA A CIASCUNA CELLA
CC
      DO 10 I=1,IDIM
	YCOORD=YMAX-DX/2.-(I-1)*DX


         DO 1 J=1,JDIM
      	XCOORD=XMIN+DX/2.+(J-1)*DX

             IPU(I,J)=-1
             IF(MASK(I,J)) THEN
                 DMIN=1.E20

                 DO 100 K=1,NSTR
CGGG   Aggiungo un IF sul valore disponibile di pioggia per eliminare 
CGGG  la stazione con dato mancante dal calcolo dei poligoni di Thissen
	              IF(XP(IMISUR,K).GE.0.) THEN 
					 D2=(XCOORD-COORDP(1,K))*(XCOORD-COORDP(1,K))+
     &(YCOORD-COORDP(2,K))*(YCOORD-COORDP(2,K))

                      IF (D2.LT.DMIN) THEN
                          DMIN=D2
                          KKK=K
                      ENDIF
	              ENDIF
100               CONTINUE
 	           IPU(I,J)=KKK
                 TH(KKK)=TH(KKK)+1
             ENDIF
CC
C     TH= AREA DEL POLIGONO DI THISSEN ASSOCIATA ALLO STRUMENTO I-ESIMO
CC
1     CONTINUE
10     CONTINUE
C
	THH=0
      DO 30 I=1,NSTR
         THH=TH(I)+THH
30    CONTINUE
      DO 40 I=1,NSTR
        TH(I)=TH(I)/THH
40    CONTINUE                     
C      CLOSE(UNIT=9)
	  WRITE(*,*) ' CHIUSURA SUBROUTINE THISSEN'
      RETURN
      END
cccccccccccc
cccccccccccccc Wribin writes the matrix M with IDIM Rows and JDIM Columns on the binary file number ILU

        SUBROUTINE WRIbin(M,IDIM,JDIM,ILU)
	  PARAMETER	(NNR=353,NNC=262)
        real M(NNR,NNC)
c        integer M(NNR,NNC)
        do 1 i=1,idim
1          write(ILU) (M(i,j),j=1,jdim)
        return
        end
C Karakoram 2.5
        SUBROUTINE WRIbin1(NV,G,IDIM,JDIM,ILU)
	  PARAMETER	(NNR=353,NNC=262)
	  real G(NNR,NNC,10) 
	    do 1 i=1,idim
1		 write(ILU) (G(i,j,NV),j=1,jdim)
        return
        end

cccccccccccc
cccccccccccccc Wribin2 writes a series of iv=1,NV  matrixes G(i,j,iv=1,..NV) with IDIM Rows and JDIM Columns on the binary file number ILU
        SUBROUTINE WRIbin2(NV,G,IDIM,JDIM,ILU)
	  PARAMETER	(NNR=353,NNC=262)
	  real G(NNR,NNC,10) 
c        integer M(NNR,NNC)
        DO IV=1,NV
	    do  i=1,idim
c           do j=1,jdim
		 write(ILU) (NINT(G(i,j,iv)),j=1,jdim)
c		 write(ILU) G(i,j,iv)
c      	    enddo
 	    ENDDO
	  ENDDO
C 
        return
        end
cccccccccccc
cccccccccccccc Wribin3 writes a series of iv=1,NV  matrixes G(i,j,iv=1,..NV) with IDIM Rows and JDIM Columns on the binary file number ILU
cccccccccccccc if G(i,j,NV) 
        SUBROUTINE WRIbin3(NV,G,luse,IDIM,JDIM,ILU)
	  PARAMETER	(NNR=353,NNC=262)
	  real G(NNR,NNC,10),wrirow(NNC) 
	  integer (kind=2) luse(NNR,NNC) 
	
c        integer M(NNR,NNC)
        DO IV=1,NV
	    do  i=1,idim
           do j=1,jdim
	
	         if(G(i,j,nv).gt.0) then
				wrirow(j)=1
			 elseif(luse(i,j).eq.1) then
	            wrirow(j)=2
	         else
	            wrirow(j)=0
			 endif
      	  enddo
			 write(ILU) (wrirow(j),j=1,jdim)
 	    ENDDO
	  ENDDO
C 
        return
        end

ccccccccccccccccccc Plottting 3D Array 
        SUBROUTINE WRIV5D(NV,G,IDIM,JDIM,ILU)
	  PARAMETER	(NNR=353,NNC=262)
	  real G(NNR,NNC,10) 
	data nt/1/
c	if (nt.eq.1) then
c       OPEN(ilu,FILE='PERC1.OUT',STATUS='unknown',FORM='FORMATTED')
c	endif  
c	write(*,*)'nt=',nt,idim,jdim
      nt=nt+1

c	  write(*,*)'ilu=',ilu
crr giovanna attenta ho cambiato	do k=1,nv
	  do  i=1,idim
	     do j=1,jdim
c		     write(*,*) (G(i,j,k),k=1,NV)
crr		     write(ILU) G(i,j,k)
	         write(ILU) G(i,j,nv)
	     enddo
	   enddo
crr	enddo
	   
	return
	end 
ccccccccccccccccccc Plottting 3D Array 
        SUBROUTINE WRIV5DNAME(NV,G,IDIM,JDIM,ILU,CHARILU)
	  PARAMETER	(NNR=353,NNC=262)
	  real G(NNR,NNC,10) 
	CHARACTER*2 CHARILU
	data nt/1/
	if (nt.eq.1) then
       OPEN(ilu,FILE=CHARILU,STATUS='unknown',FORM='BINARY')
	endif  
	write(*,*)'nt=',nt,idim,jdim
      nt=nt+1

c	  write(*,*)'ilu=',ilu
crr giovanna attenta ho cambiato	do k=1,nv
	  do  i=1,idim
	     do j=1,jdim
c		     write(*,*) (G(i,j,k),k=1,NV)
crr		     write(ILU) G(i,j,k)
	         write(ILU) G(i,j,nv)
	     enddo
	   enddo
crr	enddo
	   
	return
	end 
cccccccccccc
cccccccccccccc WriascI2 writes the matrix M of integer I2 values with IDIM Rows and JDIM Columns on the ascii file number ILU

        SUBROUTINE WRIASCI2(M,IDIM,JDIM,ILU)
	  PARAMETER	(NNR=353,NNC=262)
        integer M(NNR,NNC)
        do 1 i=1,idim
1          write(ILU,'(262I2)') (M(i,j),j=1,jdim)
        return
        end
cccccccccccc
cccccccccccccc wriascF6 writes the matrix M of real F6.0 with IDIM Rows and JDIM Columns on the ascii file number ILU

        SUBROUTINE WRIASCF6(M,IDIM,JDIM,ILU)
	  PARAMETER	(NNR=353,NNC=262)
        real M(NNR,NNC)
        do 1 i=1,idim
1          write(ILU,'(262F6.0)') (M(i,j),j=1,jdim)
        return
        end
cccccccccccccc WriASCF7 writes the matrix M of real F8.0 with IDIM Rows and JDIM Columns on the ascii file number ILU

        SUBROUTINE WRIASCF7(NV,G,IDIM,JDIM,ILU)
	  PARAMETER	(NNR=353,NNC=262)
	  real G(NNR,NNC,10) 
	        do 1 i=1,idim
1		 write(ILU,'(262F7.0)') (G(i,j,nv),j=1,jdim) 
        return
        end

ccccccccccccccccccc Plotting the 3D array G(i,j,nv) of real F6.0 values with IDIM Rows and JDIM Columns on the ascii file number ILU
        SUBROUTINE WRIVASC(NV,G,IDIM,JDIM,ILU)
	  PARAMETER	(NNR=353,NNC=262)
	  real G(NNR,NNC,10) 
	data nt/1/
c	if (nt.eq.1) then
c       OPEN(ilu,FILE='PERC1.OUT',STATUS='unknown',FORM='FORMATTED')
c	endif  
c	write(*,*)'nt=',nt,idim,jdim
      nt=nt+1

c	  write(*,*)'ilu=',ilu
crr giovanna attenta ho cambiato	do k=1,nv
	  do  i=1,idim
	     do j=1,jdim
c		     write(*,*) (G(i,j,k),k=1,NV)
crr		     write(ILU) G(i,j,k)
	         write(ILU,'(262F6.0)') G(i,j,nv)
	     enddo
	   enddo
crr	enddo
	   
	return
	end 
ccccccccccccccccccc Plotting the 3D array G(i,j,nv) of real F6.0 values with IDIM Rows and JDIM Columns on the ascii file number ILU
        SUBROUTINE WRIVASC2(NV,G,LUSE,IDIM,JDIM,ILU)
	  PARAMETER	(NNR=353,NNC=262)
	  real G(NNR,NNC,10) 
	  integer (kind=2) LUSE(NNR,NNC)
	data nt/1/
c	if (nt.eq.1) then
c       OPEN(ilu,FILE='PERC1.OUT',STATUS='unknown',FORM='FORMATTED')
c	endif  
c	write(*,*)'nt=',nt,idim,jdim
      nt=nt+1

c	  write(*,*)'ilu=',ilu
crr giovanna attenta ho cambiato	do k=1,nv
	  do  i=1,idim
	     do j=1,jdim
c		     write(*,*) (G(i,j,k),k=1,NV)
crr		     write(ILU) G(i,j,k)
			 if(G(i,j,nv).gt.0) then
				wri=1
			 elseif(luse(i,j).eq.1) then
	            wri=2
	         else
	            wri=0
			 endif
	         write(ILU,'(262F6.0)') wri
	     enddo
	   enddo
crr	enddo
	   
	return
	end 


CRRMUL
CRRMUL ********************************************************
CRRMUL SNOW.FOR
CRRMUL
$debug
C
C       SUBROUTINE SNPACK Computes water storage and heat storage 
C       for each layer of snowpack
C
C       State variables
C
C       SwST(i)  =Water saturation as percentage of porosity [% 0-1]
C       PDKEF(i) =Heat content i-layer (negative)          [J] 
C       HSNPCK   =Height of snowpack                       [m]
C       PACT(i)  =Snowpack i-layer mean temperature        [C]
C       PSNOWP   = Input; Snow net precipitation(water equivalent)  [m/SECDT]
C       DENSN(i) =Snow i-layer density                     [kg/m³]
C       SECDT     Sampling step                            [s]
C       SECDTH    Integration step                         [s]
C       STPDTH    N. of Integration steps                  []            
C       HEACOG    Snow-groud Heat flow                     [W/m2]
C
C       Input
C
C       WATINP   =Incoming Water flow                      [m/SECDT]
C       WATOU    =Outgoing Water flow                      [m/SECDTH]
C       HEAFLO(i)=Incoming Heat flow for i-th layer        [W/m2]
C
C       Other variables
C       FHI   =Porosity [% 0-1]
C       FHIe  =Effective Porosity [% 0-1]
C       SwI   =Irredusible Water saturation (% porosity)  (0.07)[% 0-1]
C       SeST  =Effective Saturation         (% porosity)        [% 0-1] 
C       ALFASN=Permeability parameter (Colbeck,WRR82) (5470000)  [m/s] 
C       FHIeKs=Permeability parameter (Colbeck, "   ) (2.788E-3) [m**2/3]
c       CS    =Specific heat of ice   (Eagleson) (2093.4)        [J/(Kg C)]
C       CJOMME=Latent heat of fusion  (80000*4.186)              [J/Kg]
C
C       EPSH    =Minimum snowpack heigt, below which snowmelt is assumed [m]
C
	SUBROUTINE SNPACK(NX,NY,SwST,PDKEF,HSNPCK,PACT,DENSN,WATINP,
     &WATOU,HEAFLO,PSNOWP,Ta,SECDT,SECDTH,STPDTH,HEACOG)
	PARAMETER (NNR=353,NNC=262)
	INTEGER STPDTH
	REAL SWST(2,NNR,NNc),PDKEF(2,NNR,NNc),PACT(2,NNR,NNc)        
	REAL HSNPCK(NNR,NNc),DENSN(2,NNR,NNc),HEAFLO(2)
Cggc        REAL P1OLD(NNR,NNc),P2OLD(NNR,NNc) 
	LOGICAL NOPACK,LAY2
      COMMON /SNOW/A0,A1,TSPRIN,ASPRIN,TRNCF,AM0,CAM,GANESN,ALFAS,
     +BETAS,STRN,CJOMME,HSNLA1,FHI,Swi,FHIeKs,BEERSN,CMELTM,SPCFCT
      COMMON /HEAT/CPw,Wden,ROhu,Cpa,RatiKh,CKa,Z0,VLh,RatiKe,Rgas,ZW,ZT
	DATA ALFASN,PI,CS,EPSH,EPST,UM1,UM2/5470000.,3.1416,2093.4,
     +0.0005,.01,0.,0./
C
C  Variables assignement
C
	WATOU=0.
	WATIN=WATINP/STPDTH
	PSNOW=PSNOWP/STPDTH
	Tafsn=Ta
	IF(HSNPCK(NX,NY).LE.EPSH) THEN
C  Total snowpack melt
	  WATOU=AMAX1(HSNPCK(NX,NY)*DENSN(1,NX,NY)/Wden+WATIN,0.)
C          P1OLD(NX,NY)=0.
C          P2OLD(NX,NY)=0.
	  NOPACK=.TRUE.
	  LAY2=.FALSE.
	  HSNST2=0.
	  HSNS21=0.
	  PDKEF2=0.
	  PACT2=0.
	  DENSN2=0.
	  SwST2=0.
	    HSNST1=0.
	    HSNS11=0.
	    PDKEF1=0.
	    PACT1=0.
	    DENSN1=0.
	    SwST1=0.
	   GOTO 1000
	ELSE
	  NOPACK=.FALSE.
	  HSNPC=HSNPCK(NX,NY)                        
	  HSNST1=AMIN1(HSNLA1,HSNPC)                
	  SwST1=SwST(1,NX,NY)                       
	  PACT1=PACT(1,NX,NY)                       
	  DENSN1=DENSN(1,NX,NY)
	  PDKEF1=PACT1*CS*DENSN1*HSNST1                     
C          if(pdkef1.lt.0.) write(*,*)'att. pdkef1 pact1',pdkef1,pact1
	  HSNST2=HSNPC-HSNST1
	  IF(HSNST2.Gt.0)THEN                   
	    LAY2=.TRUE.                             
	    SwST2=SwST(2,NX,NY)                                           
	    PACT2=PACT(2,NX,NY)                     
	    DENSN2=DENSN(2,NX,NY)                   
	    PDKEF2=PACT2*CS*DENSN2*HSNST2
	  ELSE                                      
	    LAY2=.FALSE.                            
C 2nd layer total melt
	    WATIN=WATIN+HSNST2*DENSN(2,NX,NY)/Wden
	    HSNST2=0.                               
	    HSNS21=0.                               
	    PDKEF2=0.                               
	    PACT2=0.                                
	    DENSN2=0.                               
	    SwST2=0.                                
	  ENDIF
	  DENPCK=(DENSN1*HSNST1+DENSN2*HSNST2)/HSNPC                                     
	ENDIF
C
C  Heat conducted [ PRMS Model=2*SQRT(CS*DENPCK*ThKeff/PI/SECDTH)*(PACT1-PACT2)]
C  to lower layer or ground
	ThKeff=3.2232E-06*DENPCK*DENPCK
	IF(LAY2) THEN
	   HEACON=2.*SQRT(CS*DENPCK*ThKeff/PI/SECDTH)*(PACT1-PACT2)
C   Snow ground interface heat exchange
	   HEACOG=-2.*SQRT(CS*DENPCK*ThKeff/PI/SECDTH)*PACT2
	ELSE
	   HEACON=2.*SQRT(CS*DENPCK*ThKeff/PI/SECDTH)*PACT1
	ENDIF
C  If there isn't enough free water for evaporation it doesn't occurr
C  until additional melt is supplied
C  EVMELT melt=(EVR-rain)-swst*fhi*hsn  [m]
	 IF(WATIN.LT.0) THEN
	  EVMELT=AMAX1(-WATIN-SWST1*FHI*HSNST1,0.)
	  HEAJOU=(HEAFLO(1)-HEACON)*SECDTH-EVMELT*CJOMME*Wden
	 ELSE
	  HEAJOU=(HEAFLO(1)-HEACON)*SECDTH
	 ENDIF
C
C  Melting accounting  1st Layer
C  SNMELT=Internal Snowmelt water equivalent produced during DT [m] 
C
	IF(PACT1.LT.-EPST) THEN
C               T < 0 and PDKEF < 0.
C          Cold content has to be satisfied before melt occurs
C              HEASUR=Surplus heat useful for melt [Joules]
C          Regelation of condensed water/rainfall WATIN is computed, but 
C          additional advected heat is neglected
	  IF(PSNOW.GT.0) THEN
	    HEAJOU=HEAJOU+PSNOW*(Ta-PACT1)*WDEN*CS
	    Tafsn=PACT1
	  ENDIF
	  SNREGE=-AMAX1(WATIN,0.)
	  HEAREG=SNREGE*CJOMME*WDEN
	  PDKEF1=HEAJOU+PDKEF1-HEAREG
C          IF (PDKEF1.LT.0) WRITE(*,*)'COLD REG WAT',PDKEF1,HEAREG,WATIN
	  HEASUR=AMAX1(PDKEF1,0.) 
	  SNMELT=AMIN1(HEASUR/CJOMME,HSNST1*DENSN1)/WDEN
	  PDKEF1=AMIN1(PDKEF1,0.) 
	ELSE
C  Cold content =0. then melt or regelation can occur
C  SNMELT maximum possible melt of regelation, given HEAJOU
	  IF(PSNOW.GT.0) THEN
	    HEAJOU=HEAJOU+PSNOW*Ta*WDEN*CS
	    Tafsn=PACT1
	  ENDIF
	  SNMELT=AMIN1(HEAJOU/CJOMME,HSNST1*DENSN1)/Wden
	  SNREGE=0.
	  IF(HEAJOU.LT.0) THEN         
C   Regelation SNREGElation (negative)   [m]
C   HEAREG=Heat necessary for regelation of total water content [J]
C   (Here SNMELT (negative) means maximum possible regelation)
	    SNREGE=AMIN1(AMAX1(-SWST1*FHI*HSNST1-WATIN,SNMELT),0.)
	    SNMELT=0. 
	    HEAREG=SNREGE*CJOMME*Wden
	    PDKEF1=PDKEF1+(HEAJOU-HEAREG)
	  ENDIF 
	ENDIF 
C
	WST1=SwST1-Swi
	SeST1=WST1/(1.-Swi)
C        
	IF(SeST1.GT.0) THEN
C                       Water percolation
C       UMELT = Rate of total Melt production [m]
C       CMELT = Celerity of Melt -Wave  [m/s] (Colbeck,WRR82)
C       CMELTM=Multiplication factor
C       SIGMEL=Storage constant [s]
C       WSeST =Water Storage [m]
C
C The store releases water following u=ALFASN*K*Se^3   [Colbeck (4)+]
CGG	    UFLOW=0.5*ALFASN*(FHIeKs*FHI*wst1)**3
Crrr
            UFLOW=1.*ALFASN*(FHIeKs*FHI*wst1)**3
            CMELT=CMELTM*3*ALFASN**(1./3.)*FHIeKs*UFLOW**(2./3.)
            SIGMEL=HSNST1/CMELT
            DECAY=EXP(-SECDTH/SIGMEL)
Crrr
	    UMELT=WATIN+SNMELT+SNREGE
	    WSeST=AMAX1(SeST1*FHI*HSNST1,0.)
Cgg	    UINT=UFLOW*SECDTH
Cgg	    UTOT=UMELT+UM1
Cgg	    IF(UTOT.GT.0) THEN
Cgg	      IF (UINT.GE.UTOT) THEN
Cgg		  WATOU=UTOT
Cgg		  UM1=0.
Cgg	       ELSE
Cgg		  WATOU=UINT
Cgg		  UM1=UTOT-UINT
Cgg	      ENDIF
Cgg	     ELSE
Cgg	      WATOU=0.
Cgg	    ENDIF
Cgg            WATOU=(PERC1+P1OLD(NX,NY))/2.
Cgg            P1OLD(NX,NY)=PERC1
Crrr
            IF(UMELT.GT.0) THEN
              WATOU=WSeST/SIGMEL*DECAY*SECDTH+(1-DECAY)*UMELT
            ELSE                
              WATOU=WSeST/SIGMEL*DECAY*SECDTH
            ENDIF 
Cgg	 ELSE
Cgg	   WATOU=0.
CggC           P1OLD(NX,NY)=0.
	ENDIF
C        WRITE (27,'(3(2X,F8.4))') UMELT*1000.,WATOU*1000.,UINT*1000.
C
C       Final Water storage accounting and density computations
C       If new fallen snow an additional 0-Layer is introduced
C
c      WRITE(9,*) '   HEAJOU  HEACON  SNMELT   SNREGE    PACT1',
c     &'     PDKEF1  WATIN    WATOU'
c      WRITE(9,101) HEAJOU,HEACON,SNMELT,SNREGE,PACT1,PDKEF1,WATIN,WATOU
c101   FORMAT(f11.1,f7.1,6E10.4)
	CALL ACCUMU(HSNPC,HSNST1,DENSN1,HSNS11,PSNOW,SNMELT,WATIN,WATOU)
	IF(HSNS11.GT.0) THEN
	 SwST1=AMAX1(SwST1*HSNST1/HSNS11+(SNMELT+SNREGE+WATIN-WATOU)
     &/(FHI*HSNS11),0.)
	  WEI12=HSNS11*DENSN1
	  IF(WEI12.GT.0) PACT1=PDKEF1/(WEI12*CS)
C          if(pact1.lt.-10.) write(*,*)'pact1 pdkef1',pact1,pdkef1
	ELSE
	 PACT1=0.
	 SwST1=0.
	ENDIF
C                       2nd Layer Routine
	IF(LAY2)THEN
C                       Two layers interface
	HadvMe=WATOU*Cpw*Wden*(0-PACT2)                                       
	HEAJOU=(HEACON+HEACOG+HEAFLO(2))*SECDTH+HadvMe
	WATIN=WATOU
C                          Melting accounting
C      SNMELT=Internal Snowmelt water equivalent produced during DT [m] 
	IF(PACT2.LT.-EPST) THEN
C       T < 0 and PDKEF < 0.
C       Cold content has to be satisfied before melt occurs
	  SNREGE=-AMAX1(WATIN,0.)
	  HEAREG=SNREGE*CJOMME*WDEN
	  PDKEF2=HEAJOU+PDKEF2-HEAREG
	  HEASUR=AMAX1(PDKEF2,0.) 
	  SNMELT=HEASUR/CJOMME/Wden
	  PDKEF2=AMIN1(PDKEF2,0.) 
	ELSE
C
C       Cold content =0. then melt or regelation can occur
	  SNMELT=AMIN1(HEAJOU/CJOMME,HSNST2*DENSN2)/Wden
	  SNREGE=0.
	  IF(HEAJOU.LT.0) THEN        
	    SNREGE=AMIN1(AMAX1(-SWST2*FHI*HSNST2-WATIN,SNMELT) ,0.)
	    SNMELT=0. 
	    HEAREG=SNREGE*CJOMME*Wden
	    PDKEF2=PDKEF2+(HEAJOU-HEAREG)
	  ENDIF 
	ENDIF 
C
	WSt2=SwST2-Swi 
	SeST2=WST2/(1.-Swi)
C        
	IF(SeST2.GT.0) THEN
C                               Water percolation
Cgg	   UFLOW=0.5*ALFASN*(FHIeKs*FHI*wst2)**3
Crrr
	       UFLOW=1.*ALFASN*(FHIeKs*FHI*wst2)**3
             CMELT=CMELTM*3*ALFASN**(1./3.)*FHIeKs*UFLOW**(2./3.)
              SIGMEL=HSNST2/CMELT
              DECAY=EXP(-SECDTH/SIGMEL)
Crrr
	   UMELT=WATIN+SNMELT+SNREGE
	   WSeST=AMAX1(SeST2*FHI*HSNST2,0.)
Cgg	   UINT=UFLOW*SECDTH
Cgg	   UTOT=UMELT+UM2
Cgg	   IF(UTOT.GT.0) THEN
Cgg	     IF (UINT.GT.UTOT) THEN
Cgg	       WATOU=UTOT
Cgg	       UM2=0.
Cgg	      ELSE
Cgg	       WATOU=UINT
Cgg	       UM2=UTOT-UINT
Cgg	     ENDIF
Cgg	    ELSE
Cgg	     WATOU=0.
Cgg	   ENDIF
CggC           WATOU=(PERC2+P2OLD(NX,NY))/2.
CggC           P2OLD(NX,NY)=PERC2
Crrr
           IF(UMELT.GT.0) THEN
              WATOU=WSeST/SIGMEL*DECAY*SECDTH+(1-DECAY)*UMELT
            ELSE                
              WATOU=WSeST/SIGMEL*DECAY*SECDTH
           ENDIF
Crrr
	ELSE
	  WATOU=0.
CggC          P2OLD(NX,NY)=0.
	ENDIF 
C        WRITE (28,'(3(2X,F8.4))') UMELT*1000.,WATOU*1000.,SEST2
C
C      2nd Layer Final Water storage accounting and density computations
C
c       WRITE(9,*) '   HEAJOU  HEACON  SNMELT   SNREGE    PACT2',
c     &'     PDKEF2  WATIN    WATOU'
c      WRITE(9,101) HEAJOU,HEACOG,SNMELT,SNREGE,PACT2,PDKEF2,WATIN,WATOU
	CALL ACCUMU(HSNPC,HSNST2,DENSN2,HSNS21,0.,SNMELT,WATIN,WATOU)
	IF(HSNS21.GT.0.) THEN
	 SwST2=AMAX1(SwST2*HSNST2/HSNS21+(SNMELT+SNREGE+WATIN-WATOU)
     &/(FHI*HSNS21),0.)
	  WEI22=HSNS21*DENSN2
	  IF(WEI22.GT.0) PACT2=PDKEF2/(WEI22*CS)
	ELSE
	 PACT2=0.
	 SwST2=0.
	ENDIF
C                 End of 2nd Layer computations
	ENDIF
C              2 Layers aggregations and re-computations
C       H01 =New fallen snow part of new layer 1
C       H11 =Old layer 1     part of new layer 1
C       H21 =Old layer 2     part of new layer 1
C           H02 =New fallen snow part of new layer 2
C           H12 =Old layer 1     part of new layer 2
C
1000    HSNST0=PSNOW*Wden/GANESN
C          WRITE(9,*) ' HSNST0    HSNST1   HSNS11   HSNST2   HSNS21'
C          WRITE(9,110) HSNST0,HSNST1,HSNS11,HSNST2,HSNS21
C110       FORMAT(5E10.4)
	  HSNPCK(NX,NY)=HSNST0+HSNS11+HSNS21
	  H01=AMIN1(HSNLA1,HSNST0)
	  H11=AMIN1(HSNLA1-H01,HSNS11)
	  H21=AMIN1(HSNLA1-H11-H01,HSNS21)
	  HSNST1=H01+H11+H21
	  WEIHS1=GANESN*H01+DENSN1*H11+DENSN2*H21                          
	
	IF(WEIHS1.GT.0) THEN
	   IF(PDKEF1.LT.-EPST) THEN
C   PDKEF2 will remain 0. and heat deficit is stored in PDKEF1
	     IF(PDKEF2.GT.-EPST) THEN
	       PDKEF1=(HSNST0*Tafsn*GANESN+HSNS11*PACT1*DENSN1+
     &HSNS21*PACT2*DENSN2)
	     ELSE
	       PDKEF1=(H01*Tafsn*GANESN+H11*PACT1*DENSN1+
     &H21*PACT2*DENSN2)
	     ENDIF
	   ELSE
	     PDKEF1=0.
	   ENDIF 
C             NEW Snowpack layers density is computed by mass balance            
C             Fresh snow is supposed not to have free water content              
	   DENSN(1,NX,NY)=WEIHS1/HSNST1 
	                                               
	   PACT(1,NX,NY)=PDKEF1/WEIHS1         
C           WRITE(*,*) 'PDKEF1,WEIHS1=',PDKEF1,WEIHS1
	   IF(.NOT.NOPACK) SwST(1,NX,NY)=SwST1*HSNS11/HSNST1
	 ELSE
C                          Snowpack disappeared
	   HSNPCK(NX,NY)=0.
	   PACT(1,NX,NY)=0.
	   PDKEF(1,NX,NY)=0.
	   DENSN(1,NX,NY)=0.
	   SwST(1,NX,NY)=0.
	   PACT(2,NX,NY)=0.
	   PDKEF(2,NX,NY)=0.
	   DENSN(2,NX,NY)=0.  
	   SwST(2,NX,NY)=0.
C
	   RETURN 
	 ENDIF 
C  New  2nd Layer if old layers exist or new fallen snow contributes
	    H02=AMAX1(HSNST0-H01,0.)
	    H12=HSNS11-H11
	    H22=HSNS21-H21
	  HSNST2=H02+H12+H22
	    WEIHS2=GANESN*H02+DENSN1*H12+DENSN2*H22
	IF(WEIHS2.GT.0) THEN
	   IF(PDKEF2.LT.-EPST) THEN
C   PDKEF2 will remain 0. and heat deficit is stored in PDKEF2
	      IF(PDKEF1.GT.-EPST) THEN
		PDKEF2=(HSNST0*Tafsn*GANESN+
     &HSNS11*PACT1*DENSN1+HSNS21*PACT2*DENSN2)
	      ELSE
	       PDKEF2=H02*Tafsn*GANESN+H12*PACT1*DENSN1+H22*PACT2*DENSN2
	      ENDIF
	   ELSE
	      PDKEF2=0. 
	   ENDIF
	  DENSN(2,NX,NY)=WEIHS2/HSNST2
	  PACT(2,NX,NY)=PDKEF2/WEIHS2
	  IF(.NOT.NOPACK) SwST(2,NX,NY)=SwST2*HSNS21/HSNST2
	ELSE
	  PACT(2,NX,NY)=0.
	  PDKEF(2,NX,NY)=0.  
	  DENSN(2,NX,NY)=0.
	  SwST(2,NX,NY)=0.
	  WEIHS2=0.
	ENDIF
	RETURN
	END
C
C       ACCUMU Snow accumulation Model (See NAM-II Model)
C
C       HSNPC = Input; Depth of snowpack at time t-1         [m]
C       DSTM1 = Input; Depth of snowpack layer at time t-1         [m]
C       DENSNO= Input/OUTPUT Snowpack specific weight      [Kg/mc]
C       DST   = OUTPUT;Depth of snowpack at time t           [m] 
C       PSNOW = Input; Snow net precipitation(water equivalent)  [m]
C       SNMELT  = Input; Internal Snow melt  (water equivalent)  [m]
C       WATIN   = Input; Water input         (water equivalent)  [m]
C       WATOU  = Input; Snow melt released  (water equivalent)  [m]
C
C       GANESN= New snow specific weight                   [Kg/mc]
C       Wden  = Water specific weight                      [Kg/mc]
C       ALFAS  = Parameter   (0.25 NAM-II Model)                   [m]
C       BETAS  = Parameter   (0.35)                        [ ]
C
      SUBROUTINE ACCUMU(HSNPC,DSTM1,DENSNO,DST,PSNOW,SNMELT,WATIN,WATOU)
      COMMON /SNOW/A0,A1,TSPRIN,ASPRIN,TRNCF,AM0,CAM,GANESN,ALFAS,
     +BETAS,STRN,CJOMME,HSNLA1,FHI,Swi,FHIeKs,BEERSN,CMELTM,SPCFCT
      COMMON /HEAT/CPw,Wden,ROhu,Cpa,RatiKh,CKa,Z0,VLh,RatiKe,Rgas,ZW,ZT
      IF((DSTM1.LE.0).OR.(DENSNO.LE.0)) THEN
	DST=0.
	DENSNO=0.
	RETURN
      ENDIF 
      RAPPGL=WDEN/DENSNO 
C  Reduction due to compression becouse of new fallen snow
      REDUCT=PSNOW*RAPPGL*(HSNPC/ALFAS)**BETAS
C  Reduction due to snow/water density change.Volume factor 2/PI=0.62662 or .25
      DMELT=SNMELT*RAPPGL
      DST=DSTM1-REDUCT-DMELT
      IF(DST.LE.0) THEN
	 WATOU=AMAX1(DENSNO*DSTM1/WDEN+WATIN,0.)
	 DST=0.
	 DENSNO=0.
C         WRITE(9,*)  '    DSTM1    DST       REDUCT     DMELT'
C         WRITE(9,120) DSTM1,DST,REDUCT,DMELT
C120      FORMAT(4E10.4)
	 RETURN
      ENDIF
C         WRITE(9,*)  '    DSTM1    DST       REDUCT     DMELT'
C         WRITE(9,120) DSTM1,DST,REDUCT,DMELT
	DENSNO=AMIN1(AMAX1(DENSNO*DSTM1/DST+(WATIN-WATOU)*WDEN/DST,0.),
     +1000.)
c        DENSNO=AMAX1(DENSNO*(2-DST/DSTM1)+(WATIN-WATOU)*
c     &WDEN/DSTM1,0.)

      RETURN
      END               
C
C     TSPK  Snowpack Mean temperature [.ø.C]
C
C     Input
C     HSNPCK   [m]
C     DENSN    [kg/m³]
C     PACT     [.ø.C]
      FUNCTION TSPK(NX,NY,HSNPCK,DENSN,PACT)
	PARAMETER (NNR=353,NNC=262)
	REAL HSNPCK(NNR,NNC),DENSN(2,NNR,NNC),PACT(2,NNR,NNC)        
      COMMON /SNOW/A0,A1,TSPRIN,ASPRIN,TRNCF,AM0,CAM,GANESN,ALFAS,
     +BETAS,STRN,CJOMME,HSNLA1,FHI,Swi,FHIeKs,BEERSN,CMELTM,SPCFCT
	HSNPC=HSNPCK(NX,NY)
	HS1=amin1(HSNPC,HSNLA1)
	WEIHS1=HS1*DENSN(1,NX,NY)
	IF(WEIHS1.LE.0) THEN
	  Tspk=0.
	 RETURN
	ENDIF
	WEIHS2=(HSNPC-HS1)*DENSN(2,NX,NY)
	Tspk=(WEIHS1*PACT(1,NX,NY)+WEIHS2*PACT(2,NX,NY))/
     &(WEIHS1+WEIHS2)
	RETURN
	END
C
C     WEQU    Computes snow water equivalent  [m]
C     
C     Input
C     HSNPCK   [m]
C     DENSN    [kg/m³]
      FUNCTION WEQU(NX,NY,HSNPCK,DENSN)
	PARAMETER (NNR=353,NNC=262)
	REAL HSNPCK(NNR,NNc),DENSN(2,NNR,NNc) 
      COMMON /SNOW/A0,A1,TSPRIN,ASPRIN,TRNCF,AM0,CAM,GANESN,ALFAS,
     +BETAS,STRN,CJOMME,HSNLA1,FHI,Swi,FHIeKs,BEERSN,CMELTM,SPCFCT
      COMMON /HEAT/CPw,Wden,ROhu,Cpa,RatiKh,CKa,Z0,VLh,RatiKe,Rgas,ZW,ZT
	HSNPC=HSNPCK(NX,NY)
	HS1=amin1(HSNPC,HSNLA1)
	WEQU=(HS1*DENSN(1,NX,NY)+(HSNPC-HS1)*DENSN(2,NX,NY))/WDEN
      RETURN
      END
C
C     LWC Computes liquid water content  [mm]
C     
C     Input
C     HSNPCK [m]
C     SWST   [-]
      REAL FUNCTION LWC(NX,NY,HSNPCK,SWST)
	PARAMETER (NNR=353,NNC=262)
	REAL HSNPCK(NNR,NNc),SWST(2,NNR,NNc) 
      COMMON /SNOW/A0,A1,TSPRIN,ASPRIN,TRNCF,AM0,CAM,GANESN,ALFAS,
     +BETAS,STRN,CJOMME,HSNLA1,FHI,Swi,FHIeKs,BEERSN,CMELTM,SPCFCT
      COMMON /HEAT/CPw,Wden,ROhu,Cpa,RatiKh,CKa,Z0,VLh,RatiKe,Rgas,ZW,ZT
	HSNPC=HSNPCK(NX,NY)
	HS1=amin1(HSNPC,HSNLA1)
	LWC=(HS1*SWST(1,NX,NY)*FHI+(HSNPC-HS1)*SWST(2,NX,NY)*FHI)*1000.
      RETURN
      END
CRRMUL
CRRMUL ********************************************************
CRRMUL SUB.FOR
CRRMUL

$debug
C  Spatial statistics
C	20040608  AWRITE RIDIMENSIONATA PER SOLE 2 VARIABILI, INVECE CHE 96
c20060201 melt spatial average computed here
       SUBROUTINE 
     +STATXY(WEQUI,XYMELT,MELT,XYICEM,PE,MELTICE,PEMED,PXY,PAREA,EVRXY,
     +SNOWXY, EVRMEDS, EVRMEDI, SNOWMEDS, RAINMEDS,RAINMEDI,MELTSNOW,
     +NXY,NSNOWXY,LUSE,NSNOWICE)
CCC
C       Input
C       WEQUI   Distributed water equivalent
C       XYMELT  Distributed snowmelt matrix
C       XYICEM  Distributed icemelt
C       PE	  Distributed net rainfall -(Rain-Evr)
C       PXY     Distributed net precipitation Rain-Evrpot+Snow 
C

C
C       Output
C
C       MELT    Spatial average Snowmelt
C	  MELTICE Spatial average Icemelt
C       PEMED   Spatial average net rainfall + snow/icemelt
c       PAREA   Spatial average of net precipitation
c       EVRXY   Distributed Evaporation
c       SNOWXY  Snowfall Precipitation
c       EVRMEDS Spatial Mean Evaporation over Snow 
c       EVRMEDI Spatial Mean Evaporation over snowfree area  
c       SNOWMEDS Spatial Average of snowfall over snow
c       RAINMEDS Spatial Average of rainfall over snow
c       RAINMEDI Spatial Average of rainfall over snowfree area
c       MELTSNOW Spatial Average of Snowmelt
c       NXY     # grids cells in the basin  
c       NSNOWXY # snowcover cells 
c 
c       Area Common   
c     
C       SWEM    Spatial mean Water EQUIVALENT over total basin area 
C       SWAVXY  Distributed Net shortwave radiation matrix
C       RAMUXY  Spatial mean Net shortwave radiation over area covered by snow
C       LONGM   Spatial mean Net longwave radiation over area covered by snow
C       CONVM   Spatial mean CONVECTIVE FLUX over area covered by snow
C       LATM    Spatial mean LATENT FLUX over area covered by snow
C       HEMUXY  Spatial mean HEAT BALANCE over area covered by snow
C       SWEST   Spatial standard deviation Water EQUIVALENT 
C       RASTXY  Spatial standard deviation Net shortwave radiation
C       LONGST  Spatial standard deviation Net longwave radiation
C       CONVST  Spatial standard deviation CONVECTIVE FLUX
C       LATST   Spatial standard deviation LATENT FLUX
C       HESTXY  Spatial standard deviation HEAT BALANCE 
C       MESTXY  Spatial standard deviation Melt over BASIN area
C       NXY     # grids in the basin 
C
CCC
      PARAMETER (NNR=353,NNC=262,NBAS=10)
      INTEGER (kind=2) ind(nnr,nnc),luse(nnr,nnc)
      INTEGER NB(10,7),IDSUB(nbas), IMASK(nbas)
      REAL MAX(7), MIN(7), VAL(7),DEL(7),PB(10,7),DEVIDE(10,7)
      REAL MELT(NBAS),MESTXY(NBAS),LONGM(NBAS),LONGST(NBAS),LATM(NBAS)
	REAL LATST(NBAS),MELTICE(NBAS),pemed(nbas),PAREA(NBAS)
	REAL HEMUXY(NBAS)
	REAL EVRMEDS(NBAS),EVRMEDI(NBAS),SNOWMEDS(nbas),NSNOWXY(NBAS)
	REAL RAINMEDS(NBAS),RAINMEDI(NBAS), MELTSNOW(NBAS),NSNOWICE(NBAS)
      REAL XYMELT(NNR,NNC),WEQUI(NNR,NNC),HEABXY(NNR,NNC)
	real xyicem(nnr,nnc),pe(nnr,nnc),PXY(NNR,NNC),EVRXY(NNR,NNC),
     +SNOWXY(NNR,NNC)
      REAL DEM(NNR,NNC),XNORM(3,NNR,NNC),SWST(2,NNR,NNC),
     +PDKEF(2,NNR,NNC)
      REAL PACT(2,NNR,NNC),HSNPCK(NNR,NNC),DENSN(2,NNR,NNC)
      LOGICAL MASK(NNR,NNC)       
      COMMON /DEM/ DEM,XNORM,HSENSO,NXSENS,NYSENS,NR,NC,DEMDX,DEMDY,
     +AREAA,NGRIDS
      COMMON /SNPK/ SWST,PDKEF,PACT,HSNPCK,DENSN
      COMMON /MASKD/ MASK,NXPRNT,NYPRNT,IND,NUMSUB,IDSUB,IMASK
      COMMON/RADIA/SWAVXY(nnr,nnc),RLOWXY(nnr,nnc),HCONVXY(nnr,nnc),
     +HLATXY(nnr,nnc),HPXY(nnr,nnc),SNOGRXY(nnr,nnc)
      COMMON /AGGRES/ RAMUXY(NBAS),RASTXY(NBAS),HPMUXY(NBAS),
     +HPSTXY(NBAS),LONGM,LONGST,CONVM(NBAS),CONVST(NBAS),LATM,
     +LATST,MESTXY,SWEM(NBAS),SWEST(NBAS),SVM(7)
      COMMON /MOMENT/SCSW(3),SCRL(3),SCHC(3),SCHL(3)
C
C  Be careful! NCLASS is also the dimension of the local variables
C      NB(NCLASS,IV),PB(NCLASS,IV),DEVIDE(NCLASS,IV)
C
      DATA EPSH,NCLASS,ITER/0.,10,0/  
C
      ITER=ITER+1
	DO 1 I=1,7
	   MAX(I)=-10000.
	   MIN(I)=10000.
	   SVM(I)=-99.
	   DEL(I)=0.
	   DO 1 J=1,NCLASS
	     NB(J,I)=0
	     PB(J,I)=0.
	     DEVIDE(J,I)=0. 
1       CONTINUE
	NXY=0
	DO ISS=1,NUMSUB
c20060201 melt spatial average
      MELT(ISS)=0.  
	SWEM(ISS)=0.
	MELTICE(ISS)=0.
	PEMED(ISS)=0.
	PAREA(ISS)=0.
	RAMUXY(ISS)=0.
	LONGM(ISS)=0.
	CONVM(ISS)=0.      
	LATM(ISS)=0.
	HPMUXY(ISS)=0.
	SWEST(ISS)=0.
	RASTXY(ISS)=0.
	HEMUXY(ISS)=0.
	LONGST(ISS)=0.
	CONVST(ISS)=0.
	LATST(ISS)=0.
	HPSTXY(ISS)=0.
	MESTXY(ISS)=0.
	EVRMEDS(ISS)=0.
	EVRMEDI(ISS)=0.
      SNOWMEDS(ISS)=0.
	RAINMEDS(ISS)=0.
	RAINMEDI(ISS)=0.
	NSNOWXY(ISS)=0
	NSNOWICE(ISS)=0
	ENDDO
	DO I=1,3
	   SCSW(I)=0.
	   SCRL(I)=0.
	   SCHL(I)=0.
	   SCHC(I)=0.
	END DO
	DO 100 I=1,NR
	  do 200 J=1,NC
	   IF(MASK(I,J)) THEN
C In the basin area
CGG	    IF(HSNPCK(I,J).GT.EPSH) NSNOWXY=NSNOWXY+1
C in the snowpack area (se commentata dovrebbe risultare nxy=ngrids)
	      NXY=NXY+1

	      DO ISS=1,NUMSUB
 	        NUM=IDSUB(ISS)
      		if(ind(i,j).eq.NUM) then
		  SWEM(ISS)=SWEM(ISS)+WEQUI(i,j)*1000.
	      MELTICE(ISS)=MELTICE(ISS)+XYICEM(I,J)
	      PEMED(ISS)=PEMED(ISS)+PE(I,J)
	      PAREA(ISS)=PAREA(ISS)+PXY(I,J)
c20060201 melt spatial average
            MELT(ISS)=MELT(ISS)+XYMELT(I,J)
	      RAMUXY(ISS)=RAMUXY(ISS)+SWAVXY(i,j)
	      LONGM(ISS)=LONGM(ISS)+RLOWXY(i,j)
	      CONVM(ISS)=CONVM(ISS)+HCONVXY(i,j)
	      LATM(ISS)=LATM(ISS)+HLATXY(i,j)
C HPXY is now the cold content of the first layer PDKEF(1,NX,NY)
Cgg	      write(*,*) hpmuxy,i,j,hpxy(i,j)
		  HPMUXY(ISS)=HPMUXY(ISS)+HPXY(i,j)
	    IF(HSNPCK(I,J).GT.EPSH) THEN
		   NSNOWXY(ISS)=NSNOWXY(ISS)+1
	       NSNOWICE(ISS)=NSNOWICE(ISS)+1
            EVRMEDS(ISS)=EVRMEDS(ISS)+EVRXY(I,J)
            SNOWMEDS(ISS)=SNOWMEDS(ISS)+SNOWXY(I,J)
	      RAINMEDS(ISS)=RAINMEDS(ISS)+PXY(I,J)+EVRXY(I,J)-SNOWXY(I,J)
            MELTSNOW(ISS)=MELTSNOW(ISS)+XYMELT(I,J)
	ELSE
	      IF(LUSE(I,J).EQ.1) NSNOWICE(ISS)=NSNOWICE(ISS)+1
            EVRMEDI(ISS)=EVRMEDI(ISS)+EVRXY(I,J)
	      RAINMEDI(ISS)=RAINMEDI(ISS)+PXY(I,J)+EVRXY(I,J)
	ENDIF


	      GOTO 777
	      ENDIF
	      ENDDO
777	      DO MM=1,3
		 EE=MM+1
		 SCSW(MM)=SCSW(MM)+SWAVXY(I,J)**(EE)
		 SCRL(MM)=SCRL(MM)+RLOWXY(I,J)**(EE)
		 SCHC(MM)=SCHC(MM)+HCONVXY(I,J)**(EE)
		 SCHL(MM)=SCHL(MM)+HLATXY(I,J)**(EE)
	      END DO
C
C  Minimum and maximum value
C
	      VAL(1)=WEQUI(i,j)*1000.
	      VAL(2)=SWAVXY(i,j)
	      VAL(3)=RLOWXY(i,j)
	      VAL(4)=HCONVXY(i,j)
	      VAL(5)=HLATXY(i,j)
	      VAL(6)=HPXY(i,j)
	      IF (XYMELT(I,J).LT.0.) XYMELT(I,J)=0.  
	      VAL(7)=XYMELT(i,j)
	      DO 9 IV=1,7  
	       IF(VAL(IV).GT.MAX(IV)) MAX(IV)=VAL(IV)
	       IF(VAL(IV).LT.MIN(IV)) MIN(IV)=VAL(IV)
9             CONTINUE
CC          ELSE
C             WRITE(*,*) ' NEVE SCOMPARSA NELLA CELLA',I,J
C             STOP
Cggg!!!	   ENDIF 
C Height of snowpack>0
	   ENDIF
C     Mask=.true.
C
200       continue
100     continue
C
C       WRITE(18,955) ITER,(MAX(IV),IV=1,7)
C       WRITE(19,955) ITER,(MIN(IV),IV=1,7)
C955     FORMAT(I4,1X,7F9.3)     
C
C       IF(NXY.GT.0) then
 	  DO ISS=1,NUMSUB
	  MELTICE(ISS)=MELTICE(ISS)/float(IMASK(ISS))         
	  PEMED(ISS)=PEMED(ISS)/float(IMASK(ISS))         
	  MELT(ISS)=MELT(ISS)/float(IMASK(ISS))
	  PAREA(ISS)=PAREA(ISS)/float(IMASK(ISS))
	  RAMUXY(ISS)=RAMUXY(ISS)/float(IMASK(ISS))         
	  LONGM(ISS)=LONGM(ISS)/float(IMASK(ISS))        
	  CONVM(ISS)=CONVM(ISS)/float(IMASK(ISS))         
	  LATM(ISS)=LATM(ISS)/float(IMASK(ISS)) 
   	  SWEM(ISS)=SWEM(ISS)/float(IMASK(ISS))
	  HPMUXY(ISS)=HPMUXY(ISS)/float(IMASK(ISS))         
C
        EVRMEDS(ISS)=EVRMEDS(ISS)/float(IMASK(ISS))
        SNOWMEDS(ISS)=SNOWMEDS(ISS)/float(IMASK(ISS))
        RAINMEDS(ISS)=RAINMEDS(ISS)/float(IMASK(ISS))
        MELTSNOW(ISS)=MELTSNOW(ISS)/float(IMASK(ISS))
C        IF(NSNOWXY(ISS).LT.IMASK(ISS)) THEN
	    EVRMEDI(ISS)=EVRMEDI(ISS)/float(IMASK(ISS))
	    RAINMEDI(ISS)=RAINMEDI(ISS)/float(IMASK(ISS))
C		ELSE
C		EVRMEDI(ISS)=0.
C		RAINMEDI(ISS)=0.
C	   ENDIF    
C


	  ENDDO         
	  DO MM=1,3
	     SCSW(MM)=SCSW(MM)/float(NXY)         
	     SCRL(MM)=SCRL(MM)/float(NXY)         
	     SCHC(MM)=SCHC(MM)/float(NXY)         
	     SCHL(MM)=SCHL(MM)/float(NXY)         
	  END DO
	  DO 2 IV=1,7 
	     DEL(IV)=(MAX(IV)-MIN(IV))/FLOAT(NCLASS)
	     IF (DEL(IV).LE.0.) GOTO 2 
	     DO 23 IB=1,(NCLASS-1)
23             DEVIDE(IB,IV)=MIN(IV)+FLOAT(IB)*DEL(IV)
2         CONTINUE 

      DO 300 I=1,NR
       DO 400 J=1,NC
	  IF(MASK(I,J)) THEN
C
C            IF(NXY.GT.NGRIDS) WRITE(22,*) I,J,HSNPCK(I,J),HPXY(I,J),
C     &XYMELT(I,J),SWAVXY(I,J)
C
CGG	     IF(HSNPCK(I,J).GT.EPSH) then

 	  DO ISS=1,NUMSUB
 	        NUM=IDSUB(ISS)
      		if(ind(i,j).eq.NUM) then
		RASTXY(ISS)=RASTXY(ISS)+(SWAVXY(i,j)-RAMUXY(iss))**2
		LONGST(ISS)=LONGST(ISS)+(RLOWXY(i,j)-LONGM(iss))**2
		CONVST(ISS)=CONVST(ISS)+(HCONVXY(i,j)-CONVM(iss))**2
		LATST(ISS)=LATST(ISS)+(HLATXY(i,j)-LATM(iss))**2
		HPSTXY(ISS)=HPSTXY(ISS)+(HPXY(i,j)-HPMUXY(iss))**2
          MESTXY(ISS)=MESTXY(ISS)+(XYMELT(i,j)-MELT(iss))**2
          SWEST(ISS)=SWEST(ISS)+(WEQUI(i,j)*1000.-SWEM(iss))**2
 		GOTO 779
         ENDIF
	    ENDDO
779		VAL(2)=SWAVXY(i,j)
		VAL(3)=RLOWXY(i,j)
		VAL(4)=HCONVXY(i,j)
		VAL(5)=HLATXY(i,j)
		VAL(6)=HPXY(i,j)
		VAL(7)=XYMELT(i,j)
		DO 4 IV=2,7  
		 DO 4 IB=1,(NCLASS-1)
4                   IF(VAL(IV).LE.DEVIDE(IB,IV)) NB(IB,IV)=NB(IB,IV)+1
CGG	     ENDIF
cc
c       over total basin area
cc 
C	      SWEST=SWEST+(WEQUI(i,j)*1000.-SWEM)**2
	      DO 12 IB=1,(NCLASS-1) 
12           IF((WEQUI(I,J)*1000.).LE.DEVIDE(IB,1))
     +NB(IB,1)=NB(IB,1)+1           
C
C  Mask=.true.
C
	   ENDIF
C
400        continue
300     continue
C
C        write (*,*) 'nb(ib,7)=',(nb(ib,7),ib=1,nclass)
C        write (*,*) 'devide(ib,7)=',(devide(ib,7),ib=1,nclass)
C          
C           DO 3 N=1,7
C              WRITE (*,*) 'NB=',(NB(M,N),M=1,NCLASS)
C              WRITE (*,*) 'min,max,del=',MIN(N),MAX(N),DEL(N)
C3          CONTINUE
C
	   DO 5 IV=1,7   
	      ENTROPY=0. 
	      IF(IV.EQ.1) THEN
		 NB(NCLASS,IV)=NGRIDS
c	       ELSE
		 IF (NXY.EQ.0) GOTO 5
		 NB(NCLASS,IV)=NXY
	      ENDIF  
	      IF(DEL(IV).LE.0.) THEN
		 SVM(IV)=1.
		 GOTO 5
	      ENDIF   
	       DO 6 IB=1,NCLASS
		  IF (IB.EQ.1) THEN
		      BB=FLOAT(NB(IB,IV))
		      IF (BB.LE.0.) GOTO 6
		   ELSE  
		      BB=FLOAT(NB(IB,IV)-NB(IB-1,IV))
		      IF (BB.LE.0.) GOTO 6 
		  ENDIF
		  PB(IB,IV)=BB/FLOAT(NB(NCLASS,IV)+1)
cgg momentaneamente commentato per evitare errori MAth
c		  ENTROPY=ENTROPY+PB(IB,IV)*ALOG(PB(IB,IV))
6              CONTINUE  
	   AA=FLOAT(NB(NCLASS,IV))
Cgg commentato termporaneamente perche' da' errorec
Cgg	   SVM(IV)=1.+ENTROPY/ALOG(AA)
5          CONTINUE
C
	 
C	 IF(SWEST.GT.0.) SWEST=SQRT(SWEST/float(NGRIDS))         
C	 IF(NXY.GT.1)THEN
 	  DO ISS=1,NUMSUB
	    RASTXY(ISS)=SQRT(RASTXY(ISS)/float(IMASK(ISS)))         
	    LONGST(ISS)=SQRT(LONGST(ISS)/float(IMASK(ISS)))          
	    CONVST(ISS)=SQRT(CONVST(ISS)/float(IMASK(ISS)))            
	    LATST(ISS)=SQRT(LATST(ISS)/float(IMASK(ISS)))             
	    HPSTXY(ISS)=SQRT(HPSTXY(ISS)/float(IMASK(ISS)))             
	    MESTXY(ISS)=SQRT(MESTXY(ISS)/float(IMASK(ISS)))             
	    SWEST(ISS)=SQRT(SWEST(ISS)/float(IMASK(ISS)))           
	ENDDO
C	 ENDIF 
C
C
 	RETURN
	END   
C
C
	SUBROUTINE OUTPUT(NR,NC,XYMELT,istep,IUOUT)
CC      
C    Scrive sul file TESIOUT.DAT le matrici dell'equivalente in acqua e 
C    dello scioglimento
CC
      PARAMETER (NNR=353,NNC=262)
C REAL SWST(2,20,20),PDKEF(2,20,20),PACT(2,20,20),DENSN(2,20,20)        
	REAL XYMELT(NNR,NNC),awrite(336420,2),asum(336420)
	data asum/336420*0./
	NW=0
	DO 10 I=1,NR
	   DO 11 J=1,NC
	nw=nw+1
	asum(nw)=asum(nw)+xymelt(i,j)*3.6/1000.
11         CONTINUE 
10      CONTINUE   
	if (istep.eq.25) then
	   NW=0
	   DO I=1,NR
	     DO J=1,NC
	        NW=NW+1
       	    WRITE(IUOUT,120) asum(nw)
	     enddo
	   enddo
	endif
100     FORMAT(F6.3)
120     FORMAT(48(1x,F7.2))
C
	RETURN
	END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	SUBROUTINE OUTP1(NR,NC,XYMELT,istep,IUOUT)
      PARAMETER (NNR=353,NNC=262)
	REAL XYMELT(NNR,NNC),awrite(336420,2),asum(336420)
	data asum/336420*0./
C        
	NW=0
	DO 10 I=1,NR
	   DO 11 J=1,NC
	nw=nw+1

	asum(nw)=asum(nw)+xymelt(i,j)*3.6/1000.
11         CONTINUE 
10      CONTINUE   
	if (istep.eq.25) then
	   NW=0
	   DO I=1,NR
	     DO J=1,NC
	        NW=NW+1

       	    WRITE(IUOUT,120) asum(nw)
	     enddo
	   enddo
	endif
100     FORMAT(F6.3)
120     FORMAT(48(1x,F7.2))
	RETURN
	END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	SUBROUTINE OUTP2(NR,NC,XYMELT,istep,IUOUT)
      PARAMETER (NNR=353,NNC=262)
	REAL XYMELT(NNR,NNC),awrite(336420,2),asum(336420)
	data asum/336420*0./
C        
	NW=0
	DO 10 I=1,NR
	   DO 11 J=1,NC
	nw=nw+1

	asum(nw)=asum(nw)+xymelt(i,j)*3.6/1000.
11         CONTINUE 
10      CONTINUE   
	if (istep.eq.25) then
	   NW=0
	   DO I=1,NR
	     DO J=1,NC
	        NW=NW+1

       	    WRITE(IUOUT,120) asum(nw)
	     enddo
	   enddo
	endif
100     FORMAT(F6.3)
120     FORMAT(25(1x,F7.2))
	RETURN
	END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	SUBROUTINE OUTP3(NR,NC,XYMELT,istep,IUOUT)
      PARAMETER (NNR=353,NNC=262)
	REAL XYMELT(NNR,NNC),awrite(336420,2),asum(336420)
	data asum/336420*0./
C        
	NW=0
	DO 10 I=1,NR
	   DO 11 J=1,NC
	nw=nw+1

	asum(nw)=asum(nw)+xymelt(i,j)*3.6/1000.
11         CONTINUE 
10      CONTINUE   
	if (istep.eq.25) then
	   NW=0
	   DO I=1,NR
	     DO J=1,NC
	        NW=NW+1

       	    WRITE(IUOUT,120) asum(nw)

	     enddo
	   enddo
	endif
100     FORMAT(F6.3)
120     FORMAT(25(1x,F7.2))
130     FORMAT(25(1x,F10.3))
	RETURN
	END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	SUBROUTINE OUTP4(NR,NC,XYMELT,istep,IUOUT)
      PARAMETER (NNR=353,NNC=262)
	REAL XYMELT(NNR,NNC),awrite(336420,2)
C        
	NW=0
	DO 10 I=1,NR
	   DO 11 J=1,NC
	nw=nw+1
	awrite(nw,istep)=xymelt(i,j)
11         CONTINUE 
10      CONTINUE   

	   NW=0
	   DO I=1,NR
	     DO J=1,NC
	        NW=NW+1

       	    WRITE(IUOUT,120) awrite(nw,1)
	     enddo
	   enddo
CGG	endif
100     FORMAT(F6.3)
120     FORMAT(4(1x,F7.2))
	RETURN
	END
C
C
	SUBROUTINE OUTPUS(NR,NC,HSNPCK,HEABXY)
      PARAMETER (NNR=353,NNC=262)
	REAL HSNPCK(NNR,NNC),HEABXY(NNR,NNC),TEMP(NNR)
	COMMON/RADIA/SWAVXY(nnr,nnc),RLOWXY(nnr,nnc),HCONVXY(nnr,nnc),
     +HLATXY(nnr,nnc),HPXY(nnr,nnc),SNOGRXY(nnr,nnc)
C        
c        WRITE(19,*) ' HEABXY [W/mq]' 
	DO 10 NY=1,NC
	   DO 11 I=1,NR
	      IF (HSNPCK(I,NY).LE.0) THEN
		  TEMP(I)= HEABXY(I,NY)
		  HEABXY(I,NY)=-999
	      ENDIF
11         CONTINUE 
c           WRITE(19,100) (HEABXY(NX,NY),NX=1,NR)
	   DO 12 I=1,NR
12            IF (HSNPCK(I,NY).LE.0) HEABXY(I,NY)=TEMP(I)
10      CONTINUE
c        WRITE(19,*) ' SWAVXY [W/mq]' 
	DO 20 NY=1,NC
	   DO 21 I=1,NR
	      IF (HSNPCK(I,NY).LE.0) THEN
		  TEMP(I)= SWAVXY(I,NY)
		  SWAVXY(I,NY)=-999
	      ENDIF
21         CONTINUE 
c           WRITE(19,100) (SWAVXY(NX,NY),NX=1,NR)
	   DO 22 I=1,NR
22            IF (HSNPCK(I,NY).LE.0) SWAVXY(I,NY)=TEMP(I)
20      CONTINUE
c        WRITE(19,*) ' RLOWXY [W/mq]' 
	DO 30 NY=1,NC
	   DO 31 I=1,NR
	      IF (HSNPCK(I,NY).LE.0) THEN
		  TEMP(I)= RLOWXY(I,NY)
		  RLOWXY(I,NY)=-999
	      ENDIF
31         CONTINUE 
c           WRITE(19,100) (RLOWXY(NX,NY),NX=1,NR)
	   DO 32 I=1,NR
32            IF (HSNPCK(I,NY).LE.0) RLOWXY(I,NY)=TEMP(I)
30      CONTINUE
c        WRITE(19,*) ' HCONVXY [W/mq]' 
	DO 40 NY=1,NC
	   DO 41 I=1,NR
	      IF (HSNPCK(I,NY).LE.0) THEN
		  TEMP(I)= HCONVXY(I,NY)
		  HCONVXY(I,NY)=-999
	      ENDIF
41         CONTINUE 
c           WRITE(19,100) (HCONVXY(NX,NY),NX=1,NR)
	   DO 42 I=1,NR
42            IF (HSNPCK(I,NY).LE.0) HCONVXY(I,NY)=TEMP(I)
40      CONTINUE
c        WRITE(19,*) ' HLATXY [W/mq]' 
	DO 50 NY=1,NC
	   DO 51 I=1,NR
	      IF (HSNPCK(I,NY).LE.0) THEN
		  TEMP(I)= HLATXY(I,NY)
		  HLATXY(I,NY)=-999
	      ENDIF
51         CONTINUE 
c           WRITE(19,100) (HLATXY(NX,NY),NX=1,NR)
	   DO 52 I=1,NR
52            IF (HSNPCK(I,NY).LE.0) HLATXY(I,NY)=TEMP(I)
50      CONTINUE
c        WRITE(19,*) ' HPXY [W/mq]' 
	DO 60 NY=1,NC
	   DO 61 I=1,NR
	      IF (HSNPCK(I,NY).LE.0) THEN
		  TEMP(I)= HPXY(I,NY)
		  HPXY(I,NY)=-999
	      ENDIF
61         CONTINUE 
c           WRITE(19,100) (HPXY(NX,NY),NX=1,NR)
	   DO 62 I=1,NR
62            IF (HSNPCK(I,NY).LE.0) HPXY(I,NY)=TEMP(I)
60      CONTINUE
c        WRITE(19,*) ' SNOGRXY [W/mq]' 
	DO 70 NY=1,NC
	   DO 71 I=1,NR
	      IF (HSNPCK(I,NY).LE.0) THEN
		  TEMP(I)= SNOGRXY(I,NY)
		  SNOGRXY(I,NY)=-999
	      ENDIF
71         CONTINUE 
c           WRITE(19,100) (SNOGRXY(NX,NY),NX=1,NR)
	   DO 72 I=1,NR
72            IF (HSNPCK(I,NY).LE.0) SNOGRXY(I,NY)=TEMP(I)
70      CONTINUE
100     FORMAT(20(F7.1,1X))
c        WRITE (19,*)
	RETURN
	END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
C    

    
	SUBROUTINE STATT(NSTSWA,NSTRA,SWAMUT,SWASTT,RAMUT,RASTT,
     &HEMUT,HESTT,MEMUT,MESTT,HESMUT,HESSTT,MESMUT,MESSTT)
CCC
C       OUTPUT
C
C       NSTSWA  N. of hours of presence of snow over sensor grid
C       NSTRA   N. of hours of presence of snow over the basin
C       SWAMUT  Temporal mean of net shortwave over sensor grid
C       SWASTT  Temporal standard deviation of net shortwave over sensor grid
C       RAMUT   Temporal mean of net shortwave over the basin +SNOW ON SENSOR
C       RASTT   Temporal std dev of net shortwave over the basin+SNOW ON SENSOR
C       HEMUT   Temporal mean of heat balance over the basin + SNOW ON SENS
C       HESTT   Temporal standard deviation of heat balance over the basin +S
C       MEMUT   Temporal mean of snow melt over the basin + SNOW ON SENSOR
C       MESTT   Temporal standard deviation of snow melt over the basin +S
C       HESMUT  Temporal mean of heat balance over sensor grid
C       HESSTT  Temporal std dev of heat balance over sensor grid
C       MESMUT  Temporal mean of snow melt over sensor grid
C       MESSTT  Temporal std dev of snow melt over sensor grid
CCC
C       COMMON Variables Used
C       NTIMST    N. of Time steps
C       XX(i,j,t) j=1  RAMUXY SHORTWAVE AVERAGED OVER THE BASIN AT TIME T
C                   2  HEMUXY HEAT BALANCE AVERAGED OVER THE BASIN AT T
C                   3  SWAVSG SHORTWAVE AT STATION AT TIME T
C                   4  MELT    MELT (+ RAIN) AVERAGED OVER THE BASIN AT TIME T
C                   5  HEABAS  HEAT BALANCE AT STATION
C                   6  MELTS   MELT AT STATION AT TIME T
C       H1NPCK    Presence of snow flag
CCC
      PARAMETER (NNR=353,NNC=262)
      INTEGER T 
      REAL MEMUT,MESMUT,MESTT,MESSTT
      REAL DEM(NNR,NNC),XNORM(3,NNR,NNC),SWST(2,NNR,NNC),
     +PDKEF(2,NNR,NNC)
      REAL PACT(2,NNR,NNC),HSNPCK(NNR,NNC),DENSN(2,NNR,NNC)
      COMMON /DEM/ DEM,XNORM,HSENSO,NXSENS,NYSENS,NR,NC,DEMDX,DEMDY,
     +AREAA,NGRIDS
      COMMON /SNPK/ SWST,PDKEF,PACT,HSNPCK,DENSN
      COMMON /RECORD/ XX(10,6,8784),NTIMST
CCC
C       Estimation of Temporal mean and standard deviation
CCC
	NSTSWA=0
	NSTRA=0
	SWAMUT=0.
	SWASTT=0.
	RAMUT=0.
	RASTT=0.
	HEMUT=0.
	HESTT=0.
	MEMUT=0.
	MESTT=0.
	HESMUT=0.
	HESSTT=0.
	MESMUT=0.
	MESSTT=0.
	DO 10 T=1,NTIMST
	  IF(XX(1,2,T).NE.0) then
	   IF(XX(1,5,T).NE.0) then
	      NSTSWA=NSTSWA+1
	      SWAMUT=SWAMUT+XX(1,3,t)
	      RAMUT=RAMUT+XX(1,1,t)
	      HEMUT=HEMUT+XX(1,2,t)
	      MEMUT=MEMUT+XX(1,4,t)
	      HESMUT=HESMUT+XX(1,5,t)
	      MESMUT=MESMUT+XX(1,6,t)
	   ENDIF
	   NSTRA=NSTRA+1
	  ENDIF
10      continue
CCC
C       Temporal MEAN
CCC
	IF(NSTSWA.GT.0) THEN
	  SWAMUT=SWAMUT/FLOAT(NSTSWA)  
	  HESMUT=HESMUT/FLOAT(NSTSWA)  
	  MESMUT=MESMUT/FLOAT(NSTSWA)  
	  RAMUT=RAMUT/FLOAT(NSTSWA)  
	  HEMUT=HEMUT/FLOAT(NSTSWA)  
	  MEMUT=MEMUT/FLOAT(NSTSWA)  
	ENDIF 
CC   STATNDATD DEVIATION
	DO 20 T=1,NTIMST
	   IF(NSTSWA.GT.0) THEN 
	      SWASTT=SWASTT+(XX(1,3,t)-SWAMUT)**2
	      HESSTT=HESSTT+(XX(1,5,t)-HESMUT)**2
	      MESSTT=MESSTT+(XX(1,6,t)-MESMUT)**2
	      RASTT=RASTT+(XX(1,1,t)-RAMUT)**2
	      HESTT=HESTT+(XX(1,2,t)-HEMUT)**2
	      MESTT=MESTT+(XX(1,4,t)-MEMUT)**2
	    ENDIF
20      continue
CCC
C       Temporal Standard Deviation
CCC
	IF(NSTSWA.GT.1) THEN
	  SWASTT=SQRT(SWASTT/FLOAT(NSTSWA-1))  
	  HESSTT=SQRT(HESSTT/FLOAT(NSTSWA-1))  
	  MESSTT=SQRT(MESSTT/FLOAT(NSTSWA-1))  
	  RASTT=SQRT(RASTT/FLOAT(NSTSWA-1))  
	  HESTT=SQRT(HESTT/FLOAT(NSTSWA-1))  
	  MESTT=SQRT(MESTT/FLOAT(NSTSWA-1))  
	ENDIF 
	RETURN
	END
C
      SUBROUTINE SETUP
    
C
C  [1] Parameter set up 
C  N= N. of Parameters to be optimized
C  X() Parameter to be optimized
C
      PARAMETER (NNR=353,NNC=262)
      REAL VVCFCT(12),KATM
      COMMON /PARAME/ XS(20),X(20),INDX(20),N,NPAR,STEPS(20),TESTPS(20),
     +TESSTP(20),STEP(20)
      COMMON /SOIL/ ASOIL,SLAMDA,Tmsoil,Cv(NNR,NNC),DSOIL1(NNR,NNC)
      COMMON /SNOW/A0,A1,TSPRIN,ASPRIN,TRNCF,AM0,CAM,GANESN,ALFAS,
     +BETAS,STRN,CJOMME,HSNLA1,FHI,Swi,FHIeKs,BEERSN,CMELTM,SPCFCT
      COMMON /ATMOSP/ DOX,SATM0,TRAATM,PERCLO,VVCFCT,KATM
      COMMON /HEAT/CPw,Wden,ROhu,Cpa,RatiKh,CKa,Z0,VLh,RatiKe,Rgas,ZW,ZT
      WRITE(*,*) ' SETUP'
      NZ=0
      DO 250 I=1,NPAR
      IF(INDX(I).NE.1) GOTO 250
	NZ=NZ+1
	XS(I)=X(I)
250   CONTINUE
C     Snow parameters
C   Z0,A0,A1,TSPRIN,ASPRIN,TRNCF,AM0,CAM,STRN,ALFAS,BETAS
      Z0=XS(1)
      A0=XS(2)      
      A1=XS(3)      
      TSPRIN=XS(4)
      ASPRIN=XS(5)
      TRNCF=XS(6)
      AM0=XS(7)
      CAM=XS(8)                                               
      STRN=XS(9)
      ALFAS=XS(10)
      BETAS=XS(11) 
      Swi=XS(12)
C  Soil parameters
      ASOIL=XS(13)
      SLAMDA=XS(14)
      Tmsoil=XS(15)
C  ATmosphere parameters
      SATM0=XS(16)
      PERCLO=XS(17)
      TRAATM=XS(18)
	KATM=XS(19)
C     Multiplicative factor for melt celerity
      CMELTM=XS(20)
      RETURN
      END
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C       CCCCCCCCCCCC    PROGRAMMA   S U N S E T    CCCCCCCCCCCCC
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C       IL FILE INPUT SI CHIAMA SOLE.DAT ed e` compatibile con i
C       programmi ELIOS*.FOR. Esso deve contenere i seguenti records:
C       
C       Record  1:      Titolo                                   A80
C       Record  2:      Latitudine,Longitudine,Angolo fuso orario
C                       Le coordinate geografiche vanno espresse 
C                       in gradi e primi separati dal punto.      F*
C
C
	SUBROUTINE SUNSET(XLAT,XLONG,STD,DEL,DELL,PASEC)

	REAL RADN(3)
	CHARACTER LOCAT*80
       COMMON /SUNS/ SCI,SCE,ANRAD,RADAN,PI,PI12,PI2,PASSO,SOLST,CALJOU
C
C
	WRITE(*,*)  ' SUNSET'
C       SCI=1353.
	SCE=1162.9
	ANRAD=.0174533
	RADAN=57.29578
C
C
	CALJOU=4.1868 
	SOLST=23.43*ANRAD
	PI=4.*ATAN(1.)
	PI12=12./PI
	PI2=PI/2.
C
C       **************************************************
C
	READ(2,9) LOCAT
9       FORMAT(A)
C
	READ(2,*) PASSO
	READ(2,*)(RADN(I),I=1,3)
C
C       ***  CONVERSIONE GRADI/RADIANTI DELLE COORDINATE GEOGRAF. ***
C
	XLAT=((RADN(1)-INT(RADN(1)))*1.66667+FLOAT(INT(RADN(1))))*
     #ANRAD
	XLONG = RADN(2)
	STD   = RADN(3)
	DEL=XLONG-STD                                                   
	DELL=DEL/15.                                                    
	PASEC=3600.*PASSO                                               
C
C   Input PARAMETER FILE IS CLOSED
C
c        CLOSE(2)
	RETURN
	END
C
C
C
	SUBROUTINE SOL(DATE,DAYS,ET,ALFA,R)
    
C 
C       ET     Correzione ora [ore]
C       TCMMTV Correzione tempo civile medio-tempo vero[ore]
C
      COMMON /SUNS/ SCI,SCE,ANRAD,RADAN,PI,PI12,PI2,PASSO,SOLST,CALJOU
C        
C      WRITE (*,*) 'SOL'  
      JDAY=INT(DAYS+DATE)
C
C
C       Da Price,J. C.,Photogrammetric engineering 1982 (10)
C
      TCMMTV=(8*SIN(2.*PI*(JDAY-2)/365.)+
     +10*SIN(4.*PI*(JDAY-355)/365.))/60.
      ET=-TCMMTV
      DAYX=(DAYS+DATE-1.)*6.2831853/365.256

	ALFAD=.302-22.93*COS(DAYX)-.229*COS(2.*DAYX)-.243*COS(3.*DAYX)
     $+3.851*SIN(DAYX)+.002*SIN(2.*DAYX)-.055*SIN(3.*DAYX)

	ALFA=ALFAD*.0174533
C
	R= 1.-.0167*COS(DAYX)
C
	RETURN
	END
C
	SUBROUTINE DAILY(ALFA,XLAT,TR,TS,BETAM,RISANG,PSIR)
    
C
C       *********************************************************
C       *   INPUT       ALFA  = ANGOLO DECLINAZIONE SOL.  (RAD.)*
C       *               XLAT  = LATITUDINE                (RAD.)*
C       *                                                       *  
C       *   OUTPUT      TR    = TEMPO LOC. STAND. CORR. ALBA    *
C       *               TS    =  ''    ''   ''     ''   TRAMONTO*
C       *               RISANG= ANGOLO ORARIO ALBA       (RAD.) *
C       *               PSIR  = AZIMUT ALBA (RAD.)              *
C       *               PSIS  = AZIMUT TRAMONTO (RAD.)          *
C       *               BETAM = ELEVAZ.MASSIMA DEL GIORNO (RAD.)*
C       *********************************************************
C
       COMMON /SUNS/ SCI,SCE,ANRAD,RADAN,PI,PI12,PI2,PASSO,SOLST,CALJOU
	ACOSF(X)=PI/2.-ASIN(X)
C        WRITE (*,*) 'DAILY' 
	ALFAD=ALFA*RADAN
C
	SINA=SIN(PI2-ALFA)
	COSA=COS(PI2-ALFA)
	TANA=SINA/COSA
C
	SINFI=SIN(PI2-XLAT)
	COSFI=COS(PI2-XLAT)
	TANFI=SINFI/COSFI
C
	SINBI=SIN(BETAM)
C
C       ** ORA  ALBA E TRAMONTO **
C
	CONG=ACOSF(-TAN(ALFA)*TAN(XLAT))
C       RISANG=PI-CONG
	RISANG=CONG
C
	BETAM=PI/2.-ACOSF(-SINFI*SINA+COSFI*COSA)
	TR=12.-RISANG*PI12
	TS=24.-TR
C
C       ** AZIMUT ALBA E TRAMONTO **
C
	ALAMBD=ATAN(TANA*COS(RISANG))
	TANZ=TAN(RISANG)*SIN(ALAMBD)/SIN(PI2-XLAT-ALAMBD)
	PSIR=ATAN(TANZ)
	IF(RISANG.EQ.PI2.AND.TANZ.LE.0.)THEN
	PSIR=PSIR+PI
	GO TO 9
	END IF
	IF(RISANG.GT.0.AND.RISANG.LT.PI.AND.PSIR.LE.0.)PSIR=PSIR+PI
 9      PSIS=2.*PI-PSIR
C
	RETURN
	END
C
C
	SUBROUTINE ORA(XLAT,ALFA,AST,PSIZ,THETR)
    
C       ************************************************
C       *** INPUT  ***  XLAT = LATUTUDINE              *
C       *               ALFA = DECLINAZIONE            *
C       *               AST  = TEMPO STANDARD CORRETTO *
C       *                                              *
C       *** OUTPUT ***  PSIZ = ANGOLO AZIMUTALE (RAD.) *
C       *               THETR= ANGOLO ELEVAZIONE (RAD.)*
C       ************************************************
C
	LOGICAL T1,T2,TT1,TT2,TAN1,TAN2
       COMMON /SUNS/ SCI,SCE,ANRAD,RADAN,PI,PI12,PI2,PASSO,SOLST,CALJOU
CMUL
       ACOSF(X)=1.5707693-ASIN(X)
C        WRITE (*,*) 'ORA'  
C
	SINA=SIN(PI2-ALFA)
	COSA=COS(PI2-ALFA)
	TANA=SINA/COSA
	TAND=SIN(ALFA)/COS(ALFA)
C
	SINFI=SIN(PI2-XLAT)
	COSFI=COS(PI2-XLAT)
C
C
	PT=(PI*(12.-AST)/12.)
C
C       ** ANGOLO ELEVAZIONE **
C
	COSB=SINA*SINFI*COS(PT)+COSA*COSFI
	BETA=ACOSF(COSB)
	THETR=PI2-BETA
C
C       ** ANGOLO AZIMUTALE **
C
	TT1=(PT.EQ.PI2)
	TT2=(PT.EQ.-PI2)
	IF(TT1.OR.TT2)GO TO 9
	ALAMBD=ATAN(TANA*COS(PT))
	TANZ=TAN(PT)*SIN(ALAMBD)/SIN(PI2-XLAT-ALAMBD)
	PSIZ=ATAN(TANZ)
	GO TO 19
C
 9      TANZ=1./(SINFI*TAND)
	PSIZ=ATAN(TANZ)
	TAN1=(TANZ.GT.0.)
	TAN2=(TANZ.LT.0.)
	IF(TT1.AND.TAN1)PSIZ=PSIZ
	IF(TT1.AND.TAN2)PSIZ=PSIZ+PI
	IF(TT2.AND.TAN1)PSIZ=2.*PI-ATAN(TANZ)
	IF(TT2.AND.TAN2)PSIZ=PI-PSIZ
	RETURN
C
  19    TAN1=(TANZ.GT.0.)
	TAN2=(TANZ.LT.0.)
	T1=(PT.LT.PI.AND.PT.GE.0.)
	T2=(PT.LT.0.AND.PT.GT.-PI)
	IF(T1.AND.TAN1)PSIZ=PSIZ
	IF(T1.AND.TANZ.LE.0.)PSIZ=PSIZ+PI
	IF(T2.AND.TAN1)PSIZ=PSIZ+PI
	IF(T2.AND.TAN2)PSIZ=2.*PI+PSIZ
C
	RETURN
	END
C
C       FUNCTION TLAPSE(Tasns,Hgrid,Hsenso) Determination of temperatures
C       variation with altitude
C       GAMMAT Temperature Lapse rate [.ø.C/m] 
C       (computed in subroutine dtweight)
C
	FUNCTION TLAPSE(Tasns,HGRID,Hsenso,gammat)
	TLAPSE=Tasns+GAMMAT*(HGRID-Hsenso)
	RETURN
	END
C
C       FUNCTION PRLAPSE(Pmed,Hgrid,Hsenso) Determination of precipitation
C       variation with altitude
C       PREC0 Annual Precipitation at 0 m a.s.l. [mm]
C       GAMMAPR Annual Precipitation Lapse rate [mm/m]
C

c	 DATA PREC0,GAMMAPR /717.0,0.312/ For Adamelo Glacier 
c      For Naltar We assumed no altitudinal gradient 

	FUNCTION PRLAPSE(HGRID,hspmis)
c	DATA PREC0,GAMMAPR /717.0,0.312/
c	 PRLAPSE=((PREC0+GAMMAPR*HGRID)/(PREC0+GAMMAPR*hspmis))
       PRLAPSE=1.
	RETURN
	END
C
C       FUNCTION PLAPSE(Pmed,Hgrid,Hsenso) Determination of pressure
C       variation with altitude
C       GAMMAP Pression Lapse rate [mbar/m]
C
	FUNCTION PLAPSE(Pmed,HGRID,Hsenso)
	DATA GAMMAP /0.11772/
	 PLAPSE=PMED-GAMMAP*(HGRID-Hsenso)
	RETURN
	END
C
C       FUNCTION RAINSN(Ta,P,Rain) Determination of Form of Precipitation
C              (See CEQUEAU,UBC,PRMS or Soncini Models)                 
C
C       Ta    =Input;Air Temperature          [.ø.C]
C       P     =Input;Precipitation (Water equivalent) [mm]
C       Rain  =Output; Rain precipitation    "        [mm]
C       Rainsn=Output; Snow precipitation    "        [mm] 
C
C       STRN        Threshold-type parameter  [.ø.C]
C 
	FUNCTION RAINSN(Ta,P,Rain)
      COMMON /SNOW/A0,A1,TSPRIN,ASPRIN,TRNCF,AM0,CAM,GANESN,ALFAS,
     +BETAS,STRN,CJOMME,HSNLA1,FHI,Swi,FHIeKs,BEERSN,CMELTM,SPCFCT
	TCOMP=Ta-STRN
	IF((TCOMP.GT.-1.).AND.(TCOMP.LT.1.)) THEN
	 ALFA=(TCOMP+1)/2
	 RAIN=ALFA*P
	 RAINSN=P-RAIN
	 RETURN
       ELSE IF(TCOMP.GE.1.) THEN
	 RAIN=P
	 RAINSN=0.
	 RETURN
       ELSE IF(TCOMP.LE.-1.) THEN
	 RAINSN=P
	 RAIN=0.
	 RETURN
       ENDIF
       RETURN
       END
		
C
C       RLics       Incoming Long-wave radiation under clear-sky [W/mq]
C                   Open areas (no forest)
C       Ta          Air Temperature         [.ø.C]      
C       Ea          Air Vapor pressure      [mb]
C       STEFAA      Stefan-Boltzman Law (air temperature)    [W/mq]
C
	FUNCTION RLics(Ta,EA)
	  COMMON /STEF/ STEFAA
	  RLiCS=EPSA(Ta,EA)*STEFAA
	RETURN
	END
C
C       Hp          Heat advected by rain at Temperature equal to air
C                   (Rain is freezed to 0° C Temperature   [J/mq]
C       Cpw         Specific Heat of water [J/(Kgm °C)]
C       Wden        Water density          [Kgm/mc] 
C       P           Rainfall               [mm]
C       Ta          Air Temperature        [°C]
C       Tsnow       Snow Surface Temperature       [°C]
C       
	FUNCTION HP(P,Ta,Tsnow)
      COMMON /HEAT/CPw,Wden,ROhu,Cpa,RatiKh,CKa,Z0,VLh,RatiKe,Rgas,ZW,ZT
	  HP=Cpw*P/1000*Wden*Ta
	RETURN
	END
C
C       Hc          Convective Heat       [W/mq]
C
C       Dw          Convective term       [m/s]
C       Ta          Air Temperature       [°C]
C       Tsnow       Snow Surface Temperature       [°C]
C
C       ROhu        Density of humid air  [Kgm/mc]
C       Cpa         Specific heat of air at constant pressure [J/(Kgm °C)]
C       RatiKh      Kh/Km                 [ ]
C
	FUNCTION Hc(Dw,Ta,Tsnow)
      COMMON /HEAT/CPw,Wden,ROhu,Cpa,RatiKh,CKa,Z0,VLh,RatiKe,Rgas,ZW,ZT
	  HC=ROhu*Cpa*RatiKh*Dw*(Ta-Tsnow)
	RETURN
	END
C
C       Hl          Latent Heat       [W/mq]
C
C       Evr         Evaporation rate            [Kgm/(mq s)]
C
C       VLh         Latent heat of vaporization [J/Kgm]
C
	FUNCTION Hl(Evr)
      COMMON /HEAT/CPw,Wden,ROhu,Cpa,RatiKh,CKa,Z0,VLh,RatiKe,Rgas,ZW,ZT
	  Hl=-VLh*Evr
	RETURN
	END
C
C       Ev          Evaporation rate      [Kgm/(m2 s)]
C                   Vapor pressure at snow surface is assumed equal to
C                   saturation vapor pressure at Snow Temperature (.ø.0=6.11 mb)
C
C       Dw          Convective term       [m/s]
C       EA(Ta,Rhumi) Absolute Vapor Pressure       [mb]      
C       (  EA=Rhumi/100*EAS(Ta)  )
C       Ta          Air Temperature       [.ø.C]
C       Tsnow       Snow Surface Temperature       [.ø.C]
C
C       RatiKe      Ke/Km                 [ ]
C       Rgas        Constant of perfect gas (air)   [m3 mb/(Kgm K)]
C 
	FUNCTION Ev(Dw,Ea,Ta,Tsnow)
      COMMON /HEAT/CPw,Wden,ROhu,Cpa,RatiKh,CKa,Z0,VLh,RatiKe,Rgas,ZW,ZT
c         WRITE (*,*) 'EA,TSNOW,TA=',EA,TSNOW,TA
c         WRITE (*,*) 'RATIKE,DW,RGAS=',RATIKE,DW,RGAS
c         WRITE (*,*) 'EAS=',EAS(TSNOW)
	  Ev=-0.622*RatiKe*Dw*(Ea-Eas(Tsnow))/(Rgas*(Ta+273.14))
c20060206	Ev=-0.622*RatiKe*Dw*(Ea-Eas(Ta))/(Rgas*(Ta+273.14))
	RETURN
	END
C
C       Function Potevr(DAYCUR)  Potential evaporation from soil [mm/hr]
C       DAYCUR =Current day from beginning of year
C       EVMEHO =Hourly evaporation rate (annual mean=403*1.1616 mm ) [mm/hr] 
C       PI2365 =2PIgr/365
C
	FUNCTION Potevr(DAYCUR)
	DATA EVMEHO,PI2365 /0.05344,0.017214/
	  Potevr=EVMEHO*(1.+SIN((DAYCUR-90.)*PI2365))       
	RETURN
	END
C
C
C       Dw          Convective term          [m/s]
C                   Correction for lapse-neutral-stability condition is done
C                   applying Richardson number procedure (Price,Dunne,WRR,1976)
C                   And Male-Granger formula+data(WRR,1981) 
C
C       Ta          Air  Temperature               [.ø.C]
C       Tsnow       (Snow) Surface Temperature       [.ø.C]
C       Hsnst       Height of snowpack at station [m]
C       Vw          Wind speed                    [m/s]
C
C       ZW          Height of wind speed sensor   [m]
C       K=CKa       Von Karman constant 0.41 in area COMMON /HEAT/
C       Z0surf      Surface roughness                [m]  Z0 from parameter file, or 0.01 or 0.007 for debris (Regine Hock)
C
	FUNCTION Dw(Vw,Ta,Tsnow,Hsnst,z0surf)
      COMMON /HEAT/CPw,Wden,ROhu,Cpa,RatiKh,CKa,Z0,VLh,RatiKe,Rgas,ZW,ZT
	  IF(Vw.EQ.0) THEN
	    Dw=0.
        ELSE
CGGG	   Dw=CKa*CKa*Vw*(ZW-Hsnst)/(ZT-Hsnst)/(Alog((Zw-Hsnst)/Z0))**2
C2.8	   Dw=CKa*CKa*Vw*(ZW-Hsnst)/(ZT-Hsnst)/(Alog((Zw-Hsnst)/Z0))**2
C2.8
	Dw=CKa*CKa*Vw/(Alog((Zw-Hsnst)/z0surf))**2
      ENDIF
	RETURN
	END
C                                       
C       EPSA:       Air emissivity [ ] following Satterlund eq. in WRR,1979
C       Ta          Air Temperature [.ø.C]      
C       EA(T,Rhumi) Absolute Vapor Pressure [mb]      
C       EA=Rhumi/100*EAS(Ta)
C
	FUNCTION EPSA(Ta,EA)        
	  EPSA=1.08*(1-exp(-EA**((273.14+Ta)/2016)))
	RETURN
	END
C
C       EAS:        Saturated Vapor Pressure  [mb]
C       Ta          Air Temperature [.ø.C]      
C
	FUNCTION EAS(Ta)
C          WRITE (*,*)'TA=',TA
	  if (ta.lt.-50.) ta=-50.
	  EAS=6.11*EXP(17.27*Ta/(237.3+Ta))
	RETURN
	END
C
C       RLio Long-wave radiation balance for forested/not and
C           Clouds/no conditions   [W/mq]
C       FOR   = Percentage of forested area cover 
C       Cloudi= Percentage of diffused radiation (Assumed as cloudiness)
C       Ta    =Air temperature [.ø.C]
C       Tsnow =Snow Surface Temperature       [.ø.C]
C       EA(T,Rhumi) Absolute Vapor Pressure [mb]      
C  
	FUNCTION RLiO(COVER,Cloudi,Ta,Tsnow,EA)
	 COMMON /STEF/ STEFAA
	  Emiss=(1-Cloudi)*EPsa(Ta,EA)+Cloudi
	if (emiss.lt.0.01) EMISS=0.01
C	 RLiO=(1-COVER)*EMISS*Stefaa+COVER*Stefaa-.99*STEFAN(Tsnow+273.14)
	 RLiO=(1-COVER)*EMISS*Stefaa+COVER*(EMISS-1.)*Stefaa
     +-(1.-COVER)*.99*STEFAN(Tsnow+273.14)
	RETURN
	END
C
C       GminS0   Balance of Shortwave radiation [W/Mq]
C
C       Solrad= Incoming shortwave radiation [W/mq]
C       THETR = Sun's elevation angle        [rad]
C       ALBEDO= Open areas albedo [ ]
C       FOR   = Percentage of forested area  [0,1]
C       ALBEDF= Forested areas albedo [ ]
C 
C       TRNCF = Transmission of shortwave radiation by forest (.15 by Federer)
C
	FUNCTION GminS0(Solrad,ALBEDO,FOR,ALBEDF)
      COMMON /SNOW/A0,A1,TSPRIN,ASPRIN,TRNCF,AM0,CAM,GANESN,ALFAS,
     +BETAS,STRN,CJOMME,HSNLA1,FHI,Swi,FHIeKs,BEERSN,CMELTM,SPCFCT
c	SolrAF=Solrad*TRNCF
	GminS0=FOR*(1.-ALBEDF)*Solrad+(1.-FOR)*(1.-ALBEDO)*Solrad
	RETURN
	END
C     
C       Hgroun    Heat exchange snowpack-ground [W/mq]
C       
C       Approximation 2.
C       2.   Analitical solution for heat flow through semi-infinite plane
C       without phase change with initial temperature (before snowfall) equal
C       to mean yearly air temperature and boundary condition T=0 .ø.C (snowmelt)
C       Formula (5.49)
C       See Bonacina-Cavallini-Mattarolo Trasmisione del calore Pag. 130..134
C
C       TIMENO = Hours after 15th of November (assumed first snowfall)  [hr]
C
C       ASOIL  = Soil thermal diffusivity (.001486 Eagleson-Bonacina) [mq/h]         
C       SLAMDA = Soil thermal conductivity (1.557)       [W/(m .ø.C)]
C       Tm     = Annual Mean Air temperature      (6.2)         [.ø.C]
C
	FUNCTION Hgroun(TIMENO)
	PARAMETER(NNR=353,NNC=262)
      COMMON /SOIL/ ASOIL,SLAMDA,Tmsoil,Cv(NNR,NNC),DSOIL1(NNR,NNC)
	DATA SQRPI /1.7725/
	  HGROUN=SLAMDA*Tmsoil/SQRPI/SQRT(ASOIL*TIMENO)
	RETURN
	END
C
C     This function computes the Stefan-Boltzman Law
C
C     T      Absolute Temperature [K]
C     STEFAN Global emission of a black body at temperature T [W/mq]
C
      FUNCTION STEFAN(T)
    
      DATA SIGMA /5.67/
	STEFAN=SIGMA*(T/100)**4.
      RETURN
      END
C
C
C     ALBEDO CLEAR-SKY (First days after snowfall and before metamorphosis)
C
C     T2112= Time [days] after 21st December
C
C     A0=Winter time albedo (21/12)-SPRING time albedo   (0.08)
C     A1=Parameter                    (0.16)
C     TSPRIN=Time when spring time albedo is fixed (106 gg)     
C     ASPRIN=Spring time albedo       (0.735)
C
	  FUNCTION ALBCLS(T2112)
      COMMON /SNOW/A0,A1,TSPRIN,ASPRIN,TRNCF,AM0,CAM,GANESN,ALFAS,
     +BETAS,STRN,CJOMME,HSNLA1,FHI,Swi,FHIeKs,BEERSN,CMELTM,SPCFCT
	   ALBCL=ASPRIN-A0/(1.5708)*ATAN(A1*(T2112-TSPRIN))
	   ALBCLS=AMIN1(ALBCL,1.)
	  RETURN
	  END
C
C     RADSUN    COMPUTES COMPONENTS OF SOLAR RADIATION
C
C    Input
C     SOLRAD=Solar radiation measured by sensor        [W/mq]
C     THETR = Solar altitude       [rad]
C    Output
C     SolDIF=Diffused radiation    [W/mq]
C     SolBEA=Direct beam radiation [W/mq]
C     DIFRAD=Cloudiness/Diffused radiation in excess respect theorical[0,1]
C     RCLSKY=Clear Sky radiation model
C
C     DOX   = Solar constant in 0.3-3 micron band [ca. 1300 W/mq]
C     SATM0  =Atmosphere absorbtion function at sea level (0.33) [ ]
C     TRAATM=Transmission coefficient of diffused radiation (0.2) [ ]
C     PERCLO=Percentage of maximum radiation in overcast condition 0.25 [ ]
C	KATM = Altitudinal decay function of atmospheric optical depth [m]
C
CGG	SUBROUTINE RADSUN(Solrad,THETR,SolDIF,SolBEA,DIFRAD,RCLSKY)
	SUBROUTINE RADSUN(Solrad,THETR,SolSKYDIF,SolBEA,DIFRAD,RCLSKY,
     +vfsens,asens,hsenso,SATM)
	REAL VVCFCT(12),KATM
      COMMON /ATMOSP/ DOX,SATM0,TRAATM,PERCLO,VVCFCT,KATM
C        WRITE (*,*) 'RADSUN'  
	  IF(THETR.LE.0) THEN
	    DIFRAD=1.
C	    SOLDIF=SOLRAD
	    SOLSKYDIF=SOLRAD
	    SOLBEA=0.
	    RCLSKY=0.
C
	    RETURN
	  ELSE   
c radiation model
CGG  SATM is now dependent upon sensor cell elevation, computed in the 
	    SATM=SATM0*EXP(-HSENSO/KATM)
	    EXTINC=EXP(-SATM/SIN(THETR))*SIN(THETR)
	    RCLBEA=DOX*EXTINC
CGG	    RCLDIF=TRAATM*(DOX*SIN(THETR)-RCLBEA)
	    RCLSKYDIF=TRAATM*(DOX*SIN(THETR)-RCLBEA)
	    RCLSENSDIF=((1.-vfsens)*RCLSKYDIF+asens*vfsens*RCLBEA)/
     +(1.-asens*vfsens)
CGG	    RCLSKY=RCLBEA+RCLDIF 
	    RCLSKY=RCLBEA+RCLSENSDIF 
C       Output
	   DIFRAD=(RCLSKY-Solrad)/((1-PERCLO)*RCLSKY)
	   DIFRAD=AMIN1(AMAX1(0.,DIFRAD),1.)
	  ENDIF
CGG	 SolDIF=AMIN1(Solrad,RCLDIF*(1-DIFRAD)+Solrad*DIFRAD)
	 SolSENSDIF=AMIN1(Solrad,RCLSENSDIF*(1-DIFRAD)+Solrad*DIFRAD)
	 SOLSKYDIF=(RCLSENSDIF*(1-DIFRAD)+Solrad*(DIFRAD-ASENS*VFSENS))/
     +(1-VFSENS)
CGGrrli 	 SOLDXY=SOLSKYDIF*(1.-VIEWFT(NX,NY))+
cggrrli     +ALBSMOO(NX,NY)/100.*SOLRAD*VIEWFT(NX,NY)
CGG	 SolBEA=Solrad-SolDIF
	 SolBEA=Solrad-SolSENSDIF
	RETURN
	END
C
C         ALBEDO IN PRESENCE OF DIFFUSE RADIATION and/or CLOUDS
C
C     DIFRAD= Input:Cloudiness/Diffused radiation [0,1]
C     T2112= Time [days] after 21st December
C     TEMCUM= Sum of Maximum Daily Temperatures above -10 ø C 
C             after last snowfall [øC]
C
C     ALBMEL Clear-sky albedo after snow metamorphosis
C     AM0   =Late sping albedo   (0.38)    
C     CAM  =Decay constant for aged snow albedo [1/øC]
C
       FUNCTION ALBDIF(DIFRAD,T2112,TEMCUM)
      COMMON /SNOW/A0,A1,TSPRIN,ASPRIN,TRNCF,AM0,CAM,GANESN,ALFAS,
     +BETAS,STRN,CJOMME,HSNLA1,FHI,Swi,FHIeKs,BEERSN,CMELTM,SPCFCT
C
C       Petzold-type formula
C
      VAR=0.499+0.0097*(10*DIFRAD)**3.
C
C       Aged-snow Clear-sky albedo Formula
C
      ALBMEL=AM0+(ALBCLS(T2112)-AM0)*EXP(-CAM*TEMCUM)
C
C       Aged snow albedo in presence of diffused radiation and/or Cloudiness
C
	  ALBDIF=AMIN1(ALBMEL*(1+VAR/100),1.)
	RETURN
	END
C
       FUNCTION Hg(dsoil1,julday,rain)
       real lambt,gamt,lambtp,gamtp
CGGG	DATA LAMBT,GAMT,A,omega,pi4/3.5,0.28,10.,0.017214,0.7854/
CGGG	SIBERIA
cggrrli  lambda (W/m/K) 1.676 e 0.603 vedi hpsuolo.xls dati Siberia
	DATA LAMBT,GAMT,A,omega,pi4/1.676,0.38,10.,0.017214,0.7854/
	DATA LAMBTP,GAMTP/0.603,0.63/
	if(rain.gt.0.) then
	 Hg=-sqrt(2.)*lambtp*gamtp*A*exp(-gamtp*dsoil1)
     + *sin(omega*julday-gamtp*dsoil1-pi4)
	else
	 Hg=-sqrt(2.)*lambt*gamt*A*exp(-gamt*dsoil1)
     + *sin(omega*julday-gamt*dsoil1-pi4)
	endif

	RETURN
	END
CGGG 
	FUNCTION tsc(heaWHC,Ta,Cv,zs,Tsoil,secdt,dw)
      COMMON /HEAT/CPw,Wden,ROhu,Cpa,RatiKh,CKa,Z0,VLh,RatiKe,Rgas,ZW,ZT
	Cvzsdt=Cv*zs/secdt
	Cterm=rohu*cpa*dw
	tsc=(Cvzsdt*Tsoil+Cterm*Ta+Heawhc)/(Cvzsdt+cterm)
c	tsc=((Cvzsdt-cterm)*Tsoil+Cterm*Ta+Heawhc)/Cvzsdt
c	write(27,'(5f11.3)')tsc, cvzsdt,cterm,heawhc,tsoil
	RETURN
	END
CRRMUL
CRRMUL ********************************************************
CRRMUL INTERCET.FOR
CRRMUL
$debug       
      subroutine INTERCET(T,RAIN,ETP,C,Q,EI,ET,ETA,TR,
     +rs,ra,deltag,gamma)
c parametri del metodo di Rutter-Calder: 
c p=free throughfall fraction:parte di pioggia che cade direttamente a terra,
c attraverso gli alberi
c k=drainage parameter [mm/h]
c b=drainage parameter [1/mm]
c S=maximum canopy storage capacity: capacità max di intercettazione 
c della vegetazione [mm]
c c=canopy storage [mm]      
      integer T
      real k,Kc
      DATA p,k,b,S/0.005,0.0168,1.76,0.288/
c      open(1,file='intercet.inp')
c      read(1,*)
c      open(2,file='resist.inp')
c      read(2,*)
c inizializza l'intercettazione reale e la pioggia utile
      if(t.eq.1) then
c      open(3,file='intercet.out')
c      write(3,'(a76)') 'hour     rain     ETP      c       Q      EI    
c     +   ET    ETA     TR  BALANCE'
c       read(1,'(i8,2(f8.3))') hour,rain,ETP
c       read(2,'(4(f8.3))') delta,gamma,rs,ra
       c=(1-p)*rain
       Q=(1-p)*rain
       Kc=k*(exp(b*c)-1)
       if(c.lt.S) then
        EI=min(ETP*c/S,c)
       else
        EI=min(ETP,c)
       endif
c       if(rain.eq.0) then
c        ET=ETP
c       else 
        ET=(ETP-EI)*(deltag+gamma)/(deltag+
     +   gamma*(1+rs/ra))
c       endif
       ETA=EI+ET
       TR=p*rain+Kc
       BALANCE=rain-c-EI-TR
c      write(3,'(i5,9(f8.4))') T,rain,ETP,c,Q,EI,ET,ETA,TR,BALANCE
      else       
c calcola l'intercettazione reale
c       read(1,'(i8,2(f8.3))',end=100) hour,rain,ETP
c       read(2,'(4(f8.3))') delta,gamma,rs,ra
       cold=c
       if(cold.lt.S) then
        EI=min(ETP*cold/S,cold+(1-p)*rain)
       else
        EI=min(ETP,cold+(1-p)*rain)
       endif
       Q=(1-p)*rain-EI
       Kc=k*(exp(b*cold)-1)
       c=min(cold+Q-Kc,S)
       TR=p*rain+Kc+max(cold+Q-Kc-S,0.)
c       if(rain.eq.0) then
c        ET=ETP
c       else 
        ET=(ETP-EI)*(deltag+gamma)/(deltag+
     +   gamma*(1+rs/ra))
c       endif
       ETA=EI+ET
       BALANCE=rain-(c-cold)-EI-TR
c      write(3,'(i5,9(f8.4))') T,rain,ETP,c,Q,EI,ET,ETA,TR,BALANCE
      endif
c100   continue      
      return
      end  
      
      
C      function rs(j_day,deficit)
c j_day=giorno giuliano, deficit=deficit di vapor saturo [kPa]      
C      deficit=2
C      pigreco=3.14      
C      if(deficit.lt.2.2) then
C       rs=74.5*(1-0.3*cos(2*pigreco*(j_day-222)/365))/(1-0.45*deficit)
C      else
C       rs=100000
C      endif
C      return
C      end
      
C      function ra(z,d,z0,uz)
Cc uz=velocità del vento [m/s] all'altezza z [m], hc=altezza della vegetazione [m],
c d=piano degli spostamenti nulli, z0=scabrezza
C      z=2
C      uz=2.5
C      hc=0.15
C      d=0.1*hc
C      z0=0.75*hc
C      ra=(alog((z-z0)/d))**2/(0.41**2*uz)
C      return
C      end


CRRMUL
CRRMUL ********************************************************
CRRMUL ETPSUB.FOR
CRRMUL
$debug       
C updated by Gitti and Grossi on 2006.01.24 to read zom e zov da z0 in COMMON/HEAT
********************************************************************      
C Subroutine ETPPM: calculates Potential Evapotranspiration by Penman-Monteith Equation.
      subroutine etppm(Ta,Ts,RHumi,press,Vw,energy,h0c,
     +ETP1,ETP2,ETPW,ra,deltag,gamma,hsnow)
      parameter(NNR=353,NNC=262)
      real h0c		 
      real log1,log2
      COMMON /HEAT/CPw,Wden,ROhu,Cpa,RatiKh,CKa,Z0,VLh,RatiKe,Rgas,ZW,ZT
c      COMMON /meteomo/SP10M,RH2M,TPR,RITOT,ALB,T2M,BLWB
*******************************************
C DETERMATION OF METEOROLOGICAL VARIABLES
cc
c VLhMJ = latent heat of vaporisation in MJ/kg
c gamma = psychrometric constant [kPa/°C]-4.2.28      
c ROhu = density of air [kg/m³]
c EAS = saturated vapour pressure [hPa=mbar]-4.2.2      
c EA = absolute vapour pressure [hPa=mbar]-pag.16 FAO n.24      
c deltag = gradient of saturated vapour pressure vs. temperature [kPa/°C]-4.2.3
cc      
      VLhMJ=VLh*1E-6
      gamma=0.0016286*press/10./VLhMJ
      if (h0c.le.0.5) then
	   EASkPa=EAS(Ts)/10.
C   ATTENZIONE!, usare ts anziche' ta	può comportare valori di ETPW eccessivamente
C elevati nel caso di debris perche' la differenza tra ts e ta è molto elevata 
c	   EASkPa=EAS(TA)/10.
	 else
	   EASkPa=EAS(Ta)/10.
	endif
      EAkPa=EA(EAS(ta),RHumi)/10.
      deltag=4098.*EASkPa/(237.7+Ta)**2
*************************************************      
c h0c = mean height of the crop [m]
c zom = roughness length for momentum [m]
c zov = roughness length for humidity [m]
c d = zero plane displacement for momentum and humidity [m]
c zW = height of sensor mesuring wind velocity [m]
c zt=hRH = height of sensor mesuring relative humidity [m]
cc
CGG   Former: CALCOLO NUMERO DI RICHARDSON E CORREZIONI PROFILI DI VENTO (deleted)
	if(hsnow.eq.-88) then
         zom=0.007
         zov=0.007
         d=0.
      elseif(hsnow.gt.h0c) then
cgg20060124 OLD          zom=0.003
cgg20060124 OLD         zov=0.003
cgg20060124  nel caso di neve sostituiamo z0m e z0v con il valore di z0 in disma.inp
cgg20060124  in modo da poterlo calibrare senza correggere il codice (24/01/06-GG,Gitti)
         zom=z0
         zov=z0
CMUL old version         d=0.67*h0c
CRRGGK
         d=amin1(zw-1.,hsnow)
CMUL RR:the problem is here  why ????
CMULold	   if(hsnow.gt.d) d=hsnow
         
	elseif(hsnow.gt.0.)	 then
cgg          zom=0.003
cgg         zov=0.003
cgg  nel caso di neve sostituiamo z0m e z0v con il valore di z0 in disma.inp
cgg  in modo da poterlo calibrare senza correggere il codice (24/01/06-GG,Gitti)
          zom=z0
         zov=z0
         d=0.67*h0c
CMUL RR:the problem is here. why ????
        if(hsnow.gt.d) d=amin1(zw-1.,hsnow)
      elseif(h0c.gt.0.)then
CRRGGK no snow, but vegetation
	   zom=0.123*h0c
         zov=0.1*zom
         d=0.67*h0c
    	else
CRRGK bare soil
         zom=0.01
         zov=0.01
         d=0.
	endif
Cgg!!	 CORREZIONE PER CNR-SIBERIA
CGG	zom=0.25
CGG	zov=0.25
CGG	zom=0.176
CGG	zov=0.176
CGG!!
CMUL
CMUL      write(*,*) 'zw=',zw,'d=',d,'zom=',zom
      log1=alog((zW-d)/zom)
	if(d.gt.zt) then
CMUL
CMUL	write(*,*) 'zw=',zw,'d=',d,'zom=',zom
      	log2=alog((zW-d)/zom)
	 else
      	if(d.ge.zt) then
		 log2=0.
	     else
CMUL
CMUL	write(*,*) 'zt=',zt,'d=',d,'zom=',zom
          log2=alog((zt-d)/zov)
      	endif
		if(log2.lt.0.) log2=0.
      endif
c	write(*,*)'ra,deltag,gamma=',ra,deltag,gamma      
*************************************************
C PENMAN-MONTEITH EQUATION 
      ETP1=1/VLh*(deltag*energy)/(deltag+gamma)
c    * (1+rs/ra))
c      ETP2=1/VLhMJ*(ROhu*Cpa/1000.*(EASkPa-EAkPa)/ra)/
	if(vw.eq.0.) then
	ETP2=0.
	else
      ra=(log1*log2)/(Vw*CKa**2)
CGG  CORREZIONE PROFILO VENTO CONDIZIONI STABILI/INSTABILI
Cg	if (rich.lt.0.)	then
Cggg	 condizioni instabili)
Cgg         ra=ra*(RICORF**(3./2.))
Cgg	else
Cggg	   stabili
Cgg         ra=ra*(ricorf**2.)
Cgg	endif

	  if(ra.gt.0.) then
           ETP2=1/VLh*(ROhu*Cpa*(EASkPa-EAkPa)/ra)/(deltag+gamma)
c     *(1+rs/ra))
	   else
	    etp2=0.
	  endif
	endif
      ETPW=ETP1+ETP2
C      write(*,*) 'ETP1,ETP2,ETPW=',ETP1,ETP2,ETPW
      return
      end
      

      FUNCTION EA(EAS,RHumi)
c EA=absolute vapour pressure [hPa=mbar]-pag.16 FAO n.24      
      EA=EAS*RHumi/100.
      RETURN
      END



CRRMUL
CRRMUL ********************************************************
CRRMUL END
CRRMUL


CRRMUL
CRRMUL ********************************************************
CRRMUL ETPSUB.FOR
CRRMUL
$debug       
C updated by Gitti and Grossi on 2006.01.24 to read zom e zov da z0 in COMMON/HEAT
********************************************************************      
C Subroutine ETPPM: calculates Potential Evapotranspiration by Penman-Monteith Equation.
C Revised by RR & GG to take into account snow and canopy height only at Wind speed station
C
C hsnst snow depth at station Hsnst=HSNPCK(NXSENS,NYSENS)
c h0cst canopy height at station H0c=h0c(NXSENS,NYSENS)

      subroutine etppm2(Ta,Ts,RHumi,press,Vw,energy,h0c,
     +ETP1,ETP2,ETPW,ra,deltag,gamma,hsnow,hsnst,h0cst)
      parameter(NNR=353,NNC=262)
      real h0c		 
      real log1,log2
      COMMON /HEAT/CPw,Wden,ROhu,Cpa,RatiKh,CKa,Z0,VLh,RatiKe,Rgas,ZW,ZT
*******************************************
C DETERMATION OF METEOROLOGICAL VARIABLES
cc
c VLhMJ = latent heat of vaporisation in MJ/kg
c gamma = psychrometric constant [kPa/°C]-4.2.28      
c ROhu = density of air [kg/m³]
c EAS = saturated vapour pressure [hPa=mbar]-4.2.2      
c EA = absolute vapour pressure [hPa=mbar]-pag.16 FAO n.24      
c deltag = gradient of saturated vapour pressure vs. temperature [kPa/°C]-4.2.3
cc      
      VLhMJ=VLh*1E-6
      gamma=0.0016286*press/10./VLhMJ
CRR20211226
Cc Vwresc= rescaled wind speed at sensor height in presence of snow at station. Station is assumed out of canopy as recommended
C	If(hsnst.gt.0) then
C		Vwresc=Vw*alog(Zw/z0)/alog((Zw-hsnst)/Z0)
C	else
C		Vwresc=Vw
C	endif
C_KARAKORAM_2.3 WIND SPEED NOT CORRECTED 
      Vwresc=Vw
      if (h0c.le.0.5) then
	   EASkPa=EAS(Ts)/10.
	 else
	   EASkPa=EAS(Ta)/10.
	endif
      EAkPa=EA(EAS(ta),RHumi)/10.
      deltag=4098.*EASkPa/(237.7+Ta)**2
*************************************************      
c h0c = mean height of the crop [m]
c zom = roughness length for momentum [m]
c zov = roughness length for humidity [m]
c d = zero plane displacement for momentum and humidity [m]
c zW = height of sensor mesuring wind velocity [m]
c zt=hRH = height of sensor mesuring relative humidity [m]
cc Aerodynamic coefficient at the sensor station
CC?? -88 
*************************************************
cc Aerodynamic coefficient at the site x,y

	if(hsnow.eq.-88) then
         zom=0.007
         zov=0.007
         d=0.
      elseif(hsnow.gt.h0c) then
CRR_20211226 canopy totally covered by snow
         zom=z0
         zov=z0
         d=0.
      elseif(h0c.gt.0.)then
CRR_20211226 canopy partially covered by snow

	   zom=0.123*(h0c-hsnow)
         zov=0.1*zom
         d=0.67*(h0c-hsnow)
    	else
CRR_20211226 bare soil

         zom=0.01
         zov=0.01
         d=0.
	endif
CRR_20211228 KARAKORAM 2.3
      dw=amin1(zw-1.,d)
      log1=alog((zW-dw)/zom)
      dt=amin1(zt-0.5,d)
      log2=amax1(alog((zt-dt)/zov),0.)
      	
*************************************************
C PENMAN-MONTEITH EQUATION 
      ETP1=1/VLh*(deltag*energy)/(deltag+gamma)
c    * (1+rs/ra))
c      ETP2=1/VLhMJ*(ROhu*Cpa/1000.*(EASkPa-EAkPa)/ra)/
	if(Vwresc.eq.0.) then
	ETP2=0.
	else

      ra=(log1*log2)/(Vwresc*CKa**2)

	  if(ra.gt.0.) then
           ETP2=1/VLh*(ROhu*Cpa*(EASkPa-EAkPa)/ra)/(deltag+gamma)
c     *(1+rs/ra))
	   else
	    etp2=0.
	  endif
	endif
      ETPW=ETP1+ETP2
      return
      end

