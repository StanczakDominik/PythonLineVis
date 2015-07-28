       PROGRAM MOLIBDEN

C     PROGRAM GLOWNY                                                  
C                                         
C     X-pinch
C
C     REVISED VERSION 12.01.2015
C      
C     THREE TEMPERATURE MODEL (STATIONARY)
C     WARUNEK BRZEGOWY NA GRADIENTY PREDKOSCI NIE NA CISNIENIE
c     staraya versya .AKAT. ,  4...+0.0...
c     with anomal resistivity, parabolic grid
C     vklutchena stabilizasta - lienynya vyazkost =0.03
c     vyklutchena vybros vnutrennych tochek
c     N.Yu. Orlov free paths data base
c     turn on limiters of Z and T           
C     KRAEVOY USLOVIE DLA SKOROSTEY (Z DAVLENIEM)
C     ISSKUSTVENNYA (TURBULENTANYA) LINEYNYIA VYASKOST'
c Eksi model 06.01.2015
c     w/o points remooving
       implicit double PRECISION(A-H,O-Z)
      COMMON PT(200,500,45),ndkt(200,500,18),NPDK(200,500),
     *TSS(200,500,9),NWB(200,500),WSPL(200,500,40),T
      common/nat/NBRZ(200,500,10),KPC(500),KKC(500),LL,LL1,          
     *LCZ,LKERN,KCZ,LCZ1,KCZ1,LWC,KWC,LWC2,KWC2,LMAX,KMAX
      COMMON/PARA/LWYR,KWYR,LIWYR,KIWYR,LOS,KOS
      COMMON/DANE/RO0,B0,B00,T0,ZET0,ETA0,AMJ,AME,SF,AIH,AI0,AI1,ZNJ,
     *DEI0,TF0,ALAM0,AB0,AKAA,AKII,TAUA0,TCH,OME,OMJ,AKE,AKBOL,AQ,CSW,
     *CV,AKAPE0,AKAPI0,TAU0,GAMA,PSLC,PR1,PR2,ETA1,CMIU0,DTMAX,SIG,AMAS 
      COMMON/IFPILM/sin1,cos1,pi,asek,pi2,PI6,dpi2
      COMMON/KROK/AJAW(200,500),AJAB(200,500),AJAWM,LPR,KPR
      COMMON/KRAK/TAUM,TAU77,spt,sptr
      COMMON/ELA/apt(200,500,16)
      COMMON/INITIAL/RPINCH,dpinch,adf,alfa,DEL,B01
      COMMON/ROBERT/AP(200,500)
      COMMON/JURA/dis(200,500,9),ROSMAX
      COMMON/STALA/PIW(0:73),AKBOL1,AS1,AS2,AS3,AS4,AS5,AS6
      COMMON/KAROL/SDRP,SLD,UW,UP,PS1,PS2
      COMMON/DESK/LTE,KTE,LTJ,KTJ
c      COMMON/PROM/STALA1,STALA2,STALA3
      COMMON/LIMIT/ALEP(200,500,2),OG(200,500)
      COMMON/GRZEGORZ/DPR(200,500,2)
      COMMON/ANOM/ANOM1,ANOM2,ANOM3,ANOM4,ANOM5,WAMX,LANOM,KANOM
      COMMON/VIC/CVA,CVB
      COMMON/DEG/TFF(200,500)
      COMMON/KRYT/LKRY,KKRY 
      COMMON/BASKO/EI1,EI2,EI3
      COMMON/PRZYR/ps(200,500,5)
      COMMON/ORLOV/DROS(11,11),DROP(11,11),droz(11,11),TEM(11),GES(11)
      COMMON/IRENA/LOP,KOP,LOK,KOK,LWQ,KWQ
      common/pomoc/ al(0:73),alL1(0:73),bl(0:73),bl1(0:73),cl(0:73)
      common/dif13/srza4,srza42,srz1c,srz1cc,akat1,akat01,akat02,kn5,kn7
      common/start/knwx,kntx,krok,kns
      common/srodki/as(8),nsek(8),NTX,NTX1,NBR,NBR1,rbr,rbr1 
       
          common/roznica/nges(200,500) 
       COMMON/WSTAWIONY/LSRO,KSRO,LSW,ISR(200,500)
         COMMON/WSTAWIONYBRZEG/IML1
        COMMON/NOWY/RSRO,ZSRO,LPOCZ,KPOCZ,IPOST
         common/artificial/pred(200,500,2)
         COMMON/TEST/PT1(200,500,5),LCZ5,KCZ5,LCZ6,KCZ6 
         COMMON/PETLA/LPOCZ1,KPOCZ1,IMIN,IML,NLJ,LPOCZ2,KPOCZ2
      COMMON/EWA/BRR(200,500,9),BRZ(200,500,9),BZZ(200,500,9),
     & AD(200,500,9),XA(200,500,6),AT(200,500)   
       COMMON/COLD/ACOL(6)

C
c      Common/dif41/am5,bm5,cm5,dt
        INTEGER knw,KNt,Lzap,Lmon
   
      DIMENSION POM(200,500,5),CZAS(35),den(400),PX(1200,500,3)
     ^ ,ZZ(200),NLI(400),rU(20500),zu(20500),FU(20500),PL(121,6)
     ^,nl(20500),nk(20500),DRE(11) 
       
      INTEGER IADJ(123000), IEND(20500)

      CHARACTER*15 DANE1,DANE2
      CHARACTER*75 FILSCR,FILDAT,FILDAT1,FILDAT2,SCRIPT,FILWEJ,zap,    
     ^ zap1,zap2,ZAP3,ZBIOR,FILDAT3,fild,ZAP4,FILDAT4
      CHARACTER*1 ZNAK

 101  FORMAT(///3X,1HL,3X,1HK,2X,1HR,9X,1HZ,
     *9X,2HRO,7X,4HU(R),5X,4HV(Z),6X,2HTE, 
     *7X,2HTJ,5X,4HBEFI,5X,3HZET//)   
 102  FORMAT(1X,I3,1X,I3,1X,1P9E9.2)   
 103  FORMAT(/2X,3H,I6,2X,3HDT=,1PE13.5,2X,2HL=,
     *I4,2X,2HK=,I4,2X,2HP=,1PE13.3,2X,3HRO=,1PE13.5,
     *2X,3HPS=,1PE13.3,2X,4HROS=,1PE13.5,1X,4HLSW=,I4)
 104  FORMAT(3X,3HL1=,I3,3X,3HK1=,I3,3X,
     *3HL2=,I3,3X,3HK2=,I3,3X,3HL3=,I3,
     *3X,3HK3=,I3,3X,3HL4=,I3,3X,3HK4=,I3,
     *3X,3HL5=,I3,3X,3HK5=,I3)
 105  FORMAT(//3X,5HCZAS=,1PE12.5,
     *10X,13HKROK CZASOWY=,1PE12.5//)
 107  FORMAT(3X,20I4)
 108  FORMAT(3X,12I4)
 109  FORMAT(/3X,4HLCZ=,I4,10X,5HLCZ1=,I4/)
 106  FORMAT(3X,I4,3X,1P5E11.3)
 110  FORMAT(/3X,4I4/)
 111  FORMAT(3X,I4,3X,I4,3X,1P7E11.3)
 114  FORMAT(1X,3HKN=,I6,1X,3HDT=,1PE12.3,1X,2HL=,
     *I4,1X,2HK=,I4,1X,2HP=,1PE11.3,1X,3HRO=,1PE12.5,1X,4HLSW=,I3)
 115  FORMAT(3X,I4,3X,I4,3X,1P1E11.3)
 116  FORMAT(' P_MAX(bar)',7X,I4,8X,I4,3X,1PE11.3
     *,2X,'(visc',1PE11.3,')',3X,I4)
 117  FORMAT(' Pe(bar)',2X,1PE11.3,6X,'Pi(bar)',3X,1PE11.3
     *,6X,('Prad(bar)',1PE11.3))
 118  FORMAT(' Pcold(bar)',2X,1PE11.3)
 119  FORMAT(' Pcold_max(bar)',2X,1PE11.3,3x,i4,3x,i4,3x,1PE11.3)
 216  FORMAT(4X,'Ar(MIN)= ',1PE11.4,6X,'Ar(MAX)= ',1PE11.4)
 217  FORMAT(4X,'Vr(MIN)= ',1PE11.4,6X,'Vr(MAX)= ',1PE11.4) 

      
 
       write(*,*) 'PODAJ ZBIOR WEJSCOWY   '
       read(*,*) dane1  
c       WRITE(*,*) 'PODAJ ZBIOR WYJSCIOWY DANE(XXX+1)  ' 
c       READ(*,*) DANE2
c	   OPEN(UNIT=7,FILE=DANE2,ACCESS='SEQUENTIAL',FORM=
c    *'FORMATTED')
	   
 673  continue



c      OPEN(3,FILE='e:/mhd/x/'//DANE1,STATUS='UNKNOWN',FORM=
      OPEN(3,FILE=DANE1,ACCESS='SEQUENTIAL',FORM=
     *'UNFORMATTED')

      
      TAU77=1D-12
     
      AKBOL=1.3805D-16
      GAMA=1.6666
      AME=9.1091D-28
      AQ=4.803D-10
      HKR=1.0546D-27
      rBohr=(hkr/aq)**2/ame
      CSW=3D10
      SIG=5.67D-5
      LICZ=0
      temax=2.5D-8
C      KNW=1
C      KNT=500

C     CZAS MAKSYMALNY TWWW

      TWWW=10D-8
      NZM=45
      NZMSS=9
      NZMWS=40
      TSWM=3D-9
      LSWM=48
      LWM=1
      NKZP=0
      KNN1=0
      KNBB=0
      DEB=1D7
      DEB1=1D7
      LL=35
      ll=131
      LL1=LL
      LSRO=0
      KSRO=0
      
      LPOCZ2=0
      KPOCZ2=0

      RPINCH=15D-4
      DPINCH=50D-4
      DEL=RPINCH/5.
      adf=3.
      Lmax=LL
      Kmax=350
       drca=2.*Dpinch/Dfloat(Kmax-1)
      DEL=Rpinch/5.
      U0=4D5
      pik=1.5D5
      PIK0=PIK
      Z0=U0/PIK
      

      
      knwx=knw
      kns=-1
c     SRZ1C - ROSSTYANIE NIZHE KOTOROGO UBIERAEEEM TOCHKI (DRCA -
C     PROSTANSTVIENNY SHAG (MALYI KOEFFITSENT MOZHET VESTI DO ROSTA
C     GRADIENTOV

C       srz1c=0.136*DRCA
        srz1c=0.14*DRCA
C        srz1c=0.15*DRCA   

C       srz1c=0.2*DRCA
      srz1cc=srz1c**2
c      0.3*DRCA
C      SRZC=30.*drca
C      SMAC=4.*drca
c      gama=gama
c      gama1=3.
      pslc=(gama+1.)/4.
C      pslc=3./4.
       pslc=0.03
      LICZ=0
      LCZS=0    
      LCZS1=0   
          
       lr=0
       kr=0

      SF=0.11
CC      B0=0.
      T0=1.16D4           
      T02=300.
      ZET0=0.1
      CUR=1D4       
      B0=CUR/(5.*RPINCH)
      B00=B0
      ROSMAX=RPINCH
      AMH=1.6726D-24
      AIH=2.18D-11
      AI1=1.666D-11
      ZNJ=42.
      AMAS=96.
      AMJ=AMAS*AMH
      RO0=aMj*2.5D21

      SIGEA=7D-16  
C      SIGIA=2.25E-12  
cc      TCH=3E-8
CC      AI0=4.806E-12/(1.-0.96/(ZNJ+1.)**0.257)
      CV=AKBOL/((GAMA-1.)*AMJ)
cc      ETA0=4.*SQRT(6.28*AME)*AQ**2/(3.*AKBOL**1.5)
      ETA0=AQ**2*DSQRT(AME)/(AKBOL**1.5*0.22) 
C      ETA0=4.*SQRT(ACOS(0.)*AME)*AQ**2/(3.*AKBOL**1.5)
      ETA1=16.*SIGEA*DSQRT(AKBOL*AME)/(3.*DSQRT(2.D0*3.1415D0)*AQ*AQ)
cc      AKAPE0=20.*(2./3.1415)**1.5*(AKBOL/AQ)**4/SQRT(AKBOL)
cc      AKAPE0=AKAPE0*0.43/SQRT(AME)
       AKAPE0=0.74*(AKBOL/AQ)**4/DSQRT(AKBOL)/Dsqrt(AME)
      AKAA=1./(SIGEA*dSQRT(AME))
      AKAPI0=3.28*(2./3.1415)**1.5*(AKBOL/AQ)**4/DSQRT(AKBOL)
      AKAPI0=AKAPI0/DSQRT(AMJ)
C      AKII=1./(SIGIA*SQRT(AMJ))
c      TAU0=3.*AMJ*AKBOL**1.5/(8.*SQRT(6.28*AME)*AQ**2)
c      TAU0=TAU0/AQ**2
      TAU0=0.11*(AMJ/aq)*(AKBOL/aq)/AQ**2*DSQRT(akbol/AME) 
      tau0= 0.75*(AMJ/AQ**4)/DSQRT(16.*DACOS(0D0)*ame)*AKBOL**1.5
C      TAU0=6.*ACOS(0.)/8.*AMJ*(HKR/AME)**2*(HKR/AQ**2)  
      TAUA0=3.*DSQRT(3.14D0)*AMJ/(8.*DSQRT(8.*AKBOL)*SIGEA)
      TAUA0=TAUA0/DSQRT(AME)
    
    
      PR1=SF*AKBOL**1.5/DSQRT(AME)
      PR2=AKBOL**1.5/AMJ**1.5
      CMIU0=3.*DSQRT(AMJ)*AKBOL**1.5/(4.*DSQRT(3.1415D0)*AQ**2)
      CMIU0=CMIU0/AQ**2
      DEI0=4.*3.1415*AQ*AQ/AKBOL
      DEI0=DEI0/AMJ
      TF0=(HKR/ame)*(HKR/akbol)*(3.*3.1415*3.1415)**0.6666/3.
      
c      ALAM0=(HKR/akbol)*(HKR/ame)/12.
      ALAM0=(HKR/akbol)*(HKR/ame)/2.
      AB0=AQ*AQ/(3.*AKBOL)
      OME=2.*AQ/(CSW*AMJ)
      OMJ=AQ/(CSW*AMJ)
      AKE=AKBOL*CSW/AQ
      SDRP=(GAMA-1.)*CV
      AKBOL1=AKBOL*1.16D4
      ABOR=(HKR/AQ)**2/AME
      EI1=ABOR**3
      EI1=AMJ/10.235/EI1
      EI2=0.41*AQ**2/ABOR/AKBOL
      EI3=(AQ/ABOR)**2/ABOR**2*
     ^ (3.*(2.*DACOS(0.D0))**2)**(2./3.)/5.

      ANOM1=AQ/(CSW*DSQRT(AME)*DSQRT(AMJ))
c      ANOM1=0.
	ANOM2=AME/AQ**2
cI      ANOM3=0.
      ANOM3=1d-2*DSQRT(12.566/AMJ)*AQ
      ANOM4=AKBOL/AMJ
      ANOM5=1./(12.566*AME*CSW**2)
      ANOM5=AQ/(AME*CSW)
corl      STALA1=0.4*AMAS/(1.1E7)**3
corl      STALA1=STALA1/(16.*2.*ACOS(0.)*ZNJ**4/3.**1.5)
corl      STALA2=(2.*AIH)**2*ZNJ**4
corl      STALA3=4./3**1.5
      SLD=DSQRT(AME/AMJ)*CMIU0
      AS1=1.2*AQ**2
      UW=AKBOL*DSQRT(AMJ/(8*DACOS(0.D0)))/AQ
      CVA=CV*AMJ
C      CVB=1./(GAMA-1.)
      PS1=AMAS*ZNJ
      PS2=4.*SIG/3./CSW
      PI=2.*DACOS(0.D0)
      dpi2=DSQRT(PI)/2.
      HKR=1.0546D-27
      PS2=8./15.*PI**5*AKBOL*(AKBOL/(CSW*HKR))**3/3.
corl      STALA2=1./(1.16E7)**1.493/67.02

            

      LCZ=LL1
      KCZ=KMAX
      LCZ1=LL1
      KCZ1=1
      DTMAX=1D-11
      DT=0.00001*DTMAX
      tau77=dt
      
C      RMCU=LL*DRCA*5
C      RMCU=30.*dpinch
c      RPINCH
      LWC=1
      KWC=1
      LWC2=1
      KWC2=KMAX
      LKRY=LWC
      KKRY=KWC
      LWYR=0
      KWYR=0
      LIWYR=0
      KIWYR=0
      IWYD=0
      IWYD1=0
      LMIN=1
      KMIN=1
      LRMAX=1
      KRMAX=kmax
      lb=1
      kb=1
      koniec=1
     
       
      PI=2.*DAcOS(0.D0)
     
      AseK=PI/4.
      pi2=4.*Dacos(0.D0)
c      PI6=PI/6. 
      PI6=PI/8.
c      PI6=PI/7.
      ADT=1.D0
 100   FORMAT(1P6E12.3)      
   

     

      

      DO 17 L=1,LMAX
      DO 17 K=1,KMAX
      wspl(l,k,33)=1.
        wspl(l,k,34)=1.
        wspl(l,k,35)=1.
        wspl(l,k,36)=1.
      AJAB(L,K)=1.
      NWB(L,K)=0
      NPDK(L,K)=0
      OG(L,K)=1.
      xa(l,k,3)=srz1c
       
      DO 9 N=1,NZM
   9  PT(L,K,N)=0.

      PT(L,K,3)=1.
      DO 10 N=1,10
  10  NBRZ(L,K,N)=0
      DO 11 N=1,18
  11  NDKT(L,K,N)=0
c      DO 7 N=1,2
      DO 8 N=1,NZMSS
   8  TSS(L,K,N)=0.         
      DO 18 N=1,NZMWS
  18  WSPL(L,K,N)=0.    
      APT(L,K,1)=0.
      APT(L,K,7)=1.
      WSPL(L,K,1)=1D-12
      WSPL(L,K,2)=1.
      WSPL(L,K,6)=1.
      WSPL(L,K,5)=1.
      WSPL(L,K,13)=1D10
      WSPL(L,K,21)=1D-3
      WSPL(L,K,22)=1D-3
      WSPL(L,K,4)=1.
      ALEP(L,K,1)=0.
      ALEP(L,K,2)=0.
      wspl(l,k,31)=5.*drca


      tss(l,k,7)=0.4
      apt(l,k,15)=0.     
      
  17  CONTINUE


      
     

      WRITE(*,*) ' START ODCZYTU Z DYSKU'
      READ(3) PT,NBRZ,NDKT,TSS,WSPL,T,DT,LCZ,KCZ
     *,LCZ1,KCZ1,LWC,KWC,LWC2,KWC2,NPDK,LL,KPC,KKC,APT
      CLOSE(3)
      WRITE(*,*) ' KONIEC ODCZYTU Z DYSKU DLA CZASU T=',T
      
      NPART=0.
      prmin=1d20
      

  67   continue
       ipor=0
      OPEN(36,FILE="dane2sasiedzi",ACCESS='SEQUENTIAL',FORM='FORMATTED')
 121  FORMAT(2I4,2PE13.4,2PE13.4)   
      DO 601 L=1,LL
      DO 601 K=KPC(L),KKC(L)
      if(pt(l,k,1).gt.100.) go to 601 
    

      TFf(l,k)=TF0*WSPL(L,K,13)**0.6666

      tss(l,k,6)=wspl(l,k,6)
      tss(l,k,7)=wspl(l,k,23)
      apt(l,k,6)=1.
c      PT(L,K,33)=1./3.
      apt(l,k,14)=0.
c      apt(l,k,15)=0.
      apt(l,k,9)=0.
C       ENDIF
      
      IF(PT(L,K,1).LT.100.) NPART=NPART+1
      pt(l,k,24)=0.
      pt(l,k,25)=0.
      pt(l,k,16)=0.
C      apt(l,k,13)=0
      apt(l,k,11)=1.
      AT(L,K)=1.5*AKBOL*WSPL(L,K,13)
      PT(L,K,35)=PT(L,K,13)
      pt(l,k,33)=pt(l,k,13)**0.25
      PT(L,K,43)=PT(L,K,11)
      PT(L,K,41)=PT(L,K,12)
      PT(L,K,39)=0.33333
      pt(l,k,36)=0.
      if(nbrz(l,k,1).ne.0.and.apt(l,k,2).lt.5.and.pt(l,k,1).gt.0.and.
     * abs(pt(l,k,2)).ne.dpinch) nbrz(l,k,1)=0

      wspl(l,k,31)=0.D0
        XA(L,K,1)=DSQRT(PT(L,K,11)**2+TFF(L,K)**2)
	

	DO 602 N=1,17,2

      N1=N+1
      KM=N1/2
      LI=NDKT(L,K,N)
      KI=NDKT(L,K,N1)

      RA=PT(LI,KI,1)
      ZA=PT(LI,KI,2)

      DIS(L,K,KM)=0.
      IF(PT(L,K,1).GT.100.) GO TO 602
      DIS(L,K,KM)=DSQRT((RA-PT(L,K,1))**2+(ZA-PT(L,K,2))**2)
c      write(*,*) l,k,'km=',km
c      if(l.eq.1.and.k.eq.1) then
	write(36,121) L,K, RA, ZA
c       write(*,*) 'km=',KM,'LI=',li,'KI=',ki
c      write(*,*) 'dra=',(ra-pt(l,k,1))/srz1c,'dza=',(za-pt(l,k,2))/srz1c
c      write(*,*) dis(l,k,km)/srz1c
c       write(*,*) pt(li,ki,13)**0.25/1.16d4
c       WRITE(*,*) PT(LI,KI,3)
c      endif
    
  602  CONTINUE
  601 continue
      
         KPOCZ=KPOCZ1
         LPOCZ2=LPOCZ
         KPOCZ2=KPOCZ
  

 
      STOP
      END
