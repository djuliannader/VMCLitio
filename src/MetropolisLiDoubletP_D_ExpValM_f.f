c    Program Metropolis.f
c
c    Time-stamp: <2019-02-25 12:22:34 daniel>
c
c    Lithium-like atoms (quartet state 1s2s2p)
c
c     VMC  Jan 10, 2023
c
      PROGRAM main
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c
      DOUBLE PRECISION, allocatable :: PV(:), deltaPV(:)

      character(10), allocatable :: text(:) 

      COMMON / XINPUT1 /  Bx, Z, N, NPV, ms
      COMMON / DELTA/  deltarho, deltaphi, deltaz, alambda

      COMMON / XMONITOR/ Ebest, Rbest, icall, icallmax, icallbest
c

c
c     EXECUTABLE STATEMENTS .....
c
c     
c     Variables defined to monitor the development of the calculation
      icall     = 0
      Ebest     = 100.0D0
      icallbest = icall
c
c     Reading Data ....
      read(5,*)
      read(5,*)Z    ! -> xB (in a.u);  nuclear charge
!      read(5,*)
!      read(5,*)(ms(ir), ir=1, 2)
      read(5,*)
      read(5,*)deltarho, deltaphi, deltaz  ! maximum change in rho, phi, z 
      read(5,*)
      read(5,*)N                ! Monte Carlo iterations
      read(5,*)
      read(5,*)ms               ! Monte Carlo subiterations
      read(5,*)
      read(5,*)icallmax         ! number of calls for wavefunction optimization
      read(5,*)
      read(5,*)alambda           ! lambda
      read(5,*)
      read(5,*)
      read(5,*)
      read(5,*)NPV      
      allocate(PV(NPV))
      allocate(text(NPV))
      allocate(deltaPV(NPV))
      Do i=1,NPV
         read(*,*)text(i), PV(i), deltaPV(i)
      Enddo
c
c
      Bx    = B               ! Magnetic Field in a.u. 
      sigma = sym             ! inherited definition of parity   
c     Writting data
      write(6,*)'******************************************'
      write(6,*)'Program MetroplisLiQuartetPV3_f.f   '
      write(6,*)' Ansatz A'
      write(6,*)'Nuclear Charge',Z
      write(6,*)'*****Monte Carlo parameters************'
      write(6,*)'Iterations = ',N
      write(6,*)'Subiterations = ',ms
      write(6,*)'Maximum change in rho = ',deltarho
      write(6,*)'Maximum change in phi = ',deltaphi
      write(6,*)'Maximum change in z = '  ,deltaz
      write(6,*)'**************************************'
      write(6,*)'Update parameter lambda =',alambda
      write(6,*)'******* Initial Variational Parameters**********'
      write(6,*)'name         initial value          fixed/released  '
      Do i=1,NPV
         write(*,*)text(i), PV(i), deltaPV(i)
      Enddo
      write(6,*)'***********************************'

c
c     Variational Monte Carlo performed here
      call cpu_time(start)
      call FCNE(PV,deltaPV)
      call cpu_time(finish)
c
c
c       

      write(6,*)'("Time = ',finish-start,' seconds.")'
      STOP
      END   ! END OF MAIN PROGRAM 
c
c







      
c     SUBUTINE FCN
c     Variational Energy for the 
c
c     XV: Array of Variational Parameters 
c
      SUBROUTINE FCNE(XV,DPV)
      IMPLICIT  DOUBLE PRECISION(A-H,O-Z)
c     .. COMMON BLOCK FROM INPUT DATA ...
      COMMON / XINPUT1 /  Bx, Z, N, NPV, ms
      COMMON / DELTA/  deltarho, deltaphi, deltaz, alambda
      COMMON / XMONITOR/ Ebest, Rbest, icall, icallmax, icallbest
      

      DOUBLE PRECISION   XV(NPV),DPV(NPV)
      INTEGER          NDIM
      DOUBLE PRECISION   PV(NPV)

c
            

c
c     .. COMMON BLOCK  FOR WAVE FUNCTION VARIATIONAL PARAMETERS ...

      
      write(6,*)'--------------------------------------------'
      write(6,*)'  **Optimization process begins**           '
      write(6,*)' ',icallmax,' iterations will be performed  '
      write(6,*)'--------------------------------------------'
c
      
      FVAL=0.0D0
      Do icall=1,icallmax
c
c
c     ------------------------------------------------------------------
c     ASSIGNING MINUIT VARIABLES ( XV VARIATIONAL PARAMETERS )
c
       c            =  XV(1)    !
       a            =  XV(2)    !
       alfa1        =  XV(3)    !
       alfa2        =  XV(4)
       alfa3        =  XV(5)
       alfa12       =  XV(6)    !
       c12          =  XV(7)    !
       d12          =  XV(8)    !
       alfa13       =  XV(9)    !
       alfa23       =  XV(10)   !

       
c    Saving variational parameters in a new array 
      PV=XV
       
c      Conditions for square normalizability
      ANOR1=alfa12+alfa13-alfa1
      ANOR2=alfa12+alfa23-alfa2
      ANOR3=alfa13+alfa23-alfa3
c      WRITE(*,*),'flag',ANOR1,ANOR2


      IF ((ANOR1 .GT. 0) .OR.  (ANOR2 .GT. 0) .OR. (ANOR3 .GT. 0)) THEN
      write(6,*)'-----------   Atention!!!!       --------------'
      write(6,*)'------------Nomalizability broken--------------'
      write(6,*)'-----------------------------------------------'
      write(6,*)'--- Set lambda smaller at the input file    ---'
      write(6,*)'---- or change initial set of parameters    ---'
      RETURN
      ENDIF
c
c
c
      Eprom=0.0D0
      ta=0.0D0
      tc=0.0D0
      error=0.0D0
      
c
c     Call Metropolis Routine
c     ************************************************************
      CALL METROPOLIS(PV,DPV,Eprom,ta,tc,error)
      
c     ************************************************************ 
c
c     saving proposal of variational parameters for next call
      XV=PV
c
c
c
c
      if(icall.EQ.1)then
      write(*,*)'------------------------------------------------'
      write(*,*)'Acceptation ratio ', ta
      write(*,*)'------------------------------------------------'      
      end if

c     *************************
 200  FVAL = Eprom
c     *************************
c
c     Updating E best
c     ----------------------------
c
c
c
c
c
c
c
c
      write(6,*)
      write(6,*)'****************************************************'
      write(6,*)'Li-like atoms doublet state 1s^22p (Compact Ansatz)  '
      write(6,'(A24,i8)')'Call :',icall
      write(6,*)'Evar   =',FVAL,' au'
      write(6,*)'Error  =',error,'au'
      write(6,*)'AC time=',tc,'iterations'
      write(6,*)'ratio Ms/AC t=',ms/tc,'>100 or tune Ms (input file)'
      write(6,*)
c
      write(6,111)' c       =',c
      write(6,111)' a       =',a
      write(6,111)' alfa1   =',alfa1
      write(6,111)' alfa2   =',alfa2
      write(6,111)' alfa3   =',alfa3
      write(6,111)' alfa12  =',alfa12
      write(6,111)' alfa12  =',c12
      write(6,111)' alfa12  =',d12
      write(6,111)' alfa13  =',alfa13
      write(6,111)' alfa23  =',alfa23
      write(6,*)'****************************************************'
 
      
c
 111   format(A20,2x,F12.8)
 112   format(A10,f12.8,A14,f12.8,A12,I4)
c
c
c
c
c
      write(6,*)
      enddo
c
c
c
      RETURN
      END
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c                             END FCN
c
c
      SUBROUTINE  METROPOLIS(PV,DPV,Eprom,ta,tc,error)
c     This Subroutine perform the Metropolis Algortim
      IMPLICIT  DOUBLE PRECISION(A-H,O-Z)
c
      COMMON / XINPUT1 /  Bx, Z, N, NPV, ms
      COMMON / DELTA/  deltarho, deltaphi, deltaz, alambda
c     Local Variables
      DOUBLE PRECISION coor1(3), coor2(3), coor1t(3), coor2t(3)
      DOUBLE PRECISION coor3(3), coor3t(3)
      
      DOUBLE PRECISION :: PV(NPV),DPV(NPV),UVP(NPV,3)
      DOUBLE PRECISION PR(4), PRt(4), PRH(4)
c
      double precision,dimension(2) :: psit,psi,P
      double precision,dimension(2) :: psi123,psi213,psi123t,psi213t 
      double precision,dimension(2) :: psi321,psi132,psi321t,psi132t
      double precision,dimension(2) :: psi312,psi231,psi312t,psi231t
      double precision,dimension(2) :: V123a,V213a,V321a,V132a,V312a
      double precision,dimension(2) :: V231a
      double precision,dimension(NPV) :: aK123,aK213,aK321,aK132,aK312
      double precision,dimension(NPV) :: aK231,psi123p,psi213p,psi321p
      double precision,dimension(NPV) :: psi132p,psi312p,psi231p      
      double precision,dimension(6) :: Vas,psis,psio
      double precision,dimension(6,NPV) :: aK,psip
      double precision :: Eprom, ta, tc, error
c
      

c 
c 
c      
c      
c! **************************Initial Configuration***************************      
      idum=1815  ! seed for random numbers
      icont=0
      E=0d0
      E22=0d0
      UVP(:,:)=0.0d0
      coor1(1)=0.1d0
      coor1(2)=0.2d0
      coor1(3)=0.3d0
      coor2(1)=0.2d0
      coor2(2)=0.3d0
      coor2(3)=0.4d0
      coor3(1)=0.3d0
      coor3(2)=0.4d0
      coor3(3)=0.5d0
      P(1)=0.0d0
      psi(1)=0d0
      psi(2)=0d0
      mb=N/ms

      Ebs=0d0
      E2bs=0d0      
c!*************   Metropolis   Initiates***********************************
      do i=1,mb
         Emsint=0.0D0
      do j=1,ms
c!      *********** Proposal for new coordinates*********************
         coor1t(1)=ABS(coor1(1)+deltarho*(2d0*ran2(idum)-1d0))
         coor1t(2)=ABS(coor1(2)+deltaphi*(2d0*ran2(idum)-1d0))
         coor1t(3)=coor1(3)+deltaz*(2d0*ran2(idum)-1d0)
         coor2t(1)=ABS(coor2(1)+deltarho*(2d0*ran2(idum)-1d0))
         coor2t(2)=ABS(coor2(2)+deltaphi*(2d0*ran2(idum)-1d0))
         coor2t(3)=coor2(3)+deltaz*(2d0*ran2(idum)-1d0)
         coor3t(1)=ABS(coor3(1)+deltarho*(2d0*ran2(idum)-1d0))
         coor3t(2)=ABS(coor3(2)+deltaphi*(2d0*ran2(idum)-1d0))
         coor3t(3)=coor3(3)+deltaz*(2d0*ran2(idum)-1d0)
c!      ************ Probability densities*************************
c
        ! -- Permutations
        call psiOnda(psi123t,coor1t,coor2t,coor3t,PV)
        call psiOnda(psi213t,coor2t,coor1t,coor3t,PV) 
        call psiOnda(psi321t,coor3t,coor2t,coor1t,PV)
        call psiOnda(psi132t,coor1t,coor3t,coor2t,PV) 
        call psiOnda(psi312t,coor3t,coor1t,coor2t,PV)
        call psiOnda(psi231t,coor2t,coor3t,coor1t,PV)
        
              c=PV(1)                 ! linear combination parameter

c       proyectors
        
        PRt(1)=(1/(2*sqrt(3.)))*(2*psi123t(1)+2*psi213t(1)-psi321t(1)-
     $       psi132t(1)-psi231t(1)-psi312t(1))

        PRt(2)=(1/2.0)*(psi321t(1)-
     $       psi132t(1)+psi231t(1)-psi312t(1))

        PRt(3)=(1/2.0)*(psi321t(1)-
     $       psi132t(1)-psi231t(1)+psi312t(1))

        PRt(4)=(1/(2*sqrt(3.)))*(2*psi123t(1)-2*psi213t(1)+psi321t(1)+
     $       psi132t(1)-psi231t(1)-psi312t(1))
        
        psit(1)=(PRt(1)+c*PRt(3))**2+(PRt(2)+c*PRt(4))**2 ! Density of probability

c     
        psit(2)=0.0D0 ! Imaginary part (proposal)
c
        ajact=coor1t(1)*coor2t(1)*coor3t(1)  !jacobian
c
        P(2)=ajact*( (psit(1))+(psit(2)))
        w=P(2)/P(1)  
c        !*********************************************************************
c       ! if (w > 0.00000000000000001) then
        IF (w .GE. 1D0) then    ! First conditional
           coor1=coor1t
           coor2=coor2t
           coor3=coor3t
           psi123=psi123t
           psi213=psi213t
           psi321=psi321t
           psi132=psi132t
           psi231=psi231t
           psi312=psi312t
           call VAsociado(V123a,aK123,psi123p,coor1t,coor2t,coor3t,PV)
           call VAsociado(V213a,aK213,psi213p,coor2t,coor1t,coor3t,PV) 
           call VAsociado(V321a,aK321,psi321p,coor3t,coor2t,coor1t,PV)
           call VAsociado(V132a,aK132,psi132p,coor1t,coor3t,coor2t,PV) 
           call VAsociado(V312a,aK312,psi312p,coor3t,coor1t,coor2t,PV)
           call VAsociado(V231a,aK231,psi231p,coor2t,coor3t,coor1t,PV)
           psi=psit
           P(1)=P(2)
           PR=PRt
           icont=icont+1
           call vOriginal(VV,coor1,coor2,coor3)
        ELSE 
           r=ran2(idum)
           IF (w .GT. r) THEN   ! Second conditional
           coor1=coor1t
           coor2=coor2t
           coor3=coor3t
           psi123=psi123t
           psi213=psi213t
           psi321=psi321t
           psi132=psi132t
           psi231=psi231t
           psi312=psi312t
           call VAsociado(V123a,aK123,psi123p,coor1t,coor2t,coor3t,PV)
           call VAsociado(V213a,aK213,psi213p,coor2t,coor1t,coor3t,PV) 
           call VAsociado(V321a,aK321,psi321p,coor3t,coor2t,coor1t,PV)
           call VAsociado(V132a,aK132,psi132p,coor1t,coor3t,coor2t,PV) 
           call VAsociado(V312a,aK312,psi312p,coor3t,coor1t,coor2t,PV)
           call VAsociado(V231a,aK231,psi231p,coor2t,coor3t,coor1t,PV)
           psi=psit
           P(1)=P(2)
           PR=PRt
           icont=icont+1
           call vOriginal(VV,coor1,coor2,coor3)
           END IF 
        END IF
c      !  end if
c!     ************************************************************
c     !     *************  Local energy computation********************
c      Denominator of the local energy
c      squared wave function after integrating spin
       den=(PR(1)+c*PR(3))**2+(PR(2)+c*PR(4))**2
  
c    projectors after applying Laplacian to the left
       PRH(1)=(1/(2*sqrt(3.)))*(2*psi123(1)*V123a(1)
     $       + 2*psi213(1)*V213a(1) - psi321(1)*V321a(1)
     $       - psi132(1)*V132a(1)  - psi231(1)*V231a(1)
     $       - psi312(1)*V312a(1))

       PRH(2)=(1/2.0)*(psi321(1)*V321a(1) - psi132(1)*V132a(1)
     $  + psi231(1)*V231a(1) - psi312(1)*V312a(1))

       PRH(3)=(1/2.0)*(psi321(1)*V321a(1) - psi132(1)*V132a(1)
     $  - psi231(1)*V231a(1) + psi312(1)*V312a(1))

       PRH(4)=(1/(2*sqrt(3.)))*(2*psi123(1)*V123a(1)
     $     - 2*psi213(1)*V213a(1) + psi321(1)*V321a(1)
     $     + psi132(1)*V132a(1)   - psi231(1)*V231a(1)
     $     - psi312(1)*V312a(1))

c    numerator for the local energy (potential term not included)
       anum = (PRH(1)+c*PRH(3))*(PR(1)+c*PR(3)) +
     $  (PRH(2)+c*PRH(4))*(PR(2)+c*PR(4))

c      Local energy
       Elocal=0.0*VV+1.0*anum/den
       E=E+Elocal
       E22=E22+Elocal*Elocal
       Emsint=Emsint+Elocal

c      Calling routine to calculate Dlocal and Hlocal
       Vas = (/ V123a(1),V213a(1),V321a(1),V132a(1),V312a(1),V213a(1) /)
       psio = (/ psi123(1),psi213(1),psi321(1),psi132(1),psi312(1),
     &       psi231(1) /)
       psip(1,:)=psi123p(:)
       psip(2,:)=psi213p(:)
       psip(3,:)=psi321p(:)
       psip(4,:)=psi132p(:)
       psip(5,:)=psi312p(:)
       psip(6,:)=psi213p(:)
       aK(1,:)=aK123(:)
       aK(2,:)=aK213(:)
       aK(3,:)=aK321(:)
       aK(4,:)=aK132(:)
       aK(5,:)=aK312(:)
       aK(6,:)=aK213(:)
       Call UPar(UVP,Vas,aK,psip,psio,den,PR,c,coor1,coor2,coor3,VV)

      end do
      ! sum over mb steps
      Ebs=Ebs+Emsint/ms
      E2bs=E2bs+(Emsint/ms)*(Emsint/ms)      
      end do
c      !*************** Mean energy and aceptation ratio
      Eprom=E/N 
      ta=1.0d0*icont/(1.0d0*N)

c     !*****************Error ***********************
      varb=(E2bs/mb)-(Ebs/mb)**2
      var=(E22/N)-(Eprom)**2
      tc=ms*(varb/var)
      error=Sqrt(tc*var/N)
c     *********************
c     updating variational parameters
      UVP=UVP/N
      call UpdatePar(PV,DPV,UVP,Eprom)
c      !******************************************************+
      return
      end subroutine
c      
c!     *****************************************************************************************
c!     *****************************************************************************************
c      !                  Routine which calculate trial function 
c!     ****************************************************************************************
c!     ****************************************************************************************
      subroutine  psiOnda(psi,coor1,coor2,coor3,varpars)
      IMPLICIT  DOUBLE PRECISION(A-H,O-Z)

      COMMON / XINPUT1 /  Bx, Z, N, NPV, ms 
      COMMON / DELTA/  deltarho, deltaphi, deltaz, alambda
  
c     Local variables
      double precision, dimension(3) :: coor1, coor2, coor3
      double precision, dimension(NPV) :: varpars
      double precision, dimension(2) :: psi, Va
c

c      
      pi=dacos(-1d0)
c            
c********** Cylindrical coordinates *************************************
      rho1=coor1(1)
      rho2=coor2(1)
      rho3=coor3(1)
      phi1=coor1(2)
      phi2=coor2(2)
      phi3=coor3(2)
      z1  =coor1(3)
      z2  =coor2(3)
      z3  =coor3(3)
      c=0.0D0                   ! neglecting parameter in the wave function
      a=varpars(2)
      al1=varpars(3)
      al2=varpars(4)
      al3=varpars(5)
      al12 =varpars(6)
      c12  =varpars(7)
      d12  =varpars(8)
      al13 =varpars(9)
      al23=varpars(10)
      
c     ************************ Correspondence to cartesian coordinates  ********************************************
      x1=rho1*COS(phi1)
      y1=rho1*SIN(phi1)
      x2=rho2*COS(phi2)
      y2=rho2*SIN(phi2)
      x3=rho3*COS(phi3)
      y3=rho3*SIN(phi3)
c!     ********************** Some definitions**************************************************
c

      
      x1mx2    = x1-x2
      y1my2    = y1-y2
      z1mz2    = z1-z2
      x1mx3    = x1-x3
      y1my3    = y1-y3
      z1mz3    = z1-z3
      x2mx3    = x2-x3
      y2my3    = y2-y3
      z2mz3    = z2-z3
c
c
c
      r12 =  SQRT( x1mx2*x1mx2 + y1my2*y1my2 + z1mz2*z1mz2)
      r1  =  SQRT( x1*x1 + y1*y1 + (z1)*(z1))
      r2  =  SQRT( x2*x2 + y2*y2 +  (z2)*(z2))
      r3  =  SQRT( x3*x3 + y3*y3 +  (z3)*(z3))
      r13 =  SQRT( x1mx3*x1mx3 + y1my3*y1my3 + z1mz3*z1mz3)
      r23 =  SQRT( x2mx3*x2mx3 + y2my3*y2my3 + z2mz3*z2mz3)

c     ************** Trial Function*************************************
c      Exponential factor
      exp_factor1 = exp( - al1*r1 - al2*r2 - al3*r3 
     $     + al12*r12*(1+c12*r12)/(1+d12*r12) + al13*r13 + al23*r23)
     
c     prefactor
      pre=z3*(1+a*r12)
c
c      wave function
       psi(1)=pre*exp_factor1   ! real part
       psi(2)=0.0D0             ! imaginary part

      return
      END SUBROUTINE


c      
c!     *****************************************************************************************
c!     *****************************************************************************************
c      !                  Routine which calculate the associated potential
c!     ****************************************************************************************
c!     ****************************************************************************************
      subroutine  VAsociado(Va,aK,psip,coor1,coor2,coor3,varpars)
      IMPLICIT  DOUBLE PRECISION(A-H,O-Z)

      COMMON / XINPUT1 /  Bx, Z, N, NPV, ms
      COMMON / DELTA/  deltarho, deltaphi, deltaz, alambda

c     Local variables 
      double precision, dimension(3) :: coor1, coor2, coor3
      double precision, dimension(NPV) :: varpars
      double precision, dimension(NPV) :: aK,psip
      double precision, dimension(2) ::  Va,psi
     
c
c      
      pi=dacos(-1d0)
c            
c********** Cylindrical coordinates *************************************
      rho1=coor1(1)
      rho2=coor2(1)
      rho3=coor3(1)
      phi1=coor1(2)
      phi2=coor2(2)
      phi3=coor3(2)
      z1  =coor1(3)
      z2  =coor2(3)
      z3  =coor3(3)
      c=0.0D0                   ! neglecting parameter in the wave function
      a=varpars(2)            ! neglecting parameter in the wave function
      al1=varpars(3)
      al2=varpars(4)
      al3=varpars(5)
      al12 =varpars(6)
      c12  =varpars(7)
      d12  =varpars(8)
      al13 =varpars(9)
      al23=varpars(10)

      
c     ************************ Correspondence to cartesian coordinates  ********************************************
      x1=rho1*COS(phi1)
      y1=rho1*SIN(phi1)
      x2=rho2*COS(phi2)
      y2=rho2*SIN(phi2)
      x3=rho3*COS(phi3)
      y3=rho3*SIN(phi3)
c!     ********************** Some definitions**************************************************
c

      
      x1mx2    = x1-x2
      y1my2    = y1-y2
      z1mz2    = z1-z2
      x1mx3    = x1-x3
      y1my3    = y1-y3
      z1mz3    = z1-z3
      x2mx3    = x2-x3
      y2my3    = y2-y3
      z2mz3    = z2-z3
c
c
c
      r12 =  SQRT( x1mx2*x1mx2 + y1my2*y1my2 + z1mz2*z1mz2)
      r1  =  SQRT( x1*x1 + y1*y1 + (z1)*(z1))
      r2  =  SQRT( x2*x2 + y2*y2 +  (z2)*(z2))
      r3  =  SQRT( x3*x3 + y3*y3 +  (z3)*(z3))
      r13 =  SQRT( x1mx3*x1mx3 + y1my3*y1my3 + z1mz3*z1mz3)
      r23 =  SQRT( x2mx3*x2mx3 + y2my3*y2my3 + z2mz3*z2mz3)
c     **********************More definitions ********************************   
c     ************** Trial Function*************************************
c      Exponential factor
      exp_factor1 = exp( - al1*r1 - al2*r2 - al3*r3 
     $     + al12*r12*(1+c12*r12)/(1+d12*r12) + al13*r13 + al23*r23)
     
c     prefactor
      pre=z3*(1.0+a*r12)
c
c      wave function
       psi(1)=pre*exp_factor1   ! real part
       psi(2)=0.0D0 ! imaginary part	 
      
c
c
c
c     ************  NOTE:     alfaij  are positives     **************
c

c     ---- first electron Laplacian ----------------------
       Va1=        (a*exp(-(al1*r1) - al2*r2 - al3*r3 + 
     -        (al12*r12*
     -           (1 + c12*r12))/
     -         (1 + d12*r12) + 
     -        al13*r13 + al23*r23)
     -      *(z1 - z2)**2*z3)/(r12**2)**1.5 - 
     -  (a*exp(-(al1*r1) - al2*r2 - al3*r3 + 
     -        (al12*r12*
     -           (1 + c12*r12))/
     -         (1 + d12*r12) + 
     -        al13*r13 + al23*r23)
     -      *z3)/r12 + 
     -  exp(-(al1*r1) - al2*r2 - al3*r3 + 
     -      (al12*r12*
     -         (1 + c12*r12))/
     -       (1 + d12*r12) + 
     -      al13*r13 + al23*r23)*
     -   z3*(1 + a*r12)*
     -   ((-2*al12*d12**2*(z1 - z2)**2*(1 + c12*r12))/
     -      (r12*(1 + d12*r12)**3)
     -       + (2*al12*c12*d12*(z1 - z2)**2)/
     -      (r12*(1 + d12*r12)**2)
     -       + (al12*d12*(1 + c12*r12))/
     -      (1 + d12*r12)**2 + 
     -     (al12*d12*(z1 - z2)**2*(1 + c12*r12))/
     -      ((r12**2)*(1 + d12*r12)**2) - 
     -     (al12*c12)/(1 + d12*r12) - 
     -     (al12*c12*(z1 - z2)**2)/
     -      ((r12**2)*(1 + d12*r12)) + 
     -     (al12*(z1 - z2)**2*(1 + c12*r12))/
     -      ((r12**2)**1.5*(1 + d12*r12))
     -      - (al12*(c12 + 1/r12))/
     -      (1 + d12*r12)) - 
     -  (a*exp(-(al1*r1) - al2*r2 - al3*r3 + 
     -        (al12*r12*
     -           (1 + c12*r12))/
     -         (1 + d12*r12) + 
     -        al13*r13 + al23*r23)
     -      *(z1 - z2)*z3*(-((al1*z1)/r1) + (al12*(z1 - z2)*
     -      (1 + c12*d12*(rho1**2 + rho2**2 + (z1 - z2)**2)
     - - 2*c12*d12*rho1*rho2*Cos(phi1 - phi2) + 
     -            2*c12*r12))/
     -        (r12*
     -          (1 + d12*r12)**2) + 
     -       (al13*(z1 - z3))/r13))/
     -   r12 + 
     -  (a*exp(-(al1*r1) - al2*r2 - al3*r3 + 
     -        (al12*r12*
     -           (1 + c12*r12))/
     -         (1 + d12*r12) + 
     -        al13*r13 + al23*r23)
     -      *(z1 - z2)*z3*(-((al2*z2)/r2) - (al12*(z1 - z2)*
     -      (1 + c12*d12*(rho1**2 + rho2**2 + (z1 - z2)**2)
     - - 2*c12*d12*rho1*rho2*Cos(phi1 - phi2) + 
     -            2*c12*r12))/
     -        (r12*
     -          (1 + d12*r12)**2) + 
     -       (al23*(z2 - z3))/r23))/
     -   r12 + 
     -  exp(-(al1*r1) - al2*r2 - al3*r3 + 
     -      (al12*r12*
     -         (1 + c12*r12))/
     -       (1 + d12*r12) + 
     -      al13*r13 + al23*r23)*
     -   z3*(1 + a*r12)*
     -      (-((al1*z1)/r1) + (al12*(z1 - z2)*(1 + c12*d12*(rho1**2
     - + rho2**2 + (z1 - z2)**2)-2*c12*d12*rho1*rho2*Cos(phi1 - phi2) + 
     -          2*c12*r12))/
     -      (r12*
     -        (1 + d12*r12)**2) + 
     -     (al13*(z1 - z3))/r13)*
     -      (-((al2*z2)/r2) - (al12*(z1 - z2)*(1 + c12*d12*(rho1**2
     - + rho2**2 + (z1-z2)**2) - 2*c12*d12*rho1*rho2*Cos(phi1 - phi2) + 
     -          2*c12*r12))/
     -      (r12*
     -        (1 + d12*r12)**2) + 
     -     (al23*(z2 - z3))/r23) + 
     -  Cos(phi1 - phi2)*(-((a*exp
     -           (-(al1*r1) - al2*r2 - al3*r3 + 
     -             (al12*r12*
     -                (1 + c12*r12))/
     -              (1 + d12*r12) + 
     -             al13*r13 + 
     -      al23*r23)*z3*(rho2 - rho1*Cos(phi1 - phi2))*(rho1
     -  - rho2*Cos(phi1 - phi2)))/
     -        (r12**2)**1.5) - 
     -     (a*exp(-(al1*r1) - al2*r2 - al3*r3 + 
     -           (al12*r12*
     -              (1 + c12*r12))/
     -            (1 + d12*r12) + 
     -           al13*r13 + 
     -           al23*r23)*z3*Cos(phi1 - phi2))/
     -      r12 + 
     -     exp(-(al1*r1) - al2*r2 - al3*r3 + 
     -         (al12*r12*
     -            (1 + c12*r12))/
     -          (1 + d12*r12) + 
     -         al13*r13 + al23*r23
     -         )*z3*(1 + a*r12)*
     -      ((2*al12*d12**2*(rho2 - rho1*Cos(phi1 - phi2))*(rho1
     - - rho2*Cos(phi1 - phi2))*
     -           (1 + c12*r12))/
     -         (r12*
     -           (1 + d12*r12)**3) + 
     -        (al12*d12*Cos(phi1 - phi2)*(1 + c12*r12))/
     -         (1 + d12*r12)**2 - 
     -      (al12*d12*(rho2 - rho1*Cos(phi1 - phi2))*(rho1
     -  - rho2*Cos(phi1 - phi2))*
     -           (1 + c12*r12))/
     -         ((r12**2)*(1 + d12*r12)**2)
     -         - (al12*c12*Cos(phi1 - phi2))/(1 + d12*r12) + 
     -      (al12*c12*(rho2 - rho1*Cos(phi1 - phi2))*(rho1
     -   - rho2*Cos(phi1 - phi2)))/
     -         ((r12**2)*(1 + d12*r12)) - 
     -      (al12*(rho2 - rho1*Cos(phi1 - phi2))*(rho1
     - - rho2*Cos(phi1 - phi2))*(1 + c12*r12))/
     -         ((r12**2)**1.5*
     -           (1 + d12*r12)) - 
     -        (al12*Cos(phi1 - phi2)*(1 + c12*r12))/
     -         (r12*(1 + d12*r12))
     -      + (2*al12*c12*d12*(rho2 - rho1*Cos(phi1 - phi2))*
     - (-rho1 + rho2*Cos(phi1 - phi2)))/
     -         (r12*
     -           (1 + d12*r12)**2)) + 
     -     (a*exp(-(al1*r1) - al2*r2 - al3*r3 + 
     -           (al12*r12*
     -              (1 + c12*r12))/
     -            (1 + d12*r12) + 
     -           al13*r13 + 
     -           al23*r23)*z3*(rho2 - rho1*Cos(phi1 - phi2))*
     -   (-((al1*rho1)/r1) - (al12*d12*(rho1 - rho2*Cos(phi1 - phi2))*
     -             (1 + c12*r12))/
     -           (1 + d12*r12)**2 + 
     -          (al12*(rho1 - rho2*Cos(phi1 - phi2))*(1 + c12*r12))/
     -           (r12*
     -             (1 + d12*r12)) + 
     -   (al12*c12*(rho1 - rho2*Cos(phi1 - phi2)))/(1 + d12*r12) + 
     -          (al13*(rho1 - rho3*Cos(phi1 - phi3)))/r13))/
     -      r12 + 
     -     (a*exp(-(al1*r1) - al2*r2 - al3*r3 + 
     -           (al12*r12*
     -              (1 + c12*r12))/
     -            (1 + d12*r12) + 
     -           al13*r13 + 
     -           al23*r23)*z3*(rho1 - rho2*Cos(phi1 - phi2))*
     -   (-((al2*rho2)/r2) - (al12*d12*(rho2 - rho1*Cos(phi1 - phi2))*
     -             (1 + c12*r12))/
     -           (1 + d12*r12)**2 + 
     -          (al12*(rho2 - rho1*Cos(phi1 - phi2))*(1 + c12*r12))/
     -           (r12*
     -             (1 + d12*r12)) + 
     -      (al12*c12*(rho2 - rho1*Cos(phi1 - phi2)))/(1 + d12*r12) + 
     -          (al23*(rho2 - rho3*Cos(phi2 - phi3)))/r23))/
     -      r12 + 
     -     exp(-(al1*r1) - al2*r2 - al3*r3 + 
     -         (al12*r12*
     -            (1 + c12*r12))/
     -          (1 + d12*r12) + 
     -         al13*r13 + al23*r23
     -         )*z3*(1 + a*r12)*
     -   (-((al1*rho1)/r1) - (al12*d12*(rho1 - rho2*Cos(phi1 - phi2))*
     -           (1 + c12*r12))/
     -         (1 + d12*r12)**2 + 
     -        (al12*(rho1 - rho2*Cos(phi1 - phi2))*(1 + c12*r12))/
     -         (r12*(1 + d12*r12))
     -      + (al12*c12*(rho1 - rho2*Cos(phi1 - phi2)))/(1 + d12*r12) + 
     -        (al13*(rho1 - rho3*Cos(phi1 - phi3)))/r13)*
     -   (-((al2*rho2)/r2) - (al12*d12*(rho2 - rho1*Cos(phi1 - phi2))*
     -           (1 + c12*r12))/
     -         (1 + d12*r12)**2 + 
     -        (al12*(rho2 - rho1*Cos(phi1 - phi2))*(1 + c12*r12))/
     -         (r12*(1 + d12*r12))
     -   + (al12*c12*(rho2 - rho1*Cos(phi1 - phi2)))/(1 + d12*r12) + 
     -        (al23*(rho2 - rho3*Cos(phi2 - phi3)))/r23)) + 
     -  exp(-(al1*r1) - al2*r2 - al3*r3 + 
     -      (al12*r12*
     -         (1 + c12*r12))/
     -       (1 + d12*r12) + 
     -      al13*r13 + 
     -      al23*r23)*z3*Sin(phi1 - phi2)*
     -      ((a*rho1*(rho1 - rho2*Cos(phi1 - phi2))*Sin(phi1 - phi2))
     -  /(r12**2)**1.5 - 
     -     (a*Sin(phi1 - phi2))/r12 + 
     -     (al12*(1 + a*r12)*
     -        (-d12 + (d12*rho1**2)/(r12**2) + 
     -          rho1**2/(r12**2)**1.5 + 
     -          (2*d12*rho1**2)/(r12**2) - 
     -          1/r12 + 
     -      c12*(-2 - d12**2*(rho2**2 + (z1 - z2)**2) - (3*d12*(rho2**2
     -  + (z1 - z2)**2 - 2*rho1*rho2*Cos(phi1 - phi2)))/
     -              r12) + 
     -          (rho1*rho2*Cos(phi1 - phi2)*(-1 - 3*d12*r12 + 
     -      c12*d12*(rho1**2 + rho2**2 + (z1 - z2)**2)*(-3 + d12*r12) - 
     -         2*c12*d12*rho1*rho2*Cos(phi1 - phi2)*(-3 + d12*r12)))/
     -           (r12**2)**1.5)*Sin(phi1 - phi2))/
     -      (1 + d12*r12)**3 - 
     - (a*rho1*(-((al1*rho1)/r1)-(al12*d12*(rho1-rho2*Cos(phi1 - phi2))*
     -             (1 + c12*r12))/
     -           (1 + d12*r12)**2 + 
     -      (al12*c12*(rho1 - rho2*Cos(phi1 - phi2)))/(1 + d12*r12) + 
     -          (al12*(rho1 - rho2*Cos(phi1 - phi2))*(1 + c12*r12))/
     -           (r12*
     -             (1 + d12*r12)) + 
     -   (al13*(rho1 - rho3*Cos(phi1 - phi3)))/r13)*Sin(phi1 - phi2))/
     -      r12 + 
     -      (a*(rho1 - rho2*Cos(phi1 - phi2))*(-((al12*rho1*(1
     -      + c12*d12*(rho1**2 + rho2**2 + (z1 - z2)**2)
     -  - 2*c12*d12*rho1*rho2*Cos(phi1 - phi2) + 
     -                 2*c12*r12)*Sin(phi1 - phi2))/
     -             (r12*
     -               (1 + d12*r12)**2)) + 
     -          (al23*rho3*Sin(phi2 - phi3))/r23))/
     -      r12 + 
     -     (1 + a*r12)*
     -    (-((al1*rho1)/r1) - (al12*d12*(rho1 - rho2*Cos(phi1 - phi2))*
     -           (1 + c12*r12))/
     -         (1 + d12*r12)**2 + 
     -        (al12*c12*(rho1 - rho2*Cos(phi1 - phi2)))/(1 + d12*r12) + 
     -        (al12*(rho1 - rho2*Cos(phi1 - phi2))*(1 + c12*r12))/
     -         (r12*
     -           (1 + d12*r12)) + 
     -        (al13*(rho1 - rho3*Cos(phi1 - phi3)))/r13)*
     -   (-((al12*rho1*(1 + c12*d12*(rho1**2 + rho2**2 + (z1 - z2)**2)
     -  - 2*c12*d12*rho1*rho2*Cos(phi1 - phi2) + 
     -               2*c12*r12)*Sin(phi1 - phi2))/
     -           (r12*
     -             (1 + d12*r12)**2)) + 
     -        (al23*rho3*Sin(phi2 - phi3))/r23)) + 
     -  exp(-(al1*r1) - al2*r2 - al3*r3 + 
     -      (al12*r12*
     -         (1 + c12*r12))/
     -       (1 + d12*r12) + 
     -      al13*r13 + 
     -      al23*r23)*z3*
     -      (-(Sin(phi1 - phi2)*((a*rho2*(-rho2 + rho1*Cos(phi1
     -  - phi2))*Sin(phi1 - phi2))/
     -           (r12**2)**1.5 + 
     -          (a*Sin(phi1 - phi2))/r12 + 
     -          (al12*(1 + a*r12)*
     -             (d12 - (d12*rho2**2)/(r12**2) - 
     -               rho2**2/(r12**2)**1.5 - 
     -               (2*d12*rho2**2)/(r12**2) + 
     -               1/r12 + 
     -      c12*(2 + d12**2*(rho1**2 + (z1 - z2)**2)
     -  +(3*d12*(rho1**2+(z1-z2)**2 - 2*rho1*rho2*Cos(phi1 - phi2)))/
     -                   r12) + 
     -     (rho1*rho2*Cos(phi1 - phi2)*(1 + 3*d12*r12 - 
     -     c12*d12*(rho1**2 + rho2**2 + (z1 - z2)**2)*(-3 + d12*r12) + 
     -        2*c12*d12*rho1*rho2*Cos(phi1 - phi2)*(-3 + d12*r12)))/
     -                (r12**2)**1.5)*Sin(phi1 - phi2))/
     -           (1 + d12*r12)**3 + 
     - (a*rho2*(-((al2*rho2)/r2)-(al12*d12*(rho2-rho1*Cos(phi1 - phi2))*
     -                  (1 + c12*r12))/
     -                (1 + d12*r12)**2 + 
     -      (al12*c12*(rho2 - rho1*Cos(phi1 - phi2)))/(1 + d12*r12) + 
     -     (al12*(rho2 - rho1*Cos(phi1 - phi2))*(1 + c12*r12))/
     -                (r12*
     -                  (1 + d12*r12)) + 
     -     (al23*(rho2 - rho3*Cos(phi2 - phi3)))/r23)*Sin(phi1 - phi2))/
     -           r12 + 
     -      (a*(rho2 - rho1*Cos(phi1 - phi2))*((al12*rho2*(1
     -      + c12*d12*(rho1**2 + rho2**2 + (z1 - z2)**2)
     -  - 2*c12*d12*rho1*rho2*Cos(phi1 - phi2) + 
     -                    2*c12*r12)*Sin(phi1 - phi2))/
     -                (r12*
     -                  (1 + d12*r12)**2) + 
     -               (al13*rho3*Sin(phi1 - phi3))/r13))/
     -           r12 + 
     -          (1 + a*r12)*
     -    (-((al2*rho2)/r2) - (al12*d12*(rho2 - rho1*Cos(phi1 - phi2))*
     -                (1 + c12*r12))/
     -              (1 + d12*r12)**2 + 
     -      (al12*c12*(rho2 - rho1*Cos(phi1 - phi2)))/(1 + d12*r12) + 
     -        (al12*(rho2 - rho1*Cos(phi1 - phi2))*(1 + c12*r12))/
     -              (r12*
     -                (1 + d12*r12)) + 
     -             (al23*(rho2 - rho3*Cos(phi2 - phi3)))/r23)*
     -      ((al12*rho2*(1 + c12*d12*(rho1**2 + rho2**2 + (z1 - z2)**2)
     -  - 2*c12*d12*rho1*rho2*Cos(phi1 - phi2) + 
     -                  2*c12*r12)*Sin(phi1 - phi2))/
     -              (r12*
     -                (1 + d12*r12)**2) + 
     -             (al13*rho3*Sin(phi1 - phi3))/r13))) + 
     -     Cos(phi1 - phi2)*(-((a*Cos(phi1 - phi2))/r12) + 
     -        (a*rho1*rho2*Sin(phi1 - phi2)**2)/(r12**2)**1.5 + 
     -        (al12*(1 + a*r12)*
     -  (2*c12*d12**2*rho1*rho2*Cos(phi1 - phi2)**2 + Cos(phi1 - phi2)*
     -              (-d12 - 1/r12 - 
     - c12*(2+d12**2*(rho1**2 + rho2**2 + (z1 - z2)**2) + 3*d12*r12)) + 
     -             rho1*rho2*((r12**2)**(-1.5) + 
     -                (3*d12)/(r12**2) + 
     -                c12*d12*(d12 + 3/r12))*Sin(phi1 - phi2)**2))/
     -         (1 + d12*r12)**3 - 
     -      (a*rho1*Sin(phi1 - phi2)*((al12*rho2*(1 + c12*d12*(rho1**2
     - + rho2**2+(z1 - z2)**2) - 2*c12*d12*rho1*rho2*Cos(phi1 - phi2) + 
     -                  2*c12*r12)*Sin(phi1 - phi2))/
     -              (r12*
     -                (1 + d12*r12)**2) + 
     -             (al13*rho3*Sin(phi1 - phi3))/r13))/
     -         r12 + 
     -      (a*rho2*Sin(phi1 - phi2)*(-((al12*rho1*(1 + c12*d12*(rho1**2
     -  + rho2**2+(z1 - z2)**2) - 2*c12*d12*rho1*rho2*Cos(phi1 - phi2) + 
     -                    2*c12*r12)*Sin(phi1 - phi2))/
     -                (r12*
     -                  (1 + d12*r12)**2)) + 
     -             (al23*rho3*Sin(phi2 - phi3))/r23))/
     -         r12 + 
     -        (1 + a*r12)*
     -      ((al12*rho2*(1 + c12*d12*(rho1**2 + rho2**2 + (z1 - z2)**2)
     - - 2*c12*d12*rho1*rho2*Cos(phi1 - phi2) + 
     -                2*c12*r12)*Sin(phi1 - phi2))/
     -            (r12*
     -              (1 + d12*r12)**2) + 
     -           (al13*rho3*Sin(phi1 - phi3))/r13)*
     -  (-((al12*rho1*(1 + c12*d12*(rho1**2 + rho2**2 + (z1 - z2)**2)
     -  - 2*c12*d12*rho1*rho2*Cos(phi1 - phi2) + 
     -                  2*c12*r12)*Sin(phi1 - phi2))/
     -              (r12*
     -                (1 + d12*r12)**2)) + 
     -           (al23*rho3*Sin(phi2 - phi3))/r23)))


       
c     ------------ second electron Laplacian------------


       Va2 =   -0.25*(exp(-(al1*r1) - al2*r2 - al3*r3 + 
     -        (al12*r12*
     -           (1 + c12*r12))/
     -         (1 + d12*r12) + 
     -        al13*r13 + al23*r23)
     -      *(4*al13*rho1*rho3*r1*z3*r3*r12*
     -        (1 + a*r12)*
     - (1+d12*r12)**2*(rho1**2+rho3**2 - 2*rho1*rho3*Cos(phi1 - phi3))*
     -        r23 - 
     -       4*a*rho1*rho3*r1*(z1 - z2)*r3*(1 + d12*r12)**2*
     -        (r13**2)**1.5*r23 + 
     -       4*al13*rho1*rho3*r1*z3*r3*r12*
     -        (1 + a*r12)*
     -        (1 + d12*r12)**2*Cos(phi1 - phi3)*
     -      (rho1*rho3 + (z1-z3)**2*Cos(phi1-phi3)
     -   -rho1*rho3*Cos(phi1-phi3)**2)*r23 - 
     -       4*rho1*rho3*r3*(1 + a*r12)*
     -        (r13**2)*
     -        (al13*r1*(z1 - z3)*r12*
     -           (1 + d12*r12)**2 - 
     -          al12*d12*r1*(z1 - z2)*r12*
     -           (1 + c12*r12)*r13
     -            + al12*c12*r1*(z1 - z2)*r12*
     -           (1 + d12*r12)*r13
     -            + al12*r1*(z1 - z2)*(1 + c12*r12)*
     -           (1 + d12*r12)*r13
     -            - al1*z1*r12*
     -           (1 + d12*r12)**2*
     -           r13)*r23 - 
     -      2*a*rho1*rho3*r1*z3*(rho1 - rho2*Cos(phi1 - phi2))*
     -   (1 + d12*r12)**2*
     -        Cos(phi1 - phi3)*(r13**2)*
     -        (2*al23*r3*r13*(rho3 - rho2*Cos(phi2 - phi3)) + 
     -          2*al13*r3*(rho3 - rho1*Cos(phi1 - phi3))*r23 - 
     -          2*al3*rho3*r13*r23
     -          ) - rho1*rho3*z3*(1 + a*r12)*Cos(phi1 - phi3)*
     -        r13*
     -        (2*al13*r1*r12*
     -           (1 + d12*r12)**2*(rho1 - rho3*Cos(phi1 - phi3)) - 
     -          2*al12*d12*r1*(rho1 - rho2*Cos(phi1 - phi2))*r12*
     -           (1 + c12*r12)*r13
     -            + 2*al12*c12*r1*(rho1 - rho2*Cos(phi1 - phi2))*r12*
     -           (1 + d12*r12)*r13
     -       + 2*al12*r1*(rho1 - rho2*Cos(phi1 - phi2))*(1 + c12*r12)*
     -           (1 + d12*r12)*r13
     -            - 2*al1*rho1*r12*
     -           (1 + d12*r12)**2*
     -           r13)*
     -        (2*al23*r3*r13*(rho3 - rho2*Cos(phi2 - phi3)) + 
     -          2*al13*r3*(rho3 - rho1*Cos(phi1 - phi3))*r23 - 
     -          2*al3*rho3*r13*r23
     -          ) + 4*a*rho1*rho3*r1*(z1 - z2)*z3*(1 + d12*r12)**2*
     -        (r13**2)*
     -        (-(al23*(-z2 + z3)*r3*r13) + 
     -          (al13*(z1 - z3)*r3 + al3*z3*r13)*
     -           r23) + 
     -       4*rho1*rho3*z3*(1 + a*r12)*
     -        r13*
     -        (al13*r1*(z1 - z3)*r12*
     -           (1 + d12*r12)**2 - 
     -          al12*d12*r1*(z1 - z2)*r12*
     -           (1 + c12*r12)*r13
     -            + al12*c12*r1*(z1 - z2)*r12*
     -           (1 + d12*r12)*r13
     -            + al12*r1*(z1 - z2)*(1 + c12*r12)*
     -           (1 + d12*r12)*r13
     -            - al1*z1*r12*
     -           (1 + d12*r12)**2*
     -           r13)*
     -        (-(al23*(-z2 + z3)*r3*r13) + 
     -          (al13*(z1 - z3)*r3 + al3*z3*r13)*
     -           r23) - 
     -       4*al13*rho1*rho3*r1*z3*r3*r12*
     -        (1 + a*r12)*
     - (1+d12*r12)**2*(-rho3**2-(z1-z3)**2+rho1*rho3*Cos(phi1 - phi3))*
     -        r23*Sin(phi1 - phi3)**2 + 
     -       4*a*rho1*rho3*r1*z3*r3*(rho1 - rho2*Cos(phi1 - phi2))*
     -        (1 + d12*r12)**2*(r13**2)*
     -        Sin(phi1 - phi3)*(al13*rho1*r23*Sin(phi1 - phi3) + 
     -          al23*rho2*r13*Sin(phi2 - phi3)) + 
     -       2*rho1*rho3*z3*r3*(1 + a*r12)*
     -        r13*
     -        (2*al13*r1*r12*
     -           (1 + d12*r12)**2*(rho1 - rho3*Cos(phi1 - phi3)) - 
     -          2*al12*d12*r1*(rho1 - rho2*Cos(phi1 - phi2))*r12*
     -           (1 + c12*r12)*r13
     -            + 2*al12*c12*r1*(rho1 - rho2*Cos(phi1 - phi2))*r12*
     -           (1 + d12*r12)*r13
     -      + 2*al12*r1*(rho1 - rho2*Cos(phi1 - phi2))*(1 + c12*r12)*
     -           (1 + d12*r12)*r13
     -            - 2*al1*rho1*r12*
     -           (1 + d12*r12)**2*
     -           r13)*Sin(phi1 - phi3)*
     -        (al13*rho1*r23*Sin(phi1 - phi3) + 
     -          al23*rho2*r13*Sin(phi2 - phi3)) + 
     -       2*r1*z3*(a*rho1*rho2*rho3*(1 + d12*r12)**2*
     -           (r13**2)*
     -           (2*al23*r3*r13*(rho3 - rho2*Cos(phi2 - phi3)) + 
     -             2*al13*r3*(rho3 - rho1*Cos(phi1 - phi3))*r23 - 
     -             2*al3*rho3*r13*
     -              r23)*Sin(phi1 - phi2)*Sin(phi1 - phi3) - 
     -          2*al13*rho1*rho3*r3*r12*
     -           (1 + a*r12)*
     -      (1 + d12*r12)**2*(-rho1**2 - (z1 - z3)**2
     -  + rho1*rho3*Cos(phi1 - phi3))*
     -           r23*Sin(phi1 - phi3)**2 - 
     -          rho1*rho3*(1 + a*r12)*
     -           r13*
     -           (2*al23*r3*r13*(rho3 - rho2*Cos(phi2 - phi3)) + 
     -             2*al13*r3*(rho3 - rho1*Cos(phi1 - phi3))*r23 - 
     -             2*al3*rho3*r13*
     -              r23)*Sin(phi1 - phi3)*
     -           (-(al12*rho2*r13*Sin(phi1 - phi2)) - 
     -             2*al12*c12*rho2*r12*
     -              r13*Sin(phi1 - phi2) - 
     -             al12*c12*d12*rho2*(r12**2)*
     -              r13*Sin(phi1 - phi2) - 
     -             al13*rho3*r12*Sin(phi1 - phi3) - 
     -             2*al13*d12*rho3*(r12**2)*Sin(phi1 - phi3) - 
     -             al13*d12**2*rho3*(r12**2)**1.5*Sin(phi1 - phi3)) + 
     -          2*al13*rho1*rho3*r3*r12*
     -           (1 + a*r12)*
     -           (1 + d12*r12)**2*Cos(phi1 - phi3)*
     -           r23*
     -    (Cos(phi1 - phi3)*(r13**2) - rho1*rho3*Sin(phi1 - phi3)**2) + 
     -          2*a*rho1*rho2*rho3*r3*(1 + d12*r12)**2*Cos(phi1 - phi3)*
     -           (r13**2)*Sin(phi1 - phi2)*
     -           (al13*rho1*r23*Sin(phi1 - phi3) + 
     -             al23*rho2*r13*Sin(phi2 - phi3)) - 
     -          2*rho1*rho3*r3*(1 + a*r12)*Cos(phi1 - phi3)*
     -           r13*
     -           (-(al12*rho2*r13*Sin(phi1 - phi2)) - 
     -             2*al12*c12*rho2*r12*
     -              r13*Sin(phi1 - phi2) - 
     -             al12*c12*d12*rho2*(r12**2)*
     -              r13*Sin(phi1 - phi2) - 
     -             al13*rho3*r12*Sin(phi1 - phi3) - 
     -             2*al13*d12*rho3*(r12**2)*Sin(phi1 - phi3) - 
     -             al13*d12**2*rho3*(r12**2)**1.5*Sin(phi1 - phi3))*
     -           (al13*rho1*r23*Sin(phi1 - phi3) + 
     -             al23*rho2*r13*Sin(phi2 - phi3)))))/
     -   (rho1*rho3*r1*r3*r12*
     -     (1 + d12*r12)**2*(r13**2)**1.5*
     -     r23)



            
   
c     ----------- third electron Laplacian -------

       Va3 =     (exp(-(al1*r1) - al2*r2 - al3*r3 + 
     -       (al12*r12*
     -          (1 + c12*r12))/
     -        (1 + d12*r12) + 
     -       al13*r13 + al23*r23)*
     -    (-4*al23*rho2*rho3*r2*z3*r3*r12*
     -       (1 + a*r12)*
     -       (1 + d12*r12)**2*r13*
     -      (rho2**2 + rho3**2 - 2*rho2*rho3*Cos(phi2 - phi3))
     -  - 4*a*rho2*rho3*(z1 - z2)*r2*r3*
     -       (1 + d12*r12)**2*r13*
     -       (r23**2)**1.5 + 
     -      4*rho2*rho3*r3*(1 + a*r12)*
     -       r13*(r23**2)*
     -       (al23*r2*(z2 - z3)*r12*
     -          (1 + d12*r12)**2 + 
     -         al12*d12*(z1 - z2)*r2*r12*
     -          (1 + c12*r12)*r23
     -          - al12*c12*(z1 - z2)*r2*r12*
     -          (1 + d12*r12)*r23
     -          - al12*(z1 - z2)*r2*(1 + c12*r12)*
     -          (1 + d12*r12)*r23
     -          - al2*z2*r12*
     -          (1 + d12*r12)**2*
     -          r23) + 
     -  2*a*rho2*rho3*r2*z3*(rho2-rho1*Cos(phi1-phi2))*(1+d12*r12)**2*
     -       Cos(phi2 - phi3)*(r23**2)*
     -       (2*al23*r3*r13*(rho3 - rho2*Cos(phi2 - phi3)) + 
     -         2*al13*r3*(rho3 - rho1*Cos(phi1 - phi3))*r23 - 
     -         2*al3*rho3*r13*r23)
     -        + rho2*rho3*z3*(1 + a*r12)*Cos(phi2 - phi3)*
     -       r23*
     -       (2*al23*r2*r12*
     -          (1 + d12*r12)**2*(rho2 - rho3*Cos(phi2 - phi3)) - 
     -         2*al12*d12*r2*(rho2 - rho1*Cos(phi1 - phi2))*r12*
     -          (1 + c12*r12)*r23
     -          + 2*al12*c12*r2*(rho2 - rho1*Cos(phi1 - phi2))*r12*
     -          (1 + d12*r12)*r23
     -      + 2*al12*r2*(rho2 - rho1*Cos(phi1 - phi2))*(1 + c12*r12)*
     -          (1 + d12*r12)*r23
     -          - 2*al2*rho2*r12*
     -          (1 + d12*r12)**2*
     -          r23)*
     -       (2*al23*r3*r13*(rho3 - rho2*Cos(phi2 - phi3)) + 
     -         2*al13*r3*(rho3 - rho1*Cos(phi1 - phi3))*r23 - 
     -         2*al3*rho3*r13*r23)
     -        + 4*a*rho2*rho3*(z1 - z2)*r2*z3*(1 + d12*r12)**2*
     -       (r23**2)*
     -       (-(al23*(-z2 + z3)*r3*r13) + 
     -         (al13*(z1 - z3)*r3 + al3*z3*r13)*
     -          r23) - 
     -      4*rho2*rho3*z3*(1 + a*r12)*
     -       r23*
     -       (al23*r2*(z2 - z3)*r12*
     -          (1 + d12*r12)**2 + 
     -         al12*d12*(z1 - z2)*r2*r12*
     -          (1 + c12*r12)*r23
     -          - al12*c12*(z1 - z2)*r2*r12*
     -          (1 + d12*r12)*r23
     -          - al12*(z1 - z2)*r2*(1 + c12*r12)*
     -          (1 + d12*r12)*r23
     -          - al2*z2*r12*
     -          (1 + d12*r12)**2*
     -          r23)*
     -       (-(al23*(-z2 + z3)*r3*r13) + 
     -         (al13*(z1 - z3)*r3 + al3*z3*r13)*
     -          r23) + 
     -      4*al23*rho2*rho3*r2*z3*r3*r12*
     -       (1 + a*r12)*
     -       (1 + d12*r12)**2*
     -       r13*(-rho3**2 - (z2 - z3)**2 + rho2*rho3*Cos(phi2 - phi3))*
     -      Sin(phi2 - phi3)**2-4*a*rho2*rho3*r2*z3*r3*(rho2-rho1
     -   *Cos(phi1 - phi2))*
     -       (1 + d12*r12)**2*(r23**2)*
     -       Sin(phi2 - phi3)*(al13*rho1*r23*Sin(phi1 - phi3) + 
     -         al23*rho2*r13*Sin(phi2 - phi3)) - 
     -      2*rho2*rho3*z3*r3*(1 + a*r12)*
     -       r23*
     -       (2*al23*r2*r12*
     -          (1 + d12*r12)**2*(rho2 - rho3*Cos(phi2 - phi3)) - 
     -         2*al12*d12*r2*(rho2 - rho1*Cos(phi1 - phi2))*r12*
     -          (1 + c12*r12)*r23
     -          + 2*al12*c12*r2*(rho2 - rho1*Cos(phi1 - phi2))*r12*
     -          (1 + d12*r12)*r23
     -      + 2*al12*r2*(rho2 - rho1*Cos(phi1 - phi2))*(1 + c12*r12)*
     -          (1 + d12*r12)*r23
     -          - 2*al2*rho2*r12*
     -          (1 + d12*r12)**2*
     -          r23)*Sin(phi2 - phi3)*
     -       (al13*rho1*r23*Sin(phi1 - phi3) + 
     -         al23*rho2*r13*Sin(phi2 - phi3)) + 
     -      2*al23*rho2*rho3*r2*z3*r3*r12*
     -       (1 + a*r12)*
     -       (1 + d12*r12)**2*
     -       r13*Cos(phi2 - phi3)*
     -      (-2*(z2 - z3)**2*Cos(phi2 - phi3)
     -  - 2*rho2*rho3*Sin(phi2 - phi3)**2) + 
     -      2*r2*z3*(a*rho1*rho2*rho3*(1 + d12*r12)**2*
     -          (r23**2)*
     -          (2*al23*r3*r13*(rho3 - rho2*Cos(phi2 - phi3)) + 
     -            2*al13*r3*(rho3 - rho1*Cos(phi1 - phi3))*r23 - 
     -            2*al3*rho3*r13*
     -             r23)*Sin(phi1 - phi2)*Sin(phi2 - phi3) + 
     -         2*al23*rho2*rho3*r3*r12*
     -          (1 + a*r12)*
     -          (1 + d12*r12)**2*
     -    r13*(-rho2**2 - (z2 - z3)**2 + rho2*rho3*Cos(phi2 - phi3))*
     -          Sin(phi2 - phi3)**2 - rho2*rho3*(1 + a*r12)*
     -          r23*
     -          (2*al23*r3*r13*(rho3 - rho2*Cos(phi2 - phi3)) + 
     -            2*al13*r3*(rho3 - rho1*Cos(phi1 - phi3))*r23 - 
     -            2*al3*rho3*r13*
     -             r23)*Sin(phi2 - phi3)*
     -          (-(al12*rho1*r23*Sin(phi1 - phi2)) - 
     -            2*al12*c12*rho1*r12*
     -             r23*Sin(phi1 - phi2) - 
     -            al12*c12*d12*rho1*(r12**2)*
     -             r23*Sin(phi1 - phi2) + 
     -            al23*rho3*r12*Sin(phi2 - phi3) + 
     -            2*al23*d12*rho3*(r12**2)*Sin(phi2 - phi3) + 
     -            al23*d12**2*rho3*(r12**2)**1.5*Sin(phi2 - phi3)) + 
     -         2*a*rho1*rho2*rho3*r3*(1 + d12*r12)**2*Cos(phi2 - phi3)*
     -          (r23**2)*Sin(phi1 - phi2)*
     -          (al13*rho1*r23*Sin(phi1 - phi3) + 
     -            al23*rho2*r13*Sin(phi2 - phi3)) - 
     -         2*rho2*rho3*r3*(1 + a*r12)*Cos(phi2 - phi3)*
     -          r23*
     -          (-(al12*rho1*r23*Sin(phi1 - phi2)) - 
     -            2*al12*c12*rho1*r12*
     -             r23*Sin(phi1 - phi2) - 
     -            al12*c12*d12*rho1*(r12**2)*
     -             r23*Sin(phi1 - phi2) + 
     -            al23*rho3*r12*Sin(phi2 - phi3) + 
     -            2*al23*d12*rho3*(r12**2)*Sin(phi2 - phi3) + 
     -            al23*d12**2*rho3*(r12**2)**1.5*Sin(phi2 - phi3))*
     -          (al13*rho1*r23*Sin(phi1 - phi3) + 
     -            al23*rho2*r13*Sin(phi2 - phi3)) - 
     -         2*al23*rho2*rho3*r3*r12*
     -          (1 + a*r12)*
     -          (1 + d12*r12)**2*
     -          r13*Cos(phi2 - phi3)*
     -  (Cos(phi2 - phi3)*(r23**2) - rho2*rho3*Sin(phi2 - phi3)**2))))/
     -  (4.*rho2*rho3*r2*r3*r12*
     -    (1 + d12*r12)**2*r13*
     -    (r23**2)**1.5)     
 
c    && ----- Total associated potential ----------
       
      Va(2) = 0.0D0         ! Imaginary part
      Va(1)=  (1.0*Va1+1.0*Va2+1.0*Va3)/psi(1)   ! Real part

      
      return
      END SUBROUTINE


      
c!     *****************************************************************************************
c!     *****************************************************************************************
c      !                  Routine to calculate the potential of the system
c!     ****************************************************************************************
c!     ****************************************************************************************
      SUBROUTINE  vOriginal(Vor,coor1,coor2,coor3)
      IMPLICIT  DOUBLE PRECISION(A-H,O-Z)
c
      COMMON / XINPUT1 /  Bx, Z, N, NPV, ms
      COMMON / DELTA/  deltarho, deltaphi, deltaz, alambda

c     Local variables     
      double precision, dimension(3) :: coor1,coor2,coor3
c
c     
      double precision :: pi=dacos(-1d0)
c      
c      
c!********** Cylindrical coordinates *************************************
      rho1=coor1(1)
      rho2=coor2(1)
      rho3=coor3(1)
      phi1=coor1(2)
      phi2=coor2(2)
      phi3=coor3(2)      
      z1  =coor1(3)
      z2  =coor2(3)
      z3  =coor3(3)
c!     ************************ Correspondence to cartesian coordinates  ********************************************
      x1=rho1*cos(phi1)
      y1=rho1*sin(phi1)
      x2=rho2*cos(phi2)
      y2=rho2*sin(phi2)
      x3=rho3*cos(phi3)
      y3=rho3*sin(phi3)
c!     ********************** Some definitions**************************************************
c
      x1mx2    = x1-x2
      y1my2    = y1-y2
      z1mz2    = z1-z2
      x1mx3    = x1-x3
      y1my3    = y1-y3
      z1mz3    = z1-z3
      x2mx3    = x2-x3
      y2my3    = y2-y3
      z2mz3    = z2-z3
cc
c
      r12 = sqrt( x1mx2*x1mx2 + y1my2*y1my2 + z1mz2*z1mz2)
      r1  =  sqrt( x1*x1 + y1*y1 + (z1)*(z1))
      r2  =  sqrt( x2*x2 + y2*y2 +  (z2)*(z2))
      r3  =  sqrt( x3*x3 + y3*y3 +  (z3)*(z3))
      r13 = sqrt( x1mx3*x1mx3 + y1my3*y1my3 + z1mz3*z1mz3)
      r23 = sqrt( x2mx3*x2mx3 + y2my3*y2my3 + z2mz3*z2mz3)
c!     *********************** Potential terms*********************
c
      Vor = 1.0*(1d0/r12 + 1d0/r23 + 1d0/r13) - Z/r1- Z/r2 -Z/r3
    
c           
c
      
      return
      END SUBROUTINE
c
c
c      !*****************************************************************
c      !***************   Routine to generates random numbers ***********
c      !*****************************************************************
      function ran2(idum)
      integer idum, im1, im2, imm1, ia1, ia2, iq1, iq2
      integer ir1,ir2,ntab,ndiv
      double precision ran2,am,eps,rnmx
      parameter (im1=2147483563,im2=2147483399,am=1./im1,imm1=im1-1,
     $     ia1=40014, ia2=40692, iq1=53668, iq2=52774, ir1=12211, 
     $     ir2=3791, ntab=32, ndiv=1+imm1/ntab,eps=1.2e-7,rnmx=1.-eps)
      integer idum2, j, k, iv(ntab), iy
      save iv, iy, idum2
      data idum2/123456789/, iv/ntab*0/, iy/0/
      if (idum.le.0) then
      idum=max(-idum,1)
      idum2 = idum
      do j=ntab+8,1,-1
      k = idum/iq1
      idum =ia1*(idum-k*iq1) - k*ir1
      if (idum.lt.0) idum = idum + im1
      if (j.le.ntab) iv(j) = idum
      enddo
      iy = iv(1)
      endif
      k = idum/iq1
      idum =ia1*(idum-k*iq1)-k*ir1
      if (idum.lt.0) idum = idum+im1
      k = idum2/iq2
      idum2 =ia2*(idum2-k*iq2)-k*ir2
      if (idum2.lt.0) idum2 = idum2+im2
      j = 1 + iy/ndiv
      iy = iv(j) - idum2
      iv(j) = idum
      if (iy.lt.1) iy = iy + imm1
      ran2 = min(am*iy,rnmx)
      return
      end


      Subroutine UpdatePar(PV,DPV,UVP,Eav)
      ! Subroutine update variational parameters
      IMPLICIT  DOUBLE PRECISION(A-H,O-Z)
c
      COMMON / XINPUT1 /  Bx, Z, N, NPV, ms
      COMMON / DELTA/  deltarho, deltaphi, deltaz, alambda

c     Local Variables
      double precision, dimension(NPV) :: PV, DPV
      double precision, dimension(NPV,3) :: UVP

      
      Do i=1, NPV
      PV(i) = PV(i) - (UVP(i,2)-UVP(i,3)-Eav*UVP(i,1))*DPV(i)*alambda
      ENDDO

c      print*, 'flag---->',UVP(1,:)
c      print*, 'flag---->',UVP(2,:)
c      print*, 'flag---->',-(UVP(1,2)+UVP(1,3)-Eav*UVP(1,1))
      
      end subroutine

      Subroutine UPar(UVP,Vas,aKs,psis,psi,den,PR,c,coor1,coor2,coor3,V)
c     Subroutine update variational parameters
      IMPLICIT  DOUBLE PRECISION(A-H,O-Z)
c     
      COMMON / XINPUT1 /  Bx, Z, N, NPV, ms
      COMMON / DELTA/  deltarho, deltaphi, deltaz, alambda
c     Local variables
      double precision, dimension(NPV,3) :: UVP
      double precision, dimension(6,NPV) :: aKs,psis
      double precision, dimension(6) :: Vas,psi
      double precision, dimension(4) :: PR,PRF,PRG1,PRG2
      double precision, dimension(3) :: coor1,coor2,coor3

c     variable psis receives values for Psi^prime_p and aKs for K_psi_p      

      Va123=Vas(1)
      Va213=Vas(2)
      Va321=Vas(3)
      Va132=Vas(4)
      Va312=Vas(5)
      Va231=Vas(6)

!     c=0.0D0

!     Updates for parameter c (linear combination of spin functions)
!  ***********************************************************************
      anumF = 1.0*(PR(3))*(PR(1)+c*PR(3))+
     $    1.0*(PR(4))*(PR(2)+c*PR(4))  
         
      UVP(1,1) = UVP(1,1) + anumF/den

      PRG1(3)=(1/2.0)*(psi(3)*(Va321+V) - psi(4)*(Va132+V)
     $  - psi(6)*(Va231+V) + psi(5)*(Va312+V))

       PRG1(4)=(1/(2*sqrt(3.)))*(2*psi(1)*(Va123+V)
     $     - 2*psi(2)*(Va213+V) + 1.0*psi(3)*(Va321+V)
     $     + 1.0*psi(4)*(Va132+V) - 1.0*psi(5)*(Va312+V)
     $      - 1.0*psi(6)*(Va231+V))

       anumG1 = 1.0*(PRG1(3))*(PR(1)+c*PR(3))+
     $  1.0*(PRG1(4))*(PR(2)+c*PR(4)) 
      
      UVP(1,2) = UVP(1,2) + anumG1/den
      
      UVP(1,3) = 0.0

c      print*,psi(:)


!     ***********************************************************************

!    Updates of quantities for parameters in the exponential wave function

      Do j=2,NPV

c    Formulas to calculate F for each iteration
c     projectors after applying operator to the left
         
       PRF(1)=(1/(2*sqrt(3.)))*(2*psis(1,j)
     $       + 2*psis(2,j) - psis(3,j)
     $       - psis(4,j)  - psis(5,j)
     $       - psis(6,j))

       PRF(2)=(1/2.0)*(psis(3,j) - psis(4,j)
     $  + psis(6,j) - psis(5,j))

       PRF(3)=(1/2.0)*(psis(3,j) - psis(4,j)
     $  - psis(6,j) + psis(5,j))

       PRF(4)=(1/(2*sqrt(3.)))*(2*psis(1,j)
     $     - 2*psis(2,j) + psis(3,j)
     $     + psis(4,j)   - psis(5,j)
     $      - psis(6,j))

      anumF = (PRF(1)+c*PRF(3))*(PR(1)+c*PR(3)) +
     $  (PRF(2)+c*PRF(4))*(PR(2)+c*PR(4))  
         
      UVP(j,1) = UVP(j,1) + anumF/den

c     ------------------------------------------------------------------
      
c    Formulas to calculate G part 1       
c    projectors after applying operator to the left

       PRG1(1)=(1/(2*sqrt(3.)))*(2*psis(1,j)*(Va123+V)
     $       + 2*psis(2,j)*(Va213+V) - psis(3,j)*(Va321+V)
     $       - psis(4,j)*(Va132+V)  - psis(5,j)*(Va312+V)
     $       - psis(6,j)*(Va231+V))

       PRG1(2)=(1/2.0)*(psis(3,j)*(Va321+V) - psis(4,j)*(Va132+V)
     $  + psis(6,j)*(Va231+V) - psis(5,j)*(Va312+V))

       PRG1(3)=(1/2.0)*(psis(3,j)*(Va321+V) - psis(4,j)*(Va132+V)
     $  - psis(6,j)*(Va231+V) + psis(5,j)*(Va312+V))

       PRG1(4)=(1/(2*sqrt(3.)))*(2*psis(1,j)*(Va123+V)
     $     - 2*psis(2,j)*(Va213+V) + psis(3,j)*(Va321+V)
     $     + psis(4,j)*(Va132+V) - psis(5,j)*(Va312+V)
     $      - psis(6,j)*(Va231+V))

      anumG1 = (PRG1(1)+c*PRG1(3))*(PR(1)+c*PR(3)) +
     $  (PRG1(2)+c*PRG1(4))*(PR(2)+c*PR(4)) 

      
      UVP(j,2) = UVP(j,2) + anumG1/den

c     ------------------------------------------------------------------
      
c    Formulas to calculate G part 2       
c    projectors after applying operator to the left

       PRG2(1)=(1/(2*sqrt(3.)))*(2*aKs(1,j)
     $       + 2*aKs(2,j) - aKs(3,j)
     $       - aKs(4,j) - aKs(5,j)
     $       - aKs(6,j))

       PRG2(2)=(1/2.0)*(aKs(3,j) - aKs(4,j)
     $  + aKs(6,j) - aKs(5,j))

       PRG2(3)=(1/2.0)*(aKs(3,j) - aKs(4,j)
     $  - aKs(6,j) + aKs(5,j))

       PRG2(4)=(1/(2*sqrt(3.)))*(2*aKs(1,j)
     $     - 2*aKs(2,j) + aKs(3,j)
     $     + aKs(4,j) - aKs(5,j)
     $      - aKs(6,j))

      anumG2 = (PRG2(1)+c*PRG2(3))*(PR(1)+c*PR(3)) +
     $  (PRG2(2)+c*PRG2(4))*(PR(2)+c*PR(4)) 

      
      UVP(j,3) = UVP(j,3) + anumG2/den

      enddo



      return
      end subroutine


