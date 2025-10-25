      subroutine gauss(n,a,b,x)
      implicit none
      integer n,i,j,k
      double precision a(1000,1000),s(100)
      double precision b(1000),x(1000)
      double precision smax,rmax,sum,pk,r,z
      integer p(1000)
      double precision max
c
C
      do 3 i=1,n
         p(i) = i
         smax = 0.0 
         do 2 j=1,n
            smax = max(smax,abs(a(i,j)))
 2       continue
         s(i) = smax
 3    continue
c
      do 7 k=1,n-1
         rmax = 0.0
         do 4 i=k,n
            r = abs(a(p(i),k))/s(p(i))
            if (r .gt. rmax) then
               j = i
               rmax = r
            endif
 4       continue
c
         pk = p(j)
         p(j) = p(k)
         p(k) = pk
c
         do 6 i=k+1,n      
            z = a(p(i),k)/a(p(k),k)       
            a(p(i),k) = z
            do 5 j=k+1,n    
               a(p(i),j) = a(p(i),j) - z*a(p(k),j)      
 5          continue
 6       continue  
 7    continue    
c
      do 9 k=1,n-1
         do 8 i=k+1,n
            b(p(i)) = b(p(i)) - a(p(i),k)*b(p(k))
 8       continue
 9    continue
      do 11 i=n,1,-1
         sum = b(p(i))
         do 10 j=i+1,n
            sum = sum - a(p(i),j)*x(j)
 10      continue
         x(i) = sum/a(p(i),i)
 11   continue
c
c
      return
      end 
C234567
      SUBROUTINE RES(T,Y,YPRIME,DELTA,IRES,RPAR,IPAR)
      double precision y(90),yprime(90),delta(90),rpar(1000)   
      double precision vel(93),pos(93),am(1000,1000),g(20),gp(22,22)
      double precision f(100),rl1(100),rl2(100),gi(20), rl(100)
      integer ires,ipar(1000)
      n=13
      do ii = 1, n
      pos(ii)=y(ii)  
      vel(ii)=y(ii+n)
      end do
      do i = 1,2
       rl(i) = y(2*n+i) 
      end do
      call slider(n,vel,pos,f,am,g,gi,gp,t,rl)       
      do i=1,n
       rl1(i)=0.0
      do j=1,2
        rl1(i)=rl1(i)+gp(j,i)*y(2*n+j)
      end do
      end do
      do ii = 1, n
      delta(ii)=yprime(ii)-y(n+ii)
      end do
      do i=1,n
      delta(i+n)=0.0
      do j=1,n
      delta(i+n)=delta(i+n)+am(i,j)*yprime(n+j)
      end do
      delta(i+n)=delta(i+n)-f(i)+rl1(i)
      end do
      do i=1,2
       delta(2*n+i)=0.0
      do j=1,n
        delta(2*n+i)=delta(2*n+i)+gp(i,j)*yprime(j+n)
      end do
        delta(2*n+i)=delta(2*n+i)+gi(i)
      end do
      return
      end
C
      program Slider_benchmark
      external res
      external jac
      INTEGER  NEQ, INFO(15), IDID, LRW, IWORK(100000), LIW, IPAR(1000)
      DOUBLE PRECISION x(93),pos(93),vel(93),rl(2),
     * T, Y(90), YPRIME(90), TOUT, RTOL, ATOL, RWORK(100000),
     * RPAR(1000),am(1000,1000),f(100),gp(22,22),g(22),gi(20)
      t=0.0
c     rwork(1) = 1.00
      tout=1D0
      info(1)=0
      info(2)=0
      info(3)=1
      info(4)=0
      info(5)=0
      info(6)=0
      info(7)=0
      info(8)=0
      info(9)=0
      info(10)=0
      info(11)=0
      info(12)=0
      info(13)=0
      info(14)=0
      info(15)=0
      lrw=100000
      liw=100000
      n = 13
      do ii = 1,n
       y(ii) = 0.0
       y(ii+n) = 0.0
      end do
      do ii = 1,2
      y(2*n + ii) = 0.0
      rl(ii) = 0.0
      end do
      y(3)=0.45D0
c     y(8) = 0.0001
      do ii = 1, n
      pos(ii)=y(ii)  
      vel(ii)=y(ii+n)
      end do
      call slider(n,vel,pos,f,am,g,gi,gp,t,rl)
      call gauss(n,am,f,x)
      write(*,*) (x(i),i=1,13)    
c     pause
c     do ii = 1, n
c      yprime(ii+n) = x(ii) 
c     end do
      neq=2*n+2
      rtol = 0.00000001D0 
      atol = 0.00000001D0
      open(unit=12,file='sl.dat')
      idid = 0
  30  call DDASSL (RES, NEQ, T, Y, YPRIME, TOUT, INFO, RTOL, ATOL,
     + IDID, RWORK, LRW, IWORK, LIW, RPAR, IPAR, JAC)
      write(*,*) idid
      write(12,24) t,(y(i),i=1,n)    
      write(*,24) t,(y(i),i=1,n)    
      if((idid .lt. 0) .or. (idid .gt. 3)) then 
       write(*,*) idid
       stop
       endif
      if (t .lt. tout ) go to 30
  24  format(E14.6,14G13.5)
      stop
      end

      subroutine jac
      stop 
      end
c234567
C
c$$$      subroutine slider(n,v,p,f,am,g,gi,gp,t)
c$$$      integer n
c$$$c 
c$$$c     implicit none
c$$$      integer i,j
c$$$      INTEGER    NP, NV, NL, NG, NU, LDG, IPAR(5), IFAIL
c$$$      LOGICAL    LFLAG(10)
c$$$      double precision T, P(N), V(N), U(100), RL(100), AM(100,100),
c$$$     &           GP(22,22), F(N), PDOT(100), UDOT(100), G(100),
c$$$     &           GI(100), FL(100,100), RPAR(100)
c$$$C
c$$$      np=n
c$$$      nv=n
c$$$      nl=2
c$$$      ng=2
c$$$      ldg=100
c$$$       do i=1,10
c$$$       lflag(i)=.true.
c$$$       end do
c$$$       lflag(8)=.false.
c$$$c
c$$$      rpar(1)=0.2
c$$$      rpar(2)=0.3 
c$$$      rpar(3)=9.81
c$$$c
c$$$      call FMBSE (NP, NV, NL, NG, NU, LDG, T, P, V, U,
c$$$     &                 RLAM, AM, GP, F, PDOT, UDOT, G, GI,
c$$$     &                 FL, LFLAG, RPAR, IPAR, IFAIL)
c$$$
c$$$      return
c$$$      end 
c$$$C234567
c$$$      subroutine slider(n,vel,pos,f,am,g,gp,t)
c$$$      double precision vel(n),pos(n),am(100,100),gp(22,22),
c$$$     *  t,g(22),f(n),m1,m2,m3,l1,l2,j1,j2,pi,gv,omega,tm
c$$$      pi=4.0*atan(1.0)
c$$$      m1=0.36D0
c$$$      m2=0.151104D0
c$$$      m3=0.075552D0
c$$$      l1=0.15D0
c$$$      l2=0.3D0
c$$$      j1=0.002727D0
c$$$      j2=0.0045339259D0
c$$$      gv=9.81
c$$$      omega=0.2
c$$$      tm=0.3
c$$$      am(1,1)=j1+m2*l1*l1
c$$$      am(2,1)=0.5*l1*l2*m2*cos(pos(1)-pos(2))
c$$$      am(1,2)=am(2,1)
c$$$      am(2,2)=j2
c$$$      am(3,3)=m3
c$$$      am(1,3)=0.0
c$$$      am(3,1)=0.0
c$$$      am(2,3)=0.0
c$$$      am(3,2)=0.0
c$$$      g(1)=l1*sin(pos(1))+l2*sin(pos(2))
c$$$      g(2)=pos(3)-l1*cos(pos(1))-l2*cos(pos(2))
c$$$      gp(1,1)=l1*cos(pos(1))
c$$$      gp(1,2)=l2*cos(pos(2))
c$$$      gp(1,3)=0.0
c$$$      gp(2,1)=l1*sin(pos(1))
c$$$      gp(2,2)=l2*sin(pos(2))
c$$$      gp(2,3)=1.0
c$$$      f(1)=-0.5*l1*gv*(m1+2.0*m2)*cos(pos(1))
c$$$     *  -0.5*l1*l2*m2*vel(2)*vel(2)*sin(pos(1)-pos(2))
c$$$      if (t .le. tm) then
c$$$         f(1)=f(1)+omega/tm*(1.0-cos(2*pi*t/tm))
c$$$      endif
c$$$      f(2)=-0.5*l2*gv*m2*cos(pos(2))
c$$$     *     +0.5*l1*l2*m2*vel(1)*vel(1)*sin(pos(1)-pos(2))
c$$$      f(3)=0.0
c$$$      return 
c$$$      end 
c234567
      subroutine slider(n,v,p,f,am,g,gi,gp,t,rlam)
      implicit none
      integer n
c
      integer i,j
      INTEGER    NP, NV, NL, NG, NU, LDG, IPAR(5), IFAIL
      LOGICAL    LFLAG(10)
      double precision T, P(N), V(N), U(100), RL(100), AM(1000,1000),
     &      GP(22,22), F(N), PDOT(100), UDOT(100), G(22),
     &      GI(20), FL(100,100), RPAR(9),EN, RLAM(100)
C
      np=n
      nv=n
      nl=2
      ng=2
      ldg=n+1
       do i=1,8
       lflag(i)=.true.
       end do
       lflag(10)=.true.
c
      ipar(1)=4
      rpar(1)=0.2
      rpar(2)=0.3
      rpar(3)=9.81

c     rpar(3)=0.0

      do i=4,ipar(1)+3
       rpar(i)=0.075
      end do
c      write(*,*) 'ipar',(ipar(i),i=1,5)
C      write(*,*) 'rpar',(rpar(i),i=1,ipar(1)+3)
c
      ldg = 22
      
      CALL FMBSE (NP, NV, NL, NG, NU, LDG, T, P, V, U,
     &         RLAM, AM, GP, F, PDOT, UDOT, G, GI,
     &         FL, LFLAG, RPAR, IPAR, IFAIL)

      return
      end

C
C     FMBSR:     Starrkoerpersystem
C     FMBSE:     Elastisches System
C     COMPIL ... FE Ansatz Balken fuer Verbindungsstab
C
      SUBROUTINE FMBSR (NP, NV, NL, NG, NU, LDG, T, P, V, U,
     &                 RLAM, AM, GP, F, PDOT, UDOT, G, GI,
     &                 FL, LFLAG, RPAR, IPAR, IFAIL) 
      INTEGER    NP, NV, NL, NG, NU, LDG, IPAR(1), IFAIL
      LOGICAL    LFLAG(10)
      REAL*8     T, P(NP), V(NV), U(NU), RLAM(NL), AM(100,100),
     &           GP(22,22), F(NV), PDOT(NP), UDOT(NU), G(NG),  
     &           GI(NL), FL(NV,NL), RPAR(3)
C 
C     Slider crank - Rigid system with sliding block
C     ----------------------------------------------
C  ** written by Bernd Simeon, TH Darmstadt, 05/22/95 **
C
C     FMBSR defines the rigid slider crank mechanism described in
C
C     Simeon, B.: Modelling a Flexible Slider Crank Mechanism
C     by a Mixed System of DAEs and PDEs. Technical Report
C     TH Darmstadt, to appear.
C
C     The interface (with certain extensions) is according to
C     
C     Lubich Ch., Nowak U., Engstler Ch.: MEXX - Numerical Software
C     for the Integration of Constraint Mechanical Systems.
C     Konrad-Zuse-Zentrum fuer Informationstechnik Berlin,
C     Technical Report SC 92-12.
C
C     PARAMETERS ON ENTRY:
C
C       NP       This integer defines the dimension of P.
C       NV       This integer defines the dimension of V.
C       NL=NG    This integer defines the dimension of RLAM.
C       LDG      This integer defines the leading dim. of GP.
C       T        This real variable contains the current value of the 
C                independent variable (time).
C       P(NP)    This array contains the current values of the
C                position  variables.
C                ( 3 rigid motion coordinates 
C                    phi1:  crank angle
C                    phi2:  angle of reference frame connecting rod
C                    x3  :  x-coordinate of sliding block        )
C       V(NV)    This array contains the current values of the
C                velocity variables.
C       RLAM(NL) This array contains the current values of the
C                Lagrange multipliers.
C       LFLAG(10) This logical array determines which part of the 
C                equations of motion is to be evaluated by FMBS.
C
C                if (LFLAG(1))  M(t,p) in AM
C                   (LFLAG(2))  G(t,p) in GP
C                   (LFLAG(4))  f(t,p,v,lambda,u) in F
C                   (LFLAG(5))  T(t,p)*v in PDOT
C                   (LFLAG(7))  g(t,p) in G
C                   (LFLAG(8))  gI(t,p) in GI
C                   (LFLAG(10)) z(t,p,v) in GI      
C
C     ON RETURN:
C
C       AM(3,3)  This real matrix contains the mass matrix M(t,p).
C       GP(LDG,3) This real matrix contains the constraint matrix G(t,p).
C       F(3)     This real array contains the forces f(t,p,v,lambda,u).
C       PDOT(3)  This real array contains the derivative p'=T(t,p)*v.
C       G(2)     This real array contains the constraints g(t,p).
C       GI(2)    This real array contains the constraint derivative gI(t,p).
C       IFAIL    This integer flag indicates a stop condition 
C                iff IFAIL .NE. 0.  (Unused).
C
C     SYSTEM PARAMETERS:
C
C       RPAR(1)  omega for drive torque, RPAR(2) t1 for drive torque
C       RPAR(3)  gravity constant
C
        REAL*8   GRAV, OMEGA, T1, J1, J2, L1, L2, M1, M2, M3, PI
C
        PARAMETER( M1 = 0.36D0,     M2 = 0.151104D0,
     *             M3 = 0.075552D0, 
     *             L1 = 0.15D0,     L2 = 0.30D0,
     *             J1 = 0.002727D0, J2 = 0.0045339259D0,
     *             PI = 3.1415927D0  )
C
        REAL*8   COSP1, COSP2, SINP1, SINP2, COSP12, SINP12
C
C_____________End of declaration part_________________________
C
      IFAIL = 0
C
      OMEGA = RPAR(1)
      T1    = RPAR(2)
      GRAV  = RPAR(3)
C
      COSP1  = COS(P(1))
      COSP2  = COS(P(2))
      SINP1  = SIN(P(1))
      SINP2  = SIN(P(2))
      COSP12 = COS(P(1)-P(2))
      SINP12 = SIN(P(1)-P(2))
C
      IF (LFLAG(1)) THEN 
          AM(1,1) = J1 + M2*L1*L1
          AM(2,1) = .5*L1*L2*M2*COSP12
          AM(1,2) = .5*L1*L2*M2*COSP12 
          AM(2,2) = J2
          AM(1,3) = 0.D0
          AM(2,3) = 0.D0
          AM(3,1) = 0.D0
          AM(3,2) = 0.D0
          AM(3,3) = M3
      END IF
      IF (LFLAG(2) .OR. LFLAG(3)) THEN 
            GP(1,1) = L1*COSP1
            GP(1,2) = L2*COSP2
            GP(1,3) = 0.D0
            GP(2,1) = L1*SINP1
            GP(2,2) = L2*SINP2
            GP(2,3) = 1.D0
      END IF
      IF (LFLAG(4)) THEN 
          F(1) = -.5*L1*GRAV*(M1+2.*M2)*COSP1 
     &           -.5*L1*L2*M2*V(2)*V(2)*SINP12
          IF (T .LE. T1) THEN
             F(1) =  F(1)+OMEGA/T1*(1.D0-COS(2.D0*PI*T/T1))
          END IF
          F(2) = -.5*L2*GRAV*M2*COSP2
     &           +.5*L1*L2*M2*V(1)*V(1)*SINP12
          F(3) = 0.D0
      END IF
      IF (LFLAG(5)) THEN 
          PDOT(1) = V(1)
          PDOT(2) = V(2)
          PDOT(3) = V(3)
      END IF
      IF (LFLAG(7)) THEN
            G(1) = L1*SINP1 + L2*SINP2
            G(2) = P(3) - L1*COSP1 - L2*COSP2
      END IF
      IF (LFLAG(8)) THEN
            GI(1) = 0.0D0
            GI(2) = 0.0D0
      END IF
      IF (LFLAG(10)) THEN 
            GI(1) = -L1*SINP1*V(1)*V(1) - L2*SINP2*V(2)*V(2)
            GI(2) =  L1*COSP1*V(1)*V(1) + L2*COSP2*V(2)*V(2)
      END IF
C_____________________________________________________________
C
      RETURN
      END

      SUBROUTINE FMBSE (NP, NV, NL, NG, NU, LDG, T, P, V, U,
     &                 RLAM, AM, GP, F, PDOT, UDOT, G, GI,
     &                 FL, LFLAG, RPAR, IPAR, IFAIL) 
      INTEGER    NP, NV, NL, NG, NU, LDG, IPAR(5), IFAIL
      LOGICAL    LFLAG(10)
      REAL*8     T, P(NP), V(NV), U(NU), RLAM(NL), AM(1000,1000),
     &           GP(22,22), F(NV), PDOT(NP), UDOT(NU), G(22),  
     &           GI(20), FL(NV,NL), RPAR(*)
C 
C     Slider crank - flexible model with sliding block
C     ------------------------------------------------
C  ** written by Bernd Simeon, TH Darmstadt, 05/22/95 **
C
C     FMBSE defines the flexible slider crank mechanism described in
C
C     Simeon, B.: Modelling a Flexible Slider Crank Mechanism
C     by a Mixed System of DAEs and PDEs. Technical Report
C     TH Darmstadt, to appear.
C
C     The interface (with certain extensions) is according to
C     
C     Lubich Ch., Nowak U., Engstler Ch.: MEXX - Numerical Software
C     for the Integration of Constraint Mechanical Systems.
C     Konrad-Zuse-Zentrum fuer Informationstechnik Berlin,
C     Technical Report SC 92-12.
C
C     PARAMETERS ON ENTRY:
C
C       NP       This integer defines the dimension of P.
C       NV       This integer defines the dimension of V.
C       NL=NG    This integer defines the dimension of RLAM.
C       LDG      This integer defines the leading dim. of GP.
C       T        This real variable contains the current value of the 
C                independent variable (time).
C       P(NP)    This array contains the current values of the
C                position  variables 
C                ( 3 rigid motion coordinates 
C                    phi1:  crank angle
C                    phi2:  angle of reference frame connecting rod
C                    x3  :  x-coordinate of sliding block
C                  plus NP-3 flexible motion coordinates)
C       V(NV)    This array contains the current values of the
C                velocity variables.
C       RLAM(NL) This array contains the current values of the
C                Lagrange multipliers.
C       LFLAG(10) This logical array determines which part of the 
C                equations of motion is to be evaluated by FMBS.
C
C                if (LFLAG(1))  M(t,p) in AM
C                   (LFLAG(2))  G(t,p) in GP
C                   (LFLAG(4))  f(t,p,v,lambda,u) in F
C                   (LFLAG(5))  T(t,p)*v in PDOT
C                   (LFLAG(7))  g(t,p) in G
C                   (LFLAG(8))  gI(t,p) in GI
C                   (LFLAG(10)) z(t,p,v) in GI      
C
C     ON RETURN:
C
C       AM(NV,NV) This real matrix contains the mass matrix M(t,p).
C       GP(LDG,NP) This real matrix contains the constraint matrix G(t,p).
C       F(NV)    This real array contains the forces f(t,p,v).
C       PDOT(NP) This real array contains the derivative p'=v.
C       G(NL)    This real array contains the constraints g(t,p).
C       GI(NL)   This real array contains the constraint derivative gI(t,p)
C                or second derivative terms z(t,p,v) (depending on LFLAG).
C       IFAIL    This integer flag indicates a stop condition 
C                iff IFAIL .NE. 0.  (Unused).
C
C     SYSTEM PARAMETERS:
C
C       IPAR(1)  K = Number of elements for flexible connecting rod.
C                ( NP = 3 + 2*K+2 =   3 rigid motion coordinates p
C                                  +  2*K+2 FE coordinates q        )
C       RPAR(1)  Drive torque omega
C       RPAR(2)  T1 for drive torque
C       RPAR(3)  Gravity constant
C       RPAR(4:K+3) Array of element lengths for FE ansatz
C
        INTEGER  K
        REAL*8   GRAV, OMEGA, T1, J1, J2, L1, L2, M1, M2, M3, PI,
     *           EE, NUE, BB, HH, RHO
C
C       Data set 2
C
        PARAMETER( M1 = 0.36D0,     M2 = 0.151104D0,
     *             M3 = 0.075552D0,
     *             L1 = 0.15D0,     L2 = 0.30D0,
     *             J1 = 0.002727D0, J2 = 0.0045339259D0,
     *             PI = 3.14159265358979D0,
     *             EE = .20D12,     NUE= 0.30D0,
     *             BB = 0.0080D0,   HH = 0.0080D0,
     *             RHO= 7870.0D0                         )
C
C     LOCAL VARIABLES:
C
C       Q, QD for FE coefficients and time derivatives, 
C       MQ, KQ, BQ, DQ for FE matrices.
C       Up to K=10 elements possible, i.e. max. dimension Q = 22
C
        INTEGER  NQMAX, I, J, NQ
        PARAMETER( NQMAX = 22 )
        REAL*8   Q(NQMAX), QD(NQMAX), MQ(NQMAX,NQMAX), 
     *           KQ(NQMAX,NQMAX), BQ(NQMAX), DQ(NQMAX),
     *           MQQ(NQMAX), KQQ(NQMAX), BQTQ, BQTQD,
     *           QTMQQ, QDTMQQ, DDOT,
     *           COSP1, COSP2, SINP1, SINP2, COSP12, SINP12
        SAVE     MQ, KQ, BQ, DQ
C
C       FIRST for first call - evaluation of FE matrices.
C
        LOGICAL  FIRST
        DATA     FIRST / .TRUE. /       
C
C_____________End of declaration part_________________________
C
      IFAIL = 0
      K     = IPAR(1)
      NQ    = 2*K+2
      IF (FIRST) THEN
C
C       Evaluate FE matrices and incorporate
C       boundary conditions left end of rod.
C
        CALL COMPIL ( K, NQMAX, NQ, RHO, EE, NUE, HH, BB,
     *                 RPAR(4), MQ, KQ, BQ, DQ  ) 
        DO 4 I=1,2
          DO 3 J=1,NQ
            MQ(I,J) = 0.D0
            KQ(I,J) = 0.D0
            MQ(J,I) = 0.D0
            KQ(J,I) = 0.D0
  3       CONTINUE
          MQ(I,I) = 1.0D0
          KQ(I,I) = 1.0D0
          BQ(I)   = 0.D0
          DQ(I)   = 0.D0
  4     CONTINUE
        FIRST = .FALSE.
      END IF
C
c     write(*,*)  'mq from compile'
c     write(*,988) ((mq(i,j),i=1,nq),j=1,nq)
c988  format(10G15.3)
C     Some initializations.
C
      OMEGA  = RPAR(1)
      T1     = RPAR(2)
      GRAV   = RPAR(3)
C
      COSP1  = COS(P(1))
      COSP2  = COS(P(2))
      SINP1  = SIN(P(1))
      SINP2  = SIN(P(2))
      COSP12 = COS(P(1)-P(2))
      SINP12 = SIN(P(1)-P(2))
C
      DO 6 I=1,NQ
         Q(I)  = P(3+I)
         QD(I) = V(3+I)
   6  CONTINUE
C
C     Evaluate scalar products and quadratic forms.
C
      BQTQ  = DDOT(NQ,BQ,1,Q,1)
      BQTQD = DDOT(NQ,BQ,1,QD,1)
      DO 10 I=1,NQ
         MQQ(I) = DDOT(NQ,MQ(1,I),1,Q,1)
         KQQ(I) = DDOT(NQ,KQ(1,I),1,Q,1)
  10  CONTINUE 
      QTMQQ = DDOT(NQ,Q,1,MQQ,1)
      QDTMQQ= DDOT(NQ,QD,1,MQQ,1)
      IF (LFLAG(1)) THEN 
C
C         Mass matrix - rigid motion entries.
C
          AM(1,1) = J1 + M2*L1*L1
          AM(1,2) = .5*L1*L2*M2*COSP12 
          AM(1,3) = 0.D0
          AM(2,2) = J2
          AM(2,3) = 0.D0
          AM(3,3) = M3
C
C         Superposition of flexible motion (matrix M^e).
C         
          AM(1,2) = AM(1,2) + RHO*L1*SINP12*BQTQ
          AM(2,2) = AM(2,2) + QTMQQ
C
C         Coupling terms (matrix C^T).
C
          DO 100 I=1,NQ
             AM(1,3+I) = RHO*L1*COSP12*BQ(I) 
             AM(2,3+I) = RHO*DQ(I)
             AM(3,3+I) = 0.D0
  100     CONTINUE
C
C         FE mass matrix.
C
          DO 120 I=1,NQ
             DO 110 J=1,I
                AM(3+J,3+I) = MQ(J,I)
  110        CONTINUE
  120     CONTINUE
C         
C         Lower left triangle (symmetry not exploited).
C
          DO 140 I=1,NV
             DO 130 J=I+1,NV
                AM(J,I) = AM(I,J)
  130        CONTINUE
  140     CONTINUE
      END IF
c
c     write(*,*) 'mass from fmbse'
c     write(*,989) ((am(i,j),j=1,nv),i=1,nv)
c 989  format(13G13.2)
      IF (LFLAG(2) .OR. LFLAG(3)) THEN
C
C         Constraint Jacobian.
C 
          GP(1,1) = L1*COSP1
          GP(1,2) = L2*COSP2 - Q(NQ-1)*SINP2
          GP(1,3) = 0.D0
          GP(2,1) = L1*SINP1
          GP(2,2) = L2*SINP2 + Q(NQ-1)*COSP2
          GP(2,3) = 1.D0
          DO 150 I=1,NQ
             GP(1,3+I) = 0.D0
             GP(2,3+I) = 0.D0
  150     CONTINUE
          GP(1,2+NQ) = COSP2
          GP(2,2+NQ) = SINP2
      END IF
      IF (LFLAG(4)) THEN 
C
C         Forces - rigid motion entries.
C
          F(1) = -.5*L1*GRAV*(M1+2.*M2)*COSP1 
     &           -.5*L1*L2*M2*V(2)*V(2)*SINP12
          IF (T .LE. T1) THEN
             F(1) =  F(1)+OMEGA/T1*(1.D0-COS(2.D0*PI*T/T1))
          END IF
          F(2) = -.5*L2*GRAV*M2*COSP2
     &           +.5*L1*L2*M2*V(1)*V(1)*SINP12
          F(3) = 0.d0
C
C         Superposition of flexible motion (term f^e).
C
          F(1) = F(1) + RHO*L1*V(2)*V(2)*COSP12*BQTQ
     &             - 2.*RHO*L1*V(2)*SINP12*BQTQD
          F(2) = F(2) - RHO*L1*V(1)*V(1)*COSP12*BQTQ
     &                - 2.*V(2)*QDTMQQ 
     &                + RHO*GRAV*SINP2*BQTQ
C
C         Coriolis and gravity terms flexible motion (Gamma).
C
          DO 200 I=1,NQ
             F(3+I) = V(2)*V(2)*MQQ(I) 
     &                + RHO*L1*V(1)*V(1)*SINP12*BQ(I)
     &                - RHO*GRAV*COSP2*BQ(I)
  200     CONTINUE
C
C         Stiffness term K q.
C
          DO 210 I=1,NQ
             F(3+I) = F(3+I) - KQQ(I) 
  210     CONTINUE
      END IF
      IF (LFLAG(5)) THEN 
          DO 300 I=1,NP
             PDOT(I) = V(I)
  300     CONTINUE
      END IF
      IF (LFLAG(7)) THEN
            G(1) = L1*SINP1 + L2*SINP2 + Q(NQ-1)*COSP2
            G(2) = P(3) - L1*COSP1 - L2*COSP2 + Q(NQ-1)*SINP2
      END IF
      IF (LFLAG(8)) THEN
            GI(1) = 0.0D0
            GI(2) = 0.0D0
      END IF
      IF (LFLAG(10)) THEN 
            GI(1) = -L1*SINP1*V(1)*V(1) - L2*SINP2*V(2)*V(2)
     *       -2.*SINP2*V(2)*QD(NQ-1) - COSP2*V(2)*V(2)*Q(NQ-1) 
            GI(2) =  L1*COSP1*V(1)*V(1) + L2*COSP2*V(2)*V(2)
     *       +2.*COSP2*V(2)*QD(NQ-1) - SINP2*V(2)*V(2)*Q(NQ-1)
      END IF
C _______________________________________________
C
      RETURN
      END

      SUBROUTINE COMPIL ( K, LD, NQ, RHO, EE, NU, HH, BB, DL, 
     *                    MQ, KQ, BQ, DQ  ) 
      INTEGER  NQ, K, LD
      REAL*8   RHO, EE, NU, HH, BB, DL(K), 
     *         MQ(LD,NQ), KQ(LD,NQ), BQ(NQ), DQ(NQ)
C 
C  .. Compilation of FE matrices and vectors ..
C     - cubic ansatz, 1D linear beam model -
C     - boundary conditions not included   -
C
C     ** written by Bernd Simeon, TH Darmstadt, 05/21/95 **
C
C  .. Parameters on entry ..
C
C       K         Number of elements
C       LD        Leading dimension of MQ, KQ
C       NQ        NQ = 2*K+2 Dimension of node variables
C       RHO       Mass density
C       EE        Young's modulus
C       NU        Poisson's number
C       HH        Height of beam
C       BB        Width of beam
C       DL(K)     Array of element lengths 
C        
C  .. On return ..
C
C       MQ(NQ,NQ) Mass matrix
C       KQ(NQ,NQ) Stiffness matrix
C       BQ(NQ)    Vector for integral over (psi)
C       DQ(NQ)    Vector for integral over (x*psi + y*y*psi')
C
C  .. Subroutines called ..
C
C       MASS, STIFF1, STIFF2, VEK0, VEK1, VEK4
C
C  .. Local variables ..
C
      INTEGER I, J, IK, K0
      REAL*8  ALPHA1, ALPHA2, ALPHA3, BETA1, BETA2, 
     *        SM(4,4), S1(4,4), S2(4,4), C0(4), C1(4), C4(4)
C
C_____________End of declaration part COMPIL_________________________
C
C  .. Initialize matrices and constants ..
C
      DO 20 I=1,NQ
         DO 10 J=1,NQ
            MQ(J,I) = 0.D0
            KQ(J,I) = 0.D0
 10      CONTINUE
         BQ(I) = 0.D0
         DQ(I) = 0.D0
 20   CONTINUE
      ALPHA1 = RHO*BB*HH
      ALPHA2 = RHO*BB*HH**3/12.D0
      ALPHA3 = BB*HH**3*EE/12.D0*(1.-NU)/(1.-NU-2.*NU*NU)
      BETA1  = BB*HH
      BETA2  = BB*HH**3/12.D0
C
C   .. Compilation ..
C
      CALL MASS  ( DL(1), ALPHA1, SM )
      CALL STIFF1( DL(1), ALPHA2, S1 )
      CALL STIFF2( DL(1), ALPHA3, S2 )
      CALL VEK0  ( DL(1), BETA1 , C0 )
      CALL VEK1  ( DL(1), BETA2 , C1 )
      CALL VEK4  ( DL(1), BETA1 , C4 )
      DO 40 I=1,4
         DO 30 J=1,4
            MQ(J,I) = SM(J,I) + S1(J,I)
            KQ(J,I) = S2(J,I)
 30      CONTINUE
         BQ(I) = C0(I)
         DQ(I) = C1(I) + C4(I)
 40   CONTINUE
C
      DO 100 IK=2,K
         K0 = 2*IK - 2
         CALL MASS  ( DL(IK), ALPHA1, SM )
         CALL STIFF1( DL(IK), ALPHA2, S1 )
         CALL STIFF2( DL(IK), ALPHA3, S2 )
         CALL VEK0  ( DL(IK), BETA1 , C0 )
         CALL VEK1  ( DL(IK), BETA2 , C1 )
         CALL VEK4  ( DL(IK), BETA1 , C4 )
         DO 60 I=1,4
            DO 50 J=1,4
               MQ(K0+J,K0+I) = MQ(K0+J,K0+I) + SM(J,I) + S1(J,I)
               KQ(K0+J,K0+I) = KQ(K0+J,K0+I) + S2(J,I)
 50         CONTINUE
            BQ(K0+I) = BQ(K0+I) + C0(I)
            DQ(K0+I) = DQ(K0+I) + C1(I) + C4(I)
  60      CONTINUE
100   CONTINUE
C__________________________________________________________________
C
      RETURN
      END
C
      SUBROUTINE MASS ( L, ALPHA, SM )
      REAL*8   L, ALPHA, SM(4,4)
C
C     Element mass matrix cubic ansatz.
C
      INTEGER I,J
      REAL*8  FAC
C
      FAC = ALPHA*L/420.D0
C
      SM(1,1) =   FAC*156.
      SM(1,2) =   FAC*22.*L
      SM(1,3) =   FAC*54.
      SM(1,4) = - FAC*13.*L
      SM(2,2) =   FAC*4.*L*L
      SM(2,3) =   FAC*13.*L
      SM(2,4) = - FAC*3.*L*L
      SM(3,3) =   FAC*156.
      SM(3,4) = - FAC*22.*L
      SM(4,4) =   FAC*4.*L*L
C
      DO 20 I=1,4
         DO 10 J=1,I-1
            SM(I,J) = SM(J,I)
  10     CONTINUE
  20  CONTINUE
C
      RETURN
      END
C
      SUBROUTINE STIFF1 ( L, ALPHA, S1 )
      REAL*8   L, ALPHA, S1(4,4)
C
C     Element stiffness matrix 1. derivative.
C
      INTEGER I,J
      REAL*8  FAC
C
      FAC = ALPHA/(30.D0*L)
C
      S1(1,1) =   FAC*36.
      S1(1,2) =   FAC*3.*L
      S1(1,3) = - FAC*36.
      S1(1,4) =   FAC*3.*L
      S1(2,2) =   FAC*4.*L*L
      S1(2,3) = - FAC*3.*L
      S1(2,4) = - FAC*L*L
      S1(3,3) =   FAC*36.
      S1(3,4) = - FAC*3.*L
      S1(4,4) =   FAC*4.*L*L
C
      DO 20 I=1,4
         DO 10 J=1,I-1
            S1(I,J) = S1(J,I)
  10     CONTINUE
  20  CONTINUE
C
      RETURN
      END
C
      SUBROUTINE STIFF2 ( L, ALPHA, S2 )
      REAL*8   L, ALPHA, S2(4,4)
C
C     Element stiffness matrix 2. derivative.
C
      INTEGER I,J
      REAL*8  FAC
C
      FAC = ALPHA*2.D0/L**3
C
      S2(1,1) =   FAC*6.
      S2(1,2) =   FAC*3.*L
      S2(1,3) = - FAC*6.
      S2(1,4) =   FAC*3.*L
      S2(2,2) =   FAC*2.*L*L
      S2(2,3) = - FAC*3.*L
      S2(2,4) =   FAC*L*L
      S2(3,3) =   FAC*6.
      S2(3,4) = - FAC*3.*L
      S2(4,4) =   FAC*2.*L*L
C
      DO 20 I=1,4
         DO 10 J=1,I-1
            S2(I,J) = S2(J,I)
  10     CONTINUE
  20  CONTINUE
C
      RETURN
      END
C      
      SUBROUTINE VEK0 ( L, ALPHA, C0 )
      REAL*8   L, ALPHA, C0(4)
C
C     Element vector for int(psi).
C
      REAL*8  FAC
C
      FAC = ALPHA*L/12.
C
      C0(1) =   FAC*6.
      C0(2) =   FAC*L
      C0(3) =   FAC*6.
      C0(4) = - FAC*L
C
      RETURN
      END
C
      SUBROUTINE VEK1 ( L, ALPHA, C1 )
      REAL*8   L, ALPHA, C1(4)
C
C     Element vector for int(psi').
C
      REAL*8  FAC
C
      FAC = ALPHA
C
      C1(1) = - FAC
      C1(2) =   0.0D0
      C1(3) =   FAC
      C1(4) =   0.0D0
C
      RETURN
      END
C
      SUBROUTINE VEK4 ( L, ALPHA, C4 )
      REAL*8   L, ALPHA, C4(4)
C
C     Element vector for int(psi*x).
C
      REAL*8  FAC
C
      FAC = ALPHA*L*L/10.D0
C
      C4(1) =   FAC*1.5D0
      C4(2) =   FAC*L/3.
      C4(3) =   FAC*3.5D0
      C4(4) = - FAC*L*.5
C
      RETURN
      END

