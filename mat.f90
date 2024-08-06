MODULE global
implicit none
INTEGER, PARAMETER, PRIVATE  :: dp = selected_real_kind(15,307)
integer :: i,step,nlines,n,io,nlines2
     
          
real(dp), allocatable :: x(:), V(:), dV(:)
real (dp) :: h, junk, V0, Vnp1
 
END module global
    
    
    
program one
use global
implicit none
INTEGER, PARAMETER  :: dp = selected_real_kind(15, 307)!selected_real_kind(15, 307)
integer, parameter :: timestep=500000,m=2
integer :: j,nlines3
real(dp):: VV,XX,dVV
real(dp), parameter :: epsi=1d-5
real(dp) :: t(timestep),t0,x0,dx0,q(m),dqdt(m)
real(dp) :: pi,rm,initial_sample(300)
real(dp), parameter :: dt=0.0001d0 
character(len=512) :: chr,ic,ic2(20)
real(dp) :: dat(500),dat2(500)!,dat3
 real(dp) :: x_fr_min,x_fr_max,xacc,rtsafe,rtnewt
external findroot,dvv 
character (len=512) :: initial_condition
real(dp) :: findroot_mathematica,epsilonH(timestep),rho_stop
real(dp) ::mp2,hubble2(timestep),rho(timestep),stoping_con
real(dp) :: velocity(timestep),delta_psi(timestep),eH(timestep)
real(dp):: lh(1000),tt(1000),aa
 real(dp), allocatable :: nfinal(:),x_init(:),num_der(:),num_der2(:),x_init1(:)
real(dp), allocatable :: phi(:), nn(:),x_initial2(:)
    
integer, parameter :: n_gauleg=1000
real(dp), allocatable :: h2(:), nh(:),yh2(:)
real(dp), allocatable ::ddelta_psi(:)
real(dp):: nhmin,nhmax,xh,wh(n_gauleg)
real(dp) :: yp1h,yp2h,int1,h2_new(n_gauleg),ypnh
real (dp) :: nh_new,h_res!(n_gauleg)
             
             
open(unit=99, file="potm.txt",status= "unknown") !file3.txt       
!open(unit=97, file="field.txt",status= "unknown") 
open(unit=95, file="deh.txt", status="unknown")
open(unit=96, file="res.txt", status="unknown")
!open(unit=94, file="hubble2.txt", status="unknown")
     
open(unit=93, file= "res2.txt",status="unknown")
            
pi=dacos(-1.d0)
mp2=1.d0!/ (8.d0*pi) !the potential is given in reduced mass planck unit 
nlines = 0
t0=0
         
      do 
      read(99,*,iostat=io)
      if (io/=0) exit  
      nlines = nlines + 1
      enddo
      close(99)
      
      n=nlines-2
      
      open(unit=99, file="potm.txt",status= "unknown") !for mathematica code - 3 loops corrections
      read(99,*)
      read(99,*) h, junk
      close(99)

      allocate(x(n),V(n),dv(n))
      open(unit=99, file="potm.txt",status= "unknown") ! !for mathematica file
      read(99,*) junk, V0
      do i=1,n
      read(99,*) x(i), V(i)
      enddo
      read(99,*) junk, Vnp1
      close(99)
     
     
 ! "check the potential after the interpolation: "    
      dV(1)=(V(1)-V0)/(2.0*h)
      do i=2,n-1
      dV(i)=(V(i+1)-V(i-1))/(2.0*h)
      enddo
      dV(n)=(Vnp1-V(n))/(2.0*h)
      
      !do i=1,n
      !write(98,*) x(i), V(i), dV(i)
      !enddo
      !close(98)
! end of checking the potential for interpolation      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    nlines2=120!130!170!130!130!12!130!50!1000!nlines2=
   
      allocate(nfinal(nlines2),x_init(nlines2),x_init1(nlines2),x_initial2(nlines2))
      allocate(ddelta_psi(nlines2))
      
      open(unit=10, file= "findroot.txt",status="unknown")
     read(10,*) findroot_mathematica
    ! findroot_mathematica=findroot_mathematica*1.01
     !findroot_mathematica=0.03637239479135252d0!0.03637239479135285!0.03636340307734146!0.03668279524600589!0.036100782492118635!0.035734266693000664!
       !open(unit=11, file="initial_cont.txt", status="unknown")    
       !print*, findroot_mathematica
       !read(*,*)
       do i=1,timestep
       t(i)=t0+dt*(i-1)
       enddo


!solvling nlines2 the equation of the field
DO j=2,nlines2!20,20!5,5!360!10!120!300

       lh(j)=1-j
       lh(j)=(1-10**((lh(j))/20))!30!10!40!10!30!90
       x0=findroot_mathematica*lh(j)
      
       dx0=0
      !print*,j,lh(j)
      rm=10.d0* VV(x0)
      
        !read(*,*)     
      
       q(1)=x0
       q(2)=dx0
       hubble2(1)=(rm*exp(-3*t(1))+VV(q(1)))/3!/(3*mp2)

       rho(1)= rm*exp(t(1))+VV(q(1))+0.5d0*hubble2(1)*q(2)**2.d0
       !print*, x0

         call derivs(m,t(1),rm,q,dqdt)   !derivs(n,t,q,dqdt)
         call rk4(q,dqdt,m,t(2),dt,rm)  ! rk4(y,dydx,n,x,h)
         !call odeint(q,m,t(2),t(1),epsi,dt,rad) !(ystart,nvar,x1,x2,eps,h1) 
         !print*, q(1),q(2)!,rad
        ! READ(*,*)

        !  print*, "rho=", rho(1)
       epsilonH(1)=0.5*q(2)**2.d0*8*pi
eH(1)=((3*rm*dexp(-3*t(1))-dvv(q(1))*q(2))/dsqrt(2*(rm*exp(-3*t(1))+VV(q(1)))))
!print*,eH(1),hubble2(1)
      
         ! if(j.eq.nlines2) print*, "Hcmb:", dsqrt(hubble2(1)),  eH(1),epsilonH(1)!,eH, ,
      do i=2,timestep-1 
        
       call derivs(m,t(i),rm,q,dqdt)
       call rk4(q,dqdt,m,t(i+1),dt,rm) 
    
       !call odeint(q,m,t(i+1),t(i),epsi,dt,rad)! 
   
      
     
hubble2(i)=(rm*exp(-3*t(i))+VV(q(1)))/(3*mp2)

epsilonH(i)=0.5*q(2)**2.d0*8*pi
 
eH(i)=((4*rm*dexp(-4*t(i))-dvv(q(1))*q(2))/(2*(rm*exp(-4*t(i))+VV(q(1)))))!mp2!*(mp2)! you need to check if this is coreect
delta_psi(i)=hubble2(i)*(1- dexp(-2*eH(i)*t(i)))/(8*pi**2*eH(i))

ddelta_psi(j)=delta_psi(2)

  ! write(97,*) t(i), q(1),0.5*q(2)**2.d0*8*pi

rho(i+1)= rm*exp(-3*t(i+1))+VV(q(1))
stoping_con= rm*exp(-3*t(i))+VV(q(1))+0.5*hubble2(i)*q(2)**2
velocity(i)=q(2)
        
! "terminal conditions for stopping integration of de:"
        ! if (epsilonH(i).gt.1) exit! if (eH(i).gt.1) exit
      !if(stoping_con.le.(rho(1)/100)) exit 
if( (velocity(i)**2*hubble2(i)).lt.(velocity(i-1)**2*hubble2(i-1))) exit

         
 enddo
print*, "number of efolds: ", t(i), q(1),j!x0!, j
      
     !  write(68,*) t(i), x0
        nfinal(j)=t(i)
        x_init(j)=log10(findroot_mathematica-x0)!findroot_mathematica-x0!x0!
        x_initial2(j)=(findroot_mathematica-x0)
        x_init1(j)= x0
ENDDO  
      
! "end of solvling nlines2 the equation of the field"
      
 
        do j=2,nlines2 !check the previous numbers
       write(96,*)  nfinal(j), x_init(j),x_init1(j)!,ddelta_psi(j),x0!,hubble2(i)!,x0!q(1)!j!t(i), log10(findroot_mathematica-x0)!
       !if ((nfinal(j)).gt.(nfinal(j-1))) write(96,*)  nfinal(j), x_init(j),q(1)
      ! if (j.gt.5) write(93,*)  nfinal(j),x_init1(j)    ! in order to have a better fit 
       write(93,*)  nfinal(j),x_init1(j)
        !write(95,*) nfinal(j),x_initial2(j)
!     write(94,*) nfinal(j), ((rm*exp(-4*nfinal(j))+VV(q(1)))/(3*mp2))/(4*pi**2) 
       ENDDO
       close(96)
       !close(94)
       
            allocate(num_der(nlines2),phi(nlines2),nn(nlines2))  
       
         do i =2 ,nlines2-1          
        num_der(i)= (x_init1(i+1)- x_init1(i))/(nfinal(i+1)- nfinal(i))
        if((nfinal(i+1)).gt.( nfinal(i))) then!.ne.
        write(95,*)  nfinal(i),  num_der(i)
        endif
       enddo     
       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     


    END
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!!!!!!!!!!!!        END       !!!!!!!!!!!!!!!!!!        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       

      function findroot(x) 
      implicit none
            INTEGER, PARAMETER  :: dp = selected_real_kind(15, 307)
     integer :: i
     real(dp) :: x , findroot,dvv,vv
           !open (unit=59, file="check.txt", status="unknown")
 !     
      DO I =1,100
       X=0.02 +(I-1)*0.0001D0
      findroot =dvv(x)
!    WRITE(59,*) XX, VV(xx), dVV(xx)
      ENDDO
     
!      close(59)
      
      end
       
      FUNCTION VV(xx)
      use global
      implicit none
      INTEGER, PARAMETER  :: dp1 = selected_real_kind(15, 307)
      real(dp1) :: xx, VV
      real(dp1), allocatable :: y2(:)
      allocate(y2(n))
      !open (unit=59, file="check.txt", status="unknown")
      call spline(x,V,n,dV(1),dV(n),y2)
      call splint(x,V,y2,n,xx,VV)
     
      END
       
       
      FUNCTION dVV(x)
      implicit none
      integer, parameter :: dp = selected_real_kind(15, 307)
      INTEGER :: NTAB
      REAL (dp) :: dVV,err,h,x,VV,CON,CON2,BIG,SAFE
      PARAMETER (CON=1.4,CON2=CON*CON,BIG=1.E30,NTAB=10,SAFE=2.)
      INTEGER :: i,j
      REAL (dp) :: errt,fac,hh,a(NTAB,NTAB)
      h=0.0001d0
      if(h.eq.0.) pause 'h must be nonzero in dfridr'
      hh=h
      a(1,1)=(VV(x+hh)-VV(x-hh))/(2.0*hh)
      err=BIG
      do 12 i=2,NTAB
        hh=hh/CON
        a(1,i)=(VV(x+hh)-VV(x-hh))/(2.0*hh)
        fac=CON2
        do 11 j=2,i
          a(j,i)=(a(j-1,i)*fac-a(j-1,i-1))/(fac-1.)
          fac=CON2*fac
          errt=max(abs(a(j,i)-a(j-1,i)),abs(a(j,i)-a(j-1,i-1)))
          if (errt.le.err) then
            err=errt
            dVV=a(j,i)
          endif
11      continue
        if(abs(a(i,i)-a(i-1,i-1)).ge.SAFE*err)return
12    continue
      return
      END       
       

      
      
      SUBROUTINE rk4(y,dydx,n,x,h,rad)
      implicit none
      integer, parameter :: dp = selected_real_kind(15, 307)
      INTEGER :: n
      REAL(dp) :: h,x,dydx(n),y(n),yout(n),rad
      integer, parameter :: NMAX=50
      INTEGER :: i
      REAL(dp) :: h6,hh,xh,dym(NMAX),dyt(NMAX),yt(NMAX)
      hh=h*0.5d0
      h6=h/6.d0
      xh=x+hh
      do i=1,n
       yt(i)=y(i)+hh*dydx(i)
       enddo
       call derivs(n,xh,rad,yt,dyt)
       do i=1,n
       yt(i)=y(i)+hh*dyt(i)
       enddo
       call derivs(n,xh,rad,yt,dym)
       do i=1,n
       yt(i)=y(i)+h*dym(i)
       dym(i)=dyt(i)+dym(i)
       enddo
       call derivs(n,x+h,rad,yt,dyt)
       do i=1,n
       yout(i)=y(i)+h6*(dydx(i)+dyt(i)+2.d0*dym(i))
       enddo
       do i=1,n
       y(i)=yout(i)
       enddo
       return
       END  
       
       
       
      subroutine derivs(n,t,rad,q,dqdt)
      !use global
      implicit none
      integer, parameter :: dp = selected_real_kind(15, 307)
      integer :: n
      real(dp) :: t,q(n),dqdt(n)
      real(dp) :: VV,dVV,mp,v0,x0
      real(dp) :: h2,rad,pi
      pi = dacos(-1.d0)

      mp=1.d0!/ dsqrt(8.d0*pi) !1.22d19
        !rad= 1.068351022008d-11!1.08129d-11!3.70082d-14 !3.78d-14!3.70082d-14!3.27717d-13!2.29841d-12!prosoxi den exei mpei sosta
        !print*, VV(x0)
        !read(*,*)
      !if(q(1).lt.0) then 
       !  VV(q(1))=abs(VV(q(1)))
           !endif
      h2=rad*exp(-3*t)+VV(q(1))
      h2=h2/(3*mp**2)
      
      
      dqdt(1)=q(2)
      dqdt(2)=-q(2)*(-3*rad*EXP(-3*t)+q(2)*dVV(q(1)))/ (2*(rad*exp(-3*t)+ABS(VV(q(1)))))- 3*q(2)- dVV(q(1))/h2
      
           !open (unit=59, file="check.txt", status="unknown")

     
      return
      end





      SUBROUTINE spline(x,y,n,yp1,ypn,y2)
      implicit none
      integer, parameter :: dp = selected_real_kind(15, 307)
      INTEGER :: n,NMAX
      REAL (dp) :: yp1,ypn,x(n),y(n),y2(n)
      PARAMETER (NMAX=10000)
      INTEGER :: i,k
      REAL (dp) :: p,qn,sig,un,u(NMAX)
      if (yp1.gt..99e30) then
        y2(1)=0.
        u(1)=0.
      else
        y2(1)=-0.5
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do 11 i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.
        y2(i)=(sig-1.)/p
        u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
11    continue
      if (ypn.gt..99e30) then
        qn=0.
        un=0.
      else
        qn=0.5
        un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
      do 12 k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
12    continue
      return
      END
      
      
      
      SUBROUTINE splint(xa,ya,y2a,n,x,y)
      implicit none
      integer, parameter :: dp = selected_real_kind(15, 307)
      INTEGER :: n
      REAL (dp) :: x,y,xa(n),y2a(n),ya(n)
      INTEGER :: k,khi,klo
      REAL (dp) :: a,b,h
      klo=1
      khi=n
1     if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else
          klo=k
        endif
      goto 1
      endif
      h=xa(khi)-xa(klo)
      if (h.eq.0.) pause 'bad xa input in splint'
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
      return
      END




   

      SUBROUTINE odeint(ystart,nvar,x1,x2,eps,h1,rad)
      IMPLICIT NONE
            integer, parameter :: dp = selected_real_kind(15, 307)
      INTEGER :: nbad,nok,nvar,KMAXX,MAXSTP,NMAX
      REAL(dp) ::  eps,h1,hmin,x1,x2,ystart(nvar),TINY
      EXTERNAL derivs,rkqs
      PARAMETER (MAXSTP=100000,NMAX=50,KMAXX=200,TINY=1.d-30)
      INTEGER :: i,kmax,kount,nstp
      REAL(dp) :: dxsav,h,hdid,hnext,x,xsav,dydx(NMAX),xp(KMAXX),y(NMAX),yp(NMAX,KMAXX),yscal(NMAX)
      real (dp) :: rad
      COMMON /path/ kmax,kount,dxsav,xp,yp
      hmin=0.d0
      x=x1
      h=sign(h1,x2-x1)
      nok=0
      nbad=0
      kount=0
      do 11 i=1,nvar
        y(i)=ystart(i)
11    continue
      if (kmax.gt.0) xsav=x-2.d0*dxsav
      do 16 nstp=1,MAXSTP
        call derivs(nvar,x,rad,y,dydx) !derivs(n,t,rad,q,dqdt)      
        do 12 i=1,nvar
          yscal(i)=abs(y(i))+abs(h*dydx(i))+TINY
12      continue
        if(kmax.gt.0)then
          if(abs(x-xsav).gt.abs(dxsav)) then
            if(kount.lt.kmax-1)then
              kount=kount+1
              xp(kount)=x
              do 13 i=1,nvar
                yp(i,kount)=y(i)
13            continue
              xsav=x
            endif
          endif
        endif
        if((x+h-x2)*(x+h-x1).gt.0.) h=x2-x
       call rkqs(y,dydx,nvar,x,h,eps,yscal,hdid,hnext,rad)
        if(hdid.eq.h)then
          nok=nok+1
        else
          nbad=nbad+1
        endif
        if((x-x2)*(x2-x1).ge.0.)then
          do 14 i=1,nvar
            ystart(i)=y(i)
14        continue
          if(kmax.ne.0)then
            kount=kount+1
            xp(kount)=x
            do 15 i=1,nvar
              yp(i,kount)=y(i)
15          continue
          endif
          return
        endif
        if(abs(hnext).lt.hmin) pause 'stepsize smaller than minimum in odeint'
        h=hnext
16    continue
       pause 'too many steps in odeint'
      return
      END

      SUBROUTINE rkqs(y,dydx,n,x,htry,eps,yscal,hdid,hnext,rad )
      IMPLICIT NONE
            integer, parameter :: dp = selected_real_kind(15, 307)
      INTEGER :: n,NMAX
      REAL(dp) :: eps,hdid,hnext,htry,x,dydx(n),y(n),yscal(n),rad
      EXTERNAL derivs
      PARAMETER (NMAX=50)
!CU    USES derivs,rkck
      INTEGER :: i
      REAL(dp) :: errmax,h,htemp,xnew,yerr(NMAX),ytemp(NMAX),SAFETY,PGROW,PSHRNK,ERRCON
      PARAMETER(SAFETY=0.9d0,PGROW=-0.2d0,PSHRNK=-0.25d0,ERRCON=1.89d-4)
      h=htry
1     call rkck(y,dydx,n,x,h,ytemp,yerr,rad)
      errmax=0.d0
      do 11 i=1,n
        errmax=max(errmax,abs(yerr(i)/yscal(i)))
11    continue
      errmax=errmax/eps
      if(errmax.gt.1.)then
        htemp=SAFETY*h*(errmax**PSHRNK)
        h=sign(max(abs(htemp),0.1d0*abs(h)),h)
        xnew=x+h
        if(xnew.eq.x)pause 'stepsize underflow in rkqs'
        goto 1
      else
        if(errmax.gt.ERRCON)then
          hnext=SAFETY*h*(errmax**PGROW)
        else
          hnext=5.d0*h
        endif
        hdid=h
        x=x+h
        do 12 i=1,n
          y(i)=ytemp(i)
12      continue
        return
      endif
      END


      SUBROUTINE rkck(y,dydx,n,x,h,yout,yerr,rad)
      IMPLICIT NONE
            integer, parameter :: dp = selected_real_kind(15, 307)
      INTEGER :: n,NMAX
      REAL(dp) :: h,x,dydx(n),y(n),yerr(n),yout(n),rad
      EXTERNAL derivs
      PARAMETER (NMAX=50)
!CU    USES derivs
      INTEGER :: i
      REAL(dP) :: ak2(NMAX),ak3(NMAX),ak4(NMAX),ak5(NMAX),ak6(NMAX)
     REAL(dP) ::ytemp(NMAX),A2,A3,A4,A5,A6,B21,B31,B32,B41,B42,B43,B51,B52,B53
     REAL(dP) ::B54,B61,B62,B63,B64,B65,C1,C3,C4,C6,DC1,DC3,DC4,DC5,DC6
      PARAMETER(A2=0.2d0,A3=0.3d0,A4=0.6d0,A5=1.d0,A6=0.875d0,B21=0.2d0)
     PARAMETER(B31=3.d0/40.d0,B32=9.d0/40.d0,B41=0.3d0,B42=-0.9d0,B43=1.2d0)
     PARAMETER(B51=-11.d0/54.d0,B52=2.5d0,B53=-70.d0/27.d0,B54=35.d0/27.d0)
     PARAMETER(B61=1631.d0/55296.d0,B62=175.d0/512.d0,B63=575.d0/13824.d0)
     PARAMETER(B64=44275.d0/110592.d0,B65=253.d0/4096.d0,C1=37.d0/378.d0)
     PARAMETER(C3=250.d0/621.d0,C4=125.d0/594.d0,C6=512.d0/1771.d0)
     PARAMETER(DC1=C1-2825.d0/27648.d0,DC3=C3-18575.d0/48384.d0)
     PARAMETER(DC4=C4-13525.d0/55296.d0,DC5=-277.d0/14336.d0)
     PARAMETER(DC6=C6-0.25d0)
      do 11 i=1,n
        ytemp(i)=y(i)+B21*h*dydx(i)
11    continue
      call derivs(n,x+A2*h,rad,ytemp,ak2)!derivs(n,t,rad,q,dqdt)
      do 12 i=1,n
        ytemp(i)=y(i)+h*(B31*dydx(i)+B32*ak2(i))
12    continue
      call derivs(n,x+A3*h,rad,ytemp,ak3)
      do 13 i=1,n
        ytemp(i)=y(i)+h*(B41*dydx(i)+B42*ak2(i)+B43*ak3(i))
13    continue
      call derivs(n,x+A4*h,rad,ytemp,ak4)!n,x+h,yt,dyt
      do 14 i=1,n
        ytemp(i)=y(i)+h*(B51*dydx(i)+B52*ak2(i)+B53*ak3(i)+B54*ak4(i))
14    continue
      call derivs(n,x+A5*h,rad,ytemp,ak5)
      do 15 i=1,n
        ytemp(i)=y(i)+h*(B61*dydx(i)+B62*ak2(i)+B63*ak3(i)+B64*ak4(i)+B65*ak5(i))
15    continue
      call derivs(n,x+A6*h,rad,ytemp,ak6)
      do 16 i=1,n
        yout(i)=y(i)+h*(C1*dydx(i)+C3*ak3(i)+C4*ak4(i)+C6*ak6(i))
16    continue
      do 17 i=1,n
        yerr(i)=h*(DC1*dydx(i)+DC3*ak3(i)+DC4*ak4(i)+DC5*ak5(i)+DC6*ak6(i))
17    continue
      return
      END



      FUNCTION rtsafe(funcd,x1,x2,xacc)
          integer, parameter :: dp = selected_real_kind(15, 307)
      INTEGER :: MAXIT
      REAL(dp) :: rtsafe,x1,x2,xacc
      EXTERNAL funcd
      PARAMETER (MAXIT=100)
      INTEGER :: j
      REAL(dp) :: df,dx,dxold,f,fh,fl,temp,xh,xl
      call funcd(x1,fl,df)
      call funcd(x2,fh,df)
  if((fl.gt.0..and.fh.gt.0.).or.(fl.lt.0..and.fh.lt.0.)) pause 'root must be bracketed in rtsafe'
      if(fl.eq.0.)then
        rtsafe=x1
        return
      else if(fh.eq.0.)then
        rtsafe=x2
        return
      else if(fl.lt.0.)then
        xl=x1
        xh=x2
      else
        xh=x1
        xl=x2
      endif
      rtsafe=.5*(x1+x2)
      dxold=abs(x2-x1)
      dx=dxold
      call funcd(rtsafe,f,df)
      do 11 j=1,MAXIT
        if(((rtsafe-xh)*df-f)*((rtsafe-xl)*df-f).ge.0..or. abs(2.*f).gt.abs(dxold*df) ) then
          dxold=dx
          dx=0.5*(xh-xl)
          rtsafe=xl+dx
          if(xl.eq.rtsafe)return
        else
          dxold=dx
          dx=f/df
          temp=rtsafe
          rtsafe=rtsafe-dx
          if(temp.eq.rtsafe)return
        endif
        if(abs(dx).lt.xacc) return
        call funcd(rtsafe,f,df)
        if(f.lt.0.) then
          xl=rtsafe
        else
          xh=rtsafe
        endif
11    continue
      pause 'rtsafe exceeding maximum iterations'
      return
      END


   FUNCTION rtnewt(funcd,x1,x2,xacc)
             integer, parameter :: dp = selected_real_kind(15, 307)
      INTEGER :: JMAX
      REAL(dp) :: rtnewt,x1,x2,xacc
      EXTERNAL funcd
      PARAMETER (JMAX=10000)
      INTEGER :: j
      REAL(dp) :: df,dx,f
      rtnewt=0.5d0*(x1+x2)
      do 11 j=1,JMAX
        call funcd(rtnewt,f,df)
        dx=f/df
        rtnewt=rtnewt-dx
        if((x1-rtnewt)*(rtnewt-x2).lt.0.) pause 'rtnewt jumped out of brackets'
        if(abs(dx).lt.xacc) return
11    continue
      pause 'rtnewt exceeded maximum iterations'
    END

    
     SUBROUTINE gauleg(x1,x2,x,w,n)
                  integer, parameter :: dp = selected_real_kind(15, 307)
      INTEGER:: n
        REAL(dp) :: x1,x2,x(n),w(n)
      DOUBLE PRECISION EPS
      PARAMETER (EPS=3.d-14)
      INTEGER:: i,j,m
      DOUBLE PRECISION p1,p2,p3,pp,xl,xm,z,z1
      m=(n+1)/2
      xm=0.5d0*(x2+x1)
      xl=0.5d0*(x2-x1)
      do 12 i=1,m
        z=cos(3.141592654d0*(i-.25d0)/(n+.5d0))
1       continue
          p1=1.d0
          p2=0.d0
          do 11 j=1,n
            p3=p2
            p2=p1
            p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
11        continue
          pp=n*(z*p1-p2)/(z*z-1.d0)
          z1=z
          z=z1-p1/pp
        if(abs(z-z1).gt.EPS)goto 1
        x(i)=xm-xl*z
        x(n+1-i)=xm+xl*z
        w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
        w(n+1-i)=w(i)
12    continue
      return
      END
