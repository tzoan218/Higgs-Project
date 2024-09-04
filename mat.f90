
    
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
