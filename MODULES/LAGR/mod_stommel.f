c module that contains the subroutine which returns the exact velocity of 
c a stommel flow

      module mod_stommel
      
      implicit none
      
      contains
      
      subroutine stommel_velocity(x0,y0,ii,jj,L,eps,u,v)
        
        implicit none
        
        integer ii,jj
        real*8 x,y,u,v,eps,x0,y0,tao0,L,scale
        
        REAL*8, PARAMETER :: Pi = 3.1415927
        
        x = x0/dfloat(ii)
        y = y0/dfloat(jj)
        
        scale = L/dfloat(ii)
       
        
        u = pi**2*(1-x-dexp(-x/eps))*dcos(pi*y)
        v= -pi*(-1+dexp(-x/eps)/eps)*dsin(pi*y)
        
        u = u!/dfloat(ii)
        v = v!/dfloat(ii)
        
        return
      
      end subroutine stommel_velocity
      
      subroutine rk4_stommel(x0,y0,ii,jj,eps,L,dt,x_diff,y_diff)
      
      implicit none
      
      real*8 x0,y0,eps,dt,x_diff,y_diff,dt2,L,xt,yt
      real*8 u1,v1,u2,v2,u3,v3,u4,v4
      integer ii,jj
      
      dt2 = dt/2.
      
      call stommel_velocity(x0,y0,ii,jj,L,eps,u1,v1)
      
      xt = x0 + dt2*u1
      yt = y0 + dt2*v1
      
      call stommel_velocity(xt,yt,ii,jj,L,eps,u2,v2)
      
      xt = x0 + dt2*u2
      yt = y0 + dt2*v2
      
      call stommel_velocity(xt,yt,ii,jj,L,eps,u3,v3)
      
      xt = x0 + dt*u3
      yt = y0 + dt*v3
      
      call stommel_velocity(xt,yt,ii,jj,L,eps,u4,v4)
      
      x_diff = (dt/6.)*(u1+2*(u2+u3)+u4)
      y_diff = (dt/6.)*(v1+2*(v2+v3)+v4)
      
      
      
      end subroutine rk4_stommel
      
      end module mod_stommel
