c module that contains the subroutine to advect particles using the exact velocity
c that describes two equal and opposite plane waves

      module mod_twowave
      
      implicit none
      
      contains
      
      subroutine twowave_velocity(amp,cff1,cff2,cff3,cff4,x,y,u,v)
c subroutine that returns the exact velocity of two plane waves
c amp is the wave amplitude
c cff are the dimensionless wavenumbers
c x,y is the particle location
c u,v are the velocity components

      implicit none
      
      integer ii,jj
      real*8 amp,cff1,cff2,cff3,cff4,x,y,u,v,pi
      
      data pi/3.14159265358979323846D0/
      
      u = amp*(cff2*dsin(cff1*x+cff2*y) + cff4*(dsin(cff3*x + cff4*y)))
      v = -amp*(cff1*dsin(cff1*x+cff2*y) + cff3*(dsin(cff3*x + cff4*y)))
      
      end subroutine twowave_velocity
c -----------------------------------------------------------------------
c ----------------------------------------------------------------------
      
      subroutine propagtating_twowave_velocity(amp,cff1,cff2,cff3,cff4
     &      ,x,y,time,u,v)
c subroutine that returns the exact velocity of two propagating plane waves
c amp is the wave amplitude
c time is the time (multipled by factor of 2pi and divided by the period)
c cff are the dimensionless wavenumbers
c x,y is the particle location
c u,v are the velocity components

      implicit none
      
      integer ii,jj
      real*8 amp,cff1,cff2,cff3,cff4,x,y,u,v,pi,time
      
      data pi/3.14159265358979323846D0/
      
      u = amp*(cff2*dsin(cff1*x+cff2*y-time) 
     &      + cff4*(dsin(cff3*x + cff4*y-time)))
      v = -amp*(cff1*dsin(cff1*x+cff2*y-time)
     &  + cff3*(dsin(cff3*x + cff4*y-time)))
      
      end subroutine propagtating_twowave_velocity
      
c -------------------------------------------------------------------------
c -------------------------------------------------------------------------
      
      subroutine rk4_twowave(x,y,cff1,cff2,cff3,cff4
     &  ,amp,dt,x_diff,y_diff)
c subroutine that advects a particle at (x,y) by velocity using rk4
c for a time t

c cff are dimensionless wave numbers 
c (x,y) is particle locadion
c amp is amplitude of wave
c x_diff,y_diff are the x,y displacements

      implicit none
      
      real*8 x,y,cff1,cff2,cff3,cff4,amp,x_diff,y_diff,time_old
     & ,time_half,time_new,dt 
      
      real*8 u1,v1,u2,v2,u3,v3,u4,v4,xt,yt,dt2
      
      dt2 = dt/2.
      
      
c perform rk4
      
      call twowave_velocity(amp,cff1,cff2,cff3,cff4,x,y
     &      ,u1,v1)
      
      xt = x + dt2*u1
      yt = y + dt2*v1
      
      call twowave_velocity(amp,cff1,cff2,cff3,cff4,xt,yt
     &      ,u2,v2)
      
      xt = x + dt2*u2
      yt = y + dt2*v2
      
      call twowave_velocity(amp,cff1,cff2,cff3,cff4,xt,yt
     &      ,u3,v3)
      
      xt = x + dt*u3
      yt = y + dt*v3
      
      call twowave_velocity(amp,cff1,cff2,cff3,cff4,xt,yt
     &      ,u4,v4)  
      
      x_diff = (dt/6.)*(u1+2*(u2+u3)+u4)
      y_diff = (dt/6.)*(v1+2*(v2+v3)+v4)   
      
      end subroutine
      
c ------------------------------------------------------------------------
c ------------------------------------------------------------------------
      
      subroutine exact_twowave(x,y,cff1,cff2,cff3,cff4
     &  ,amp,dt,x_diff,y_diff)
     
c subroutine that advects particles using two plane waves
c uses no approximation - completely accurate

      implicit none
      
      real*8 x,y,cff1,cff2,cff3,cff4,amp,x_diff,y_diff,dt
      
      real*8 u,v
      
      call twowave_velocity(amp,cff1,cff2,cff3,cff4,x,y,u,v)
      
      x_diff = u*dt
      y_diff = v*dt      
      
      end subroutine exact_twowave

      end module mod_twowave
