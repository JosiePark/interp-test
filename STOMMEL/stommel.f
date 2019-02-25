c
c Code that caclulates trajectories of particles due to a stommel flow
c Testing NETCDF

c isolve=0 --- 2dcubic interpolation
c isolve=1 --- bicubic interpolation
c isolve=2 --- no interpolation
c irk4=0 no rk4
c irk4=1 rk4
c set itraj = 1 to store points at different times so can determine trajectories        
        
        program main
        
        use mod_2dcubic
        use mod_bicubic
        use mod_laplace
        use mod_random
        use mod_traj_netcdf
        use mod_stommel
        use mod_rk4
c

      
      
      implicit none
      
      
      integer ii,jj,i,j
      
    
      REAL*8, PARAMETER :: Pi = 3.1415927
      
      parameter(ii=512,jj=ii)
      

      real*8 psi(ii,jj),x(ii),y(jj),L,eps,start,finish
      
      integer iseed,npoints,n,nrec,npoints_sqrt,k,nt,isolve,irk4
      
      parameter(npoints_sqrt=500,npoints = npoints_sqrt**2)
      
      real*8 x0(npoints),y0(npoints),a_mat(ii,jj,4,4),xt,yt,dt,dt05,dt6
      real*8 u_point0,v_point0,u_point1,v_point1,u_point2,v_point2
      real*8 u_point3,v_point3,x_diff,y_diff
      real*8 a(ii,jj),b(ii,jj),c(ii,jj),d(ii,jj),psi_interp_x(4)
      real*8 basinscale,scale,uscale,tscale,time_day,x_mean,y_mean
      real*8 time_day_halfstep, time_half,time_new,time_old
      real*8 u(ii,jj), v(ii,jj),tao0,dt_day

      character(100) :: file_name_x,file_name_y,file_name
c      parameter(file_name = 'stommel_traj.nc')
      
      integer l_tot_day,max_time1,max_time,l_day,k_o,k_t,nrec_t,k_out
      integer k_p,k_traj,k_save
        real *8, dimension(:,:),allocatable ::x_traj,y_traj
        real*8 k_s

      
      
      parameter(max_time = 100)
      
      parameter(isolve=0,irk4=0)
      
      parameter(basinscale = 520.D5)
      

      
c      data dt/86400./
c      data dt/43200./
c      data dt/21600./
c      data dt/10800./
      data dt/8640./
c      data dt/4320./
c      data dt/2160./
c      data dt/1080./
      
      dt_day = dt
      
      if (irk4.eq.0) then
      
      if (isolve.eq.0) then
         write(file_name,'("stommel_bicubic_dt_",i5,".nc")')int(dt_day)


      elseif(isolve.eq.1) then
        write(file_name,'("stommel_2Dcubic_dt_",i5,".nc")')int(dt_day)


      elseif(isolve.eq.2) then
         write (file_name, '("stommel_exact_dt_",i5,".nc")')int(dt_day)

      endif
      
      elseif(irk4.eq.1) then
      
      if (isolve.eq.0) then
      write(file_name,'("stommel_rk4_bicubic_dt_",i5,".nc")')int(dt_day)


      elseif(isolve.eq.1) then
      write(file_name,'("stommel_rk4_2Dcubic_dt_",i5,".nc")')int(dt_day)


      elseif(isolve.eq.2) then
      write (file_name, '("stommel_rk4_exact_dt_",i5,".nc")')int(dt_day)

      endif
      endif
      
      
      
      
c
c--- nondimensionalization
c
      uscale=1.
      scale=basinscale/dfloat(jj)
      tscale=scale/uscale
      
c--- counters
c
      dt=dt/tscale
      dt05=.5*dt
      dt6=dt/6.0

      k_out=1     ! diagnostics and output every k_out days
      k_traj=1     ! saving trajectories every k_traj days
      !k_save = 50
      


      l_tot_day=int(86400./(dt*tscale)+0.001)  ! time steps in one day
      print*, 'total time steps in a day',l_tot_day
      max_time1=max_time*l_tot_day             ! total number of time steps

      time_day=0.
      l_day=0
      k_o=0
      nrec=0

      k_t=0
      nrec_t=0
      
      !allocate(x_traj(npoints,max_time1),y_traj(npoints,max_time1))
      
      call create_NETCDF_file(file_name,npoints)

      
      do i = 1,ii
        
        x(i) = dfloat(i-1)/dfloat(ii)
        
      enddo
      
      do j = 1,jj
      
        y(j) = dfloat(j-1)/dfloat(jj)
        
      enddo
      
      !print*,'x,y=',x,y
      
      eps=0.04
      !tao0 = 100
      do i = 1,ii
        do j = 1,jj
        
            
        psi(i,j) = pi*(1- x(i) 
     &        - dexp(-x(i)/(eps)))*dsin(pi*y(j))
        u(i,j) = -pi**2*(1-x(i)-
     &   dexp(-x(i)/eps))*dcos(pi*y(j))
        v(i,j) = pi*(-1
     &   +dexp(-x(i)/eps)/eps)*dsin(pi*y(j))
        enddo
      enddo
      
      psi = psi*dfloat(ii)
      !u = u/dfloat(jj)
      !v = v/dfloat(ii)
      
      
      open(11, file= 'stommel_u.dat',form='unformatted')
      open(12, file= 'stommel_v.dat',form='unformatted')
      write(11) u
      write(12) v
      close(11)
      close(12)
      
      
      open(10,file= 'stommel.dat',form='unformatted')
      write(10) psi
      write(*,*) 'stommel.dat written'
      
      close(10)
      
      
      
      
      

      
      
      !call create_NETCDF(traj_name,npoints)

c ------------ LAGRANGIAN PARTICLES ------------------

        iseed = 102
        n=0
        nrec = 0
        do j = 1,npoints_sqrt
          do i = 1,npoints_sqrt
            n = n+1
            !x0(n) = ran1(iseed)*dfloat(ii)
            !y0(n) = ran1(iseed)*dfloat(jj)
            x0(n) = 150
            y0(n) = 250
            
          enddo
        enddo

      call cpu_time(start)
        
c      if (isolve.eq.0) then
        
c      call A_matrix(ii,jj,psi,A_mat)
      
c      else if (isolve.eq.1) then
      
c      call cubic_coeff_x(ii,jj,psi
c     &,a,b,c,d)
     
c      endif
      
      do k = 1,max_time1
  
         time_old=time_day
         time_day=time_day+dt*tscale/86400.
         time_day_halfstep=time_day-0.5D0*dt*tscale/86400.

         l_day=l_day+1

         if(l_day.eq.l_tot_day)then
           l_day=0
           k_o=k_o+1
           k_p=k_p+1
           k_t=k_t+1
         endif
         
      if (isolve.eq.0) then
        
      call A_matrix(ii,jj,psi,A_mat)
      
      else if (isolve.eq.1) then
      
      call cubic_coeff_x(ii,jj,psi
     &,a,b,c,d)
     
      endif
      
      if (irk4.eq.1) then
      if (isolve.eq.0) then
        
      call A_matrix(ii,jj,psi,A_mat)
      
      else if (isolve.eq.1) then
      
      call cubic_coeff_x(ii,jj,psi
     &,a,b,c,d)
     
      endif
      
      endif

      if(irk4.eq.1) then
      ! perform runge kutta forth order
      
      do n = 1,npoints
      
      if (isolve.eq.0) then
      
      call rk4_bicubic(ii,jj,x0(n),y0(n),dt,dfloat(0),a_mat,a_mat,a_mat
     & ,x_diff,y_diff)
      
      elseif (isolve.eq.1) then
      
      call rk4_2dcubic(ii,jj,x0(n),y0(n),dt,dfloat(0),a
     & ,b,c,d,a,b,c,d,a,b
     & ,c,d,x_diff,y_diff)
      
      elseif(isolve.eq.2) then
      
      call rk4_stommel(x0(n),y0(n),ii,jj,eps,L,dt,x_diff,y_diff)
      endif
      
      x0(n)=x0(n)+x_diff
      y0(n)=y0(n)+y_diff
       
       enddo
       
       elseif (irk4.eq.0) then
       
       do n = 1,npoints
      
      
            !write(*,*)'x,y=',x0(n), y0(n)
            
            if (isolve.eq.0) then
      call bicubic(ii,jj,a_mat,x0(n),y0(n),u_point0,v_point0)
      
      elseif (isolve.eq.1) then
      
      call cubic_poly_x(ii,jj,x0(n),y0(n)
     & ,a,b,c,d,psi_interp_x)

      call vel(ii,jj,psi_interp_x
     & ,a,b,c,d,x0(n),y0(n),u_point0,v_point0)
     
        
      
      elseif (isolve.eq.2) then
      call stommel_velocity(x0(n),y0(n)
     & ,ii,jj,eps,basinscale,u_point0,v_point0)
      endif
      !write(*,*) 'velocity calculated'
      
      !write(*,*)'first time step',u_point0,v_point0
      
      x0(n) = x0(n) + dt*(u_point0)
      y0(n) = y0(n) + dt*(v_point0)
      
     
       
       enddo
       
       endif
       
c--- Lagrangian diagnostics
c
         if(k_t.eq.k_traj)then
           write(*,*)'Diagnostics at time=',time_day
           k_o=0
           nrec=nrec+1
           k_t = 0
           
           x_mean=0.
           y_mean=0.
           do n=1,npoints
              x_mean=x_mean+x0(n)
              y_mean=y_mean+y0(n)
           enddo
           x_mean=x_mean/dfloat(npoints)
           y_mean=y_mean/dfloat(npoints)

           write(*,'(A8,F6.1,A10,F6.1)')'x1_mean=',x_mean
     &                               ,'; y1_mean=',y_mean
           
        
        
            call write_NETCDF_file(file_name,npoints,x0,y0
     & ,time_day
     & ,nrec)
     
        endif

        enddo
       
       call cpu_time(finish)
       
c       print*,file_name_x, file_name_y
       
       write(*,*) 'Time taken to run = ',finish-start
       print*, 'Time steps=',max_time1
      
      
      end program
