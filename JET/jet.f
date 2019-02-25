c code that advects particles by a snapshot of a jet
c no exact interpolation method possible
c no rk4 method used

      program jet
      
        use mod_2dcubic
        use mod_bicubic
        use mod_laplace
        use mod_random
        use mod_traj_netcdf
        use mod_rk4

      implicit none
      
      integer ii,jj,npoints,npoints_sqrt,iseed,max_time,isolve,max_time1
     &, l_tot_day,kk,l_day,n,step_tim,i,j,k_t,k_traj,nrec,irk4 
      real*8 U0,dt,tscale,uscale,basinscale,scale,time_day,new_tim
     &, finish, start  
      
      character*(*) qg_name
      character*100 file_name
      
      parameter(qg_name = 'QG.nc') ! name of file with jet stream function
      parameter(ii=512,jj=ii,U0=6.D0,npoints_sqrt = 500
     &, npoints = npoints_sqrt**2)
      parameter(max_time = 100) ! run time in days
      parameter(isolve = 1)
      parameter(basinscale = 520.D5)
      parameter(irk4=0)
     
      real*8 psi1(ii,jj),psi2(ii,jj),M(ii,jj,4,4),a(ii,jj),b(ii,jj)
     & ,c(ii,jj),d(ii,jj)
     & ,psi1_x(4),psi2_x(4),u,v
     & ,x1(npoints),x2(npoints),y1(npoints),y2(npoints)
     & ,x_diff,y_diff
       
      integer x1_coord(npoints),x2_coord(npoints),y1_coord(npoints)
     & ,y2_coord(npoints)
     
     
      data dt/8640./ ! time step in seconds
c      data dt/4320./
c      data dt/2160./
c      data dt/1080./
       
       if(irk4.eq.0) then
      if(isolve.eq.0) then
      write(file_name,'("jet_2Dcubic__dt_",i5,".nc")')int(dt)
      elseif(isolve.eq.1)then
      write(file_name,'("jet_bicubic_dt_",i5,".nc")')int(dt)
      endif
      elseif(irk4.eq.1) then
      if(isolve.eq.0) then
      write(file_name,'("jet_rk4_2Dcubic_dt_",i5,".nc")')
     &      int(dt)
      elseif(isolve.eq.1)then
      write(file_name,'("jet_rk4_bicubic_dt_",i5,".nc")')
     &      int(dt)
      endif      
      endif
      
      call read_netcdf(qg_name,psi1,psi2,ii,jj,dfloat(0)
     & ,new_tim,step_tim)
      call create_NETCDF_file(file_name,npoints)
      
c      n = 0
c      iseed=102
c      do j=1,npoints_sqrt
c      do i=1,npoints_sqrt
c      n = n+1
c      x1(n) = ran1(iseed)*dfloat(ii)
c      y1(n) = ran1(iseed)*dfloat(jj)
c      enddo
c      enddo

      x1(1) = 100.
      y1(1) = 220.
      

      

      
      
      uscale = 1
      scale = basinscale/dfloat(jj)
      tscale = scale/uscale
      
      dt = dt/tscale
      
      k_traj = 1 ! save every k_traj days
      k_t = 0
      nrec = 0
      
      l_tot_day = int(86400./(dt*tscale)+0.001) ! time steps in one day
      max_time1 = max_time*l_tot_day ! total number of time steps
      
      print*,'total number of time steps =',max_time1
      
      time_day = 0.
      
      call cpu_time(start)
      
      do kk = 1,max_time1
      
      time_day = time_day + dt*tscale/86400.
      !print*,'time_day=',time_day
      l_day = l_day + 1
      
      if(l_day .eq. l_tot_day) then
       l_day = 0
       k_t = k_t + 1
      endif
      !print*, 'k_t,k_traj=',k_t,k_traj
      
      if(isolve.eq.0) then
      call cubic_coeff_x(ii,jj,psi1,a,b,c,d)
      elseif(isolve.eq.1) then
      call A_matrix(ii,jj,psi1,M)
      endif
      
      if (irk4.eq.1) then
      if(isolve.eq.0) then
      call cubic_coeff_x(ii,jj,psi1,a,b,c,d)
      elseif(isolve.eq.1) then
      call A_matrix(ii,jj,psi1,M)
      endif
      endif
c ----- ADVECT PARTICLES

        if(irk4.eq.0) then

        if(isolve.eq.0) then
        do n = 1,npoints
        call cubic_poly_x(ii,jj,x1(n),y1(n),a,b,c,d,psi1_x)
        call vel(ii,jj,psi1_x,a,b,c,d,x1(n),y1(n),u,v)
        enddo
        
        elseif(isolve.eq.1) then
        do n = 1,npoints
        call bicubic(ii,jj,M,x1(n),y1(n),u,v)
        enddo
        
        endif
        
        do n = 1,npoints
        x1(n) = x1(n) + dt*(u+U0)
        y1(n) = y1(n) + dt*v
        
            if(x1(n).lt.0.)then
              x1(n)=dfloat(ii)+x1(n)
            endif
            if(x1(n).gt.dfloat(ii))then
              x1(n)=x1(n)-dfloat(ii)
            endif
            if(y1(n).lt.0.)then
              y1(n)=dfloat(jj)+y1(n)
            endif
            if(y1(n).gt.dfloat(jj))then
              y1(n)=y1(n)-dfloat(jj)
            endif
        enddo
        
        elseif ( irk4.eq.1) then
        do n = 1,npoints
        if (isolve.eq.0) then
        call rk4_2dcubic(ii,jj,x1(n),y1(n),dt,U0,a,b,c
     & ,d,a,b,c,d,a,b,c,d
     & ,x_diff,y_diff)
        elseif (isolve.eq.1) then
        call rk4_bicubic(ii,jj,x1(n),y1(n),dt,U0,M,M,M
     & ,x_diff,y_diff)
        endif
        enddo
        do n = 1,npoints
        x1(n) = x1(n) + x_diff
        y1(n) = y1(n) + y_diff
        
            if(x1(n).lt.0.)then
              x1(n)=dfloat(ii)+x1(n)
            endif
            if(x1(n).gt.dfloat(ii))then
              x1(n)=x1(n)-dfloat(ii)
            endif
            if(y1(n).lt.0.)then
              y1(n)=dfloat(jj)+y1(n)
            endif
            if(y1(n).gt.dfloat(jj))then
              y1(n)=y1(n)-dfloat(jj)
            endif
        enddo
        endif
        if(k_t.eq.k_traj)then
           k_t=0
           nrec=nrec+1
c           write(11) x1_save,y1_save

        print*,'Time (day) = ',time_day

           call write_netcdf_file(file_name,npoints,x1,y1
     & ,time_day,nrec)

         endif    
      enddo
      
      call cpu_time(finish)
      
      print*, 'Time taken to run was',finish-start,'seconds'
      
      stop
      
      end program jet 
