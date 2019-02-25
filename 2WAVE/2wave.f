c
c Code that caclulates trajectories of particles due to two plane waves
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
        use mod_twowave
        use mod_rk4
c
      implicit none
      integer ips,itraj,isolve,max_av,ii,jj,max_time,max_time_rec
     & ,npoints,j,i,n,k_out,k_traj,k_o
     & ,l_tot_day,l_day,k_t,max_time1,k_phase,k_p,iseed,kk,nrec,nrec_t
     & ,npoints_sqrt,irk4
      real*8 basinscale,H1,H2,Rd,U1,V1,T_p,beta,visc
     & ,step
     & ,dt,dt2,dt6,time_day_halfstep,time_day
     & ,amp,k1,l1,k2,l2,cff1,cff2,cff3,cff4
     & ,time_old,time_half,time_new
     & ,xt,yt,u0_point,v0_point,u1_point,v1_point
     & ,u2_point,v2_point,u3_point,v3_point
     & ,x_diff,y_diff,x1_mean,y1_mean,read_tim

      parameter(ips=1,itraj=1,isolve=0,max_av=20,irk4=1
     & ,basinscale=520.D5,H1=1.D5,H2=3.D5
c     & ,ii=256,jj=ii
     & ,ii=512,jj=ii
c     & ,ii=1024,jj=ii
c     & ,ii=2048,jj=ii
     & ,Rd=40.D5,visc=20.D4,T_p=50.D0
     & ,U1=0.D0,V1=0.D0
c     & ,max_time=1*int(T_p)
     & ,max_time=100
     & ,max_time_rec=max_time/10
     & ,npoints_sqrt=500,npoints=npoints_sqrt**2)

      real*8 x_c(ii),y_c(jj)
     & ,psi_real(ii,jj), psi_imag(ii,jj)
     & ,u(ii,jj),v(ii,jj)


     & ,x0(npoints),y0(npoints)
     & ,x1(npoints),y1(npoints)
     & ,x1_traj(npoints),y1_traj(npoints)
     & ,xt1,yt1,phase

      real*8 disp_1(3,max_time_rec)
     & ,x1_save(npoints),y1_save(npoints)

      
      real*8 S1,S2,SS
     & ,uscale,scale,tscale,visc_nondim,beta_nondim
     & ,pi
     
      real*8 a_real(ii,jj),b_real(ii,jj)
     & ,c_real(ii,jj), d_real(ii,jj)
     & ,a_imag(ii,jj),b_imag(ii,jj)
     & ,c_imag(ii,jj), d_imag(ii,jj)
     & ,psi_x_real(4), psi_x_imag(4)
     & ,u0_real,u0_imag, v0_real,v0_imag
     & ,u1_real,u1_imag,v1_real,v1_imag
     & ,u2_real,u2_imag,v2_real,v2_imag
     & ,u3_real,u3_imag,v3_real,v3_imag
     & ,matrix_real(ii,jj,4,4)
     & ,matrix_imag(ii,jj,4,4)
     & ,factor
     & ,u_real,v_real
     
        integer dt_day
     
      real start,finish

     
      character*100 FILE_NAME 


c      data dt/86400./
c      data dt/43200./
c      data dt/21600./
c      data dt/10800./
      data dt/8640./
c      data dt/4320./
c      data dt/2160./
c      data dt/1080./
c       data dt/540./

      data pi/3.14159265358979323846D0/
      data beta/2.D-13/

    
      common /ONE/ x0,y0,x1,y1,x1_traj,y1_traj
      
      dt_day = int(dt)
      


      if(irk4.eq.0) then
       if (isolve.eq.0) then
         write(file_name,'("2wave_2Dcubic_dt_",i5,".nc")')dt_day


      elseif(isolve.eq.1) then
        write(file_name,'("2wave_bicubic_dt_",i5,".nc")')dt_day


      elseif(isolve.eq.2) then
         write (file_name, '("2wave_exact_dt_",i5,".nc")')dt_day
         print*,file_name

      endif
      
      elseif(irk4.eq.1) then
             if (isolve.eq.0) then
         write(file_name,'("2wave_rk4_2Dcubic_dt_",i5,".nc")')dt_day


      elseif(isolve.eq.1) then
        write(file_name,'("2wave_rk4_bicubic_dt_",i5,".nc")')dt_day


      elseif(isolve.eq.2) then
         write (file_name, '("2wave_rk4_exact_dt_",i5,".nc")')dt_day
         print*,file_name

      endif
      endif
      
      call create_NETCDF_file(file_name,npoints)



c
c--- nondimensionalization
c
      uscale=1.
      scale=basinscale/dfloat(jj)
      tscale=scale/uscale

c--- simple velocity fields
c
      do i=1,ii
         x_c(i)=dfloat(i-1)
      enddo

      do j=1,jj
         y_c(j)=dfloat(j-1)
      enddo
c
c--- amp=25 is threshold for linear interpolation with dt=43200
c
      amp=200.D0

      k1=-2.D0 !
      l1=+2.D0 ! number of periods per basin width
      k2=-2.D0 !
      l2=-2.D0 !

      cff1=2.D0*pi*k1/(dfloat(ii))
      cff2=2.D0*pi*l1/(dfloat(jj))
      cff3=2.D0*pi*k2/(dfloat(ii))
      cff4=2.D0*pi*l2/(dfloat(jj))
      
      factor= abs(cff2)+abs(cff4)
      

      do j=1,jj
         do i=1,ii
         
            psi_real(i,j) = amp*(dcos(cff1*x_c(i) + cff2*y_c(j))
     &                      + dcos(cff3*x_c(i) + cff4*y_c(j)))
     
           !psi_real(i,j) = psi_real(i,j)/factor
        
        
            
            u(i,j) = -amp*(cff2*dsin(cff1*x_c(i) + cff2*y_c(j))
     &                      +cff4*dsin(cff3*x_c(i) + cff4*y_c(j)))
            v(i,j) = amp*(cff1*dsin(cff1*x_c(i) + cff2*y_c(j))
     &                      +cff3*dsin(cff3*x_c(i) + cff4*y_c(j)))
            !u(i,j) = u(i,j)/factor
            !v(i,j) = v(i,j)/factor
            

    
         
         enddo
      enddo

      open(11,file = '2wave.dat',form = 'unformatted')
      write(11) psi_real
      close(11)
      
      open(12, file='2wave_u.dat',form='unformatted')
      write(12) u
      close(12)
      
      open(13,file='2wave_v.dat',form='unformatted')
      write(13) v
      close(13)
      

      



      
c
c---
c
c
c############################## LAGRANGIAN PARTICLES ##############################
c
      iseed=102
      step=dfloat(jj)/dfloat(npoints_sqrt)
      n=0
      do j=1,npoints_sqrt
         do i=1,npoints_sqrt
            n=n+1
c            x0(n)=step*dfloat(i-1)
c            y0(n)=step*dfloat(j-1)

            x0(n)=ran1(iseed)*dfloat(ii)
            y0(n)=ran1(iseed)*dfloat(jj)
         enddo
      enddo

c        x0(1) = 100.
c        y0(1) = 120.

      do n=1,npoints
         if(x0(n).lt.0.D0 .or. x0(n).gt.dfloat(ii))then
           write(*,*)n,x0(n)
           stop
         endif
         if(y0(n).lt.0.D0 .or. y0(n).gt.dfloat(jj))then
           write(*,*)n,y0(n)
           stop
         endif

         x1(n)=x0(n)
         y1(n)=y0(n)
         x1_traj(n) = x0(n)
         y1_traj(n) = y0(n)

      enddo
c
c--- counters
c
      dt=dt/tscale
      dt2=.5*dt
      dt6=dt/6.0

      k_out=10     ! diagnostics and output every k_out days
      k_traj=1     ! saving trajectories every k_traj days

c      k_phase=int(T_p)  ! randomization time parameter
      k_phase=75
      k_p=0
      phase=0.

      l_tot_day=int(86400./(dt*tscale)+0.001)  ! time steps in one day
      max_time1=max_time*l_tot_day             ! total number of time steps

      time_day=0.
      l_day=0
      k_o=0
      nrec=0

      k_t=0
      nrec_t=0
c
c#################### MAIN CYCLE
c
c      open(11,file='particles_traj_qg2.d',form='unformatted')
c      if(itraj.eq.1) then
c            call create_traj_file(file_name,npoints)
c      end if
        call cpu_time(start)
      do 1000 kk=1,max_time1
      
      if (isolve.eq.0) then
      call cubic_coeff_x(ii,jj,psi_real,a_real,b_real,c_real,d_real)

      else if (isolve.eq.1) then
      
      call A_matrix(ii,jj,psi_real,matrix_real)
      
      endif
      
      if (irk4.eq.1) then
            if (isolve.eq.0) then
      call cubic_coeff_x(ii,jj,psi_real,a_real,b_real,c_real,d_real)

      else if (isolve.eq.1) then
      
      call A_matrix(ii,jj,psi_real,matrix_real)
      
      endif
      endif


      

c         write(*,*)'kk=',kk
         time_old=2.*pi*time_day
         time_day=time_day+dt*tscale/86400.
         time_day_halfstep=time_day-0.5D0*dt*tscale/86400.

         l_day=l_day+1

         if(l_day.eq.l_tot_day)then
           l_day=0
           k_o=k_o+1
           k_p=k_p+1
           k_t=k_t+1
         endif

         time_new=2.*pi*time_day
         time_half=2.*pi*time_day_halfstep

c---insert
c         phase=0.
c         time_old=0.
c         time_new=0.
c         time_half=0.
c---end
c
c--- advect Lagrangian points: L1

        if(irk4.eq.1) then
c
         do n=1,npoints
         
            if (isolve.eq.0) then
            call rk4_2dcubic(ii,jj,x1(n),y1(n),dt,dfloat(0),a_real
     & ,b_real,c_real,d_real,a_real,b_real,c_real,d_real,a_real,b_real
     & ,c_real,d_real,x_diff,y_diff)
     
            elseif(isolve.eq.1) then
            
            call rk4_bicubic(ii,jj,x1(n),y1(n),dt,dfloat(0),matrix_real
     & ,matrix_real,matrix_real,x_diff,y_diff)
            elseif(isolve.eq.2) then
            
            call rk4_twowave(x1(n),y1(n),cff1,cff2,cff3,cff4
     &  ,amp,dt,x_diff,y_diff)
     
            endif
            
            x1(n)=x1(n)+x_diff
            y1(n)=y1(n)+y_diff
            x1_traj(n)=x1_traj(n)+x_diff
            y1_traj(n)=y1_traj(n)+y_diff

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
            !write(*,*),'at final step, xt,yt=',x1(n),y1(n)
         enddo
         
         elseif(irk4.eq.0) then
         
         do n=1,npoints
         
            if (isolve.eq.0) then
            call cubic_poly_x(ii,jj,x1(n),y1(n)
     & ,a_real,b_real,c_real,d_real,psi_x_real)

      call vel(ii,jj,psi_x_real
     & ,a_real,b_real,c_real,d_real
     & ,x1(n),y1(n),u_real,v_real)
     
            elseif(isolve.eq.1) then
            
            call bicubic(ii,jj,matrix_real,x1(n),y1(n),u_real,v_real)
     
            elseif(isolve.eq.2) then
            
            call twowave_velocity(amp,cff1,cff2,cff3,cff4
     &       ,x1(n),y1(n),u_real,v_real)
     
            endif
            
            x1(n)=x1(n)+dt*u_real
            y1(n)=y1(n)+dt*v_real

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
            !write(*,*),'at final step, xt,yt=',x1(n),y1(n)
         enddo
         
         endif
         
         
         
         

c
c--- Lagrangian diagnostics
c
         if(k_o.eq.k_out)then
           write(*,*)'Diagnostics at time=',time_day
           k_o=0
           nrec=nrec+1

           x1_mean=0.
           y1_mean=0.
           do n=1,npoints
              x1_mean=x1_mean+x1_traj(n)
              y1_mean=y1_mean+y1_traj(n)
           enddo
           x1_mean=x1_mean/dfloat(npoints)
           y1_mean=y1_mean/dfloat(npoints)

           write(*,'(A8,F6.1,A10,F6.1)')'x1_mean=',x1_mean
     &                               ,'; y1_mean=',y1_mean

           do i=1,3
              disp_1(i,nrec)=0.
           enddo
           do n=1,npoints
              disp_1(1,nrec)=disp_1(1,nrec)+(x1_traj(n)-x0(n))**2
              disp_1(2,nrec)=disp_1(2,nrec)
     &                          +(x1_traj(n)-x0(n))*(y1_traj(n)-y0(n))
              disp_1(3,nrec)=disp_1(3,nrec)+(y1_traj(n)-y0(n))**2

           enddo
           do i=1,3
              disp_1(i,nrec)=disp_1(i,nrec)/dfloat(npoints)
           enddo

           write(*,'(I5,A7,F7.0,A7,F7.0,A7,F7.0)')nrec
     &                                        ,'  D_xx=',disp_1(1,nrec)
     &                                        ,'; D_xy=',disp_1(2,nrec)
     &                                        ,'; D_yy=',disp_1(3,nrec)


c           open(22,file='particles_qg2.d',form='unformatted')
c           write(22) x1,y1
c           close(22)
         endif
c
c--- trajectories diagnostics
c
         if(itraj.eq.1 .and. k_t.eq.k_traj)then
           k_t=0
           nrec_t=nrec_t+1
           do n=1,npoints
              x1_save(n)=x1(n)
              y1_save(n)=y1(n)
           enddo
c           write(11) x1_save,y1_save

        print*,'Time (day) = ',time_day

           call write_netcdf_file(file_name,npoints,x1_save,y1_save
     & ,time_day,nrec_t)

         endif
c
c--- change phase
c
         if(k_p.eq.k_phase)then
           k_p=0
           phase=2.*pi*ran1(iseed)
         endif

1000  continue
c
c#################################################
c
c      close(11)
c      write(*,*)'File particles_traj_qg2.d was written...'

c      open(11,file='particles_qg2.d',form='unformatted')
c      write(11) x1,y1
c      close(11)
c      write(*,*)'File particles_qg2.d was written...'

      call cpu_time(finish)
      
      print*, 'Time taken to run was', finish-start ,'seconds'

      stop
      end
