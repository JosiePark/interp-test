      module MOD_traj_NETCDF

      contains
      subroutine create_traj_file(file_name,npoints)
c Creates a file called file_name, and takes the number of lagrangian particles as the input
      implicit none
      include 'netcdf.inc'
      character*(*) file_name
      integer npoints, ncid, retval,nlvls
      parameter(nlvls = 1)

      
c Define names of variables
      character*(*) t_name, loc_name, part_name,traj_name
      character*(*) coor_name,lvl_name
      parameter (t_name = 'Time',loc_name ='Location of Particle')
      parameter(part_name = 'Particle Number')
      parameter(traj_name = 'Trajectories')
      parameter(coor_name = 'Domain Coordinate')
      parameter(lvl_name = 'Layer')

      
c Define units
      character*(*) t_units, x_units, y_units, units
      parameter (t_units = 'Days',x_units = 'cms',y_units = 'cms')
      parameter (units = 'units')
      
c initilisae dimension ids
      integer t_dimid, loc_dimid, part_dimid, traj_dimid, coor_dimid
      integer lvl_dimid

c initialise variable ids
      integer t_varid, loc_varid,part_varid, traj_varid,coor_varid

      
c number of dimensions in space(2: x & y)
      integer spatial_dims
      parameter (spatial_dims = 2)
c number of dimensions in the particle location variables : space, particle number, time, layer
      integer ndims
      parameter (ndims = 4)
      integer dimid(ndims)

      
c create NETCDF file
      retval = nf_create(file_name,nf_clobber,ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
c Define dimensions in the NETCDF files
      retval = nf_def_dim(ncid,part_name,npoints,part_dimid) ! defines particle numbers
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_def_dim(ncid,t_name,nf_unlimited,t_dimid) ! time, allowed to go to infinity
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_def_dim(ncid,loc_name,spatial_dims,loc_dimid) ! spatial dimensions for the location
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_def_dim(ncid,lvl_name,nlvls,lvl_dimid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
c Define variables
      retval = nf_def_var(ncid,t_name,nf_double,1,t_dimid,t_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_put_att_text(ncid,t_varid,units,len(t_units)
     & ,t_units)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      dimid(1) = loc_dimid
      dimid(2) = part_dimid
      dimid(3) = lvl_dimid
      dimid(4) = t_dimid
      
      
c Define Trajectory matrix, which stores locations of all npoints lagrangian particles at time t

      retval = nf_def_var(ncid,traj_name,nf_double,ndims,dimid
     & ,traj_varid)
      if (retval .ne. nf_noerr) then 
      print *, 'traj not created'
      call handle_err(retval)
      endif
      
      retval = nf_def_var(ncid,coor_name,nf_int,ndims,dimid,
     & coor_varid)
      if (retval .ne. nf_noerr) call handle_err(retval) 
      
      retval = nf_put_att_text(ncid,traj_varid,units
     & ,len(x_units),x_units)
      if (retval .ne. nf_noerr) call handle_err(retval)
    
      
      retval = nf_enddef(ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      
     
      
      retval = nf_close(ncid)                                         !! Close file
      if (retval .ne. nf_noerr) call handle_err(retval)
   
      print *,'!!! Netcdf file created i.e. ', FILE_NAME, '!!!'
      print *,'!!!                                         !!!'

      end
      
c -------------------------------------------------------------------------------------- 
c -------------------- SUBROUTINE TO WRITE TO A NETCDF FILE ------------------------------
      subroutine write_traj_file(file_name,npoints,x1,y1,x2,y2
     & ,x1_coord
     & ,y1_coord,x2_coord,y2_coord,save_time
     & ,step_time)
c writes the new particle locations x,y to the NETCDF file at the save_time
      implicit none
      include 'netcdf.inc'
      
      character*(*) file_name
      integer npoints,retval,ncid,i,j, step_time,traj_varid,eddy_varid
     & ,x1_coord(npoints),y1_coord(npoints),coor_varid,coor_eddy_varid
     & ,x2_coord(npoints),y2_coord(npoints)
      real*8 x1(npoints),y1(npoints), save_time,x2(npoints),y2(npoints)
      
c Temporary trajectory variable
      real*8 traj_temp(2,npoints,2)
      
      integer coor(2,npoints,2)
      
      integer nlvls
      parameter(nlvls = 1)
      
c Define names
      character*(*) traj_name, t_name,coor_name,lvl_name
      parameter(traj_name = 'Trajectories')
      parameter(t_name = 'Time')
      parameter(coor_name = 'Domain Coordinate')
      parameter(lvl_name = 'Layer')

      
c number of dimensions in space(2: x & y)
      integer spatial_dims, t_varid, lvl_dimid
      parameter (spatial_dims = 2)
c number of dimensions in the particle location variables : location, coordinate, particle number, time
      integer ndims
      parameter (ndims = 4)
      integer dimid(ndims)
      integer start(ndims), count(ndims)
      
c open file in read mode
      retval = nf_open(file_name, nf_write,ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
c request the variable ids of the variables in the file
      retval = nf_inq_varid(ncid,traj_name,traj_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_inq_varid(ncid,t_name,t_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)  
      
      retval = nf_inq_varid(ncid,coor_name,coor_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      

c Write x, y to the temporary trajectory variable
      do j = 1,npoints
              traj_temp(1,j,1) = x1(j)
              traj_temp(2,j,1) = y1(j)
              traj_temp(1,j,2) = x2(j)
              traj_temp(2,j,2) = y2(j)
              coor(1,j,1) = x1_coord(j)
              coor(2,j,1) = y1_coord(j)
              coor(1,j,2) = x2_coord(j)
              coor(2,j,2) = y2_coord(j)
      end do
      
c Now tell the program where to write the data in the NETCDF file
      count(1) = spatial_dims
      count(2) = npoints
      count(3) = nlvls
      count(4) = 1
      start(1) = 1
      start(2) = 1
      start(3) = 1
      start(4) = step_time ! tells the program to start writing from step time, in the time variable, but writes in all columns in the location and particle numbers

c Update data
      retval = nf_put_vara_double(ncid,traj_varid,start,count,traj_temp)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_vara_int(ncid,coor_varid,start,count,coor)
      if (retval .ne. nf_noerr) call handle_err(retval)
      

c Update time
      retval = nf_put_vara_double(ncid, t_varid, step_time, 1,
     +       save_time)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_close(ncid)                                         !! Close NetCDF file                
      if (retval .ne. nf_noerr) call handle_err(retval)
        
      print *,'!!! Netcdf file modified i.e. ', FILE_NAME, '!!!'
      print *,'!!!                                         !!!'
      end 
      
      
c --------------------------------------------------------------
c SUBROUTINE TO READ TRAJECTORY FILE FROM READ TIME
        subroutine read_traj_file(file_name,npoints,x1,y1,x2,y2
     &   ,x1_coord,y1_coord,x2_coord,y2_coord 
     & ,read_tim,new_tim,step_tim)
     
        implicit none
        include 'netcdf.inc'
        
        integer npoints,x1_coord(npoints),y1_coord(npoints),step_tim
        integer x2_coord(npoints),y2_coord(npoints)
        character*(*) file_name
        real*8 x1(npoints),y1(npoints),read_tim,new_tim
        real*8 x2(npoints),y2(npoints)
        integer retval,ncid, traj_varid, t_varid, traj_dimid, t_dimid
     & ,coor_varid,coor_dimid,t_len,i
        real *8, dimension (:), allocatable :: tp_time
        
        character*(*) traj_name, t_name,coor_name,lvl_name
        parameter(traj_name = 'Trajectories')
        parameter(t_name = 'Time')
        parameter(coor_name = 'Domain Coordinate')
        parameter(lvl_name = 'Layer')
        
c number of dimensions in space(2: x & y)
      integer spatial_dims
      parameter (spatial_dims = 2)
      integer nlvls
      parameter (nlvls = 1)
c number of dimensions in the particle location variables : location, coordinate, particle number, time
      integer ndims
      parameter (ndims = 4)
      integer dimid(ndims)
      integer start(ndims), count(ndims)
        
        retval = nf_open(file_name,nf_nowrite,ncid) !open file in read mode
        if (retval .ne. nf_noerr) call handle_err(retval)
        
        retval = nf_inq_varid(ncid,traj_name,traj_varid)
        if (retval .ne. nf_noerr) then
        print *, 'trajectories not found'
        call handle_err(retval)
        endif
        retval = nf_inq_varid(ncid,t_name,t_varid)
        if (retval .ne. nf_noerr) then 
        print *, 'time not found'
        call handle_err(retval)
        endif
        retval = nf_inq_varid(ncid,coor_name,coor_varid)
        if (retval .ne. nf_noerr) then 
        print *, 'coordinate not found'
        call handle_err(retval)
        endif
        
        retval = nf_inq_dimid(ncid,traj_name,traj_dimid)
        if (retval .ne. nf_noerr) then 
        print *, 'traj dimid not found'
        call handle_err(retval)
        endif
        retval = nf_inq_dimid(ncid,t_name,t_dimid)
        if (retval .ne. nf_noerr)then
        print *, 'time dimid not found'
         call handle_err(retval)
         endif
        retval = nf_inq_dimid(ncid,coor_name,coor_dimid)
        if (retval .ne. nf_noerr) then 
        print *,'coordinate dimid not found'
        call handle_err(retval)
        endif
        
        allocate(tp_time(t_len))
        retval = nf_get_vara_double(ncid,t_varid,1,t_len,tp_time)
        if (retval .ne. nf_noerr) then 
        print *, 'time not read'
        call handle_err(retval)
        endif
        
        if(read_tim.ne.0.) then
          if(read_tim.gt.tp_time(t_len)) then
            print *,'Maximum time saved at t = ',tp_time(t_len),' Days'
            print *,'Starting simulation from this time .......' 
            new_tim = tp_time(t_len)
            step_tim = t_len
          else
            do i = 1, t_len
              if(tp_time(i).ge.read_tim) then
                step_tim = i
                new_tim = tp_time(i)
                print *,'Starting simulation from t = ',new_tim,' Days'
                stop
              end if
            end do
          end if
        else
          new_tim = tp_time(t_len)
          step_tim = t_len
          print *,'Starting simulation from t = ',new_tim,' Days'
        end if  
        
      count(1) = 1
      count(2) = npoints
      count(3) = 1
      count(4) = 1
      start(1) = 1
      start(2) = 1
      start(3) = 1
      start(4) = step_tim ! tells the program to start writing from step time, in the time variable, but writes in all columns in the location and particle numbers
      
      retval = nf_get_vara_double(ncid,traj_varid,start,count,x1)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_get_vara_int(ncid,coor_varid,start,count,x1_coord)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      
      start(1) = 2
      retval = nf_get_vara_double(ncid,traj_varid,start,count,y1)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_get_vara_int(ncid,coor_varid,start,count,y1_coord)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      start(1) = 1
      start(3) = 2
      
      retval = nf_get_vara_double(ncid,traj_varid,start,count,x2)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_get_vara_int(ncid,coor_varid,start,count,x2_coord)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      start(1) = 2
      
      retval = nf_get_vara_double(ncid,traj_varid,start,count,y2)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_get_vara_int(ncid,coor_varid,start,count,y2_coord)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      
        retval = nf_close(ncid)                                         !! Close file
        if (retval .ne. nf_noerr) call handle_err(retval)
        
        deallocate (tp_time)
        
       end
c -----------------------------------------------------------------------
c -----------------------------------------------------------------------
       
              subroutine handle_err(errcode)
        implicit none
        include 'netcdf.inc'
        integer errcode

        print *, 'Error: ', nf_strerror(errcode)
        stop 2
       end
       
c ---------------------------------------------------------------------------
c ---------------------------------------------------------------------------
c SUBROUTINE THAT CREATES NETCDF FILES THAT STORE THE LOCATIONS OF LAGRANGIAN PARTICLES AT DIFFERENT RELEASE TIMES

      subroutine create_release_file(file_name,npoints,release_no)
c npoints is the number of lagrangian particles per release
c release_no is the number of releases

      implicit none
      include 'netcdf.inc'
      
      character*(*) file_name
      
      integer npoints,release_no,retval,ncid
      
c Define names of variables
      character*(*) t_name, loc_name, part_name,traj_name
      character*(*) coor_name,lvl_name,rel_name,rec_name
      parameter (t_name = 'Time',loc_name ='Location of Particle')
      parameter(part_name = 'Particle Number')
      parameter(traj_name = 'Trajectories')
      parameter(coor_name = 'Domain Coordinate')
      parameter(lvl_name = 'Layer')
      parameter(rel_name = 'Release Number')
      parameter(rec_name = 'Record Time')
      
c Define units
      character*(*) t_units, x_units, y_units, units
      parameter (t_units = 'Days',x_units = 'grid point'
     &      ,y_units = 'grid point')
      parameter (units = 'units')
c initilisae dimension ids
      integer t_dimid, loc_dimid, part_dimid, traj_dimid, coor_dimid
      integer lvl_dimid, rel_dimid,rec_dimid

c initialise variable ids
      integer t_varid,loc_varid,part_varid,traj_varid,coor_varid
     & ,rel_varid,rec_varid

      
c number of dimensions in space(2: x & y)
      integer spatial_dims
      parameter (spatial_dims = 2)
c number of dimensions in the particle location variables : space, particle number, realisation number, time, layer
      integer ndims
      parameter (ndims = 5)
      integer dimid(ndims)
c number of layers
      integer nlvls
      parameter (nlvls=2)
c dimension of the time array 
      integer tdimid(2)
      

      
c create NETCDF file
      retval = nf_create(file_name,nf_clobber,ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
c Define dimensions in the NETCDF files
      retval = nf_def_dim(ncid,part_name,npoints,part_dimid) ! defines particle numbers
      if (retval .ne. nf_noerr) call handle_err(retval)
      
c      retval = nf_def_dim(ncid,t_name,nf_unlimited,t_dimid) ! time, allowed to go to infinity
c      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_def_dim(ncid,loc_name,spatial_dims,loc_dimid) ! spatial dimensions for the location
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_def_dim(ncid,lvl_name,nlvls,lvl_dimid)
      if (retval .ne. nf_noerr) call handle_err(retval)  
      
      retval = nf_def_dim(ncid,rel_name,release_no,rel_dimid)
      if (retval .ne. nf_noerr) call handle_err(retval) 
      
      retval = nf_def_dim(ncid,rec_name,nf_unlimited,rec_dimid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
c Define variables
c Time -->
c Need time to be size of number of releases so that it stores 
c time record for each realisation
      
      dimid(1) = loc_dimid
      dimid(2) = part_dimid
      dimid(3) = rel_dimid
      dimid(4) = lvl_dimid
      dimid(5) = rec_dimid
      
      tdimid(1) = rel_dimid
      tdimid(2) = rec_dimid
      
      
c Define Trajectory matrix, which stores locations of all npoints lagrangian particles at time t

      retval = nf_def_var(ncid,traj_name,nf_double,ndims,dimid
     & ,traj_varid)
      if (retval .ne. nf_noerr) then 
      print *, 'traj not created'
      call handle_err(retval)
      endif
      
      retval = nf_def_var(ncid,coor_name,nf_int,ndims,dimid,
     & coor_varid)
      if (retval .ne. nf_noerr) call handle_err(retval) 
      
      retval = nf_put_att_text(ncid,traj_varid,units
     & ,len(x_units),x_units)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_def_var(ncid,t_name,nf_double,2,tdimid,t_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_put_att_text(ncid,t_varid,units,len(t_units)
     & ,t_units)
      if (retval .ne. nf_noerr) call handle_err(retval)
    
      
      retval = nf_enddef(ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      
      retval = nf_close(ncid)                                         !! Close file
      if (retval .ne. nf_noerr) call handle_err(retval)
   
      print *,'!!! Netcdf file created i.e. ', FILE_NAME, '!!!'
      print *,'!!!                                         !!!'

      end subroutine create_release_file  
      
c ----- WRITE TO THE RELEASE FILE -------------------
c -----------------------------------------------------
c ---------------------------------------
      
      subroutine write_release_file(file_name,npoints,release_no,x1,y1
     & ,x2,y2,x1_coord,y1_coord,x2_coord,y2_coord
     & ,release,save_time,step_time)

c release_no is the number of realisations
c release is the realisation you wish to write to
      
      implicit none
      include 'netcdf.inc'
      
      character*(*) file_name
      integer npoints,retval,ncid,i,j, step_time,traj_varid
     & ,x1_coord(npoints),y1_coord(npoints),coor_varid
     & ,x2_coord(npoints),y2_coord(npoints)
      integer release,release_no
      real*8 x1(npoints),y1(npoints), save_time,x2(npoints),y2(npoints)
      
c Temporary trajectory variable
      real*8 traj_temp(2,npoints,2)
      
      integer coor(2,npoints,2)
      
      integer nlvls
      parameter(nlvls = 2)
      
c Define names
      character*(*) traj_name, t_name,coor_name,lvl_name, rel_name
     & ,rec_name      
      parameter(traj_name = 'Trajectories')
      parameter(t_name = 'Time')
      parameter(coor_name = 'Domain Coordinate')
      parameter(lvl_name = 'Layer')
      parameter(rel_name = 'Release Number') 
      parameter(rec_name = 'Record Time')

      
c number of dimensions in space(2: x & y)
      integer spatial_dims, t_varid, lvl_dimid
      parameter (spatial_dims = 2)
c number of dimensions in the particle location variables : location, coordinate, particle number, time
      integer ndims
      parameter (ndims = 5)
      integer dimid(ndims)
      integer start(ndims), count(ndims)
      
      integer tdimid(2)
      integer start1(2),count1(2)
      
c open file in read mode
      retval = nf_open(file_name, nf_write,ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
c request the variable ids of the variables in the file
      retval = nf_inq_varid(ncid,traj_name,traj_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_inq_varid(ncid,t_name,t_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)  
      
      retval = nf_inq_varid(ncid,coor_name,coor_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      

c Write x, y to the temporary trajectory variable
      do j = 1,npoints
              traj_temp(1,j,1) = x1(j)
              traj_temp(2,j,1) = y1(j)
              traj_temp(1,j,2) = x2(j)
              traj_temp(2,j,2) = y2(j)
              coor(1,j,1) = x1_coord(j)
              coor(2,j,1) = y1_coord(j)
              coor(1,j,2) = x2_coord(j)
              coor(2,j,2) = y2_coord(j)
      end do
      
c Now tell the program where to write the data in the NETCDF file
      count(1) = spatial_dims
      count(2) = npoints
      count(3) = 1
      count(4) = nlvls
      count(5) = 1
      start(1) = 1
      start(2) = 1
      start(3) = release ! tells the program to start writing at the release number
      start(4) = 1
      start(5) = step_time ! tells the program to start writing from step time, in the time variable, but writes in all columns in the location and particle numbers

c Update data
      retval = nf_put_vara_double(ncid,traj_varid,start,count,traj_temp)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_vara_int(ncid,coor_varid,start,count,coor)
      if (retval .ne. nf_noerr) call handle_err(retval)
      

c Update time

      count1(1) = 1
      count1(2) = 1
      start1(1) = release
      start1(2) = step_time
      
      retval = nf_put_vara_double(ncid, t_varid, start1, count1,
     +       save_time)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_close(ncid)                                         !! Close NetCDF file                
      if (retval .ne. nf_noerr) call handle_err(retval)
        
      print *,'!!! Netcdf file modified i.e. ', FILE_NAME, '!!!'
      
      
      end subroutine write_release_file
      
      subroutine write_NETCDF_file(file_name,npoints,x,y,save_time
     & ,step_time)
c writes the new particle locations x,y to the NETCDF file at the save_time
      implicit none
      include 'netcdf.inc'
      
      character*(*) file_name
      integer npoints,retval,ncid,i, step_time,traj_varid
      real*8 x(npoints),y(npoints), save_time
c Temporary trajectory variable
      real*8 traj_temp(2,npoints)
      
c Define names
      character*(*) traj_name, t_name
      parameter(traj_name = 'Trajectories')
      parameter(t_name = 'Time')
      
c number of dimensions in space(2: x & y)
      integer spatial_dims, t_varid
      parameter (spatial_dims = 2)
c number of dimensions in the particle location variables : space, particle number, time
      integer ndims
      parameter (ndims = 3)
      integer dimid(ndims)
      integer start(ndims), count(ndims)
      
c open file in read mode
      retval = nf_open(file_name, nf_write,ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
c request the variable ids of the variables in the file
      retval = nf_inq_varid(ncid,traj_name,traj_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)

      
      retval = nf_inq_varid(ncid,t_name,t_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)  

c Write x, y to the temporary trajectory variable
      do i = 1,npoints
            traj_temp(1,i) = x(i)
            traj_temp(2,i) = y(i)
      end do
      
c Now tell the program where to write the data in the NETCDF file
      count(1) = spatial_dims
      count(2) = npoints
      count(3) = 1
      start(1) = 1
      start(2) = 1
      start(3) = step_time ! tells the program to start writing from step time, in the time variable, but writes in all columns in the location and particle numbers

c Update data
      retval = nf_put_vara_double(ncid,traj_varid,start,count,traj_temp)
      if (retval .ne. nf_noerr) call handle_err(retval)


c Update time
      retval = nf_put_vara_double(ncid, t_varid, step_time, 1,
     +       save_time)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_close(ncid)                                         !! Close NetCDF file                
      if (retval .ne. nf_noerr) call handle_err(retval)
        
      print *,'!!! Netcdf file modified i.e. ', FILE_NAME, '!!!'
      print *,'!!!                                         !!!'
      end
      
      subroutine create_NETCDF_file(file_name,npoints)
c Creates a file called file_name, and takes the number of lagrangian particles as the input
      implicit none
      include 'netcdf.inc'
      character*(*) file_name
      integer npoints, ncid, retval

      
c Define names of variables
      character*(*) t_name, loc_name, part_name,traj_name
      parameter (t_name = 'Time',loc_name ='Location of Particle')
      parameter(part_name = 'Particle Number')
      parameter(traj_name = 'Trajectories')
      
c Define units
      character*(*) t_units, x_units, y_units, units
      parameter (t_units = 'Days',x_units = 'cms',y_units = 'cms')
      parameter (units = 'units')
      
c initilisae dimension ids
      integer t_dimid, loc_dimid, part_dimid, traj_dimid
c initialise variable ids
      integer t_varid, loc_varid,part_varid, traj_varid
      
c number of dimensions in space(2: x & y)
      integer spatial_dims
      parameter (spatial_dims = 2)
c number of dimensions in the particle location variables : space, particle number, time
      integer ndims
      parameter (ndims = 3)
      integer dimid(ndims)
      
c create NETCDF file
      retval = nf_create(file_name,nf_clobber,ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
c Define dimensions in the NETCDF files
      retval = nf_def_dim(ncid,part_name,npoints,part_dimid) ! defines particle numbers
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_def_dim(ncid,t_name,nf_unlimited,t_dimid) ! time, allowed to go to infinity
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_def_dim(ncid,loc_name,spatial_dims,loc_dimid) ! spatial dimensions for the location
      if (retval .ne. nf_noerr) call handle_err(retval)
      
c Define variables
      retval = nf_def_var(ncid,t_name,nf_double,1,t_dimid,t_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
c      retval = nf_def_var(ncid,loc_name,nf_double,1,loc_dimid,loc_varid)
c      if (retval .ne. nf_noerr) call handle_err(retval)
      
c Assign units to varibles
c      retval = nf_put_att_text(ncid,loc_varid,units,len(x_units),
c     & ,x_units)
c      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_put_att_text(ncid,t_varid,units,len(t_units)
     & ,t_units)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      dimid(1) = loc_dimid
      dimid(2) = part_dimid
      dimid(3) = t_dimid
      
c Define Trajectory matrix, which stores locations of all npoints lagrangian particles at time t

      retval = nf_def_var(ncid,traj_name,nf_double,ndims,dimid
     & ,traj_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_put_att_text(ncid,traj_varid,units
     & ,len(x_units),x_units)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_enddef(ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_close(ncid)                                         !! Close file
      if (retval .ne. nf_noerr) call handle_err(retval)
   
      print *,'!!! Netcdf file created i.e. ', FILE_NAME, '!!!'
      print *,'!!!                                         !!!'

      end
      
       subroutine read_netcdf(FILE_NAME,psi1,psi2,Num_x,Num_y,read_tim,
     +  new_tim, step_tim)
       
        implicit none
        include 'netcdf.inc'
        
        integer Num_x, Num_y, step_tim
        character*(*) FILE_NAME
        real*8 psi1(Num_x,Num_y), psi2(Num_x,Num_y), read_tim, new_tim
        
        integer NDIMS, NLVLS, ETYPE, i, j
        parameter (NDIMS = 4, ETYPE = 4, NLVLS=2)
        integer ncid, retval, dimids(NDIMS), dimid(2),
     +   x_dimid, y_dimid, rec_dimid, lvl_dimid, Et_dimid, 
     +   start(NDIMS), count(NDIMS), start1(2), count1(2)
     
        character*(*) Psi_NAME, REC_NAME                                       
        parameter(Psi_NAME='Stream Function',REC_NAME='Time')
        integer psi_varid, t_varid, t_len
        real *8 tp_psi(Num_x,Num_y,NLVLS)
        real *8, dimension (:), allocatable :: tp_time
        
        retval = nf_open(FILE_NAME, nf_nowrite, ncid)                   !! Open file in read mode
        if (retval .ne. nf_noerr) call handle_err(retval)
        
        retval = nf_inq_varid(ncid, Psi_NAME, psi_varid)                !! Request varid
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_inq_varid(ncid, REC_NAME, t_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        
        retval = nf_inq_dimid(ncid, REC_NAME, rec_dimid)                !! Request dimension length
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_inq_dimlen(ncid, rec_dimid, t_len)
        if (retval .ne. nf_noerr) call handle_err(retval)
        
        allocate(tp_time(t_len))
        retval = nf_get_vara_double(ncid, t_varid, 1, t_len, tp_time)   ! Read time first
        if (retval .ne. nf_noerr) call handle_err(retval)
        
c        print *, t_len, tp_time(t_len)
        
        if(read_tim.ne.0.) then
          if(read_tim.gt.tp_time(t_len)) then
            print *,'Maximum time saved at t = ',tp_time(t_len),' Days'
            print *,'Starting simulation from this time .......' 
            new_tim = tp_time(t_len)
            step_tim = t_len
          else
            do i = 1, t_len
              if(tp_time(i).ge.read_tim) then
                step_tim = i
                new_tim = tp_time(i)
                print *,'Starting simulation from t = ',new_tim,' Days'
                exit
              end if
            end do
          end if
        else
          new_tim = tp_time(t_len)
          step_tim = t_len
          print *,'Starting simulation from t = ',new_tim,' Days'
        end if
        
        count(1) = Num_x                                                                     
        count(2) = Num_y
        count(3) = 1                                                    !! Start, count tell the program where to  
        count(4) = 1                                                    !! write data in NetCDF file 
        start(1) = 1
        start(2) = 1
        start(3) = 1
        start(4) = step_tim
        
        retval = nf_get_vara_double(ncid, psi_varid, start, count,      !! Read Psi
     +        psi1)
        if (retval .ne. nf_noerr) call handle_err(retval)
        
        start(3) = 2
        retval = nf_get_vara_double(ncid, psi_varid, start, count,
     +        psi2)
        if (retval .ne. nf_noerr) call handle_err(retval)        
     
        retval = nf_close(ncid)                                         !! Close file
        if (retval .ne. nf_noerr) call handle_err(retval)
        
        deallocate (tp_time)
        
       end
       
      end module MOD_traj_NETCDF
      

