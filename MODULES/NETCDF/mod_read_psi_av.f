        
        module MOD_read_psi_av
        contains
        subroutine psi_av(file_name,ii,jj,psi1_av,psi2_av)
      
      implicit none
      include 'netcdf.inc'
      
      character*(*) file_name
      integer ii,jj
      real*8 psi1_av(ii,jj), psi2_av(ii,jj), time_day
      
      integer ncid,retval,ndims,nlvls
      parameter(nlvls = 2,ndims=3)
      
      character*(*) psiav_name,t_name
      parameter(psiav_name = 'Time-averaged Stream Function')
      parameter(t_name = 'Time')
      
      integer t_varid, t_dimid,psiav_varid
      integer count(ndims), start(ndims)
      
        retval = nf_open(FILE_NAME, nf_nowrite, ncid)                   !! Open file in read mode
        if (retval .ne. nf_noerr) call handle_err(retval)
        
        retval = nf_inq_varid(ncid, Psiav_NAME, psiav_varid)                !! Request varid
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_inq_varid(ncid, t_NAME, t_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        
        retval = nf_inq_dimid(ncid, t_NAME, t_dimid)                !! Request dimension length
        if (retval .ne. nf_noerr) call handle_err(retval)
        
        retval = nf_get_vara_double(ncid, t_varid, 1, 1, time_day)   ! Read time first
        if (retval .ne. nf_noerr) call handle_err(retval)
        
        count(1) = ii
        count(2) = jj
        count(3) = 1
        start(1) = 1
        start(2) = 1
        start(3) = 1
        
        retval = nf_get_vara_double(ncid, psiav_varid, start, count,      !! Read Psi
     +        psi1_av)
        if (retval .ne. nf_noerr) call handle_err(retval)
        
        start(3) = 2
        retval = nf_get_vara_double(ncid, psiav_varid, start, count,
     +        psi2_av)
        if (retval .ne. nf_noerr) call handle_err(retval)        
     
        retval = nf_close(ncid)                                         !! Close file
        if (retval .ne. nf_noerr) call handle_err(retval)
        
        print*,'psi_av read'
        
       
        end subroutine psi_av
        
        subroutine handle_err(errcode)
        implicit none
        include 'netcdf.inc'
        integer errcode

        print *, 'Error: ', nf_strerror(errcode)
        stop 
        end subroutine
        
        end module MOD_read_psi_av
