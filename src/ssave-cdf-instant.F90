!  CVS: $Id: ssave-cdf.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     =================
      subroutine SSAVEINS
#include <def-undef.h>
!     =================
!     output in NETcdf format
!     written by liu hai long 2001 jun
use param_mod
use pconst_mod
use dyn_mod
use tracer_mod
use forc_mod
use work_mod
use buf_mod, only:t_cpl,s_cpl,u_cpl,v_cpl,dhdx,dhdy,qheat      
use domain
use gather_scatter
use distribution
use blocks
!use mct_mod
!use esmf
!use seq_flds_mod
!use seq_cdata_mod
!use seq_infodata_mod
!use seq_timemgr_mod

  implicit none


      logical       :: write_restart     ! restart now
      character (len=18) :: fname
      integer :: klevel, iiday, nwmf, mpierr
      integer(kind(1))   :: curr_ymd     ! Current date YYYYMMDD
      integer :: fh
      type (block) ::this_block
      integer :: subarray_size(2), subarray_start(2),psize(2)
      integer :: lsize(2),lstarts(2),local_rank(2)
      integer :: global_dims(2)
      integer(kind=4) :: file_type, buf_type
      integer(kind=MPI_OFFSET_KIND) :: offset

         number_day = iday +1 !LPF 20120816
         !number_month = month
      !if (mod(number_month,12)==0) then
         if ( number_day  > imd ) then
           if (mod(number_month,12)==0) then
            !DO LPF 20200401
            iyfm=int(number_month/12)+1
            !iyfm=iyfm+1 !20200401
            !number_month=1
            number_month=number_month+1
            !DO LPF 20200401 
            number_day=1
           else
             number_day = 1
             number_month = number_month + 1
           !if (mod(number_month,12)==0) then
           ! iyfm=iyfm+1
           end if
         end if
         nwmf= iyfm
         !DO LPF 20200401
         !mon0=number_month
         mon0=number_month-(iyfm-1)*12
         !DO LPF 20200401
         !mon0 = mod(number_month-1,12) + 1
         !mon0=number_month
         iiday = nnmonth(mon0) + iday + 1+ (nwmf-1)*365

!
     write_restart = .false.
     if ( rest_freq < 0) then
        if (number_day == 1) write_restart = .true.
     else if ( mod(iiday-1,rest_freq) == 0  ) then
        write_restart = .true.
     end if
     if(mytid==0) write(*,*)'in instant,iday=,imd=,mon0=, nwmf=, rest_freq= ',iday,imd,mon0,nwmf,rest_freq
     if(mytid==0) write(*,*)'in instant,write_restart', write_restart
!
!    if ( mod(iday,rest_freq) == 0 .or. iday == imd .or. iday == 10 .or. iday ==20) then
!    if ( mod(iday,rest_freq) == 0 .or. iday == 1 ) then
!    if ( mod(iiday,rest_freq) == 1  ) then

        if ( write_restart  ) then
!
!       if ( mytid == 0 ) then
         fname(1:8)='fort.22.'
         fname(13:13)='-'
         fname(16:16)='-'
         write(fname(14:15),'(i2.2)')mon0
         write(fname(9:12),'(i4.4)')nwmf
         write(fname(17:18),'(i2.2)') number_day
 
!        if ( number_day ==1 .and. mon0 < 12) then
!            write(fname(14:15),'(i2.2)')mon0+1
!        end if
!
!        if ( number_day ==1 .and. mon0 == 12) then
!            write(fname(14:15),'(i2.2)')mon0-11
!            write(fname(9:12),'(i4.4)')nwmf+1
!        end if
        if ( mytid == 0 ) then
         open (17, file="rpointer.ocn", form='formatted')
         write(17,'(a18)') fname
         close(17)
         write(*,*)'fname=',fname
#ifdef SEQ_IO
         open(22,file=trim(out_dir)//fname,form='unformatted') !jenny change
#endif
       end if
!
#ifdef SEQ_IO
        allocate(buffer(imt_global,jmt_global))
#else
       call mpi_barrier(mpi_comm_world,ierr)
#endif

#ifndef SEQ_IO

        this_block = get_block(mytid+1,1)
        psize(1) = (imt_global-1)/licom_BlockSizeX + 1
        psize(2) = (jmt_global-1)/licom_BlockSizeY + 1

        global_dims = (/imt_global, jmt_global/)
        offset = 0 
        subarray_size = (/BLCKX, BLCKY/)
        if (mytid/psize(1) == (psize(2)-1)) then
                subarray_size(2) = NJMT - (psize(2)-1)*BLCKY
        endif

        local_rank(1) = mod(mytid,psize(1)) 
        local_rank(2) = mytid / psize(1) 
        
        subarray_start(1) = local_rank(1) * BLCKX
        subarray_start(2) = local_rank(2) * BLCKY

        call MPI_Type_create_subarray(2, global_dims, subarray_size, &
        subarray_start, MPI_ORDER_FORTRAN, MPI_DBL, file_type, ierr)
        call MPI_Type_commit(file_type, ierr)

        !create buffer type
        lsize = (/imt,jmt/)
        lstarts = (/nghost,nghost/)
    
        call MPI_Type_create_subarray(2, lsize, subarray_size, &
        lstarts, MPI_ORDER_FORTRAN, MPI_DBL, buf_type, ierr)
        call MPI_Type_commit(buf_type, ierr)

        call MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_WRONLY + & 
        MPI_MODE_CREATE, MPI_INFO_NULL, fh, ierr)

        call MPI_File_set_view(fh, offset, MPI_DBL, file_type, &
        'native', MPI_INFO_NULL, ierr)
   
        where (vit(:,:,1,:) < 0.5D0) h0 = 0.0_r8 
        call MPI_File_write_all(fh, h0,1,buf_type,MPI_STATUS_IGNORE, ierr)

        where (viv < 0.5D0) u = 0.0_r8
        call MPI_File_write_all(fh, u, km, buf_type,MPI_STATUS_IGNORE, ierr)

        where (viv < 0.5D0) v = 0.0_r8
        call MPI_File_write_all(fh, v, km, buf_type,MPI_STATUS_IGNORE, ierr)

        where (vit < 0.5D0) at(:,:,:,1,:) = 0.0_r8
        where (vit < 0.5D0) at(:,:,:,2,:) = 0.0_r8
        call MPI_File_write_all(fh, at, km*2, buf_type,MPI_STATUS_IGNORE, ierr)

        where (vit < 0.5D0) ws(:,:,1:km,:) = 0.0_r8
        call MPI_File_write_all(fh, ws, km, buf_type,MPI_STATUS_IGNORE, ierr)

        where (viv(:,:,1,:) < 0.5D0) su = 0.0_r8
        call MPI_File_write_all(fh, su, 1, buf_type,MPI_STATUS_IGNORE, ierr)

        where (viv(:,:,1,:) < 0.5D0) sv = 0.0_r8
        call MPI_File_write_all(fh, sv, 1, buf_type,MPI_STATUS_IGNORE, ierr)

        where (vit(:,:,1,:) < 0.5D0) swv = 0.0_r8
        call MPI_File_write_all(fh, swv, 1, buf_type,MPI_STATUS_IGNORE, ierr)

        where (vit(:,:,1,:) < 0.5D0) sshf = 0.0_r8
        call MPI_File_write_all(fh, sshf, 1, buf_type,MPI_STATUS_IGNORE, ierr)

        where (vit(:,:,1,:) < 0.5D0) lthf = 0.0_r8
        call MPI_File_write_all(fh, lthf, 1, buf_type,MPI_STATUS_IGNORE, ierr)

        where (vit(:,:,1,:) < 0.5D0) fresh = 0.0_r8
        call MPI_File_write_all(fh, fresh, 1, buf_type,MPI_STATUS_IGNORE, ierr)

#ifdef COUP
        where (vit(:,:,1,:) < 0.5D0) t_cpl = 0.0_r8
        call MPI_File_write_all(fh, t_cpl, 1, buf_type,MPI_STATUS_IGNORE, ierr)

        where (vit(:,:,1,:) < 0.5D0) s_cpl = 0.0_r8
        call MPI_File_write_all(fh, s_cpl, 1, buf_type,MPI_STATUS_IGNORE, ierr)

        where (vit(:,:,1,:) < 0.5D0) u_cpl = 0.0_r8
        call MPI_File_write_all(fh, u_cpl, 1, buf_type,MPI_STATUS_IGNORE, ierr)

        where (vit(:,:,1,:) < 0.5D0) v_cpl = 0.0_r8
        call MPI_File_write_all(fh, v_cpl, 1, buf_type,MPI_STATUS_IGNORE, ierr)
        
        where (vit(:,:,1,:) < 0.5D0) dhdx = 0.0_r8
        call MPI_File_write_all(fh, dhdx, 1, buf_type,MPI_STATUS_IGNORE, ierr)

        where (vit(:,:,1,:) < 0.5D0) dhdy = 0.0_r8
        call MPI_File_write_all(fh, dhdy, 1, buf_type,MPI_STATUS_IGNORE, ierr)

        where (vit(:,:,1,:) < 0.5D0) qheat = 0.0_r8
        call MPI_File_write_all(fh, qheat, 1, buf_type,MPI_STATUS_IGNORE, ierr)
#endif

        call MPI_File_close(fh, ierr)
        call MPI_Type_free(file_type,ierr)
        call MPI_Type_free(buf_type,ierr)

      if(mytid==0) then
        !open(22,file=trim(out_dir)//fname,form='unformatted',access = 'append',convert='little_endian')
        open(22,file=trim(out_dir)//fname,form='unformatted',status='old',position='append',convert='little_endian')
        write(22) number_month, number_day
        close(22)
      end if

#else
! seq io
       where (vit(:,:,1,:) < 0.5D0) h0 = 0.0_r8
       call gather_global(buffer, h0, master_task,distrb_clinic)
!
       if (mytid==0) then
          WRITE (22)buffer
       end if
!
      where (viv < 0.5D0) u = 0.0_r8
      do klevel=1,km
          call gather_global(buffer, u(:,:,klevel,:), master_task,distrb_clinic)
          if (mytid==0) then
             WRITE (22)buffer
          end if
       end do
         !write(*,*)'finish U'
!        
      where (viv < 0.5D0) v = 0.0_r8
      do klevel=1,km
          call gather_global(buffer, v(:,:,klevel,:), master_task,distrb_clinic)
          if (mytid==0) then
             WRITE (22)buffer
          end if
       end do
         !write(*,*)'finish V'
!
      where (vit < 0.5D0) at(:,:,:,1,:) = 0.0_r8
      do klevel=1,km
         call gather_global(buffer, at(:,:,klevel,1,:), master_task,distrb_clinic)
         if (mytid==0) then
             WRITE (22)buffer
         end if
      end do
         !write(*,*)'finish at1'
!
      where (vit < 0.5D0) at(:,:,:,2,:) = 0.0_r8
      do klevel=1,km
         call gather_global(buffer, at(:,:,klevel,2,:), master_task,distrb_clinic)
         if (mytid==0) then
           WRITE (22)buffer
         end if
      end do
         !write(*,*)'finish at2'
!lhl20110728
      do klevel=1,km
         where (vit(:,:,klevel,:) < 0.5D0) ws(:,:,klevel,:) = 0.0_r8
         call gather_global(buffer, ws(:,:,klevel,:), master_task,distrb_clinic)
         if (mytid==0) then
             WRITE (22)buffer
         end if
      end do
         !write(*,*)'finish at2'
!        
      where (viv(:,:,1,:) < 0.5D0) su = 0.0_r8
      call gather_global(buffer, su, master_task,distrb_clinic)
      if (mytid==0) then
         WRITE (22)buffer
      end if
         !write(*,*)'finish su'

      where (viv(:,:,1,:) < 0.5D0) sv = 0.0_r8
      call gather_global(buffer, sv, master_task,distrb_clinic)
      if (mytid==0) then
         WRITE (22)buffer
      end if

         !write(*,*)'finish sv'
      where (vit(:,:,1,:) < 0.5D0) swv = 0.0_r8
     call gather_global(buffer, swv, master_task,distrb_clinic)
     if (mytid==0) then
         WRITE (22)buffer
     end if
         !write(*,*)'finish lwv'
      where (vit(:,:,1,:) < 0.5D0) sshf = 0.0_r8
     call gather_global(buffer, sshf, master_task,distrb_clinic)
     if (mytid==0) then
         WRITE (22)buffer
     end if
         !write(*,*)'finish sshf'
      where (vit(:,:,1,:) < 0.5D0) lthf = 0.0_r8
     call gather_global(buffer, lthf, master_task,distrb_clinic)
     if (mytid==0) then
         WRITE (22)buffer
     end if
         !write(*,*)'finish lthf'
      where (vit(:,:,1,:) < 0.5D0) fresh = 0.0_r8
     call gather_global(buffer, fresh, master_task,distrb_clinic)
     if (mytid==0) then
         WRITE (22)buffer
     end if
         !write(*,*)'finish fresh'
     if (mytid==0) then
         write(22) number_month, number_day
     end if
!
#ifdef COUP
      where (vit(:,:,1,:) < 0.5D0) t_cpl = 0.0_r8
     call gather_global(buffer, t_cpl, master_task,distrb_clinic)
     if (mytid==0) then
         WRITE (22)buffer
     end if
!         write(*,*)'finish t_cpl'
      where (vit(:,:,1,:) < 0.5D0) s_cpl = 0.0_r8
     call gather_global(buffer, s_cpl, master_task,distrb_clinic)
     if (mytid==0) then
         WRITE (22)buffer
     end if
!         write(*,*)'finish s_cpl'
      where (vit(:,:,1,:) < 0.5D0) u_cpl = 0.0_r8
     call gather_global(buffer, u_cpl, master_task,distrb_clinic)
     if (mytid==0) then
         WRITE (22)buffer
     end if
!         write(*,*)'finish u_cpl'
      where (vit(:,:,1,:) < 0.5D0) v_cpl = 0.0_r8
     call gather_global(buffer, v_cpl, master_task,distrb_clinic)
     if (mytid==0) then
         WRITE (22)buffer
     end if
!         write(*,*)'finish v_cpl'
      where (vit(:,:,1,:) < 0.5D0) dhdx = 0.0_r8
     call gather_global(buffer, dhdx, master_task,distrb_clinic)
     if (mytid==0) then
         WRITE (22)buffer
     end if
!         write(*,*)'finish dhdx'
      where (vit(:,:,1,:) < 0.5D0) dhdy = 0.0_r8
     call gather_global(buffer, dhdy, master_task,distrb_clinic)
     if (mytid==0) then
         WRITE (22)buffer
     end if
!         write(*,*)'finish dhdy'
      where (vit(:,:,1,:) < 0.5D0) qheat = 0.0_r8
     call gather_global(buffer, qheat, master_task,distrb_clinic)
     if (mytid==0) then
         WRITE (22)buffer
     end if
!         write(*,*)'finish q'
#endif
! end coup
     deallocate( buffer)
      if(mytid==0) then
        close(22)
      end if
#endif
!end seq io
  end if

      return
      end

#if (defined NETCDF) || (defined ALL)
      SUBROUTINE check_err (iret)
#include <netcdf.inc>
      INTEGER :: iret
      IF (iret /= NF_NOERR) THEN
         PRINT *, nf_strerror (iret)
         STOP
      END IF
      END SUBROUTINE check_err
#endif

