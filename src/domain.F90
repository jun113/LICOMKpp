!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module domain

!BOP
! !MODULE: domain
!
! !DESCRIPTION:
!  This module contains the model domain and routines for initializing
!  the domain.  It also initializes the decompositions and
!  distributions across processors/threads by calling relevent
!  routines in the block, distribution modules.
!
! !REVISION HISTORY:
!  SVN:$Id: domain.F90 12674 2008-10-31 22:21:32Z njn01 $

! !USES:


   use precision_mod
   use param_mod
   use msg_mod
   use POP_BlocksMod
   use POP_DistributionMod
   use POP_HaloMod
   use POP_GridHorzMod
   use control_mod
   !use control_licom_mod
   use LICOM_Error_mod
   use broadcast
   use blocks
   use distribution

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS

   public  :: init_domain_blocks ,&
              init_domain_distribution

! !PUBLIC DATA MEMBERS:

   integer (int_kind), public :: &
      nblocks_clinic   ,&! actual number of blocks on this processor
      nblocks_tropic     !   in each distribution

   integer (int_kind), dimension(:), pointer, public :: &
      blocks_clinic    ,&! block ids for local blocks in baroclinic dist
      blocks_tropic      ! block ids for local blocks in barotropic dist

   type (POP_distrb), public :: & !  block distribution info
      POP_distrbClinic    ,&! block distribution for baroclinic part
      POP_distrbTropic      ! block distribution for barotropic part
   !ALLOCATABLE
   integer(i4), dimension(:), pointer, public :: &
   ! integer(i4), dimension(:), allocatable, public :: &
      POP_distrbClinic_blockLocation, &
      POP_distrbClinic_blockLocalID,  &
      POP_distrbClinic_blockGlobalID

   type (distrb), public :: & !  block distribution info
      distrb_clinic    ,&! block distribution for baroclinic part
      distrb_tropic      ! block distribution for barotropic part

!------------------------------------------------------------
! Lets keep track of the land blocks for parallel IO reasons
!------------------------------------------------------------
   integer(int_kind), public :: &
      nblocks_land    ! acount number of land blocks assigned to processor

   integer(int_kind), dimension(:), pointer, public :: &
      blocks_land     ! blocks ids for land block

   type (distrb), public :: &! block distribution info for land
      distrb_land

   type (POP_halo), public :: &!  ghost cell update info
      POP_haloClinic         ,&! halo information for baroclinic part
      POP_haloTropic           ! halo information for barotropic part

   logical (log_kind), public :: &!
      ltripole_grid        ! flag to signal use of tripole grid


!EOP
!BOC
!-----------------------------------------------------------------------
!
!   module private variables - for the most part these appear as
!   module variables to facilitate sharing info between init_domain1
!   and init_domain2.
!
!-----------------------------------------------------------------------

    integer (i4) ::          &
       clinicDistributionMethod, &! method for distributing blocks
       tropicDistributionMethod   ! method for distributing blocks

    character (char_len) ::      &
       clinic_distribution_type, &! method to use for distributing
       tropic_distribution_type, &!    blocks in each case
       ew_boundary_type,         &! type of domain bndy in each logical
       ns_boundary_type           !    direction (ew is i, ns is j)

    integer (int_kind) :: &! decomposition info
       nprocs_clinic     ,&! num of processors in baroclinic dist
       nprocs_tropic       ! num of processors in barotropic dist

    logical, public :: profile_barrier
!EOC
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: init_domain_blocks
! !INTERFACE:

 subroutine init_domain_blocks

! !DESCRIPTION:
!  This routine reads in domain information and calls the routine
!  to set up the block decomposition.
!
! !REVISION HISTORY:
!  same as module
!EOP
!BOC
!----------------------------------------------------------------------
!
!  local variables
!
!----------------------------------------------------------------------

   integer (int_kind) :: &
      errorCode,         &! returned error code
      nml_error           ! namelist read error flag

!----------------------------------------------------------------------
!
!  input namelists
!
!----------------------------------------------------------------------

   namelist /domain_nml/ nprocs_clinic, nprocs_tropic, &
                         clinic_distribution_type, &
                         tropic_distribution_type, &
                         ew_boundary_type,         &
                         ns_boundary_type,         &
                         profile_barrier

!----------------------------------------------------------------------
!
!  read domain information from namelist input
!
!----------------------------------------------------------------------

   stdout  = 6
   errorCode = 0

   nprocs_clinic = -1
   nprocs_tropic = -1
   clinic_distribution_type = 'balanced'
   tropic_distribution_type = 'cartesian'
   ew_boundary_type = 'cyclic'
   ns_boundary_type = 'closed'
   profile_barrier  = .false.

   if (my_task == master_task) then
      open (11, file='ocn.parm', status='old',iostat=nml_error)
      if (nml_error /= 0) then
         nml_error = -1
      else
         nml_error =  1
      endif
      do while (nml_error > 0)
         read(11, nml=domain_nml,iostat=nml_error)
      end do
      if (nml_error == 0) close(11)
   endif

   call broadcast_scalar(nml_error, master_task)
   if (nml_error /= 0) then
       write (6,*) "can not read the namelist domain_nml"
       call flush(6)
       call MPI_ABORT(MPI_COMM_OCN,errorCode, ierr)
   endif

   call broadcast_scalar(nprocs_clinic,            master_task)
   call broadcast_scalar(nprocs_tropic,            master_task)
   call broadcast_scalar(clinic_distribution_type, master_task)
   call broadcast_scalar(tropic_distribution_type, master_task)
   call broadcast_scalar(ew_boundary_type,         master_task)
   call broadcast_scalar(ns_boundary_type,         master_task)
   call broadcast_scalar(profile_barrier,          master_task)

   select case (trim(clinic_distribution_type))
   case ('cartesian','Cartesian','CARTESIAN')
      clinicDistributionMethod = POP_distribMethodCartesian
   case ('balanced','Balanced','BALANCED')
      clinicDistributionMethod = POP_distribMethodRake
   case ('spacecurve','Spacecurve','SPACECURVE')
      clinicDistributionMethod = POP_distribMethodSpacecurve
   case default
      call LICOM_ErrorSet(errorCode, &
         'POP_DomainInit: unknown clinic distribution type')
      return
   end select

   select case (trim(tropic_distribution_type))
   case ('cartesian','Cartesian','CARTESIAN')
      tropicDistributionMethod = POP_distribMethodCartesian
   case ('balanced','Balanced','BALANCED')
      tropicDistributionMethod = POP_distribMethodRake
   case ('spacecurve','Spacecurve','SPACECURVE')
      tropicDistributionMethod = POP_distribMethodSpacecurve
   case default
      call LICOM_ErrorSet(errorCode, &
         'POP_DomainInit: unknown tropic distribution type')
      return
   end select

!----------------------------------------------------------------------
!
!  perform some basic checks on domain
!
!----------------------------------------------------------------------

   if (trim(ns_boundary_type) == 'tripole') then
      ltripole_grid = .true.
   else
      ltripole_grid = .false.
   endif

   if (imt_global < 1 .or. jmt_global < 1 .or. km < 1) then
      !***
      !*** domain size zero or negative
      !***
      call exit_licom(sigAbort,'Invalid domain: size < 1') ! no domain
   else if (ntra < 2) then
      !***
      !*** ntra must be at least 2 to hold temp,salinitiy
      !***
      call exit_licom(sigAbort,'Invalid tracer number: nt < 2')
   else if (nghost < 2) then
      !***
      !*** must have at least 2 layers of ghost cells
      !***
      call exit_licom(sigAbort,'Not enough ghost cells allocated')
   endif

!----------------------------------------------------------------------
!
!  compute block decomposition and details
!
!----------------------------------------------------------------------

   call create_blocks(imt_global, jmt_global, trim(ew_boundary_type), &
                                            trim(ns_boundary_type))
   call POP_BlocksCreate(imt_global, jmt_global,   &
                         trim(ew_boundary_type), &
                         trim(ns_boundary_type), errorCode)
   if (errorCode /= 0) then
      call LICOM_ErrorSet(errorCode, &
               'init_domain_blocks: error creating blocks')
      return
   endif

!----------------------------------------------------------------------
!
!  Now we need grid info before proceeding further
!  Print some domain information
!
!----------------------------------------------------------------------

   if (my_task == master_task) then
     write(stdout,'(a18)') 'Domain Information'
     write(stdout,'(a26,i6)') '  Horizontal domain: nx = ',imt_global
     write(stdout,'(a26,i6)') '                     ny = ',jmt_global
     write(stdout,'(a26,i6)') '  Vertical   domain: km = ',km
     write(stdout,'(a26,i6)') '  Number of tracers: ntra = ',ntra
     write(stdout,'(a26,i6)') '  Block size:  nx_block = ',nx_block
     write(stdout,'(a26,i6)') '               ny_block = ',ny_block
     write(stdout,'(a29,i6)') '  Processors for baroclinic: ', &
                                 nprocs_clinic
     write(stdout,'(a29,i6)') '  Processors for barotropic: ', &
                                 nprocs_tropic
     write(stdout,'(a31,a10)') '  Distribution for baroclinic: ', &
                                 trim(clinic_distribution_type)
     write(stdout,'(a31,a10)') '  Distribution for barotropic: ', &
                                 trim(tropic_distribution_type)
     write(stdout,'(a25,i2)') '  Number of ghost cells: ', nghost
   endif

!----------------------------------------------------------------------
!EOC

 end subroutine init_domain_blocks

!***********************************************************************
!BOP
! !IROUTINE: init_domain_distribution
! !INTERFACE:

 subroutine init_domain_distribution(KMTG)

! !DESCRIPTION:
!  This routine calls appropriate setup routines to distribute blocks
!  across processors and defines arrays with block ids for any local
!  blocks. Information about ghost cell update routines is also
!  initialized here through calls to the appropriate boundary routines.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), dimension(imt_global,jmt_global), intent(in) :: &
      KMTG             ! global KMT (topography) field

   integer (i4) :: &
      errorCode

!EOP
!BOC
!----------------------------------------------------------------------
!
!  local variables
!
!----------------------------------------------------------------------

   character (char_len) :: outstring

   integer (int_kind), parameter :: &
      max_work_unit=10   ! quantize the work into values from 1,max

   integer (int_kind) :: &
      i,j,k,n              ,&! dummy loop indices
      count1, count2       ,&! dummy counters
      work_unit            ,&! size of quantized work unit
      nblocks_tmp          ,&! temporary value of nblocks
      nblocks_tmp_clinic   ,&! num blocks on proc for clinic
      nblocks_max_clinic   ,&! max blocks on proc for clinic
      nblocks_tmp_tropic   ,&! num blocks on proc for tropic
      nblocks_max_tropic     ! max blocks on proc for tropic

   integer (int_kind), dimension(:), allocatable :: &
      nocn               ,&! number of ocean points per block
      work_per_block       ! number of work units per block

   type (block) :: &
      this_block         ! block information for current block

   integer (int_kind) :: jblock
!----------------------------------------------------------------------
!
!  estimate the amount of work per processor using the topography
!
!----------------------------------------------------------------------

   allocate(nocn(nblocks_tot))

   nocn = 0
   do n=1,nblocks_tot
      this_block = get_block(n,n)
      jblock = this_block%jblock
      do j=this_block%jb,this_block%je
         if (this_block%j_glob(j) > 0) then
            do i=this_block%ib,this_block%ie
               if (this_block%i_glob(i) > 0) then
                  if(KMTG(this_block%i_glob(i), this_block%j_glob(j)) > 0)  &
                                nocn(n) = nocn(n) + 1
               endif
            end do
         endif
      end do

      !*** with array syntax, we actually do work on non-ocean
      !*** points, so where the block is not completely land,
      !*** reset nocn to be the full size of the block
      ! wjl 20240224
      ! if ( my_task == master_task) write(*,*) "Points in Blocks", n, nocn(n)
      if (nocn(n) > 0) nocn(n) = nx_block*ny_block
   end do

   work_unit = maxval(nocn)/max_work_unit + 1

   !*** find number of work units per block

   allocate(work_per_block(nblocks_tot))

   where (nocn > 0)
     work_per_block = nocn/work_unit + 1
   elsewhere
     work_per_block = 0
   end where
   deallocate(nocn)
   if(my_task == master_task) then 
      write(stdout,'(a22,i6)') ' Active Ocean blocks: ',count(work_per_block > 0)
   endif



!----------------------------------------------------------------------
!
!  determine the distribution of blocks across processors
!
!----------------------------------------------------------------------

   distrb_tropic = create_distribution(tropic_distribution_type, &
                                       nprocs_tropic, work_per_block)

   distrb_clinic = create_distribution(clinic_distribution_type, &
                                       nprocs_clinic, work_per_block)

   POP_distrbClinic = POP_DistributionCreate(clinicDistributionMethod, &
                         nprocs_clinic, work_per_block, errorCode)

   if (errorCode /= 0) then
      call LICOM_ErrorSet(errorCode, &
         'POP_DomainInitDistrb: error creating clinic distrb')
      return
   endif

!  POP_distrbTropic = POP_DistributionCreate(tropicDistributionMethod, &
!                        nprocs_tropic, work_per_block, errorCode)

!  if (errorCode /= 0) then
!     call LICOM_ErrorSet(errorCode, &
!        'POP_DomainInitDistrb: error creating tropic distrb')
!     return
!  endif

   deallocate(work_per_block)

!----------------------------------------------------------------------
!
!  allocate and determine block id for any local blocks in each
!  distribution.
!
!----------------------------------------------------------------------

   call create_local_block_ids(blocks_clinic, distrb_clinic)
   call create_local_block_ids(blocks_tropic, distrb_tropic)

   if (associated(blocks_clinic)) then
      nblocks_clinic = size(blocks_clinic)
   else
      nblocks_clinic = 0
   endif
   nblocks_max_clinic = 0
   do n=0,distrb_clinic%nprocs - 1
     nblocks_tmp_clinic = nblocks_clinic
     call broadcast_scalar(nblocks_tmp_clinic, n)
     nblocks_max_clinic = max(nblocks_max_clinic,nblocks_tmp_clinic)
   end do

   if (nblocks_max_clinic > max_blocks_clinic) then
     write(outstring,*) 'clinic blocks exceed max: increase max to',&
                         nblocks_max_clinic
     call exit_LICOM(sigAbort,trim(outstring))
   else if (nblocks_max_clinic < max_blocks_clinic) then
     write(outstring,*) 'clinic blocks too large: decrease max to',&
                         nblocks_max_clinic
     if (my_task == master_task) write(stdout,*) trim(outstring)
   endif

!  if (my_task < distrb_tropic%nprocs .and. &
!      associated(blocks_tropic)) then
!    nblocks_tropic = size(blocks_tropic)
!  else
!    nblocks_tropic = 0
!  endif

!  nblocks_max_tropic = 0
!  do n=0,distrb_tropic%nprocs - 1
!    nblocks_tmp = nblocks_tropic
!    call broadcast_scalar(nblocks_tmp, n)
!    nblocks_max_tropic = max(nblocks_max_tropic,nblocks_tmp)
!  end do

!  if (nblocks_max_tropic > max_blocks_tropic) then
!    write(outstring,*) 'tropic blocks exceed max: increase max to',&
!                        nblocks_max_tropic
!    call exit_LICOM(sigAbort,trim(outstring))
!  else if (nblocks_max_tropic < max_blocks_tropic) then
!    write(outstring,*) 'tropic blocks too large: decrease max to',&
!                        nblocks_max_tropic
!    if (my_task == master_task) write(stdout,*) trim(outstring)
     !call exit_LICOM(sigAbort,trim(outstring))
!  endif

!----------------------------------------------------------------------
!
!  set up ghost cell updates for each distribution
!  Boundary types are cyclic, closed, or tripole
!
!----------------------------------------------------------------------

   POP_haloClinic = POP_HaloCreate(POP_distrbClinic,   &
                               trim(ns_boundary_type), &
                               trim(ew_boundary_type), &
                               imt_global, errorCode)

   if (errorCode /= 0) then
      call LICOM_ErrorSet(errorCode, &
         'POP_DomainInitDistrb: error creating clinic halo')
      return
   endif

!  POP_haloTropic = POP_HaloCreate(POP_distrbTropic,   &
!                              trim(ns_boundary_type), &
!                              trim(ew_boundary_type), &
!                              imt_global, errorCode)

!  if (errorCode /= 0) then
!     call LICOM_ErrorSet(errorCode, &
!        'POP_DomainInitDistrb: error creating tropic halo')
!     return
!  endif

!----------------------------------------------------------------------
!EOC

 end subroutine init_domain_distribution

!***********************************************************************

 end module domain

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
