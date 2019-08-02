! This module was core to my PhD work on Radiative transfer modelling of
! stellar environments. This module defines the 3D astrophysical environment
! in which the radiative transfer simulation occurs. 

module model_definition_module
use amr_module
use amrray_module
use rtglobal_module
use constants_module
!use camera_module
use dust_module
use lines_module
use stars_module
use quantum_module
use montecarlo_module
use namelist_module

!------------------------------------------------------------------------
! Here you can define your own variables, arrays etc.
!------------------------------------------------------------------------

integer :: userdef_nx,userdef_ny,userdef_nz
doubleprecision :: userdef_diskInR,userdef_diskOutR,userdef_h0,userdef_r0
doubleprecision :: userdef_beta,userdef_alpha,userdef_diskMass
...
contains

!------------------------------------------------------------------------
! This subroutine allows you to specify defaults values for your
! variables.
!
! WARNING: In this subroutine you are not allowed to use write(stdo,*),
!          because the stdo is not yet set. Reason: The defaults are
!          set before the command-line options are interpreted (so that
!          the defaults can be overwritten by command-line options),
!          but the stdo depends on whether the user calls RADMC-3D as
!          a child or not, which is given by the command-line options.
!------------------------------------------------------------------------


!subroutine userdef_paramlimits()
!#  implicit none
!#
!#  open(unit=1,file="params.py")
!#
!#  write(1,*), "distance = [100,1500]"
!  write(1,*), "diskInR = [5,10]"


 ! close(1)

!end subroutine userdef_paramlimits

subroutine userdef_defaults()
  implicit none
  !
  userdef_nx = 80
  userdef_ny = 60
  userdef_nz = 1
!  userdef_diskInR = 9.d0 * au
! userdef_diskOutR = 500.d0 *au
  userdef_h0 = 17.d0 !astronomical units
  userdef_r0 = 100.d0 !astronomical units
  userdef_beta = 1.02d0
  userdef_alpha = 2.4d0
  userdef_diskMass = 9.d-6/2.d0
  userdef_sfx = 1.25d0
  userdef_sfy = 1.14d0
  userdef_sf = 1.14d0
  userdef_starLum = 10000.0
  userdef_starTemp = 35000.0
  userdef_starLum_2 = 0.0
  userdef_starTemp_2 = 0.0
  userdef_amin = 0.05
  userdef_amax = 10.0
  userdef_apow = 3.5
end subroutine userdef_defaults

....

!------------------------------------------------------------------------
! Here you can interpret your own command-line options.
!------------------------------------------------------------------------
subroutine userdef_commandline(buffer,numarg,iarg,fromstdi,gotit)
  implicit none
  character*100 :: buffer
  integer :: iarg,numarg
  logical :: gotit,fromstdi
  !
  !

  if(buffer(1:8).eq.'diskMass') then
      if(iarg.gt.numarg) then
          write(stdo,*) 'Please enter an appropriate disk mass'
          stop
      endif
      call ggetarg(iarg,buffer,fromstdi)
      iarg = iarg+1
      read(buffer,*) userdef_diskMass
      gotit = .true.


  elseif(buffer(1:4).eq.'amax') then
      if(iarg.gt.numarg) then
          write(stdo,*) 'Please enter an appropriate amax'
          stop
      endif
      call ggetarg(iarg,buffer,fromstdi)
      iarg = iarg+1
      read(buffer,*) userdef_amax
      gotit = .true.


....



!------------------------------------------------------------------------
! Define (base) grid setup here,
! read your own frequency grid or set up your own stellar sources.
! No example given here, because it would interfere with basic operations.
!------------------------------------------------------------------------
subroutine userdef_prep_model()
  implicit none
  doubleprecision, allocatable :: xi(:),yi(:),zi(:)
  doubleprecision :: dr1, dr2
  integer :: i,ntot
  !
  ! Tell the code which grid type is used (AMR)
  !
  igrid_type      = 0             ! Regular grid
  igrid_coord     = 100             ! Spherical
  !
  ! Communicate these values to the AMR module
  !
  amr_style       = igrid_type
  amr_coordsystem = igrid_coord
  !
  ! Set up the regular grid coordinates in temporary arrays
  !
  allocate(xi(userdef_nx+1),yi(userdef_ny+1),zi(userdef_nz+1))
  !xi = (/0.1*au,10.0*au,50.0*au,400.0*au/)
  !yi = (/0.d0,0.3*pi,0.5*pi/)
  zi = (/0.D0,twopi/)
  ! calculat xi with logspace
  !!


  dr1 = (userdef_diskOutR*au-userdef_diskInR*au)*(userdef_sfx-1.d0)/(userdef_sfx**userdef_nx-1.d0)
  do i=1,userdef_nx+1
      xi(i) = userdef_diskInR*au + dr1 * (userdef_sfx**(i-1.d0)-1.d0)/(userdef_sfx-1.d0)
  enddo
  !!! mirror symmetry 0.5*pi

  dr2 = (0.5d0*pi)*(userdef_sfy-1.d0)/(userdef_sfy**userdef_ny-1.d0)

  do i=1,userdef_ny+1
    yi(i) = dr2 * (userdef_sfy**(i-1.d0)-1.d0)/(userdef_sfy-1.d0)

  enddo
  yi = 0.5d0*pi - yi(userdef_ny+1:1:-1)
  yi(1) = 0.d0

  !!0.5D0*pi as we consider mirror symmetry of north to south hemisphere


  ! Compute the total number of cells, so that we can tell the AMR
  ! initializer routine below how much memory it has to reserve. Note
  ! that if you want to refine lateron (while building up the model
  ! for instance) you must make sure here to reserve enough extra
  ! space for all the extra branches and leafs to fit in. For models
  ! with deep refinement levels the nr of branches can be up to (at
  ! most) 15% larger than the nr of leafs (0.125+0.125^2+0.125^3....
  ! = 0.142857).
  !
  ntot = userdef_nx*userdef_ny*userdef_nz
  !
  ! Now set up the grid, for now just a regular grid
  !
  call amr_initialize(.true.,.true.,.false.,            &
       userdef_nx,userdef_ny,userdef_nz,xi,yi,zi,1,1,  &
       ntot,ntot,.false.,.false.,.false.)
  !
  ! Set up the index lists so that we know how to connect each AMR grid
  ! cell to a data memory location.
  !
  call amr_compute_list_all()
  !
  ! Deallocate temporary arrays
  !write(stdo,*), xi,yi
  !
  deallocate(xi,yi,zi)
  !
  ! Done
  !
end subroutine userdef_prep_model


!------------------------------------------------------------------------
! Define physical model here
!------------------------------------------------------------------------
subroutine userdef_setup_model()
  implicit none
  integer :: ix,iy,iz,ierr,index
  integer, parameter :: out_unit=20
  doubleprecision :: xc,yc,zc
  type(amr_branch), pointer :: b
  !
  ! Make sure that the dust data is read, because we will need
  ! that in this model
  !
  call read_dustdata(1)
  if(dust_nr_species.ne.1) stop 991
  !
  ! Allocate the dust/hp density array
  !
  allocate(dustdens(1:dust_nr_species,1:nrcells),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR: Could not allocate dust density array'
     stop
  endif
  allocate(dust_massdust(1:dust_nr_species),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR: Could not allocate dust mass array'
     stop
  endif
  !d
  do iz=1,amr_grid_nz
     do iy=1,amr_grid_ny
        do ix=1,amr_grid_nx
           b     => amr_grid_branch(ix,iy,iz)%link
           index = b%leafindex
           xc    = amr_finegrid_xc(b%ixyzf(1),1,b%level)
           yc    = amr_finegrid_xc(b%ixyzf(2),2,b%level)
           zc    = amr_finegrid_xc(b%ixyzf(3),3,b%level)
           dustdens(1,index) = userdef_dustdens(xc,yc,zc)
        enddo
     enddo
  enddo

  call compute_dust_mass()
  !write(stdo,*),dust_massdusttot/1.99d33


  !Normalise for disk mass
  do iz=1,amr_grid_nz
      do iy=1,amr_grid_ny
         do ix=1,amr_grid_nx
            b     => amr_grid_branch(ix,iy,iz)%link
            index = b%leafindex
            xc    = amr_finegrid_xc(b%ixyzf(1),1,b%level)
            yc    = amr_finegrid_xc(b%ixyzf(2),2,b%level)
            zc    = amr_finegrid_xc(b%ixyzf(3),3,b%level)
            ! 0.5 factor need when using equatorial symmetric disk model

            dustdens(1,index) = dustdens(1,index)*(0.5d0*userdef_diskMass*1.988435d33/dust_massdusttot)
         enddo
      enddo
   enddo
   open (unit=out_unit,file="dust_density.txt",action="write",status="replace")
   write (out_unit,*) MAXVAL(dustdens)
   write (out_unit,*) dustdens
   close (out_unit)

   call compute_dust_mass()
   write(stdo,*),dust_massdusttot/1.988435d33
   !  write(stdo,*) 'xc'
   !write(stdo,*),xc
   ! write(stdo,*) 'yc'
  !write(stdo,*),yc
   !   write(stdo,*) 'zc'
!write(stdo,*),zc


end subroutine userdef_setup_model
!-----------------------------------------------------------------------
!                           Density Function
!-----------------------------------------------------------------------
function userdef_dustdens(xc,yc,zc)
    use ioput_module
    implicit none
    doubleprecision :: userdef_dustdens
    doubleprecision :: xc,yc,zc
    doubleprecision :: rr,zz,rcyl
    doubleprecision :: hp,sig
    integer, parameter :: out_unit=20


    rcyl = xc*sin(yc)
    zz = xc*cos(yc)

    hp = userdef_h0*au*(rcyl/(userdef_r0*au))**(userdef_beta)
    sig = (rcyl/(userdef_r0*au))**(-userdef_alpha)

    if (rcyl.ge.userdef_diskInR*au .and. rcyl.le.userdef_diskOutR*au) then
        userdef_dustdens = sig*exp(-0.5D0*(zz/hp)**2)
        userdef_dustdens = max(userdef_dustdens,0.d0)
    else
        userdef_dustdens = 0.d0
    endif
!dust density
    !open (unit=out_unit,file="dust_density.txt",action="write",status="replace")
    !write (out_unit,*) MAXVAL(userdef_dustdens)
    !write (out_unit,*) userdef_dustdens
    !close (out_unit)

    return



end function userdef_dustdens


end module userdef_module
