! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download



! This is the interface between the Rutgers version of the ROMS ocean model and DART.
! The required public interfaces arguments CANNOT be changed.
!
! Unlike the UCLA version, the interface does not use/require precomputed 
! forward operator values output directly from ROMS (i.e. s4dvar.in:MODname).
! A model interpolation routine is implemented to find the expected value
! of the observation along any vertical location. The interpolation is 
! implemented using the quad utility functions.   
!
! In a barotropic mode, the Rutgers ROMS applies a LF-AM3 with FB Feedback 
! numerical stepping scheme (myroms.org/wiki/Numerical_Solution_Technique)
! It also applies an AB3 scheme for 3D momenta
! LF : Leapfrog;
! AM3: Three-time Adams-Moulton corrector;
! FB : Forward-Backward;
! AB3: Third-order Adams-Bashforth step
!
! As such, several variables in the restart files have more than one
! time level. Instead of dealing with multi-level variables, we opted in 
! this interface to use the ROMS history files as ensemble members. All
! variables in these files have single time levels, making interpolation 
! and other computations a lot more straight-forward. 
!
! The conversion from the sigma terrain-following (S-level) vertical 
! coordinate for physical depth (Z-level) is done: (1) either by 
! reading the z-level coordinates directly from the input file if 
! they are available, or (2) by actually computing them. The 
! computation here is done in a similar fashion to ROMS.



module model_mod

! Modules that are absolutely required for use are listed
use             types_mod, only : r4, r8, digits12, SECPERDAY, DEG2RAD, rad2deg, PI,   &
                                  MISSING_I, MISSING_R8, i4, i8, vtablenamelength
use      time_manager_mod, only : time_type, set_time, set_date, get_date, get_time,   &
                                  print_time, print_date,                              &
                                  set_calendar_type, get_calendar_type,                &
                                  operator(*),  operator(+), operator(-),              &
                                  operator(>),  operator(<), operator(/),              &
                                  operator(/=), operator(<=)
use          location_mod, only : location_type, set_location, get_location,           &
                                  write_location, set_location_missing,                &
                                  loc_get_close_obs => get_close_obs,                  &
                                  loc_get_close_state => get_close_state,              &
                                  convert_vertical_obs, convert_vertical_state,        &
                                  VERTISHEIGHT, VERTISSURFACE, query_location,         &
                                  set_vertical, get_close_type
use         utilities_mod, only : register_module, error_handler, do_nml_term,         &
                                  E_MSG, logfileunit, nmlfileunit, file_to_text,       &
                                  get_unit, do_output, to_upper, do_nml_file,          &
                                  find_namelist_in_file, check_namelist_read,          &
                                  open_file, file_exist, find_textfile_dims, E_ERR,    &
                                  close_file, string_to_real, string_to_logical
use        quad_utils_mod, only : quad_interp_handle, set_quad_coords,                 &
                                  init_quad_interp, finalize_quad_interp,              &
                                  quad_lon_lat_locate, quad_lon_lat_evaluate,          &
                                  GRID_QUAD_FULLY_IRREGULAR, QUAD_LOCATED_LON_EDGES,   &
                                  QUAD_LOCATED_CELL_CENTERS, QUAD_LOCATED_LAT_EDGES
use          obs_kind_mod, only : QTY_TEMPERATURE, QTY_SALINITY, QTY_DRY_LAND,         &
                                  QTY_U_CURRENT_COMPONENT, QTY_V_CURRENT_COMPONENT,    &
                                  QTY_SEA_SURFACE_HEIGHT, get_index_for_quantity
use   state_structure_mod, only : get_dart_vector_index, get_model_variable_indices,   &
                                  get_domain_size, get_varid_from_kind,  add_domain,   &
                                  state_structure_info, get_kind_index, get_num_dims
use  netcdf_utilities_mod, only : nc_add_global_attribute, nc_synchronize_file,        &    
                                  nc_add_global_creation_time, nc_begin_define_mode,   &
                                  nc_end_define_mode, nc_open_file_readonly,           &    
                                  nc_close_file, nc_get_global_attribute,              &    
                                  nc_get_dimension_size, nc_get_variable, nc_check,    &
                                  nc_get_attribute_from_variable, nc_synchronize_file, &
                                  nc_variable_exists, nc_dimension_exists 
use       location_io_mod, only : nc_write_location_atts, nc_get_location_varids,      &
                                  nc_write_location
use     default_model_mod, only : nc_write_model_vars, init_conditions, init_time,     &
                                  adv_1step
use     mpi_utilities_mod, only : my_task_id
use        random_seq_mod, only : random_seq_type, init_random_seq, random_gaussian
use  ensemble_manager_mod, only : ensemble_type
use distributed_state_mod, only : get_state, get_state_array
use netcdf

implicit none
private

! these routines must be public and you cannot change
! the arguments - they will be called *from* the DART code.

! routines in this list have code in this module
public :: get_model_size,                      &
          get_state_meta_data,                 &
          model_interpolate,                   &
          shortest_time_between_assimilations, &
          static_init_model,                   &
          end_model,                           &
          nc_write_model_atts,                 &
          write_model_time,                    &
          read_model_time

! code for these routines are in other modules
public :: nc_write_model_vars,                 &
          pert_model_copies,                   &
          adv_1step,                           &
          init_time,                           &
          init_conditions,                     &
          convert_vertical_obs,                &
          convert_vertical_state,              &
          get_close_obs,                       &
          get_close_state

! not required interfaces but useful for utility programs
public :: get_time_information
          
! version controlled file description for error handling, do not edit
character(len=256), parameter :: source  = 'ROMS/model_mod.f90'

character(len=512) :: string1, string2, string3
logical, save      :: module_initialized = .false.

! things which can/should be in the model_nml
logical  :: output_state_vector         = .false.
integer  :: assimilation_period_days    = 1
integer  :: assimilation_period_seconds = 0
integer  :: vert_localization_coord     = VERTISHEIGHT
integer  :: debug                       = 0                ! turn up for more and more debug messages
character(len=256) :: roms_filename     = 'roms_input.nc'

! DART contents are specified in the input.nml:&model_nml namelist.
integer, parameter              :: MAX_STATE_VARIABLES                          = 8
integer, parameter              :: table_columns                                = 5
character(len=vtablenamelength) :: variables(MAX_STATE_VARIABLES*table_columns) = ' '
character(len=vtablenamelength) :: var_names(MAX_STATE_VARIABLES)               = ' '
logical                         :: update_list(MAX_STATE_VARIABLES)             = .FALSE.
integer                         :: kind_list(MAX_STATE_VARIABLES)               = MISSING_I
real(r8)                        :: clamp_vals(MAX_STATE_VARIABLES, 2)           = MISSING_R8

namelist /model_nml/            &
   assimilation_period_days,    &
   assimilation_period_seconds, &
   roms_filename,               &
   vert_localization_coord,     &
   debug,                       &
   variables

integer :: nfields ! This is the number of variables in the DART state vector
integer :: domid   ! global variable for state_structure_mod routines

! Interpolation grid handles 
type(quad_interp_handle) :: interp_t_grid, &
                            interp_u_grid, &
                            interp_v_grid

! Grid parameters - the values will be read from a
! standard ROMS namelist and filled in here
! Nx, Ny and Nz are the size of the rho grids
integer  :: Nx = -1, Ny = -1, Nz = -1
integer  :: Nu = -1, Nv = -1, Nw = -1
integer  :: Nc = 4          ! number of corners of the quad for interpolation
integer  :: Nd = 3          ! 3D location for the obs
real(r8) :: Zm = -5000.0_r8 ! a masking factor to account for land

! model_interpolate failure codes
integer, parameter :: QUAD_LOCATE_FAILED   = 13
integer, parameter :: QUAD_EVALUATE_FAILED = 21
integer, parameter :: SSH_QUAD_EVAL_FAILED = 34
integer, parameter :: QUAD_MAYBE_ON_LAND   = 55
integer, parameter :: OBS_TOO_DEEP         = 89

! ROMS related grid variables
real(r8), allocatable :: ULAT(:,:), ULON(:,:), UDEP(:,:,:), &
                         TLAT(:,:), TLON(:,:), TDEP(:,:,:), &
                         VLAT(:,:), VLON(:,:), VDEP(:,:,:)
logical, allocatable  :: TMSK(:,:), UMSK(:,:), VMSK(:,:)
real(r8), allocatable :: h(:,:), Cr(:), sr(:)  
real(r8)              :: hc         ! critical depth (m)
integer               :: Vt         ! transformation formula from ROMS
integer               :: ix, iy, ik

type(time_type) :: model_timestep
integer         :: model_size       ! the state vector length

contains


!-----------------------------------------------------------------------
! Called to do one time initialization of the model.
! In this case, it reads in the grid information, the namelist
! containing the variables of interest, where to get them, their size,
! their associated DART quantity, etc.
!
! In addition to harvesting the model metadata (grid,
! desired model advance step, etc.), it also fills a structure
! containing information about what variables are where in the DART
! framework.

subroutine static_init_model()

integer :: iunit, io
integer :: ss, dd
integer :: ncid

character(len=32) :: calendar
character(len=*), parameter :: routine = 'static_init_model'

type(time_type) :: model_time

if ( module_initialized ) return

! * read in the grid sizes from grid file
! * allocate space, and read in actual grid values
! * figure out model timestep
! * Compute the model size

! Print module information to log file
call register_module(source)

module_initialized = .true.

! Read the DART namelist for this model
call find_namelist_in_file('input.nml', 'model_nml', iunit)
read(iunit, nml = model_nml, iostat = io)
call check_namelist_read(iunit, io, 'model_nml')

! Record the namelist values used for the run
if (do_nml_file()) write(nmlfileunit, nml=model_nml)
if (do_nml_term()) write(     *     , nml=model_nml)

model_timestep = set_model_time_step()
call get_time(model_timestep, ss, dd)

write(string1, *) 'Assimilation period is ', dd, ' days ', ss, ' seconds'
call error_handler(E_MSG, routine, string1, source)

ncid = nc_open_file_readonly(roms_filename, routine)

call get_time_information(roms_filename, ncid, 'ocean_time', 'ocean_time', &
                          calendar=calendar, last_time=model_time)

call set_calendar_type(trim(calendar))

! Get the ROMS grid -- sizes and variables.
call get_grid_dimensions()
call get_grid()

! parse_variable_input() fills var_names, kind_list, clamp_vals, update_list
call parse_variable_input(variables, nfields)

domid = add_domain(roms_filename, & 
                         nfields, &
                       var_names, &
                       kind_list, &
                      clamp_vals, &
                     update_list)

if (debug > 0) call state_structure_info(domid)

call nc_close_file(ncid, routine)

model_size = get_domain_size(domid)

call setup_interpolation()

end subroutine static_init_model


!-----------------------------------------------------------------------
! Returns the size of the DART state vector (i.e. model) as an integer.
! Required for all applications.

function get_model_size()

integer(i8) :: get_model_size

if (.not. module_initialized) call static_init_model

get_model_size = model_size

end function get_model_size


!-----------------------------------------------------------------------
! Given an integer index into the state vector structure, returns the
! associated location. A second intent(out) optional argument quantity
! can be returned if the model has more than one type of field (for
! instance temperature and zonal wind component). This interface is
! required for all filter applications as it is required for computing
! the distance between observations and state variables.
!
! index_in: the index into the DART state vector
! location: the location at that index
! qty:      the DART QUANTITY at that index

subroutine get_state_meta_data(index_in, location, qty)

integer(i8),         intent(in)  :: index_in
type(location_type), intent(out) :: location
integer, optional,   intent(out) :: qty

! Local variables
integer  :: iloc, jloc, kloc
integer  :: myvarid, myqty

if (.not. module_initialized) call static_init_model

call get_model_variable_indices(index_in, iloc, jloc, kloc, var_id=myvarid)

myqty = get_kind_index(domid, myvarid)

if (myqty == QTY_U_CURRENT_COMPONENT) then
   location = set_location(ULON(iloc,jloc), ULAT(iloc,jloc), UDEP(iloc,jloc,kloc), VERTISHEIGHT)

elseif (myqty == QTY_V_CURRENT_COMPONENT) then
   location = set_location(VLON(iloc,jloc), VLAT(iloc,jloc), VDEP(iloc,jloc,kloc), VERTISHEIGHT)

elseif (myqty == QTY_SEA_SURFACE_HEIGHT) then
   ! TODO: So, because we are also doing waves and tides SSH (or zeta) could have a negative depth 
   ! Need to check this down the line
   location = set_location(TLON(iloc,jloc), TLAT(iloc,jloc), 0.0_r8, VERTISSURFACE)

else  ! Everything else is assumed to be on the rho points
   location = set_location(TLON(iloc,jloc), TLAT(iloc,jloc), TDEP(iloc,jloc,kloc), VERTISHEIGHT)

endif

! return state quantity for this index if requested
if (present(qty)) then 
   qty = myqty
   if (point_on_land(qty, iloc, jloc)) qty = QTY_DRY_LAND
endif

end subroutine get_state_meta_data


!-----------------------------------------------------------------------
! Model interpolate will interpolate any DART state variable
! (i.e. T, S, U, V, zeta) to the given location given a state vector.
! The type of the variable being interpolated is obs_type since
! normally this is used to find the expected value of an observation
! at some location. The interpolated value is returned in interp_vals
! and istatus is 0 for success. NOTE: This is a workhorse routine and is
! the basis for all the forward observation operator code.
!
! state_handle: DART ensemble handle
! ens_size:     DART ensemble size
! location:     Location of interest
! qty:          DART quantity
! expected_obs: Estimated value of the DART state at the location
!               of interest (the interpolated value).
! istatus:      Interpolation status ... 0 == success, /=0 is a failure

subroutine model_interpolate(state_handle, ens_size, location, qty, expected_obs, istatus)

type(ensemble_type), intent(in)  :: state_handle
integer,             intent(in)  :: ens_size
type(location_type), intent(in)  :: location
integer,             intent(in)  :: qty
real(r8),            intent(out) :: expected_obs(:)
integer,             intent(out) :: istatus(:)

! ---- Local
integer     :: varid, lstatus
integer     :: qstatus, vstatus 
integer     :: sshid, salid, i
integer     :: lon_c(Nc), lat_c(Nc)  ! lon and lat indices of the 4 quad corners
real(r8)    :: pdbar(ens_size)       ! Depth of input point converted from (m) to (dbars)
real(r8)    :: expected_T(ens_size)  ! Ensemble of expected values for potential temperature
real(r8)    :: expected_S(ens_size)  ! Ensemble of expected values for salinity
real(r8)    :: lon_lat_vrt(Nd)       ! lon, lat, vert of the point to interpolate
real(r8)    :: corners(Nc, ens_size) ! State values at the quad corners
real(r8)    :: SSHcorn(Nc, ens_size) ! SSH values at the quad corners
integer(i8) :: dartidx               ! Index into the DART state
logical     :: on_land              


type(quad_interp_handle) :: interp

if (.not. module_initialized) call static_init_model

! Successful istatus is 0
expected_obs = MISSING_R8
istatus = 99

! Get the id of the state variable for interpolation 
! Numbered according to the table list in the input.nml
! For instance: 
! T -> 1; S -> 2; U -> 3; V -> 4, Z -> 5
! If it's 0 then, it's not in the state and 
! we should just fail
varid = get_varid_from_kind(domid, qty)
if (varid < 0) return 

! Unpack the location of the obs 
lon_lat_vrt = get_location(location)

! Query the grid 
! Should have 3 options: 
! T grid, U grid, and V grid 
interp = get_interp_handle(qty)

! Locate the quad horizontal and get 
! the indices of its four corners
call quad_lon_lat_locate(interp, lon_lat_vrt(1), lon_lat_vrt(2), &
                         lon_c, lat_c, lstatus) 
if (lstatus /= 0) then  
   istatus = QUAD_LOCATE_FAILED
   return 
endif

if (debug > 0) then 
   do i = 1, Nc
      write(*, '(A, i1, A, i4, A, i4, A, f10.6, A, f10.6, A)') &
               'Corner #', i, ': [', lon_c(i), ',', lat_c(i), '] -> [', &
               TLON(lon_c(i), lat_c(i)),  ',', TLAT(lon_c(i), lat_c(i)), ']'
   enddo
endif

! Any part of the quad on land? 
call quad_on_land(qty, lon_c, lat_c, on_land)  
if (on_land) then 
   istatus = QUAD_MAYBE_ON_LAND
   return 
endif

! We always need to get the value of SSH at the 4 corners
! because we need to compute the depth using it.
! Recall that SSH can be both -ve and +ve
sshid = get_varid_from_kind(domid, QTY_SEA_SURFACE_HEIGHT)
do i = 1, Nc
   dartidx       = get_dart_vector_index(lon_c(i), lat_c(i), 1, domid, sshid)
   SSHcorn(i, :) = get_state(dartidx, state_handle)
enddo

! Find the state at the four corners of the quad
! No vertical interpolation for SSH
if (qty == QTY_SEA_SURFACE_HEIGHT) then
   ! TODO: Any other surface obs? Maybe anomalies?  
   ! For now, the assumption is that we're working with history files
   ! so each variable has a single time level. If for some reason we
   ! want to go back to using the restart files, here we need to query
   ! the number of dimensions for each variable and interpolate along
   ! the proper time level. 
   call quad_lon_lat_evaluate(interp, lon_lat_vrt(1), lon_lat_vrt(2), &
                              lon_c, lat_c, ens_size, SSHcorn,        & 
                              expected_obs, qstatus)
   istatus = 0
   if (qstatus /= 0) istatus = SSH_QUAD_EVAL_FAILED
   return
endif

! If the qty is not at the surface, then 
! we need to do some vertical interpolation 
call vert_interp(varid, ens_size, lon_lat_vrt, lon_c, lat_c, &
                 state_handle, SSHcorn, corners, vstatus)

if (vstatus /= 0) then 
   istatus = vstatus
   return
endif

! Do the interpolation
call quad_lon_lat_evaluate(interp, lon_lat_vrt(1), lon_lat_vrt(2), &
                           lon_c, lat_c, ens_size, corners,        & 
                           expected_obs, qstatus)

if (qstatus /= 0) then
   istatus = QUAD_EVALUATE_FAILED
   return
endif

if(qty == QTY_TEMPERATURE) then
  ! Set the potential temperature ensemble values 
  expected_T = expected_obs

  ! Need Salinity to compute in-situ temperature
  salid = get_varid_from_kind(domid, QTY_SALINITY)

  ! Pressure in decibars 
  pdbar = 0.59808_r8*(exp(-0.025_r8*lon_lat_vrt(3)) - 1.0_r8) + &
          1.00766_r8*lon_lat_vrt(3) + 2.28405e-6_r8*lon_lat_vrt(3)**2

  ! Vertical interp for salt
  call vert_interp(salid, ens_size, lon_lat_vrt, lon_c, lat_c, &
                 state_handle, SSHcorn, corners, vstatus)
  
  call quad_lon_lat_evaluate(interp, lon_lat_vrt(1), lon_lat_vrt(2), &
                             lon_c, lat_c, ens_size, corners,        &
                             expected_S, qstatus)

  if (qstatus /= 0) then
     istatus = QUAD_EVALUATE_FAILED
     return
  endif

  ! Deduce the in-situ temperature values
  expected_obs = sensible_temp(expected_T, expected_S, pdbar) 
endif

istatus = 0

end subroutine model_interpolate


!-----------------------------------------------------------------------
! Returns the the time step of the model; the smallest increment in
! time that the model is capable of advancing the ROMS state.

function shortest_time_between_assimilations()

type(time_type) :: shortest_time_between_assimilations

if (.not. module_initialized) call static_init_model

shortest_time_between_assimilations = model_timestep

end function shortest_time_between_assimilations


!-----------------------------------------------------------------------
! Writes the model-specific attributes to a DART 'diagnostic' netCDF file.
! This includes coordinate variables and some metadata, but NOT the
! actual DART state.
!
! ncid the netCDF handle of the DART diagnostic file opened by
!                 assim_model_mod:init_diag_output

subroutine nc_write_model_atts(ncid, domain_id)

integer, intent(in) :: ncid      ! netCDF file identifier
integer, intent(in) :: domain_id

integer :: nDimensions, nVariables, nAttributes, unlimitedDimID

! for the dimensions and coordinate variables
integer :: nxirhoDimID, nxiuDimID, nxivDimID
integer :: netarhoDimID, netauDimID, netavDimID
integer :: nsrhoDimID, nswDimID
integer :: VarID

! local variables
character(len=256) :: filename

if (.not. module_initialized) call static_init_model

! Write Global Attributes
call nc_begin_define_mode(ncid)
call nc_add_global_creation_time(ncid)
call nc_add_global_attribute(ncid, "model_source", source)
call nc_add_global_attribute(ncid, "model", "ROMS")
call nc_end_define_mode(ncid)

! Flush the buffer and leave netCDF file open
call nc_synchronize_file(ncid)

end subroutine nc_write_model_atts


!-----------------------------------------------------------------------
! Using a pre-defined factor, create ensemble perturbations around 
! a single state.  

subroutine pert_model_copies(state_ens_handle, ens_size, perturbation_amplitude, interf_provided)

type(ensemble_type), intent(inout) :: state_ens_handle
integer,             intent(in)    :: ens_size
real(r8),            intent(in)    :: perturbation_amplitude
logical,             intent(out)   :: interf_provided

integer             :: var_type
integer             :: j,i 
integer(i8)         :: dart_index
type(location_type) :: location

! Storage for a random sequence for perturbing a single initial state
type(random_seq_type) :: random_seq

if (.not. module_initialized) call static_init_model

interf_provided = .true.

! Initialize random number sequence
call init_random_seq(random_seq, my_task_id())

! Perturb the salinity and temperature fields
do i = 1, state_ens_handle%my_num_vars
   dart_index = state_ens_handle%my_vars(i)
   call get_state_meta_data(dart_index, location, var_type)

   if (var_type == QTY_TEMPERATURE .or. var_type == QTY_SALINITY) then
      do j = 1, ens_size
         state_ens_handle%copies(j,i) = random_gaussian(random_seq, &
            state_ens_handle%copies(j,i), perturbation_amplitude)
      enddo
   endif

enddo

end subroutine pert_model_copies


!-----------------------------------------------------------------------
! writes the time of the current state and (optionally) the time
! to be conveyed to ROMS to dictate the length of the forecast.
! This file is then used by scripts to modify the ROMS run.
! The format in the time information is totally at your discretion.
!
! ncfile_out name of the file
! model_time the current time of the model state
! adv_to_time the time in the future of the next assimilation.

subroutine write_model_time(ncid, model_time, adv_to_time)

character(len=*), parameter :: routine = 'write_model_time'

integer,         intent(in)           :: ncid
type(time_type), intent(in)           :: model_time
type(time_type), intent(in), optional :: adv_to_time

integer         :: io, varid, seconds, days
type(time_type) :: origin_time, deltatime
real(digits12)  :: run_duration

if (.not. module_initialized) call static_init_model

if (present(adv_to_time)) then
   string3 = time_to_string(adv_to_time)
   
   write(string1, *) 'ROMS/DART not configured to advance ROMS.'
   write(string2, *) 'called with optional advance_to_time of'
   
   call error_handler(E_ERR, routine, string1, &
              source, text2=string2,text3=string3)
endif

! If the ocean_time variable exists, we are updating a ROMS file,
! if not ... must be updating a DART diagnostic file.

io = nf90_inq_varid(ncid, 'ocean_time', varid)
if (io == NF90_NOERR) then
   call get_time_information('unknown', ncid, 'ocean_time', 'ocean_time', &
                myvarid=varid, origin_time=origin_time)
   deltatime = model_time - origin_time
   call get_time(deltatime, seconds, days)
   run_duration = real(days,digits12)*86400.0_digits12 + real(seconds,digits12)
   call nc_check(nf90_put_var(ncid, varid, run_duration), routine, 'put_var')
   return
endif

io = nf90_inq_varid(ncid, 'time', varid)
if (io == NF90_NOERR) then
   call get_time_information('unknown', ncid, 'time', 'time', &
                myvarid=varid, origin_time=origin_time)
   deltatime = model_time - origin_time
   call get_time(deltatime, seconds, days)
   run_duration = real(days,digits12)*86400.0_digits12 + real(seconds,digits12)
   call nc_check(nf90_put_var(ncid, varid, run_duration), routine, 'put_var')
   return
endif

end subroutine write_model_time


!--------------------------------------------------------------------
! read the time from the input file
! filename name of file that contains the time

function read_model_time(filename)

character(len=*), parameter :: routine = 'read_model_time'

character(len=*), intent(in) :: filename
type(time_type)              :: read_model_time

integer :: ncid

if (.not. module_initialized) call static_init_model

if (.not. file_exist(filename)) then
   write(string1,*) 'cannot open file ', trim(filename),' for reading.'
   call error_handler(E_ERR, routine, string1, source)
endif

ncid = nc_open_file_readonly(trim(filename), routine)

call get_time_information(filename, ncid, 'ocean_time', 'ocean_time', last_time=read_model_time)
call nc_close_file(ncid, routine)

end function read_model_time


!-----------------------------------------------------------------------
! Set the desired minimum model advance time. This is generally NOT the
! dynamical timestep of the model, but rather the shortest forecast length
! you are willing to make. This impacts how frequently the observations
! may be assimilated.

function set_model_time_step()

type(time_type) :: set_model_time_step

! assimilation_period_seconds, assimilation_period_days are from the namelist
set_model_time_step = set_time(assimilation_period_seconds, assimilation_period_days)

end function set_model_time_step


!-----------------------------------------------------------------------
! Read the grid dimensions from the ROMS grid netcdf file.
! By reading the dimensions first, we can use them in variable
! declarations later - which is faster than using allocatable arrays.

subroutine get_grid_dimensions()

character(len=*), parameter :: routine = 'get_grid_dimensions'

integer :: ncid
integer :: Nxi_rho, Nxi_u, Nxi_v
integer :: Neta_rho, Neta_u, Neta_v
integer :: Ns_rho, Ns_w

! Read the (static) grid dimensions from the ROMS grid file.
ncid = nc_open_file_readonly(roms_filename, routine)

Nxi_rho  = nc_get_dimension_size(ncid, 'xi_rho',  routine)
Nxi_u    = nc_get_dimension_size(ncid, 'xi_u',    routine)
Nxi_v    = nc_get_dimension_size(ncid, 'xi_v',    routine)
Neta_rho = nc_get_dimension_size(ncid, 'eta_rho', routine)
Neta_u   = nc_get_dimension_size(ncid, 'eta_u',   routine)
Neta_v   = nc_get_dimension_size(ncid, 'eta_v',   routine)

! Read the vertical dimensions 
Ns_rho   = nc_get_dimension_size(ncid, 's_rho', routine)
Ns_w     = nc_get_dimension_size(ncid, 's_w'  , routine)

call nc_close_file(ncid, routine)

Nx = Nxi_rho  ! Setting the nominal value of the 'global' variables
Ny = Neta_rho ! Setting the nominal value of the 'global' variables
Nz = Ns_rho   ! Setting the nominal value of the 'global' variables

Nu = Nxi_u
Nv = Neta_v
Nw = Ns_w

end subroutine get_grid_dimensions


!-----------------------------------------------------------------------
! Read the actual grid values from the ROMS netcdf file.
! We might need to adapt this to allow for leapfrog time stepping
! The assumption is that we're reading from history files 
! which only have one time-level. 

subroutine get_grid()

integer                     :: ncid
character(len=*), parameter :: routine = 'get_grid'

real(r8), allocatable       :: zeta(:,:), zeta_lf(:,:,:) 
real(r8), allocatable       :: mask(:,:)   ! Land mask: 0 land & 1 water

allocate(ULAT(Nu, Ny), ULON(Nu, Ny), UMSK(Nu, Ny), UDEP(Nu, Ny, Nz))
allocate(VLAT(Nx, Nv), VLON(Nx, Nv), VMSK(Nx, Nv), VDEP(Nx, Nv, Nz))
allocate(TLAT(Nx, Ny), TLON(Nx, Ny), TMSK(Nx, Ny), TDEP(Nx, Ny, Nz))

allocate(h(Nx,Ny), Cr(Nz), sr(Nz))

! s-coordinates at RHO-points
!do ik = 1, Nz
!   sr(ik) = (ik - 0.5_r8 - Nz) / Nz
!enddo
! s-coordinates at W-points
!do ik = 0, Nz
!   sw(ik+1) = (real(ik) - Nz) / Nz
!enddo

! Read the vertical information 
ncid = nc_open_file_readonly(roms_filename, routine)

! Varibles needed to calculate depth given the free surface
call nc_get_variable(ncid, 'h'         , h , routine)
call nc_get_variable(ncid, 'hc'        , hc, routine)
call nc_get_variable(ncid, 'Cs_r'      , Cr, routine)
call nc_get_variable(ncid, 's_rho'     , sr, routine)
call nc_get_variable(ncid, 'Vtransform', Vt, routine)

! Check whether physical depth data exist in the restart file
if (nc_variable_exists(ncid, 'z_rho')) then 

   ! It exists, so grab all of them
   call nc_get_variable(ncid, 'z_rho', TDEP, routine)
   call nc_get_variable(ncid, 'z_u'  , UDEP, routine)
   call nc_get_variable(ncid, 'z_v'  , VDEP, routine)
   TDEP = -TDEP
   UDEP = -UDEP
   VDEP = -VDEP
else
   
   ! Need to compute the physical depth
   ! Check whether zeta is coming from a history or a restart file
   ! The 3rd dimension shouldn't exist in a history file
   allocate(zeta(Nx, Ny))

   if (nc_dimension_exists(ncid, 'three')) then 
      allocate(zeta_lf(Nx, Ny, 3))
      zeta = zeta_lf(:, :, 1)
   else 
      call nc_get_variable(ncid, 'zeta', zeta, routine)
   endif 

   ! Find depth at the center of the cells
   do ix = 1, Nx
      do iy = 1, Ny 
         call compute_physical_depth(h(ix, iy), zeta(ix, iy), TDEP(ix, iy, :))
      enddo
   enddo

   ! Average in xi and eta directions
   call U_V_depths()
endif

! Read in the land mask
TMSK = .false.
UMSK = .false.
VMSK = .false.

allocate(mask(Nx, Ny))
call nc_get_variable(ncid, 'mask_rho', mask, routine)
where(mask < 1.0_r8) TMSK = .true. ! mask is active where land is
deallocate(mask)

allocate(mask(Nu, Ny))
call nc_get_variable(ncid, 'mask_u', mask, routine)
where(mask < 1.0_r8) UMSK = .true.
deallocate(mask)

allocate(mask(Nx, Nv))
call nc_get_variable(ncid, 'mask_v', mask, routine)
where(mask < 1.0_r8) VMSK = .true.
deallocate(mask)

! Read the rest of the grid information from the traditional grid file
call nc_get_variable(ncid, 'lon_rho', TLON, routine)
call nc_get_variable(ncid, 'lat_rho', TLAT, routine)
call nc_get_variable(ncid, 'lon_u'  , ULON, routine)
call nc_get_variable(ncid, 'lat_u'  , ULAT, routine)
call nc_get_variable(ncid, 'lon_v'  , VLON, routine)
call nc_get_variable(ncid, 'lat_v'  , VLAT, routine)

where (TLON < 0.0_r8) TLON = TLON + 360.0_r8
where (ULON < 0.0_r8) ULON = ULON + 360.0_r8
where (VLON < 0.0_r8) VLON = VLON + 360.0_r8

! Be aware that all the depths are negative values.
! The surface of the ocean is 0.0, the deepest is a big negative value.
if (do_output() .and. debug > 0) then
    write(string1,*) '    min/max ULON ', minval(ULON), maxval(ULON)
    write(string2,*)     'min/max ULAT ', minval(ULAT), maxval(ULAT)
    write(string3,*)     'min/max UDEP ', minval(UDEP, MASK=UDEP>Zm), maxval(UDEP)
    call error_handler(E_MSG, routine, string1, text2=string2, text3=string3)

    write(string1,*) '    min/max VLON ', minval(VLON), maxval(VLON)
    write(string2,*)     'min/max VLAT ', minval(VLAT), maxval(VLAT)
    write(string3,*)     'min/max VDEP ', minval(VDEP, MASK=VDEP>Zm), maxval(VDEP)
    call error_handler(E_MSG, routine, string1, text2=string2, text3=string3)

    write(string1,*) '    min/max TLON ', minval(TLON), maxval(TLON)
    write(string2,*)     'min/max TLAT ', minval(TLAT), maxval(TLAT)
    write(string3,*)     'min/max TDEP ', minval(TDEP, MASK=TDEP>Zm), maxval(TDEP)
    call error_handler(E_MSG, routine, string1, text2=string2, text3=string3)
endif

call nc_close_file(ncid, routine)

end subroutine get_grid


!-----------------------------------------------------------------------
! Get physical depth z_rho using ROMS formulations:
! refer to set_depth.F
subroutine compute_physical_depth(b, z, d)

character(len=*), parameter :: routine = 'compute_physical_depth'

real(r8), intent(in)  :: b, z   ! bathymetry and free surface
real(r8), intent(out) :: d(Nz)  ! depth at all levels
real(r8)              :: z0

if (debug > 100) then 
   do ik = 1, Nz
      write(*, '(A, i3, A, f10.7)') 'Layer: ', ik, ', S_RHO: ', sr(ik) 
   enddo
endif

! Compute z at RHO points 
do ik = 1, Nz
   if (Vt == 1) then ! Original transformation
      z0    = hc*(sr(ik) - Cr(ik)) + Cr(ik)*b 
      d(ik) = z0 + z * (1.0_r8 + z0/b)
   elseif (Vt == 2) then ! New transformation
      z0    = (hc*sr(ik) + Cr(ik)*b)/(hc+b) 
      d(ik) = z + z0*(z+b) 
   else
      string1 = 'Unsupported Vtransform'
      call error_handler(E_ERR, routine, string1)
   endif 
enddo

! Reverse the sign: 
! Below the surface is +ve depth
! Surface is 0
! Above the surface is -ve depth
d = -d

end subroutine compute_physical_depth


!-----------------------------------------------------------------------
! Using z_rho, compute z_u and z_v
subroutine U_V_depths()

! Compute z_u; average in xi direction 
do ix = 1, Nx-1
   do iy = 1, Ny
      UDEP(ix, iy, :) = 0.5_r8 * (TDEP(ix, iy, :) + TDEP(ix+1, iy, :))
   enddo
enddo

! Compute z_v; average in eta direction 
do ix = 1, Nx
   do iy = 1, Ny-1
      VDEP(ix, iy, :) = 0.5_r8 * (TDEP(ix, iy, :) + TDEP(ix, iy+1, :))
   enddo
enddo

end subroutine U_V_depths 


!-----------------------------------------------------------------------
! Apply vertical interpolation using the four corners of the quad

subroutine vert_interp(varid, ens_size, lon_lat_vrt, lon_c, lat_c, &
                       state_handle, SSHcorn, corners, dstatus)

integer, intent(in)             :: varid
integer, intent(in)             :: ens_size 
real(r8), intent(in)            :: lon_lat_vrt(Nd)   
integer, intent(in)             :: lon_c(Nc)
integer, intent(in)             :: lat_c(Nc)
type(ensemble_type), intent(in) :: state_handle
real(r8), intent(in)            :: SSHcorn(Nc, ens_size) 
real(r8), intent(out)           :: corners(Nc, ens_size) ! State values at the quad corners

integer     :: i, ie, il
integer     :: dstatus, Nl, lz
integer     :: lev(ens_size, 2)      ! Bottom (1) and Top (2) level indices that bounds point of interest 
integer     :: levels(2*ens_size)    ! Unique vertical levels across the ensemble
real(r8)    :: depths(Nz)            ! Depths for each ensemble member
real(r8)    :: lev_frc(ens_size)     ! Fractional distances between bot and top 
real(r8)    :: Zvals(Nz, ens_size)   ! All ensemble values to be interpolated
real(r8)    :: tops, bots            ! Ensemble level values
real(r8)    :: bath                  ! Bathymetry at a single point
integer(i8) :: dartidx     

levels = 0

do i = 1, Nc
   bath = h(lon_c(i), lat_c(i))
  
   ! Physical depth is a function of SSH and differs 
   ! across the ensemble
   do ie = 1, ens_size
      call compute_physical_depth(bath, SSHcorn(i, ie), depths) 
      call depth_bounds(lon_lat_vrt(3), depths, lev(ie, 1), lev(ie, 2), lev_frc(ie), dstatus)
     
      if (dstatus /= 0) return
   enddo

   ! Now, get the state. First, find unique vertical levels 
   ! used in the ensemble of vertical interpolation
   ! This is to avoid computing/calling the state for different 
   ! members at the same level. So, we only get the state 
   ! ensemble at each of the interpolating levels once!
   call unique_levels(lev(:, 1), lev(:, 2), ens_size, levels, Nl)

   ! Get the state at all requested depths
   do il = 1, Nl
      lz = levels(il)
      dartidx      = get_dart_vector_index(lon_c(i), lat_c(i), lz, domid, varid)
      Zvals(lz, :) = get_state(dartidx, state_handle)

      !write(*, '(A, i3, A, i3, A, i12, A, i2, A, 3F10.6)') & 
      !         '  il: ', il, ', lvl: ', lz, ', DART indx: ', dartidx, ', Zvals(', lz, ', :): ', Zvals(lz, :)
   enddo

   ! Do the vertical interpolation using the 
   ! state at all depth array (i.e., Zvals)
   do ie = 1, ens_size
      bots = Zvals(lev(ie, 1), ie)
      tops = Zvals(lev(ie, 2), ie)

      corners(i, ie) = bots + lev_frc(ie)*(tops-bots)
   enddo
enddo

end subroutine vert_interp


!------------------------------------------------------------------
! Any model specific distance calcualtion can be done here
subroutine get_close_obs(gc, base_loc, base_type, locs, loc_qtys, loc_types, &
                         num_close, close_ind, dist, ens_handle)

type(get_close_type),          intent(in)    :: gc            ! handle to a get_close structure
integer,                       intent(in)    :: base_type     ! observation TYPE
type(location_type),           intent(inout) :: base_loc      ! location of interest
type(location_type),           intent(inout) :: locs(:)       ! obs locations
integer,                       intent(in)    :: loc_qtys(:)   ! QTYS for obs
integer,                       intent(in)    :: loc_types(:)  ! TYPES for obs
integer,                       intent(out)   :: num_close     ! how many are close
integer,                       intent(out)   :: close_ind(:)  ! incidies into the locs array
real(r8),            optional, intent(out)   :: dist(:)       ! distances in radians
type(ensemble_type), optional, intent(in)    :: ens_handle

character(len=*), parameter :: routine = 'get_close_obs'

call loc_get_close_obs(gc, base_loc, base_type, locs, loc_qtys, loc_types, &
                          num_close, close_ind, dist, ens_handle)

end subroutine get_close_obs


!------------------------------------------------------------------
! Any model specific distance calcualtion can be done here
subroutine get_close_state(gc, base_loc, base_type, locs, loc_qtys, loc_indx, &
                           num_close, close_ind, dist, ens_handle)

type(get_close_type),          intent(in)    :: gc           ! handle to a get_close structure
type(location_type),           intent(inout) :: base_loc     ! location of interest
integer,                       intent(in)    :: base_type    ! observation TYPE
type(location_type),           intent(inout) :: locs(:)      ! state locations
integer,                       intent(in)    :: loc_qtys(:)  ! QTYs for state
integer(i8),                   intent(in)    :: loc_indx(:)  ! indices into DART state vector
integer,                       intent(out)   :: num_close    ! how many are close
integer,                       intent(out)   :: close_ind(:) ! indices into the locs array
real(r8),            optional, intent(out)   :: dist(:)      ! distances in radians
type(ensemble_type), optional, intent(in)    :: ens_handle

character(len=*), parameter :: routine = 'get_close_state'

integer  :: ii

call loc_get_close_state(gc, base_loc, base_type, locs, loc_qtys, loc_indx, &
                            num_close, close_ind, dist, ens_handle)

if (.not. present(dist)) return

! Put any land points very far away
! so they are not updated by obs
do ii = 1, num_close
  if(loc_qtys(close_ind(ii)) == QTY_DRY_LAND) dist(ii) = 1.0e9_r8
enddo

end subroutine get_close_state


!------------------------------------------------------------
! calculate sensible (in-situ) temperature from 
! local pressure, salinity, and potential temperature

elemental function sensible_temp(potemp, s, lpres)

real(r8), intent(in)  :: potemp ! potential temperature in C
real(r8), intent(in)  :: s      ! salinity Practical Salinity Scale 1978 (PSS-78)
real(r8), intent(in)  :: lpres  ! pressure in decibars
real(r8) :: sensible_temp ! in-situ (sensible) temperature (C)

integer  :: i,j,n
real(r8) :: dp,p,q,r1,r2,r3,r4,r5,s1,t,x

s1 = s - 35.0_r8
p  = 0.0_r8
t  = potemp

dp = lpres - p
n  = int (abs(dp)/1000.0_r8) + 1
dp = dp/n

do i=1,n
   do j=1,4

      r1 = ((-2.1687e-16_r8 * t + 1.8676e-14_r8) * t - 4.6206e-13_r8) * p
      r2 = (2.7759e-12_r8*t - 1.1351e-10_r8) * s1
      r3 = ((-5.4481e-14_r8 * t + 8.733e-12_r8) * t - 6.7795e-10_r8) * t
      r4 = (r1 + (r2 + r3 + 1.8741e-8_r8)) * p + (-4.2393e-8_r8 * t+1.8932e-6_r8) * s1
      r5 = r4 + ((6.6228e-10_r8 * t-6.836e-8_r8) * t + 8.5258e-6_r8) * t + 3.5803e-5_r8

      x  = dp*r5

      if (j == 1) then
         t = t + 0.5_r8 * x
         q = x
         p = p + 0.5_r8 * dp
      
      else if (j == 2) then
         t = t + 0.29298322_r8 * (x-q)
         q = 0.58578644_r8 * x + 0.121320344_r8 * q

      else if (j == 3) then
         t = t + 1.707106781_r8 * (x-q)
         q = 3.414213562_r8*x - 4.121320344_r8*q
         p = p + 0.5_r8*dp

      else ! j must == 4
         t = t + (x - 2.0_r8 * q) / 6.0_r8

      endif

   enddo ! j loop
enddo ! i loop

sensible_temp = t

end function sensible_temp


!-----------------------------------------------------------------------
! Fill the array of requested variables, dart kinds, possible min/max
! values and whether or not to update the field in the output file.
!
! state_variables the list of variables and kinds from model_mod_nml
! ngood the number of variable/KIND pairs specified

subroutine parse_variable_input(state_variables, ngood)

character(len=*), parameter   :: routine = 'parse_variable_input'
character(len=*), intent(in)  :: state_variables(:)
integer,          intent(out) :: ngood

integer :: i
character(len=NF90_MAX_NAME) :: varname       ! column 1
character(len=NF90_MAX_NAME) :: dartstr       ! column 2
character(len=NF90_MAX_NAME) :: minvalstring  ! column 3
character(len=NF90_MAX_NAME) :: maxvalstring  ! column 4
character(len=NF90_MAX_NAME) :: state_or_aux  ! column 5   change to updateable

ngood = 0
MyLoop : do i = 1, MAX_STATE_VARIABLES

   varname      = trim(state_variables(table_columns*i-4))
   dartstr      = trim(state_variables(table_columns*i-3))
   minvalstring = trim(state_variables(table_columns*i-2))
   maxvalstring = trim(state_variables(table_columns*i-1))
   state_or_aux = trim(state_variables(table_columns*i  ))

   if ( varname == ' ' .and. dartstr == ' ' ) exit MyLoop ! Found end of list.

   if ( varname == ' ' .or. dartstr == ' ' ) then
      string1 = 'model_nml:model "variables" not fully specified'
      call error_handler(E_ERR, routine, string1, source)
   endif

   ! Make sure DART kind is valid
   if( get_index_for_quantity(dartstr) < 0 ) then
      write(string1,'(''there is no quantity <'',a,''> in obs_kind_mod.f90'')') trim(dartstr)
      call error_handler(E_ERR, routine, string1, source)
   endif

   call to_upper(minvalstring)
   call to_upper(maxvalstring)
   call to_upper(state_or_aux)

   var_names(   i) = varname
   kind_list(   i) = get_index_for_quantity(dartstr)
   clamp_vals(i,1) = string_to_real(minvalstring)
   clamp_vals(i,2) = string_to_real(maxvalstring)
   update_list( i) = string_to_logical(state_or_aux, 'UPDATE')

   ngood = ngood + 1

enddo MyLoop

if (ngood == MAX_STATE_VARIABLES) then
   string1 = 'WARNING: There is a possibility you need to increase ''MAX_STATE_VARIABLES'''
   write(string2,'(''WARNING: you have specified at least '',i4,'' perhaps more.'')')ngood
   call error_handler(E_MSG, routine, string1, source, text2=string2)
endif

end subroutine parse_variable_input


!-----------------------------------------------------------------------
! Find the named variable (often 'ocean_time') in a ROMS netCDF file.
! If it is not found, it is a fatal error.
!
! filename the name of the ROMS netCDF file
!                 (used to generate useful error messages).
! ncid the netCDF handle to the ROMS netCDF file.
! variable name which contains the time
! calendar the character string indicating the calendar in use
! last_time_index the value of the last time dimension
! last_time the time/date of the last time
! origin_time the base time other times are relative to
! all_times an array of all times in the variable

subroutine get_time_information(filename, ncid, var_name, dim_name, myvarid, &
               calendar, last_time_index, last_time, origin_time, all_times)

character(len=*),            intent(in)  :: filename
integer,                     intent(in)  :: ncid
character(len=*),            intent(in)  :: var_name
character(len=*),            intent(in)  :: dim_name
integer,           optional, intent(out) :: myvarid
character(len=32), optional, intent(out) :: calendar
integer,           optional, intent(out) :: last_time_index
type(time_type),   optional, intent(out) :: last_time
type(time_type),   optional, intent(out) :: origin_time
type(time_type),   optional, intent(out) :: all_times(:)

character(len=*), parameter   :: routine = 'get_time_information'

integer           :: ios, DimID, VarID, dimlen, i
character(len=64) :: unitstring
character(len=32) :: calendarstring

integer :: year, month, day, hour, minute, second, rc
real(digits12), allocatable :: these_times(:)
type(time_type) :: time_offset, base_time

logical :: offset_in_seconds  ! if .false., assuming offset in days

integer :: original_calendar_type

call nc_check(nf90_inq_dimid(ncid,dim_name,dimid=DimID), &
       routine,'cannot find "'//trim(dim_name)//'" dimension in '//trim(filename))

call nc_check(nf90_inquire_dimension(ncid, DimID, len=dimlen), &
       routine, 'inquire_dimension '//trim(dim_name)//' from '//trim(filename))
if (present(last_time_index)) last_time_index = dimlen

call nc_check(nf90_inq_varid(ncid, var_name, VarID), &
       routine, 'inq_varid '//trim(var_name)//' from '//trim(filename))
if (present(myvarid)) myvarid = VarID

! assume gregorian calendar unless there's a calendar attribute saying elsewise
rc = nf90_get_att(ncid, VarID, 'calendar', calendarstring)
if (rc /= nf90_noerr) calendarstring = 'gregorian'

if (index(calendarstring, 'gregorian') == 0) then
   write(string1, *) 'expecting '//trim(var_name)//' calendar of "gregorian"'
   write(string2, *) 'got '//trim(calendarstring)
   write(string3, *) 'from file "'//trim(filename)//'"'
   call error_handler(E_MSG,routine, string1, &
             source, text2=string2, text3=string3)
else
   ! coerce all forms of gregorian to the one DART supports
   ! 'gregorian_proleptic' needs to be changed, for example.
   calendarstring = 'gregorian'
endif

if (present(calendar)) calendar = trim(calendarstring)

if (present(last_time) .or. present(origin_time) .or. present(all_times)) then

   ! May need to put the calendar back to some original value
   original_calendar_type = get_calendar_type()

   ! We need to set the calendar to interpret the time values
   ! do we need to preserve the original calendar setting if there is one?
   call set_calendar_type( trim(calendarstring) )

   ! Make sure the calendar is expected form
   ! var_name:units    = "seconds since 1999-01-01 00:00:00" ;
   !                      1234567890123
   !   OR
   ! var_name:units    = "days since 1999-01-01 00:00:00" ;
   !                      1234567890

   call nc_check(nf90_get_att(ncid, VarID, 'units', unitstring), &
          routine, 'get_att '//trim(var_name)//' units '//trim(filename))

   ! decode the start time of the time variable - expecting time to be coded
   ! as an offset to some base
   if (unitstring(1:13) == 'seconds since') then
      read(unitstring,'(14x,i4,5(1x,i2))',iostat=ios)year,month,day,hour,minute,second
      if (ios /= 0) then
         write(string1, *) 'Unable to read time variable units. Error status was ',ios
         write(string2, *) 'expected "seconds since YYYY-MM-DD HH:MM:SS"'
         write(string3, *) 'was      "'//trim(unitstring)//'"'
         call error_handler(E_ERR, routine, string1, &
                source, text2=string2, text3=string3)
      endif
      offset_in_seconds = .true.

   else if (unitstring(1:10) == 'days since') then
      read(unitstring,'(11x,i4,5(1x,i2))',iostat=ios)year,month,day,hour,minute,second
      if (ios /= 0) then
         write(string1,*)'Unable to read time variable units. Error status was ',ios
         write(string2,*)'expected "days since YYYY-MM-DD HH:MM:SS"'
         write(string3,*)'was      "'//trim(unitstring)//'"'
         call error_handler(E_ERR, routine, string1, &
                source, text2=string2, text3=string3)
      endif
      offset_in_seconds = .false.

   else
      write(string1,*)'expecting time attribute units of "seconds since ..." -OR-'
      write(string2,*)'                              "days since ..."'
      write(string3,*)'got "'//trim(unitstring)//'"'
      call error_handler(E_ERR,routine, string1, &
                source, text2=string2, text3=string3)
   endif

   base_time = set_date(year, month, day, hour, minute, second)

   if (present(origin_time)) origin_time = base_time

   if (present(last_time) .or. present(all_times)) then

      ! big_integer may overflow a 32bit integer, so declare it 64bit
      ! and parse it into an integer number of days and seconds, both
      ! of which can be 32bit. Our set_time, set_date routines need 32bit integers.

      allocate(these_times(dimlen))

      call nc_check(nf90_get_var( ncid, VarID, these_times), &
             routine, 'get_var '//trim(var_name)//' from '//trim(filename))

      if (present(last_time)) then
         time_offset = convert_to_time_offset(these_times(dimlen), offset_in_seconds)
         last_time = base_time + time_offset
      endif

      if (present(all_times)) then
         do i=1, dimlen
            time_offset = convert_to_time_offset(these_times(i), offset_in_seconds)
            all_times(i) = base_time + time_offset
         enddo
      endif

      if (do_output() .and. debug > 0 .and. present(last_time)) then
         call print_time(last_time, str='last roms time is ',iunit=logfileunit)
         call print_time(last_time, str='last roms time is ')
         call print_date(last_time, str='last roms date is ',iunit=logfileunit)
         call print_date(last_time, str='last roms date is ')
      endif

      deallocate(these_times)

   endif

   call set_calendar_type(original_calendar_type)

endif

end subroutine get_time_information


!-----------------------------------------------------------------------
! convert a fractional day to a dart time type

function convert_to_time_offset(offset, offset_in_seconds)

real(digits12), intent(in) :: offset
logical,        intent(in) :: offset_in_seconds
type(time_type) :: convert_to_time_offset

integer(i8) :: big_integer
integer :: some_seconds, some_days

if (offset_in_seconds) then
   big_integer  = int(offset,i8)
   some_days    = big_integer / (24*60*60)
   big_integer  = big_integer - (some_days * (24*60*60))
   some_seconds = int(big_integer,i4)
else
   ! offset in fractional days
   some_days    = int(offset)
   some_seconds = (offset - some_days) * (24*60*60)
endif

convert_to_time_offset = set_time(some_seconds, some_days)

end function convert_to_time_offset


!-----------------------------------------------------------------------
! convert DART time type into a character string with the
! format of YYYYMMDDhh ... or DDhh
!
! time_to_string the character string containing the time
! t:        the time
! interval: logical flag describing if the time is to be
!           interpreted as a calendar date or a time increment.
!           If the flag is merely present, the time is to be
!           interpreted as an increment and the format is simply
!           DDhh. If the flag is not present, the time is a full
!           calendar (Gregorian) date and will be renedered with
!           the YYYYMMDDhh format.

function time_to_string(t, interval)

character(len=*), parameter   :: routine = 'time_to_string'

character(len=19)             :: time_to_string
type(time_type),   intent(in) :: t
logical, optional, intent(in) :: interval

! local variables
integer :: iyear, imonth, iday, ihour, imin, isec
integer :: ndays, nsecs
logical :: dointerval

if (present(interval)) then
   dointerval = interval
else
   dointerval = .false.
endif

! for interval output, output the number of days, then hours, mins, secs
! for date output, use the calendar routine to get the year/month/day hour:min:sec
if (dointerval) then
   call get_time(t, nsecs, ndays)
   if (ndays > 99) then
      write(string1, *) 'interval number of days is ', ndays
      call error_handler(E_ERR, routine, 'interval days cannot be > 99', &
                         source, text2=string1)
   endif
   ihour = nsecs / 3600
   nsecs = nsecs - (ihour * 3600)
   imin  = nsecs / 60
   nsecs = nsecs - (imin * 60)
   isec  = nsecs
   
   write(time_to_string, '(I2.2,I2.2)') ndays, ihour
else
   call get_date(t, iyear, imonth, iday, ihour, imin, isec)
   write(time_to_string, '(I4.4,5(A1,I2.2))') &
          iyear, '-', imonth, '-', iday, ' ', ihour, ':', imin, ':', isec
endif

end function time_to_string


!-----------------------------------------------------------------------
! Initialize the quad interpolation tools

subroutine setup_interpolation

! temp
call init_quad_interp(GRID_QUAD_FULLY_IRREGULAR, Nx, Ny, &
                      QUAD_LOCATED_CELL_CENTERS,         &
                      .false., .false., .false.,         & 
                      interp_t_grid)
call set_quad_coords(interp_t_grid, TLON, TLAT, TMSK)

! u-momentum
call init_quad_interp(GRID_QUAD_FULLY_IRREGULAR, Nu, Ny, &
                      QUAD_LOCATED_LON_EDGES,            &    
                      .false., .false., .false.,         &    
                      interp_u_grid)
call set_quad_coords(interp_u_grid, ULON, ULAT, UMSK)

! v-momentum
call init_quad_interp(GRID_QUAD_FULLY_IRREGULAR, Nx, Nv, &
                      QUAD_LOCATED_LAT_EDGES,            &    
                      .false., .false., .false.,         &    
                      interp_v_grid)
call set_quad_coords(interp_v_grid, VLON, VLAT, VMSK)

end subroutine setup_interpolation


!------------------------------------------------------------
! Returns the appropriate quad_interp handle

function get_interp_handle(qty)

type(quad_interp_handle) :: get_interp_handle
integer, intent(in)      :: qty

if (on_u_grid(qty)) then
  get_interp_handle = interp_u_grid
elseif (on_v_grid(qty)) then
  get_interp_handle = interp_v_grid
else
  get_interp_handle = interp_t_grid
endif

end function


!---------------------------------------------------------
! Check if we're on the U grid

function on_u_grid(qty)

integer, intent(in) :: qty
logical             :: on_u_grid

on_u_grid = .false.
if (qty == QTY_U_CURRENT_COMPONENT) on_u_grid = .true.

end function on_u_grid


!------------------------------------------------------------
! Check if we're on the V grid

function on_v_grid(qty)

integer, intent(in) :: qty
logical             :: on_v_grid

on_v_grid = .false.
if (qty == QTY_V_CURRENT_COMPONENT) on_v_grid = .true.

end function on_v_grid


!----------------------------------------------------------
! Check if we're on the T grid

function on_t_grid(qty)

integer, intent(in) :: qty
logical             :: on_t_grid

on_t_grid = .true.
if (qty == QTY_U_CURRENT_COMPONENT .or. & 
    qty == QTY_V_CURRENT_COMPONENT) on_t_grid = .false.

end function on_t_grid


!-----------------------------------------------------------
! Figure out whether a point is on land or not

function point_on_land(qty, ilon, ilat)

integer :: qty, ilon, ilat
logical :: point_on_land

if (on_u_grid(qty)) then 
   point_on_land = UMSK(ilon, ilat) 
elseif (on_v_grid(qty)) then 
   point_on_land = VMSK(ilon, ilat) 
else
   point_on_land = TMSK(ilon, ilat) 
endif

end function point_on_land


!-----------------------------------------------------------
! Figure out whether any part of the quad is on land 

subroutine quad_on_land(qty, lons, lats, lstatus)

integer, intent(in)  :: qty, lons(:), lats(:)
logical, intent(out) :: lstatus

lstatus = .false.

if (on_u_grid(qty)) then 
   lstatus = UMSK(lons(1), lats(1)) .or. &
             UMSK(lons(2), lats(2)) .or. &
             UMSK(lons(3), lats(3)) .or. &
             UMSK(lons(4), lats(4))
elseif (on_v_grid(qty)) then 
   lstatus = VMSK(lons(1), lats(1)) .or. &
             VMSK(lons(2), lats(2)) .or. &
             VMSK(lons(3), lats(3)) .or. &
             VMSK(lons(4), lats(4))
else 
   lstatus = TMSK(lons(1), lats(1)) .or. &
             TMSK(lons(2), lats(2)) .or. &
             TMSK(lons(3), lats(3)) .or. &
             TMSK(lons(4), lats(4))
endif   

end subroutine quad_on_land


!-----------------------------------------------------------------------
! Does any shutdown and clean-up needed for model.

subroutine depth_bounds(in_depth, all_depth, bot, top, frc, istatus)

real(r8), intent(in)  :: in_depth
real(r8), intent(in)  :: all_depth(Nz)
integer,  intent(out) :: bot, top
real(r8), intent(out) :: frc
integer,  intent(out) :: istatus

integer :: lo, hi, mid

istatus = 0

! Too shallow
if (in_depth <= all_depth(Nz)) then
   bot = Nz - 1
   top = Nz
   frc = 1.0_r8
   return
endif

! Too deep
if (in_depth > all_depth(1)) then
   bot = -1
   top = -1
   frc = 0.0_r8

   istatus = OBS_TOO_DEEP
   return
endif

! Binary search: find bot such that all_depth(bot) < in_depth <= all_depth(top)
lo = 1
hi = Nz

do while (hi - lo > 1)
   mid = (lo + hi) / 2
   if (in_depth > all_depth(mid)) then
      hi = mid
   else
      lo = mid
   endif
enddo

bot = lo
top = lo + 1

if (all_depth(top) == all_depth(bot)) then
   frc = 0.0_r8
else
   ! linear interpolation 
   ! This will not divide by zero
   frc = (in_depth - all_depth(bot)) / (all_depth(top) - all_depth(bot))
endif

end subroutine depth_bounds


!-----------------------------------------------------------------------
! Given a set of bottom and top levels for each ensemble member, 
! what is the list of unique vertical levels across the ensemble?

subroutine unique_levels(val1, val2, N, val_un, N_un)
  
  integer, intent(in)  :: N
  integer, intent(in)  :: val1(N), val2(N)
  integer, intent(out) :: val_un(2*N)
  integer, intent(out) :: N_un

  integer :: val(2*N)
  integer :: i, j
  logical :: is_new
  integer :: count

  val_un = 0

  ! Combine inputs
  val(1:N) = val1
  val(N+1:2*N) = val2

  count = 0
  do i = 1, 2*N
     is_new = .true.
     ! Check if val(i) is already in val_un
     do j = 1, count
        if (val(i) == val_un(j)) then
           is_new = .false.
           exit
        endif
     enddo
     ! Add to val_un if it's new
     if (is_new) then
        count = count + 1
        val_un(count) = val(i)
     endif
  enddo

  N_un = count
end subroutine unique_levels


!-----------------------------------------------------------------------
! Does any shutdown and clean-up needed for model.

subroutine end_model()

deallocate(ULAT, ULON, UDEP, UMSK)
deallocate(VLAT, VLON, VDEP, VMSK)
deallocate(TLAT, TLON, TDEP, TMSK)
deallocate(h, Cr, sr)  

end subroutine end_model


!===================================================================
! End of model_mod
!===================================================================
end module model_mod
