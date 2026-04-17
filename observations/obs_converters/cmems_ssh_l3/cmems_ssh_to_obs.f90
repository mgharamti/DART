! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

! Copernicus SSH converter. The data is based on the 
! "Global Ocean Along Track L3 Sea Surface Heights NRT" 
! available at: https://data.marine.copernicus.eu/product/SEALEVEL_GLO_PHY_L3_NRT_008_044/description

! The converter reads in the raw csv file and grabs the sla_filtered
! which is the low-pass filtered sea level anomaly, with DAC, ocean tide, 
! internal tide, and long-wavelength-error corrections applied. 
! This obs should be smoother, less noisy, and closer to the modeled
! mesoscale signal than its "sla_unfiltered" counterpart.

! The observation error standard deviation is parameterized as function 
! of depth such that it's small in the deep ocean and large over the shelf/
! shallow water. This is because in shallow/coastal regions, data tends 
! to be more contaminated from land, unresolved tides, internal tides, 
! wetting/drying effects, stronger small-scale variability, and 
! interpolation mismatch. To this end, the bathymetry is read in from 
! the model's grid file and used in the obs error function. 

! SSH data is in meters. 


program cmems_ssh_to_obs

use types_mod,            only : r8, MISSING_I, MISSING_R8
use time_manager_mod,     only : time_type, set_calendar_type, GREGORIAN, get_time,   &
                                 set_date, print_date, operator(+), operator(-)
use utilities_mod,        only : initialize_utilities, find_namelist_in_file,         &
                                 nmlfileunit, error_handler, do_nml_term, E_ERR,      &
                                 finalize_utilities, do_nml_file, get_next_filename,  &
                                 find_textfile_dims, file_exist, E_MSG, to_upper
use location_mod,         only : VERTISSURFACE
use obs_sequence_mod,     only : obs_type, obs_sequence_type, init_obs, get_num_obs,  &
                                 static_init_obs_sequence, init_obs_sequence,         &
                                 set_copy_meta_data, set_qc_meta_data, write_obs_seq, &
                                 destroy_obs_sequence, destroy_obs
use obs_utilities_mod,    only : create_3d_obs, add_obs_to_seq
use obs_kind_mod,         only : SATELLITE_SSH
use read_csv_mod,         only : csv_file_type, csv_get_nrows, csv_get_field,         &
                                 csv_open, csv_close, csv_print_header
use netcdf_utilities_mod, only : nc_check, nc_open_file_readonly, nc_close_file,      &
                                 nc_get_variable, nc_get_dimension_size

implicit none

character(len=*), parameter :: source = 'cmems_ssh_to_obs'

character(len=*), parameter :: SSH_OBS           = 'SLA_FILTERED'
integer,          parameter :: INFORMATION_LINES = 5

character(len=512) :: string1

! File variables
character(len=256) :: next_infile
integer            :: io, iunit, filenum
integer            :: num_new_obs, nfiles
integer            :: num_valid_obs
logical            :: first_obs = .true.

! Obs sequence
type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs, prev_obs
type(time_type)         :: prev_time
integer                 :: num_copies = 1, &     ! number of copies in sequence
                           num_qc     = 1        ! number of QC entries

! Data arrays
real(r8), allocatable   :: lon(:), lat(:) ! Location 
real(r8), allocatable   :: val(:)         ! Obs value
integer,  allocatable   :: vqc(:)         ! Obs QC 

real(r8), allocatable   :: glon(:, :), glat(:, :)  ! Ocean grid 
real(r8), allocatable   :: bath(:, :), mask(:, :)  ! Bathymetry and land mask
real(r8), allocatable   :: lon_1d(:) , lat_1d(:)

logical :: roms_fast_lookup

character(len=512), allocatable :: dat(:), par(:)

! SSH data structure
type altimetry_obs 
    type(time_type) :: dat 
    real(r8) :: lat 
    real(r8) :: lon 
    real(r8) :: obs
    real(r8) :: err
    real(r8) :: oqc 
end type altimetry_obs

type(altimetry_obs), allocatable :: sat(:)

! csv obs file 
type(csv_file_type) :: cf

!------------------------------------------------------------------------
!  Declare namelist parameters
character(len=256) :: file_list         = ''                 ! List of incoming obs csv files
character(len=256) :: ocean_in          = 'roms_restart.nc'  ! Ocean file for reading bathymetry
character(len=256) :: file_out          = 'obs_seq.ssh'      ! Resulting obs seq file
integer            :: avg_obs_per_file  = 500000             ! Average number of obs in file
real(r8)           :: obs_error_min     = 0.04_r8            ! Obs err sd [m]; minimum for open ocean conditions
real(r8)           :: obs_error_max     = 0.08_r8            ! Obs err sd [m]; maximum for shallow water
real(r8)           :: transition_depth  = 500.0_r8           ! Bathymetric transition depth
logical            :: debug             = .true.             ! verbose output


namelist /cmems_ssh_to_obs_nml/ file_list,        &
                                ocean_in,         &
                                file_out,         &
                                obs_error_min,    &
                                obs_error_max,    & 
                                transition_depth, &
                                avg_obs_per_file, &
                                debug

! Start Converter
call initialize_utilities()

! Read the namelist options
call find_namelist_in_file('input.nml', 'cmems_ssh_to_obs_nml', iunit)
read(iunit, nml = cmems_ssh_to_obs_nml, iostat = io)

if (do_nml_file()) write(nmlfileunit, nml=cmems_ssh_to_obs_nml)
if (do_nml_term()) write(     *     , nml=cmems_ssh_to_obs_nml)

! Set the calendar kind
call set_calendar_type(GREGORIAN)

! Get number of observations
num_new_obs = avg_obs_per_file

call find_textfile_dims(file_list, nfiles)
num_new_obs = avg_obs_per_file * nfiles

! Initialize obs seq file
call static_init_obs_sequence()
call init_obs(obs, num_copies, num_qc)
call init_obs(prev_obs, num_copies, num_qc)

if (file_exist(file_out)) then
   write(*, '(/, A)') 'Output file: '//trim(adjustl(file_out))//' exists. Replacing it ...'
else
   write(*, '(/, A)') 'Creating "'//trim(adjustl(file_out))//'" file.'
endif

call init_obs_sequence(obs_seq, num_copies, num_qc, num_new_obs)
call set_copy_meta_data(obs_seq, num_copies, 'SSH observation')
call set_qc_meta_data(obs_seq, num_qc, 'SSH QC')

! Read the ocean grid, mask, and bathymetry
call read_ocean(ocean_in)

! Loop over the obs files 
do filenum = 1, nfiles
   ! Retrieve the obs file 
   next_infile = get_next_filename(file_list, filenum)
   if (len_trim(next_infile) == 0) exit 
   if (debug) write(*, '(/, 2X, A, i0, X, A)') 'Input file: #', filenum, trim(next_infile)

   ! Read data from file
   call get_satellite_data(next_infile)

   ! Add obs to sequence
   call fill_obs(obs, prev_obs, prev_time) 

   ! Release memory
   call cleanup()
enddo

call finish_obs()
call finalize_utilities()


contains


!-----------------------------------------------
! Read SSH, date, and metadata from file  
subroutine get_satellite_data(filename)

character(len=*), intent(in) :: filename
character(len=*), parameter  :: routine = 'get_sat_data'

integer :: k, i, nobs

character(len=64) :: obs_name

! Open csv file and get dims
call csv_open(filename, cf, skip_lines=INFORMATION_LINES, context=routine)
nobs = csv_get_nrows(cf)

if (debug) call csv_print_header(cf)

allocate(dat(nobs), sat(nobs))
allocate(lat(nobs), lon(nobs))
allocate(val(nobs), vqc(nobs), par(nobs))

! Read the data
call csv_get_field(cf, 'time'     , dat, routine)
call csv_get_field(cf, 'longitude', lon, routine)
call csv_get_field(cf, 'latitude' , lat, routine)
call csv_get_field(cf, 'value'    , val, routine)
call csv_get_field(cf, 'valueQc'  , vqc, routine)
call csv_get_field(cf, 'parameter', par, routine)

call csv_close(cf)

k = 0

! Cleanup the data 
do i = 1, nobs

   if (lat(i) == MISSING_R8 .or. &
       lon(i) == MISSING_R8 .or. &
       val(i) == MISSING_R8 .or. &
       vqc(i) == MISSING_I) cycle
   
   ! Get rid of any non-ssh/bad obs
   obs_name = trim(par(i))
   call to_upper(obs_name)

   if (obs_name /= SSH_OBS .or. vqc(i) /= 0) cycle 

   if (lon(i) < 0.0_r8) lon(i) = lon(i) + 360.0_r8

   k = k + 1

   sat(k)%lon = lon(i)
   sat(k)%lat = lat(i)
   sat(k)%dat = parse_time(dat(i))
   sat(k)%obs = val(i)
   sat(k)%err = ssh_obs_error(lon(i), lat(i))
   sat(k)%oqc = vqc(i)

enddo
num_valid_obs = k

! Get rid of unused entries
if (num_valid_obs < nobs) sat = sat(1:num_valid_obs)

if (debug) then
   ! print some of the obs to check
   print * 
   do i = 1, min(num_valid_obs, 10)
      write(string1, '(4X, A, i6, 4(A,f10.4), A)') '* obs #', i, ', lat:'  , sat(i)%lat, &
                                          ', lon:', sat(i)%lon , ', SSH:'  , sat(i)%obs, &
                                          ', QC:' , sat(i)%oqc , ', date: '
      call print_date(sat(i)%dat, str=string1)
   enddo
endif

end subroutine get_satellite_data  


!------------------------------------------------------
! Add the metadata in the obs sequence  
subroutine fill_obs(obs, prev_obs, prev_time)

type(obs_type),        intent(inout) :: obs, prev_obs
type(time_type),       intent(inout) :: prev_time

type(time_type) :: odat
real(r8)        :: olon, olat, oval, oerr
integer         :: iobs, osec, oday

do iobs = 1, num_valid_obs
   olon = sat(iobs)%lon
   olat = sat(iobs)%lat
   odat = sat(iobs)%dat
   oval = sat(iobs)%obs
   oerr = sat(iobs)%err

   call get_time(odat, osec, oday)
   call create_3d_obs(olat, olon, 0.0_r8, VERTISSURFACE, oval, SATELLITE_SSH, oerr, oday, osec, 0.0_r8, obs)
   call add_obs_to_seq(obs_seq, obs, odat, prev_obs, prev_time, first_obs)
enddo

end subroutine fill_obs


!---------------------------------------------
! Use local bathymetry depth and define an obs 
! error sd that increases as water gets shallower
! deep ocean          -> smaller SSH error
! shelf/shallow water -> larger  SSH error
!
! sigma = sigma_min + (sigma_max - sigma_min) * exp(-h/h0)
! where:
!    - sigma    : error sd
!    - sigma_min: open ocean minimum error
!    - sigma_max: shallow-water maximum error
!    - h        : local bathymetry
!    - h0       : transition depth scale 
function ssh_obs_error(lon, lat) result(error)

real(r8), intent(in) :: lon, lat
real(r8)             :: sigma_min, sigma_max, h0
real(r8)             :: error, h_obs

sigma_min = obs_error_min
sigma_max = obs_error_max
h0        = transition_depth

! Get bathymetry
h_obs = depth_at_obs(lon, lat)

if (h_obs /= h_obs .or. h_obs <= 0.0_r8) then
   error = sigma_max
   return
endif

error = sigma_min + (sigma_max - sigma_min) * exp(-h_obs / h0)

if (1>2) then 
   write(*, '(4(A, F13.6))') 'lon:', lon  , ', lat:', lat, & 
                    ', bathymetry:', h_obs, ', err:', error 
endif

end function ssh_obs_error


!------------------------------------------------------------------
! Given an obs location, get the nearest estimate of the bathymetry
function depth_at_obs(lon_obs, lat_obs) result(h_obs)

real(r8), intent(in) :: lon_obs, lat_obs

real(r8) :: h_obs

! Fast lookup for ROMS's regular grid 
if (roms_fast_lookup) then
   h_obs = depth_roms_fast(lon_obs, lat_obs)
else
   ! Other grids: brute force
   h_obs = depth_bruteforce(lon_obs, lat_obs)
endif

end function depth_at_obs


!-----------------------------------------------
! Find depth using a pair of 1D location arrays
function depth_roms_fast(lon_obs, lat_obs) result(h_obs)

real(r8), intent(in) :: lon_obs, lat_obs

integer  :: i0, j0, i, j, i1, i2, j1, j2
real(r8) :: dlon, dlat, d2, d2min, h_obs

! find nearest indices in 1D
i0 = nearest_index_1d(lon_obs, lon_1d)
j0 = nearest_index_1d(lat_obs, lat_1d)

! local search window (within 2 cells)
! Idea is that instead of searching the 
! entire grid, I'll look only in a local
! neighborhood around the candidate index, i.e., i0, j0
! My window: A 5x5 patch (25 cells) + Clamp the edges 
! i0-2 -> i0+2
! j0-2 -> j0+2
i1 = max(1, i0-2)
i2 = min(size(glon,1), i0+2)
j1 = max(1, j0-2)
j2 = min(size(glon,2), j0+2)

d2min = huge(1.0_r8)
h_obs = MISSING_R8

do j = j1, j2
   do i = i1, i2
      if (mask(i,j) == 0.0_r8) cycle

      dlon = abs(glon(i,j) - lon_obs)
      dlon = min(dlon, 360.0_r8 - dlon) ! in case of a wrap-around

      dlat = glat(i,j) - lat_obs
      d2   = dlon*dlon + dlat*dlat

      ! Compare current distance with the best so far "d2min"
      ! If smaller, we need to update the min and store the 
      ! corresponding depth "h_obs"
      if (d2 == d2 .and. d2 < d2min) then
         d2min = d2
         h_obs = bath(i,j)
      endif
   enddo
enddo

end function depth_roms_fast


!-----------------------------------------------
! Scan the entire grid to find approximate depth
function depth_bruteforce(lon_obs, lat_obs) result(h_obs)

real(r8), intent(in) :: lon_obs, lat_obs

integer  :: i, j
real(r8) :: dlon, dlat, d2, d2min, h_obs

d2min = huge(1.0_r8)
h_obs = MISSING_R8

do j = 1, size(glon,2)
   do i = 1, size(glon,1)
      if (mask(i,j) == 0.0_r8) cycle

      dlon = abs(glon(i,j) - lon_obs)
      dlon = min(dlon, 360.0_r8 - dlon)

      dlat = glat(i,j) - lat_obs
      d2   = dlon*dlon + dlat*dlat

      if (d2 == d2 .and. d2 < d2min) then
         d2min = d2
         h_obs = bath(i,j)
      endif
   enddo
enddo

end function depth_bruteforce


!----------------------------------------
! Nearest index in a 1D array to a number
integer function nearest_index_1d(x, arr) result(idx)

real(r8), intent(in) :: x
real(r8), intent(in) :: arr(:)

integer  :: k
real(r8) :: d, dmin

dmin = huge(1.0_r8)
idx  = 1

do k = 1, size(arr)
   d = abs(arr(k) - x)
   if (d < dmin) then
      dmin = d
      idx = k
   endif
enddo

end function nearest_index_1d


!-------------------------------------------
! Check whether lon/lat arrays are monotonic 
logical function is_monotonic(arr)

real(r8), intent(in) :: arr(:)
integer :: k

is_monotonic = .true.

do k = 2, size(arr)
   if (arr(k) <= arr(k-1)) then
      is_monotonic = .false.
      return
   endif
enddo

end function is_monotonic


!-----------------------------------------
! Using the date string, set the DART time
function parse_time(datestr) result(obs_time)
   
character(len=*), parameter :: routine = 'parse_time'
   
character(len=*), intent(in) :: datestr
type(time_type)              :: obs_time

integer :: year, month, day, hour, minute, second
   
!example: 2026-04-16T08:43:54.000Z                      
read(datestr(1:4),  *) year
read(datestr(6:7),  *) month
read(datestr(9:10), *) day
read(datestr(12:13),*) hour
read(datestr(15:16),*) minute
read(datestr(18:19),*) second

obs_time = set_date(year, month, day, hour, minute, second)
                                            
end function parse_time 


!--------------------------------------------------
! Read ocean grid and bathymetry from the relavant 
! ocean model 
subroutine read_ocean(file)

character(len=*), parameter :: routine = 'read_ocean'

character(len=*), intent(in) :: file
integer :: ncid, nx, ny

! open Ocean file 
ncid = nc_open_file_readonly(file, routine)

! We need to read the bathymetry from the ocean file. 
! In this case, ROMS is used. The following code can 
! be adjusted for other ocean configurations. 
nx = nc_get_dimension_size(ncid, 'xi_rho',  routine)
ny = nc_get_dimension_size(ncid, 'eta_rho', routine)

allocate(bath(nx, ny), glon(nx, ny), &
         glat(nx, ny), mask(nx, ny), &
         lon_1d(nx)  , lat_1d(ny))

call nc_get_variable(ncid, 'h'       , bath, routine)
call nc_get_variable(ncid, 'mask_rho', mask, routine)
call nc_get_variable(ncid, 'lon_rho' , glon, routine)
call nc_get_variable(ncid, 'lat_rho' , glat, routine)

! Only possible because the ROMS grid in this case
! is regular but not equally spaced. 
where (glon < 0.0_r8) glon = glon + 360.0_r8

lon_1d = glon(:,1)
lat_1d = glat(1,:)

roms_fast_lookup = is_monotonic(lon_1d) .and. is_monotonic(lat_1d)

call nc_close_file(ncid, source)  

end subroutine read_ocean


!------------------------------------------------------------
! All collected obs (if any) are written in the seq file. 
subroutine finish_obs()

integer :: obs_num

obs_num = get_num_obs(obs_seq)

if (obs_num > 0) then
   if (debug) write(*, '(/, A, i0, A)') '> Ready to write ', obs_num, ' observations:'

   call write_obs_seq(obs_seq, file_out)
   call destroy_obs(obs)
   call destroy_obs_sequence(obs_seq)
else
   string1 = 'No obs were converted.'
   call error_handler(E_ERR, source, string1)
endif

deallocate(glon, glat, bath, mask)
deallocate(lon_1d, lat_1d)

call error_handler(E_MSG, source, 'Finished successfully.')

end subroutine finish_obs


!------------------------------------------------------------
! Clear up memory 
subroutine cleanup()

deallocate(lon, lat, val, vqc)
deallocate(dat, par, sat)

end subroutine cleanup

end program cmems_ssh_to_obs
