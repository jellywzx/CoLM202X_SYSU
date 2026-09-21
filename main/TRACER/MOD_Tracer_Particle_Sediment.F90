#include <define.h>

! Sediment is a river-routing tracer: it hard-USEs MOD_Grid_RiverLakeNetwork and
! other GridRiverLakeFlow-only modules, so it can only be compiled when routing
! is compiled in. Guarding on TRACER alone left a dangling USE when
! GridRiverLakeFlow was off (e.g. SinglePoint); the register manifest is guarded
! to match (see include/tracer_lifecycle_providers.inc).
#if (defined TRACER) && (defined GridRiverLakeFlow)
MODULE MOD_Tracer_Particle_Sediment
!-------------------------------------------------------------------------------------
! DESCRIPTION:
!
!   Sediment transport module for GridRiverLakeFlow.
!   Ported from CaMa-Flood sediment module (CoLM-sed-master).
!
!   Physics:
!     - Suspended sediment advection (upstream scheme)
!     - Bedload transport (Ashida-Michiue shear-velocity form; often grouped with
!       Meyer-Peter-Mueller-type bedload relations)
!     - Suspension-deposition exchange (Uchida & Fukuoka 2019)
!     - Hillslope erosion (precipitation-driven, Sunada & Hasegawa 1993)
!     - Vertical bed layer redistribution
!
!   References:
!     - Uchida & Fukuoka (2019) suspension velocity formula
!     - Egiazoroff equation for mixed-size critical shear
!     - Shields curve for critical shear stress
!
!-------------------------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_Namelist, only: DEF_UnitCatchment_file, DEF_USE_BIFURCATION, DEF_USE_LEVEE, &
      DEF_hist_vars
   USE MOD_Vars_Global, only: spval
   USE, INTRINSIC :: IEEE_ARITHMETIC, only: ieee_is_finite
   USE, INTRINSIC :: IEEE_EXCEPTIONS, only: ieee_get_halting_mode, ieee_set_halting_mode, &
      ieee_invalid, ieee_divide_by_zero, ieee_overflow
   IMPLICIT NONE
   PRIVATE

   !-------------------------------------------------------------------------------------
   ! Module Parameters
   !-------------------------------------------------------------------------------------
   integer,  save :: nsed           ! Number of sediment size classes
   integer,  save :: totlyrnum      ! Number of deposition layers
   integer,  save :: nlfp_sed       ! Number of floodplain layers for slope

   real(r8), save :: lambda         ! Porosity [-]
   real(r8), save :: lyrdph         ! Active layer depth [m]
   real(r8), save :: psedD          ! Internal sediment density [g/cm3]
   real(r8), save :: pwatD          ! Internal water density [g/cm3]
   real(r8), save :: visKin         ! Kinematic viscosity [m2/s]
   real(r8), save :: vonKar         ! Von Karman coefficient [-]
   real(r8), save :: pset           ! Settling velocity multiplier [-]

   ! Sediment yield parameters
   real(r8), save :: pyld           ! Yield coefficient
   real(r8), save :: pyldc          ! Slope exponent
   real(r8), save :: pyldpc         ! Precipitation exponent
   real(r8), save :: dsylunit       ! Unit conversion factor
   real(r8), save :: sed_ignore_dph ! Minimum water depth for active suspended-sediment processes [m]
   real(r8), save :: sed_cfl_adv    ! CFL factor for suspended-sediment advection [-]
   real(r8), save :: sed_dt_max     ! Maximum sediment substep [s]
   real(r8), save :: sed_bed_depth  ! Initial named deposit-bed depth [m]
   real(r8), save :: sed_max_conc   ! Maximum total suspended sediment volume concentration [-]

   real(r8), parameter :: SED_DEFAULT_LAMBDA = 0.4_r8
   real(r8), parameter :: SED_DEFAULT_LYRDPH = 0.05_r8
   real(r8), parameter :: SED_DEFAULT_DENSITY = 2.65_r8
   real(r8), parameter :: SED_DEFAULT_WATER_DENSITY = 1.0_r8
   real(r8), parameter :: SED_DEFAULT_VISKIN = 1.0e-6_r8
   real(r8), parameter :: SED_DEFAULT_VONKAR = 0.4_r8
   real(r8), parameter :: SED_DEFAULT_PSET = 1.0_r8
   integer,  parameter :: SED_DEFAULT_TOTLYRNUM = 5
   real(r8), parameter :: SED_DEFAULT_CFL_ADV = 0.5_r8
   real(r8), parameter :: SED_DEFAULT_IGNORE_DPH = 0.05_r8
   real(r8), parameter :: SED_DEFAULT_DT_MAX = 3600._r8
   integer,  parameter :: SED_MAX_ADV_SUBSTEPS = 100000
   integer,  parameter :: SED_RESTART_SCHEMA_VERSION = 1
   real(r8), parameter :: SED_DEFAULT_BED_DEPTH = 10._r8
   character(len=*), parameter :: SED_DEFAULT_DIAMETER = '0.0002,0.002,0.02'
   real(r8), parameter :: SED_DEFAULT_PYLD = 0.01_r8
   real(r8), parameter :: SED_DEFAULT_PYLDC = 2.0_r8
   real(r8), parameter :: SED_DEFAULT_PYLDPC = 2.0_r8
   real(r8), parameter :: SED_DEFAULT_DSYLUNIT = 1.0e-6_r8

   real(r8), parameter :: SED_DEFAULT_MAX_CONC = 0.1_r8  ! Maximum sediment concentration (10% by volume, matches CoLM-sed-master)
   real(r8), parameter :: SED_BEDLOAD_COEFF = 17._r8
   real(r8), parameter :: SED_PRECIP_THRESHOLD_MM_DAY = 2._r8
   real(r8), parameter :: EXCH_SHEARVEL_MIN = 1.e-4_r8
   real(r8), parameter :: EXCH_SHEARVEL_BLEND = 2._r8 * EXCH_SHEARVEL_MIN
   real(r8), parameter :: EXCH_ZD_MAX = 100._r8
   real(r8), parameter :: SED_BALANCE_ABS_TOL = 1.e-10_r8
   real(r8), parameter :: SED_BALANCE_REL_TOL = 1.e-10_r8
   ! Numerical near-dry threshold for sediment carrier transport.
   ! Shallower water is treated as a residual numerical film and must not
   ! generate sediment advection CFL constraints.
   real(r8), parameter :: SED_NEAR_DRY_DEPTH = 1.e-4_r8

#ifdef CoLMDEBUG
   ! Fixed Amazon benchmark stations from the HYBAM reference network.
   ! Coordinates are decimal degrees (north/east positive).  The routing cell
   ! match is resolved once at initialization, then reused every routing period.
   integer, parameter :: SED_N_DIAG_STATIONS = 4
   character(len=16), parameter :: SED_DIAG_STATION_NAMES(SED_N_DIAG_STATIONS) = (/ &
      'Serrinha        ', 'Porto Velho     ', 'Manacapuru      ', 'Obidos          ' /)
   real(r8), parameter :: SED_DIAG_STATION_LAT(SED_N_DIAG_STATIONS) = (/ &
      -0.4500_r8, -8.7370_r8, -3.3122_r8, -1.9470_r8 /)
   real(r8), parameter :: SED_DIAG_STATION_LON(SED_N_DIAG_STATIONS) = (/ &
      -64.8300_r8, -63.9200_r8, -60.6303_r8, -55.5110_r8 /)
   integer, save :: sed_diag_station_ucid(SED_N_DIAG_STATIONS) = 0
   integer, allocatable, save :: sed_diag_station_local_i(:)
   real(r8), save :: sed_diag_station_model_lon(SED_N_DIAG_STATIONS) = 0._r8
   real(r8), save :: sed_diag_station_model_lat(SED_N_DIAG_STATIONS) = 0._r8
   real(r8), save :: sed_diag_station_distance_km(SED_N_DIAG_STATIONS) = 0._r8
#endif

   !-------------------------------------------------------------------------------------
   ! Static Data (read from DEF_UnitCatchment_file)
   !-------------------------------------------------------------------------------------
   real(r8), allocatable :: sed_frc   (:,:)    ! Sediment fraction [nsed, numucat]
   real(r8), allocatable :: sed_slope (:,:)    ! Floodplain slope [nlfp_sed, numucat]
   real(r8), allocatable :: sDiam     (:)      ! Grain diameter [nsed]
   real(r8), allocatable :: setvel    (:)      ! Settling velocity [nsed]
   real(r8), allocatable :: sDiam_from_param(:)! Grain diameters supplied by DEF_TRACER_PARAM_FILES

   !-------------------------------------------------------------------------------------
   ! State Variables
   !-------------------------------------------------------------------------------------
   real(r8), allocatable :: sedcon  (:,:)      ! Suspended sediment concentration [nsed, numucat]
   real(r8), allocatable :: sedsto  (:,:)      ! Suspended solid volume [m3, nsed, numucat]
   real(r8), allocatable :: layer   (:,:)      ! Active layer storage [nsed, numucat]
   real(r8), allocatable :: seddep  (:,:,:)    ! Deposition layer storage [nsed, totlyrnum, numucat]

   !-------------------------------------------------------------------------------------
   ! Diagnostic Variables
   !-------------------------------------------------------------------------------------
   real(r8), allocatable :: sedout  (:,:)      ! Suspended sediment outflow [nsed, numucat]
   real(r8), allocatable :: bedout  (:,:)      ! Bedload solid-volume flux [nsed, numucat]
   real(r8), allocatable :: sedinp  (:,:)      ! Erosion input [nsed, numucat]
   real(r8), allocatable :: netflw  (:,:)      ! Net bed-water exchange flux [nsed, numucat]
                                                ! Includes suspension-deposition exchange (Es-D)
                                                ! and shallow-cell erosion input deposited directly
                                                ! into bed layer.  Positive = net entrainment.
   real(r8), allocatable :: exch_es_raw(:,:)   ! Raw entrainment flux Es from exchange formula [nsed, numucat]
   real(r8), allocatable :: exch_d_raw (:,:)   ! Raw deposition flux D from exchange formula [nsed, numucat]
   real(r8), allocatable :: exch_es_eff(:,:)   ! Effective entrainment flux applied after limits [nsed, numucat]
   real(r8), allocatable :: exch_d_eff (:,:)   ! Effective deposition flux applied after limits [nsed, numucat]
   real(r8), allocatable :: netflw_adv_step(:,:) ! Current CFL-substep cap/dry deposition [nsed, numucat]
   real(r8), allocatable :: exch_d_adv_step(:,:) ! Current CFL-substep effective deposition [nsed, numucat]
   real(r8), allocatable :: shearvel(:)        ! Shear velocity [numucat]
   real(r8), allocatable :: critshearvel(:,:)  ! Critical shear velocity [nsed, numucat]
   real(r8), allocatable :: susvel  (:,:)      ! Suspension velocity [nsed, numucat]

   !-------------------------------------------------------------------------------------
   ! Accumulated Variables for Sediment Time-stepping (per-cell)
   !-------------------------------------------------------------------------------------
   real(r8), allocatable :: sed_acc_time (:)   ! Accumulated time [numucat]
   real(r8), allocatable :: sed_acc_v2   (:)   ! Accumulated velocity**2 * dt [numucat]
   real(r8), allocatable :: sed_acc_wdsrf(:)   ! Accumulated water depth*dt [numucat]
   real(r8), allocatable :: sed_acc_rivsto(:)  ! Accumulated routed water storage*dt [m3 s]
   real(r8), allocatable :: sed_acc_rivout(:)  ! Accumulated discharge*dt [numucat]
   real(r8), allocatable :: sed_acc_abs_rivout(:) ! Accumulated abs(discharge)*dt [numucat]
   real(r8), allocatable :: sed_acc_floodarea(:) ! Accumulated flood area*dt [numucat]
   real(r8), allocatable :: sed_acc_carrier_time(:) ! Time with depth above sediment carrier threshold [s]
   real(r8), allocatable :: sed_acc_wdsrf_min(:), sed_acc_wdsrf_max(:) ! Period min/max water depth [m]
   real(r8), allocatable :: sed_acc_rivsto_min(:), sed_acc_rivsto_max(:) ! Period min/max carrier storage [m3]
   real(r8), allocatable :: sed_acc_rivout_min(:), sed_acc_rivout_max(:) ! Period min/max face discharge [m3/s]
   real(r8), allocatable :: sed_acc_pos_rivout(:), sed_acc_neg_rivout(:) ! Time-integrated +/- discharge [m3]
   real(r8), allocatable :: sed_acc_near_dry_abs_rivout(:) ! Abs discharge filtered by instantaneous near-dry gate [m3]
   real(r8), allocatable :: sed_precip(:)      ! Accumulated precipitation [mm, for diagnostics]
   real(r8), allocatable :: sed_precip_yield(:) ! Accumulated (rate_mm_hr)^pyldpc * dt [numucat]
                                                ! Pre-computed per forcing step to avoid Jensen bias
   real(r8), save        :: sed_precip_time    ! Accumulated precipitation time [s]

   !-------------------------------------------------------------------------------------
   ! Accumulated Variables for History Output
   !-------------------------------------------------------------------------------------
   real(r8), save        :: sed_hist_acctime   ! Module-private total time for history averaging [s]
   real(r8), allocatable :: a_sedcon  (:,:)    ! Accumulated sedcon
   real(r8), allocatable :: a_sedout  (:,:)    ! Accumulated sedout
   real(r8), allocatable :: a_bedout  (:,:)    ! Accumulated bedload solid-volume flux
   real(r8), allocatable :: a_sedinp  (:,:)    ! Accumulated sedinp
   real(r8), allocatable :: a_netflw  (:,:)    ! Accumulated netflw
   real(r8), allocatable :: a_layer   (:,:)    ! Accumulated layer
   real(r8), allocatable :: a_shearvel(:)      ! Accumulated shearvel

   !-------------------------------------------------------------------------------------
   ! Public Subroutines
   !-------------------------------------------------------------------------------------
   PUBLIC :: register_sediment_tracer_provider

   integer, parameter :: MAX_SED_PARAM_CLASSES = 100
   type :: sediment_parameter_type
      integer  :: nsed = -1
      real(r8) :: grain_diameter(MAX_SED_PARAM_CLASSES) = -1._r8
      real(r8) :: grain_density = -1._r8
      real(r8) :: water_density = -1._r8
      real(r8) :: porosity = -1._r8
      integer  :: ndeposit_layers = -1
      real(r8) :: ignore_depth_m = -1._r8
      real(r8) :: active_layer_depth = -1._r8
      real(r8) :: viscosity = -1._r8
      real(r8) :: von_karman = -1._r8
      real(r8) :: settling_multiplier = -1._r8
      real(r8) :: yield_coefficient = -1._r8
      real(r8) :: slope_exponent = -1._r8
      real(r8) :: precipitation_exponent = -1._r8
      real(r8) :: unit_conversion = -1._r8
      real(r8) :: cfl_adv = -1._r8
      real(r8) :: max_timestep_s = -1._r8
      real(r8) :: bed_depth = -1._r8
      real(r8) :: max_concentration = -1._r8
   end type sediment_parameter_type

   integer, save :: sediment_itrc = 0
CONTAINS

   !-------------------------------------------------------------------------------------
   logical FUNCTION sediment_particle_enabled()
   !-------------------------------------------------------------------------------------
      IMPLICIT NONE
      sediment_particle_enabled = sediment_itrc > 0
   END FUNCTION sediment_particle_enabled

   !-------------------------------------------------------------------------------------
   integer FUNCTION sediment_tracer_index()
   !-------------------------------------------------------------------------------------
   IMPLICIT NONE
      sediment_tracer_index = sediment_itrc

   END FUNCTION sediment_tracer_index

   !-------------------------------------------------------------------------------------
   SUBROUTINE register_sediment_tracer_provider()
   !-------------------------------------------------------------------------------------
      USE MOD_Tracer_Defs, only: FAMILY_PARTICLE, STATE_OWNER_PROVIDER, REACTION_NONE
      USE MOD_Tracer_Lifecycle, only: tracer_lifecycle_hooks_type, register_tracer_provider
      IMPLICIT NONE
      type(tracer_lifecycle_hooks_type) :: hooks

      hooks = tracer_lifecycle_hooks_type()
      hooks%route_init            => grid_sediment_init
      hooks%route_read_restart    => read_sediment_restart
      hooks%route_forcing_put     => sediment_forcing_put
      hooks%route_diag_accumulate => sediment_diag_accumulate
      hooks%route_calc            => grid_sediment_calc
      hooks%route_history         => write_sediment_history
      hooks%route_flush_history   => flush_sediment_history
      hooks%route_write_restart   => write_sediment_restart
      hooks%route_final           => grid_sediment_final

      CALL register_tracer_provider('SEDIMENT', 'SED', 'sediment', &
         FAMILY_PARTICLE, STATE_OWNER_PROVIDER, REACTION_NONE, hooks, sediment_itrc)

   END SUBROUTINE register_sediment_tracer_provider

#ifdef CoLMDEBUG
   !-------------------------------------------------------------------------------------
   SUBROUTINE initialize_sediment_diag_stations()
   ! Resolve each fixed HYBAM station to one deterministic nearest routing unit.
   !-------------------------------------------------------------------------------------
   USE MOD_Grid_RiverLakeNetwork, only: numucat, ucat_ucid, x_ucat, y_ucat, griducat
   USE MOD_Vars_Global, only: pi
   IMPLICIT NONE

   integer :: ista, i, best_i, candidate_ucid
   real(r8) :: cell_lon, cell_lat, dlon, dlat, d2, best_d2, global_d2
   real(r8) :: model_lon_local, model_lat_local
   real(r8) :: tie_tol

      IF (.not. p_is_worker) RETURN

      IF (allocated(sed_diag_station_local_i)) deallocate(sed_diag_station_local_i)
      allocate(sed_diag_station_local_i(SED_N_DIAG_STATIONS))
      sed_diag_station_local_i = 0
      sed_diag_station_ucid = 0
      sed_diag_station_model_lon = 0._r8
      sed_diag_station_model_lat = 0._r8
      sed_diag_station_distance_km = 0._r8

      DO ista = 1, SED_N_DIAG_STATIONS
         best_d2 = huge(1._r8)
         best_i = 0

         DO i = 1, numucat
            cell_lon = -180._r8 + (real(x_ucat(i),r8) - 0.5_r8) &
               * 360._r8 / real(griducat%nlon,r8)
            cell_lat = 90._r8 - (real(y_ucat(i),r8) - 0.5_r8) &
               * 180._r8 / real(griducat%nlat,r8)
            dlon = (cell_lon - SED_DIAG_STATION_LON(ista)) &
               * cos(SED_DIAG_STATION_LAT(ista) * pi / 180._r8)
            dlat = cell_lat - SED_DIAG_STATION_LAT(ista)
            d2 = dlon * dlon + dlat * dlat
            IF (d2 < best_d2) THEN
               best_d2 = d2
               best_i = i
            ELSE IF (d2 == best_d2 .and. best_i > 0) THEN
               IF (ucat_ucid(i) < ucat_ucid(best_i)) THEN
                  best_d2 = d2
                  best_i = i
               ENDIF
            ENDIF
         ENDDO

         global_d2 = best_d2
#ifdef USEMPI
         CALL mpi_allreduce(MPI_IN_PLACE, global_d2, 1, MPI_REAL8, MPI_MIN, p_comm_worker, p_err)
#endif

         candidate_ucid = huge(1)
         tie_tol = max(1.e-12_r8, 1.e-10_r8 * max(global_d2, 1._r8))
         IF (best_i > 0 .and. abs(best_d2 - global_d2) <= tie_tol) THEN
            candidate_ucid = ucat_ucid(best_i)
         ENDIF
#ifdef USEMPI
         CALL mpi_allreduce(MPI_IN_PLACE, candidate_ucid, 1, MPI_INTEGER, MPI_MIN, p_comm_worker, p_err)
#endif
         sed_diag_station_ucid(ista) = candidate_ucid

         DO i = 1, numucat
            IF (ucat_ucid(i) == candidate_ucid) THEN
               sed_diag_station_local_i(ista) = i
               EXIT
            ENDIF
         ENDDO

         model_lon_local = 0._r8
         model_lat_local = 0._r8
         IF (sed_diag_station_local_i(ista) > 0) THEN
            i = sed_diag_station_local_i(ista)
            model_lon_local = -180._r8 + (real(x_ucat(i),r8) - 0.5_r8) &
               * 360._r8 / real(griducat%nlon,r8)
            model_lat_local = 90._r8 - (real(y_ucat(i),r8) - 0.5_r8) &
               * 180._r8 / real(griducat%nlat,r8)
         ENDIF
#ifdef USEMPI
         CALL mpi_allreduce(MPI_IN_PLACE, model_lon_local, 1, MPI_REAL8, MPI_SUM, p_comm_worker, p_err)
         CALL mpi_allreduce(MPI_IN_PLACE, model_lat_local, 1, MPI_REAL8, MPI_SUM, p_comm_worker, p_err)
#endif
         sed_diag_station_model_lon(ista) = model_lon_local
         sed_diag_station_model_lat(ista) = model_lat_local
         sed_diag_station_distance_km(ista) = sqrt(max(global_d2,0._r8)) * 111.195_r8
      ENDDO

      IF (p_iam_worker == 0) THEN
         WRITE(*,'(A)') 'Sediment Amazon benchmark station mapping:'
         DO ista = 1, SED_N_DIAG_STATIONS
            WRITE(*,'(2X,A16,A,2(F10.4,1X),A,I0,A,2(F10.4,1X),A,F9.3)') &
               trim(SED_DIAG_STATION_NAMES(ista)), ' target(lat lon)=', &
               SED_DIAG_STATION_LAT(ista), SED_DIAG_STATION_LON(ista), &
               ' ucat=', sed_diag_station_ucid(ista), ' model(lat lon)=', &
               sed_diag_station_model_lat(ista), sed_diag_station_model_lon(ista), &
               ' distance_km=', sed_diag_station_distance_km(ista)
         ENDDO
      ENDIF

   END SUBROUTINE initialize_sediment_diag_stations
#endif

   !-------------------------------------------------------------------------------------
   SUBROUTINE grid_sediment_init()
   !-------------------------------------------------------------------------------------
   USE netcdf
   USE MOD_NetCDFSerial
   USE MOD_Grid_RiverLakeNetwork, only: numucat, totalnumucat, &
      ucat_data_address, topo_rivwth, topo_rivlen
   IMPLICIT NONE

   character(len=256) :: parafile
   integer :: ncid, dimid, ierr
   integer :: i

      IF (.not. sediment_particle_enabled()) RETURN

      IF (p_is_io) THEN
         WRITE(*,*) 'Initializing sediment module...'
      ENDIF

      ! Set species-local defaults.  Overrides come from the tracer parameter
      ! file selected with DEF_TRACER_PARAM_FILES = 'SEDIMENT:<path>'.
      lambda    = SED_DEFAULT_LAMBDA
      lyrdph    = SED_DEFAULT_LYRDPH
      psedD     = SED_DEFAULT_DENSITY
      pwatD     = SED_DEFAULT_WATER_DENSITY
      visKin    = SED_DEFAULT_VISKIN
      vonKar    = SED_DEFAULT_VONKAR
      pset      = SED_DEFAULT_PSET
      totlyrnum = SED_DEFAULT_TOTLYRNUM
      pyld      = SED_DEFAULT_PYLD
      pyldc     = SED_DEFAULT_PYLDC
      pyldpc    = SED_DEFAULT_PYLDPC
      dsylunit  = SED_DEFAULT_DSYLUNIT
      sed_ignore_dph = SED_DEFAULT_IGNORE_DPH
      sed_cfl_adv = SED_DEFAULT_CFL_ADV
      sed_dt_max = SED_DEFAULT_DT_MAX
      sed_bed_depth = SED_DEFAULT_BED_DEPTH
      sed_max_conc = SED_DEFAULT_MAX_CONC

      parafile = DEF_UnitCatchment_file

      ! Read dimensions directly from NetCDF dimension names
      IF (p_is_master) THEN
         ierr = nf90_open(trim(parafile), NF90_NOWRITE, ncid)
         IF (ierr /= NF90_NOERR) THEN
            WRITE(*,*) 'ERROR: Cannot open UnitCatchment file: ', trim(parafile)
            CALL CoLM_stop()
         ENDIF

         ierr = nf90_inq_dimid(ncid, 'sed_n', dimid)
         IF (ierr /= NF90_NOERR) THEN
            WRITE(*,*) 'ERROR: Dimension sed_n not found in ', trim(parafile)
            CALL CoLM_stop()
         ENDIF
         ierr = nf90_inquire_dimension(ncid, dimid, len=nsed)

         ierr = nf90_inq_dimid(ncid, 'slope_layers', dimid)
         IF (ierr /= NF90_NOERR) THEN
            WRITE(*,*) 'ERROR: Dimension slope_layers not found in ', trim(parafile)
            CALL CoLM_stop()
         ENDIF
         ierr = nf90_inquire_dimension(ncid, dimid, len=nlfp_sed)

         ierr = nf90_close(ncid)
      ENDIF

#ifdef USEMPI
      CALL mpi_bcast(nsed, 1, MPI_INTEGER, p_address_master, p_comm_glb, p_err)
      CALL mpi_bcast(nlfp_sed, 1, MPI_INTEGER, p_address_master, p_comm_glb, p_err)
#endif

      IF (p_is_io) THEN
         WRITE(*,*) 'Sediment module: nsed=', nsed, ' nlfp_sed=', nlfp_sed
      ENDIF

      CALL read_sediment_parameter_file()
      CALL validate_sediment_parameters()

      CALL parse_grain_diameters()
      CALL calc_settling_velocities()
      CALL validate_sediment_parameters()

      ! Print the final effective parameters after all overrides/defaults
      ! and derived settling velocities have been resolved.
      CALL print_sediment_runtime_parameters()

      CALL read_sediment_static_data(parafile)
      CALL allocate_sediment_vars()
      CALL initialize_sediment_state()
#ifdef CoLMDEBUG
      IF (p_is_worker) CALL initialize_sediment_diag_stations()
#endif

      IF (p_is_io) THEN
         WRITE(*,*) 'Sediment module initialized successfully.'
      ENDIF

   END SUBROUTINE grid_sediment_init

   !-------------------------------------------------------------------------------------
   SUBROUTINE read_sediment_parameter_file()
   !-------------------------------------------------------------------------------------
   USE MOD_Tracer_Defs, only: tracer_param_file_for_index
   IMPLICIT NONE

   type(sediment_parameter_type) :: DEF_SEDIMENT
   character(len=512) :: file_param
   logical :: found, fexists
   integer :: ierr, unit_nml, ised
   character(len=512) :: iomsg
   namelist /nl_colm_sediment_parameter/ DEF_SEDIMENT

      IF (sediment_itrc <= 0) sediment_itrc = sediment_tracer_index()
      CALL tracer_param_file_for_index(sediment_itrc, 'SEDIMENT,SED', file_param, found)
      IF (.not. found) THEN
         IF (p_is_io) WRITE(*,'(A)') &
            'ERROR: missing DEF_TRACER_PARAM_FILES mapping for the active SEDIMENT tracer.'
         CALL CoLM_stop()
      ENDIF

      INQUIRE(file=trim(file_param), exist=fexists)
      IF (.not. fexists) THEN
         IF (p_is_io) WRITE(*,'(A,A)') 'ERROR: missing sediment parameter file: ', trim(file_param)
         CALL CoLM_stop()
      ENDIF

      open(newunit=unit_nml, status='OLD', file=trim(file_param), form='FORMATTED')
      iomsg = ''
      read(unit_nml, nml=nl_colm_sediment_parameter, iostat=ierr, iomsg=iomsg)
      close(unit_nml)
      IF (ierr /= 0) THEN
         IF (p_is_io) THEN
            WRITE(*,'(A,A)') &
            'ERROR: invalid &nl_colm_sediment_parameter in ', trim(file_param)
            WRITE(*,'(A)') TRIM(iomsg)
         ENDIF
         CALL CoLM_stop()
      ENDIF

      ! Ordered comparisons do not reject NaN: every comparison with NaN is
      ! false, so an invalid optional value could otherwise be mistaken for an
      ! omitted sentinel and silently fall back to a default.  Validate the
      ! parsed record itself before applying any optional override.
      IF (.not. ieee_is_finite(DEF_SEDIMENT%grain_density) .or. &
          .not. ieee_is_finite(DEF_SEDIMENT%water_density) .or. &
          .not. ieee_is_finite(DEF_SEDIMENT%porosity) .or. &
          .not. ieee_is_finite(DEF_SEDIMENT%ignore_depth_m) .or. &
          .not. ieee_is_finite(DEF_SEDIMENT%active_layer_depth) .or. &
          .not. ieee_is_finite(DEF_SEDIMENT%viscosity) .or. &
          .not. ieee_is_finite(DEF_SEDIMENT%von_karman) .or. &
          .not. ieee_is_finite(DEF_SEDIMENT%settling_multiplier) .or. &
          .not. ieee_is_finite(DEF_SEDIMENT%yield_coefficient) .or. &
          .not. ieee_is_finite(DEF_SEDIMENT%slope_exponent) .or. &
          .not. ieee_is_finite(DEF_SEDIMENT%precipitation_exponent) .or. &
          .not. ieee_is_finite(DEF_SEDIMENT%unit_conversion) .or. &
          .not. ieee_is_finite(DEF_SEDIMENT%cfl_adv) .or. &
          .not. ieee_is_finite(DEF_SEDIMENT%max_timestep_s) .or. &
          .not. ieee_is_finite(DEF_SEDIMENT%bed_depth) .or. &
          .not. ieee_is_finite(DEF_SEDIMENT%max_concentration) .or. &
          any(.not. ieee_is_finite(DEF_SEDIMENT%grain_diameter))) THEN
         IF (p_is_io) WRITE(*,'(A)') &
            'ERROR: sediment parameter file contains NaN or infinite values.'
         CALL CoLM_stop()
      ENDIF

      IF (DEF_SEDIMENT%nsed > 0 .and. DEF_SEDIMENT%nsed /= nsed) THEN
         IF (p_is_io) WRITE(*,'(A,I0,A,I0)') &
            'ERROR: sediment parameter nsed=', DEF_SEDIMENT%nsed, &
            ' does not match UnitCatchment sed_n=', nsed
         CALL CoLM_stop()
      ENDIF
      IF (DEF_SEDIMENT%grain_density >= 0._r8) THEN
         IF (DEF_SEDIMENT%grain_density < 100._r8) THEN
            IF (p_is_io) WRITE(*,'(A,ES12.4)') &
               'ERROR: sediment grain_density must be supplied in kg/m3, got ', &
               DEF_SEDIMENT%grain_density
            CALL CoLM_stop()
         ENDIF
         psedD = DEF_SEDIMENT%grain_density / 1000._r8
      ENDIF
      IF (DEF_SEDIMENT%water_density >= 0._r8) THEN
         IF (DEF_SEDIMENT%water_density < 100._r8) THEN
            IF (p_is_io) WRITE(*,'(A,ES12.4)') &
               'ERROR: sediment water_density must be supplied in kg/m3, got ', &
               DEF_SEDIMENT%water_density
            CALL CoLM_stop()
         ENDIF
         pwatD = DEF_SEDIMENT%water_density / 1000._r8
      ENDIF
      IF (DEF_SEDIMENT%porosity >= 0._r8) lambda = DEF_SEDIMENT%porosity
      IF (DEF_SEDIMENT%ndeposit_layers > 0) totlyrnum = DEF_SEDIMENT%ndeposit_layers
      IF (DEF_SEDIMENT%ignore_depth_m >= 0._r8) sed_ignore_dph = DEF_SEDIMENT%ignore_depth_m
      IF (DEF_SEDIMENT%active_layer_depth > 0._r8) lyrdph = DEF_SEDIMENT%active_layer_depth
      IF (DEF_SEDIMENT%viscosity > 0._r8) visKin = DEF_SEDIMENT%viscosity
      IF (DEF_SEDIMENT%von_karman > 0._r8) vonKar = DEF_SEDIMENT%von_karman
      IF (DEF_SEDIMENT%settling_multiplier > 0._r8) pset = DEF_SEDIMENT%settling_multiplier
      IF (DEF_SEDIMENT%yield_coefficient >= 0._r8) pyld = DEF_SEDIMENT%yield_coefficient
      IF (DEF_SEDIMENT%slope_exponent >= 0._r8) pyldc = DEF_SEDIMENT%slope_exponent
      IF (DEF_SEDIMENT%precipitation_exponent >= 0._r8) pyldpc = DEF_SEDIMENT%precipitation_exponent
      IF (DEF_SEDIMENT%unit_conversion >= 0._r8) dsylunit = DEF_SEDIMENT%unit_conversion
      IF (DEF_SEDIMENT%cfl_adv >= 0._r8) sed_cfl_adv = DEF_SEDIMENT%cfl_adv
      IF (DEF_SEDIMENT%max_timestep_s > 0._r8) sed_dt_max = DEF_SEDIMENT%max_timestep_s
      IF (DEF_SEDIMENT%bed_depth > 0._r8) sed_bed_depth = DEF_SEDIMENT%bed_depth
      IF (DEF_SEDIMENT%max_concentration >= 0._r8) sed_max_conc = DEF_SEDIMENT%max_concentration

      IF (DEF_SEDIMENT%grain_diameter(1) > 0._r8) THEN
         IF (nsed > MAX_SED_PARAM_CLASSES) THEN
            IF (p_is_io) WRITE(*,'(A,I0,A,I0)') &
               'ERROR: sediment parameter file supports at most ', MAX_SED_PARAM_CLASSES, &
               ' grain classes, got ', nsed
            CALL CoLM_stop()
         ENDIF
         IF (allocated(sDiam_from_param)) deallocate(sDiam_from_param)
         allocate(sDiam_from_param(nsed))
         DO ised = 1, nsed
            sDiam_from_param(ised) = DEF_SEDIMENT%grain_diameter(ised)
            IF (.not. ieee_is_finite(sDiam_from_param(ised)) .or. &
                sDiam_from_param(ised) <= 0._r8) THEN
               IF (p_is_io) WRITE(*,'(A,I0)') 'ERROR: missing/invalid sediment grain_diameter class ', ised
               CALL CoLM_stop()
            ENDIF
         ENDDO
      ENDIF

      IF (p_is_io) WRITE(*,'(A,A)') 'Sediment parameters loaded from ', trim(file_param)

   END SUBROUTINE read_sediment_parameter_file

   !-------------------------------------------------------------------------------------
   SUBROUTINE validate_sediment_parameters()
   !-------------------------------------------------------------------------------------
   IMPLICIT NONE

      IF (DEF_USE_BIFURCATION) THEN
         IF (p_is_io) WRITE(*,'(A)') &
            'ERROR: sediment bifurcation transport is not yet implemented; disable DEF_USE_BIFURCATION or remove the SEDIMENT particle tracer.'
         CALL CoLM_stop()
      ENDIF
      IF (DEF_USE_LEVEE) THEN
         IF (p_is_io) WRITE(*,'(A)') &
            'ERROR: sediment levee transport is not yet implemented; disable DEF_USE_LEVEE or remove the SEDIMENT particle tracer.'
         CALL CoLM_stop()
      ENDIF
      IF (nsed <= 0) THEN
         IF (p_is_io) WRITE(*,*) 'ERROR: sediment sed_n must be > 0, got ', nsed
         CALL CoLM_stop()
      ENDIF
      IF (nlfp_sed <= 0) THEN
         IF (p_is_io) WRITE(*,*) 'ERROR: sediment slope_layers must be > 0, got ', nlfp_sed
         CALL CoLM_stop()
      ENDIF
      IF (.not. ieee_is_finite(lambda) .or. &
          .not. ieee_is_finite(lyrdph) .or. &
          .not. ieee_is_finite(psedD) .or. &
          .not. ieee_is_finite(pwatD) .or. &
          .not. ieee_is_finite(visKin) .or. &
          .not. ieee_is_finite(vonKar) .or. &
          .not. ieee_is_finite(pset) .or. &
          .not. ieee_is_finite(pyld) .or. &
          .not. ieee_is_finite(pyldc) .or. &
          .not. ieee_is_finite(pyldpc) .or. &
          .not. ieee_is_finite(dsylunit) .or. &
          .not. ieee_is_finite(sed_ignore_dph) .or. &
          .not. ieee_is_finite(sed_cfl_adv) .or. &
          .not. ieee_is_finite(sed_dt_max) .or. &
          .not. ieee_is_finite(sed_bed_depth) .or. &
          .not. ieee_is_finite(sed_max_conc)) THEN

         IF (p_is_io) WRITE(*,'(A)') &
            'ERROR: sediment scalar configuration contains NaN or infinite values.'
         CALL CoLM_stop()
      ENDIF
      IF (lambda < 0._r8 .or. lambda >= 1._r8) THEN
         IF (p_is_io) WRITE(*,*) 'ERROR: sediment porosity must satisfy 0 <= lambda < 1, got ', lambda
         CALL CoLM_stop()
      ENDIF
      IF (lyrdph <= 0._r8) THEN
         IF (p_is_io) WRITE(*,*) 'ERROR: sediment active_layer_depth must be > 0, got ', lyrdph
         CALL CoLM_stop()
      ENDIF
      IF (psedD <= pwatD) THEN
         IF (p_is_io) WRITE(*,*) 'ERROR: psedD <= pwatD is non-physical for sediment settling/bedload:', psedD, pwatD
         CALL CoLM_stop()
      ENDIF
      IF (pwatD <= 0._r8 .or. visKin <= 0._r8 .or. vonKar <= 0._r8 .or. pset <= 0._r8) THEN
         IF (p_is_io) WRITE(*,*) 'ERROR: sediment density/viscosity/vonKarman/settling multiplier must be positive.'
         CALL CoLM_stop()
      ENDIF
      IF (totlyrnum <= 0) THEN
         IF (p_is_io) WRITE(*,*) 'ERROR: sediment ndeposit_layers must be > 0, got ', totlyrnum
         CALL CoLM_stop()
      ENDIF
      IF (sed_cfl_adv <= 0._r8 .or. sed_cfl_adv > 1._r8) THEN
         IF (p_is_io) WRITE(*,*) 'ERROR: sediment cfl_adv must satisfy 0 < cfl_adv <= 1, got ', sed_cfl_adv
         CALL CoLM_stop()
      ENDIF
      IF (sed_dt_max <= 0._r8) THEN
         IF (p_is_io) WRITE(*,*) 'ERROR: sediment max_timestep_s must be > 0, got ', sed_dt_max
         CALL CoLM_stop()
      ENDIF
      IF (sed_max_conc <= 0._r8 .or. sed_max_conc > 1._r8) THEN
           IF (p_is_io) WRITE(*,*) 'ERROR: sediment max_concentration must satisfy 0 < max_concentration <= 1, got ', sed_max_conc
           CALL CoLM_stop()
      ENDIF
      IF (sed_bed_depth <= 0._r8) THEN
         IF (p_is_io) WRITE(*,*) 'ERROR: sediment bed_depth must be > 0, got ', sed_bed_depth
         CALL CoLM_stop()
      ENDIF
      IF (sed_bed_depth < lyrdph * real(totlyrnum, r8)) THEN
         IF (p_is_io) WRITE(*,*) &
            'ERROR: sediment bed_depth must cover active plus named layers; got ', &
            sed_bed_depth, ' minimum ', lyrdph * real(totlyrnum, r8)
         CALL CoLM_stop()
      ENDIF
      IF (sed_ignore_dph < 0._r8) THEN
         IF (p_is_io) WRITE(*,*) 'ERROR: sediment ignore_depth_m must be >= 0, got ', sed_ignore_dph
         CALL CoLM_stop()
      ENDIF
      IF (pyld < 0._r8 .or. pyldc < 0._r8 .or. pyldpc < 0._r8 .or. dsylunit < 0._r8) THEN
         IF (p_is_io) WRITE(*,*) 'ERROR: sediment yield parameters must be non-negative.'
         CALL CoLM_stop()
      ENDIF
      IF (allocated(sDiam)) THEN
         IF (any(.not. ieee_is_finite(sDiam)) .or. any(sDiam <= 0._r8)) THEN
            IF (p_is_io) WRITE(*,'(A)') &
               'ERROR: sediment grain diameters must be finite and positive.'
            CALL CoLM_stop()
         ENDIF
            IF (lyrdph < maxval(sDiam)) THEN
               IF (p_is_io) WRITE(*,'(A,ES12.4,A,ES12.4)') &
                  'WARNING: sediment active_layer_depth is smaller than max grain diameter; keep user value. lyrdph=', &
                  lyrdph, ', max_sDiam=', maxval(sDiam)
               ! 不再执行：lyrdph = maxval(sDiam)
            ENDIF
      ENDIF
      IF (allocated(setvel)) THEN
         IF (any(.not. ieee_is_finite(setvel)) .or. any(setvel < 0._r8)) THEN
            IF (p_is_io) WRITE(*,'(A)') &
               'ERROR: sediment settling velocities must be finite and non-negative.'
            CALL CoLM_stop()
         ENDIF
      ENDIF
      IF (p_is_worker .and. allocated(sed_frc)) THEN
         IF (any(.not. ieee_is_finite(sed_frc)) .or. any(sed_frc < 0._r8)) THEN
            IF (p_is_io) WRITE(*,'(A)') &
               'ERROR: sediment fractions must be finite and non-negative.'
            CALL CoLM_stop()
         ENDIF
      ENDIF
      IF (p_is_worker .and. allocated(sed_slope)) THEN
         IF (any(.not. ieee_is_finite(sed_slope)) .or. any(sed_slope < 0._r8) .or. &
             any(abs(sed_slope) > 0.5_r8 * abs(spval))) THEN
            IF (p_is_io) WRITE(*,'(A)') &
               'ERROR: sediment slopes must be finite, non-negative, and not missing values.'
            CALL CoLM_stop()
         ENDIF
      ENDIF

   END SUBROUTINE validate_sediment_parameters

   !-------------------------------------------------------------------------------------
   SUBROUTINE grid_sediment_calc(deltime)
   ! Main sediment calculation. Called from MOD_Grid_RiverLakeFlow after water routing.
   !-------------------------------------------------------------------------------------
   USE MOD_Grid_RiverLakeNetwork, only: numucat, topo_rivwth, topo_rivlen, &
      topo_rivman, topo_area, ucat_next, ucat_ucid, x_ucat, y_ucat, griducat
   USE MOD_Const_Physical, only: grav
   IMPLICIT NONE

   real(r8), intent(in) :: deltime

   real(r8) :: sed_time_remaining, dt_morph, dt_adv, dt_adv_remaining
   real(r8) :: avg_v2, avg_wdsrf, avg_rivsto, avg_rivout, avg_abs_rivout
   real(r8) :: sed_flow_cancel_ratio, sed_carrier_wet_fraction
   real(r8), allocatable :: rivsto(:), rivout(:), rivout_abs(:), bed_area(:), fldfrc(:)
   logical,  allocatable :: wet_seen(:), shallow_seen(:), source_seen(:)
   logical,  allocatable :: susp_seen(:), bed_seen(:), exch_pos_seen(:), exch_neg_seen(:)
   logical,  allocatable :: es_raw_seen(:), d_raw_seen(:), es_eff_seen(:), d_eff_seen(:)
   real(r8) :: precip_time_local
   integer  :: i, ised, iter_sed, iter_adv
   integer  :: clk_total_start, clk_total_end
   integer  :: clk_phase_start, clk_phase_end, clk_rate
   real(r8) :: t_total, t_yield, t_adv, t_input, t_exchange, t_layer, t_diag
   real(r8) :: max_sed_precip_local, max_precip_rate_local, max_slope_local
   real(r8) :: max_sedcon_local, max_sedout_local, max_bedout_local
   real(r8) :: max_sedinp_local, max_netflw_local, max_shearvel_local
   real(r8) :: sum_layer_local, sum_seddep_local, sum_sedsto_local
   real(r8) :: sum_sedinp_local, sum_sedout_down_local, sum_sedout_up_local
   real(r8) :: sum_sedout_abs_local, sum_netflw_pos_local, sum_netflw_neg_local
   real(r8) :: sum_bed_bulk_local, sum_bed_solid_local
   real(r8) :: sum_rivout_signed_local, sum_rivout_abs_local
   real(r8) :: max_flow_cancel_local
   real(r8) :: sum_es_raw_local, sum_d_raw_local
   real(r8) :: max_es_raw_local, max_d_raw_local
   real(r8) :: sum_es_eff_local, sum_d_eff_local
   real(r8) :: max_es_eff_local, max_d_eff_local
   real(r8) :: sum_inst_near_dry_abs_q_local
   real(r8) :: sum_period_near_dry_abs_q_local, sum_period_near_dry_signed_q_local
   real(r8) :: precip_diag_global(3), diag_max_global(11), diag_sum_global(17)
   real(r8) :: carrier_filter_diag_global(3)
   real(r8) :: dt_cfl_local, dt_cfl_global, dt_cell
   real(r8) :: sedout_val, bedout_val, netflw_val, d_eff_val
   integer  :: n_wet_local, n_shallow_local, n_source_local
   integer  :: n_susp_local, n_bed_local
   integer  :: n_exchange_pos_local, n_exchange_neg_local
   integer  :: n_es_raw_local, n_d_raw_local, n_es_eff_local, n_d_eff_local
   integer  :: n_flow_cancel_local, n_period_near_dry_local
   integer  :: diag_count_global(12), carrier_filter_count_global(1)
   logical  :: invalid_sediment_found
   integer  :: extreme_meta_local(2,4), extreme_meta_global(2,4)
   integer  :: extreme_cell_global(2), iextreme
   real(r8) :: extreme_state_local(2,10), extreme_state_global(2,10)
   real(r8), allocatable :: extreme_sedcon_local(:,:), extreme_sedcon_global(:,:)
   real(r8), allocatable :: extreme_sedout_local(:,:), extreme_sedout_global(:,:)
   real(r8), allocatable :: extreme_sedinp_local(:,:), extreme_sedinp_global(:,:)
   real(r8), allocatable :: extreme_netflw_local(:,:), extreme_netflw_global(:,:)
   real(r8) :: extreme_lon, extreme_lat, extreme_trigger_value
   real(r8) :: sedcon_total_diag, sedout_total_diag, sedout_abs_total_diag
   real(r8) :: ssc_mg_l_diag, ssl_t_day_diag, ssl_abs_t_day_diag
#ifdef CoLMDEBUG
   integer  :: ista, station_i
   real(r8) :: station_state_local(SED_N_DIAG_STATIONS,4)
   real(r8) :: station_state_global(SED_N_DIAG_STATIONS,4)
   real(r8), allocatable :: station_sedcon_local(:,:), station_sedcon_global(:,:)
   real(r8), allocatable :: station_sedout_local(:,:), station_sedout_global(:,:)
   real(r8) :: station_sedcon_total, station_sedout_total, station_sedout_abs_total
   real(r8) :: station_ssc_mg_l, station_ssl_t_day, station_ssl_abs_t_day
#endif
   real(r8) :: timing_local(8)
   real(r8) :: timing_min(8)
   real(r8) :: timing_max(8)
   real(r8) :: timing_sum(8)
   real(r8) :: timing_mean(8)
   real(r8), parameter :: CFL_RIVOUT_EPS = 1.e-12_r8


      IF (.not. sediment_particle_enabled()) RETURN
      IF (.not. p_is_worker) RETURN

#ifdef CoLMDEBUG
      CALL system_clock(clk_total_start, clk_rate)
#endif

      allocate(rivsto(numucat))
      allocate(rivout(numucat))
      allocate(rivout_abs(numucat))
      allocate(bed_area(numucat))
      allocate(fldfrc(numucat))
      allocate(wet_seen(numucat), shallow_seen(numucat), source_seen(numucat))
      allocate(susp_seen(numucat), bed_seen(numucat), exch_pos_seen(numucat), exch_neg_seen(numucat))
      allocate(es_raw_seen(numucat), d_raw_seen(numucat), es_eff_seen(numucat), d_eff_seen(numucat))
#ifdef CoLMDEBUG
      allocate(extreme_sedcon_local(nsed,2), extreme_sedcon_global(nsed,2))
      allocate(extreme_sedout_local(nsed,2), extreme_sedout_global(nsed,2))
      allocate(extreme_sedinp_local(nsed,2), extreme_sedinp_global(nsed,2))
      allocate(extreme_netflw_local(nsed,2), extreme_netflw_global(nsed,2))
      allocate(station_sedcon_local(nsed,SED_N_DIAG_STATIONS))
      allocate(station_sedcon_global(nsed,SED_N_DIAG_STATIONS))
      allocate(station_sedout_local(nsed,SED_N_DIAG_STATIONS))
      allocate(station_sedout_global(nsed,SED_N_DIAG_STATIONS))
#endif

      ! Compute flooded fraction from current routing period accumulators only,
      ! not from history-period averages.  This preserves the flood-exposure
      ! time-series at routing-period resolution for hillslope erosion.
      DO i = 1, numucat
         IF (sed_acc_time(i) > 0._r8 .and. topo_area(i) > 0._r8) THEN
            fldfrc(i) = sed_acc_floodarea(i) / sed_acc_time(i) / topo_area(i)
         ELSE
            fldfrc(i) = 0._r8
         ENDIF
         fldfrc(i) = min(max(fldfrc(i), 0._r8), 1._r8)
      ENDDO

      iter_sed = 0
      iter_adv = 0
      t_yield = 0._r8
      t_adv = 0._r8
      t_input = 0._r8
      t_exchange = 0._r8
      t_layer = 0._r8
      t_diag = 0._r8
      sum_sedinp_local = 0._r8
      sum_sedout_down_local = 0._r8
      sum_sedout_up_local = 0._r8
      sum_sedout_abs_local = 0._r8
      sum_netflw_pos_local = 0._r8
      sum_netflw_neg_local = 0._r8
      sum_rivout_signed_local = 0._r8
      sum_rivout_abs_local = 0._r8
      max_flow_cancel_local = 0._r8
      sum_inst_near_dry_abs_q_local = 0._r8
      sum_period_near_dry_abs_q_local = 0._r8
      sum_period_near_dry_signed_q_local = 0._r8
      sum_es_raw_local = 0._r8
      sum_d_raw_local = 0._r8
      sum_es_eff_local = 0._r8
      sum_d_eff_local = 0._r8
      n_flow_cancel_local = 0
      n_period_near_dry_local = 0
      wet_seen = .false.
      shallow_seen = .false.
      source_seen = .false.
      susp_seen = .false.
      bed_seen = .false.
      exch_pos_seen = .false.
      exch_neg_seen = .false.
      es_raw_seen = .false.
      d_raw_seen = .false.
      es_eff_seen = .false.
      d_eff_seen = .false.

#ifdef CoLMDEBUG
      ! Start the wall-clock timer before any sediment-period diagnostics or
      ! sediment operators.  The corresponding stop and MPI aggregation are
      ! performed only after the full sediment calculation has completed.
      CALL system_clock(clk_total_start, clk_rate)
#endif

      ! Store precipitation averaging time before reset
      precip_time_local = sed_precip_time

#ifdef CoLMDEBUG
      max_sed_precip_local = 0._r8
      max_precip_rate_local = 0._r8
      max_slope_local = 0._r8
      max_sedcon_local = 0._r8
      max_sedout_local = 0._r8
      max_bedout_local = 0._r8
      max_sedinp_local = 0._r8
      max_netflw_local = 0._r8
      max_shearvel_local = 0._r8
      max_es_raw_local = 0._r8
      max_d_raw_local = 0._r8
      max_es_eff_local = 0._r8
      max_d_eff_local = 0._r8
      extreme_meta_local = 0
      extreme_meta_global = 0
      extreme_cell_global = huge(1)
      extreme_state_local = 0._r8
      extreme_state_global = 0._r8
      extreme_sedcon_local = 0._r8
      extreme_sedcon_global = 0._r8
      extreme_sedout_local = 0._r8
      extreme_sedout_global = 0._r8
      extreme_sedinp_local = 0._r8
      extreme_sedinp_global = 0._r8
      extreme_netflw_local = 0._r8
      extreme_netflw_global = 0._r8
      IF (numucat > 0) THEN
         max_sed_precip_local = maxval(sed_precip)
         max_precip_rate_local = max_sed_precip_local / max(precip_time_local, 1.e-20_r8)
         max_slope_local = maxval(sed_slope)
      ENDIF
      precip_diag_global = (/ max_sed_precip_local, max_precip_rate_local, max_slope_local /)
#ifdef USEMPI
      CALL mpi_allreduce(MPI_IN_PLACE, precip_diag_global, size(precip_diag_global), &
         MPI_REAL8, MPI_MAX, p_comm_worker, p_err)
#endif

      ! Diagnostic: check precipitation forcing reaching sediment module
      IF (p_iam_worker == 0) THEN
         WRITE(*,'(A,ES10.3,A,ES10.3,A,ES10.3,A,ES10.3)') &
            'Sediment precip diag: prcp_time=', precip_time_local, &
            ', max_sed_precip=', precip_diag_global(1), &
            ', max_precip_rate[mm/s]=', precip_diag_global(2), &
            ', max_slope=', precip_diag_global(3)
      ENDIF
#endif

      sed_time_remaining = deltime

      DO WHILE (sed_time_remaining > 0._r8)
         iter_sed = iter_sed + 1
         dt_morph = min(sed_time_remaining, sed_dt_max)

         ! Calculate average water flow variables from per-cell accumulators
         DO i = 1, numucat
            IF (sed_acc_time(i) > 0._r8) THEN
               avg_v2     = sed_acc_v2(i)     / sed_acc_time(i)
               avg_wdsrf  = sed_acc_wdsrf(i)  / sed_acc_time(i)
               avg_rivsto = sed_acc_rivsto(i) / sed_acc_time(i)
               avg_rivout = sed_acc_rivout(i) / sed_acc_time(i)
               avg_abs_rivout = sed_acc_abs_rivout(i) / sed_acc_time(i)
               IF (avg_abs_rivout > CFL_RIVOUT_EPS) THEN
                  sed_flow_cancel_ratio = 1._r8 - min(abs(avg_rivout) / avg_abs_rivout, 1._r8)
               ELSE
                  sed_flow_cancel_ratio = 0._r8
               ENDIF
               sed_carrier_wet_fraction = sed_acc_carrier_time(i) / sed_acc_time(i)
               IF (iter_sed == 1) THEN
                  sum_inst_near_dry_abs_q_local = sum_inst_near_dry_abs_q_local &
                     + sed_acc_near_dry_abs_rivout(i) / sed_acc_time(i)
               ENDIF

               ! -------------------------------------------------------------
               ! TEMP DEBUG: diagnose cells producing very small sediment CFL
               ! Print cells that would require more than 10000 advection steps.
               ! Remove after diagnosing the CFL problem.
               ! -------------------------------------------------------------
               IF (avg_rivsto > 0._r8 .and. avg_abs_rivout > CFL_RIVOUT_EPS) THEN

                  dt_cell = sed_cfl_adv * avg_rivsto / avg_abs_rivout

                  IF (dt_cell < dt_morph / 10000._r8) THEN

                     WRITE(*,'(A)') '========== SED_CFL_DEBUG =========='
                     WRITE(*,'(A,I0)')      'worker             = ', p_iam_worker
                     WRITE(*,'(A,I0)')      'cell i             = ', i
                     WRITE(*,'(A,I0)')      'ucat_next          = ', ucat_next(i)

                     WRITE(*,'(A,ES20.10)') 'avg_rivsto [m3]    = ', avg_rivsto
                     WRITE(*,'(A,ES20.10)') 'avg_rivout [m3/s]  = ', avg_rivout
                     WRITE(*,'(A,ES20.10)') 'avg_abs_rivout     = ', avg_abs_rivout

                     WRITE(*,'(A,ES20.10)') 'avg_wdsrf [m]      = ', avg_wdsrf
                     WRITE(*,'(A,ES20.10)') 'avg_v2             = ', avg_v2
                     WRITE(*,'(A,ES20.10)') 'sed_acc_time [s]   = ', sed_acc_time(i)
                     WRITE(*,'(A,ES20.10)') 'carrier_wet_frac   = ', sed_carrier_wet_fraction
                     WRITE(*,'(A,ES20.10)') 'min_wdsrf [m]      = ', sed_acc_wdsrf_min(i)
                     WRITE(*,'(A,ES20.10)') 'max_wdsrf [m]      = ', sed_acc_wdsrf_max(i)
                     WRITE(*,'(A,ES20.10)') 'min_rivsto [m3]    = ', sed_acc_rivsto_min(i)
                     WRITE(*,'(A,ES20.10)') 'max_rivsto [m3]    = ', sed_acc_rivsto_max(i)
                     WRITE(*,'(A,ES20.10)') 'min_rivout [m3/s]  = ', sed_acc_rivout_min(i)
                     WRITE(*,'(A,ES20.10)') 'max_rivout [m3/s]  = ', sed_acc_rivout_max(i)
                     WRITE(*,'(A,ES20.10)') 'int_pos_Q [m3]     = ', sed_acc_pos_rivout(i)
                     WRITE(*,'(A,ES20.10)') 'int_neg_Q_abs [m3] = ', sed_acc_neg_rivout(i)
                     WRITE(*,'(A,ES20.10)') 'inst_filtered_abs_Q[m3] = ', &
                        sed_acc_near_dry_abs_rivout(i)

                     WRITE(*,'(A,ES20.10)') 'dt_cell [s]        = ', dt_cell
                     WRITE(*,'(A,ES20.10)') 'required_substeps  = ', &
                        dt_morph / dt_cell

                     WRITE(*,'(A,ES20.10)') 'flow_cancel_ratio  = ', sed_flow_cancel_ratio

                     IF (avg_rivout > 0._r8) THEN
                        WRITE(*,'(A)') 'flow_direction     = FORWARD'
                     ELSEIF (avg_rivout < 0._r8) THEN
                        WRITE(*,'(A)') 'flow_direction     = REVERSE'
                     ELSE
                        WRITE(*,'(A)') 'flow_direction     = ZERO-MEAN/OSCILLATORY'
                     ENDIF

                     WRITE(*,'(A)') '==================================='

                  ENDIF
               ENDIF

               ! -------------------------------------------------------------
               ! Defensive cleanup of an inconsistent near-dry carrier state.
               !
               ! Keep the DEBUG block above this point so that the raw state
               ! that would otherwise control the CFL remains visible.
               !
               ! This catches cases such as the Amazon failure:
               !   avg_wdsrf ~ 0
               !   avg_rivout > 0
               ! where the routing-period mean state is below the sediment
               ! carrier threshold but a short wet pulse still contributed
               ! forward-dominant carrier flux. Filter only the carrier flux:
               ! avg_rivsto and avg_v2 still describe the period-mean water
               ! state used by concentration, shear, and exchange diagnostics.
               ! -------------------------------------------------------------
               IF (avg_wdsrf <= SED_NEAR_DRY_DEPTH .and. &
                   avg_abs_rivout > CFL_RIVOUT_EPS .and. avg_rivout > 0._r8 .and. &
                   sed_flow_cancel_ratio <= 0.5_r8) THEN

                  IF (iter_sed == 1) THEN
                     sum_period_near_dry_abs_q_local = sum_period_near_dry_abs_q_local &
                        + avg_abs_rivout
                     sum_period_near_dry_signed_q_local = sum_period_near_dry_signed_q_local &
                        + avg_rivout
                     n_period_near_dry_local = n_period_near_dry_local + 1
                  ENDIF

                  avg_rivout = 0._r8
                  avg_abs_rivout = 0._r8

               ENDIF


               IF (iter_sed == 1) THEN
                  sum_rivout_signed_local = sum_rivout_signed_local + avg_rivout
                  sum_rivout_abs_local = sum_rivout_abs_local + avg_abs_rivout
                  IF (avg_abs_rivout > CFL_RIVOUT_EPS) THEN
                     max_flow_cancel_local = max(max_flow_cancel_local, sed_flow_cancel_ratio)
                     IF (sed_flow_cancel_ratio > 0.5_r8) n_flow_cancel_local = n_flow_cancel_local + 1
                  ENDIF
               ENDIF

               ! Shear velocity from RMS velocity: u* = sqrt(g * n^2 * <v^2> * d^(-1/3))
               ! Using <v^2> (mean of squared velocity) avoids sign cancellation
               ! when flow direction oscillates (tidal/backwater areas).
               ! HYDRO sets velocity to zero only while reservoir water is
               ! stationary. Keying on live flow preserves river shear before
               ! a scheduled reservoir is actually built.
               IF (avg_wdsrf > 0._r8 .and. avg_v2 > 0._r8) THEN
                  shearvel(i) = sqrt(grav * topo_rivman(i)**2 * avg_v2 &
                     * avg_wdsrf**(-1._r8/3._r8))
               ELSE
                  shearvel(i) = 0._r8
               ENDIF
               CALL calc_critical_shear_egiazoroff(i, shearvel(i), critshearvel(:,i))
               CALL calc_suspend_velocity(critshearvel(:,i), shearvel(i), susvel(:,i))

               ! Use the HYDRO-owned water storage so reservoirs/lakes keep
               ! sediment concentration and CFL consistent with water routing.
               rivsto(i) = max(avg_rivsto, 0._r8)
               rivout(i) = avg_rivout
               rivout_abs(i) = avg_abs_rivout
               bed_area(i) = topo_rivwth(i) * topo_rivlen(i)
               IF (avg_v2 <= 0._r8) THEN
                  bed_area(i) = max(bed_area(i), sed_acc_floodarea(i) / sed_acc_time(i))
               ENDIF
            ELSE
               shearvel(i) = 0._r8
               critshearvel(:,i) = 1.e20_r8
               susvel(:,i) = 0._r8
               rivsto(i) = 0._r8
               rivout(i) = 0._r8
               rivout_abs(i) = 0._r8
               bed_area(i) = topo_rivwth(i) * topo_rivlen(i)
            ENDIF
         ENDDO

         ! Derive concentration exactly once at the routing-period boundary.
         ! Every operator below mutates the shared canonical solid volume;
         ! none may reconstruct mass from concentration and HYDRO storage.
         IF (iter_sed == 1) CALL begin_suspended_period(rivsto)

         dt_cfl_local = dt_morph
         DO i = 1, numucat
            IF (rivsto(i) <= 0._r8) CYCLE
            ! Advection transports every wet cell, including cells below the
            ! morphology/exchange depth threshold.  Therefore every wet carrier
            ! must constrain the same CFL step; excluding shallow cells here
            ! would silently bypass sed_cfl_adv for exactly those cells.
            IF (rivout_abs(i) <= CFL_RIVOUT_EPS) CYCLE
            dt_cell = sed_cfl_adv * rivsto(i) / rivout_abs(i)
            dt_cfl_local = min(dt_cfl_local, dt_cell)
         ENDDO
#ifdef USEMPI
         dt_cfl_global = dt_cfl_local
         CALL mpi_allreduce(MPI_IN_PLACE, dt_cfl_global, 1, MPI_REAL8, MPI_MIN, p_comm_worker, p_err)
#else
         dt_cfl_global = dt_cfl_local
#endif
         ! This form also rejects NaN because every ordered comparison with
         ! NaN is false. Never enter the advection loop with a zero step.
         IF (.not. (dt_cfl_global > 0._r8)) THEN
            IF (p_iam_worker == 0) THEN
               WRITE(*,'(A,ES12.4,A,ES12.4)') &
                  'ERROR sediment: invalid CFL timestep=', dt_cfl_global, &
                  ' s; morphology interval=', dt_morph
            ENDIF
            CALL CoLM_stop('sediment advection CFL timestep must be valid and positive')
         ENDIF
         IF (dt_cfl_global < dt_morph / real(SED_MAX_ADV_SUBSTEPS, r8)) THEN
            IF (p_iam_worker == 0) THEN
               WRITE(*,'(A,ES12.4,A,ES12.4,A,ES12.4,A,I0)') &
                  'ERROR sediment: pathological CFL timestep=', dt_cfl_global, &
                  ' s; morphology interval=', dt_morph, &
                  ' s; required advection substeps=', dt_morph / dt_cfl_global, &
                  ' exceeds limit=', SED_MAX_ADV_SUBSTEPS
            ENDIF
            CALL CoLM_stop('sediment advection CFL requires too many substeps')
         ENDIF

         CALL system_clock(clk_phase_start)
         CALL calc_sediment_yield(fldfrc, topo_area, precip_time_local)
         CALL system_clock(clk_phase_end)
         IF (clk_rate > 0) t_yield = t_yield + real(clk_phase_end - clk_phase_start, r8) / real(clk_rate, r8)

         CALL system_clock(clk_phase_start)
         CALL calc_sediment_exchange(dt_morph, rivsto, bed_area)
         CALL system_clock(clk_phase_end)
         IF (clk_rate > 0) t_exchange = t_exchange + real(clk_phase_end - clk_phase_start, r8) / real(clk_rate, r8)

         CALL system_clock(clk_phase_start)
         CALL apply_sediment_input(dt_morph, rivsto, bed_area)
         CALL system_clock(clk_phase_end)
         IF (clk_rate > 0) t_input = t_input + real(clk_phase_end - clk_phase_start, r8) / real(clk_rate, r8)

         dt_adv_remaining = dt_morph
         DO WHILE (dt_adv_remaining > 0._r8)
            iter_adv = iter_adv + 1
            dt_adv = min(dt_adv_remaining, dt_cfl_global)

            CALL system_clock(clk_phase_start)
            CALL calc_sediment_advection(dt_adv, rivout, rivout_abs, rivsto)
            CALL system_clock(clk_phase_end)
            IF (clk_rate > 0) t_adv = t_adv + real(clk_phase_end - clk_phase_start, r8) / real(clk_rate, r8)

            CALL debug_check_sediment_fields('after advection', iter_sed, iter_adv, &
               dt_morph, dt_adv, dt_cfl_global, rivout, rivout_abs, rivsto, invalid_sediment_found)
            IF (invalid_sediment_found) THEN
               CALL CoLM_stop('invalid finite value in sediment fields after advection')
            ENDIF

            CALL system_clock(clk_phase_start)
            CALL accumulate_sediment_output(dt_adv)
            CALL system_clock(clk_phase_end)
            IF (clk_rate > 0) t_diag = t_diag + real(clk_phase_end - clk_phase_start, r8) / real(clk_rate, r8)

            CALL debug_check_sediment_fields('after accumulate_sediment_output', iter_sed, iter_adv, &
               dt_morph, dt_adv, dt_cfl_global, rivout, rivout_abs, rivsto, invalid_sediment_found)
            IF (invalid_sediment_found) THEN
               CALL CoLM_stop('invalid finite value in sediment fields after accumulation')
            ENDIF

            IF (numucat > 0) THEN
               DO i = 1, numucat
                  wet_seen(i) = wet_seen(i) .or. (rivsto(i) > 0._r8)
                  shallow_seen(i) = shallow_seen(i) .or. &
                     (rivsto(i) > 0._r8 .and. rivsto(i) < topo_rivwth(i) * topo_rivlen(i) * sed_ignore_dph)
                  max_shearvel_local = max(max_shearvel_local, shearvel(i))

                  DO ised = 1, nsed
                     sedout_val = sedout(ised,i)
                     bedout_val = bedout(ised,i)
                     netflw_val = netflw(ised,i) + netflw_adv_step(ised,i)
                     d_eff_val = exch_d_eff(ised,i) + exch_d_adv_step(ised,i)

#ifdef CoLMDEBUG
                     IF (sedcon(ised,i) > max_sedcon_local .or. &
                         (sedcon(ised,i) == max_sedcon_local .and. sedcon(ised,i) > 0._r8 .and. &
                          (extreme_meta_local(1,1) == 0 .or. ucat_ucid(i) < extreme_meta_local(1,1)))) THEN
                        max_sedcon_local = sedcon(ised,i)
                        extreme_meta_local(1,:) = (/ ucat_ucid(i), x_ucat(i), y_ucat(i), ised /)
                        IF (sed_acc_time(i) > 0._r8) THEN
                           extreme_state_local(1,1) = sed_acc_wdsrf(i) / sed_acc_time(i)
                           extreme_state_local(1,2) = sed_acc_wdsrf_min(i)
                           extreme_state_local(1,3) = sed_acc_wdsrf_max(i)
                        ELSE
                           extreme_state_local(1,1:3) = 0._r8
                        ENDIF
                        extreme_state_local(1,4) = rivsto(i)
                        extreme_state_local(1,5) = rivout(i)
                        extreme_state_local(1,6) = rivout_abs(i)
                        extreme_state_local(1,7) = sed_acc_rivout_min(i)
                        extreme_state_local(1,8) = sed_acc_rivout_max(i)
                        extreme_state_local(1,9) = shearvel(i)
                        extreme_state_local(1,10) = bed_area(i)
                        extreme_sedcon_local(:,1) = sedcon(:,i)
                        extreme_sedout_local(:,1) = sedout(:,i)
                        extreme_sedinp_local(:,1) = sedinp(:,i)
                        extreme_netflw_local(:,1) = netflw(:,i) + netflw_adv_step(:,i)
                     ENDIF
                     IF (abs(sedout_val) > max_sedout_local .or. &
                         (abs(sedout_val) == max_sedout_local .and. abs(sedout_val) > 0._r8 .and. &
                          (extreme_meta_local(2,1) == 0 .or. ucat_ucid(i) < extreme_meta_local(2,1)))) THEN
                        max_sedout_local = abs(sedout_val)
                        extreme_meta_local(2,:) = (/ ucat_ucid(i), x_ucat(i), y_ucat(i), ised /)
                        IF (sed_acc_time(i) > 0._r8) THEN
                           extreme_state_local(2,1) = sed_acc_wdsrf(i) / sed_acc_time(i)
                           extreme_state_local(2,2) = sed_acc_wdsrf_min(i)
                           extreme_state_local(2,3) = sed_acc_wdsrf_max(i)
                        ELSE
                           extreme_state_local(2,1:3) = 0._r8
                        ENDIF
                        extreme_state_local(2,4) = rivsto(i)
                        extreme_state_local(2,5) = rivout(i)
                        extreme_state_local(2,6) = rivout_abs(i)
                        extreme_state_local(2,7) = sed_acc_rivout_min(i)
                        extreme_state_local(2,8) = sed_acc_rivout_max(i)
                        extreme_state_local(2,9) = shearvel(i)
                        extreme_state_local(2,10) = bed_area(i)
                        extreme_sedcon_local(:,2) = sedcon(:,i)
                        extreme_sedout_local(:,2) = sedout(:,i)
                        extreme_sedinp_local(:,2) = sedinp(:,i)
                        extreme_netflw_local(:,2) = netflw(:,i) + netflw_adv_step(:,i)
                     ENDIF
#else
                     max_sedcon_local = max(max_sedcon_local, sedcon(ised,i))
                     max_sedout_local = max(max_sedout_local, abs(sedout_val))
#endif
                     max_bedout_local = max(max_bedout_local, abs(bedout_val))
                     max_sedinp_local = max(max_sedinp_local, sedinp(ised,i))
                     max_netflw_local = max(max_netflw_local, abs(netflw_val))
                     max_es_raw_local = max(max_es_raw_local, exch_es_raw(ised,i))
                     max_d_raw_local = max(max_d_raw_local, exch_d_raw(ised,i))
                     max_es_eff_local = max(max_es_eff_local, exch_es_eff(ised,i))
                     max_d_eff_local = max(max_d_eff_local, d_eff_val)

                     sum_sedinp_local = sum_sedinp_local + sedinp(ised,i) * dt_adv
                     sum_sedout_abs_local = sum_sedout_abs_local + abs(sedout_val) * dt_adv
                     sum_es_raw_local = sum_es_raw_local + exch_es_raw(ised,i) * dt_adv
                     sum_d_raw_local = sum_d_raw_local + exch_d_raw(ised,i) * dt_adv
                     sum_es_eff_local = sum_es_eff_local + exch_es_eff(ised,i) * dt_adv
                     sum_d_eff_local = sum_d_eff_local + d_eff_val * dt_adv

                     IF (sedout_val > 0._r8) THEN
                        sum_sedout_down_local = sum_sedout_down_local + sedout_val * dt_adv
                        susp_seen(i) = .true.
                     ELSEIF (sedout_val < 0._r8) THEN
                        sum_sedout_up_local = sum_sedout_up_local - sedout_val * dt_adv
                        susp_seen(i) = .true.
                     ENDIF
                     IF (bedout_val /= 0._r8) bed_seen(i) = .true.
                     IF (sedinp(ised,i) > 0._r8) source_seen(i) = .true.
                     IF (netflw_val > 0._r8) THEN
                        sum_netflw_pos_local = sum_netflw_pos_local + netflw_val * dt_adv
                        exch_pos_seen(i) = .true.
                     ELSEIF (netflw_val < 0._r8) THEN
                        sum_netflw_neg_local = sum_netflw_neg_local - netflw_val * dt_adv
                        exch_neg_seen(i) = .true.
                     ENDIF
                     IF (exch_es_raw(ised,i) > 0._r8) es_raw_seen(i) = .true.
                     IF (exch_d_raw(ised,i) > 0._r8) d_raw_seen(i) = .true.
                     IF (exch_es_eff(ised,i) > 0._r8) es_eff_seen(i) = .true.
                     IF (d_eff_val > 0._r8) d_eff_seen(i) = .true.
                  ENDDO
               ENDDO
            ENDIF

            dt_adv_remaining = dt_adv_remaining - dt_adv
         ENDDO

         CALL system_clock(clk_phase_start)
         CALL calc_layer_redistribution(bed_area)
         CALL system_clock(clk_phase_end)
         IF (clk_rate > 0) t_layer = t_layer + real(clk_phase_end - clk_phase_start, r8) / real(clk_rate, r8)

         sed_time_remaining = sed_time_remaining - dt_morph
      ENDDO

      ! Publish the diagnostic concentration only after every operator has
      ! completed; canonical suspended solid volume is already updated.
      CALL commit_suspended_period(rivsto)

      ! Accumulate total time for history output averaging
      sed_hist_acctime = sed_hist_acctime + deltime

      ! Reset accumulation variables
      sed_acc_time(:)      = 0._r8
      sed_acc_v2(:)        = 0._r8
      sed_acc_wdsrf(:)     = 0._r8
      sed_acc_rivsto(:)    = 0._r8
      sed_acc_rivout(:)    = 0._r8
      sed_acc_abs_rivout(:)= 0._r8
      sed_acc_floodarea(:) = 0._r8
      sed_acc_carrier_time(:) = 0._r8
      sed_acc_wdsrf_min(:) = huge(1._r8)
      sed_acc_wdsrf_max(:) = 0._r8
      sed_acc_rivsto_min(:) = huge(1._r8)
      sed_acc_rivsto_max(:) = 0._r8
      sed_acc_rivout_min(:) = huge(1._r8)
      sed_acc_rivout_max(:) = -huge(1._r8)
      sed_acc_pos_rivout(:) = 0._r8
      sed_acc_neg_rivout(:) = 0._r8
      sed_acc_near_dry_abs_rivout(:) = 0._r8
      sed_precip(:)        = 0._r8
      sed_precip_yield(:)  = 0._r8
      sed_precip_time      = 0._r8

#ifdef CoLMDEBUG
      ! Stop timing only after the complete sediment calculation for this
      ! routing period has finished.  Aggregate worker timings here so
      ! SED_PERF reports actual min/mean/max costs rather than pre-compute zeros.
      CALL system_clock(clk_total_end)
      IF (clk_rate > 0) THEN
         t_total = real(clk_total_end - clk_total_start, r8) / real(clk_rate, r8)
      ELSE
         t_total = -1._r8
      ENDIF

      ! timing_local:
      ! 1 total
      ! 2 yield
      ! 3 advection
      ! 4 input
      ! 5 exchange
      ! 6 layer redistribution
      ! 7 diagnostics
      ! 8 other/uninstrumented
      timing_local(1) = t_total
      timing_local(2) = t_yield
      timing_local(3) = t_adv
      timing_local(4) = t_input
      timing_local(5) = t_exchange
      timing_local(6) = t_layer
      timing_local(7) = t_diag
      timing_local(8) = t_total - &
         (t_yield + t_adv + t_input + t_exchange + t_layer + t_diag)

      timing_min = timing_local
      timing_max = timing_local
      timing_sum = timing_local

#ifdef USEMPI
      CALL mpi_allreduce(MPI_IN_PLACE, timing_min, size(timing_min), &
         MPI_REAL8, MPI_MIN, p_comm_worker, p_err)
      CALL mpi_allreduce(MPI_IN_PLACE, timing_max, size(timing_max), &
         MPI_REAL8, MPI_MAX, p_comm_worker, p_err)
      CALL mpi_allreduce(MPI_IN_PLACE, timing_sum, size(timing_sum), &
         MPI_REAL8, MPI_SUM, p_comm_worker, p_err)
#endif

      timing_mean = timing_sum / real(p_np_worker, r8)

      IF (p_iam_worker == 0) THEN
         ! Keep the legacy timing lines for backward-compatible log parsing.
         WRITE(*,'(A,I6,A,I6,A,F12.3,A,F12.3,A)') 'Sediment timing: morph_substeps=', iter_sed, &
            ', adv_substeps=', iter_adv, ', total=', t_total, ' s, routing_dt=', deltime, ' s'
         WRITE(*,'(A,6(F10.3,A))') 'Sediment timing detail [s]: yield=', t_yield, &
            ', adv=', t_adv, ', input=', t_input, ', exch=', t_exchange, &
            ', layer=', t_layer, ', diag=', t_diag

         WRITE(*,'(A,I0,A,I0,A,I0,A,F12.3)') &
            'SED_PERF substeps: morph=', iter_sed, &
            ', adv=', iter_adv, &
            ', workers=', p_np_worker, &
            ', routing_dt=', deltime
         WRITE(*,'(A,3F12.3)') &
            'SED_PERF total  min/mean/max [s] = ', &
            timing_min(1), timing_mean(1), timing_max(1)
         WRITE(*,'(A,3F12.3)') &
            'SED_PERF adv    min/mean/max [s] = ', &
            timing_min(3), timing_mean(3), timing_max(3)
         WRITE(*,'(A,3F12.3)') &
            'SED_PERF exch   min/mean/max [s] = ', &
            timing_min(5), timing_mean(5), timing_max(5)
         WRITE(*,'(A,3F12.3)') &
            'SED_PERF diag   min/mean/max [s] = ', &
            timing_min(7), timing_mean(7), timing_max(7)
         WRITE(*,'(A,3F12.3)') &
            'SED_PERF other  min/mean/max [s] = ', &
            timing_min(8), timing_mean(8), timing_max(8)
      ENDIF
#endif

      sum_layer_local = 0._r8
      sum_seddep_local = 0._r8
      sum_sedsto_local = 0._r8
      sum_bed_bulk_local = 0._r8
      sum_bed_solid_local = 0._r8
      n_wet_local = 0
      n_shallow_local = 0
      n_source_local = 0
      n_susp_local = 0
      n_bed_local = 0
      n_exchange_pos_local = 0
      n_exchange_neg_local = 0
      n_es_raw_local = 0
      n_d_raw_local = 0
      n_es_eff_local = 0
      n_d_eff_local = 0
      IF (numucat > 0) THEN
         sum_layer_local = sum(layer)
         sum_seddep_local = sum(seddep)
         sum_bed_bulk_local = sum_layer_local + sum_seddep_local
         sum_bed_solid_local = (1._r8 - lambda) * sum_bed_bulk_local
         sum_sedsto_local = sum(sedsto)
         IF (deltime > 0._r8) THEN
            sum_sedinp_local = sum_sedinp_local / deltime
            sum_sedout_down_local = sum_sedout_down_local / deltime
            sum_sedout_up_local = sum_sedout_up_local / deltime
            sum_sedout_abs_local = sum_sedout_abs_local / deltime
            sum_netflw_pos_local = sum_netflw_pos_local / deltime
            sum_netflw_neg_local = sum_netflw_neg_local / deltime
            sum_es_raw_local = sum_es_raw_local / deltime
            sum_d_raw_local = sum_d_raw_local / deltime
            sum_es_eff_local = sum_es_eff_local / deltime
            sum_d_eff_local = sum_d_eff_local / deltime
         ENDIF
         n_wet_local = count(wet_seen)
         n_shallow_local = count(shallow_seen)
         n_source_local = count(source_seen)
         n_susp_local = count(susp_seen)
         n_bed_local = count(bed_seen)
         n_exchange_pos_local = count(exch_pos_seen)
         n_exchange_neg_local = count(exch_neg_seen)
         n_es_raw_local = count(es_raw_seen)
         n_d_raw_local = count(d_raw_seen)
         n_es_eff_local = count(es_eff_seen)
         n_d_eff_local = count(d_eff_seen)
      ENDIF
      diag_max_global = (/ max_sedcon_local, max_sedout_local, max_bedout_local, &
         max_sedinp_local, max_netflw_local, max_shearvel_local, max_flow_cancel_local, &
         max_es_raw_local, max_d_raw_local, max_es_eff_local, max_d_eff_local /)
      diag_sum_global = (/ sum_layer_local, sum_seddep_local, sum_sedsto_local, &
         sum_sedinp_local, sum_sedout_down_local, sum_sedout_up_local, &
         sum_sedout_abs_local, sum_netflw_pos_local, sum_netflw_neg_local, &
         sum_bed_bulk_local, sum_bed_solid_local, sum_rivout_signed_local, &
         sum_rivout_abs_local, sum_es_raw_local, sum_d_raw_local, &
         sum_es_eff_local, sum_d_eff_local /)
      diag_count_global = (/ n_wet_local, n_shallow_local, n_source_local, &
         n_susp_local, n_bed_local, n_exchange_pos_local, n_exchange_neg_local, &
         n_es_raw_local, n_d_raw_local, n_es_eff_local, n_d_eff_local, &
         n_flow_cancel_local /)
      carrier_filter_diag_global = (/ sum_inst_near_dry_abs_q_local, &
         sum_period_near_dry_abs_q_local, sum_period_near_dry_signed_q_local /)
      carrier_filter_count_global = (/ n_period_near_dry_local /)
#ifdef CoLMDEBUG
      station_state_local = 0._r8
      station_state_global = 0._r8
      station_sedcon_local = 0._r8
      station_sedcon_global = 0._r8
      station_sedout_local = 0._r8
      station_sedout_global = 0._r8
      DO ista = 1, SED_N_DIAG_STATIONS
         station_i = 0
         IF (allocated(sed_diag_station_local_i)) station_i = sed_diag_station_local_i(ista)
         IF (station_i > 0) THEN
            station_state_local(ista,1) = rivout(station_i)
            station_state_local(ista,2) = rivout_abs(station_i)
            station_state_local(ista,3) = rivsto(station_i)
            station_state_local(ista,4) = shearvel(station_i)
            station_sedcon_local(:,ista) = sedcon(:,station_i)
            station_sedout_local(:,ista) = sedout(:,station_i)
         ENDIF
      ENDDO
#ifdef USEMPI
      CALL mpi_allreduce(MPI_IN_PLACE, diag_max_global, size(diag_max_global), &
         MPI_REAL8, MPI_MAX, p_comm_worker, p_err)
#endif

      ! Select one deterministic owner for each global extreme. If the same
      ! maximum occurs in multiple cells, use the smallest global ucat ID.
      extreme_cell_global = huge(1)
      IF (extreme_meta_local(1,1) > 0 .and. max_sedcon_local == diag_max_global(1)) &
         extreme_cell_global(1) = extreme_meta_local(1,1)
      IF (extreme_meta_local(2,1) > 0 .and. max_sedout_local == diag_max_global(2)) &
         extreme_cell_global(2) = extreme_meta_local(2,1)
#ifdef USEMPI
      CALL mpi_allreduce(MPI_IN_PLACE, extreme_cell_global, size(extreme_cell_global), &
         MPI_INTEGER, MPI_MIN, p_comm_worker, p_err)
#endif

      DO iextreme = 1, 2
         IF (extreme_meta_local(iextreme,1) /= extreme_cell_global(iextreme)) THEN
            extreme_meta_local(iextreme,:) = 0
            extreme_state_local(iextreme,:) = 0._r8
            extreme_sedcon_local(:,iextreme) = 0._r8
            extreme_sedout_local(:,iextreme) = 0._r8
            extreme_sedinp_local(:,iextreme) = 0._r8
            extreme_netflw_local(:,iextreme) = 0._r8
         ENDIF
      ENDDO

      extreme_meta_global = extreme_meta_local
      extreme_state_global = extreme_state_local
      extreme_sedcon_global = extreme_sedcon_local
      extreme_sedout_global = extreme_sedout_local
      extreme_sedinp_global = extreme_sedinp_local
      extreme_netflw_global = extreme_netflw_local
#ifdef USEMPI
      CALL mpi_allreduce(MPI_IN_PLACE, extreme_meta_global, size(extreme_meta_global), &
         MPI_INTEGER, MPI_SUM, p_comm_worker, p_err)
      CALL mpi_allreduce(MPI_IN_PLACE, extreme_state_global, size(extreme_state_global), &
         MPI_REAL8, MPI_SUM, p_comm_worker, p_err)
      CALL mpi_allreduce(MPI_IN_PLACE, extreme_sedcon_global, size(extreme_sedcon_global), &
         MPI_REAL8, MPI_SUM, p_comm_worker, p_err)
      CALL mpi_allreduce(MPI_IN_PLACE, extreme_sedout_global, size(extreme_sedout_global), &
         MPI_REAL8, MPI_SUM, p_comm_worker, p_err)
      CALL mpi_allreduce(MPI_IN_PLACE, extreme_sedinp_global, size(extreme_sedinp_global), &
         MPI_REAL8, MPI_SUM, p_comm_worker, p_err)
      CALL mpi_allreduce(MPI_IN_PLACE, extreme_netflw_global, size(extreme_netflw_global), &
         MPI_REAL8, MPI_SUM, p_comm_worker, p_err)
      station_state_global = station_state_local
      station_sedcon_global = station_sedcon_local
      station_sedout_global = station_sedout_local
      CALL mpi_allreduce(MPI_IN_PLACE, station_state_global, size(station_state_global), &
         MPI_REAL8, MPI_SUM, p_comm_worker, p_err)
      CALL mpi_allreduce(MPI_IN_PLACE, station_sedcon_global, size(station_sedcon_global), &
         MPI_REAL8, MPI_SUM, p_comm_worker, p_err)
      CALL mpi_allreduce(MPI_IN_PLACE, station_sedout_global, size(station_sedout_global), &
         MPI_REAL8, MPI_SUM, p_comm_worker, p_err)
      CALL mpi_allreduce(MPI_IN_PLACE, diag_sum_global, size(diag_sum_global), &
         MPI_REAL8, MPI_SUM, p_comm_worker, p_err)
      CALL mpi_allreduce(MPI_IN_PLACE, diag_count_global, size(diag_count_global), &
         MPI_INTEGER, MPI_SUM, p_comm_worker, p_err)
      CALL mpi_allreduce(MPI_IN_PLACE, carrier_filter_diag_global, &
         size(carrier_filter_diag_global), MPI_REAL8, MPI_SUM, p_comm_worker, p_err)
      CALL mpi_allreduce(MPI_IN_PLACE, carrier_filter_count_global, &
         size(carrier_filter_count_global), MPI_INTEGER, MPI_SUM, p_comm_worker, p_err)
#else
      station_state_global = station_state_local
      station_sedcon_global = station_sedcon_local
      station_sedout_global = station_sedout_local
#endif

      ! Diagnostic summary of global sediment state (worker 0 only)
      IF (p_iam_worker == 0) THEN
         WRITE(*,'(A,ES10.3,A,ES10.3,A,ES10.3)') &
            'Sediment diag: max_sedcon=', diag_max_global(1), &
            ', max_sedout=', diag_max_global(2), ', max_bedout=', diag_max_global(3)
         WRITE(*,'(A,ES10.3,A,ES10.3,A,ES10.3)') &
            'Sediment diag: max_sedinp=', diag_max_global(4), &
            ', max_netflw=', diag_max_global(5), ', max_shearvel=', diag_max_global(6)
         WRITE(*,'(A,ES10.3,A,ES10.3,A,ES10.3)') &
            'Sediment diag: sum_layer=', diag_sum_global(1), &
            ', sum_seddep=', diag_sum_global(2), ', sum_sedsto=', diag_sum_global(3)
         WRITE(*,'(A,ES10.3,A,ES10.3)') &
            'Sediment diag: sum_bed_bulk=', diag_sum_global(10), &
            ', sum_bed_solid=', diag_sum_global(11)
         WRITE(*,'(A,ES10.3,A,ES10.3,A,ES10.3,A,I9)') &
            'Sediment flow diag: sum_rivout_signed=', diag_sum_global(12), &
            ', sum_rivout_abs=', diag_sum_global(13), &
            ', max_cancel=', diag_max_global(7), ', n_flow_cancel=', diag_count_global(12)
         WRITE(*,'(A,ES10.3,A,ES10.3,A,ES10.3,A,I9)') &
            'Sediment carrier filter diag: inst_near_dry_abs_q=', &
            carrier_filter_diag_global(1), ', period_near_dry_abs_q=', &
            carrier_filter_diag_global(2), ', period_near_dry_signed_q=', &
            carrier_filter_diag_global(3), ', n_period_near_dry=', &
            carrier_filter_count_global(1)
         WRITE(*,'(A,ES10.3,A,ES10.3,A,ES10.3,A,ES10.3,A,ES10.3)') &
            'Sediment diag: sum_sedinp=', diag_sum_global(4), &
            ', sum_sedout_down=', diag_sum_global(5), ', sum_sedout_up=', diag_sum_global(6), &
            ', sum_sedout_abs=', diag_sum_global(7), ', sum_netflw_pos=', diag_sum_global(8)
         WRITE(*,'(A,ES10.3)') 'Sediment diag: sum_netflw_neg=' , diag_sum_global(9)
         WRITE(*,'(A,ES10.3,A,ES10.3,A,ES10.3,A,ES10.3)') &
            'Sediment exchange diag raw: sum_Es=', diag_sum_global(14), ', sum_D=', diag_sum_global(15), &
            ', max_Es=', diag_max_global(8), ', max_D=', diag_max_global(9)
         WRITE(*,'(A,ES10.3,A,ES10.3,A,ES10.3,A,ES10.3)') &
            'Sediment exchange diag eff: sum_Es=', diag_sum_global(16), ', sum_D=', diag_sum_global(17), &
            ', max_Es=', diag_max_global(10), ', max_D=', diag_max_global(11)
         WRITE(*,'(A,I9,A,I9,A,I9,A,I9,A,I9,A,I9,A,I9)') &
            'Sediment counts: wet=', diag_count_global(1), ', shallow_wet=', diag_count_global(2), &
            ', source=', diag_count_global(3), ', susp=', diag_count_global(4), ', bed=', diag_count_global(5), &
            ', exch_pos=', diag_count_global(6), ', exch_neg=', diag_count_global(7)
         WRITE(*,'(A,I9,A,I9,A,I9,A,I9)') 'Sediment exchange counts raw: Es=', diag_count_global(8), &
            ', D=', diag_count_global(9), ', eff_Es=', diag_count_global(10), ', eff_D=', diag_count_global(11)

         DO iextreme = 1, 2
            IF (extreme_meta_global(iextreme,1) <= 0) CYCLE
            extreme_lon = -180._r8 + (real(extreme_meta_global(iextreme,2),r8) - 0.5_r8) &
               * 360._r8 / real(griducat%nlon,r8)
            extreme_lat = 90._r8 - (real(extreme_meta_global(iextreme,3),r8) - 0.5_r8) &
               * 180._r8 / real(griducat%nlat,r8)
            sedcon_total_diag = sum(extreme_sedcon_global(:,iextreme))
            sedout_total_diag = sum(extreme_sedout_global(:,iextreme))
            sedout_abs_total_diag = sum(abs(extreme_sedout_global(:,iextreme)))
            ssc_mg_l_diag = sedcon_total_diag * psedD * 1.e6_r8
            ssl_t_day_diag = sedout_total_diag * psedD * 86400._r8
            ssl_abs_t_day_diag = sedout_abs_total_diag * psedD * 86400._r8

            IF (iextreme == 1) THEN
               extreme_trigger_value = diag_max_global(1)
               WRITE(*,'(A)') 'Sediment extreme MAX_SEDCON:'
               WRITE(*,'(A,ES12.4,A,I0)') '  trigger_value[m3/m3]=', extreme_trigger_value, &
                  ', trigger_class=', extreme_meta_global(iextreme,4)
            ELSE
               extreme_trigger_value = diag_max_global(2)
               WRITE(*,'(A)') 'Sediment extreme MAX_SEDOUT:'
               WRITE(*,'(A,ES12.4,A,I0)') '  trigger_value[m3/s]=', extreme_trigger_value, &
                  ', trigger_class=', extreme_meta_global(iextreme,4)
            ENDIF
            WRITE(*,'(A,I0,A,I0,A,I0,A,F11.5,A,F10.5)') &
               '  ucat=', extreme_meta_global(iextreme,1), &
               ', x=', extreme_meta_global(iextreme,2), ', y=', extreme_meta_global(iextreme,3), &
               ', lon=', extreme_lon, ', lat=', extreme_lat
            WRITE(*,'(A,3(ES12.4,A))') '  depth avg/min/max [m]=', &
               extreme_state_global(iextreme,1), ' ', extreme_state_global(iextreme,2), ' ', &
               extreme_state_global(iextreme,3), ''
            WRITE(*,'(A,3(ES12.4,A))') '  rivsto[m3], rivout[m3/s], abs_rivout[m3/s]=', &
               extreme_state_global(iextreme,4), ' ', extreme_state_global(iextreme,5), ' ', &
               extreme_state_global(iextreme,6), ''
            WRITE(*,'(A,4(ES12.4,A))') '  rivout_min/max[m3/s], shearvel[m/s], bed_area[m2]=', &
               extreme_state_global(iextreme,7), ' ', extreme_state_global(iextreme,8), ' ', &
               extreme_state_global(iextreme,9), ' ', extreme_state_global(iextreme,10), ''
            WRITE(*,'(A)',advance='no') '  sedcon_by_class[m3/m3]='
            DO ised = 1, nsed
               WRITE(*,'(1X,ES12.4)',advance='no') extreme_sedcon_global(ised,iextreme)
            ENDDO
            WRITE(*,*)
            WRITE(*,'(A)',advance='no') '  sedout_by_class[m3/s]='
            DO ised = 1, nsed
               WRITE(*,'(1X,ES12.4)',advance='no') extreme_sedout_global(ised,iextreme)
            ENDDO
            WRITE(*,*)
            WRITE(*,'(A)',advance='no') '  sedinp_by_class[m3/s]='
            DO ised = 1, nsed
               WRITE(*,'(1X,ES12.4)',advance='no') extreme_sedinp_global(ised,iextreme)
            ENDDO
            WRITE(*,*)
            WRITE(*,'(A)',advance='no') '  netflw_by_class[m3/s]='
            DO ised = 1, nsed
               WRITE(*,'(1X,ES12.4)',advance='no') extreme_netflw_global(ised,iextreme)
            ENDDO
            WRITE(*,*)
            WRITE(*,'(A,ES12.4,A,ES12.4)') '  sedcon_total[m3/m3]=', sedcon_total_diag, &
               ', SSC_total[mg/L]=', ssc_mg_l_diag
            WRITE(*,'(A,ES12.4,A,ES12.4,A,ES12.4,A,ES12.4)') &
               '  sedout_total_signed[m3/s]=', sedout_total_diag, &
               ', sedout_total_abs[m3/s]=', sedout_abs_total_diag, &
               ', SSL_signed[t/day]=', ssl_t_day_diag, ', SSL_abs[t/day]=', ssl_abs_t_day_diag
         ENDDO

         WRITE(*,'(A)') 'Sediment Amazon benchmark station diag:'
         DO ista = 1, SED_N_DIAG_STATIONS
            station_sedcon_total = sum(station_sedcon_global(:,ista))
            station_sedout_total = sum(station_sedout_global(:,ista))
            station_sedout_abs_total = sum(abs(station_sedout_global(:,ista)))
            station_ssc_mg_l = station_sedcon_total * psedD * 1.e6_r8
            station_ssl_t_day = station_sedout_total * psedD * 86400._r8
            station_ssl_abs_t_day = station_sedout_abs_total * psedD * 86400._r8

            WRITE(*,'(2X,A16,A,I0,A,F8.3,A,2(F10.4,1X))') &
               trim(SED_DIAG_STATION_NAMES(ista)), ' ucat=', sed_diag_station_ucid(ista), &
               ' match_km=', sed_diag_station_distance_km(ista), ' model(lat lon)=', &
               sed_diag_station_model_lat(ista), sed_diag_station_model_lon(ista)
            WRITE(*,'(A,ES12.4,A,ES12.4,A,ES12.4,A,ES12.4)') &
               '    Q_signed[m3/s]=', station_state_global(ista,1), &
               ', Q_abs[m3/s]=', station_state_global(ista,2), &
               ', rivsto[m3]=', station_state_global(ista,3), &
               ', shearvel[m/s]=', station_state_global(ista,4)
            WRITE(*,'(A)',advance='no') '    sedcon_by_class[m3/m3]='
            DO ised = 1, nsed
               WRITE(*,'(1X,ES12.4)',advance='no') station_sedcon_global(ised,ista)
            ENDDO
            WRITE(*,*)
            WRITE(*,'(A)',advance='no') '    sedout_by_class[m3/s]='
            DO ised = 1, nsed
               WRITE(*,'(1X,ES12.4)',advance='no') station_sedout_global(ised,ista)
            ENDDO
            WRITE(*,*)
            WRITE(*,'(A,ES12.4,A,ES12.4,A,ES12.4)') &
               '    sedcon_total[m3/m3]=', station_sedcon_total, &
               ', SSC_total[mg/L]=', station_ssc_mg_l, &
               ', sedout_total[m3/s]=', station_sedout_total
            WRITE(*,'(A,ES12.4,A,ES12.4,A,ES12.4)') &
               '    sedout_abs_total[m3/s]=', station_sedout_abs_total, &
               ', SSL_signed[t/day]=', station_ssl_t_day, &
               ', SSL_abs[t/day]=', station_ssl_abs_t_day
         ENDDO
      ENDIF
#endif

#ifdef CoLMDEBUG
      deallocate(extreme_sedcon_local, extreme_sedcon_global, extreme_sedout_local, extreme_sedout_global, &
         extreme_sedinp_local, extreme_sedinp_global, extreme_netflw_local, extreme_netflw_global, &
         station_sedcon_local, station_sedcon_global, station_sedout_local, station_sedout_global)
#endif
      deallocate(rivsto, rivout, rivout_abs, bed_area, fldfrc, wet_seen, shallow_seen, source_seen, &
         susp_seen, bed_seen, exch_pos_seen, exch_neg_seen, es_raw_seen, &
         d_raw_seen, es_eff_seen, d_eff_seen)

   END SUBROUTINE grid_sediment_calc

   !-------------------------------------------------------------------------------------
   SUBROUTINE accumulate_sediment_output(dt)
   !-------------------------------------------------------------------------------------
   USE MOD_Grid_RiverLakeNetwork, only: numucat
   IMPLICIT NONE
   real(r8), intent(in) :: dt
   integer :: i

      IF (.not. p_is_worker) RETURN
      IF (numucat <= 0) RETURN

      DO i = 1, numucat
         a_sedcon(:,i) = a_sedcon(:,i) + sedcon(:,i) * dt
         a_sedout(:,i) = a_sedout(:,i) + sedout(:,i) * dt
         a_bedout(:,i) = a_bedout(:,i) + bedout(:,i) * dt
         a_sedinp(:,i) = a_sedinp(:,i) + sedinp(:,i) * dt
         a_netflw(:,i) = a_netflw(:,i) + &
            (netflw(:,i) + netflw_adv_step(:,i)) * dt
         a_layer(:,i)  = a_layer(:,i)  + layer(:,i)  * dt
         a_shearvel(i) = a_shearvel(i) + shearvel(i) * dt
      ENDDO

   END SUBROUTINE accumulate_sediment_output

   !-------------------------------------------------------------------------------------
   SUBROUTINE debug_check_sediment_fields(context, iter_sed, iter_adv, dt_morph, dt_adv, &
      dt_cfl_global, rivout, rivout_abs, rivsto, found)
   !-------------------------------------------------------------------------------------
   IMPLICIT NONE
   character(len=*), intent(in) :: context
   integer, intent(in) :: iter_sed, iter_adv
   real(r8), intent(in) :: dt_morph, dt_adv, dt_cfl_global
   real(r8), intent(in) :: rivout(:), rivout_abs(:), rivsto(:)
   logical, intent(out) :: found
   logical :: halt_invalid, halt_zero, halt_overflow

      found = .false.

      CALL ieee_get_halting_mode(ieee_invalid, halt_invalid)
      CALL ieee_get_halting_mode(ieee_divide_by_zero, halt_zero)
      CALL ieee_get_halting_mode(ieee_overflow, halt_overflow)
      CALL ieee_set_halting_mode(ieee_invalid, .false.)
      CALL ieee_set_halting_mode(ieee_divide_by_zero, .false.)
      CALL ieee_set_halting_mode(ieee_overflow, .false.)

      CALL debug_check_sediment_field(context, 'sedout', sedout, iter_sed, iter_adv, &
         dt_morph, dt_adv, dt_cfl_global, rivout, rivout_abs, rivsto, found)
      CALL debug_check_sediment_field(context, 'bedout', bedout, iter_sed, iter_adv, &
         dt_morph, dt_adv, dt_cfl_global, rivout, rivout_abs, rivsto, found)
      CALL debug_check_sediment_field(context, 'sedcon', sedcon, iter_sed, iter_adv, &
         dt_morph, dt_adv, dt_cfl_global, rivout, rivout_abs, rivsto, found)
      CALL debug_check_sediment_field(context, 'sedsto', sedsto, iter_sed, iter_adv, &
         dt_morph, dt_adv, dt_cfl_global, rivout, rivout_abs, rivsto, found)
      CALL debug_check_sediment_field(context, 'layer', layer, iter_sed, iter_adv, &
         dt_morph, dt_adv, dt_cfl_global, rivout, rivout_abs, rivsto, found)
      CALL debug_check_sediment_field(context, 'netflw', netflw, iter_sed, iter_adv, &
         dt_morph, dt_adv, dt_cfl_global, rivout, rivout_abs, rivsto, found)

      CALL ieee_set_halting_mode(ieee_invalid, halt_invalid)
      CALL ieee_set_halting_mode(ieee_divide_by_zero, halt_zero)
      CALL ieee_set_halting_mode(ieee_overflow, halt_overflow)

   END SUBROUTINE debug_check_sediment_fields

   !-------------------------------------------------------------------------------------
   SUBROUTINE debug_check_sediment_field(context, field_name, field, iter_sed, iter_adv, &
      dt_morph, dt_adv, dt_cfl_global, rivout, rivout_abs, rivsto, found)
   !-------------------------------------------------------------------------------------
   USE MOD_Grid_RiverLakeNetwork, only: numucat, ucat_next
   IMPLICIT NONE
   character(len=*), intent(in) :: context, field_name
   real(r8), intent(in) :: field(:,:)
   integer, intent(in) :: iter_sed, iter_adv
   real(r8), intent(in) :: dt_morph, dt_adv, dt_cfl_global
   real(r8), intent(in) :: rivout(:), rivout_abs(:), rivsto(:)
   logical, intent(inout) :: found
   integer :: i, ised

      IF (found) RETURN
      DO i = 1, numucat
         DO ised = 1, nsed
            IF (.not. ieee_is_finite(field(ised,i))) THEN
               found = .true.
               WRITE(*,'(A)') '========== SEDIMENT_FINITE_DEBUG =========='
               WRITE(*,'(A,A)') 'context            = ', trim(context)
               WRITE(*,'(A,A)') 'bad field          = ', trim(field_name)
               WRITE(*,'(A,I0)') 'worker             = ', p_iam_worker
               WRITE(*,'(A,I0)') 'cell i             = ', i
               WRITE(*,'(A,I0)') 'sediment class     = ', ised
               WRITE(*,'(A,I0)') 'ucat_next          = ', ucat_next(i)
               WRITE(*,'(A,I0)') 'iter_sed           = ', iter_sed
               WRITE(*,'(A,I0)') 'iter_adv           = ', iter_adv
               WRITE(*,'(A,ES20.10)') 'dt_morph [s]       = ', dt_morph
               WRITE(*,'(A,ES20.10)') 'dt_adv [s]         = ', dt_adv
               WRITE(*,'(A,ES20.10)') 'dt_cfl_global [s]  = ', dt_cfl_global
               WRITE(*,'(A,ES20.10)') 'bad value          = ', field(ised,i)
               WRITE(*,'(A,ES20.10)') 'sedout             = ', sedout(ised,i)
               WRITE(*,'(A,ES20.10)') 'bedout             = ', bedout(ised,i)
               WRITE(*,'(A,ES20.10)') 'sedcon             = ', sedcon(ised,i)
               WRITE(*,'(A,ES20.10)') 'sedsto [m3]        = ', sedsto(ised,i)
               WRITE(*,'(A,ES20.10)') 'layer              = ', layer(ised,i)
               WRITE(*,'(A,ES20.10)') 'netflw             = ', netflw(ised,i)
               WRITE(*,'(A,ES20.10)') 'netflw_adv_step    = ', netflw_adv_step(ised,i)
               WRITE(*,'(A,ES20.10)') 'exch_d_adv_step    = ', exch_d_adv_step(ised,i)
               WRITE(*,'(A,ES20.10)') 'rivout [m3/s]      = ', rivout(i)
               WRITE(*,'(A,ES20.10)') 'rivout_abs [m3/s]  = ', rivout_abs(i)
               WRITE(*,'(A,ES20.10)') 'rivsto [m3]        = ', rivsto(i)
               WRITE(*,'(A,ES20.10)') 'shearvel           = ', shearvel(i)
               WRITE(*,'(A)') '=========================================='
               RETURN
            ENDIF
         ENDDO
      ENDDO

   END SUBROUTINE debug_check_sediment_field

   !-------------------------------------------------------------------------------------
   SUBROUTINE sediment_diag_accumulate(dt_all, irivsys, ucatfilter, veloc, wdsrf, rivsto_input, rivout_fc, floodarea)
   ! Accumulate water flow variables for sediment calculation.
   ! Called once per water sub-step with full arrays (not per-cell).
   !-------------------------------------------------------------------------------------
   USE MOD_Grid_RiverLakeNetwork, only: numucat
   IMPLICIT NONE
   real(r8), intent(in) :: dt_all(:)       ! Time step per river system [numrivsys]
   integer,  intent(in) :: irivsys(:)      ! Cell -> river system mapping [numucat]
   logical,  intent(in) :: ucatfilter(:)   ! Active cell mask [numucat]
   real(r8), intent(in) :: veloc(:)        ! River velocity [numucat]
   real(r8), intent(in) :: wdsrf(:)        ! Water depth [numucat]
   real(r8), intent(in) :: rivsto_input(:) ! HYDRO water storage [m3, numucat]
   real(r8), intent(in) :: rivout_fc(:)    ! Downstream face flux [numucat]
   real(r8), intent(in) :: floodarea(:)    ! Flooded area [m^2, numucat]
   integer  :: i
   real(r8) :: dt

      IF (.not. sediment_particle_enabled()) RETURN
      IF (.not. p_is_worker) RETURN
      IF (numucat <= 0) RETURN

      DO i = 1, numucat
         IF (.not. ucatfilter(i)) CYCLE
         IF (irivsys(i) < 1 .or. irivsys(i) > size(dt_all)) CYCLE
         dt = dt_all(irivsys(i))
         sed_acc_time(i)      = sed_acc_time(i)      + dt
         sed_acc_wdsrf(i)     = sed_acc_wdsrf(i)     + wdsrf(i)       * dt
         sed_acc_rivsto(i)    = sed_acc_rivsto(i)    + rivsto_input(i)* dt
         sed_acc_floodarea(i) = sed_acc_floodarea(i) + floodarea(i)   * dt
         sed_acc_wdsrf_min(i) = min(sed_acc_wdsrf_min(i), wdsrf(i))
         sed_acc_wdsrf_max(i) = max(sed_acc_wdsrf_max(i), wdsrf(i))
         sed_acc_rivsto_min(i) = min(sed_acc_rivsto_min(i), rivsto_input(i))
         sed_acc_rivsto_max(i) = max(sed_acc_rivsto_max(i), rivsto_input(i))
         sed_acc_rivout_min(i) = min(sed_acc_rivout_min(i), rivout_fc(i))
         sed_acc_rivout_max(i) = max(sed_acc_rivout_max(i), rivout_fc(i))
         sed_acc_pos_rivout(i) = sed_acc_pos_rivout(i) + max(rivout_fc(i), 0._r8) * dt
         sed_acc_neg_rivout(i) = sed_acc_neg_rivout(i) + max(-rivout_fc(i), 0._r8) * dt
         ! Do not let a residual numerical water film carry sediment.
         ! The Amazon failure showed depths of O(1e-5--1e-16 m) paired
         ! with finite face fluxes, which drives V/|Q| -> 0 and destroys
         ! the explicit sediment-advection CFL timestep.
         IF (wdsrf(i) > SED_NEAR_DRY_DEPTH) THEN
            sed_acc_carrier_time(i) = sed_acc_carrier_time(i) + dt
            sed_acc_v2(i)         = sed_acc_v2(i)         + veloc(i)**2 * dt
            sed_acc_rivout(i)     = sed_acc_rivout(i)     + rivout_fc(i) * dt
            sed_acc_abs_rivout(i) = sed_acc_abs_rivout(i) + abs(rivout_fc(i)) * dt
         ELSE
            sed_acc_near_dry_abs_rivout(i) = sed_acc_near_dry_abs_rivout(i) &
               + abs(rivout_fc(i)) * dt
         ENDIF
      ENDDO

   END SUBROUTINE sediment_diag_accumulate

   !-------------------------------------------------------------------------------------
   SUBROUTINE sediment_forcing_put(precip, dt)
   !-------------------------------------------------------------------------------------
   ! Accumulate precipitation forcing for sediment yield calculation.
   ! The yield power-law term (rate_mm_hr)^pyldpc is accumulated per forcing step
   ! to avoid Jensen's inequality bias: <P^p> >= <P>^p for convex p>1.
   ! The precipitation threshold is NOT applied here; it is evaluated later using
   ! the routing-period mean rain rate so the trigger definition stays unchanged.
   !-------------------------------------------------------------------------------------
   USE MOD_Grid_RiverLakeNetwork, only: numucat
   IMPLICIT NONE
   real(r8), intent(in) :: precip(:)   ! precipitation rate [mm/s]
   real(r8), intent(in) :: dt          ! forcing time step [s]
   integer  :: i
   real(r8) :: rate_mm_hr

      IF (.not. sediment_particle_enabled()) RETURN
      IF (.not. p_is_worker) RETURN
      IF (numucat <= 0) RETURN

      DO i = 1, numucat
         sed_precip(i) = sed_precip(i) + max(0._r8, precip(i)) * dt      ! for diagnostics

         ! Accumulate yield power-law term per step; threshold is checked later
         ! from the routing-period mean rain rate.
         rate_mm_hr = max(0._r8, precip(i)) * 3600._r8
         sed_precip_yield(i) = sed_precip_yield(i) + rate_mm_hr**pyldpc * dt
      ENDDO
      sed_precip_time = sed_precip_time + dt

   END SUBROUTINE sediment_forcing_put

   !-------------------------------------------------------------------------------------
   SUBROUTINE parse_grain_diameters()
   !-------------------------------------------------------------------------------------
   IMPLICIT NONE
   character(len=256) :: str
   integer :: i, j, k, n
   integer :: iostat

      allocate(sDiam(nsed))
      IF (allocated(sDiam_from_param)) THEN
         sDiam(:) = sDiam_from_param(:)
         IF (p_is_io) WRITE(*,*) 'Grain diameters (m) from sediment tracer parameter file:', sDiam
         RETURN
      ENDIF
      str = trim(adjustl(SED_DEFAULT_DIAMETER))

      n = 0
      j = 1
      DO i = 1, len_trim(str)
         IF (str(i:i) == ',' .or. i == len_trim(str)) THEN
            n = n + 1
            IF (i == len_trim(str)) THEN
               k = i
            ELSE
               k = i - 1
            ENDIF
            IF (n <= nsed) THEN
               read(str(j:k), *, iostat=iostat) sDiam(n)
               IF (iostat /= 0) THEN
                  IF (p_is_io) THEN
                     WRITE(*,*) 'ERROR: Failed to parse grain diameter at position', n
                  ENDIF
                  CALL CoLM_stop()
               ENDIF
            ENDIF
            j = i + 1
         ENDIF
      ENDDO

      IF (n /= nsed) THEN
         IF (p_is_io) THEN
            WRITE(*,*) 'ERROR: Number of diameters does not match nsed:', n, nsed
         ENDIF
         CALL CoLM_stop()
      ENDIF

      DO i = 1, nsed
         IF (.not. ieee_is_finite(sDiam(i)) .or. sDiam(i) <= 0._r8) THEN
            IF (p_is_io) WRITE(*,*) 'ERROR: Grain diameter must be positive, class:', i
            CALL CoLM_stop()
         ENDIF
      ENDDO

      IF (p_is_io) WRITE(*,*) 'Grain diameters (m):', sDiam

   END SUBROUTINE parse_grain_diameters

   !-------------------------------------------------------------------------------------
   SUBROUTINE calc_settling_velocities()
   !-------------------------------------------------------------------------------------
   USE MOD_Const_Physical, only: grav
   IMPLICIT NONE
   real(r8) :: sTmp
   integer :: i

      allocate(setvel(nsed))

      DO i = 1, nsed
         sTmp = 6.0_r8 * visKin / sDiam(i)
         setvel(i) = pset * (sqrt(2.0_r8/3.0_r8 * (psedD-pwatD)/pwatD * grav * sDiam(i) &
                    + sTmp*sTmp) - sTmp)
      ENDDO

      IF (p_is_io) WRITE(*,*) 'Settling velocities (m/s):', setvel

   END SUBROUTINE calc_settling_velocities

   !-------------------------------------------------------------------------------------
   SUBROUTINE read_sediment_static_data(parafile)
   !-------------------------------------------------------------------------------------
   USE MOD_Grid_RiverLakeNetwork, only: numucat, readin_riverlake_parameter
   IMPLICIT NONE
   character(len=*), intent(in) :: parafile

      CALL readin_riverlake_parameter(parafile, 'sed_frc', rdata2d=sed_frc)
      CALL readin_riverlake_parameter(parafile, 'sed_slope', rdata2d=sed_slope)

      ! Validate dimensions
      IF (p_is_worker .and. numucat > 0) THEN
         IF (size(sed_frc,1) /= nsed) THEN
            IF (p_is_io) WRITE(*,*) 'ERROR: sed_frc dim1 =', size(sed_frc,1), ' expected nsed =', nsed
            CALL CoLM_stop()
         ENDIF
         IF (size(sed_slope,1) /= nlfp_sed) THEN
            IF (p_is_io) WRITE(*,*) 'ERROR: sed_slope dim1 =', size(sed_slope,1), ' expected nlfp_sed =', nlfp_sed
            CALL CoLM_stop()
         ENDIF
      ENDIF

      ! Validate distributed values before MAX/normalization can mask invalid
      ! negative or non-finite input from the static-data file.
      CALL validate_sediment_parameters()
      IF (allocated(sed_slope)) sed_slope = max(sed_slope, 0._r8)
      CALL normalize_sed_frc()

      IF (p_is_io) WRITE(*,*) 'Sediment static data read successfully.'

   END SUBROUTINE read_sediment_static_data

   !-------------------------------------------------------------------------------------
   SUBROUTINE normalize_sed_frc()
   !-------------------------------------------------------------------------------------
   USE MOD_Grid_RiverLakeNetwork, only: numucat
   IMPLICIT NONE
   integer :: i
   real(r8) :: frc_sum

      IF (.not. p_is_worker) RETURN
      IF (numucat <= 0) RETURN

      sed_frc = max(sed_frc, 0._r8)
      DO i = 1, numucat
         frc_sum = sum(sed_frc(:,i))
         IF (frc_sum > 0._r8) THEN
            sed_frc(:,i) = sed_frc(:,i) / frc_sum
         ELSE
            sed_frc(:,i) = 1._r8 / real(nsed, r8)
         ENDIF
      ENDDO

   END SUBROUTINE normalize_sed_frc

   !-------------------------------------------------------------------------------------
   SUBROUTINE allocate_sediment_vars()
   !-------------------------------------------------------------------------------------
   USE MOD_Grid_RiverLakeNetwork, only: numucat
   IMPLICIT NONE

      IF (.not. p_is_worker) RETURN

      ! Allocate on ALL workers, including numucat=0 (zero-length arrays).
      ! This is required because calc_sediment_advection passes module arrays
      ! as array sections to worker_push_data, which all workers must call.
      allocate(sedcon(nsed, numucat))
      allocate(sedsto(nsed, numucat))
      allocate(layer (nsed, numucat))
      allocate(seddep(nsed, totlyrnum, numucat))

      allocate(sedout      (nsed, numucat))
      allocate(bedout      (nsed, numucat))
      allocate(sedinp      (nsed, numucat))
      allocate(netflw      (nsed, numucat))
      allocate(exch_es_raw (nsed, numucat))
      allocate(exch_d_raw  (nsed, numucat))
      allocate(exch_es_eff (nsed, numucat))
      allocate(exch_d_eff  (nsed, numucat))
      allocate(netflw_adv_step(nsed, numucat))
      allocate(exch_d_adv_step(nsed, numucat))
      allocate(shearvel    (numucat))
      allocate(critshearvel(nsed, numucat))
      allocate(susvel      (nsed, numucat))

      allocate(sed_acc_time     (numucat))
      allocate(sed_acc_v2       (numucat))
      allocate(sed_acc_wdsrf    (numucat))
      allocate(sed_acc_rivsto   (numucat))
      allocate(sed_acc_rivout   (numucat))
      allocate(sed_acc_abs_rivout(numucat))
      allocate(sed_acc_floodarea(numucat))
      allocate(sed_acc_carrier_time(numucat))
      allocate(sed_acc_wdsrf_min(numucat))
      allocate(sed_acc_wdsrf_max(numucat))
      allocate(sed_acc_rivsto_min(numucat))
      allocate(sed_acc_rivsto_max(numucat))
      allocate(sed_acc_rivout_min(numucat))
      allocate(sed_acc_rivout_max(numucat))
      allocate(sed_acc_pos_rivout(numucat))
      allocate(sed_acc_neg_rivout(numucat))
      allocate(sed_acc_near_dry_abs_rivout(numucat))
      allocate(sed_precip       (numucat))
      allocate(sed_precip_yield (numucat))

      allocate(a_sedcon  (nsed, numucat))
      allocate(a_sedout  (nsed, numucat))
      allocate(a_bedout  (nsed, numucat))
      allocate(a_sedinp  (nsed, numucat))
      allocate(a_netflw  (nsed, numucat))
      allocate(a_layer   (nsed, numucat))
      allocate(a_shearvel(numucat))

      sedcon       = 0._r8;  sedsto       = 0._r8
      layer        = 0._r8;  seddep       = 0._r8
      sedout       = 0._r8;  bedout       = 0._r8;  sedinp       = 0._r8
      netflw       = 0._r8
      exch_es_raw  = 0._r8;  exch_d_raw   = 0._r8
      exch_es_eff  = 0._r8;  exch_d_eff   = 0._r8
      netflw_adv_step = 0._r8; exch_d_adv_step = 0._r8
      shearvel     = 0._r8;  critshearvel = 0._r8
      susvel       = 0._r8
      sed_acc_time  = 0._r8;  sed_acc_v2        = 0._r8
      sed_acc_wdsrf = 0._r8;  sed_acc_rivsto    = 0._r8
      sed_acc_rivout = 0._r8
      sed_acc_abs_rivout = 0._r8
      sed_acc_floodarea = 0._r8
      sed_acc_carrier_time = 0._r8
      sed_acc_wdsrf_min = huge(1._r8)
      sed_acc_wdsrf_max = 0._r8
      sed_acc_rivsto_min = huge(1._r8)
      sed_acc_rivsto_max = 0._r8
      sed_acc_rivout_min = huge(1._r8)
      sed_acc_rivout_max = -huge(1._r8)
      sed_acc_pos_rivout = 0._r8
      sed_acc_neg_rivout = 0._r8
      sed_acc_near_dry_abs_rivout = 0._r8
      sed_precip    = 0._r8;  sed_precip_yield = 0._r8
      sed_precip_time = 0._r8
      sed_hist_acctime = 0._r8
      a_sedcon     = 0._r8;  a_sedout     = 0._r8;  a_bedout     = 0._r8
      a_sedinp     = 0._r8;  a_netflw     = 0._r8;  a_layer      = 0._r8
      a_shearvel   = 0._r8

   END SUBROUTINE allocate_sediment_vars

   !-------------------------------------------------------------------------------------
   SUBROUTINE initialize_sediment_state()
   !-------------------------------------------------------------------------------------
   USE MOD_Grid_RiverLakeNetwork, only: numucat, topo_rivwth, topo_rivlen
   IMPLICIT NONE
   integer :: i, ilyr

      IF (.not. p_is_worker) RETURN
      IF (numucat <= 0) RETURN

      DO i = 1, numucat
         layer(:,i) = lyrdph * topo_rivwth(i) * topo_rivlen(i) * sed_frc(:,i)
         DO ilyr = 1, totlyrnum - 1
            seddep(:,ilyr,i) = layer(:,i)
         ENDDO
         seddep(:,totlyrnum,i) = max(sed_bed_depth - lyrdph*totlyrnum, 0._r8) &
            * topo_rivwth(i) * topo_rivlen(i) * sed_frc(:,i)
      ENDDO

   END SUBROUTINE initialize_sediment_state

   !-------------------------------------------------------------------------------------
   FUNCTION calc_critical_shear_vel_sq(diam) RESULT(csvel_sq)
   ! Return the SQUARE of critical shear velocity in (cm/s)^2.
   ! Callers apply sqrt() and * 0.01 to obtain the velocity in m/s.
   ! Matches CaMa-Flood's calc_criticalShearVelocity (see its comment: "[(cm/s)^2]").
   !-------------------------------------------------------------------------------------
   IMPLICIT NONE
   real(r8), intent(in) :: diam
   real(r8) :: csvel_sq
   real(r8) :: cA, cB

      cB = 1._r8
      IF (diam >= 0.00303_r8) THEN
         cA = 80.9_r8
      ELSEIF (diam >= 0.00118_r8) THEN
         cA = 134.6_r8;  cB = 31._r8 / 32._r8
      ELSEIF (diam >= 0.000565_r8) THEN
         cA = 55._r8
      ELSEIF (diam >= 0.000065_r8) THEN
         cA = 8.41_r8;   cB = 11._r8 / 32._r8
      ELSE
         cA = 226._r8
      ENDIF

      csvel_sq = cA * (diam * 100._r8) ** cB

   END FUNCTION calc_critical_shear_vel_sq

   !-------------------------------------------------------------------------------------
   SUBROUTINE calc_critical_shear_egiazoroff(i, svel, csvel_out)
   !-------------------------------------------------------------------------------------
   IMPLICIT NONE
   integer,  intent(in)  :: i
   real(r8), intent(in)  :: svel
   real(r8), intent(out) :: csvel_out(nsed)
   real(r8) :: dMean, csVel0_sq, layer_sum
   integer  :: ised

      layer_sum = sum(layer(:,i))
      IF (layer_sum <= 0._r8) THEN
         csvel_out(:) = 1.e20_r8
         RETURN
      ENDIF

      dMean = 0._r8
      DO ised = 1, nsed
         dMean = dMean + sDiam(ised) * layer(ised,i) / layer_sum
      ENDDO

      csVel0_sq = calc_critical_shear_vel_sq(dMean)   ! (cm/s)^2

      DO ised = 1, nsed
         IF (sDiam(ised) / dMean >= 0.4_r8) THEN
            csvel_out(ised) = sqrt(csVel0_sq * sDiam(ised) / dMean) * &
               (log10(19._r8) / log10(19._r8 * sDiam(ised) / dMean)) * 0.01_r8
         ELSE
            csvel_out(ised) = sqrt(0.85_r8 * csVel0_sq) * 0.01_r8
         ENDIF
      ENDDO

   END SUBROUTINE calc_critical_shear_egiazoroff

   !-------------------------------------------------------------------------------------
   SUBROUTINE calc_suspend_velocity(csvel, svel, susvel_out)
   !-------------------------------------------------------------------------------------
   IMPLICIT NONE
   real(r8), intent(in)  :: csvel(nsed)
   real(r8), intent(in)  :: svel
   real(r8), intent(out) :: susvel_out(nsed)
   real(r8) :: alpha, a, cB, sTmp
   integer  :: ised

      alpha = vonKar / 6._r8
      a = 0.08_r8
      cB = 1._r8 - lambda
      susvel_out(:) = 0._r8

      DO ised = 1, nsed
         IF (csvel(ised) > svel) CYCLE
         IF (svel <= 0._r8) CYCLE
         sTmp = setvel(ised) / alpha / svel
         susvel_out(ised) = max(setvel(ised) * cB / (1._r8 + sTmp) * &
            (1._r8 - a*sTmp) / (1._r8 + (1._r8-a)*sTmp), 0._r8)
      ENDDO

   END SUBROUTINE calc_suspend_velocity

   !-------------------------------------------------------------------------------------
   SUBROUTINE begin_suspended_period(rivsto)
   ! Derive concentration once from canonical suspended solid volume and current
   ! HYDRO water storage.  Clean-water volume changes therefore dilute/concentrate
   ! particles without changing their volume.
   !-------------------------------------------------------------------------------------
   USE MOD_Grid_RiverLakeNetwork, only: numucat
   IMPLICIT NONE

   real(r8), intent(in) :: rivsto(:)
   integer :: i

      IF (.not. p_is_worker) RETURN

      IF (any(.not. (sedsto >= 0._r8))) THEN
         CALL CoLM_stop('invalid suspended sediment solid volume at period boundary')
      ENDIF

      DO i = 1, numucat
         IF (rivsto(i) > 0._r8) THEN
            sedcon(:,i) = sedsto(:,i) / rivsto(i)
         ELSE
            ! The shared mass remains intact until the first advection substep,
            ! which deposits it to the bed and books that transfer exactly once.
            sedcon(:,i) = 0._r8
         ENDIF
      ENDDO

   END SUBROUTINE begin_suspended_period

   !-------------------------------------------------------------------------------------
   SUBROUTINE commit_suspended_period(rivsto)
   ! Publish the diagnostic concentration derived from canonical solid volume.
   !-------------------------------------------------------------------------------------
   USE MOD_Grid_RiverLakeNetwork, only: numucat
   IMPLICIT NONE

   real(r8), intent(in) :: rivsto(:)
   integer :: i

      IF (.not. p_is_worker) RETURN

      DO i = 1, numucat
         IF (rivsto(i) > 0._r8) THEN
            sedcon(:,i) = sedsto(:,i) / rivsto(i)
         ELSE
            IF (any(sedsto(:,i) > 0._r8)) THEN
               CALL CoLM_stop('dry cell retained suspended sediment after advection')
            ENDIF
            sedcon(:,i) = 0._r8
         ENDIF
      ENDDO
   END SUBROUTINE commit_suspended_period

   !-------------------------------------------------------------------------------------
   SUBROUTINE calc_sediment_advection(dt, rivout_signed, rivout_abs, rivsto)
   ! Preserve forward and reverse transport volumes when the routing-period
   ! signed mean cancels.  The two directional means integrate to the supplied
   ! signed and absolute discharge diagnostics.
   IMPLICIT NONE

   real(r8), intent(in) :: dt
   real(r8), intent(in) :: rivout_signed(:), rivout_abs(:), rivsto(:)
   real(r8), allocatable :: rivout_forward(:), rivout_reverse(:)
   real(r8), allocatable :: sedout_first(:,:), bedout_first(:,:)
   real(r8), allocatable :: netflw_adv_first(:,:), exch_d_adv_first(:,:)

      allocate(rivout_forward(size(rivout_signed)), rivout_reverse(size(rivout_signed)))
      allocate(sedout_first(nsed, size(rivout_signed)), bedout_first(nsed, size(rivout_signed)))
      allocate(netflw_adv_first(nsed, size(rivout_signed)), exch_d_adv_first(nsed, size(rivout_signed)))

      rivout_forward = max(0._r8, 0.5_r8 * (rivout_abs + rivout_signed))
      rivout_reverse = min(0._r8, 0.5_r8 * (rivout_signed - rivout_abs))

      CALL calc_sediment_advection_one_direction(dt, rivout_forward, rivsto)
      sedout_first = sedout
      bedout_first = bedout
      netflw_adv_first = netflw_adv_step
      exch_d_adv_first = exch_d_adv_step

      CALL calc_sediment_advection_one_direction(dt, rivout_reverse, rivsto)
      sedout = sedout + sedout_first
      bedout = bedout + bedout_first
      netflw_adv_step = netflw_adv_step + netflw_adv_first
      exch_d_adv_step = exch_d_adv_step + exch_d_adv_first

      deallocate(rivout_forward, rivout_reverse, sedout_first, bedout_first, &
         netflw_adv_first, exch_d_adv_first)
   END SUBROUTINE calc_sediment_advection

   !-------------------------------------------------------------------------------------
   SUBROUTINE calc_sediment_advection_one_direction(dt, rivout, rivsto)
   ! Flux-based advection scheme: each cell computes its downstream face flux,
   ! then push_ups2ucat gathers upstream fluxes. Each cell updates its own storage.
   ! This correctly handles cross-MPI transport including reverse flow.
   !
   ! Sign convention for sedout/bedout:
   !   positive = sediment flows downstream (cell loses mass)
   !   negative = sediment flows upstream   (cell gains mass via reverse flow)
   !-------------------------------------------------------------------------------------
   USE MOD_Const_Physical, only: grav
   USE MOD_Grid_RiverLakeNetwork, only: numucat, ucat_next, topo_rivwth, &
      push_next2ucat, push_ups2ucat
   USE MOD_WorkerPushData
   IMPLICIT NONE

   real(r8), intent(in) :: dt
   real(r8), intent(in) :: rivout(:)
   real(r8), intent(in) :: rivsto(:)

   real(r8), allocatable :: sedcon_next(:,:), layer_next(:,:), critshearvel_next(:,:)
   real(r8), allocatable :: sed_ups(:,:), bed_ups(:,:)
   real(r8), allocatable :: avail_sto(:,:), avail_bed_solid(:,:)
   real(r8), allocatable :: shearvel_next(:), rivwth_next(:)
   real(r8), allocatable :: cell_mass_before(:)

   integer  :: i, ised
   real(r8) :: plusVel, minusVel, layer_sum, sedsto_sum
   real(r8) :: sedsto_neg_tol, layer_neg_tol
      real(r8) :: dTmp(nsed)

      IF (.not. p_is_worker) RETURN

      ! `netflw`/`exch_d_eff` are morphology-interval rates and are integrated
      ! over every CFL substep.  Cap and dry-cell deposition below belongs only
      ! to this CFL substep, so reset its separate contribution on entry.  If it
      ! were added to the persistent base rate, the first substep's deposition
      ! would be integrated again by every later substep.
      netflw_adv_step(:,:) = 0._r8
      exch_d_adv_step(:,:) = 0._r8

      allocate(sedcon_next (nsed, numucat))
      allocate(layer_next (nsed, numucat))
      allocate(critshearvel_next(nsed, numucat))
      allocate(sed_ups    (nsed, numucat))
      allocate(bed_ups    (nsed, numucat))
      allocate(shearvel_next(numucat))
      allocate(rivwth_next(numucat))
      allocate(cell_mass_before(numucat))

      ! Get downstream/source-cell state needed for reverse-flow upwind transport.
      DO ised = 1, nsed
         CALL worker_push_data(push_next2ucat, sedcon(ised,:), sedcon_next(ised,:), &
            fillvalue = 0._r8)
         CALL worker_push_data(push_next2ucat, layer(ised,:), layer_next(ised,:), &
            fillvalue = 0._r8)
         CALL worker_push_data(push_next2ucat, critshearvel(ised,:), critshearvel_next(ised,:), &
            fillvalue = 1.e20_r8)
      ENDDO
      CALL worker_push_data(push_next2ucat, shearvel, shearvel_next, fillvalue = 0._r8)
      CALL worker_push_data(push_next2ucat, topo_rivwth, rivwth_next, fillvalue = 0._r8)

      ! --- Step 1: Compute flux at each cell's downstream face ---
      DO i = 1, numucat
         ! Suspended sediment flux (upstream scheme)
         IF (rivout(i) >= 0._r8) THEN
            ! Forward flow: use own concentration
            sedout(:,i) = sedcon(:,i) * rivout(i)
         ELSE
            ! Reverse flow: use downstream cell's concentration
            sedout(:,i) = sedcon_next(:,i) * rivout(i)   ! negative
         ENDIF

         ! Bedload solid-volume flux using the upstream/source cell of the face.
         ! This uses an Ashida-Michiue-style shear-velocity form with coefficient 17,
         ! not the classic Meyer-Peter-Mueller coefficient-8 expression.
         ! Forward flow uses local bed state; reverse flow uses downstream bed state.
         bedout(:,i) = 0._r8
         IF (rivout(i) > 0._r8) THEN
            layer_sum = sum(layer(:,i))
            IF (.not. all(critshearvel(:,i) >= shearvel(i)) .and. layer_sum > 0._r8) THEN
               DO ised = 1, nsed
                  IF (critshearvel(ised,i) >= shearvel(i) .or. layer(ised,i) <= 0._r8) CYCLE
                  plusVel  = shearvel(i) + critshearvel(ised,i)
                  minusVel = shearvel(i) - critshearvel(ised,i)
                  bedout(ised,i) = SED_BEDLOAD_COEFF * topo_rivwth(i) * plusVel * minusVel * minusVel &
                     / ((psedD-pwatD)/pwatD) / grav * layer(ised,i) / layer_sum
               ENDDO
            ENDIF
         ELSEIF (rivout(i) < 0._r8) THEN
            layer_sum = sum(layer_next(:,i))
            IF (.not. all(critshearvel_next(:,i) >= shearvel_next(i)) .and. layer_sum > 0._r8) THEN
               DO ised = 1, nsed
                  IF (critshearvel_next(ised,i) >= shearvel_next(i) .or. layer_next(ised,i) <= 0._r8) CYCLE
                  plusVel  = shearvel_next(i) + critshearvel_next(ised,i)
                  minusVel = shearvel_next(i) - critshearvel_next(ised,i)
                  bedout(ised,i) = -SED_BEDLOAD_COEFF * rivwth_next(i) * plusVel * minusVel * minusVel &
                     / ((psedD-pwatD)/pwatD) / grav * layer_next(ised,i) / layer_sum
               ENDDO
            ENDIF
         ENDIF
      ENDDO

      ! --- Step 2a: Rate-limit FORWARD outflow (source = self, one edge per cell) ---
      ! Forward outflow is committed first; Step 2b will use the remaining storage
      ! for reverse extraction ("forward committed first" strategy).
      DO i = 1, numucat
         DO ised = 1, nsed
            IF (sedout(ised,i) > 0._r8) THEN
               IF (sedsto(ised,i) > 0._r8) THEN
                  sedout(ised,i) = min(sedout(ised,i), sedsto(ised,i) / dt)
               ELSE
                  sedout(ised,i) = 0._r8
               ENDIF
            ENDIF
            IF (bedout(ised,i) > 0._r8) THEN
               IF (layer(ised,i) > 0._r8) THEN
                  bedout(ised,i) = min(bedout(ised,i), (1._r8 - lambda) * layer(ised,i) / dt)
               ELSE
                  bedout(ised,i) = 0._r8
               ENDIF
            ENDIF
         ENDDO
      ENDDO

      ! --- Step 2b: Rate-limit REVERSE outflow (source = downstream cell) ---
      ! Multiple upstream cells may reverse-drain the same downstream source.
      ! Gather total reverse demand per source cell, compute scale factor, distribute back.
      ! All workers must call limit_reverse_flux (MPI communication inside).
      !
      ! "Forward committed first" strategy: the available supply at each source cell
      ! is local storage minus the forward outflow already committed in Step 2a.
      allocate(avail_sto(nsed, numucat))
      allocate(avail_bed_solid(nsed, numucat))
      DO i = 1, numucat
         DO ised = 1, nsed
            avail_sto(ised,i) = max(sedsto(ised,i) - max(sedout(ised,i), 0._r8) * dt, 0._r8)
            avail_bed_solid(ised,i) = max((1._r8 - lambda) * layer(ised,i) - max(bedout(ised,i), 0._r8) * dt, 0._r8)
         ENDDO
      ENDDO
      CALL limit_reverse_flux(sedout, avail_sto, dt)
      CALL limit_reverse_flux(bedout, avail_bed_solid, dt)
      deallocate(avail_sto, avail_bed_solid)

      ! --- Step 3: Gather upstream fluxes via MPI-safe communication ---
      DO ised = 1, nsed
         CALL worker_push_data(push_ups2ucat, sedout(ised,:), sed_ups(ised,:), &
            fillvalue = 0._r8, mode = 'sum')
         CALL worker_push_data(push_ups2ucat, bedout(ised,:), bed_ups(ised,:), &
            fillvalue = 0._r8, mode = 'sum')
      ENDDO

      ! --- Step 4: Update each cell's storage ---
      ! Net change = - own_downstream_flux + sum_of_upstream_fluxes
      cell_mass_before = sum(sedsto, dim=1) + (1._r8 - lambda) * sum(layer, dim=1)
      DO i = 1, numucat
         DO ised = 1, nsed
            sedsto(ised,i) = sedsto(ised,i) - sedout(ised,i) * dt + sed_ups(ised,i) * dt
            layer(ised,i) = layer(ised,i) + (-bedout(ised,i) + bed_ups(ised,i)) * dt / (1._r8 - lambda)
         ENDDO
      ENDDO

      ! --- Step 5: Safety clamp ---
      DO i = 1, numucat
         DO ised = 1, nsed

            sedsto_neg_tol = SED_BALANCE_ABS_TOL + SED_BALANCE_REL_TOL * &
               max(1._r8, abs(sedsto(ised,i)), abs(sedout(ised,i) * dt), &
               abs(sed_ups(ised,i) * dt))
            layer_neg_tol = SED_BALANCE_ABS_TOL + SED_BALANCE_REL_TOL * &
               max(1._r8, abs(layer(ised,i)), &
               abs(bedout(ised,i) * dt / (1._r8 - lambda)), &
               abs(bed_ups(ised,i) * dt / (1._r8 - lambda)))
            IF (sedsto(ised,i) < -sedsto_neg_tol .or. &
                layer(ised,i) < -layer_neg_tol) THEN
               WRITE(*,'(A)') '========== SEDIMENT_NEGATIVE_STORAGE =========='
               WRITE(*,'(A,I0)')       'worker          = ', p_iam_worker
               WRITE(*,'(A,I0)')       'cell i          = ', i
               WRITE(*,'(A,I0)')       'ucat_next       = ', ucat_next(i)
               WRITE(*,'(A,I0)')       'sediment class  = ', ised
               WRITE(*,'(A,ES24.14)')  'dt              = ', dt

               WRITE(*,'(A,ES24.14)')  'sedsto before clamp = ', sedsto(ised,i)
               WRITE(*,'(A,ES24.14)')  'layer before clamp  = ', layer(ised,i)

               WRITE(*,'(A,ES24.14)')  'sedout          = ', sedout(ised,i)
               WRITE(*,'(A,ES24.14)')  'sed_ups         = ', sed_ups(ised,i)
               WRITE(*,'(A,ES24.14)')  'bedout          = ', bedout(ised,i)
               WRITE(*,'(A,ES24.14)')  'bed_ups         = ', bed_ups(ised,i)

               WRITE(*,'(A,ES24.14)')  'rivsto          = ', rivsto(i)
               WRITE(*,'(A,ES24.14)')  'rivout          = ', rivout(i)

               WRITE(*,'(A)') '==============================================='
               FLUSH(6)
            ENDIF

            sedsto(ised,i) = max(sedsto(ised,i), 0._r8)
            layer(ised,i)  = max(layer(ised,i), 0._r8)

         ENDDO
      ENDDO

      ! --- Step 6: Update concentration; deposit stranded sediment in dry cells ---
      DO i = 1, numucat
         IF (rivsto(i) > 0._r8) THEN
            sedsto_sum = sum(sedsto(:,i))
            IF (sedsto_sum > rivsto(i) * sed_max_conc) THEN
               dTmp(:) = (sedsto_sum - rivsto(i) * sed_max_conc) * &
                  sedsto(:,i) / max(sedsto_sum, 1.e-20_r8)
               dTmp(:) = min(dTmp(:), sedsto(:,i))
               ! SEDIMENT_DRY_CAP_DEPOSIT_CREDIT: cap-induced deposition is a
               ! real suspended-to-bed transfer and must enter netflw/history.
               netflw_adv_step(:,i) = netflw_adv_step(:,i) - dTmp(:) / dt
               exch_d_adv_step(:,i) = exch_d_adv_step(:,i) + dTmp(:) / dt
               sedsto(:,i) = sedsto(:,i) - dTmp(:)
               layer(:,i) = layer(:,i) + dTmp(:) / (1._r8 - lambda)
            ENDIF
            sedcon(:,i) = sedsto(:,i) / rivsto(i)
         ELSE
            IF (sum(sedsto(:,i)) > 0._r8) THEN
               ! SEDIMENT_DRY_CAP_DEPOSIT_CREDIT: stranded dry-cell suspended
               ! material is deposited into the bed and credited diagnostically.
               netflw_adv_step(:,i) = netflw_adv_step(:,i) - sedsto(:,i) / dt
               exch_d_adv_step(:,i) = exch_d_adv_step(:,i) + sedsto(:,i) / dt
               layer(:,i) = layer(:,i) + sedsto(:,i) / (1._r8 - lambda)
            ENDIF
            sedcon(:,i) = 0._r8
            sedsto(:,i) = 0._r8
         ENDIF
         CALL assert_sediment_mass_balance('advection', i, cell_mass_before(i), &
            sum(sedsto(:,i)) + (1._r8 - lambda) * sum(layer(:,i)), &
            dt * sum(-sedout(:,i) + sed_ups(:,i) - bedout(:,i) + bed_ups(:,i)), &
            rivsto_i = rivsto(i), &
            rivout_i = rivout(i))
      ENDDO

      deallocate(sedcon_next, layer_next, critshearvel_next, sed_ups, bed_ups, &
         shearvel_next, rivwth_next, cell_mass_before)

   END SUBROUTINE calc_sediment_advection_one_direction

   !-------------------------------------------------------------------------------------
   SUBROUTINE limit_reverse_flux(flux, storage, dt)
   ! Unified source-cell scaling for reverse flow.
   !
   ! Problem: multiple upstream cells may reverse-drain the same downstream source.
   ! Each edge's |flux| is the demand; the source cell's storage is the supply.
   !
   ! Algorithm:
   !   1. Extract per-edge reverse demand: rev_demand(i) = max(-flux(i), 0) * dt
   !   2. push_ups2ucat: total_demand(j) = sum of rev_demand from all upstream edges
   !   3. Compute rate(j) = min(storage(j) / total_demand(j), 1) at each source cell
   !   4. push_next2ucat: distribute rate(j) back to each upstream cell as rate_edge(i)
   !   5. Scale: flux(i) *= rate_edge(i) for reverse edges
   !-------------------------------------------------------------------------------------
   USE MOD_Grid_RiverLakeNetwork, only: numucat, push_next2ucat, push_ups2ucat
   USE MOD_WorkerPushData
   IMPLICIT NONE

   real(r8), intent(inout) :: flux(:,:)      ! (nsed, numucat)
   real(r8), intent(in)    :: storage(:,:)   ! (nsed, numucat)
   real(r8), intent(in)    :: dt

   real(r8), allocatable :: rev_demand(:)    ! per-edge reverse demand for one grain class
   real(r8), allocatable :: total_demand(:)  ! total demand at each source cell
   real(r8), allocatable :: rate_src(:)      ! scale factor at source cell
   real(r8), allocatable :: rate_edge(:)     ! scale factor distributed to edges
   integer :: ised, i

      allocate(rev_demand  (numucat))
      allocate(total_demand(numucat))
      allocate(rate_src    (numucat))
      allocate(rate_edge   (numucat))

      DO ised = 1, nsed

         ! Step 1: extract reverse demand per edge
         DO i = 1, numucat
            rev_demand(i) = max(-flux(ised,i), 0._r8) * dt
         ENDDO

         ! Step 2: gather total demand at each source cell (downstream cell)
         CALL worker_push_data(push_ups2ucat, rev_demand, total_demand, &
            fillvalue = 0._r8, mode = 'sum')

         ! Step 3: compute scale factor at each source cell
         DO i = 1, numucat
            IF (total_demand(i) > 1.e-20_r8) THEN
               rate_src(i) = min(storage(ised,i) / total_demand(i), 1._r8)
            ELSE
               rate_src(i) = 1._r8
            ENDIF
         ENDDO

         ! Step 4: distribute scale factor back to upstream edges
         CALL worker_push_data(push_next2ucat, rate_src, rate_edge, fillvalue = 1._r8)

         ! Step 5: apply scale to reverse edges
         DO i = 1, numucat
            IF (flux(ised,i) < 0._r8) THEN
               flux(ised,i) = flux(ised,i) * rate_edge(i)
            ENDIF
         ENDDO

      ENDDO

      deallocate(rev_demand, total_demand, rate_src, rate_edge)

   END SUBROUTINE limit_reverse_flux

   !-------------------------------------------------------------------------------------
   SUBROUTINE apply_sediment_input(dt, rivsto, bed_area)
   ! Apply hillslope erosion input after exchange, following CoLM-sed-master more closely.
   ! Add input to suspended storage when enough water is present, then apply a single
   ! sed_max_conc cap. For shallow/dry cells, deposit directly into the bed layer.
   !-------------------------------------------------------------------------------------
   USE MOD_Grid_RiverLakeNetwork, only: numucat
   IMPLICIT NONE

   real(r8), intent(in) :: dt
   real(r8), intent(in) :: rivsto(:), bed_area(:)

      real(r8) :: sedsto_sum, dTmp(nsed), mass_before
      integer :: i

      IF (.not. p_is_worker) RETURN
      IF (numucat <= 0) RETURN

      DO i = 1, numucat
         IF (sum(sedinp(:,i)) <= 0._r8) CYCLE
         mass_before = sum(sedsto(:,i)) + (1._r8 - lambda) * sum(layer(:,i))

         IF (rivsto(i) >= bed_area(i) * sed_ignore_dph) THEN
            sedsto(:,i) = sedsto(:,i) + sedinp(:,i) * dt
            sedsto_sum = sum(sedsto(:,i))
            IF (sedsto_sum > rivsto(i) * sed_max_conc) THEN
               dTmp(:) = (sedsto_sum - rivsto(i) * sed_max_conc) &
                  * sedsto(:,i) / max(sedsto_sum, 1.e-20_r8)
               dTmp(:) = min(dTmp(:), sedsto(:,i))
               netflw(:,i) = netflw(:,i) - dTmp(:) / dt
               exch_d_eff(:,i) = exch_d_eff(:,i) + dTmp(:) / dt
               sedsto(:,i) = sedsto(:,i) - dTmp(:)
               layer(:,i) = layer(:,i) + dTmp(:) / (1._r8 - lambda)
            ENDIF
            sedcon(:,i) = sedsto(:,i) / rivsto(i)
         ELSE
            ! Shallow/dry cell: deposit erosion input directly into the bed layer.
            ! Record as negative netflw (net deposition) so that diagnostics capture
            ! this pathway in the mass balance.
            layer(:,i) = layer(:,i) + sedinp(:,i) * dt / (1._r8 - lambda)
            netflw(:,i) = netflw(:,i) - sedinp(:,i)
            exch_d_eff(:,i) = exch_d_eff(:,i) + sedinp(:,i)
         ENDIF
         CALL assert_sediment_mass_balance('hillslope input', i, mass_before, &
            sum(sedsto(:,i)) + (1._r8 - lambda) * sum(layer(:,i)), &
            sum(sedinp(:,i)) * dt, &
            rivsto_i = rivsto(i))
      ENDDO

   END SUBROUTINE apply_sediment_input

   !-------------------------------------------------------------------------------------
   SUBROUTINE calc_sediment_exchange(dt, rivsto, bed_area)
   !-------------------------------------------------------------------------------------
   USE MOD_Grid_RiverLakeNetwork, only: numucat
   IMPLICIT NONE

   real(r8), intent(in) :: dt
   real(r8), intent(in) :: rivsto(:), bed_area(:)

   real(r8) :: Es(nsed), D(nsed), Zd(nsed)
   real(r8) :: dTmp(nsed)
   real(r8) :: layer_sum, area, dTmp1, sedsto_sum, shear_eff, d_raw, mass_before
   real(r8) :: rouse_factor, profile_factor, transition_weight
   integer  :: i, ised
      IF (.not. p_is_worker) RETURN
      IF (numucat <= 0) RETURN

      exch_es_raw(:,:) = 0._r8
      exch_d_raw(:,:) = 0._r8
      exch_es_eff(:,:) = 0._r8
      exch_d_eff(:,:) = 0._r8

      DO i = 1, numucat
         IF (rivsto(i) < bed_area(i) * sed_ignore_dph) THEN
            netflw(:,i) = 0._r8
            CYCLE
         ENDIF
         mass_before = sum(sedsto(:,i)) + (1._r8 - lambda) * sum(layer(:,i))

         layer_sum = sum(layer(:,i))
         area = bed_area(i)

         IF (layer_sum <= 0._r8 .or. all(susvel(:,i) <= 0._r8)) THEN
            Es(:) = 0._r8
         ELSE
            Es(:) = susvel(:,i) * (1._r8 - lambda) * area * layer(:,i) / layer_sum
            Es(:) = max(Es(:), 0._r8)
         ENDIF

         IF (all(setvel(:) <= 0._r8)) THEN
            D(:) = 0._r8
         ELSE
            ! STATIC-WATER SETTLING uses depth-mean concentration: particles
            ! settle even when shear is zero. Above that limit, the Rouse
            ! profile converts depth-mean to near-bed concentration. A
            ! smoothstep blend avoids a discontinuity where the closures meet.
            shear_eff = max(shearvel(i), EXCH_SHEARVEL_MIN)
            DO ised = 1, nsed
               IF (shearvel(i) <= EXCH_SHEARVEL_MIN) THEN
                  rouse_factor = 1._r8
               ELSE
                  Zd(ised) = min(6._r8 * setvel(ised) / vonKar / shear_eff, EXCH_ZD_MAX)
                  IF (abs(Zd(ised)) < 1.0e-8_r8) THEN
                     profile_factor = 1._r8
                  ELSE
                     profile_factor = Zd(ised) / (1._r8 - exp(-Zd(ised)))
                  ENDIF
                  transition_weight = min(max((shearvel(i) - EXCH_SHEARVEL_MIN) / &
                     (EXCH_SHEARVEL_BLEND - EXCH_SHEARVEL_MIN), 0._r8), 1._r8)
                  transition_weight = transition_weight * transition_weight * &
                     (3._r8 - 2._r8 * transition_weight)
                  rouse_factor = 1._r8 + transition_weight * (profile_factor - 1._r8)
               ENDIF
               d_raw = setvel(ised) * area * sedcon(ised,i) * rouse_factor
               D(ised) = min(max(d_raw, 0._r8), sedsto(ised,i) / dt)
            ENDDO
         ENDIF

         exch_es_raw(:,i) = Es(:)
         exch_d_raw(:,i) = D(:)
         netflw(:,i) = Es(:) - D(:)

         DO ised = 1, nsed
            IF (abs(netflw(ised,i)) < 1.e-20_r8) THEN
               CYCLE
            ELSEIF (netflw(ised,i) > 0._r8) THEN
               dTmp1 = netflw(ised,i) * dt / (1._r8 - lambda)
               IF (dTmp1 < layer(ised,i)) THEN
                  layer(ised,i) = layer(ised,i) - dTmp1
               ELSE
                  netflw(ised,i) = layer(ised,i) * (1._r8 - lambda) / dt
                  layer(ised,i) = 0._r8
               ENDIF
               sedsto(ised,i) = sedsto(ised,i) + netflw(ised,i) * dt
            ELSE
               IF (abs(netflw(ised,i)) * dt < sedsto(ised,i)) THEN
                  sedsto(ised,i) = max(sedsto(ised,i) - abs(netflw(ised,i)) * dt, 0._r8)
               ELSE
                  netflw(ised,i) = -sedsto(ised,i) / dt
                  sedsto(ised,i) = 0._r8
               ENDIF
               layer(ised,i) = layer(ised,i) + abs(netflw(ised,i)) * dt / (1._r8 - lambda)
            ENDIF
         ENDDO

         DO ised = 1, nsed
            IF (Es(ised) >= D(ised)) THEN
               exch_d_eff(ised,i) = D(ised)
               exch_es_eff(ised,i) = D(ised) + max(netflw(ised,i), 0._r8)
            ELSE
               exch_es_eff(ised,i) = Es(ised)
               exch_d_eff(ised,i) = Es(ised) + max(-netflw(ised,i), 0._r8)
            ENDIF
         ENDDO

         ! Enforce concentration cap after exchange (matches CaMa's unconditional cap).
         ! Without this, strong entrainment (Es >> D) could push concentration above
         ! sed_max_conc indefinitely when there is no erosion input to trigger the
         ! cap in apply_sediment_input.
         sedsto_sum = sum(sedsto(:,i))
         IF (rivsto(i) > 0._r8 .and. sedsto_sum > rivsto(i) * sed_max_conc) THEN
            dTmp(:) = (sedsto_sum - rivsto(i) * sed_max_conc) &
               * sedsto(:,i) / max(sedsto_sum, 1.e-20_r8)
            dTmp(:) = min(dTmp(:), sedsto(:,i))
            netflw(:,i) = netflw(:,i) - dTmp(:) / dt
            exch_d_eff(:,i) = exch_d_eff(:,i) + dTmp(:) / dt
            sedsto(:,i) = sedsto(:,i) - dTmp(:)
            layer(:,i) = layer(:,i) + dTmp(:) / (1._r8 - lambda)
         ENDIF

         IF (rivsto(i) > 0._r8) THEN
            sedcon(:,i) = sedsto(:,i) / rivsto(i)
         ENDIF
         CALL assert_sediment_mass_balance('exchange', i, mass_before, &
            sum(sedsto(:,i)) + (1._r8 - lambda) * sum(layer(:,i)), 0._r8, &
            rivsto_i = rivsto(i))
      ENDDO

   END SUBROUTINE calc_sediment_exchange

   !-------------------------------------------------------------------------------------
   SUBROUTINE calc_layer_redistribution(bed_area)
   ! Bug fix: seddepP now uses (nsed, totlyrnum+1) to match seddep layout (nsed, totlyrnum, numucat)
   !-------------------------------------------------------------------------------------
   USE MOD_Grid_RiverLakeNetwork, only: numucat
   IMPLICIT NONE

   real(r8), intent(in) :: bed_area(:)

   real(r8) :: lyrvol, diff, mass_before
   real(r8) :: layerP(nsed), seddepP(nsed, totlyrnum+1), tmp(nsed)
   real(r8) :: tmpsum
   integer  :: i, ilyr, jlyr
   integer  :: slyr

      IF (.not. p_is_worker) RETURN
      IF (numucat <= 0) RETURN

      DO i = 1, numucat
         lyrvol = lyrdph * bed_area(i)
         mass_before = (1._r8 - lambda) * &
            (sum(layer(:,i)) + sum(seddep(:,:,i)))
         slyr = 0

         layer(:,i) = max(layer(:,i), 0._r8)
         seddep(:,:,i) = max(seddep(:,:,i), 0._r8)

         IF (sum(layer(:,i)) + sum(seddep(:,:,i)) <= lyrvol) THEN
            layer(:,i) = layer(:,i) + sum(seddep(:,:,i), dim=2)
            seddep(:,:,i) = 0._r8
            CALL assert_sediment_mass_balance('layer redistribution', i, mass_before, &
               (1._r8 - lambda) * (sum(layer(:,i)) + sum(seddep(:,:,i))), 0._r8)
            CYCLE
         ENDIF

         layerP(:) = layer(:,i)
         IF (sum(layerP(:)) >= lyrvol) THEN
            layer(:,i) = layerP(:) * min(lyrvol / max(sum(layerP(:)), 1.e-20_r8), 1._r8)
            layerP(:) = max(layerP(:) - layer(:,i), 0._r8)
            slyr = 0
         ELSEIF (sum(seddep(:,:,i)) > 0._r8) THEN
            layerP(:) = 0._r8
            DO ilyr = 1, totlyrnum
               diff = lyrvol - sum(layer(:,i))
               IF (diff <= 0._r8) EXIT
               tmpsum = sum(seddep(:,ilyr,i))
               IF (tmpsum <= diff) THEN
                  layer(:,i) = layer(:,i) + seddep(:,ilyr,i)
                  seddep(:,ilyr,i) = 0._r8
                  slyr = ilyr + 1
               ELSE
                  IF (tmpsum > 1.e-20_r8) THEN
                     tmp(:) = diff * seddep(:,ilyr,i) / tmpsum
                  ELSE
                     tmp(:) = 0._r8
                  ENDIF
                  layer(:,i) = layer(:,i) + tmp(:)
                  seddep(:,ilyr,i) = max(seddep(:,ilyr,i) - tmp(:), 0._r8)
                  slyr = ilyr
                  EXIT
               ENDIF
            ENDDO
         ELSE
            seddep(:,:,i) = 0._r8
            CALL assert_sediment_mass_balance('layer redistribution', i, mass_before, &
               (1._r8 - lambda) * (sum(layer(:,i)) + sum(seddep(:,:,i))), 0._r8)
            CYCLE
         ENDIF

         ! If the active layer was compressed, layerP stores excess material that
         ! still needs to be pushed into the bed even when the existing bed is empty.
         IF (sum(seddep(:,:,i)) <= 0._r8 .and. sum(layerP(:)) <= 0._r8) THEN
            CALL assert_sediment_mass_balance('layer redistribution', i, mass_before, &
               (1._r8 - lambda) * (sum(layer(:,i)) + sum(seddep(:,:,i))), 0._r8)
            CYCLE
         ENDIF

         ! seddepP: (nsed, totlyrnum+1) -- slot 1 = excess from layer, slots 2: = bed layers
         seddepP(:,1) = layerP(:)
         seddepP(:,2:totlyrnum+1) = seddep(:,1:totlyrnum,i)
         seddep(:,:,i) = 0._r8

         DO ilyr = 1, totlyrnum - 1
            IF (sum(seddep(:,ilyr,i)) >= lyrvol) CYCLE
            DO jlyr = slyr + 1, totlyrnum + 1
               diff = lyrvol - sum(seddep(:,ilyr,i))
               IF (diff <= 0._r8) EXIT
               tmpsum = sum(seddepP(:,jlyr))
               IF (tmpsum <= diff) THEN
                  seddep(:,ilyr,i) = seddep(:,ilyr,i) + seddepP(:,jlyr)
                  seddepP(:,jlyr) = 0._r8
               ELSE
                  IF (tmpsum > 1.e-20_r8) THEN
                     tmp(:) = diff * seddepP(:,jlyr) / tmpsum
                  ELSE
                     tmp(:) = 0._r8
                  ENDIF
                  seddep(:,ilyr,i) = seddep(:,ilyr,i) + tmp(:)
                  seddepP(:,jlyr) = max(seddepP(:,jlyr) - tmp(:), 0._r8)
                  EXIT
               ENDIF
            ENDDO
         ENDDO

         IF (sum(seddepP) > 0._r8) THEN
            seddep(:,totlyrnum,i) = seddep(:,totlyrnum,i) + sum(seddepP, dim=2)
         ENDIF
         CALL assert_sediment_mass_balance('layer redistribution', i, mass_before, &
            (1._r8 - lambda) * (sum(layer(:,i)) + sum(seddep(:,:,i))), 0._r8)
      ENDDO

   END SUBROUTINE calc_layer_redistribution


   !-------------------------------------------------------------------------------------
   SUBROUTINE assert_sediment_mass_balance(context, index, before, after, expected_change, &
                                        rivsto_i, rivout_i)
   USE MOD_Grid_RiverLakeNetwork, only: ucat_next
   IMPLICIT NONE

   character(len=*), intent(in) :: context
   integer, intent(in) :: index
   real(r8), intent(in) :: before, after, expected_change
   real(r8), intent(in), optional :: rivsto_i, rivout_i

   real(r8) :: residual, scale, tolerance
   integer :: ised, ilyr

   residual = after - before - expected_change
   scale = max(abs(before), abs(after), abs(expected_change))
   tolerance = SED_BALANCE_ABS_TOL + SED_BALANCE_REL_TOL * scale

   IF (.not. ieee_is_finite(residual) .or. abs(residual) > tolerance) THEN

      ! IMPORTANT:
      ! Do NOT use "IF (p_is_io)" here.
      ! Sediment calculation runs on worker ranks, so the failing worker
      ! must print its own diagnostics.

      WRITE(*,'(A)') ' '
      WRITE(*,'(A)') '============================================================'
      WRITE(*,'(A)') 'SEDIMENT MASS BALANCE FAILURE'
      WRITE(*,'(A)') '============================================================'

      WRITE(*,'(A,A)')        'context          = ', trim(context)
      WRITE(*,'(A,I0)')       'worker           = ', p_iam_worker
      WRITE(*,'(A,I0)')       'local cell i     = ', index
      WRITE(*,'(A,I0)')       'ucat_next        = ', ucat_next(index)

      WRITE(*,'(A,ES24.14)')  'mass before      = ', before
      WRITE(*,'(A,ES24.14)')  'mass after       = ', after
      WRITE(*,'(A,ES24.14)')  'expected change  = ', expected_change
      WRITE(*,'(A,ES24.14)')  'residual         = ', residual
      WRITE(*,'(A,ES24.14)')  'tolerance        = ', tolerance
      WRITE(*,'(A,ES24.14)')  'scale            = ', scale

      IF (present(rivsto_i)) THEN
         WRITE(*,'(A,ES24.14)') 'rivsto [m3]      = ', rivsto_i
      ELSE
         WRITE(*,'(A)')         'rivsto [m3]      = not supplied'
      ENDIF

      IF (present(rivout_i)) THEN
         WRITE(*,'(A,ES24.14)') 'rivout [m3/s]    = ', rivout_i
      ELSE
         WRITE(*,'(A)')         'rivout [m3/s]    = not supplied'
      ENDIF

      WRITE(*,'(A)') ' '
      WRITE(*,'(A)') '--- sediment classes ---'

      DO ised = 1, nsed
         WRITE(*,'(A,I0)') 'sediment class = ', ised

         WRITE(*,'(A,ES24.14)') '  sedsto           = ', sedsto(ised,index)
         WRITE(*,'(A,ES24.14)') '  sedcon           = ', sedcon(ised,index)
         WRITE(*,'(A,ES24.14)') '  layer            = ', layer(ised,index)
         WRITE(*,'(A,ES24.14)') '  sedout           = ', sedout(ised,index)
         WRITE(*,'(A,ES24.14)') '  bedout           = ', bedout(ised,index)
         WRITE(*,'(A,ES24.14)') '  sedinp           = ', sedinp(ised,index)
         WRITE(*,'(A,ES24.14)') '  netflw           = ', netflw(ised,index)
         WRITE(*,'(A,ES24.14)') '  exch_es_raw      = ', exch_es_raw(ised,index)
         WRITE(*,'(A,ES24.14)') '  exch_d_raw       = ', exch_d_raw(ised,index)
         WRITE(*,'(A,ES24.14)') '  exch_es_eff      = ', exch_es_eff(ised,index)
         WRITE(*,'(A,ES24.14)') '  exch_d_eff       = ', exch_d_eff(ised,index)
         WRITE(*,'(A,ES24.14)') '  netflw_adv_step  = ', netflw_adv_step(ised,index)
         WRITE(*,'(A,ES24.14)') '  exch_d_adv_step  = ', exch_d_adv_step(ised,index)

         DO ilyr = 1, totlyrnum
            WRITE(*,'(A,I0,A,ES24.14)') &
               '  seddep layer ', ilyr, ' = ', seddep(ised,ilyr,index)
         ENDDO
      ENDDO

      WRITE(*,'(A)') ' '
      WRITE(*,'(A,ES24.14)') 'sum sedsto       = ', sum(sedsto(:,index))
      WRITE(*,'(A,ES24.14)') 'sum layer        = ', sum(layer(:,index))
      WRITE(*,'(A,ES24.14)') 'sum seddep       = ', sum(seddep(:,:,index))
      WRITE(*,'(A,ES24.14)') 'bed solid layer  = ', &
         (1._r8-lambda) * sum(layer(:,index))
      WRITE(*,'(A,ES24.14)') 'bed solid total  = ', &
         (1._r8-lambda) * (sum(layer(:,index)) + sum(seddep(:,:,index)))

      WRITE(*,'(A)') '============================================================'
      WRITE(*,'(A)') ' '

      FLUSH(6)

      CALL CoLM_stop('sediment mass balance failure')
   ENDIF

   END SUBROUTINE assert_sediment_mass_balance
   !-------------------------------------------------------------------------------------

   !-------------------------------------------------------------------------------------
   SUBROUTINE calc_sediment_yield(fldfrc, grarea, prcp_time)
   !-------------------------------------------------------------------------------------
   ! Compute hillslope erosion input using pre-accumulated yield power-law term.
   ! sed_precip_yield stores sum[ (rate_mm_hr)^pyldpc * dt ] over forcing steps.
   ! Dividing by prcp_time gives the time-averaged <(rate_mm_hr)^pyldpc>, which
   ! correctly preserves the high-intensity contribution (no Jensen bias).
   ! The precipitation threshold is still evaluated from the routing-period mean
   ! rain rate, preserving the original trigger definition.
   !-------------------------------------------------------------------------------------
   USE MOD_Grid_RiverLakeNetwork, only: numucat
   IMPLICIT NONE

   real(r8), intent(in) :: fldfrc(:), grarea(:)
   real(r8), intent(in) :: prcp_time    ! Precipitation accumulation time [s]

   real(r8) :: precip_yield_avg, precip_rate_avg, precip_mm_day_avg
   integer  :: i, ilyr

      IF (.not. p_is_worker) RETURN
      IF (numucat <= 0) RETURN

      sedinp(:,:) = 0._r8

      IF (prcp_time <= 0._r8) RETURN

      DO i = 1, numucat
         precip_rate_avg = sed_precip(i) / prcp_time
         precip_mm_day_avg = precip_rate_avg * 86400._r8
         IF (precip_mm_day_avg <= SED_PRECIP_THRESHOLD_MM_DAY) CYCLE
         IF (sed_precip_yield(i) <= 0._r8) CYCLE

         ! Time-averaged yield power-law term: <(rate_mm_hr)^pyldpc>
         precip_yield_avg = sed_precip_yield(i) / prcp_time

         DO ilyr = 1, nlfp_sed
            IF (fldfrc(i) * nlfp_sed > real(ilyr, r8)) CYCLE

            sedinp(:,i) = sedinp(:,i) + &
               pyld * precip_yield_avg * sed_slope(ilyr,i)**pyldc / 3600._r8 &
               * grarea(i) * min(real(ilyr, r8)/real(nlfp_sed, r8) - fldfrc(i), 1._r8/real(nlfp_sed, r8)) &
               * dsylunit * sed_frc(:,i)
         ENDDO
      ENDDO

   END SUBROUTINE calc_sediment_yield

   !-------------------------------------------------------------------------------------
   SUBROUTINE read_sediment_restart(file_restart)
   ! Read sediment state from restart file using a temp buffer (vector_read_and_scatter
   ! requires allocatable argument, so we cannot pass array slices directly).
   !-------------------------------------------------------------------------------------
   USE netcdf
   USE MOD_NetCDFSerial
   USE MOD_Vector_ReadWrite
   USE MOD_Grid_RiverLakeNetwork, only: numucat, totalnumucat, ucat_data_address, &
      topo_rivwth, topo_rivlen
   IMPLICIT NONE

   character(len=*), intent(in) :: file_restart
      real(r8), allocatable :: buf(:)
      integer :: ised, ilyr, ncid, varid, ierr
   character(len=16) :: cised, cilyr
   character(len=64) :: vname
      logical :: file_ok, has_schema, has_complete, has_canonical_mass
      logical :: has_any_sediment, legacy_nonzero
      logical :: legacy_state_nonzero
      integer :: nread
      logical :: meta_bad

      ! All processes must participate (MPI collective calls inside).
      ! Do NOT return early on non-workers before MPI calls.
      IF (.not. sediment_particle_enabled()) RETURN

      ! Detect any recognizable sediment field, not just sedcon_1.  A partial
      ! transaction must fail strict reads rather than masquerade as a cold start.
      file_ok = .false.
      has_schema = .false.
      has_complete = .false.
      has_canonical_mass = .false.
      has_any_sediment = .false.
      IF (p_is_master) THEN
         inquire(file=trim(file_restart), exist=file_ok)
         IF (file_ok) THEN
            ierr = nf90_open(trim(file_restart), NF90_NOWRITE, ncid)
            IF (ierr == NF90_NOERR) THEN
               has_schema = (nf90_inq_varid(ncid, 'sed_restart_schema_meta', varid) == NF90_NOERR)
               has_complete = (nf90_inq_varid(ncid, 'sed_restart_complete_meta', varid) == NF90_NOERR)
               has_canonical_mass = (nf90_inq_varid(ncid, 'sedsto_1', varid) == NF90_NOERR)
               has_any_sediment = has_schema .or. has_complete .or. has_canonical_mass
               has_any_sediment = has_any_sediment .or. &
                  (nf90_inq_varid(ncid, 'sed_n_meta', varid) == NF90_NOERR)
               has_any_sediment = has_any_sediment .or. &
                  (nf90_inq_varid(ncid, 'sedcon_1', varid) == NF90_NOERR)
               has_any_sediment = has_any_sediment .or. &
                  (nf90_inq_varid(ncid, 'layer_1', varid) == NF90_NOERR)
               has_any_sediment = has_any_sediment .or. &
                  (nf90_inq_varid(ncid, 'seddep_1_1', varid) == NF90_NOERR)
               has_any_sediment = has_any_sediment .or. &
                  (nf90_inq_varid(ncid, 'sed_acc_time', varid) == NF90_NOERR)
               has_any_sediment = has_any_sediment .or. &
                  (nf90_inq_varid(ncid, 'a_sedcon_1', varid) == NF90_NOERR)
               file_ok = has_any_sediment
               IF (.not. has_any_sediment) ierr = nf90_close(ncid)
            ELSE
               file_ok = .false.
            ENDIF
         ENDIF
      ENDIF
#ifdef USEMPI
      CALL mpi_bcast(file_ok, 1, MPI_LOGICAL, p_address_master, p_comm_glb, p_err)
      CALL mpi_bcast(has_schema, 1, MPI_LOGICAL, p_address_master, p_comm_glb, p_err)
      CALL mpi_bcast(has_complete, 1, MPI_LOGICAL, p_address_master, p_comm_glb, p_err)
      CALL mpi_bcast(has_canonical_mass, 1, MPI_LOGICAL, p_address_master, p_comm_glb, p_err)
#endif

      IF (.not. file_ok) THEN
         IF (p_is_io) WRITE(*,*) 'Sediment restart variables not found, using initial state.'
         RETURN
      ENDIF

      IF (.not. has_schema .and. (has_complete .or. has_canonical_mass)) THEN
         IF (p_is_io) WRITE(*,'(A)') &
            'ERROR: partial sediment transaction has canonical fields but no schema marker.'
         IF (p_is_master) ierr = nf90_close(ncid)
         CALL CoLM_stop()
      ENDIF

      nread = 0
      meta_bad = .false.

      IF (has_schema) THEN
         CALL check_sediment_restart_scalar(file_restart, ncid, &
            'sed_restart_schema_meta', real(SED_RESTART_SCHEMA_VERSION, r8), &
            numucat, ucat_data_address, meta_bad)
         CALL check_sediment_restart_scalar(file_restart, ncid, &
            'sed_restart_complete_meta', 1._r8, &
            numucat, ucat_data_address, meta_bad)
         CALL check_sediment_restart_scalar(file_restart, ncid, 'sed_nlfp_meta', &
            real(nlfp_sed, r8), numucat, ucat_data_address, meta_bad)
         CALL check_sediment_restart_scalar(file_restart, ncid, 'sed_lambda_meta', &
            lambda, numucat, ucat_data_address, meta_bad)
         CALL check_sediment_restart_scalar(file_restart, ncid, 'sed_lyrdph_meta', &
            lyrdph, numucat, ucat_data_address, meta_bad)
         CALL check_sediment_restart_scalar(file_restart, ncid, 'sed_psedd_meta', &
            psedD, numucat, ucat_data_address, meta_bad)
         CALL check_sediment_restart_scalar(file_restart, ncid, 'sed_pwatd_meta', &
            pwatD, numucat, ucat_data_address, meta_bad)
         CALL check_sediment_restart_scalar(file_restart, ncid, 'sed_viskin_meta', &
            visKin, numucat, ucat_data_address, meta_bad)
         CALL check_sediment_restart_scalar(file_restart, ncid, 'sed_vonkar_meta', &
            vonKar, numucat, ucat_data_address, meta_bad)
         CALL check_sediment_restart_scalar(file_restart, ncid, 'sed_pset_meta', &
            pset, numucat, ucat_data_address, meta_bad)
         CALL check_sediment_restart_scalar(file_restart, ncid, 'sed_ignore_dph_meta', &
            sed_ignore_dph, numucat, ucat_data_address, meta_bad)
         CALL check_sediment_restart_scalar(file_restart, ncid, 'sed_cfl_adv_meta', &
            sed_cfl_adv, numucat, ucat_data_address, meta_bad)
         CALL check_sediment_restart_scalar(file_restart, ncid, 'sed_dt_max_meta', &
            sed_dt_max, numucat, ucat_data_address, meta_bad)
         CALL check_sediment_restart_scalar(file_restart, ncid, 'sed_bed_depth_meta', &
            sed_bed_depth, numucat, ucat_data_address, meta_bad)
         CALL check_sediment_restart_scalar(file_restart, ncid, 'sed_pyld_meta', &
            pyld, numucat, ucat_data_address, meta_bad)
         CALL check_sediment_restart_scalar(file_restart, ncid, 'sed_pyldc_meta', &
            pyldc, numucat, ucat_data_address, meta_bad)
         CALL check_sediment_restart_scalar(file_restart, ncid, 'sed_pyldpc_meta', &
            pyldpc, numucat, ucat_data_address, meta_bad)
         CALL check_sediment_restart_scalar(file_restart, ncid, 'sed_dsylunit_meta', &
            dsylunit, numucat, ucat_data_address, meta_bad)
         CALL check_sediment_restart_scalar(file_restart, ncid, 'sed_max_conc_meta', &
            sed_max_conc, numucat, ucat_data_address, meta_bad)
         CALL check_sediment_restart_scalar(file_restart, ncid, 'sed_bedload_coeff_meta', &
            SED_BEDLOAD_COEFF, numucat, ucat_data_address, meta_bad)
         CALL check_sediment_restart_scalar(file_restart, ncid, 'sed_precip_threshold_meta', &
            SED_PRECIP_THRESHOLD_MM_DAY, numucat, ucat_data_address, meta_bad)
         CALL check_sediment_restart_scalar(file_restart, ncid, 'sed_exch_shear_min_meta', &
            EXCH_SHEARVEL_MIN, numucat, ucat_data_address, meta_bad)
         CALL check_sediment_restart_scalar(file_restart, ncid, 'sed_exch_shear_blend_meta', &
            EXCH_SHEARVEL_BLEND, numucat, ucat_data_address, meta_bad)
         CALL check_sediment_restart_scalar(file_restart, ncid, 'sed_exch_zd_max_meta', &
            EXCH_ZD_MAX, numucat, ucat_data_address, meta_bad)
         DO ised = 1, nsed
            WRITE(cised, '(I0)') ised
            CALL check_sediment_restart_scalar(file_restart, ncid, &
               'sed_diam_meta_' // trim(cised), sDiam(ised), numucat, &
               ucat_data_address, meta_bad)
            CALL check_sediment_restart_scalar(file_restart, ncid, &
               'sed_setvel_meta_' // trim(cised), setvel(ised), numucat, &
               ucat_data_address, meta_bad)
         ENDDO
      ENDIF

      CALL require_sediment_restart_var(file_restart, ncid, 'sed_n_meta', buf, numucat, ucat_data_address)
      IF (allocated(buf)) THEN
         IF (p_is_worker .and. numucat > 0) meta_bad = meta_bad .or. any(nint(buf(:)) /= nsed)
         deallocate(buf)
      ENDIF
      CALL require_sediment_restart_var(file_restart, ncid, 'sed_totlyrnum_meta', buf, numucat, ucat_data_address)
      IF (allocated(buf)) THEN
         IF (p_is_worker .and. numucat > 0) meta_bad = meta_bad .or. any(nint(buf(:)) /= totlyrnum)
         deallocate(buf)
      ENDIF
      IF (has_schema) THEN
         DO ised = 1, nsed
            WRITE(cised, '(I0)') ised
            vname = 'sed_frc_meta_' // trim(cised)
            CALL require_sediment_restart_var(file_restart, ncid, vname, buf, &
               numucat, ucat_data_address)
            IF (allocated(buf)) THEN
               IF (p_is_worker .and. numucat > 0) meta_bad = meta_bad .or. &
                  any(.not. sediment_restart_real_matches(buf, sed_frc(ised,:)))
               deallocate(buf)
            ENDIF
         ENDDO
         DO ilyr = 1, nlfp_sed
            WRITE(cilyr, '(I0)') ilyr
            vname = 'sed_slope_meta_' // trim(cilyr)
            CALL require_sediment_restart_var(file_restart, ncid, vname, buf, &
               numucat, ucat_data_address)
            IF (allocated(buf)) THEN
               IF (p_is_worker .and. numucat > 0) meta_bad = meta_bad .or. &
                  any(.not. sediment_restart_real_matches(buf, sed_slope(ilyr,:)))
               deallocate(buf)
            ENDIF
         ENDDO
         CALL require_sediment_restart_var(file_restart, ncid, 'sed_rivwth_meta', &
            buf, numucat, ucat_data_address)
         IF (allocated(buf)) THEN
            IF (p_is_worker .and. numucat > 0) meta_bad = meta_bad .or. &
               any(.not. sediment_restart_real_matches(buf, topo_rivwth))
            deallocate(buf)
         ENDIF
         CALL require_sediment_restart_var(file_restart, ncid, 'sed_rivlen_meta', &
            buf, numucat, ucat_data_address)
         IF (allocated(buf)) THEN
            IF (p_is_worker .and. numucat > 0) meta_bad = meta_bad .or. &
               any(.not. sediment_restart_real_matches(buf, topo_rivlen))
            deallocate(buf)
         ENDIF
      ENDIF
#ifdef USEMPI
      CALL mpi_allreduce(MPI_IN_PLACE, meta_bad, 1, MPI_LOGICAL, MPI_LOR, p_comm_glb, p_err)
#endif
      IF (meta_bad) THEN
         IF (p_is_io) WRITE(*,'(A)') &
            'ERROR: sediment restart schema/configuration does not match the active sediment model.'
         IF (p_is_master) ierr = nf90_close(ncid)
         CALL CoLM_stop()
      ENDIF

      DO ised = 1, nsed
         WRITE(cised, '(I0)') ised
         vname = 'sedcon_' // trim(cised)
         CALL require_sediment_restart_var(file_restart, ncid, vname, buf, numucat, ucat_data_address)
         IF (allocated(sedcon) .and. p_is_worker .and. numucat > 0) THEN
            sedcon(ised,:) = buf(:)
            nread = nread + 1
         ENDIF
         IF (allocated(buf)) deallocate(buf)
      ENDDO

      IF (has_schema) THEN
         DO ised = 1, nsed
            WRITE(cised, '(I0)') ised
            vname = 'sedsto_' // trim(cised)
            CALL require_sediment_restart_var(file_restart, ncid, vname, buf, &
               numucat, ucat_data_address)
            IF (allocated(sedsto) .and. p_is_worker .and. numucat > 0) THEN
               sedsto(ised,:) = buf(:)
               nread = nread + 1
            ENDIF
            IF (allocated(buf)) deallocate(buf)
         ENDDO
      ELSE
         ! Legacy files stored concentration without its carrier water volume,
         ! so nonzero suspended mass cannot be recovered exactly.  Zero is the
         ! sole lossless migration and initializes the canonical mass to zero.
         legacy_nonzero = .false.
         IF (p_is_worker .and. allocated(sedcon)) THEN
            legacy_nonzero = any(.not. (sedcon == 0._r8))
         ENDIF
#ifdef USEMPI
         CALL mpi_allreduce(MPI_IN_PLACE, legacy_nonzero, 1, MPI_LOGICAL, &
            MPI_LOR, p_comm_glb, p_err)
#endif
         IF (legacy_nonzero) THEN
            IF (p_is_io) WRITE(*,'(A)') &
               'ERROR: legacy sediment restart has nonzero concentration but no exact suspended mass.'
            IF (p_is_master) ierr = nf90_close(ncid)
            CALL CoLM_stop()
         ENDIF
         IF (allocated(sedsto)) sedsto(:,:) = 0._r8
      ENDIF

      DO ised = 1, nsed
         WRITE(cised, '(I0)') ised
         vname = 'layer_' // trim(cised)
         CALL require_sediment_restart_var(file_restart, ncid, vname, buf, numucat, ucat_data_address)
         IF (allocated(layer) .and. p_is_worker .and. numucat > 0) THEN
            layer(ised,:) = buf(:)
            nread = nread + 1
         ENDIF
         IF (allocated(buf)) deallocate(buf)
      ENDDO

      DO ised = 1, nsed
         WRITE(cised, '(I0)') ised
         DO ilyr = 1, totlyrnum
            WRITE(cilyr, '(I0)') ilyr
            vname = 'seddep_' // trim(cised) // '_' // trim(cilyr)
            CALL require_sediment_restart_var(file_restart, ncid, vname, buf, numucat, ucat_data_address)
            IF (allocated(seddep) .and. p_is_worker .and. numucat > 0) THEN
               seddep(ised,ilyr,:) = buf(:)
               nread = nread + 1
            ENDIF
            IF (allocated(buf)) deallocate(buf)
         ENDDO
      ENDDO

      CALL require_sediment_restart_var(file_restart, ncid, 'sed_acc_time', buf, numucat, ucat_data_address)
      IF (allocated(sed_acc_time) .and. p_is_worker .and. numucat > 0) sed_acc_time(:) = buf(:)
      IF (allocated(buf)) deallocate(buf)

      CALL require_sediment_restart_var(file_restart, ncid, 'sed_acc_v2', buf, numucat, ucat_data_address)
      IF (allocated(sed_acc_v2) .and. p_is_worker .and. numucat > 0) sed_acc_v2(:) = buf(:)
      IF (allocated(buf)) deallocate(buf)

      CALL require_sediment_restart_var(file_restart, ncid, 'sed_acc_wdsrf', buf, numucat, ucat_data_address)
      IF (allocated(sed_acc_wdsrf) .and. p_is_worker .and. numucat > 0) sed_acc_wdsrf(:) = buf(:)
      IF (allocated(buf)) deallocate(buf)

      CALL require_sediment_restart_var(file_restart, ncid, 'sed_acc_rivsto', buf, numucat, ucat_data_address)
      IF (allocated(sed_acc_rivsto) .and. p_is_worker .and. numucat > 0) sed_acc_rivsto(:) = buf(:)
      IF (allocated(buf)) deallocate(buf)

      CALL require_sediment_restart_var(file_restart, ncid, 'sed_acc_rivout', buf, numucat, ucat_data_address)
      IF (allocated(sed_acc_rivout) .and. p_is_worker .and. numucat > 0) sed_acc_rivout(:) = buf(:)
      IF (allocated(buf)) deallocate(buf)

      CALL require_sediment_restart_var(file_restart, ncid, 'sed_acc_abs_rivout', buf, numucat, ucat_data_address)
      IF (allocated(sed_acc_abs_rivout) .and. p_is_worker .and. numucat > 0) sed_acc_abs_rivout(:) = buf(:)
      IF (allocated(buf)) deallocate(buf)

      CALL require_sediment_restart_var(file_restart, ncid, 'sed_acc_floodarea', buf, numucat, ucat_data_address)
      IF (allocated(sed_acc_floodarea) .and. p_is_worker .and. numucat > 0) sed_acc_floodarea(:) = buf(:)
      IF (allocated(buf)) deallocate(buf)

      CALL require_sediment_restart_var(file_restart, ncid, 'sed_precip', buf, numucat, ucat_data_address)
      IF (allocated(sed_precip) .and. p_is_worker .and. numucat > 0) sed_precip(:) = buf(:)
      IF (allocated(buf)) deallocate(buf)

      CALL require_sediment_restart_var(file_restart, ncid, 'sed_precip_yield', buf, numucat, ucat_data_address)
      IF (allocated(sed_precip_yield) .and. p_is_worker .and. numucat > 0) sed_precip_yield(:) = buf(:)
      IF (allocated(buf)) deallocate(buf)

      CALL require_sediment_restart_var(file_restart, ncid, 'sed_precip_time_vec', buf, numucat, ucat_data_address)
      IF (allocated(buf)) THEN
         IF (p_is_worker .and. numucat > 0) THEN
            sed_precip_time = buf(1)
         ELSE
            sed_precip_time = 0._r8
         ENDIF
      ENDIF
      IF (allocated(buf)) deallocate(buf)

      DO ised = 1, nsed
         WRITE(cised, '(I0)') ised
         vname = 'a_sedcon_' // trim(cised)
         CALL require_sediment_restart_var(file_restart, ncid, vname, buf, numucat, ucat_data_address)
         IF (allocated(a_sedcon) .and. p_is_worker .and. numucat > 0) a_sedcon(ised,:) = buf(:)
         IF (allocated(buf)) deallocate(buf)

         vname = 'a_sedout_' // trim(cised)
         CALL require_sediment_restart_var(file_restart, ncid, vname, buf, numucat, ucat_data_address)
         IF (allocated(a_sedout) .and. p_is_worker .and. numucat > 0) a_sedout(ised,:) = buf(:)
         IF (allocated(buf)) deallocate(buf)

         vname = 'a_bedout_' // trim(cised)
         CALL require_sediment_restart_var(file_restart, ncid, vname, buf, numucat, ucat_data_address)
         IF (allocated(a_bedout) .and. p_is_worker .and. numucat > 0) a_bedout(ised,:) = buf(:)
         IF (allocated(buf)) deallocate(buf)

         vname = 'a_sedinp_' // trim(cised)
         CALL require_sediment_restart_var(file_restart, ncid, vname, buf, numucat, ucat_data_address)
         IF (allocated(a_sedinp) .and. p_is_worker .and. numucat > 0) a_sedinp(ised,:) = buf(:)
         IF (allocated(buf)) deallocate(buf)

         vname = 'a_netflw_' // trim(cised)
         CALL require_sediment_restart_var(file_restart, ncid, vname, buf, numucat, ucat_data_address)
         IF (allocated(a_netflw) .and. p_is_worker .and. numucat > 0) a_netflw(ised,:) = buf(:)
         IF (allocated(buf)) deallocate(buf)

         vname = 'a_layer_' // trim(cised)
         CALL require_sediment_restart_var(file_restart, ncid, vname, buf, numucat, ucat_data_address)
         IF (allocated(a_layer) .and. p_is_worker .and. numucat > 0) a_layer(ised,:) = buf(:)
         IF (allocated(buf)) deallocate(buf)
      ENDDO

      CALL require_sediment_restart_var(file_restart, ncid, 'a_shearvel', buf, numucat, ucat_data_address)
      IF (allocated(a_shearvel) .and. p_is_worker .and. numucat > 0) a_shearvel(:) = buf(:)
      IF (allocated(buf)) deallocate(buf)

      CALL require_sediment_restart_var(file_restart, ncid, 'sed_hist_acctime_vec', buf, numucat, ucat_data_address)
      IF (allocated(buf)) THEN
         IF (p_is_worker .and. numucat > 0) THEN
            sed_hist_acctime = buf(1)
         ELSE
            sed_hist_acctime = 0._r8
         ENDIF
      ENDIF
      IF (allocated(buf)) deallocate(buf)

      IF (.not. has_schema) THEN
         ! Without a descriptor the old bed/queued/history state cannot be
         ! interpreted under possibly changed porosity, grain, or yield
         ! parameters.  Migration is lossless only when every such value is zero.
         legacy_state_nonzero = .false.
         IF (p_is_worker) THEN
            legacy_state_nonzero = any(.not. (sedcon == 0._r8)) .or. &
               any(.not. (layer == 0._r8)) .or. any(.not. (seddep == 0._r8)) .or. &
               any(.not. (sed_acc_time == 0._r8)) .or. &
               any(.not. (sed_acc_v2 == 0._r8)) .or. &
               any(.not. (sed_acc_wdsrf == 0._r8)) .or. &
               any(.not. (sed_acc_rivsto == 0._r8)) .or. &
               any(.not. (sed_acc_rivout == 0._r8)) .or. &
               any(.not. (sed_acc_abs_rivout == 0._r8)) .or. &
               any(.not. (sed_acc_floodarea == 0._r8)) .or. &
               any(.not. (sed_precip == 0._r8)) .or. &
               any(.not. (sed_precip_yield == 0._r8)) .or. &
               any(.not. (a_sedcon == 0._r8)) .or. any(.not. (a_sedout == 0._r8)) .or. &
               any(.not. (a_bedout == 0._r8)) .or. any(.not. (a_sedinp == 0._r8)) .or. &
               any(.not. (a_netflw == 0._r8)) .or. any(.not. (a_layer == 0._r8)) .or. &
               any(.not. (a_shearvel == 0._r8)) .or. &
               .not. (sed_precip_time == 0._r8) .or. .not. (sed_hist_acctime == 0._r8)
         ENDIF
#ifdef USEMPI
         CALL mpi_allreduce(MPI_IN_PLACE, legacy_state_nonzero, 1, MPI_LOGICAL, &
            MPI_LOR, p_comm_glb, p_err)
#endif
         IF (legacy_state_nonzero) THEN
            IF (p_is_io) WRITE(*,'(A)') &
               'ERROR: legacy sediment restart lacks a descriptor for nonzero bed or accumulated state.'
            IF (p_is_master) ierr = nf90_close(ncid)
            CALL CoLM_stop()
         ENDIF
         IF (p_is_io) WRITE(*,'(A)') &
            'Sediment restart: migrated provably empty legacy state to canonical schema.'
      ENDIF

      CALL validate_sediment_checkpoint_state('read')

      IF (p_is_master) ierr = nf90_close(ncid)

      IF (p_is_io) WRITE(*,*) 'Sediment restart: read', nread, 'prognostic variables and strict accumulators.'

   END SUBROUTINE read_sediment_restart

   !-------------------------------------------------------------------------------------
   SUBROUTINE try_read_restart_var(file_restart, ncid, vname, buf, numucat, ucat_data_address, var_ok)
   ! Check if variable exists in restart file; if yes, read and scatter; if no, skip.
   !-------------------------------------------------------------------------------------
   USE netcdf
   USE MOD_Vector_ReadWrite
   USE MOD_DataType
   IMPLICIT NONE

   character(len=*), intent(in) :: file_restart, vname
   integer, intent(in) :: ncid, numucat
   type(pointer_int32_1d), intent(in) :: ucat_data_address(0:)
   real(r8), allocatable, intent(inout) :: buf(:)
   logical, intent(out) :: var_ok

   integer :: varid

      var_ok = .false.
      IF (p_is_master) THEN
         var_ok = (nf90_inq_varid(ncid, trim(vname), varid) == NF90_NOERR)
      ENDIF
#ifdef USEMPI
      CALL mpi_bcast(var_ok, 1, MPI_LOGICAL, p_address_master, p_comm_glb, p_err)
#endif

      IF (var_ok) THEN
         CALL vector_read_and_scatter(file_restart, buf, numucat, trim(vname), ucat_data_address)
      ELSE
         IF (p_is_io) WRITE(*,*) '  Sediment restart: variable "' // trim(vname) // '" not found, skipped.'
      ENDIF

   END SUBROUTINE try_read_restart_var

   !-------------------------------------------------------------------------------------
   SUBROUTINE require_sediment_restart_var(file_restart, ncid, vname, buf, numucat, ucat_data_address)
   ! Strict sediment restart read: once any sediment restart is present, all
   ! prognostic, queued, and history-continuity variables must be present.
   !-------------------------------------------------------------------------------------
   USE MOD_DataType
   IMPLICIT NONE

   character(len=*), intent(in) :: file_restart, vname
   integer, intent(in) :: ncid, numucat
   type(pointer_int32_1d), intent(in) :: ucat_data_address(0:)
   real(r8), allocatable, intent(inout) :: buf(:)
   logical :: var_ok

      CALL try_read_restart_var(file_restart, ncid, vname, buf, numucat, ucat_data_address, var_ok)
      IF (.not. var_ok) THEN
         IF (p_is_io) WRITE(*,'(A,A,A)') &
            'ERROR: incomplete sediment restart; required variable "', trim(vname), '" is missing.'
         CALL CoLM_stop()
      ENDIF

   END SUBROUTINE require_sediment_restart_var

   !-------------------------------------------------------------------------------------
   ELEMENTAL LOGICAL FUNCTION sediment_restart_real_matches(actual, expected)
   ! Restart metadata is written and read as r8.  The tolerance permits only
   ! roundoff-level serialization noise and deliberately rejects NaN/Inf.
   !-------------------------------------------------------------------------------------
   IMPLICIT NONE
   real(r8), intent(in) :: actual, expected
   real(r8) :: tol

      tol = 64._r8 * epsilon(1._r8) * max(1._r8, abs(expected))
      sediment_restart_real_matches = (actual == actual) .and. &
         (abs(actual) <= huge(actual)) .and. (abs(actual - expected) <= tol)

   END FUNCTION sediment_restart_real_matches

   !-------------------------------------------------------------------------------------
   ELEMENTAL LOGICAL FUNCTION sediment_restart_value_finite(value)
   !-------------------------------------------------------------------------------------
   IMPLICIT NONE
   real(r8), intent(in) :: value

      sediment_restart_value_finite = (value == value) .and. (abs(value) <= huge(value))

   END FUNCTION sediment_restart_value_finite

   !-------------------------------------------------------------------------------------
   SUBROUTINE validate_sediment_checkpoint_state(context)
   ! One collective validator is shared by restart read and write.  Signed
   ! transport accumulators need only be finite; mass, bed, carrier-time, and
   ! nonnegative diagnostic accumulators must also be >= 0.
   !-------------------------------------------------------------------------------------
   IMPLICIT NONE
   character(len=*), intent(in) :: context
   logical :: state_bad

      state_bad = .false.
      IF (p_is_worker) THEN
         state_bad = any(.not. sediment_restart_value_finite(sedcon)) .or. &
            any(sedcon < 0._r8) .or. &
            any(.not. sediment_restart_value_finite(sedsto)) .or. &
            any(sedsto < 0._r8) .or. &
            any(.not. sediment_restart_value_finite(layer)) .or. any(layer < 0._r8) .or. &
            any(.not. sediment_restart_value_finite(seddep)) .or. any(seddep < 0._r8) .or. &
            any(.not. sediment_restart_value_finite(sed_acc_time)) .or. any(sed_acc_time < 0._r8) .or. &
            any(.not. sediment_restart_value_finite(sed_acc_v2)) .or. any(sed_acc_v2 < 0._r8) .or. &
            any(.not. sediment_restart_value_finite(sed_acc_wdsrf)) .or. any(sed_acc_wdsrf < 0._r8) .or. &
            any(.not. sediment_restart_value_finite(sed_acc_rivsto)) .or. any(sed_acc_rivsto < 0._r8) .or. &
            any(.not. sediment_restart_value_finite(sed_acc_rivout)) .or. &
            any(.not. sediment_restart_value_finite(sed_acc_abs_rivout)) .or. &
            any(sed_acc_abs_rivout < 0._r8) .or. &
            any(.not. sediment_restart_value_finite(sed_acc_floodarea)) .or. &
            any(sed_acc_floodarea < 0._r8) .or. &
            any(.not. sediment_restart_value_finite(sed_precip)) .or. any(sed_precip < 0._r8) .or. &
            any(.not. sediment_restart_value_finite(sed_precip_yield)) .or. &
            any(sed_precip_yield < 0._r8) .or. &
            any(.not. sediment_restart_value_finite(a_sedcon)) .or. any(a_sedcon < 0._r8) .or. &
            any(.not. sediment_restart_value_finite(a_sedout)) .or. &
            any(.not. sediment_restart_value_finite(a_bedout)) .or. &
            any(.not. sediment_restart_value_finite(a_sedinp)) .or. any(a_sedinp < 0._r8) .or. &
            any(.not. sediment_restart_value_finite(a_netflw)) .or. &
            any(.not. sediment_restart_value_finite(a_layer)) .or. any(a_layer < 0._r8) .or. &
            any(.not. sediment_restart_value_finite(a_shearvel)) .or. any(a_shearvel < 0._r8) .or. &
            .not. sediment_restart_value_finite(sed_precip_time) .or. sed_precip_time < 0._r8 .or. &
            .not. sediment_restart_value_finite(sed_hist_acctime) .or. sed_hist_acctime < 0._r8
      ENDIF
#ifdef USEMPI
      CALL mpi_allreduce(MPI_IN_PLACE, state_bad, 1, MPI_LOGICAL, MPI_LOR, p_comm_glb, p_err)
#endif
      IF (state_bad) THEN
         IF (p_is_io) WRITE(*,'(A,A,A)') &
            'ERROR: sediment checkpoint ', trim(context), &
            ' contains non-finite or physically negative state.'
         CALL CoLM_stop()
      ENDIF

   END SUBROUTINE validate_sediment_checkpoint_state

   !-------------------------------------------------------------------------------------
   SUBROUTINE check_sediment_restart_scalar(file_restart, ncid, vname, expected, &
      numucat, ucat_data_address, meta_bad)
   !-------------------------------------------------------------------------------------
   USE MOD_DataType
   IMPLICIT NONE

   character(len=*), intent(in) :: file_restart, vname
   integer, intent(in) :: ncid, numucat
   real(r8), intent(in) :: expected
   type(pointer_int32_1d), intent(in) :: ucat_data_address(0:)
   logical, intent(inout) :: meta_bad
   real(r8), allocatable :: buf(:)

      CALL require_sediment_restart_var(file_restart, ncid, vname, buf, &
         numucat, ucat_data_address)
      IF (allocated(buf)) THEN
         IF (p_is_worker .and. numucat > 0) THEN
            meta_bad = meta_bad .or. any(.not. sediment_restart_real_matches(buf, expected))
         ENDIF
         deallocate(buf)
      ENDIF

   END SUBROUTINE check_sediment_restart_scalar

   !-------------------------------------------------------------------------------------
   SUBROUTINE write_sediment_scalar_meta(file_restart, vname, value)
   !-------------------------------------------------------------------------------------
   USE MOD_Vector_ReadWrite
   USE MOD_Grid_RiverLakeNetwork, only: numucat, totalnumucat, ucat_data_address
   IMPLICIT NONE

   character(len=*), intent(in) :: file_restart, vname
   real(r8), intent(in) :: value
   real(r8), allocatable :: scalar_vec(:)

      IF (p_is_worker) THEN
         allocate(scalar_vec(numucat))
         scalar_vec(:) = value
      ELSE
         allocate(scalar_vec(0))
      ENDIF
      CALL vector_gather_and_write(scalar_vec, size(scalar_vec), totalnumucat, &
         ucat_data_address, file_restart, vname, 'ucatch')
      deallocate(scalar_vec)

   END SUBROUTINE write_sediment_scalar_meta

   !-------------------------------------------------------------------------------------
   SUBROUTINE write_sediment_restart(file_restart)
   !-------------------------------------------------------------------------------------
   USE MOD_Vector_ReadWrite
   USE MOD_Grid_RiverLakeNetwork, only: numucat, totalnumucat, ucat_data_address, &
      topo_rivwth, topo_rivlen
   IMPLICIT NONE

   character(len=*), intent(in) :: file_restart
   integer :: ised, ilyr
   character(len=16) :: cised, cilyr
   real(r8) :: dummy_sed(1)   ! Dummy for non-worker processes (vlen=0, never accessed)
   real(r8), allocatable :: scalar_vec(:)

      ! All processes must participate (MPI collective calls inside vector_gather_and_write).
      IF (.not. sediment_particle_enabled()) RETURN
      CALL validate_sediment_checkpoint_state('write')

      ! Invalidate the transaction before writing any payload.  The final
      ! schema-valued completion marker is written only after every field.
      CALL write_sediment_scalar_meta(file_restart, 'sed_restart_schema_meta', &
         real(SED_RESTART_SCHEMA_VERSION, r8))
      CALL write_sediment_scalar_meta(file_restart, 'sed_restart_complete_meta', 0._r8)
      CALL write_sediment_scalar_meta(file_restart, 'sed_n_meta', real(nsed, r8))
      CALL write_sediment_scalar_meta(file_restart, 'sed_totlyrnum_meta', real(totlyrnum, r8))
      CALL write_sediment_scalar_meta(file_restart, 'sed_nlfp_meta', real(nlfp_sed, r8))
      CALL write_sediment_scalar_meta(file_restart, 'sed_lambda_meta', lambda)
      CALL write_sediment_scalar_meta(file_restart, 'sed_lyrdph_meta', lyrdph)
      CALL write_sediment_scalar_meta(file_restart, 'sed_psedd_meta', psedD)
      CALL write_sediment_scalar_meta(file_restart, 'sed_pwatd_meta', pwatD)
      CALL write_sediment_scalar_meta(file_restart, 'sed_viskin_meta', visKin)
      CALL write_sediment_scalar_meta(file_restart, 'sed_vonkar_meta', vonKar)
      CALL write_sediment_scalar_meta(file_restart, 'sed_pset_meta', pset)
      CALL write_sediment_scalar_meta(file_restart, 'sed_ignore_dph_meta', sed_ignore_dph)
      CALL write_sediment_scalar_meta(file_restart, 'sed_cfl_adv_meta', sed_cfl_adv)
      CALL write_sediment_scalar_meta(file_restart, 'sed_dt_max_meta', sed_dt_max)
      CALL write_sediment_scalar_meta(file_restart, 'sed_bed_depth_meta', sed_bed_depth)
      CALL write_sediment_scalar_meta(file_restart, 'sed_pyld_meta', pyld)
      CALL write_sediment_scalar_meta(file_restart, 'sed_pyldc_meta', pyldc)
      CALL write_sediment_scalar_meta(file_restart, 'sed_pyldpc_meta', pyldpc)
      CALL write_sediment_scalar_meta(file_restart, 'sed_dsylunit_meta', dsylunit)
      CALL write_sediment_scalar_meta(file_restart, 'sed_max_conc_meta', sed_max_conc)
      CALL write_sediment_scalar_meta(file_restart, 'sed_bedload_coeff_meta', SED_BEDLOAD_COEFF)
      CALL write_sediment_scalar_meta(file_restart, 'sed_precip_threshold_meta', &
         SED_PRECIP_THRESHOLD_MM_DAY)
      CALL write_sediment_scalar_meta(file_restart, 'sed_exch_shear_min_meta', EXCH_SHEARVEL_MIN)
      CALL write_sediment_scalar_meta(file_restart, 'sed_exch_shear_blend_meta', EXCH_SHEARVEL_BLEND)
      CALL write_sediment_scalar_meta(file_restart, 'sed_exch_zd_max_meta', EXCH_ZD_MAX)

      DO ised = 1, nsed
         WRITE(cised, '(I0)') ised
         CALL write_sediment_scalar_meta(file_restart, &
            'sed_diam_meta_' // trim(cised), sDiam(ised))
         CALL write_sediment_scalar_meta(file_restart, &
            'sed_setvel_meta_' // trim(cised), setvel(ised))
         IF (p_is_worker) THEN
            CALL vector_gather_and_write(sed_frc(ised,:), numucat, totalnumucat, &
               ucat_data_address, file_restart, 'sed_frc_meta_' // trim(cised), 'ucatch')
         ELSE
            CALL vector_gather_and_write(dummy_sed, 0, totalnumucat, ucat_data_address, &
               file_restart, 'sed_frc_meta_' // trim(cised), 'ucatch')
         ENDIF
      ENDDO
      DO ilyr = 1, nlfp_sed
         WRITE(cilyr, '(I0)') ilyr
         IF (p_is_worker) THEN
            CALL vector_gather_and_write(sed_slope(ilyr,:), numucat, totalnumucat, &
               ucat_data_address, file_restart, 'sed_slope_meta_' // trim(cilyr), 'ucatch')
         ELSE
            CALL vector_gather_and_write(dummy_sed, 0, totalnumucat, ucat_data_address, &
               file_restart, 'sed_slope_meta_' // trim(cilyr), 'ucatch')
         ENDIF
      ENDDO
      IF (p_is_worker) THEN
         CALL vector_gather_and_write(topo_rivwth, numucat, totalnumucat, &
            ucat_data_address, file_restart, 'sed_rivwth_meta', 'ucatch')
         CALL vector_gather_and_write(topo_rivlen, numucat, totalnumucat, &
            ucat_data_address, file_restart, 'sed_rivlen_meta', 'ucatch')
      ELSE
         CALL vector_gather_and_write(dummy_sed, 0, totalnumucat, ucat_data_address, &
            file_restart, 'sed_rivwth_meta', 'ucatch')
         CALL vector_gather_and_write(dummy_sed, 0, totalnumucat, ucat_data_address, &
            file_restart, 'sed_rivlen_meta', 'ucatch')
      ENDIF

      DO ised = 1, nsed
         WRITE(cised, '(I0)') ised
         IF (allocated(sedcon)) THEN
            CALL vector_gather_and_write (&
               sedcon(ised,:), numucat, totalnumucat, ucat_data_address, file_restart, &
               'sedcon_' // trim(cised), 'ucatch')
         ELSE
            CALL vector_gather_and_write (&
               dummy_sed, 0, totalnumucat, ucat_data_address, file_restart, &
               'sedcon_' // trim(cised), 'ucatch')
         ENDIF
      ENDDO

      DO ised = 1, nsed
         WRITE(cised, '(I0)') ised
         IF (allocated(sedsto)) THEN
            CALL vector_gather_and_write (&
               sedsto(ised,:), numucat, totalnumucat, ucat_data_address, file_restart, &
               'sedsto_' // trim(cised), 'ucatch')
         ELSE
            CALL vector_gather_and_write (&
               dummy_sed, 0, totalnumucat, ucat_data_address, file_restart, &
               'sedsto_' // trim(cised), 'ucatch')
         ENDIF
      ENDDO

      DO ised = 1, nsed
         WRITE(cised, '(I0)') ised
         IF (allocated(layer)) THEN
            CALL vector_gather_and_write (&
               layer(ised,:), numucat, totalnumucat, ucat_data_address, file_restart, &
               'layer_' // trim(cised), 'ucatch')
         ELSE
            CALL vector_gather_and_write (&
               dummy_sed, 0, totalnumucat, ucat_data_address, file_restart, &
               'layer_' // trim(cised), 'ucatch')
         ENDIF
      ENDDO

      DO ised = 1, nsed
         WRITE(cised, '(I0)') ised
         DO ilyr = 1, totlyrnum
            WRITE(cilyr, '(I0)') ilyr
            IF (allocated(seddep)) THEN
               CALL vector_gather_and_write (&
                  seddep(ised,ilyr,:), numucat, totalnumucat, ucat_data_address, file_restart, &
                  'seddep_' // trim(cised) // '_' // trim(cilyr), 'ucatch')
            ELSE
               CALL vector_gather_and_write (&
                  dummy_sed, 0, totalnumucat, ucat_data_address, file_restart, &
                  'seddep_' // trim(cised) // '_' // trim(cilyr), 'ucatch')
            ENDIF
         ENDDO
      ENDDO

      IF (allocated(sed_acc_time)) THEN
         CALL vector_gather_and_write(sed_acc_time, numucat, totalnumucat, ucat_data_address, &
            file_restart, 'sed_acc_time', 'ucatch')
      ELSE
         CALL vector_gather_and_write(dummy_sed, 0, totalnumucat, ucat_data_address, &
            file_restart, 'sed_acc_time', 'ucatch')
      ENDIF

      IF (allocated(sed_acc_v2)) THEN
         CALL vector_gather_and_write(sed_acc_v2, numucat, totalnumucat, ucat_data_address, &
            file_restart, 'sed_acc_v2', 'ucatch')
      ELSE
         CALL vector_gather_and_write(dummy_sed, 0, totalnumucat, ucat_data_address, &
            file_restart, 'sed_acc_v2', 'ucatch')
      ENDIF

      IF (allocated(sed_acc_wdsrf)) THEN
         CALL vector_gather_and_write(sed_acc_wdsrf, numucat, totalnumucat, ucat_data_address, &
            file_restart, 'sed_acc_wdsrf', 'ucatch')
      ELSE
         CALL vector_gather_and_write(dummy_sed, 0, totalnumucat, ucat_data_address, &
            file_restart, 'sed_acc_wdsrf', 'ucatch')
      ENDIF

      IF (allocated(sed_acc_rivsto)) THEN
         CALL vector_gather_and_write(sed_acc_rivsto, numucat, totalnumucat, ucat_data_address, &
            file_restart, 'sed_acc_rivsto', 'ucatch')
      ELSE
         CALL vector_gather_and_write(dummy_sed, 0, totalnumucat, ucat_data_address, &
            file_restart, 'sed_acc_rivsto', 'ucatch')
      ENDIF

      IF (allocated(sed_acc_rivout)) THEN
         CALL vector_gather_and_write(sed_acc_rivout, numucat, totalnumucat, ucat_data_address, &
            file_restart, 'sed_acc_rivout', 'ucatch')
      ELSE
         CALL vector_gather_and_write(dummy_sed, 0, totalnumucat, ucat_data_address, &
            file_restart, 'sed_acc_rivout', 'ucatch')
      ENDIF

      IF (allocated(sed_acc_abs_rivout)) THEN
         CALL vector_gather_and_write(sed_acc_abs_rivout, numucat, totalnumucat, ucat_data_address, &
            file_restart, 'sed_acc_abs_rivout', 'ucatch')
      ELSE
         CALL vector_gather_and_write(dummy_sed, 0, totalnumucat, ucat_data_address, &
            file_restart, 'sed_acc_abs_rivout', 'ucatch')
      ENDIF

      IF (allocated(sed_acc_floodarea)) THEN
         CALL vector_gather_and_write(sed_acc_floodarea, numucat, totalnumucat, ucat_data_address, &
            file_restart, 'sed_acc_floodarea', 'ucatch')
      ELSE
         CALL vector_gather_and_write(dummy_sed, 0, totalnumucat, ucat_data_address, &
            file_restart, 'sed_acc_floodarea', 'ucatch')
      ENDIF

      IF (allocated(sed_precip)) THEN
         CALL vector_gather_and_write(sed_precip, numucat, totalnumucat, ucat_data_address, &
            file_restart, 'sed_precip', 'ucatch')
      ELSE
         CALL vector_gather_and_write(dummy_sed, 0, totalnumucat, ucat_data_address, &
            file_restart, 'sed_precip', 'ucatch')
      ENDIF

      IF (allocated(sed_precip_yield)) THEN
         CALL vector_gather_and_write(sed_precip_yield, numucat, totalnumucat, ucat_data_address, &
            file_restart, 'sed_precip_yield', 'ucatch')
      ELSE
         CALL vector_gather_and_write(dummy_sed, 0, totalnumucat, ucat_data_address, &
            file_restart, 'sed_precip_yield', 'ucatch')
      ENDIF

      IF (p_is_worker) THEN
         allocate(scalar_vec(numucat))
         scalar_vec(:) = sed_precip_time
      ELSE
         allocate(scalar_vec(0))
      ENDIF
      CALL vector_gather_and_write(scalar_vec, size(scalar_vec), totalnumucat, ucat_data_address, &
         file_restart, 'sed_precip_time_vec', 'ucatch')
      IF (allocated(scalar_vec)) deallocate(scalar_vec)

      DO ised = 1, nsed
         WRITE(cised, '(I0)') ised
         IF (allocated(a_sedcon)) THEN
            CALL vector_gather_and_write(a_sedcon(ised,:), numucat, totalnumucat, ucat_data_address, &
               file_restart, 'a_sedcon_' // trim(cised), 'ucatch')
            CALL vector_gather_and_write(a_sedout(ised,:), numucat, totalnumucat, ucat_data_address, &
               file_restart, 'a_sedout_' // trim(cised), 'ucatch')
            CALL vector_gather_and_write(a_bedout(ised,:), numucat, totalnumucat, ucat_data_address, &
               file_restart, 'a_bedout_' // trim(cised), 'ucatch')
            CALL vector_gather_and_write(a_sedinp(ised,:), numucat, totalnumucat, ucat_data_address, &
               file_restart, 'a_sedinp_' // trim(cised), 'ucatch')
            CALL vector_gather_and_write(a_netflw(ised,:), numucat, totalnumucat, ucat_data_address, &
               file_restart, 'a_netflw_' // trim(cised), 'ucatch')
            CALL vector_gather_and_write(a_layer(ised,:), numucat, totalnumucat, ucat_data_address, &
               file_restart, 'a_layer_' // trim(cised), 'ucatch')
         ELSE
            CALL vector_gather_and_write(dummy_sed, 0, totalnumucat, ucat_data_address, &
               file_restart, 'a_sedcon_' // trim(cised), 'ucatch')
            CALL vector_gather_and_write(dummy_sed, 0, totalnumucat, ucat_data_address, &
               file_restart, 'a_sedout_' // trim(cised), 'ucatch')
            CALL vector_gather_and_write(dummy_sed, 0, totalnumucat, ucat_data_address, &
               file_restart, 'a_bedout_' // trim(cised), 'ucatch')
            CALL vector_gather_and_write(dummy_sed, 0, totalnumucat, ucat_data_address, &
               file_restart, 'a_sedinp_' // trim(cised), 'ucatch')
            CALL vector_gather_and_write(dummy_sed, 0, totalnumucat, ucat_data_address, &
               file_restart, 'a_netflw_' // trim(cised), 'ucatch')
            CALL vector_gather_and_write(dummy_sed, 0, totalnumucat, ucat_data_address, &
               file_restart, 'a_layer_' // trim(cised), 'ucatch')
         ENDIF
      ENDDO

      IF (allocated(a_shearvel)) THEN
         CALL vector_gather_and_write(a_shearvel, numucat, totalnumucat, ucat_data_address, &
            file_restart, 'a_shearvel', 'ucatch')
      ELSE
         CALL vector_gather_and_write(dummy_sed, 0, totalnumucat, ucat_data_address, &
            file_restart, 'a_shearvel', 'ucatch')
      ENDIF

      IF (p_is_worker) THEN
         allocate(scalar_vec(numucat))
         scalar_vec(:) = sed_hist_acctime
      ELSE
         allocate(scalar_vec(0))
      ENDIF
      CALL vector_gather_and_write(scalar_vec, size(scalar_vec), totalnumucat, ucat_data_address, &
         file_restart, 'sed_hist_acctime_vec', 'ucatch')
      IF (allocated(scalar_vec)) deallocate(scalar_vec)

      CALL write_sediment_scalar_meta(file_restart, 'sed_restart_complete_meta', &
         1._r8)

   END SUBROUTINE write_sediment_restart

   !-------------------------------------------------------------------------------------
   SUBROUTINE write_sediment_history (file_hist_ucat, itime_in_file_ucat)
   !-------------------------------------------------------------------------------------
   USE MOD_Grid_RiverLakeNetwork, only: numucat
   ! Route history is dispatched, not written directly: the same call must
   ! land in the single file under DEF_HIST_mode='one' and in this IO group's
   ! shard under 'block'.
   USE MOD_Grid_RiverLakeHistRoute, only: route_hist_write_ucat
   IMPLICIT NONE

   character(len=*), intent(in) :: file_hist_ucat
   integer, intent(in) :: itime_in_file_ucat

   real(r8), allocatable :: a_sedcon_avg(:,:)
   real(r8), allocatable :: a_sedout_avg(:,:)
   real(r8), allocatable :: a_bedout_avg(:,:)
   real(r8), allocatable :: a_sedinp_avg(:,:)
   real(r8), allocatable :: a_netflw_avg(:,:)
   real(r8), allocatable :: a_layer_avg(:,:)
   real(r8), allocatable :: a_shearvel_avg(:)
   integer :: ised
   character(len=16) :: cised

      IF (.not. sediment_particle_enabled()) RETURN

      ! Allocate on ALL processes (zero-size on non-workers) to avoid
      ! passing unallocated arrays to the route-history dispatcher.
      IF (p_is_worker .and. numucat > 0) THEN
         allocate (a_sedcon_avg  (nsed, numucat))
         allocate (a_sedout_avg  (nsed, numucat))
         allocate (a_bedout_avg  (nsed, numucat))
         allocate (a_sedinp_avg  (nsed, numucat))
         allocate (a_netflw_avg  (nsed, numucat))
         allocate (a_layer_avg   (nsed, numucat))
         allocate (a_shearvel_avg(numucat))

         IF (sed_hist_acctime > 0._r8) THEN
            a_shearvel_avg = a_shearvel / sed_hist_acctime
            DO ised = 1, nsed
               a_sedcon_avg(ised,:) = a_sedcon(ised,:) / sed_hist_acctime
               a_sedout_avg(ised,:) = a_sedout(ised,:) / sed_hist_acctime
               a_bedout_avg(ised,:) = a_bedout(ised,:) / sed_hist_acctime
               a_sedinp_avg(ised,:) = a_sedinp(ised,:) / sed_hist_acctime
               a_netflw_avg(ised,:) = a_netflw(ised,:) / sed_hist_acctime
               a_layer_avg(ised,:)  = a_layer(ised,:)  / sed_hist_acctime
            ENDDO
         ELSE
            a_shearvel_avg = 0.
            a_sedcon_avg   = 0.
            a_sedout_avg   = 0.
            a_bedout_avg   = 0.
            a_sedinp_avg   = 0.
            a_netflw_avg   = 0.
            a_layer_avg    = 0.
         ENDIF
      ELSE
         ! Allocate with nsed in first dim so a_xxx_avg(ised,:) is a valid zero-length slice.
         allocate (a_sedcon_avg  (nsed, 0))
         allocate (a_sedout_avg  (nsed, 0))
         allocate (a_bedout_avg  (nsed, 0))
         allocate (a_sedinp_avg  (nsed, 0))
         allocate (a_netflw_avg  (nsed, 0))
         allocate (a_layer_avg   (nsed, 0))
         allocate (a_shearvel_avg(0))
      ENDIF

      IF (DEF_hist_vars%sedcon) THEN
         DO ised = 1, nsed
            WRITE(cised, '(I0)') ised
            CALL route_hist_write_ucat (a_sedcon_avg(ised,:), 'f_sedcon_' // trim(cised), &
               longname = 'suspended sediment concentration, size class ' // trim(cised), &
               units = 'm^3/m^3')
         ENDDO
      ENDIF

      IF (DEF_hist_vars%sedout) THEN
         DO ised = 1, nsed
            WRITE(cised, '(I0)') ised
            CALL route_hist_write_ucat (a_sedout_avg(ised,:), 'f_sedout_' // trim(cised), &
               longname = 'suspended sediment flux, size class ' // trim(cised), &
               units = 'm^3/s')
         ENDDO
      ENDIF

      IF (DEF_hist_vars%bedout) THEN
         DO ised = 1, nsed
            WRITE(cised, '(I0)') ised
            CALL route_hist_write_ucat (a_bedout_avg(ised,:), 'f_bedout_' // trim(cised), &
               longname = 'bedload solid-volume flux, size class ' // trim(cised), &
               units = 'm^3/s')
         ENDDO
      ENDIF

      IF (DEF_hist_vars%sedinp) THEN
         DO ised = 1, nsed
            WRITE(cised, '(I0)') ised
            CALL route_hist_write_ucat (a_sedinp_avg(ised,:), 'f_sedinp_' // trim(cised), &
               longname = 'sediment erosion input, size class ' // trim(cised), &
               units = 'm^3/s')
         ENDDO
      ENDIF

      IF (DEF_hist_vars%netflw) THEN
         DO ised = 1, nsed
            WRITE(cised, '(I0)') ised
            CALL route_hist_write_ucat (a_netflw_avg(ised,:), 'f_netflw_' // trim(cised), &
               longname = 'net bed-water exchange flux (incl. shallow deposit), size class ' // trim(cised), &
               units = 'm^3/s')
         ENDDO
      ENDIF

      IF (DEF_hist_vars%sedlayer) THEN
         DO ised = 1, nsed
            WRITE(cised, '(I0)') ised
            CALL route_hist_write_ucat (a_layer_avg(ised,:), 'f_layer_' // trim(cised), &
               longname = 'active layer storage, size class ' // trim(cised), &
               units = 'm^3')
         ENDDO
      ENDIF

      IF (DEF_hist_vars%shearvel) THEN
         CALL route_hist_write_ucat (a_shearvel_avg, 'f_shearvel', &
            longname = 'shear velocity', &
            units = 'm/s')
      ENDIF

      IF (allocated(a_sedcon_avg  )) deallocate (a_sedcon_avg  )
      IF (allocated(a_sedout_avg  )) deallocate (a_sedout_avg  )
      IF (allocated(a_bedout_avg  )) deallocate (a_bedout_avg  )
      IF (allocated(a_sedinp_avg  )) deallocate (a_sedinp_avg  )
      IF (allocated(a_netflw_avg  )) deallocate (a_netflw_avg  )
      IF (allocated(a_layer_avg   )) deallocate (a_layer_avg   )
      IF (allocated(a_shearvel_avg)) deallocate (a_shearvel_avg)

   END SUBROUTINE write_sediment_history

   !-------------------------------------------------------------------------------------
   SUBROUTINE flush_sediment_history()
   !-------------------------------------------------------------------------------------
   IMPLICIT NONE

      IF (.not. sediment_particle_enabled()) RETURN
      IF (.not. allocated(a_sedcon)) RETURN

      a_sedcon   = 0.
      a_sedout   = 0.
      a_bedout   = 0.
      a_sedinp   = 0.
      a_netflw   = 0.
      a_layer    = 0.
      a_shearvel = 0.
      sed_hist_acctime = 0.

   END SUBROUTINE flush_sediment_history

   !-------------------------------------------------------------------------------------
   SUBROUTINE grid_sediment_final()
   !-------------------------------------------------------------------------------------
   IMPLICIT NONE
      IF (allocated(sed_frc      )) deallocate(sed_frc      )
      IF (allocated(sed_slope    )) deallocate(sed_slope    )
      IF (allocated(sDiam        )) deallocate(sDiam        )
      IF (allocated(sDiam_from_param)) deallocate(sDiam_from_param)
      IF (allocated(setvel       )) deallocate(setvel       )
      IF (allocated(sedcon       )) deallocate(sedcon       )
      IF (allocated(sedsto       )) deallocate(sedsto       )
      IF (allocated(layer        )) deallocate(layer        )
      IF (allocated(seddep       )) deallocate(seddep       )
      IF (allocated(sedout       )) deallocate(sedout       )
      IF (allocated(bedout       )) deallocate(bedout       )
      IF (allocated(sedinp       )) deallocate(sedinp       )
      IF (allocated(netflw       )) deallocate(netflw       )
      IF (allocated(exch_es_raw  )) deallocate(exch_es_raw  )
      IF (allocated(exch_d_raw   )) deallocate(exch_d_raw   )
      IF (allocated(exch_es_eff  )) deallocate(exch_es_eff  )
      IF (allocated(exch_d_eff   )) deallocate(exch_d_eff   )
      IF (allocated(netflw_adv_step)) deallocate(netflw_adv_step)
      IF (allocated(exch_d_adv_step)) deallocate(exch_d_adv_step)
      IF (allocated(shearvel     )) deallocate(shearvel     )
      IF (allocated(critshearvel )) deallocate(critshearvel )
      IF (allocated(susvel       )) deallocate(susvel       )
      IF (allocated(sed_acc_time )) deallocate(sed_acc_time )
      IF (allocated(sed_acc_v2   )) deallocate(sed_acc_v2   )
      IF (allocated(sed_acc_wdsrf)) deallocate(sed_acc_wdsrf)
      IF (allocated(sed_acc_rivsto)) deallocate(sed_acc_rivsto)
      IF (allocated(sed_acc_rivout)) deallocate(sed_acc_rivout)
      IF (allocated(sed_acc_abs_rivout)) deallocate(sed_acc_abs_rivout)
      IF (allocated(sed_acc_floodarea)) deallocate(sed_acc_floodarea)
      IF (allocated(sed_acc_carrier_time)) deallocate(sed_acc_carrier_time)
      IF (allocated(sed_acc_wdsrf_min)) deallocate(sed_acc_wdsrf_min)
      IF (allocated(sed_acc_wdsrf_max)) deallocate(sed_acc_wdsrf_max)
      IF (allocated(sed_acc_rivsto_min)) deallocate(sed_acc_rivsto_min)
      IF (allocated(sed_acc_rivsto_max)) deallocate(sed_acc_rivsto_max)
      IF (allocated(sed_acc_rivout_min)) deallocate(sed_acc_rivout_min)
      IF (allocated(sed_acc_rivout_max)) deallocate(sed_acc_rivout_max)
      IF (allocated(sed_acc_pos_rivout)) deallocate(sed_acc_pos_rivout)
      IF (allocated(sed_acc_neg_rivout)) deallocate(sed_acc_neg_rivout)
      IF (allocated(sed_acc_near_dry_abs_rivout)) deallocate(sed_acc_near_dry_abs_rivout)
      IF (allocated(sed_precip   )) deallocate(sed_precip   )
      IF (allocated(sed_precip_yield)) deallocate(sed_precip_yield)
      IF (allocated(a_sedcon     )) deallocate(a_sedcon     )
      IF (allocated(a_sedout     )) deallocate(a_sedout     )
      IF (allocated(a_bedout     )) deallocate(a_bedout     )
      IF (allocated(a_sedinp     )) deallocate(a_sedinp     )
      IF (allocated(a_netflw     )) deallocate(a_netflw     )
      IF (allocated(a_layer      )) deallocate(a_layer      )
      IF (allocated(a_shearvel   )) deallocate(a_shearvel   )
#ifdef CoLMDEBUG
      IF (allocated(sed_diag_station_local_i)) deallocate(sed_diag_station_local_i)
#endif
   END SUBROUTINE grid_sediment_final

   !-------------------------------------------------------------------------------------
   SUBROUTINE print_sediment_runtime_parameters()
   !-------------------------------------------------------------------------------------
   IMPLICIT NONE

   integer :: ised

   ! Print exactly once.
   IF (.not. p_is_master) RETURN

   WRITE(*,'(A)') ' '
   WRITE(*,'(A)') '============================================================'
   WRITE(*,'(A)') 'SEDIMENT EFFECTIVE PARAMETERS'
   WRITE(*,'(A)') '============================================================'

   WRITE(*,'(A,I0)')       'nsed                         = ', nsed
   WRITE(*,'(A,I0)')       'nlfp_sed                     = ', nlfp_sed
   WRITE(*,'(A,I0)')       'ndeposit_layers              = ', totlyrnum

   WRITE(*,'(A,ES20.10)')  'grain_density [kg/m3]        = ', psedD * 1000._r8
   WRITE(*,'(A,ES20.10)')  'water_density [kg/m3]        = ', pwatD * 1000._r8
   WRITE(*,'(A,ES20.10)')  'porosity                     = ', lambda

   WRITE(*,'(A,ES20.10)')  'ignore_depth_m [m]           = ', sed_ignore_dph
   WRITE(*,'(A,ES20.10)')  'active_layer_depth [m]       = ', lyrdph
   WRITE(*,'(A,ES20.10)')  'bed_depth [m]                = ', sed_bed_depth

   WRITE(*,'(A,ES20.10)')  'viscosity [m2/s]             = ', visKin
   WRITE(*,'(A,ES20.10)')  'von_karman                   = ', vonKar
   WRITE(*,'(A,ES20.10)')  'settling_multiplier          = ', pset

   WRITE(*,'(A,ES20.10)')  'yield_coefficient            = ', pyld
   WRITE(*,'(A,ES20.10)')  'slope_exponent               = ', pyldc
   WRITE(*,'(A,ES20.10)')  'precipitation_exponent       = ', pyldpc
   WRITE(*,'(A,ES20.10)')  'unit_conversion              = ', dsylunit

   WRITE(*,'(A,ES20.10)')  'cfl_adv                      = ', sed_cfl_adv
   WRITE(*,'(A,ES20.10)')  'max_timestep_s [s]           = ', sed_dt_max

   WRITE(*,'(A)') ' '
   WRITE(*,'(A)') '--- Grain classes ---'

   DO ised = 1, nsed
      WRITE(*,'(A,I0,A,ES20.10)') &
         'grain_diameter(', ised, ') [m] = ', sDiam(ised)

      WRITE(*,'(A,I0,A,ES20.10)') &
         'settling_velocity(', ised, ') [m/s] = ', setvel(ised)
   ENDDO

   WRITE(*,'(A)') ' '
   WRITE(*,'(A)') '--- Compiled sediment numerical constants ---'

   WRITE(*,'(A,ES20.10)') 'SED_NEAR_DRY_DEPTH [m]       = ', SED_NEAR_DRY_DEPTH
   WRITE(*,'(A,I0)')      'SED_MAX_ADV_SUBSTEPS          = ', SED_MAX_ADV_SUBSTEPS
   WRITE(*,'(A,ES20.10)') 'SED_BALANCE_ABS_TOL           = ', SED_BALANCE_ABS_TOL
   WRITE(*,'(A,ES20.10)') 'SED_BALANCE_REL_TOL           = ', SED_BALANCE_REL_TOL
   WRITE(*,'(A,ES20.10)') 'SED_BEDLOAD_COEFF             = ', SED_BEDLOAD_COEFF
   WRITE(*,'(A,ES20.10)') 'SED_PRECIP_THRESHOLD_MM_DAY   = ', &
      SED_PRECIP_THRESHOLD_MM_DAY

   WRITE(*,'(A)') '============================================================'
   WRITE(*,'(A)') ' '

   END SUBROUTINE print_sediment_runtime_parameters
   !-------------------------------------------------------------------------------------

END MODULE MOD_Tracer_Particle_Sediment
#endif
