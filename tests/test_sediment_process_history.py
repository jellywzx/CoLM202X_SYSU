"""Compile shipped sediment operators and verify their diagnostic NetCDF records."""

from pathlib import Path
import re
import subprocess

import numpy as np
import pytest
from netCDF4 import Dataset

from fortran_test_support import netcdf_fortran_flags, require_runnable_fortran_compiler


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "main/TRACER/MOD_Tracer_Particle_Sediment.F90"


def routine(source, name):
    match = re.search(
        rf"^\s*(?:SUBROUTINE|ELEMENTAL LOGICAL FUNCTION) {name}\b.*?"
        rf"^\s*END (?:SUBROUTINE|FUNCTION) {name}\b", source, re.M | re.S | re.I,
    )
    assert match is not None, name
    return match.group()


STUBS = r"""
module MOD_Precision
 integer, parameter :: r8=kind(1.d0)
end module
module MOD_SPMD_Task
 logical :: p_is_worker=.true., p_is_master=.true., p_is_io=.true.
contains
 subroutine CoLM_stop(message)
  character(*), optional, intent(in) :: message
  if(present(message)) print *, message
  error stop 99
 end subroutine
end module
module MOD_Namelist
 logical :: DEF_USE_BIFURCATION=.false., DEF_USE_LEVEE=.false.
 integer :: DEF_REST_CompressLevel=0
 type history_flags
  logical :: sedcon=.false.,sedout=.false.,bedout=.false.,sedinp=.false.
  logical :: netflw=.false.,sedlayer=.false.,shearvel=.false.,seddiag=.true.
 end type
 type(history_flags) :: DEF_hist_vars
contains
 function get_unitcatchment_file() result(name)
  character(len=1) :: name
  name=' '
 end function
end module
module MOD_Vars_Global
 use MOD_Precision
 real(r8), parameter :: spval=-9999._r8
end module
module MOD_DataType
 type pointer_int32_1d
  integer, pointer :: val(:)=>null()
 end type
end module
module MOD_Grid_RiverLakeNetwork
 use MOD_Precision
 use MOD_DataType
 integer :: numucat=2,totalnumucat=2,npthlev_bif=0,npthout_local=0
 integer :: push_next2ucat=0,push_ups2ucat=0
 real(r8) :: topo_rivwth(2)=1._r8
 type(pointer_int32_1d) :: ucat_data_address(0:0)
end module
module MOD_Const_Physical
 use MOD_Precision
 real(r8), parameter :: grav=9.81_r8
end module
module MOD_WorkerPushData
 use MOD_Precision
contains
 subroutine worker_push_data(mapping,source,dest,fillvalue,mode)
  integer,intent(in) :: mapping
  real(r8),intent(in) :: source(:),fillvalue
  real(r8),intent(out) :: dest(:)
  character(*),optional,intent(in) :: mode
  dest=fillvalue
 end subroutine
end module
module probe_netcdf
 use netcdf
 use MOD_Precision
 use MOD_Vars_Global
contains
 subroutine check(status)
  integer,intent(in) :: status
  if(status/=nf90_noerr) then
   print *, nf90_strerror(status)
   error stop 98
  endif
 end subroutine
 subroutine write_row(file,name,row,units,longname)
  character(*),intent(in) :: file,name
  real(r8),intent(in) :: row(:)
  character(*),optional,intent(in) :: units,longname
  integer :: ncid,dimid,varid,status
  logical :: exists
  if(size(row)==0) return
  inquire(file=file,exist=exists)
  if(exists) then
   call check(nf90_open(file,nf90_write,ncid))
   call check(nf90_redef(ncid))
   call check(nf90_inq_dimid(ncid,'ucatch',dimid))
  else
   call check(nf90_create(file,nf90_clobber,ncid))
   call check(nf90_def_dim(ncid,'ucatch',size(row),dimid))
  endif
  status=nf90_inq_varid(ncid,name,varid)
  if(status/=nf90_noerr) then
   call check(nf90_def_var(ncid,name,nf90_double,[dimid],varid))
   call check(nf90_put_att(ncid,varid,'_FillValue',spval))
  endif
  if(present(units)) call check(nf90_put_att(ncid,varid,'units',units))
  if(present(longname)) call check(nf90_put_att(ncid,varid,'long_name',longname))
  call check(nf90_enddef(ncid))
  call check(nf90_put_var(ncid,varid,row))
  call check(nf90_close(ncid))
 end subroutine
end module
module MOD_Vector_ReadWrite
 use MOD_Precision
 use MOD_DataType
 use probe_netcdf
contains
 subroutine vector_gather_and_write(row,nlocal,ntotal,address,file,name,dimname)
  real(r8),intent(in) :: row(:)
  integer,intent(in) :: nlocal,ntotal
  type(pointer_int32_1d),intent(in) :: address(0:)
  character(*),intent(in) :: file,name,dimname
  call write_row(file,name,row)
 end subroutine
 subroutine vector_read_and_scatter(file,row,nlocal,name,address)
  character(*),intent(in) :: file,name
  real(r8),allocatable,intent(out) :: row(:)
  integer,intent(in) :: nlocal
  type(pointer_int32_1d),intent(in) :: address(0:)
  integer :: ncid,varid
  allocate(row(nlocal))
  if(nlocal==0) return
  call check(nf90_open(file,nf90_nowrite,ncid))
  call check(nf90_inq_varid(ncid,name,varid))
  call check(nf90_get_var(ncid,varid,row))
  call check(nf90_close(ncid))
 end subroutine
end module
module MOD_Grid_RiverLakeHistRoute
 use MOD_Precision
 use probe_netcdf
 integer :: dispatches=0
 character(len=256) :: history_file='history.nc'
contains
 subroutine route_hist_write_ucat(vector,varname,longname,units)
  real(r8),intent(in) :: vector(:)
  character(*),intent(in) :: varname
  character(*),optional,intent(in) :: longname,units
  dispatches=dispatches+1
  call write_row(trim(history_file),varname,vector,units,longname)
 end subroutine
end module
"""


DRIVER = r"""
program main
 use sediment_probe
 use MOD_Grid_RiverLakeNetwork
 use MOD_Grid_RiverLakeHistRoute
 use probe_netcdf
 implicit none
 character(len=32) :: mode
 real(r8) :: water(2),area(2),before(2),deep(2,2),saved(2,2,10),peak(2)
 integer :: ncid
 call get_command_argument(1,mode)
 nsed=2;totlyrnum=2;lambda=.4_r8;psedD=2.65_r8;pwatD=1._r8
 nlfp_sed=1;lyrdph=.05_r8;sed_ignore_dph=.05_r8;sed_max_conc=.01_r8;vonKar=.4_r8
 if(mode=='empty_worker' .or. mode=='nonworker') then
  numucat=0
  p_is_worker=mode=='empty_worker'
  call write_sediment_history('unused',1)
  if(dispatches/=21) error stop 'empty rank missed diagnostic collectives'
  print *, 'EMPTY_OK'
  stop
 endif
 call allocate_sediment_vars()
 allocate(setvel(nsed));setvel=0._r8
 water=[100._r8,0._r8];area=1._r8
 if(mode=='settling') then
  sedsto(:,1)=[.2_r8,.4_r8];sedcon(:,1)=sedsto(:,1)/water(1)
  layer(:,1)=.1_r8;susvel(:,1)=[2._r8,4._r8];setvel=1000._r8
  call calc_sediment_exchange(10._r8,water,area)
  call calc_layer_redistribution(area)
  call accumulate_sediment_output(10._r8)
  sed_hist_acctime=10._r8
  call write_sediment_history('unused',1)
  stop
 else if(mode=='other') then
  deep(:,1)=[1._r8,2._r8];deep(:,2)=[.1_r8,.2_r8]
  call apply_bif_sediment_credits(10._r8,water,deep,0._r8*deep,0._r8*deep,0._r8*deep)
  sedsto_protected(:,1)=[2._r8,1._r8]
  deep=sedsto_protected
  call transfer_levee_sediment(10._r8,10._r8,water,[50._r8,0._r8],[50._r8,0._r8],deep,.false.)
  call accumulate_sediment_output(10._r8)
  sed_hist_acctime=10._r8
  call write_sediment_history('unused',1)
  stop
 endif
 sedsto(:,1)=[.1_r8,.2_r8];sedcon(:,1)=sedsto(:,1)/water(1)
 layer(:,1)=5._r8;susvel(:,1)=[.1_r8,.2_r8]
 sedinp(:,1)=[.02_r8,.01_r8];sedinp(:,2)=[.04_r8,.06_r8]
 call calc_sediment_exchange(10._r8,water,area)
 if(abs(sum(a_seddiag(:,1,SED_DIAG_CAP_EXCH))-.2_r8)>1.e-12_r8) error stop 'exchange cap'
 call apply_sediment_input(10._r8,water,area)
 if(abs(sum(a_seddiag(:,1,SED_DIAG_CAP_HILL))-.3_r8)>1.e-12_r8) error stop 'input cap'
 if(sum(a_seddiag(:,2,SED_DIAG_CAP_HILL))/=0._r8) error stop 'dry input counted as cap'
 before=(1._r8-lambda)*(sum(layer,dim=1)+sum(sum(seddep,dim=1),dim=1))
 call calc_layer_redistribution(area)
 deep=(1._r8-lambda)*sum(seddep,dim=2)
 if(maxval(abs(before-(1._r8-lambda)*sum(layer,dim=1)-sum(deep,dim=1)))>1.e-10_r8) &
  error stop 'layer mass'
 call accumulate_sediment_output(10._r8)
 if(maxval(abs(a_seddiag(:,:,SED_DIAG_DEP)/10._r8-deep))>1.e-12_r8) error stop 'bed sampling'
 sed_hist_acctime=10._r8
 if(mode=='save') then
  call sediment_process_history_checkpoint('checkpoint.nc',.true.,-1)
  stop
 else if(mode=='sidecar') then
  saved=a_seddiag;peak=a_ssc_precap_max
  call sediment_history_acc_sidecar('checkpoint.nc',.true.)
  call flush_sediment_history()
  call sediment_history_acc_sidecar('checkpoint.nc',.false.)
  if(maxval(abs(a_seddiag-saved))>1.e-12_r8 .or. any(a_ssc_precap_max/=peak)) error stop 'sidecar values'
  if(sed_diag_acctime/=10._r8 .or. sed_hist_acctime/=10._r8) error stop 'sidecar clocks'
 else if(mode=='restore' .or. mode=='legacy') then
  saved=a_seddiag;peak=a_ssc_precap_max
  call flush_sediment_history()
  sed_hist_acctime=1000._r8
  call check(nf90_open('checkpoint.nc',nf90_nowrite,ncid))
  call sediment_process_history_checkpoint('checkpoint.nc',.false.,ncid)
  call check(nf90_close(ncid))
  if(mode=='restore') then
   if(maxval(abs(a_seddiag-saved))>1.e-12_r8 .or. any(a_ssc_precap_max/=peak)) error stop 'restart values'
   if(sed_diag_acctime/=10._r8) error stop 'restart clock'
  else
   if(any(a_seddiag/=0._r8) .or. any(a_ssc_precap_max/=spval)) error stop 'legacy values'
   if(sed_diag_acctime/=0._r8) error stop 'legacy clock'
  endif
 endif
 susvel=0._r8;sedinp=0._r8
 call calc_sediment_exchange(30._r8,water,area)
 call calc_layer_redistribution(area)
 call accumulate_sediment_output(30._r8)
 sed_hist_acctime=sed_hist_acctime+30._r8
 call write_sediment_history('unused',1)
 if(dispatches/=21) error stop 'history variables'
 call flush_sediment_history()
 if(sed_diag_acctime/=0._r8 .or. any(a_seddiag/=0._r8) .or. &
  any(a_ssc_precap_max/=spval)) error stop 'history flush'
 ! Different configured grain density must affect the SSC conversion.
 psedD=3._r8
 call record_sediment_precap(1,[.1_r8,.2_r8],10._r8)
 if(abs(a_ssc_precap_max(1)-90000._r8)>1.e-9_r8) error stop 'configured density'
 ! Two advection caps in one interval: 2 -> 1 -> .5 m3 suspended solid.
 call flush_sediment_history()
 psedD=2.65_r8;sedsto=0._r8;layer=0._r8;sedcon=0._r8
 sedsto(:,1)=[1._r8,1._r8];sedsto(:,2)=[.1_r8,.2_r8]
 critshearvel=1.e20_r8;shearvel=0._r8;sedinp=0._r8
 call calc_sediment_advection(10._r8,[0._r8,0._r8],[0._r8,0._r8],water,[50._r8,0._r8])
 if(abs(sum(a_seddiag(:,1,SED_DIAG_CAP_ADV))-1.5_r8)>1.e-12_r8) error stop 'two advection caps'
 if(any(a_seddiag(:,2,SED_DIAG_CAP_ADV)/=0._r8)) error stop 'dry advection counted as cap'
 exch_es_raw=0._r8;exch_d_raw=0._r8;exch_es_eff=0._r8;exch_d_eff=0._r8
 call accumulate_sediment_output(10._r8)
 history_file='advection.nc'
 call write_sediment_process_history()
end program
"""


@pytest.fixture(scope="module")
def executable(tmp_path_factory):
    folder = tmp_path_factory.mktemp("sediment_process_history")
    compiler = require_runnable_fortran_compiler(folder)
    source = SOURCE.read_text()
    declarations = source.split("MODULE MOD_Tracer_Particle_Sediment", 1)[1].split("\nCONTAINS", 1)[0]
    declarations = re.sub(r"^\s*(?:PRIVATE|PUBLIC ::.*)$", "", declarations, flags=re.M)
    names = (
        "allocate_sediment_vars", "record_sediment_precap", "record_sediment_cap_deposition",
        "calc_sediment_exchange", "apply_sediment_input", "calc_layer_redistribution",
        "accumulate_sediment_output", "write_sediment_history", "write_sediment_process_history",
        "flush_sediment_history", "sediment_process_history_checkpoint", "validate_sediment_process_history",
        "sediment_restart_value_finite", "try_read_restart_var", "require_sediment_restart_var",
        "calc_sediment_advection", "calc_sediment_advection_one_direction", "limit_reverse_flux",
        "sediment_history_acc_sidecar", "apply_bif_sediment_credits", "transfer_levee_sediment",
    )
    helpers = r"""
 logical function sediment_particle_enabled()
  sediment_particle_enabled=.true.
 end function
 subroutine assert_sediment_mass_balance(context,index,before,after,expected_change,rivsto_i,rivout_i)
  character(*),intent(in) :: context
  integer,intent(in) :: index
  real(r8),intent(in) :: before,after,expected_change
  real(r8),optional,intent(in) :: rivsto_i,rivout_i
  if(abs(after-before-expected_change)>1.e-9_r8) error stop 'sediment mass balance'
 end subroutine
"""
    code = STUBS + "\nmodule sediment_probe" + declarations + "\ncontains\n" + helpers
    code += "\n".join(routine(source, name) for name in names) + "\nend module\n" + DRIVER
    path = folder / "probe.F90"
    path.write_text(code)
    includes, libs = netcdf_fortran_flags()
    built = subprocess.run(
        [compiler, "-cpp", "-ffree-line-length-0", "-fcheck=all", "-ffpe-trap=invalid,zero,overflow",
         *includes, str(path), *libs, "-o", "probe"],
        cwd=folder, text=True, capture_output=True, timeout=120,
    )
    assert built.returncode == 0, built.stdout + built.stderr
    return folder / "probe"


def run(executable, folder, mode="normal", success=True):
    result = subprocess.run([str(executable), mode], cwd=folder, text=True, capture_output=True, timeout=30)
    if success:
        assert result.returncode == 0, result.stdout + result.stderr
    else:
        assert result.returncode != 0
    return result


def read(path):
    with Dataset(path) as ds:
        return {name: np.ma.filled(var[:], -9999.0) for name, var in ds.variables.items()}


def test_exchange_input_caps_and_time_weighted_storage(executable, tmp_path):
    run(executable, tmp_path)
    data = read(tmp_path / "history.nc")
    assert len(data) == 21
    assert data["f_es_raw_1"][0] == pytest.approx(0.0075)
    assert data["f_es_eff_2"][0] == pytest.approx(0.015)
    assert sum(data[f"f_capdep_exch_{k}"][0] for k in (1, 2)) == pytest.approx(0.005)
    assert sum(data[f"f_capdep_hill_{k}"][0] for k in (1, 2)) == pytest.approx(0.0075)
    assert sum(data[f"f_sedsto_{k}"][0] for k in (1, 2)) == pytest.approx(1.0)
    assert data["f_ssc_precap_max"][0] == pytest.approx(34450.0)
    assert data["f_ssc_precap_max"][1] == -9999.0
    assert sum(data[f"f_capdep_hill_{k}"][1] for k in (1, 2)) == 0.0
    with Dataset(tmp_path / "history.nc") as ds:
        assert ds["f_ssc_precap_max"].units == "mg L-1"
        assert ds["f_seddep_solid_1"].units == "m^3"
        assert ds["f_es_raw_1"].units == "m^3/s"
    adv = read(tmp_path / "advection.nc")
    assert sum(adv[f"f_capdep_adv_{k}"][0] for k in (1, 2)) == pytest.approx(0.15)
    assert sum(adv[f"f_capdep_adv_{k}"][1] for k in (1, 2)) == 0.0


def test_checkpoint_continuation_preserves_numerators_and_peak(executable, tmp_path):
    run(executable, tmp_path, "save")
    run(executable, tmp_path, "restore")
    data = read(tmp_path / "history.nc")
    assert data["f_es_raw_1"][0] == pytest.approx(0.0075)
    assert data["f_ssc_precap_max"][0] == pytest.approx(34450.0)
    with Dataset(tmp_path / "checkpoint.nc") as ds:
        assert len(ds.variables) == 22  # 10 diagnostics * 2 classes + peak + clock


def test_old_checkpoint_starts_new_diagnostic_clock_without_dilution(executable, tmp_path):
    with Dataset(tmp_path / "checkpoint.nc", "w") as ds:
        ds.createDimension("ucatch", 2)
    run(executable, tmp_path, "legacy")
    data = read(tmp_path / "history.nc")
    assert data["f_es_raw_1"][0] == 0.0
    assert sum(data[f"f_sedsto_{k}"][0] for k in (1, 2)) == pytest.approx(1.0)
    assert data["f_ssc_precap_max"][0] == pytest.approx(26500.0)


@pytest.mark.parametrize("mode", ["empty_worker", "nonworker"])
def test_empty_ranks_dispatch_allocated_zero_size_vectors(executable, tmp_path, mode):
    assert "EMPTY_OK" in run(executable, tmp_path, mode).stdout


def test_new_checkpoint_with_missing_diagnostic_fails(executable, tmp_path):
    with Dataset(tmp_path / "checkpoint.nc", "w") as ds:
        ds.createDimension("ucatch", 2)
        ds.createVariable("sed_diag_acctime_vec", "f8", ("ucatch",))[:] = 10.0
    result = run(executable, tmp_path, "restore", success=False)
    assert "required variable" in result.stdout and "a_es_raw_1" in result.stdout


def test_stock_limited_exchange_retains_physical_settling(executable, tmp_path):
    run(executable, tmp_path, "settling")
    data = read(tmp_path / "history.nc")
    assert data["f_es_raw_1"][0] == pytest.approx(0.6)
    assert data["f_d_raw_1"][0] == pytest.approx(0.02)
    assert data["f_es_eff_1"][0] == pytest.approx(0.026)
    assert data["f_d_eff_1"][0] == pytest.approx(0.02)
    assert data["f_capdep_exch_1"][0] == 0.0


def test_bif_levee_caps_include_protected_deposition_but_visible_peak(executable, tmp_path):
    run(executable, tmp_path, "other")
    data = read(tmp_path / "history.nc")
    assert sum(data[f"f_capdep_other_{k}"][0] for k in (1, 2)) == pytest.approx(0.45)
    assert sum(data[f"f_capdep_other_{k}"][1] for k in (1, 2)) == 0.0
    assert sum(data[f"f_d_eff_{k}"][0] for k in (1, 2)) == pytest.approx(0.45)
    assert data["f_ssc_precap_max"][0] == pytest.approx(79500.0)


def test_partial_history_sidecar_roundtrip(executable, tmp_path):
    run(executable, tmp_path, "sidecar")
    data = read(tmp_path / "history.nc")
    assert data["f_es_raw_1"][0] == pytest.approx(0.0075)
    assert data["f_ssc_precap_max"][0] == pytest.approx(34450.0)
