"""Real NetCDF/MPI round-trip checks for sediment process history diagnostics."""

from pathlib import Path
import os
import re
import shlex
import shutil
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "main/MATERIAL/FAMILY_PARTICLE/MOD_MT_Sediment.F90"


def mpi_fortran_flags():
    """Return runtime-check flags supported by the MPI wrapper's compiler."""
    try:
        version = subprocess.check_output(
            ["mpif90", "--version"], text=True, stderr=subprocess.STDOUT
        )
    except subprocess.CalledProcessError as error:
        version = error.output
    lowered = version.lower()
    if "gnu fortran" in lowered:
        return ["-ffree-line-length-none", "-fcheck=all"], version.splitlines()[0]
    if "intel" in lowered or "ifort" in lowered or "ifx" in lowered:
        return ["-extend-source", "512", "-check", "all"], version.splitlines()[0]
    return [], version.splitlines()[0] if version.splitlines() else "unknown compiler"


def extract(source: str, name: str) -> str:
    pattern = (
        rf"^\s*(?:(?:ELEMENTAL\s+)?(?:LOGICAL\s+)?FUNCTION|SUBROUTINE)\s+{name}\b.*?"
        rf"^\s*END\s+(?:FUNCTION|SUBROUTINE)\s+{name}\b"
    )
    match = re.search(pattern, source, re.MULTILINE | re.DOTALL | re.IGNORECASE)
    if match is None:
        raise AssertionError(f"missing routine {name}")
    return match.group()


class SedimentProcessHistoryMPI(unittest.TestCase):
    def test_real_mpi_netcdf_history_and_restart_roundtrip(self):
        required = ("mpif90", "mpiexec", "nf-config", "nc-config")
        missing = [tool for tool in required if shutil.which(tool) is None]
        if missing:
            self.skipTest("real MPI/NetCDF test requires: " + ", ".join(missing))

        source = SOURCE.read_text(encoding="utf-8")
        routines = "\n".join(
            extract(source, name)
            for name in (
                "try_read_restart_var",
                "require_sediment_restart_var",
                "sediment_restart_value_finite",
                "mt_sediment_snapshot_hist_acc",
                "sediment_process_history_checkpoint",
                "validate_sediment_process_history",
                "write_sediment_process_history",
                "flush_sediment_history",
            )
        )
        program = MPI_PROGRAM.replace("@PRODUCTION_ROUTINES@", routines)

        try:
            flags = shlex.split(subprocess.check_output(["nf-config", "--fflags"], text=True))
            libs = shlex.split(subprocess.check_output(["nf-config", "--flibs"], text=True))
            libs += shlex.split(subprocess.check_output(["nc-config", "--libs"], text=True))
        except subprocess.CalledProcessError as error:
            self.skipTest(f"NetCDF configuration command failed: {error}")
        compiler_flags, compiler_name = mpi_fortran_flags()

        with tempfile.TemporaryDirectory(prefix="mt-seddiag-mpi-") as tmp:
            folder = Path(tmp)
            probe_source = folder / "probe.F90"
            probe = folder / "probe"
            probe_source.write_text(
                "program probe\n use mpi\n use netcdf\n print *, nf90_noerr\nend program\n",
                encoding="utf-8",
            )
            probed = subprocess.run(
                [
                    "mpif90", "-cpp", *compiler_flags, *flags,
                    str(probe_source), *libs, "-o", str(probe),
                ],
                cwd=folder,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                timeout=60,
            )
            if probed.returncode != 0:
                detail = " ".join(probed.stdout.splitlines()[:3])
                self.skipTest(
                    f"{compiler_name} cannot use the NetCDF-Fortran module selected by "
                    f"nf-config: {detail}"
                )

            test_source = folder / "test.F90"
            executable = folder / "test"
            test_source.write_text(program, encoding="utf-8")
            built = subprocess.run(
                [
                    "mpif90", "-cpp", "-DUSEMPI", *compiler_flags,
                    *flags, str(test_source), *libs, "-o", str(executable),
                ],
                cwd=folder,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                timeout=120,
            )
            self.assertEqual(built.returncode, 0, built.stdout)
            environment = dict(os.environ)
            environment.setdefault("OMPI_ALLOW_RUN_AS_ROOT", "1")
            environment.setdefault("OMPI_ALLOW_RUN_AS_ROOT_CONFIRM", "1")
            ran = subprocess.run(
                ["mpiexec", "-n", "4", str(executable)],
                cwd=folder,
                env=environment,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                timeout=60,
            )
            self.assertEqual(ran.returncode, 0, ran.stdout)
            self.assertIn("SED_DIAG_MPI_NETCDF_OK", ran.stdout)


MPI_PROGRAM = r"""
module MOD_Precision
 integer,parameter :: r8=kind(1.d0)
end module

module MOD_DataType
 type pointer_int32_1d
  integer,allocatable :: val(:)
 end type
end module

module MOD_Grid_RiverLakeNetwork
 use MOD_DataType
 integer :: numucat=0,totalnumucat=2
 type(pointer_int32_1d),allocatable :: ucat_data_address(:)
end module

module MOD_Test_IO
 use mpi
 use netcdf
 use MOD_Precision
 implicit none
contains
 subroutine check_nc(status)
  integer,intent(in) :: status
  integer :: ierr
  if(status/=nf90_noerr)then
   print *,trim(nf90_strerror(status))
   call MPI_Abort(MPI_COMM_WORLD,91,ierr)
  endif
 end subroutine

 subroutine gather_write(local,nlocal,ntotal,file,name,dimname,longname,units)
  real(r8),intent(in) :: local(:)
  integer,intent(in) :: nlocal,ntotal
  character(*),intent(in) :: file,name,dimname
  character(*),optional,intent(in) :: longname,units
  integer :: rank,nproc,ierr,ncid,dimid,varid,status
  integer,allocatable :: counts(:),displs(:)
  real(r8),allocatable :: global(:)
  logical :: exists
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nproc,ierr)
  allocate(counts(nproc),displs(nproc));counts=0;displs=0
  call MPI_Gather(nlocal,1,MPI_INTEGER,counts,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  if(rank==0)then
   do status=2,nproc
    displs(status)=displs(status-1)+counts(status-1)
   enddo
   allocate(global(max(ntotal,1)));global=0._r8
  else
   allocate(global(1));global=0._r8
  endif
  call MPI_Gatherv(local,nlocal,MPI_DOUBLE_PRECISION,global,counts,displs, &
       MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  if(rank==0)then
   inquire(file=trim(file),exist=exists)
   if(exists)then
    call check_nc(nf90_open(trim(file),NF90_WRITE,ncid))
    call check_nc(nf90_redef(ncid))
    status=nf90_inq_dimid(ncid,trim(dimname),dimid)
    if(status/=nf90_noerr)call check_nc(nf90_def_dim(ncid,trim(dimname),ntotal,dimid))
   else
    call check_nc(nf90_create(trim(file),NF90_CLOBBER,ncid))
    call check_nc(nf90_def_dim(ncid,trim(dimname),ntotal,dimid))
   endif
   status=nf90_inq_varid(ncid,trim(name),varid)
   if(status/=nf90_noerr)call check_nc(nf90_def_var(ncid,trim(name),NF90_DOUBLE,[dimid],varid))
   if(present(longname))call check_nc(nf90_put_att(ncid,varid,'long_name',trim(longname)))
   if(present(units))call check_nc(nf90_put_att(ncid,varid,'units',trim(units)))
   call check_nc(nf90_enddef(ncid))
   call check_nc(nf90_put_var(ncid,varid,global(1:ntotal)))
   call check_nc(nf90_close(ncid))
  endif
  deallocate(counts,displs,global)
 end subroutine

 subroutine scatter_read(file,name,local,nlocal)
  character(*),intent(in) :: file,name
  real(r8),allocatable,intent(inout) :: local(:)
  integer,intent(in) :: nlocal
  integer :: rank,nproc,ierr,ncid,varid,i,total
  integer,allocatable :: counts(:),displs(:)
  real(r8),allocatable :: global(:)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nproc,ierr)
  allocate(counts(nproc),displs(nproc));counts=0;displs=0
  call MPI_Gather(nlocal,1,MPI_INTEGER,counts,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  total=1
  if(rank==0)then
   do i=2,nproc
    displs(i)=displs(i-1)+counts(i-1)
   enddo
   total=max(sum(counts),1)
  endif
  allocate(global(total));global=0._r8
  if(rank==0)then
   call check_nc(nf90_open(trim(file),NF90_NOWRITE,ncid))
   call check_nc(nf90_inq_varid(ncid,trim(name),varid))
   call check_nc(nf90_get_var(ncid,varid,global(1:sum(counts))))
   call check_nc(nf90_close(ncid))
  endif
  if(allocated(local))deallocate(local)
  allocate(local(nlocal))
  call MPI_Scatterv(global,counts,displs,MPI_DOUBLE_PRECISION,local,nlocal, &
       MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  deallocate(counts,displs,global)
 end subroutine
end module

module MOD_Vector_ReadWrite
 use MOD_Precision
 use MOD_DataType
 use MOD_Test_IO
contains
 subroutine vector_gather_and_write(row,nlocal,ntotal,address,file,name,dimname)
  real(r8),intent(in) :: row(:)
  integer,intent(in) :: nlocal,ntotal
  type(pointer_int32_1d),intent(in) :: address(0:)
  character(*),intent(in) :: file,name,dimname
  call gather_write(row,nlocal,ntotal,file,name,dimname)
 end subroutine
 subroutine vector_read_and_scatter(file,row,nlocal,name,address)
  character(*),intent(in) :: file,name
  real(r8),allocatable,intent(inout) :: row(:)
  integer,intent(in) :: nlocal
  type(pointer_int32_1d),intent(in) :: address(0:)
  call scatter_read(file,name,row,nlocal)
 end subroutine
end module

module MOD_Grid_RiverLakeHistRoute
 use MOD_Precision
 use MOD_Grid_RiverLakeNetwork,only:numucat,totalnumucat
 use MOD_Test_IO
 character(len=256) :: history_file='history.nc'
contains
 subroutine route_hist_write_ucat(row,name,longname,units)
  real(r8),intent(in) :: row(:)
  character(*),intent(in) :: name
  character(*),optional,intent(in) :: longname,units
  character(len=256) :: lname,u
  lname='';u=''
  if(present(longname))lname=longname
  if(present(units))u=units
  call gather_write(row,size(row),totalnumucat,history_file,name,'ucatch',trim(lname),trim(u))
 end subroutine
end module

module sediment_test
 use mpi
 use netcdf
 use MOD_Precision
 implicit none
 integer,parameter :: SED_N_HISTORY_DIAG=10
 integer,parameter :: SED_DIAG_ES_RAW=1,SED_DIAG_D_RAW=2
 integer,parameter :: SED_DIAG_ES_EFF=3,SED_DIAG_D_EFF=4
 integer,parameter :: SED_DIAG_CAP_ADV=5,SED_DIAG_CAP_EXCH=6,SED_DIAG_CAP_HILL=7
 integer,parameter :: SED_DIAG_STO=8,SED_DIAG_DEP=9,SED_DIAG_CAP_OTHER=10
 character(len=16),parameter :: SED_HISTORY_DIAG_NAMES(SED_N_HISTORY_DIAG)= &
  [character(len=16)::'es_raw','d_raw','es_eff','d_eff','capdep_adv', &
   'capdep_exch','capdep_hill','sedsto','seddep_solid','capdep_other']
 character(len=128),parameter :: SED_HISTORY_DIAG_LONGNAMES(SED_N_HISTORY_DIAG)= &
  [character(len=128)::'raw bed entrainment solid-volume flux', &
   'physical settling solid-volume flux','effective bed entrainment solid-volume flux', &
   'effective deposition solid-volume flux','advection cap deposition', &
   'exchange cap deposition','hillslope cap deposition','suspended storage', &
   'below-active-layer storage','other cap deposition']
 real(r8),parameter :: spval=-1.d36
 integer :: nsed=2,p_comm_glb=MPI_COMM_WORLD,p_address_master=0,p_err
 logical :: p_is_worker=.false.,p_is_master=.false.,p_is_io=.false.
 real(r8) :: sed_diag_acctime=0._r8,sed_hist_acctime=0._r8
 logical :: sedhistsnap_valid=.false.
 real(r8) :: sedhistsnap_acctime=0._r8,sedhistsnap_diag_acctime=0._r8
 real(r8),allocatable :: a_seddiag(:,:,:),a_ssc_precap_max(:)
 real(r8),allocatable :: a_sedcon(:,:),a_sedout(:,:),a_bedout(:,:),a_sedinp(:,:)
 real(r8),allocatable :: a_netflw(:,:),a_layer(:,:),a_shearvel(:)
 real(r8),allocatable :: sedhistsnap_con(:,:),sedhistsnap_out(:,:),sedhistsnap_bed(:,:)
 real(r8),allocatable :: sedhistsnap_inp(:,:),sedhistsnap_net(:,:),sedhistsnap_lay(:,:)
 real(r8),allocatable :: sedhistsnap_shear(:),sedhistsnap_diag(:,:,:),sedhistsnap_ssc_precap(:)
contains
 logical function sediment_particle_enabled()
  sediment_particle_enabled=.true.
 end function
 subroutine CoLM_stop(message)
  character(*),optional,intent(in) :: message
  integer :: ierr
  if(present(message))print *,trim(message)
  call MPI_Abort(MPI_COMM_WORLD,92,ierr)
 end subroutine
@PRODUCTION_ROUTINES@
end module

program main
 use mpi
 use netcdf
 use MOD_Precision
 use MOD_Grid_RiverLakeNetwork
 use MOD_Grid_RiverLakeHistRoute
 use MOD_Test_IO,only:check_nc
 use sediment_test
 implicit none
 integer :: rank,ierr,nuc,i,j,k,ncid,nvars,varid
 real(r8),allocatable :: saved(:,:,:),peak(:),values(:)
 character(len=32) :: units
 call MPI_Init(ierr)
 p_err=ierr
 call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
 p_is_master=rank==0;p_is_io=rank==0
 p_is_worker=rank==1.or.rank==2
 numucat=0
 if(rank==1)numucat=2
 nuc=0
 if(p_is_worker)nuc=numucat
 allocate(ucat_data_address(0:0))
 allocate(a_seddiag(nsed,nuc,SED_N_HISTORY_DIAG),a_ssc_precap_max(nuc))
 allocate(a_sedcon(nsed,nuc),a_sedout(nsed,nuc),a_bedout(nsed,nuc),a_sedinp(nsed,nuc))
 allocate(a_netflw(nsed,nuc),a_layer(nsed,nuc),a_shearvel(nuc))
 a_seddiag=0._r8;a_ssc_precap_max=spval
 a_sedcon=1._r8;a_sedout=2._r8;a_bedout=3._r8;a_sedinp=4._r8
 a_netflw=5._r8;a_layer=6._r8;a_shearvel=7._r8
 if(rank==1)then
  do k=1,SED_N_HISTORY_DIAG
   do j=1,numucat
    do i=1,nsed
     a_seddiag(i,j,k)=10._r8*real(100*k+10*i+j,r8)
    enddo
   enddo
  enddo
  a_ssc_precap_max=[111._r8,222._r8]
  sed_diag_acctime=10._r8;sed_hist_acctime=10._r8
 endif
 allocate(saved(nsed,nuc,SED_N_HISTORY_DIAG),peak(nuc))
 saved=a_seddiag;peak=a_ssc_precap_max
 history_file='history.nc'
 call write_sediment_process_history()
 call sediment_process_history_checkpoint('checkpoint.nc',.true.,-1)
 call mt_sediment_snapshot_hist_acc()
 call flush_sediment_history()
 if(p_is_worker)then
  if(any(a_seddiag/=0._r8).or.sed_diag_acctime/=0._r8)call CoLM_stop('flush failed')
  if(any(sedhistsnap_diag/=saved).or.any(sedhistsnap_ssc_precap/=peak)) &
   call CoLM_stop('terminal snapshot failed')
 endif
 ncid=-1
 if(rank==0)call check_nc(nf90_open('checkpoint.nc',NF90_NOWRITE,ncid))
 call sediment_process_history_checkpoint('checkpoint.nc',.false.,ncid)
 if(rank==0)call check_nc(nf90_close(ncid))
 if(p_is_worker)then
  if(any(a_seddiag/=saved).or.any(a_ssc_precap_max/=peak))call CoLM_stop('restart values failed')
 endif
 if(rank==1)then
  a_seddiag=a_seddiag+3._r8*saved
  sed_diag_acctime=40._r8
 endif
 history_file='continued.nc'
 call write_sediment_process_history()
 call MPI_Barrier(MPI_COMM_WORLD,ierr)
 if(rank==0)then
  call check_nc(nf90_open('history.nc',NF90_NOWRITE,ncid))
  call check_nc(nf90_inquire(ncid,nVariables=nvars))
  if(nvars/=21)call CoLM_stop('history variable count failed')
  call check_nc(nf90_inq_varid(ncid,'f_es_raw_1',varid))
  units='';call check_nc(nf90_get_att(ncid,varid,'units',units))
  if(trim(units)/='m^3/s')call CoLM_stop('flux units failed')
  allocate(values(2));call check_nc(nf90_get_var(ncid,varid,values))
  if(any(abs(values-[111._r8,112._r8])>1.d-12))call CoLM_stop('history values failed')
  deallocate(values)
  call check_nc(nf90_inq_varid(ncid,'f_sedsto_1',varid))
  units='';call check_nc(nf90_get_att(ncid,varid,'units',units))
  if(trim(units)/='m^3')call CoLM_stop('storage units failed')
  call check_nc(nf90_inq_varid(ncid,'f_ssc_precap_max',varid))
  units='';call check_nc(nf90_get_att(ncid,varid,'units',units))
  if(trim(units)/='mg L-1')call CoLM_stop('peak units failed')
  call check_nc(nf90_close(ncid))
  call check_nc(nf90_open('continued.nc',NF90_NOWRITE,ncid))
  call check_nc(nf90_inq_varid(ncid,'f_es_raw_1',varid))
  allocate(values(2));call check_nc(nf90_get_var(ncid,varid,values))
  if(any(abs(values-[111._r8,112._r8])>1.d-12))call CoLM_stop('restart continuation failed')
  deallocate(values)
  call check_nc(nf90_close(ncid))
  call check_nc(nf90_open('checkpoint.nc',NF90_NOWRITE,ncid))
  call check_nc(nf90_inquire(ncid,nVariables=nvars))
  if(nvars/=22)call CoLM_stop('checkpoint variable count failed')
  call check_nc(nf90_close(ncid))
  call check_nc(nf90_create('legacy.nc',NF90_CLOBBER,ncid))
  call check_nc(nf90_def_dim(ncid,'ucatch',2,varid))
  call check_nc(nf90_close(ncid))
 endif
 call MPI_Barrier(MPI_COMM_WORLD,ierr)
 sed_hist_acctime=123._r8
 ncid=-1
 if(rank==0)call check_nc(nf90_open('legacy.nc',NF90_NOWRITE,ncid))
 call sediment_process_history_checkpoint('legacy.nc',.false.,ncid)
 if(rank==0)call check_nc(nf90_close(ncid))
 if(p_is_worker)then
  if(any(a_seddiag/=0._r8).or.any(a_ssc_precap_max/=spval))call CoLM_stop('legacy reset failed')
  if(sed_diag_acctime/=0._r8.or.sed_hist_acctime/=123._r8)call CoLM_stop('legacy clocks failed')
 endif
 if(rank==0)print *,'SED_DIAG_MPI_NETCDF_OK'
 call MPI_Finalize(ierr)
end program
"""


if __name__ == "__main__":
    unittest.main()
