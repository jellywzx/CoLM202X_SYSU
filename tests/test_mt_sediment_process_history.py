"""Targeted checks for Material Transport sediment process diagnostics."""

from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "main/MATERIAL/FAMILY_PARTICLE/MOD_MT_Sediment.F90"
NAMELIST = ROOT / "share/MOD_Namelist.F90"


def routine(source: str, name: str) -> str:
    match = re.search(
        rf"^\s*SUBROUTINE {name}\b.*?^\s*END SUBROUTINE {name}\b",
        source,
        re.MULTILINE | re.DOTALL | re.IGNORECASE,
    )
    if match is None:
        raise AssertionError(f"missing routine {name}")
    return match.group()


class SedimentProcessHistoryTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source = SOURCE.read_text(encoding="utf-8")
        cls.namelist = NAMELIST.read_text(encoding="utf-8")

    def test_history_contract_and_window_lifecycle(self):
        for name in (
            "es_raw", "d_raw", "es_eff", "d_eff", "capdep_adv",
            "capdep_exch", "capdep_hill", "sedsto", "seddep_solid",
            "capdep_other",
        ):
            self.assertIn(f"'{name}'", self.source)
        self.assertIn("'f_ssc_precap_max'", self.source)
        self.assertIn("units = 'mg L-1'", self.source)
        self.assertIn("units = 'm^3'", self.source)
        self.assertIn("units = 'm^3/s'", self.source)
        self.assertIn("logical :: seddiag", self.namelist)
        self.assertIn("sync_hist_vars_one (DEF_hist_vars%seddiag", self.namelist)

        grid = routine(self.source, "grid_sediment_calc")
        self.assertLess(
            grid.index("CALL calc_layer_redistribution(bed_area)"),
            grid.index("CALL accumulate_sediment_output(dt_morph)"),
        )
        flush = routine(self.source, "flush_sediment_history")
        drop = routine(self.source, "mt_sediment_drop_hist_window")
        for body in (flush, drop):
            self.assertIn("a_seddiag = 0._r8", body)
            self.assertIn("a_ssc_precap_max = spval", body)
            self.assertIn("sed_diag_acctime = 0._r8", body)

    def test_restart_and_empty_rank_contract(self):
        checkpoint = routine(self.source, "sediment_process_history_checkpoint")
        for field in (
            "sed_diag_acctime_vec", "a_ssc_precap_max",
            "a_' // trim(SED_HISTORY_DIAG_NAMES",
        ):
            self.assertIn(field, checkpoint)
        self.assertIn("IF (.not. found) THEN", checkpoint)
        self.assertIn("allocate(row(0))", checkpoint)
        snapshot = routine(self.source, "mt_sediment_snapshot_hist_acc")
        self.assertIn("sedhistsnap_diag = a_seddiag", snapshot)
        self.assertIn("sedhistsnap_ssc_precap = a_ssc_precap_max", snapshot)
        self.assertIn("sedhistsnap_diag_acctime = sed_diag_acctime", snapshot)

    def test_compiled_dt_units_and_cap_accounting(self):
        compiler = shutil.which("gfortran")
        if compiler is None:
            self.skipTest("gfortran is unavailable")

        extracted = "\n".join(
            routine(self.source, name)
            for name in (
                "accumulate_sediment_output",
                "record_sediment_precap",
                "record_sediment_cap_deposition",
            )
        )
        program = f"""
module MOD_Precision
 integer, parameter :: r8=kind(1.d0)
end module
module MOD_Grid_RiverLakeNetwork
 integer :: numucat=2
end module
module sediment_probe
 use MOD_Precision
 implicit none
 real(r8), parameter :: spval=-1.d36
 logical :: p_is_worker=.true.
 integer :: nsed=2
 real(r8) :: lambda=.4_r8, psedD=2.65_r8
 integer, parameter :: SED_N_HISTORY_DIAG=10
 integer, parameter :: SED_DIAG_ES_RAW=1,SED_DIAG_D_RAW=2
 integer, parameter :: SED_DIAG_ES_EFF=3,SED_DIAG_D_EFF=4
 integer, parameter :: SED_DIAG_CAP_ADV=5,SED_DIAG_CAP_EXCH=6,SED_DIAG_CAP_HILL=7
 integer, parameter :: SED_DIAG_STO=8,SED_DIAG_DEP=9,SED_DIAG_CAP_OTHER=10
 real(r8), allocatable :: sedcon(:,:),sedout(:,:),bedout(:,:),sedinp(:,:),netflw(:,:)
 real(r8), allocatable :: layer(:,:),shearvel(:),netflw_adv_step(:,:)
 real(r8), allocatable :: exch_es_raw(:,:),exch_d_raw(:,:),exch_es_eff(:,:),exch_d_eff(:,:)
 real(r8), allocatable :: exch_d_adv_step(:,:),sedsto(:,:),seddep(:,:,:)
 real(r8), allocatable :: a_sedcon(:,:),a_sedout(:,:),a_bedout(:,:),a_sedinp(:,:)
 real(r8), allocatable :: a_netflw(:,:),a_layer(:,:),a_shearvel(:)
 real(r8), allocatable :: a_seddiag(:,:,:),a_ssc_precap_max(:)
 real(r8) :: sed_diag_acctime=0._r8
contains
 subroutine CoLM_stop(message)
  character(*),intent(in) :: message
  error stop message
 end subroutine
{extracted}
end module
program main
 use sediment_probe
 implicit none
 integer :: n
 n=2
 allocate(sedcon(n,2),sedout(n,2),bedout(n,2),sedinp(n,2),netflw(n,2))
 allocate(layer(n,2),shearvel(2),netflw_adv_step(n,2),exch_es_raw(n,2),exch_d_raw(n,2))
 allocate(exch_es_eff(n,2),exch_d_eff(n,2),exch_d_adv_step(n,2),sedsto(n,2),seddep(n,2,2))
 allocate(a_sedcon(n,2),a_sedout(n,2),a_bedout(n,2),a_sedinp(n,2),a_netflw(n,2))
 allocate(a_layer(n,2),a_shearvel(2),a_seddiag(n,2,10),a_ssc_precap_max(2))
 sedcon=0;sedout=0;bedout=0;sedinp=0;netflw=0;layer=0;shearvel=0;netflw_adv_step=0
 exch_es_raw=0;exch_d_raw=0;exch_es_eff=0;exch_d_eff=0;exch_d_adv_step=0
 sedsto=0;seddep=0;a_sedcon=0;a_sedout=0;a_bedout=0;a_sedinp=0;a_netflw=0
 a_layer=0;a_shearvel=0;a_seddiag=0;a_ssc_precap_max=spval
 exch_es_raw(:,1)=[.1_r8,.2_r8];exch_d_eff(:,1)=[.3_r8,.4_r8]
 exch_d_adv_step(:,1)=[.05_r8,.06_r8];sedsto(:,1)=[2._r8,3._r8]
 seddep(:,1,1)=[4._r8,5._r8];seddep(:,2,1)=[6._r8,7._r8]
 call record_sediment_cap_deposition(1,'adv',[.2_r8,.3_r8])
 call record_sediment_precap(1,[.1_r8,.2_r8],10._r8)
 call accumulate_sediment_output(10._r8)
 if(any(abs(a_seddiag(:,1,SED_DIAG_ES_RAW)-[1._r8,2._r8])>1.d-12)) error stop 1
 if(any(abs(a_seddiag(:,1,SED_DIAG_D_EFF)-[3.5_r8,4.6_r8])>1.d-12)) error stop 2
 if(any(abs(a_seddiag(:,1,SED_DIAG_STO)-[20._r8,30._r8])>1.d-12)) error stop 3
 if(any(abs(a_seddiag(:,1,SED_DIAG_DEP)-[60._r8,72._r8])>1.d-12)) error stop 4
 if(any(abs(a_seddiag(:,1,SED_DIAG_CAP_ADV)-[.2_r8,.3_r8])>1.d-12)) error stop 5
 if(abs(a_ssc_precap_max(1)-79500._r8)>1.d-9) error stop 6
 if(sed_diag_acctime/=10._r8 .or. a_ssc_precap_max(2)/=spval) error stop 7
 print '(A)', 'SED_DIAG_OK'
end program
"""
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            source = folder / "probe.F90"
            executable = folder / "probe"
            source.write_text(program, encoding="utf-8")
            built = subprocess.run(
                [compiler, "-cpp", "-ffree-line-length-0", "-fcheck=all", str(source), "-o", str(executable)],
                universal_newlines=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                cwd=str(folder),
                timeout=60,
            )
            self.assertEqual(built.returncode, 0, built.stdout + built.stderr)
            ran = subprocess.run(
                [str(executable)], universal_newlines=True, stdout=subprocess.PIPE,
                stderr=subprocess.PIPE, timeout=30
            )
            self.assertEqual(ran.returncode, 0, ran.stdout + ran.stderr)
            self.assertIn("SED_DIAG_OK", ran.stdout)


if __name__ == "__main__":
    unittest.main()
