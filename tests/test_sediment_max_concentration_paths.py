"""Keep every river sediment concentration cap tied to max_concentration.

The sediment bed/suspended balance is managed by the existing process routines.
These assertions guard against accidentally re-introducing a fixed 0.01 cap
into a subset of the advection, levee or bifurcation pathways.
"""
from pathlib import Path
import re

SOURCE = (
    Path(__file__).resolve().parents[1]
    / "main/TRACER/MOD_Tracer_Particle_Sediment.F90"
).read_text(encoding="utf-8")


def subroutine(name: str) -> str:
    start = re.search(rf"\bSUBROUTINE\s+{name}\s*\(", SOURCE, re.I)
    assert start is not None, f"{name} missing"
    end = re.search(rf"\bEND\s+SUBROUTINE\s+{name}\b", SOURCE[start.end():], re.I)
    assert end is not None, f"{name} unterminated"
    return SOURCE[start.start(): start.end() + end.end()]


def test_default_override_validation_and_restart_use_one_limit():
    assert "SED_DEFAULT_MAX_CONC = 0.01_r8" in SOURCE
    assert "sed_max_conc = SED_DEFAULT_MAX_CONC" in SOURCE
    assert "DEF_SEDIMENT%max_concentration >= 0._r8" in SOURCE
    assert "sed_max_conc <= 0._r8 .or. sed_max_conc > 1._r8" in SOURCE
    assert "'sed_max_conc_meta'" in subroutine("read_sediment_restart")
    assert "'sed_max_conc_meta'" in subroutine("write_sediment_restart")


def test_all_cap_paths_use_configured_ceiling():
    assert not re.search(r"\bMAX_SED_CONC\b", SOURCE), (
        "legacy hard-coded maximum is hiding the configured limit"
    )
    # Uses include comparisons and the corresponding deposited excess.
    expected = {
        "transfer_levee_sediment": 2,
        "apply_bif_sediment_credits": 1,
        "calc_sediment_advection": 4,
        "apply_sediment_input": 2,
        "calc_sediment_exchange": 2,
    }
    for name, count in expected.items():
        body = subroutine(name)
        uses = re.findall(r"\*\s*sed_max_conc\b", body)
        assert len(uses) == count, (name, len(uses), count)
        assert "record_sediment_cap_deposition" in body, name
        assert "/ (1._r8 - lambda)" in body, name


def test_near_dry_filter_diagnostics_have_correct_type_and_grid_scope():
    calc = subroutine("grid_sediment_calc")
    assert re.search(
        r"real\(r8\)\s*::\s*carrier_filter_diag_global\(3\)", calc, re.I
    )
    assert "carrier_filter_count_global(1)" in calc
    assert "USE MOD_Grid_RiverLakeNetwork, only: x_ucat, y_ucat, griducat" in calc
    assert "CALL mpi_allreduce(MPI_IN_PLACE, carrier_filter_diag_global" in calc
    assert "Sediment carrier filter diag: inst_near_dry_abs_q=" in calc
    assert all(
        f"carrier_filter_diag_global({i})" in calc for i in (1, 2, 3)
    )
