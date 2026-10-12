"""Keep generated interception/thermal build switches consistent with TRACER.

The extended modules expose tracer phase-change outputs and include optional
flood-evaporation arguments only with CANOPY_FLOOD_FEEDBACK.  Test actual
create_defineh.bash output rather than just searching source text.
"""

from pathlib import Path
import shutil
import subprocess

import pytest


ROOT = Path(__file__).resolve().parents[1]
GENERATOR = ROOT / ".github" / "workflows" / "create_defineh.bash"


@pytest.mark.parametrize(
    ("domain", "cama", "tracer", "expected_impl", "expected_flood"),
    [
        ("GRID", "CaMaOFF", "TRACERON", "BUILD_EXTENDED", "FLOOD_ON"),
        ("GRID", "CaMaON", "TRACERON", "BUILD_EXTENDED", "FLOOD_ON"),
        ("GRID", "CaMaOFF", "TRACEROFF", "BUILD_STANDARD", "FLOOD_OFF"),
        ("GRID", "CaMaON", "TRACEROFF", "BUILD_STANDARD", "FLOOD_OFF"),
        ("SinglePoint", "CaMaOFF", "TRACERON", "BUILD_EXTENDED", "FLOOD_OFF"),
        ("SinglePoint", "CaMaOFF", "TRACEROFF", "BUILD_STANDARD", "FLOOD_OFF"),
    ],
)
def test_generated_interception_and_flood_switches(
    tmp_path: Path, domain: str, cama: str, tracer: str,
    expected_impl: str, expected_flood: str,
):
    if shutil.which("cpp") is None:
        pytest.skip("requires the C preprocessor used by the Makefile")
    (tmp_path / "include").mkdir()
    subprocess.run(
        ["bash", str(GENERATOR), domain, "LULC_IGBP_PC", "URBANOFF",
         "vanGenu", cama, "BGCOFF", "CROPOFF", tracer],
        cwd=tmp_path,
        check=True,
        capture_output=True,
        text=True,
    )

    probe = """#include <define.h>
#ifdef TRACER
TRACER_ON
#else
TRACER_OFF
#endif
#ifdef extend_interception
BUILD_EXTENDED
#else
BUILD_STANDARD
#endif
#ifdef CANOPY_FLOOD_FEEDBACK
FLOOD_ON
#else
FLOOD_OFF
#endif
"""
    proc = subprocess.run(
        ["cpp", "-P", "-I", str(tmp_path / "include"), "-"],
        input=probe,
        text=True,
        capture_output=True,
        check=True,
    )
    flags = {line.strip() for line in proc.stdout.splitlines()}
    assert ("TRACER_ON" in flags) == (tracer == "TRACERON")
    assert expected_impl in flags
    assert expected_flood in flags
