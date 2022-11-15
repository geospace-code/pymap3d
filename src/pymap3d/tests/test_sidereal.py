from datetime import datetime
from math import radians

import pymap3d.haversine as pmh
import pymap3d.sidereal as pmd
import pytest
from pytest import approx

lon = -148
t0 = datetime(2014, 4, 6, 8)
sra = 2.90658
ha = 45.482789587392013


@pytest.mark.parametrize("time", [t0, [t0]])
def test_sidereal(time: datetime) -> None:
    # http://www.jgiesen.de/astro/astroJS/siderealClock/
    tsr = pmd.datetime2sidereal(time, radians(lon))
    if isinstance(tsr, list):
        tsr = tsr[0]
    assert tsr == approx(sra, rel=1e-5)


def test_anglesep() -> None:
    pytest.importorskip("astropy")
    assert pmh.anglesep(35, 23, 84, 20) == approx(ha)


def test_anglesep_list() -> None:
    # %% compare with astropy
    np = pytest.importorskip("numpy")
    ha1 = pmh.anglesep([radians(35)], [radians(23)], [radians(84)], [radians(20)], False)
    assert np.isclose(ha1, radians(ha)).all()


def test_anglesep_meeus() -> None:
    # %% compare with astropy
    assert pmh.anglesep_meeus(35, 23, 84, 20) == approx(ha)


def test_anglesep_meeus_list() -> None:
    # %% compare with astropy
    np = pytest.importorskip("numpy")
    ha1 = pmh.anglesep_meeus([radians(35)], [radians(23)], [radians(84)], [radians(20)], False)
    assert np.isclose(ha1, radians(ha)).all()
