import pytest
from pytest import approx
from math import radians
from datetime import datetime

import pymap3d.sidereal as pmd
import pymap3d.haversine as pmh

lon = -148
t0 = datetime(2014, 4, 6, 8)
sra = 2.90658
ha = 45.482789587392013


@pytest.mark.parametrize(
    "time, use_astropy", [(t0, False), (t0, True), ([t0], False), ([t0], True)]
)
def test_sidereal(time, use_astropy):
    if use_astropy:
        pytest.importorskip("astropy")
    # http://www.jgiesen.de/astro/astroJS/siderealClock/
    tsr = pmd.datetime2sidereal(time, radians(lon), use_astropy=use_astropy)
    if isinstance(tsr, list):
        tsr = tsr[0]
    assert tsr == approx(sra, rel=1e-5)


def test_anglesep():
    pytest.importorskip("astropy")
    assert pmh.anglesep(35, 23, 84, 20) == approx(ha)


def test_anglesep_meeus():
    # %% compare with astropy
    assert pmh.anglesep_meeus(35, 23, 84, 20) == approx(ha)
