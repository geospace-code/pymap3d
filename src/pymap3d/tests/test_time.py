from datetime import datetime

import pymap3d.sidereal as pms
import pytest
from pymap3d.timeconv import str2dt
from pytest import approx

t0 = datetime(2014, 4, 6, 8)


def test_juliantime():
    assert pms.juliandate(t0) == approx(2.456753833333e6)


def test_types():
    np = pytest.importorskip("numpy")
    assert str2dt(t0) == t0  # passthrough
    assert str2dt("2014-04-06T08:00:00") == t0
    ti = str2dt("2014-04-06T08:00:00")
    assert ti == t0

    assert (np.asarray(str2dt(t0)) == t0).all()


def test_datetime64():
    np = pytest.importorskip("numpy")
    t1 = np.datetime64(t0)
    assert str2dt(t1) == t0

    t1 = np.array([np.datetime64(t0), np.datetime64(t0)])
    assert (str2dt(t1) == t0).all()
