from __future__ import annotations

from datetime import datetime
from importlib.util import find_spec
from typing import cast

try:
    from numpy import datetime64
    from numpy.typing import NDArray
except ImportError:
    pass

import pymap3d.sidereal as pms
import pytest
from pymap3d.timeconv import str2dt
from pytest import approx

t0 = datetime(2014, 4, 6, 8)
t1 = datetime(2014, 4, 6, 8, 1, 2)
t0_str = "2014-04-06T08:00:00"
t1_str = "2014-04-06T08:01:02"


def test_juliantime() -> None:
    assert pms.juliandate(t0) == approx(2.456753833333e6)


def test_types() -> None:
    np = pytest.importorskip("numpy")
    assert str2dt(t0) == t0  # passthrough

    if find_spec("dateutil"):
        assert str2dt(t0_str) == t0
    else:
        with pytest.raises(ImportError) as excinfo:
            str2dt(t0_str)
        assert str(excinfo.value) == "pip install python-dateutil"

    if find_spec("dateutil"):
        ti = str2dt([t0_str, t1_str])
        to = [t0, t1]
        assert ti == to
    else:
        with pytest.raises(ImportError) as excinfo:
            str2dt([t0_str, t1_str])
        assert str(excinfo.value) == "pip install python-dateutil"

    t2 = [t0, t0]
    assert (np.asarray(str2dt(t2)) == t0).all()


def test_datetime64() -> None:
    np = pytest.importorskip("numpy")
    t1: datetime64 | NDArray[datetime64]
    t1 = cast(datetime64, np.datetime64(t0))
    assert str2dt(t1) == t0

    t1 = cast(NDArray[datetime64], np.array([np.datetime64(t0), np.datetime64(t0)]))
    assert (str2dt(t1) == t0).all()


def test_xarray_time() -> None:
    xarray = pytest.importorskip("xarray")

    t = {"time": t0}
    ds = xarray.Dataset(t)
    assert str2dt(ds["time"]) == t0

    t2 = {"time": [t0, t0]}
    ds = xarray.Dataset(t2)
    assert (str2dt(ds["time"]) == t0).all()


def test_pandas_time() -> None:
    pandas = pytest.importorskip("pandas")

    t = pandas.Series(t0)
    assert (str2dt(t) == t0).all()

    t = pandas.Series([t0, t0])
    assert (str2dt(t) == t0).all()
