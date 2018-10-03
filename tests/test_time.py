#!/usr/bin/env python
import pytest
from pytest import approx
from pymap3d.timeconv import str2dt
import pymap3d.sidereal as pms
from datetime import datetime
import numpy as np

t0 = datetime(2014, 4, 6, 8)


def test_juliantime():
    assert pms.juliandate(t0) == approx(2.456753833333e6)


def test_types():
    assert str2dt(t0) == t0  # passthrough
    assert str2dt('2014-04-06T08:00:00') == t0
    ti = [str2dt('2014-04-06T08:00:00'), str2dt('2014-04-06T08:01:02')]
    to = [t0, datetime(2014, 4, 6, 8, 1, 2)]
    assert ti == to   # even though ti is numpy array of datetime and to is list of datetime

    t1 = [t0, t0]
    assert (np.asarray(str2dt(t1)) == t0).all()


def test_datetime64():
    t1 = np.datetime64(t0)
    assert str2dt(t1) == t0

    t1 = np.array([np.datetime64(t0), np.datetime64(t0)])
    assert (str2dt(t1) == t0).all()


def test_xarray_time():
    xarray = pytest.importorskip('xarray')

    t = {'time': t0}
    ds = xarray.Dataset(t)
    assert str2dt(ds['time']) == t0

    t2 = {'time': [t0, t0]}
    ds = xarray.Dataset(t2)
    assert (str2dt(ds['time']) == t0).all()


def test_pandas_time():
    pandas = pytest.importorskip('pandas')

    t = pandas.Series(t0)
    assert (str2dt(t) == t0).all()

    t = pandas.Series([t0, t0])
    assert (str2dt(t) == t0).all()


if __name__ == '__main__':
    pytest.main(['-x', __file__])
