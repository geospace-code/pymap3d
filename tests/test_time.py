#!/usr/bin/env python
import pytest
from pymap3d.timeconv import str2dt
from datetime import datetime

t0 = datetime(2014, 4, 6, 8)


def test_str2dt():
    assert str2dt(t0) == t0  # passthrough
    assert str2dt('2014-04-06T08:00:00') == t0
    ti = [str2dt('2014-04-06T08:00:00'), str2dt('2014-04-06T08:01:02')]
    to = [t0, datetime(2014, 4, 6, 8, 1, 2)]
    assert ti == to   # even though ti is numpy array of datetime and to is list of datetime


def test_xarray_time():
    xarray = pytest.importorskip('xarray')

    t = {'time': t0}

    ds = xarray.Dataset(t)
    assert str2dt(ds['time']) == t0


def test_pandas_time():
    pandas = pytest.importorskip('pandas')

    t = pandas.Series(t0)
    assert str2dt(t) == t0


if __name__ == '__main__':
    pytest.main([__file__])
