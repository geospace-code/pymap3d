# Copyright (c) 2014-2018 Michael Hirsch, Ph.D.
""" convert strings to datetime """

from __future__ import annotations

from datetime import datetime

try:
    import dateutil.parser
except ImportError:
    pass

__all__ = ["str2dt"]


def str2dt(time: str | datetime) -> datetime:
    """
    Converts times in string or list of strings to datetime(s)

    Parameters
    ----------

    time : str or datetime.datetime or numpy.datetime64

    Results
    -------

    t : datetime.datetime

    """
    if isinstance(time, datetime):
        return time
    elif isinstance(time, str):
        try:
            return dateutil.parser.parse(time)
        except NameError:
            raise ImportError("pip install dateutil")

    # some sort of iterable
    try:
        if isinstance(time[0], datetime):
            return time
        elif isinstance(time[0], str):
            return [dateutil.parser.parse(t) for t in time]
    except IndexError:
        pass
    except NameError:
        raise ImportError("pip install dateutil")

    # pandas/xarray
    try:
        return time.values.astype("datetime64[us]").astype(datetime)
    except AttributeError:
        pass

    # Numpy.datetime64
    try:
        return time.astype(datetime)
    except AttributeError:
        pass

    return time
