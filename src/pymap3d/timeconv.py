"""convert strings to datetime"""

from __future__ import annotations

from datetime import datetime

from ._typing import DatetimeArray, DatetimeLike


__all__ = ["str2dt"]


def str2dt(time: DatetimeLike) -> DatetimeArray:
    """
    Converts times in string or list of strings to datetime(s)

    Parameters
    ----------

    time : str or datetime.datetime or numpy.datetime64

    Results
    -------

    t : datetime.datetime or numpy.datetime64
    """

    if isinstance(time, datetime):
        return time
    elif isinstance(time, str):
        return datetime.fromisoformat(time)
    else:
        return time.astype(datetime)
