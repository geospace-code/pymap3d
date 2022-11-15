# Copyright (c) 2014-2018 Michael Hirsch, Ph.D.
""" convert strings to datetime """

from __future__ import annotations

from datetime import datetime
from typing import Any, List, cast, overload

try:
    from numpy import datetime64
    from numpy.typing import NDArray
except ImportError:
    pass

try:
    import dateutil.parser
except ImportError:
    pass

__all__ = ["str2dt"]


@overload
def str2dt(time: str | datetime | datetime64) -> datetime:
    pass


@overload
def str2dt(time: list[str] | list[datetime]) -> list[datetime]:
    pass


@overload
def str2dt(time: NDArray[datetime64]) -> NDArray[Any]:
    pass


def str2dt(
    time: str | datetime | datetime64 | list[str] | list[datetime] | NDArray[datetime64],
) -> datetime | list[datetime] | NDArray[Any]:
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
    if isinstance(time, list):
        try:
            if isinstance(time[0], datetime):
                return cast(List[datetime], time)
            elif isinstance(time[0], str):
                return [dateutil.parser.parse(cast(str, t)) for t in time]
        except NameError:
            raise ImportError("pip install dateutil")

    # pandas/xarray
    try:
        return time.values.astype("datetime64[us]").astype(datetime)  # type: ignore[no-any-return, union-attr]
    except AttributeError:
        pass

    # Numpy.datetime64
    try:
        return time.astype(datetime)
    except AttributeError:
        pass

    return time  # type: ignore[return-value]
