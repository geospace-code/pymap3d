# Copyright (c) 2014-2018 Michael Hirsch, Ph.D.
""" convert strings to datetime """
from datetime import datetime

try:
    from dateutil.parser import parse
except ImportError:
    parse = None
try:
    import numpy as np
except ImportError:
    np = None


__all__ = ["str2dt"]


def str2dt(time: datetime) -> datetime:
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
        if parse is None:
            raise TypeError("expected datetime")
        return parse(time)
    elif np is not None and isinstance(time, np.datetime64):
        return time.astype(datetime)
    else:  # some sort of iterable
        try:
            if isinstance(time[0], datetime):
                return time
            elif np is not None and isinstance(time[0], np.datetime64):
                return time.astype(datetime)
            elif isinstance(time[0], str):
                if parse is None:
                    raise TypeError("expected datetime")
                return [parse(t) for t in time]
        except (IndexError, TypeError):
            pass

        # last resort--assume pandas/xarray

        return time.values.astype("datetime64[us]").astype(datetime)
