# Copyright (c) 2014-2018 Michael Hirsch, Ph.D.
""" convert strings to datetime """
from datetime import datetime
from dateutil.parser import parse
import numpy as np

__all__ = ['str2dt']


def str2dt(time: datetime) -> np.ndarray:
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
        return parse(time)
    elif isinstance(time, np.datetime64):
        return time.astype(datetime)
    else:  # some sort of iterable
        try:
            if isinstance(time[0], datetime):
                return time
            elif isinstance(time[0], np.datetime64):
                return time.astype(datetime)
            elif isinstance(time[0], str):
                return [parse(t) for t in time]
        except (IndexError, TypeError):
            pass

        # last resort--assume pandas/xarray

        return time.values.astype('datetime64[us]').astype(datetime)
