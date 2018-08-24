# Copyright (c) 2014-2018 Michael Hirsch, Ph.D.
from typing import Union, List
from datetime import datetime
from dateutil.parser import parse
import numpy as np


def str2dt(time: Union[str, datetime, list, tuple, np.ndarray]) -> Union[datetime, List[datetime], np.ndarray]:
    """
    Converts times in string or list of strings to datetime(s)

    output: datetime
    """
    if isinstance(time, (float, int)) or (isinstance(time, (tuple, list, np.ndarray)) and isinstance(time[0], (float, int))):
        return time  # assuming Unix epoch time or radians

    if isinstance(time, datetime):
        return time
    elif isinstance(time, str):
        time = parse(time)
    elif time is not None:  # iterable
        time = [parse(t) for t in time]
    else:
        raise TypeError('unknown time spec type {}'.format(type(time)))

    return time
