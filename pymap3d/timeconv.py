import numpy as np
from six import string_types
from dateutil.parser import parse
from datetime import datetime

def str2dt(t):
    """
    output: datetime
    """

    t = np.atleast_1d(t)
    if isinstance(t[0], string_types):
        t = np.array([parse(T) for T in t])

    assert isinstance(t[0], datetime), 'did not convert {} to datetime'.format(type(t[0]))

    return t[()]   # [()] to pop scalar from 0d array while being OK with ndim>0