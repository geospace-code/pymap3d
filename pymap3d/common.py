from dateutil.parser import parse
from datetime import datetime
from numpy import atleast_1d
from six import string_types

def str2dt(t):
    t = atleast_1d(t)
    if isinstance(t[0],string_types):
        t = [parse(T) for T in t]

    assert isinstance(t[0],datetime)

    return t