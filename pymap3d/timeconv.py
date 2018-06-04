# Copyright (c) 2014-2018 Michael Hirsch, Ph.D.
from typing import Union
from datetime import datetime
from dateutil.parser import parse


def str2dt(t: Union[str, datetime]):
    """
    Converts times in string or list of strings to datetime(s)

    output: datetime
    """
    if isinstance(t, datetime):
        return t
    elif isinstance(t, str):
        t = parse(t)
    elif t is not None:
        t = [parse(T) for T in t]

    return t
