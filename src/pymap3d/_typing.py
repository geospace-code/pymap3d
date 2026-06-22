"""Typing helpers for pymap3d.

NumPy is strictly optional at runtime.
"""

from __future__ import annotations

from datetime import datetime
from typing import TYPE_CHECKING, Any, Union

# typing.Union is for Python 3.9. Python 3.10+ can use the | operator for unions.

if TYPE_CHECKING:
    # Only loaded by type checkers (mypy, pyright, etc.)
    import numpy as np
    import numpy.typing as npt

    ArrayLike = npt.ArrayLike
    NDArray = npt.NDArray

    DatetimeLike = Union[str, datetime, np.datetime64]
    FloatLike = Union[float, np.floating]

    FloatArray = npt.NDArray[np.floating]
    DatetimeArray = Union[npt.NDArray[np.datetime64], np.datetime64, datetime]
else:
    ArrayLike = Any
    NDArray = Any

    DatetimeLike = Union[str, datetime]
    FloatLike = float

    FloatArray = float
    DatetimeArray = datetime


__all__ = ["ArrayLike", "NDArray", "FloatArray", "FloatLike", "DatetimeArray", "DatetimeLike"]
