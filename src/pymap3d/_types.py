from typing import Any, Sequence, Union

try:
    from numpy.typing import NDArray
except ImportError:
    pass

ArrayLike = Union[Sequence[float], "NDArray[Any]"]
