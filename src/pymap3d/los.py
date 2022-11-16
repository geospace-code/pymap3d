""" Line of sight intersection of space observer to ellipsoid """

from __future__ import annotations

from math import nan, pi
from typing import Any, Sequence, overload

try:
    from numpy import asarray
    from numpy.typing import NDArray
except ImportError:
    pass

from ._types import ArrayLike
from .ecef import ecef2geodetic, enu2uvw, geodetic2ecef
from .ellipsoid import Ellipsoid
from .enu import aer2enu
from .mathfun import sqrt

__all__ = ["lookAtSpheroid"]


@overload
def lookAtSpheroid(
    lat0: float,
    lon0: float,
    h0: float,
    az: float,
    tilt: float,
    ell: Ellipsoid | None = None,
    deg: bool = True,
) -> tuple[NDArray[Any], NDArray[Any], NDArray[Any]]:
    pass


@overload
def lookAtSpheroid(
    lat0: ArrayLike,
    lon0: ArrayLike,
    h0: ArrayLike,
    az: ArrayLike,
    tilt: ArrayLike,
    ell: Ellipsoid | None = None,
    deg: bool = True,
) -> tuple[NDArray[Any], NDArray[Any], NDArray[Any]]:
    pass


def lookAtSpheroid(
    lat0: float | ArrayLike,
    lon0: float | ArrayLike,
    h0: float | ArrayLike,
    az: float | ArrayLike,
    tilt: float | ArrayLike,
    ell: Ellipsoid | None = None,
    deg: bool = True,
) -> tuple[float, float, float] | tuple[NDArray[Any], NDArray[Any], NDArray[Any]]:
    """
    Calculates line-of-sight intersection with Earth (or other ellipsoid) surface from above surface / orbit

    Parameters
    ----------

    lat0 : float
           observer geodetic latitude
    lon0 : float
           observer geodetic longitude
    h0 : float
        observer altitude (meters)  Must be non-negative since this function doesn't consider terrain
    az : float
        azimuth angle of line-of-sight, clockwise from North
    tilt : float
        tilt angle of line-of-sight with respect to local vertical (nadir = 0)
    ell : Ellipsoid, optional
          reference ellipsoid
    deg : bool, optional
          degrees input/output  (False: radians in/out)

    Results
    -------

    lat : float
           geodetic latitude where the line-of-sight intersects with the Earth ellipsoid
    lon : float
           geodetic longitude where the line-of-sight intersects with the Earth ellipsoid
    d : float
        slant range (meters) from starting point to intersect point

    Values will be NaN if the line of sight does not intersect.

    Algorithm based on https://medium.com/@stephenhartzell/satellite-line-of-sight-intersection-with-earth-d786b4a6a9b6 Stephen Hartzell
    """

    if ell is None:
        ell = Ellipsoid.from_name("wgs84")

    try:
        lat0 = asarray(lat0)
        lon0 = asarray(lon0)
        h0 = asarray(h0)
        az = asarray(az)
        tilt = asarray(tilt)
        if (h0 < 0).any():
            raise ValueError("Intersection calculation requires altitude  [0, Infinity)")
    except NameError:
        if h0 < 0:  # type: ignore[operator]
            raise ValueError("Intersection calculation requires altitude  [0, Infinity)")

    assert not isinstance(tilt, Sequence)

    a = ell.semimajor_axis
    b = ell.semimajor_axis
    c = ell.semiminor_axis

    el = tilt - 90.0 if deg else tilt - pi / 2

    e, n, u = aer2enu(az, el, srange=1.0, deg=deg)  # type: ignore[arg-type]
    # fixed 1 km slant range

    u, v, w = enu2uvw(e, n, u, lat0, lon0, deg=deg)  # type: ignore[arg-type]
    x, y, z = geodetic2ecef(lat0, lon0, h0, deg=deg)  # type: ignore[arg-type]

    value = -(a**2) * b**2 * w * z - a**2 * c**2 * v * y - b**2 * c**2 * u * x
    radical = (
        a**2 * b**2 * w**2
        + a**2 * c**2 * v**2
        - a**2 * v**2 * z**2
        + 2 * a**2 * v * w * y * z
        - a**2 * w**2 * y**2
        + b**2 * c**2 * u**2
        - b**2 * u**2 * z**2
        + 2 * b**2 * u * w * x * z
        - b**2 * w**2 * x**2
        - c**2 * u**2 * y**2
        + 2 * c**2 * u * v * x * y
        - c**2 * v**2 * x**2
    )

    magnitude = a**2 * b**2 * w**2 + a**2 * c**2 * v**2 + b**2 * c**2 * u**2

    # %%   Return nan if radical < 0 or d < 0 because LOS vector does not point towards Earth
    try:
        radical[radical < 0] = nan  # type: ignore[index]
    except TypeError:
        if radical < 0:
            radical = nan

    d = (value - a * b * c * sqrt(radical)) / magnitude

    try:
        d[d < 0] = nan
    except TypeError:
        if d < 0:
            d = nan

    # %% cartesian to ellipsodal
    lat, lon, _ = ecef2geodetic(x + d * u, y + d * v, z + d * w, deg=deg)

    try:
        return lat.squeeze()[()], lon.squeeze()[()], d.squeeze()[()]
    except AttributeError:
        return lat, lon, d
