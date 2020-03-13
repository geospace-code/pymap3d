""" transforms involving ECI earth-centered inertial """

from datetime import datetime
import numpy as np
import typing

from .sidereal import datetime2sidereal

try:
    from astropy.time import Time
except ImportError:
    Time = None

__all__ = ["eci2ecef", "ecef2eci"]

if typing.TYPE_CHECKING:
    from numpy import ndarray


def eci2ecef(
    x: "ndarray", y: "ndarray", z: "ndarray" = None, time: datetime = None, *, useastropy: bool = True
) -> typing.Tuple["ndarray", "ndarray", "ndarray"]:
    """
    Observer => Point  ECI  =>  ECEF

    defaults to J2000 frame

    Parameters
    ----------
    x : "ndarray"
        ECI x-location [meters]
    y : "ndarray"
        ECI y-location [meters]
    z : "ndarray"
        ECI z-location [meters]
    time : datetime.datetime
        time of obsevation (UTC)
    useastropy : bool, optional
        use AstroPy for conversion

    Results
    -------
    x : "ndarray"
        target x ECEF coordinate
    y : "ndarray"
        target y ECEF coordinate
    z : "ndarray"
        target z ECEF coordinate
    """
    # %%
    x = np.atleast_1d(x)
    y = np.atleast_1d(y)
    z = np.atleast_1d(z)
    if not x.shape == y.shape == z.shape:
        raise ValueError("shapes of ECI x,y,z must be identical")

    useastropy = useastropy and Time is not None

    if useastropy:
        gst = Time(time).sidereal_time("apparent", "greenwich").radian
    else:
        gst = datetime2sidereal(time, 0.0)

    gst = np.atleast_1d(gst)
    if gst.ndim != 1:
        raise ValueError("GST must be scalar or vector in radians")
    if gst.size == 1:
        gst *= np.ones(x.size)
    if gst.size != x.size:
        raise ValueError("GST must be scalar or same length as positions")

    eci = np.column_stack((x.ravel(), y.ravel(), z.ravel()))
    ecef = np.empty((x.size, 3))
    for i in range(eci.shape[0]):
        ecef[i, :] = _rottrip(gst[i]) @ eci[i, :]

    xecef = ecef[:, 0].reshape(x.shape)
    yecef = ecef[:, 1].reshape(x.shape)
    zecef = ecef[:, 2].reshape(x.shape)

    return xecef, yecef, zecef


def ecef2eci(
    x: "ndarray", y: "ndarray", z: "ndarray" = None, time: datetime = None, *, useastropy: bool = True
) -> typing.Tuple["ndarray", "ndarray", "ndarray"]:
    """
    Point => Point   ECEF => ECI

    Parameters
    ----------

    x : "ndarray"
        target x ECEF coordinate
    y : "ndarray"
        target y ECEF coordinate
    z : "ndarray"
        target z ECEF coordinate
    time : datetime.datetime
        time of observation
    useastropy : bool, optional
        use AstroPy for conversion


    Results
    -------
    x : "ndarray"
        target x ECI coordinate
    y : "ndarray"
        target y ECI coordinate
    z : "ndarray"
        target z ECI coordinate
    """
    # %%
    x = np.atleast_1d(x)
    y = np.atleast_1d(y)
    z = np.atleast_1d(z)
    if not x.shape == y.shape == z.shape:
        raise ValueError("shapes of ECI x,y,z must be identical")

    useastropy = useastropy and Time is not None

    if useastropy:
        gst = Time(time).sidereal_time("apparent", "greenwich").radian
    else:
        gst = datetime2sidereal(time, 0.0)

    gst = np.atleast_1d(gst)
    if gst.ndim != 1:
        raise ValueError("GST must be scalar or vector in radians")
    if gst.size == 1:
        gst *= np.ones(x.size)
    if gst.size != x.size:
        raise ValueError("GST must be scalar or same length as positions")

    ecef = np.column_stack((x.ravel(), y.ravel(), z.ravel()))
    eci = np.empty((x.size, 3))
    for i in range(x.size):
        eci[i, :] = _rottrip(gst[i]).T @ ecef[i, :]  # this one is transposed

    xeci = eci[:, 0].reshape(x.shape)
    yeci = eci[:, 1].reshape(x.shape)
    zeci = eci[:, 2].reshape(x.shape)

    return xeci, yeci, zeci


def _rottrip(ang: float) -> np.ndarray:
    """
    transformation matrix

    Parameters
    ----------

    ang : float
        angle to transform (radians)

    Returns
    -------
    T : numpy.ndarray of float
        3 x 3 transformation matrix
    """
    return np.array([[np.cos(ang), np.sin(ang), 0], [-np.sin(ang), np.cos(ang), 0], [0, 0, 1]])
