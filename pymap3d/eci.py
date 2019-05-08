""" transforms involving ECI earth-centered inertial """

from datetime import datetime
import numpy as np
from typing import Tuple

from .sidereal import datetime2sidereal
try:
    from astropy.time import Time
except ImportError:
    Time = None


def eci2ecef(x: float, y: float, z: float = None,
             time: datetime = None, *,
             useastropy: bool = True) -> Tuple[float, float, float]:
    """
    Observer => Point  ECI  =>  ECEF

    defaults to J2000 frame

    Parameters
    ----------
    x : float
        ECI x-location [meters]
    y : float
        ECI y-location [meters]
    z : float
        ECI z-location [meters]
    time : datetime.datetime
        time of obsevation (UTC)
    useastropy : bool, optional
        use AstroPy for conversion

    Results
    -------
    x : float
        target x ECEF coordinate
    y : float
        target y ECEF coordinate
    z : float
        target z ECEF coordinate
    """
# %%
    # FIXME: temporary old API, which was a single N x 3 numpy.ndarray
    if z is None and isinstance(y, (str, datetime)):
        time = y
        z = x[:, 2]
        y = x[:, 1]
        x = x[:, 0]
# %%
    x = np.atleast_1d(x)
    y = np.atleast_1d(y)
    z = np.atleast_1d(z)
    if not x.shape == y.shape == z.shape:
        raise ValueError('shapes of ECI x,y,z must be identical')

    useastropy = useastropy and Time is not None

    if useastropy:
        gst = Time(time).sidereal_time('apparent', 'greenwich').radian
    else:
        gst = datetime2sidereal(time, 0.)

    gst = np.atleast_1d(gst)
    if gst.ndim != 1 or not isinstance(gst[0], float):
        raise ValueError('GST must be vector in radians')
    if gst.size == 1:
        gst *= np.ones(x.size)
    if gst.size != x.size:
        raise ValueError('GST must be scalar or same length as positions')

    eci = np.column_stack((x.ravel(), y.ravel(), z.ravel()))
    ecef = np.empty((x.size, 3))
    for i in range(eci.shape[0]):
        ecef[i, :] = _rottrip(gst[i]) @ eci[i, :]

    xecef = ecef[:, 0].reshape(x.shape)
    yecef = ecef[:, 1].reshape(x.shape)
    zecef = ecef[:, 2].reshape(x.shape)

    return xecef, yecef, zecef


def ecef2eci(x: float, y: float, z: float = None,
             time: datetime = None, *,
             useastropy: bool = True) -> Tuple[float, float, float]:
    """
    Point => Point   ECEF => ECI

    Parameters
    ----------

    x : float
        target x ECEF coordinate
    y : float
        target y ECEF coordinate
    z : float
        target z ECEF coordinate
    time : datetime.datetime
        time of observation
    useastropy : bool, optional
        use AstroPy for conversion


    Results
    -------
    x : float
        target x ECI coordinate
    y : float
        target y ECI coordinate
    z : float
        target z ECI coordinate
    """
# %%
    # FIXME: temporary old API, which was a single N x 3 numpy.ndarray
    if z is None and isinstance(y, (str, datetime)):
        time = y
        z = x[:, 2]
        y = x[:, 1]
        x = x[:, 0]
# %%
    x = np.atleast_1d(x)
    y = np.atleast_1d(y)
    z = np.atleast_1d(z)
    if not x.shape == y.shape == z.shape:
        raise ValueError('shapes of ECI x,y,z must be identical')

    useastropy = useastropy and Time is not None

    if useastropy:
        gst = Time(time).sidereal_time('apparent', 'greenwich').radian
    else:
        gst = datetime2sidereal(time, 0.)

    gst = np.atleast_1d(gst)
    if gst.ndim != 1 or not isinstance(gst[0], float):
        raise ValueError('GST must be vector in radians')
    if gst.size == 1:
        gst *= np.ones(x.size)
    if gst.size != x.size:
        raise ValueError('GST must be scalar or same length as positions')

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
    return np.array([[np.cos(ang), np.sin(ang), 0],
                     [-np.sin(ang), np.cos(ang), 0],
                     [0, 0, 1]])
