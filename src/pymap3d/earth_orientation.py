"""Internal Earth orientation helpers for pure-Python inertial/terrestrial transforms."""

from __future__ import annotations

from datetime import datetime, timedelta, timezone
from math import cos, pi, sin, tau

ARCSEC2RAD = pi / 648000.0

# UTC effective dates and cumulative TAI-UTC leap seconds.
# Sources: IERS/ITU leap second history through 2017-01-01.
_LEAP_SECONDS = (
    (datetime(1972, 1, 1, tzinfo=timezone.utc), 10),
    (datetime(1972, 7, 1, tzinfo=timezone.utc), 11),
    (datetime(1973, 1, 1, tzinfo=timezone.utc), 12),
    (datetime(1974, 1, 1, tzinfo=timezone.utc), 13),
    (datetime(1975, 1, 1, tzinfo=timezone.utc), 14),
    (datetime(1976, 1, 1, tzinfo=timezone.utc), 15),
    (datetime(1977, 1, 1, tzinfo=timezone.utc), 16),
    (datetime(1978, 1, 1, tzinfo=timezone.utc), 17),
    (datetime(1979, 1, 1, tzinfo=timezone.utc), 18),
    (datetime(1980, 1, 1, tzinfo=timezone.utc), 19),
    (datetime(1981, 7, 1, tzinfo=timezone.utc), 20),
    (datetime(1982, 7, 1, tzinfo=timezone.utc), 21),
    (datetime(1983, 7, 1, tzinfo=timezone.utc), 22),
    (datetime(1985, 7, 1, tzinfo=timezone.utc), 23),
    (datetime(1988, 1, 1, tzinfo=timezone.utc), 24),
    (datetime(1990, 1, 1, tzinfo=timezone.utc), 25),
    (datetime(1991, 1, 1, tzinfo=timezone.utc), 26),
    (datetime(1992, 7, 1, tzinfo=timezone.utc), 27),
    (datetime(1993, 7, 1, tzinfo=timezone.utc), 28),
    (datetime(1994, 7, 1, tzinfo=timezone.utc), 29),
    (datetime(1996, 1, 1, tzinfo=timezone.utc), 30),
    (datetime(1997, 7, 1, tzinfo=timezone.utc), 31),
    (datetime(1999, 1, 1, tzinfo=timezone.utc), 32),
    (datetime(2006, 1, 1, tzinfo=timezone.utc), 33),
    (datetime(2009, 1, 1, tzinfo=timezone.utc), 34),
    (datetime(2012, 7, 1, tzinfo=timezone.utc), 35),
    (datetime(2015, 7, 1, tzinfo=timezone.utc), 36),
    (datetime(2017, 1, 1, tzinfo=timezone.utc), 37),
)


def ensure_utc(time: datetime) -> datetime:
    """Normalize datetimes to UTC while preserving naive UTC behavior."""

    if time.tzinfo is None:
        return time.replace(tzinfo=timezone.utc)

    return time.astimezone(timezone.utc)


def utc_with_offset(time: datetime, seconds: float) -> datetime:
    """Apply a seconds offset to a UTC-normalized datetime."""

    return ensure_utc(time) + timedelta(seconds=seconds)


def leap_seconds(time: datetime) -> int:
    """TAI-UTC for the supplied UTC epoch."""

    utc = ensure_utc(time)

    leaps = _LEAP_SECONDS[0][1]
    for effective, total in _LEAP_SECONDS:
        if utc >= effective:
            leaps = total
        else:
            break

    return leaps


def utc_to_tt(time: datetime) -> datetime:
    """Convert UTC to TT using the embedded leap-second table."""

    return utc_with_offset(time, leap_seconds(time) + 32.184)


def juliandate(time: datetime) -> float:
    """
    Python datetime to Julian date.

    Adapted from Vallado/Meeus using UTC-normalized inputs.
    """

    time = ensure_utc(time)

    if time.month < 3:
        year = time.year - 1
        month = time.month + 12
    else:
        year = time.year
        month = time.month

    A = int(year / 100.0)
    B = 2 - A + int(A / 4.0)
    C = (
        (
            (
                (
                    time.second
                    + time.microsecond / 1e6
                )
                / 60.0
                + time.minute
            )
            / 60.0
            + time.hour
        )
        / 24.0
    )

    return (
        int(365.25 * (year + 4716))
        + int(30.6001 * (month + 1))
        + time.day
        + B
        - 1524.5
        + C
    )


def julian_centuries(jd: float) -> float:
    """Julian centuries from J2000."""

    return (jd - 2451545.0) / 36525.0


def greenwich_mean_sidereal_time(jd_ut1: float) -> float:
    """Greenwich mean sidereal time in radians."""

    t_ut1 = julian_centuries(jd_ut1)
    gmst_sec = (
        67310.54841
        + (876600 * 3600 + 8640184.812866) * t_ut1
        + 0.093104 * t_ut1**2
        - 6.2e-6 * t_ut1**3
    )

    return gmst_sec * tau / 86400.0 % tau


def earth_rotation_angle(jd_ut1: float) -> float:
    """Earth rotation angle in radians, per IAU 2000 formulation."""

    return tau * (0.7790572732640 + 1.00273781191135448 * (jd_ut1 - 2451545.0)) % tau


def mean_obliquity(ttt: float) -> float:
    """Mean obliquity of the ecliptic in radians."""

    eps_arcsec = 84381.448 - 46.8150 * ttt - 0.00059 * ttt**2 + 0.001813 * ttt**3
    return eps_arcsec * ARCSEC2RAD


def nutation_short(ttt: float) -> tuple[float, float]:
    """
    Short-form nutation model in radians.

    Uses the dominant Meeus/Vallado terms to capture most of the nutation signal
    without bringing in the full IAU 2000B coefficient set.
    """

    L = (280.4665 + 36000.7698 * ttt) % 360.0
    Lp = (218.3165 + 481267.8813 * ttt) % 360.0
    omega = (
        125.04452 - 1934.136261 * ttt + 0.0020708 * ttt**2 + ttt**3 / 450000.0
    ) % 360.0

    L = L * pi / 180.0
    Lp = Lp * pi / 180.0
    omega = omega * pi / 180.0

    dpsi_arcsec = (
        -17.20 * sin(omega)
        - 1.32 * sin(2.0 * L)
        - 0.23 * sin(2.0 * Lp)
        + 0.21 * sin(2.0 * omega)
    )
    deps_arcsec = (
        9.20 * cos(omega)
        + 0.57 * cos(2.0 * L)
        + 0.10 * cos(2.0 * Lp)
        - 0.09 * cos(2.0 * omega)
    )

    return dpsi_arcsec * ARCSEC2RAD, deps_arcsec * ARCSEC2RAD


def equation_of_equinoxes(ttt: float) -> float:
    """Approximate equation of the equinoxes in radians."""

    dpsi, deps = nutation_short(ttt)
    return dpsi * cos(mean_obliquity(ttt) + deps)


def greenwich_apparent_sidereal_time(jd_ut1: float, jd_tt: float) -> float:
    """Greenwich apparent sidereal time in radians."""

    ttt = julian_centuries(jd_tt)
    gast = greenwich_mean_sidereal_time(jd_ut1) + equation_of_equinoxes(ttt)
    return gast % tau


def rotation_x(angle: float) -> tuple[tuple[float, float, float], ...]:
    c = cos(angle)
    s = sin(angle)
    return ((1.0, 0.0, 0.0), (0.0, c, s), (0.0, -s, c))


def rotation_y(angle: float) -> tuple[tuple[float, float, float], ...]:
    c = cos(angle)
    s = sin(angle)
    return ((c, 0.0, -s), (0.0, 1.0, 0.0), (s, 0.0, c))


def rotation_z(angle: float) -> tuple[tuple[float, float, float], ...]:
    c = cos(angle)
    s = sin(angle)
    return ((c, s, 0.0), (-s, c, 0.0), (0.0, 0.0, 1.0))


def matmul3(a, b):
    """3x3 matrix multiply."""

    return tuple(
        tuple(sum(a[i][k] * b[k][j] for k in range(3)) for j in range(3)) for i in range(3)
    )


def matvec3(a, v):
    """3x3 matrix by 3-vector."""

    return tuple(sum(a[i][k] * v[k] for k in range(3)) for i in range(3))


def transpose3(a):
    """Transpose a 3x3 matrix."""

    return tuple(tuple(a[j][i] for j in range(3)) for i in range(3))


def precession_matrix(ttt: float):
    """IAU 1976 precession matrix, J2000 mean equator/equinox to mean-of-date."""

    zeta = (
        2306.2181 * ttt + 0.30188 * ttt**2 + 0.017998 * ttt**3
    ) * ARCSEC2RAD
    theta = (
        2004.3109 * ttt - 0.42665 * ttt**2 - 0.041833 * ttt**3
    ) * ARCSEC2RAD
    zed = (2306.2181 * ttt + 1.09468 * ttt**2 + 0.018203 * ttt**3) * ARCSEC2RAD

    return matmul3(rotation_z(-zed), matmul3(rotation_y(theta), rotation_z(-zeta)))


def nutation_matrix(ttt: float):
    """Nutation matrix, mean-of-date to true-of-date."""

    dpsi, deps = nutation_short(ttt)
    eps0 = mean_obliquity(ttt)
    eps = eps0 + deps

    return matmul3(rotation_x(-eps), matmul3(rotation_z(-dpsi), rotation_x(eps0)))


def polar_motion_matrix(xp_arcsec: float = 0.0, yp_arcsec: float = 0.0):
    """Polar motion matrix, true Earth-fixed to ITRF-like ECEF."""

    xp = xp_arcsec * ARCSEC2RAD
    yp = yp_arcsec * ARCSEC2RAD
    return matmul3(rotation_y(xp), rotation_x(yp))


def eci_to_ecef_matrix(
    time: datetime,
    delta_ut1: float = 0.0,
    xp: float = 0.0,
    yp: float = 0.0,
):
    """
    Rotation matrix for the pure-Python ECI->ECEF path.

    The model uses:
    - UTC -> UT1 with user-supplied ``delta_ut1`` seconds
    - UTC -> TT via embedded leap seconds
    - IAU 1976 precession
    - truncated nutation / apparent sidereal time
    - optional polar motion ``xp`` / ``yp`` in arcseconds
    """

    jd_ut1 = juliandate(utc_with_offset(time, delta_ut1))
    jd_tt = juliandate(utc_to_tt(time))
    ttt = julian_centuries(jd_tt)

    p = precession_matrix(ttt)
    n = nutation_matrix(ttt)
    st = rotation_z(greenwich_apparent_sidereal_time(jd_ut1, jd_tt))
    pm = polar_motion_matrix(xp, yp)

    return matmul3(pm, matmul3(st, matmul3(n, p)))
