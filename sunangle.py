import numpy

def get_day_ephemeris(year, month, day):
    return 367 * year - (7 * (year + (month + 9) / 12) / 4) + (275 * month / 9) + day - 730531.5

def get_zenith_cos_angle(time, lats, lons):
    """Calculate the solar Zenith angle -- the angle of the sun above the horizon -- and return
    its cosine.  Time is a Python datetime object (ensure the timezone is set correctly, or use UTC).
    Lats and lons are Numpy arrays of size M and N, with latitude and longitude in degrees north
    and east.  The result is an array of size MxN."""

    PI = numpy.pi
    TWO_PI = 2 * numpy.pi

    MAX_ITERATIONS = 5

    day_ephemeris = get_day_ephemeris(time.year, time.month, time.day)
    utnew = numpy.zeros_like(lons)
    lats = numpy.radians([lats]).transpose()
    sinPhi = numpy.sin(lats)
    cosPhi = numpy.cos(lats)
    lons = numpy.radians(lons)
    for ct in range(MAX_ITERATIONS):
        utold = utnew
        days = day_ephemeris + utold / TWO_PI
        t = days / 36525
        l = 4.8949504201433 + 628.331969753199 * t
        g = 6.2400408 + 628.3019501 * t
        ec = 0.033423 * numpy.sin(g) + 0.00034907 * numpy.sin(2 * g)
        lam = l + ec
        e = -ec + 0.0430398 * numpy.sin(2 * lam) - 0.00092502 * numpy.sin(4 * lam)
        utnew = PI - e - lons

    obl = 0.409093 - 0.0002269 * t
    delta = numpy.sin(obl) * numpy.sin(lam)
    delta = numpy.arcsin(delta)

    hour = (time.hour + time.minute/60.0 + time.second/3600.0)/24.0*TWO_PI
    h = hour - utnew
    s = sinPhi * numpy.sin(delta) + cosPhi * numpy.cos(delta) * numpy.cos(h)
    return s


if __name__ == '__main__':
    import pytz, datetime
    zone = pytz.timezone("UTC")
    now = datetime.datetime.now(zone)
    lat, long = -41, 174   # Wellington
    zenith_angle = numpy.degrees(numpy.arccos(get_zenith_cos_angle(now, lat, long)))
    print 'Sun angle:', zenith_angle
