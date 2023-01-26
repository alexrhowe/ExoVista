import numpy as np

# converts orbital elements to cartesian coordinates
def cartesian(GM, kepcoords):
    epsilon = 1.e-10
    
    a = kepcoords[0]
    ecc = kepcoords[1]
    i = kepcoords[2]
    longnode = kepcoords[3]
    argperi = kepcoords[4]
    meananom = kepcoords[5]
    
    # compute eccentric anomaly
    E0 = meananom * 1. # need a deep copy
    E1 = meananom + ecc * np.sin(E0)
    E2 = meananom + ecc * np.sin(E1)
    
    while np.max(np.abs(E0-E2)) > epsilon:
        E1 = meananom + ecc * np.sin(E0)
        E2 = meananom + ecc * np.sin(E1)
        den = E2 - 2.*E1 + E0
        j = np.where(np.abs(den) > epsilon)[0]
        E0[j] = E0[j] - (E1[j]-E0[j])*(E1[j]-E0[j])/den[j]
        j = np.where(np.abs(den) < epsilon)[0]
        E0[j] = E2[j]
        E2[j] = E1[j]
    cosE = np.cos(E0)
    sinE = np.sin(E0)
    
    # compute unrotated positions and velocities
    foo = np.sqrt(1.0 - ecc*ecc)
    meanmotion = np.sqrt(GM/(a*a*a))
    x = a * (cosE - ecc)
    y = foo * (a * sinE)
    z = np.zeros(len(y))
    denom = 1. / (1.0 - ecc * cosE)
    
    xd = (-a * meanmotion * sinE) * denom
    yd = foo * (a * meanmotion * cosE * denom)
    zd = np.zeros(len(yd))

    # rotate by argument of perihelion in the orbit plane
    cosw = np.cos(argperi)
    sinw = np.sin(argperi)
    xp = x * cosw - y * sinw
    yp = x * sinw + y * cosw
    zp = z
    xdp = xd * cosw - yd * sinw
    ydp = xd * sinw + yd * cosw
    zdp = zd
    
    # rotate by inclination about x axis
    cosi = np.cos(i)
    sini = np.sin(i)
    x = xp
    y = yp * cosi - zp * sini
    z = yp * sini + zp * cosi
    xd = xdp
    yd = ydp * cosi - zdp * sini
    zd = ydp * sini + zdp * cosi
    
    # rotate by longitude of node about z axis
    cosnode = np.cos(longnode)
    sinnode = np.sin(longnode)
    xf = x * cosnode - y * sinnode
    yf = x * sinnode + y * cosnode
    zf = z
    
    vx = xd * cosnode - yd * sinnode
    vy = xd * sinnode + yd * cosnode
    vz = zd
    
    x = xf
    y = yf
    z = zf

    return [x, y, z, vx, vy, vz]


# converts cartesian coordinates to orbital elements
def keplerian(GM, cartcoords):
    # NOTE: this does not correctly calculate the
    # meananom (or eccanom) for a hyperbolic orbit
    # (It relies on e<1.)

    x = cartcoords[0]
    y = cartcoords[1]
    z = cartcoords[2]
    vx = cartcoords[3]
    vy = cartcoords[4]
    vz = cartcoords[5]
    
    # find direction of angular momentum vector
    rxv_x = y * vz - z * vy
    rxv_y = z * vx - x * vz
    rxv_z = x * vy - y * vx
    hs = rxv_x * rxv_x + rxv_y * rxv_y + rxv_z * rxv_z
    
    h = np.sqrt(hs)
    
    r = np.sqrt(x * x + y * y + z * z)
    vs = vx * vx + vy * vy + vz * vz
    rdotv = x * vx + y * vy + z * vz
    rdot = rdotv / r
    parameter = hs / GM
    
    i = np.arccos(rxv_z / h)
    
    longnode = i * 0.0 # zero by default
    j = [jj for jj in range(0,len(longnode)) if rxv_x[jj] != 0 or rxv_y[jj] != 0]
    longnode[j] = np.arctan2(rxv_x[j], -rxv_y[j])

    a = 1.0 / (2.0 / r - vs / GM)
    
    ecostrueanom = parameter / r - 1.0
    esintrueanom = rdot * h / GM
    e = np.sqrt(ecostrueanom * ecostrueanom + esintrueanom * esintrueanom)
    
    trueanom = i * 0.0 # zero by default
    j = [jj for jj in range(0,len(esintrueanom)) if esintrueanom[jj] != 0 or ecostrueanom[jj] != 0]
    trueanom[j] = np.arctan2(esintrueanom[j], ecostrueanom[j])
    
    cosnode = np.cos(longnode)
    sinnode = np.sin(longnode)
    
    # u is the argument of latitude
    rcosu = x * cosnode + y * sinnode
    rsinu = (y * cosnode - x * sinnode)/np.cos(i)
    
    u = i * 0.0 # zero by default
    j = [j for j in range(0,len(rsinu)) if rsinu[j] != 0 or rcosu != 0]
    u[j] = np.arctan2(rsinu[j], rcosu[j])
    
    argperi = (u - trueanom) % (2.*np.pi)
    
    eccanom = 2.0 * np.arctan(np.sqrt((1.0 - e)/(1.0 + e)) * np.tan(trueanom/2.0))
    meananom = eccanom - e * np.sin(eccanom)
    
    return [a, e, i, longnode, argperi, meananom]
