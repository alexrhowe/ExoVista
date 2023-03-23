import numpy as np
from src import coordinates as coord
from src import wrapIntegrator

max_n_planets = 100
max_n_equations = 6*max_n_planets
epsilon = 1.e-12

def nbody(cartin, GM, R, istar, curr_time, desired_time):
    n_planets = len(GM)
    n_equations = 6*n_planets
    tlist = []
    transmaster = []
    for i in range(1,n_planets): transmaster.append([])
    
    x = cartin[0]
    y = cartin[1]
    z = cartin[2]
    vx = cartin[3]
    vy = cartin[4]
    vz = cartin[5]

    r1 = np.sqrt((x[1]-x[0])**2 + (y[1]-y[0])**2 + (z[1]-z[0])**2)
    v1 = np.sqrt((vx[1]-vx[0])**2 + (vy[1]-vy[0])**2 + (vz[1]-vz[0])**2)
    guess = r1/v1 * 0.01  # initial guess is 0.01/2pi inner planet orbits
    
    state = np.zeros(n_equations)
    for i in range(0,n_planets):
        c1 = 6*i
        state[c1] = x[i]
        state[c1+1] = y[i]
        state[c1+2] = z[i]
        state[c1+3] = vx[i]
        state[c1+4] = vy[i]
        state[c1+5] = vz[i]

    t = curr_time
    deltat = guess
    
    integrator = wrapIntegrator.PyIntegrator()
    integrator.SetupIntegrator(GM, state, t, deltat)

    step = 0
    
    while t < desired_time:
        integrator.integrate()
        state = integrator.getState()
        t = integrator.getTime()
        deltat = integrator.getDeltaT()
        
        mintrans, translist = detect_transits(R, istar, state)
        tlist.append(t)
        for i in range(0,n_planets-1): transmaster[i].append(translist[i])
        
        if mintrans < 10*deltat:
            deltat *= mintrans / 10*deltat
            integrator.setDeltaT(deltat)
        step += 1
            
    zero_deltat = desired_time - t

    mintrans, translist = detect_transits(R, istar, state)
    tlist.append(t)
    for i in range(0,n_planets-1): transmaster[i].append(translist[i])
    
    while abs(zero_deltat) > epsilon:
        integrator.setDeltaT(zero_deltat)
        
        integrator.integrate()
        state = integrator.getState()
        t = integrator.getTime()
        deltat = integrator.getDeltaT()
        
        mintrans, translist = detect_transits(R, istar, state)
        tlist.append(t)
        for i in range(0,n_planets-1): transmaster[i].append(translist[i])
        
        if mintrans < 10*deltat:
            deltat *= mintrans / 10*deltat
            integrator.setDeltaT(deltat)
        zero_deltat = desired_time - t
        step += 1
        
    for i in range(0,n_planets):
        c1 = 6*i
        x[i] = state[c1]
        y[i] = state[c1+1]
        z[i] = state[c1+2]
        vx[i] = state[c1+3]
        vy[i] = state[c1+4]
        vz[i] = state[c1+5]

    cartout = [x,y,z,vx,vy,vz]
    del integrator
    
    return cartout, tlist, transmaster


def detect_transits(R, istar, state):
    n_planets = len(R)
    sini = np.sin(istar*np.pi/180.)
    cosi = np.cos(istar*np.pi/180.)
    
    sep = np.zeros(n_planets-1)
    ttrans = np.zeros(n_planets-1)
    xstar = state[0]
    ystar = state[1] * cosi - state[2] * sini
    translist = np.zeros(n_planets-1)
    
    for i_pl  in range(1,n_planets):
        xpl = state[6*i_pl]
        ypl = state[6*i_pl+1] * cosi - state[6*i_pl+2] * sini
        zpl = state[6*i_pl+1] * sini - state[6*i_pl+2] * cosi
        vxpl = state[6*i_pl+3]
        vypl = state[6*i_pl+4] * cosi - state[6*i_pl+5] * sini
        
        sep[i_pl-1] = max(1, np.sqrt((xpl-xstar)**2 + (ypl-ystar)**2) / (R[0]+R[i_pl]))
        vxy = np.sqrt(vxpl**2 + vypl**2) / (R[0]+R[i_pl])
        ttrans[i_pl-1] = sep[i_pl-1]/vxy
        
        if sep[i_pl-1]==1 and zpl>0: translist[i_pl-1] = 1
        elif sep[i_pl-1]==1 and zpl<0: translist[i_pl-1] = -1
        
    return np.min(ttrans), translist
