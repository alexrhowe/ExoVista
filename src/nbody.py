import numpy as np
from src import coordinates as coord

max_n_planets = 100
max_n_equations = 6*max_n_planets
epsilon = 1.e-12
G = 4.*np.pi*np.pi # gravitational constant in AU^3*SM^-1*yr^-2

def nbody(cartin, GM, curr_time, desired_time):
    n_planets = len(GM)
    n_equations = 6*n_planets

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
    
    while t < desired_time:
        state, t, deltat = integrate(GM, state, t, deltat)
    zero_deltat = desired_time - t
    while abs(zero_deltat) > epsilon:
        state, t, zero_deltat = integrate(GM, state, t, zero_deltat)
        zero_deltat = desired_time - t
    for i in range(0,n_planets):
        c1 = 6*i
        x[i] = state[c1]
        y[i] = state[c1+1]
        z[i] = state[c1+2]
        vx[i] = state[c1+3]
        vy[i] = state[c1+4]
        vz[i] = state[c1+5]

    cartout = [x,y,z,vx,vy,vz]
    return cartout

# Bulirsch-Stoer integrator
# Note: this is a Python translation of a C translation of a FORTRAN code.
# Original FORTRAN by Michael Henon and Jack Wisdom
def integrate(GM, state, t, deltat):
    
    lt_BS = [1,2,3,4,6,8,12,16,24,32]
    temp_BS = np.zeros((max_n_equations,12))
    d_BS = np.zeros(6)
    n_equations = len(state)
    
    delta_t_BS = deltat
    xa_BS = t
    
    deriv = equations_of_motion(GM, state)

    for i_BS in range(0,n_equations):
        temp_BS[i_BS,1] = abs(state[i_BS])
        if temp_BS[i_BS,1] < epsilon: temp_BS[i_BS,1] = epsilon
        temp_BS[i_BS,4] = deriv[i_BS]
        temp_BS[i_BS,0] = state[i_BS]

    done = False
    while not done:
        
        xb_BS = delta_t_BS + xa_BS
    
        for m_BS in range(0,10):
            flt_BS = lt_BS[m_BS]
            varm_BS = 0.
            m1_BS = min(m_BS,6)
            if m1_BS != 0:
                for k_BS in range(0,m1_BS):
                    d_BS[k_BS] = (flt_BS/lt_BS[m_BS-k_BS-1])**2
            h_BS = delta_t_BS / flt_BS
            hd_BS = 0.5*h_BS
            
            for i_BS in range(0,n_equations):
                temp_BS[i_BS,3] = temp_BS[i_BS,0]
                state[i_BS] = temp_BS[i_BS,0] + hd_BS*temp_BS[i_BS,4]
                
            i1max_BS = 2*flt_BS - 1
            t = xa_BS
            
            for i1_BS in range(0,i1max_BS):
                t += hd_BS
                deriv = equations_of_motion(GM,state)
                
                for i_BS in range(0,n_equations):
                    foo_BS = state[i_BS]
                    temp_BS[i_BS,1] = max(temp_BS[i_BS,1], abs(foo_BS))
                    if temp_BS[i_BS,1] == 0.:
                        print('FAIL: temp_BS[i_BS.1] = 1')
                        return
                    eta2_BS = temp_BS[i_BS,3] + h_BS*deriv[i_BS]
                    temp_BS[i_BS,3] = state[i_BS]
                    state[i_BS] = eta2_BS
                    
            t = xb_BS
            deriv = equations_of_motion(GM,state)
            
            for i_BS in range(0,n_equations):
                dta_BS = temp_BS[i_BS,11]
                yb_BS = (temp_BS[i_BS,3] + state[i_BS] + hd_BS*deriv[i_BS]) / 2.
                c_BS = yb_BS
                temp_BS[i_BS,11] = yb_BS
                
                if m1_BS != 0:
                    for k_BS in range(0,m1_BS):
                        b1_BS = d_BS[k_BS] * dta_BS
                        den_BS = b1_BS - c_BS
                        dtn_BS = dta_BS
                        if den_BS != 0:
                            b_BS = (c_BS - dta_BS) / den_BS
                            dtn_BS = c_BS * b_BS
                            c_BS = b1_BS * b_BS
                        dta_BS = temp_BS[i_BS,11-k_BS-1]
                        temp_BS[i_BS,11-k_BS-1] = dtn_BS
                        yb_BS += dtn_BS
                    var_BS = abs(temp_BS[i_BS,2]-yb_BS) / temp_BS[i_BS,1]
                    varm_BS = max(varm_BS, var_BS)
                temp_BS[i_BS,2] = yb_BS
                
            if m_BS < 3: varma_BS = varm_BS
            elif varm_BS <= epsilon:
                done = True
                break
            elif varm_BS >= varma_BS: break
            
        if not done: delta_t_BS /= 2.

    t = xb_BS
    for i_BS in range(0,n_equations):
        state[i_BS] = temp_BS[i_BS,2]
    delta_t = delta_t_BS * 1.5 * 0.6**float(m_BS-m1_BS)

    return state, t, delta_t

def equations_of_motion(GM,state):
    n_planets = int(len(state)/6)
    deriv = np.zeros(len(state))

    # velocity from state
    for i in range(0,n_planets):
        c1 = 6*i
        deriv[c1]   = state[c1+3]
        deriv[c1+1] = state[c1+4]
        deriv[c1+2] = state[c1+5]
        deriv[c1+3] = 0.
        deriv[c1+4] = 0.
        deriv[c1+5] = 0.

    # pairwise accelerations
    for i in range(0,n_planets):
        c1 = 6*i
        for j in range(0,i):
            c2 = 6*j
            dx = state[c1  ] - state[c2]
            dy = state[c1+1] - state[c2+1]
            dz = state[c1+2] - state[c2+2]
            r2 = dx*dx + dy*dy + dz*dz
            fij = 1./(r2**1.5)
            deriv[c1+3] += -GM[j]*fij*dx
            deriv[c1+4] += -GM[j]*fij*dy
            deriv[c1+5] += -GM[j]*fij*dz
            deriv[c2+3] += GM[i]*fij*dx
            deriv[c2+4] += GM[i]*fij*dy
            deriv[c2+5] += GM[i]*fij*dz
            
    return deriv
