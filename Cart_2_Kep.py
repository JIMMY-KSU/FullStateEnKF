
# coding: utf-8

# In[ ]:

import numpy as np
import math

def RV2COE(R, V, mu):

    (h_vector, h) = compute_angular_momentum(R, V)
    (r, v) = compute_rv_magnitude(R, V)
    E = compute_specific_energy(r, v, mu)
    a = compute_semimajor_axis(mu, E)
    e = compute_eccentricity(h, a, mu)
    i = compute_inclination(h_vector, h)
    cap_omega = compute_RAAN(h_vector)
    arg_of_latitude = compute_arg_of_latitude(R, cap_omega, i)
    p = compute_period(a, e)
    nu = compute_true_anomaly(R, V, r, p, mu)
    omega = compute_arg_of_periapse(arg_of_latitude, nu)
    
    if (e < .0001): #omega = undefined
        print('Eccentricity is ~0 and omega is undefined. Use Argument of Latitude instead.')
        print('Argument of Latitude = ', arg_of_latitude)   
    if (i < .0001): #cap_omega = undefined
        print('The inclination is ~zero and RAAN is undefined.')
    if (e < .0001) and (i > -.0001 and i < .0001):
        print('The eccentricity & inclination is ~zero.')
        
    return (a, e, i, cap_omega, omega, nu)


# 1.Compute specific angular momentum
def compute_angular_momentum(R, V):
    h_vector = np.cross(R, V)
    h = math.sqrt(math.pow(h_vector[0], 2) + math.pow(h_vector[1], 2) + math.pow(h_vector[2], 2))
    #print('h = ', h)
    return (h_vector, h)

# 2.Compute the radius
def compute_rv_magnitude(R,V):
    r = math.sqrt(math.pow(R[0], 2) + math.pow(R[1], 2) + math.pow(R[2], 2))
    #print('r = ', r)
    v = math.sqrt(math.pow(V[0], 2) + math.pow(V[1], 2) + math.pow(V[2], 2))
    #print('v = ', v)
    return (r, v)

# 3.Compute specific energy & verify elliptical motion
def compute_specific_energy(r, v, mu):
    E = (math.pow(v, 2)/2) - (mu/r)
    #print('E = ', E)
    return E

# 4.Compute semi-major axis
def compute_semimajor_axis(mu, E):
    a = -mu/(2*E)
    #print('a = ', a)
    return a

# 5.Compute eccentricity
def compute_eccentricity(h, a, mu):
    e = math.sqrt(1 - math.pow(h, 2)/(a*mu))
    #print('e = ', e)
    return e

# 6.Compute inclination [0, 180]
def compute_inclination(h_vector, h):
    i = math.acos(h_vector[2]/h)
    #print('i = ', i)
    return i

# 7.Compute RAAN (cap omega) [0, 360]
def compute_RAAN(h_vector):
    cap_omega = math.atan2(h_vector[0], -h_vector[1])
    if (cap_omega < 0):
        cap_omega = cap_omega + 2*math.pi
    #print('cap_omega = ', cap_omega)
    return cap_omega

# 8.Compute argument of latitude, omega + v [0, 360]
def compute_arg_of_latitude(R, cap_omega, i):
    param2 = R[0]*math.cos(cap_omega) + R[1]*math.sin(cap_omega)
    arg_of_latitude = math.atan2(R[2]/math.sin(i), param2)
    #print('arg_of_latitude = ', arg_of_latitude)
    return arg_of_latitude

def compute_period(a, e):
    p = a*(1-math.pow(e, 2))
    return p

# 9.Compute true anomaly (nu) [0, 360]
def compute_true_anomaly(R, V, r, p, mu):
    V_dot_R = V[0]*R[0] + V[1]*R[1] + V[2]*R[2]
    first_term = math.sqrt(p/mu)*(V_dot_R)
    nu = math.atan2(first_term, p-r)
    if (nu < 0):
        nu = nu + 2*math.pi
    #print('nu (true anomaly) = ', nu)
    return nu

# 10.Compute argument of periapse, omega [0, 360]
def compute_arg_of_periapse(arg_of_latitude, nu):
    omega = arg_of_latitude - nu
    if (omega < 0):
        omega = omega + 2*math.pi
    #print('omega = ', omega)
    return omega

