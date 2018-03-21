
# coding: utf-8

# In[ ]:

# ***Question 2 (Keplerian Orbital Elements to Cartesian)
import numpy as np
import math

def COE2RV(a, e, i, cap_omega, omega, nu, mu):
    radius = compute_radius(a, e, nu)
    h = compute_angular_momentum(a, e, mu)
    r_eci = compute_position(radius, i, cap_omega, omega, nu)
    p = compute_period(a, e)
    v_eci = compute_velocity(r_eci, h, e, radius, p, nu, omega, i, cap_omega)
    return (r_eci, v_eci)
    
def compute_period(a, e):
    p = a*(1-math.pow(e, 2))
    return p

# 4. Compute radius
def compute_radius(a, e, nu):
    p = a*(1 - math.pow(e, 2))
    denominator = 1 + e*math.cos(nu)
    radius = p/denominator
    #print('radius = ', radius)
    return radius

# 5. Compute th specific angular momentum, h
def compute_angular_momentum(a, e, mu):
    term = 1 - math.pow(e, 2)
    h = math.sqrt(mu*a*term)
    #print('h = ', h)
    return h

# 6.Compute position components (ECI)
def compute_position(radius, i, cap_omega, omega, nu):
    x_eci = radius*(math.cos(cap_omega)*math.cos(omega + nu) -                     math.sin(cap_omega)*math.sin(omega + nu)*math.cos(i))
    y_eci = radius*(math.sin(cap_omega)*math.cos(omega + nu) +                     math.cos(cap_omega)*math.sin(omega + nu)*math.cos(i))
    z_eci = radius*(math.sin(i)*math.sin(omega + nu))
    r_eci = [x_eci, y_eci, z_eci]
    #print('R_ECI = (', r_eci[0], ',', r_eci[1], ',', r_eci[2], ')')
    return r_eci

# 7.Compute velocity components (ECI)
def compute_velocity(r_eci, h, e, radius, p, nu, omega, i, cap_omega):
    vx_eci = (r_eci[0]*h*e/(radius*p))*math.sin(nu) -                 (h/radius)*(math.cos(cap_omega)*math.sin(omega + nu) +                 math.sin(cap_omega)*math.cos(omega + nu)*math.cos(i))
    vy_eci = (r_eci[1]*h*e/(radius*p))*math.sin(nu) -                 (h/radius)*(math.sin(cap_omega)*math.sin(omega + nu) -                 math.cos(cap_omega)*math.cos(omega + nu)*math.cos(i))
    vz_eci = (r_eci[2]*h*e/(radius*p))*math.sin(nu) + (h/radius)*math.sin(i)*math.cos(omega + nu)
    v_eci = [vx_eci, vy_eci, vz_eci]
    #print('V_ECI = (', v_eci[0], ',', v_eci[1], ',', v_eci[2], ')')
    return v_eci

