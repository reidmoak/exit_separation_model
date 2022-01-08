#!/usr/bin/env python3

import numpy as np

# CONSTANTS
g = 9.81

class Skydiver:

    def __init__(self, m):
        self.m = float(m)

    def compute_u(self, t, Va, Q):
        u = 2*self.m*Va/(Q*Va*t + 2*self.m)
        return u
    
    def compute_w(self, t, Va, Q):
        w = np.sqrt(2*self.m*g/Q)*np.tanh(-1*t*np.sqrt(Q*g/(2*self.m)))
        return w

    def compute_x(self, t, x_offset, Va, Q):
        x = 2*self.m/Q * np.log(Q*Va/(2*self.m)*t + 1) + x_offset
        return x

    def compute_z(self, t, z0, Va, Q):
        z = -2*self.m/Q * np.log(np.cosh(t*np.sqrt(g*Q/(2*self.m)))) + z0
        return z
