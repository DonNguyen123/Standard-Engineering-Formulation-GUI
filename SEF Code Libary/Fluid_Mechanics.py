#Don's Standard Engineering Formulation

import math
'''
=====================================
Standard Fluid Equations
=====================================
'''

"""
Calculates the discharge (Q) and velocity (v) using Manning's Equation.

Parameters:
k (float): Manning's roughness coefficient
n (float): Roughness coefficient
a (float): Cross-sectional area of flow (ft^2 or m^2)
rh (float): Hydraulic radius (ft or m)
s (float): Slope (ft/ft or m/m)

Returns:
q (float): Discharge (ft^3/sec or m^3/s)
v (float): Velocity (ft/sec or m/s)
"""

def manning_equation(k, n, a, rh, s):
    q = (k / n) * a * (rh**(2/3)) * (s**(1/2))
    v = (k / n) * (rh**(2/3)) * (s**(1/2))

    # Create a formatted string representation of the results
    class FormattedDict(dict):
        def __str__(self):
            return '\n'.join(f'{key}: {value:_.4g}' if isinstance(value, (int, float)) 
                           else f'{key}: {value}'
                           for key, value in self.items())
        
        def __repr__(self):
            return self.__str__()

    return FormattedDict({
        'Discharge': q,
        'Velocity': v,
    })

"""
Calculates the discharge (Q) and velocity (v) using Hazen-Williams Equation.

Parameters:
k1 (float): Coefficient based on units (0.849 for SI units, 1.318 for USCS units)
c (float): Roughness coefficient, as tabulated in the Civil Engineering section
A (float): Cross sectional area (ft^2 or m^2)
S (float): Slope (ft/ft or m/m)
rh (float): Hydraulic radius (ft or m)

Returns:
q (float): Discharge (ft^3/sec or m^3/s)
v (float): Velocity (ft/sec or m/s)
"""

def hazen_williams_equation(k1, c, S, A, rh):
    v = k1 * c * (rh**0.63) * (S**0.54)
    q = k1 * c * A * (rh**0.63) * (S**0.54)

    # Create a formatted string representation of the results
    class FormattedDict(dict):
        def __str__(self):
            return '\n'.join(f'{key}: {value:_.4g}' if isinstance(value, (int, float)) 
                           else f'{key}: {value}'
                           for key, value in self.items())
        
        def __repr__(self):
            return self.__str__()

    return FormattedDict({
        'Discharge': q,
        'Velocity': v,
    })

"""
Calculates the discharge (Q) from a submerged orifice operating under steady-flow conditions.

Parameters:
c (float): Coefficient of discharge of the orifice
a (float): Cross-sectional area of the orifice (ft^2 or m^2)
g (float): Acceleration due to gravity (ft/s^2 or m/s^2)
h1 (float): Upstream fluid depth (ft or m)
h2 (float): Downstream fluid depth (ft or m)

Returns:
q (float): Discharge (ft^3/sec or m^3/s)
"""
    
def submerged_orifice_discharge(c, a, g, h1, h2):
    q = c * a * math.sqrt(2 * g * (h1 - h2))
    return q



"""
Calculates the discharge (Q) from an orifice discharging freely into the atmosphere.

Parameters:
c (float): Coefficient of discharge of the orifice
a0 (float): Cross-sectional area of the orifice opening (ft^2 or m^2)
g (float): Acceleration due to gravity (ft/s^2 or m/s^2)
h (float): Height of fluid above the orifice (ft or m)

Returns:
q (float): Volumetric flow (ft^3/sec or m^3/s)
"""
    
def free_orifice_discharge(c, a0, g, h):
    q = c * a0 * math.sqrt(2 * g * h)
    return q



"""
Calculates the time required to drain a tank.

Parameters:
a (float): Cross-sectional area of the tank (ft^2 or m^2)
a0 (float): Cross-sectional area of the orifice opening (ft^2 or m^2)
h1 (float): Initial height of fluid in the tank (ft or m)
h2 (float): Final height of fluid in the tank (ft or m)
g (float): Gravity (m/s^2 or ft/s^2)

Returns:
dt (float): Time required to drain the tank (seconds)
"""
    
def tank_drain_time(a, a0, h1, h2, g):
    dt = ((2 * (a / a0)) / (math.sqrt(2 * g))) * (h1**(1/2) - h2**(1/2))
    return dt


