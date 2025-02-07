import math
import numpy as np
import matplotlib.pyplot as plt

'''
=====================================
Standard Steady State Heat Transfer Functions 
=====================================
'''

"""
Calculate heat transfer rate using Fourier's Law of Conduction

Parameters:
k (float): thermal conductivity [W/(m·K)]
A (float): surface area perpendicular to direction of heat transfer (m²)
dT (float): temperature difference (K)
dx (float): thickness or distance (m)

Returns:
Q_dot (float): Rate of heat transfer (W)
"""
    
def conduction_rate(k, A, dT, dx):
    Q_dot = -k * A * (dT/dx)
    return Q_dot



"""
Calculate convective heat transfer rate using Newton's Law of Cooling

Parameters:
h (float): convection heat-transfer coefficient [W/(m²·K)]
A (float): convection surface area (m²)
T_wall (float): wall surface temperature (K)
T_fluid (float): bulk fluid temperature (K)

Returns:
Q_dot (float): Rate of heat transfer (W)
"""
    
def convection_rate(h, A, T_wall, T_fluid):
    Q_dot = h * A * (T_wall - T_fluid)
    return Q_dot



"""
Calculate radiation heat transfer rate

Parameters:
emissivity (float): emissivity of the body
A (float): body surface area (m²)
T (float): absolute temperature (K)
Unit (str): metric or imperial

Returns:
Q_dot (float): Rate of heat transfer (W)
"""

def radiation_rate(emmisstivity, A, T, unit):
    if unit == "metric":
        sigma = 5.67*  (10**(-8))  # Stefan-Boltzmann constant [W/(m²·K⁴)]
    if unit == "imperial":
        sigma = 1.713 * (10**(-9))  # Stefan-Boltzmann constant [BTU/(hr·ft²·R⁴)]
    Q_dot = emmisstivity * sigma * A * (T**4)
    return Q_dot












'''
=====================================
Standard Conduction Functions 
=====================================
'''

"""
Calculate heat transfer rate through a plane wall

Parameters:
k (float): thermal conductivity [W/(m·K)]
A (float): wall surface area normal to heat flow (m²)
L (float): wall thickness (m)
T1 (float): temperature of one surface (K)
T2 (float): temperature of other surface (K)

Returns:
Q_dot (float): Rate of heat transfer (W)
"""

def plane_wall_conduction(k, A, L, T1, T2):
    Q_dot = (k * A * (T1 - T2)) / L
    return Q_dot



"""
Calculate heat transfer rate through a cylindrical wall

Parameters:
k (float): thermal conductivity [W/(m·K)]
L (float): cylinder length (m)
r1 (float): inner radius (m)
r2 (float): outer radius (m)
T1 (float): inner surface temperature (K)
T2 (float): outer surface temperature (K)

Returns:
Q_dot (float): Rate of heat transfer (W)
"""
    
def cylindrical_wall_conduction(k, L, r1, r2, T1, T2):
    Q_dot = (2 * math.pi() * k * L * (T1 - T2)) / (math.log(r2/r1))
    return Q_dot



"""
Calculate critical insulation radius

Parameters:
k_insulation (float): thermal conductivity of insulation [W/(m·K)]
h_infinity (float): convection heat transfer coefficient [W/(m²·K)]

Returns:
radius (float): Critical insulation radius (m)
"""
    
def critical_insulation_radius(k_insulation, h_infinity):
    radius = k_insulation / h_infinity
    return radius












'''
=====================================
Thermal Resistance Functions 
=====================================
'''

"""
Calculate thermal resistance for a plane wall

Parameters:
L (float): wall thickness (m)
k (float): thermal conductivity [W/(m·K)]
A (float): surface area (m²)

Returns:
R (float): Thermal resistance (K/W)
"""
    
def thermal_resistance_plane_wall(L, k, A):
    R = L / (k * A)
    return R



"""
Calculate thermal resistance for a cylindrical wall

Parameters:
r1 (float): inner radius (m)
r2 (float): outer radius (m)
k (float): thermal conductivity [W/(m·K)]
L (float): cylinder length (m)

Returns:
R (float): Thermal resistance (K/W)
"""
    
def thermal_resistance_cylinder(r1, r2, k, L):
    R = (math.log(r2/r1)) / (2 * math.pi() * k * L)
    return R



"""
Calculate thermal resistance for convection

Parameters:
h (float): convection heat transfer coefficient [W/(m²·K)]
A (float): surface area (m²)

Returns:
R (float): Thermal resistance (K/W)
"""
    
def convection_resistance(h, A):
    R = 1 / (h * A)
    return R












'''
=====================================
Lumped Capactience Functions 
=====================================
'''

"""
Calculate temperature variation with time using lumped capacitance model

Parameters:
T_initial (float): initial body temperature (K)
T_fluid (float): fluid temperature (K)
h (float): convection heat transfer coefficient [W/(m²·K)]
A_s (float): surface area of body (m²)
density (float): density of body (kg/m³)
volumne (float): volume of body (m³)
c_p (float): heat capacity [J/(kg·K)]
t (float): time (s)
Bi (float): Biot number (unitless)

Returns:
T (float): Temperature at time t (K)
"""
    
def lumped_capacitance_temperature_varianece(T_initial, T_fluid, h, A_s, 
    density, volumne, c_p, t, Bi):
    if Bi <= 0.1:
        beta = h * A_s / (density * volumne * c_p)
        T = (T_initial - T_fluid) * math.exp(-beta * t)
        return T
    else:
        print("Error. Bi value is beneath 0.1, cannot used lumped capcitance model.")



"""
Calculate total heat transferred using lumped capacitance model

Parameters:
density (float): density of body (kg/m³)
volumne (float): volume of body (m³)
c_p (float): heat capacity [J/(kg·K)]
T_initial (float): initial temperature (K)
T_final (float): final temperature (K)
Bi (float): Biot number (unitless)

Returns:
Q (float): Total heat transferred (J)
"""
    
def lumped_capacitance_total_heat_transferred(density, volumne, c_p, T_initial, T_final, Bi):
    if Bi <= 0.1:
        Q = density * volumne * c_p * (T_initial - T_final)
        return Q
    else:
        print("Error. Bi value is beneath 0.1, cannot used lumped capcitance model.")
        
        










'''
=====================================
Sudden Convection Functions 
=====================================
'''

"""
Calculates the temperature at any location within a solid with sudden convection.

Parameters:
k (float): Thermal conductivity of the solid (W/mK)
rho (float): Density of the solid (kg/m^3)
c (float): Specific heat of the solid (J/kgK)

Returns:
alpha (float): Thermal Diffustivity Value (m^2 / s)
"""
    
def calculate_thermal_diffusivity(k, rho, c):
    alpha = k / (rho * c)
    return alpha



"""
Calculates the temperature at any location within a solid with sudden convection.

Parameters:
L (float): Half-thickness of the plane wall, or radius of the cylinder/sphere (m)
t (float): Time (s)
x (float): Distance from the midplane of the wall (m)
r (float): Radial distance from the centerline of the cylinder/centerpoint of the sphere (m)
T_internal (float): Internal center temperture of object (K)
T_external (float): External temperture of object (k)
T_surrounding (float): External temperture of surrounding
xi (float): Dimentionless reference values
C (float): Dimentionless reference values
Fo (float): Fourier_number
type (text): Can be cylcinder, sphere, or wall

Returns:
T (float): Temperature at the specified location and time (K)
"""
    
def solid_sudden_convection(T_internal, T_external, T_surrounding, Fo, a, xi, C, type, L = None, r = None):
    if Fo >= 0.2:
        if type == "plane_wall":
            T = (((np.log((((T_external-T_surrounding)/(T_internal-T_surrounding))/ C))) / (xi**2)) * (L**2)) / a
            return T
        elif type == "cylinder" or type == "sphere":
            T = (((np.log((((T_external-T_surrounding)/(T_internal-T_surrounding))/ C))) / (xi**2)) * (r**2)) / a
            return T
        else:
            print("Type not supported. Pick cyclidner, sphere, or plane_wall")
    else:
        print("Fourier value is not over 0.2, cannot get time for sudden convection.")
    











'''
=====================================
Radition Between Bodies Functions 
=====================================
'''
"""
Calculates the net heat transfer rate from a small body to its surroundings.

Parameters:
emissivity (float): Emissivity of the body surface
surface_area (float): Surface area of the body (m^2)
body_temp (float): Absolute temperature of the body surface (K)
surrounding_temp (float): Absolute temperature of the surroundings (K)
unit: metric or imperial

Returns:
Q (float): Net heat transfer rate from the body (W)
"""
    
def net_heat_transfer_small_body_to_enviroment(emissivity, surface_area, body_temp, surrounding_temp, unit):
    if unit == "metric":
        sigma = 5.67*  (10**(-8))  # Stefan-Boltzmann constant [W/(m²·K⁴)]
    if unit == "imperial":
        sigma = 1.713 * (10**(-9))  # Stefan-Boltzmann constant [BTU/(hr·ft²·R⁴)]
    Q = emissivity * sigma * surface_area * (body_temp**4 - surrounding_temp**4)
    return Q



"""
Calculates the net energy exchange by radiation between two black bodies.

Parameters:
internal_surface_area (float): Surface area of the first black body (m^2)
view_factor (float): View factor between the two bodies
internal_body_temp (float): Absolute temperature of the first black body (K)
enclosuing_body_temp (float): Absolute temperature of the second black body (K)
unit: metric or imperial

Returns:
Q (float): Net heat transfer rate between the two black bodies (W)
"""
    
def net_heat_transfer_two_black_bodies(internal_surface_area, view_factor, internal_body_temp, enclosing_body_temp, unit):
    if unit == "metric":
        sigma = 5.67*  (10**(-8))  # Stefan-Boltzmann constant [W/(m²·K⁴)]
    if unit == "imperial":
        sigma = 1.713 * (10**(-9))  # Stefan-Boltzmann constant [BTU/(hr·ft²·R⁴)]
    Q = internal_surface_area * view_factor * sigma * (internal_body_temp**4 - enclosing_body_temp**4)
    return Q



"""
Calculates the net energy exchange by radiation between two diffuse-gray surfaces.

Parameters:
internal_area (float): Surface area of the first diffuse-gray surface (m^2)
external_area (float): Surface area of the second diffuse-gray surface (m^2)
internal_emissivity (float): Emissivity of the first diffuse-gray surface
external_emissivity (float): Emissivity of the second diffuse-gray surface
internal_temp (float): Absolute temperature of the first diffuse-gray surface (K)
external_temp (float): Absolute temperature of the second diffuse-gray surface (K)
view_factor (float): View factor between the two diffuse-gray surfaces
unit: metric or imperial

Returns:
Q (float): Net heat transfer rate between the two diffuse-gray surfaces (W)
"""
    
def net_heat_transfer_between_enclosed_gray_surfaces(internal_area, external_area, internal_emissivity, 
    external_emissivity, internal_temp, external_temp, view_factor, unit):
    if unit == "metric":
        sigma = 5.67*  (10**(-8))  # Stefan-Boltzmann constant [W/(m²·K⁴)]
    if unit == "imperial":
        sigma = 1.713 * (10**(-9))  # Stefan-Boltzmann constant [BTU/(hr·ft²·R⁴)]
    Q = (sigma * ((internal_temp**4)-(external_temp**4))) / \
        (((1-internal_emissivity) / (internal_emissivity * internal_area)) + 
         (1 / (internal_area * view_factor)) + 
         ((1 - external_emissivity) / (external_emissivity*external_area)))
    return Q



"""
    Calculates the net energy exchange by radiation between two parallel plates with a thin low-emissivity shield.
    
    Parameters:
    area_1 (float): Surface area of the first parallel plate (m^2)
    area_2 (float): Surface area of the second parallel plate (m^2)
    area_radiation_shield (float): Surface area of the radiation shield (m^2)
    emissivity_1 (float): Emissivity of the first parallel plate
    emissivity_2 (float): Emissivity of the second parallel plate
    emissivity_shield_to_a1 (float): Emissivity of the radiation shield towards the first parallel plate
    emissivity_shield_to_a2 (float): Emissivity of the radiation shield towards the second parallel plate
    temp_1 (float): Absolute temperature of the first parallel plate (K)
    temp_2 (float): Absolute temperature of the second parallel plate (K)
    temp_radition_shield (float): Absolute temperature of the radiation shield (K)
    view_factor_a1_to_shield (float): View factor from the first parallel plate to the radiation shield
    view_factor_a3_to_shield (float): View factor from the radiation shield to the second parallel plate
    unit: metric or imperial
    
    Returns:
    Q (float): Net heat transfer rate between the two parallel plates with the radiation shield (W)
    """
    
def net_heat_transfer_between_two_grey_bodies_with_radiation_shield(area_1, area_2, area_radiation_shield, emissivity_1, emissivity_2, 
    emissivity_shield_to_a1, emissivity_shield_to_a2, temp_1, temp_2, temp_3, view_factor_a1_to_shield, view_factor_a2_to_shield, unit):
    if unit == "metric":
        sigma = 5.67*  (10**(-8))  # Stefan-Boltzmann constant [W/(m²·K⁴)]
    if unit == "imperial":
        sigma = 1.713 * (10**(-9))  # Stefan-Boltzmann constant [BTU/(hr·ft²·R⁴)]
        
    Q = (sigma * (temp_1**4)-(temp_2**4)) / \
        (((1- emissivity_1) / (emissivity_1 *area_1)) + 
         (1 / (area_1 * view_factor_a1_to_shield))+ ((1- emissivity_shield_to_a1)/ (emissivity_shield_to_a1 * area_radiation_shield)) +
         ((1- emissivity_shield_to_a2)/ (emissivity_shield_to_a2 * area_radiation_shield)) +
         (1 / (area_radiation_shield * view_factor_a2_to_shield)) +((1- emissivity_2) / (emissivity_2 *area_2)))
    return Q



"""
Calculates the net energy exchange by radiation between two surfaces with a reradiating surface.

Parameters:
area_1 (float): Surface area of the first surface (m^2)
area_2 (float): Surface area of the second surface (m^2)
area_r (float): Surface area of the reradiating surface (m^2)
emissivity_1 (float): Emissivity of the first surface
emissivity_2 (float): Emissivity of the second surface
emissivity_r (float): Emissivity of the reradiating surface
temp_1 (float): Absolute temperature of the first surface (K)
temp_2 (float): Absolute temperature of the second surface (K)
temp_r (float): Absolute temperature of the reradiating surface (K)
view_factor_1_2 (float): View factor from the first surface to the second surface
view_factor_1_r (float): View factor from the first surface to the reradiating surface
view_factor_2_r (float): View factor from the second surface to the reradiating surface
unit (str): metric of imperial

Returns:
Q (float): Net heat transfer rate between the two surfaces with the reradiating surface (W)
"""
    
def net_heat_transfer_with_reradiating_surface(area_1, area_2, area_r, emissivity_1, emissivity_2, emissivity_r, temp_1, temp_2, temp_r, 
    view_factor_1_2, view_factor_1_r, view_factor_2_r, unit):
    
    if unit == "metric":
        sigma = 5.67*  (10**(-8))  # Stefan-Boltzmann constant [W/(m²·K⁴)]
    if unit == "imperial":
        sigma = 1.713 * (10**(-9))  # Stefan-Boltzmann constant [BTU/(hr·ft²·R⁴)]
        
    return (sigma * (temp_1**4)-(temp_2**4)) / \
        (((1- emissivity_1) / (emissivity_1 *area_1)) + ((1- emissivity_2) / (emissivity_2 *area_2))
         + (1 / ((area_1 * view_factor_1_r) + (((1/(area_1 * view_factor_1_r)) + (1/(area_2 * view_factor_2_r)))**(-1)))))