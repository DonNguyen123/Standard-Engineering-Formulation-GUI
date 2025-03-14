import math
import numpy as np


    
'''
=====================================
Standard Material Properties
=====================================
'''

"""
Calculate the strength and elastic modulus of a reinforced plastic composite.

Parameters:
sigma_f (float): Strength of the fiber (MPa)
sigma_m (float): Strength of the matrix (MPa)
E_f (float): Elastic modulus of the fiber (GPa)
E_m (float): Elastic modulus of the matrix (GPa)
x (float): Volume fraction of fibers in the composite (dimensionless)

Returns:
dict: Dictionary containing:
    - sigma_c: Strength of the composite (MPa)
    - E_c: Elastic modulus of the composite (GPa)
    - P_f_ratio: Fraction of the total load carried by the fibers
"""

def composite_properties(sigma_f, sigma_m, E_f, E_m, x):
    # Calculate strength of the composite
    sigma_c = (x * sigma_f) + ((1 - x) * sigma_m)
    
    # Calculate elastic modulus of the composite
    E_c = x * E_f + ((1 - x) * E_m)
    
    # Calculate fraction of the total load carried by the fibers
    P_f_ratio = (x * E_f) / ((x * E_f) + ((1 - x) * E_m))
    
    class FormattedDict(dict):
        def __str__(self):
            return '\n'.join(f'{key}: {value:_.4g}' if isinstance(value, (int, float)) 
                        else f'{key}: {value}'
                        for key, value in self.items())
        
        def __repr__(self):
            return self.__str__()
    
    return FormattedDict({
        'sigma_c [MPa]': sigma_c,
        'E_c [GPa]': E_c,
        'P_f_ratio': P_f_ratio
    })



"""
Calculate the thermal conductivity of a ceramic material considering porosity.

Parameters:
k0 (float): Thermal conductivity of the ceramic at zero porosity (W/m·K)
P (float): Volume fraction of pores in the ceramic (dimensionless)

Returns:
float: Thermal conductivity of the ceramic (W/m·K)
"""

def thermal_conductivity_ceramic(k0, P):
    return k0 * (1 - P)



"""
Calculate the tensile strength and modulus of elasticity of a polycrystalline ceramic considering porosity.

Parameters:
UTS0 (float): Tensile strength at zero porosity (MPa)
E0 (float): Modulus of elasticity at zero porosity (GPa)
P (float): Volume fraction of pores in the ceramic (dimensionless)
n (float): Exponent ranging between 4 and 7 (dimensionless)

Returns:
dict: Dictionary containing:
    - UTS: Tensile strength of the ceramic (MPa)
    - E: Modulus of elasticity of the ceramic (GPa)
"""

def ceramic_properties(UTS0, E0, P, n):
    # Calculate tensile strength of the ceramic
    UTS = UTS0 * np.exp(-n * P)
    
    # Calculate modulus of elasticity of the ceramic
    E = E0 * (1 - 1.9 * P + 0.9 * P**2)
    
    class FormattedDict(dict):
        def __str__(self):
            return '\n'.join(f'{key}: {value:_.4g}' if isinstance(value, (int, float)) 
                        else f'{key}: {value}'
                        for key, value in self.items())
        
        def __repr__(self):
            return self.__str__()
    
    return FormattedDict({
        'UTS [MPa]': UTS,
        'E [GPa]': E
    })












'''
=====================================
Standard Rolling Properties
=====================================
'''

"""
Calculate the maximum possible draft in the flat-rolling process.

Parameters:
mu (float): Coefficient of friction (dimensionless)
R (float): Roll radius (meters)

Returns:
float: Maximum possible draft (meters)
"""

def maximum_draft(mu, R):
    return mu**2 * R



"""
Calculate the roll force in the flat-rolling process.

Parameters:
L (float): Roll-strip contact length (meters)
w (float): Width of the strip (meters)
Y_avg (float): Average true stress of the strip in the roll gap (Pa)
include_friction (bool): Whether to include friction (default=True)

Returns:
float: Roll force (Newtons)
"""

def roll_force(L, w, Y_avg, include_friction=True):
    F = L * w * Y_avg
    if include_friction:
        F = F * 1.2  # Increase force by 20% to account for friction
    return F



"""
Calculate the power required for the flat-rolling process.

Parameters:
F (float): Roll force (Newtons)
L (float): Roll-strip contact length (meters)
N (float): Revolutions per minute of the roll (rpm)
units (str): Units for power ('SI' for kW, 'English' for hp)

Returns:
float: Power required (kW or hp)
"""

def rolling_power(F, L, N, units='SI'):
    if units.lower() == 'si':
        power = (2 * np.pi * F * L * N) / 60000  # Power in kW
    elif units.lower() == 'english':
        power = (2 * np.pi * F * L * N) / 33000  # Power in hp
    else:
        raise ValueError("Units must be 'SI' or 'English'")
    return power



"""
Calculate the forward slip in the flat-rolling process.

Parameters:
V_f (float): Exit velocity of the strip (meters/second)
V_r (float): Surface speed of the roll (meters/second)

Returns:
float: Forward slip (dimensionless)
"""

def forward_slip(V_f, V_r):
    return (V_f - V_r) / V_r












'''
=====================================
Standard Forging, Extrusion, Drawing, Hubbing Properties
=====================================
'''

"""
Calculate the extrusion force for a given billet material and process conditions.

Parameters:
Y (float): Yield stress of the billet material (Pa)
alpha (float): Die angle (radians)
mu (float): Coefficient of friction (dimensionless)
R (float): Extrusion ratio (dimensionless)
A_o (float): Initial cross-sectional area of the billet (m²)
A_f (float): Final cross-sectional area of the billet (m²)
k (float): Extrusion Constant
alternate (str): Choice for alternate method of force

Returns:
float: Extrusion force (Newtons)
"""

def extrusion_force(Y, alpha, mu, R, A_o,A_f= None, alternate = "False", k= None):
    if alternate == "False":
        p = Y * (1 + (np.tan(alpha) / mu)) * ((R**(mu / (np.tan(alpha)))) - 1)
        F = p * A_o
        return F
    else:
        F = A_o*k*math.log(A_o/A_f)
        return F        



"""
Calculate the cold extrusion force.

Parameters:
A_o (float): Initial cross-sectional area of the slug (m²)
Y_avg (float): Average true stress of the material (Pa)
epsilon (float): True Strain (dimensionless)

Returns:
float: Cold extrusion force (Newtons)
"""

def cold_extrusion_force(A_o, Y_avg, epsilon):
    F = 1100 * A_o * Y_avg * epsilon
    return F



"""
Calculate the drawing force for wire drawing.

Parameters:
Y_avg (float): Average true stress of the material in the die gap (Pa)
A_f (float): Final cross-sectional area of the wire (m²)
A_o (float): Initial cross-sectional area of the wire (m²)
mu (float): Coefficient of friction (dimensionless)
alpha (float): Die angle (radians)

Returns:
float: Drawing force (Newtons)
"""

def drawing_force(Y_avg, A_f, A_o, mu, alpha):
    F = Y_avg * A_f * (((1 + mu / alpha) * np.log(A_o / A_f)) + ((2 / 3) * alpha))
    return F



"""
Calculate the hubbing force for creating a die cavity.

Parameters:
UTS (float): Ultimate tensile strength of the material (Pa)
A (float): Projected area of the impression (m²)

Returns:
float: Hubbing force (Newtons)
"""

def hubbing_force(UTS, A):
    F = 3 * UTS * A
    return F



"""
Calculate the forging force for open-die forging of a solid cylindrical workpiece.

Parameters:
Y_f (float): Flow stress of the material (Pa)
r (float): Instantaneous radius of the workpiece (m)
b (float): Instantaneous height of the workpiece (m)
mu (float): Coefficient of friction between the workpiece and the die (dimensionless)

Returns:
float: Forging force (Newtons)
"""

def forging_force(Y_f, r, b, mu):
    F = Y_f * np.pi * (r**2) * (1 + ((2 * mu * r) / (3 * b)))
    return F



"""
Calculate the die pressure distribution in upsetting of a cylinder (friction hill).

Parameters:
Y (float): Flow stress of the material (Pa)
r_o (float): Outer radius of the cylinder (m)
r (float): Radial position where pressure is calculated (m)
b (float): Height of the cylinder (m)
mu (float): Coefficient of friction (dimensionless)

Returns:
float: Die pressure at the given radial position (Pa)
"""

def die_pressure(Y, r_o, r, b, mu):
    p = Y * np.exp(2 * mu * (r_o - r) / b)
    return p












'''
=====================================
Standard Punching, Drawing, Forming Properties
=====================================
'''

"""
Calculate the punch force required to punch out a blank.

Parameters:
T (float): Sheet thickness (meters)
L (float): Total length being sheared (meters)
UTS (float): Ultimate tensile strength of the material (Pa)

Returns:
float: Punch force (Newtons)
"""

def punch_force(T, L, UTS):
    return 0.7 * T * L * UTS



"""
Calculate the dent resistance of a sheet metal panel.

Parameters:
Y (float): Yield stress of the material (Pa)
T (float): Sheet thickness (meters)
S (float): Panel stiffness (Pa·m³)

Returns:
float: Dent resistance (dimensionless)
"""

def dent_resistance(Y, T, S):
    return ((Y**2) * (T**4)) / S



"""
Calculate the bend allowance for a sheet metal part.

Parameters:
alpha (float): Bend angle (radians)
R (float): Bend radius (meters)
T (float): Sheet thickness (meters)
k (float): Constant (typically 0.33 for R < 2T, 0.5 for R > 2T)

Returns:
float: Bend allowance (meters)
"""

def bend_allowance(alpha, R, T, k):
    return alpha * (R + (k * T))



"""
Calculate the minimum bend radius to avoid cracking.

Parameters:
T (float): Sheet thickness (meters)
r (float): Tensile reduction of area (dimensionless)

Returns:
float: Minimum bend radius (meters)
"""

def minimum_bend_radius(T, r):
    return T * (50 / r - 1)



"""
Calculate the maximum bending force for sheets and plates.

Parameters:
k (float): Factor (0.3 for wiping die, 0.7 for U-die, 1.3 for V-die)
Y (float): Yield stress of the material (Pa)
L (float): Length of the bend (meters)
T (float): Sheet thickness (meters)
W (float): Die opening (meters)

Returns:
float: Bending force (Newtons)
"""

def maximum_bending_force(k, Y, L, T, W):
    return (k * Y * L * (T**2)) / W



"""
Calculate the limiting drawing ratio (LDR) for deep drawing.

Parameters:
D_o (float): Maximum blank diameter (meters)
D_p (float): Punch diameter (meters)

Returns:
float: Limiting drawing ratio (dimensionless)
"""

def limiting_drawing_ratio(D_o, D_p):
    return D_o / D_p



"""
Calculate the normal anisotropy (R) of a sheet metal.

Parameters:
epsilon_w (float): Width strain (dimensionless)
epsilon_t (float): Thickness strain (dimensionless)

Returns:
float: Normal anisotropy (dimensionless)
"""

def normal_anisotropy(epsilon_w, epsilon_t):
    return epsilon_w / epsilon_t



"""
Calculate the planar anisotropy (ΔR) of a sheet metal.

Parameters:
R_0 (float): R value at 0 degrees (dimensionless)
R_45 (float): R value at 45 degrees (dimensionless)
R_90 (float): R value at 90 degrees (dimensionless)

Returns:
float: Planar anisotropy (dimensionless)
"""

def planar_anisotropy(R_0, R_45, R_90):
    return (R_0 - 2 * R_45 + R_90) / 2












'''
=====================================
Standard Compaction and Machining Properties
=====================================
'''

"""
Calculate the cutting ratio (chip-thickness ratio) in machining.

Parameters:
t_o (float): Depth of cut (meters)
t_c (float): Chip thickness (meters)

Returns:
float: Cutting ratio (dimensionless)
"""

def cutting_ratio(t_o, t_c):
    return t_o / t_c



"""
Calculate the shear angle in machining.

Parameters:
alpha (float): Rake angle (degrees)
beta (float): Friction angle (degrees)

Returns:
float: Shear angle (degrees)
"""

def shear_angle(alpha, beta):
    phi = 45 + (alpha / 2) - (beta / 2)
    return phi



"""
Calculate the chip velocity in machining.

Parameters:
V (float): Cutting speed (meters/second)
r (float): Cutting ratio (dimensionless)

Returns:
float: Chip velocity (meters/second)
"""

def chip_velocity(V, r):
    return V * r



"""
Calculate the effective rake angle in oblique cutting.

Parameters:
i (float): Inclination angle (degrees)
alpha_n (float): Normal rake angle (degrees)

Returns:
float: Effective rake angle (degrees)
"""

def effective_rake_angle(i, alpha_n):
    alpha_e = np.arcsin((np.sin(np.radians(i))**2) + ((np.cos(np.radians(i))**2) * (np.sin(np.radians(alpha_n)))))
    return np.degrees(alpha_e)



"""
Calculate the shear force in machining.

Parameters:
F_c (float): Cutting force (Newtons)
F_t (float): Thrust force (Newtons)
phi (float): Shear angle (degrees)

Returns:
float: Shear force (Newtons)
"""

def shear_force(F_c, F_t, phi):
    return (F_c * np.cos(np.radians(phi))) - (F_t * np.sin(np.radians(phi)))



"""
Calculate the normal force in machining.

Parameters:
F_c (float): Cutting force (Newtons)
F_t (float): Thrust force (Newtons)
phi (float): Shear angle (degrees)

Returns:
float: Normal force (Newtons)
"""

def normal_force(F_c, F_t, phi):
    return (F_c * np.sin(np.radians(phi))) + (F_t * np.cos(np.radians(phi)))



"""
Calculate the coefficient of friction at the tool-chip interface.

Parameters:
F_c (float): Cutting force (Newtons)
F_t (float): Thrust force (Newtons)
alpha (float): Rake angle (degrees)

Returns:
float: Coefficient of friction (dimensionless)
"""

def coefficient_of_friction(F_c, F_t, alpha):
    return (F_t + (F_c * np.tan(np.radians(alpha)))) / (F_c - (F_t * np.tan(np.radians(alpha))))



"""
Calculate the power input in cutting.

Parameters:
F_c (float): Cutting force (Newtons)
V (float): Cutting speed (meters/second)

Returns:
float: Power input (Watts)
"""

def cutting_power(F_c, V):
    return F_c * V



"""
Calculate the specific energy for shearing in machining.

Parameters:
F_s (float): Shear force (Newtons)
V_s (float): Shear velocity (meters/second)
w (float): Width of cut (meters)
t_o (float): Depth of cut (meters)
V (float): Cutting speed (meters/second)

Returns:
float: Specific energy for shearing (Joules/m³)
"""

def specific_energy_shearing(F_s, V_s, w, t_o, V):
    return (F_s * V_s) / (w * t_o * V)



"""
Calculate the specific energy for friction in machining.

Parameters:
F (float): Friction force (Newtons)
r (float): cutting ratio (unitless)
w (float): Width of cut (meters)
t_o (float): Depth of cut (meters)

Returns:
float: Specific energy for friction (Joules/m³)
"""

def specific_energy_friction(F, r, w, t_o):
    return (F * r) / (w * t_o)



"""
Calculate the total specific energy in machining.

Parameters:
u_s (float): Specific energy for shearing (Joules/m³)
u_f (float): Specific energy for friction (Joules/m³)

Returns:
float: Total specific energy (Joules/m³)
"""

def total_specific_energy(u_s, u_f):
    return u_s + u_f



"""
Calculate the peak pressure generated in explosive forming.

Parameters:
K (float): Constant (dimensionless)
W (float): Explosive charge weight (kilograms)
R (float): Distance from explosive to workpiece (meters)
a (float): Exponent (dimensionless)

Returns:
float: Peak pressure (Pascals)
"""

def explosive_forming_pressure(K, W, R, a):
    return K * (((W**(1/3)) / R)**a)



"""
Calculate the pressure distribution during compaction.

Parameters:
p_o (float): Pressure at the punch (Pascals)
mu (float): Coefficient of friction (dimensionless)
k (float): Interparticle friction factor (dimensionless)
x (float): Distance along the compact (meters)
D (float): Compact diameter (meters)

Returns:
float: Pressure at distance x (Pascals)
"""

def compaction_pressure(p_o, mu, k, x, D):
    return p_o * np.exp(-4 * mu * k * x / D)












'''
=====================================
Standard Welding, Grinding, Milling Properties
=====================================
'''

"""
Calculate the heat input in arc welding.

Parameters:
V (float): Voltage applied (Volts)
I (float): Current (Amperes)
nu (float): Welding speed (meters/second)
e (float): Efficiency of the process (dimensionless)

Returns:
float: Heat input (Joules/length)
"""

def heat_input_welding(V, I, nu, e):
    return e * V * I / nu



"""
Calculate the welding speed in arc welding.

Parameters:
V (float): Voltage applied (Volts)
I (float): Current (Amperes)
u (float): Specific energy required for melting (Joules/m³)
A (float): Cross-section of the weld (m²)
e (float): Efficiency of the process (dimensionless)

Returns:
float: Welding speed (meters/second)
"""

def welding_speed(V, I, u, A, e):
    return e * V * I / (u * A)



"""
Calculate the cutting depth in laser cutting.

Parameters:
C (float): Constant for the process (dimensionless)
P (float): Power input (Watts)
nu (float): Cutting speed (meters/second)
d (float): Laser-spot diameter (meters)

Returns:
float: Cutting depth (meters)
"""

def laser_cutting_depth(C, P, nu, d):
    return C * P / (nu * d)



"""
Calculate the material removal rate in electrical-discharge machining (EDM).

Parameters:
I (float): Current (Amperes)
T_w (float): Melting point of the workpiece (°C)

Returns:
float: Material removal rate (mm³/min)
"""

def edm_material_removal_rate(I, T_w):
    return 4*(10**4) * I *(T_w**(-1.23))



"""
Calculate the undeformed chip length in surface grinding.

Parameters:
D (float): Wheel diameter (meters)
d (float): Wheel depth of cut (meters)

Returns:
float: Undeformed chip length (meters)
"""

def undeformed_chip_length(D, d):
    return np.sqrt(D * d)



"""
Calculate the undeformed chip thickness in surface grinding.

Parameters:
v (float): Workpiece velocity (meters/second)
V (float): Wheel tangential velocity (meters/second)
C (float): Number of cutting points per unit area (1/m²)
r (float): Ratio of chip width to average undeformed chip thickness (dimensionless)
d (float): Wheel depth of cut (meters)
D (float): Wheel diameter (meters)

Returns:
float: Undeformed chip thickness (meters)
"""

def undeformed_chip_thickness(v, V, C, r, d, D):
    return np.sqrt(((4 * v) / (V * C * r)) * (math.sqrt(d / D)))



"""
Calculate the cutting speed in peripheral milling.

Parameters:
D (float): Cutter diameter (meters)
N (float): Rotational speed of the cutter (revolutions/second)

Returns:
float: Cutting speed (meters/second)
"""

def milling_cutting_speed(D, N):
    return np.pi * D * N



"""
Calculate the undeformed chip thickness in slab milling.

Parameters:
f (float): Feed per tooth of the cutter (meters/tooth)
d (float): Depth of cut (meters)
D (float): Cutter diameter (meters)

Returns:
float: Undeformed chip thickness (meters)
"""

def milling_undeformed_chip_thickness(f, d, D):
    return 2 * f * np.sqrt(d / D)



"""
Calculate the feed per tooth in milling.

Parameters:
nu (float): Linear speed of the workpiece (meters/second)
N (float): Rotational speed of the cutter (revolutions/second)
n (float): Number of teeth on the cutter (dimensionless)

Returns:
float: Feed per tooth (meters/tooth)
"""

def milling_feed_per_tooth(nu, N, n):
    return nu / (N * n)



"""
Calculate the cutting time in milling.

Parameters:
l (float): Length of the workpiece (meters)
l_c (float): Horizontal extent of the cutter's first contact with the workpiece (meters)
nu (float): Linear speed of the workpiece (meters/second)

Returns:
float: Cutting time (seconds)
"""

def milling_cutting_time(l, l_c, nu):
    return (l + l_c) / nu



"""
Calculate the material removal rate in milling.

Parameters:
w (float): Width of the workpiece (meters)
d (float): Depth of cut (meters)
nu (float): Linear speed of the workpiece (meters/second)

Returns:
float: Material removal rate (m³/second)
"""
def milling_material_removal_rate(w, d, nu):
    return w * d * nu












'''
=====================================
Standard Milling Economic Properties
=====================================
'''

"""
Calculate the optimum cutting speed in machining economics.

Parameters:
C (float): Constant (dimensionless)
n (float): Taylor's tool life exponent (dimensionless)
T_c (float): Tool change time (seconds)
m (float): Machine cost rate (dollars/second)
T_i (float): Tool indexing time (seconds)

Returns:
float: Optimum cutting speed (meters/second)
"""

def optimum_cutting_speed(C, n, T_c, m, T_i):
    return C / ((1/n - 1) * (T_c / m + T_i))**n



"""
Calculate the optimum tool life in machining economics.

Parameters:
n (float): Taylor's tool life exponent (dimensionless)
T_c (float): Tool change time (seconds)
m (float): Machine cost rate (dollars/second)
T_i (float): Tool indexing time (seconds)

Returns:
float: Optimum tool life (seconds)
"""

def optimum_tool_life(n, T_c, m, T_i):
    return (1/n - 1) * (T_c / m + T_i)