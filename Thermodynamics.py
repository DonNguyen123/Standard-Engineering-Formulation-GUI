import pyromat as pm
import numpy as np
from typing import Union, List, Dict, Tuple

import pyromat as pm
print(pm.__path__)
'''
=====================================
Standard Thermodynamic Properities
=====================================
'''

"""
Calculate thermodynamic properties of various substances using PYroMat library.

Parameters:
substance (str): Substance identifier. See PYroMat documentation for full list of available substances. Examples:
- 'mp.H2O' for water
- 'mp.C2H6' for ethane
- 'mp.N2' for nitrogen
- 'mp.O2' for oxygen
- 'mp.CO2' for carbon dioxide

T (float): Temperature in Celsius
P (float): Pressure in kPa (optional, defaults to 101.325 kPa)

Returns:
dict: Dictionary containing calculated properties:
- 'T': Temperature (°C)
- 'P': Pressure (kPa)
- 'h': Specific enthalpy (kJ/kg)
- 'h_vap': Specific vapor enthalpy (kJ/kg)
- 'v': Specific volume (m³/kg)
- 'v_vap': Specific vapor volume (m³/kg)
- 's': Specific entropy (kJ/(kg·K))
- 's_vap': Specific vapor entropy (kJ/(kg·K))
- 'cp': Specific heat capacity (kJ/(kg·K))
- 'cv': Specific heat capacity at constant volume (kJ/(kg·K))
- 'rho': Density (kg/m³)
- 'phase': Phase of substance ('liquid', 'vapor', or 'supercritical')
- 'M': Molar mass (kg/kmol)
- 'Z': Compressibility factor
"""
    
def thermo_properties(substance, T=None, P=None):
    pm.config['unit_pressure'] = 'kPa'    # Ensure pressure is in kPa
    pm.config['unit_temperature'] = 'K'   # Use Kelvin internally
    try:
        # Initialize substance from PYroMat
        subst = pm.get(substance)
    except Exception as e:
        raise ValueError(f"Invalid substance identifier: {substance}. Error: {str(e)}")
    
    # Set default pressure if not provided
    if P is None:
        P = 101.325
        
    # Get critical properties if available
    try:
        T_crit = subst.critical('T') - 273.15  # Convert to Celsius
        P_crit = subst.critical('P')
    except:
        T_crit = None
        P_crit = None
        
    
    if T is None:
        raise ValueError("Temperature must be provided for forward mode")
        
    # Convert temperature to Kelvin for PYroMat
    T_K = T + 273.15
    
    try:
        # Calculate basic properties
        h = subst.h(T_K, P)  # kJ/kg
        s = subst.s(T_K, P)  # kJ/kg·K
        v = 1/subst.d(T_K, P)  # m³/kg
        cp = subst.cp(T_K, P)  # kJ/kg·K
        cv = subst.cv(T_K, P)  # kJ/kg·K
        rho = subst.d(T_K, P)  # kg/m³
        
        # Get molar mass
        M = subst.mw()  # kg/kmol
        
        # Calculate compressibility factor
        Z = (P * v * M) / (8.314 * T_K)
        
        # Determine phase if critical points are available
        if T_crit is not None and P_crit is not None:
            if T > T_crit or P > P_crit:
                phase = 'supercritical'
            else:
                try:
                    # Use saturation pressure to determine phase
                    P_sat = subst.ps(T_K)
                    if abs(P - P_sat) < 1e-3:
                        phase = 'saturated'
                    elif P > P_sat:
                        phase = 'liquid'
                    else:
                        phase = 'vapor'
                except:
                    phase = 'unknown'
        else:
            phase = 'unknown'
                
    except Exception as e:
        raise ValueError(f"Error calculating properties: {str(e)}")
    
    
    try:
        # Saturation temperature (liquid and vapor are defined at the same temperature)
        if hasattr(subst, 'Ts'):
            T_sat = T + 273.15  # Use input temperature as baseline
            pm.config['unit_temperature'] = 'K' 
            # Saturated vapor properties
            h_vap = subst.h(T_sat, x=1)
            s_vap = subst.s(T_sat, x=1)
            v_vap = 1 / subst.d(T_sat, x=1)

        else:
            h_vap = s_vap = v_vap = None

    except Exception as e:
        raise ValueError(f"Error calculating saturated properties: {str(e)}")
    
    # Create a formatted string representation of the results
    class FormattedDict(dict):
        def __str__(self):
            return '\n'.join(f'{key}: {value:_.4g}' if isinstance(value, (int, float)) 
                           else f'{key}: {value}'
                           for key, value in self.items())
        
        def __repr__(self):
            return self.__str__()

    return FormattedDict({
        'Temperature [°C]': T,
        'Pressure [kPa]': P,
        'Specific Enthalpy [kJ/kg]': h,
        'Saturated Vapor Enthalpy [kJ/kg]': h_vap,
        'Specific Entropy [kJ/kg·K]': s,
        'Saturated Vapor Entropy [kJ/kg·K]': s_vap,
        'Specific Volume [m³/kg]': v,
        'Saturated Vapor Volume [m³/kg]': v_vap,
        'Specific Heat cp [kJ/kg·K]': cp,
        'Specific Heat cv [kJ/kg·K]': cv,
        'Density [kg/m³]': rho,
        'Phase': phase,
        'Molar Mass [kg/kmol]': M,
        'Compressibility Factor': Z
    })
    
    

"""
Calculate thermal efficiency of a heat engine.

Parameters:
work_out (float): Work output of the system (kJ)
heat_in (float): Heat input to the system (kJ)

Returns:
(float): Thermal efficiency as a decimal
"""
def calculate_efficiency(work_out, heat_in):
    return work_out / heat_in



"""
Calculate the theoretical maximum (Carnot) efficiency of a heat engine.

Parameters:
T_hot (float): Temperature of hot reservoir
T_cold (float): Temperature of cold reservoir
temp_unit (str): Temperature unit ('C', 'K', or 'F')

Returns:
(float): Carnot efficiency as a decimal
"""
def carnot_efficiency(T_hot, T_cold, temp_unit):
    # Convert temperatures to Kelvin
    if temp_unit.upper() == 'C':
        T_hot_K = T_hot + 273.15
        T_cold_K = T_cold + 273.15
    elif temp_unit.upper() == 'F':
        T_hot_K = (T_hot - 32) * 5/9 + 273.15
        T_cold_K = (T_cold - 32) * 5/9 + 273.15
    elif temp_unit.upper() == 'K':
        T_hot_K = T_hot
        T_cold_K = T_cold
    else:
        raise ValueError("Temperature unit must be 'C', 'K', or 'F'")
    
    return 1 - (T_cold_K / T_hot_K)



"""
Calculate isentropic efficiency of a turbine or compressor.

Parameters:
actual_work (float): Actual work done (kJ)
ideal_work (float): Ideal (isentropic) work (kJ)

Returns:
(float): Isentropic efficiency as a decimal
"""

def isentropic_efficiency(actual_work, ideal_work):
    return actual_work / ideal_work



"""
Calculate the effectiveness of a heat exchanger considering heat capacity rates.

Parameters:
m_dot_hot (float): Mass flow rate of hot fluid (kg/s)
cp_hot (float): Specific heat capacity of hot fluid (kJ/kg·K)
T_hot_in (float): Hot fluid inlet temperature (°C)
T_hot_out (float): Hot fluid outlet temperature (°C)
m_dot_cold (float): Mass flow rate of cold fluid (kg/s)
cp_cold (float): Specific heat capacity of cold fluid (kJ/kg·K)
T_cold_in (float): Cold fluid inlet temperature (°C)
T_cold_out (float): Cold fluid outlet temperature (°C)
exchanger_type (str): Type of heat exchanger ('parallel', 'counter', 'shell-tube')

Returns:
dict: Dictionary containing:
    - effectiveness: Heat exchanger effectiveness
    - C_min: Minimum heat capacity rate
    - C_max: Maximum heat capacity rate
    - C_ratio: Heat capacity rate ratio
    - NTU: Number of transfer units (if UA is provided)
    - actual_heat_transfer: Actual heat transfer rate
    - max_possible_heat_transfer: Maximum possible heat transfer rate
"""

def heat_exchanger_effectiveness(m_dot_hot,cp_hot,T_hot_in,T_hot_out,m_dot_cold,cp_cold,T_cold_in,T_cold_out, exchanger_type):
    # Calculate heat capacity rates
    C_hot = m_dot_hot * cp_hot    # kW/K
    C_cold = m_dot_cold * cp_cold  # kW/K
    
    # Determine C_min and C_max
    C_min = min(C_hot, C_cold)
    C_max = max(C_hot, C_cold)
    C_ratio = C_min / C_max
    
    # Calculate actual heat transfer rate
    Q_actual_hot = C_hot * (T_hot_in - T_hot_out)
    Q_actual_cold = C_cold * (T_cold_out - T_cold_in)
    
    # Use average of hot and cold side calculations to account for measurement uncertainties
    Q_actual = (Q_actual_hot + Q_actual_cold) / 2
    
    # Calculate maximum possible heat transfer rate
    Q_max = C_min * (T_hot_in - T_cold_in)
    
    # Calculate effectiveness
    effectiveness = Q_actual / Q_max
    
    # Create result dictionary
    results = {
        'effectiveness': effectiveness,
        'C_min': C_min,
        'C_max': C_max,
        'C_ratio': C_ratio,
        'actual_heat_transfer': Q_actual,
        'max_possible_heat_transfer': Q_max
    }
    
    # Additional calculations based on exchanger type
    if exchanger_type.lower() == 'parallel':
        # Theoretical effectiveness for parallel flow
        # ε = (1 - exp(-NTU*(1 + C_ratio))) / (1 + C_ratio)
        results['theoretical_max_effectiveness'] = 1 / (1 + C_ratio)
    
    elif exchanger_type.lower() == 'counter':
        # Theoretical effectiveness for counter flow
        # ε = (1 - exp(-NTU*(1 - C_ratio))) / (1 - C_ratio*exp(-NTU*(1 - C_ratio)))
        if C_ratio < 1:
            results['theoretical_max_effectiveness'] = 1
        else:
            results['theoretical_max_effectiveness'] = C_ratio
    
    elif exchanger_type.lower() == 'shell-tube':
        # For shell and tube heat exchangers, effectiveness depends on number of shell passes
        # This is a simplified approximation for one shell pass
        results['theoretical_max_effectiveness'] = 0.8  # Typical maximum for one shell pass

    class FormattedDict(dict):
        def __str__(self):
            return '\n'.join(f'{key}: {value:_.4g}' if isinstance(value, (int, float)) 
                           else f'{key}: {value}'
                           for key, value in self.items())
        
        def __repr__(self):
            return self.__str__()
        
    return FormattedDict(results)



"""
Calculate pressure ratio between two states.

Parameters:
P2 (float): Final pressure (kPa)
P1 (float): Initial pressure (kPa)

Returns:
(float): Pressure ratio (dimensionless)
"""
def pressure_ratio(P2, P1):
    return P2 / P1



"""
Calculate specific work required for compression process.

Parameters:
substance (str): Substance identifier (e.g., 'mp.Air')
T1 (float): Initial temperature (°C)
P1 (float): Initial pressure (kPa)
P2 (float): Final pressure (kPa)
eta_isentropic (float): Isentropic efficiency (default=1.0 for ideal case)

Returns:
dict: Dictionary containing work and final temperature
"""

def specific_work_compressor(substance,T1,P1,P2,eta_isentropic):

    # Get initial state properties
    state1 = thermo_properties(substance, T1, P1)
    h1 = state1['Specific Enthalpy [kJ/kg]']
    s1 = state1['Specific Entropy [kJ/kg·K]']
    
    # Find temperature at state 2s (isentropic process)
    def entropy_diff(T2s):
        state2s = thermo_properties(substance, T2s, P2)
        return state2s['Specific Entropy [kJ/kg·K]'] - s1
    
    # Use numerical method to find T2s where entropy matches
    from scipy.optimize import fsolve
    T2s = fsolve(entropy_diff, T1)[0]
    
    # Calculate actual final state
    state2s = thermo_properties(substance, T2s, P2)
    h2s = state2s['Specific Entropy [kJ/kg·K]']
    
    # Calculate actual enthalpy change using isentropic efficiency
    delta_h_actual = (h2s - h1) / eta_isentropic
    
    # Calculate actual final temperature
    h2_actual = h1 + delta_h_actual
    
    def enthalpy_diff(T2):
        state2 = thermo_properties(substance, T2, P2)
        return state2['Specific Enthalpy [kJ/kg]'] - h2_actual
    
    T2_actual = fsolve(enthalpy_diff, T2s)[0]
    
        # Create a formatted string representation of the results
    class FormattedDict(dict):
        def __str__(self):
            return '\n'.join(f'{key}: {value:_.4g}' if isinstance(value, (int, float)) 
                           else f'{key}: {value}'
                           for key, value in self.items())
        
        def __repr__(self):
            return self.__str__()
    
    return FormattedDict({
        'specific_work': delta_h_actual,
        'final_temperature': T2_actual,
        'isentropic_final_temperature': T2s
    })



"""
Calculate available energy (exergy) of a substance.

Parameters:
substance (str): Substance identifier
T (float): Temperature (°C)
P (float): Pressure (kPa)
T_ambient (float): Ambient temperature (°C)
P_ambient (float): Ambient pressure (kPa)

Returns:
float: Specific available energy (kJ/kg)
"""

def available_energy(substance,T,P,T_ambient,P_ambient):
    # Get properties at given state
    state = thermo_properties(substance, T, P)
    h = state['Specific Enthalpy [kJ/kg]']
    s = state['Specific Entropy [kJ/kg·K]']
    
    # Get properties at ambient conditions
    state_ambient = thermo_properties(substance, T_ambient, P_ambient)
    h0 = state_ambient['Specific Enthalpy [kJ/kg]']
    s0 = state_ambient['Specific Entropy [kJ/kg·K]']
    
    # Calculate specific available energy
    T0_K = T_ambient + 273.15
    return (h - h0) - (T0_K * (s - s0))