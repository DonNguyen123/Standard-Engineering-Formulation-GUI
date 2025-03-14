import numpy as np
import math


'''
=====================================
Standard Electrical Formulas
=====================================
'''

"""
Calculates the voltage (V), current (I), and resistance (R) using Ohm's Law.

Parameters:
voltage (float, optional): Voltage (V)
current (float, optional): Current (A)
resistance (float, optional): Resistance (Ω)

Returns:
Dictionary with calculated values based on provided parameters
"""

def ohms_law(voltage=None, current=None, resistance=None):
    # Validate inputs - exactly two parameters must be provided
    params = [voltage, current, resistance]
    if sum(p is not None for p in params) != 2:
        raise ValueError("Exactly two parameters must be provided")
    
    # Calculate the missing parameter
    if voltage is None:
        voltage = current * resistance
    elif current is None:
        current = voltage / resistance
    elif resistance is None:
        resistance = voltage / current
    
    # Create a formatted string representation of the results
    class FormattedDict(dict):
        def __str__(self):
            return '\n'.join(f'{key}: {value:_.4g}' if isinstance(value, (int, float)) 
                           else f'{key}: {value}'
                           for key, value in self.items())
        
        def __repr__(self):
            return self.__str__()

    return FormattedDict({
        'Voltage': voltage,
        'Current': current,
        'Resistance': resistance,
    })



"""
Calculates the power (P) in an electrical circuit.

Parameters:
voltage (float, optional): Voltage (V)
current (float, optional): Current (A)
resistance (float, optional): Resistance (Ω)
power (float, optional): Power (W)

Returns:
Dictionary with calculated values based on provided parameters
"""

def power_calculation(voltage=None, current=None, resistance=None, power=None):
    # Count provided parameters
    params = [voltage, current, resistance, power]
    provided = sum(p is not None for p in params)
    
    if provided < 1:
        raise ValueError("At least one parameter must be provided")
    
    # Calculate based on available parameters
    if power is None:
        if voltage is not None and current is not None:
            power = voltage * current
        elif voltage is not None and resistance is not None:
            power = (voltage ** 2) / resistance
        elif current is not None and resistance is not None:
            power = (current ** 2) * resistance
        else:
            raise ValueError("Insufficient parameters to calculate power")
    
    # Calculate other parameters if missing
    if voltage is None and current is not None and resistance is not None:
        voltage = current * resistance
    elif voltage is None and current is not None and power is not None:
        voltage = power / current
    elif voltage is None and resistance is not None and power is not None:
        voltage = math.sqrt(power * resistance)
        
    if current is None and voltage is not None and resistance is not None:
        current = voltage / resistance
    elif current is None and voltage is not None and power is not None:
        current = power / voltage
    elif current is None and resistance is not None and power is not None:
        current = math.sqrt(power / resistance)
        
    if resistance is None and voltage is not None and current is not None:
        resistance = voltage / current
    elif resistance is None and voltage is not None and power is not None:
        resistance = (voltage ** 2) / power
    elif resistance is None and current is not None and power is not None:
        resistance = power / (current ** 2)
    
    # Create a formatted string representation of the results
    class FormattedDict(dict):
        def __str__(self):
            return '\n'.join(f'{key}: {value:_.4g}' if isinstance(value, (int, float)) 
                           else f'{key}: {value}'
                           for key, value in self.items())
        
        def __repr__(self):
            return self.__str__()

    return FormattedDict({
        'Power': power,
        'Voltage': voltage,
        'Current': current,
        'Resistance': resistance,
    })




"""
Calculates parameters for resistors in series or parallel

Parameters:
resistors (list): List of resistance values (Ω)

Returns:
total_resistance (float): Total resistance (Ω)
current (float, optional): Current if voltage is provided (A)
voltage_drops (list, optional): Voltage drops across each resistor if voltage is provided (V)
configuration (str): 'series' or 'parallel'
"""

def resistence_network(resistors, voltage=None, configeration = None):
    if configeration == 'series':
        if not resistors:
            raise ValueError("At least one resistor must be provided")
        
        total_resistance = sum(resistors)
        result = {'Total Resistance': total_resistance}
        
        # If voltage is provided, calculate current and voltage drops
        if voltage is not None:
            current = voltage / total_resistance
            result['Current'] = current
            
            voltage_drops = [current * r for r in resistors]
            result['Voltage Drops'] = voltage_drops
        
        # Create a formatted string representation of the results
        class FormattedDict(dict):
            def __str__(self):
                return '\n'.join(f'{key}: {value:_.4g}' if isinstance(value, (int, float)) 
                            else f'{key}: {value}'
                            for key, value in self.items())
            
            def __repr__(self):
                return self.__str__()
        
        return FormattedDict(result)

    if configeration == 'parallel':
        if not resistors:
            raise ValueError("At least one resistor must be provided")
        
        # Calculate reciprocal of total resistance
        reciprocal_sum = sum(1/r for r in resistors)
        total_resistance = 1 / reciprocal_sum
        
        result = {'Total Resistance': total_resistance}
        
        # If voltage is provided, calculate currents through each resistor
        if voltage is not None:
            currents = [voltage / r for r in resistors]
            total_current = sum(currents)
            
            result['Individual Currents'] = currents
            result['Total Current'] = total_current
        
        # Create a formatted string representation of the results
        class FormattedDict(dict):
            def __str__(self):
                return '\n'.join(f'{key}: {value:_.4g}' if isinstance(value, (int, float)) 
                            else f'{key}: {value}'
                            for key, value in self.items())
            
            def __repr__(self):
                return self.__str__()
        
        return FormattedDict(result)



"""
Calculates capacitance in series and parallel configurations.

Parameters:
capacitors (list): List of capacitance values (F)
configuration (str): 'series' or 'parallel'

Returns:
total_capacitance (float): Total capacitance (F)
"""

def capacitor_network(capacitors, configuration):
    if not capacitors:
        raise ValueError("At least one capacitor must be provided")
    
    if configuration.lower() == 'series':
        # For series capacitors, reciprocal of sum of reciprocals
        reciprocal_sum = sum(1/c for c in capacitors)
        total_capacitance = 1 / reciprocal_sum
    elif configuration.lower() == 'parallel':
        # For parallel capacitors, sum of capacitances
        total_capacitance = sum(capacitors)
    else:
        raise ValueError("Configuration must be 'series' or 'parallel'")
    
    # Create a formatted string representation of the results
    class FormattedDict(dict):
        def __str__(self):
            return '\n'.join(f'{key}: {value:_.4g}' if isinstance(value, (int, float)) 
                           else f'{key}: {value}'
                           for key, value in self.items())
        
        def __repr__(self):
            return self.__str__()
    
    return FormattedDict({
        'Configuration': configuration,
        'Total Capacitance': total_capacitance
    })



"""
Calculates inductance in series and parallel configurations.

Parameters:
inductors (list): List of inductance values (H)
configuration (str): 'series' or 'parallel'

Returns:
total_inductance (float): Total inductance (H)
"""

def inductor_network(inductors, configuration):
    if not inductors:
        raise ValueError("At least one inductor must be provided")
    
    if configuration.lower() == 'series':
        # For series inductors, sum of inductances
        total_inductance = sum(inductors)
    elif configuration.lower() == 'parallel':
        # For parallel inductors, reciprocal of sum of reciprocals
        reciprocal_sum = sum(1/l for l in inductors)
        total_inductance = 1 / reciprocal_sum
    else:
        raise ValueError("Configuration must be 'series' or 'parallel'")
    
    # Create a formatted string representation of the results
    class FormattedDict(dict):
        def __str__(self):
            return '\n'.join(f'{key}: {value:_.4g}' if isinstance(value, (int, float)) 
                           else f'{key}: {value}'
                           for key, value in self.items())
        
        def __repr__(self):
            return self.__str__()
    
    return FormattedDict({
        'Configuration': configuration,
        'Total Inductance': total_inductance
    })



"""
Calculates parameters for an RC circuit.

Parameters:
resistance (float): Resistance (Ω)
capacitance (float): Capacitance (F)
frequency (float, optional): Frequency (Hz)

Returns:
time_constant (float): Time constant (s)
impedance (float, optional): Impedance if frequency is provided (Ω)
phase_angle (float, optional): Phase angle if frequency is provided (rad)
"""

def rc_circuit(resistance, capacitance, frequency=None):
    time_constant = resistance * capacitance
    result = {'Time Constant': time_constant}
    
    if frequency is not None:
        angular_frequency = 2 * math.pi * frequency
        capacitive_reactance = 1 / (angular_frequency * capacitance)
        impedance = math.sqrt(resistance**2 + capacitive_reactance**2)
        phase_angle = math.atan(capacitive_reactance / resistance)
        
        result['Capacitive Reactance'] = capacitive_reactance
        result['Impedance'] = impedance
        result['Phase Angle (rad)'] = phase_angle
        result['Phase Angle (deg)'] = math.degrees(phase_angle)
    
    # Create a formatted string representation of the results
    class FormattedDict(dict):
        def __str__(self):
            return '\n'.join(f'{key}: {value:_.4g}' if isinstance(value, (int, float)) 
                           else f'{key}: {value}'
                           for key, value in self.items())
        
        def __repr__(self):
            return self.__str__()
    
    return FormattedDict(result)



"""
Calculates parameters for an RL circuit.

Parameters:
resistance (float): Resistance (Ω)
inductance (float): Inductance (H)
frequency (float, optional): Frequency (Hz)

Returns:
time_constant (float): Time constant (s)
impedance (float, optional): Impedance if frequency is provided (Ω)
phase_angle (float, optional): Phase angle if frequency is provided (rad)
"""

def rl_circuit(resistance, inductance, frequency=None):
    time_constant = inductance / resistance
    result = {'Time Constant': time_constant}
    
    if frequency is not None:
        angular_frequency = 2 * math.pi * frequency
        inductive_reactance = angular_frequency * inductance
        impedance = math.sqrt(resistance**2 + inductive_reactance**2)
        phase_angle = math.atan(inductive_reactance / resistance)
        
        result['Inductive Reactance'] = inductive_reactance
        result['Impedance'] = impedance
        result['Phase Angle (rad)'] = phase_angle
        result['Phase Angle (deg)'] = math.degrees(phase_angle)
    
    # Create a formatted string representation of the results
    class FormattedDict(dict):
        def __str__(self):
            return '\n'.join(f'{key}: {value:_.4g}' if isinstance(value, (int, float)) 
                           else f'{key}: {value}'
                           for key, value in self.items())
        
        def __repr__(self):
            return self.__str__()
    
    return FormattedDict(result)



"""
Calculates parameters for an RLC circuit.

Parameters:
resistance (float): Resistance (Ω)
inductance (float): Inductance (H)
capacitance (float): Capacitance (F)
frequency (float, optional): Frequency (Hz)

Returns:
resonant_frequency (float): Resonant frequency (Hz)
quality_factor (float): Quality factor (dimensionless)
impedance (float, optional): Impedance if frequency is provided (Ω)
phase_angle (float, optional): Phase angle if frequency is provided (rad)
"""

def rlc_circuit(resistance, inductance, capacitance, frequency=None):
    resonant_frequency = 1 / (2 * math.pi * math.sqrt(inductance * capacitance))
    damping_factor = resistance / (2 * math.sqrt(inductance / capacitance))
    quality_factor = 1 / (2 * damping_factor)
    
    result = {
        'Resonant Frequency': resonant_frequency,
        'Damping Factor': damping_factor,
        'Quality Factor': quality_factor
    }
    
    if frequency is not None:
        angular_frequency = 2 * math.pi * frequency
        inductive_reactance = angular_frequency * inductance
        capacitive_reactance = 1 / (angular_frequency * capacitance)
        reactance = inductive_reactance - capacitive_reactance
        impedance = math.sqrt(resistance**2 + reactance**2)
        phase_angle = math.atan(reactance / resistance)
        
        result['Inductive Reactance'] = inductive_reactance
        result['Capacitive Reactance'] = capacitive_reactance
        result['Net Reactance'] = reactance
        result['Impedance'] = impedance
        result['Phase Angle (rad)'] = phase_angle
        result['Phase Angle (deg)'] = math.degrees(phase_angle)
    
    # Create a formatted string representation of the results
    class FormattedDict(dict):
        def __str__(self):
            return '\n'.join(f'{key}: {value:_.4g}' if isinstance(value, (int, float)) 
                           else f'{key}: {value}'
                           for key, value in self.items())
        
        def __repr__(self):
            return self.__str__()
    
    return FormattedDict(result)