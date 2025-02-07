import numpy as np
from sympy import Symbol
from sympy.abc import s, t
import sympy as sp
from typing import Union, List, Tuple

from scipy import signal
import matplotlib.pyplot as plt
import control
from control import matlab
import re


'''
=====================================
Standard Control Functions 
=====================================
'''

"""
Analyze system stability using Routh-Hurwitz criterion.

Parameters: List of coefficients in descending order of s. Example: s^3 + 2s^2 + 3s + 4 would be [1, 2, 3, 4]. Can be numerical or symbolic (using SymPy symbols).

Returns Tuple containing:
- bool: True if system is stable, False otherwise
- str: Explanation of stability analysis
- ndarray: The computed Routh array

Example usage:
if __name__ == "__main__":
    # Numerical example
    coeffs_numerical = [1, 2, 3, 4]  # s^3 + 2s^2 + 3s + 4
    stable, explanation, array = routh_hurwitz(coeffs_numerical)
    print("Numerical Example:")
    print(format_routh_array(array))
    print("\nStability Analysis:")
    print(explanation)

    # Symbolic example
    from sympy.abc import a, b, c
    coeffs_symbolic = [1, a, b, c]  # s^3 + as^2 + bs + c
    stable, explanation, array = routh_hurwitz(coeffs_symbolic)
    print("\nSymbolic Example:")
    print(format_routh_array(array))
    print("\nStability Analysis:")
    print(explanation)
    """

def routh_hurwitz(coefficients: List[Union[float, Symbol]]) -> Tuple[bool, str, np.ndarray]:
    order = len(coefficients) - 1
    n = order + 1
    
    # Initialize Routh array
    from sympy import symbols
    array = np.zeros(((n), (n + 1) // 2), dtype=object)
    
    # Fill first two rows
    for i in range(0, (n + 1) // 2):
        if 2 * i < n:
            array[0, i] = coefficients[2 * i]
        if 2 * i + 1 < n:
            array[1, i] = coefficients[2 * i + 1]
    
    # Calculate remaining rows
    explanation = []
    for i in range(2, n):
        for j in range((n + 1) // 2 - 1):
            if array[i-1, 0] == 0:
                # Special case: when first element is zero
                degree = n - i + 1
                array[i-1, 0] = degree * array[i-2, 0]
                explanation.append(f"Row {i-1}: First element was zero, using auxiliary polynomial")
            
            determinant = (array[i-1, 0] * array[i-2, j+1] - 
                         array[i-2, 0] * array[i-1, j+1])
            array[i, j] = determinant / array[i-1, 0]
    
    # Check stability by analyzing first column
    first_column = array[:, 0]
    sign_changes = 0
    prev_sign = None
    
    for i, value in enumerate(first_column):
        if value == 0:
            explanation.append(f"Row {i} has a zero in first column - system may have pure imaginary roots")
            return False, "\n".join(explanation), array
        
        try:
            # Handle both numeric and symbolic cases
            if isinstance(value, (int, float)):
                current_sign = 1 if value > 0 else -1
            else:
                # For symbolic expressions, we can only note that stability depends on the sign
                explanation.append(f"Row {i} contains symbolic term {value} - stability depends on its sign")
                continue
                
            if prev_sign is not None and current_sign != prev_sign:
                sign_changes += 1
                explanation.append(f"Sign change detected at row {i}: from {prev_sign} to {current_sign}")
            prev_sign = current_sign
        except Exception as e:
            explanation.append(f"Could not determine sign of term in row {i}: {value}")
            return False, "\n".join(explanation), array
    
    is_stable = sign_changes == 0
    if is_stable:
        explanation.append("System is stable - no sign changes in first column")
    else:
        explanation.append(f"System is unstable - {sign_changes} sign changes detected")
    
    return is_stable, "\n".join(explanation), array

def format_routh_array(array: np.ndarray) -> str:
    """Helper function to format the Routh array for display"""
    order = len(array) - 1
    result = []
    for i, row in enumerate(array):
        row_str = f"s^{order - i}: "
        row_str += " ".join(str(x) for x in row if x != 0)
        result.append(row_str)
    return "\n".join(result)




"""
Simulate and plot PID system response to various input types.

Parameters:
-----------
Kp (float): Proportional gain
Ki (float): Intergal gain
Kd (float): Derivative gain
input_type (str): Type of input signal ('step', 'ramp', 'impulse', or 'points')
t_span (float): Simulation time span
num_points (int): Number of simulation points
custom_points (list): List of custom input points (used when input_type='points')
amplitude (float): Amplitude of the input signal
    
Returns:
t (array): Time points
y (array): System response
input_signal (array): Input signal values
    
Example usage:
if __name__ == "__main__":
    # Step response
    t, y, input_signal = simulate_pid(Kp=10, Ki=0.1, Kd=1, input_type='step')

    # Ramp response
    # t, y, input_signal = simulate_pid(Kp=1.0, Ki=0.1, Kd=0.05, input_type='ramp')

    # Impulse response
    #t, y, input_signal = simulate_pid(Kp=1.0, Ki=0.1, Kd=0.05, input_type='impulse')

    # Custom points response
    # points = [0, 0.5, 1.0, 0.8, 0.6, 0.4]
    # t, y, input_signal = simulate_pid(Kp=1.0, Ki=0.1, Kd=0.05, 
    #                                  input_type='points', custom_points=points)

    """
    
def simulate_pid(Kp, Ki, Kd, input_type='step', t_span=30, num_points=1000, 
                custom_points=None, amplitude=1.0):
     # Initialize time array
    t = np.linspace(0, t_span, num_points)
    dt = t_span / (num_points - 1)
    
    # Define input signal based on type
    if input_type == 'step':
        input_values = np.where(t > 0, amplitude, 0)
    elif input_type == 'ramp':
        input_values = amplitude * t
    elif input_type == 'impulse':
        input_values = np.where(t < 0.1, amplitude, 0)
    elif input_type == 'points':
        if custom_points is None:
            raise ValueError("Must provide custom_points for 'points' input type")
        input_values = np.interp(t, np.linspace(0, t_span, len(custom_points)), custom_points)
    else:
        raise ValueError("Invalid input_type. Must be 'step', 'ramp', 'impulse', or 'points'")

    def pid_system(state, t):
        """
        State space representation of the system
        state[0] = position
        state[1] = velocity
        state[2] = integral of error
        """
        position, velocity, error_integral = state
        
        # Get desired position at current time
        target = np.interp(t, np.linspace(0, t_span, len(input_values)), input_values)
        
        # Calculate error
        error = target - position
        
        # PID control law
        control = (Kp * error +                  # Proportional term
                  Ki * error_integral +          # Integral term
                  Kd * (-velocity))              # Derivative term
        
        # State derivatives
        d_position = velocity
        d_velocity = control  # Assuming unit mass
        d_error_integral = error
        
        return [d_position, d_velocity, d_error_integral]

    # Initial conditions: [position, velocity, error_integral]
    initial_state = [0, 0, 0]
    
    # Simulate system response
    from scipy.integrate import odeint
    solution = odeint(pid_system, initial_state, t)
    position = solution[:, 0]
    
    # Plot results
    plt.figure(figsize=(12, 6))
    plt.plot(t, input_values, 'r--', label='Input (Setpoint)')
    plt.plot(t, position, 'b-', label='System Response')
    plt.grid(True)
    plt.xlabel('Time')
    plt.ylabel('Amplitude')
    plt.title(f'PID System Response (Kp={Kp}, Ki={Ki}, Kd={Kd})')
    plt.legend()
    plt.show()
    
    return t, position, input_values



"""
Convert s-domain expression to time domain, analyze, and plot the result.

Parameters:
expr_str (str): String representation of the s-domain expression (e.g., "1/(s^2 + 2*s + 1)")
t_start (float): Start time for plotting
t_end (float): End time for plotting
num_points (int): Number of points to plot

Returns:
dict: Dictionary containing the analysis results

# Example usage
if __name__ == "__main__":
    # Example 1: Second order system
    print("\nExample 1: Second order system")
    result1 = analyze_and_plot_laplace("1/(s**2 + 2*s + 1)")
    print(f"Time domain expression: {result1['time_domain_expr']}")
    print(f"Poles: {result1['poles']}")
    print(f"Zeros: {result1['zeros']}")
    
    # Example 2: System with zeros
    print("\nExample 2: System with zeros")
    result2 = analyze_and_plot_laplace("(s + 1)/(s**2 + 2*s + 2)")
    print(f"Time domain expression: {result2['time_domain_expr']}")
    print(f"Poles: {result2['poles']}")
    print(f"Zeros: {result2['zeros']}")
    """
    
def analyze_and_plot_time_from_laplace(expr_str, t_start=0, t_end=10, num_points=1000):
    try:
        # Initialize printing
        from sympy import init_printing
        init_printing(use_unicode=True)
        
        # Convert string to sympy expression
        s_expr = sp.sympify(expr_str)
        
        # Get the time domain expression
        from sympy import inverse_laplace_transform
        time_expr = inverse_laplace_transform(s_expr, s, t)
        
        # Find poles
        # Convert to rational function first
        num, den = sp.fraction(s_expr)
        poles = sp.solve(den, s)
        zeros = sp.solve(num, s)
        
        # Create lambda function for numerical evaluation
        time_fn = sp.lambdify(t, time_expr, modules=['numpy'])
        
        # Generate time points for plotting
        t_points = np.linspace(t_start, t_end, num_points)
        
        try:
            # Evaluate the function
            y_points = time_fn(t_points)
            
            # Create the plot
            plt.figure(figsize=(12, 8))
            
            # Plot the time response
            plt.subplot(2, 1, 1)
            plt.plot(t_points, y_points, 'b-', label='Time Response')
            plt.grid(True)
            plt.xlabel('Time (s)')
            plt.ylabel('Amplitude')
            plt.title(f'Time Domain Response: {time_expr}')
            plt.legend()
            
            plt.tight_layout()
            plt.show()
            
        except Exception as e:
            print(f"Error in plotting: {str(e)}")
        
        return {
            "time_domain_expr": time_expr,
            "poles": poles,
            "zeros": zeros
        }
        
    except Exception as e:
        return f"Error processing expression: {str(e)}"
    
    
    
"""
Plot the root locus of a transfer function given as a string expression.

Parameters:
expression (str): Transfer function expression (e.g., "1/(s^2 + 2*s + 1)")
title (str): Optional plot title
xlim (tuple): Optional x-axis limits (min, max)
ylim (tuple): Optional y-axis limits (min, max)

Example usage:
plot_root_locus("1/(s^2 + 2*s + 1)")
plot_root_locus("(s + 1)/(s^3 + 4*s^2 + 5*s + 2)")
"""
    
def plot_root_locus(expression, title=None, xlim=None, ylim=None):
    # Clean the expression and standardize the format
    expression = expression.replace(' ', '')
    expression = expression.replace('**', '^')
    
    # Split into numerator and denominator
    if '/' in expression:
        num_str, den_str = expression.split('/')
        num_str = num_str.strip('()')
        den_str = den_str.strip('()')
    else:
        num_str = expression
        den_str = '1'
    
    def parse_polynomial(poly_str):
        """Parse a polynomial string into coefficient list."""
        if poly_str == '1': return [1]
        if poly_str == '0': return [0]
        
        # Initialize coefficient dictionary
        coeffs = {}
        max_power = 0
        
        # Replace minus signs for easier splitting
        poly_str = poly_str.replace('-', '+-')
        if poly_str.startswith('+'): poly_str = poly_str[1:]
        
        # Split into terms
        terms = poly_str.split('+')
        
        for term in terms:
            if not term: continue
            
            # Handle negative terms
            sign = 1
            if term.startswith('-'):
                sign = -1
                term = term[1:]
            
            # Parse the term
            if term == 's':
                coeffs[1] = sign * 1
                max_power = max(max_power, 1)
            elif 's' not in term:
                # Constant term
                coeffs[0] = sign * float(term)
            else:
                # Term with s
                if '^' in term:
                    # Has power
                    parts = term.split('s^')
                    power = int(parts[1])
                    coeff = parts[0]
                else:
                    # Linear term
                    parts = term.split('s')
                    power = 1
                    coeff = parts[0]
                
                if coeff == '':
                    coeffs[power] = sign * 1
                elif coeff == '*':
                    coeffs[power] = sign * 1
                else:
                    coeff = coeff.strip('*')
                    coeffs[power] = sign * float(coeff)
                
                max_power = max(max_power, power)
        
        # Convert to coefficient list
        coeff_list = [0] * (max_power + 1)
        for power, coeff in coeffs.items():
            coeff_list[max_power - power] = coeff
        
        return coeff_list
    
    # Get coefficients and create transfer function
    num = parse_polynomial(num_str)
    den = parse_polynomial(den_str)
    sys = control.TransferFunction(num, den)
    
    # Create figure with adjusted size and get rid of excess white space
    plt.figure(figsize=(12, 8))  # Increased figure size
    
    # Plot root locus with adjusted parameters
    control.root_locus(sys, grid=True, plot=True)
    
    # Get the current axis
    ax = plt.gca()
    
    # Make the plot more compact
    plt.tight_layout(pad=2.0)  # Reduce padding around the plot
    
    # Set title and labels with adjusted font sizes
    if title is None:
        title = f"Root Locus of G(s) = {expression}"
    plt.title(title, pad=2, fontsize=12)  # Adjusted title padding and font size
    plt.xlabel('Real Axis', fontsize=10)
    plt.ylabel('Imaginary Axis', fontsize=10)
    
    # Set axis limits if provided, otherwise auto-adjust
    if xlim:
        plt.xlim(xlim)
    if ylim:
        plt.ylim(ylim)
    
    # Add grid with better visibility
    plt.grid(True, linestyle='--', alpha=0.7)
    
    # Make plot lines thicker for better visibility
    for line in ax.get_lines():
        line.set_linewidth(1.5)
    
    # Adjust layout to fill space
    plt.tight_layout()
    
    plt.show()
    

