#Don's Standard Properties Formulation

import math
import numpy as np
import matplotlib.pyplot as plt

'''
=====================================
Standard Proprties Functions 
=====================================
'''

"""
Calculate the moment of inertia about the x-axis and y-axis passing through the centroid/axis for various 2-d shapes.

Parameters:
shape (str): The name of the geometric shape.
bottom_base (float): The length of the bottom base for shapes like triangles, rectangles, trapezoids.
top_base (float): The length of the top base for trapezoids.
height (float): The height of the shape.
radius (float): The radius for circular shapes.
chord_length (float): The length of the chord for circular segments.
inner_radius (float): The inner radius for hollow circles.

Returns:
(float): The moment of inertia about the x-axis and y-axis passing through the centroid/axis.
"""

def moment_of_inertia_x_center(shape, bottom_base=None, top_base=None, height=None, radius=None, chord_length=None, inner_radius=None):
    
    if shape == "right_triangle":
        if bottom_base is None or height is None:
            raise ValueError("Bottom base and height are required for a right triangle.")
        return (bottom_base * height**3) / 36

    elif shape == "rectangle":
        if bottom_base is None or height is None:
            raise ValueError("Bottom base and height are required for a rectangle.")
        return (bottom_base * height**3) / 12

    elif shape == "trapezoid":
        if bottom_base is None or top_base is None or height is None:
            raise ValueError("Bottom base, top base, and height are required for a trapezoid.")
        return (height**3 * (bottom_base**2 + 4*bottom_base*top_base + top_base**2)) / (36 * (bottom_base + top_base))

    elif shape == "circle":
        if radius is None:
            raise ValueError("Radius is required for a circle.")
        return (math.pi * radius**4) / 4

    elif shape == "hollow_circle":
        if radius is None or inner_radius is None:
            raise ValueError("Outer radius and inner radius are required for a hollow circle.")
        return (math.pi / 4) * (radius**4 - inner_radius**4)

    elif shape == "semicircle":
        if radius is None:
            raise ValueError("Radius is required for a semicircle.")
        return ((((9*((math.pi)**2)))-64) * (radius**4)) / (72 * math.pi)

    elif shape == "parabola":
        if bottom_base is None or height is None:
            raise ValueError("Base and height are required for a parabola.")
        return (8/15) * bottom_base * (height**3)

    elif shape == "half_parabola":
        if bottom_base is None or height is None:
            raise ValueError("Base and height are required for a half parabola.")
        return (8/175) * bottom_base * (height**3)

    else:
        raise ValueError("Unsupported shape. Check the function documentation for supported shapes.")

def moment_of_inertia_y_center(shape, bottom_base=None, top_base=None, height=None, radius=None, chord_length=None, inner_radius=None):
    
    if shape == "right_triangle":
        if bottom_base is None or height is None:
            raise ValueError("Bottom base and height are required for a right triangle.")
        return (bottom_base**3 * height) / 48

    elif shape == "rectangle":
        if bottom_base is None or height is None:
            raise ValueError("Bottom base and height are required for a rectangle.")
        return (height * bottom_base**3) / 12

    elif shape == "trapezoid":
        if bottom_base is None or top_base is None or height is None:
            raise ValueError("Bottom base, top base, and height are required for a trapezoid.")
        return (height * ((2*top_base)+bottom_base)) / (3*(top_base+bottom_base))

    elif shape == "circle":
        if radius is None:
            raise ValueError("Radius is required for a circle.")
        return (math.pi * radius**4) / 4

    elif shape == "hollow_circle":
        if radius is None or inner_radius is None:
            raise ValueError("Outer radius and inner radius are required for a hollow circle.")
        return (math.pi / 4) * (radius**4 - inner_radius**4)

    elif shape == "semicircle":
        if radius is None:
            raise ValueError("Radius is required for a semicircle.")
        return (1/8) * math.pi * (radius**4)

    elif shape == "parabola":
        if bottom_base is None or height is None:
            raise ValueError("Base and height are required for a parabola.")
        return (16/175) * (bottom_base**3) * height

    elif shape == "half_parabola":
        if bottom_base is None or height is None:
            raise ValueError("Base and height are required for a half parabola.")
        return (19/480) * (bottom_base**3) * height

    else:
        raise ValueError("Unsupported shape. Check the function documentation for supported shapes.")
    
def moment_of_inertia_x_axis(shape, bottom_base=None, top_base=None, height=None, radius=None, chord_length=None, inner_radius=None):

    if shape == "right_triangle":
        if bottom_base is None or height is None:
            raise ValueError("Bottom base and height are required for a right triangle.")
        return (bottom_base * height**3) / 12

    elif shape == "rectangle":
        if bottom_base is None or height is None:
            raise ValueError("Bottom base and height are required for a rectangle.")
        return (bottom_base * height**3) / 3

    elif shape == "trapezoid":
        if bottom_base is None or top_base is None or height is None:
            raise ValueError("Bottom base, top base, and height are required for a trapezoid.")
        return (height**3 * (bottom_base + (3* top_base))) / 12

    elif shape == "circle":
        if radius is None:
            raise ValueError("Radius is required for a circle.")
        return (5*(math.pi * radius**4)) / 4

    elif shape == "hollow_circle":
        if radius is None or inner_radius is None:
            raise ValueError("Outer radius and inner radius are required for a hollow circle.")
        return ((5*math.pi*(radius**4)) / 4) - (math.pi*((radius**2) * (inner_radius**2))) - ((math.pi*(inner_radius**4)) / 4)

    elif shape == "semicircle":
        if radius is None:
            raise ValueError("Radius is required for a semicircle.")
        return (math.pi * radius**4) / 8 

    elif shape == "parabola":
        if bottom_base is None or height is None:
            raise ValueError("Base and height are required for a parabola.")
        return (4/15) * bottom_base * height**3

    elif shape == "half_parabola":
        if bottom_base is None or height is None:
            raise ValueError("Base and height are required for a half parabola.")
        return (2/15) * bottom_base * height**3

    else:
        raise ValueError("Unsupported shape. Check the function documentation for supported shapes.")

def moment_of_inertia_y_axis(shape, bottom_base=None, top_base=None, height=None, radius=None, chord_length=None, inner_radius=None, offset= None):
    
    if shape == "right_triangle":
        if bottom_base is None or height is None:
            raise ValueError("Bottom base and height are required for a right triangle.")
        return (bottom_base**3 * height) / 4

    elif shape == "rectangle":
        if bottom_base is None or height is None:
            raise ValueError("Bottom base and height are required for a rectangle.")
        return (height * bottom_base**3) / 3

    elif shape == "trapezoid":
        if bottom_base is None or top_base is None or height is None:
            raise ValueError("Bottom base, top base, and height are required for a trapezoid.")
        return (height / 12) * (top_base**3 + 3*top_base*offset**2 + 3*top_base**2*offset +
                               bottom_base**3 + bottom_base*offset**2 + top_base*bottom_base**2 +
                               top_base*bottom_base*offset + bottom_base*top_base**2)
    elif shape == "circle":
        if radius is None:
            raise ValueError("Radius is required for a circle.")
        return (5*(math.pi * radius**4)) / 4

    elif shape == "hollow_circle":
        if radius is None or inner_radius is None:
            raise ValueError("Outer radius and inner radius are required for a hollow circle.")
        return ((5*math.pi*(radius**4)) / 4) - (math.pi*((radius**2) * (inner_radius**2))) - ((math.pi*(inner_radius**4)) / 4)
    
    elif shape == "semicircle":
        if radius is None:
            raise ValueError("Radius is required for a semicircle.")
        return (5*(math.pi * radius**4)) / 8

    elif shape == "parabola":
        if bottom_base is None or height is None:
            raise ValueError("Base and height are required for a parabola.")
        return (4/7) * bottom_base**3 * height

    elif shape == "half_parabola":
        if bottom_base is None or height is None:
            raise ValueError("Base and height are required for a half parabola.")
        return (2/7) * bottom_base**3 * height

    else:
        raise ValueError("Unsupported shape. Check the function documentation for supported shapes.")
    
    

"""
Calculate the area for various shapes.

Parameters:
shape (str): The name of the geometric shape.
bottom_base (float): The length of the bottom base for shapes like triangles, rectangles, trapezoids.
top_base (float): The length of the top base for trapezoids.
height (float): The height of the shape.
radius (float): The radius for circular shapes.
chord_length (float): The length of the chord for circular segments.
inner_radius (float): The inner radius for hollow circles.

Returns:
(float): The area of the shape.
"""

def calculate_area(shape, bottom_base=None, top_base=None, height=None, radius=None, chord_length=None, inner_radius=None):
    
    if shape == "right_triangle":
        if bottom_base is None or height is None:
            raise ValueError("Bottom base and height are required for a right triangle.")
        return 0.5 * bottom_base * height

    elif shape == "rectangle":
        if bottom_base is None or height is None:
            raise ValueError("Bottom base and height are required for a rectangle.")
        return bottom_base * height

    elif shape == "trapezoid":
        if bottom_base is None or top_base is None or height is None:
            raise ValueError("Bottom base, top base, and height are required for a trapezoid.")
        return 0.5 * (bottom_base + top_base) * height

    elif shape == "circle":
        if radius is None:
            raise ValueError("Radius is required for a circle.")
        return math.pi * radius**2

    elif shape == "hollow_circle":
        if radius is None or inner_radius is None:
            raise ValueError("Outer radius and inner radius are required for a hollow circle.")
        return math.pi * (radius**2 - inner_radius**2)

    elif shape == "semicircle":
        if radius is None:
            raise ValueError("Radius is required for a semicircle.")
        return 0.5 * math.pi * radius**2

    elif shape == "parabola":
        if bottom_base is None or height is None:
            raise ValueError("Base and height are required for a parabola.")
        return (2/3) * bottom_base * height

    elif shape == "half_parabola":
        if bottom_base is None or height is None:
            raise ValueError("Base and height are required for a half parabola.")
        return (1/3) * bottom_base * height

    else:
        raise ValueError("Unsupported shape. Check the function documentation for supported shapes.")
    


"""
Calculate the mass moment of inertia about the x, y and z axis for various 3-d shapes. Also calculates the moment of interia about the centroid.
WARNING: Orientation of axis may be identical to centroid axis orientation. Please check orientation of commonplace textbooks

Parameters:
mass (float): Mass of the object in kg
shape (str): Shape type ('rod', 'hollow_disk', 'thin_disk', etc.)
radius_outer (float, optional): Outer radius in meters
radius_inner (float, optional): Inner radius in meters for hollow shapes
length (float, optional): Length in meters
width (float, optional): Width in meters
height (float, optional): Height in meters

Returns:
(float): Moment of inertia about x-axis
"""
    
def calculate_moment_inertia_x_axis(mass, shape, radius_outer=None, radius_inner= None, length=None, width= None, height=None):
    
    if shape == 'rod':
        I_xx = 0
    elif shape == 'hollow_disk':
        if radius_outer is None and radius_inner is None:
            raise ValueError("Radius required for hollow disk")
        I_xx = 1.5 * mass * (((radius_inner+radius_outer)/2)**2)
    elif shape == 'thin_disk':
        if radius_outer is None:
            raise ValueError("Radius required for thin disk")
        I_xx = 0.25 * mass * (radius_outer**2)
    elif shape == 'cylinder':
        if radius_outer is None:
            raise ValueError("Radius required for cylinder")
        I_xx = (1/12) * mass * ((3* (radius_outer**2)) + (4*(height*2)))
    elif shape == 'hollow_cylinder':
        if radius_outer is None and radius_inner is None:
            raise ValueError("Radius required for hollow cylinder")
        I_xx = (1/12) * mass * ((3* (radius_outer**2)) + (3* (radius_inner**2)) + (4*(height*2)))
    elif shape == 'sphere':
        if radius_outer is None:
            raise ValueError("Radius required for sphere")
        I_xx = 0.4 * mass * (radius_outer**2)
    elif shape == 'hemisphere':
        if radius_outer is None:
            raise ValueError("Radius required for hemisphere")
        I_xx = 0.4 * mass * (radius_outer**2)
    elif shape == 'rectangular_plate':
        if width is None:
            raise ValueError("Height required for rectangular plate")
        I_xx = (1/12) * mass * (width**2)
    elif shape == 'cone':
        if radius_outer is None or height is None:
            raise ValueError("Both radius and height required for cone")
        I_xx = 0.05 * mass * ((3 * (radius_outer**2)) + (2*(height**2)))
    else:
        raise ValueError("Unsupported shape specified")
    
    # Print results
    print(f"Mass Moment of Inertia about x-axis - {shape.capitalize()}:")
    print(f"Ixx: {I_xx:.3e}")
    
    return I_xx
    
def calculate_moment_inertia_y_axis(mass, shape, radius_outer=None, radius_inner= None, length=None, width= None, height=None):
    
    if shape == 'rod':
        I_yy = (mass * (length**3))/12
    elif shape == 'hollow_disk':
        if radius_outer is None and radius_inner is None:
            raise ValueError("Radius required for hollow disk")
        I_yy = 1.5 * mass * (((radius_inner+radius_outer)/2)**2)
    elif shape == 'thin_disk':
        if radius_outer is None:
            raise ValueError("Radius required for thin disk")
        I_yy = 0.25 * mass * (radius_outer**2)
    elif shape == 'cylinder':
        if radius_outer is None:
            raise ValueError("Radius required for cylinder")
        I_yy = 0.5 * mass * ((radius_outer**2))
    elif shape == 'hollow_cylinder':
        if radius_outer is None and radius_inner is None:
            raise ValueError("Radius required for hollow cylinder")
        I_yy = 0.5 * mass * ((radius_outer**2) + (radius_inner**2))
    elif shape == 'sphere':
        if radius_outer is None:
            raise ValueError("Radius required for sphere")
        I_yy = 0.4 * mass * (radius_outer**2)
    elif shape == 'hemisphere':
        if radius_outer is None:
            raise ValueError("Radius required for hemisphere")
        I_yy = 0.4 * mass * (radius_outer**2)
    elif shape == 'rectangular_plate':
        if length is None:
            raise ValueError("Length required for rectangular plate")
        I_yy = (1/12) * mass * (length**2)
    elif shape == 'cone':
        if radius_outer is None or height is None:
            raise ValueError("Both radius and height required for cone")
        I_yy = 0.05 * mass * ((3 * (radius_outer**2)) + (2*(height**2)))
    else:
        raise ValueError("Unsupported shape specified")
    
    # Print results
    print(f"Mass Moment of Inertia about x-axis - {shape.capitalize()}:")
    print(f"Ixx: {I_yy:.3e}")
    
    return I_yy
    
def calculate_moment_inertia_z_axis(mass, shape, radius_outer=None, radius_inner= None, length=None, width= None, height=None):
    
    if shape == 'rod':
        I_zz = (mass * (length**3))/3
    elif shape == 'hollow_disk':
        if radius_outer is None and radius_inner is None:
            raise ValueError("Radius required for hollow disk")
        I_zz = 3 * mass * (((radius_inner+radius_outer)/2)**2)
    elif shape == 'thin_disk':
        if radius_outer is None:
            raise ValueError("Radius required for thin disk")
        I_zz = 0.5 * mass * (radius_outer**2)
    elif shape == 'cylinder':
        if radius_outer is None:
            raise ValueError("Radius required for cylinder")
        I_zz = (1/12) * mass * ((3* (radius_outer**2)) + (4*(height*2)))
    elif shape == 'hollow_cylinder':
        if radius_outer is None and radius_inner is None:
            raise ValueError("Radius required for hollow cylinder")
        I_zz = (1/12) * mass * ((3* (radius_outer**2)) + (3* (radius_inner**2)) + (4*(height*2)))
    elif shape == 'sphere':
        if radius_outer is None:
            raise ValueError("Radius required for sphere")
        I_zz = 0.4 * mass * (radius_outer**2)
    elif shape == 'hemisphere':
        if radius_outer is None:
            raise ValueError("Radius required for hemisphere")
        I_zz = 0.4 * mass * (radius_outer**2)
    elif shape == 'rectangular_plate':
        if width is None:
            raise ValueError("Height required for rectangular plate")
        I_zz = (1/12) * mass * ((width**2) + (length**2))
    elif shape == 'cone':
        if radius_outer is None or height is None:
            raise ValueError("Both radius and height required for cone")
        I_zz = 0.3 * (mass * (length**3))
    else:
        raise ValueError("Unsupported shape specified")
    
    # Print results
    print(f"Mass Moment of Inertia about z-axis - {shape.capitalize()}:")
    print(f"Ixx: {I_zz:.3e}")
    
    return I_zz
    
def calculate_moment_inertia_x_centroid(mass, shape, radius_outer=None, radius_inner= None, length=None, width= None, height=None):
    
    if shape == 'rod':
        I_xx_c = 0
    elif shape == 'hollow_disk':
        if radius_outer is None and radius_inner is None:
            raise ValueError("Radius required for hollow disk")
        I_xx_c = mass * (((radius_inner+radius_outer)/2)**2)
    elif shape == 'thin_disk':
        if radius_outer is None:
            raise ValueError("Radius required for thin disk")
        I_xx_c = 0.25 * mass * (radius_outer**2)
    elif shape == 'cylinder':
        if radius_outer is None:
            raise ValueError("Radius required for cylinder")
        I_xx_c = (1/12) * mass * ((3* (radius_outer**2)) + ((height*2)))
    elif shape == 'hollow_cylinder':
        if radius_outer is None and radius_inner is None:
            raise ValueError("Radius required for hollow cylinder")
        I_xx_c = (1/12) * mass * ((3* (radius_outer**2)) + (3* (radius_inner**2)) + ((height*2)))
    elif shape == 'sphere':
        if radius_outer is None:
            raise ValueError("Radius required for sphere")
        I_xx_c = 0.4 * mass * (radius_outer**2)
    elif shape == 'hemisphere':
        if radius_outer is None:
            raise ValueError("Radius required for hemisphere")
        I_xx_c = (83/320) * mass * (radius_outer**2)
    elif shape == 'rectangular_plate':
        if width is None:
            raise ValueError("Height required for rectangular plate")
        I_xx_c = (1/12) * mass * (width**2)
    elif shape == 'cone':
        if radius_outer is None or height is None:
            raise ValueError("Both radius and height required for cone")
        I_xx_c = (3/80) * mass * ((4 * (radius_outer**2)) + ((height**2)))
    else:
        raise ValueError("Unsupported shape specified")
    
    # Print results
    print(f"Mass Moment of Inertia about x-centroid - {shape.capitalize()}:")
    print(f"Ixx: {I_xx_c:.3e}")
    
    return I_xx_c
    
def calculate_moment_inertia_y_centroid(mass, shape, radius_outer=None, radius_inner= None, length=None, width= None, height=None):
    
    if shape == 'rod':
        I_yy_c = (mass * (length**3))/12
    elif shape == 'hollow_disk':
        if radius_outer is None and radius_inner is None:
            raise ValueError("Radius required for hollow disk")
        I_yy_c = 0.5 * mass * (((radius_inner+radius_outer)/2)**2)
    elif shape == 'thin_disk':
        if radius_outer is None:
            raise ValueError("Radius required for thin disk")
        I_yy_c = 0.25 * mass * (radius_outer**2)
    elif shape == 'cylinder':
        if radius_outer is None:
            raise ValueError("Radius required for cylinder")
        I_yy_c = 0.5 * mass * ((radius_outer**2))
    elif shape == 'hollow_cylinder':
        if radius_outer is None and radius_inner is None:
            raise ValueError("Radius required for hollow cylinder")
        I_yy_c = 0.5 * mass * ((radius_outer**2) + (radius_inner**2))
    elif shape == 'sphere':
        if radius_outer is None:
            raise ValueError("Radius required for sphere")
        I_yy_c = 0.4 * mass * (radius_outer**2)
    elif shape == 'hemisphere':
        if radius_outer is None:
            raise ValueError("Radius required for hemisphere")
        I_yy_c = (83/320) * mass * (radius_outer**2)
    elif shape == 'rectangular_plate':
        if length is None:
            raise ValueError("Length required for rectangular plate")
        I_yy_c = (1/12) * mass * (length**2)
    elif shape == 'cone':
        if radius_outer is None or height is None:
            raise ValueError("Both radius and height required for cone")
        I_yy_c = (3/80) * mass * ((3 * (radius_outer**2)) + ((height**2)))
    else:
        raise ValueError("Unsupported shape specified")
    
    # Print results
    print(f"Mass Moment of Inertia about y-centroid - {shape.capitalize()}:")
    print(f"I_yy_cc: {I_yy_c:.3e}")
    
    return I_yy_c
    
def calculate_moment_inertia_z_centroid(mass, shape, radius_outer=None, radius_inner= None, length=None, width= None, height=None):
    
    if shape == 'rod':
        I_zz_c = (mass * (length**3))/12
    elif shape == 'hollow_disk':
        if radius_outer is None and radius_inner is None:
            raise ValueError("Radius required for hollow disk")
        I_zz_c = mass * (((radius_inner+radius_outer)/2)**2)
    elif shape == 'thin_disk':
        if radius_outer is None:
            raise ValueError("Radius required for thin disk")
        I_zz_c = 1.5 * mass * (radius_outer**2)
    elif shape == 'cylinder':
        if radius_outer is None:
            raise ValueError("Radius required for cylinder")
        I_zz_c = (1/12) * mass * ((3* (radius_outer**2)) + ((height*2)))
    elif shape == 'hollow_cylinder':
        if radius_outer is None and radius_inner is None:
            raise ValueError("Radius required for hollow cylinder")
        I_zz_c = (1/12) * mass * ((3* (radius_outer**2)) + (3* (radius_inner**2)) + ((height*2)))
    elif shape == 'sphere':
        if radius_outer is None:
            raise ValueError("Radius required for sphere")
        I_zz_c = 0.4 * mass * (radius_outer**2)
    elif shape == 'hemisphere':
        if radius_outer is None:
            raise ValueError("Radius required for hemisphere")
        I_zz_c = 0.4 * mass * (radius_outer**2)
    elif shape == 'rectangular_plate':
        if width is None:
            raise ValueError("Height required for rectangular plate")
        I_zz_c = (1/12) * mass * ((width**2) + (length**2))
    elif shape == 'cone':
        if radius_outer is None or height is None:
            raise ValueError("Both radius and height required for cone")
        I_zz_c = 0.3 * (mass * (length**3))
    else:
        raise ValueError("Unsupported shape specified")
    
    # Print results
    print(f"Mass Moment of Inertia about z-centroid - {shape.capitalize()}:")
    print(f"Ixx: {I_zz_c:.3e}")
    
    return I_zz_c



"""
Calculate the volume of various geometric shapes.

Parameters:
shape (str): Shape type ('cylinder', 'cone', 'hemisphere', 'sphere', 'hollow_cylinder', etc.)
radius (float, optional): Outer radius in meters
height (float, optional): Height in meters
length (float, optional): Length in meters
width (float, optional): Width in meters
radius_inner (float, optional): Inner radius for hollow shapes in meters

Returns:
(float): Volume in cubic meters

Raises:
ValueError: If required parameters are missing or shape is not supported
"""

def calculate_volume(shape, radius=None, height=None, length=None, width=None, radius_inner=None):
    
    if shape == 'cylinder':
        if radius is None or height is None:
            raise ValueError("Both radius and height required for cylinder")
        volume = math.pi * (radius ** 2) * height
        
    elif shape == 'hollow_cylinder':
        if radius is None or height is None or radius_inner is None:
            raise ValueError("Outer radius, inner radius, and height required for hollow cylinder")
        if radius_inner >= radius:
            raise ValueError("Inner radius must be smaller than outer radius")
        volume = math.pi * height * (radius**2 - radius_inner**2)
        
    elif shape == 'cone':
        if radius is None or height is None:
            raise ValueError("Both radius and height required for cone")
        volume = (1/3) * math.pi * (radius**2) * height
        
    elif shape == 'hemisphere':
        if radius is None:
            raise ValueError("Radius required for hemisphere")
        volume = (2/3) * math.pi * (radius**3)
        
    elif shape == 'sphere':
        if radius is None:
            raise ValueError("Radius required for sphere")
        volume = (4/3) * math.pi * (radius**3)
        
    elif shape == 'rectangular_prism':
        if length is None or width is None or height is None:
            raise ValueError("Length, width, and height required for rectangular prism")
        volume = length * width * height
        
    else:
        raise ValueError(f"Unsupported shape: {shape}")
    
    # Print results for verification
    print(f"Volume of {shape}: {volume:.3e} m³")
    return volume

"""
Calculate the radius of gyration of a cross-section.The radius of gyration (r) is a geometric property of a cross-section that describes how the area is distributed relative to its centroidal axis. It is used in the analysis of columns and beams, particularly for stability calculations.

Parameters
I (float) = area moment of inertia (second moment of area)
A (float) = cross-sectional area

Returns:
r (float) = radius of gyration
"""

def calculate_radius_of_gyration (I, A):
    r = (I/A)**(1/2)
    return r



"""
Calculate various dimensionless numbers used in fluid mechanics, heat transfer, and material properties.

Parameters:
density (float): Fluid density (kg/m³)
velocity (float): Flow velocity (m/s)
char_length (float): Characteristic length (m)
dyn_viscosity (float): Dynamic viscosity (Pa·s)
gravity (float): Acceleration due to gravity (m/s²), defaults to 9.81
surface_tension (float): Surface tension (N/m)
pressure_diff (float): Pressure difference (Pa)
specific_heat_capacity (float): Specific heat capacity (J/kg·K)
thermal_conductivity (float): Thermal conductivity (W/m·K)
thermal_diffustivity (float): Theremal diffustivity (J/m^3 K)
heat_transfer_coeff (float): Heat transfer coefficient (W/m²·K)
beta (float): Thermal expansion coefficient (1/K)
temp_diff (float): Temperature difference (K)
kin_viscosity (float): Kinematic viscosity (m²/s)
grashof (float): Grashof number
prandtl (float): Prandtl number
mass_diffusivity (float): Mass diffusivity (m²/s)
mass_transfer_coeff (float): Mass transfer coefficient (m/s)
transverse_strain (float): Transverse strain
axial_strain (float): Axial strain
time (float): Time (s)
length (float): length (s)

number_type (str): Type of dimensionless number to calculate:
'reynolds' - (ρVL)/μ: Ratio of inertial forces to viscous forces
'froude' - V/√(gL): Ratio of inertial forces to gravitational forces
'weber' - (ρV²L)/σ: Ratio of inertial forces to surface tension forces
'euler' - ΔP/(ρV²): Ratio of pressure forces to inertial forces
'prandtl' - (Cpμ)/k: Ratio of momentum diffusivity to thermal diffusivity
'nusselt' - (hL)/k: Ratio of convective to conductive heat transfer
'grashof' - (gβ(Ts-T∞)L³)/ν²: Ratio of buoyancy to viscous forces
'rayleigh' - Gr × Pr: Product of Grashof and Prandtl numbers
'schmidt' - ν/D: Ratio of momentum diffusivity to mass diffusivity
'sherwood' - (kL)/D: Ratio of convective to diffusive mass transfer
'poisson' - -εtransverse/εaxial: Ratio of transverse to axial strain
'biot' - (hL)/k: Ratio of heat transfer resistances
'fourier' - (at)/L²: Diffusion time scale for heat transfer

Returns:
(float): Calculated dimensionless number
"""
    
def calculate_dimensionless_number(number_type, density = None, velocity = None, char_length = None, 
    dyn_viscosity = None, surface_tension = None, pressure_difference = None, specific_heat_capacity = None,
    thermal_conductivity = None, heat_transfer_coeff = None, beta = None, temp_diff = None, 
    kin_viscosity = None, mass_diffusivity = None, mass_transfer_coeff = None, transverse_strain = None,
    axial_strain = None, gravity = None, thermal_diffustivity = None, time = None, length = None):
    if number_type == 'reynolds':
        Re = (density * velocity * char_length) / dyn_viscosity
        return Re
    
    elif number_type == 'froude':
        Fr = velocity / ((gravity * char_length)**0.5)
        return Fr
        
    elif number_type == 'weber':
        We = (density * (velocity**2) * char_length) / surface_tension
        return We
    
    elif number_type == 'euler':
        Eu = pressure_difference / (density * (velocity**2))
        return Eu
        
    elif number_type == 'prandtl':
        Pr = (specific_heat_capacity * dyn_viscosity) / thermal_conductivity
        return Pr
    
    elif number_type == 'nusselt':
        Nu = (heat_transfer_coeff *char_length) / thermal_conductivity
        return Nu
    
    elif number_type == 'grashof':
        Gr = (gravity * beta * temp_diff * char_length**3) / (kin_viscosity**2)
        return Gr
    
    elif number_type == 'rayleigh':
        Ra = Gr * Pr
        return Ra
    
    elif number_type == 'schmidt':
        Sc = kin_viscosity / mass_diffusivity
        return Sc
    
    elif number_type == 'sherwood':
        Sh = (mass_transfer_coeff * char_length) / mass_diffusivity
        return Sh
    
    elif number_type == 'poisson':
        v = transverse_strain / axial_strain
        return v
    
    elif number_type == 'biot':
        Bi = (heat_transfer_coeff * char_length) / thermal_conductivity
        return Bi
    
    elif number_type == 'fourier':
        Fo = (thermal_diffustivity * time) / (length**2)
        return Fo
    
    else:
        raise ValueError(f"Unknown dimensionless number type: {number_type}")








