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
float: Volume in cubic meters

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
Calculate the radius of gyration of a cross-section. The radius of gyration (r) is a geometric property of a cross-section that describes how the area is distributed relative to its centroidal axis. It is used in the analysis of columns and beams, particularly for stability calculations.

Parameters:
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

Returns:
(float): Calculated dimensionless number
"""
    
def calculate_dimensionless_number(number_type, density = None, velocity = None, char_length = None, 
    dyn_viscosity = None, surface_tension = None, pressure_difference = None, specific_heat_capacity = None,
    thermal_conductivity = None, heat_transfer_coeff = None, beta = None, temp_diff = None, 
    kin_viscosity = None, mass_diffusivity = None, mass_transfer_coeff = None, transverse_strain = None,
    axial_strain = None, gravity = None):
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
    
    else:
        raise ValueError(f"Unknown dimensionless number type: {number_type}")










'''
=====================================
Standard Mechancial Functions 
=====================================
'''

"""
Calculate and print the principal stresses and their components.

Parameters:
sigma_x (float): Normal stress in the x-direction.
sigma_y (float): Normal stress in the y-direction.
tau_xy (float): Shear stress in the xy-plane.

Returns:
(tuple): A tuple containing the maximum principal stress, minimum principal stress, and maximum shear stress.
"""

def calculate_principal_stresses(sigma_x, sigma_y, tau_xy):
    # Calculate average normal stress
    sigma_avg = (sigma_x + sigma_y) / 2

    # Calculate radius of Mohr's circle
    R = math.sqrt(((sigma_x - sigma_y) / 2)**2 + tau_xy**2)

    # Calculate principal stresses
    sigma_1 = sigma_avg + R
    sigma_2 = sigma_avg - R

    # Calculate maximum shear stress
    tau_max = R

    # Calculate angle of principal plane
    theta_p = 0.5 * math.atan2(2 * tau_xy, sigma_x - sigma_y)

    # Print results
    print(f"Principal Stresses and Components:")
    print(f"Maximum Principal Stress (σ1): {sigma_1:.2f}")
    print(f"Minimum Principal Stress (σ2): {sigma_2:.2f}")
    print(f"Maximum Shear Stress (τmax): {tau_max:.2f}")
    print(f"Average Normal Stress (σavg): {sigma_avg:.2f}")
    print(f"Angle of Principal Plane (θp): {math.degrees(theta_p):.2f}°")
    print(f"Normal Stress in x-direction (σx): {sigma_x:.2f}")
    print(f"Normal Stress in y-direction (σy): {sigma_y:.2f}")
    print(f"Shear Stress in xy-plane (τxy): {tau_xy:.2f}")

    return sigma_1, sigma_2, tau_max



"""
Calculate the bulk modulus of elasticity using Young's modulus and Poisson's ratio.

Parameters:
E (float): Young's modulus (Pa)
nu (float): Poisson's ratio

Returns:
(float): Bulk modulus of elasticity (Pa)
"""
    
def calculate_bulk_modulus(E, nu):
    # Check if Poisson's ratio is valid (must be less than 0.5 for physical materials)
    if nu >= 0.5:
        raise ValueError("Poisson's ratio must be less than 0.5 for stable materials")
    
    # Calculate bulk modulus (K = E/(3(1-2ν)))
    K = E / (3 * (1 - (2*nu)))
    
    # Print results
    print(f"Bulk Modulus (K): {K:.2e} Pa")
    
    return K



"""
Calculate 3D strain components from stress using Hooke's law.

Parameters:
sigma_x (float): Normal stress in x-direction (Pa)
sigma_y (float): Normal stress in y-direction (Pa)
sigma_z (float): Normal stress in z-direction (Pa)
E (float): Young's modulus (Pa)
nu (float): Poisson's ratio

Returns:
(tuple): Strain components (epsilon_x, epsilon_y, epsilon_z)
"""
    
def calculate_3d_strain_from_stress(sigma_x, sigma_y, sigma_z, E, nu):
    # Calculate strain components
    epsilon_x = (1/E) * (sigma_x - nu*(sigma_y + sigma_z))
    epsilon_y = (1/E) * (sigma_y - nu*(sigma_x + sigma_z))
    epsilon_z = (1/E) * (sigma_z - nu*(sigma_x + sigma_y))
    
    # Print results
    print(f"3D Strain Components from Stress:")
    print(f"Strain in x-direction (εx): {epsilon_x:.2e}")
    print(f"Strain in y-direction (εy): {epsilon_y:.2e}")
    print(f"Strain in z-direction (εz): {epsilon_z:.2e}")
    
    return epsilon_x, epsilon_y, epsilon_z



"""
Calculate 3D stress components from strain using Hooke's law.

Parameters:
epsilon_x (float): Strain in x-direction
epsilon_y (float): Strain in y-direction
epsilon_z (float): Strain in z-direction
E (float): Young's modulus (Pa)
nu (float): Poisson's ratio

Returns:
(tuple): Stress components (sigma_x, sigma_y, sigma_z)
"""
    
def calculate_3d_stress_from_strain(epsilon_x, epsilon_y, epsilon_z, E, nu):
    # Calculate factor for cleaner equations
    factor = E/((1 + nu)*(1 - (2*nu)))
    
    # Calculate stress components
    sigma_x = factor * (((1-nu)*epsilon_x) + (nu*(epsilon_y + epsilon_z)))
    sigma_y = factor * (((1-nu)*epsilon_y) + (nu*(epsilon_x + epsilon_z)))
    sigma_z = factor * (((1-nu)*epsilon_z) + (nu*(epsilon_x + epsilon_y)))
    
    # Print results
    print(f"3D Stress Components from Strain:")
    print(f"Stress in x-direction (σx): {sigma_x:.2e} Pa")
    print(f"Stress in y-direction (σy): {sigma_y:.2e} Pa")
    print(f"Stress in z-direction (σz): {sigma_z:.2e} Pa")
    
    return sigma_x, sigma_y, sigma_z



"""
Calculate strain energy density from stress components using equation 7-57b.

Parameters:
sigma_x (float): Normal stress in x-direction (Pa)
sigma_y (float): Normal stress in y-direction (Pa)
sigma_z (float): Normal stress in z-direction (Pa)
E (float): Young's modulus (Pa)
nu (float): Poisson's ratio

Returns:
(float): Strain energy density (J/m³)
"""
    
def calculate_strain_energy_density_from_stress(sigma_x, sigma_y, sigma_z, E, nu):
    # Calculate first term (σx² + σy² + σz²)
    stress_squared_sum = (sigma_x**2 + sigma_y**2 + sigma_z**2)
    
    # Calculate second term (σxσy + σyσz + σzσx)
    stress_products = (sigma_x*sigma_y + sigma_y*sigma_z + sigma_z*sigma_x)
    
    # Calculate strain energy density
    u = ((1/(2*E)) * stress_squared_sum) - ((nu/E) * stress_products)
    
    # Print results
    print(f"Strain Energy Density from Stress Components:")
    print(f"Strain Energy Density (u): {u:.2e} J/m³")
    
    return u



"""
Calculate strain energy density from strain components using equation 7-57c.

Parameters:
epsilon_x (float): Strain in x-direction
epsilon_y (float): Strain in y-direction
epsilon_z (float): Strain in z-direction
E (float): Young's modulus (Pa)
nu (float): Poisson's ratio

Returns:
float: Strain energy density (J/m³)
"""
    
def calculate_strain_energy_density_from_strain(epsilon_x, epsilon_y, epsilon_z, E, nu):
    # Calculate factor outside brackets
    factor = E / (2 * (1 + nu) * (1 - 2*nu))
    
    # Calculate first term ((1-ν)(εx² + εy² + εz²))
    strain_squared_term = (1 - nu) * (epsilon_x**2 + epsilon_y**2 + epsilon_z**2)
    
    # Calculate second term (2ν(εxεy + εyεz + εzεx))
    strain_products_term = 2 * nu * (epsilon_x*epsilon_y + epsilon_y*epsilon_z + epsilon_z*epsilon_x)
    
    # Calculate strain energy density
    u = factor * (strain_squared_term + strain_products_term)
    
    # Print results
    print(f"Strain Energy Density from Strain Components:")
    print(f"Strain Energy Density (u): {u:.2e} J/m³")

    return u



"""
Calculate the Modified-Mohr Effective Stress for 3D stress state.

Parameters:
sigma1 (float): First principal stress
sigma2 (float): Second principal stress
sigma3 (float): Third principal stress
S_ut (float): Ultimate tensile strength
S_yc (float): Yield strength in compression

Returns:
(float): Modified-Mohr effective stress
"""
    
def calculate_modified_mohr_stress(sigma1, sigma2, sigma3, S_ut, S_yc):
    # Calculate the term (2S_ut - |S_yc|)/|S_yc| which is used in all C calculations
    stress_ratio_term = (2 * S_ut - abs(S_yc)) / abs(S_yc)
    
    # Calculate C1, C2, and C3 according to the formulas
    C1 = 0.5 * (abs(sigma1 - sigma2) + (stress_ratio_term * (sigma1 + sigma2)))
    C2 = 0.5 * (abs(sigma2 - sigma3) + (stress_ratio_term * (sigma2 + sigma3)))
    C3 = 0.5 * (abs(sigma3 - sigma1) + (stress_ratio_term * (sigma3 + sigma1)))
    
    # Find the maximum value among C1, C2, C3, sigma1, sigma2, sigma3
    max_value = max(C1, C2, C3, sigma1, sigma2, sigma3)
    
    # Apply the condition: if MAX < 0, then σ̄ = 0
    effective_stress = max(0, max_value)
    
    return effective_stress








'''
=====================================
Standard Fatigue Functions
=====================================
'''

"""
Calculate the uncorrected fatigue strength estimate based on material and ultimate strength.

Parameters:
material (str): one of 'steel', 'iron', 'aluminum', or 'copper'
ultimate_strength (float): ultimate strength in MPa

Return: 
(float): uncorrected fatigue strength estimate in MPa
"""

def calculate_uncorrected_fatigue_strength(material, ultimate_strength ,units):
    if units == "MPa":
        
        if material == 'steel':
            if ultimate_strength < 1400:
                return 0.5 * ultimate_strength
            else:
                return 700
            
        elif material == 'iron':
            if ultimate_strength < 400:
                return 0.4 * ultimate_strength
            else:
                return 160
            
        elif material == 'aluminum':
            if ultimate_strength < 330:
                return 0.4 * ultimate_strength
            else:
                return 130
            
        elif material == 'copper':
            if ultimate_strength < 280:
                return 0.4 * ultimate_strength
            else:
                return 100
        
        else:
            raise ValueError("Invalid material. Choose 'steel', 'iron', 'aluminum', or 'copper'.")
        
    if units == "kpsi":
        
        if material == 'steel':
            if ultimate_strength < 200:
                return 0.5 * ultimate_strength
            else:
                return 100
            
        elif material == 'iron':
            if ultimate_strength < 60:
                return 0.4 * ultimate_strength
            else:
                return 24
            
        elif material == 'aluminum':
            if ultimate_strength < 48:
                return 0.4 * ultimate_strength
            else:
                return 19
            
        elif material == 'copper':
            if ultimate_strength < 40:
                return 0.4 * ultimate_strength
            else:
                return 14
        
        else:
            raise ValueError("Invalid material. Choose 'steel', 'iron', 'aluminum', or 'copper'.")




"""
Calculate the temperature correction factor (C_temp).

Parameters:
temperature (float): temperature in Celsius or Fahrenheit
units (str): 'C' for Celsius or 'F' for Fahrenheit

Return:
(float): temperature correction factor
"""
    
def temperature_correction_factor(temperature, units='C'):
    if units == 'F':
        if temperature <= 840:
            return 1
        elif 840 < temperature <= 1020:
            return 1 - 0.0032 * (temperature - 840)
    
    elif units == 'C':
        if temperature <= 450:
            return 1
        elif 450 < temperature <= 550:
            return 1 - 0.0058 * (temperature - 450)
    
    else:
        raise ValueError("Invalid units. Use 'C' for Celsius or 'F' for Fahrenheit.")
    
    raise ValueError(f"Temperature {temperature}°{units} is out of the supported range.")



"""
Calculate the load correction factor (C_load).

Parameters:
load_type (str): 'bending' or 'axial'

Return: 
(float): load correction factor
"""
    
def load_correction_factor(load_type):
    if load_type == 'bending':
        return 1
    
    elif load_type == 'axial':
        return 0.70
    
    else:
        raise ValueError("Invalid load type. Use 'bending' or 'axial'.")



"""
Calculate the surface finish correction factor (C_surf).

Parameters:
ultimate_strength (float): ultimate strength in MPa or psi
surface_finish (str): 'ground', 'machined', 'hot-rolled', or 'as-forged'
units (str): 'MPa' or 'psi'

Return:
(float): surface finish correction factor
"""

def surface_finish_correction_factor(ultimate_strength, surface_finish, units='MPa'):
    coefficients = {
        'ground': (1.58, -0.085, 2.411, -0.085),
        'machined': (4.51, -0.265, 16.841, -0.265),
        'hot-rolled': (57.7, -0.718, 2052.9, -0.718),
        'as-forged': (272, -0.995, 38545.0, -0.995)
    }
    
    if surface_finish not in coefficients:
        raise ValueError("Invalid surface finish. Choose 'ground', 'machined', 'hot-rolled', or 'as-forged'.")
    
    if units == 'MPa':
        A, b = coefficients[surface_finish][:2]
    elif units == 'psi':
        A, b = coefficients[surface_finish][2:]
    else:
        raise ValueError("Invalid units. Use 'MPa' or 'psi'.")
    
    C_surf = A * (ultimate_strength ** b)
    return min(C_surf, 1.0)



"""
Calculate the equivalent diameter for various cross-sections based on the 95% stressed area.

Parameters:
shape (str), type of cross-section ('circle', 'rectangle', 'channel', or 'I-beam')
diameter (float): diameter for circular cross-sections (mm)
base (float): width/base of the cross-section (mm)
height (float): height of the cross-section (mm)
thickness (float): thickness for channel and I-beam cross-sections (mm)
distance_from_top (float): distance from the top of the channel section to the point of interest (mm)
rotation (str): specifies if the circular section is rotating ('yes' or 'no')

axis (str): specifies the axis of interest for channel and I-beam sections
For channel: 'vertical_center' or 'horizontal'
For I-beam: 'horizontal_center' or 'vertical_center'

Return: 
(float): equivalent diameter based on the 95% stressed area (mm)

Notes:
- For circular sections, specify 'diameter' and 'rotation'.
- For rectangular sections, specify 'base' and 'height'.
- For channel sections, specify 'base', 'height', 'thickness', 'distance_from_top', and 'axis'.
- For I-beam sections, specify 'base', 'height', 'thickness', and 'axis'.
"""
    
def calculate_cross_section_diameter(shape=None, diameter=None, base=None, height=None, thickness=None, distance_from_top=None, rotation=None, axis=None):
    def calculate_equivalent_diameter(A95):
        return math.sqrt(A95 / 0.0766)

    if shape == "circle" and diameter is not None:
        if rotation == "yes":
            A95 = 0.0766 * diameter**2
        elif rotation == "no":
            A95 = 0.010462 * diameter**2
        else:
            raise ValueError("For circular cross-sections, rotation must be 'yes' or 'no'")
        return calculate_equivalent_diameter(A95)

    elif shape == "rectangle" and base is not None and height is not None:
        A95 = 0.05 * base * height
        return calculate_equivalent_diameter(A95)

    elif shape == "channel" and base is not None and height is not None and thickness is not None and distance_from_top is not None:
        if axis == "vertical_center":
            A95 = 0.05 * base * height
        elif axis == "horizontal_noncenter":
            A95 = 0.05 * base * thickness + 0.10 * thickness * (height - thickness)
            if distance_from_top > thickness:
                A95 += 0.05 * base * distance_from_top + 0.10 * distance_from_top * (height - distance_from_top)
        else:
            raise ValueError("For channel cross-sections, axis must be 'vertical_center' or 'horizontal'")
        return calculate_equivalent_diameter(A95)

    elif shape == "I-beam" and base is not None and height is not None and thickness is not None:
        if axis == "horizontal_center":
            A95 = 0.05 * base * height
        elif axis == "vertical_center":
            A95 = 0.1 * base * thickness
        else:
            raise ValueError("For I-beam cross-sections, axis must be 'horizontal_center' or 'vertical_center'")
        return calculate_equivalent_diameter(A95)

    else:
        raise ValueError("Invalid combination of parameters or shape not recognized.")



"""
Calculate the size correction factor (C_size).

Parameters:
diameter (float): diameter in mm or inches
units (str): 'mm' or 'in'

Return: 
(float): size correction factor
"""

def size_correction_factor(diameter, units='mm'):
    
    if units == 'mm':
        if diameter <= 8:
            return 1
        elif 8 < diameter <= 250:
            return 1.189 * (diameter ** -0.097)
    
    elif units == 'in':
        if diameter <= 0.3:
            return 1
        elif 0.3 < diameter <= 10:
            return 0.869 * (diameter ** -0.097)
    
    else:
        raise ValueError("Invalid units. Use 'mm' or 'in'.")
    
    raise ValueError(f"Diameter {diameter} {units} is out of the supported range.")



"""
Calculate the reliability correction factor (C_reliab) based on the desired reliability percentage.

Parameters:
reliability_percent (float): desired reliability percentage (50 <= reliability_percent <= 99.9999)

Return:
(float): reliability correction factor
"""
    
def reliability_correction_factor(reliability_percent):
    reliability_factors = {
        50: 1.000,
        90: 0.897,
        95: 0.868,
        99: 0.814,
        99.9: 0.753,
        99.99: 0.702,
        99.999: 0.659,
        99.9999: 0.620
    }
    
    if reliability_percent not in reliability_factors:
        raise ValueError("Invalid reliability percentage. Choose from: 50, 90, 95, 99, 99.9, 99.99, 99.999, or 99.9999")
    
    return reliability_factors[reliability_percent]



"""
Calculate the corrected fatigue strength considering all correction factors.

Parameters:
material (str): one of 'steel', 'iron', 'aluminum', or 'copper'
ultimate_strength (float): ultimate strength in MPa or psi
material_strength_units (str): 'MPa' or 'psi' (default is 'MPa')
value_temp (float): temperature value (default is 20)
temperature_units (str): 'C' or 'F' (default is 'C')
load_type (str): 'bending' or 'axial' (default is 'bending')
surface_finish (str): 'ground', 'machined', 'hot-rolled', or 'as-forged' (default is 'ground')
diameter (float): diameter of the component (default is 7)
diameter_units (str): 'mm' or 'in' (default is 'mm')
percent_reliability (float): desired reliability percentage (default is 50)

Return: 
(float): corrected fatigue strength in MPa
"""
    
def calculate_corrected_fatigue_strength_MPA(material,ultimate_strength,material_strength_units = "Mpa", temperture =20 , temperture_units = "C",
    load_type = "axial", surface_finish = "machined", surface_finish_units = "MPa", diameter_units = "mm", percent_reliability = 99, 
    shape = "circle", base = None, height = None, thickness = None, distance_from_top =None, axis = "horizontal" ):
    
   Se = calculate_uncorrected_fatigue_strength(material, ultimate_strength ,material_strength_units)
   Ct = temperature_correction_factor (temperture, temperture_units)
   Cl = load_correction_factor (load_type)
   Csr = surface_finish_correction_factor (ultimate_strength, surface_finish, surface_finish_units)
   Cr = reliability_correction_factor (percent_reliability)
   
   diameter = calculate_cross_section_diameter(shape, base, height, thickness, distance_from_top, axis="horizontal")
   Cs = size_correction_factor (diameter, diameter_units)
   
   Corrected_fatigure_strength_MPa = Se*Ct*Cl*Csr*Cs*Cr
   Corrected_fatigure_strength_kpsi = Corrected_fatigure_strength_MPa * 0.14503773800722
   return Corrected_fatigure_strength_MPa, Corrected_fatigure_strength_kpsi












'''
=====================================
Standard Beam and Shaft Functions
=====================================
'''

"""
Calculate various types of stresses in a beam based on loading conditions

Parameters:
load_type (str): Type of stress calculation to perform:
- "normal_axial": Normal stress due to axial loading
- "direct_shear": Average shear stress
- "bending_axial": Normal stress due to bending
- "bending_shear": Shear stress due to bending
- "torsion_shear": Shear stress due to torsion

P (float): Axial or shear force 
A (float): Cross-sectional area 
M (float): Bending moment 
y (float): Distance from neutral axis 
I (float): Second moment of area (moment of inertia) 
T (float): Applied torque 
r (float): Radial distance from center 
J (float): Polar moment of inertia 
V (float): Vertical shear force 
Q (float): First moment of area 
b (float): Width of section at point of interest
c_i (float): Inner radius curvature 
r_i (float): Inner radius 
Q_val: Section factor for torsion

Returns:
(float): Calculated stress 
"""
    
def standard_stress_calcualtion_on_beam(load_type, condition = None, P = None, A = None, M = None, 
    y = None, I = None, T = None, r = None, J = None, V = None, Q = None, b = None, e = None,
    c_i = None, r_i = None, Q_val = None):
    
    if load_type == "normal_axial":
        sigma = P/A
    
    elif load_type == "direct_shear":
        sigma = P/A
    
    elif load_type == "bending_axial":
        sigma = (M * y)/I
    
    elif load_type == "bending_shear":
        sigma = (V * Q) / (I * b)
        
    elif load_type == "torsional_shear":
        sigma = (T * r) / (J)
    
    elif load_type == "transverse_shear" and condition == "curved_beam":
        sigma = (M * c_i) / (e * A * r_i)
        
    elif load_type == "transverse_shear" and condition == "rectangular_beam":
        sigma = 1.5 * (V / A)
        
    elif load_type == "transverse_shear" and condition == "circular_beam":
        sigma = (4/3) * (V / A)
        
    elif load_type == "torsional_shear" and condition == "non_circular_beam":
        sigma = (T / Q_val)
        
    else:
        print("Error, load type not avaliable")
        
    return sigma



"""
Calculate torsional constant Q_val for various cross-sectional shapes.

Parameters:
shape (str): Cross-section shape type:
- "solid_square"
- "hollow_square"
- "solid_rectangle"
- "hollow_rectangle"
- "solid_ellipse"
- "hollow_ellipse"
- "circular_tube"
- "arbitrary"

length (float): Side length for squares, width for rectangles, or semi-major axis for ellipses
width (float): Width for rectangles or semi-minor axis for ellipses
thickness (float): Wall thickness for hollow sections
radius (float): Radius for circular tube
median_length (float): Length of median line for arbitrary shape

Returns:
(float): Torsional constant Q_val [length^3]
"""
    
def calculate_torsional_Q_val(shape, length=None, width=None, thickness=None, radius=None, median_length_line=None):
    if shape == "solid_square":
        if length is None:
            raise ValueError("Square side length required")
        return (1/0.6) * ((length/2)**3)
        
    elif shape == "hollow_square":
        if length is None or thickness is None:
            raise ValueError("Square side length and thickness required")
        return (2 * (thickness) * ((length - thickness)**2))
        
    elif shape == "solid_rectangle":
        if length is None or width is None:
            raise ValueError("Rectangle length and width required")
        return (8 * ((length/2)**2) * ((width/2)**2)) / ((3 *(length/2)) + (1.8*(width/2)))
        
    elif shape == "hollow_rectangle":
        if length is None or width is None or thickness is None:
            raise ValueError("Rectangle length, width, and thickness required")
        return (2 * (thickness) * (((length/2) - thickness)) * (((width/2) - thickness))) 
        
    elif shape == "solid_ellipse":
        if length is None or width is None:  # length = a (semi-major), width = b (semi-minor)
            raise ValueError("Ellipse semi-major and semi-minor axes required")
        return (math.pi * ((length/2)) * ((width/2)**2)) / 2
        
    elif shape == "hollow_ellipse":
        if length is None or width is None or thickness is None:
            raise ValueError("Ellipse semi-major, semi-minor axes, and thickness required")
        return ((math.pi * ((length/2)) * ((width/2)**2))/2) * (1 - ((1 - (thickness/(length/2)))**4))
        
    elif shape == "circular_tube":
        if radius is None or thickness is None:
            raise ValueError("Tube radius and thickness required")
        if thickness >= radius:
            raise ValueError("Thickness must be less than radius")
        return ((4*(math.pi**2))* ((radius)**2) * ((thickness)**2)) / ((6 * math.pi * radius) + (1.8*thickness))
        
    elif shape == "open_arbitrary":
        if thickness is None or median_length_line is None:
            raise ValueError("Thickness and median line length required")
        return ((thickness**2) * (median_length_line**2)) / ((1.8 * thickness) + (3 * median_length_line))
        
    else:
        raise ValueError(f"Invalid shape: {shape}")



"""
Shaft Torsional Deflection (Section 10.9)

Parameters:
T (float): Torque
l (float): Length
G (float): Shear modulus
J (float): Polar moment of inertia
K (float): Torsional constant for non-round shafts

Return: 
theta (float): Angular deflection
"""
    
def shaft_torsional_deflection(type, T, l, G, J, K):
    if type == "round":
        deflection = (T * l) / (G * J)
    
    elif type == "non_round":
        deflection = (T * l) / (G * K)
        
    else:
        print ("Error, type not avaliable")
        
    return deflection



"""
Calculate torsional constant K for various cross-sectional shapes.

Parameters:
shape (str): Cross-section shape type:
- "solid_square"
- "hollow_square"
- "solid_rectangle"
- "hollow_rectangle"
- "solid_ellipse"
- "hollow_ellipse"
- "circular_tube"
- "arbitrary"

length (float): Side length for squares, width for rectangles, or semi-major axis for ellipses
width (float): Width for rectangles or semi-minor axis for ellipses
thickness (float): Wall thickness for hollow sections
radius (float): Radius for circular tube
median_length (float): Length of median line for arbitrary shape

Returns:
(float): Torsional constant K [length^4]
"""
    
def calculate_torsional_K(shape, length=None, width=None, thickness=None, radius=None, median_length_line=None):
    if shape == "solid_square":
        if length is None:
            raise ValueError("Square side length required")
        return 2.25 * ((length/2)**4)
        
    elif shape == "hollow_square":
        if length is None or thickness is None:
            raise ValueError("Square side length and thickness required")
        return (2 * (thickness**2) * (length - thickness)**4) / ((2 * length * thickness) - (2 * (thickness**2)))
        
    elif shape == "solid_rectangle":
        if length is None or width is None:
            raise ValueError("Rectangle length and width required")
        return (length/2) * ((width/2)**3) * ((16/3) - (3.36 * ((width/2)/(length/2)) * (1 - (((width/2)**4)/(12*((length/2)**4))))))
        
    elif shape == "hollow_rectangle":
        if length is None or width is None or thickness is None:
            raise ValueError("Rectangle length, width, and thickness required")
        return (2 * (thickness**2) * (((length/2) - thickness)**2) * (((width/2) - thickness)**2)) / \
               (((length/2)*thickness) + ((width/2)*thickness) - (2*(thickness**2)))
        
    elif shape == "solid_ellipse":
        if length is None or width is None:  # length = a (semi-major), width = b (semi-minor)
            raise ValueError("Ellipse semi-major and semi-minor axes required")
        return (math.pi * ((length/2)**3) * ((width/2)**3)) / (((length/2)**2) + ((width/2)**2))
        
    elif shape == "hollow_ellipse":
        if length is None or width is None or thickness is None:
            raise ValueError("Ellipse semi-major, semi-minor axes, and thickness required")
        return (math.pi * ((length/2)**3) * ((width/2)**3)) / (((length/2)**2) + ((width/2)**2)) * (1 - ((1 - (thickness/(length/2)))**4))
        
    elif shape == "circular_tube":
        if radius is None or thickness is None:
            raise ValueError("Tube radius and thickness required")
        if thickness >= radius:
            raise ValueError("Thickness must be less than radius")
        return (math.pi * (2/3)) * (thickness**3) * radius
        
    elif shape == "open_arbitrary":
        if thickness is None or median_length_line is None:
            raise ValueError("Thickness and median line length required")
        return (1/3) * (thickness**3) * median_length_line
        
    else:
        raise ValueError(f"Invalid shape: {shape}")



"""
General Shaft-Design Equation

Parameters:
N_f (float): Safety factor
K_f (float): Fatigue stress-concentration factor for bending
M_a (float): Alternating bending moment
K_fa (float): Fatigue stress-concentration factor for alternating torsion
T_a (float): Alternating torsional moment
K_fm (float): Fatigue stress-concentration factor for mean bending
M_m (float): Mean bending moment
K_fsm (float): Fatigue stress-concentration factor for mean torsion
T_m (float): Mean torsional moment
S_f (float): Fatigue strength
S_ut (float): Ultimate tensile strength

Return:
d (float): Shaft diameter
"""
    
def general_shaft_design(N_f, K_f, M_a, K_fa, T_a, K_fm, M_m, K_fsm, T_m, S_f, S_ut):
    term1 = ((((K_f * M_a)**2) + (0.75 * ((K_fa * T_a)**2)))**0.5) / S_f
    term2 = ((((K_fm * M_m)**2) + (0.75 * ((K_fsm * T_m)**2)))**0.5) / S_ut
    
    d = (((32 * N_f) / math.pi) * (term1 + term2))**(1/3)
    return d



"""
Pressure Generated by an Interference Fit 

Parameters:
delta (float): Radial interference
r_o (float): Outer radius of hub
r_i (float): Inner radius of shaft
r (float): Radius of fit
E_o (float): Young's modulus of hub material
v_o (float): Poisson's ratio of hub material
E_i (float): Young's modulus of shaft material
v_i (float): Poisson's ratio of shaft material

Return: 
p (float): Pressure
"""
    
def pressure_interference_fit(delta, r_o, r_i, r, E_o, v_o, E_i, v_i):
    term1 = (r_o**2 + r**2) / (E_o * (r_o**2 - r**2)) + (v_o / E_o)
    term2 = (r**2 + r_i**2) / (E_i * (r**2 - r_i**2)) - (v_i / E_i)
    p = (0.5 * delta) / (r * (term1 + term2))
    return p



"""
Tangential Stresses in Shaft and Hub of an Interference Fit (Section 10.11)

Parameters:
p (float): Pressure
r (float): Radius of fit
r_i (float): Inner radius of shaft

Return: 
sigma_t_shaft (float): Tangential stress in shaft
"""
    
def tangential_stress_interference_fit(p, r, r_i): 
    sigma_t_shaft = -p * (r**2 + r_i**2) / (r**2 - r_i**2)
    return sigma_t_shaft












'''
=====================================
Standard Material Functions 
=====================================
'''

"""
Calculate the fracture toughness (KIC) for a material.

Parameters:
stress (float): Applied engineering stress (σ) in MPa or appropriate units
crack_length (float): Crack length (a) in meters or appropriate units
geometry_factor (float, optional): Geometrical factor (Y), defaults to 1.0
crack_type (str, optional): Type of crack ("interior" or "exterior"), 
                            affects default geometry factor

Returns:
(float): Fracture toughness (KIC) in MPa⋅√m or appropriate units

Note: For exterior cracks, Y ≈ 1.1; for interior cracks, Y ≈ 1.0
"""
    
def calculate_fracture_toughness(stress, crack_length, geometry_factor=1.0, crack_type="interior"):
    if crack_length <= 0:
        raise ValueError("Crack length must be positive")
    if stress <= 0:
        raise ValueError("Stress must be positive")
    
    # Adjust geometry factor based on crack type if not explicitly provided
    if geometry_factor == 1.0:  # Only adjust if not explicitly set
        if crack_type.lower() == "exterior":
            geometry_factor = 1.1
        elif crack_type.lower() == "interior":
            geometry_factor = 1.0
        else:
            raise ValueError("Invalid crack type. Must be 'interior' or 'exterior'")
    
    # Calculate fracture toughness: KIC = Y⋅σ⋅√(π⋅a)
    fracture_toughness = geometry_factor * stress * math.sqrt(math.pi * crack_length)
    
    return fracture_toughness



"""
Calculate the steady-state creep rate (dε/dt).

Parameters:
stress (float): Applied stress (σ) in MPa or appropriate units
temperature (float): Absolute temperature (T) in Kelvin
activation_energy (float): Activation energy for creep (Q) in J/mol
pre_exp_constant (float): Pre-exponential constant (A)
stress_sensitivity (float): Stress sensitivity exponent (n)
gas_constant (float, optional): Ideal gas constant (R) in J/(mol⋅K), defaults to 8.314 J/(mol⋅K)

Returns:
(float): Creep rate (dε/dt) in appropriate units
"""
    
def calculate_creep_rate(stress, temperature, activation_energy, 
                        pre_exp_constant, stress_sensitivity, gas_constant=8.314):
    if temperature <= 0:
        raise ValueError("Temperature must be positive (in Kelvin)")
    if stress <= 0:
        raise ValueError("Stress must be positive")
    if activation_energy <= 0:
        raise ValueError("Activation energy must be positive")
    
    # Calculate creep rate: dε/dt = A⋅σ^n⋅e^(-Q/RT)
    creep_rate = (pre_exp_constant * 
                 (stress ** stress_sensitivity) * 
                 math.exp(-activation_energy / (gas_constant * temperature)))
    
    return creep_rate












'''
=====================================
Standard Column Functions
=====================================
'''

"""
Calculates the effective slenderness ratio of a column based on its length, radius of gyration, and end conditions.

Parameters:
L (float): Actual length of the column
r (float): Radius of gyration of the column cross-section
end_condition (str): Type of end supports. Options are:
- "rounded_rounded": Both ends rounded
- "pinned_pinned": Both ends pinned
- "fixed_free": One end fixed, other end free (cantilever)
- "fixed_pinned": One end fixed, other end pinned
- "fixed_fixed": Both ends fixed

Returns:
(float): Effective slenderness ratio (S_r)
"""
    
def calculate_effective_slenderness_ratio(L, r, end_condition = None):
    S_r = L/r
    
    if end_condition == "rounded_rounded":
        S_r == S_r
    elif end_condition == "pinned_pinned":
        S_r == S_r
    elif end_condition == "fixed_free":
        S_r == 2.1 * S_r
    elif end_condition == "fixed_pinned":
        S_r == 0.8 * S_r
    elif end_condition == "fixed_fixed":
        S_r == 0.65 * S_r
        
    return S_r



"""
Calculates the transition point between Euler and Johnson buckling formulas.

Parameters:
E (float): Modulus of elasticity of the material
Sr (float): Effective slenderness ratio

Returns:
tuple: (S_rD, S_yc):
S_rD (float): Slenderness ratio at transition point
S_yc (float): Compressive yield stress at transition point
"""
    
def calculate_euler_johnson_cutoff (E, Sr):
    S_yc = ((2*(math.pi**2)) * E) / (Sr**2)
    S_rD = (math.pi) * (((2*E) / S_yc) ** (1/2))
    return S_rD, S_yc



"""
Calculates the critical buckling load for a column using either Euler or Johnson formula based on the slenderness ratio.

Parameters:
L (float): Length of column
r (float): Radius of gyration
end_condition (str): Type of end supports (see calculate_effective_slenderness_ratio)
E (float): Modulus of elasticity
A (float): Cross-sectional area

Returns:
(float): Critical buckling load (P_cr)
"""

def calcualte_critical_column_buckling_load (L, r, end_condition, E, A):
    S_r = calculate_effective_slenderness_ratio (L, r, end_condition)
    S_rD, S_yc = calculate_euler_johnson_cutoff (E, S_r)
    
    if S_r < S_rD:
        P_cr = A * (S_yc - ((1/E)*(((S_yc * S_r)/(2 * math.pi))**2)))
        return P_cr
    
    if S_r >= S_rD:
        P_cr = ((math.pi**2) * E * A)/(S_r**2)
        return P_cr












'''
=====================================
Standard Pressurized Cylinder Functions
=====================================
'''

"""
Calculate stresses in a pressurized cylinder at a given radius.

Parameters:
r (float): Radius at which to calculate stresses
ri (float): Inner radius of cylinder
ro (float): Outer radius of cylinder
pi (float): Internal pressure
po (float): External pressure

Returns:
(dict): Dictionary containing radial, tangential, and axial stresses
"""
    
def calculate_cylinder_stresses(type, r, ri, ro, pi, po):
    # Calculate common terms to avoid repetition
    ri_squared = ri**2
    ro_squared = ro**2
    r_squared = r**2
    
    # Calculate denominator term used in multiple equations
    denom = ro_squared - ri_squared
    
    # Tangential (hoop) stress
    if type == "tangental":
        sigma_t = (pi * ri_squared - po * ro_squared) / denom + \
        (ri_squared * ro_squared * (pi - po)) / (r_squared * denom)
    
    # Radial stress
    if type == "radial":
        sigma_r = (pi * ri_squared - po * ro_squared) / denom - \
        (ri_squared * ro_squared * (pi - po)) / (r_squared * denom)
    
    # Axial stress - equation 4.47c
    if type == "axial":
        sigma_a = (pi * ri_squared - po * ro_squared) / denom
    
    return {
        'tangential_stress': sigma_t,
        'radial_stress': sigma_r,
        'axial_stress': sigma_a
    }
    


'''
=====================================
Standard Singularity function
=====================================
'''

"""
Analyze beam and create diagrams for various loading conditions.

Parameters:
L (float): Length of beam in meters
loads (dict): Dictionary containing load definitions:
'point_loads': list of tuples [(magnitude in kN, position in m), ...]
'distributed_loads': list of tuples [(intensity in kN/m, start_pos in m, end_pos in m), ...]
'moment_loads': list of tuples [(magnitude in kN·m, position in m), ...] (positive = CCW)
'supports': list of support positions [left_support, right_support] in m
EI (float): Flexural rigidity (E*I) in kN·m², default 1e6 kN·m²

Returns:
figure: matplotlib Figure object containing all diagrams
results (dict): Dictionary containing calculation results:
    'reactions': tuple (R1, R2) - support reactions in kN
    'shear': array - shear force diagram values in kN
    'moment': array - bending moment diagram values in kN·m
    'slope': array - slope diagram values in radians
    'deflection': array - deflection diagram values in meters
    'x': array - position values along beam in meters

Notes:
- Positive forces point upward
- Positive moments are counterclockwise (CCW)
- Generates 5 diagrams: loading, shear, moment, slope, and deflection
- All positions are measured from left support (x=0)

Example Useage:
if __name__ == "__main__":
beam_length = 10
loading_conditions = {
'point_loads': [(20, 7)],
'distributed_loads': [(15, 2, 5)],
'moment_loads': [(30, 4), (-20, 8)],
'supports': [0, 10]
}

# Analyze beam
EI = 2e8
fig, results = analyze_beam(beam_length, loading_conditions, EI)
plt.show(block=True)  # Force the window to stay open
"""

def analyze_simply_supported_beam(L, loads, EI=1e6, beam_type = None):
    # Create position array for calculations
    x = np.linspace(0, L, 1000)
    dx = x[1] - x[0]
    
    def calculate_reactions():
        """Calculate support reactions including moment loads"""
        # For a simply supported beam
        moment_sum = 0
        force_sum = 0
        
        # Add point loads contribution
        for P, a in loads.get('point_loads', []):
            moment_sum += P * a
            force_sum += P
            
        # Add distributed loads contribution
        for w, start, end in loads.get('distributed_loads', []):
            force = w * (end - start)
            center = (start + end) / 2
            moment_sum += force * center
            force_sum += force
            
        # Add moment loads contribution
        for M, a in loads.get('moment_loads', []):
            moment_sum += M  # Direct moment contribution
            
        R2 = moment_sum / L     # Right reaction
        R1 = force_sum - R2     # Left reaction
        
        return R1, R2
    
    def calculate_shear(x, R1, R2):
        """Calculate shear force at position x"""
        V = np.zeros_like(x)
        
        # Add support reactions
        V += R1 * (x >= 0)  # Left support
        V += R2 * (x >= L)  # Right support
        
        # Add point loads
        for P, a in loads.get('point_loads', []):
            V -= P * (x >= a)
            
        # Add distributed loads
        for w, start, end in loads.get('distributed_loads', []):
            V -= w * np.maximum(0, np.minimum(x - start, end - start))
        
        return V
        
    def calculate_moment(x, V):
        """Calculate bending moment including moment loads"""
        M = -np.cumsum(V) * dx
        
        # Create a copy of M to avoid modifying the original array
        M_total = M.copy()
        
        # Add moment loads (step function)
        for M_load, a in loads.get('moment_loads', []):
            M_total += M_load * (x >= a)
        
        # Adjust to make moment zero at supports
        M_total -= np.interp(x, [0, L], [M_total[0], M_total[-1]])
        
        return M_total
    
    def calculate_slope_deflection(x, M):
        """Calculate slope and deflection from moment"""
        slope = np.cumsum(M) * dx / EI
        slope -= slope[0]
        
        deflection = np.cumsum(slope) * dx
        support_deflection = np.interp(x, [0, L], [deflection[0], deflection[-1]])
        deflection -= support_deflection
        
        return slope, deflection
    
    def plot_moment_load(ax, pos, magnitude, size=1.5):
        """Plot moment load symbol"""
        circle = plt.Circle((pos, 0), size, fill=False)
        ax.add_artist(circle)
        
        # Add arrow to show direction (positive = CCW)
        if magnitude > 0:
            arrow_start = pos - size
            arrow_end = pos + size
        else:
            arrow_start = pos + size
            arrow_end = pos - size
            
        arrow_y = size
        ax.arrow(arrow_start, arrow_y, (arrow_end - arrow_start) * 0.3, 0,
                head_width=0.3, head_length=0.3, fc='r', ec='r')
        
        # Add magnitude label
        ax.text(pos, 2*size, f'{abs(magnitude)}kN·m', ha='center')
    
    # Calculate all beam responses
    R1, R2 = calculate_reactions()
    V = calculate_shear(x, R1, R2)
    M = calculate_moment(x, V)
    slope, deflection = calculate_slope_deflection(x, M)
    
    # Plot diagrams
    plt.ion()  # Turn on interactive mode
    fig, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(5, 1, figsize=(12, 15))
    
    # Loading diagram
    ax1.plot([0, L], [0, 0], 'k-', linewidth=2)
    ax1.set_title('Loading Diagram')
    # Plot support reactions
    ax1.arrow(0, 0, 0, R1/5, head_width=0.2, head_length=R1/25, fc='b', ec='b')
    ax1.arrow(L, 0, 0, R2/5, head_width=0.2, head_length=R2/25, fc='b', ec='b')
    ax1.text(0, R1/4, f'R1={R1:.1f}kN', ha='right')
    ax1.text(L, R2/4, f'R2={R2:.1f}kN', ha='left')
    
    # Add loads to diagram
    for P, a in loads.get('point_loads', []):
        ax1.arrow(a, 0, 0, -P/2, head_width=0.1, head_length=2, fc='r', ec='r')
    for w, start, end in loads.get('distributed_loads', []):
        ax1.plot([start, end], [-w, -w], 'r-')
        ax1.fill([start, start, end, end], [0, -w, -w, 0], alpha=0.1, color='r')
    for M_load, a in loads.get('moment_loads', []):
        plot_moment_load(ax1, a, M_load)
        
    # Set consistent y-limits for loading diagram
    ax1.set_ylim([-max([w for w, _, _ in loads.get('distributed_loads', [(0, 0, 0)])]) * 1.5,
                  max(R1, R2) * 0.4])
    
    # Shear diagram
    ax2.plot(x, V, 'b-')
    ax2.fill_between(x, 0, V, alpha=0.1)
    ax2.set_title('Shear Diagram')
    ax2.grid(True)
    
    # Moment diagram
    ax3.plot(x, M, 'g-')
    ax3.fill_between(x, 0, M, alpha=0.1)
    ax3.set_title('Moment Diagram')
    ax3.grid(True)
    
    # Slope diagram
    ax4.plot(x, slope, 'r-')
    ax4.fill_between(x, 0, slope, alpha=0.1)
    ax4.set_title('Slope Diagram')
    ax4.grid(True)
    
    # Deflection diagram
    ax5.plot(x, deflection, 'm-')
    ax5.fill_between(x, 0, deflection, alpha=0.1)
    ax5.set_title('Deflection Diagram')
    ax5.grid(True)
    
    plt.tight_layout()
    return fig, {
        'reactions': (R1, R2),
        'shear': V,
        'moment': M,
        'slope': slope,
        'deflection': deflection,
        'x': x
    }



"""
Analyze cantilever beam and create diagrams for various loading conditions. The fixed support is at x=0, and the free end is at x=L.

Parameters:
L (float): Length of beam in meters
loads (dict): Dictionary containing load definitions:
'point_loads': list of tuples [(magnitude in kN, position in m), ...]
'distributed_loads': list of tuples [(intensity in kN/m, start_pos in m, end_pos in m), ...]
'moment_loads': list of tuples [(magnitude in kN·m, position in m), ...] (positive = CCW)
EI (float): Flexural rigidity (E*I) in kN·m², default 1e6 kN·m²

Returns:
figure: matplotlib Figure object containing all diagrams
results (dict): Dictionary containing calculation results:
'reactions': tuple (R, M) - support reaction force (kN) and moment (kN·m)
'shear': array - shear force diagram values in kN
'moment': array - bending moment diagram values in kN·m
'slope': array - slope diagram values in radians
'deflection': array - deflection diagram values in meters
'x': array - position values along beam in meters

Example Useage:
if __name__ == "__main__":
# Example usage for a cantilever beam
beam_length = 6  # meters
loading_conditions = {
'point_loads': [(10, 4)],  # 10 kN at x=4m
'distributed_loads': [(5, 1, 3)],  # 5 kN/m from x=1m to x=3m
'moment_loads': [(15, 2)]  # 15 kN·m at x=2m
}

{'point_loads': [(10, 4)],'distributed_loads': [(5, 1, 3)], 'moment_loads': [(15, 2)]}

# Analyze beam
EI = 2e8  # kN·m²
fig, results = analyze_cantilever_beam(beam_length, loading_conditions, EI)
plt.show(block=True)
"""
    
def analyze_cantilever_beam(L, loads, EI=1e6):
    # Create position array for calculations
    x = np.linspace(0, L, 1000)
    dx = x[1] - x[0]
    
    def calculate_reactions():
        """Calculate fixed support reactions"""
        force_sum = 0
        moment_sum = 0
        
        # Add point loads contribution
        for P, a in loads.get('point_loads', []):
            force_sum += P
            moment_sum += P * a
            
        # Add distributed loads contribution
        for w, start, end in loads.get('distributed_loads', []):
            force = w * (end - start)
            center = (start + end) / 2
            force_sum += force
            moment_sum += force * center
            
        # Add moment loads contribution
        for M, a in loads.get('moment_loads', []):
            moment_sum += M
            
        R = force_sum  # Reaction force at support
        M = moment_sum  # Reaction moment at support
        
        return R, M
    
    def calculate_shear(x, R):
        """Calculate shear force at position x"""
        V = np.zeros_like(x)
        
        # Add support reaction
        V += R * (x >= 0)
        
        # Add point loads
        for P, a in loads.get('point_loads', []):
            V -= P * (x >= a)
            
        # Add distributed loads
        for w, start, end in loads.get('distributed_loads', []):
            V -= w * np.maximum(0, np.minimum(x - start, end - start))
        
        return V
        
    def calculate_moment(x, V, M_support):
        """Calculate bending moment"""
        # Initialize with support moment
        M = np.zeros_like(x)
        M += M_support
        
        # Add moment due to shear forces
        M -= np.cumsum(V) * dx
        
        # Add moment loads
        for M_load, a in loads.get('moment_loads', []):
            M += M_load * (x >= a)
        
        return M
    
    def calculate_slope_deflection(x, M):
        """Calculate slope and deflection from moment"""
        # For cantilever: slope and deflection are zero at fixed support (x=0)
        slope = np.cumsum(M) * dx / EI
        deflection = np.cumsum(slope) * dx
        
        # No need to adjust for support conditions since cantilever is fixed at x=0
        return slope, deflection
    
    def plot_moment_load(ax, pos, magnitude, size=1.5):
        """Plot moment load symbol"""
        circle = plt.Circle((pos, 0), size, fill=False)
        ax.add_artist(circle)
        
        # Add arrow to show direction (positive = CCW)
        if magnitude > 0:
            arrow_start = pos - size
            arrow_end = pos + size
        else:
            arrow_start = pos + size
            arrow_end = pos - size
            
        arrow_y = size
        ax.arrow(arrow_start, arrow_y, (arrow_end - arrow_start) * 0.3, 0,
                head_width=0.3, head_length=0.3, fc='r', ec='r')
        
        ax.text(pos, 2*size, f'{abs(magnitude)}kN·m', ha='center')
    
    # Calculate all beam responses
    R, M_support = calculate_reactions()
    V = calculate_shear(x, R)
    M = calculate_moment(x, V, M_support)
    slope, deflection = calculate_slope_deflection(x, M)
    
    # Plot diagrams
    plt.ion()
    fig, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(5, 1, figsize=(12, 15))
    
    # Loading diagram
    ax1.plot([0, L], [0, 0], 'k-', linewidth=2)
    ax1.set_title('Loading Diagram')
    
    # Plot fixed support
    support_height = 4
    ax1.plot([0, 0], [-support_height, support_height], 'k-', linewidth=3)
    ax1.fill([-0.2, 0, 0, -0.2], [-support_height, -support_height, support_height, support_height],
            'gray', alpha=0.3)
    
    # Plot support reactions
    ax1.arrow(0, 0, 0, R/5, head_width=0.2, head_length=R/25, fc='b', ec='b')
    ax1.text(0.5, R/4, f'R={R:.1f}kN', ha='left')
    
    # Plot support moment
    if abs(M_support) > 0:
        plot_moment_load(ax1, 0, M_support)
        ax1.text(0.5, -3, f'M={M_support:.1f}kN·m', ha='left')
    
    # Add loads to diagram
    for P, a in loads.get('point_loads', []):
        ax1.arrow(a, 0, 0, -P/2, head_width=0.1, head_length=2, fc='r', ec='r')
    for w, start, end in loads.get('distributed_loads', []):
        ax1.plot([start, end], [-w, -w], 'r-')
        ax1.fill([start, start, end, end], [0, -w, -w, 0], alpha=0.1, color='r')
    for M_load, a in loads.get('moment_loads', []):
        plot_moment_load(ax1, a, M_load)
    
    # Set consistent y-limits for loading diagram
    max_distributed = max([w for w, _, _ in loads.get('distributed_loads', [(0, 0, 0)])], default=0)
    ax1.set_ylim([-max_distributed * 1.5, max(abs(R) * 0.4, support_height)])
    
    # Shear diagram
    ax2.plot(x, V, 'b-')
    ax2.fill_between(x, 0, V, alpha=0.1)
    ax2.set_title('Shear Diagram')
    ax2.grid(True)
    
    # Moment diagram
    ax3.plot(x, M, 'g-')
    ax3.fill_between(x, 0, M, alpha=0.1)
    ax3.set_title('Moment Diagram')
    ax3.grid(True)
    
    # Slope diagram
    ax4.plot(x, slope, 'r-')
    ax4.fill_between(x, 0, slope, alpha=0.1)
    ax4.set_title('Slope Diagram')
    ax4.grid(True)
    
    # Deflection diagram
    ax5.plot(x, deflection, 'm-')
    ax5.fill_between(x, 0, deflection, alpha=0.1)
    ax5.set_title('Deflection Diagram')
    ax5.grid(True)
    
    plt.tight_layout()
    return fig, {
        'reactions': (R, M_support),
        'shear': V,
        'moment': M,
        'slope': slope,
        'deflection': deflection,
        'x': x
    }

