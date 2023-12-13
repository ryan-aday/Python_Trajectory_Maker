'''
runTraj.py
Written by Ryan Aday
Copyright @2023
V4: 12/12/2023- Added ENU/ECEF solver, commented out the flat LLA solver

The current ENU/ECEF algorithm first converts the emplacement location from geodetic/LLA into ENU coordinates. 
The XYZ trajectory data is then added to the ENU coordinates, after which the sums are converted
into geodetic/LLA coordinates. These are converted to geocentric coordinates (WGS84 currently). 

The Hughes transform matrix performs a coordinate transformation by rotating and translating the emplacement
coordinate system into that of the globe model. 
This approach is incorrect, as it assumes incorrectly that the earth's surface is flat as it plots trajectories.

The LLA flat earth algorithm first converts the emplacement location from geodetic/LLA into flat-earth XYZ coordinates. 
The XYZ trajectory data is then added to the flat-earth XYZ coordinates, after which the sums are converted
into geodetic/LLA coordinates. These are converted to geocentric coordinates (WGS84 currently). 
This approach is incorrect, as it assumes incorrectly that the x, y, and z changes are linear. They are not, 
which is what the ENU solution resolves. 
'''

# Import libraries
import math, numpy

'''
#############################################################################################
User Constant Inumpyuts
#############################################################################################
'''

#Defines whether this is a global or relative reference to make the .trj file
isGlobal = True # boolean, True/False

# Define filename
filename = 'traj_check.txt'

# WGS-84 Constants (needed for a global reference)
wgs84_lat = 32.416 # degrees
wgs84_lon = -106.279 # degrees
wgs84_alt = 1206.0 # meters

# Target Constants
distance_from_emplacement = 80000.0 # meters
search_sector_angle = 120.0 # degrees
bearing_N = 0.0 # 6133.3 # mills (Note: 1 degree = 160/9 milla)
speed = 210.0 # m/s
altitude = 80000.0 # meters
start_time = 0.0 # seconds
end_time = 5000.0 # seconds

# Extra Constants
pitch = 0.0 # degrees, negative is down
yaw = 0.0 # RCS degrees

'''
#############################################################################################
Conversion Functions ***DO NOT MODIFY
#############################################################################################
'''

### Hughes transform functions ###
def Rz(z_angle):
    return numpy.array([\
    [math.cos(z_angle), math.sin(z_angle), 0],\
    [-math.sin(z_angle), math.cos(z_angle), 0],\
    [0, 0, 1],])
'''
def Ry(y_angle):
    return numpy.array([\
    [math.cos(y_angle), 0, -math.sin(y_angle)],\
    [0, 1, 0],\
    [math.sin(y_angle), 0, math.cos(y_angle)],])
'''    
def Rx(x_angle):
    return numpy.array([\
    [1, 0, 0],\
    [0, math.cos(x_angle), math.sin(x_angle)],\
    [0, -math.sin(x_angle), math.cos(x_angle)],])
'''   
def Rt(x_angle, y_angle, z_angle):
        return Rz(z_angle) * Ry(y_angle) * Rx(x_angle)
        
def H(Rt, x_ECEF, y_ECEF, z_ECEF):
    Hf = numpy.zeros((4,4))
    
    # Rotation matrix entries
    Hf[0][0] = Rt[0][0]
    Hf[0][1] = Rt[0][1]    
    Hf[0][2] = Rt[0][2]
    Hf[1][0] = Rt[1][0]
    Hf[1][1] = Rt[1][1]    
    Hf[1][2] = Rt[1][2]   
    Hf[2][0] = Rt[2][0]
    Hf[2][1] = Rt[2][1]    
    Hf[2][2] = Rt[2][2]
    
    # Translation matrix entries
    Hf[0][3] = x_ECEF
    Hf[1][3] = y_ECEF
    Hf[2][3] = z_ECEF
    Hf[3][3] = 1
    
    return Hf

def resM(Ht, x, y, z):
    ref = numpy.array([x, y, z, 1])
    # X, Y, Z should be the absolute ECEF coordinates
    return numpy.dot(Ht, ref)

def angle(v1, v2):
    # Assumes v1 and v2 are both 1x3 vectors, outputs in rads
    return numpy.arccos(numpy.dot(v1, v2)/(numpy.linalg.norm(v1) * numpy.linalg.norm(v2)))
'''
### Flat LLA solution functions ###
'''
def geodetic_to_geocentric(ellps, lat, lon, h):
    a, rf = ellps
    
    b = a * (1 - 1/rf)
    e = math.sqrt(1- (b**2/a**2))
    
    lat_rad = math.radians(lat)
    lon_rad = math.radians(lon)
    
    N = a / math.sqrt(1 - e**2*math.sin(lat_rad)**2)
    
    X = ((N + h) * math.cos(lat_rad) * math.cos(lon_rad))
    Y = ((N + h) * math.cos(lat_rad) * math.sin(lon_rad))    
    Z = ((1 - e) ** 2 * N + h) * math.sin(lat_rad)

    return X, Y, Z
    
def ECEF_flat_to_geodetic(ellps, x, y, z):
    R, rf = ellps
    r = math.sqrt(x**2 + z**2)
    
    lat = math.degrees(math.pi/2 - r/R)
    lon = math.degrees(math.acos(z/r) * sign(x))
    h = y
    
    return lat, lon, h
    
def geodetic_to_ECEF_flat(ellps, lat, lon, h):
    R, rf = ellps
    
    lat_rad = math.radians(lat)
    lon_rad = math.radians(lon)
    
    A = (math.pi / 2 - lat_rad) * R
    x = A * math.sin(lon_rad)
    y = h
    z = A * math.cos(lon_rad)
    return x, y, z
    
    
def sign(x):
    return 1.0 if x>= 0 else -1.0
'''

### ENU solution functions ###
def geodetic_to_ENU(ellps, lat, lon, x, y, z):
    lat_rad = math.radians(lat)
    lon_rad = math.radians(lon)

    R1 = Rx(math.pi/2 - lat_rad)
    R3 = Rz(math.pi/2 + lon_rad)
    xyz = numpy.array([x, y, z])
    ENU = numpy.dot(numpy.dot(R1, R3), xyz)

    easting = ENU[0]
    northing = ENU[1]
    up = ENU[2]

    return easting, northing, up

def ENU_to_geodetic(ellps, lat, lon, easting, northing, up):
    lat_rad = math.radians(lat)    
    lon_rad = math.radians(lon)

    R1 = Rx(-(math.pi/2 - lat_rad))
    R3 = Rz(-(math.pi/2 + lon_rad))
    ENU = numpy.array([x, y, z])
    ENU = numpy.dot(numpy.dot(R1, R3), xyz)
    xyz = numpy.dot(numpy.dot(R3, R1), ENU)

    x = xyz[0]
    y = xyz[1]
    z = xyz[2]

    return x, y, z

'''
#############################################################################################
Functions
#############################################################################################
'''    

# Nodel Params
#grs80 = (6378137, 298.257222100882711)
wgs84 = (6378137, 298.257223563)

# ECEF Reference Vectors
'''
x_ECEF_ref = [1, 0, 0]
y_ECEF_ref = [0, 1, 0]
z_ECEF_ref = [0, 0, 1]
'''

# ECEF Conversion of emplacement location
x_emplacement_ECEF, y_emplacement_ECEF, z_emplacement_ECEF = geodetic_to_geocentric(wgs84, wgs84_lat, wgs84_lon, wgs84_alt)
emplacement_ECEF = [x_emplacement_ECEF, y_emplacement_ECEF, z_emplacement_ECEF]

# Derive H Transform Matrix
'''
x_angle = angle(emplacement_ECEF, x_ECEF_ref)
y_angle = angle(emplacement_ECEF, y_ECEF_ref)
z_angle = angle(emplacement_ECEF, z_ECEF_ref)

R = Rt(x_angle, y_angle, z_angle)
Ht = H(R, x_emplacement_ECEF, y_emplacement_ECEF, z_emplacement_ECEF)
'''

# Set up flat XYZ references for the flat globe solution
# x_flat_ref, y_flat_ref, z_flat_ref = geodetic_to_ECEF_flat(wgs84, wgs84_lat, wgs84_lon, wgs84_alt)

# Set up ENU references for the ENU/ECEF solution
easting_ref, northing_ref, up_ref = geodetic_to_ENU(wgs84, wgs84_lat, wgs84_lon, wgs84_alt, x_emplacement_ECEF, y_emplacement_ECEF, z_emplacement_ECEF)

# Set the time interval to be 0.1 seconds
time_interval = 0.1 # seconds

# Calculate the start angle offset(assumed to be from the aft)
bearing_N_deg = bearing_N * 9 / 160
start_angle_offset = -math.radians(search_sector_angle / 2 + bearing_N_deg)

# Calculate the total time range
total_time = end_time - start_time

# Calculate angular velocity (radians per second)
angular_velocity = math.radians(search_sector_angle) / total_time

# Initialize arrays to store
t_intervals = []
x_coordinates = []
y_coordinates = []
z_coordinates = []
x_ECEF_coordinates = []
y_ECEF_coordinates = []
z_ECEF_coordinates = []
pitch_intervals = []
yaw_intervals = []

# Calculate coordinates for each time interval
for t in numpy.arange(start_time, end_time + time_interval, time_interval):
    '''
    #############################################################################################
    Start of user-defined X, Y, Z trajectory characteristics
    #############################################################################################
    '''
    
    # Calculate X, Y, Z coordinates
    x = distance_from_emplacement * math.sin(start_angle_offset + angular_velocity * t) # km
    y = altitude # km
    z = distance_from_emplacement * math.cos(start_angle_offset + angular_velocity * t) # km
          
    '''
    #############################################################################################
    End of user-defined X, Y, Z trajectory characteristics
    #############################################################################################
    '''   

    # Hughes transform solution
    '''
    #sol = resM(Ht, x_ECEF, y_ECEF, z_ECEF)
    sol = resM(Ht, x, y, z)

    x_ECEF = sol[0]
    y_ECEF = sol[1]
    z_ECEF = sol[2]
    '''
    
    # Convert to ECEF using the flat LLA solution
    '''
    x_flat = x + x_flat_ref
    y_flat = y + y_flat_ref    
    z_flat = z + z_flat_ref

    lat, lon, h = ECEF_flat_to_geodetic(wgs84, x_flat, y_flat, z_flat)
    x_ECEF, y_ECEF, z_ECEF = geodetic_to_geocentric(wgs84, lat, lon, h)
    '''

    # Convert to ECEF using the ENU solution
    x_ECEF, y_ECEF, z_ECEF = ENU_to_geocentric(wgs84, wgs84_lat, wgs84_lon, easting_ref + x, northing_ref + z, up_ref + y)
    
    # Store data in arrays
    t_intervals.append(t) # seconds
    x_coordinates.append(x/1000) # km
    y_coordinates.append(y/1000) # km
    z_coordinates.append(z/1000) # km
    x_ECEF_coordinates.append(x_ECEF) # meters
    y_ECEF_coordinates.append(y_ECEF) # meters
    z_ECEF_coordinates.append(z_ECEF) # meters
    pitch_intervals.append(pitch) # degrees
    yaw_intervals.append(yaw) # degrees
    
# Create a list of tuples
trajectory_data = [(t, x, y, z, xE, yE, zE, p, yw) for t, x, y, z, xE, yE, zE, p, yw in \
zip(t_intervals, x_coordinates, y_coordinates, z_coordinates,\
x_ECEF_coordinates, y_ECEF_coordinates, z_ECEF_coordinates, pitch_intervals, yaw_intervals)]

# Write the data to a .txt file
with open(filename, 'w') as file:
    if (isGlobal):
        print("\nCreating ECEF trajectory.")
        for data_point in trajectory_data:
            file.write("%s %s %s %s %s %s\n" % (format(data_point[0], '.1f'), format(data_point[4], '.6f'),\
            format(data_point[5], '.6f'), format(data_point[6], '.6f'), format(data_point[7], '.4f'), format(data_point[8], '.4f')))
    else:
        print("\n Creating emplacement-relative trajectory.")
        for data_point in trajectory_data:
            file.write("%s %s %s %s %s %s\n" % (format(data_point[0], '.1f'), format(data_point[1], '.6f'),\
            format(data_point[2], '.6f'), format(data_point[3], '.6f'), format(data_point[7], '.4f'), format(data_point[8], '.4f')))
    file.write("%s\n" % (-1))
print("Trajectory data saved to '%s'\n" % (filename))
