'''
runTraj.py
Written by Ryan Aday
Copyright @2023
V2: 12/06/2023- Added Hughes transform matrix capability, revised relative to ECEF algorithm

The Hughes transform matrix performs a coordinate transformation by rotating and translating the emplacement coordinate system into that of the globe model. 
This approach is correct.

The existing relative to ECEF algorithm converts the relative Cartesian coordinates into the relative latitude/longitude angles & altitude.
This is then added to the radar latitude/longitude angles and the altitude to get the correct ECEF coordinates of a trajectory.
This approach is wrong.
'''

# Import libraries
import math, numpy

'''
#############################################################################################
User Constant Inputs
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
distance_from_radar = 80000.0 # meters
search_sector_angle = 120.0 # degrees
bearing_N = 0.0 # 6133.3 # mills (Note: 1 degree = 160/9 milla)
speed = 210.0 # m/s
altitude = 40000.0 # meters
start_time = 0.0 # seconds
end_time = 5000.0 # seconds

# Extra Constants
pitch = 0.0 # degrees, negative is down
RCS = 10.0 # dB, for nominal frequency
low_RCS = 10.0 # dB, for subnominal frequency
high_RCS = 10.0 # dB, for supernominal frequency
yaw = 0.0 # RCS degrees

'''
#############################################################################################
Conversion Functions ***DO NOT MODIFY
#############################################################################################
'''

# Uncomment for Hughes transform
def Rz(z_angle):
    return numpy.array([\
    [math.cos(z_angle), math.sin(z_angle), 0],\
    [-math.sin(z_angle), math.cos(z_angle), 0],\
    [0, 0, 1],])

def Ry(y_angle):
    return numpy.array([\
    [math.cos(y_angle, 0, math.sin(y_angle)],\
    [0, 1, 0],\
    [-math.sin(y_angle), 0, math.cos(y_angle)],])
    
def Rx(x_angle):
    return numpy.array([\
    [1, 0, 0],\
    [0, -math.cos(x_angle), -math.sin(z.angle)],\
    [0, math.sin(x_angle), math.cos(x_angle)],])
    
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

def resM(H, x, y, z):
    ref = numpy.array([[x], [y], [z], [1]])
    # X, Y, Z should be the absolute ECEF coordinates
    return H * ref

def angle(v1, v2):
    # Assumes v1 and v2 are both 1x3 vectors, outputs in rads
    return numpy.arccos(numpy.dot(v1, v2)/(numpy.linalg.norm(v1) * numpy.linalg.norm(v2)))

def geodetic_to_geocentric(ellps, lat, lon, h):
    a, rf = ellps
    lat_rad = math.radians(lat)
    lon_rad = math.radians(lon)
    N = a / math.sqrt(1 - (1 - (1 - 1 / rf) ** 2) * (math.sin(lat_rad))**2)
    X = ((N + h) * math.cos(lat_rad) * math.cos(lon_rad))
    Y = ((N + h) * math.cos(lat_rad) * math.sin(lon_rad))    
    Z = ((1 - 1 / rf) ** 2 * N + h) * math.sin(lat_rad)
    
    return X, Y, Z
    
def relative_to_geodetic(ellps, x, y, z):
    R, rf = ellps
    
    lon = math.degrees(math.atan2(y, x))
    p = math.sqrt(x**2 + y**2)
    lat = math.degrees(math.atan2(p, z))
    h = y
    
    return lat, lon, h
    
'''
#############################################################################################
Functions
#############################################################################################
'''    

# Nodel Params
#grs80 = (6378137, 298.257222100882711)
wgs84 = (6378137, 298.257223563)

# ECEF Reference Vectors
X_ECEF_ref = [1, 0, 0]
Y_ECEF_ref = [0, 1, 0]
Z_ECEF_ref = [0, 0, 1]

# ECEF Conversion of emplacement location
ECEF_X, ECEF_Y, ECEF_Z = geodetic_to_geocentric(wgs84, wgs84_lat, wgs84_lon, wgs_84_alt)
radar_ECEF = [ECEF_X, ECEF_Y, ECEF_Z]

# Derive H Transform Matrix
x_angle = angle(radar_ECEF, X_ECEF_ref)
y_angle = angle(radar_ECEF, Y_ECEF_ref)
z_angle = angle(radar_ECEF, Z_ECEF_ref)

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
RCS_intervals = []
low_RCS_intervals = []
high_RCS_intervals = []
yaw_intervals = []

# Calculate coordinates for each time interval
for t in numpy.arange(start_time, end_time + time_interval, time_interval):
    '''
    #############################################################################################
    Start of user-defined X, Y, Z trajectory characteristics
    #############################################################################################
    '''
    
    # Calculate X, Y, Z coordinates
    x = distance_from_radar * math.sin(start_angle_offset + angular_velocity * t) # km
    y = altitude # km
    z = distance_from_radar * math.cos(start_angle_offset + angular_velocity * t) # km
    
    # LLA solution
    '''
    lat_sim, lon_sim, h_sim = relative_to_geodetic(wgs84, x, y, z)
    x_ECEF, y_ECEF, z_ECEF = geodetic_to_geocentric(wgs84, wgs84_lat + lat_sim,\
    wgs84_lon + lon_sim, wgs84_alt + h_sim)
    '''
    
    # Hughes transform solution
    sol = resM(H, x, y, z)
    x_ECEF = sol[0]
    y_ECEF = sol[1]
    z_ECEF = sol[2]

    '''
    #############################################################################################
    End of user-defined X, Y, Z trajectory characteristics
    #############################################################################################
    '''        
    
    t_intervals.append(t) # seconds
    x_coordinates.append(x/1000) # km
    y_coordinates.append(y/1000) # km
    z_coordinates.append(z/1000) # km
    x_ECEF_coordinates.append(x_ECEF) # meters
    y_ECEF_coordinates.append(y_ECEF) # meters
    z_ECEF_coordinates.append(z_ECEF) # meters
    pitch_intervals.append(pitch) # degrees
    RCS_intervals.append(RCS) # dB
    low_RCS_intervals.append(low_RCS) # dB
    high_RCS_intervals.append(high_RCS) # dB
    yaw_intervals.append(yaw) # degrees
    
# Create a list of tuples
trajectory_data = [(t, x, y, z, xE, yE, zE, p, R, lR, hR, yw) for t, x, y, z, xE, yE, zE, p, R, lR, hR, yw in \
x_ECEF_coordinates, y_ECEF_coordinates, z_ECEF_coordinates,\
pitch_intervals, RCS_intervals, low_RCS_intervals, high_RCS_intervals, yaw_intervals)]

# Write the data to a .txt file
with open(filename, 'w') as file:
    if (isGlobal):
        print("\n Creating ECEF trajectory.")
        for data_point in trajectory_data:
            file.write("%s %s %s %s %s %s %s %s %s\n" % (format(data_point[0], '.1f'), format(data_point[4], '.6f'),\
            format(data_point[5], '.6f'), format(data_point[6], '.6f'), format(data_point[7], '.4f'), format(data_point[8], '.4f'),\
            format(data_point[9], '.6f'), format(data_point[10], '.6f'), format(data_point[11], '.6f')))
    else:
        print("\n Creating emplacement-relative trajectory.")
        for data_point in trajectory_data:
            file.write("%s %s %s %s %s %s %s %s %s\n" % (format(data_point[0], '.1f'), format(data_point[1], '.6f'),\
            format(data_point[2], '.6f'), format(data_point[3], '.6f'), format(data_point[7], '.4f'), format(data_point[8], '.4f'),\
            format(data_point[9], '.6f'), format(data_point[10], '.6f'), format(data_point[11], '.6f')))
    file.write("%s\n" % (-1))
print("Trajectory data saved to '%s'\n" % (filename))
