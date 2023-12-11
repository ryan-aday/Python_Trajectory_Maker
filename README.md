# Python_Trajectory_Maker
Creates files based off of a trajectory moving about a fixed emplacement.
Note: Currently, this script creates a trajectory circling the emplacement. You can modify the .py script to modify the trajectories for your own use. 

# References
[Unit Converter](https://www.unitconverters.net/angle/degree-to-mil.html)

[Eccentricity Calcs](https://www.johndcook.com/blog/2022/10/14/eccentricity-flatness-aspect/)

[Geometric to Geodetic Relationships](https://www.oc.nps.edu/oc2902w/coord/coordcvt.pdf)

[Coordinate Transforms](https://x-lumin.com/wp-content/uploads/2020/09/Coordinate_Transforms.pdf)

[Geodetic, Geocentric, to Flat Earth Relationships](http://walter.bislins.ch/bloge/index.asp?page=Globe+and+Flat+Earth+Transformations+and+Mappings)

# Libraries
numpy, math

# To run this script: 

1. Type gedit runTraj. 

2. Scroll to the User Constant Inputs section. Input units are described in the adjacent comments. 
  * If the .trj will use a global WGS84 reference, set isGlobal to True. Otherwise, set isGlobal to False. 
filename will be set to the character array you feed it.  
  * I.e. 'traj_check.trj' will be the file name assigned to the output file. 
  * wgs84_lat, wgs84_lon, and wgs84_alt will be assigned to the latitude, longitude, and altitude in wgs-84 format. Please use standard wgs-84 units when assigning these values. 

3. Target Constants are self-explanatory. 
  * search_sector_angle is the full angle of the radar search sector, including the alternate search sector (typically 120 degrees) 
  * bearing_N is the bearing of the radar relative to true north. 
  * speed is the fixed average speed of the SRN. This may need to be modified if the speed is variable. 
  * altitude is the average altitude of the SRN. This may need to be modified if the altitude is variable. 
  * start_time and end_time are the start and end times of the scenario, respectively. 

4. Extra Constants are also self explanatory. These populate columns 5-9 of the .trj file, and are currently configured to be static. 
  * pitch and yaw are the pitch and yaw of the trajectory per instance, respectively. 

5. Scroll to the user defined X, Y, and Z trajectory characteristics settings. 
  * Note that this is currently set up for a trajectory encircling a stationary radar at a fixed altitude. This can be modified. 
  * X is the cross range, Y is the altitude, Z is the downrange. 

6. Save the .py file and close it. 

7. Type python runTraj.py and press Enter. A message displaying what you saved your file as will pop up. 
