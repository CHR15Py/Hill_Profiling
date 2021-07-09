# Hill_Profiling

# The aim of the project was to analyse the shape of a hill's profile i.e. how it appeared on the horizon from a given point
# in order to investigate whether the hill's Old English name was a description of it's shape.

# I was given the DTM function which provided an approximation of altitude at any integer coordinates.

# First step was to write a function that determined whether or not two points were intervisible (intervis).
# I did this in two ways. The best way (in terms of speed and accuracy) was to create a grid with lines 5m apart and then find which grid squares 
# were directly underneath the line between the two points (grid_squ). Then I used the function interp_surf to create a surface function for each grid square and used 
# line_surf_inter to dtermine whether the line between the two points intercepted the surface of the grid square. If the line did not intercept any of the grid suqare surfaces, 
# the two points were intervisible.

# The other method recursively checked whether the height of the midpoint was above the surface, but this turned out to be slower and less accurate given I could only check
# the altitude at integer coordinates.

# The hor_angle function gave the angle to the horizon by checking the visibility of point at regularly reduced angles, until a point wasn't visible, and then used
# the bisection method to increase accuracy.

# Using hor_angle in multiple directions, I could get the profile of the hill and then analyse the shape.
