# SKYSURF-Completeness-GUI
Completeness calculator to find the 50% completeness limits and completeness percentages


Magnitude Calculation:
Given the object size, background, and exposure time for a filter, this tool calculates the 50% completeness magnitude.


Probability Calculation:
Given the object size, background, and exposure time for a filter, this tool calculates the completeness percentage at a specific object magnitude. We accomplish this by calculating the 50% completeness magnitude for the object and comparing it with the measured 50% completeness extrapolated to the object's size. Subtracting these two completeness magnitudes gives the delta which we shift the calculated completeness magnatudes to follow the data for non 50% completeness values. We arrive at an estimation for completeness percentage that is relatively accurate.


Packages: Requires sys, PyQt5, numpy, and pandas
