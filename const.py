#!/usr/bin/env python3

# CONVERSION CONSTANTS
KT_TO_MPS = 0.51444444444444   
LB_TO_KG = 0.45359237          
M_TO_FT = 3.280839895          
MPS_TO_MPH = 2.2369362920544025

# PHYSICS CONSTANTS
g = 9.81                # Gravitational acceleration constant in m/s^2

# AERODYNAMIC CONSTANTS
CD_HD = 0.8             # Drag Coefficient of a Head Down flyer (see https://tinyurl.com/ya99r6b4)                
CD_BELLY = 1.0          # Estimated drag coefficient of a belly flyer                                             
A_HD = 0.35             # Head Down surface area (see http://labman.phys.utk.edu/phys221core/modules/m3/drag.html)
A_BELLY = 16.0*68.0 / 144.0 / ((M_TO_FT)**2) * 0.7 # Estimated Human body surface area - see NOTE 1               

# MISC CONSTANTS
BLAZE_IT = 420

# NOTE 1:                                                                      
# Human body surface area -- Estimated by taking my width*height in inches,    
# converting to feet, then converting to meters and multiplying by an estimate 
# of the percentage of the rectangular surface that would be filled by a human 
# body falling with an arch                                                    
