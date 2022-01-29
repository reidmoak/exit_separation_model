#!/usr/bin/env python3

# CONVERSION CONSTANTS
KT_TO_FPS = 1.6878098571012
FPS_TO_MPH = 0.68181818181818

# PHYSICS CONSTANTS
g = 32.17404855643      # Gravitational acceleration constant in ft/s^2

# AERODYNAMIC CONSTANTS
#CD_HD = 1.0             # Drag Coefficient of a Head Down flyer (see https://tinyurl.com/ya99r6b4)                
CD_HD = 0.5             # Drag Coefficient of a Head Down flyer (see https://tinyurl.com/ya99r6b4)                
CD_BELLY = 1.0          # Estimated drag coefficient of a belly flyer                                             
#A_HD = 18.0*25.0 / 144.0 # HD Surface area in ft^2
A_HD = 45.0*25.0 / 144.0 # HD Surface area in ft^2
A_BELLY = 16.0*68.0 / 144.0 * 0.7 # Estimated Human body surface area - see NOTE 1               

# MISC CONSTANTS
BLAZE_IT = 420

# NOTE 1:                                                                      
# Human body surface area -- Estimated by taking my width*height in inches,    
# converting to feet, then converting to meters and multiplying by an estimate 
# of the percentage of the rectangular surface that would be filled by a human 
# body falling with an arch                                                    
