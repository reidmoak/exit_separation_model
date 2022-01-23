#!/usr/bin/env python3

import const
import os
import sys
from termcolor import colored
from types import ModuleType

import winds

class Params:

    def __init__(self):
        self.EXIT_ALT        = 13500       # Exit altitude in feet
        self.BREAKOFF_ALT    = 5000        # Breakoff altitude in feet
        self.PULL_ALT        = 3500        # Pull altitude in feet
        self.IDEAL_SEP       = 1000        # Ideal exit separation (1000 feet)
        self.V_upper         = 10          # Uppers in knots
        self.V_air           = 70          # Jump run airspeed in knots 
        self.jump_run        = 270         # Jump run direction in degrees
        self.t_sep           = 10          # Exit separation in seconds
        self.num_rw_groups   = 3           # Number of belly groups on the load
        self.num_ff_groups   = 4           # Number of freefly groups on the load
        self.weight          = 160         # Average jumper exit weight in pounds
        
        self.setup()


    # TODO: Have this print out units, explanations of the variables, etc.
    def show(self):
        # NOTE: Unfortunately, not sure how I would make this not be hard-coded
        var_help = {}
        var_help['EXIT_ALT'] = "Exit altitude in feet"
        var_help['BREAKOFF_ALT'] = "Breakoff altitude in feet"
        var_help['PULL_ALT'] = "Pull altitude in feet"
        var_help['IDEAL_SEP'] = "Ideal exit separation in feet (default: 1000)"
        var_help['weight'] = "Average jumper exit weight in pounds"
        var_help['V_upper'] = "Average winds aloft uppers in knots"
        var_help['V_air'] = "Jump run aircraft airspeed in knots"
        var_help['jump_run'] = "Jump run direction in degrees"
        var_help['t_sep'] = "Exit separation in seconds"
        var_help['num_rw_groups'] = "Number of belly groups on the load"
        var_help['num_ff_groups'] = "Number of freefly groups on the load"
        print("")
        for i, key in enumerate(self.__dict__):
            print(f"{str(i+1) + ') ' + key:<20}{self.__dict__[key]:<10}{var_help[key]:<30}")
        print("")

    def setup(self):
        os.system('clear')
        keys = list(self.__dict__.keys())
        while(True):
            self.show()
            ans = input("Enter number of variable to modify, or \'q\' to quit: ") 
            if ans == 'q':
                break
            else:
                print("")
                try:
                    index = int(ans)-1
                    ans = input("Enter new value for " + keys[index] + ": ")
                    while(True):
                        try:
                            self.__dict__[keys[index]] = int(ans)
                            break
                        except ValueError:
                            os.system('clear')
                            print(colored("Invalid entry, input must be a number.\n", 'red'))
                            ans = input("Enter new value for " + keys[index] + ": ")
                except ValueError:
                    os.system('clear')
                    print(colored("Invalid entry, input must be a number.", 'red'))
                    continue
                except IndexError:
                    os.system('clear')
                    print(colored("Invalid entry, input must be between 1 and " + 
                          str(len(keys)) + ".", 'red'))
                    continue
            os.system('clear')
