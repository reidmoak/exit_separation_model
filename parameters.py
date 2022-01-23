#!/usr/bin/env python3

import const
import os
from termcolor import colored
from types import ModuleType

class Params:

    def __init__(self):
        self.EXIT_ALT        = 13500       # Exit altitude in feet
        self.BREAKOFF_ALT    = 5000        # Breakoff altitude in feet
        self.PULL_ALT        = 3500        # Pull altitude in feet
        self.IDEAL_SEP       = 1000        # Ideal exit separation (1000 feet)

        self.SIM_TIME        = 120         # Simulation time in seconds

        self.weight          = 160         # Jumper exit weight in pounds

        self.V_upper         = 10          # Uppers in knots
        self.V_air           = 70          # Jump run airspeed in knots 

        self.t_sep           = 10          # Exit separation in seconds
        self.num_groups      = 7           # Number of groups on the load

        self.setup()


    # TODO: Have this print out units, explanations of the variables, etc.
    def show(self):
        print("")
        for i, key in enumerate(self.__dict__):
            print(f"{str(i+1) + ') ' + key:<30}{self.__dict__[key]:<40}")
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
