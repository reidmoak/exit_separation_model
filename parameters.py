#!/usr/bin/env python3

import os
from termcolor import colored

class Params:

    def __init__(self):
        self.EXIT_ALT        = 13500       # Exit altitude in feet
        self.BREAKOFF_ALT    = 5000        # Breakoff altitude in feet
        self.PULL_ALT        = 3500        # Pull altitude in feet
        self.IDEAL_SEP       = 1000        # Ideal exit separation (1000 feet)
        self.V_upper         = 10          # Uppers in knots
        self.weight          = 160         # Average jumper exit weight in pounds
        self.num_rw_groups   = 3           # Number of belly groups on the load
        self.num_ff_groups   = 4           # Number of freefly groups on the load
        self.aircraft        = "Caravan"   # Type of aircraft
        self.jump_run        = 270         # Jump run direction in degrees
        self.t_sep           = 10          # Exit separation in seconds

    def show(self, numbered):
        # Variable descriptions (kind of have to be hardcoded...)
        var_help = {}
        var_help['EXIT_ALT'] = "Exit altitude in feet"
        var_help['BREAKOFF_ALT'] = "Breakoff altitude in feet"
        var_help['PULL_ALT'] = "Pull altitude in feet"
        var_help['IDEAL_SEP'] = "Ideal exit separation in feet (default: 1000)"
        var_help['weight'] = "Average jumper exit weight in pounds"
        var_help['V_upper'] = "Average winds aloft uppers in knots - only used when simple_winds is True"
        var_help['aircraft'] = "Type of aircraft (affects ground speed at exit altitude)"
        var_help['jump_run'] = "Jump run direction in degrees"
        var_help['t_sep'] = "Exit separation in seconds"
        var_help['num_rw_groups'] = "Number of belly groups on the load"
        var_help['num_ff_groups'] = "Number of freefly groups on the load"

        if numbered is True:
            for i, key in enumerate(self.__dict__):
                if 'jump_run' in key or 't_sep' in key:
                    print(colored(f"{str(i+1) + ') ' + key:<20}{self.__dict__[key]:<10}{var_help[key]:<30}", 'cyan'))
                else:
                    print(f"{str(i+1) + ') ' + key:<20}{self.__dict__[key]:<10}{var_help[key]:<30}")
        else:
            for key in self.__dict__:
                if 'jump_run' in key or 't_sep' in key:
                    print(colored(f"{key:<20}{self.__dict__[key]:<10}{var_help[key]:<30}", 'cyan'))
                else:
                    print(f"{key:<20}{self.__dict__[key]:<10}{var_help[key]:<30}")
        print("")

    def setup(self):
        aircraft_list = {
            1: "Caravan",
            2: "Otter",
            3: "Skyvan",
            4: "Cessna 182"
        }
        def print_aircrafts():
            print("Aircraft options:")
            for index in aircraft_list:
                print("\t" + str(index) + ") " + aircraft_list[index])
            print("")

        os.system('clear')
        keys = list(self.__dict__.keys())

        while(True):
            print("")
            self.show(True)
            ans = input("Enter number of variable to modify, or \'q\' to quit: ") 
            if ans == 'q':
                break
            else:
                print("")
                try:
                    index = int(ans)-1
                    # Aircraft is special case, due to needing specific options
                    if "aircraft" in keys[index]:
                        print_aircrafts()
                    ans = input("Enter new value for " + keys[index] + ": ")
                    while(True):
                        try:
                            # Aircraft is special case (str versus int)
                            if "aircraft" in keys[index]:
                                if int(ans) < 1 or int(ans) > 4:
                                    os.system('clear')
                                    print(colored("Invalid aircraft entry.\n", 'red'))
                                    print_aircrafts()
                                    ans = input("Enter new value for " + keys[index] + ": ")
                                else:
                                    print(aircraft_list.get(ans))
                                    self.__dict__[keys[index]] = aircraft_list.get(int(ans))
                                    break
                            else:
                                self.__dict__[keys[index]] = int(ans)
                                break
                        except ValueError:
                            os.system('clear')
                            print(colored("Invalid entry, input must be a number.\n", 'red'))
                            if "aircraft" in keys[index]:
                                print_aircrafts()
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
        
        # For printing purposes, mod jump_run by 360 
        if self.jump_run >= 360 or self.jump_run < 0:
            self.jump_run = self.jump_run % 360
