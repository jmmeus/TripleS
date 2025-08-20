
import numpy as np
from satellite import Satellite
from body import Body

class System(Body):

    def __init__(self, info : float | list[Satellite] | bool) -> None:
        """Initializes the System with either a total mass, a list of satellites, or a boolean for random satellite masses.
        Args:
            info (float, list, bool): 
                - If float, it represents the total mass of the system.
                - If list, it should contain Satellite objects.
                - If bool, it should be True to generate random satellite masses.
        Raises:
            TypeError: If info is not a float, list, or boolean.
            ValueError: If info is a list and empty.
        """
        
        self.info = info
         
        if isinstance(info, float):
            # if info is a float, then it is the total mass of the system
            self.total_mass = info
            self.satellites: list = self.satellites_from_total_mass(info)

        elif isinstance(info, list):
            # if info is a list, then it is a list of satellites
            if not all(isinstance(sat, Satellite) for sat in info):
                raise TypeError("All elements in the list must be Satellite objects")
            
            if len(info) == 0:
                raise ValueError("List of satellites cannot be empty")

            self.satellites = info
            self.total_mass = np.sum(self.satellites)

        elif isinstance(info, bool):
            # if info is a boolean, then it is True and we assign random satellite masses
            if info:
                self.satellites = self.random_mayhem()
                self.total_mass = np.sum(self.satellites)
            else:
                raise TypeError("Boolean input must be True in order to assign random satellite masses")
        
        else:
            # Throw informative error to explain what each type is for
            raise TypeError("Input must be a float (total mass), a list (satellites), or a boolean (True for random satellites)")
                
    @staticmethod
    def satellites_from_total_mass(total_mass: float) -> list:
        #creates a list of satellites that add up to the total mass inputted value
        satellite_list = []
        tot_mass = total_mass - np.sum(satellite_list)
        satellite = np.random(0, tot_mass)

        satellite_list.append(satellite)
        
        return satellite_list

    @staticmethod
    def random_mayhem() -> list:
        #random draw from a gaussian distribution (copy across from my code)

        #gaussian parameters
        mean = 0.0       # center of the gaussian
        std_dev = 1.0    # standard deviation
        num_sats = 1000  # number of satellites in system
        
        satellite_list = np.abs(np.random.normal(loc=mean, scale=std_dev, size=num_sats))
       
        return satellite_list 