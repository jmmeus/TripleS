
import numpy as np

class System(Body):


    def __init__(self, info):
         self.info = info
         
         if not isinstance(info, float):
            
            if not isinstance(info, list):

                if info == True:
                        raise TypeError("Boolean input must be True in order to assign random satellite masses")

                if not isinstance(info, bool):

                    raise TypeError("Input parameter must be a total mass (float), a list of satellties (list) or True (boolean) to assign random satellites)")
                    
                

         
    def total_mass(self):
        #creates a list of satellites that add up to the total mass inputted value 
        satellite_list = []
        tot_mass = info - np.sum(satellite_list)
        satellite = np.random(0, tot_mass)

        satellite_list.append(satellite)
        
        return satellite_list


    def sat_list(self):

        satellite_list = info 

        return satellite_list
    

    def random_mayhem(self):
        #random draw from a gaussian distribution (copy across from my code) 

        #gaussian parameters
        mean = 0.0       # center of the gaussian
        std_dev = 1.0    # standard deviation
        num_sats = 1000  # number of satellites in system
        
        satellite_list = np.abs(np.random.normal(loc=mean, scale=std_dev, size=num_sats))
       
        return satellite_list 
    


    if info == float:
        satellite_list = total_mass(self)

    if info == list:
        satellite_list == sat_list(self)

    if info == True:
        satellite_list = random_mayhem(self)

    

    
    
    
  