import numpy as np
from satellite import Satellite
from body import Body

class System(Body):

    def __init__(self, info : float | list[Satellite] | bool, dt: float = 1.0):

        """Initializes the System with either a total mass, a list of satellites, or a boolean for random satellite masses.
        Args:
            info (float, list, bool): 
                - If float, it represents the total mass of the system.
                - If list, it should contain Satellite objects.
                - If bool, it should be True to generate random satellite masses.
            dt (float): Time step for the simulation, default is 1.0.
        Raises:
            TypeError: If info is not a float, list, or boolean.
            ValueError: If info is a list and empty.
        """
        self.dt = dt
        self.info = info
         
        if isinstance(info, float):
            # if info is a float, then it is the total mass of the system
            self.total_mass = info
            self.satellites: list[Satellite] = self.satellites_from_total_mass(info)

        elif isinstance(info, list):
            # if info is a list, then it is a list of satellites
            if not all(isinstance(sat, Satellite) for sat in info):
                raise TypeError("All elements in the list must be Satellite objects")
            
            if len(info) == 0:
                raise ValueError("List of satellites cannot be empty")

            self.satellites: list[Satellite] = info
            self.total_mass = np.sum([sat.mass for sat in info])

        elif isinstance(info, bool):
            # if info is a boolean, then it is True and we assign random satellite masses
            if info:
                self.satellites: list[Satellite] = self.random_mayhem()
                self.total_mass = np.sum([sat.mass for sat in self.satellites])
            else:
                raise TypeError("Boolean input must be True in order to assign random satellite masses")
        
        else:
            # Throw informative error to explain what each type is for
            raise TypeError("Input must be a float (total mass), a list (satellites), or a boolean (True for random satellites)")
                
    @staticmethod
    def satellites_from_total_mass(total_mass: float | int) -> list[Satellite]:
        # creates a list of satellites that add up to the total mass inputted value
        satellite_list: list[Satellite] = []
        if total_mass <= 0:
            raise ValueError("Total mass must be greater than zero")
        if not isinstance(total_mass, float | int):
            raise TypeError("Total mass must be a float or int")

        # Generate random satellite masses until they sum to total_mass
        while np.sum([sat.mass for sat in satellite_list]) < total_mass:
            remaining_mass: float = total_mass - np.sum([sat.mass for sat in satellite_list])
            rand_mass: float = np.random.uniform(0, total_mass)

            # Ensure the random mass does not exceed the remaining total mass
            if rand_mass > remaining_mass:
                rand_mass = remaining_mass

            # Need a radius too, let's assume a random radius between 1 and 10 for simplicity
            # TODO: Implement a more sophisticated radius generation // maybe based on mass?
            rand_radius: float = np.random.uniform(1, 10)

            satellite_list.append(Satellite(mass=rand_mass, radius=rand_radius))

        # Ensure the total mass is exactly equal to the input value
        if np.sum([sat.mass for sat in satellite_list]) != total_mass:
            raise ValueError("Generated satellites do not sum to the total mass")
        
        return satellite_list

    @staticmethod
    def random_mayhem() -> list[Satellite]:
        #random draw from a gaussian distribution (copy across from my code)

        #gaussian parameters
        mean = 0.0       # center of the gaussian
        std_dev = 1.0    # standard deviation
        num_sats = np.random.randint(1, 30)  # number of satellites in system

        satellite_list: list[Satellite] = []
        for _ in range(num_sats):
            mass = np.abs(np.random.normal(loc=mean, scale=std_dev))
            radius = np.random.uniform(1, 10)  # Random radius between 1 and 10
            satellite_list.append(Satellite(mass=mass, radius=radius))

        return satellite_list


if __name__ == '__main__':
    # Initialize system with random satellites using random_mayhem
    print("Initializing system with random satellites...")
    system = System(True)  # Use random_mayhem
    
    # Display system information
    print(f"\nGenerated {len(system.satellites)} satellites")
    print(f"Total system mass: {system.total_mass:.3f}")
    print(f"Time step (dt): {system.dt}")
    
    # Show sample satellites (first 3 for brevity)
    print("\nSample satellites:")
    for i, sat in enumerate(system.satellites[:3]):
        print(f"  Satellite {i}: mass={sat.mass:.3f}, radius={sat.radius:.3f}")
    
    if len(system.satellites) > 3:
        print(f"  ... and {len(system.satellites) - 3} more satellites")
    
    # Time simulation setup
    t = 0.0
    num_steps = 5
    print(f"\nRunning simulation for {num_steps} time steps:")
    
    # Time iteration
    for step in range(num_steps):
        t += system.dt
        print(f"\n--- Time: {t:.1f} ---")
        
        # Get coordinates for each satellite (showing first 3 for clarity)
        for i, sat in enumerate(system.satellites[:3]):
            coords = sat.get_coords(t)  # Expecting (x, y, z) tuple
            print(f"  Satellite {i}: coordinates {coords}")
        
        if len(system.satellites) > 3:
            print(f"  ... coordinates for {len(system.satellites) - 3} more satellites")
    
    print("\nSimulation complete!")