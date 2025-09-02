from body import Body
import math

class Sun(Body):
    """Class containing basic structure components for the sun"""

    def __init__(self, name : str, temp : float, luminosity : float) -> None:
        self.name = name
        self.temp = temp
        self.luminosity = luminosity


    def get_coords(self):
        # Method that calculates and returns the coordinates of the sun. Defaults to 0,0,0
        coords = (0,0,0)
        return coords
    
    def get_temp(self):
        # Return temperature

        return self.temp
    
    def get_radius(self):
        # Use the Stefan-Boltzmann law to calculate radius based on luminosity and temperature

        sbc = 5.670374419e-8

        radius = math.sqrt(self.luminosity/(4*math.pi*sbc*(self.temp)**4))

        return radius


    def has_hat_on(self) -> bool:
        # Determines whether the sun has its hat on or not

        if self.temp > 1000:
            status = True
        else:
            status = False
        return status
    

if __name__ == "__main__":
    sun = Sun("Sunny",5778,3.828e26)

    radius = sun.get_radius()
    print(f"The radius of the sun is {radius}m")