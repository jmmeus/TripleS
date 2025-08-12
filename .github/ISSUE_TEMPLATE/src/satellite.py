# satellite.py
# MEW 12th Aug 2025

from dataclasses import dataclass
from .body import Body

@dataclass
class SurfaceProp:
    pass

@dataclass
class OrbitalProp:
    pass


class Satellite(Body):
    
    def __init__(self, mass, radius):

        self.mass = mass
        self.radius = radius

    def get_temp(self):

        # Placeholder for temperature calculation logic
        return 0.0

    def get_coords(self):

        # Placeholder for coordinates calculation logic
        return (0.0, 0.0, 0.0)