# satellite.py
# MEW 12th Aug 2025

from dataclasses import dataclass
import param
import numpy as np
from body import Body

@dataclass
class OrbitalProp:
    "A class for the orbital elements, using the Keplerian Parameterisation"
    eccentricity:   float                                                           # eccentricity takes values in (0,1) for closed elliptical trajectories, 1 for parabolic trajectories and (1, inf) for hyperbolic orbits
    semiMajorAxis:  float                                                           # semi-major axis takes values >0
    periapsis:      float                                                           # periapsis takes values >0
    inclination:    float = param.Number(default=0, bounds=(0, np.pi / 2))          # inclination angle in radians, takes values in [0, pi/2]
    rightAscension: float = param.Number(default=0, bounds=(0, 2 * np.pi))          # right ascension, angle in radians, in (0, 2 pi)
    true_anomaly:   float = param.Number(default=0, bounds=(0, 2 * np.pi))          # True anomaly, in (0, 2 * pi)

@dataclass
class SurfaceProp:
    "A class storing properties of the surface of satellites."
    albedo:     float = param.Number(default=0.5, bounds=(0,1))                  # Albedo value should be between 0 (perfect absorber) and 1 (perfect reflector)
    emisivity:  float = param.Number(default=0.5, bounds=(0,1))                  # Emisivity value also between 0 (perfect reflector) and 1 (perfect emitter)
    

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