# satellite.py
# MEW 12th Aug 2025

from dataclasses import dataclass
import numpy as np
from body import Body

@dataclass
class OrbitalProp:
    "A class for the orbital elements, using the Keplerian Parameterisation"
    semi_major_axis:  float          # semi-major axis takes values >0
    periapsis:      float          # periapsis takes values >0
    eccentricity:   float = 0         # eccentricity takes values in (0,1) for closed elliptical trajectories, 1 for parabolic trajectories and (1, inf) for hyperbolic orbits
    inclination:    float = 0          # inclination angle in radians, takes values in [0, pi/2]
    right_ascension: float = 0          # right ascension, angle in radians, in (0, 2 pi)
    true_anomaly:    float = 0          # True anomaly, in (0, 2 * pi)

    def __post_init__(self):
        if (self.semi_major_axis <= 0):
            raise ValueError(f"Invalid semi-major axis: {self.semi_major_axis}. Must be greater than 0.")
        if (self.periapsis <= 0):
            raise ValueError(f"Invalid periapsis: {self.periapsis}. Must be greater than 0.")
        if not (0 <= self.inclination <= np.pi / 2):
            raise ValueError(f"Invalid inclination: {self.inclination}. Must be in [0, pi/2].")
        if not (0 <= self.right_ascension < 2 * np.pi):
            raise ValueError(f"Invalid right ascension: {self.right_ascension}. Must be in [0, 2*pi).")
        if not (0 <= self.true_anomaly < 2 * np.pi):
            raise ValueError(f"Invalid true anomaly: {self.true_anomaly}. Must be in [0, 2*pi).")

@dataclass
class SurfaceProp:
    "A class storing properties of the surface of satellites."
    albedo:     float = 0.5                 # Albedo value should be between 0 (perfect absorber) and 1 (perfect reflector)
    emisivity:  float = 0.5                 # Emisivity value also between 0 (perfect reflector) and 1 (perfect emitter)

    def __post_init__(self):
        #check that the values are between (0,1)
        if not (0 <= self.albedo <= 1):
            raise ValueError(f"Invalid albedo value: {self.albedo}. Must be between 0 and 1.")
        if not (0 <= self.emisivity <= 1):
            raise ValueError(f"Invalid emisivity value: {self.emisivity}. Must be between 0 and 1.")
    

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