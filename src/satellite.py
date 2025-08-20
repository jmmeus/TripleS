# satellite.py
# MEW 12th Aug 2025

from dataclasses import dataclass
from typing import Tuple
import numpy as np
from body import Body
import keplerian_orbitals as ko
from keplerian_orbitals import Vector3D, Radians, Eccentricity, TrueAnomaly, EccentricAnomaly, HyperbolicAnomaly, MeanAnomaly
import math

@dataclass
class OrbitalProp:
    "A class for the orbital elements, using the Keplerian Parameterisation"
    eccentricity:   float                                                           # eccentricity takes values in (0,1) for closed elliptical trajectories, 1 for parabolic trajectories and (1, inf) for hyperbolic orbits
    semiMajorAxis:  float                                                           # semi-major axis takes values >0
    periapsis:      float                                                           # periapsis takes values >0
    inclination:    float = 0.0                                                     # inclination angle in degrees, takes values in [0, 90]
    rightAscension: float = 0.0                                                     # right ascension, angle in degrees, in (0, 360)
    true_anomaly:   float = 0.0                                                     # True anomaly, in degrees (0, 360)

@dataclass
class SurfaceProp:
    "A class storing properties of the surface of satellites."
    albedo:     float = 0.5                                                         # Albedo value should be between 0 (perfect absorber) and 1 (perfect reflector)
    emissivity: float = 0.5                                                         # Emissivity value also between 0 (perfect reflector) and 1 (perfect emitter)
    

class Satellite(Body):

    def __init__(self, mass: float, radius: float) -> None:
        """Initializes a satellite with its mass and radius.
        Args:
            mass (float): Mass of the satellite in kg.
            radius (float): Radius of the satellite in km.
        """

        self.mass: float = mass
        self.radius: float = radius
        self.mu: float = self._gravitational_parameter_from_mass()

        self.orbital_prop: OrbitalProp = OrbitalProp(eccentricity=0.0, semiMajorAxis=1.0, periapsis=1.0)
        self.surface_prop: SurfaceProp = SurfaceProp()

    def _gravitational_parameter_from_mass(self) -> float:
        """Calculates the gravitational parameter (mu) from the mass of the satellite."""
        # Should probably add G to a constants file?
        return self.mass * 6.67430e-11  # Gravitational constant in m^3 kg^-1 s^-2

    def get_temp(self) -> float:
        """Calculate and return the temperature of the satellite.
        
        Returns:
            Temperature in Kelvin (placeholder implementation).
        """
        # Placeholder for temperature calculation logic
        return 0.0

    def _get_orbital_elements_radians(self) -> Tuple[float, float, Radians, Radians, Radians, TrueAnomaly]:
        """Convert orbital elements to radians and return as tuple."""
        a = self.orbital_prop.semiMajorAxis
        e = self.orbital_prop.eccentricity
        i = Radians(math.radians(self.orbital_prop.inclination))
        raan = Radians(math.radians(self.orbital_prop.rightAscension))
        omega = Radians(math.radians(self.orbital_prop.periapsis))
        nu0 = TrueAnomaly(math.radians(self.orbital_prop.true_anomaly))
        return a, e, i, raan, omega, nu0

    def _propagate_true_anomaly(self, t: float) -> Tuple[float, float]:
        """Propagate true anomaly and calculate orbital radius at time t."""
        a, e, _, _, _, nu0 = self._get_orbital_elements_radians()
        
        if e == 1.0:
            raise NotImplementedError("Parabolic orbits not supported.")

        if e < 1.0:
            # Elliptical orbit
            n = math.sqrt(self.mu / a**3)
            E0: EccentricAnomaly = ko.eccentric_anomaly_from_true(nu0, Eccentricity(e))
            m0 = E0 - e*math.sin(E0)
            m = m0 + n * t
            m = (m + math.pi) % (2*math.pi) - math.pi
            E = ko.solve_kepler_elliptic(MeanAnomaly(m), Eccentricity(e))
            cosE, sinE = math.cos(E), math.sin(E)
            denom = 1 - e*cosE
            sin_nu = math.sqrt(1 - e*e) * sinE / denom
            cos_nu = (cosE - e) / denom
            nu = math.atan2(sin_nu, cos_nu)
        else:
            # Hyperbolic orbit
            n = math.sqrt(self.mu / abs(a)**3)
            H0: HyperbolicAnomaly = ko.hyperbolic_anomaly_from_true(nu0, Eccentricity(e))
            m0 = e*math.sinh(H0) - H0
            m = m0 + n*t
            H = ko.solve_kepler_hyperbolic(MeanAnomaly(m), Eccentricity(e))
            coshH, sinhH = math.cosh(H), math.sinh(H)
            denom = e*coshH - 1
            sin_nu = math.sqrt(e*e - 1)*sinhH / denom
            cos_nu = (e - coshH) / denom
            nu = math.atan2(sin_nu, cos_nu)
        
        # Calculate orbital radius
        p = a * (1 - e*e)
        r_mag = p / (1 + e*math.cos(nu))
        
        return nu, r_mag

    def _calculate_perifocal_position(self, nu: float, r_mag: float) -> Vector3D:
        """Calculate position vector in perifocal coordinate frame."""
        return (r_mag*math.cos(nu), r_mag*math.sin(nu), 0.0)

    def _calculate_perifocal_velocity(self, nu: float, p: float, e: float) -> Vector3D:
        """Calculate velocity vector in perifocal coordinate frame."""
        return (
            -math.sqrt(self.mu/p) * math.sin(nu),
            math.sqrt(self.mu/p) * (e + math.cos(nu)),
            0.0
        )

    def _get_transformation_matrix(self) -> ko.Matrix3x3:
        """Get transformation matrix from perifocal to ECI coordinates."""
        _, _, i, raan, omega, _ = self._get_orbital_elements_radians()
        return ko.matmul33(ko.matmul33(ko.R3(raan), ko.R1(i)), ko.R3(omega))

    def _transform_to_eci(self, perifocal_vector: Vector3D) -> Vector3D:
        """Transform vector from perifocal to ECI coordinate frame."""
        R_total = self._get_transformation_matrix()
        return ko.matmul3(R_total, perifocal_vector)

    def get_coords(self, t: float = 0.0) -> Vector3D:
        """Return position (km) at time t after epoch."""
        nu, r_mag = self._propagate_true_anomaly(t)
        r_pf = self._calculate_perifocal_position(nu, r_mag)
        return self._transform_to_eci(r_pf)

    def get_velocity(self, t: float = 0.0) -> Vector3D:
        """Return velocity (km/s) at time t after epoch."""
        a, e, _, _, _, _ = self._get_orbital_elements_radians()
        nu, _ = self._propagate_true_anomaly(t)
        p = a * (1 - e*e)
        v_pf = self._calculate_perifocal_velocity(nu, p, e)
        return self._transform_to_eci(v_pf)


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    
    # Earth's gravitational parameter (km^3/s^2)
    MU_EARTH = 398600.4418
    
    print("Testing Satellite get_coords method with three different orbit types...")
    
    # ========================================
    # Test 1: Elliptical Orbit (e < 1)
    # ========================================
    print("\n1. Creating elliptical orbit satellite...")
    sat_elliptical = Satellite(mass=500.0, radius=1.0)  # 500 kg, 1 km radius
    sat_elliptical.mu = MU_EARTH  # Use Earth's mu instead of satellite's tiny mu
    sat_elliptical.orbital_prop.semiMajorAxis = 7000.0  # km (LEO)
    sat_elliptical.orbital_prop.eccentricity = 0.5      # moderately elliptical
    sat_elliptical.orbital_prop.inclination = 30.0      # degrees
    sat_elliptical.orbital_prop.rightAscension = 45.0   # degrees
    sat_elliptical.orbital_prop.periapsis = 90.0        # degrees (argument of periapsis)
    sat_elliptical.orbital_prop.true_anomaly = 0.0      # start at periapsis
    
    # Calculate orbital period for elliptical orbit
    period_elliptical = 2 * math.pi * math.sqrt(sat_elliptical.orbital_prop.semiMajorAxis**3 / MU_EARTH)
    print(f"   Orbital period: {period_elliptical/60:.1f} minutes")
    
    # ========================================
    # Test 2: Nearly Circular Orbit (e ≈ 0)
    # ========================================
    print("\n2. Creating circular orbit satellite...")
    sat_circular = Satellite(mass=1000.0, radius=2.0)  # 1000 kg, 2 km radius
    sat_circular.mu = MU_EARTH
    sat_circular.orbital_prop.semiMajorAxis = 42164.0   # km (geostationary altitude)
    sat_circular.orbital_prop.eccentricity = 0.001      # nearly circular
    sat_circular.orbital_prop.inclination = 0.0         # equatorial
    sat_circular.orbital_prop.rightAscension = 0.0      # degrees
    sat_circular.orbital_prop.periapsis = 0.0           # degrees
    sat_circular.orbital_prop.true_anomaly = 0.0        # start position
    
    # Calculate orbital period for circular orbit
    period_circular = 2 * math.pi * math.sqrt(sat_circular.orbital_prop.semiMajorAxis**3 / MU_EARTH)
    print(f"   Orbital period: {period_circular/3600:.1f} hours")
    
    # ========================================
    # Test 3: Hyperbolic Orbit (e > 1)
    # ========================================
    print("\n3. Creating hyperbolic orbit satellite...")
    sat_hyperbolic = Satellite(mass=100.0, radius=0.5)  # 100 kg, 0.5 km radius
    sat_hyperbolic.mu = MU_EARTH
    sat_hyperbolic.orbital_prop.semiMajorAxis = -10000.0  # km (negative for hyperbolic)
    sat_hyperbolic.orbital_prop.eccentricity = 1.5        # hyperbolic
    sat_hyperbolic.orbital_prop.inclination = 45.0        # degrees
    sat_hyperbolic.orbital_prop.rightAscension = 180.0    # degrees
    sat_hyperbolic.orbital_prop.periapsis = 270.0         # degrees
    sat_hyperbolic.orbital_prop.true_anomaly = 0.0        # start at periapsis
    
    # Calculate periapsis distance
    periapsis_distance = abs(sat_hyperbolic.orbital_prop.semiMajorAxis) * (sat_hyperbolic.orbital_prop.eccentricity - 1)
    print(f"   Periapsis distance: {periapsis_distance:.1f} km")
    
    # ========================================
    # Propagate orbits and collect positions
    # ========================================
    print("\n4. Propagating orbits...")
    
    # Time arrays for each orbit type
    t_elliptical = np.linspace(0, period_elliptical, 100)
    t_circular = np.linspace(0, period_circular, 100)
    t_hyperbolic = np.linspace(-3600, 3600, 100)  # ±1 hour around periapsis
    
    # Storage for positions
    positions_elliptical = []
    positions_circular = []
    positions_hyperbolic = []
    
    # Propagate elliptical orbit
    for t in t_elliptical:
        pos = sat_elliptical.get_coords(t)
        positions_elliptical.append(pos)
    
    # Propagate circular orbit
    for t in t_circular:
        pos = sat_circular.get_coords(t)
        positions_circular.append(pos)
    
    # Propagate hyperbolic orbit
    for t in t_hyperbolic:
        try:
            pos = sat_hyperbolic.get_coords(t)
            positions_hyperbolic.append(pos)
        except Exception as e:
            print(f"   Warning: Error at t={t:.1f}s: {e}")
            continue
    
    # ========================================
    # Create 3D visualization
    # ========================================
    print("\n5. Creating 3D visualization...")
    
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')
    
    # Plot elliptical orbit
    if positions_elliptical:
        pos_array = np.array(positions_elliptical)
        ax.plot(pos_array[:, 0], pos_array[:, 1], pos_array[:, 2], 
                'b-', linewidth=2, label=f'Elliptical (e={sat_elliptical.orbital_prop.eccentricity})')
        ax.scatter(pos_array[0, 0], pos_array[0, 1], pos_array[0, 2], 
                  c='blue', s=50, marker='o', label='Elliptical start')
    
    # Plot circular orbit
    if positions_circular:
        pos_array = np.array(positions_circular)
        ax.plot(pos_array[:, 0], pos_array[:, 1], pos_array[:, 2], 
                'g-', linewidth=2, label=f'Circular (e={sat_circular.orbital_prop.eccentricity})')
        ax.scatter(pos_array[0, 0], pos_array[0, 1], pos_array[0, 2], 
                  c='green', s=50, marker='s', label='Circular start')
    
    # Plot hyperbolic orbit
    if positions_hyperbolic:
        pos_array = np.array(positions_hyperbolic)
        ax.plot(pos_array[:, 0], pos_array[:, 1], pos_array[:, 2], 
                'r-', linewidth=2, label=f'Hyperbolic (e={sat_hyperbolic.orbital_prop.eccentricity})')
        ax.scatter(pos_array[0, 0], pos_array[0, 1], pos_array[0, 2], 
                  c='red', s=50, marker='^', label='Hyperbolic start')
    
    # Add Earth at origin (scaled for visibility)
    ax.scatter(0, 0, 0, c='cyan', s=100, marker='o', label='Earth')
    
    # Set labels and title
    ax.set_xlabel('X (km)')
    ax.set_ylabel('Y (km)')
    ax.set_zlabel('Z (km)')
    ax.set_title('Satellite Orbital Trajectories')
    ax.legend()
    
    # Set equal aspect ratio
    max_range = 0
    for positions in [positions_elliptical, positions_circular, positions_hyperbolic]:
        if positions:
            pos_array = np.array(positions)
            max_range = max(max_range, np.max(np.abs(pos_array)))
    
    ax.set_xlim([-max_range, max_range])
    ax.set_ylim([-max_range, max_range])
    ax.set_zlim([-max_range, max_range])
    
    # ========================================
    # Print validation information
    # ========================================
    print("\n6. Validation results:")
    
    # Check circular orbit radius consistency
    if positions_circular:
        pos_array = np.array(positions_circular)
        radii = np.sqrt(np.sum(pos_array**2, axis=1))
        mean_radius = np.mean(radii)
        radius_std = np.std(radii)
        print(f"   Circular orbit: mean radius = {mean_radius:.1f} km, std = {radius_std:.3f} km")
        print(f"   Expected radius: {sat_circular.orbital_prop.semiMajorAxis:.1f} km")
    
    # Check elliptical orbit periapsis/apoapsis
    if positions_elliptical:
        pos_array = np.array(positions_elliptical)
        radii = np.sqrt(np.sum(pos_array**2, axis=1))
        min_radius = np.min(radii)
        max_radius = np.max(radii)
        a = sat_elliptical.orbital_prop.semiMajorAxis
        e = sat_elliptical.orbital_prop.eccentricity
        expected_periapsis = a * (1 - e)
        expected_apoapsis = a * (1 + e)
        print(f"   Elliptical orbit: periapsis = {min_radius:.1f} km (expected {expected_periapsis:.1f})")
        print(f"   Elliptical orbit: apoapsis = {max_radius:.1f} km (expected {expected_apoapsis:.1f})")
    
    # Check hyperbolic orbit periapsis
    if positions_hyperbolic:
        pos_array = np.array(positions_hyperbolic)
        radii = np.sqrt(np.sum(pos_array**2, axis=1))
        min_radius = np.min(radii)
        print(f"   Hyperbolic orbit: minimum radius = {min_radius:.1f} km (expected ~{periapsis_distance:.1f})")
    
    print("\nTest complete! Close the plot window to exit.")
    plt.show()