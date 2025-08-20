from __future__ import annotations

import math
from typing import Tuple, NewType

# Type aliases for better readability and maintainability
Vector3D = Tuple[float, float, float]
Matrix3x3 = Tuple[Tuple[float, float, float], Tuple[float, float, float], Tuple[float, float, float]]

# Semantic types for orbital mechanics
Radians = NewType('Radians', float)
Eccentricity = NewType('Eccentricity', float)
MeanAnomaly = NewType('MeanAnomaly', float)
EccentricAnomaly = NewType('EccentricAnomaly', float)
HyperbolicAnomaly = NewType('HyperbolicAnomaly', float)
TrueAnomaly = NewType('TrueAnomaly', float)

# Constants
TWO_PI = 2.0 * math.pi
ATANH_CLAMP_LIMIT = 0.999999999999  # Numerical stability limit for atanh

# ---------------------------
# Angle conversions
# ---------------------------

def eccentric_anomaly_from_true(nu: TrueAnomaly, e: Eccentricity) -> EccentricAnomaly:
    """
    Convert true anomaly to eccentric anomaly for elliptical orbits.
    
    Args:
        nu: True anomaly in radians
        e: Eccentricity (must be < 1 for elliptical orbits)
        
    Returns:
        Eccentric anomaly in radians, normalized to range [0, 2π)
        
    Note:
        This function is valid only for elliptical orbits (e < 1).
        For parabolic (e = 1) or hyperbolic (e > 1) orbits, use other methods.
        
        Range Convention: The returned eccentric anomaly is normalized to [0, 2π).
        Some libraries expect angles in (-π, π]. If needed, convert using:
        E_symmetric = E - 2π if E > π else E
    """
    cos_nu, sin_nu = math.cos(nu), math.sin(nu)
    sqrt_1_e2 = math.sqrt(max(0.0, 1 - e*e))
    denom = 1 + e * cos_nu
    cos_eccentric = (e + cos_nu) / denom
    sin_eccentric = (sqrt_1_e2 * sin_nu) / denom
    eccentric = math.atan2(sin_eccentric, cos_eccentric)
    return EccentricAnomaly(eccentric if eccentric >= 0 else eccentric + TWO_PI)

def hyperbolic_anomaly_from_true(nu: TrueAnomaly, e: Eccentricity) -> HyperbolicAnomaly:
    """
    Convert true anomaly to hyperbolic anomaly for hyperbolic orbits.
    
    Args:
        nu: True anomaly in radians
        e: Eccentricity (must be > 1 for hyperbolic orbits)
        
    Returns:
        Hyperbolic anomaly
        
    Note:
        This function is valid only for hyperbolic orbits (e > 1).
        The result is clamped to maintain numerical stability.
    """
    k = math.sqrt((e - 1) / (e + 1))
    t = math.tan(nu / 2.0)
    x = max(min(k * t, ATANH_CLAMP_LIMIT), -ATANH_CLAMP_LIMIT)
    return HyperbolicAnomaly(2.0 * math.atanh(x))

# ---------------------------
# Kepler equation solvers
# ---------------------------

def solve_kepler_elliptic(
    mean_anomaly: MeanAnomaly, 
    eccentricity: Eccentricity, 
    tolerance: float = 1e-12, 
    max_iterations: int = 60
) -> EccentricAnomaly:
    """
    Solve Kepler's equation for elliptical orbits using Newton-Raphson method.
    
    Solves: M = E - e*sin(E) for E
    
    Args:
        mean_anomaly: Mean anomaly M in radians
        eccentricity: Orbital eccentricity (must be < 1)
        tolerance: Convergence tolerance for Newton-Raphson iteration
        max_iterations: Maximum number of iterations
        
    Returns:
        Eccentric anomaly E in radians
        
    Note:
        Uses Newton-Raphson method with smart initial guess based on eccentricity.
    """
    # Smart initial guess: use M for low eccentricity, π for high eccentricity
    eccentric = mean_anomaly if eccentricity < 0.8 else math.pi
    
    for _ in range(max_iterations):
        f = eccentric - eccentricity * math.sin(eccentric) - mean_anomaly
        fp = 1 - eccentricity * math.cos(eccentric)
        delta_eccentric = -f / fp
        eccentric += delta_eccentric
        if abs(delta_eccentric) < tolerance:
            break
    return EccentricAnomaly(eccentric)

def solve_kepler_hyperbolic(
    mean_anomaly: MeanAnomaly, 
    eccentricity: Eccentricity, 
    tolerance: float = 1e-12, 
    max_iterations: int = 80
) -> HyperbolicAnomaly:
    """
    Solve Kepler's equation for hyperbolic orbits using Newton-Raphson method.
    
    Solves: M = e*sinh(H) - H for H
    
    Args:
        mean_anomaly: Mean anomaly M
        eccentricity: Orbital eccentricity (must be > 1)
        tolerance: Convergence tolerance for Newton-Raphson iteration
        max_iterations: Maximum number of iterations
        
    Returns:
        Hyperbolic anomaly H
        
    Note:
        Uses Newton-Raphson method with initial guess from inverse hyperbolic sine.
    """
    # Initial guess using inverse hyperbolic sine, with fallback for M=0
    hyperbolic = math.asinh(mean_anomaly / eccentricity) if mean_anomaly != 0 else 1e-3
    
    for _ in range(max_iterations):
        f = eccentricity * math.sinh(hyperbolic) - hyperbolic - mean_anomaly
        fp = eccentricity * math.cosh(hyperbolic) - 1
        delta_hyperbolic = -f / fp
        hyperbolic += delta_hyperbolic
        if abs(delta_hyperbolic) < tolerance:
            break
    return HyperbolicAnomaly(hyperbolic)

# ---------------------------
# Rotation utilities
# ---------------------------

def R3(theta: Radians) -> Matrix3x3:
    """
    Create a rotation matrix for rotation about the Z-axis (right-hand rule).
    
    Args:
        theta: Rotation angle in radians (positive for counter-clockwise rotation)
        
    Returns:
        3x3 rotation matrix for Z-axis rotation
        
    Example:
        >>> R3(math.pi/2)  # 90 degree counter-clockwise rotation about Z-axis
    """
    c, s = math.cos(theta), math.sin(theta)
    return (
        (c, -s, 0.0),
        (s,  c, 0.0),
        (0.0, 0.0, 1.0)
    )

def R1(theta: Radians) -> Matrix3x3:
    """
    Create a rotation matrix for rotation about the X-axis (right-hand rule).
    
    Args:
        theta: Rotation angle in radians (positive for counter-clockwise rotation)
        
    Returns:
        3x3 rotation matrix for X-axis rotation
        
    Example:
        >>> R1(math.pi/2)  # 90 degree counter-clockwise rotation about X-axis
    """
    c, s = math.cos(theta), math.sin(theta)
    return (
        (1.0, 0.0, 0.0),
        (0.0, c, -s),
        (0.0, s,  c)
    )

def matmul3(matrix: Matrix3x3, vector: Vector3D) -> Vector3D:
    """
    Multiply a 3x3 matrix by a 3D vector.
    
    Args:
        matrix: 3x3 transformation matrix
        vector: 3D vector to transform
        
    Returns:
        Transformed 3D vector
        
    Note:
        Performs matrix-vector multiplication: result = matrix @ vector
    """
    return (
        matrix[0][0]*vector[0] + matrix[0][1]*vector[1] + matrix[0][2]*vector[2],
        matrix[1][0]*vector[0] + matrix[1][1]*vector[1] + matrix[1][2]*vector[2],
        matrix[2][0]*vector[0] + matrix[2][1]*vector[1] + matrix[2][2]*vector[2],
    )

def matmul33(matrix_a: Matrix3x3, matrix_b: Matrix3x3) -> Matrix3x3:
    """
    Multiply two 3x3 matrices.
    
    Args:
        matrix_a: First 3x3 matrix (left operand)
        matrix_b: Second 3x3 matrix (right operand)
        
    Returns:
        Result of matrix_a @ matrix_b
        
    Note:
        Order matters in matrix multiplication: AB ≠ BA in general
    """
    return (
        (
            sum(matrix_a[0][k] * matrix_b[k][0] for k in range(3)),
            sum(matrix_a[0][k] * matrix_b[k][1] for k in range(3)),
            sum(matrix_a[0][k] * matrix_b[k][2] for k in range(3))
        ),
        (
            sum(matrix_a[1][k] * matrix_b[k][0] for k in range(3)),
            sum(matrix_a[1][k] * matrix_b[k][1] for k in range(3)),
            sum(matrix_a[1][k] * matrix_b[k][2] for k in range(3))
        ),
        (
            sum(matrix_a[2][k] * matrix_b[k][0] for k in range(3)),
            sum(matrix_a[2][k] * matrix_b[k][1] for k in range(3)),
            sum(matrix_a[2][k] * matrix_b[k][2] for k in range(3))
        )
    )
