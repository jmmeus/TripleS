# Abstract base classes for the solar system project

from abc import ABC, abstractmethod

class body(ABC):

    def __init__(self, mass: float, radius: float):
        self.mass = mass
        self.radius = radius

    @abstractmethod
    def get_temp(self) -> float:
        pass

    @abstractmethod
    def get_coords(self) -> tuple:
        pass