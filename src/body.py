from abc import ABC, abstractmethod

class Body(ABC):

    def __init__(self, mass: float, radius: float):
        self.mass = mass
        self.radius = radius

    @abstractmethod
    def get_temp(self) -> float:
        pass

    @abstractmethod
    def get_coords(self) -> tuple:
        pass
