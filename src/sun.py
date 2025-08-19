from body import Body

class Sun(Body):
    def __init__(self, name : str, temp : float) -> None:
        self.name = name
        self.temp = temp

    def hasHatOn(self) -> bool:


        if self.temp > 1000:
            status = True
        else:
            status = False
        return status
    

if __name__ == "__main__":
    sun = Sun("Sunny",1100)
