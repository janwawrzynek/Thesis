class CarPark:
    def _init_(self):
        self.cars = []  # Initialize an empty list to hold cars

    def add_car(self, car):
        self.cars.append(car)  # Add a car to the park

    def show_cars(self):
        for car in self.cars:
            print(f"Car make: {car.make}, Colour: {car.colour}, Speed: {car.speed} mph")
    




class Vehicle:
    # properties (attributes)
    def _init_(self, make_var):
        self.colour = "red"
        self.speed = 69
        self.make = make_var


    # methods (functions)
    def fizzle(self):
        print("fizlly!!!!")
        print("I am of make", self.make)
    
    def accelerate(self):
        self.speed += 10
        print("Accelerating! New speed is", self.speed)
    
    def go(self):
        print("Vroom! The", self.make, "is going at", self.speed, "mph.")

    def foo(self):
        self.speed = 100


class Car(Vehicle):
    # Inherits from Vehicle
    def _init_(self, make_var):
        super()._init_(make_var)  # Call the parent class constructor
    
    def beep(self):
        print("Beep! Beep! This is a", self.make, "car.")

class Motorbike(Vehicle):
    # Inherits from Vehicle
    def _init_(self, make_var):
        super()._init_(make_var)  # Call the parent class constructor


car1 = Car("Toyota") # object instantiation
car1 = Car._init_("Toyota") # another way to instantiate
        # sets make to "Toyota", colour to "red", and speed to 69 
        # the only thing the user provides is the make
car2 = Car("Honda") # another object instantiation
motorbike1 = Motorbike("Yamaha") # motorbike object instantiation

car1.fizzle() # calling a method on car1
car2.fizzle() # calling a method on car2

car1.accelerate() # calling a method on car1
car1.go() # calling another method on car1


# Create a car park instance
car_park = CarPark()
car_park.add_car(car1)  # Add car1 to the car park
car_park.add_car(car2)  # Add car2 to the car park
car_park.add_car(motorbike1)  # Add motorbike1 to the car park
# Show all cars in the car park
car_park.show_cars()  # Display the cars in the car park


car1.beep()  # Calling the beep method on car1
motorbike1.beep()  # Calling the beep method on motorbike1


car1.foo()