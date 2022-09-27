from labware import plates

test_plate = plates.Plate384Echo()
test_plate.fill_plate("Coke", 3000)
test_plate.fill_well("B1", "Fanta", 5000)
print(test_plate)
print(test_plate)
test_plate.empty_well("C1")
print(test_plate.well("C1"))
print(test_plate)
test_plate.to_csv("test.csv", save_volumes=True)
empty_plate = plates.Plate384Echo()
empty_plate.from_csv("test.csv", vol=10000)
print(empty_plate)
empty_dict = empty_plate.to_dict()
assert type(empty_dict) is dict
print(empty_dict)
test_plate.fill_span("E5", "H8", "M5", 2000)
