import neuralfoil as nf
import numpy as np
import airfoil
import solver
import time
from random import random


foil = airfoil.Airfoil(0.02, 0.4, 0.12)
points = []
for x, y in zip(foil.x_points[::-1], foil.y_points[::-1]):
    points.append((x, y))
array = np.array(points)

start = time.time_ns()
for _ in range(500):
    foil = airfoil.Airfoil(random(), random(), random())
    foil.create_points()
    aero = nf.get_aero_from_coordinates(  # You can use xy airfoil coordinates as an entry point
        coordinates=foil.array,
        alpha=5,  # Vectorize your evaluations across `alpha` and `Re`
        Re=5e6,
    )
    print(float(aero["CL"][0])/float(aero["CD"][0]))
print(time.time_ns() - start)
print(solver.calculateLD(foil.x_points, foil.y_points, 5, 5e6)[0])