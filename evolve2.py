#!/usr/env python3

# Evolving take 2
# Random: random NACA airfoils, alpha, Re_chord
# Eval:   solver L/D
# Breed:  average points, alpha, Re_chord

# IMPORTS
from airfoil import Airfoil
import neuralfoil as nf
import random
import numpy as np


def evolve(population_count, generations, elite_count, breed_count,
           random_function, eval_function, breed_function):
    population = [random_function() for _ in range(population_count)]
    new_population = [x for x in population]
    try:
        for generation in range(generations):
            results = []
            for x in population:
                results.append(eval_function(x))
                print("Result made")
            copy_results = [x for x in results]
            for _ in range(elite_count):
                index = results.index(max(results))
                new_population.append(population[index])
                results[index] = -9999999
            breeders = []
            for _ in range(breed_count):
                index = copy_results.index(max(copy_results))
                breeders.append(population[index])
                copy_results[index] = -9999999
            while len(new_population) < population_count:
                new_population.append(
                    breed_function(
                        random.choice(breeders),
                        random.choice(breeders)
                    )
                )
            population = new_population
            new_population = []
            # print(population[0])
            print(eval_function(population[0]))
            print(generation)
    except KeyboardInterrupt:
        print("Cancelling")
        results = [eval_function(x) for x in population]
    else:
        print("Final analysis")
        results = [eval_function(x) for x in population]
    print(max(results))
    return population[results.index(max(results))]


def randomFunction():
    while True:
        try:
            foil = Airfoil(random.random(), random.random(), random.random())
            break
        except ZeroDivisionError:
            pass
    foil.create_points()
    return (foil, random.randint(-20, 20), random.randint(int(1e5), int(1e7)))


def evalFunction(system):
    foil = system[0]
    aero = nf.get_aero_from_coordinates(
        coordinates=foil.array,
        alpha=system[1],
        Re=system[2],
    )
    return (aero["CL"]/aero["CD"])[0]


def breedFunction(s1, s2):
    f1 = s1[0]
    f2 = s2[0]
    x1 = f1.x_points
    x2 = f2.x_points
    y1 = f1.y_points
    y2 = f2.y_points
    x3 = [(x1[i]+x2[i])/2 for i in range(len(x1))]
    y3 = [(y1[i]+y2[i])/2 for i in range(len(y1))]
    for i in range(len(x3)):
        while random.random() < 0.05:
            x3[i] += (0.5 - random.random())/100
        while random.random() < 0.05:
            y3[i] += (0.5 - random.random())/100
    f3 = Airfoil()
    f3.x_points = np.array(x3)
    f3.y_points = np.array(y3)
    f3.create_points()
    alpha_deg = (s1[1] + s2[1]) / 2
    while random.random() < 0.05:
        alpha_deg += random.random() - 0.5
    Re_chord = (s1[2] + s2[2]) / 2
    while random.random() < 0.05:
        Re_chord += random.randint(int(-1e5), int(1e5))
    return (f3, alpha_deg, Re_chord)


# MAIN
def main() -> None:
    best = evolve(500, 150, 7, 250, randomFunction, evalFunction, breedFunction)
    print("Got best")
    best[0].name = "points"
    print("Renamed best")
    print(best[0])
    print("plotted best")
    best[0].savePlot("take2")
    print("saved diagram")
    best[0].savePoints("take2")
    print(best[1], best[2])
    with open("take2/summary.txt", 'w') as f:
        f.write(f"""Angle of attack = {best[1]} degrees
Reynolds number = {best[2]}""")


# RUN
if __name__ == "__main__":
    main()
