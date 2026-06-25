#!/usr/env python3

# Evolving take 4 
# Random: random NACA airfoils, alpha, Re_chord
# Eval:   solver L/D
# Breed:  average max_camber_dev, loc, thickness, alpha, Re_chord
# Record: max, minimum, mean

# IMPORTS
from airfoil import Airfoil
import neuralfoil as nf
import random


def evolve(population_count, generations, elite_count, breed_count,
           random_function, eval_function, breed_function):
    population = [random_function() for _ in range(population_count)]
    new_population = []
    print("started")
    try:
        for generation in range(generations):
            results = []
            for x in population:
                results.append(eval_function(x))
                # print("Result made")
            with open("take4/record.csv", 'a') as f:
                to_write = ','.join([f"{num:30.28f}" for num in results]) + '\n'
                f.write(to_write)
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
    with open("take4/record.csv", 'a') as f:
        to_write = ','.join([f"{num:30.28f}" for num in results]) + '\n'
        f.write(to_write)
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
    if aero["CD"] == 0:
        return 0
    return float(aero["CL"][0])/float(aero["CD"][0])


def breedFunction(s1, s2):
    f1 = s1[0]
    f2 = s2[0]
    avg_dev = (f1.max_camber_dev + f2.max_camber_dev) / 2
    avg_loc = (f1.max_camber_loc + f2.max_camber_loc) / 2
    avg_thick = (f1.thickness + f2.thickness) / 2
    while random.random() < 0.1:
        avg_dev += (0.5 - random.random()) / 15
    while random.random() < 0.1:
        avg_loc += (0.5 - random.random()) / 15
    while random.random() < 0.1:
        avg_thick += (0.5 - random.random()) / 15
    f3 = Airfoil(avg_dev, avg_loc, avg_thick)
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
    best = evolve(10, 175, 1, 5, randomFunction, evalFunction, breedFunction)
    print("Got best")
    best[0].name = "points"
    print("Renamed best")
    print(best[0])
    print("plotted best")
    best[0].savePlot("take4")
    print("saved diagram")
    best[0].savePoints("take4")
    print(best[1], best[2])
    with open("take4/summary.txt", 'w') as f:
        f.write(f"""Angle of attack = {best[1]} degrees
Reynolds number = {best[2]}
Max camber deviation = {best[0].max_camber_dev}
Max camber location = {best[0].max_camber_loc}
Thickness = {best[0].thickness}""")


# RUN
if __name__ == "__main__":
    main()
