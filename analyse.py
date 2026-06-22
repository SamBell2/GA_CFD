import solver
import numpy as np

folder = input("Enter folder: ")
x_points = []
y_points = []
with open(folder + '/' + "points.csv") as f:
    for line in f.readlines():
        splitline = line.split(',')
        x_points.append(float(splitline[0]))
        y_points.append(float(splitline[1]))
x = np.array(x_points)
y = np.array(y_points)
with open(folder + '/' + "summary.txt") as f:
    lines = f.readlines()
    alpha_deg = int(lines[0].split()[-2])
    Re_chord = int(lines[1].split()[-1])
print(solver.calculateLD(x, y, alpha_deg, Re_chord)[0])
