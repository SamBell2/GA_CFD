import neuralfoil as nf
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
    alpha_deg = float(lines[0].split()[-2])
    Re_chord = float(lines[1].split()[-1])
points = []
for x, y in zip(x_points[::-1], y_points[::-1]):
    points.append((x, y))
array = np.array(points)
aero = nf.get_aero_from_coordinates(
    coordinates=array,
    alpha=alpha_deg,
    Re=Re_chord,
)
print((aero["CL"]/aero["CD"])[0])
