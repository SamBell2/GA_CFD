import matplotlib.pyplot as plt

results = []
with open("record.csv") as f:
    for i, line in enumerate(f.readlines()):
        splitline = line.split(',')
        line_results = sorted([float(x) for x in splitline])
        length = len(line_results)
        if length % 2 == 0:
            median = line_results[length//2 - 1] + line_results[length//2]
            median /= 2
            if length % 4 == 0:
                q1 = line_results[length//4 - 1] + line_results[length//4]
                q1 /= 2
                q3 = (
                    line_results[(length//4)*3 - 1]
                    + line_results[length//4 * 3]
                )
                q3 /= 2
            else:
                q1 = line_results[length//4]
                q3 = line_results[length//4*3]
        else:
            median = line_results[length//2]
            if (length//2) % 2 == 0:
                q1 = (
                    line_results[:length//2][length//4]
                    + line_results[:length//2][length//4 - 1]
                )
                q1 /= 2
                q3 = (
                    line_results[length//2+1:][length//4]
                    + line_results[length//2+1:][length//4 - 1]
                )
                q3 /= 2
            else:
                q1 = line_results[length//4]
                q3 = line_results[length//4*3]
        results.append({
            "generation": i,
            "results": line_results,
            "max": max(line_results),
            "min": min(line_results),
            "mean": sum(line_results)/len(line_results),
            "q1": q1,
            "median": median,
            "q3": q3
        })
mean_x_points = [x["mean"] for x in results]
max_x_points = [x["max"] for x in results]
min_x_points = [x["min"] for x in results]
q1_x_points = [x["q1"] for x in results]
median_x_points = [x["median"] for x in results]
q3_x_points = [x["q3"] for x in results]
y_points = [x["generation"] for x in results]
plt.xlabel("Generation")
plt.ylabel("L/D")
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.plot(y_points, mean_x_points, label="Mean")
plt.plot(y_points, max_x_points, label="Maximum")
# plt.plot(y_points, min_x_points, label="Minimum")
plt.plot(y_points, median_x_points, label="Median")
plt.fill_between(
    y_points,
    q1_x_points,
    q3_x_points,
    alpha=0.5,
    label="Interquartile range"
)
plt.legend()
plt.show()
