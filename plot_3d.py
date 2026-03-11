import csv
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

data = []

# Read CSV
with open("suzanne.csv", "r") as f:
    reader = csv.reader(f)
    for row in reader:
        data.append([float(x) for x in row])

n = len(data[0])

x = [row[0] for row in data]
y = [row[1] for row in data]
z = [row[2] for row in data]

fig = plt.figure()
ax = fig.add_subplot(projection="3d")

if n == 4:
    c = [row[3] for row in data]
    sc = ax.scatter(x, y, z, c=c, cmap="viridis")
    plt.colorbar(sc, label="Mass?")
else:
    ax.scatter(x, y, z)

ax.set_xlabel("Aspect Ratio")
ax.set_ylabel("Range")
ax.set_zlabel("Initial Angle")

plt.show()