import csv
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

data = []
filename = ""

choice = int(input("Enter choice of plane: "))

if choice == 0:
    filename = "suzanne.csv"
elif choice == 1:
    filename = "alkonost.csv"
elif choice == 2:
    filename = "super_dart.csv"
elif choice == 3:
    filename = "chinese.csv"

# Read CSV
with open(filename, "r") as f:
    reader = csv.reader(f)
    for row in reader:
        data.append([float(x) for x in row])

n = len(data[0])

x = [row[0] for row in data]
z = [row[1] for row in data]
y = [row[2] for row in data]

fig = plt.figure()
ax = fig.add_subplot(projection="3d")

if n == 4:
    c = [row[3] for row in data]
    sc = ax.scatter(x, y, z, c=c, cmap="viridis")
    plt.colorbar(sc, label="Mass?")
else:
    ax.scatter(x, y, z)

ax.set_xlabel("Aspect Ratio")
ax.set_zlabel("Range")
ax.set_ylabel("Initial Angle")

plt.show()