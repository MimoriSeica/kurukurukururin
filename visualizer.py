import matplotlib.pyplot as plt
import matplotlib.patches as patches

fig = plt.figure()
ax = plt.axes()

while True:
    str = input()
    if str == "END":
        break

    if str == "SEGMENT":
        pass

    if str == "LINE":
        pass

    if str == "CYCLE":
        X, Y, R = map(float, input().split())
        ax.add_patch(patches.Circle(xy=(X, Y), radius=R, ec='b', fill=False))

    if str == "FAN-SHAPED":
        X, Y, R, theta1, theta2 = map(float, input().split())
        ax.add_patch(patches.Wedge(center=(X, Y), r=R ,theta1 = theta1, theta2 = theta2, ec='b', fill=False))

    if str == "POLYGON":
        n = int(input())
        vertexs = tuple()
        for i in range(n):
            xy = map(float, input().split())
            vertexs += tuple(xy)

        ax.add_patch(plt.Polygon(vertexs, ))

    if str == "SAVE":
        fileName = input()
        plt.axis('scaled')
        plt.savefig(fileName)
        plt.cla()

plt.axis('scaled')
plt.show()
