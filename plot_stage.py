import matplotlib.pyplot as plt

L, r = map(int, input().split())

sx, sy = map(int, input().split())
plt.plot(sx, sy, "o")
plt.plot([sx + L, sx - L], [sy, sy], 'k-')
gx, gy = map(int, input().split())
plt.plot(gx, gy, "o")

n = int(input())
for i in range(n):
    x1, y1, x2, y2 = map(int, input().split())
    plt.plot([x1, x2], [y1, y2], 'k-', color = "red")

plt.show()
