import matplotlib.pyplot as plt

n = int(input())
for _ in range(n):
    f, ax = plt.subplots(1, 1, figsize=(8, 6))
    m = int(input())
    for i in range(m):
        x, y, flag = map(float, input().split())
        if flag == 1:
            ax.plot(x, y, "o", color = "red")
        else:
            ax.plot(x, y, "o", color = "green")

    m = int(input())
    for i in range(m):
        x1, y1, x2, y2 = map(float, input().split())
        ax.plot([x1, x2], [y1, y2], linestyle= "solid", color = "blue")

    m = int(input())
    for i in range(m):
        x1, y1, x2, y2 = map(float, input().split())
        ax.plot([x1, x2], [y1, y2], linestyle= "dashed", color = "cyan")

    m = int(input())
    for i in range(m):
        x1, y1, x2, y2 = map(float, input().split())
        ax.plot([x1, x2], [y1, y2], linestyle= "dotted", color = "crimson")

    plt.show();
