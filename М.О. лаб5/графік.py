import matplotlib.pyplot as plt
import numpy as np

def f(x):
    return 1.75 * (x + 2.3) * (x - 1.05) * (x - 1.7) * (x - 3.0) + 1.0

x_pts = np.array([-1.35, -1.20, -1.05])
f_pts = f(x_pts)

a1 = (f_pts[1] - f_pts[0]) / (x_pts[1] - x_pts[0])
a2 = ((f_pts[2] - f_pts[0]) / (x_pts[2] - x_pts[0]) - a1) / (x_pts[2] - x_pts[1])

def parabola(x):
    return f_pts[0] + a1 * (x - x_pts[0]) + a2 * (x - x_pts[0]) * (x - x_pts[1])

x_vals = np.linspace(-1.5, -0.9, 400)

plt.figure(figsize=(10, 6))
plt.plot(x_vals, f(x_vals), label='цільова функція f(x)', color='blue', linewidth=2)
plt.plot(x_vals, parabola(x_vals), label='апроміксація параболи y(x)', color='orange', linestyle='--')
plt.scatter(x_pts, f_pts, color='red', s=80, zorder=5, label='початкові точки x1, x2, x3')

plt.title('метод пауелла: квадратична апроміксація першої ітерації')
plt.xlabel('x')
plt.ylabel('f(x)')
plt.grid(True, linestyle=':')
plt.legend()
plt.show()