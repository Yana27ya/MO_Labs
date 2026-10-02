import matplotlib.pyplot as plt
import numpy as np

#варіант 7
def f(x):
    return 1.75 * (x + 2.3) * (x - 1.05) * (x - 1.7) * (x - 3.0) + 1.0

#генерація точок для плавного графіка на інтервалі Свенна
x_vals = np.linspace(-1.4, -1.0, 400)
y_vals = f(x_vals)

plt.figure(figsize=(9, 5))
plt.plot(x_vals, y_vals, color='#8B7355', linewidth=2, label='f(x)')

points_x = [-1.16459, -1.23541, -1.27918, -1.30623, -1.32295]
points_y = [f(px) for px in points_x]

plt.scatter(points_x, points_y, color='#B22222', s=60, zorder=5, label='Пробні точки (x1, x2)')

plt.title('Метод золотого перерізу', color='#4A3B32', pad=15)
plt.xlabel('x', color='#4A3B32')
plt.ylabel('f(x)', color='#4A3B32')
plt.grid(color='#E6DCC5', linestyle='--')
plt.legend()

plt.show()