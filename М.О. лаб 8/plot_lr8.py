import numpy as np
import matplotlib.pyplot as plt

def F(x1, x2):
    return 3*x1**2 + x2**2 + 2*x1*x2 + 16*x1 + 4*x2

def grad(x):
    return np.array([6*x[0] + 2*x[1] + 16, 2*x[1] + 2*x[0] + 4], dtype=float)

def direction(g, metric):
    if metric == "sphere":
        return -g / np.linalg.norm(g)
    elif metric == "cube":
        return -np.sign(g)
    else:
        j = np.argmax(np.abs(g))
        d = np.zeros_like(g)
        d[j] = -np.sign(g[j])
        return d

def get_trajectory(x0, h0, eps, metric):
    x = np.array(x0, dtype=float)
    Fx = F(x[0], x[1])
    h = h0
    pts = [x.copy()]
    for _ in range(5000):
        g = grad(x)
        if np.linalg.norm(g) <= eps:
            break
        d = direction(g, metric)
        while True:
            xb = x + h * d
            Fb = F(xb[0], xb[1])
            if Fb < Fx:
                break
            h /= 2.0
        x, Fx = xb, Fb
        pts.append(x.copy())
    return np.array(pts)

x0 = [3.0, -3.0]
t_sph = get_trajectory(x0, 1.0, 0.01, "sphere")
t_cub = get_trajectory(x0, 1.0, 0.01, "cube")
t_oct = get_trajectory(x0, 1.0, 0.01, "octa")

X1, X2 = np.meshgrid(np.linspace(-6, 5, 400), np.linspace(-5, 5, 400))
Z = F(X1, X2)

plt.figure(figsize=(10, 7))
levels = [-20, -15, -5, 10, 30, 54, 85]
cs = plt.contour(X1, X2, Z, levels=levels, colors="gray", linewidths=0.9, alpha=0.7)
plt.clabel(cs, inline=True, fontsize=8, fmt="F=%.0f")

plt.plot(t_sph[:, 0], t_sph[:, 1], "o-", color="#1f77b4", lw=2, ms=4, label=f"Сферична ({len(t_sph)-1} ітер.)")
plt.plot(t_cub[:, 0], t_cub[:, 1], "s-", color="#ff7f0e", lw=2, ms=4, label=f"Кубічна ({len(t_cub)-1} ітер.)")
plt.plot(t_oct[:, 0], t_oct[:, 1], "^-", color="#2ca02c", lw=2, ms=4, label=f"Октаедрична ({len(t_oct)-1} ітер.)")

plt.scatter([3], [-3], color="black", s=80, zorder=6, label="Старт $x^0 = (3; -3)$")
plt.scatter([-3], [1], color="red", s=90, zorder=6, label="Мінімум $x^* = (-3; 1)$")

plt.title("Градієнтні методи: лінії рівня та траєкторії (Варіант 7, $h_0 = 1$)", fontsize=13)
plt.xlabel("$x_1$", fontsize=12)
plt.ylabel("$x_2$", fontsize=12)
plt.grid(True, alpha=0.3)
plt.legend(fontsize=10)
plt.tight_layout()
plt.savefig("lr8_trajectories_var7.png", dpi=300)
plt.show()