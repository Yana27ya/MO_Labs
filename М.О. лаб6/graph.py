import numpy as np
import matplotlib.pyplot as plt

def df(x):
    return 7.0 * x**3 - 18.1125 * x**2 - 11.165 * x + 31.019625

def d2f(x):
    return 21.0 * x**2 - 36.225 * x - 11.165

x0 = -1.35
d1_0 = df(x0)
d2_0 = d2f(x0)

x1 = -1.295531
d1_1 = df(x1)
d2_1 = d2f(x1)

x2 = -1.293605

x_vals = np.linspace(-1.36, -1.27, 400)
phi_vals = df(x_vals)

plt.figure(figsize=(10, 6))
plt.plot(x_vals, phi_vals, color="#2ca02c", lw=2.5, label=r"$\varphi(x) = f'(x)$")
plt.axhline(0, color="black", lw=1)

x_tan0 = np.linspace(-1.355, -1.29, 100)
y_tan0 = d1_0 + d2_0 * (x_tan0 - x0)
plt.plot(x_tan0, y_tan0, "--", color="#1f77b4", lw=2, label=f"Дотична в $x_0 = {x0:.4f}$")
plt.scatter([x0], [d1_0], color="#1f77b4", s=60, zorder=5)
plt.scatter([x1], [0], color="black", s=40, zorder=5)
plt.vlines(x0, 0, d1_0, linestyles="dotted", color="gray")
plt.annotate(f"$x_0={x0:.4f}$", (x0, 0.2), ha="center", fontsize=10)

x_tan1 = np.linspace(-1.30, -1.285, 100)
y_tan1 = d1_1 + d2_1 * (x_tan1 - x1)
plt.plot(x_tan1, y_tan1, "--", color="#ff7f0e", lw=2, label=f"Дотична в $x_1 = {x1:.4f}$")
plt.scatter([x1], [d1_1], color="#ff7f0e", s=60, zorder=5)
plt.scatter([x2], [0], color="black", s=40, zorder=5)
plt.vlines(x1, 0, d1_1, linestyles="dotted", color="gray")
plt.annotate(f"$x_1={x1:.4f}$", (x1, 0.15), ha="right", fontsize=10)

plt.scatter([x2], [0], color="#d62728", s=80, zorder=6, label=f"Корінь $x^* = {x2:.6f}$")

plt.title("Метод Ньютона (Варіант 7)", fontsize=14)
plt.xlabel("$x$", fontsize=12)
plt.ylabel(r"$\varphi(x)$", fontsize=12)
plt.grid(True, alpha=0.3)
plt.legend(fontsize=10)
plt.tight_layout()

plt.savefig("newton_tangents_var7.png", dpi=300)
plt.show()