import numpy as np
import matplotlib.pyplot as plt

def df(x):
    return 7.0 * x**3 - 18.1125 * x**2 - 11.165 * x + 31.019625


x0 = -1.35
x1 = -1.05
x2 = -1.283969 # Результат 1-ї ітерації
x3 = -1.293485 # Результат 2-ї ітерації

y0 = df(x0)
y1 = df(x1)
y2 = df(x2)
y3 = df(x3)

x_vals = np.linspace(-1.36, -1.04, 400)
phi_vals = df(x_vals)

plt.figure(figsize=(10, 6))
plt.plot(x_vals, phi_vals, color="#2ca02c", lw=2.5, label=r"$\varphi(x) = f'(x)$")
plt.axhline(0, color="black", lw=1)

# Перша січна (через x0 та x1)
plt.plot([x0, x1], [y0, y1], "--", color="#1f77b4", lw=2, label=f"Січна 1 (через $x_0$, $x_1$)")
plt.scatter([x0, x1], [y0, y1], color="#1f77b4", s=60, zorder=5)
plt.scatter([x2], [0], color="black", s=40, zorder=5) # Перетин 1
plt.vlines(x2, 0, y2, linestyles="dotted", color="gray")

# Друга січна (через x0 та x2, бо ліва точка прилипла)
plt.plot([x0, x2], [y0, y2], "--", color="#ff7f0e", lw=2, label=f"Січна 2 (через $x_0$, $x_2$)")
plt.scatter([x2], [y2], color="#ff7f0e", s=60, zorder=5)
plt.scatter([x3], [0], color="black", s=40, zorder=5) # Перетин 2

# Фінальна точка
plt.scatter([x3], [0], color="#d62728", s=80, zorder=6, label=f"Наближення кореня $x^*$")

# Оформлення
plt.title("Метод січних (Варіант 7)", fontsize=14)
plt.xlabel("$x$", fontsize=12)
plt.ylabel(r"$\varphi(x) = f'(x)$", fontsize=12)
plt.grid(True, alpha=0.3)
plt.legend(fontsize=10)
plt.tight_layout()
plt.savefig("secants_var7.png", dpi=300)
plt.show()