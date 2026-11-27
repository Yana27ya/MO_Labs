import numpy as np
import matplotlib.pyplot as plt

def F(x1, x2):
    return 3*x1**2 + x2**2 + 2*x1*x2 + 16*x1 + 4*x2

def grad(x):
    return np.array([6*x[0] + 2*x[1] + 16, 2*x[0] + 2*x[1] + 4], dtype=float)

def line_search(x, d, delta0=0.01, eps=1e-6):
    phi = lambda h: F(x[0] + h*d[0], x[1] + h*d[1])
    h_prev, f_prev = 0.0, phi(0.0)
    step = delta0
    h_curr, f_curr = step, phi(step)
    if f_curr >= f_prev:
        while f_curr >= f_prev and step > 1e-12:
            step /= 2.0
            h_curr, f_curr = step, phi(step)
        a, b = 0.0, h_curr * 2.0
    else:
        while True:
            step *= 2.0
            h_next = h_curr + step
            f_next = phi(h_next)
            if f_next >= f_curr:
                a, b = h_prev, h_next
                break
            h_prev, f_prev, h_curr, f_curr = h_curr, f_curr, h_next, f_next
    tau = (np.sqrt(5.0) - 1.0) / 2.0
    x1, x2 = b - tau*(b - a), a + tau*(b - a)
    f1, f2 = phi(x1), phi(x2)
    while (b - a) > eps:
        if f1 < f2:
            b, x2, f2 = x2, x1, f1
            x1 = b - tau*(b - a); f1 = phi(x1)
        else:
            a, x1, f1 = x1, x2, f2
            x2 = a + tau*(b - a); f2 = phi(x2)
    return 0.5 * (a + b)

def traj_sphere_lr8(x0, h0=1.0, eps=0.01):
    x = np.array(x0, float); Fx = F(x[0], x[1]); h = h0
    pts = [x.copy()]
    for _ in range(500):
        g = grad(x)
        if np.linalg.norm(g) <= eps: break
        d = -g / np.linalg.norm(g)
        while True:
            xb = x + h*d; Fb = F(xb[0], xb[1])
            if Fb < Fx: break
            h /= 2.0
        x, Fx = xb, Fb; pts.append(x.copy())
    return np.array(pts)

def traj_fr(x0, eps=0.01, steepest=False):
    x = np.array(x0, float); g = grad(x); d = -g.copy()
    pts = [x.copy()]
    for k in range(500):
        if np.linalg.norm(g) <= eps: break
        h = line_search(x, d)
        x = x + h*d; pts.append(x.copy())
        g_new = grad(x)
        beta = 0.0 if (steepest or (k + 1) % 2 == 0) else np.dot(g_new, g_new) / np.dot(g, g)
        d = -g_new + beta*d; g = g_new
    return np.array(pts)

x0 = [3.0, -3.0]
t_sph = traj_sphere_lr8(x0)
t_sd  = traj_fr(x0, steepest=True)
t_fr  = traj_fr(x0, steepest=False)

X1, X2 = np.meshgrid(np.linspace(-6, 5, 400), np.linspace(-5, 5, 400))
Z = F(X1, X2)

plt.figure(figsize=(10, 7))
cs = plt.contour(X1, X2, Z, levels=[-20, -15, -7.7, 5, 25, 54, 85], colors="gray", linewidths=0.9, alpha=0.7)
plt.clabel(cs, inline=True, fontsize=8, fmt="F=%.0f")

plt.plot(t_sph[:,0], t_sph[:,1], "o-", color="#1f77b4", lw=1.8, ms=4, label=f"Сферична ЛР 8 ({len(t_sph)-1} ітер.)")
plt.plot(t_sd[:,0],  t_sd[:,1],  "s-", color="#9467bd", lw=2.0, ms=4, label=f"Антиградієнт + точний пошук ({len(t_sd)-1} ітер.)")
plt.plot(t_fr[:,0],  t_fr[:,1],  "^-", color="#2ca02c", lw=2.5, ms=7, label=f"Флетчер–Рівс ({len(t_fr)-1} ітер.)")

plt.scatter([3], [-3], color="black", s=80, zorder=6, label="Старт $x_{(0)} = (3; -3)$")
plt.scatter([-3], [1], color="red", s=90, zorder=6, label="Мінімум $x^* = (-3; 1)$")

plt.title("Порівняння траєкторій (Варіант 7, ЛР 9)", fontsize=13)
plt.xlabel("$x_1$", fontsize=12); plt.ylabel("$x_2$", fontsize=12)
plt.grid(True, alpha=0.3); plt.legend(fontsize=10); plt.tight_layout()
plt.savefig("lr9_trajectories_var7.png", dpi=300)
plt.show()