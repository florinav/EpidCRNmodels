# %% imports and parameters
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

al = 0.8   # alpha
be = 1.2   # beta

def may_leonard_rhs(t, x):
    x1, x2, x3 = x
    return [
        x1*(1 - x1 - al*x2 - be*x3),
        x2*(1 - be*x1 - x2 - al*x3),
        x3*(1 - al*x1 - be*x2 - x3),
    ]

t_span = (0.0, 2000.0)
t_eval = np.linspace(t_span[0], t_span[1], 20001)

# several initial conditions in the simplex
inits = [
    [0.7, 0.2, 0.1],
    [0.6, 0.3, 0.1],
    [0.55, 0.25, 0.20],
    [0.50, 0.40, 0.10],
]

sims = [solve_ivp(may_leonard_rhs, t_span, x0, t_eval=t_eval,
                  rtol=1e-8, atol=1e-10)
        for x0 in inits]

# %% Figure 1: log10 min{x1,x2,x3} vs time for several inits
plt.figure(figsize=(8,4))
for k, sol in enumerate(sims, start=1):
    mins = np.min(sol.y, axis=0)
    plt.plot(sol.t, np.log10(mins), label=f"init {k}")
plt.xlabel("time")
plt.ylabel(r"$\log_{10}\min\{x_1,x_2,x_3\}$")
plt.legend()
plt.tight_layout()
# plt.savefig("may_leonard_min_log.png", dpi=300)

# %% Figure 2: simplex trajectory colored by log distance to boundary
sol = sims[0]                   # use first trajectory
x1, x2, x3 = sol.y

# barycentric coordinates -> 2D simplex
v1 = np.array([0.0, 0.0])
v2 = np.array([1.0, 0.0])
v3 = np.array([0.5, np.sqrt(3)/2])
pts = np.outer(x1, v1) + np.outer(x2, v2) + np.outer(x3, v3)

vals = np.log10(np.min(sol.y, axis=0))

plt.figure(figsize=(6,6))
sc = plt.scatter(pts[:,0], pts[:,1], c=vals, s=2)
cbar = plt.colorbar(sc)
cbar.set_label(r"$\log_{10}\min\{x_1,x_2,x_3\}$")
plt.axis("equal")
plt.axis("off")
plt.tight_layout()
# plt.savefig("may_leonard_simplex_colored.png", dpi=300)

# %% Figure 3: long-horizon time series for one trajectory
plt.figure(figsize=(8,4))
labels = ["x1", "x2", "x3"]
for i in range(3):
    plt.plot(sol.t, sol.y[i], label=labels[i])
plt.xlabel("time")
plt.ylabel("population")
plt.legend()
plt.tight_layout()
# plt.savefig("may_leonard_timeseries_long.png", dpi=300)
