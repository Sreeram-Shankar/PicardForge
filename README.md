# 📦 PicardForge
### *A fixed-step ODE & PDE time-integration library featuring Picard–Gauss–Seidel implicit solvers and classical explicit methods*

---

## ✨ Overview

**PicardForge** is a lightweight numerical solver library implementing a complete suite of **fixed-step ODE solvers**, designed especially for **semi-discretized PDEs** (heat equation, diffusion, conduction models, etc.).

Unlike traditional solver libraries that depend on expensive Newton iterations, **PicardForge uses Picard fixed-point iteration with Gauss–Seidel relaxation** to solve all implicit systems. This drastically reduces cost for diffusion-type problems while maintaining full A-stability and high accuracy for stiff systems.

The library includes:

- **Explicit Runge–Kutta (RK1–RK6)**
- **Adams–Bashforth multistep (AB2–AB5)**
- **Adams–Moulton multistep (AM2–AM5)**
- **Backward Differentiation Formulas (BDF1–BDF6)**
- **SDIRK (2nd–4th order)**
- **Fully implicit Gauss, Radau, and Lobatto collocation IRK (s = 1–5)**

All implicit families are solved using **Picard–Gauss–Seidel**, giving a robust, Jacobian-free integration strategy well-suited for large systems arising from PDE discretization.

---

## 🚀 Features

### ✔ Full suite of classic numerical integrators

| Family | Methods | Notes |
|-------|---------|-------|
| **Explicit RK** | RK1–RK6 | Standard explicit Butcher tables |
| **Adams–Bashforth** | AB2–AB5 | Explicit multistep |
| **Adams–Moulton** | AM2–AM5 | Implicit multistep (Picard solved) |
| **BDF** | BDF1–BDF6 | All implicit, solved with Picard–GS |
| **SDIRK** | SDIRK2–SDIRK4 | Diagonally implicit RK |
| **Gauss–Legendre IRK** | s = 1–5 | A-stable, symplectic |
| **Radau IIA IRK** | s = 1–5 | L-stable; great for stiff PDEs |
| **Lobatto IIIC IRK** | s = 2–5 | Symmetric, stiffly accurate |

### ✔ Picard–Gauss–Seidel nonlinear iteration
- No Jacobians needed  
- No Newton linear solves  
- Stage-by-stage relaxation  
- Ideal for diffusion-dominated PDEs  
- Stable and efficient for mildly nonlinear systems  

### ✔ Exact Butcher tableau generation
The `generation/` folder contains high-precision quadrature-based generators for:
- Gauss–Legendre IRK
- Radau IIA IRK
- Lobatto IIIC IRK

Adapted from **MethodForge** and simplified to return only numeric tables.

### ✔ PDE-ready design
Tailored for time integration of systems from:
- radial heat conduction  
- parabolic PDE discretization  
- finite differences / volumes  

### ✔ Minimal dependencies
Only:
- `numpy`
- `mpmath` (table generation only)

---

## 📁 Repository Structure

```
PicardForge/
│
├── solvers/
│   ├── rk.py              # RK1–RK6 explicit solvers
│   ├── ab.py              # AB2–AB5 explicit multistep
│   ├── am.py              # AM2–AM5 implicit multistep (Picard)
│   ├── bdf.py             # BDF1–BDF6 implicit multistep (Picard)
│   ├── irk.py             # Gauss/Radau/Lobatto collocation IRK (Picard)
│   ├── sdirk.py           # SDIRK2–SDIRK4 implicit RK
│   └── __init__.py        # Unified solver_map
│
├── generation/
│   ├── gauss_legendre.py
│   ├── radau.py
│   └── lobatto.py
│
└── README.md
```

---

## 🧠 How the Picard–Gauss–Seidel nonlinear iteration works

All implicit solvers need to solve the collocation equation:

\[
Y = y_n + h A F(Y)
\]

Rather than forming a full Jacobian matrix and solving a nonlinear system with Newton’s method, we use a **Picard fixed-point iteration** enhanced with **Gauss–Seidel relaxation**.

### **Picard iteration**
\[
Y^{(k+1)} = y_n + h A F\Big(Y^{(k)}\Big)
\]

### **Gauss–Seidel sweep**

```
for k in 1..sweeps:
    for i in 1..s:
        Y[i] = y + h * Σ_j A[i,j] * f(t + c[j]*h, Y[j])
```

This approach is:

- matrix-free  
- stable for diffusion problems  
- much faster than Newton  
- convergent for nearly-linear PDE operators  

For your use case (diffusion/conduction planetary models), Picard–GS is ideal.

---

## 🔧 Usage Examples

### **Simple ODE example**

```python
import numpy as np
from solvers import solver_map

def f(t, y):
    return -0.01 * y

y0 = np.array([1.0])
t, Y = solver_map["R5"](f, (0, 10.0), y0, h=0.1)

print("Final value:", Y[-1])
```

### **Semi-discretized PDE example**

```python
from solvers import solver_map
from my_pde_model import rhs  # your PDE's right-hand side

T0 = initial_temperature_profile()
t, T = solver_map["BDF4"](rhs, (0, 1e7), T0, h=1e5)
```

---

## 📈 Stability Notes

- **Explicit solvers (RK, AB)**  
  Not suited for stiff PDEs — limited by von Neumann stability constraints.

- **Implicit multistep (AM, BDF)**  
  A-stable (AM2) or strongly A-stable-ish (higher orders). Good for stiff systems.

- **Collocation IRK (Gauss, Radau, Lobatto)**  
  - Gauss → A-stable, symplectic  
  - Radau → L-stable (best for stiff PDEs)  
  - Lobatto → symmetric, good for reversible problems  

- **Picard–GS iteration**  
  Converges extremely well for heat/diffusion equations.

---

## 🛠 Requirements

- Python **3.8+**
- NumPy
- mpmath (for generating IRK tables)

## 📜 License

Released under the **MIT License**.

---

## 🤝 Contributing

Contributions are welcome!  
Ideas include:

- implementing adaptive IRK  
- embedded error estimators  
- Jacobian-free Newton–Krylov  
- multigrid-enhanced Picard  
- GPU-accelerated f(t, y) evaluation  

Enjoy exploring implicit and explicit time integration with **PicardForge**!
