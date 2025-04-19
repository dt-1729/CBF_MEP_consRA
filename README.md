# CBF-Based Optimization for Constrained Resource Allocation

This repository provides a modular implementation of several optimization approaches for solving **constrained resource allocation** and **clustering** problems. The primary focus is our proposed method that uses **Control Barrier Functions (CBFs)** to guide a continuous-time optimization process to KKT points while respecting all constraints.

We compare this approach with:
- Classical Deterministic Annealing (DA) [Kenneth Rose, 1998](https://doi.org/10.1109/5.726788)
- Penalty-based DA for soft constraint handling [IEEE ICC, 2022](https://ieeexplore.ieee.org/abstract/document/10093253)
- Safe Gradient Flow (SGF) using barrier-guided dynamics [TAC, 2022](https://doi.org/10.1109/TAC.2022.3200517)
- Standard SLSQP optimization from [SciPy](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html)

---

## 🔧 Features

- CBF-based control-theoretic optimization with constraint satisfaction
- Unconstrained and penalty-based deterministic annealing (DA)
- Gradient-based safety-aware flow (SGF)
- Built-in SciPy SLSQP baseline
- Visualizations including trajectories, allocation costs, and capacity bounds
- Utilities for reshaping, projection, and synthetic problem generation

---

## 📁 Repository Structure

```
📦your-repo/
├── class_flp.py               # Core FLP class: all methods for clustering, DA, SGF, CBF
├── utils.py                   # Helper functions: reshaping, projection, test case generation
├── clustering_comparison.ipynb         # Solves one example using all methods, compares results
├── clustering_compare_annealing.ipynb  # Compares methods at a fixed beta (annealing stage)
├── README.md                  # This file
```

---

## 🚀 Getting Started

### 1. Install dependencies

```bash
pip install numpy scipy matplotlib cvxpy
```

### 2. Run an example

Open the Jupyter notebooks:

- `clustering_comparison.ipynb` — End-to-end example solving a problem with all methods
- `clustering_compare_annealing.ipynb` — Compares performance at a fixed beta value

---

## 🧠 Methods Implemented

| Method | Description |
|--------|-------------|
| **CBF-based (proposed)** | Designs a control system over the optimization variable that guarantees convergence to a KKT point while satisfying constraints throughout. [[Paper](https://arxiv.org/abs/2504.01378)] |
| **Penalty DA** | Introduces constraint penalties into deterministic annealing. Based on: _"Inequality Constraints in Facility Location and Related Problems," IEEE ICC, 2022_. [[Paper](https://ieeexplore.ieee.org/abstract/document/10093253)] |
| **Classical DA** | Unconstrained deterministic annealing for clustering. Based on: _Kenneth Rose, "Deterministic Annealing for Clustering...," Proc. IEEE, 1998_. [[Paper](https://doi.org/10.1109/5.726788)] |
| **SGF (Safe Gradient Flow)** | Gradient-based control using barrier functions to enforce constraints. Based on: _Allibhoy & Cortés, "Control-Barrier-Function-Based Design of Gradient Flows for Constrained Nonlinear Programming," TAC, 2022_. [[Paper](https://doi.org/10.1109/TAC.2022.3200517)] |
| **SLSQP (SciPy)** | Standard constrained optimization using `scipy.optimize.minimize(method='SLSQP')`. [[Docs](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html)] |

---

## 📝 Citation

If you use this codebase in your work, please cite:

> Alisina Bayati and Dhananjay Tiwary. *"A Control Barrier Function Approach to Constrained Resource Allocation Problems."* arXiv preprint, 2024. [[arXiv:2504.01378]](https://arxiv.org/abs/2504.01378)

---

## 👨‍💻 Authors

- **Dhananjay Tiwary** ([dt-1729](https://github.com/dt-1729)) — main contributor and core developer
- **Alisina Bayati** ([alisina75](https://github.com/alisina75)) — research lead and co-author

---

## 📜 License

MIT License  
© 2025 Dhananjay Tiwary (dt-1729) and Alisina Bayati (alisina75)



