# CBF-Based Optimization for Constrained Resource Allocation

This repository implements a **control design approach to optimization** proposed in [our paper](https://arxiv.org/abs/2504.01378). The method is showcased for a broad class of **resource allocation problems that can be formulated as generalized facility location problems (FLPs)**, where the assignment costs and problem interpretations may vary depending on the application, but still adhere to capacity constraints, where control barrier functions (CBFs) are used to ensure feasibility while guiding the solution toward KKT points.

The repository includes an implementation of the **classical Deterministic Annealing (DA)** algorithm [Kenneth Rose, 1998](https://doi.org/10.1109/5.726788), which solves the unconstrained FLP based on the Maximum Entropy Principle (MEP).

In our **MEP-based formulation** of **constrained** FLP, the algorithm requires solving an internal **constrained optimization problem** at each iteration. We compare our control-theoretic solution of this subproblem against several alternative approaches:

- **Penalty-based DA**, which incorporates inequality constraints through soft penalties [IEEE ICC, 2022](https://ieeexplore.ieee.org/abstract/document/10093253)
- **Safe Gradient Flow (SGF)**, a method based on barrier-function-guided dynamics [TAC, 2022](https://doi.org/10.1109/TAC.2022.3200517)
- **Sequential Least Squares Programming (SLSQP)**, a general-purpose solver from [SciPy](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html)

---

## 🔧 Features

- CBF-based control-theoretic optimization with constraint satisfaction
- Unconstrained and penalty-based deterministic annealing (DA)
- Gradient-based safety-aware flow (SGF)
- Built-in SciPy SLSQP baseline
- Visualizations including trajectories, allocation costs, and capacity bounds
- Utilities for reshaping, projection, and synthetic problem generation

---

# CBF-Based Optimization for Constrained Resource Allocation

This repository implements a **control-theoretic optimization approach** proposed in [our paper](https://arxiv.org/abs/2504.01378), targeting a broad class of **resource allocation problems** formulated as generalized **facility location problems** with flexible assignment costs and application-specific interpretations. Our approach uses **Control Barrier Functions (CBFs)** to ensure feasibility while steering the system toward KKT-optimality.

We also include the **classical Deterministic Annealing (DA)** method [Kenneth Rose, 1998](https://doi.org/10.1109/5.726788) for solving the unconstrained version of the facility location problem using the **Maximum Entropy Principle (MEP)**.

In the **MEP-based formulation** of constrained problems, each iteration requires solving a **constrained optimization subproblem**. We compare our control-based solution to several alternatives:

- **Penalty-based DA** – handles constraints via soft penalties [IEEE ICC, 2022](https://ieeexplore.ieee.org/abstract/document/10093253)
- **Safe Gradient Flow (SGF)** – uses barrier-function-based dynamics to converge safely to KKT points [TAC, 2022](https://ieeexplore.ieee.org/document/10224270)
- **Sequential Least Squares Programming (SLSQP)** – classical gradient-based method from [SciPy](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html)

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

| Method                       | Description                                                                                                                                                                                                                                                           |
| ---------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **CBF-based (proposed)**     | Designs a control system over the optimization variable that guarantees convergence to a KKT point while satisfying constraints throughout. [[Paper](https://arxiv.org/abs/2504.01378)]                                                                               |
| **Penalty DA**               | Introduces constraint penalties into deterministic annealing. Based on: *"Inequality Constraints in Facility Location and Related Problems," IEEE ICC, 2022*. [[Paper](https://ieeexplore.ieee.org/abstract/document/10093253)]                                       |
| **Classical DA**             | Unconstrained deterministic annealing for clustering. Based on: *Kenneth Rose, "Deterministic Annealing for Clustering...," Proc. IEEE, 1998*. [[Paper](https://doi.org/10.1109/5.726788)]                                                                            |
| **SGF (Safe Gradient Flow)** | Gradient-based control using barrier functions to enforce constraints. Based on: *Allibhoy & Cortés, "Control-Barrier-Function-Based Design of Gradient Flows for Constrained Nonlinear Programming," TAC, 2022*. [[Paper](https://doi.org/10.1109/TAC.2022.3200517)] |
| **SLSQP (SciPy)**            | Standard constrained optimization using `scipy.optimize.minimize(method='SLSQP')`. [[Docs](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html)]                                                                                        |

---

## 📝 Citation

If you use this codebase in your work, please cite:

> Alisina Bayati and Dhananjay Tiwary. *"A Control Barrier Function Approach to Constrained Resource Allocation Problems."* arXiv preprint, 2024. [[arXiv:2504.01378\]](https://arxiv.org/abs/2504.01378)

---

## 👨‍💻 Authors

- **Dhananjay Tiwary** ([dt-1729](https://github.com/dt-1729)) — main contributor and core developer
- **Alisina Bayati** ([alisina75](https://github.com/alisina75)) — research lead and co-author

---

## 📜 License

MIT License\
© 2025 Dhananjay Tiwary (dt-1729) and Alisina Bayati (alisina75)










