# CBF-Based Optimization for Constrained Resource Allocation

This repository provides an implementation of a **control-theoretic optimization approach** proposed in [our paper](https://arxiv.org/abs/2504.01378). The method targets a broad class of **resource allocation problems** that can be formulated as generalized **facility location problems**—with flexible assignment cost structures and application-specific interpretations—subject to capacity constraints. Our approach uses **Control Barrier Functions (CBFs)** to ensure constraint satisfaction while steering the solution toward KKT-optimality.

We also include an implementation of the **classical Deterministic Annealing (DA)** algorithm [Kenneth Rose, 1998](https://doi.org/10.1109/5.726788), which solves the **unconstrained** version of the facility location problem based on the **Maximum Entropy Principle (MEP)**.

In our **MEP-based formulation** of constrained resource allocation, each iteration involves solving an internal **constrained optimization problem**. We compare our control-based solution of this subproblem with several alternatives:

- **Penalty-based DA** – adds soft penalties to handle constraints [IEEE ICC, 2022](https://ieeexplore.ieee.org/abstract/document/10093253)
- **Safe Gradient Flow (SGF)** – uses control-barrier-function-based dynamics to guide the system to KKT points [TAC, 2022](https://ieeexplore.ieee.org/document/10224270)
- **Sequential Least Squares Programming (SLSQP)** – a classical gradient-based solver from [SciPy](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html)

---

## 🔧 Features

- Generates and simulates **synthetic scenarios** for constrained facility location problems
- Uses a **Maximum Entropy Principle (MEP)** formulation to structure the optimization process
- Solves the internal constrained optimization subproblem using multiple methods:
  - ✅ Our proposed **CBF-based control-theoretic approach**
  - ✅ **Safe Gradient Flow (SGF)** with barrier-function-guided dynamics
  - ✅ **SciPy’s SLSQP** (Sequential Least Squares Programming)
  - ✅ **Penalty-based Deterministic Annealing (DA)** with soft constraints
- Includes classical **unconstrained DA** as a baseline for comparison
- Provides visualization tools for trajectories, cost trends, and capacity constraints
- Includes utility functions for reshaping variables, projecting onto the simplex, and generating test cases

---

## 📁 Repository Structure

```
📦CBF_MEP_based_constrained_Resource_Allocation/
├── class_flp.py               # Core FLP class: defines the MEP-based formulation, 
│                              # including unconstrained DA, penalty-based DA, SGF, and CBF methods
├── utils.py                   # Helper functions: reshaping, projection, and synthetic test case generation
├── clustering_comparison.ipynb         # Solves a full example using all methods; performs full annealing over a range of beta values
├── clustering_compare_annealing.ipynb  # Compares different constrained optimization approaches at a fixed beta value
├── README.md                  # Project overview and documentation
```


## 🚀 Getting Started

### 1. Install dependencies

```bash
pip install numpy scipy matplotlib cvxpy
```

### 2. Run an example

Open the Jupyter notebooks:

- `clustering_comparison.ipynb` — Runs a complete example of the constrained facility location problem using all methods (CBF, SGF, SLSQP, penalty-based DA), performs annealing over a range of β values, and compares their performance.
- `clustering_compare_annealing.ipynb` — Evaluates and benchmarks the constrained optimization methods at a **fixed β**, highlighting differences in convergence, cost, and feasibility.

---

## 🧠 Methods Implemented

| Method | Description |
|--------|-------------|
| **CBF-based (proposed)** | Designs a control system over the optimization variable that guarantees convergence to a KKT point while satisfying constraints throughout. [[Paper](https://arxiv.org/abs/2504.01378)] |
| **Penalty DA** | Introduces constraint penalties into deterministic annealing. Based on: _"Inequality Constraints in Facility Location and Related Problems," IEEE ICC, 2022_. [[Paper](https://ieeexplore.ieee.org/abstract/document/10093253)] |
| **Classical DA** | Unconstrained deterministic annealing for clustering. Based on: _Kenneth Rose, "Deterministic Annealing for Clustering...," Proc. IEEE, 1998_. [[Paper](https://doi.org/10.1109/5.726788)] |
| **SGF (Safe Gradient Flow)** | Gradient-based control using barrier functions to enforce constraints. Based on: _Allibhoy & Cortés, "Control-Barrier-Function-Based Design of Gradient Flows for Constrained Nonlinear Programming," TAC, 2022_. [[Paper](https://ieeexplore.ieee.org/document/10224270)] |
| **SLSQP (SciPy)** | Standard constrained optimization using `scipy.optimize.minimize(method='SLSQP')`. [[Docs](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html)] |

---

## 📝 Citation

If you use this codebase in your work, please cite:

> Alisina Bayati and Dhananjay Tiwary. *"A Control Barrier Function Approach to Constrained Resource Allocation Problems."* arXiv preprint, 2024. [[arXiv:2504.01378]](https://arxiv.org/abs/2504.01378)

---

## 👨‍💻 Authors

- **Dhananjay Tiwary** ([dt-1729](https://github.com/dt-1729)) — core developer and implementation lead  
- **Alisina Bayati** ([alisina75](https://github.com/alisina75)) — led theoretical development and contributed to implementation
- 
---

## 📜 License

MIT License  
© 2025 Dhananjay Tiwary (dt-1729) and Alisina Bayati (alisina75)








