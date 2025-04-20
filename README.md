# A Control Barrier Function Approach to Constrained Resource Allocation Problems in a Maximum Entropy Principle Framework

This repository implements the **control-theoretic optimization framework** proposed in [our paper](https://arxiv.org/abs/2504.01378), targeting a class of **constrained resource allocation problems** modeled as generalized **facility location problems (FLPs)**. The method uses **Control Barrier Functions (CBFs)** and a **Control Lyapunov Function (CLF)** to ensure constraint satisfaction and convergence to KKT points.

We also include the **classical Deterministic Annealing (DA)** [Rose, 1998](https://doi.org/10.1109/5.726788) for the unconstrained case, and compare our CBF-based approach with:

- **Penalty-based DA** [IEEE ICC, 2022](https://ieeexplore.ieee.org/document/10093253)
- **Safe Gradient Flow (SGF)** [TAC, 2022](https://ieeexplore.ieee.org/document/10224270)
- **SciPy’s SLSQP solver** [SciPy Docs](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html)

For full methodology and problem formulation, please refer to the [paper](https://arxiv.org/abs/2504.01378).

---

## 📁 Repository Structure

```
📦CBF_MEP_based_constrained_Resource_Allocation/
├── class_flp.py               # Core FLP class: all methods for MEP-based FLP, including classic DA, CBF (our proposed approach), SGF, and penalty-based DA
├── utils.py                   # Helper functions: reshaping, projection, test case generation
├── clustering_compare_annealing.ipynb  # Solves one example completely using all methods, compares results
├── clustering_comparison.ipynb         # Compares methods at a fixed beta (annealing stage)
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

- `clustering_compare_annealing.ipynb` — End-to-end example solving a problem with all methods
- `clustering_comparison.ipynb` — Compares performance at a fixed beta value

---

## 🧠 Methods Implemented & Visualized

| Method                       | Description                                                                                                                                                                                                                                                           |
| ---------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **CBF-based (proposed)**     | Designs a control system over the optimization variable that guarantees convergence to a KKT point while satisfying constraints throughout. [[Paper](https://arxiv.org/abs/2504.01378)]                                                                               |
| **Penalty DA**               | Introduces constraint penalties into deterministic annealing. Based on: *"Inequality Constraints in Facility Location and Related Problems," IEEE ICC, 2022*. [[Paper](https://ieeexplore.ieee.org/abstract/document/10093253)]                                       |
| **Classical DA**             | Unconstrained deterministic annealing for clustering. Based on: *Kenneth Rose, "Deterministic Annealing for Clustering...," Proc. IEEE, 1998*. [[Paper](https://doi.org/10.1109/5.726788)]                                                                            |
| **SGF (Safe Gradient Flow)** | Gradient-based control using barrier functions to enforce constraints. Based on: *Allibhoy & Cortés, "Control-Barrier-Function-Based Design of Gradient Flows for Constrained Nonlinear Programming," TAC, 2022*. [[Paper](https://doi.org/10.1109/TAC.2022.3200517)] |
| **SLSQP (SciPy)**            | Standard constrained optimization using `scipy.optimize.minimize(method='SLSQP')`. [[Docs](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html)]    

---

> 🎥 **Solution evolution during annealing across different methods:**

![Solution Animation](files/solution_comparison%20(1).gif)

---

## 📝 Citation

If you use this codebase in your work, please cite:

> Alisina Bayati and Dhananjay Tiwary. *"A Control Barrier Function Approach to Constrained Resource Allocation Problems in a Maximum Entropy Principle Framework."* arXiv preprint, 2024. [[arXiv:2504.01378]](https://arxiv.org/abs/2504.01378)

---

## 👨‍💻 Authors

- **Dhananjay Tiwary** ([dt-1729](https://github.com/dt-1729)) — core developer and implementation lead  
- **Alisina Bayati** ([alisina75](https://github.com/alisina75)) — led theoretical development and contributed to implementation

---

## 📜 License

MIT License  
© 2025 Dhananjay Tiwary (dt-1729) and Alisina Bayati (alisina75)














