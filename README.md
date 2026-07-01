# Numerical Solvers for Second-Order ODE Boundary Value Problems (BVP)

This repository contains a professional numerical analysis tool written in C for solving second-order ordinary differential equations (ODEs) with boundary conditions (Boundary Value Problems - BVP). The project implements multiple mathematical approaches to solve equations in the general form:
$$u'' + p(x)u' + q(x)u = f(x), \quad x \in (a,b)$$
with Dirichlet boundary conditions: $u(a) = u_a, \ e.g., \ u(b) = u_b$.

*Note: The source code, analytical comments, and mathematical functions inside this repository are written in Polish.*

---

## 🧠 Mathematical Methods Implemented

The project demonstrates two distinct paradigms for solving Boundary Value Problems:

1. **Finite Difference Method (Metoda Różnicowa):** * Discretizes the continuous differential equation into a system of linear equations.
   * Leverages an optimized **Thomas Algorithm (Algorytm Trójdiagonalny Thomasa)** for $O(n)$ time complexity matrix inversion, avoiding heavy Gaussian elimination.
2. **Shooting Method (Metoda Strzałów):**
   * Transforms the BVP into an Initial Value Problem (IVP).
   * Implements a custom **Root-Finding Bisection Algorithm (Bisekcja)** paired with a dynamic boundary-bracketing routine (`uchwyc_konce`) to find the correct initial slope that satisfies the target boundary condition $u(b) = u_b$.

---

## 📊 Visualization

The solver is integrated with the **DISLIN** scientific plotting library. Upon computing the numerical grid, the program automatically generates structural visualization charts showcasing the convergence behavior, axis grids, and custom background themes directly inside a graphical user interface window.

---

## 📁 Repository Structure

* `ode2_dislin.c` — Main entry point containing the configuration of the differential functions $p(x), q(x), f(x)$ and the execution pipeline.
* `metoda_strzalow.c` — Core algorithmic backend executing both the Finite Difference grid approximations and the Shooting/Bisection loops.
* `ode_dislin.dev` — Project configuration file for the **Embarcadero Dev-C++** Integrated Development Environment.
* `ludcmp.c` / `lubksb.c` — Supporting LU Decomposition routines for advanced matrix operations.

---

## 🛠️ Environment & Prerequisites

* **IDE:** Embarcadero Dev-C++ (configured with Project Compiler Settings for DISLIN linking).
* **Compiler Flags:** The project requires active linking against the DISLIN engine and Windows GUI API. The following linker flags are embedded in the `.dev` project file:
  ```text
  -ldismg -lgdi32 -luser32 -lopengl32
