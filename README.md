# e8-parity-bisection

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22962066.svg)](https://doi.org/10.5281/zenodo.22962066)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Standard: C++20](https://img.shields.io/badge/C%2B%2B-20-blue.svg)](https://en.cppreference.com/w/cpp/20)
[![OEIS: A004010](https://img.shields.io/badge/OEIS-A004010-green.svg)](https://oeis.org/A004010)

Hardware-native, integer-only discrete lattice point enumeration of the 8-dimensional exceptional Lie algebra root lattice ($E_8$) collapsed to quasi-quadratic $\mathcal{O}(r^2)$ complexity via Diophantine scalar dilation and dimension-paired bisection.

This repository provides the reference C++20 implementation, benchmark harness, and formal verification suite for **Paper V** of the *Discrete Lattice Research Suite*.

---

## Theoretical Foundation

The classical 8-dimensional Gosset lattice $E_8$ (the densest known packing in $\mathbb{R}^8$ with kissing number 240) is standardly defined as the union of the root lattice $D_8$ and a shifted half-integer coset:

```math
E_8 = D_8 \cup \left( D_8 + \left(\tfrac{1}{2}\right)^8 \right) = \left\{ x \in \mathbb{Z}^8 \cup \left(\mathbb{Z} + \tfrac{1}{2}\right)^8 \;\middle|\; \sum_{i=1}^8 x_i \equiv 0 \pmod 2 \right\}
```

Evaluating high-dimensional volumes on this continuous representation traditionally requires floating-point operations, transcendental modular forms, or exponential $\mathcal{O}(r^8)$ coordinate walks.

### 1. The $2E_8$ Diophantine Dilation
To eliminate fractional half-integer registers and restore integer-native ALU execution, we apply a uniform global scalar dilation $y = 2x$. Under this mapping:
* The continuous half-integer coordinates map strictly into $\mathbb{Z}^8$.
* The squared Euclidean norm scales by a factor of 4:
  $$\Vert{}y\Vert{}^2 = 4\Vert{}x\Vert{}^2$$
* Because all minimal shell vectors in $E_8$ satisfy $\Vert{}x\Vert{}^2 = 2k$ ($k \in \mathbb{N}$), the dilated norm space evaluates as:
  $$\Vert{}y\Vert{}^2 = 8k \quad \Longrightarrow \quad \Vert{}y\Vert{}^2 \equiv 0 \pmod 8$$

### 2. All-Even and All-Odd Sublattice Partitioning
Under the $y = 2x$ dilation, the geometry partitions cleanly into two disjoint Diophantine parity manifolds governed by an invariant modulo-4 sum constraint:

1. **All-Even Parity Sublattice ($2D_8$):**
   $$y_i \equiv 0 \pmod 2 \quad \forall i, \quad \sum_{i=1}^8 y_i \equiv 0 \pmod 4$$
2. **All-Odd Parity Sublattice ($2 [ D_8 + (\tfrac{1}{2})^8 ]$):** $y_i \equiv 1 \pmod 2 \quad \forall i, \quad \sum_{i=1}^8 y_i \equiv 0 \pmod 4$
   $$y_i \equiv 1 \pmod 2 \quad \forall i, \quad \sum_{i=1}^8 y_i \equiv 0 \pmod 4$$

### 3. Orthogonal Dimension-Paired Bisection ($4\text{D} \times 4\text{D}$)
Because squared Euclidean distance is additively separable:

$$\Vert{}y\Vert{}^2 = \sum_{i=1}^4 y_i^2 + \sum_{j=5}^8 y_j^2 = S_A + S_B$$

The 8-dimensional space decouples into two orthogonal 4D submanifolds ($\mathbb{Z}^8 \cong \mathbb{Z}^4 \times \mathbb{Z}^4$). By pre-filtering 4D profiles into modulo-4 parity classes and executing sparse-index sweeps, the terminal volume integration collapses from $\mathcal{O}(r^8)$ to a single-pass 1D dot product in strict $\mathcal{O}(r^2)$.

---

## Empirical Verification & Benchmarks

The benchmark suite verifies exact lattice shell counts against the classical modular theta series divisor sum sieve:

$$\Theta_{E_8}(q) = 1 + 240 \sum_{k=1}^\infty \sigma_3(k) q^{2k} \quad (\text{OEIS A004010})$$

where $\sigma_3(k) = \sum_{d\vert{}k} d^3$.

### Benchmark Milestone ($n = 10{,}000$)
* **Target Radial Shell:** $n = 10{,}000$ (dilated boundary $\Vert{}y\Vert{}^2 \le 80{,}000$)
* **Accumulated Point Count:** `649,533,725,496,097,441`
* **Execution Time:** ~606 ms on consumer x86_64 hardware (single core)
* **Mathematical Accuracy:** 100% bit-exact parity match across all evaluated shells ($k = 0 \dots 10{,}000$) with zero floating-point emulation.

---

## Quickstart

### Prerequisites
* A 64-bit C++20 compliant compiler (`g++` $\ge 11$, `clang++` $\ge 13$, or MSVC $\ge 2019$).
* Target architecture with 64-bit integer registers (`uint64_t` / `__int128_t`).

### 1. Compilation
Clone the repository and compile with maximum optimization:

```bash
git clone [https://github.com/Andrew-J-Morris/e8-parity-bisection.git](https://github.com/Andrew-J-Morris/e8-parity-bisection.git)
cd e8-parity-bisection
g++ -O3 -std=c++20 -march=native E8_Parity_Bisection.cpp -o e8_benchmark
