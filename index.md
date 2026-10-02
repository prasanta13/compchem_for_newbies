---
layout: default
title: CompChem for Newbies
---

**Status:** this site is a work in progress. The introduction below, the Basis Sets page and the HF chapter are written. The MP2, CCSD and DFT chapter pages are placeholders with a one-line description and links; their content has not been written yet.

# Foreword

Often I think about how could I have learned this subject better. I believe there may be numerous other students who also think the same. As a kinesthetic learner, I often wonder about the equation(s), how do they look and how they shaped the modern **Quantum Chemistry** from **Quantum Mechanics**. This is an ongoing project which I want to continue as long as I can. This is also a testament of my learning. As I learn, I write more to learn more. 

At first, I will try to cover the **Hartree-Fock** approximation. It is imperative that all of us, students of computational chemistry to grasp this particular method. Then I will try to cover the basis-sets, different integrals and implement them with simple python codes.

### Quantum Mechanics

The concept of finding approximate solution to the **Schrödinger equation** is central: Here, it also means that we try to find the solution of non-relativistic time-independent **Schrödinger equation**.

$$
\hat{H} \Psi = E \Psi
$$

where:  
- $\hat{H}$  is the Hamiltonian operator for a system of nuclei and electrons having position vector $\mathbf{R}_A$ and $\mathbf{r}_i$ respectively.
- $\Psi$ is the wavefunction  
- $E$ is the energy

The electronic $\hat{H}$ for N electrons and M nuclei, in the Born-Oppenheimer approximation and in atomic units, is

$$
\hat{H}_{\text{elec}} = -\sum_{i=1}^{N} \frac{1}{2}\nabla^{2}_{i} - \sum_{i=1}^{N}\sum_{A=1}^{M} \frac{Z_A}{r_{iA}} + \sum_{i=1}^{N}\sum_{j>i}^{N} \frac{1}{r_{ij}}
$$

### Where to go next

The chapters build on this introduction. Start with [Basis sets](basis_set.md), then [Hartree-Fock](hf.md).

