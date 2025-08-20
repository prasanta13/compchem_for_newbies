---
title: "Introduction to Computational Chemistry"
---

# Introduction to Computational Chemistry

Computational chemistry allows us to study molecules and reactions using computer simulations.  
We can calculate energies, molecular geometries, vibrational frequencies, and more.

## Key Concepts

### Quantum Mechanics

The **Schrödinger equation** is central:

$$
\hat{H} \Psi = E \Psi
$$

where:  
- \(\hat{H}\) is the Hamiltonian operator  
- \(\Psi\) is the wavefunction  
- \(E\) is the energy

### Basis Sets

Basis sets approximate the molecular orbitals as linear combinations of known functions:

\[
\phi_i(\mathbf{r}) = \sum_\mu C_{\mu i} \chi_\mu(\mathbf{r})
\]

- \(\chi_\mu(\mathbf{r})\) are atomic orbitals  
- \(C_{\mu i}\) are coefficients

### Hartree-Fock Approximation

The **Hartree-Fock equation** in matrix form:

\[
\mathbf{F}\mathbf{C} = \mathbf{S}\mathbf{C}\mathbf{\varepsilon}
\]

- \(\mathbf{F}\) is the Fock matrix  
- \(\mathbf{S}\) is the overlap matrix  
- \(\mathbf{C}\) contains MO coefficients  
- \(\mathbf{\varepsilon}\) contains orbital energies

### Inline Math Example

Energy of hydrogen atom: \(E_n = -\frac{13.6\,\text{eV}}{n^2}\)

---

You can add figures, tables, and code blocks as usual in Markdown. For example:

```python
from pyscf import gto, scf

mol = gto.M(atom='H 0 0 0; H 0 0 0.74', basis='sto-3g')
mf = scf.RHF(mol).run()
print(mf.energy_tot())

