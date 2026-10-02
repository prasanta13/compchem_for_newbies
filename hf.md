---
layout: default
title: Hartree-Fock (HF)
---

# Hartree–Fock (HF)

The Hartree–Fock method is the foundation of most ab initio quantum chemistry.  
It approximates the **N-electron wavefunction** as a single Slater determinant.

In this chapter I want to explain what that sentence means, starting from zero. I assume you know some algebra, what a derivative and an integral are, and what a matrix is. I do not assume any physical chemistry. If you have not read the page on [basis sets](basis_set.md) yet, read it first, because every practical Hartree-Fock calculation is built on a basis set.

The plan is:

1. the problem we are trying to solve;
2. the mean-field idea behind Hartree-Fock;
3. how to build a wavefunction out of orbitals;
4. the variational principle;
5. deriving the Hartree-Fock energy;
6. deriving the Hartree-Fock equations;
7. bringing in the basis set (the Roothaan-Hall equations);
8. the integrals, one at a time;
9. the Fock matrix;
10. the self-consistent field procedure;
11. doing all of this in PySCF, including writing our own Hartree-Fock program;
12. what Hartree-Fock gets wrong.

## 1. The problem we are trying to solve

On the [home page](index.md) I wrote down the Schrödinger equation, $\hat{H}\Psi = E\Psi$, and the electronic Hamiltonian:

$$
\hat{H}_{\text{elec}} = -\sum_{i=1}^{N} \frac{1}{2}\nabla^{2}_{i}
 - \sum_{i=1}^{N}\sum_{A=1}^{M} \frac{Z_A}{r_{iA}}
 + \sum_{i=1}^{N}\sum_{j>i}^{N} \frac{1}{r_{ij}}
$$

Let me read it term by term. There are $N$ electrons, labelled $i$ and $j$, and $M$ nuclei, labelled $A$.

1. $-\frac{1}{2}\nabla^{2}_{i}$ is the **kinetic energy** of electron $i$. The symbol $\nabla^2$ ("del squared") measures how curved a function is; the more sharply the wavefunction bends, the higher the kinetic energy.
2. $-Z_A / r_{iA}$ is the **attraction** between electron $i$ and nucleus $A$. Here $Z_A$ is the charge of the nucleus (1 for hydrogen, 8 for oxygen) and $r_{iA}$ is the distance between them. It is negative because opposite charges attract.
3. $1/r_{ij}$ is the **repulsion** between electrons $i$ and $j$, where $r_{ij}$ is the distance between them. The condition $j > i$ makes sure each pair is counted once.

Two remarks about this equation.

**Atomic units.** There are no constants such as the electron mass or the charge of the electron in the equation, because we measure everything in atomic units, where those constants are all equal to 1. Lengths are measured in bohr (1 bohr = 0.529 Å) and energies in Hartree (1 Hartree = 27.211 eV = 2625.5 kJ/mol).

**The Born-Oppenheimer approximation.** A nucleus is at least about 1800 times heavier than an electron, so on the time scale of electronic motion the nuclei look frozen. We therefore keep the nuclei fixed, solve for the electrons, and add the repulsion between the nuclei at the end as a constant:

$$
E_{\text{nuc}} = \sum_{A=1}^{M}\sum_{B>A}^{M} \frac{Z_A Z_B}{R_{AB}}
$$

**Why is this hard?** The first two terms each involve one electron at a time. If only those existed, the problem would fall apart into $N$ separate one-electron problems, each of which is easy. The third term spoils this: it couples every pair of electrons, so where electron 1 is likely to be depends on where electron 2 is, and so on. The wavefunction $\Psi$ is a function of the positions of all the electrons at once, which is $3N$ coordinates. For water, with 10 electrons, that is a function of 30 variables. Storing it on a grid with only 10 points along each of those 30 directions would need $10^{30}$ numbers. We need an approximation.

## 2. The mean-field idea

Think of walking through a crowded railway station. You do not track every single person around you; you react to the crowd as a whole, avoiding the busy areas. Hartree's idea was the same: let each electron move in the **average** electric field created by all the other electrons, instead of reacting to each of them individually.

This turns one impossible $N$-electron problem into $N$ one-electron problems. Each electron gets its own **orbital**, a function of just three coordinates.

There is a catch. The average field felt by one electron depends on the orbitals of all the others, and those are exactly what we are trying to find. So we have to guess the orbitals, compute the field, find new orbitals in that field, compute a new field, and repeat until nothing changes. When the orbitals produce a field that reproduces the same orbitals, the solution is **self-consistent**. This is why Hartree-Fock is also called the **self-consistent field** (SCF) method.

What do we lose? Real electrons avoid each other instantly: if one electron is on the left of a bond, the other one is a little more likely to be on the right. A mean field cannot describe this, because each electron only sees an average. The energy we miss is called the **correlation energy**. Hartree-Fock typically recovers more than 99% of the total energy of a molecule, but the missing part is often about the same size as the energies chemists care about, such as reaction energies. Recovering it is the job of the later chapters: [MP2](mp2.md), [CCSD](ccsd.md) and [DFT](dft.md).

## 3. Building a wavefunction out of orbitals

### Spin

Electrons have a property called spin, which can be "up" or "down". We write the two possibilities as spin functions $\alpha(\omega)$ and $\beta(\omega)$, where $\omega$ is a spin coordinate. A **spin orbital** is a spatial orbital $\phi(\mathbf{r})$ multiplied by a spin function:

$$
\psi(\mathbf{x}) = \phi(\mathbf{r})\,\alpha(\omega) \quad \text{or} \quad \psi(\mathbf{x}) = \phi(\mathbf{r})\,\beta(\omega)
$$

where $\mathbf{x} = (\mathbf{r}, \omega)$ collects the position and the spin of an electron. The spin functions are orthonormal: integrating over $\omega$ gives $\langle \alpha \vert \alpha \rangle = \langle \beta \vert \beta \rangle = 1$ and $\langle \alpha \vert \beta \rangle = 0$.

A note on notation, which I keep throughout: $\Psi$ (capital) is the full wavefunction, $\psi_a$ (small) is a spin orbital, $\phi_i$ is a spatial orbital and $\chi_\mu$ is a basis function.

### The simplest guess, and why it fails

The simplest way to combine orbitals is to multiply them: electron 1 in spin orbital $\psi_a$, electron 2 in $\psi_b$,

$$
\Psi(\mathbf{x}_1, \mathbf{x}_2) = \psi_a(\mathbf{x}_1)\, \psi_b(\mathbf{x}_2)
$$

This is called a Hartree product, and it breaks a basic rule of nature. Electrons are identical, so swapping two of them must not change anything we can measure, and quantum mechanics requires more: the wavefunction must change sign when two electrons are swapped. This is the **antisymmetry principle**, and the Pauli exclusion principle you may know from general chemistry follows from it. The Hartree product does not change sign when you swap $\mathbf{x}_1$ and $\mathbf{x}_2$.

### The Slater determinant

The fix is to subtract the swapped version. For two electrons:

$$
\Psi(\mathbf{x}_1, \mathbf{x}_2)
= \frac{1}{\sqrt{2}}
\begin{vmatrix}
\psi_a(\mathbf{x}_1) & \psi_b(\mathbf{x}_1) \\
\psi_a(\mathbf{x}_2) & \psi_b(\mathbf{x}_2)
\end{vmatrix}
= \frac{1}{\sqrt{2}} \Big[ \psi_a(\mathbf{x}_1)\psi_b(\mathbf{x}_2) - \psi_b(\mathbf{x}_1)\psi_a(\mathbf{x}_2) \Big]
$$

Check the two properties for yourself:

- swap $\mathbf{x}_1$ and $\mathbf{x}_2$ and the bracket changes sign, as required;
- put both electrons in the same spin orbital, $\psi_a = \psi_b$, and the bracket is zero. Two electrons cannot occupy the same spin orbital: that is the Pauli exclusion principle.

The factor $1/\sqrt{2}$ makes the wavefunction normalised. For $N$ electrons we write an $N \times N$ determinant with prefactor $1/\sqrt{N!}$, called a **Slater determinant**.

**The Hartree-Fock approximation** is: the wavefunction is a single Slater determinant, and we look for the best possible spin orbitals to put in it.

## 4. The variational principle

What does "best" mean? The answer comes from the **variational principle**: for any trial wavefunction $\Psi$, the energy it predicts is never lower than the true ground-state energy $E_0$,

$$
E[\Psi] = \frac{\langle \Psi \vert \hat{H} \vert \Psi \rangle}{\langle \Psi \vert \Psi \rangle} \ge E_0
$$

The angle brackets are shorthand for integrals over the coordinates of all electrons:

$$
\langle \Psi \vert \hat{H} \vert \Psi \rangle = \int \Psi^{\ast}\, \hat{H}\, \Psi \; d\mathbf{x}_1 \cdots d\mathbf{x}_N,
\qquad
\langle \Psi \vert \Psi \rangle = \int \Psi^{\ast}\, \Psi \; d\mathbf{x}_1 \cdots d\mathbf{x}_N
$$

(The star means complex conjugate; for the real functions we use from Section 5 onwards it does nothing.) So the lower the energy, the closer we are to the truth. The best Slater determinant is the one with the **lowest energy**. To find it we need two things: a formula for the energy of a Slater determinant, and a way to minimise it.

## 5. Deriving the Hartree-Fock energy

First, split the Hamiltonian into one-electron and two-electron parts:

$$
\hat{H}_{\text{elec}} = \sum_{i=1}^{N} \hat{h}(i) + \sum_{i=1}^{N}\sum_{j>i}^{N} \frac{1}{r_{ij}},
\qquad
\hat{h}(i) = -\frac{1}{2}\nabla^{2}_{i} - \sum_{A=1}^{M} \frac{Z_A}{r_{iA}}
$$

$\hat{h}$ is called the **core Hamiltonian**: the energy of one electron moving among the bare nuclei, as if the other electrons were not there.

### Two electrons, step by step

Take the two-electron determinant from Section 3, with orthonormal spin orbitals: $\langle \psi_a \vert \psi_a \rangle = \langle \psi_b \vert \psi_b \rangle = 1$ and $\langle \psi_a \vert \psi_b \rangle = 0$. Write $\vert pq \rangle$ for the product $\psi_p(\mathbf{x}_1)\psi_q(\mathbf{x}_2)$, so that $\Psi = \frac{1}{\sqrt{2}}\left( \vert ab \rangle - \vert ba \rangle \right)$.

**The one-electron part.** $\hat{h}(1)$ only acts on electron 1, so for electron 2 the integral is just an overlap:

$$
\langle pq \vert \hat{h}(1) \vert rs \rangle = \langle \psi_p \vert \hat{h} \vert \psi_r \rangle \, \langle \psi_q \vert \psi_s \rangle
$$

Expanding the determinant gives four terms:

$$
\begin{aligned}
\langle \Psi \vert \hat{h}(1) \vert \Psi \rangle
&= \tfrac{1}{2} \Big[ \langle ab \vert \hat{h}(1) \vert ab \rangle - \langle ab \vert \hat{h}(1) \vert ba \rangle - \langle ba \vert \hat{h}(1) \vert ab \rangle + \langle ba \vert \hat{h}(1) \vert ba \rangle \Big] \\
&= \tfrac{1}{2} \Big[ h_{aa}\cdot 1 - h_{ab}\cdot 0 - h_{ba}\cdot 0 + h_{bb}\cdot 1 \Big]
= \tfrac{1}{2} \left( h_{aa} + h_{bb} \right)
\end{aligned}
$$

where $h_{aa} = \langle \psi_a \vert \hat{h} \vert \psi_a \rangle$. The cross terms vanish because the orbitals are orthogonal. Electron 2 gives the same result, so the total one-electron energy is $h_{aa} + h_{bb}$: simply the core energy of each occupied spin orbital.

**The two-electron part.** Define

$$
\langle pq \vert rs \rangle = \iint \psi_p^{\ast}(\mathbf{x}_1)\, \psi_q^{\ast}(\mathbf{x}_2)\, \frac{1}{r_{12}}\, \psi_r(\mathbf{x}_1)\, \psi_s(\mathbf{x}_2)\; d\mathbf{x}_1\, d\mathbf{x}_2
$$

Then

$$
\langle \Psi \vert \tfrac{1}{r_{12}} \vert \Psi \rangle
= \tfrac{1}{2} \Big[ \langle ab \vert ab \rangle - \langle ab \vert ba \rangle - \langle ba \vert ab \rangle + \langle ba \vert ba \rangle \Big]
$$

Renaming the integration variables $\mathbf{x}_1 \leftrightarrow \mathbf{x}_2$ does not change an integral, and it shows that $\langle ba \vert ba \rangle = \langle ab \vert ab \rangle$ and $\langle ba \vert ab \rangle = \langle ab \vert ba \rangle$. So

$$
\langle \Psi \vert \tfrac{1}{r_{12}} \vert \Psi \rangle = \langle ab \vert ab \rangle - \langle ab \vert ba \rangle = J_{ab} - K_{ab}
$$

and the total energy of the two-electron determinant is

$$
E = h_{aa} + h_{bb} + J_{ab} - K_{ab}
$$

The two new quantities are worth understanding properly.

- The **Coulomb integral** $$J_{ab} = \iint \lvert \psi_a(\mathbf{x}_1) \rvert^2 \frac{1}{r_{12}} \lvert \psi_b(\mathbf{x}_2) \rvert^2 \, d\mathbf{x}_1 d\mathbf{x}_2$$ is ordinary electrostatics: the repulsion between the charge cloud of electron 1 and the charge cloud of electron 2.
- The **exchange integral** $K_{ab} = \langle ab \vert ba \rangle$ has no classical counterpart. It appears only because the wavefunction is antisymmetric. Integrating over spin shows that it is zero unless $\psi_a$ and $\psi_b$ have the same spin (because $\langle \alpha \vert \beta \rangle = 0$). It enters with a minus sign and lowers the energy: antisymmetry keeps electrons of the same spin away from each other, so they repel each other less.

### Any number of electrons

The same bookkeeping works for $N$ electrons (the general version is known as the Slater-Condon rules). Every electron contributes its core energy, and every pair of electrons contributes a Coulomb term minus an exchange term:

$$
E_{\text{HF}} = \sum_{a=1}^{N} h_{aa} + \frac{1}{2} \sum_{a=1}^{N} \sum_{b=1}^{N} \left( J_{ab} - K_{ab} \right)
$$

The $\frac{1}{2}$ is there because the double sum counts each pair twice. The double sum also includes $a = b$, which would be an electron repelling itself; but $J_{aa} = K_{aa}$, so these terms cancel exactly. Hartree-Fock has no self-interaction, which is a nice property that some approximate DFT methods lack.

### Closed shells and spatial orbitals

Most stable molecules have an even number of electrons, with each spatial orbital $\phi_i$ holding one spin-up and one spin-down electron. This is called **restricted Hartree-Fock** (RHF), and from now on I only treat this case. There are $N/2$ occupied spatial orbitals, and the orbitals can be chosen to be real functions, so from here on I drop the complex conjugates.

Integrating out the spin, each pair of spatial orbitals $i, j$ gives four spin combinations ($\alpha\alpha$, $\alpha\beta$, $\beta\alpha$, $\beta\beta$). All four contribute a Coulomb term, but only the two same-spin combinations contribute exchange. With the $\frac{1}{2}$ in front, this gives $\frac{1}{2}(4J_{ij} - 2K_{ij}) = 2J_{ij} - K_{ij}$:

$$
E_{\text{RHF}} = 2\sum_{i=1}^{N/2} h_{ii} + \sum_{i=1}^{N/2}\sum_{j=1}^{N/2} \left( 2J_{ij} - K_{ij} \right)
$$

Here everything is now written with spatial orbitals:

$$
h_{ii} = \int \phi_i(\mathbf{r})\, \hat{h}\, \phi_i(\mathbf{r})\, d\mathbf{r},
\qquad
J_{ij} = (ii \vert jj),
\qquad
K_{ij} = (ij \vert ji)
$$

using the **chemists' notation** for two-electron integrals:

$$
(ij \vert kl) = \iint \phi_i(\mathbf{r}_1)\, \phi_j(\mathbf{r}_1)\, \frac{1}{r_{12}}\, \phi_k(\mathbf{r}_2)\, \phi_l(\mathbf{r}_2)\; d\mathbf{r}_1\, d\mathbf{r}_2
$$

In chemists' notation, the left pair of indices belongs to electron 1 and the right pair to electron 2.

A quick check with H₂: there is one occupied orbital, so $E = 2h_{11} + 2J_{11} - K_{11} = 2h_{11} + J_{11}$. Two electrons, each with core energy $h_{11}$, repelling each other once. That is exactly what we expect.

## 6. Deriving the Hartree-Fock equations

Now we minimise $E_{\text{RHF}}$ by changing the orbitals, with one condition: the orbitals must stay orthonormal, $\langle \phi_i \vert \phi_j \rangle = \delta_{ij}$ (where $\delta_{ij}$ is 1 if $i = j$ and 0 otherwise). The standard tool for minimising with conditions attached is the method of **Lagrange multipliers**: add each condition to the quantity being minimised, multiplied by an unknown number, and minimise the combination freely. Here the conditions are labelled by pairs $i, j$, so the multipliers form a symmetric matrix $\varepsilon_{ij}$:

$$
\mathcal{L} = E_{\text{RHF}} - 2\sum_{i=1}^{N/2}\sum_{j=1}^{N/2} \varepsilon_{ij} \left( \langle \phi_i \vert \phi_j \rangle - \delta_{ij} \right)
$$

(The factor 2 is only there to make the result look tidy.) At the minimum, a small change $\phi_i \to \phi_i + \delta\phi_i$ in any orbital must leave $\mathcal{L}$ unchanged to first order. Let us collect the first-order changes term by term, for one particular orbital $\phi_i$.

- $h_{ii}$ contains $\phi_i$ twice, so $\delta h_{ii} = 2\langle \delta\phi_i \vert \hat{h} \vert \phi_i \rangle$, and $\delta \big( 2\sum h \big) = 4\langle \delta\phi_i \vert \hat{h} \vert \phi_i \rangle$.
- The Coulomb sum $\sum_{kl} 2J_{kl}$ contains $\phi_i$ in $J_{il}$ and in $J_{ki}$, twice each. Collecting these gives $8 \sum_j \langle \delta\phi_i \vert \hat{J}_j \vert \phi_i \rangle$.
- The exchange sum $\sum_{kl} K_{kl}$ similarly gives $4 \sum_j \langle \delta\phi_i \vert \hat{K}_j \vert \phi_i \rangle$.
- The constraint term gives $-4\sum_j \varepsilon_{ij} \langle \delta\phi_i \vert \phi_j \rangle$.

Here I have introduced the **Coulomb operator** $\hat{J}_j$ and the **exchange operator** $\hat{K}_j$, defined by what they do to a function $\phi$:

$$
\hat{J}_j\, \phi(\mathbf{r}_1) = \left[ \int \frac{\phi_j(\mathbf{r}_2)^2}{r_{12}}\, d\mathbf{r}_2 \right] \phi(\mathbf{r}_1),
\qquad
\hat{K}_j\, \phi(\mathbf{r}_1) = \left[ \int \frac{\phi_j(\mathbf{r}_2)\, \phi(\mathbf{r}_2)}{r_{12}}\, d\mathbf{r}_2 \right] \phi_j(\mathbf{r}_1)
$$

$\hat{J}_j$ is easy to picture: the bracket is the electrostatic potential created by the electron cloud in orbital $j$. $\hat{K}_j$ is stranger: it swaps $\phi$ and $\phi_j$ between the two positions, which is the operator form of the exchange integral.

Adding up all the first-order changes and setting the total to zero:

$$
4 \left\langle \delta\phi_i \,\middle\vert\, \hat{h}\phi_i + \sum_{j} \left( 2\hat{J}_j - \hat{K}_j \right)\phi_i - \sum_j \varepsilon_{ij}\phi_j \right\rangle = 0
$$

Since $\delta\phi_i$ can be any small change at all, the function on the right of the bar must itself be zero. Defining the **Fock operator**

$$
\hat{f} = \hat{h} + \sum_{j=1}^{N/2} \left( 2\hat{J}_j - \hat{K}_j \right)
$$

we get $\hat{f}\,\phi_i = \sum_j \varepsilon_{ij}\, \phi_j$.

One last simplification. If you mix the occupied orbitals among themselves with a rotation (a unitary transformation), the Slater determinant changes at most by a sign, so the energy and $\hat{f}$ do not change. We can therefore choose the rotation that makes the matrix $\varepsilon_{ij}$ diagonal. This gives the **Hartree-Fock equations**:

$$
\hat{f}\, \phi_i = \varepsilon_i\, \phi_i
$$

This looks just like a one-electron Schrödinger equation, with the Fock operator in place of the Hamiltonian. The $\hat{h}$ part is an electron among the nuclei, and the $2\hat{J}_j - \hat{K}_j$ part is the mean field of all the other electrons, as promised in Section 2. The eigenvalue $\varepsilon_i$ is called the **orbital energy**.

And here is the catch from Section 2 in mathematical form: $\hat{f}$ is built from the occupied orbitals through $\hat{J}_j$ and $\hat{K}_j$, so the equation we are solving depends on its own solution. We must solve it iteratively.

Two useful facts about orbital energies:

- **The total energy is not the sum of orbital energies.** Each $\varepsilon_i = h_{ii} + \sum_j (2J_{ij} - K_{ij})$ includes the repulsion with all the other electrons, so summing them counts every repulsion twice. The correct relation is $E_{\text{RHF}} = \sum_{i} \left( h_{ii} + \varepsilon_i \right)$, which you can check by substituting $\varepsilon_i$.
- **Koopmans' theorem.** Minus the energy of the highest occupied orbital (HOMO) is an estimate of the first ionisation energy of the molecule. I use this in Section 11.

## 7. Bringing in the basis set: the Roothaan-Hall equations

The Hartree-Fock equation is an equation for unknown functions, and solving it directly on a grid is only practical for atoms. In 1951 Roothaan and Hall independently showed how to turn it into a matrix problem: expand each orbital in a basis set, exactly as on the [basis sets](basis_set.md) page,

$$
\phi_i(\mathbf{r}) = \sum_{\nu=1}^{K} C_{\nu i}\, \chi_\nu(\mathbf{r})
$$

Substitute this into $\hat{f}\phi_i = \varepsilon_i\phi_i$:

$$
\sum_{\nu} C_{\nu i}\, \hat{f}\chi_\nu = \varepsilon_i \sum_{\nu} C_{\nu i}\, \chi_\nu
$$

Now multiply both sides by $\chi_\mu(\mathbf{r})$ and integrate over $\mathbf{r}$. This turns the functions into numbers:

$$
\sum_{\nu} F_{\mu\nu}\, C_{\nu i} = \varepsilon_i \sum_{\nu} S_{\mu\nu}\, C_{\nu i},
\qquad
F_{\mu\nu} = \int \chi_\mu\, \hat{f}\, \chi_\nu \, d\mathbf{r},
\qquad
S_{\mu\nu} = \int \chi_\mu\, \chi_\nu \, d\mathbf{r}
$$

Writing this for all orbitals at once gives the **Roothaan-Hall equations**:

$$
\mathbf{F}\mathbf{C} = \mathbf{S}\mathbf{C}\boldsymbol{\varepsilon}
$$

- $\mathbf{F}$ is the **Fock matrix** ($K \times K$);
- $\mathbf{S}$ is the **overlap matrix** ($K \times K$);
- column $i$ of $\mathbf{C}$ holds the coefficients of orbital $i$;
- $\boldsymbol{\varepsilon}$ is a diagonal matrix of orbital energies.

If $\mathbf{S}$ were the identity matrix this would be an ordinary eigenvalue problem, $\mathbf{F}\mathbf{C} = \mathbf{C}\boldsymbol{\varepsilon}$. It is not, because basis functions on neighbouring atoms overlap. This is called a **generalised eigenvalue problem**. The standard trick is to build $\mathbf{X} = \mathbf{S}^{-1/2}$, solve the ordinary problem $\mathbf{F}^{\prime}\mathbf{C}^{\prime} = \mathbf{C}^{\prime}\boldsymbol{\varepsilon}$ with $\mathbf{F}^{\prime} = \mathbf{X}\mathbf{F}\mathbf{X}$, and transform back with $\mathbf{C} = \mathbf{X}\mathbf{C}^{\prime}$. In practice, `scipy.linalg.eigh(F, S)` does all of this for you.

The infinite problem has become a $K \times K$ matrix problem. All that is left is to compute the matrix elements, which are integrals over basis functions.

## 8. The integrals, one at a time

Every integral we need involves only the basis functions, which are fixed. So the program computes them all once at the start and reuses them. For each one I give its meaning and the name PySCF uses for it.

**Overlap**, `int1e_ovlp`:

$$
S_{\mu\nu} = \int \chi_\mu(\mathbf{r})\, \chi_\nu(\mathbf{r})\, d\mathbf{r}
$$

How much two basis functions occupy the same region of space. The diagonal is 1, because the functions are normalised. For two functions on neighbouring atoms it is between 0 and 1, and larger when the atoms are closer.

**Kinetic energy**, `int1e_kin`:

$$
T_{\mu\nu} = \int \chi_\mu(\mathbf{r}) \left( -\frac{1}{2}\nabla^2 \right) \chi_\nu(\mathbf{r})\, d\mathbf{r}
$$

Tight, sharply curved functions have high kinetic energy. This is the quantum mechanical reason an electron does not simply fall into the nucleus: squeezing it into a smaller space raises its kinetic energy.

**Nuclear attraction**, `int1e_nuc`:

$$
V_{\mu\nu} = \int \chi_\mu(\mathbf{r}) \left( -\sum_{A=1}^{M} \frac{Z_A}{\lvert \mathbf{r} - \mathbf{R}_A \rvert} \right) \chi_\nu(\mathbf{r})\, d\mathbf{r}
$$

The attraction to all the nuclei. It is negative. Together, $\mathbf{H}^{\text{core}} = \mathbf{T} + \mathbf{V}$ is the matrix of the core Hamiltonian $\hat{h}$.

**Two-electron repulsion integrals**, `int2e`:

$$
(\mu\nu \vert \lambda\sigma) = \iint \chi_\mu(\mathbf{r}_1)\, \chi_\nu(\mathbf{r}_1)\, \frac{1}{r_{12}}\, \chi_\lambda(\mathbf{r}_2)\, \chi_\sigma(\mathbf{r}_2)\; d\mathbf{r}_1\, d\mathbf{r}_2
$$

Read it as the electrostatic repulsion between two charge clouds: $\chi_\mu\chi_\nu$ for electron 1 and $\chi_\lambda\chi_\sigma$ for electron 2. These integrals are the expensive part of Hartree-Fock. With four indices there are $K^4$ of them. Swapping $\mu \leftrightarrow \nu$, swapping $\lambda \leftrightarrow \sigma$, or swapping the two pairs does not change the value, which cuts the number of distinct integrals by about a factor of 8. For water in cc-pVDZ, $K = 24$, so there are $24^4 = 331{,}776$ integrals, of which 45,150 are distinct. This $K^4$ growth is why Hartree-Fock is said to scale formally as the fourth power of the basis size, and why the [Gaussian product theorem](basis_set.md) matters so much.

## 9. The Fock matrix

To build $F_{\mu\nu}$ we need the Coulomb and exchange operators in the basis. Substituting $\phi_j = \sum_\lambda C_{\lambda j}\chi_\lambda$ into the definitions of $\hat{J}_j$ and $\hat{K}_j$ from Section 6:

$$
\int \chi_\mu\, \hat{J}_j\, \chi_\nu \, d\mathbf{r} = \sum_{\lambda\sigma} C_{\lambda j} C_{\sigma j}\, (\mu\nu \vert \lambda\sigma),
\qquad
\int \chi_\mu\, \hat{K}_j\, \chi_\nu \, d\mathbf{r} = \sum_{\lambda\sigma} C_{\lambda j} C_{\sigma j}\, (\mu\lambda \vert \nu\sigma)
$$

The coefficients of the occupied orbitals always appear in the same combination, so we give it a name, the **density matrix**:

$$
D_{\lambda\sigma} = 2 \sum_{j=1}^{N/2} C_{\lambda j}\, C_{\sigma j}
$$

The factor 2 is the two electrons in each orbital. The density matrix contains everything about where the electrons are: the electron density is $\rho(\mathbf{r}) = \sum_{\lambda\sigma} D_{\lambda\sigma}\, \chi_\lambda(\mathbf{r})\chi_\sigma(\mathbf{r})$, and integrating it gives the number of electrons, $\sum_{\lambda\sigma} D_{\lambda\sigma} S_{\lambda\sigma} = N$.

Summing over the occupied orbitals $j$ with the factors $2\hat{J}_j - \hat{K}_j$ then gives the Fock matrix:

$$
F_{\mu\nu} = H^{\text{core}}_{\mu\nu} + \sum_{\lambda\sigma} D_{\lambda\sigma} \left[ (\mu\nu \vert \lambda\sigma) - \frac{1}{2} (\mu\lambda \vert \nu\sigma) \right]
$$

or in short, $\mathbf{F} = \mathbf{H}^{\text{core}} + \mathbf{J} - \frac{1}{2}\mathbf{K}$, with $J_{\mu\nu} = \sum_{\lambda\sigma} D_{\lambda\sigma}(\mu\nu \vert \lambda\sigma)$ and $K_{\mu\nu} = \sum_{\lambda\sigma} D_{\lambda\sigma}(\mu\lambda \vert \nu\sigma)$. These two lines become two lines of Python in Section 11.

Finally, the energy. Starting from $E_{\text{RHF}} = \sum_i (h_{ii} + \varepsilon_i)$ and writing $h_{ii}$ and $\varepsilon_i$ in the basis:

$$
E_{\text{elec}} = \frac{1}{2} \sum_{\mu\nu} D_{\mu\nu} \left( H^{\text{core}}_{\mu\nu} + F_{\mu\nu} \right),
\qquad
E_{\text{total}} = E_{\text{elec}} + E_{\text{nuc}}
$$

## 10. The self-consistent field procedure

Everything is now in place. The Hartree-Fock algorithm is:

1. Choose the molecule (the positions of the nuclei) and a basis set.
2. Compute $\mathbf{S}$, $\mathbf{H}^{\text{core}}$ and all $(\mu\nu \vert \lambda\sigma)$, once.
3. Guess a density matrix $\mathbf{D}$. The simplest guess ignores electron repulsion completely: solve $\mathbf{H}^{\text{core}}\mathbf{C} = \mathbf{S}\mathbf{C}\boldsymbol{\varepsilon}$ and build $\mathbf{D}$ from it. (PySCF uses a better default guess, built from the densities of the free atoms.)
4. Build $\mathbf{F}$ from $\mathbf{D}$.
5. Solve $\mathbf{F}\mathbf{C} = \mathbf{S}\mathbf{C}\boldsymbol{\varepsilon}$.
6. Build a new $\mathbf{D}$ from the $N/2$ orbitals with the lowest energies (filling orbitals from the bottom up, the aufbau principle).
7. Compute the energy. If the energy and $\mathbf{D}$ have stopped changing, stop: the solution is self-consistent. Otherwise go back to step 4.

Sometimes this simple loop oscillates instead of settling down. Real programs, including PySCF, speed up and stabilise convergence with a technique called DIIS, which builds the next Fock matrix from a clever combination of the previous ones.

## 11. Hartree-Fock with PySCF

All the examples below were run with PySCF 2.14.0, and the outputs shown are the ones I got. Your last digits may differ slightly with other versions.

### A first calculation

```python
from pyscf import gto, scf

mol = gto.M(atom='H 0 0 0; H 0 0 0.74', basis='sto-3g')
mf = scf.RHF(mol)
energy = mf.kernel()

print('converged:', mf.converged)
print('total energy (Hartree):', energy)
print('orbital energies (Hartree):', mf.mo_energy)
print('MO coefficients (one orbital per column):')
print(mf.mo_coeff)
```

`gto.M` builds the molecule: two hydrogen atoms 0.74 Å apart, with the STO-3G basis. `scf.RHF(mol)` sets up a restricted Hartree-Fock calculation, and `kernel()` runs the SCF loop of Section 10 and returns the total energy:

```text
converged SCF energy = -1.11675930739643
converged: True
total energy (Hartree): -1.1167593073964255
orbital energies (Hartree): [-0.57855386  0.67114349]
MO coefficients (one orbital per column):
[[ 0.54884228 -1.21245192]
 [ 0.54884228  1.21245192]]
```

STO-3G gives each hydrogen one 1s function, so $K = 2$ and there are two molecular orbitals. Look at the coefficients. The first orbital has a negative energy, is occupied by both electrons, and has coefficients of the same sign on both atoms: this is the **bonding** orbital, $\chi_{1s}^{A} + \chi_{1s}^{B}$. The second has a positive energy, is empty, and has coefficients of opposite sign: the **antibonding** orbital, $\chi_{1s}^{A} - \chi_{1s}^{B}$. (The overall sign of each column is arbitrary, so do not worry if both of yours are flipped.)

### Looking at the integrals

```python
import numpy as np
from pyscf import gto

mol = gto.M(atom='H 0 0 0; H 0 0 0.74', basis='sto-3g')

S = mol.intor('int1e_ovlp')   # overlap
T = mol.intor('int1e_kin')    # kinetic energy
V = mol.intor('int1e_nuc')    # nuclear attraction
eri = mol.intor('int2e')      # two-electron integrals (mu nu|lambda sigma)

np.set_printoptions(precision=4, suppress=True)
print('S =\n', S)
print('T =\n', T)
print('V =\n', V)
print('eri has shape', eri.shape)
print('(11|11) =', eri[0, 0, 0, 0])
print('(11|22) =', eri[0, 0, 1, 1])
print('nuclear repulsion =', mol.energy_nuc())
```

The output:

```text
S =
 [[1.     0.6599]
 [0.6599 1.    ]]
T =
 [[0.76  0.237]
 [0.237 0.76 ]]
V =
 [[-1.881  -1.1963]
 [-1.1963 -1.881 ]]
eri has shape (2, 2, 2, 2)
(11|11) = 0.7746059439198978
(11|22) = 0.5699948822432332
nuclear repulsion = 0.7151043390810812
```

Python counts from 0, so `eri[0, 0, 1, 1]` is $(11 \vert 22)$. Things to notice, and to connect with Section 8:

- $\mathbf{S}$ has 1 on the diagonal, and the two 1s functions overlap by about two thirds;
- $\mathbf{T}$ is positive and $\mathbf{V}$ is negative;
- `eri` is a $2 \times 2 \times 2 \times 2$ array, stored in chemists' notation;
- $(11 \vert 22)$ is smaller than $(11 \vert 11)$: two charge clouds on different atoms are further apart than a cloud with itself, so they repel less.

### Writing our own Hartree-Fock program

This is the part I think teaches the most. The script below implements Sections 9 and 10 directly, for water, and compares the answer with PySCF. Every line maps onto an equation above.

```python
import numpy as np
from scipy.linalg import eigh
from pyscf import gto, scf

mol = gto.M(
    atom='''
    O   0.000   0.000   0.000
    H   0.000   0.757   0.587
    H   0.000  -0.757   0.587
    ''',
    basis='sto-3g',
)

# Step 2: all the integrals, computed once
S = mol.intor('int1e_ovlp')
Hcore = mol.intor('int1e_kin') + mol.intor('int1e_nuc')
eri = mol.intor('int2e')
E_nuc = mol.energy_nuc()
nocc = mol.nelectron // 2          # number of doubly occupied orbitals

def density(C):
    """D = 2 * sum over occupied orbitals of C C^T"""
    C_occ = C[:, :nocc]
    return 2 * C_occ @ C_occ.T

# Step 3: core guess, ignoring electron repulsion
eps, C = eigh(Hcore, S)
D = density(C)

E_old = 0.0
for iteration in range(1, 101):
    # Step 4: Fock matrix, F = Hcore + J - K/2
    J = np.einsum('pqrs,rs->pq', eri, D)   # J_pq = sum_rs (pq|rs) D_rs
    K = np.einsum('prqs,rs->pq', eri, D)   # K_pq = sum_rs (pr|qs) D_rs
    F = Hcore + J - 0.5 * K

    # Step 7: energy
    E_elec = 0.5 * np.sum(D * (Hcore + F))
    E_total = E_elec + E_nuc
    print(f'iteration {iteration:3d}   E = {E_total:.10f}')
    if abs(E_total - E_old) < 1e-10:
        break
    E_old = E_total

    # Steps 5 and 6: solve FC = SCe, build the new density
    eps, C = eigh(F, S)
    D = density(C)

mf = scf.RHF(mol)
mf.kernel()
print('our energy:  ', E_total)
print('PySCF energy:', mf.e_tot)
```

The output:

```text
iteration   1   E = -73.2327077545
iteration   2   E = -74.9457239958
iteration   3   E = -74.9622044928
...
iteration  13   E = -74.9630631297
iteration  14   E = -74.9630631297
converged SCF energy = -74.9630631297276
our energy:   -74.9630631297234
PySCF energy: -74.9630631297276
```

Our thirty-odd lines of Python agree with PySCF to about $10^{-11}$ Hartree. Our loop has no DIIS, so it needs more iterations than PySCF (14 here, against 6 for PySCF), and for harder molecules it may not converge at all. For simplicity it also only checks the energy for convergence, while real programs also check that $\mathbf{D}$ has stopped changing. Once it works, try a bigger basis set, or stretch one of the O-H bonds and watch how the number of iterations changes.

### Orbital energies and Koopmans' theorem

```python
from pyscf import gto, scf

mol = gto.M(
    atom='''
    O   0.000   0.000   0.000
    H   0.000   0.757   0.587
    H   0.000  -0.757   0.587
    ''',
    basis='cc-pvdz',
)
mf = scf.RHF(mol).run()

hartree_to_ev = 27.211386
nocc = mol.nelectron // 2
homo = mf.mo_energy[nocc - 1]
lumo = mf.mo_energy[nocc]

print(f'HOMO energy: {homo:.4f} Hartree = {homo * hartree_to_ev:.2f} eV')
print(f'LUMO energy: {lumo:.4f} Hartree = {lumo * hartree_to_ev:.2f} eV')
print(f'Koopmans estimate of the ionisation energy: {-homo * hartree_to_ev:.2f} eV')
```

The output:

```text
converged SCF energy = -76.0267656731175
HOMO energy: -0.4931 Hartree = -13.42 eV
LUMO energy: 0.1854 Hartree = 5.05 eV
Koopmans estimate of the ionisation energy: 13.42 eV
```

Compare the Koopmans estimate, 13.4 eV, with the experimental first ionisation energy of water, about 12.6 eV. It is in the right region, but not exact. Two errors are hidden in it, and they partly cancel: the other orbitals would relax if an electron were really removed (which Koopmans' theorem ignores), and Hartree-Fock lacks correlation energy.

## 12. What Hartree-Fock gets wrong

Hartree-Fock gives reasonable molecular geometries and a qualitatively correct picture of orbitals and bonding. Its errors all come from the mean-field approximation of Section 2.

- **Correlation energy.** The difference between the exact energy and the Hartree-Fock limit is the correlation energy. It is a small fraction of the total energy, but comparable in size to reaction energies and barrier heights, so Hartree-Fock reaction energies are often not accurate enough.
- **Dispersion.** The weak attraction between neutral molecules that holds, for example, noble-gas atoms together in a liquid (London dispersion) is a pure correlation effect. Hartree-Fock misses it completely.
- **Breaking bonds.** RHF forces both electrons of a bond into the same spatial orbital. Stretch H₂ to a large distance and RHF still keeps both electrons together, giving a wavefunction that is half "H⁺ plus H⁻". The energy at dissociation is therefore much too high.

Going beyond Hartree-Fock means adding correlation. [MP2](mp2.md) and [CCSD](ccsd.md) start from the Hartree-Fock orbitals and correct the wavefunction; [DFT](dft.md) takes a different route through the electron density. All of them reuse the integrals, the basis sets and much of the machinery described here.

## Summary

- Hartree-Fock approximates the wavefunction as a single Slater determinant, which makes it antisymmetric.
- Each electron moves in the average field of the others; the Fock operator is the Hamiltonian for that situation.
- Expanding the orbitals in a basis turns the Hartree-Fock equations into the matrix equation $\mathbf{F}\mathbf{C} = \mathbf{S}\mathbf{C}\boldsymbol{\varepsilon}$.
- The matrices are built from four kinds of integrals: overlap, kinetic, nuclear attraction and two-electron repulsion.
- Because $\mathbf{F}$ depends on $\mathbf{C}$, the equations are solved iteratively until self-consistent.

## Further reading

- A. Szabo and N. S. Ostlund, *Modern Quantum Chemistry*, Dover, chapters 2 and 3. The classic, careful derivation; the notation on this page follows it closely.
- The [Crawford group programming projects](https://github.com/CrawfordGroup/ProgrammingProjects): project 3 walks you through writing a Hartree-Fock program, step by step.

👉 More on HF: [Wikipedia](https://en.wikipedia.org/wiki/Hartree%E2%80%93Fock_method)  
👉 PySCF docs: [PySCF Hartree–Fock](https://pyscf.org/user/hf.html)
