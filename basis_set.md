---
layout: default
title: Basis Sets
---

# Basis Sets

Before I can explain Hartree-Fock properly, I need to explain basis sets, because every practical Hartree-Fock calculation is built on one. If you have never studied physical chemistry, this is a good place to start. The idea is much simpler than the name suggests.

## The idea in one sentence

A basis set is a fixed collection of simple, known functions, and we build every unknown orbital as a weighted sum of them.

## An analogy: mixing paint

Suppose you need to match the colour of a wall, but you only own red, yellow and blue paint. You cannot buy the exact colour, so you mix: some red, a little yellow, more blue. What you actually have to work out is the *recipe*, the amount of each paint. If you own more colours of paint you can match the wall more closely, but working out the recipe takes more effort.

A basis set works in the same way:

- the colour of the wall is the true molecular orbital, which we do not know;
- the paints are the basis functions, which we choose before the calculation starts;
- the recipe is a list of numbers called coefficients, which the calculation finds for us;
- more basis functions give a better orbital, and a more expensive calculation.

## What is an orbital, mathematically?

An orbital is a function. You give it a point in space, $\mathbf{r} = (x, y, z)$, and it gives you back a number, $\phi(\mathbf{r})$ (say, atomic orbital). From the LCAO (Linear Combination of Atomic Orbitals), perspective, an MO (Molecular Orbital) can be formed from the linear combination of the atomic orbitals (refer to the recipe, here again). Finding the true MO exactly would mean finding its value at infinitely many points, which no computer can do. So we restrict ourselves to orbitals of a particular form:

$$
\phi_i(\mathbf{r}) = \sum_{\mu=1}^{K} C_{\mu i}\thinspace \chi_\mu(\mathbf{r})
$$

Here:

- $\chi_\mu$ (the Greek letter "chi") is basis function number $\mu$, and there are $K$ of them;
- $C_{\mu i}$ is how much of basis function $\mu$ goes into orbital $i$;
- the impossible problem "find a function" has become the finite problem "find $K$ numbers for each orbital".

This is called the **linear combination of atomic orbitals** (LCAO) approximation, because the basis functions are centred on the atoms and look like atomic orbitals. A simple example: in the H₂ molecule, the lowest orbital is roughly a 1s function on atom A plus a 1s function on atom B, mixed in equal amounts:

$$
\phi_1 \approx C_{A}\thinspace \chi_{1s}^{A} + C_{B}\thinspace \chi_{1s}^{B}, \qquad C_A = C_B
$$

## Slater functions and Gaussian functions

What should the building blocks look like? For the hydrogen atom we know the exact 1s orbital. In atomic units (lengths measured in bohr, 1 bohr = 0.529 Angstroem) it is

$$
\psi_{1s}(r) = \frac{1}{\sqrt{\pi}}\thinspace e^{-r}
$$

where $r$ is the distance from the nucleus. Functions with the shape $e^{-\zeta r}$ are called **Slater-type orbitals** (STOs). The number $\zeta$ ("zeta") controls how tight or spread out the function is. Slater functions have the right shape: a sharp point (a "cusp") at the nucleus, and a slow decay far away from it. Unfortunately, they have one big practical problem. The Hartree-Fock or other methods method needs integrals that involve four basis functions at once, possibly sitting on four different atoms (I explain these in the [Hartree-Fock chapter](hf.md)). With Slater functions those integrals are very hard to compute.

**Gaussian-type orbitals** (GTOs) use the shape $e^{-\alpha r^2}$ instead (compare with Gaussian or a normal distribution curve). On their own they have the wrong shape: they are flat at the nucleus (no cusp) and die away too quickly at large distance. But they have one enormous advantage, the **Gaussian product theorem**: the product of two Gaussians centred on two different atoms is another Gaussian, centred at a point between them.

$$
e^{-\alpha \lvert \mathbf{r}-\mathbf{A} \rvert^2}\mskip5mu e^{-\beta \lvert \mathbf{r}-\mathbf{B} \rvert^2}
= K_{AB}\mskip5mu e^{-(\alpha+\beta) \lvert \mathbf{r}-\mathbf{P} \rvert^2},
\qquad
\mathbf{P} = \frac{\alpha \mathbf{A} + \beta \mathbf{B}}{\alpha + \beta},
\qquad
K_{AB} = e^{-\frac{\alpha\beta}{\alpha+\beta} \lvert \mathbf{A}-\mathbf{B} \rvert^2}
$$

Because of this, an integral over four Gaussians on four atoms collapses into an integral over two Gaussians, and those have closed formulas that a computer can evaluate very quickly. This is why nearly every quantum chemistry program, including PySCF, uses Gaussian functions.

## Contracted Gaussians: fixing the shape

To fix the wrong shape, we add up several Gaussians of different widths with a fixed recipe. Each individual Gaussian is called a **primitive**, and the sum is called a **contracted Gaussian**:

$$
\chi(\mathbf{r}) = \sum_{k=1}^{L} d_k\thinspace N_k\thinspace e^{-\alpha_k r^2}
$$

- $\alpha_k$ are the **exponents** (how tight each primitive is);
- $d_k$ are the **contraction coefficients** (how much of each primitive);
- $N_k$ are normalisation constants, which make each primitive "size one".

The important point is that $\alpha_k$ and $d_k$ are fixed by the people who designed the basis set. They never change during your calculation. The calculation only finds the coefficients $C_{\mu i}$ that mix the contracted functions into molecular orbitals.

The smallest common basis set, **STO-3G**, approximates each Slater function by a contraction of 3 Gaussians. The name literally means "Slater-type orbital, approximated with 3 Gaussians".

### Try it: how good is STO-3G?

The script below compares a Slater 1s function with its STO-3G imitation for hydrogen. The value $\zeta = 1.24$ is the one used for hydrogen in STO-3G (a hydrogen atom inside a molecule is a little tighter than a free one). The exponents and coefficients are the standard STO-3G values; the first line also prints them straight from PySCF's basis-set library so you can compare.

```python
import numpy as np
import matplotlib.pyplot as plt
from pyscf import gto

# STO-3G for hydrogen, as stored in PySCF:
# [angular momentum, [exponent, coefficient], [exponent, coefficient], ...]
print(gto.basis.load('sto-3g', 'H'))

exps = np.array([3.42525091, 0.62391373, 0.16885540])
coeffs = np.array([0.15432897, 0.53532814, 0.44463454])

r = np.linspace(0, 4, 400)  # distance from the nucleus, in bohr

# Slater 1s function with zeta = 1.24, normalised
zeta = 1.24
slater = np.sqrt(zeta**3 / np.pi) * np.exp(-zeta * r)

# STO-3G: three normalised Gaussian primitives, mixed with fixed coefficients
norms = (2 * exps / np.pi) ** 0.75
sto3g = np.sum(coeffs[:, None] * norms[:, None] * np.exp(-exps[:, None] * r**2), axis=0)

plt.plot(r, slater, label='Slater 1s (zeta = 1.24)')
plt.plot(r, sto3g, '--', label='STO-3G')
plt.xlabel('distance from nucleus (bohr)')
plt.ylabel('value of the function')
plt.legend()
plt.show()
```

In the printed list, the leading `0` is the angular momentum (0 means an s function), followed by three exponent and coefficient pairs. In the plot, the two curves are almost indistinguishable beyond about 0.3 bohr. The clear difference is right at the nucleus: the Slater function comes to a sharp point (about 0.78 at $r = 0$), while the Gaussian sum is rounded and flat (about 0.63). The tails differ too, because Gaussians die away faster than $e^{-\zeta r}$, but only far out where both functions are already tiny.

## Shells: s, p, d and f functions

So far I have only talked about spherical s functions. Real basis sets also contain functions with a direction, which look like the p, d and f orbitals you may have seen in general chemistry:

- an **s** shell contains 1 function;
- a **p** shell contains 3 functions ($p_x$, $p_y$, $p_z$);
- a **d** shell contains 5 functions;
- an **f** shell contains 7 functions.

A shell is a group of functions that share the same exponents and contraction coefficients and only differ in direction. (Using cartesian reprentation, one would have 6 Cartesian d functions instead of 5. PySCF uses 5 by default; setting `cart=True` in `gto.M` switches to 6. If you use Gaussian suite of programs, using *6D 10F* in the route allow to use cartesian basis sets)

## Reading basis-set names

Basis-set names look cryptic, but each part means something.

**Minimal basis (STO-3G).** One function for each orbital that is occupied in the free atom. Hydrogen gets a 1s function: 1 function. Oxygen gets 1s, 2s, $2p_x$, $2p_y$, $2p_z$: 5 functions. So water, H₂O, has $5 + 1 + 1 = 7$ basis functions.

**Split valence (6-31G).** Read the name as "6 / 3 1". Each core orbital (such as the oxygen 1s) is one contraction of 6 primitives. Each valence orbital is split into two functions: an inner one made of 3 primitives and an outer one made of 1 primitive. Why split? Atoms change size when they form bonds. With one tight and one loose function per valence orbital, the calculation can make an orbital bigger or smaller by mixing them in different amounts. Water in 6-31G has $9 + 2 + 2 = 13$ functions.

**Polarisation functions (6-31G\* or 6-31G(d)).** These add d functions on heavy atoms; a second star, or (d,p), also adds p functions on hydrogen. In a bond, the electron cloud of an atom is pulled towards its neighbour. An s function on its own cannot shift off-centre, but mixing in a little p function can; in the same way, mixing d into p lets p orbitals bend.

**Diffuse functions (6-31+G).** A plus sign adds very spread-out s and p functions on heavy atoms; a double plus also adds them on hydrogen. They matter for anions, excited states and weak interactions between molecules, where electrons spend time far from the nuclei.

**Correlation-consistent sets (cc-pVDZ, cc-pVTZ, cc-pVQZ).** These were designed by Dunning so that results improve in a steady, predictable way as you go up the series. The name reads "correlation-consistent, polarised valence, double / triple / quadruple zeta". "Double zeta" simply means two functions per valence orbital, "triple zeta" three, and so on (the name comes from the Slater exponent $\zeta$). An `aug-` prefix, as in aug-cc-pVDZ, adds diffuse functions.

| Basis set | Kind | Basis functions for water |
|-----------|------|---------------------------|
| STO-3G | minimal | 7 |
| 6-31G | split valence | 13 |
| cc-pVDZ | double zeta + polarisation | 24 |
| cc-pVTZ | triple zeta + polarisation | 58 |

For cc-pVDZ, oxygen has 3 s, 2 p and 1 d shells ($3 + 6 + 5 = 14$ functions) and each hydrogen has 2 s and 1 p shells ($2 + 3 = 5$), giving $14 + 5 + 5 = 24$.

## Inspecting a basis set in PySCF

```python
from pyscf import gto

mol = gto.M(
    atom='''
    O   0.000   0.000   0.000
    H   0.000   0.757   0.587
    H   0.000  -0.757   0.587
    ''',
    basis='sto-3g',
)

print('number of basis functions:', mol.nao_nr())
for label in mol.ao_labels():
    print(label)
```

The coordinates are in Angstroem by default. The labels tell you which atom each function sits on and what type it is, for example the oxygen 1s, 2s and three 2p functions followed by one 1s function on each hydrogen. Change `basis='sto-3g'` to `basis='cc-pvdz'` and run it again: you should get 24 functions, including d functions on oxygen and p functions on the hydrogens.

You can also use different basis sets on different elements, for example `basis={'O': 'cc-pvdz', 'H': 'sto-3g'}`.

## More functions, better answer: the basis-set limit

The script below runs a Hartree-Fock calculation on water with bigger and bigger basis sets.

```python
from pyscf import gto, scf

water = '''
O   0.000   0.000   0.000
H   0.000   0.757   0.587
H   0.000  -0.757   0.587
'''

for basis in ['sto-3g', '6-31g', 'cc-pvdz', 'cc-pvtz', 'cc-pvqz']:
    mol = gto.M(atom=water, basis=basis, verbose=0)
    mf = scf.RHF(mol)
    energy = mf.kernel()
    print(f'{basis:8s} {mol.nao_nr():4d} functions   E = {energy:.6f} Hartree')
```

With PySCF 2.14.0 I get:

```text
sto-3g      7 functions   E = -74.963063 Hartree
6-31g      13 functions   E = -75.983948 Hartree
cc-pvdz    24 functions   E = -76.026766 Hartree
cc-pvtz    58 functions   E = -76.057114 Hartree
cc-pvqz   115 functions   E = -76.064777 Hartree
```

There are two things to see. First, the energy goes down as the basis grows: the energy drops a lot at first, then by smaller and smaller amounts. This is the **variational principle** at work (see the [Hartree-Fock chapter](hf.md)): a more flexible recipe can only get closer to the true orbitals, never further away, as long as the bigger basis contains everything the smaller one could describe. The value the energy settles towards is called the **Hartree-Fock limit**. Second, the calculations get slower. The number of two-electron integrals grows roughly as $K^4$, so doubling the number of basis functions makes that part of the work about 16 times bigger.

Choosing a basis set is therefore always a compromise between accuracy and cost.

## Summary

- A basis set is a fixed set of functions; orbitals are built as weighted sums of them.
- The calculation only finds the mixing coefficients $C_{\mu i}$; the basis functions themselves are fixed.
- Gaussian functions have the wrong shape but make integrals cheap, so we glue several together into contracted functions.
- Bigger basis sets give lower, more accurate energies and cost much more.

Next: [Hartree-Fock](hf.md), where these basis functions are put to work.

## Further reading

- A. Szabo and N. S. Ostlund, *Modern Quantum Chemistry*, Dover, section 3.6 on basis sets.
- [Basis Set Exchange](https://www.basissetexchange.org/): download almost any basis set and see its exponents and coefficients.
- PySCF docs: [Molecular structure and basis sets](https://pyscf.org/user/gto.html)
