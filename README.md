# Directed Scale-Free Network Generator (Fortran)

This repository implements a directed scale-free graph model based on

> Bollobás, Borgs, Chayes, Riordan  
> *Directed Scale-Free Graphs* (2003)

The model generates a growing directed network with preferential attachment
depending on in-degree and out-degree.

---

## Model Description

At each discrete time step, one directed edge is added.

With probabilities:

- **(A)** α : Add a new vertex `v` and an edge `v → w`
  - `w` chosen proportional to `d_in(w) + δ_in`

- **(B)** β : Add an edge `v → w` between existing vertices
  - `v` chosen proportional to `d_out(v) + δ_out`
  - `w` chosen proportional to `d_in(w) + δ_in`

- **(C)** γ : Add a new vertex `w` and an edge `v → w`
  - `v` chosen proportional to `d_out(v) + δ_out`

with

α + β + γ = 1

---

## Features

- Directed graph
- No self-loops
- No multiple edges
- Adjacency matrix output
- User-defined final number of vertices and edges
- Initial condition:
  - 2 vertices
  - 1 edge: 1 → 2

---

## Parameters

| Parameter | Meaning |
|-----------|---------|
| α | Probability of rule (A) |
| β | Probability of rule (B) |
| γ | Probability of rule (C) |
| δ_in | Initial attractiveness for in-degree |
| δ_out | Initial attractiveness for out-degree |
| N_final | Final number of vertices |
| E_final | Final number of edges |

---

## Compilation

Using Intel ifx:

```bash
ifx -O src/DIRSF.f90 -o DIRSF.exe -mcmodel=medium

