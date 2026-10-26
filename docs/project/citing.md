# Citing ADAM

If you use ADAM in work that leads to a scientific publication, please cite the following paper:

> S. Zaghi, F. Salvadore, A. Di Mascio, G. Rossi —
> **Efficient GPU parallelization of adaptive mesh refinement technique for high-order compressible solver with immersed boundary** —
> *Computers and Fluids*, 266 (2023) 106040.
> DOI: [10.1016/j.compfluid.2023.106040](https://doi.org/10.1016/j.compfluid.2023.106040)

The paper describes the ADAM framework architecture, the AMR/IB coupling strategy, the GPU parallelization approach (CUDA Fortran), and demonstrates strong scaling on a shock–sphere interaction benchmark. A preprint is available in [docs/papers/zaghi-2023-computer_fluids.pdf](../papers/zaghi-2023-computer_fluids.pdf).

BibTeX entry:

```bibtex
@article{zaghi2023adam,
  author  = {Zaghi, S. and Salvadore, F. and {Di Mascio}, A. and Rossi, G.},
  title   = {Efficient {GPU} parallelization of adaptive mesh refinement technique
             for high-order compressible solver with immersed boundary},
  journal = {Computers \& Fluids},
  volume  = {266},
  pages   = {106040},
  year    = {2023},
  doi     = {10.1016/j.compfluid.2023.106040},
}
```

## Authors

- Andrea di Mascio — [andrea.dimascio@univaq.it](mailto:andrea.dimascio@univaq.it)
- Federico Negro — [federico.negro.01@gmail.com](mailto:federico.negro.01@gmail.com)
- Giacomo Rossi — [giacomo.rossi@amd.com](mailto:giacomo.rossi@amd.com)
- Francesco Salvadore — [f.salvadore@cineca.it](mailto:f.salvadore@cineca.it)
- Stefano Zaghi — [stefano.zaghi@cnr.it](mailto:stefano.zaghi@cnr.it)

ADAM is dual-licensed under the [MIT License](/project/license#mit) and the [GNU Lesser General Public License v3.0](/project/license#lgplv3) (LGPLv3); you may choose either.

> Copyright (C) Andrea Di Mascio, Federico Negro, Giacomo Rossi, Francesco Salvadore, Stefano Zaghi.
