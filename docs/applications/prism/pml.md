# Absorbing layers

Two layers absorb the waves that reach the domain faces: the **perfectly matched layer** (`[PML]`), an auxiliary-equation
PML of four formulations, and the **far-wave layer** (`[fWLayer]`), a damping of the outgoing wave applied after every
step. Each is set face by face and works on top of the face's [boundary condition](./boundary-conditions). The keys are
in the [input reference](./input#pml).

## Perfectly matched layer

### Layer geometry and storage

`[PML] width` is a physical width; on each active face (`x_minus_layer` … `z_plus_layer`) the layer covers the cells
within it. The PML state lives only on the boundary blocks that carry a layer, compacted face by face: each face stores
four split variables, electric first and magnetic second, e.g. on an $x$ face
$(\psi_{E_y}, \psi_{E_z}, \psi_{H_y}, \psi_{H_z})$, over the layer cells of each block.

### Profiles

With $d$ the distance of a cell from the boundary cell and $L$ the largest such distance inside the layer among the blocks of the rank (no MPI reduction: a defect when the layer spans several ranks, [#108](https://github.com/szaghi/adam/issues/108)), the depth
$\xi = \operatorname{clamp}(1 - d/L, 0, 1)$ is 1 at the boundary and 0 at the inner edge of the layer:

| `PML_type` | Damping $\gamma$ | $\alpha$, $\kappa$ | Parameters |
|---|---|---|---|
| `CLASSIC` | $\gamma_{\max}\,\xi^{m}$ | $\alpha = 0$, $\kappa = 1$ | `gamma_max`, `gamma_exponent` ($m$) |
| `CFS` (CPML) | $\gamma_{\max}\,\xi^{2}$ | $\alpha = \alpha_{\max}(1-\xi)^2$, $\kappa = 1 + (\kappa_{\max}-1)\,\xi^2$ | `gamma_max`, `alfa_max`, `k_max` |
| `BERMUDEZ` | $\beta / (W d/L + \epsilon)^{m}$, $\epsilon = 10^{-6}$ | $\alpha = 0$, $\kappa = 1$ | `beta`, `gamma_exponent` |
| `CLASSIC_DIRECT` | $\gamma_{\max}\,\xi^{m}$ | — | `gamma_max`, `gamma_exponent` |

$W$ is the layer width. The Bermúdez profile is the unbounded one, singular at the boundary.

### Auxiliary equations

`CLASSIC`, `CFS` and `BERMUDEZ` are auxiliary-differential-equation (ADE) PMLs. On an $x$ face each split variable obeys

$$\frac{\partial \psi}{\partial t} = \frac{\gamma}{\kappa^2}\,\frac{\partial F}{\partial x} - \Big(\alpha + \frac{\gamma}{\kappa}\Big)\psi,$$

with $F$ the tangential field it is paired with ($D_y/\varepsilon_0$, $D_z/\varepsilon_0$, $B_y/\mu_0$, $B_z/\mu_0$), and the
Maxwell residual of the tangential components in the layer gains the stretched derivative and the split variable:

$$R \leftarrow R + \Big(\frac{1}{\kappa} - 1\Big)\big(\pm\partial_x F'\big) + \psi',$$

$F'$ and $\psi'$ being the partners of the curl term. The split variables advance with the field, in the same SSP
Runge–Kutta stages (`prism_rk_pml_object`).

**`CLASSIC_DIRECT` is not an ADE PML.** It skips the auxiliary variables and their Runge–Kutta integration and damps the
tangential components directly, $R \leftarrow R - \gamma\,q$. It exists as a testing instrument.

### Time step

A PML caps the time step at $\Delta t \le C_\gamma/\gamma_{\max}^{\rm eff}$, with $\gamma_{\max}^{\rm eff}$ the largest
damping over the layer cells and `C_gamma` (default 0.1) the safety factor.

### Constraints

| Backend | Accepted |
|---|---|
| CPU | `CLASSIC`, `CLASSIC_DIRECT`, `BERMUDEZ`, `CFS`; `fd_centered`; an SSP Runge–Kutta scheme; no divergence correction, no constrained transport |
| FNL | `fd_centered`; `runge_kutta` with an SSP scheme; no hyperbolic correction, no constrained transport |

Anything else is fatal at initialization, and so is a section with every face flag `.false.` or a missing parameter
of the chosen type.

## Far-wave layer

The far-wave layer damps the outgoing part of the wave in a band of cells along the chosen faces. `[fWLayer] width` is
physical: the cells per face are $C = \min(n, \lceil W/\Delta s\rceil)$. After every step (and once at t0) each pair of
tangential components $(\alpha_D, \beta_B)$, $(\beta_D, \alpha_B)$ of a layer cell is corrected towards the outgoing
characteristic,

$$D'_{\alpha} = \frac{1}{2\sqrt{\mu_0}}\Big[s_2(f-1)\,B_{\beta}\sqrt{\varepsilon_0} + (f+1)\,D_{\alpha}\sqrt{\mu_0}\Big],$$

and likewise for the other three, with $s_2 = \pm1$ the face sign. The factor $f$ is evaluated per cell, from the
distance $d$ of the cell from the boundary cell over the layer extent $L$:

$$f = \frac{1}{f_i}\log_{10}\!\Big[r\,(10^{f_i} - 1) + 1\Big], \quad r = \operatorname{clamp}(d/L, 0, 1), \quad
f_i = \begin{cases}\frac{1}{150}(-7C^2 + 255C + 250) & C < 40\\ 25 & C \ge 40\end{cases}$$

so $f = 0$ at the boundary and $f = 1$ (no correction) at the inner edge. $W = 0$ disables the layer. One function
(`compute_fwl_factor`) serves both backends. The removed key `C` (a cell count) is rejected by the removed-keys lint.
