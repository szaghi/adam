---
title: adam_prism_parameters
---

# adam_prism_parameters

> PRISM, general parameters.

**Source**: `src/app/prism/common/adam_prism_parameters.F90`

**Dependencies**

```mermaid
graph LR
  adam_prism_parameters["adam_prism_parameters"] --> penf["penf"]
```

## Variables

| Name | Type | Attributes | Description |
|------|------|------------|-------------|
| `NV_MAX` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | parameter | Maximum number of variables for static arrays dimensioning. |
| `MU0` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | parameter | Vacuum magnetic permeability. |
| `EPS0` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | parameter | Vacuum dielectric constant. |
| `C0` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | parameter | Vacuum speed of light. |
| `E_CHARGE` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | parameter | Elettron charge value. |
| `E_MASS` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | parameter | Elettron mass value. |
| `Q_OVER_M` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | parameter | Q/m value for elettron. |
| `K_B` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | parameter | Boltzmann coefficient. |
| `T_E` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | parameter | Plasma Temperature. |
| `VTH` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | parameter | Thermal Velocity value. |
| `PI` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | parameter | Pi value. |
| `MU0_SQ` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | parameter | Mu0 square root |
| `MU0_SQ_I2` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | parameter | Half of MU0 square root inverse |
| `EPS0_SQ` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | parameter | EPS0 square root |
| `EPS0_SQ_I2` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | parameter | Half of EPS0 square root inverse |
| `EV` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | target | Eigenvalues. |
| `ER` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | target | Right eigenvectors. |
| `EL` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | target | Right eigenvectors. |
| `IEV` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | target | Identity eigenvalues. |
| `IEV_D` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | target | Identity eigenvalues with D divergence cleaning. |
| `IEV_B` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | target | Identity eigenvalues with B divergence cleaning. |
| `IEV_D_B` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | target | Identity eigenvalues with D & B divergence cleaning. |
| `IERL` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | target | Identity eigenvectors. |
| `IERL_D` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | target | Identity eigenvectors, D divergence cleaning. |
| `IERL_B` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | target | Identity eigenvectors, B divergence cleaning. |
| `IERL_D_B` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | target | Identity eigenvectors, D & B divergence cleaning. |
