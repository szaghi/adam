---
title: adam_chase_parameters
---

# adam_chase_parameters

> CHASE, general parameters.

**Source**: `src/app/chase/common/adam_chase_parameters.F90`

**Dependencies**

```mermaid
graph LR
  adam_chase_parameters["adam_chase_parameters"] --> penf["penf"]
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
| `EV` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | target | Eigenvalues. |
| `ER` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | target | Right eigenvectors. |
| `EL` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | target | Right eigenvectors. |
| `IEV` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | target | Identity eigenvalues. |
| `IERL` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | target | Identity eigenvectors. |
