---
title: fundal_laplace_dev_routine_external
---

# fundal_laplace_dev_routine_external

> External laplace routine module.

**Source**: `src/third_party/FUNDAL/src/tests/laplace/fundal_laplace_dev_routine.F90`

**Dependencies**

```mermaid
graph LR
  fundal_laplace_dev_routine_external["fundal_laplace_dev_routine_external"] --> iso_fortran_env["iso_fortran_env"]
```

## Contents

- [laplace](#laplace)

## Subroutines

### laplace

```fortran
subroutine laplace(n, m, A, Anew, error)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `n` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  |  |
| `m` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  |  |
| `A` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | inout | target |  |
| `Anew` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | inout | target |  |
| `error` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | inout |  |  |
