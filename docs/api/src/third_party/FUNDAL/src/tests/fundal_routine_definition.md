---
title: fundal_routine_definition
---

# fundal_routine_definition

> FUNDAL, use device memory in external routine test.
 External routine module definition.

**Source**: `src/third_party/FUNDAL/src/tests/fundal_external_routine_test.F90`

**Dependencies**

```mermaid
graph LR
  fundal_routine_definition["fundal_routine_definition"] --> iso_fortran_env["iso_fortran_env"]
```

## Contents

- [do_work_ptr](#do-work-ptr)
- [do_work_unstr](#do-work-unstr)

## Subroutines

### do_work_ptr

Do work on device memory, pointer approach.

```fortran
subroutine do_work_ptr(n, a)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `n` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Array dimension. |
| `a` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | inout | target | Array. |

### do_work_unstr

Do work on device memory, unstructured approach.

```fortran
subroutine do_work_unstr(n, a)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `n` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Array dimension. |
| `a` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | inout | target | Array. |
