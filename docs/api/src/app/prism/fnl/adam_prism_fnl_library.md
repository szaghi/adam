---
title: adam_prism_fnl_library
---

# adam_prism_fnl_library

> ADAM PRISM FNL library, entry for all PRISM FNL (and common) classes and libraries.

**Source**: `src/app/prism/fnl/adam_prism_fnl_library.F90`

**Dependencies**

```mermaid
graph LR
  adam_prism_fnl_library["adam_prism_fnl_library"] --> adam_fnl_library["adam_fnl_library"]
  adam_prism_fnl_library["adam_prism_fnl_library"] --> adam_prism_common_library["adam_prism_common_library"]
  adam_prism_fnl_library["adam_prism_fnl_library"] --> adam_prism_fnl_coil_object["adam_prism_fnl_coil_object"]
  adam_prism_fnl_library["adam_prism_fnl_library"] --> adam_prism_fnl_external_fields_kernels["adam_prism_fnl_external_fields_kernels"]
  adam_prism_fnl_library["adam_prism_fnl_library"] --> adam_prism_fnl_fWLayer_object["adam_prism_fnl_fWLayer_object"]
```
