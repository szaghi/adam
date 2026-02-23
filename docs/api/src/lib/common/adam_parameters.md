---
title: adam_parameters
---

# adam_parameters

> ADAM, general parameters.

**Source**: `src/lib/common/adam_parameters.f90`

**Dependencies**

```mermaid
graph LR
  adam_parameters["adam_parameters"] --> penf["penf"]
```

## Variables

| Name | Type | Attributes | Description |
|------|------|------------|-------------|
| `BC_PERIODIC` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | parameter | Flag (reserved) for periodical boundary conditions. |
| `TO_BE_REFINED` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | parameter | Flag for node/block to be refined. |
| `TO_BE_DEREFINED` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | parameter | Flag for node/block to be derefined. |
| `TO_NOT_TOUCH` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | parameter | Flag for node/block to be untouched. |
| `FEC_1_6_ARRAY` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | parameter | Mapping fec1-26 to fec1-6 for boundaries. |
| `FEC_TO_DELTA` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | parameter | Neighor map. |
| `DELTA_TO_FEC` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | parameter | Delta to fec map. |
