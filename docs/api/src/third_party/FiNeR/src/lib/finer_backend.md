---
title: finer_backend
---

# finer_backend

> FiNeR backends: globals definition.

**Source**: `src/third_party/FiNeR/src/lib/finer_backend.f90`

**Dependencies**

```mermaid
graph LR
  finer_backend["finer_backend"] --> penf["penf"]
```

## Variables

| Name | Type | Attributes | Description |
|------|------|------------|-------------|
| `ERR_OPTION_NAME` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | parameter | Error flag for trapping errors in option name. |
| `ERR_OPTION_VALS` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | parameter | Error flag for trapping errors in option values. |
| `ERR_OPTION` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | parameter | Error flag for trapping errors in option. |
| `ERR_SECTION_NAME` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | parameter | Error flag for trapping errors in section name. |
| `ERR_SECTION_OPTIONS` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | parameter | Error flag for trapping errors in section options. |
| `ERR_SECTION` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | parameter | Error flag for trapping errors in section. |
| `ERR_SOURCE_MISSING` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | parameter | Error flag for trapping errors in file when source is missing. |
| `DEF_OPT_SEP` | character(len=1) | parameter | Default separator of option name/value. |
| `COMMENTS` | character(len=*) | parameter | Characters used for defining a comment line. |
| `INLINE_COMMENT` | character(len=1) | parameter | Inline comment delimiter. |
