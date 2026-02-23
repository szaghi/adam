---
title: finer_option_t
---

# finer_option_t

> Option class definition.

**Source**: `src/third_party/FiNeR/src/lib/finer_option_t.F90`

**Dependencies**

```mermaid
graph LR
  finer_option_t["finer_option_t"] --> finer_backend["finer_backend"]
  finer_option_t["finer_option_t"] --> penf["penf"]
  finer_option_t["finer_option_t"] --> stringifor["stringifor"]
```

## Contents

- [option](#option)
- [option](#option)
- [free](#free)
- [get_pairs](#get-pairs)
- [parse](#parse)
- [get_option](#get-option)
- [get_a_option](#get-a-option)
- [parse_comment](#parse-comment)
- [parse_name](#parse-name)
- [parse_value](#parse-value)
- [print_option](#print-option)
- [set_option](#set-option)
- [set_a_option](#set-a-option)
- [save_option](#save-option)
- [assign_option](#assign-option)
- [count_values](#count-values)
- [name_len](#name-len)
- [values_len](#values-len)
- [option_eq_string](#option-eq-string)
- [option_eq_character](#option-eq-character)
- [new_option](#new-option)

## Derived Types

### option

Option data of sections.

#### Components

| Name | Type | Attributes | Description |
|------|------|------------|-------------|
| `oname` | type([string](/api/src/third_party/StringiFor/src/lib/stringifor_string_t#string)) |  | Option name. |
| `ovals` | type([string](/api/src/third_party/StringiFor/src/lib/stringifor_string_t#string)) |  | Option values. |
| `ocomm` | type([string](/api/src/third_party/StringiFor/src/lib/stringifor_string_t#string)) |  | Eventual option inline comment. |

#### Type-Bound Procedures

| Name | Attributes | Description |
|------|------------|-------------|
| `count_values` | pass(self) | Counting option value(s). |
| `free` | pass(self) | Free dynamic memory. |
| `get` |  | Get option value (scalar). |
| `get_pairs` | pass(self) | Return option name/values pairs. |
| `name_len` | pass(self) | Return option name length. |
| `parse` | pass(self) | Parse option data. |
| `print` | pass(self) | Pretty print data. |
| `save` | pass(self) | Save data. |
| `set` |  | Set option value (scalar). |
| `values_len` | pass(self) | Return option values length. |
| `assignment(=)` |  | Assignment overloading. |
| `operator(==)` |  | Equal operator overloading. |
| `get_option` | pass(self) | Get option value (scalar). |
| `get_a_option` | pass(self) | Get option value (array). |
| `parse_comment` | pass(self) | Parse option inline comment. |
| `parse_name` | pass(self) | Parse option name. |
| `parse_value` | pass(self) | Parse option values. |
| `set_option` | pass(self) | Set option value (scalar). |
| `set_a_option` | pass(self) | Set option value (array). |
| `assign_option` | pass(lhs) | Assignment overloading. |
| `option_eq_string` | pass(lhs) | Equal to string logical operator. |
| `option_eq_character` | pass(lhs) | Equal to character logical operator. |

## Interfaces

### option

Overload `option` name with a function returning a new (initiliazed) option instance.

**Module procedures**: [`new_option`](/api/src/third_party/FiNeR/src/lib/finer_option_t#new-option)

## Subroutines

### free

Free dynamic memory.

**Attributes**: elemental

```fortran
subroutine free(self)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `self` | class([option](/api/src/third_party/FiNeR/src/lib/finer_option_t#option)) | inout |  | Option data. |

### get_pairs

Return option name/values pairs.

**Attributes**: pure

```fortran
subroutine get_pairs(self, pairs)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `self` | class([option](/api/src/third_party/FiNeR/src/lib/finer_option_t#option)) | in |  | Option data. |
| `pairs` | character(len=:) | out | allocatable | Option name/values pairs. |

**Call graph**

```mermaid
flowchart TD
  loop["loop"] --> get_pairs["get_pairs"]
  option_pairs["option_pairs"] --> get_pairs["get_pairs"]
  get_pairs["get_pairs"] --> chars["chars"]
  get_pairs["get_pairs"] --> is_allocated["is_allocated"]
  style get_pairs fill:#3e63dd,stroke:#99b,stroke-width:2px
```

### parse

Parse option data from a source string.

**Attributes**: elemental

```fortran
subroutine parse(self, sep, source, error)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `self` | class([option](/api/src/third_party/FiNeR/src/lib/finer_option_t#option)) | inout |  | Option data. |
| `sep` | character(len=*) | in |  | Separator of option name/value. |
| `source` | type([string](/api/src/third_party/StringiFor/src/lib/stringifor_string_t#string)) | inout |  | String containing option data. |
| `error` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | out |  | Error code. |

**Call graph**

```mermaid
flowchart TD
  load["load"] --> parse["parse"]
  parse["parse"] --> parse["parse"]
  parse_file["parse_file"] --> parse["parse"]
  parse_options["parse_options"] --> parse["parse"]
  search["search"] --> parse["parse"]
  parse["parse"] --> parse_comment["parse_comment"]
  parse["parse"] --> parse_name["parse_name"]
  parse["parse"] --> parse_value["parse_value"]
  style parse fill:#3e63dd,stroke:#99b,stroke-width:2px
```

### get_option

for getting option data value (scalar).

```fortran
subroutine get_option(self, val, error)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `self` | class([option](/api/src/third_party/FiNeR/src/lib/finer_option_t#option)) | in |  | Option data. |
| `val` | class(*) | inout |  | Value. |
| `error` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | out | optional | Error code. |

**Call graph**

```mermaid
flowchart TD
  get_option["get_option"] --> chars["chars"]
  get_option["get_option"] --> is_allocated["is_allocated"]
  get_option["get_option"] --> to_number["to_number"]
  style get_option fill:#3e63dd,stroke:#99b,stroke-width:2px
```

### get_a_option

Get option data values (array).

```fortran
subroutine get_a_option(self, val, delimiter, error)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `self` | class([option](/api/src/third_party/FiNeR/src/lib/finer_option_t#option)) | in |  | Option data. |
| `val` | class(*) | inout |  | Value. |
| `delimiter` | character(len=*) | in | optional | Delimiter used for separating values. |
| `error` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | out | optional | Error code. |

**Call graph**

```mermaid
flowchart TD
  get_a_option["get_a_option"] --> chars["chars"]
  get_a_option["get_a_option"] --> is_allocated["is_allocated"]
  get_a_option["get_a_option"] --> split["split"]
  get_a_option["get_a_option"] --> to_number["to_number"]
  style get_a_option fill:#3e63dd,stroke:#99b,stroke-width:2px
```

### parse_comment

Parse option inline comment trimming it out from pure value string.

**Attributes**: elemental

```fortran
subroutine parse_comment(self)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `self` | class([option](/api/src/third_party/FiNeR/src/lib/finer_option_t#option)) | inout |  | Option data. |

**Call graph**

```mermaid
flowchart TD
  parse["parse"] --> parse_comment["parse_comment"]
  parse_comment["parse_comment"] --> is_allocated["is_allocated"]
  parse_comment["parse_comment"] --> slice["slice"]
  style parse_comment fill:#3e63dd,stroke:#99b,stroke-width:2px
```

### parse_name

Parse option name from a source string.

**Attributes**: elemental

```fortran
subroutine parse_name(self, sep, source, error)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `self` | class([option](/api/src/third_party/FiNeR/src/lib/finer_option_t#option)) | inout |  | Option data. |
| `sep` | character(len=*) | in |  | Separator of option name/value. |
| `source` | type([string](/api/src/third_party/StringiFor/src/lib/stringifor_string_t#string)) | in |  | String containing option data. |
| `error` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | out |  | Error code. |

**Call graph**

```mermaid
flowchart TD
  parse["parse"] --> parse_name["parse_name"]
  parse["parse"] --> parse_name["parse_name"]
  parse_name["parse_name"] --> slice["slice"]
  style parse_name fill:#3e63dd,stroke:#99b,stroke-width:2px
```

### parse_value

Parse option value from a source string.

**Attributes**: elemental

```fortran
subroutine parse_value(self, sep, source, error)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `self` | class([option](/api/src/third_party/FiNeR/src/lib/finer_option_t#option)) | inout |  | Option data. |
| `sep` | character(len=*) | in |  | Separator of option name/value. |
| `source` | type([string](/api/src/third_party/StringiFor/src/lib/stringifor_string_t#string)) | in |  | String containing option data. |
| `error` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | out |  | Error code. |

**Call graph**

```mermaid
flowchart TD
  parse["parse"] --> parse_value["parse_value"]
  parse_value["parse_value"] --> slice["slice"]
  style parse_value fill:#3e63dd,stroke:#99b,stroke-width:2px
```

### print_option

Print data with a pretty format.

```fortran
subroutine print_option(self, unit, retain_comments, pref, iostat, iomsg)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `self` | class([option](/api/src/third_party/FiNeR/src/lib/finer_option_t#option)) | in |  | Option data. |
| `unit` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Logic unit. |
| `retain_comments` | logical | in |  | Flag for retaining eventual comments. |
| `pref` | character(len=*) | in | optional | Prefixing string. |
| `iostat` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | out | optional | IO error. |
| `iomsg` | character(len=*) | out | optional | IO error message. |

**Call graph**

```mermaid
flowchart TD
  print_option["print_option"] --> is_allocated["is_allocated"]
  style print_option fill:#3e63dd,stroke:#99b,stroke-width:2px
```

### set_option

Set option data value (scalar).

**Attributes**: pure

```fortran
subroutine set_option(self, val)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `self` | class([option](/api/src/third_party/FiNeR/src/lib/finer_option_t#option)) | inout |  | Option data. |
| `val` | class(*) | in |  | Value. |

**Call graph**

```mermaid
flowchart TD
  set_option["set_option"] --> str["str"]
  style set_option fill:#3e63dd,stroke:#99b,stroke-width:2px
```

### set_a_option

Set option data value (array).

**Attributes**: pure

```fortran
subroutine set_a_option(self, val, delimiter)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `self` | class([option](/api/src/third_party/FiNeR/src/lib/finer_option_t#option)) | inout |  | Option data. |
| `val` | class(*) | in |  | Value. |
| `delimiter` | character(len=*) | in | optional | Delimiter used for separating values. |

**Call graph**

```mermaid
flowchart TD
  set_a_option["set_a_option"] --> str["str"]
  set_a_option["set_a_option"] --> strip["strip"]
  style set_a_option fill:#3e63dd,stroke:#99b,stroke-width:2px
```

### save_option

Save data.

```fortran
subroutine save_option(self, unit, retain_comments, iostat, iomsg)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `self` | class([option](/api/src/third_party/FiNeR/src/lib/finer_option_t#option)) | in |  | Option data. |
| `unit` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Logic unit. |
| `retain_comments` | logical | in |  | Flag for retaining eventual comments. |
| `iostat` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | out | optional | IO error. |
| `iomsg` | character(len=*) | out | optional | IO error message. |

**Call graph**

```mermaid
flowchart TD
  save_option["save_option"] --> is_allocated["is_allocated"]
  style save_option fill:#3e63dd,stroke:#99b,stroke-width:2px
```

### assign_option

Assignment between two options.

**Attributes**: elemental

```fortran
subroutine assign_option(lhs, rhs)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `lhs` | class([option](/api/src/third_party/FiNeR/src/lib/finer_option_t#option)) | inout |  | Left hand side. |
| `rhs` | type([option](/api/src/third_party/FiNeR/src/lib/finer_option_t#option)) | in |  | Rigth hand side. |

**Call graph**

```mermaid
flowchart TD
  assign_option["assign_option"] --> is_allocated["is_allocated"]
  style assign_option fill:#3e63dd,stroke:#99b,stroke-width:2px
```

## Functions

### count_values

Get the number of values of option data.

**Attributes**: elemental

**Returns**: integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables))

```fortran
function count_values(self, delimiter) result(Nv)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `self` | class([option](/api/src/third_party/FiNeR/src/lib/finer_option_t#option)) | in |  | Option data. |
| `delimiter` | character(len=*) | in | optional | Delimiter used for separating values. |

**Call graph**

```mermaid
flowchart TD
  count_values["count_values"] --> count_values["count_values"]
  count_values["count_values"] --> count_values["count_values"]
  file_ini_autotest["file_ini_autotest"] --> count_values["count_values"]
  count_values["count_values"] --> is_allocated["is_allocated"]
  style count_values fill:#3e63dd,stroke:#99b,stroke-width:2px
```

### name_len

Return option name length.

**Attributes**: elemental

**Returns**: `integer`

```fortran
function name_len(self) result(length)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `self` | class([option](/api/src/third_party/FiNeR/src/lib/finer_option_t#option)) | in |  | Option data. |

**Call graph**

```mermaid
flowchart TD
  max_chars_len["max_chars_len"] --> name_len["name_len"]
  name_len["name_len"] --> is_allocated["is_allocated"]
  style name_len fill:#3e63dd,stroke:#99b,stroke-width:2px
```

### values_len

Return option values length.

**Attributes**: elemental

**Returns**: `integer`

```fortran
function values_len(self) result(length)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `self` | class([option](/api/src/third_party/FiNeR/src/lib/finer_option_t#option)) | in |  | Option data. |

**Call graph**

```mermaid
flowchart TD
  max_chars_len["max_chars_len"] --> values_len["values_len"]
  values_len["values_len"] --> is_allocated["is_allocated"]
  style values_len fill:#3e63dd,stroke:#99b,stroke-width:2px
```

### option_eq_string

Equal to string logical operator.

**Attributes**: elemental

**Returns**: `logical`

```fortran
function option_eq_string(lhs, rhs) result(is_it)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `lhs` | class([option](/api/src/third_party/FiNeR/src/lib/finer_option_t#option)) | in |  | Left hand side. |
| `rhs` | type([string](/api/src/third_party/StringiFor/src/lib/stringifor_string_t#string)) | in |  | Right hand side. |

### option_eq_character

Equal to character logical operator.

**Attributes**: elemental

**Returns**: `logical`

```fortran
function option_eq_character(lhs, rhs) result(is_it)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `lhs` | class([option](/api/src/third_party/FiNeR/src/lib/finer_option_t#option)) | in |  | Left hand side. |
| `rhs` | character(kind=[CK](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables), len=*) | in |  | Right hand side. |

### new_option

Return a new (initiliazed) option instance.

**Attributes**: elemental

**Returns**: type([option](/api/src/third_party/FiNeR/src/lib/finer_option_t#option))

```fortran
function new_option(option_name, option_values, option_comment)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `option_name` | character(len=*) | in | optional | Option name. |
| `option_values` | character(len=*) | in | optional | Option values. |
| `option_comment` | character(len=*) | in | optional | Option comment. |
