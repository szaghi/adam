# ASCOT

**ASCOT** — **A**DAM **S**li**c**es **C**onverter — is a minimal standalone utility that converts ADAM slice binary files into Tecplot-compatible ASCII format for post-processing and visualisation.

## Source Layout

```
src/app/ascot/
└── ascot.F90    # Single-file program (no modules)
```

ASCOT has no common or backend subdirectories. The entire tool is a single self-contained Fortran program.

## File Formats

### Input — ADAM Slice Binary (`.mat`)

Fortran unformatted stream binary with the following layout:

| Section | Content | Fortran type |
|---------|---------|-------------|
| Header | `[ni, nj, nk, nv]` — grid dimensions and variable count | `integer(I4P)(4)` |
| Data | `nv` values per cell, looped in k → j → i order | `real(R8P)(nv)` per cell |

The data loop is equivalent to:

```fortran
do k = 1, nk
  do j = 1, nj
    do i = 1, ni
      read(file_unit_input) vars(1:nv)
    enddo
  enddo
enddo
```

### Output — Tecplot ASCII (`.dat`)

Standard Tecplot point-data format:

```
VARIABLES = "x" "y" "Z" "v4" "v5" ...
ZONE T="ADAM slice", I=ni, J=nj, K=nk
v1 v2 v3 ... vnv
v1 v2 v3 ... vnv
...
```

One data line is written per cell, with values separated by spaces.

## Command-Line Interface

```
ascot -i <input_file> [-v "<variable names>"] [-o <output_file>]
```

At minimum, `-i` and its argument must be supplied (the program prints usage and stops if fewer than two arguments are given).

| Flag | Argument | Default | Description |
|------|----------|---------|-------------|
| `-i` | `<file>` | — | Input ADAM slice binary file (**required**) |
| `-o` | `<file>` | `slice.dat` | Output Tecplot ASCII file |
| `-v` | `"<names>"` | see below | Space-separated variable name string for the Tecplot `VARIABLES` line |

### Default Variable Names

When `-v` is omitted, ASCOT generates a default `VARIABLES` header that assumes the first three stored values are spatial coordinates and names subsequent variables sequentially:

```
VARIABLES = "x" "y" "Z" "v4" "v5" ...
```

For a file with `nv` variables, names `v4` through `vnv` are appended automatically.

### Usage Examples

```bash
# Full specification: custom variable names and output file
ascot -i slice.mat -v "rho rhou rhov rhow rhoe" -o slice.dat

# Custom output file, default variable names
ascot -i slice.mat -o slice.dat

# Minimum: default output file (slice.dat) and default variable names
ascot -i slice.mat
```

## Building

```bash
FoBiS.py build -mode ascot-nvf-cuda
```

The executable is written to `exe/ascot`.

## License

ASCOT is part of the ADAM framework, released under the [GNU Lesser General Public License v3.0](https://www.gnu.org/licenses/lgpl-3.0.html) (LGPLv3).

> Copyright (C) Andrea Di Mascio, Federico Negro, Giacomo Rossi, Francesco Salvadore, Stefano Zaghi.
