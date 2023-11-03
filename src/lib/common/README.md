<a name="top"></a>

# ADAM library common sources

> This subdirectory contains the ADAM library common sources, namely the sources that are common to all backend.

The main documentation of these sources is contained in the following sections:

| [Copyrights](#copyrights) | [API Documentation](#api-documentation) |

Go to [Top](#top)

# Copyrights

ADAM is currently a closed project:

> Copyright (C) Di Mascio/Rossi/Salvadore/Zaghi, Inc - All Rights Reserved.
>
> Unauthorized copying of these source files, via any medium is strictly prohibited, proprietary and confidential.
> Written by Andrea di Mascio, Giacomo Rossi, Francesco Salvadore and Stefano Zaghi, September 2023.

Future versions could be released with a more Free Open Source Software (FOSS) licence.

Go to [Top](#top)

# API Documentation

Currently, the following sources compose the subdirectory:

+ `adam_adam_object.F90` is a **composer** for correctly handling TREE and FIELD classes, see [ADAM object API](https://szaghi.github.io/adam/type/adam_object.html) for more details;
+ `adam_amr_object.F90` is an **AMR handler** for AMR markers, see [AMR object API](https://szaghi.github.io/adam/type/amr_object.html) for more details;
+ `adam_field_object.F90` is the **FIELD** class, see [FIELD object API](https://szaghi.github.io/adam/type/field_object.html) for more details;
+ `adam_grid_object.F90` is the **GRID** class, see [GRID object API](https://szaghi.github.io/adam/type/grid_object.html) for more details;
+ `adam_ib_object.F90` is an **IB handler** for immersed boundary solids markers, see [IB object API](https://szaghi.github.io/adam/type/ib_object.html) for more details;
+ `adam_memory_library.F90` is an **memory library** for handling memory (CPU) allocation, see [memory library API](https://szaghi.github.io/adam/module/adam_memory_lib.html) for more details;
+ `adam_mpih_object.F90` is an **MPI handler** incorporated into all ADAM classes, see [MPIH object API](https://szaghi.github.io/adam/type/mpih_object.html) for more details;
+ `adam_parameters.f90` is a **global parameter module**, see [ADAM parameters module API](https://szaghi.github.io/adam/module/adam_parameters.html) for more details;
+ `adam_slices_object.F90` is an **slices handler** for easy create IO domain slices, see [SLICES object API](https://szaghi.github.io/adam/type/slices_object.html) for more details;
+ `adam_tree_bucket_object.f90` is a **bucket class** on top of which TREE class is built, see [TREE BUCKET object API](https://szaghi.github.io/adam/type/tree_bucket_object.html) for more details;
+ `adam_tree_node_object.f90` is a **node class** on top of which TREE class is built, see [TREE NODE object API](https://szaghi.github.io/adam/type/tree_node_object.html) for more details;
+ `adam_tree_object.F90` is the **TREE** class, see [TREE object API](https://szaghi.github.io/adam/type/tree_object.html) for more details;

Go to [Top](#top)
