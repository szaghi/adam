<a name="top"></a>

# ADAM library NVF sources

> This subdirectory contains the ADAM library kernels of CUDAFortran NVF GPU-offloading backend.

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

+ `adam_base_nvf_object.F90` is the **base backend** class for NVF acceleration, see [BASE NVF object API](https://szaghi.github.io/adam/type/base_nvf_object.html) for more details;
+ `adam_ib_nvf_kernels.F90` is an **kernels library** for handling IB computation accelerated in NVF backend, see [memory library API](https://szaghi.github.io/adam/module/adam_ib_nvf_kernels.html) for more details;
+ `adam_memory_nvf_lib.F90` is an **memory library** for handling memory (GPU) allocation in NVF backend, see [memory library API](https://szaghi.github.io/adam/module/adam_memory_nvf_lib.html) for more details.

Go to [Top](#top)

