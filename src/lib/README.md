<a name="top"></a>

# ADAM library

> ADAM library, sources of ADAM library, part of ADAM framework.

This subdirectory contains the sources of ADAM library, e.g. the AMR library on top of which the ADAM framework is built.

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

These sources are organized in 3 sudirectories:

+ [`common`](./common/README.md): this subdirectory contains the ADAM library common sources, namely the sources that are common to all backend;
+ [`gmp`](./gmp/README.md): this subdirectory contains the ADAM library kernels of OpenMP GPU-offloading backend;
+ [`nvf`](./nvf/README.md): this subdirectory contains the ADAM library kernels of CUDAFortran NVF GPU-offloading backend.

Refer to the documentation contained into each sources subdirectory for more details.

Go to [Top](#top)
