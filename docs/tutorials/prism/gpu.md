# 4. The same on the GPU

The pulse of the first page on PRISM's GPU backend (FNL: OpenACC with nvfortran): the same input, the same outputs,
and the same fields to round-off.

## Build

The FNL modes need nvfortran on `PATH` and a varset that names your GPU's compute capability and HDF5 build
([install](/start/install), [build modes](/reference/build-modes)):

```bash
fobis build --mode prism-fnl-nvf --varset local_nvf    # exe/adam_prism_fnl; local_nvf: the cc89 workstation
```

## Run

The input of [page 1](./pulse), unchanged, one GPU per rank.

**On a WSL2 workstation, set the [GPU5](/cookbook/gpu#gpu5-run-on-a-development-box-with-wsl2-not-on-a-cluster) settings before the command below, or it aborts at the first ghost
exchange.**

```bash
mkdir -p run/pulse-fnl && cp src/tests/prism/regression/periodic-pulse-blocks/input.ini run/pulse-fnl/
cd run/pulse-fnl && mpirun -np 2 ../../exe/adam_prism_fnl input.ini > log.txt
```

The state lives on the device for the whole run and comes back to the host only for the outputs; the files are the
CPU run's, in the same format.

## Compare with the CPU run

The last rows of the two energy histories:

<<< @/snippets/tutorials/prism-gpu-energy.txt

and the fields of step 192, block by block (`compare_saves.py` with the FNL directory as the second run):

```bash
../../exe/.regression-venv/bin/python ../../src/tests/tutorials/compare_saves.py ../pulse periodic_pulse_blocks 192 192 .
```

<<< @/snippets/tutorials/prism-gpu-compare.txt

The two backends run the same algorithm; they differ only in the order of floating-point operations, so they agree
to round-off, not bitwise.

::: tip What you learned
- The GPU build is another build mode of the same sources; the input and the outputs do not change.
- CPU and GPU agree to round-off; compare them field by field (`compare_saves.py` prints every field of the save), not only by integrals.
- WSL2 checks correctness; performance is measured on a cluster.
:::

Previous: [3. A perfectly matched layer](./pml). This is the end of the PRISM track; the PRISM pages
([models](/applications/prism/models), [numerics](/applications/prism/numerics), [PML](/applications/prism/pml),
[PIC](/applications/prism/pic)) go further.
