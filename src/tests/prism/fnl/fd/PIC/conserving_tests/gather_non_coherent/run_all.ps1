# Launch the FNL conserving tests of this directory from PowerShell through WSL (sequential run, see run_all.sh).
# Usage:  .\run_all.ps1                     -> all the tests of this directory
#         .\run_all.ps1 totorica londrillo  -> only the listed tests
wsl.exe -e bash /home/federico/fortran/adam/src/tests/prism/fnl/fd/PIC/conserving_tests/gather_non_coherent/run_all.sh @args
