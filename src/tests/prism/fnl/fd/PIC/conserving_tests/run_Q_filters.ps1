# Launch the FNL TSC esirkepov_Q<N> conserving tests from PowerShell through WSL (sequential run, see run_Q_filters.sh).
# Usage:  .\run_Q_filters.ps1          -> Q5, Q6, Q8
#         .\run_Q_filters.ps1 Q5 Q8    -> only the listed tests
wsl.exe -e bash /home/federico/fortran/adam/src/tests/prism/fnl/fd/PIC/conserving_tests/run_Q_filters.sh @args
