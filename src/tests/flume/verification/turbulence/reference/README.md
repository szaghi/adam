# Reference data of the turbulence verification (issue #119)

`spectral_Re1600_512.gdiag`: the 512^3 pseudo-spectral DNS of the Taylor-Green vortex at Re 1600, incompressible, of the High-Order CFD Workshop
(case C3.5). Columns: time, kinetic energy, its dissipation rate -dE/dt, enstrophy, all per unit volume. Supplied by the owner on 2026-11-02, as
received except the line ends, CRLF converted to LF (the values are untouched). Read by `tgv_dns.py` and `make_doc_figures.py` (`fig_tgv_iles`, `fig_tgv_models`, `fig_tgv_hybrid`).
