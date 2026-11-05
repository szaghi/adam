# Reference data of the turbulence verification (issue #119)

`spectral_Re1600_512.gdiag`: the 512^3 pseudo-spectral DNS of the Taylor-Green vortex at Re 1600, incompressible, of the High-Order CFD Workshop
(case C3.5). Columns: time, kinetic energy, its dissipation rate -dE/dt, enstrophy, all per unit volume. Supplied by the owner on 2026-11-02, as
received except the line ends, CRLF converted to LF (the values are untouched). Read by `tgv_dns.py` and `make_doc_figures.py` (`fig_tgv_iles`, `fig_tgv_models`, `fig_tgv_hybrid`).

`M2p5_Stat.txt`, `M6Tw076_Stat.txt`: wall-normal profiles at one station of the compressible zero-pressure-gradient boundary-layer DNS of Zhang, Duan &
Choudhari (2018, AIAA J 56(11), doi:10.2514/1.J057296), cases Mach 2.5 adiabatic (T_w/T_r = 1, Re_tau 510) and Mach 5.86 (T_w/T_r = 0.76, Re_tau 453), as
downloaded on 2026-11-05 from the NASA Turbulence Modeling Resource (`tmbwg.github.io/turbmodels/Other_DNS_Data/supersonic_hypersonic_flatplate.html`, there `<case>_Stat.dat`),
untouched but the extension (`.dat` is ignored under src/tests). `tbl_m2p5_profile.txt` is the first converted by `tbl_oracle.py profile` into the inflow table of the digital-filter face (issue #119, T3.7).
