# Compressible boundary layer

| ID | Milestone | Type | Legs | Backends | Figure |
|---|---|---|---|---|---|
| TB | [#119](https://github.com/szaghi/adam/issues/119) · T3.7, [#121](https://github.com/szaghi/adam/issues/121) · S5 | P+C | tb1–tb3 | CPU, FNL | yes |

## Case and scope

T3.7 validates the wall model on spatially developing zero-pressure-gradient boundary layers: the DNS of Zhang, Duan &
Choudhari (2018, AIAA J 56(11)) at Mach 2.5 (adiabatic wall) and Mach 5.86 (cold wall). The station profiles come from the
NASA Turbulence Modeling Resource (`turbulence/reference/`). The layer is fed through the face `inflow-digital-filter`
([turbulence](../../turbulence#digital-filter-inflow)): the DNS station profile scaled to the inflow thickness
$\delta_i$, with fluctuations that have the DNS Reynolds stresses. $C_f$ is the wall model's own $\tau_w$, written by
`[diagnostics] wall_history_save` (leg tb2). The legs check the inflow generator and the wall history; the campaign
measures $C_f$, $B_q$, the shape factor and the profiles against the DNS.

## Legs

`src/tests/flume/verification/turbulence/check.sh`, CPU and FNL.

| Leg | What it verifies | Measured |
|---|---|---|
| tb1 | digital-filter inflow: the generator's stresses, temporal and spanwise correlations, periodic span, plane restart (`test_inflow_df`); a Mach 2.5 box: the inflow ghosts carry the table; one rank against 2; restart with the `.dff` against straight, control without it | stresses 0.4–0.7 % (bound 3 %), correlation exact, integral length 5 % short (bound 10 %); $\langle u\rangle$ at $\delta/2$ 0.877 against 0.891; $2\cdot10^{-16}$; fields and histories bit for bit (CPU and FNL), the control differs; FNL against CPU $10^{-15}$ |
| tb2 | wall-model history: $\tau_w$ of wm5's laminar channel against the driving force; one rank against 2 on the Mach 2.5 box; CPU against GPU | $3.9\cdot10^{-4}$ (bound $10^{-3}$); $4\cdot10^{-16}$; $1.1\cdot10^{-13}$ |
| tb3 | recycling-rescaling inflow ([turbulence](../../turbulence#recycling-rescaling-inflow)): blend, damping, plane, zero spanwise means and the rescaling on a synthetic block (`test_inflow_rr`); the Mach 2.5 box, 120 digital-filter steps then the recycling face (station 2.5, ratios 1): the ghosts rebuilt independently by `tbl_oracle.py` from the saved state on one rank and 2; restart at step 122 against straight; a rescaled run's inflow mean and fluctuations | to round-off; $6.7\cdot10^{-16}$ over 768 ghosts on both, one rank against 2 $1.6\cdot10^{-16}$; bit for bit (CPU and FNL); $\langle u\rangle$ at $\delta/2$ 0.888 against 0.891, $v_{rms}$, $w_{rms}$ 0.010, 0.011 (floor $10^{-3}$); FNL against CPU $2.7\cdot10^{-15}$ |

## Mach 2.5

Case M2p5 of the DNS: Mach 2.5, adiabatic wall ($T_w/T_r=1$), $Re_\theta=2835$ and $Re_\tau=510$ at its station.

| Setup | |
|---|---|
| box, cells | $90\times6\times6\,\delta_i$, $1080\times72\times72$, cubic $\delta_i/12$; periodic span, extrapolated top and outflow |
| gas | $Re=U_\infty\delta_i/\nu_\infty=19370$ (the DNS's 4 mm $\delta_i$), $\mu\propto T^{0.76}$ (the DNS uses Sutherland), $Pr=0.71$ |
| wall | `wall-modelled`, $T_w=2.104\,T_\infty$, the 3rd cell, semi-local damping, its eddy viscosity below $h$ (the setting that matched the [Mach 3 channel](./channel#mach-3)) |
| inflow | integral lengths $0.8$, $0.3$, $0.3\,\delta_i$ (streamwise as a time scale $L_x/U_\infty$, wall-normal, spanwise) |
| scheme | `hybrid-kep-weno`, no subgrid model |
| run | `tbl_cases_t37_l90.txt`, to $t=360\,\delta_i/U_\infty$ (4 flow-throughs); the wall history every 20 steps and 34 checkpoints of $180<t\le360$, averaged over the span and in time (`tbl_oracle.py compare`); 0.40 s per step on 4 A30 |

![The Mach 2.5 boundary layer against the DNS](/flume/tbl-m2p5.png)

*Left: $C_f$ from the wall model's $\tau_w$ (blue) and from the momentum integral $2\,d\theta/dx$ (grey), and $Re_\theta$ (right axis) against
$x$; the square is the DNS at its $Re_\theta$. Middle: $u^+_{VD}$ at the station where $Re_\theta$ is the DNS's, scaled with the run's $u_\tau$. Right:
the rms over $u_\tau$ at the same station; dashed: the DNS.*

| At $Re_\theta=2836$ ($x=55\,\delta_i$) | LES | DNS |
|---|---|---|
| $C_f$, the wall model's $\tau_w$ | **2.21e-3** (the two halves of the window 2.19e-3, 2.22e-3) | 2.31e-3 |
| $C_f$, momentum integral | 1.56e-3 | |
| $H$, against the edge state ([the shape factor](#the-shape-factor-and-the-free-stream)) | 4.08 | 4.10 (its profile: 4.22) |
| $u^+_{VD}$ at $y^+=30$, $100$, $200$ | 13.1, 17.5, 19.8 | 13.6, 16.7, 18.8 |
| $u_{rms}$, $v_{rms}$, $w_{rms}$ over $u_\tau$ at $y/\delta=0.2$ | 2.06, 0.81, 1.22 | 1.49, 0.93, 1.12 |
| $-\langle u'v'\rangle/u_\tau^2$ at $y/\delta=0.2$, $0.5$ | 0.80, 0.50 | 0.62, 0.37 |

**What it shows.**
- **$C_f$ within 4.4 % of the DNS** at its $Re_\theta$, with the setting that matched the Mach 3 channel. The two halves
  of the averaging window differ by 1.1 %. Most of the 4.4 % is the inflow's recovery, not the plateau: farther
  downstream, at $Re_\theta=3000$–$3350$, $C_f$ is within 1.6–0.4 % of the DNS's value scaled as $Re_\theta^{-1/4}$,
  and at $\delta_i/16$ it is 0.2 % at the station ([resolution](#resolution-mach-2-5-at-delta-i-16)).
- **The inflow recovers by $x\approx50\,\delta_i$.** The layer reaches the DNS's $Re_\theta$ at $x=55$, against 53 in
  the DNS. $C_f$ falls from $2.8\cdot10^{-3}$ at the inlet to $1.8\cdot10^{-3}$ at $x\approx25$, the adaptation of the
  synthetic eddies, then rises to a plateau of $2.2\cdot10^{-3}$ from $x\approx50$. That matches the $40\,\delta_i$
  recovery length of a digital filter at this Mach number (Huang, Duan & Choudhari 2022).
- **The momentum integral does not measure $C_f$ here.** Fitted over $\pm3\,\delta_i$, $2\,d\theta/dx$ swings between
  $0.9$ and $4\cdot10^{-3}$ along $x$ on the same averaged data, and is $-33$ % at the station.
- **The log layer sits 0.9–1.0 above the DNS in $u^+_{VD}$**, and $u_{rms}$ peaks too high, the near-wall streaks of a
  coarse LES. The resolved $-\langle u'v'\rangle$ in the outer layer is above the DNS's. At $\delta_i/16$ the log layer,
  $u_{rms}$ and $-\langle u'v'\rangle$ move toward the DNS.
- **$H$ is the DNS's within 1 %** (4.08), taken against the local free stream as the DNS takes it.

## The flow

![Instantaneous fields of the two boundary layers](/flume/tbl-fields.png)

*The last checkpoint of the Mach 2.5 run (top block, $90\,\delta_i$ box, $t=360$) and of the [Mach 5.86](#mach-5-86-cold-wall) cold-wall run (bottom block). Per run: the
streamwise velocity in the $x$-$y$ plane at mid-span over the whole box, with the 99 % line of the spanwise average (black); the wall-parallel plane
at $y=0.29\,\delta_i$ (its own colour range); a $20\,\delta_i$ window of the $x$-$y$ plane and the temperature in the cross-flow plane at its centre.
The figure is `make_doc_figures.py --only tbl-fields`.*

![The Mach 2.5 boundary layer and the interaction in 3-D](/flume/tbl-3d.png)

*The same planes as the faces of a cutaway, $u/U_\infty$, $y$ drawn stretched twice: the Mach 2.5 layer over $15<x<45$, the cross-flow face at
$x=45$ (top), and the Mach 2.28 [shock/boundary-layer interaction](../streams/sbli) over $35<x<60$, the cross-flow face at
the impingement $x=60$ (bottom). `make_doc_figures.py --only tbl-3d`.*

- **The streaks at $y=0.29\,\delta_i$ ($0.18\,\delta$) are outer-scaled.** The first minimum of the spanwise correlation
  of $u$ puts them $1.5$–$1.7\,\delta_i$ apart at Mach 2.5 and $1.8$–$2.0\,\delta_i$ at Mach 5.86, about one $\delta$,
  the same near the inflow ($x<15$) and downstream ($30<x<80$). At $\delta_i/12$ the first cell is $0.05\,\delta$, about
  26 wall units at $Re_\tau=510$; the near-wall cycle is the wall model's, as in the channel.
- **The edge of the layer is corrugated by bulges** of about $\delta$ in length. The cold-wall layer's cross-section
  shows the hot fluid near the wall (up to $5.5\,T_\infty$) lifted in ejections; at Mach 2.5 the near-wall temperature is
  the adiabatic $2.05\,T_\infty$.

## Resolution: Mach 2.5 at $\delta_i/16$ {#resolution-mach-2-5-at-delta-i-16}

| Setup | |
|---|---|
| grid | the same $90\,\delta_i$ box at $\delta_i/16$: $1440\times96\times96$ cells |
| wall | the matching point still the 3rd cell, so $h$ shrinks with the grid, as in Kawai & Larsson (2012) |
| inflow | the filter's node spacing the new $\Delta y$; everything else as at $\delta_i/12$ |
| run | `tbl_cases_t37_r16.txt`, 35264 steps, 30 checkpoints of $180<t\le360$; 0.80 s per step on 4 A30 |

![The Mach 2.5 boundary layer at two resolutions against the DNS](/flume/tbl-m2p5-res.png)

*Left: $C_f$ from the wall model's $\tau_w$ at $\delta_i/12$ (dashed) and $\delta_i/16$ (solid); the square is the DNS at its station. Middle:
$u^+_{VD}$ at the station of each run. Right: the rms and $-\langle u'v'\rangle$ over $u_\tau$; dotted: the DNS.*

| | $\delta_i/12$ | $\delta_i/16$ | DNS |
|---|---|---|---|
| $C_f$ at $Re_\theta=2835$ | 2.208e-3 (halves 2.195e-3, 2.220e-3) | **2.306e-3** (halves 2.299e-3, 2.312e-3) | 2.310e-3 |
| $C_f$ against the DNS scaled as $Re_\theta^{-1/4}$, at $Re_\theta=2600$, 2835, 3000, 3200, 3350 | −12.0, −3.9, −1.6, −0.9, −0.4 % | −3.6, −0.3, +1.0, +0.7, −0.5 % | |
| $H$, against the edge state | 4.08 | 4.23 | 4.10 (its profile: 4.22) |
| $u^+_{VD}$ at $y^+=30$, $100$, $200$ | 13.1, 17.5, 19.8 | 12.8, 17.3, 19.2 | 13.6, 16.7, 18.8 |
| $u_{rms}$, $v_{rms}$, $w_{rms}$ over $u_\tau$ at $y/\delta=0.2$ | 2.06, 0.81, 1.22 | 1.79, 0.83, 1.16 | 1.49, 0.93, 1.12 |
| $-\langle u'v'\rangle/u_\tau^2$ at $y/\delta=0.2$, $0.5$ | 0.80, 0.50 | 0.66, 0.43 | 0.62, 0.37 |

**What it shows.**
- **$C_f$ is grid-converged once the layer has recovered.** Downstream of $x\approx60$ ($Re_\theta\ge3000$) the two grids
  agree with each other and with the DNS's $C_f$ scaled as $Re_\theta^{-1/4}$ within 1–2 %. The scaling is an
  approximation (the DNS file holds one station), so these figures carry its error too.
- **The finer grid shortens the inflow's recovery.** The $C_f$ dip after the inlet ends by $x\approx40$ at $\delta_i/16$,
  against $x\approx55$–$60$ at $\delta_i/12$: the $-4.4$ % at the DNS station of the $\delta_i/12$ run is mostly the
  layer still recovering. At $\delta_i/16$ the station is past it, and $C_f$ is $-0.2$ %; the halves of the window differ
  by 0.6 %.
- **The near-wall and log-layer statistics converge toward the DNS**: the $u_{rms}$ peak at $y/\delta=0.2$ from 2.06 to
  1.79, $-\langle u'v'\rangle$ from 0.80 to 0.66, the log-layer offset in $u^+_{VD}$ at $y^+=200$ from +1.0 to +0.4.
- **$H$ is the DNS's on both grids** (4.08, 4.23, against 4.10) taken against the local free stream
  ([the shape factor](#the-shape-factor-and-the-free-stream)).

## The shape factor and the free stream

The DNS's profiles are normalised by its free stream at the station, and are 1 there. The run's free stream is not the
nominal one: the layer thickens downstream and turns the supersonic stream over it, a weak compression of order
$\gamma M^2 (d\delta^*/dx)/\sqrt{M^2-1}$ (linear supersonic theory). `tbl_oracle.py compare` therefore takes $\delta^*$ and
$\theta$ against the local edge state: $\rho_e$ and the Favre $u_e$ are the means over $1.5$–$2.5\,\delta$, the integrals
run to $1.5\,\delta$, and the DNS's own profile is integrated the same way (4.22 at Mach 2.5, against the header's rounded
4.1). Taken against the nominal $\rho_\infty$, $U_\infty$ and integrated to the top of the box, $H$ is 14 % low at Mach
2.5 and 30 % low at Mach 5.86.

![H against the edge state and against the nominal free stream](/flume/tbl-edge.png)

*Left and right: $H$ of the Mach 2.5 baseline and the Mach 5.86 run against $x$, against the edge state (solid) and the nominal free stream
(dashed); the DNS's header and its profile's at its station. Middle: the edge pressure and density over the nominal ones against $x$.*

**What it shows.**
- **The free stream is compressed.** At Mach 2.5 $p_e/p_\infty$ rises to 1.027 at $x\approx15$ and is 1.013 at the
  station, $\rho_e/\rho_\infty$ 1.009; at Mach 5.86 1.08 and 1.042. Integrated to the top of the box ($6\,\delta_i$, more
  than $3\,\delta$), the denser stream over the layer takes 0.5 off $H$ at Mach 2.5 and 4 at Mach 5.86. The dip of the
  nominal $H$ to 2.8 at $x\approx15$ is the compression's peak, not the layer.
- **Against the edge state $H$ is the DNS's**: 4.08 at Mach 2.5 (4.10, its profile 4.22), 13.4 at Mach 5.86 (13.6, its
  profile 14.9). Taken over $2.5$–$3.5$, $3$–$5$ or $4$–$5.5\,\delta_i$, the edge moves it by 0.04 at Mach 2.5.

## Sensitivity: five settings at $\delta_i/12$

The Mach 2.5 case with one setting changed in each run. Same averaging, 33–35 checkpoints of $180<t\le360$.

| Run | Change from the baseline | Cases file |
|---|---|---|
| lengths | the filter's integral lengths from $0.8, 0.3, 0.3$ to $1.2, 0.5, 0.5\,\delta_i$ | `tbl_cases_t37_hdef.txt` |
| span | span doubled to $12\,\delta_i$, $1080\times72\times144$ cells | `tbl_cases_t37_hdef.txt` |
| WALE, Vreman | `[turbulence] sgs_model = wale`, `vreman`, default constants | `tbl_cases_t37_sgs.txt` |

![Runs of the Mach 2.5 layer](/flume/tbl-m2p5-hdef.png)

*Four of the runs (the span of $12\,\delta_i$ in the table only). Left: $H$ against the edge state, against $x$; the square is the DNS at its station. Middle: $u^+_{VD}$ at the station of each run. Right:
$-\langle u'v'\rangle/u_\tau^2$ against $y/\delta$; black, the DNS.*

| At $Re_\theta=2835$ | baseline | lengths $1.2, 0.5, 0.5\,\delta_i$ | span $12\,\delta_i$ | WALE | Vreman | DNS |
|---|---|---|---|---|---|---|
| station $x/\delta_i$ | 55.2 | 51.4 | 54.6 | 54.7 | 57.5 | 53 |
| $H$, against the edge state | 4.08 | 4.08 | 4.11 | 4.05 | 4.11 | 4.10 (its profile: 4.22) |
| $C_f$, the wall model's $\tau_w$ | 2.208e-3 (−4.4 %) | 2.248e-3 (−2.7 %) | 2.205e-3 (−4.6 %) | 2.235e-3 (−3.3 %) | 2.206e-3 (−4.5 %) | 2.310e-3 |
| $u^+_{VD}$ at $y^+=30$, $100$, $200$ | 13.1, 17.5, 19.8 | 13.0, 17.7, 19.8 | 13.0, 17.5, 19.9 | 13.2, 17.4, 19.7 | 13.0, 17.5, 19.7 | 13.6, 16.7, 18.8 |
| $u_{rms}$, $v_{rms}$, $w_{rms}$ over $u_\tau$ at $y/\delta=0.2$ | 2.06, 0.81, 1.22 | 1.97, 0.75, 1.15 | 2.08, 0.79, 1.22 | 1.97, 0.76, 1.19 | 1.94, 0.74, 1.27 | 1.49, 0.93, 1.12 |
| $-\langle u'v'\rangle/u_\tau^2$ at $y/\delta=0.2$, $0.5$ | 0.80, 0.50 | 0.66, 0.57 | 0.74, 0.56 | 0.72, 0.55 | 0.66, 0.49 | 0.62, 0.37 |

**What it shows.**
- **None of the settings changes the layer much.** $H$ stays within 0.06 of the DNS's header value, $C_f$ within 2 % of
  the baseline's; longer integral lengths shorten the inflow's recovery a little (the layer reaches $Re_\theta=2835$ at
  $x=51$, $C_f$ −2.7 %), as $\delta_i/16$ does.
- **The subgrid model is active, and does not matter at this grid.** WALE's eddy viscosity, computed from a checkpoint,
  is about 1.5 times the molecular one over $0.3$–$0.75\,\delta$ and 0.4 at the edge.
- **What is left is in the stresses**, not in the mean: the outer $-\langle u'v'\rangle$ is 0.49–0.57 at $0.5\,\delta$
  against 0.37 in the DNS, and the near-wall $u_{rms}$ peak 1.9–2.1 against 1.5. At $\delta_i/16$ both move toward the
  DNS ([resolution](#resolution-mach-2-5-at-delta-i-16)).

## Mach 5.86, cold wall

Case M6Tw076 of the same DNS: Mach 5.86, a cold wall ($T_w/T_r=0.76$, $T_w=5.45\,T_\infty$), $Re_\theta=9455$ and
$Re_\tau=453$ at its station ($x=54\,\delta_i$), wall heat-flux parameter $B_q=q_w/(\rho_w c_p u_\tau T_w)=-0.02$ (heat
into the wall). The wall history gives $q_w$ as well as $\tau_w$.

| Setup | |
|---|---|
| box, cells | $120\times6\times6\,\delta_i$, $1440\times72\times72$, cubic $\delta_i/12$: the filter's recovery length grows with the Mach number (Huang, Duan & Choudhari 2022) |
| gas | $\mu\propto T^{0.76}$ with $\mu_w$ matched to the DNS's Sutherland value at the wall, where the model acts ($\mu_w/\mu_\infty=5.13$ against $5.45^{0.76}=3.63$): $Re=97120$, the DNS's $137350$ per $\delta_i$ over $1.414$; $Re_\theta$ is reported with $137350$ (`tbl_oracle.py compare --re-inf`) |
| gas, Sutherland run | `viscosity_law = sutherland` ([#121](https://github.com/szaghi/adam/issues/121) S1), $S=110.4$ K, the DNS's $Re=137350$ per $\delta_i$ (`tbl_cases_t37_hdef.txt`, 29 checkpoints) |
| wall, inflow, scheme | as at Mach 2.5 |
| run | `tbl_cases_t37_m6.txt`, to $t=360$ (21604 steps); 29 checkpoints of $180<t\le360$ and the wall history every 20 steps; 0.54 s per step on 4 A30 |

![The Mach 5.86 cold-wall boundary layer against the DNS](/flume/tbl-m6.png)

*The panels of the Mach 2.5 figure; the title adds $B_q$ at the station.*

| At $Re_\theta=9453$ ($x=55\,\delta_i$) | LES | LES, Sutherland | DNS |
|---|---|---|---|
| $C_f$, the wall model's $\tau_w$ | **9.47e-4** (the two halves of the window 9.70e-4, 9.25e-4); $x=60$–$110$: 9.78e-4 | 9.51e-4; $x=60$–$110$: 9.80e-4 | 9.84e-4 |
| $B_q$, the wall model's $q_w$ | **−0.0189** (halves −0.0191, −0.0187); $x=60$–$110$: −0.0194 | −0.0190 | −0.02 (one significant figure) |
| $\delta/\delta_i$ | 1.61 | 1.61 | 1.72 |
| $H$, against the edge state | 13.41 | 13.71 | 13.6 (its profile: 14.9) |
| $u^+_{VD}$ at $y^+=30$, $100$, $200$ | 12.4, 17.9, 20.8 | 12.4, 17.9, 20.6 | 13.9, 17.2, 19.5 |
| $u_{rms}$, $v_{rms}$, $w_{rms}$ over $u_\tau$ at $y/\delta=0.2$ | 2.21, 0.56, 0.94 | 2.24, 0.55, 0.90 | 1.33, 0.82, 0.99 |
| $-\langle u'v'\rangle/u_\tau^2$ at $y/\delta=0.2$, $0.5$ | 0.64, 0.39 | 0.65, 0.35 | 0.49, 0.25 |

**What it shows.**
- **$C_f$ within 3.8 % of the DNS** at its $Re_\theta$, with no change of the wall-model setting from Mach 2.5; past
  $x=60$ the plateau is 0.7 % below.
- **$B_q$ agrees to the precision of the reference.** The DNS file gives $-B_q=0.02$, one significant figure, and the run
  gives 0.0189 at the station, 0.0194 on the plateau.
- **$C_f$ drifts by 4.6 % between the two halves of the window; $q_w/\tau_w$ by 0.4 %.** The layer itself changes slowly
  over the window, and the wall model's heat flux follows its friction. 180 time units are 1.5 flow-throughs of this box,
  so the $C_f$ figure carries about $\pm2$ %.
- **The inflow recovers by $x\approx55$–$60\,\delta_i$**, farther than at Mach 2.5 (50): $C_f$ falls from
  $1.4\cdot10^{-3}$ at the inlet to $0.8\cdot10^{-3}$ at $x\approx33$, then rises to its plateau. The layer reaches the
  DNS's $Re_\theta$ at $x=55$, against 54 in the DNS.
- **$H$ is the DNS's** against the edge state (13.4; the header 13.6, its profile 14.9): the free stream over this layer
  is compressed by 4 % ([the shape factor](#the-shape-factor-and-the-free-stream)).
- **The stresses have the Mach 2.5 faults, larger.** $\delta$ is 6 % thin; $u_{rms}$ peaks at 1.7 times the DNS's at
  $y/\delta=0.2$, and $u^+_{VD}$ is 1.5 below the DNS at $y^+=30$ (the first cells under the matching point) and 0.7–1.3
  above it at $y^+=100$–$200$. At $\delta_i/12$ the layer of $Re_\tau=453$ has about 20 cells across $\delta$; at Mach 2.5
  a finer grid moves the log layer and the stresses toward the DNS ([resolution](#resolution-mach-2-5-at-delta-i-16)).
- **Sutherland's law with the DNS's Reynolds number changes little.** The power-law viscosity, matched at the wall, is
  1.41 times the DNS's in the free stream; Sutherland's law removes that: $H$ 13.71 against 13.41, $C_f$ −3.4 % against
  −3.8 % at the station, the plateau 9.80e-4, $B_q$ −0.0190.
- **The momentum integral does not measure $C_f$ here either** ($-62$ % at the station, negative in places).

## The wall model fed with the resolved convection

The zero-gradient check of the convection-fed model before it is used in the shock interaction
([turbulence](../../turbulence#fed-with-the-resolved-convection)): fed with the instantaneous resolved convection it
raised the $C_f$ of the interaction's incoming layer by 5 to 12 % ([interaction](../streams/sbli#the-run-s4h)).

| Setup | |
|---|---|
| runs | the Mach 2.5 case above (`m2p5-wm-l90`, the baseline, statistics $180<t\le360$), restarted at $t=360$ to 540 with `wall_model_pressure_gradient` (filter 10), `wall_model_convection` and `wall_model_resolved_stress`; statistics over the 27 checkpoints of $400<t\le540$ (`tbl_oracle.py compare --t0 400`); one node, 4 A30, about 3 h each |
| rs | the resolved shear stress removed from the forcing, the velocity scale of the eddy viscosity $\sqrt{(|\tau_w|+|S|)/\rho}$ (an intermediate build) |
| rs2 | the key's model: the resolved stress removed and the velocity scale $\sqrt{|\tau_w|/\rho}$ |
| a priori | `wmf_apriori.py` on the filtered averages of rs's `.wmf` (77760 wall columns): the model's variants on the same samples |
| f50, f100 | the slow forcing: rs restarted at $t=540$ to 900 with `wall_model_velocity_scale = sum` (the scale $\sqrt{(\lvert\tau_w\rvert+\lvert S\rvert)/\rho}$, the resolved stress removed) and the forcing filter 50 or 100 (`wall_model_pressure_gradient_filter`, against rs's 10); statistics over the 34 checkpoints of $720<t\le900$; one prod-gn job, a run per node, 6 h |

![The convection-fed wall model on the Mach 2.5 zero-gradient layer](/flume/tbl-convection.png)

*Left: $C_f$ from the wall model's $\tau_w$, the DNS at the baseline's station; right: the ratio to the baseline,
the stretch $45\le x\le80$ shaded; rs, rs2 and the slow forcing f50, f100. `make_doc_figures.py --only tbl-convection`.*

| | baseline | rs | rs2 | f50 | f100 |
|---|---|---|---|---|---|
| $C_f$ at the station ($Re_\theta=2835$) against the DNS | $-4.3$ % | $+2.5$ % | **$-4.4$ %** | $-1.7$ % | $-2.4$ % |
| $C_f$ over $45\le x\le80$ against the baseline | | $+6.6$ % | **$-0.3$ %** | $+2.2$ % | **$+1.9$ %** |
| the same, range along $x$ | | $+3.8$ to $+11.3$ % | $-4.3$ to $+3.3$ % | $-1.0$ to $+3.8$ % | $-1.9$ to $+5.1$ % |
| $C_f$ over $20\le x<40$ (the recovery from the inflow) against the baseline | | $+0.8$ % | $-8.6$ % | $-5.6$ % | $-6.2$ % |

| a priori, mean $\tau_w$ against the equilibrium model | |
|---|---|
| the pressure gradient alone | $+3.9$ % |
| fed with the convection | $+22.9$ % |
| and the resolved stress removed | $+12.4$ % |
| and the velocity scale from $\tau_w$ | $-1.5$ % |
| the resolved stress removed, fed the span-mean forcing | $+6.7$ % |

- **Removing the resolved stress is not enough.** The filtered forcing at the sample fluctuates as much as $\tau_w$ (rms
  $1.0\cdot10^{-3}$ against a mean $\tau_w$ of $0.97\cdot10^{-3}$), and in the velocity scale $|S|$ averages far above the
  mean forcing: the eddy viscosity grows and $\tau_w$ with it (rs, +6.6 % in the run, +12.4 % a priori).
- **With the forcing out of the eddy viscosity the layer is the baseline's.** rs2's $C_f$ is the baseline's within
  0.3 % over the developed stretch and $-4.4$ % against the DNS at the station (the baseline $-4.3$ %), as the a priori
  $-1.5$ % predicted: the gate passes.
- **Where the layer recovers from the inflow it is lower.** Over $20\le x<40$, behind the digital-filter inflow, rs2's
  $C_f$ is 8.6 % below the baseline's: there the resolved convection does not vanish, and the model answers to it. The
  DNS station lies beyond, so which of the two is right there is not measured.
- **The slow forcing keeps the layer within about 2 %.** With the scale $|\tau_w|+|S|$ (which keeps the eddy viscosity under an
  adverse gradient, [interaction](../streams/sbli#a-priori-the-velocity-scale-of-the-eddy-viscosity)) and a forcing filtered over 50
  or 100 instead of 10, the developed stretch sits $+2.2$ and $+1.9$ % above the baseline, against rs's $+6.6$ %: the filter does
  most of the work between 10 and 50, little beyond (the a priori limit of an infinite filter, $+1.7$ and $-0.7$ %, is near). With
  the criterion of 2 % fixed before the runs, f100 passes and f50 falls 0.2 % short, a difference within the scatter along $x$ (the
  first 20 checkpoints of the same window gave $+0.4$ and $+0.8$ %). Both move the station toward the DNS ($-1.7$, $-2.4$ % against
  the baseline's $-4.3$ %).

## The viscosity coupling

The zero-gradient gate of `wall_model_coupling = viscosity` ([turbulence](../../turbulence#viscosity-coupling),
[#126](https://github.com/szaghi/adam/issues/126)): the no-slip condition kept and the face stress the larger of the first cell's
molecular no-slip stress and the model's, along the first cell's wall-parallel velocity. In the shock interaction it separates the
wall stress in the mean ([S4l](../streams/sbli#the-viscosity-coupling-s4l)). Criteria fixed before the runs: $C_f$ over
$45\le x\le80$ within 2 % of the flux-coupled control, and the station against the DNS no worse than the baseline's $-4.3$ %.

| Setup | |
|---|---|
| runs | rs's state at $t=540$ (the section above) restarted to 900 with the plain equilibrium model and its eddy viscosity below $h$ (the baseline's setting, $h=5/24\,\delta_i$, 2.5 cells); vc with `wall_model_coupling = viscosity`, ctrl with the default `flux`; statistics over the 33 and 34 checkpoints of $720<t\le900$ (`tbl_oracle.py compare --t0 720`); one prod-gn job, a run per node (4 A30), 4 h |
| grid | the baseline's: $\delta_i/12$, $1080\times72\times72$ over $90\times6\times6\,\delta_i$, the first cell centre at $y^+\approx10$ |

![The viscosity coupling on the Mach 2.5 zero-gradient layer](/flume/tbl-coupling.png)

*Top: $C_f$ (the stress applied) against $x$, the dots at each run's station; vc against ctrl, the developed stretch shaded and the
gate's band. Bottom: the van Driest velocity at the station against the DNS; $C_f$ over $45\le x\le80$ against time.
`make_doc_figures.py --only tbl-coupling`.*

| | baseline | ctrl, flux | vc, viscosity |
|---|---|---|---|
| $C_f$ at the station ($Re_\theta=2835$) against the DNS | $-4.3$ % | $-4.3$ % | **$-11.3$ %** |
| $C_f$ over $45\le x\le80$ against ctrl (range along $x$) | | | **$-8.7$ %** ($-11.8$ to $-7.7$ %) |
| $C_f$ over $20\le x<40$ against ctrl | | | $-12.9$ % |
| $H$ at the station (DNS 4.10) | | 4.10 | 4.22 |
| $u^+_{VD}$ of the first three cells, $y^+\approx10$, 30, 50 (DNS 8.8, 13.6, 15.0) | | 11.3, 13.1, 14.7 | 8.9, 11.9, 14.1 |
| molecular no-slip stress of the first cell against the applied, $u_1^+/y_1^+$ | | 1.08 | 0.88 |

- **The gate fails: the coupling lowers the zero-gradient stress by 9 %.** Over the developed stretch vc's $C_f$ is 8.7 % below
  ctrl's, and $-11.3$ % against the DNS at the station; ctrl reproduces the baseline ($-4.3$ %). The drop settles within about 60
  time units of the restart and holds over the window (bottom right), so it is the coupling's, not a transient.
- **How: the maximum drags the first cell down to its no-slip velocity, and the sample with it.** Under the flux coupling the first
  cell runs at $u^+=11.3$ at $y^+\approx10$, faster than the DNS's 8.9 there: its no-slip stress would be 1.08 times the model's. Switched to
  the viscosity coupling, the maximum takes that larger stress (vc's $C_f$ first rises above ctrl's, bottom right), the first cell slows
  to $u^+=8.9$, the DNS's value within 2 %, and the cells above it slow too ($u^+=11.9$ at $y^+\approx30$ against the DNS's 13.6). The model reads
  the slower velocity at $h$ and returns a lower $\tau_w$, which the face then applies wherever it exceeds the molecular one.
- **De Vanna et al. (2022) reproduce the model's stress exactly with the resolved gradient** (an effective wall viscosity
  $\mu_{eff}=\tau_w^{wm}/(\partial u/\partial y)$, no lower bound stated); the lower bound $\mu_w^{sgs}\ge0$ here, which keeps the
  molecular stress where the wall is resolved, is what lets the first cell's no-slip gradient win. The variant without it, the
  model's magnitude along the first cell's velocity, is not tested.

## Recycling inflow

The recovery gate of the `inflow-recycling` face ([turbulence](../../turbulence#recycling-rescaling-inflow),
[#121](https://github.com/szaghi/adam/issues/121) S5): the digital filter's layer needs about $50\,\delta_i$ to reach its
developed $C_f$ (the dip in [the flow](#the-flow)); recycling should shorten it. The measure, fixed before the run: $C_f$ over
van Driest II at the run's own $Re_\theta(x)$, and the first $x$ beyond which that ratio stays within 3 % of its mean over
$45\le x\le75$ (the outflow end excluded). The reference is #126's ctrl ([the viscosity coupling](#the-viscosity-coupling)):
the same restart, model, window and grid, the digital filter at the face.

| Setup | |
|---|---|
| run | rs's state at $t=540$ restarted to 900 as ctrl (the plain equilibrium model, its eddy viscosity below $h$), the face `inflow-recycling` at station $x_r=20\,\delta_i$ (240 cells; `make_tbl.py --inflow recycling --rr-x 20`: thickness ratio 1.3305, $\beta=1.0343$, span shift one half); statistics over the 34 checkpoints of $720<t\le900$ (`tbl_oracle.py compare --t0 720`); one prod-gn node (4 A30), 4 h 36 min |
| grid | the baseline's: $\delta_i/12$, $1080\times72\times72$ over $90\times6\times6\,\delta_i$ |

![The recycling inflow against the digital filter on the Mach 2.5 layer](/flume/tbl-recycling.png)

*Left: $C_f$ (the stress applied) against $x$, the dots at each run's station. Middle: $C_f$ over van Driest II at the run's
$Re_\theta$, the band 3 % around each run's mean over $45\le x\le75$ and the recovery length (dotted); the recycling station
dashed. Right: the van Driest velocity at the station against the DNS. `make_doc_figures.py --only tbl-recycling`.*

| | ctrl, digital filter | rr, recycling |
|---|---|---|
| recovery length, 3 % of the run's own mean over $45\le x\le75$ (the gate's measure) | $50\,\delta_i$ (baseline 52) | **$26\,\delta_i$** |
| deepest departure beyond $x=20$ from the common level (1.004, the digital filter runs over $60\le x\le80$) | $-23.9$ % at $x=27$ | $+3.8$ % at $x=50$, $-2.0$ % at $x=22$ |
| first $x$ within 4 %, 3 % of the common level | 53, 55 | **10**, 55 |
| $C_f$ over $60\le x\le80$ against van Driest II | 1.004 | 1.005 |
| station ($Re_\theta=2835$): $x$, $C_f$ against the DNS, the ratio against the common level | 55.9, $-4.3$ %, $-2.8$ % | 43.3, $+1.4$ %, $+3.0$ % |
| $H$ at the station (DNS 4.10); $Re_\theta$ at the face | 4.10; 1469 | 3.98; 1518 |
| $u^+_{VD}$ at $y^+=30$, 100, 200 (DNS 13.6, 16.7, 18.8) | 13.0, 17.5, 19.9 | 13.1, 17.6, 19.8 |
| $u_{rms}$, $v_{rms}$, $w_{rms}$ over $u_\tau$ at $y/\delta=0.2$ (DNS 1.49, 0.93, 1.12) | 1.91, 0.81, 1.23 | 1.90, 0.74, 1.09 |
| wall time per step (whole run, one node; ctrl ran beside vc on a second node) | 0.61 s | 0.63 s (+3 %) |

- **Recycling removes the digital filter's dip.** Beyond the station the ratio never falls more than 2 % below the level both
  digital-filter runs reach downstream, against 24 % at $x\approx27$ under the digital filter; it is within 4 % of that level
  from $x=10$, the digital filter from $x=53$. The gate's measure gives $26\,\delta_i$ against 50, half.
- **The 3 % criterion against the common level does not separate them: both enter at $x=55$.** The recycled layer runs up to
  3.8 % *above* the downstream level (more than 3 % over $43<x<55$, more than 2 % over $42<x<59$) before it joins it; the
  gate's measure, which takes each run's own mean (1.020 for rr), absorbs that bump into its band. Which of the two references
  is right is a choice of criterion; the departure (4 % against 24 %) is not.
- **The station comparison is not a better wall model.** rr's $+1.4$ % and ctrl's $-4.3$ % against the DNS are each run's
  departure from the same downstream level ($+3.0$ and $-2.8$ %), taken where $Re_\theta=2835$ falls: inside rr's bump
  ($x=43$), at the end of ctrl's recovery ($x=56$). The level itself, carried to $Re_\theta=2835$ by van Driest II, is $-1.6$ %
  against the DNS for both. The profiles at the station agree within 0.1 in $u^+_{VD}$; the recycled layer's $v_{rms}$ and
  $w_{rms}$ at $0.2\,\delta$ are 9 and 11 % lower than ctrl's.
- **The cause of the bump is not identified.** The face's layer is 4 % thicker than the table's ($\delta$ 1.04 against 1.00,
  $Re_\theta$ 1518 against 1469): the correlation's ratios, fixed for the run, against the layer's own growth to the station are
  the first suspect, not tested.
- **The cost is 3 %**, the whole run's wall time per step against ctrl's on another node of the same kind; the first minutes ran
  at 0.82 s per step, a start-up transient.

## Not done

- The matching height at fixed grid; the outer-layer and near-wall stresses at $\delta_i/12$.
- $\delta_i/16$ at Mach 5.86; a longer averaging window at Mach 5.86.
- $Re_\theta$ and the station against the edge state: they use the nominal free stream.
- Mean profiles against the DNS beyond the station.
- The viscosity coupling without the lower bound on the wall's subgrid viscosity (the model's stress along the first cell's velocity).
- Recycling: the bump above the downstream level (ratios measured from the run against the correlation's), the station's
  distance, the recycling frequency $U/x_r$ in the incoming layer, a longer window.
