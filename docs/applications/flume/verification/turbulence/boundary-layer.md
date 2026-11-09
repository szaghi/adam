# Compressible boundary layer

| ID | Milestone | Type | Legs | Backends | Figure |
|---|---|---|---|---|---|
| TB | [#119](https://github.com/szaghi/adam/issues/119) · T3.7 | P+C | tb1–tb2 | CPU, FNL | yes |

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
  and at $\delta_i/16$ it is 0.3 % at the station ([resolution](#resolution-mach-2-5-at-delta-i-16)).
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

## Not done

- The matching height at fixed grid; the outer-layer and near-wall stresses at $\delta_i/12$.
- $\delta_i/16$ at Mach 5.86; a longer averaging window at Mach 5.86.
- $Re_\theta$ and the station against the edge state: they use the nominal free stream.
- Mean profiles against the DNS beyond the station.
