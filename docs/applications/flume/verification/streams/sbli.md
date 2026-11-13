# Shock/boundary-layer interaction

| ID | Milestone | Type | Legs | Backends | Figure |
|---|---|---|---|---|---|
| — (campaign `tbl_cases_s4.txt`) | [#121](https://github.com/szaghi/adam/issues/121) · S4b–S4h | P | none (the pressure-gradient key: wm12) | FNL | yes |

## Case and scope

The interaction of the IUSTI experiment (Dupont, Haddad & Debiève 2006, JFM 559) as computed by the STREAmS DNS
(Bernardini, Modesti, Salvadore & Pirozzoli 2021, CPC, Fig. 7): Mach 2.28, an oblique shock of $8^\circ$ flow deflection
on an adiabatic wall, $Re_\theta\approx5000$. The layer is the [T3.7 setup](../turbulence/boundary-layer); the shock comes
in through the top face, the [shock generator](./shock-generator) of S3. Nothing is fitted to the interaction: the wall
model is the one of the zero-pressure-gradient layers. The oracle is the interaction length, the separation from $C_f$
and the wall pressure against the DNS and the experiment.

## Setup

| Setup | |
|---|---|
| box, cells | $100\times12\times6\,\delta_i$, $1200\times144\times72$, cubic $\delta_i/12$ |
| shock | $\beta=32.685^\circ$, entering the top at $x=41.30$ so that the inviscid shock reaches the wall at $x_{imp}=60$, past the inflow's recovery; exact regular reflection: $p_2/p_1=1.602$, $p_3/p_1=2.458$, reflected at $29.77^\circ$ |
| gas | $Re=U_\infty\delta_i/\nu_\infty=35500$, Sutherland's law (S1) with $S/T_\infty=0.776$ (110.4 K over the experiment's 142 K), $Pr=0.71$ |
| wall, inflow, scheme | `wall-modelled` at $T_w=1.918\,T_\infty$ (the recovery temperature), the 3rd cell, semi-local damping, eddy viscosity below $h$; the M2p5 DNS profile with its temperature excess scaled by $(2.28/2.5)^2$ (`reference/tbl_m2p28_profile.txt`); `hybrid-kep-weno` |
| run | `tbl_cases_s4.txt`, to $t=1200$ on 4 A30, 0.74 s per step; statistics: the span and the 50 checkpoints of $310\le t\le934$ (`sbli_oracle.py compare --t0 300`) |
| references | the DNS wall pressure and the IUSTI points of STREAmS Fig. 7a, whose $x^*$ is measured from the reflected shock ([below](#where-the-compression-starts)); the STREAmS DNS at $Re_\theta=6882$ of Bernardini, Della Posta, Salvadore & Martelli (2022, arXiv:2204.10531) Fig. 4a; Volpiani, Bernardini & Larsson (2018, PRF 3) Fig. 12; digitised from the vector paths (`reference/sbli_m2p28_*`, error in the headers) |

`sbli_oracle.py` takes $\delta_0$, $\delta^*$ and $\theta$ of the undisturbed layer from a linear fit over
$[x_{imp}-15,\,x_{imp}-6]$ extrapolated to $x_{imp}$; $C_f$ from the wall model's $\tau_w$; the reflected shock from the
mid-pressure crossings of the rows $2.5<y/\delta_0<5$, a line extrapolated to the wall: its apparent origin $x_{ref}$
gives $L_{int}=x_{imp}-x_{ref}$, the length STREAmS scales $x$ with. It also reports the separation of the resolved flow, where the
time-mean $u$ of the first cell is negative, which $\tau_w$ misses when the reversed layer lies below the sample
([below](#the-reversed-flow-the-dns-s-bubble-and-the-runs-below-the-sample)).

## Results

![The interaction: wall pressure and skin friction](/flume/sbli.png)

*Left: the mean wall pressure $(p_w-p_\infty)/(p_2-p_\infty)$ against $x^*=(x-x_{imp})/L_{int}$, the DNS of STREAmS (line), the DNS at
$Re_\theta=6882$ (dashed) and the experiment (triangles).
Right: $C_f$; dotted the apparent origin of the reflected shock ($x^*=-1$), dashed the inviscid impingement. `make_doc_figures.py --only sbli`.*

| | FLUME, wall model | Reference |
|---|---|---|
| $Re_\theta$ at $x_{imp}$, nominal / against the edge state | 5247 / 4667 | about 5000 (DNS, experiment) |
| $\delta_0$, $\delta^*_0$, $H$ at $x_{imp}-10$ | $1.807$, $0.531\,\delta_i$, 3.58 | |
| reflected shock: angle to the wall | $28.95^\circ$ | $29.77^\circ$ inviscid |
| $L_{int}/\delta_0$ | **3.15** | 3.22, Volpiani et al. (DNS, $Re_\theta=1050$); 3.30, Bernardini et al. (2022, DNS, $Re_\theta=6882$) |
| separation: $L_{sep}/\delta_0$ | **0**, $C_f$ never negative, minimum $0.22\cdot10^{-3}$ at $x^*\approx-0.25$ | 1.62, Volpiani et al., $Re_\theta=1050$; 2.16, Bernardini et al. (2022), $Re_\theta=6882$; the experiment has a bubble |
| $(p_w-p_\infty)/(p_2-p_\infty)$ at $x^*=-1$, $-0.5$, $0$, $0.5$, $1$, $2$ | 0.02, 1.12, 1.72, 2.09, 2.24, 2.37 | DNS: 0.64, 1.75, 2.07, 2.22, 2.29, 2.36; at $Re_\theta=6882$: 0.87, 1.65, 2.01, 2.18, 2.27, 2.35 |
| Souverein $L^*=L_{int}/(\delta^*g_3)$, separation parameter $Se$ | 1.93, 1.20 | |

## The flow

![The interaction's flow field](/flume/sbli-fields.png)

*One checkpoint ($t=934$). From the top: numerical schlieren of the mid-span plane; the spanwise-averaged pressure with the sonic line (white), the
$\bar u=0$ line (red, absent) and the oracle's shock fits (orange, its crossings as dots); $u$ in the wall-parallel plane $y=0.29$
($0.16\,\delta_0$; dashed the reflected shock's apparent origin, dotted $x_{imp}$); the cross-flow plane at $x_{imp}$ and spanwise-averaged
profiles. `make_doc_figures.py --only sbli-fields`.*

![The Mach 2.5 boundary layer and the interaction in 3-D](/flume/tbl-3d.png)

*A cutaway, $u/U_\infty$, $y$ drawn stretched twice: the Mach 2.5 [boundary layer](../turbulence/boundary-layer) over $15<x<45$ (top), and
this interaction over $35<x<60$, the cross-flow face at the impingement $x=60$ (bottom). `make_doc_figures.py --only tbl-3d`.*

## What it shows

- **The interaction length is right.** $L_{int}/\delta_0=3.15$ against 3.22 in the DNS of Volpiani et al. (2018, at
  $Re_\theta=1050$), and the reflected shock leaves at $28.95^\circ$, $0.8^\circ$ under the inviscid angle. The
  downstream pressure plateau is the DNS's within 4 % ($x^*=2$).
- **The mean flow does not separate.** $C_f$ dips to $0.22\cdot10^{-3}$ and recovers; the DNS, the experiment and
  Volpiani et al. have a bubble. It comes close: in one checkpoint the spanwise-averaged $u$ of the first cell falls to
  $0.003\,U_\infty$ between $x=55$ and 65, and the wall-parallel plane at $0.16\,\delta_0$ has patches of $u<0$ (down to
  $-0.14\,U_\infty$) between $x=57.7$ and 61.9. The wall model's $\tau_w$ follows the velocity at its matching height,
  the 3rd cell, whose mean stays positive.
- **The wall pressure rises late.** It starts at $x^*\approx-0.75$ and is half-way to $p_2$ at $-0.62$; the DNS are there
  at $-1.25$ and $-1.09$, in the separation that the DNS have and this run has not. Downstream of the impingement the
  pressure is the DNS's within 4 %. (The first version of this page read the reference's $x^*$ from $x_{imp}$, as the
  STREAmS paper's text says, and found the rise early: the figure measures it from the reflected shock,
  [below](#where-the-compression-starts).)
- **Two candidate causes of the missing separation, tested since.** The equilibrium model has no pressure-gradient term,
  and it computes $\tau_w$ from the velocity at the 3rd cell, $0.21\,\delta_i$ or $0.12\,\delta_0$ above the wall, which
  the near-wall deceleration under the adverse gradient reaches late. [S4c](#a-refined-wall-grid-under-the-interaction)
  halves the matching height: the dip of $C_f$ halves, the mean does not separate.
  [S4d](#the-wall-model-with-the-pressure-gradient) adds the pressure gradient: a bubble of $0.6\,\delta_0$, the
  compression unchanged.

## A refined wall grid under the interaction

Run S4c is the S4b case with one more level, 2:1, under the interaction (`tbl_cases_s4c.txt`). The wall model and its matching
point (the 3rd cell) are unchanged, so in the refined box the matching height halves with the cells. A bubble in this
run points at the near-wall grid; none points at the equilibrium model.

| Setup | |
|---|---|
| grid | S4b's box and base grid ($\delta_i/12$) with a static 2:1 box: `make_tbl.py --refine-box 30,75,2.5` refines a block when its centre is inside $30\le x\le75$, $y\le2.5$, which takes the blocks of $25\le x\le75$, $y<3$, every $z$, to $\delta_i/24$; 120 blocks (64 fine, 56 coarse), about 23 M cells |
| matching height | the 3rd cell: $0.104\,\delta_i$ in the box against S4b's $0.208\,\delta_i$ ($0.06$ against $0.12\,\delta_0$) |
| seams | the box's upstream face is at $x=25$, 35 $\delta_i$ ahead of the impingement and 29 ahead of S4b's reflected-shock origin; its top at $y=3$, $1.7\,\delta_0$, inside the layer (FV reflux at the 2:1 faces) |
| start | from $t=0$ with the box in place, not from S4b's state: moving a developed state onto a refined grid needs runtime regridding along a modelled wall, which is not verified (a static box is what the wall-model leg `wm3` verifies). The cost is S4b's $\approx300$ time units of transient again |
| run | to $t=550$ on 2 nodes (8 A30), a checkpoint every 2000 steps ($\approx14$ time units); statistics over $t>300$, where S4b's had converged (S4b changed by at most 1 % between the windows $310$–$934$ and $310$–$1200$) |
| oracle | `sbli_oracle.py compare --t0 300`: `tbl_oracle.averages` maps the two levels onto the fine grid (coarse cells repeated); $C_f$ comes from the wall model's history, 1200 bins of equal width weighted by area, so the levels need no special case. `field_slices.py` cuts the planes on the fine grid the same way |

### Cost on two nodes

The first FNL run of FLUME across two medooza nodes. S4b ran 64 blocks on 4 A30 (one node) at 0.74 s per step; S4c
runs 120 blocks of the same size on 8 A30 at 0.95 s per step at the start and 0.98 at $t\approx400$ (the time step
falls as the layer develops; the steps cost the same). Per step the busiest GPU advances 4.1 M cells per second against
S4b's 4.2 M: the inter-node halo exchange and the 2:1 seams (reflux included) cost about 2–3 %.

The step was slower because of the partition, not the network. The 120 blocks landed 20/13/19/17/13/19/19/0 on the 8
ranks: one GPU held no block, and the busiest held 20 against an even 15, which made each step 1.33 times slower than
an even split. `mpi_redistribute` then walked the blocks greedily and kept the 8 siblings of a refined parent together;
a model of the routine reproduces that split exactly ([#123](https://github.com/szaghi/adam/issues/123)). A static grid
is now balanced by block count (`[amr] balance_static`, on by default: 15 per rank on this tree), and a runtime regrid
by a min-max split of the sibling sets (busiest 18). Evenly split, the step would take about 0.71 s, so the run would
need about 13.5 hours instead of 18.

### Results

![S4c against S4b: wall pressure and skin friction](/flume/sbli-refined.png)

*Left: the mean wall pressure of S4b and S4c against the DNS and the experiment. Right: $C_f$ from the wall model's
$\tau_w$ against $x$, the 2:1 box of S4c shaded. `make_doc_figures.py --only sbli-refined`.*

The statistics are the span and the 21 checkpoints of $307.5\le t\le550$ (`sbli_oracle.py compare --t0 300`); S4b's
column repeats the table above.

| | S4c, 2:1 box | S4b |
|---|---|---|
| $Re_\theta$ at $x_{imp}$, nominal; $\delta_0$ | 5314; $1.824$ | 5247; $1.807$ |
| reflected shock: angle to the wall | $28.63^\circ$ | $28.95^\circ$ |
| $L_{int}/\delta_0$ | **3.04** | 3.15 (DNS of Volpiani et al.: 3.22) |
| separation: $L_{sep}/\delta_0$ | **0**, $C_f$ never negative, minimum $0.10\cdot10^{-3}$ at $x=58.4$ | 0, minimum $0.22\cdot10^{-3}$ at $x=58.6$ |
| $(p_w-p_\infty)/(p_2-p_\infty)$ at $x^*=-0.5$, $0$, $0.5$, $1$, $2$ | 0.89, 1.76, 2.14, 2.26, 2.39 | 1.12, 1.72, 2.09, 2.24, 2.37 (DNS: 1.75, 2.07, 2.22, 2.29, 2.36) |
| $C_f$ upstream of the box ($x\le24$) | S4b's within 0.4 % | |
| $C_f$ in the box, before the interaction | $-11$ % at $x=26$–$30$, behind the 2:1 face; $+6$ % at $x=35$; $+10$ % at $x=45$–$55$ | |
| $C_f$ downstream of the box ($x=76$–$90$) | $-7$ % | |

![The flow of S4c](/flume/sbli-refined-fields.png)

*The panels of [the flow](#the-flow) for S4c at $t=550$, on the fine grid of the box (the coarse cells repeated
outside it, the box outlined in orange in the schlieren). `make_doc_figures.py --only sbli-refined-fields`.*

- **The finer near-wall grid does not open a mean bubble.** Halving the matching height halves the dip of $C_f$, at the
  same place, and the mean stays positive. In single checkpoints the flow at the wall reverses in patches: at $t=550$
  the spanwise-averaged $u$ of the first cells is negative over a short stretch near $x=58.5$ (the red $u=0$ line), and
  the cross-flow plane at $x_{imp}$ has reversed cells along the wall. The reversal is intermittent; the mean has none.
- **The interaction length shrinks slightly,** 3.04 against 3.15, and the rise is shallower ($0.89$ at $x^*=-0.5$
  against S4b's $1.12$ and the DNS's $1.75$), its onset unchanged ([below](#where-the-compression-starts)): the shape of
  the interaction is S4b's.
- **The zero-pressure-gradient layer depends on the matching height.** Outside the box $C_f$ is S4b's; inside it, after
  an adjustment over about $10\,\delta_i$ behind the 2:1 face, it is 10 % higher before the interaction. The
  equilibrium model's $\tau_w$ follows the velocity at its matching point, and that point moved from $0.12$ to
  $0.06\,\delta_0$ (the log-layer mismatch of Kawai & Larsson 2012). A matching height fixed in $\delta_i$, not in
  cells, would keep the two levels consistent.
- **What remains is the model.** The equilibrium model integrates a constant total stress across its layer, so its
  $\tau_w$ turns negative only where the velocity at the matching height reverses on average. With the matching point
  at $0.06\,\delta_0$ it does not; the pressure-gradient term, which lets $\tau_w$ reverse under a strong adverse
  gradient while the matching velocity is still positive, was tested next
  ([below](#the-wall-model-with-the-pressure-gradient)): a bubble of $0.6\,\delta_0$, a quarter of the DNS's.

## The wall model with the pressure gradient

The equilibrium model has its $\tau_w$ from a constant-stress layer, so it reverses only where the mean velocity at the
matching height does. With `wall_model_pressure_gradient` ([turbulence](../../turbulence#pressure-gradient)) the layer
carries $\tau_w+y\,dp/dx$, and an adverse gradient strong enough turns $\tau_w$ negative while the sample still moves
downstream.

### A priori, on the mean fields

Before any run, the model was applied to the span and time averages of S4b and S4c: per $x$, the mean sample at the
run's matching point and the gradient of the mean wall pressure, in a Python port of the model (`sbli_oracle.py means`
and `apriori`; the port reproduces the Fortran solver to $10^{-14}$).

| | S4b, $h=0.208\,\delta_i$ | S4c, $h=0.104\,\delta_i$ in the box |
|---|---|---|
| the model without $dp/dx$ on the mean sample, against the run's $C_f$ ($x=35$–$50$) | $-1.7$ % | $-2.0$ % |
| without $dp/dx$: $C_f<0$ | nowhere (minimum $0.10\cdot10^{-3}$) | nowhere (minimum $0.03\cdot10^{-3}$) |
| **with $dp/dx$: $C_f<0$** | $x=56.1$–$61.2$, $L_{sep}/\delta_0=2.81$ | $x=56.1$–$60.7$, $L_{sep}/\delta_0=2.54$ |
| with $dp/dx$: $C_f$ minimum | $-1.58\cdot10^{-3}$ at $x=56.5$ | $-1.28\cdot10^{-3}$ at $x=58.9$ |

![The equilibrium model with and without the pressure gradient on the mean samples](/flume/sbli-apriori.png)

*Top: $C_f$ of the runs (the wall model's history), of the model on the mean sample, and of the model with the mean
$dp/dx$, the reversed stretch shaded; bottom: the gradient of the mean wall pressure. `make_doc_figures.py --only
sbli-apriori`.*

- **The term opens a bubble on both grids,** 2.5–2.8 $\delta_0$ long, against 1.62 in the DNS of Volpiani et al.
  (at $Re_\theta=1050$). The model without it reproduces the runs' $C_f$ to 2 %, so the port is the run's model.
- **The term dominates there:** $(dp/dx)\,h/\tau_w$ reaches about 70 in the bubble, far from the equilibrium layer
  the model assumes. The model keeps $dp/dx$ and drops the convective terms, which partly balance it in the wall layer:
  it can over-respond.
- **The mean field is the attached one.** With a bubble the pressure rise spreads and $dp/dx$ falls, so the run
  decides.
- **The instantaneous gradient cannot be used as it is.** On a single checkpoint of S4b upstream of the interaction
  ($x=35$–$50$) the turbulent wall-pressure gradient has an rms of $1.1\cdot10^{-2}$ against a mean of
  $8\cdot10^{-4}$, so $|dp/dx|\,h/\tau_w$ has a median of 1.4, and the model with it raises the mean $\tau_w$ of the
  zero-pressure-gradient layer by 13 %. The run filters the gradient in time (`wall_model_pressure_gradient_filter`).

### The run

| Setup | |
|---|---|
| start | S4b's last state ($t=1200$), the developed layer and interaction; its block layout kept (`[amr] balance_static = .false.`) |
| model | `wall_model_pressure_gradient = .true.`, `wall_model_pressure_gradient_filter = 10` (about $5\,\delta_0/U_\infty$, far below the period of the interaction's low-frequency motion); the rest as S4b (the 3rd cell, $h=0.208\,\delta_i$, semi-local damping, the model's eddy viscosity below $h$) |
| run | to $t=1450$ on one node (4 A30), 0.74 s per step, S4b's cost; statistics over the 24 checkpoints of $t\ge1300$ (`sbli_oracle.py compare --t0 1300`); the input is S4b's `make_tbl.py` line with `--restart --time-max 1450 --wm-pressure-gradient --wm-pg-filter 10` (no case file: a one-off restart) |

![S4d against S4b and S4c: wall pressure and skin friction](/flume/sbli-pg.png)

*Left: the mean wall pressure of S4b, S4c and S4d against the DNS and the experiment; right: $C_f$ from the wall
model's $\tau_w$, the reversed stretch of S4d shaded. `make_doc_figures.py --only sbli-pg`.*

| | S4d, with $dp/dx$ | S4b | S4c, 2:1 box |
|---|---|---|---|
| **separation: $L_{sep}/\delta_0$** | **0.60**, $C_f<0$ from $x=56.6$ to $57.7$, minimum $-1.9\cdot10^{-3}$ | 0 | 0 |
| $L_{int}/\delta_0$ (DNS of Volpiani et al.: 3.22) | 3.12 | 3.15 | 3.04 |
| reflected shock: angle to the wall | $28.2^\circ$ | $28.95^\circ$ | $28.6^\circ$ |
| $(p_w-p_\infty)/(p_2-p_\infty)$ at $x^*=-0.5$ (DNS 1.75; at $Re_\theta=6882$ 1.65) | 0.95 | 1.12 | 0.89 |
| $Re_\theta$ at $x_{imp}$; $\delta_0$ | 5288; 1.794 | 5247; 1.807 | 5314; 1.824 |
| $C_f$ upstream ($x=35$–$50$) against S4b | $-1.7$ % (S4b's own windows differ by up to 1 %; the unfiltered gradient gave $+13$ % a priori) | | |

![The flow of S4d](/flume/sbli-pg-fields.png)

*The panels of [the flow](#the-flow) for S4d at $t=1450$. `make_doc_figures.py --only sbli-pg-fields`.*

- **The mean flow separates,** for the first time in these runs: $C_f$ is negative over $0.60\,\delta_0$ at the foot
  of the reflected shock's compression, where the filtered gradient is largest, and reattaches before the
  impingement. The a priori test on the mean fields gave 2.5–2.8 $\delta_0$ for the same pressure field; the run's
  is a quarter of that, the gradient being filtered and the pressure rise unchanged.
- **The bubble is a quarter of the DNS's.** The DNS at $Re_\theta=6882$ separates from $x^*\approx-0.8$ to $-0.15$,
  $L_{sep}/\delta_0\approx2.16$, ahead of the impingement as here; this one spans $x^*=-0.60$ to $-0.41$, inside it.
- **The interaction does not change.** The wall pressure still rises at $x^*\approx-0.7$, half-way at $-0.56$ against
  the DNS's $-1.09$, and $L_{int}$ is S4b's: the model's stress reverses, but the resolved flow above the matching
  height, which sets the pressure, does not.
- **What remains.** The model keeps $dp/dx$ and drops the convective terms (the known inconsistency): a
  non-equilibrium layer (the thin-boundary-layer equations with convection, or an integral model) is the next test of
  the wall model. Its target is the DNS's separation: from $x^*\approx-0.8$, $2.2\,\delta_0$ long, and the wall
  pressure half-way at $x^*\approx-1.1$ ([below](#where-the-compression-starts)).

## Where the compression starts

The earlier versions of this page found the wall pressure rising too early: from $x^*\approx-0.7$, where the DNS of
STREAmS Fig. 7a was flat until $-0.3$. The comparison was wrong, not the run. The STREAmS paper defines
$x^*=(x-x_{imp})/L$, but its figure follows the IUSTI convention of Dupont et al. (2006), where the origin is the
mean position of the reflected shock (Bernardini et al. 2022, the footnote to their $L_{int}$). Shifted by one,
$x^*-1$, the Fig. 7a curve falls on the STREAmS DNS at $Re_\theta=6882$ of Bernardini et al. (2022, Fig. 4a, $x^*$
from $x_{imp}$). Over $-1.8\le x^*\le2.5$ the rms difference is $0.034\,p_\infty$, and the best shift is 1.02. A
reading of the curve itself agrees: at $x^*=0$ the printed DNS is half-way up the rise ($0.64$), the value at the foot
of a shock, not on the plateau a separated interaction has at its impingement. The two DNS and the experiment then
agree with each other.

![Where the compression starts](/flume/sbli-onset.png)

*Left: the references as printed (dotted, open triangles) and shifted (solid, filled), the DNS at $Re_\theta=6882$
(dashed), and the runs. Right: the mean Mach number of the incoming layer at $x_{imp}-5\,\delta_0$ against the DNS
profile the inflow is made from. `sbli_oracle.py fields` and `compression`; `make_doc_figures.py --only sbli-onset`.*

| | S4b | S4c, 2:1 box | S4d, with $dp/dx$ | DNS, $Re_\theta=6882$ |
|---|---|---|---|---|
| $L_{int}/\delta_0$ | 3.16 | 3.05 | 3.13 | 3.30 |
| wall pressure at 5 % of $p_2-p_\infty$: $x^*$ ($\delta_0$ ahead of $x_{imp}$) | $-0.74$ (2.34) | $-0.76$ (2.33) | $-0.68$ (2.12) | $-1.25$ (4.1) |
| at 50 %: $x^*$ ($\delta_0$ ahead of $x_{imp}$) | $-0.62$ (1.96) | $-0.66$ (2.00) | $-0.56$ (1.77) | $-1.09$ (3.6) |
| $L_{sep}/\delta_0$ | 0 | 0 | 0.60 | 2.16 ($x^*=-0.8$ to $-0.15$) |
| incoming layer at $x_{imp}-5\,\delta_0$: mean sonic height | $0.028\,\delta_0$ | $0.031\,\delta_0$ (1.4 cells) | $0.041\,\delta_0$ | $0.036\,\delta$ (the inflow's DNS profile) |

$L_{int}$ here is from the time-mean field (`fields`, `compression`), 0.01 above `compare`'s in the tables above.

- **The incoming layer is not the cause.** Its mean Mach number is the DNS profile's to 0.03 above $0.3\,\delta$ and
  to 0.15 at the 2nd cell ($0.13\,\delta$, lower: the log-layer mismatch of the wall model), and its mean sonic height,
  $0.03$–$0.04\,\delta_0$, is the DNS profile's $0.036\,\delta$. The first cell, below that height, runs at
  $M\approx0.95$ against the DNS's 0.85: the wall model's slip, one cell thick.
- **The compression is late by about $0.5\,L_{int}$, $1.6$–$1.8\,\delta_0$, on every grid and model.** Below
  $0.5\,\delta_0$ it crosses the layer as a vertical front at $x^*\approx-0.7$ to $-0.6$; above, it is the incident
  shock. The reflected shock's apparent origin is right ($L_{int}$ within 5 %), the plateau behind the impingement too.
- **What is missing is the separation.** In the DNS a bubble of $2.2\,\delta_0$ starting at $x^*\approx-0.8$ lifts
  the flow and carries the compression to $x^*\approx-1.2$. The equilibrium model opens none (S4b, S4c) or a quarter of
  it (S4d), and the late compression follows. Halving the near-wall cells (S4c) leaves the onset where it was. The pressure-gradient term moved the onset by $0.06\,L_{int}$
  *downstream*. S4e qualifies this bullet and the first: without the model's eddy viscosity below $h$ the near-wall
  layer slows and the onset moves $1.2\,\delta_0$ upstream with no bubble
  ([below](#without-the-eddy-viscosity-below-the-matching-height)).
- **For the next wall model,** non-equilibrium: the target is the DNS's bubble ($x^*\approx-0.8$ to $-0.15$) and its
  wall pressure, half-way at $x^*\approx-1.1$. The run's wall pressure is printed by `sbli_oracle.py compare` (at fixed
  $x^*$, against both DNS) and `compression` (where it reaches 5, 10 and 50 %); the oracle reports the run's $L_{sep}$,
  the DNS's is quoted from Bernardini et al. (2022).

## Wall models for this interaction in the literature

Before a new wall model, what the literature says about wall-modelled LES of separating interactions, and the test it
suggests on these runs.

| Source | Case | Finding |
|---|---|---|
| Hickel, Touber, Bodart & Larsson (2012, CTR Summer Program) | this interaction (IUSTI, $8^\circ$), a priori on the wall-resolved LES of Touber & Sandham (2009); an adverse-pressure-gradient layer | Above the viscous layer the convective term balances the pressure gradient (Bernoulli), so a model that keeps $dp/dx$ and drops convection, as here, is inconsistent: on the adverse-gradient layer it separates much too early, and the model of Duprat et al. (2011), whose velocity scale S4d uses, over-predicts the wall shear of the interaction. Convection parametrised as $-\frac{1}{\rho}\frac{dp}{dx}\min(y/y_{pg},1)$, $y_{pg}=4\,(\mu_w^2/(\rho_w\lvert dp/dx\rvert))^{1/3}$, gives a skin friction close to the reference |
| Larsson, Kawai, Bodart & Bermejo-Moreno (2016, Mech. Eng. Rev. 3) | review | Equilibrium models can be accurate in these interactions "provided that the outer part of the boundary layer is well resolved" (a duct interaction near incipient separation; a ramp at $Re_\theta\approx50\,000$, where the equilibrium model fell within the experimental uncertainty as the non-equilibrium one did): no conclusive demonstration that a non-equilibrium model is better yet. Outer-layer grid: $\Delta x\approx0.08\,\delta$, $\Delta z\approx0.05\,\delta$, $\Delta y$ from $0.02\,\delta$ at $y=0.2\,\delta$ to $0.05\,\delta$ at $\delta$ |
| Fu, Karp, Bose, Moin & Urzay (2021, JFM) | Mach 6, laminar inflow, transitional interaction | The equilibrium-model LES delays the separation and advances the reattachment: a shorter bubble than the DNS's (the symptom here, in a different regime) |

**Against these runs.** S4b's grid ($0.046\,\delta_0$, cubic) is finer than the guideline in $x$, equal in $z$ and 2.3
times coarser in $y$ at $0.2\,\delta$; S4c's box ($0.023\,\delta_0$) meets it everywhere, and does not separate either.
The Bernoulli-consistent model a priori on the mean fields (`sbli_oracle.py apriori --ypg 4`, the port of the
pressure-gradient model with $y$ replaced by $g(y)=y-y^2/(2y_{pg})$ below $y_{pg}$ and $y_{pg}/2$ above):

| | model, no $dp/dx$ | with $dp/dx$ (S4d's model) | Bernoulli-consistent, $y^*_{pg}=4$ | DNS, $Re_\theta=6882$ |
|---|---|---|---|---|
| S4b: $L_{sep}/\delta_0$ | 0 | 2.81 | 1.01 ($x=57.8$–$59.6$) | 2.16 |
| S4c: $L_{sep}/\delta_0$ | 0 | 2.54 | 1.53 ($x=57.4$–$60.2$) | |

- **The a priori numbers bound the model, not the run.** S4d's run kept a quarter of its a priori bubble (0.60 against
  2.81), and its strongly negative $\tau_w$ left the interaction where it was: the resolved flow above the matching
  height decides. A Bernoulli-consistent model alone is unlikely to reach $2.2\,\delta_0$.
- **The suspect is the resolved near-wall flow.** The first cell runs at $M\approx0.95$–$0.98$ against the DNS profile's
  0.85 at that height ([above](#where-the-compression-starts)): extra momentum at the wall that the adverse gradient has
  to stop. The model's eddy viscosity below the matching height (`wall_model_eddy_viscosity`) shapes that layer, so
  S4e takes it out.

### Without the eddy viscosity below the matching height

| Setup | |
|---|---|
| start | S4b's last state ($t=1200$); its block layout kept (`[amr] balance_static = .false.`) |
| model | S4b's without `wall_model_eddy_viscosity` (the cells below $h$ get the SGS model's viscosity only); no $dp/dx$ |
| run | to $t=1450$ on one node (4 A30), 3 h 24 min; statistics over the 23 checkpoints of $t\ge1300$ (`sbli_oracle.py compare --t0 1300`, `fields`, `compression`); the input is S4b's line of `tbl_cases_s4.txt` without `--wm-eddy-viscosity`, with `--restart --time-max 1450` (no case file: a one-off restart) |

![S4e against S4b and S4d](/flume/sbli-nonu.png)

*Left: the mean wall pressure of S4b, S4d and S4e against the DNS and the experiment; middle: $C_f$ from the wall
model's $\tau_w$, the reversed stretch of S4d shaded; right: the mean Mach number of the incoming layer at
$x_{imp}-5\,\delta_0$ against the DNS profile the inflow is made from. `make_doc_figures.py --only sbli-nonu`.*

| | S4e, no eddy viscosity below $h$ | S4b | DNS, $Re_\theta=6882$ |
|---|---|---|---|
| $L_{int}/\delta_0$ | **3.51** | 3.16 | 3.30 |
| wall pressure at 5 % of $p_2-p_\infty$: distance ahead of $x_{imp}$ | $3.55\,\delta_0$ | $2.34\,\delta_0$ | $4.13\,\delta_0$ ($x^*=-1.25$) |
| the same at 50 % | $2.64\,\delta_0$ | $1.96\,\delta_0$ | $3.60\,\delta_0$ ($x^*=-1.09$) |
| separation: $L_{sep}/\delta_0$ | **0**, $C_f$ minimum $0.27\cdot10^{-3}$ at $x=58.5$ | 0, minimum $0.21\cdot10^{-3}$ at $x=58.6$ | 2.16 |
| incoming layer at $x_{imp}-5\,\delta_0$: first-cell Mach number | 0.75 (DNS profile 0.82 there) | 0.98 (DNS profile 0.83) | |
| the same: second and third cells | 1.08, 1.30 (DNS 1.23, 1.37) | 1.13, 1.26 (DNS 1.24, 1.37) | |
| the same: sonic height | $0.061\,\delta$ | $0.030\,\delta$ | $0.036\,\delta$ |
| incoming layer: $\delta_0$, $\delta^*$, $Re_\theta$ | 1.85, 0.572, 5463 | 1.81, 0.531, 5271 | |
| $C_f$ at $x=35$–$45$ | $+4.8$ % on S4b | | |

- **The eddy viscosity moves the onset, not the separation.** Without it the first cell slows from $M=0.98$ to 0.75,
  the sonic line rises from $0.030$ to $0.061\,\delta$, and the compression starts $1.2\,\delta_0$ earlier: about
  two thirds of the way to the DNS at 5 % of the rise, 40 % of the way at 50 %. $L_{int}$ overshoots the DNS by 6 %
  (S4b undershoots by 4 %). The mean flow still does not separate, and the minimum of $C_f$ stays near S4b's. The
  pressure-gradient term (S4d) did the opposite: a bubble, the onset moved by $0.06\,L_{int}$ downstream.
- **The incoming layer is not the DNS's either way.** S4b's first cell is too fast and its next two too slow; S4e's first
  three are all too slow, and the sonic height goes from 15 % below the DNS's to 70 % above. $\delta^*$ grows by 8 % and
  $C_f$ upstream by 5 %: S4e changes the layer that meets the shock, so it is a diagnostic, not a candidate setting.
  The run stayed stable for its 250 time units (the eddy viscosity was introduced against the first-cell drain of the
  Mach 3 channel, [turbulence](../../turbulence)); longer runs were not made.
- **The mean is short.** 150 time units against S4b's 890; the interim mean over $t=1301$–$1363$ gave
  $L_{int}/\delta_0=3.52$, the full one 3.51.
- **What it points to.** The two features have different levers: the onset answers to the momentum of the resolved
  near-wall layer, the separation to the sign of $\tau_w$. A Bernoulli-consistent model changes $\tau_w$ only (a
  priori $1.0$–$1.5\,\delta_0$, above), so it is expected to open a bubble without moving the onset, as S4d did; the
  onset needs an incoming near-wall layer closer to the DNS profile, which is the matching between the model and the
  resolved cells, not the model's momentum balance.

## The near-wall cells of the incoming layer

Before changing how the model meets the resolved cells: how far the cells of the incoming layer are from the DNS profile
the inflow is made from, on both grids and with and without the eddy viscosity below $h$.

| Setup | |
|---|---|
| data | the span and time averages of `sbli_oracle.py fields` at $x_{imp}-5\,\delta_0$ (S4b $t=310$–1200, S4c 308–550, S4d 1305–1450, S4e 1301–1450); S4c on the fine grid of its box |
| reference | `reference/tbl_m2p28_profile.txt` at the cell centres, $y/\delta$ with the run's local $\delta$, $u/u_e$ |
| figure | `make_doc_figures.py --only sbli-cells` |

![The near-wall cells of the incoming layer](/flume/sbli-cells.png)

*Left: the mean velocity of the first 8 cells of S4b, S4c, S4d and S4e against the DNS profile; right: the difference of
each cell from the DNS at its centre. A cell holds the mean over its height: in the first cell, where the profile bends,
the DNS's cell mean is below its centre value, so the excess of the runs is a lower bound.*

| $u/u_e$, cells 1 – 4 | 1 | 2 | 3 | 4 |
|---|---|---|---|---|
| DNS at S4b's centres ($y/\delta$ 0.025, 0.075, 0.12, 0.17) | 0.478 | 0.669 | 0.722 | 0.757 |
| S4b | 0.558 (+17 %) | 0.625 ($-7$ %) | 0.684 ($-5$ %) | 0.758 (0 %) |
| S4d, with $dp/dx$ | +14 % | $-9$ % | $-7$ % | 0 % |
| S4e, no eddy viscosity below $h$ | 0.441 ($-7$ %) | 0.602 ($-10$ %) | 0.698 ($-3$ %) | $-1.5$ % |
| S4c, the box ($\Delta y/2$; $y/\delta$ 0.012, 0.037, 0.062, 0.087) | +70 % | +2 % | $-2$ % | +2 % |

The mean stress through the face between cells 1 and 2, in units of $\tau_w$ (the time-mean strain times the mean
viscosities; the resolved part $-\bar\rho\,(\overline{uv}-\bar u\bar v)$):

| | viscous | resolved | the model's eddy viscosity | sum |
|---|---|---|---|---|
| S4b | 0.03 | 0.20 | 0.59 | 0.82 |
| S4e | 0.08 | 0.33 | — | 0.41 |

- **The first cell slips with the eddy viscosity and lags without it.** With the model's eddy viscosity below $h$ it
  is 14–17 % faster than the DNS on S4b's grid and 70 % on S4c's finer one; without it 7 % slower.
- **Cells 2 and 3 are 5–10 % slow on S4b's grid, within 2 % on S4c's.** That is the log-layer mismatch of a coarse
  near-wall grid (Kawai & Larsson 2012), and S4c's $\Delta y=0.025\,\delta$ removes it. From the 4th cell up every run
  is within 2.5 % of the DNS.
- **The budget does not close on the means.** No subgrid model runs (the hybrid KEP-WENO scheme is the implicit LES),
  so with the eddy viscosity the explicit stresses carry 0.82 of $\tau_w$ through the first face, without it 0.41. The
  rest is the scheme's dissipation, or the correlation of the eddy viscosity with the strain, which the means cannot
  separate: whether the resolved and the modelled stress overlap below $h$ is not measurable from these fields.
- **A DNS-like incoming layer does not move the onset.** S4c's layer is the DNS's from the second cell up, and its
  onset is S4b's ([above](#a-refined-wall-grid-under-the-interaction)). S4e moves the onset with a layer slower than
  the DNS's. So matching the model to the resolved cells better, by itself, is not the lever for the interaction. The
  DNS's onset comes with its bubble, and the separation is what the runs lack (S4d opened a quarter of it,
  [above](#the-wall-model-with-the-pressure-gradient)).

## The pressure gradient without the convection, with it, and without the eddy viscosity

Two more restarts of S4b at $t=1200$, run to $t=1450$, statistics over $t\ge1300$ (`sbli_oracle.py compare --t0 1300`,
`fields`, `compression`), as S4d and S4e:

| Setup | S4g | S4f |
|---|---|---|
| model | S4d's (`wall_model_pressure_gradient`, filter 10) with `wall_model_pressure_gradient_ypg = 4`: the Bernoulli-consistent variant ([turbulence](../../turbulence#bernoulli-consistent-variant)) | S4d's, without `wall_model_eddy_viscosity` |
| question | does the consistent model keep a bubble? (a priori $1.0\,\delta_0$ against S4d's 2.8) | does S4d's bubble grow when the near-wall flow is slower (S4e's layer)? |
| run | one node (4 A30), 0.97 s per step (S4d 0.92, S4e 0.68), 24 checkpoints | one node (4 A30), 3 h 52 min, 22 checkpoints |

![The pressure-gradient runs against S4b and S4e](/flume/sbli-pg-runs.png)

*Left: the mean wall pressure of S4b, S4e, S4d, S4g and S4f against the DNS and the experiment; right: $C_f$ from the
wall model's $\tau_w$, the reversed stretches shaded. The legend's $L_{int}$ is `compare`'s, 0.01 below the table's.
`make_doc_figures.py --only sbli-pg-runs`.*

| | S4b | S4e, no $\nu_t$ below $h$ | S4d, $dp/dx$ | S4g, Bernoulli-consistent | S4f, $dp/dx$, no $\nu_t$ below $h$ | DNS, $Re_\theta=6882$ |
|---|---|---|---|---|---|---|
| $L_{int}/\delta_0$ (the time-mean field) | 3.16 | 3.51 | 3.13 | 3.24 | 3.37 | 3.30 |
| $L_{sep}/\delta_0$ | 0 | 0 | 0.60 | **0** | 0.14 ($x=58.4$–58.6) | 2.16 |
| $C_f$ minimum ($10^{-3}$) | 0.21 | 0.27 | $-1.87$ | 0.15 | $-0.02$ | |
| 5 % of the rise, ahead of $x_{imp}$ | $2.34\,\delta_0$ | $3.55\,\delta_0$ | $2.12\,\delta_0$ | $2.44\,\delta_0$ | $3.07\,\delta_0$ | $4.13\,\delta_0$ |
| 50 % of the rise | $1.96\,\delta_0$ | $2.64\,\delta_0$ | $1.77\,\delta_0$ | $2.04\,\delta_0$ | $2.42\,\delta_0$ | $3.60\,\delta_0$ |
| first-cell Mach number at $x_{imp}-5\,\delta_0$ | 0.98 | 0.75 | 0.94 | 0.95 | 0.70 | 0.83 |

- **The consistent model does not separate.** With convection balancing the pressure gradient above $y_{pg}$ (a few
  viscous lengths of the gradient, far below $h$), $\tau_w$ stays positive: S4g's $C_f$ dips to $0.15\cdot10^{-3}$, S4b's
  0.21, and its interaction is S4b's (onset within $0.1\,\delta_0$, $L_{int}$ within 0.08). S4d's bubble was the
  pressure gradient acting over the whole layer, the inconsistency Hickel et al. (2012) point out; a priori the
  consistent model kept $1.0\,\delta_0$ of it, the run none.
- **A slower near-wall flow does not grow S4d's bubble.** Without the eddy viscosity (S4f) the first cell slows to
  $M=0.70$, as in S4e, and the bubble shrinks from 0.60 to $0.14\,\delta_0$ and moves $0.7\,\delta_0$ downstream, toward the
  impingement ($x=58.4$–58.6 against S4d's 56.6–57.7); its $C_f$ also dips to $0.15\cdot10^{-3}$ at $x\approx55.4$, where
  its compression starts, without reversing; the onset moves upstream as in S4e, by a little less ($3.07\,\delta_0$ against 3.55 at
  5 %). The two levers do not add: each moves its own feature, and together they move both less.
- **What none of the five runs does.** The DNS separates over $2.16\,\delta_0$ from $x^*\approx-0.8$; the runs open at
  most $0.6\,\delta_0$ and only with the inconsistent model. The sample sits at $h\approx0.12\,\delta$ and $y_{pg}$
  about ten times lower ($\approx0.02\,\delta_i$ under the interaction's gradient); where the reversed layer is thinner
  than $h$ the sample still moves downstream, and an equilibrium-type model, which takes $\tau_w$ from that sample and a
  consistently balanced gradient, cannot see the reversal. The DNS's bubble is taller than $h$ in its large state, and the runs' own mean
  reversed layer lies below $h$ ([below](#the-reversed-flow-the-dns-s-bubble-and-the-runs-below-the-sample)). The convection the reversal needs is the flow's own, which a model fed with the resolved
  field (the integral model of Yang et al. 2015, the thin-boundary-layer equations with convection) carries and these
  do not.

## The reversed flow: the DNS's bubble and the runs' below the sample

Whether a reversed layer can hide below the matching height: how tall the DNS's bubble is, and where the runs' own
reversed flow sits.

**The DNS.** Bernardini et al. (2022) give the bubble's length ($2.16\,\delta_0$ in the mean, from $x^*\approx-0.8$
to $-0.15$) and its time history, but no height. Their Fig. 10 shows two instantaneous states on a vertical slice, in
$y^*=y/L_{int}$: in the large one the $u=0$ line reaches $y^*\approx0.14$–0.16, about $0.5\,\delta_0$, over
$x^*\approx-0.85$ to 0.1; in the other the bubble is almost absent. The area of the bubble swings between two
probable values with intermittent large excursions (their Fig. 9). Our sample sits at $h=0.208\,\delta_i=0.115\,\delta_0
=0.035\,L_{int}$: the large state is about $4h$ thick, ten of S4b's cells.

**The runs.** `sbli_oracle.py reversal` reads each checkpoint in $x_{imp}-4.5\,\delta_0<x<x_{imp}+\delta_0$, $y<1.5\,\delta_0$:

| | S4b | S4c, 2:1 box | S4d, $dp/dx$ | S4e, no $\nu_t$ below $h$ | S4f, $dp/dx$, no $\nu_t$ | S4g, Bernoulli-consistent |
|---|---|---|---|---|---|---|
| checkpoints | 71 | 21 | 24 | 23 | 22 | 24 |
| area with $u<0$ per unit span, mean (max), $\delta_0^2$ | 0.068 (0.142) | 0.069 (0.084) | 0.049 (0.086) | 0.140 (0.200) | 0.113 (0.191) | 0.077 (0.106) |
| highest reversed cell of a column, 99th percentile, mean (max), $\delta_0$ | 0.31 (0.44) | 0.30 (0.40) | 0.31 (0.44) | 0.38 (0.52) | 0.38 (0.52) | 0.34 (0.49) |
| first cells reversed | 9.5 % | 14.5 % | 6.1 % | 23.3 % | 17.9 % | 9.8 % |
| third cells (the sample) reversed | 5.1 % | 8.7 % | 3.9 % | 8.2 % | 6.6 % | 5.9 % |
| **time-mean $u<0$ in the first cell** | none | $x=57.8$–58.9, $0.66\,\delta_0$ | none | $x=56.8$–58.8, $1.13\,\delta_0$ | $x=57.5$–58.9, $0.77\,\delta_0$ | none |
| time-mean $u<0$ in cells 2, 3 | none | cell 2 $0.32\,\delta_0$; cell 3 none | none | none | none | none |
| $L_{sep}$ from the wall model's $\tau_w$ (`compare`) | 0 | 0 | 0.60 | 0 | 0.14 | 0 |

![The reversed flow of the runs](/flume/sbli-reversal.png)

*Left: the area with $u<0$ per unit span per checkpoint; right: the fraction of the first and of the third cells (the
wall model's sample) reversed in the window. `sbli_oracle.py reversal`, `make_doc_figures.py --only sbli-reversal`.
The time-mean reversal of the cells is printed by `sbli_oracle.py compression`.*

- **The runs reverse as high as the DNS's large state, but rarely.** In every run some columns reverse up to
  $0.4$–$0.5\,\delta_0$; the span-averaged area stays at 0.02–$0.2\,\delta_0^2$, within about twice its mean in each run (S4b's 71
  checkpoints: 0.04–0.14 around 0.068), where the DNS's bubble area has a bimodal distribution with a long tail of
  large excursions (their Fig. 9; its normalisation is not given, and a single slice is not a span average, so the
  areas themselves are not compared). The sample sees these reversals: its cell is reversed 35–64 % as often as the first cell.
- **In the mean, three runs separate below the sample.** S4c's resolved flow reverses in the mean over $0.66\,\delta_0$
  in its first cell and $0.32\,\delta_0$ in its second; S4e's over $1.13\,\delta_0$ ($x^*=-0.49$ to $-0.19$, against the
  DNS's $-0.8$ to $-0.15$), S4f's over $0.77\,\delta_0$; the third cell, the sample, never. The wall model's $\tau_w$, which
  takes the sample, stays positive there: `compare` reports no separation for S4c and S4e, and $0.14\,\delta_0$ for S4f.
- **Under that layer the wall stress has the wrong sign.** A positive $\tau_w$ takes streamwise momentum out of a first
  cell whose mean flow is already reversed, where the wall would push it forward. The model is consistent only while the
  reversed layer reaches $h$.
- **So the bubble is partly there, and the oracle's $L_{sep}$ misses it.** The mean reversal of the first cell is the
  separation a wall-resolved run would report; on it, S4e is at half the DNS's length and the onset of its reversal is at
  $x^*=-0.49$, against $-0.8$. What the equilibrium-type models lack is a wall stress consistent with a reversed layer
  below $h$: a model fed with the flow's own convection between the wall and $h$ (thin-boundary-layer equations, an
  integral model) has it by construction.

### A priori: the wall model fed with the resolved convection

The thin-boundary-layer balance between the wall and $h$ keeps the convection with the pressure gradient (Wang &
Moin 2002; Park & Moin 2014): $(\mu+\mu_t)\,dU/dy=\tau_w+\int_0^y F\,dy'$, $F=dp/dx+\partial(\rho uu)/\partial x+
\partial(\rho uv)/\partial y$. Fed with $F$ of the resolved cells below the sample (linear in $y$ between the cell
centres, $F=dp/dx$ at the wall, where the convection vanishes) and applied to the mean field of a run
(`sbli_oracle.py convective`):

| $C_f<0$ a priori, $\delta_0$ | equilibrium | $dp/dx$ alone (S4d's) | Bernoulli-consistent (S4g's) | fed with the resolved convection | mean $u<0$ of the first cell |
|---|---|---|---|---|---|
| S4b's mean field | 0 | 2.86 | 1.01 | 3.04 ($x=56.5$–61.9) | none |
| S4d's | 0 | 2.88 | 0.74 | 2.79 | none |
| S4e's | 0 | 2.57 | 1.31 | 3.33 ($x=56.4$–62.5) | 1.08 |
| S4f's | 0 | 3.11 | 1.31 | 2.98 | 0.72 |

![The wall model a priori on the mean fields](/flume/sbli-convective.png)

*$C_f$ of the four models on the mean fields of S4b (left) and S4e (right), the mean $u$ of the first cell dotted.
`sbli_oracle.py convective`, `make_doc_figures.py --only sbli-convective`.*

- **Where the resolved flow obeys Bernoulli, the model follows it.** At the steepest gradient ($x\approx56.5$) the
  resolved convection of cells 1–6 balances 0.9–1.4 of $-dp/dx$, and the model fed with it loses the spike of the
  pressure gradient alone ($C_f=-1.58\cdot10^{-3}$ there for S4b), as the Bernoulli-consistent model does.
- **Where the resolved near-wall flow is slow or reversed, the convection is small** (0 to 0.3 of $-dp/dx$ in cells
  1–3 over $x=57.5$–58.5) and the model behaves as the pressure gradient alone: $C_f<0$ over about $3\,\delta_0$, from
  $x^*\approx-0.6$, downstream into the reattachment.
- **The a priori does not tell this model from S4d's** (2.8–3.3 against 2.6–3.1), and S4d's run kept a fifth of its
  a priori bubble. Unlike S4d's model, this one changes with the resolved flow: its $\tau_w$ turns negative where the
  resolved convection stops balancing the pressure gradient, which is where S4c, S4e and S4f reverse below $h$. Only a run says
  whether the feedback opens the DNS's bubble.

### The run: S4h

S4b restarted at $t=1200$ with S4d's model fed with the resolved convection (`wall_model_pressure_gradient`, filter 10,
`wall_model_convection`; [turbulence](../../turbulence#fed-with-the-resolved-convection)), the eddy viscosity below $h$
kept, to $t=1450$ on one node (4 A30, 5 h 27 min), statistics over the 24 checkpoints of $t\ge1300$; it is the black
line of [the figure above](#the-pressure-gradient-without-the-convection-with-it-and-without-the-eddy-viscosity).

| | S4h, fed with the resolved convection | S4d, $dp/dx$ | S4g, Bernoulli-consistent | S4b | DNS, $Re_\theta=6882$ |
|---|---|---|---|---|---|
| $L_{int}/\delta_0$ (the time-mean field) | **3.36** | 3.13 | 3.24 | 3.16 | 3.30 |
| $L_{sep}/\delta_0$ from $\tau_w$ | 0.18 ($x=55.96$–56.29, $x^*=-0.66$ to $-0.61$) | 0.60 ($x^*\approx-0.6$ to $-0.4$) | 0 | 0 | 2.16 ($x^*=-0.8$ to $-0.15$) |
| $C_f$ minimum ($10^{-3}$) | $-0.30$ at $x=56.1$ | $-1.87$ | 0.15 | 0.21 | |
| time-mean $u<0$ in the first cell | $x=58.3$–59.0, $0.41\,\delta_0$ | none | none | none | |
| 5 % of the rise, ahead of $x_{imp}$ | $2.63\,\delta_0$ | $2.12\,\delta_0$ | $2.44\,\delta_0$ | $2.34\,\delta_0$ | $4.13\,\delta_0$ |
| 50 % of the rise | $2.23\,\delta_0$ | $1.77\,\delta_0$ | $2.04\,\delta_0$ | $1.96\,\delta_0$ | $3.60\,\delta_0$ |
| first-cell Mach number at $x_{imp}-5\,\delta_0$ | 0.91 | 0.94 | 0.95 | 0.98 | 0.83 |
| reversed area per unit span, mean (`reversal`), $\delta_0^2$ | 0.091 | 0.049 | 0.077 | 0.068 | |
| **$C_f$ upstream against S4b**: $x=15$–25, 25–35, 35–45, 45–52 | **+5, +6, +9.5, +12 %** | +1, +1, −2, −2 % | 0, 0, −3, −3 % | | |

- **The interaction moves toward the DNS.** $L_{int}$ is the DNS's within 2 % (the closest of the runs), the onset moves
  upstream by $0.3\,\delta_0$ at 5 % and at 50 % of the rise, the first cell slows from $M=0.98$ to 0.91 (the DNS profile
  0.83), and the wall stress reverses at $x^*\approx-0.63$, where the DNS's bubble begins. The bubble is short, $0.18$
  against $2.16\,\delta_0$; the resolved first cell reverses in the mean over another $0.41\,\delta_0$ further down.
- **The zero-gradient layer moves too, which it should not.** Upstream of the interaction $C_f$ rises by 5 % near the
  inflow and by 12 % ahead of the interaction, where S4d and S4g stay within 3 % (S4e and S4f, without the
  eddy viscosity, within 7 %), and the recovery behind the
  reattachment is higher still. The convection fed to the model is the instantaneous resolved one, and its time mean
  holds the gradient of the resolved shear stress, $\partial\langle\rho u'v'\rangle/\partial y$, large in the cells below $h$
  ([the budget above](#the-near-wall-cells-of-the-incoming-layer): 0.2 of $\tau_w$ through the first face, about 1 by the
  third). The model's own eddy viscosity already carries that stress: it is counted twice, and $\tau_w$ grows. Park &
  Moin (2014) meet the same double count in their thin-boundary-layer model fed with the LES and reduce the model's
  eddy viscosity dynamically, so that modelled and resolved stresses at $h$ add up to the total.
- **So the convection is a lever, and this form of it is not yet a model.** It moves both the onset and $L_{int}$ the
  right way, which no equilibrium-type variant did, but a model that changes the zero-gradient layer by 10 % fails the
  layer before the interaction. The correction is the reduced eddy viscosity (or the convection of the mean flow only),
  checked on the zero-gradient boundary layer before the interaction.

## Not done

- A separation of the right length and place. The equilibrium-type models do not give it: the pressure gradient alone opens $0.6\,\delta_0$ by an inconsistent balance, the consistent one none, the slower near-wall layer $0.14\,\delta_0$ ([above](#the-pressure-gradient-without-the-convection-with-it-and-without-the-eddy-viscosity)). The model fed with the resolved convection (S4h) moves the onset and $L_{int}$ the right way and opens $0.18\,\delta_0$ where the DNS's bubble begins, but counts the resolved stress twice ($C_f$ of the zero-gradient layer +5 to +12 %): left, its reduced eddy viscosity (Park & Moin 2014), checked on the zero-gradient layer first; and a matching height fixed in $\delta_i$ across the levels.
- The late compression ($x^*\approx-0.6$ half-way against the DNS's $-1.1$), unchanged by the near-wall grid and by the pressure-gradient term; moved 40 % of the way by taking out the eddy viscosity below $h$ (S4e), which also moves the incoming layer off the DNS's: the matching between the model and the resolved near-wall cells.
- The wall-pressure fluctuations and the low-frequency unsteadiness.
- The recycling–rescaling inflow (S5).
