# Shock/boundary-layer interaction

| ID | Milestone | Type | Legs | Backends | Figure |
|---|---|---|---|---|---|
| — (campaign `tbl_cases_s4.txt`) | [#121](https://github.com/szaghi/adam/issues/121) · S4b–S4d | P | none (the pressure-gradient key: wm12) | FNL | yes |

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
gives $L_{int}=x_{imp}-x_{ref}$, the length STREAmS scales $x$ with.

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
  *downstream*.
- **For the next wall model,** non-equilibrium: the target is the DNS's bubble ($x^*\approx-0.8$ to $-0.15$) and its
  wall pressure, half-way at $x^*\approx-1.1$. The run's wall pressure is printed by `sbli_oracle.py compare` (at fixed
  $x^*$, against both DNS) and `compression` (where it reaches 5, 10 and 50 %); the oracle reports the run's $L_{sep}$,
  the DNS's is quoted from Bernardini et al. (2022).

## Not done

- A separation of the right length and place: the pressure-gradient term opens a short bubble ahead of the impingement; a non-equilibrium wall model (convection), and a matching height fixed in $\delta_i$ across the levels.
- The late compression ($x^*\approx-0.6$ half-way against the DNS's $-1.1$), unchanged by the near-wall grid and by the pressure-gradient term: the missing separation, above.
- The wall-pressure fluctuations and the low-frequency unsteadiness.
- The recycling–rescaling inflow (S5).
