# SOTA: telescoping fully-conservative schemes at 2:1 AMR jumps with collocated cell-centered variables (FD + FV)

**Date:** 2026-07-02 · **Method:** deep-research harness — 5 search angles (foundations, high-order FV/FD, production frameworks, cell-centered div-B, contrarian/limitations), 20 sources fetched, 82 claims extracted, top 25 adversarially verified (3 independent refutation votes each). **Result: 25/25 confirmed, 0 refuted.**
**Scope:** collocated cell-centered FD and FV only. Staggered/Yee CT, face-centered B, edge-E registers, mimetic, DG/FEM excluded by construction.

---

## 1. The one-paragraph answer

Discrete telescoping conservation (Σ fine-face fluxes = coarse-face flux) at a 2:1 jump is **not emergent from any cell-centered discretization; it is imposed** — by Berger–Colella refluxing (flux registers) or its same-timestep equivalent (direct flux replacement). This is exact **by construction, at any flux order**, for FV. For **cell-centered FD it is achievable only when the scheme is written in flux (conservative) form**, coarse/fine patches are conformed at **flux points** (not solution points), and the reflux correction is applied; the 2025 state of the art (Wang et al., arXiv:2511.08335) demonstrates this and documents its price: **once flux correction activates, global conservation-error convergence drops to 2nd order** ("any means to preserve the conservation will reduce the schemes to no more than second order"). For div-B with cell-centered storage and no staggering, the surveyed SOTA is **GLM/Dedner cleaning** — the induction equation stays in flux form (so reflux applies to B), but div-B is **transported and damped, never machine-zero**. No surveyed method achieves machine-precision solenoidality at a seam with purely cell-centered storage.

---

## 2. Verified mechanism landscape

### 2.1 FV cell-centered — Berger–Colella reflux (canonical, exact, production-standard) — HIGH confidence, 3-0

Fine-level fluxes are accumulated on the coarse-fine faces during (subcycled) fine advances; at sync, the coarse cells adjacent to the interface receive the coarse-level divergence of the mismatch `F_coarse − Σ(area/time-weighted F_fine)`. The effective coarse-face flux then *equals* the fine sum — telescoping is restored identically, independent of the flux's spatial order.

- **AMReX** `FluxRegister`: "accumulates and ultimately stores the net difference in fluxes between the coarse grid and fine grid advance over each face"; defect = "area-weighted fluxes from the fine grid advance do not in general match the underlying flux from the coarse grid face"; remedy modifies "coarse cells immediately next to the coarse-fine interface" via "the coarse-level divergence of the flux register data." `YAFluxRegister` lifecycle: `reset → CrseAdd → FineAdd → Reflux`.
  https://amrex-codes.github.io/amrex/docs_html/AmrCore.html · https://amrex-codes.github.io/amrex/doxygen/classamrex_1_1YAFluxRegister.html
- **Origins:** Berger & Colella 1989 (JCP 82:64) §4; reflux-divergence form Martin & Colella 2000.
- **PLUTO-AMR** (Mignone et al. 2012, ApJS 198:7, **Eqs. 48–49**): coarse flux at a shared edge replaced by the space-time average of fine fluxes, coarse cells corrected accordingly.
  https://iopscience.iop.org/article/10.1088/0067-0049/198/1/7
- **Athena++** (Stone, Tomida, White & Felker 2020, ApJS 249:4): **same-timestep flux replacement** — no subcycling, "the flux used to update the coarse cell on the face that overlaps with the fine cells is simply replaced with the area-weighted sum of the fluxes from these four fine cells" so that "the area integral of the fluxes … on cell faces at the boundaries between MeshBlocks on different levels must be exactly equal." Identical guarantee, different cadence. **NB: this covers the cell-centered conserved variables only — Athena++'s B is staggered CT, out of scope.**
  https://iopscience.iop.org/article/10.3847/1538-4365/ab929b

### 2.2 High-order FV at the jump — McCorquodale–Colella + the constrained ghost-fill — HIGH, 3-0

- **McCorquodale & Colella 2011** (Comm. Appl. Math. Comput. Sci. 6(1)): 4th-order space+time FV on locally refined Cartesian grids; flux registers on the l/l+1 boundary faces carry conservation.
- **Load-bearing subtlety** (Guzik dissertation, escholarship qt8sd6j4fc): nominally 4th-order polynomial ghost-fill is **NOT sufficient** for observed 4th-order convergence — the interpolation must be **constrained** so that Σ⟨U⟩ over the interpolated fine cells = (n_ref)^D ⟨U⟩ of the underlying coarse cell (constrained least-squares cubic). I.e. even the *ghost-fill* needs a conservation constraint to keep the order, separate from the reflux.
  https://escholarship.org/content/qt8sd6j4fc/qt8sd6j4fc_noSplash_d30ae33354338cef900c92f1c624319f.pdf
- **Mapped/curvilinear grids** (Guzik–McCorquodale–Colella, OSTI 1034480): freestream preservation at the jump additionally requires **metric telescoping** — the coarse-face area vector N_s overwritten by the exact sum of fine-face N_s on shared edges, with both tangential metric components derived from that single value. A second, purely geometric telescoping identity layered on top of flux reflux.
  https://www.osti.gov/servlets/purl/1034480

### 2.3 Cell-centered FD — possible, but only in flux form, and with a hard order cap — HIGH, 3-0

**Wang, Zeng, Deiterding, Yang, Liang (Nov 2025), arXiv:2511.08335** — cell-centered FD-WENO on AMReX, ratios 2 and 4. The sharpest primary source for the whole question:

1. "any high-order interpolation in the finite-difference framework is **inherently non-conservative**" → ghost-fill (accuracy) and reflux (conservation) are **decoupled concerns**;
2. "solution points are placed at the cell centers … the domain boundaries are **flux points rather than solution points**" — coarse/fine patches conformed at flux points → "inter-patch conservation to be rigorously enforced";
3. "the coarse grid solution is corrected by an equivalent flux difference δF_{I+1/2} at the end of a coarse step"; under explicit RK "the flux differences at coarse-fine interfaces at each stage must be **scaled and accumulated**";
4. ghost-fill is **hybrid**: non-conservative 5th-order WENO interpolation in smooth regions, 2nd-order **conservative** interpolation at shocks, switched by a scale-irrelevant approximate-Riemann troubled-cell detector — because pointwise replacement "would … violate the mass conservation on the coarse grid" with shocks at the interface;
5. **THE LIMITATION:** enforcing conservation via flux correction "limits the global conservation error convergence to **second order**" — "any means to preserve the conservation will reduce the schemes to no more than second order."

Demonstrated (Table 6 convergence + §6.2 shock-interface tests), not asserted. *Caveat: preprint, not yet peer-reviewed; its "flux points" are auxiliary interface storage for the register — ordinary flux-register machinery, not Yee staggering.*
https://arxiv.org/abs/2511.08335

**Corroborating limitation statement** (contrarian-angle source, preview-level): for point-value FD-WENO on AMR it is "difficult to design a scheme that is BOTH conservative AND high order (above second order)" due to coarse/fine mass inconsistency in the FD formulation. And the AFD-WENO line (arXiv:2403.01266 angle) shows pointwise FD *can* be written in exact flux-conservation form — flux form is the enabling property, not cell-centering per se.

### 2.4 The register-free alternative — Freret–Groth heterogeneous blocks — MEDIUM, 3-0

**Freret, Ivan, De Sterck & Groth 2019** (J. Sci. Comput. 79:176–208), high-order CENO cell-centered FV MHD at 2:1 jumps: "the ghost cells of a block are stored directly at the resolution of the neighboring blocks" (hanging nodes), no high-order restriction/prolongation operators, and "there is no need for an additional correction to enforce the flux conservation properties … automatically satisfied by the non-uniform treatment." **Medium confidence:** the discrete proof is deferred to a prior Freret–Groth reference, single-group, not independently corroborated. Conceptually: conservation by making both sides literally evaluate the same face flux at the same resolution, instead of correcting a mismatch after the fact.
https://link.springer.com/article/10.1007/s10915-018-0844-1

### 2.5 div-B / curl-type equations with cell-centered storage, no staggering — HIGH, 3-0

- **SOTA = GLM/Dedner cleaning.** Mignone, Tzeferacos & Bodo 2010 (JCP 229, arXiv:1001.2832): 3rd/5th-order **conservative FD** GLM-MHD, "all primary flow variables discretized at the zone center"; div-B via "a generalized Lagrange multiplier yielding a mixed hyperbolic/parabolic correction"; explicitly avoids elliptic cleaning AND "the additional complexities required by staggered mesh algorithms." Mignone & Tzeferacos 2010 (JCP 229:2117, arXiv:0911.3410) — same for CTU FV: "cell-centered spatial collocation (no staggered mesh) of all flow variables, including the magnetic field … fully conservative in mass, momentum, magnetic induction and energy."
  https://arxiv.org/abs/1001.2832 · https://arxiv.org/pdf/0911.3410
- **The distinction that matters:** "conservative in magnetic induction" = the B-update is in flux form (⇒ refluxing applies to the B rows like any conserved row). It does **NOT** mean div-B = 0 to machine precision. GLM **transports and damps** divergence errors; it never zeroes them. No surveyed cell-centered method achieves machine-zero div-B at a seam.
- Freret et al. 2019 independently choose GLM for the same reasons (cell-centered, AMR-compatible).
- A non-staggered *alternative* exists at preview level (arXiv:2409.14992, Lagrangian-relaxation/entropy approach) but is not a solenoidality-to-round-off mechanism either.

---

## 3. Direct answers to the posed sub-questions

**What guarantees exact telescoping, and at what order?** The reflux correction itself — it is definitional, not emergent. FV: exact at any flux order (the correction is exact algebra on face-integrated fluxes; the *solution* order at the interface is a separate ghost-fill concern, cf. McCorquodale–Colella's constrained interpolation). FD: exact only for flux-form schemes conformed at flux points; the *conservation-error* convergence is then capped at 2nd order (Wang 2025).

**Can centered, non-flux-form FD be made telescoping at the seam?** **No published treatment found.** Every in-scope source presupposes flux/divergence form; the pointwise-FD lines that do achieve exact conservation (AFD-WENO) do it precisely by exhibiting an exact flux form. Flux form is the de-facto prerequisite. This is an absence of construction, **not an impossibility theorem** — none was found.

**Impossibility/limitation results (the strongest found):**
1. Wang 2025: conservation enforcement ⇒ conservation-error order ≤ 2 for cell-centered FD (empirical/constructive, not a theorem).
2. High-order interpolation for FD ghost-fill is inherently non-conservative (Wang 2025; consistent with the classical under-determination results for pointwise div-free prolongation).
3. McCorquodale–Colella/Guzik: unconstrained 4th-order ghost interpolation loses the observed order — the conservation constraint on the interpolant is necessary.
4. Mapped grids add a *geometric* telescoping requirement (metric sums) on top of flux telescoping.

**Framework map (verified subset):** AMReX = FluxRegister/YAFluxRegister (subcycled accumulate + reflux). Chombo/Castro = same family. PLUTO-AMR = reflux Eqs. 48–49 + GLM for B. Athena++ = same-step flux replacement for cell-centered conserved vars (B is staggered CT — outside this scope). Flash-X/Spark, BATSRUS cell-centered mode: named in scope but **not directly sourced this pass** (open question), though Spark's flux-register usage is documented in Couch 2021 (prior #13 review).

---

## 4. Open questions (from the verification pass)

1. Is the 2nd-order conservation-error cap fundamental for collocated FD, or specific to Wang's reflux construction? Any construction with strict telescoping AND >2nd-order conservation error?
2. What is the actual discrete proof behind Freret–Groth's "automatic" register-free conservation across hanging-node faces?
3. Any genuinely centered (non-flux-form) FD scheme made telescoping at a seam? (All surveyed sources assume flux form.)
4. How exactly do BATSRUS cell-centered mode / Flash-X Spark / PLUTO reflux the GLM ψ and induction fluxes at seams, and what div-B bound do they achieve there?

---

## 5. Implications for ADAM/PRISM #13 Phase B — *my analysis, not verified findings*

1. **The rmf-amr-fd seam div(B) leak is exactly what the literature predicts.** Nothing in the cell-centered SOTA makes a centered-FD `div_h(curl_h)` identity survive a 2:1 jump; the matched-stencil cancellation is a uniform-grid property, and no surveyed mechanism restores it algebraically at the seam.
2. **The Phase-B acceptance "pointwise ∇·B ≤ 1e-13 at seam cells" is beyond the surveyed cell-centered SOTA.** Production cell-centered MHD (PLUTO, Mignone FD-GLM, Freret CENO) accepts truncation-order seam div-B and *cleans* it (GLM). Machine-zero at the seam is delivered only by staggered CT — excluded by design here. **The #13 §2.1 premise ("cell-centered + interface machinery ⇒ algebraic ∇·B=0, demonstrated by Flash-X/BATSRUS/Athena++") is not supported by the verified sources**: Athena++'s telescoping FluxCorrection covers the hydro conserved variables while its B is staggered CT; the genuinely cell-centered production codes use cleaning, which is not algebraic zero.
3. **On the §6 decision:** the literature's answer for the *conservation* half is option (a)'s shape — flux-form seam operator + reflux (already largely built: M0–M4) — with the caveat that conservation enforcement caps that error's order at 2. For the *divergence* half, the literature operates in option (b)'s frame: treat seam div-B as truncation, verify convergence under refinement, and control it with cleaning (PRISM already has `divergence_correction = hyperbolic` = GLM-style). The unrun convergence-under-refinement measurement is precisely the test the field's framing demands.
4. **Two concrete mechanisms worth evaluating against PRISM's architecture** if a stronger-than-cleaning seam treatment is still wanted: (i) Wang-2025's flux-point conforming + δF for the FD path (makes the FD seam conservative, does not zero div-B); (ii) Freret–Groth native-resolution ghosts (both sides evaluate the same face flux — removes the mismatch at the source; heavier restructuring of the ghost machinery, overlaps Phase C).

---

*All quoted claims passed 3-0 adversarial verification against the cited primary sources. arXiv:2511.08335 is a preprint — treat its quantitative claims as strong-but-provisional.*
