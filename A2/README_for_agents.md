# README for AI agents

You are helping Group 22 with assignment A2 of a Classical Mechanics course.
Read this file fully before changing anything in `A2/`.

## Context

- Assignment spec: `../Project01student.md` (A1 spec, describes A2 goals).
- Group's A1 design report: `../力學_project1_G22.md`. It defines Experiments 1–4
  and the Δt-halving check. **The experiment design is the group's decision; do not
  redesign it.** You may add diagnostics that support it.
- Course toolkit: `../numerical_tools_student_release/student_release/`.
  `mechanics_integrators.py` is the primary source for how the integrators work.

## Hard rules

1. **Never edit** `mechanics_integrators.py` or anything under
   `numerical_tools_student_release/`. Import it via `common.py` (which adds it to
   `sys.path`).
2. Use the toolkit's `rk4` and `velocity_verlet` wherever their interfaces apply.
   The toolkit's VV accepts only `a(t, x)` — no velocity dependence. For damped
   cases, use the group's own variants in `common.py`:
   `vv_damped_lagged` (naive, lagged velocity) and `vv_damped_implicit`
   (velocity update solved exactly for linear damping). Keep both: the contrast
   between them is a reported finding ("option 3": the toolkit limit is a result).
3. All physical parameters and default Δt live in `common.py`. Do not hard-code
   different values inside experiment scripts unless the experiment needs a
   deliberately different value; if so, name it and print it in the summary.
4. Every experiment script must:
   - be runnable standalone (`python expN_*.py`) and via `run_all.py`;
   - use the non-interactive `Agg` backend and save PNGs to `figures/`;
   - write a plain-text summary to `results/expN_summary.txt` (numbers the group
     will quote in the report);
   - run at Δt and Δt/2 and report the ratio of numerical errors (A1 section 4).
5. Figure text in English (Chinese fonts may be missing). Series colors: use
   `common.COLORS` in fixed order (blue, orange, aqua); exact/analytic references
   are dashed neutral gray. Linewidth ≈ 1.5–2. One y-axis per panel; never dual axes.
6. Report results honestly. If a result contradicts a prediction in the A1 report,
   say so in the summary text — do not tune parameters to force the prediction.
7. No new dependencies beyond numpy and matplotlib (no scipy). The complete
   elliptic integral for the exact nonlinear period is computed with the AGM in
   `common.py`.
8. Do not commit or push unless a group member asks.

## Conventions

- State: `x` = θ (rad), `v` = ω (rad/s); per unit mass; g = 9.8, L = 1.
- Acceleration signatures follow the toolkit: `a(t, x, v)` for RK4, `a(t, x)` for VV.
- Energy: `common.energy(theta, omega, model)` with model `"nonlinear"` or `"linear"`.
- Periods: `common.measure_periods(t, theta)` (upward zero crossings, linear interp).

## File map

| File | Role |
|---|---|
| `common.py` | params, models, energy, exact solutions, period tools, damped VV variants, plotting style |
| `exp1_amplitude.py` | Exp 1 — amplitude effect on period, no damping, RK4 |
| `exp2_energy.py` | Exp 2 — energy error / phase space, VV vs RK4, no damping |
| `exp3_damped_accuracy.py` | Exp 3 — damped linear vs analytic; toolkit limit; convergence |
| `exp4_damped_nonlinear.py` | Exp 4 — damped, linear vs nonlinear, RK4 |
| `run_all.py` | runs all experiments |

## Verifying a change

Run `python run_all.py` from `A2/`. It must finish without errors and regenerate
every file in `figures/` and `results/`. Check convergence ratios in the summaries:
≈4 for second-order methods, ≈16 for RK4 (when not at round-off level).
