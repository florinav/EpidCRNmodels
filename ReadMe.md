# EpidCRN Package - ReadMe

## Coding Conventions

**Variable Naming:**
- Use global names consistently across notebooks: RN, var, alp, bet, gam, jac, eig, par
- Never use numbered variants (RN1, RN2, etc.) - not necessary since we don't run multiple networks simultaneously
- Always use 3-letter abbreviations when printing: "RN=", "gam=", "eig="
- Always print variable name before value: `Print["RN=",RN," gam=",gam," eig=",eig]`
- Never use Print statements without variable names

**Function Outputs:**
- All important functions must output lists (not associations)
- Usage statements appear only in loader (EpidCRN.wl) and only once per function

**Testing:**
- Each new function/script must have uncommented test at end of subpackage file (before End[])
- Tests run automatically when package loads
- Follow pattern in Core.wl for ODE2RN test
- Provide self-contained tests for each new submodule

## Key Functions

**Boundary Analysis:**
- `bdFp[RHS, var, mSi]` - Computes boundary equilibria on siphon facets. For each minimal siphon in mSi: sets siphon variables to 0, solves the resulting system, filters out DFE and negative solutions using onlyNN. Returns list where `bdFp[[j]]` is either positive solutions (list of rule lists) or `{equations, variables}` for siphon j when no positive solutions exist.

**Bifurcation Analysis:**
- `scan[RHS, var, par, persRule, plotInd, mSi, gridRes, steadyTol, stabTol, chopTol, R01, R02, R12, R21, rangeExtension]` - Scans parameter space and classifies equilibria.
  - **Reproduction numbers**: R01 = R0 for strain 1 at DFE, R02 = R0 for strain 2 at DFE, R12 = invasion number for strain 1 into endemic E2, R21 = invasion number for strain 2 into endemic E1
  - **Classification logic**:
    - DFE: R01 < 1 and R02 < 1
    - E1: R01 > 1 and R02 < 1 (strain 1 only)
    - E2: R02 > 1 and R01 < 1 (strain 2 only)
    - Coexistence region: R01 > 1 and R02 > 1
      - If R12 < 1 and R21 < 1: competitive exclusion (winner determined by R01 vs R02)
      - If R12 > 1 and R21 < 1: strain 1 wins (E1)
      - If R12 < 1 and R21 > 1: strain 2 wins (E2)
      - If R12 > 1 and R21 > 1: potential coexistence - solve for equilibrium and check eigenvalues
  - **Performance optimization**: Stability (eigenvalue) analysis performed only in coexistence region (R01>1, R02>1, R12>1, R21>1) to reduce computation time. Other regions classified by reproduction numbers only.
  - **rangeExtension**: Parameter controlling plot range extension beyond intersection point (default 2.5 = 250%)

## Recent Work: Two-Reduction Regime Analysis (December 2024)

### Motivation
Integrated Feliu-Lax-Walcher-Wiuf (FLWW) asymptotic reduction theory with Saez-Feliu-Wiuf (SFW) algebraic elimination theory to validate QSS reductions via singular perturbation analysis.

### Key Theoretical Contributions

**FLWW Paper ("Quasi-steady state and singular perturbation reduction for reaction networks"):**
- Provides necessary/sufficient conditions for Tikhonov-Fenichel reduction validity
- **Blanket condition (i):** B₀(x) invertible ⟺ spanning trees exist in G_U (Lemma 11)
- **Blanket condition (ii):** Critical manifold stationary ⟺ cycle conditions (Lemma 12)
- **Corollary 13:** Setting inflow rates (∗→U reactions) to zero guarantees TFPV
- **Key insight:** Same eliminable set can be valid in different parameter regimes

**Comparison: SFW vs FLWW**
| Aspect | SFW (2017) | FLWW |
|--------|------------|------|
| Purpose | Identify eliminable sets | Validate QSS asymptotically |
| Conditions | Noninteracting, U-linear, spanning tree | + Blanket conditions (i), (ii) |
| Output | Reduced network structure | Parameter regimes where valid |
| Validation | Algebraic (steady state) | Asymptotic (Tikhonov-Fenichel) |

### Implementation: John.nb

**File:** `C:\Users\flori\Dropbox\EpidCRNmodels\multiStrain\John.nb`

Comprehensive single-file analysis of John's two-strain SIR model demonstrating:

**Model:** 5D system {s, i1, i2, r1, r2} with:
- Cross-immunity (σ₁, σ₂ parameters)
- Demographic turnover (Λ, μ rates)
- Strain-specific infection/recovery (β₁, β₂, γ₁, γ₂)
- Waning immunity (θ₁, θ₂)

**Maximal Eliminable Sets Found:**
1. {i1, i2} → reduces to 3D system in {s, r1, r2}
2. {r2, r1, s} → reduces to 2D system in {i1, i2}
3. {i1, r1} → reduces to 3D system

**Two Reduction Regimes Illustrated:**

*Regime 1: Fast Infection Dynamics*
- **Parameters:** β₁=2.0, β₂=1.5, γ₁=1.0, γ₂=0.8 (fast) vs Λ=0.05, μ=0.002, θ₁=0.05, θ₂=0.04 (slow)
- **Elimination:** {i1, i2} via QSS assumption (i₁' = 0, i₂' = 0)
- **Reduced system:** 3D in {s, r1, r2}
- **Reconstruction:** i₁ = φ₁(s, r1, r2), i₂ = φ₂(s, r1, r2)
- **Biological interpretation:** Epidemic on fast timescale → track slow demographic changes

*Regime 2: Fast Demographic Turnover*
- **Parameters:** Λ=2.0, μ=0.8, θ₁=1.0, θ₂=0.8 (fast) vs β₁=0.2, β₂=0.15, γ₁=0.1, γ₂=0.08 (slow)
- **Elimination:** {s, r1, r2} via QSS assumption (s' = 0, r₁' = 0, r₂' = 0)
- **Reduced system:** 2D in {i1, i2}
- **Reconstruction:** s = ψ₁(i1, i2), r₁ = ψ₂(i1, i2), r₂ = ψ₃(i1, i2)
- **Biological interpretation:** Demography on fast timescale → track slow epidemic spread

**File Structure:**
```mathematica
(* Section 1: MODEL DEFINITION *)
RHS = {...};
var = {s, i1, i2, r1, r2};
{RN, rts, spe, alp, bet, gam} = ODE2RN[RHS, var];

(* Section 2: BOUNDARY ANALYSIS *)
{RHS1, var, par, cp, mSi, Jx, Jy, cDFE, E0, K, R0A, infVars, ng} = bdAn[RN, rts, var];
{RHS, var, par, cp, mSi, Jx, Jy, cDFE, E0, K, R0A, infVars, EA} = bdCo[RN, rts, var];

(* Section 3: REDUCTION ANALYSIS *)
allMaxElim = maxElim[RN, rts];

(* Section 4: GRAPHICAL ILLUSTRATION - automatic if ≥2 sets *)
If[Length[allMaxElim] >= 2,
  U1 = allMaxElim[[1]];
  U2 = allMaxElim[[2]];
  (* Regime 1 comparison plots *)
  (* Regime 2 comparison plots *)
];
```

**Design Philosophy:**
- ONE file per model (named after the model)
- User changes RHS/var at top → everything updates automatically
- Automatically finds all maximal eliminable sets
- If ≥2 sets exist, automatically illustrates both regimes
- Verifies FSW conditions for each set
- Creates comparison plots: Blue = Full system, Red = Reduced system

**Graphical Output:**
- Slow variables: Full vs Reduced trajectories overlay
- Fast variables: Reconstructed via QSS vs Full system
- Validates asymptotic agreement in appropriate parameter regimes

### Related Files

**Analysis summaries:**
- `TwoReductionsSummary.md` - Detailed comparison of both regimes
- `illustrateReductions.nb` - Generic tool (works for any CRN)
- `illustrateTwoReductions.nb` - Earlier working prototype
- `compareReductionsFinal.nb` - Development version

**Recommendation for New Models:**
1. Copy John.nb to NewModel.nb
2. Change RHS and var definitions
3. Run → automatic two-regime analysis if multiple eliminable sets exist

---

## Current Plan

Integration of FLWW asymptotic theory with SFW elimination theory complete for John's model. Framework ready for application to other multi-strain/multi-pathway epidemic models.
