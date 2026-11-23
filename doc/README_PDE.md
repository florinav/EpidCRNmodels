# EpidCRN graph-diffusion extension

Three self-contained `.wl` files extend EpidCRN to multi-patch
(graph-Laplacian discretization of partly diffusive PDE) epidemic models.

## Files

| File | Purpose |
|------|---------|
| `GraphDiffuseRN.wl` | Defines `GraphDiffuseRN`; demonstrates 2-patch SIS (Task 4 test) |
| `NeumannLaplacian1D.wl` | Defines `NeumannLaplacian1D`; prints eigenvalues |
| `PatchNGM.wl` | Defines `PatchNGM`; verifies constant-mode equality for SIS and FLS |
| `tests/SEIRGraphTest.wl` | 2-patch SEIR: Task 4 with and without I-diffusion |
| `tests/FLSGraphTest.wl` | 3-patch FLS/cholera (B diffusing): Task 4 |
| `tests/HeterogeneousSISTest.wl` | Task 5: heterogeneous beta counterexample, limiting cases |

## Functions

### `GraphDiffuseRN[RHS, var, diff, L]`

Extends a local (one-patch) ODE to a multi-patch system on the graph with
Laplacian `L`.

- `RHS` — local RHS vector, length `nS`
- `var` — local variable list, length `nS`
- `diff` — `Association{species -> d_coeff}`, absent species get `d = 0`
- `L` — `nP × nP` graph Laplacian (row sums = 0, diagonal ≤ 0)

Returns `{varPatch, RHSPatch}`.

Variable naming: species `s` in patch `p` → `ToExpression["s" <> ToString[p]]`,
e.g. `s1`, `s2` for a 2-patch model with local variable `s`.

### `NeumannLaplacian1D[n, h]`

Returns the `n × n` finite-difference Laplacian for the 1D diffusion equation
`u_t = d u_{xx}` on a uniform grid with spacing `h` and homogeneous Neumann
(zero-flux) boundary conditions.

Row sums are 0; eigenvalues are `0, -2/h², -2/h²(1-cos(π/(n-1))), …`.

### `PatchNGM[FF, VV, infVarsLocal, diff, L]`

Constructs the patch next-generation matrix via Kronecker products:

```
F_patch = I_nP ⊗ FF
V_patch = I_nP ⊗ VV − L ⊗ D_inf
K_patch = F_patch · V_patch⁻¹
```

where `D_inf = Diagonal(d_1, …, d_k)` for the `k` infectious variables.

## Mathematical results verified

**Constant-mode theorem** (Task 4): For homogeneous coefficients and any
graph Laplacian `L`, the spectral radius of `K_patch` equals `R_0^local`.

This is the discrete analogue of `R_0^PDE = R_0^ODE` under homogeneous
Neumann boundary conditions.

**Heterogeneous counterexample** (Task 5): For a 2-patch SIS model with
`β₁ ≠ β₂`:

- `d → 0`:  `R_0^patch → max(R_0^local,1, R_0^local,2)`
- `d → ∞`:  `R_0^patch → (β₁ + β₂) / (2(μ + γ))` (arithmetic mean)
- Intermediate `d`:  `R_0^patch` is strictly between the two local values

so `R_0^patch ≠ R_0^local,p` for either patch.

## Usage in a model file

```mathematica
(* Load the function *)
<<GraphDiffuseRN`;   (* or copy the Module definition directly *)

var = {s, i};
RHS = {mu - mu*s - be*s*i, be*s*i - (mu+ga)*i};
L3  = NeumannLaplacian1D[3, 1];
{varPatch, RHSPatch} = GraphDiffuseRN[RHS, var, <|i -> d|>, L3];

(* Standard EpidCRN pipeline *)
{RN, rts, spe, alp, bet, gam, rnR} = ODE2RNs[RHSPatch, varPatch, prF];
{..., FF, VV, K, R0A, infV} = NGMRN[RN, rts, varPatch];
```
