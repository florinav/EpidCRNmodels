# ReadMeEpi - epi.wl Package Documentation

## Core Principle: Species Format

**CRITICAL INVARIANT**: `spe == ToString[var]` must always be `True`

- `var`: List of symbol variables, e.g., `{m0, m1, m2, mc, M1, M2}`
- `spe`: List of string species names, e.g., `{"m0", "m1", "m2", "mc", "M1", "M2"}`
- **NEVER use ToLowerCase** - it destroys the correspondence (M1 → "m1" conflicts with m1 → "m1")

## Reaction Network Format (RN)

RN is a **list of Rules** where each element has format: `lhs -> rhs`

### Components of lhs and rhs:
- **Strings**: Species names like `"m0"`, `"M1"`, `"mc"`
- **Times**: Coefficients like `2*"m0"` (2 molecules of m0)
- **Plus**: Sums like `"m0" + "m1"` (m0 and m1 together)
- **Zero**: `0` (empty side - for birth/death reactions)

### Examples:
```mathematica
0 -> "m0"              (* Birth: nothing → m0 *)
"m0" -> 0              (* Death: m0 → nothing *)
"m0" + "m1" -> "m1"    (* Catalysis: m0 consumed, m1 catalyzes *)
"m0" + "m1" -> 2*"m1"  (* Infection: m0 infected by m1 *)
"m0" -> "M1"           (* Activation: m0 becomes M1 *)
```

## ODE2RN Function

**Input**: `ODE2RN[RHS, var, prF]`
- `RHS`: List of RHS expressions for ODEs
- `var`: List of variable symbols
- `prF`: Optional formatting function (default: Identity)

**Output**: `{RN, rts, spe, alp, bet, gam, rnRed}`
1. `RN`: Reaction network (list of rules)
2. `rts`: Reaction rates (list of expressions)
3. `spe`: Species names (list of strings) - **MUST equal ToString[var]**
4. `alp`: Alpha matrix (reactant stoichiometry)
5. `bet`: Beta matrix (product stoichiometry)
6. `gam`: Gamma matrix (net stoichiometry = bet - alp)
7. `rnRed`: Reduced reactions with catalysts displayed

**Side effect**: Prints reduced reactions with rates via `prF`

## Matrix Formats (alp, bet, gam)

- **Dimensions**: n × m (n species, m reactions)
- **alp[i,j]**: Stoichiometric coefficient of species i as reactant in reaction j
- **bet[i,j]**: Stoichiometric coefficient of species i as product in reaction j
- **gam[i,j]**: Net change = bet[i,j] - alp[i,j]

**Example**: For `"m0" + 2*"m1" -> "m2"` as reaction 1 with species order `{"m0", "m1", "m2"}`:
```
alp[[All, 1]] = {1, 2, 0}  (* m0: 1, m1: 2, m2: 0 *)
bet[[All, 1]] = {0, 0, 1}  (* m0: 0, m1: 0, m2: 1 *)
gam[[All, 1]] = {-1, -2, 1} (* net change *)
```

## minSiph Function

**Usage**: `minSiph[var, RN]`
- `var`: Variable list (symbols) - e.g., `{m0, m1, m2, mc, M1, M2}`
- `RN`: Reaction network (list of rules)

**Output**: `{minimal, nonMinimal}`
- `minimal`: List of minimal siphons (lists of species strings)
- `nonMinimal`: List of non-minimal siphons

**How minSiph Works**:

1. **Initialize species list**: Converts vars to strings preserving case using `toStr` helper
   ```mathematica
   toStr[s_String] := s;
   toStr[s_Symbol] := SymbolName[s];
   species = toStr /@ vars;
   ```
   **CRITICAL**: Uses `SymbolName` NOT `ToString` to avoid Format definition issues

2. **Parse reactions**: Uses `asoRea[RN]` to extract substrates and products
   ```mathematica
   reactionsAsso = asoRea[reactions];
   subs = Lookup[reactionsAsso[[j]], "Substrates", {}];
   prods = Lookup[reactionsAsso[[j]], "Products", {}];
   ```

3. **Build binary stoichiometric matrices** (n species × m reactions):
   - `alpha[i,j] = 1` if species i appears in reactants of reaction j
   - `beta[i,j] = 1` if species i appears in products of reaction j
   - Simple presence check - coefficients are irrelevant for siphon detection

4. **Find all siphons**: A subset W ⊆ species is a siphon if:
   - For all reactions where W produces something (β_W > 0)
   - W must also consume something (α_W > 0)
   - Checked via `isSiph[W, species, alpha, beta]`
   - **Species with only birth reactions (0 -> "m0") are correctly excluded**

5. **Extract minimal siphons**: Keep only siphons that are not supersets of other siphons

**CRITICAL FIXES**:
- ODE2RN: Changed from `ToLowerCase` to `ToString /@ var` to preserve case
- minSiph: Changed from `ToString` to `SymbolName` (via toStr) to avoid Format definitions
- minSiph: Now uses `asoRea` for parsing (matches Siphons.wl implementation)
- Result: Species like M1 and m1 are now correctly distinguished

## Common Issues and Debugging

### Issue 1: minSiph returns {}
**Cause**: Species mismatch between RN and spe
**Solution**: Ensure `spe == ToString[var]` (no ToLowerCase!)

### Issue 2: Species with capital letters disappear
**Cause**: ToLowerCase converts M1 → "m1", conflicting with m1
**Solution**: Never use ToLowerCase on species names

### Issue 3: bdAn hangs in infinite loop
**Cause**: bdAn was calling `minSiph[spe, RN]` where spe comes from extMat with ToLowerCase → duplicates → Subsets generates massive redundant candidates
**Solution**: Changed bdAn to call `minSiph[var, RN]` to use original variable names (epi.wl:537)
**Fixed**: AlizE.nb with uppercase species (M1, M2) now works correctly

## Helper Functions

### toStr - Symbol to String Converter

**Purpose**: Converts symbols to strings preserving case, avoiding Format definition issues

**Definitions**:
```mathematica
toStr[s_String] := s;
toStr[s_Symbol] := SymbolName[s];
toStr[s_] := SymbolName[s];
```

**Why SymbolName not ToString**:
- `ToString` respects Format definitions which corrupts strings
- Example: If `Format[i1] := Subscript[i,1]`, then `ToString[i1]` gives formatted string
- `SymbolName[i1]` always gives `"i1"` regardless of Format
- See CLAUDE.md for full explanation

**Usage in epi.wl**: minSiph (line 408) to convert var to species strings

### comp2Asso - Complex to Association Parser

**Purpose**: Parses reaction side expressions (lhs or rhs of RN rules) into associations

**Input**: Expression like `"m0" + 2*"m1"` or `0` or `"M1"`

**Output**: Association like `<|"m0" -> 1, "m1" -> 2|>` or `<||>` (empty)

**Usage in epi.wl** (3 locations):

1. **Line 31-32: `extSpe`** - Extract species names from reactions
   - Why needed: Gets all unique species by extracting keys from comp2Asso output

2. **Line 238-239: `defi`** - Build stoichiometric matrices for deficiency calculation
   - Why needed: Deficiency δ = Nc - l - s requires stoichiometric rank s = rank(γ)
   - Must have actual coefficients (not just 0/1) to compute matrix rank
   - Called by `extMat` at line 322

3. **Line 282-283: `extMat`** - Build stoichiometric matrices with coefficients
   - Why needed: Constructs full α, β, γ matrices with stoichiometric coefficients
   - Used for building RHS = γ · Rv (mass action kinetics)

**Note**: `minSiph` was simplified to NOT use comp2Asso - it now uses `Cases[expr, _String, Infinity]` to extract species directly since it only needs presence/absence detection

**How it works**:
- Expands expression into terms (Plus → list)
- For each term:
  - String alone → coefficient 1
  - `coeff * "species"` → extract coefficient
  - Converts species to lowercase (for consistency)
- Returns Association: species → coefficient

**Example**:
```mathematica
comp2Asso["m0" + 2*"M1"]  (* Returns <|"m0" -> 1, "m1" -> 2|> *)
comp2Asso[0]              (* Returns <||> (empty) *)
comp2Asso["M1"]           (* Returns <|"m1" -> 1|> *)
```

**CRITICAL**: comp2Asso converts to lowercase! This is OK because:
- It's only used internally within functions
- Final `spe` output from ODE2RN preserves original case

## NGM Function - Next Generation Matrix

**Purpose**: Computes the Next Generation Matrix K for epidemic models at Disease-Free Equilibrium (DFE)

**Three versions exist:**

### Version 1: Simple (Boundary.wl lines 230-258, commented out)
```mathematica
NGM[mod_, infVars_] := Module[{...},
  (* Simple F computation *)
  V1 = -Jx /. Thread[X[[infc]] -> 0];
  F1 = Jx + V1 /. Thread[X[[inf]] -> 0];
  F = posM[F1];  (* Just apply posM to F1 *)
  V = F - Jx;

  (* No singularity check - directly inverts V *)
  K = (F . Inverse[V]) /. Thread[X[[inf]] -> 0] // FullSimplify;
  {Jx, F, V, K, Jy, Jxy, Jyx, Kd}
]
```
**Issues**: No singularity check → Inverse::sing errors

### Version 2: With Fuser (Boundary.wl lines 261-291, active)
```mathematica
NGM[mod_, infVars_, Fuser_:{}] := Module[{...},
  If[Fuser === {},
    F = posM[F1];
    V = F - Jx;,
    (* User provides F - check if Inverse[V] is positive *)
    V = Fuser - Jx;
    Vinv = Inverse[V] /. Thread[X[[inf]] -> 0];  (* NO SINGULARITY CHECK *)
    allPos = AllTrue[Flatten[Vinv], onlyP];
    If[allPos,
      F = Fuser;  (* Accept user F *),
      F = posM[F1];  (* Reject, use posM *)
      V = F - Jx;
    ];
  ];
  K = (F . Inverse[V]) /. Thread[X[[inf]] -> 0] // Simplify;  (* NO CHECK *)
  {Jx, F, V, K, Jy, Jxy, Jyx, Kd}
]
```
**Issues**: No singularity check before Inverse[V] → Inverse::sing errors

### Version 3: With singularity checks (epi.wl lines 492-541, active)
```mathematica
NGM[mod_, infVars_, Fuser_:{}] := Module[{...},
  If[Fuser === {},
    F = posM[F1];
    V = F - Jx;,
    V = Fuser - Jx;
    (* CHECK V SINGULARITY BEFORE INVERTING *)
    Module[{Vdfe = V /. Thread[X[[inf]] -> 0]},
      If[MatrixQ[Vdfe] && MatrixRank[Vdfe] == Length[Vdfe],
        Vinv = Inverse[V] /. Thread[X[[inf]] -> 0];
        allPos = AllTrue[Flatten[Vinv], onlyP];
        If[allPos, F = Fuser;, F = posM[F1]; V = F - Jx;],
        (* V singular - fallback to posM *)
        F = posM[F1]; V = F - Jx;
      ]
    ];
  ];
  (* CHECK V SINGULARITY BEFORE COMPUTING K *)
  Module[{Vdfe = V /. Thread[X[[inf]] -> 0]},
    If[MatrixQ[Vdfe] && MatrixRank[Vdfe] == Length[Vdfe],
      K = (F . Inverse[V]) /. Thread[X[[inf]] -> 0] // Simplify;
      Kd = (Inverse[V] . F) /. Thread[X[[inf]] -> 0] // Simplify;,
      K = 0; Kd = 0;  (* Singular - return 0 *)
    ]
  ];
  {Jx, F, V, K, Jy, Jxy, Jyx, Kd}
]
```
**Fix**: Singularity checks prevent Inverse::sing errors, returns K=0 when V singular

**Current status**: epi.wl uses Version 3 (with singularity checks)

### asoRea - Association Reaction Parser

**Status**: **USED BY minSiph** in epi.wl

**Definition**: Line 338
**Usage**: minSiph uses it to parse reactions into substrates/products

**What it does**: Converts RN to list of associations with "Substrates" and "Products" keys

**Example**:
```mathematica
asoRea[{"m0" + "m1" -> 2*"m1"}]
(* Returns {<|"Substrates" -> {"m0", "m1"}, "Products" -> {"m1"}|>} *)
```

**Why needed**: minSiph needs to extract substrate and product species from each reaction for building alpha/beta matrices
