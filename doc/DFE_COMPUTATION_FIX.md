# Disease-Free Equilibrium (DFE) Computation - Critical Fix

## The Problem

When computing Disease-Free Equilibria using minimal siphons, the substitution fails because of a type mismatch:

**User's Code:**
```mathematica
mSi = minSiph[spe, asoRea[RN]]
cInf = Thread[Union[mSi//Flatten]->0]
eqD = Thread[(RHS/.cInf)==0]
```

**What Goes Wrong:**
1. `minSiph` returns siphons as **lists of strings** (e.g., `{{"v1", "v2", "v3"}, {}}`)
2. `cInf` creates **string-based rules**: `{"v1" -> 0, "v2" -> 0, "v3" -> 0}`
3. `RHS` contains **symbols** (`v1`, `v2`, `v3`), not strings
4. The substitution `RHS/.cInf` **fails silently** because strings don't match symbols
5. Result: `eqD` still contains `v1`, `v2`, `v3` instead of setting them to 0

**Example Output Showing the Bug:**
```mathematica
cInf = {"v1" -> 0, "v2" -> 0, "v3" -> 0}  (* String rules *)

eqD = {k1 v1 - k4 v1 x1 - k7 v1 z == 0,   (* v1 NOT substituted! *)
       k2 v2 - k5 v2 x2 - k8 v2 z == 0,   (* v2 NOT substituted! *)
       k3 v3 - k6 v3 x3 - k9 v3 z == 0,   (* v3 NOT substituted! *)
       ...}
```

## The Solution

**Convert strings to symbols using `ToExpression` before creating the replacement rules:**

```mathematica
cInf = Thread[ToExpression[Union[mSi//Flatten]] -> 0]
```

**Complete Corrected Workflow:**
```mathematica
(* Extract minimal siphons *)
mSi = minSiph[spe, asoRea[RN]]

(* CORRECT: Convert strings to symbols *)
cInf = Thread[ToExpression[Union[mSi//Flatten]] -> 0]

(* Now substitution works *)
eqD = Thread[(RHS/.cInf) == 0]

(* Solve for DFE *)
var = ToExpression[spe]
cDFE = Solve[eqD, var]
```

## Why This Happens

The EpidCRN package uses **strings** for species names throughout the computation pipeline:
- `extSpe[RN]` returns strings like `{"v1", "v2", "v3", "x1", ...}`
- `minSiph` returns siphons as lists of strings
- But `RHS` (the right-hand side of the ODE system) contains **symbols** for computation

To bridge this gap, use `ToExpression` to convert strings → symbols when doing substitutions.

## Key Rule

**Strings for storage and indexing → Symbols for mathematics**

- Use **strings** when working with: `spe`, `mSi`, species names
- Convert to **symbols** with `ToExpression` when doing: substitutions, solving, calculus

## Testing

The file `test_fixes.wl` (Test 6) demonstrates both the wrong and correct approaches.
