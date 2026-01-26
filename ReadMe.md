# EpidCRN Package - ReadMe

## Standard Notebook Startup (for RN-input notebooks)

```mathematica
ClearAll["Global`*"]
SetDirectory[NotebookDirectory[]]; SetOptions[$FrontEndSession, NotebookAutoSave -> True];
NotebookSave[];
AppendTo[$Path, FileNameJoin[{$HomeDirectory, "Dropbox", "EpidCRNmodels"}]];
<<EpidCRN`;
(* 0: Latex dictionary *)
prF[ex_] := ex /. {
  x1->Subscript[x,1], x2->Subscript[x,2], x3->Subscript[x,3],
  S1->Subscript[S,1], S2->Subscript[S,2],
  B1->Subscript[B,1], B2->Subscript[B,2],
  et -> η, de->δ, be->β, la -> λ, La->Λ, mu -> μ,
  mu1 -> Subscript[μ, 1], mu2 -> Subscript[μ, 2],
  ga1 -> Subscript[γ, 1], ga2 -> Subscript[γ, 2],
  a1 -> Subscript[α, 1], a2 -> Subscript[α, 2], al -> α,
  be1 -> Subscript[β, 1], be2 -> Subscript[β, 2], be3 -> Subscript[β, 3],
  ep1 -> Subscript[ε, 1], ep2 -> Subscript[ε, 2],
  om -> ω};

(* Then define RN and var *)
RN = {...};
var = {...};
```

**Key points:**
- Always reload package with `<<EpidCRN\`` (no caching issues)
- prF is model-specific LaTeX dictionary for pretty printing
- Modify prF rules for each model's variables

---

## CRITICAL RULES FOR AI AGENTS

1. **⚠️ NEVER MODIFY ANY PACKAGE WITHOUT EXPLAINING WHAT YOU DO (EXCEPT ClaudeFuckups.wl) ⚠️**
   - All subpackages (epi.wl, Core.wl, Equilibria.wl, etc.) contain trusted functions
   - Before modifying any package file, you must explain what you plan to change and why
   - ClaudeFuckups.wl is the ONLY exception - you can freely modify it for experimental versions
   - If you need to create a new version of a function, create it in ClaudeFuckups.wl first

2. **⚠️ ALWAYS TEST ClaudeFuckups.wl COMPILES BEFORE PRESENTING ⚠️**
   - Before offering any code in ClaudeFuckups.wl to the user, verify the package loads without errors
   - Test by loading: `<<EpidCRN\`` and checking for compilation errors
   - Only modify ClaudeFuckups.wl (not other subpackages) until user approves promotion

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

## Package Structure

**Subpackages:**
- **epi.wl** - Foundational functions: minSiph, ODE2RN, NGM, etc. **NEVER modify without permission.**
- **Siphons.wl** - Siphon analysis functions
- **CRNT.wl** - Chemical Reaction Network Theory
- **Utils.wl** - Utility functions
- **Boundary.wl** - Boundary equilibria analysis: Hur3M, Hur4M, Hur5M
- **Equilibria.wl** - Equilibria classification: clsEq, clsEqb, sipLat
- **ClaudeFuckups.wl** - Experimental/untested versions of functions. Use for testing new versions before promoting to main subpackages.
- **Core.wl** - (Legacy, separate from epi.wl)

**Known Issue - minSiph Usage Display:**
- When clicking on `minSiph` in the ?EpidCRN`* list, the full code definition is displayed instead of just the usage statement
- This occurs even though the usage statement is correctly defined in EpidCRN.wl (line 9)
- All other functions display only their usage statement as expected
- This is a package structure issue that needs investigation
- **Workaround:** Use `?minSiph` directly to see the usage statement
- **TODO:** Fix package structure to ensure minSiph displays usage like other functions

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

---

\subsection{Equilibria Classification Module (January 2025)}

\subsubsection{Motivation}
After \texttt{bdAn} computes minimal siphons, DFE, reproduction functions and $R_0$, we need systematic classification of equilibria into: rational, non-rational with explicit RUR, and time-outs.

\subsubsection{Algorithm Overview}
\textbf{Input:} Output from \texttt{bdAn}

\begin{verbatim}
{RHS, var, par, cp, mSi, Jx, Jy, cDFE, E0, F, V, K, R0A, infVars}
\end{verbatim}

\textbf{Steps:}
\begin{enumerate}
\item Construct lattice of minimal siphons, order by decreasing length
\item Remove cDFE rule (DFE already found by \texttt{bdAn}, pruned from future results)
\item Process rules by length level (longest first)
\item For each siphon: attempt \texttt{Solve} with 300s timeout and 1GB memory limit
\item Classify solutions: rational vs non-rational (detecting \texttt{Sqrt}, \texttt{Root}, \texttt{AlgebraicNumber})
\item Filter DFE from stored solutions, but include in counts
\item Count RUR as equations (not solutions): even number of solutions counted as half (likely $\pm$ pairs)
\end{enumerate}

\textbf{Output:} Per-siphon list: \texttt{\{siphon, nRat, nRur, \{ratSols, rurSols\}\}} or \texttt{"limit"} if Sp/Time limit exceeded

\subsubsection{Implementation: Equilibria.wl}

\textbf{Location:} \texttt{C:\textbackslash Users\textbackslash flori\textbackslash Dropbox\textbackslash EpidCRNmodels\textbackslash Equilibria.wl}

\textbf{Function Hierarchy:}

\textbf{Level 1:} \texttt{sipLat[mSi, vars]} -- Builds siphon lattice
\begin{itemize}
\item Input: minimal siphons from \texttt{bdAn}, variable list
\item Output: list of siphon rules ordered by decreasing length
\item Removes cDFE rule (union of all minimal siphons)
\end{itemize}

\textbf{Level 2:} \texttt{slvLvl[eqs, vars, rules, solved, dfe]} -- Solves one length level
\begin{itemize}
\item For each rule: \texttt{Solve[eqs /. rule, vars]} with 300s timeout
\item Filters: DFE, duplicates from solved list
\item Separates rational (polynomial ratios) from non-rational (contains \texttt{Sqrt}, \texttt{Root})
\item Returns: \texttt{\{ratSols, nonRatSols, timeOuts\}}
\end{itemize}

\textbf{Level 3:} \texttt{slvAll[eqs, vars, sipRules, dfe]} -- Iterates all levels
\begin{itemize}
\item Processes longest siphons first
\item Maintains solved accumulator to avoid duplicates
\item Returns: \texttt{\{ratSols, rurSols, timeOuts\}}
\end{itemize}

\textbf{Level 4:} \texttt{clsEqb[RHS, var, mSi, E0, prF]} -- Main interface
\begin{itemize}
\item Creates equilibrium equations: \texttt{Thread[RHS == 0]}
\item Solves for each siphon separately
\item \texttt{prF}: Optional print flag (default \texttt{False}). If \texttt{True}, prints results per siphon with rational solutions tracked to avoid reprinting
\item Returns per-siphon classification: list of \texttt{\{siphon, nRat, nRur, \{ratSols, rurSols\}\}}
\item \textbf{Counting logic:}
  \begin{itemize}
  \item nRat: Total rational solutions (including DFE if present)
  \item nRur: Number of distinct RUR equations (not solutions). Heuristic: if even number of solutions, likely $\pm$ pairs, count as half
  \item Stored solutions: Exclude DFE from storage
  \end{itemize}
\item \textbf{RUR handling:} If solutions are extremely long ($>$5000 chars), treats as RUR and solves reduced system keeping $x_1$ (substrate) as parameter
\end{itemize}

\subsubsection{Usage Example}

\begin{verbatim}
<<EpidCRN`;

(* Define model *)
RHS = {La - be1*i1*s - be2*i2*s - mu*s,
       be1*i1*s - (ga1 + mu)*i1,
       be2*i2*s - (ga2 + mu)*i2};
var = {s, i1, i2};

(* Convert to reaction network *)
{RN, rts, spe, alp, bet, gam, rnRed} = ODE2RN[RHS, var];

(* Boundary analysis *)
{RHS, var, par, cp, mSi, Jx, Jy, cDFE, E0, F, V, K, R0A, infVars} =
  bdAn[RN, rts, var];

(* Classify equilibria per siphon *)
perSiphon = clsEqb[RHS, var, mSi, E0];

(* With printing enabled *)
perSiphon = clsEqb[RHS, var, mSi, E0, True];

(* Each entry: {siphon, nRat, nRur, {ratSols, rurSols}} *)
(* If prF=False, manually print results *)
Do[
  {sip, nRat, nRur, sols} = perSiphon[[i]];
  Print["Siphon ", i, ": ", sip, " rat=", nRat, " rur=", nRur];
, {i, Length[perSiphon]}];
\end{verbatim}

\subsubsection{Test Results}

\textbf{Two-strain SIS model:}
\begin{itemize}
\item Variables: \texttt{\{s, i1, i2\}}
\item Minimal siphons: \texttt{\{\{i1\}, \{i2\}\}}
\item Per-siphon results:
  \begin{itemize}
  \item Siphon 1 \texttt{\{i2\}}: \texttt{rat=2, rur=0} (1 DFE + 1 endemic with strain 1 only)
  \item Siphon 2 \texttt{\{i1\}}: \texttt{rat=2, rur=0} (1 DFE + 1 endemic with strain 2 only)
  \end{itemize}
\end{itemize}

\textbf{OSN model (Online Social Network with botnets):}
\begin{itemize}
\item Variables: \texttt{\{x1, x2, x3, B1, S1, B2, S2, R\}} (8 variables)
\item Minimal siphons: \texttt{\{\{x2\}, \{B1, S1\}, \{B2, S2\}\}}
\item Siphon lattice: 6 siphons (excluding cDFE)
\item Per-siphon results:
  \begin{itemize}
  \item Siphon 1 \texttt{\{B1,S1,B2,S2\}}: \texttt{rat=3, rur=0} (1 DFE + 2 endemic, stores 2)
  \item Siphons 2,3 (len=3): \texttt{rat=1, rur=0} (DFE only)
  \item Siphons 4,5 \texttt{\{B2,S2\}}, \texttt{\{B1,S1\}}: \texttt{rat=3, rur=1} (1 DFE + 2 endemic + 1 RUR equation with 2 roots)
  \item Siphon 6 \texttt{\{x2\}}: \texttt{rat=1, rur=0} (DFE only)
  \end{itemize}
\end{itemize}

\textbf{Key Implementation Details:}

\textbf{1. Detection of non-rational solutions:}
\begin{verbatim}
isRational[s_] := FreeQ[s, Root | RootOf | AlgebraicNumber] &&
  FreeQ[s, _Power?((!IntegerQ[#[[2]]] && Head[#[[2]]] =!= Symbol)&)]
\end{verbatim}
This correctly identifies \texttt{Sqrt[x]} as \texttt{Power[x, 1/2]} (non-integer exponent) and classifies it as non-rational.

\textbf{2. RUR equation counting:}
RUR count represents number of distinct algebraic equations, not individual roots. When \texttt{Solve} returns multiple solutions differing by $\pm$ signs (e.g., $x = a \pm \sqrt{b}$), these are roots of a single equation. Heuristic: if even number of RUR solutions, count as half (likely $\pm$ pairs). For epidemic models, only positive roots are biologically valid.

\textbf{3. RUR safeguard - keeping substrate as parameter:}
CRITICAL: If a solution passes \texttt{isRational} but has extremely long expressions ($>$5000 characters), it is likely irrational. For RUR cases, the proper approach is:
\begin{itemize}
\item Keep $x_1$ (substrate, first variable in \texttt{var}) as a free parameter
\item Delete the first equation (equation for $x_1$)
\item Solve the reduced system for remaining variables
\item Result: RUR solutions expressed in terms of $x_1$
\end{itemize}
This prevents spurious "rational" solutions with extremely long denominators that are actually irrational. The code automatically tries this reduced approach when RUR is detected and uses it if it produces cleaner output.

\subsubsection{Bug Fix: Complementarity Filtering (January 2025)}

\textbf{Problem:} Siphon 1 was not printing its boundary equilibria, and Siphons 4-5 were reprinting solutions already shown in Siphon 1.

\textbf{Root Cause:} The \texttt{hasOtherSiphonZero} helper function builds \texttt{antiSiph} -- the union of OTHER minimal siphons. A solution on one siphon (e.g., $\{B_2, S_2\} = 0$) should have variables from other minimal siphons non-zero (e.g., $\{B_1, S_1, x_2\} \neq 0$). The bug occurred because:
\begin{enumerate}
\item Initially used \texttt{Complement[Rest[vars], siphonVars]} which included ALL non-siphon variables (like R which isn't in any minimal siphon)
\item This incorrectly filtered out valid boundary equilibria with $R=0$
\item Later fix used \texttt{SubsetQ[siphonVars, mSi\_element]} to find other siphons, but \texttt{SubsetQ} failed because mSi elements might be strings while siphonVars contains symbols
\end{enumerate}

\textbf{Solution:} Convert mSi elements to symbols before SubsetQ test:
\begin{verbatim}
hasOtherSiphonZero[s_] := Module[{otherSiphons, antiSiph, mSiSymbols},
  (* Convert mSi elements to symbols first *)
  mSiSymbols = Map[
    Map[If[StringQ[#], ToExpression[#], #]&, #]&,
    mSi
  ];
  (* Find minimal siphons NOT contained in current siphon *)
  otherSiphons = Select[mSiSymbols, !SubsetQ[siphonVars, #]&];
  (* Union of other minimal siphons *)
  antiSiph = Flatten[otherSiphons, 1];
  antiSiph = DeleteDuplicates[antiSiph];
  (* Check if any are zero *)
  !nonZat[s, antiSiph]
];
\end{verbatim}

\textbf{Example (OSN model):}
\begin{itemize}
\item \texttt{mSi = \{\{x2\}, \{B1, S1\}, \{B2, S2\}\}}
\item For Siphon 1 with \texttt{siphonVars = \{B1, S1, B2, S2\}}:
  \begin{itemize}
  \item Before fix: \texttt{antiSiph = \{x2, x3, R\}} (WRONG - includes R which isn't in any minimal siphon)
  \item After fix: \texttt{antiSiph = \{x2\}} (CORRECT - only the other minimal siphon)
  \item Solutions with $x_2 > 0$ now correctly pass filter and print
  \end{itemize}
\end{itemize}

\textbf{Related fixes:}
\begin{itemize}
\item \texttt{prF} parameter added to clsEqb signature: formatting function applied with \texttt{//prF}
\item Return structure changed from \texttt{\{siphonVars, nRat, nRur, sols\}} to \texttt{\{nRat, nRur, sols\}}
\item \texttt{printedRat} list tracks already-printed solutions to prevent duplication
\item RUR filtering uses \texttt{Join[solRed[[i]], rule]} to complete solution before checking complementarity
\end{itemize}

\subsubsection{Integration}
\begin{enumerate}
\item File added: \texttt{Equilibria.wl} in root directory
\item Loader updated: \texttt{Get[FileNameJoin[\{root, "Equilibria.wl"\}]]} in \texttt{EpidCRN.wl} (line 219)
\item Usage statements added to \texttt{EpidCRN.wl} (lines 69-72)
\item Test files: \texttt{testEquilibria2.wl}, \texttt{testOSNfull.wl}, \texttt{testSqrt.wl}
\end{enumerate}

\subsection{clsEq: Endemic Equilibrium Solver with Timeout (2026-01-14)}

\subsubsection{Purpose}
New function \texttt{clsEq[RHS, var, mSi, E0, prF, tim]} specifically tackles endemic equilibrium (no siphon variables zero) with configurable timeout and fallback to RUR.

\subsubsection{Algorithm}
\begin{enumerate}
\item First run \texttt{clsEqb[RHS, var, mSi, E0, prF]} and print all siphon results
\item Then attempt unrestricted \texttt{Solve[eqs, vars]} with timeout \texttt{tim} for endemic equilibrium
\item If timeout: Attempt RUR approach (drop first equation, solve reduced system, compute polynomial in x1)
\item If RUR also times out: Print "TIME-OUT" and return empty list
\item Filter: Remove DFE, solutions with any siphon zero, all-zero solutions
\end{enumerate}

\subsubsection{Output}
\begin{itemize}
\item Success (Solve): Prints \texttt{nSol=N} followed by \texttt{E1=Short[...]} through \texttt{EN=Short[...]}
\item Success (RUR): Prints \texttt{nRUR=N} followed by \texttt{ratS=Short[...]} and \texttt{RUR=Short[...]} for each polynomial
\item Failure: Prints "TIME-OUT"
\item Returns: List of solution rules (empty if timeout or no solutions)
\end{itemize}

\subsubsection{Usage}
\texttt{endemic = clsEq[RHS, var, mSi, E0, prF, 60]}
Test file: \texttt{testClsEq.wl} with 60-second timeout.

## Reserved Variable Names

**IMPORTANT**: The following variable names are reserved by the package and should NOT be used in model RHS definitions:

- **E0**: Equilibrium values (numeric/symbolic values of variables at equilibrium)
- **cE0**: Equilibrium rule (list of rules var -> value at DFE, returned by minSiph)
- **cDFE**: Disease-Free Equilibrium rule (siphon variables -> 0)
- **V, F, K**: NGM matrices (Transition, New infections, Next Generation Matrix)
- **Jx, Jy, Jxy, Jyx, Kd**: Jacobian matrices and blocks
- **RN, rts, var, par**: Reaction network, rates, variables, parameters
- **alp, bet, gam, spe**: Stoichiometric matrices (alpha, beta, gamma, species)
- **mSi, inf, infc, mod, cp, sol**: Minimal siphons, infection vars, complements, model, constraints, solutions
- **I, E**: Reserved to avoid conflict with built-in Mathematica symbols (imaginary unit, exponential)

Use the \texttt{chk[RHS]} function immediately after defining RHS to detect conflicts:
```mathematica
RHS = {...};
chk[RHS];  (* Warns if any variables conflict with reserved names *)
```

**Function Signatures**: All package functions MUST use reserved names in their input/output signatures (e.g., return `V` not `VV`). Internal implementations may use alternative names (e.g., `VV` internally) but must map back to reserved names before returning.

---

