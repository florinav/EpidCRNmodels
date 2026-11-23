#!/usr/bin/env python3
"""
yaml2wl.py — convert a flepiMoP YAML config to an EpidCRN-compatible .wl file.

Usage (Windows Command Prompt or PowerShell):
    python yaml2wl.py  path\to\config.yml
    python yaml2wl.py  path\to\config.yml  output.wl

If the second argument is omitted the output file is placed next to the YAML
with the same stem and extension .wl.

Requires:  Python >= 3.9,  PyYAML  (pip install pyyaml)
"""

import sys
import re
from pathlib import Path

# ── dependency check ──────────────────────────────────────────────────────────
try:
    import yaml
except ImportError:
    sys.exit(
        "\nPyYAML is not installed.\n"
        "Open Command Prompt and run:\n\n"
        "    pip install pyyaml\n\n"
        "Then try again.\n"
    )

# ── helpers ───────────────────────────────────────────────────────────────────

def safe(name: str) -> str:
    """
    Lowercase a compartment name so it is safe in Mathematica.
    SEIR → seir;  S → s;  E → e  (avoids E = Euler, I = imaginary unit).
    """
    return name.lower()

def to_mma(expr: str) -> str:
    """Convert a Python-style expression string to Mathematica syntax."""
    expr = str(expr).strip()
    expr = re.sub(r'\s*\*\*\s*', '^', expr)   # Python ** → Mathematica ^
    expr = re.sub(r'\s+', '', expr)            # remove all spaces
    return expr

def monomial(rate: str, prop_to: list, prop_exp: list) -> str:
    """
    Build  rate * sp1^e1 * sp2^e2 * ...
    prop_to  : list; each element is either a list of species names
               or a bare string (unstratified case).
    prop_exp : list of exponent strings, parallel to prop_to.
    For stratified models only the first species in each group is used;
    full stratification support is a planned extension.
    """
    factors = [rate]
    for group, exp in zip(prop_to, prop_exp):
        sp = safe(group[0] if isinstance(group, list) else str(group))
        e  = to_mma(str(exp))
        factors.append(sp if e == '1' else f'{sp}^{e}')
    return '*'.join(factors)

# ── Greek/display map for prF ─────────────────────────────────────────────────
GREEK = {
    'sigma': r'\[Sigma]',
    'gamma': r'\[Gamma]',
    'ro':    r'Subscript[R,0]',
    'beta':  r'\[Beta]',
    'mu':    r'\[Mu]',
    'delta': r'\[Delta]',
    'alpha': r'\[Alpha]',
    'eta':   r'\[Eta]',
    'lam':   r'\[Lambda]',
    'nu':    r'\[Nu]',
    'kappa': r'\[Kappa]',
    'rho':   r'\[Rho]',
}

# ── main conversion ───────────────────────────────────────────────────────────

def yaml2wl(yaml_path: str, wl_path: str = None):
    src = Path(yaml_path)
    if not src.exists():
        sys.exit(f"File not found: {src}")

    dst = Path(wl_path) if wl_path else src.with_suffix('.wl')

    with open(src, encoding='utf-8') as fh:
        cfg = yaml.safe_load(fh)

    # ── compartments ──────────────────────────────────────────────────────────
    comp   = cfg.get('compartments', {})
    stages = comp.get('infection_stage', [])
    if not stages:
        sys.exit("No compartments.infection_stage found in the YAML.")

    var_orig = list(stages)               # e.g. ["S","E","I","R"]
    var_mma  = [safe(s) for s in stages]  # e.g. ["s","e","i","r"]

    # ── parameters ────────────────────────────────────────────────────────────
    seir   = cfg.get('seir', {})
    params = seir.get('parameters', {})

    # ── transitions → ODE accumulator ────────────────────────────────────────
    transitions = seir.get('transitions', [])
    if not transitions:
        sys.exit("No seir.transitions found in the YAML.")

    ode = {v: [] for v in var_mma}

    for tr in transitions:
        # rate: YAML stores it as a one-element list
        rate_raw = tr['rate']
        rate = to_mma(rate_raw[0] if isinstance(rate_raw, list) else rate_raw)

        # proportional_to: normalise to list-of-lists
        prop_raw = tr.get('proportional_to', [])
        prop_to  = [p if isinstance(p, list) else [p] for p in prop_raw]
        prop_exp = list(tr.get('proportion_exponent', ['1'] * len(prop_to)))

        flux = monomial(rate, prop_to, prop_exp)

        # source loses flux; destination gains flux
        for field, sign in [('source', '-'), ('destination', '+')]:
            names = tr[field]
            if not isinstance(names, list):
                names = [names]
            for name in names:
                key = safe(name)
                if key in ode:
                    ode[key].append(sign + flux)

    # ── prF substitution rules ────────────────────────────────────────────────
    prf_rules = []
    for pname in params:
        g = GREEK.get(pname.lower())
        if g:
            prf_rules.append(f'  {pname} -> {g}')
    for v_mma, v_orig in zip(var_mma, var_orig):
        prf_rules.append(f'  {v_mma} -> Subscript[{v_orig}, ""]')

    # ── format ODE entries ────────────────────────────────────────────────────
    def fmt(v):
        terms = ode[v]
        if not terms:
            return '  0'
        return '  ' + ''.join(terms).lstrip('+')

    # ── assemble .wl text ─────────────────────────────────────────────────────
    model = cfg.get('name', src.stem)
    L = []   # lines

    L += [
        '(* ::Package:: *)',
        '',
        'ClearAll["Global`*"];',
        'SetDirectory[NotebookDirectory[]];',
        'SetOptions[$FrontEndSession, NotebookAutoSave -> True];',
        'NotebookSave[];',
        'AppendTo[$Path, FileNameJoin[{$HomeDirectory, "Dropbox", "EpidCRNmodels"}]];',
        '<<EpidCRN`;',
        '',
        f'(* Model : {model} *)',
        f'(* Source: {src.name} — generated by yaml2wl.py *)',
        '',
        '(* Parameters — kept symbolic; values from YAML are noted below *)',
    ]
    for pname, pval in params.items():
        v = pval.get('value', '?') if isinstance(pval, dict) else pval
        L.append(f'(* {pname} = {v} *)')

    L += [
        '',
        f'var = {{{", ".join(var_mma)}}};',
        '',
        'RHS = {',
        ',\n'.join(fmt(v) for v in var_mma),
        '};',
        '',
    ]

    if prf_rules:
        L += [
            'prF[ex_] := ex /. {',
            ',\n'.join(prf_rules),
            '};',
        ]
    else:
        L.append('prF[ex_] := ex;  (* add display rules here *)')

    L += [
        '',
        '(* ── CRN decomposition ───────────────────────────────────────── *)',
        '{RN, rts, spe, alp, bet, gam, rnR} = ODE2RNs[RHS, var, prF];',
        'Print["verif RHS=", RHS - gam . rts // FullSimplify];',
        '',
        '(* ── Minimal siphons ─────────────────────────────────────────── *)',
        'allSi = minSiph[var, RN];',
        'Print["minSiphs=", allSi // prF];',
        '',
        '(* ── NGM and R0 ──────────────────────────────────────────────── *)',
        '{RHS2, par, cp, mSi, Jx, Jy, Jxy, fx, cDFE, cE0, FF, VV, K, R0A, infV} =',
        '    NGMRN[RN, rts, var];',
        'Print["E0=", cE0 // prF, "  infVars=", infV, "  R0=", R0A // prF];',
        'Print["K=", K // prF // MatrixForm];',
        '',
        '(* ── Child-selection and oscillation exclusion ───────────────── *)',
        '{n, m} = {Length[var], Length[RN]};',
        'RM = symRM[alp, n, m];',
        'co = Rest @ Reverse @ CoefficientList[',
        '         Expand[CharacteristicPolynomial[gam . RM, lam]], lam];',
        'nCS = negCS[co, gam, n];',
        'upf = Select[nCS, upfQ];',
        '{gH, edges} = hasseCS[nCS];',
        'Print["negCS: n=", n, " m=", m, "  #nCS=", Length[nCS],',
        '      "  #UPF=", Length[upf], "  #Hasse-edges=", Length[edges]];',
        'riCores  = oscRI[nCS, gH, spe, rts];',
        'riiCores = oscRII[nCS, gH, spe, rts];',
        'Print["oscRI #cores=", Length[riCores],',
        '      "  oscRII #cores=", Length[riiCores]];',
    ]

    dst.write_text('\n'.join(L), encoding='utf-8')
    print(f"Written: {dst}")

# ── entry point ───────────────────────────────────────────────────────────────

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(0)
    yaml2wl(*sys.argv[1:3])
