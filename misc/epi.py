"""
epi.py - Minimal epidemic modeling package
Python translation of epi.wl
Contains: ODE2RN, bdAn, and all dependencies
"""

import numpy as np
import sympy as sp
from sympy import symbols, Symbol, expand, simplify, solve, Poly, fraction, Together
from sympy.matrices import Matrix
from itertools import combinations
from typing import List, Dict, Tuple, Union, Any
import re

# =========================================================================
# FROM UTILS.WL
# =========================================================================

def posM(expr):
    """Keep syntactically positive terms"""
    if isinstance(expr, list):
        # Matrix case
        if isinstance(expr[0], list):
            result = []
            for row in expr:
                new_row = []
                for e in row:
                    expanded = sp.expand(e)
                    if expanded.is_Add:
                        terms = [t for t in expanded.as_ordered_terms()
                                if not (t.is_number and t < 0) and
                                not (t.is_Mul and any(a.is_number and a < 0 for a in t.as_ordered_factors()))]
                        new_row.append(sum(terms) if terms else 0)
                    elif expanded.is_number and expanded < 0:
                        new_row.append(0)
                    else:
                        new_row.append(expanded)
                result.append(new_row)
            return result
        else:
            # Vector case
            result = []
            for e in expr:
                expanded = sp.expand(e)
                if expanded.is_Add:
                    terms = [t for t in expanded.as_ordered_terms()
                            if not (t.is_number and t < 0) and
                            not (t.is_Mul and any(a.is_number and a < 0 for a in t.as_ordered_factors()))]
                    result.append(sum(terms) if terms else 0)
                elif expanded.is_number and expanded < 0:
                    result.append(0)
                else:
                    result.append(expanded)
            return result
    else:
        # Single expression
        expanded = sp.expand(expr)
        if expanded.is_Add:
            terms = [t for t in expanded.as_ordered_terms()
                    if not (t.is_number and t < 0) and
                    not (t.is_Mul and any(a.is_number and a < 0 for a in t.as_ordered_factors()))]
            return sum(terms) if terms else 0
        elif expanded.is_number and expanded < 0:
            return 0
        else:
            return expanded

def isNNe(m):
    """Check if expression is non-negative"""
    # Simplified version - just checks if coefficients are non-negative
    return True  # Placeholder - full implementation needs coefficient checking

# =========================================================================
# FROM CORE.WL
# =========================================================================

def Par(RHS, X):
    """Extract parameters from dynamics"""
    all_vars = set()
    for expr in RHS:
        all_vars.update(expr.free_symbols)
    X_set = set(X)
    return sorted(list(all_vars - X_set), key=str)

def extSpe(reactions):
    """Extract species names from reactions"""
    all_species = set()

    for rxn in reactions:
        if isinstance(rxn, tuple):
            left, right = rxn
            # Extract from left side
            if isinstance(left, str) and left not in ['0', '']:
                all_species.add(left.lower())
            elif hasattr(left, 'free_symbols'):
                for sym in left.free_symbols:
                    if isinstance(sym, Symbol):
                        all_species.add(str(sym).lower())

            # Extract from right side
            if isinstance(right, str) and right not in ['0', '']:
                all_species.add(right.lower())
            elif hasattr(right, 'free_symbols'):
                for sym in right.free_symbols:
                    if isinstance(sym, Symbol):
                        all_species.add(str(sym).lower())

    return sorted(list(all_species))

def allT(expr):
    """Extract all terms from expanded polynomial"""
    expanded = sp.expand(expr)
    if expanded.is_Add:
        return list(expanded.as_ordered_terms())
    else:
        return [expanded]

def isN(term):
    """Detect if a term has negative sign"""
    if term.is_number and term < 0:
        return True
    if term.is_Mul:
        factors = term.as_ordered_factors()
        if -1 in factors:
            return True
        if any(f.is_number and f < 0 for f in factors):
            return True
    return False

def comp2Asso(expr):
    """Convert compound expression to association (dict)"""
    if expr == 0 or expr is None or expr == '0':
        return {}

    terms = allT(expr) if hasattr(expr, 'is_Add') else [expr]
    result = {}

    for term in terms:
        if term.is_Mul and len(term.as_ordered_factors()) == 2:
            factors = term.as_ordered_factors()
            # Find which is coefficient and which is species
            if factors[0].is_number:
                coeff = factors[0]
                species = str(factors[1]).lower()
            else:
                coeff = factors[1] if factors[1].is_number else 1
                species = str(factors[0]).lower()
            result[species] = result.get(species, 0) + coeff
        elif isinstance(term, Symbol):
            species = str(term).lower()
            result[species] = result.get(species, 0) + 1
        elif term.is_number:
            continue
        else:
            species = str(term).lower()
            result[species] = result.get(species, 0) + 1

    return result

def asoRea(RN):
    """Convert RN to association (dict) format"""
    result = []
    for rxn in RN:
        if isinstance(rxn, tuple):
            substrates = comp2Asso(rxn[0])
            products = comp2Asso(rxn[1])
        else:
            substrates = comp2Asso(rxn)
            products = {}

        result.append({
            "Substrates": substrates,
            "Products": products
        })

    return result

def extMat(reactions, speOrder=None):
    """Extract stoichiometric matrices"""
    # Extract species
    if speOrder:
        spe = [str(s).lower() for s in speOrder]
    else:
        spe = extSpe(reactions)

    numSpecies = len(spe)
    numReactions = len(reactions)

    # Initialize matrices
    al = [[0 for _ in range(numReactions)] for _ in range(numSpecies)]
    be = [[0 for _ in range(numReactions)] for _ in range(numSpecies)]

    # Build stoichiometric matrices
    for j in range(numReactions):
        rxn = reactions[j]
        if isinstance(rxn, tuple):
            reactants = comp2Asso(rxn[0])
            products = comp2Asso(rxn[1])
        else:
            reactants = comp2Asso(rxn)
            products = {}

        for i in range(numSpecies):
            if spe[i] in reactants:
                al[i][j] = reactants[spe[i]]
            if spe[i] in products:
                be[i][j] = products[spe[i]]

    # Calculate gamma
    gamma = [[be[i][j] - al[i][j] for j in range(numReactions)] for i in range(numSpecies)]

    # Convert species to variables
    var = [Symbol(s) for s in spe]

    # Calculate reaction rates (mass action kinetics)
    rv = []
    for j in range(numReactions):
        rate = 1
        for i in range(numSpecies):
            if al[i][j] != 0:
                rate *= var[i] ** al[i][j]
        rv.append(rate)

    # Rate constants
    tk = [Symbol(f'k{j+1}') for j in range(numReactions)]
    Rv = [tk[j] * rv[j] for j in range(numReactions)]

    # Compute RHS
    RHS = []
    for i in range(numSpecies):
        rhs_i = sum(gamma[i][j] * Rv[j] for j in range(numReactions))
        RHS.append(sp.expand(rhs_i))

    # Deficiency (simplified)
    defResult = "deficiency_info"

    return (spe, al, be, gamma, Rv, RHS, defResult)

def ODE2RN(RHS_list, var_list):
    """Convert ODE system to reaction network"""
    n = len(var_list)
    spe = [str(v).lower() for v in var_list]

    # Extract terms
    allTermsByEq = [allT(RHS_list[i]) for i in range(n)]
    allTermsFlat = [term for eq_terms in allTermsByEq for term in eq_terms]

    # Separate sources and sinks
    sources = [t for t in allTermsFlat if not isN(t)]
    sinks = [t for t in allTermsFlat if isN(t)]
    sources = list(dict.fromkeys(sources))  # Remove duplicates
    sinks = list(dict.fromkeys(sinks))

    # Build rtsRaw
    rtsRaw = sources + [-s for s in sinks]
    rtsRaw = list(dict.fromkeys([r for r in rtsRaw if r != 0]))

    # Build gamma matrix
    gamRaw = []
    for i in range(n):
        eqTerms = allT(RHS_list[i])
        row = []
        for j in range(len(rtsRaw)):
            coeff = 0
            for term in eqTerms:
                termAbs = -term if isN(term) else term
                if termAbs == rtsRaw[j]:
                    coeff = -1 if isN(term) else 1
                    break
            row.append(coeff)
        gamRaw.append(row)

    # Build alpha matrix
    alpRaw = []
    for i in range(n):
        row = []
        for j in range(len(rtsRaw)):
            rate = rtsRaw[j]
            if var_list[i] not in rate.free_symbols:
                row.append(0)
            else:
                row.append(sp.degree(rate, var_list[i]))
        alpRaw.append(row)

    # Merge reactions with same source and gamma
    alpGamPairs = [tuple(alpRaw[i][j] for i in range(n)) + tuple(gamRaw[i][j] for i in range(n))
                   for j in range(len(rtsRaw))]
    uniquePairs = list(dict.fromkeys(alpGamPairs))

    rts = []
    gam = [[0 for _ in range(len(uniquePairs))] for _ in range(n)]
    alp = [[0 for _ in range(len(uniquePairs))] for _ in range(n)]

    for j, pair in enumerate(uniquePairs):
        # Find indices with this pair
        mergeIdx = [k for k, p in enumerate(alpGamPairs) if p == pair]
        # Sum rates
        rts.append(sum(rtsRaw[k] for k in mergeIdx))
        # Set alp and gam
        for i in range(n):
            alp[i][j] = pair[i]
            gam[i][j] = pair[n + i]

    # Build beta
    bet = [[gam[i][j] + alp[i][j] for j in range(len(rts))] for i in range(n)]

    # Build reaction network
    RN = []
    for j in range(len(rts)):
        # Left side
        leftTerms = []
        for i in range(n):
            if alp[i][j] == 1:
                leftTerms.append(Symbol(spe[i]))
            elif alp[i][j] > 1:
                leftTerms.append(alp[i][j] * Symbol(spe[i]))
        left = sum(leftTerms) if leftTerms else 0

        # Right side
        rightTerms = []
        for i in range(n):
            if bet[i][j] == 1:
                rightTerms.append(Symbol(spe[i]))
            elif bet[i][j] > 1:
                rightTerms.append(bet[i][j] * Symbol(spe[i]))
        right = sum(rightTerms) if rightTerms else 0

        RN.append((left, right))

    return (RN, rts, spe, alp, bet, gam)

# =========================================================================
# FROM SIPHONS.WL
# =========================================================================

def pos(v):
    """Check if vector is positive"""
    return all(x >= 0 for x in v) and sum(v) > 0

def isSiph(W, species, alpha, beta):
    """Check if set W is a siphon"""
    # Convert species names to indices
    indices = []
    for w in W:
        try:
            indices.append(species.index(w))
        except ValueError:
            return False

    if not indices:
        return False

    nReac = len(alpha[0])
    alphaW = [[alpha[i][j] for j in range(nReac)] for i in indices]
    betaW = [[beta[i][j] for j in range(nReac)] for i in indices]

    # Find producing reactions
    producingReactions = []
    for j in range(nReac):
        col = [betaW[i][j] for i in range(len(indices))]
        if pos(col):
            producingReactions.append(j)

    # If no producing reactions, vacuously true
    if not producingReactions:
        return True

    # All producing reactions must also consume
    for j in producingReactions:
        col = [alphaW[i][j] for i in range(len(indices))]
        if not pos(col):
            return False

    return True

def minSiph(vars_input, reactions):
    """Find minimal siphons"""
    # Convert to species strings
    if all(isinstance(v, str) for v in vars_input):
        species = vars_input
    else:
        species = [str(v).lower() for v in vars_input]

    reactionsAsso = asoRea(reactions)

    n = len(species)
    m = len(reactions)

    # Build alpha and beta matrices
    alpha = [[0 for _ in range(m)] for _ in range(n)]
    beta = [[0 for _ in range(m)] for _ in range(n)]

    for j in range(m):
        subs = reactionsAsso[j]["Substrates"]
        prods = reactionsAsso[j]["Products"]

        for i in range(n):
            if species[i] in subs:
                alpha[i][j] = 1
            if species[i] in prods:
                beta[i][j] = 1

    # Find all siphons
    siphons = []
    for size in range(1, n + 1):
        for subset in combinations(species, size):
            if isSiph(list(subset), species, alpha, beta):
                siphons.append(list(subset))

    # Find minimal siphons
    minimal = []
    for s in siphons:
        is_minimal = True
        for t in siphons:
            if t != s and set(t).issubset(set(s)):
                is_minimal = False
                break
        if is_minimal:
            minimal.append(s)

    nonMinimal = [s for s in siphons if s not in minimal]

    return (minimal, nonMinimal)

# =========================================================================
# FROM BOUNDARY.WL
# =========================================================================

def NGM(mod, infVars, Fuser=None):
    """Compute Next Generation Matrix"""
    RHS = mod[0]
    X = mod[1]
    par = mod[2] if len(mod) >= 3 else []

    # Convert infVars to indices
    inf = []
    for v in infVars:
        try:
            inf.append(X.index(v))
        except ValueError:
            continue

    n = len(X)
    infc = [i for i in range(n) if i not in inf]

    # Compute Jacobian blocks
    Jx = [[sp.diff(RHS[inf[i]], X[inf[j]]) for j in range(len(inf))] for i in range(len(inf))]

    if infc:
        Jy = [[sp.diff(RHS[infc[i]], X[infc[j]]) for j in range(len(infc))] for i in range(len(infc))]
        Jxy = [[sp.diff(RHS[inf[i]], X[infc[j]]) for j in range(len(infc))] for i in range(len(inf))]
        Jyx = [[sp.diff(RHS[infc[i]], X[inf[j]]) for j in range(len(inf))] for i in range(len(infc))]
    else:
        Jy = []
        Jxy = []
        Jyx = []

    # Compute V and F
    if infc:
        subs_infc = [(X[i], 0) for i in infc]
        V1 = [[-Jx[i][j] for j in range(len(inf))] for i in range(len(inf))]
        for i in range(len(inf)):
            for j in range(len(inf)):
                V1[i][j] = V1[i][j].subs(subs_infc)
    else:
        V1 = [[-Jx[i][j] for j in range(len(inf))] for i in range(len(inf))]

    subs_inf = [(X[i], 0) for i in inf]
    F1 = [[Jx[i][j] + V1[i][j] for j in range(len(inf))] for i in range(len(inf))]
    for i in range(len(inf)):
        for j in range(len(inf)):
            F1[i][j] = F1[i][j].subs(subs_inf)

    if Fuser is None:
        F = posM(F1)
        V = [[F[i][j] - Jx[i][j] for j in range(len(inf))] for i in range(len(inf))]
    else:
        V = [[Fuser[i][j] - Jx[i][j] for j in range(len(inf))] for i in range(len(inf))]
        F = Fuser

    # Compute K = F * V^-1
    V_mat = Matrix(V)
    F_mat = Matrix(F)

    K_mat = F_mat * V_mat.inv()
    for i in range(len(inf)):
        for j in range(len(inf)):
            K_mat[i, j] = sp.simplify(K_mat[i, j].subs(subs_inf))

    K = [[K_mat[i, j] for j in range(len(inf))] for i in range(len(inf))]

    Kd = list(K_mat.eigenvals().keys())

    return (Jx, F, V, K, Jy, Jxy, Jyx, Kd)

def bdAn(RN, rts, var):
    """Boundary analysis"""
    # Extract stoichiometric matrices
    spe, alp, bet, gam, Rv, RHS_dummy, def_info = extMat(RN, var)

    # Compute RHS
    RHS = []
    for i in range(len(gam)):
        rhs_i = sum(gam[i][j] * rts[j] for j in range(len(rts)))
        RHS.append(sp.expand(rhs_i))

    # Extract parameters
    par = Par(RHS, var)
    cp = [p > 0 for p in par]

    # Find minimal siphons
    mS = minSiph(spe, RN)
    mSi = mS[0]

    # infVars as strings
    infVars = sorted(list(set([s for siphon in mSi for s in siphon])))

    # Convert to symbols for substitution
    infVars_syms = [Symbol(s) for s in infVars]

    # Disease-free equilibrium
    cDFE = [(v, 0) for v in infVars_syms]
    RDFE = [rhs.subs(cDFE) for rhs in RHS]
    eq0 = [sp.Eq(r, 0) for r in RDFE]
    var0 = [v for v in var if v not in infVars_syms]

    try:
        E0_partial = solve(eq0, var0, dict=True)
        if E0_partial:
            E0 = list(E0_partial[0].items()) + cDFE
        else:
            E0 = cDFE
    except:
        E0 = cDFE

    # NGM analysis
    mod = [RHS, var, par]
    ngm = NGM(mod, infVars_syms)
    Jx = ngm[0]
    Jy = ngm[4]
    K = ngm[3]

    # Eigenvalues
    K_mat = Matrix(K)
    eigenvals = list(K_mat.eigenvals().keys())
    R0A = [ev for ev in eigenvals if ev != 0]

    return {
        'RHS': RHS,
        'var': var,
        'par': par,
        'cp': cp,
        'mSi': mSi,
        'Jx': Jx,
        'Jy': Jy,
        'cDFE': cDFE,
        'E0': E0,
        'K': K,
        'R0A': R0A,
        'infVars': infVars,
        'ngm': ngm
    }

# =========================================================================
# SUBSCRIPT CONVERTER
# =========================================================================

def subsCon(s):
    """Convert variable names to subscripted form (for display)"""
    if isinstance(s, Symbol):
        s = str(s)

    match = re.match(r'^([A-Za-z]+)([0-9]+)$', s)
    if match:
        base = match.group(1)
        digits = match.group(2)
        # Return formatted string (LaTeX or Unicode subscripts)
        subscript_map = str.maketrans('0123456789', '₀₁₂₃₄₅₆₇₈₉')
        return base + digits.translate(subscript_map)
    else:
        return s

# TEST - Commented out
"""
from sympy import symbols
s, i1, i2, i12, be1, be2, mu0, mu1, mu2, mu3, R1, R2, R3, s0, et1, et2, ga1, ga2 = symbols(
    's i1 i2 i12 be1 be2 mu0 mu1 mu2 mu3 R1 R2 R3 s0 et1 et2 ga1 ga2')

RHS1 = [-be1*i12*s - be2*i12*s - mu0*s - (i1*mu1*R1*s)/s0 - (i2*mu2*R2*s)/s0 - (i12*mu3*R3*s)/s0 + mu0*s0,
        -et1*i1*i12 - ga1*i1*i2 - i1*mu1 + be1*i12*s + (i1*mu1*R1*s)/s0,
        -ga2*i1*i2 - et2*i12*i2 - i2*mu2 + be2*i12*s + (i2*mu2*R2*s)/s0,
        et1*i1*i12 + ga1*i1*i2 + ga2*i1*i2 + et2*i12*i2 - i12*mu3 + (i12*mu3*R3*s)/s0]
var = [s, i1, i2, i12]

RN, rts, spe, alp, bet, gam = ODE2RN(RHS1, var)
result = bdAn(RN, rts, var)
RHS = result['RHS']

diff = [simplify(RHS[i] - RHS1[i]) for i in range(len(RHS1))]
print("RHS - RHS1 =", diff)  # Should be [0, 0, 0, 0]
"""
