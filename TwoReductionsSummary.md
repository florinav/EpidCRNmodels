# Two Reduction Strategies for John's Two-Strain SIR Model

## Summary

Successfully illustrated both reduction regimes for the two-strain SIR model, validating the Feliu-Lax-Walcher-Wiuf (FLWW) asymptotic reduction theory.

## Model

**Full 5D System:**
- Variables: {s, i1, i2, r1, r2}
- 14 reactions (from ODE2RN conversion)

**Two Maximal Eliminable Sets Found:**
1. {i1, i2} → reduces to 3D system in {s, r1, r2}
2. {s, r1, r2} → reduces to 2D system in {i1, i2}

Both satisfy FSW conditions (noninteracting, U-linear, spanning tree).

## Regime 1: Fast Infection Dynamics

**Assumption:** Infection/recovery rates >> demographic turnover

**Parameters:**
- Fast: β1=2.0, β2=1.5, γ1=1.0, γ2=0.8 (infection/recovery)
- Slow: λ=0.05, μ=0.002, θ1=0.05, θ2=0.04 (demography/waning)

**Reduction Strategy:**
- Eliminate {i1, i2} via QSS assumption (i1' = 0, i2' = 0)
- Reduced 3D system in slow variables {s, r1, r2}
- Fast variables reconstructed: i1 = φ1(s, r1, r2), i2 = φ2(s, r1, r2)

**Result:**
- Blue (Full 5D) and Red (Reduced 3D) curves overlap for s(t), r1(t), r2(t)
- Reconstructed i1(t), i2(t) match full system
- ✓ Confirms asymptotic validity in fast-infection regime

## Regime 2: Fast Demographic Turnover

**Assumption:** Birth/death/waning rates >> infection/recovery

**Parameters:**
- Fast: λ=2.0, μ=0.8, θ1=1.0, θ2=0.8 (demography/waning)
- Slow: β1=0.2, β2=0.15, γ1=0.1, γ2=0.08 (infection/recovery)

**Reduction Strategy:**
- Eliminate {s, r1, r2} via QSS assumption (s' = 0, r1' = 0, r2' = 0)
- Reduced 2D system in slow variables {i1, i2}
- Fast variables reconstructed: s = ψ1(i1, i2), r1 = ψ2(i1, i2), r2 = ψ3(i1, i2)

**Result:**
- Blue (Full 5D) and Red (Reduced 2D) curves overlap for i1(t), i2(t)
- Reconstructed s(t), r1(t), r2(t) match full system
- ✓ Confirms asymptotic validity in fast-demography regime

## Key Insights from FLWW Paper

1. **Parameter-dependent validity:** Same elimination set valid in different parameter regimes
   - Not just structural (FSW conditions) but also quantitative (blanket conditions i, ii)

2. **Tikhonov-Fenichel parameter values (TFPV):**
   - Blanket condition (i): B₀(x) invertible ⟺ spanning trees exist (Lemma 11)
   - Blanket condition (ii): Critical manifold is stationary ⟺ cycle conditions (Lemma 12)

3. **Graphical criteria:**
   - G_U multidigraph encodes valid parameter regimes
   - Corollary 13: Setting rates of ∗→U reactions to zero guarantees TFPV

4. **Biological interpretation:**
   - Regime 1: "Epidemic on fast timescale" → track slow demographic changes
   - Regime 2: "Demography on fast timescale" → track slow epidemic spread

## Comparison to Saez-Feliu-Wiuf

| Aspect | SFW (2017) | FLWW (This paper) |
|--------|------------|-------------------|
| **What it does** | Identifies eliminable sets algebraically | Validates QSS asymptotically |
| **Conditions** | Noninteracting, U-linear, spanning tree | + Blanket conditions (i), (ii) |
| **Output** | Reduced network structure | Parameter regimes where reduction valid |
| **Validation** | Algebraic (steady state) | Asymptotic (Tikhonov-Fenichel) |

## Files Created

1. `illustrateTwoReductions.nb` - Main comparison script
2. `testBothElim.nb` - Verifies both {i1,i2} and {s,r1,r2} satisfy FSW
3. `testJohnMaxElim.nb` - Finds all maximal eliminable sets

## Conclusion

The visualizations confirm that:
- Different reduction strategies are valid in different parameter regimes
- Choice depends on **biological timescale separation assumption**
- FLWW theory provides rigorous mathematical foundation for QSS reduction
- Both reductions preserve slow dynamics while approximating fast dynamics

**Practical recommendation:**
- Use Regime 1 reduction when studying long-term demographic/waning effects
- Use Regime 2 reduction when studying disease invasion/outbreak dynamics
