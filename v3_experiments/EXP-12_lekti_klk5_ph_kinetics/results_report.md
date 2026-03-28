# EXP-12: LEKTI-KLK5 pH-Dependent Kinetics — Results Report

## Abstract

This report documents the outcome of EXP-12, which aimed to compute the pH-dependent inhibition constant (Ki) for the LEKTI-KLK5 interaction across the physiologically relevant pH range (5.5–7.4). The experiment is classified as **INCONCLUSIVE** due to two independent blocking factors: (1) no LEKTI-KLK5 complex structure exists in the PDB, and (2) the pipeline lacks constant-pH molecular dynamics methodology required to model pH-dependent protonation equilibria during simulation.

## Introduction / Background

The LEKTI-KLK5 interaction is pH-dependent, with a dramatic shift in inhibitory potency between acidic (pH 5.5, stratum corneum) and neutral (pH 7.4, deeper epidermis) conditions. Deraison et al. (2007) demonstrated that Ki shifts by 10–100× between pH 5.5 and pH 7.4, providing a biochemical mechanism for the spatial regulation of KLK5 activity across skin layers.

This pH dependence arises from the protonation states of key histidine residues at the LEKTI-KLK5 interface and in the KLK5 active site. At low pH, histidine protonation alters the electrostatic complementarity of the interface, weakening the interaction and releasing KLK5 to perform its desquamation function.

Modeling this pH-dependent behavior computationally requires either constant-pH MD — which dynamically adjusts protonation states during simulation — or carefully designed comparative simulations at different fixed protonation states with prior pKa prediction.

## Hypothesis

**H-12:** The pipeline can reproduce the experimentally observed 10–100× shift in Ki for LEKTI-KLK5 binding between pH 5.5 and pH 7.4.

## Methods

**Planned approach:** Perform umbrella sampling calculations of ΔG_bind at both pH 5.5 and pH 7.4 conditions, either using constant-pH MD or comparative fixed-protonation-state simulations with pKa-predicted protonation states at each pH.

**Pipeline components:** OpenMM 8.1, AMBER14/SB force field, TIP3P water model, PDBFixer (pH 7.0 protonation only).

**Status:** Methods could not be executed due to (1) missing LEKTI-KLK5 complex structure and (2) absence of constant-pH MD methodology.

## Controls

Not applicable — experiment was not executed.

## Results

**Status: INCONCLUSIVE — experiment could not be executed.**

This experiment is blocked by two independent factors:

1. **Missing structural data:** No LEKTI-KLK5 complex structure exists in the PDB (same structural limitation as EXP-03). Without a validated starting complex, no binding free energy calculations can be performed at any pH.

2. **Methodology limitation:** The pipeline uses fixed protonation states set during PDBFixer preparation at pH 7.0. pH-dependent simulations require either:
   - **Constant-pH MD:** Dynamically titrates residues during simulation, sampling protonation equilibria alongside conformational degrees of freedom. Not implemented in the current pipeline.
   - **Fixed-protonation comparative approach:** Separate simulations at different fixed protonation states, guided by pKa prediction (e.g., PROPKA or H++). While conceptually simpler, this approach requires careful validation and introduces systematic uncertainty from static protonation assumptions.

   Neither approach is within the current pipeline's validated capabilities.

**Feature F-12 (pH-dependent Ki for LEKTI-KLK5):** Not measured.

**Benchmark comparison:** Not possible.

## Discussion

pH-dependent protein-protein interactions present a significant computational challenge because they couple protonation equilibria (quantum mechanical in nature) with conformational dynamics (classical). The standard MD approach of assigning fixed protonation states at a single pH is fundamentally inadequate for capturing the continuous modulation of binding affinity across a pH gradient.

**Constant-pH MD** methods have been developed in several MD packages:
- **Discrete constant-pH MD** (Mongan et al., 2004): Uses Monte Carlo protonation state sampling at intervals during MD, implemented in Amber.
- **Continuous constant-pH MD** (Wallace & Shen, 2011): Uses λ-dynamics to treat protonation as a continuous variable, available in CHARMM.
- **GPU-accelerated constant-pH MD** (Harris et al., 2019): Recent implementations in OpenMM/Amber that make constant-pH simulations tractable for larger systems.

The current pipeline does not implement any of these methods. PDBFixer assigns protonation states at a single pH (7.0) and these remain fixed throughout the simulation. To model pH 5.5 behavior, one would need to re-protonate all titratable residues (particularly histidines) according to their predicted pKa values at pH 5.5, which requires pKa prediction tools (PROPKA, H++) not currently integrated into the pipeline.

The biological significance of the pH-dependent LEKTI-KLK5 interaction makes this an important target for future pipeline extensions, but it requires both structural data and methodological developments beyond the current scope.

## Conclusions

EXP-12 is INCONCLUSIVE. pH-dependent LEKTI-KLK5 kinetics cannot be computed because (a) no LEKTI-KLK5 complex structure exists in the PDB, and (b) the pipeline lacks constant-pH MD methodology required to model pH-dependent protonation equilibria. This experiment requires both the structural determination of a LEKTI-KLK5 complex and the implementation of constant-pH MD or equivalent methodology.

## Figures

No figures generated — experiment not executed due to missing structural data.

## References

1. Deraison, C., et al. (2007). LEKTI fragments specifically inhibit KLK5, KLK7, and KLK14 and control desquamation through a pH-dependent interaction. *Molecular Biology of the Cell*, 18(9), 3607–3619.
2. Mongan, J., Case, D. A., & McCammon, J. A. (2004). Constant pH molecular dynamics in generalized Born implicit solvent. *Journal of Computational Chemistry*, 25(16), 2038–2048.
3. Wallace, J. A., & Shen, J. K. (2011). Continuous constant pH molecular dynamics in explicit solvent with pH-based replica exchange. *Journal of Chemical Theory and Computation*, 7(8), 2617–2629.
4. Harris, R. C., Shen, J., & Case, D. A. (2019). GPU-accelerated implementation of continuous constant pH molecular dynamics in Amber. *Journal of Chemical Theory and Computation*, 15(11), 6284–6298.

---

**Author:** Ryan Kamp
**Affiliation:** Dept. of Computer Science, University of Cincinnati
**Email:** kamprj@mail.uc.edu
**GitHub:** ryanjosephkamp
