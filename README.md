# Replication Package: Digital Media Use and Cognitive Ability

This repository contains the replication code and data preparation pipeline for the empirical analysis in:

> **Linarte Teran, J., Campbell, A.C., Xu, J., Li, S.** (2026). *Beyond the Screen: Disentangling Causal Effects and Familial Confounding in Adolescent Screen Time and Cognition*.

---

## 1. Project Overview

This study investigates whether screen time, across 6 distinct domains, exerts a causal effect on cognitive ability, or whether observed associations are driven by unobserved familial and genetic confounding.

Using monozygotic (MZ) twin pairs from the German **TwinLife** panel study, we implement the **ICE FALCON** (*Inference on Causality with Exchangeability: Familial Association and Lower-level Conditioning*) framework. By conditioning individual cognitive outcomes on co-twin digital media exposures (and vice versa for reverse causality), the model decomposes observational associations into direct causal parameters versus shared familial/genetic background.

---

## 2. Sample & Data Access

### Analytical Sample
* **Sample Size:** $N = 1{,}942$ individuals ($971$ intact monozygotic twin pairs).
* **Target Population:** Adolescent and young adult twins from **Cohorts 2, 3, and 4** (aged 10–25 at Wave 1).

### Obtaining the Raw Data
In compliance with data protection laws and the GESIS user agreement, raw TwinLife microdata cannot be hosted directly in this repository. Researchers can obtain access free of charge for scientific replication:

1. Register at the **GESIS Data Archive**
2. Request the **TwinLife Scientific Use File (SUF)**
3. Download the following required Stata files:
   * `ZA6701_family_wide_wid1_v9-0-0.dta` (Wave 1)
   * `ZA6701_family_wide_wid3_v9-0-0.dta` (Wave 3)
   * `ZA6701_family_wide_wid7_v9-0-0.dta` (Wave 7)
   * `ZA6701_weights_v9-0-0.dta` (Survey weights)

---

## 3. Software Requirements & System Specifications

* **Statistical Software:** Stata 17 or higher (tested on Stata/SE 17 and Stata/MP 18).
* **Dependencies:** None.
* **Hardware & Runtime Notice:**
  * Estimating 12 models across 6 screen time domains (Forward and Reverse specifications) with **10,000 cluster-bootstrap resamples** is computationally intensive.
  * Expected runtime: **~2 to 4 hours** on a standard multi-core desktop.
  * For testing purposes, you may adjust `local NumberOfBootstrapResamples = 50` in `02_ICEFALCON_Analysis.do` to verify execution in under 2 minutes.
* **Reproducibility:** The cluster-bootstrap pseudo-random number generator is initialized with specific seeds to replicate exact empirical standard errors and $p$-values.

---
