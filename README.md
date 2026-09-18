# Cross-Validated Regression Analysis of Functional Brain Network Activity and Dimensional Psychopathology

R pipeline for leave-one-out cross-validated (LOOCV) regression of brain network activity (ToM & Pain) from dimensional psychopathology (ASR) (simulated data because of sensitive patient data).

--- 
# Overview

The primary objective of this workflow is to estimate the out-of-sample predictive performance ($R^2_{\text{CV}}$) of psychopathological symptom dimensions on neural network activity, avoiding optimistic in-sample fit estimates.

Two predictor sets are evaluated against functional activation in two distinct networks:
1. **ASR DSM-Oriented Scales**: Depressive Problems, Anxiety Problems, Somatic Problems, Avoidant Personality Problems, Inattention, Hyperactivity, and Antisocial Personality Problems. *(Note: The combined ADH Problems scale is excluded to prevent perfect collinearity with the Inattention and Hyperactivity subscales).*
2. **ASR Broad-Band Overscales**: Internalizing, Externalizing, and Thought Problems.

Target outcomes:
* **ToM Network Activity** (Region-of-Interest contrast estimates / beta weights)
* **Pain Network Activity** (Region-of-Interest contrast estimates / beta weights)

---

# Statistical Methodology

Models are fit using Generalised Linear Models (`glm` with Gaussian family / OLS). Validation is carried out via **Leave-One-Out Cross-Validation (LOOCV)** using the `boot` package.

Out-of-sample performance is quantified as:

$$R^2_{\text{CV}} = 1 - \frac{\text{MSE}_{\text{CV}}}{\text{Var}(Y)}$$

Where:
* $\text{MSE}_{\text{CV}}$ represents the raw mean squared prediction error across all $N$ validation folds (`delta[1]`).
* $\text{Var}(Y)$ represents the sample variance of the respective target variable within the complete-case analysis cohort.

*Note on $R^2_{\text{CV}}$:* Whilst in-sample $R^2$ are positive, cross-validated $R^2$ values can be negative if model predictions on unseen folds perform worse than the sample mean baseline ($\bar{Y}$), indicating absence of generalisable signal or overfitting.

---

# Reproduction & Usage

The script is self-contained. To ensure direct reproducibility without distributing sensitive patient data, the pipeline automatically generates a synthetic sample dataset ($N = 150$) that matches the exact distributional properties and schema expected by the models.

# Prerequisites

* **R version**: $\ge 4.1.0$
* **Required packages**: `dplyr`, `boot`, `tibble`

Install missing dependencies via R:
```r
install.packages(c("dplyr", "boot", "tibble"))
