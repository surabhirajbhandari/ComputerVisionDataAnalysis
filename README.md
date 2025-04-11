
# Computer Vision Target Detection — Statistical Data Analysis in MATLAB

This project analyzes the statistical performance of a machine vision system for target detection using MATLAB.  
It combines ROC curve analysis, hypothesis testing, and signal processing to distinguish between “target present” and “target absent” scenarios — with the goal of identifying optimal thresholds and improving system accuracy.

---

##  Project Overview

In real-world sensor systems (SONAR, RF, IR, LIDAR), uncertainty from noise and environmental factors often causes overlap between "target present" and "target absent" data. This project explores how statistical analysis and signal processing techniques can be used to improve target detection accuracy.

The analysis is divided into three parts:

###  Part 1: ROC Analysis  
- Empirical probability densities estimated for both target states  
- Optimal threshold identified via Youden’s Index  
- ROC curve plotted and AUC computed to measure classifier performance  
- Confusion matrix, error rate, and PPV evaluated at multiple thresholds

###  Part 2: Hypothesis Testing & Bootstrapping  
- Best-fit probability distributions identified using Chi-Squared tests  
- Gamma (for absent) and Rician (for present) distributions selected  
- Parametric and non-parametric bootstrapping performed (5000 trials)  
- Confidence intervals for AUC computed  
- Theoretical vs. empirical ROC curves visualized

###  Part 3: Signal Processing & Performance Enhancement  
- Synthetic datasets generated from fitted distributions  
- Applied signal processing techniques:
  - Arithmetic mean  
  - Geometric mean  
  - Maximum  
- ROC analysis repeated over 100 trials for each technique  
- Final performance metrics:
  - Mean AUC  
  - Error rate  
  - PPV  
  - Custom Performance Index (PI)

---

## Key Outcomes

| Metric                        | Baseline     | After Signal Processing |
|------------------------------|--------------|--------------------------|
| **Initial AUC**              | 0.906        | **0.985** (Arithmetic)  |
| **Error Rate**               | 21 / 130     | **7 / 130**             |
| **Positive Predictive Value**| 0.80         | **0.934**               |
| **Best Fit Distributions**   | Gamma (H₀), Rician (H₁) | ✅ Confirmed via Bootstrapping |

---
---

## Tools & Techniques

- **Language**: MATLAB
- **Stats**: ROC Curves, AUC, Confusion Matrices, PPV, Error Rates
- **Distributions**: Gamma, Rician, Lognormal, Rayleigh, Weibull, etc.
- **Analysis Methods**:
  - Chi-Squared Goodness-of-Fit
  - Parametric & Non-Parametric Bootstrapping
  - Kernel Density Estimation
  - Neyman-Pearson Euclidean Optimization
- **Signal Processing Techniques**: Arithmetic Mean, Maximum, Geometric Mean

---

## How to Run

1. Open MATLAB and clone this repo or download the `.zip`  
2. Open the matlab files  
3. Ensure `subu.xlsx` is in the same directory  - This is the dataset
4. Run the script to generate:
   - ROC Curves  
   - Confusion & Transition Matrices  
   - Density plots with shaded error regions  
   - Summary performance metrics in command window  

---



---

