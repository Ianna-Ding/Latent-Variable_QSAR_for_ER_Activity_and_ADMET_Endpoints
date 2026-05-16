# Mathematical Latent-Variable Analysis of ERα Activity and ADMET Endpoints

> A reproducible statistical-bioinformatics case study on high-dimensional biomedical descriptor data, focused on PCA, PLS, O2PLS-inspired data integration, goodness-of-fit diagnostics, and interpretable latent-variable modelling.

This repository revisits the breast-cancer-related ERα activity and ADMET dataset from Problem D of the 2021 China Post-Graduate Mathematical Contest in Modeling.  

The central question is:

> How much of ERα activity and ADMET behaviour can be explained by low-dimensional, interpretable latent structure in a high-dimensional molecular descriptor matrix?

The analysis uses molecular descriptors as a controlled case study for problems that also arise in omics integration:

- high-dimensional feature spaces with strong correlation and redundancy;
- latent biological or chemical factors hidden behind many observed variables;
- the need to separate shared structure, block-specific structure, and noise;
- model assessment beyond predictive accuracy;
- interpretable dimension reduction for biomedical data.

Although the dataset is QSAR rather than true multi-omics, the statistical structure is close to a two-block data-integration problem. This makes it a useful small-scale test bed for PCA, PLS, O2PLS-style decomposition, and future extensions toward probabilistic two-way PLS and bidimensionally linked matrix factorization.

---

## 1. Methodological Themes

This project is designed around four methodological themes.

### 1.1 High-dimensional biomedical data

The training set contains

$$n=1974 \quad \text{compounds}, \qquad p=729 \quad \text{molecular descriptors}.$$

The response variables are:

- continuous ERα activity: `IC50_nM`, `pIC50`;
- five binary ADMET endpoints: `Caco-2`, `CYP3A4`, `hERG`, `HOB`, `MN`.

The main descriptor matrix is denoted by
$$
X \in \mathbb{R}^{n \times p},
$$
where rows are compounds and columns are molecular descriptors. After removing constant descriptors, the primary modelling matrix has

$$
X_c \in \mathbb{R}^{1974 \times 504}.
$$
The response used for regression is

$$
y = \mathrm{pIC50} = 9 - \log_{10}(\mathrm{IC50}_{\mathrm{nM}}),
$$
so that larger values correspond to stronger ERα inhibitory activity.

### 1.2 Latent-variable modelling

The main methods are PCA, PLS regression, and O2PLS-inspired block integration. These methods are useful because they represent the data through latent scores:

$$
X \approx T P^\top,
$$
where `T` contains compound scores and `P` contains descriptor loadings. This gives a low-dimensional representation that can be interpreted geometrically and chemically.

The project also includes nonlinear predictive models in the last notebook as **complementary benchmarks**. They are included to put the latent-variable results in context and to check whether flexible prediction methods capture additional structure in the data.

> Do nonlinear models explain systematic residual structure left by a transparent latent-variable model?

The emphasis remains on interpretable modelling, diagnostics, and mathematical structure, while using predictive benchmarks as a useful reference point.

### 1.3 Model assessment and goodness-of-fit

The project does not stop at fitted values. It explicitly studies:

- cross-validated prediction error;
- Q² and RMSE;
- residual calibration;
- leverage and score-space outliers;
- PCA Hotelling $T^2$ and SPE/Q residuals;
- residual patterns across the observed activity range;
- sensitivity to unusual compounds;
- comparison between linear latent-variable models and nonlinear benchmarks.

### 1.4 Data integration as a bridge to omics

The descriptor matrix is split into chemically meaningful blocks:

$$
X = [X^{(1)}, X^{(2)}, X^{(r)}]
$$
where

- $X^{(1)}$: local electronic / atom-environment descriptors;
- $X^{(2)}$: topological connectivity / structural descriptors;
- $X^{(r)}$: remaining descriptors not assigned to the two main blocks.

This block structure is used to mimic a two-omics integration problem:

$$
X^{(1)} \leftrightarrow X^{(2)}
$$
The goal is to identify:

- latent structure shared between blocks;
- block-specific structure;
- residual variation;
- relationships between these components and ERα/ADMET endpoints.

---

## 2. Dataset overview

Raw files:

```text
data/raw/Molecular_Descriptor.xlsx
data/raw/ERα_activity.xlsx
data/raw/ADMET.xlsx
data/raw/molecular_descriptors_summary.xlsx
```

Training data:

| Object | Dimension / content |
|---|---:|
| Molecular descriptors | 1,974 compounds × 729 descriptors |
| ERα activity | IC50 and pIC50 |
| ADMET endpoints | 5 binary endpoints |
| Test compounds | 50 compounds |

ADMET endpoints:

| Endpoint | Interpretation |
|---|---|
| `Caco-2` | intestinal epithelial permeability / absorption |
| `CYP3A4` | metabolism-related cytochrome P450 endpoint |
| `hERG` | cardiac safety / cardiotoxicity-related endpoint |
| `HOB` | human oral bioavailability |
| `MN` | micronucleus genetic-toxicity endpoint |

This repository uses the training set for statistical modelling and interpretation. The test set is not treated as a source of biological validation because the emphasis is methodological.

---

## 3. Repository structure

```text
Latent-Variable_QSAR_for_ER_Activity_and_ADMET_Endpoints/
├── README.md
├── data/
│   └── raw/
│       ├── Molecular_Descriptor.xlsx
│       ├── ERα_activity.xlsx
│       ├── ADMET.xlsx
│       └── molecular_descriptors_summary.xlsx
├── code/
│   ├── 01_data_audit_and_problem_formulation.ipynb
│   ├── 02_pca_descriptor_space.ipynb
│   ├── 03_pls_er_activity.ipynb
│   ├── 04_model_diagnostics.ipynb
│   ├── 05_o2pls_inspired_descriptor_block_integration.ipynb
│   └── 06_nonlinear_benchmark_and_model_comparison.ipynb
├── figures/
├── results/
└── requirements.txt
```

---

## 4. Notebook map

| Notebook | Purpose | Statistical focus |
|---|---|---|
| `01_data_audit_and_problem_formulation.ipynb` | Audit raw data and define the modelling problem | alignment, missingness, descriptor screening, block definition |
| `02_pca_descriptor_space.ipynb` | Analyse descriptor-space geometry | PCA/SVD, explained variance, loadings, score-space outliers |
| `03_pls_er_activity.ipynb` | Model ERα pIC50 using latent components | PLS regression, component selection, VIP/loadings |
| `04_model_diagnostics.ipynb` | Assess model reliability | residuals, calibration, leverage, outlier sensitivity |
| `05_o2pls_inspired_descriptor_block_integration.ipynb` | Study shared and block-specific descriptor structure | two-block integration, CCA-style shared scores, block-wise PLS |
| `06_nonlinear_benchmark_and_model_comparison.ipynb` | Benchmark remaining nonlinear signal | model-misspecification check |

---

## 5. Mathematical formulation

### 5.1 Preprocessing and descriptor audit

Let

$$
X_{raw} \in \mathbb{R}^{1974 \times 729}
$$
be the original descriptor matrix. The audit step checks:

- row alignment across descriptor, activity, and ADMET files;
- missing values;
- duplicate SMILES;
- descriptor type;
- zero-variance descriptors;
- near-zero-variance descriptors;
- descriptor-category mapping.

Empirical audit summary:

| Quantity | Value |
|---|---:|
| Training compounds | 1,974 |
| Test compounds | 50 |
| Raw descriptors | 729 |
| Constant descriptors | 225 |
| Primary non-constant descriptors | 504 |
| Near-constant descriptors | 184 |
| Binary-like descriptors | 241 |
| Integer low-cardinality descriptors | 102 |
| Continuous float-like descriptors | 359 |

The primary matrix used for PCA/PLS is obtained by removing constant descriptors and standardising each remaining descriptor:

$$
\tilde X_{ij} = \frac{X_{ij} - \bar X_j}{s_j}
$$
This gives each descriptor comparable scale before latent-variable analysis.

---

### 5.2 PCA: unsupervised descriptor-space geometry

PCA approximates the standardised descriptor matrix by a rank-$k$ matrix:

$$
\tilde X \approx T_k P_k^\top,
$$
where

$$
T_k = \tilde X V_k.
$$
Equivalently, PCA solves

$$
\min_{\mathrm{rank}(M) \le k} \|\tilde X - M\|_F^2.
$$
The principal axes are obtained from the singular value decomposition

$$
\tilde X = U D V^\top.
$$
The score matrix $T_k = U_k D_k$ represents compounds in a low-dimensional latent space; the loading matrix $P_k=V_k$ indicates which descriptors drive each direction.

#### PCA results

For the primary 504-descriptor matrix:

| PCA summary | Explained variance |
|---|---:|
| PC1 | 19.1% |
| PC1–PC2 | 28.8% |
| PC1–PC5 | 46.1% |
| PC1–PC10 | 58.5% |
| PC1–PC20 | 72.9% |
| PC1–PC50 | 91.6% |

This confirms that the descriptor space is strongly low-dimensional relative to its original 504 non-constant descriptors.

#### PCA outlier diagnostics

Two complementary diagnostics are used:

1. **Hotelling score distance**

$$
T_i^2 = t_i^\top \Lambda_k^{-1} t_i
$$

where $t_i$ is the score vector for compound $i$.

2. **Squared prediction error / Q residual**

$$
\mathrm{SPE}_i = \|\tilde x_i - \hat x_i\|_2^2
$$

Empirical result:

| Outlier category | Count |
|---|---:|
| Not flagged | 1,939 |
| $T^2$-only | 15 |
| SPE-only | 15 |
| Both $T^2$ and SPE | 5 |
| Total PCA-flagged compounds | 35 |

This separates compounds that are extreme **within** the principal-component subspace from compounds that are poorly reconstructed by the low-rank PCA model.

---

### 5.3 PLS regression: supervised latent variables for ERα activity

PLS regression searches for latent score directions that explain both descriptor variation and the response $y$. For a univariate response, the first PLS weight vector can be written as

$$
w_1
= \arg\max_{\|w\|=1} \operatorname{Cov}(\tilde X w, y)^2
= \arg\max_{\|w\|=1} w^\top \tilde X^\top yy^\top \tilde X w
$$
The latent score is

$$
t_1 = \tilde X w_1
$$
After deflation, additional components are extracted. The final model is

$$
y = T q + e
$$
where $T$ contains PLS scores and $q$ contains regression coefficients on the latent scores.

The number of components is selected using cross-validation. Predictive performance is summarised by

$$
Q^2 = 1 - \frac{\sum_i (y_i - \hat y_i^{CV})^2}{\sum_i (y_i - \bar y)^2}
$$

#### PLS results

For the primary 504-descriptor matrix:

| Quantity | Value |
|---|---:|
| Best mean-fold RMSE | 43 components |
| One-SE selected model | 20 components |
| Cross-validated RMSE | 0.856 |
| Cross-validated MAE | 0.641 |
| Cross-validated Q² | 0.638 |
| Training R² | 0.737 |

The one-SE rule is preferred because it gives a more parsimonious model while staying close to the best cross-validated error.

#### Descriptor importance

Variable importance in projection is computed as

$$
\mathrm{VIP}_j
=
\sqrt{
 p \cdot
 \frac{
 \sum_{a=1}^{A} SSY_a \frac{w_{ja}^2}{\|w_a\|^2}
 }{
 \sum_{a=1}^{A} SSY_a
 }
}
$$
where $SSY_a$ is the amount of response variation explained by component $a$.

Top descriptors include hydroxyl-related descriptors, lipophilicity descriptors, ring descriptors, hydrogen-bond descriptors, and molecular distance edge descriptors, for example:

- `maxsOH`, `minsOH`, `maxHsOH`, `minHsOH`;
- `MDEC-23`;
- `MLogP`, `CrippenLogP`, `LipoaffinityIndex`;
- `MLFER_A`;
- `nT6Ring`, `n6Ring`;
- `nHBAcc`, `maxHBd`.

The interpretation is not that these individual variables are causal, but that the corresponding descriptor families define latent directions associated with ERα activity.

---

### 5.4 Model diagnostics and goodness-of-fit

A latent-variable model should be assessed both globally and locally.

#### Residual structure

The cross-validated residual is

$$
e_i^{CV} = y_i - \hat y_i^{CV}.
$$
A large-residual threshold is defined by the 95th percentile of $|e_i^{CV}|$:

$$
|e_i^{CV}| > 1.670.
$$
Empirical diagnostics:

| Diagnostic | Value |
|---|---:|
| Large CV residuals | 99 compounds |
| Proportion large residuals | 5.0% |
| Median absolute CV residual | 0.489 |
| CV RMSE | 0.856 |
| CV Q² | 0.638 |
| Calibration slope | 0.918 |

Residuals are not uniformly distributed across the response range. The model tends to over-predict weakly active compounds and under-predict strongly active compounds, which is typical shrinkage behaviour for regularised or latent low-rank regression.

#### Outlier sensitivity

Removing PCA-flagged compounds gives:

| Dataset | CV RMSE | CV Q² |
|---|---:|---:|
| All compounds | 0.856 | 0.638 |
| Excluding PCA outliers | 0.832 | 0.655 |

This suggests that unusual descriptor-space compounds contribute to prediction difficulty but do not fully explain model error.

---

### 5.5 O2PLS-inspired descriptor-block integration

The two main descriptor blocks are:

| Block | Description | Non-constant descriptors |
|---|---|---:|
| $X^{(1)}$ | local electronic / atom-environment descriptors | 335 |
| $X^{(2)}$ | topological connectivity / structural descriptors | 129 |

This repository does not claim to perform genuine multi-omics O2PLS, because both blocks are derived from the same molecular descriptor table. Instead, it uses an O2PLS-inspired formulation as a transparent statistical exercise.

A classical two-block latent decomposition can be written schematically as

$$
X^{(1)} = T W^\top + T_\perp W_\perp^\top + E
$$

$$
X^{(2)} = U C^\top + U_\perp C_\perp^\top + F
$$

where

- $T, U$: joint latent scores shared between blocks;
- $T_\perp, U_\perp$: block-specific latent scores;
- $E,F$: residual noise.

In this project, block PCA scores are first computed, then cross-block association is studied through canonical-correlation-style analysis and residualisation.

Let

$$
Z_1 = \mathrm{PCA}(X^{(1)}), \qquad Z_2 = \mathrm{PCA}(X^{(2)})
$$
where each contains the leading block-specific PC scores. Shared structure is estimated by directions $a,b$ solving

$$
\max_{a,b} \operatorname{Corr}(Z_1a, Z_2b)
$$
Block-specific residual structure is then approximated by

$$
R_1 = Z_1 - \hat Z_1(Z_2),
\qquad
R_2 = Z_2 - \hat Z_2(Z_1)
$$
This gives a practical approximation to the conceptual O2PLS decomposition:

$$
\text{observed block variation}
=
\text{shared variation}
+
\text{block-specific variation}
+
\text{residual noise}
$$

#### Block PCA results

| Block | PC1 | PC1–PC2 | PC1–PC5 | PC1–PC10 |
|---|---:|---:|---:|---:|
| Local electronic / atom-environment block | 14.0% | 23.9% | 41.4% | 56.0% |
| Topological connectivity block | 28.9% | 40.8% | 58.6% | 73.0% |

The topological block is more compressible than the local electronic block, while the electronic block retains more dispersed variation.

#### Cross-block association

The leading canonical correlations between block-PC spaces are:

| Component | Canonical correlation |
|---|---:|
| 1 | 0.984 |
| 2 | 0.944 |
| 3 | 0.932 |
| 4 | 0.814 |
| 5 | 0.795 |

Cross-block prediction in retained PC space gives:

| Direction | $R^2$ in retained PC space | Residual fraction |
|---|---:|---:|
| Predict local electronic block from topological block | 0.619 | 0.381 |
| Predict topological block from local electronic block | 0.656 | 0.344 |

This indicates substantial shared structure, but also a non-negligible block-specific component.

#### Relationship to ERα and ADMET endpoints

Selected latent-score associations:

| Latent score | Association with pIC50 |
|---|---:|
| Local electronic PC3 | Pearson \(-0.454\) |
| Joint average component 1 | Pearson \(+0.446\) |
| Joint average component 3 | Pearson \(-0.398\) |
| Topological PC1 | Pearson \(+0.347\) |

Selected ADMET associations:

| Endpoint / latent score | Association |
|---|---:|
| `CYP3A4` vs topological PC1 | \(+0.614\) |
| `Caco-2` vs joint average component 1 | \(-0.566\) |
| `hERG` vs joint average component 1 | \(+0.558\) |
| `CYP3A4` vs joint average component 1 | \(+0.558\) |

These results suggest that the shared descriptor structure is not merely technical redundancy; it is also related to biological activity and ADMET behaviour.

---

### 5.6 Block-wise PLS comparison

PLS models are fitted using different descriptor sets.

| Feature set | CV RMSE | CV Q² | Comment |
|---|---:|---:|---|
| Local electronic + topological blocks | 0.834 | 0.656 | best minimum-RMSE block model |
| Primary full non-constant matrix | 0.838 | 0.653 | similar to two-block model |
| Local electronic block only | 0.865 | 0.630 | stronger than topological alone |
| Topological block only | 0.898 | 0.602 | useful but less complete |
| Remaining unassigned descriptors | 1.025 | 0.481 | weaker independent signal |

This comparison supports the block-integration view: both local electronic and topological descriptors contribute to activity prediction, and their combination captures more information than either block alone.

---

### 5.7 Nonlinear benchmark

The final notebook compares the linear latent-variable model with nonlinear benchmark models such as random forests and histogram gradient boosting.

This comparison helps quantify how much predictive performance is gained when more flexible response surfaces are used, while keeping PCA/PLS/O2PLS-style scores and loadings as the primary objects for interpretation.

| Model | Feature representation | CV RMSE | CV Q² |
|---|---|---:|---:|
| Linear PLS baseline | primary descriptors | 0.856 | 0.638 |
| Histogram gradient boosting | primary descriptors | 0.689 | 0.766 |
| Random forest | primary descriptors | 0.707 | 0.752 |
| Histogram gradient boosting | PLS scores | 0.752 | 0.720 |
| RBF SVR | PLS scores | 0.761 | 0.714 |

Interpretation:

- nonlinear models improve prediction;
- therefore, some nonlinear signal remains after linear PLS;
- PCA/PLS/O2PLS-style models remain the important scientific objects because they provide loadings, scores, decompositions, and diagnostics;
- the benchmark results motivate future work on interpretable nonlinear or probabilistic latent-variable extensions.

---

## 6. Relation to probabilistic O2PLS and bidimensional integration

This project can be extended in two natural directions.

### 6.1 Probabilistic two-way PLS

The O2PLS-inspired analysis here is deterministic and exploratory. A probabilistic version would introduce explicit latent variables and noise models, for example:

$$
X^{(1)} = T W^\top + T_\perp W_\perp^\top + E
$$

$$
X^{(2)} = U C^\top + U_\perp C_\perp^\top + F
$$

with distributional assumptions such as

$$
E \sim \mathcal{N}(0, \sigma_1^2 I),
\qquad
F \sim \mathcal{N}(0, \sigma_2^2 I)
$$
This would allow likelihood-based model comparison, uncertainty quantification, and more formal goodness-of-fit testing.

Possible next steps:

- implement an explicit probabilistic two-block latent-variable model;
- compare deterministic O2PLS-style scores with probabilistic posterior scores;
- study component selection by likelihood, information criteria, cross-validation, and permutation testing;
- add simulation studies where the true joint and block-specific components are known.

### 6.2 Bidimensionally linked matrices

A future version can generalise from two descriptor blocks to a collection of matrices indexed by two dimensions, for example:

$$
X_{ij}, \qquad i = 1,\ldots,I, \quad j = 1,\ldots,J.
$$
A bidimensional factorisation can be written schematically as

$$
X_{ij}
=
G_{ij}
+
R_{ij}
+
C_{ij}
+
I_{ij}
+
E_{ij}
$$
where

- $G_{ij}$: globally shared structure;
- $R_{ij}$: row-shared structure;
- $C_{ij}$: column-shared structure;
- $I_{ij}$: individual matrix-specific structure;
- $E_{ij}$: residual noise.

For biomedical data, such a design could represent multiple omics platforms measured across multiple tissues, time points, disease groups, or experimental conditions.

This repository currently provides the first layer of that direction: careful matrix auditing, PCA/PLS score construction, two-block decomposition, residual diagnostics, and endpoint association analysis.

---

## 7. Main empirical findings

| Theme | Finding |
|---|---|
| Descriptor audit | 225 of 729 descriptors are constant; the primary non-constant matrix has 504 descriptors. |
| PCA | 20 PCs explain about 72.9% of descriptor variance. |
| PCA diagnostics | 35 compounds are flagged by score-distance and/or reconstruction-error diagnostics. |
| PLS | A 20-component one-SE PLS model gives CV RMSE 0.856 and Q² 0.638. |
| Descriptor interpretation | Hydroxyl, lipophilicity, ring, hydrogen-bond, and molecular-distance descriptors are prominent. |
| Residual diagnostics | Errors are larger at activity extremes, suggesting shrinkage and possible nonlinear structure. |
| Block integration | Local electronic and topological blocks share strong latent structure but also retain block-specific variation. |
| ADMET association | Shared latent scores are related to ADMET endpoints such as Caco-2, CYP3A4, and hERG. |
| Nonlinear benchmark | Nonlinear models improve Q², indicating remaining nonlinear signal. |

---

## 8. Reproducibility

The analysis is organised to be rerunnable from raw files.

Recommended execution order:

```bash
01_data_audit_and_problem_formulation.ipynb
02_pca_descriptor_space.ipynb
03_pls_er_activity.ipynb
04_model_diagnostics.ipynb
05_o2pls_inspired_descriptor_block_integration.ipynb
06_nonlinear_benchmark_and_model_comparison.ipynb
```

Typical Python dependencies:

```text
numpy
pandas
scipy
scikit-learn
matplotlib
openpyxl
```

The current implementation is in Python. A useful future step is to re-implement the core PCA/PLS/O2PLS-style functions in R, add unit tests, and package the workflow more formally for statistical-bioinformatics reuse.

---

## 9. Limitations

This project is deliberately cautious about interpretation.

1. The dataset is not true multi-omics data. The two-block integration analysis is an O2PLS-inspired descriptor-block analysis.
2. Descriptor importance should not be interpreted as causal chemical evidence.
3. ERα/ADMET predictions are computational and do not validate drug efficacy or safety.
5. More formal probabilistic modelling is needed for likelihood-based inference and uncertainty quantification.

---

## 10. Planned extensions

The next improvements are methodological rather than purely predictive:

- implement a closer O2PLS / probabilistic O2PLS model;
- compare deterministic and probabilistic latent scores;
- add formal goodness-of-fit tests for low-rank latent-variable models;
- construct simulation studies with known joint and block-specific structure;
- extend the two-block design toward bidimensionally linked matrix factorization;
- port the core workflow to R and package reusable functions;
- apply the same pipeline to a genuine public multi-omics breast-cancer dataset.

---

## 11. Keywords

```text
statistical bioinformatics
latent-variable modelling
high-dimensional data integration
PCA
SVD
PLS regression
O2PLS
probabilistic PLS
bidimensionally linked matrices
model diagnostics
goodness-of-fit
QSAR
ERα
ADMET
breast cancer
interpretable modelling
```

---

## 12. Disclaimer

This repository is for statistical modelling and research-training purposes. The results should not be interpreted as experimental validation of any compound's biological, pharmacological, or clinical effectiveness.
