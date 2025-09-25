# ANOVAtoolbox

This repository contains a collection of functions for performing various statistical analyses and generating visualizations. The functions are designed to work with different types of data and provide comprehensive outputs for data analysis.

## Functions Included
The functions are listed in alphabetical order:

### 1. `ancova_analysis` - ANCOVA Analysis
This function performs an ANCOVA analysis by comparing a numeric dependent variable across groups while controlling for a covariate. It checks the homogeneity of regression slopes, assesses the normality of residuals (via Shapiro-Wilk test and Q-Q plot), tests the homogeneity of variances using Levene's test, calculates partial eta squared effect sizes, computes estimated marginal means (EMMs) with 95% confidence intervals, and generates several diagnostic plots.

### 2. `anova_bin` - One-way ANOVA for Binary Response Variable
This function performs a one-way ANOVA analysis for a binary response variable. It fits a logistic regression model with the specified response variable and group variable(s) and performs deviance ANOVA to assess the significance of the group effect. The function also computes the effect size (odds ratios) and confidence intervals, and generates a bar plot showing the proportion of the response variable by group.

### 3. `anova_count` - Deviance ANOVA for Count Data
This function performs a deviance analysis of variance (ANOVA) for count data using Poisson regression. It takes a data frame containing the variables of interest and performs the ANOVA analysis. The function provides the deviance ANOVA table, post-hoc test results, effect size (incidence rate ratios), confidence intervals, and a bar plot of the response by group.

### 4. `anova_glm` - GLM ANOVA
This function performs a Generalized Linear Model (GLM) ANOVA. It fits the GLM model to the data and allows the user to review diagnostics such as residuals, ANOVA test results, effect sizes, confidence intervals, and plots of the effects.

### 5. `anova_kw` - Kruskal-Wallis Test and Post-Hoc Analysis with Boxplot Visualization
This function conducts a Kruskal-Wallis test and performs post-hoc analysis for comparing the distribution of a numeric response variable across different groups defined by one or more categorical variables. It also generates a boxplot to visualize the distribution of the response variable by group.

### 6. `anova_rm` - Repeated Measures ANOVA with Diagnostics and Effect Size Calculation
This function fits a repeated measures analysis of variance (ANOVA) model to data in long format, using one or more within-subject (and optionally between-subject) factors. It uses `afex::aov_ez` for model fitting, assesses model assumptions (normality and sphericity), computes effect sizes (partial eta²), and generates diagnostic plots (residuals vs index, normal Q-Q, and LM diagnostic plots). Optionally, if an estimated marginal means specification is provided, the function computes and plots estimated marginal means.

### 7. `anova_welch` - Welch's ANOVA and Post-hoc Tests
This function performs Welch's Analysis of Variance (ANOVA) and post-hoc tests for group comparisons. It calculates means, standard deviations, confidence intervals, effect sizes, and generates a bar plot with means and confidence intervals.

### 8. `manova_analysis` - MANOVA Analysis
This function performs a multivariate analysis of variance (MANOVA) for comparing multiple dependent variables across groups. It conducts assumption tests (including Mardia's multivariate normality test and Box's M test for homogeneity of covariance matrices) fits the MANOVA model, runs follow-up univariate ANOVAs, computes partial eta squared effect sizes, calculates estimated marginal means with confidence intervals, and generates diagnostic plots and a canonical discriminant analysis plot.

## Use Table

| Function            | When to Use                                                                                                                   | Assumptions                                                                                                                                         |
|---------------------|-------------------------------------------------------------------------------------------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------|
| **ancova_analysis** | When comparing a numeric dependent variable across groups while controlling for a covariate, with verification of ANCOVA assumptions. | - Dependent variable and covariate must be numeric.<br>- Independent variable must have at least two levels.<br>- Assumptions on slopes, normality, and homogeneity of variances are evaluated. |
| **anova_bin**       | When analyzing binary response variables and comparing group effects using logistic regression.                                | - Response variable is binary (0/1).<br>- Observations are independent.<br>- Assumes logistic model conditions (e.g., linearity in the log odds). |
| **anova_count**     | When working with count data to assess group differences via Poisson regression.                                             | - Response consists of non-negative integer counts.<br>- Data follow a Poisson distribution (variance ≈ mean).<br>- Observations are independent. |
| **anova_glm**       | When performing an ANOVA within a generalized linear model framework for various data distributions.                         | - Assumptions depend on the chosen GLM family and link function (e.g., normality for Gaussian models).<br>- Observations are independent.         |
| **anova_kw**        | When the numeric response does not meet normality assumptions; ideal for comparing medians across groups using a non-parametric approach. | - Data are at least ordinal.<br>- Samples are independent.<br>- Group distributions should have similar shapes aside from median differences.    |
| **anova_rm**        | When data are in long format and you need to perform a repeated measures ANOVA analysis with diagnostics and effect size calculations. | - Data should be in long format.<br>- Within-subject factors must have more than one level (for sphericity checks if applicable).<br>- Observations are independent. |
| **anova_welch**     | When comparing group means in the presence of unequal variances; uses Welch’s ANOVA with post-hoc tests.                        | - Within-group data are approximately normally distributed.<br>- Observations are independent.<br>- Does not require equal variances across groups. |
| **manova_analysis** | When analyzing multivariate data with multiple dependent variables across groups, requiring follow-up univariate tests, effect sizes, and diagnostic evaluations. | - Multiple dependent variables should exhibit multivariate normality (or robustness is considered).<br>- Homogeneity of covariance matrices and independence are assumed. |

## Notes
Each function is accompanied by a detailed description of its parameters and the returned results. Please refer to the individual function documentation within the repository for more detailed information on each analysis and usage examples.
