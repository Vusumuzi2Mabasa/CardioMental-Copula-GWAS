"""Generate a research manuscript from audited simulation outputs."""
import csv,json
from prepare import V,ROOT,dump,sha
from model import SCENARIOS
from analyse import LABELS,METHODS

def read(name):
    with (V/'results'/name).open() as f:return list(csv.DictReader(f))
def pct(x):return f'{100*float(x):.1f}%'
def mdtable(headers,rows):return '\n'.join(['| '+' | '.join(headers)+' |','| '+' | '.join(['---']*len(headers))+' |']+['| '+' | '.join(map(str,r))+' |' for r in rows])
def main():
    audit=json.loads((V/'results/audit.json').read_text());assert audit['passed'] and audit['summary_checks_passed']
    agg=read('aggregate_rates.csv');pairs=list(csv.DictReader((V/'reference/causal_pairs.csv').open()));phen=read('phenotypes.csv');cop=read('copula_diagnostic.csv');truth=read('generating_genetic_covariance.csv')
    def get(reg,sc,method):return next(r for r in agg if r['regime']==reg and r['scenario']==sc and r['method']==method)
    def rate(reg,sc,method):return pct(get(reg,sc,method)['rate'])
    def interval(reg,sc,method):
        r=get(reg,sc,method);p=float(r['rate']);n=int(r['attempted'])
        if p==0:return f"0.0% (evaluation upper bound {pct(1-.025**(1/n))})"
        if p==1:return f"100.0% (evaluation lower bound {pct(.025**(1/n))})"
        return f"{pct(r['rate'])} ({pct(r['bootstrap_lower'])} to {pct(r['bootstrap_upper'])})"
    def meanphen(reg,key):
        rows=[r for r in phen if r['regime']==reg];return sum(float(r[key]) for r in rows)/len(rows)
    causal_table=mdtable(['Population','Region GRCh37','Reference donors','Study donors','Distinct pair study r²'],[[r['population'],['22:20.0–20.1 Mb','22:30.0–30.1 Mb','22:40.0–40.1 Mb'][int(r['region'])],40,68 if r['population']=='YRI' else 59,f"{float(r['study_r2_high']):.3f}"] for r in pairs])
    phenotype_table=mdtable(['Phenotype summary','Ideal','Realistic'],[
        ['Complete observations',f"{meanphen('ideal','n_complete'):.0f}",f"{meanphen('realistic','n_complete'):.0f}"],
        ['Observed depression proxy',pct(meanphen('ideal','prevalence')),pct(meanphen('realistic','prevalence'))],
        ['Pressure mean in illustrative mmHg',f"{meanphen('ideal','pressure_mean'):.1f}",f"{meanphen('realistic','pressure_mean'):.1f}"],
        ['Pressure SD in illustrative mmHg',f"{meanphen('ideal','pressure_sd'):.1f}",f"{meanphen('realistic','pressure_sd'):.1f}"],
        ['Treated fraction',pct(meanphen('ideal','treated_fraction')),pct(meanphen('realistic','treated_fraction'))],
        ['Missing fraction',pct(meanphen('ideal','missing_fraction')),pct(meanphen('realistic','missing_fraction'))]])
    def result_table(reg):
        return mdtable(['Generating architecture','Joint association','Tagging conjunction','Coloc PP4 > 0.8','Calibrated coloc'],[[LABELS[s]]+[rate(reg,s,m) for m in [METHODS[0],METHODS[1],METHODS[3],METHODS[4]]] for s in SCENARIOS])
    ci_table=mdtable(['Realistic architecture','PP4 > 0.8 with 95% interval','Calibrated coloc with 95% interval'],[[LABELS[s],interval('realistic',s,METHODS[3]),interval('realistic',s,METHODS[4])] for s in ['shared_same','shared_opposing','distinct_low','distinct_high']])
    coprows=[]
    for s in ['null','cardio_only','shared_same','distinct_high']:
        rows=[r for r in cop if r['regime']=='realistic' and r['scenario']==s]
        coprows.append([LABELS[s],sum(int(r['attempted']) for r in rows),pct(sum(float(r['copula_nominal_rate']) for r in rows)/6),pct(sum(float(r['marginal_anchor_joint_rate']) for r in rows)/6)])
    cop_table=mdtable(['Realistic architecture','Datasets','Copula anchor p < 0.05','Marginal joint anchor p < 0.05'],coprows)
    truth_table=mdtable(['Architecture','Shared direct variants','Genetic predictor covariance','Same variant component'],[[LABELS[s],int(float(next(r for r in truth if r['scenario']==s)['shared_direct_variants'])),f"{sum(float(r['genetic_predictor_covariance']) for r in truth if r['scenario']==s)/6:.4f}",f"{sum(float(r['same_variant_diagonal']) for r in truth if r['scenario']==s)/6:.4f}"] for s in ['null','shared_same','shared_opposing','distinct_low','distinct_high']])
    text=r'''# Shared genetic effects on cardiometabolic and mental health traits under African linkage disequilibrium

**Vusumuzi Mabasa**

**Simulation study.** Every phenotype, disease indicator and genetic effect is simulated. Public 1000 Genomes haplotypes supply local LD only. This work contains no patient discoveries, AWI-Gen analysis or claim that shared genetic causes have been established in a real population.

[Word manuscript](report/CardioMental_Copula_GWAS_Report_v3.docx) · [Previous oracle study](https://github.com/Vusumuzi2Mabasa/CardioMental-Copula-GWAS/blob/a648ebd/README.md) · [Frozen version three protocol](v3/protocol/design.json) · [Complete results](v3/results/aggregate_rates.csv)

<!-- manuscript -->
## Abstract

Cardiometabolic and mental disorders can co-occur because of shared genetic effects, distinct correlated variants, environmental dependence or combinations of these mechanisms. We evaluated whether joint association, conjunction testing and Bayesian colocalisation distinguish these explanations under African local LD. Public Yoruba and Luhya haplotypes generated 31,200 synthetic datasets across three chromosome 22 regions, eight genetic architectures and two phenotype regimes. The study used 12,000 calibration datasets and 19,200 independent evaluation datasets, each with 3,200 generated participants. Quantitative blood pressure and binary depression proxies incorporated treatment, measurement error, misclassification and missing observations in the realistic regime. Primary colocalisation used separate participant halves. Under realistic concordant sharing, PP4 > 0.8 declared sharing in @SHARED_RAW@ of datasets and composite-null calibration did so in @SHARED_CAL@. Under distinct high-LD causes, these rates were @DISTINCT_RAW@ and @DISTINCT_CAL@. Two study reference pools contained perfectly correlated causal candidates, making shared and distinct causal configurations observationally equivalent. The simulations support a distinction between detecting association with multiple phenotypes and identifying a shared causal variant. They establish methodological performance under known generating mechanisms, rather than shared causation of real clinical comorbidity.

## Introduction

The co-occurrence of cardiometabolic disease and mental ill-health motivates genetic studies of multimorbidity. A common phenotype correlation can nevertheless arise without any variant affecting both traits. Conversely, one causal variant may affect both phenotypes even when their total correlation is modest or their effects have opposing directions. Association with two diseases, colocalisation of genetic signals and a causal effect of one disease on another are therefore separate questions.

This distinction is relevant to research on shared cardiometabolic and mood-disorder genetics in African populations, including the programme described in the Wits announcement of Dr Vivien Chebii's Wellcome award [1]. The present independent simulation is thematically aligned with that question; it is not part of that grant and does not use its data. The 1000 Genomes Project provides public population-specific haplotypes for investigating one component of the problem, local linkage disequilibrium [2,3]. YRI and LWK are two sampled populations and do not represent all African genomic diversity.

The preceding version of this repository tested Gaussian/probit copula association procedures with independent common markers. Its oracle experiment showed conditional gains from correctly modelling dependence, without a general advantage or a resolved benefit from knowing nuisance parameters. Those results are preserved. The current extension changes the scientific target: it evaluates evidence for shared causal variants when nearby distinct causes can produce apparently shared associations.

## Methods

### Study design and preservation

The design was frozen in local commit 5999591 before the main simulations. Two populations, three regions and two phenotype regimes define 12 strata. Each stratum contains 200 independent evaluation replicates for each of eight architectures, giving 19,200 evaluation datasets. Five no-sharing architectures each contribute 200 independent calibration datasets per stratum, giving 12,000 calibration datasets. Calibration and evaluation use separate seed namespaces, and cutoffs were saved before evaluation began. All 31,200 datasets are retained in the published record. The earlier 92,000-dataset oracle experiment and the initial study remain separate; their denominators are not pooled with this extension.

### African reference haplotypes and causal configurations

We verified the complete public 1000 Genomes Phase 3 v5b chromosome 22 VCF and population-panel files against their SHA256 checksums [2,3]. Within YRI and LWK, a seeded split reserved 40 donors for external reference characterization and used the remaining 68 YRI and 59 LWK donors for generating study haplotypes. Donor pools were disjoint. Public outputs omit donor identifiers but retain derived haplotypes, population labels, split-list digests, source URLs and file checksums. These derived public data remain attributable to the original resource; omitting identifiers is not a claim of anonymisation.

We retained PASS, completely phased, biallelic, non-palindromic SNPs with minor-allele frequency at least 0.10 in each population and donor pool. We selected 64 coordinate-ordered common markers per 100 kb region. Sampling two whole local haplotypes independently with replacement generates each synthetic individual's genotype vector:

$$
G_i=h_{a_i}+h_{b_i},\qquad a_i,b_i\sim\operatorname{Uniform}\{1,\ldots,H\}.
$$

This preserves empirical within-region haplotype dependence and finite-panel allele frequencies. It does not generate genome-wide ancestry, pedigrees, migration or new recombination within the window. Independently simulated participants can reuse the same finite set of haplotypes. Donors are reference material, not the 3,200 simulated participants.

External reference LD selected an anchor and a distinct high-LD variant with the largest reference r² below 0.98 and at least 0.50. A second distinct variant minimised reference r² with the anchor and had r² below 0.05. These choices used no phenotypes. Study-pool LD was reported separately and was not used to revise the causal pair after selection. In the third region, the chosen high-LD pairs were perfectly correlated in both study pools; this prespecified consequence of finite reference sampling was retained as an identifiability stress test.

The eight architectures are: no effects; a cardiometabolic effect only; a mental-health effect only; one concordantly shared cause; one oppositely directed shared cause; distinct low-LD causes; distinct high-LD causes; and a shared cause plus an additional distinct mental-health cause. Active continuous effects are 0.18 residual-scale units per ALT allele and active mental-health effects are 0.35 probit-liability units, with a negative mental-health effect for opposing sharing. The last architecture intentionally violates the single-causal assumption of standard coloc. Effects were fixed before the main results were examined.

The known generating genetic covariance between the continuous predictor and the mental-health liability predictor is

$$
\operatorname{Cov}(g^\mathsf{T}\beta_C,g^\mathsf{T}\beta_M)=\beta_C^\mathsf{T}\Sigma_G\beta_M,\qquad\Sigma_G=2\operatorname{Cov}(h).
$$

The same-variant diagonal contribution is the sum of beta Cj times beta Mj times Sigma Gjj; off-diagonal terms arise from LD between different causal variants. This analytic truth is reported separately from inferred sharing. It is covariance between genetic predictors before treatment, thresholding and misclassification, rather than genetic covariance of the observed diagnoses. A nonzero total genetic covariance does not require any directly shared causal variant.

### Phenotype generation

Let A denote standard-normal age clipped to ±2.5, S a centred Bernoulli sex indicator, E a standard-normal measured exposure, and g the genotype centred using the generating haplotype pool. The cardiometabolic predictor and mental-health liability are

$$
Y^*=0.25A+0.15S+0.20E+0.08A^2+g^\mathsf{T}\beta_C+\epsilon_C,
$$

$$
L=-1.15+0.20A+0.15S+0.20E+0.05A^2+g^\mathsf{T}\beta_M+\epsilon_M,\qquad D^*=\mathbf{1}(L>0).
$$

In the ideal regime, epsilon C is Z1 and epsilon M is rho Z1 plus sqrt(1-rho²) Z2, where independent Z1 and Z2 are standard normal and rho = 0.35. The liability intercept is fixed across architectures; it is not retuned to force equal prevalence under alternatives. Dependence can produce comorbidity even in the architecture with no genetic effects.

In the realistic regime, dependence varies with age, and the continuous residual includes heavy tails and asymmetry:

$$
\rho(A)=\tanh\{\operatorname{atanh}(0.30)+0.35A\}.
$$

$$
\epsilon_C=0.70Z_1+0.45T_5/\sqrt{5/3}+0.15(Z_1^2-1).
$$

Here T5 is an independent Student t variable. The mental-health residual remains rho(A) Z1 plus sqrt(1-rho(A)²) Z2. The function rho is the Z1 coupling parameter; after adding the other continuous residual terms it is not the Pearson correlation between the observed traits.

Treatment B has probability expit(-1.8 + 0.7Y* + 0.3A), and the measured cardiometabolic outcome is Y = Y* - 0.6B + 0.25 Z3. For presentation, pressure equals 120 + 12Y in illustrative mmHg; the GWAS uses Y. Treatment thus lowers the proxy pressure by 7.2 mmHg. It depends on the pressure predictor and age, not mental diagnosis, avoiding a medication pathway that would silently turn a mental-only genetic architecture into an effect on both outcomes. The binary observation has sensitivity 0.90 and specificity 0.98 for D*. A complete observation is unavailable with probability expit(-2.7 + 0.25A + 0.20E). Availability is missing at random conditional on fitted covariates. Complete-case analysis was used. These values are plausible stress-test choices, not estimates from a clinical cohort.

Descriptive comorbidity is the co-occurrence of observed D = 1 and proxy pressure at least 140. Its odds ratio uses a 0.5 cell correction. The genetic analysis retains continuous pressure; it does not test two clinical diagnoses or the comorbidity indicator directly.

### Marginal GWAS and joint association

Each retained marker was fitted separately using linear regression for Y and logistic regression for observed D, with intercept, A, S, E and A². Logistic effects are estimated log odds ratios; generating mental effects are probit-liability coefficients and are not numerically equated to them. A batched Newton algorithm maximised the logistic likelihood, with bounded step increments and an explicit convergence requirement. If s and H denote the score and information at an iterate,

$$
\theta_{t+1}=\theta_t+H_t^{-1}s_t.
$$

$$
\widehat{\operatorname{Var}}_{HC0}(\widehat\theta)=H^{-1}\left(\sum_i s_i s_i^\mathsf{T}\right)H^{-1}.
$$

Step damping is a solver device; it does not replace convergence checks. HC0 sandwich variances supply marginal Wald p-values. Cross-products of the two regression influence functions estimate cross-trait covariance. For the two estimated marker effects b and their joint sandwich covariance V, the omnibus statistic is

$$
T=b^\mathsf{T}V^{-1}b,\qquad p_{joint}=\Pr(\chi^2_2\ge T).
$$

Its null states that both marginal coefficients are zero. Rejection establishes association with at least one phenotype, and a Bonferroni correction across 64 markers gives a region-level joint association p-value. This does not test a shared causal variant.

### Explicit conjunction of genetic associations

For a marker j, the conjunction null states that at least one phenotype has no marginal association:

$$
H_{0,j}^{conj}:\{\beta_{C,j}=0\}\ \operatorname{or}\ \{\beta_{M,j}=0\},\qquad
p_{IUT,j}=\max(p_{C,j},p_{M,j}),\quad
p_{region}=\min\{1,64\min_j p_{IUT,j}\}.
$$

This intersection-union test requires evidence for both associations and does not require independence between their p-values. The Bonferroni correction controls marker search under the corresponding region-wide conjunction null when each marginal p-value is valid. Distinct causal variants in LD can make a tag genuinely associated with both phenotypes: a conjunction rejection in that architecture is therefore not a false marginal-association finding, but would be a false attribution if interpreted as one shared cause. We report this distinction rather than treating conjunction as proof of pleiotropy.

### Shared causal variant inference

We used the standard single-causal-variant coloc approximate Bayes factor framework [4,5]. The primary continuous GWAS uses the first half of complete participants and the binary GWAS uses the second half. The partition depends only on participant index. This removes cross-trait sampling-error dependence between the two primary summary datasets, while reducing effective sample size. It estimates sharing from separate subsamples of the same generating population, rather than an individual-level copula likelihood. A whole-sample overlapping analysis is a diagnostic sensitivity with both a sample-size advantage and possible correlated errors; their effects cannot be separated by that contrast.

For estimated marker effect b, variance V and prior effect variance W, the Wakefield log Bayes factor is

$$
\log BF=\tfrac12\{\log(1-r)+rz^2\},\qquad r=\frac{W}{V+W},\quad z=\frac{b}{\sqrt V}.
$$

We used model-based GWAS variances for standard ABFs, a continuous prior SD of 0.15 in the fitted outcome units and a logistic prior SD of 0.20 log-odds units. These prior scales differ from generating effect sizes. With per-marker BF vectors BC and BM, unnormalised hypothesis weights are

$$
w_0=1,\quad w_1=p_1\sum_j BC_j,\quad w_2=p_2\sum_j BM_j.
$$

$$
w_3=p_1p_2\sum_{j\ne k}BC_jBM_k,\quad w_4=p_{12}\sum_j BC_jBM_j.
$$

Normalising the five weights gives posterior probabilities for neither trait associated, cardiometabolic association only, mental association only, distinct causes, and one shared cause. The primary priors are p1 = p2 = 10^-4 and p12 = 10^-5, with p12 sensitivity at 10^-6 and 10^-4. PP4 > 0.8 is the conventional decision rule. These weights and posteriors were checked against the established R coloc package, version 5.2.3. Colocalisation is conditional on the candidate set, effect priors and single-causal model; it cannot establish a biological mechanism or the causal direction between diseases.

### Composite null calibration and uncertainty

The no-shared-cause calibration set comprises no effects, cardiometabolic only, mental health only, distinct low-LD causes and distinct high-LD causes. For each population, region and phenotype regime, we independently calibrated the split and overlapping PP4 scores. For 200 calibration replicates per architecture, the rank allowance is k = floor(0.05 × 201) = 10. For an upper-tail score Q, the cutoff is the largest 191st ascending score across the five no-sharing architectures:

$$
c=\max_{s\in\mathcal N}Q_{s,(191)},\qquad\text{declare sharing only if }Q>c.
$$

For each individual calibration architecture, exchangeability gives average rejection probability at most 10/201 with the strict order-statistic rule in the absence of ties; taking the maximum is conservative for all included architectures. This guarantee is averaged over calibration banks, not conditional on one realised bank, and applies only to the generating mechanisms included in calibration. PP4 is an evidence score in this frequentist calibration; the resulting threshold does not create a universal posterior probability of causation. A separately calibrated tagging test uses the smallest 10th ascending region p-value across the no-effect and trait-specific banks; distinct-cause scenarios are not its marginal-association null.

Per-cell rates include all 200 attempted evaluation datasets, with Wilson intervals. Aggregate rates weight the six population-region cells equally within each phenotype regime and architecture. Four hundred paired stratified bootstrap draws independently resample calibration and evaluation replicates, recompute all cutoffs and preserve method pairing. Reported aggregate percentile intervals include both sources of simulation uncertainty. An ordinary bootstrap degenerates if every evaluation decision is zero or one. For boundary entries, Table 6 instead displays a finite-evaluation binomial bound conditional on the realised cutoffs; this bound does not include uncertainty from unseen calibration-bank realisations. Raw bootstrap quantiles remain published and a zero-width bootstrap interval must not be interpreted as certainty about a population rate. Comparisons are descriptive and do not provide multiplicity-adjusted superiority claims. Any failed main fit is retained as a non-discovery; a failed calibration fit blocks cutoff generation.

### Copula diagnostics and validation

The original Gaussian/probit code is retained unchanged. The first 20 evaluation replicates in every stratum and architecture form a fixed 1,920-dataset diagnostic panel. We fitted varying dependence with intercept and age and evaluated the sandwich score at the genotype-selected anchor. Both marginal models use the same five covariates as the new GWAS layer. The Gaussian/probit likelihood is correctly motivated by the ideal data and is deliberately stressed by heavy tails, treatment and observed diagnostic error in the realistic regime. These nominal anchor tests have a different search burden from region-level tests; they are not used for a power ranking against the full-panel methods or for claims of shared causation.

Independent controls compare the batched logistic solver against a separate optimiser, validate OLS and HC0 calculations, flip allele coding, replay seeded datasets and demonstrate exact non-identifiability when causal genotype columns are equal. R lm and glm independently check full-panel estimates and variances in prespecified fixtures; R coloc checks posterior enumeration. The source identity, raw chunk hashes and original-result checksums are audited before publication.

## Results

### Reference LD and identifiability

The retained panels contain 192 marker positions and population-specific empirical LD. Distinct high-LD causal pairs have study r² from @RMIN@ to 1.000. The third-region pair is perfectly correlated in both study pools. Under the same random stream, swapping the mental effect from the shared anchor to the distinct perfectly correlated variant leaves all observed genotypes and phenotypes unchanged. This exact counterexample shows why an association procedure cannot always resolve sharing even when both trait associations are strong.

@CAUSAL_TABLE@

Table 1. Public reference and study donor pools. The 64 markers per region are shared across populations; causal index choices are population-specific and were selected using external reference LD.

![African haplotype LD across the three regions](v3/figures/african_ld.png)

Figure 1. Study-pool squared haplotype correlations for YRI and LWK. Population labels describe the public reference samples, not simulated diagnoses or inferred ancestry in a real cohort. Diagonal blocks and perfectly correlated tags limit causal resolution.

The generating covariance decomposition further separates sharing from LD. Distinct high-LD causes can contribute covariance between the two genetic predictors while having zero same-variant contribution. Opposing shared effects instead contribute negative covariance. These are known properties of the generating mechanisms, not estimated discoveries.

@TRUTH_TABLE@

Table 2. Analytic generating genetic-predictor covariance averaged equally over six population-region cells. Same-variant sharing and total genetic covariance are different quantities.

### Phenotype characteristics

The realistic regime reduces complete observations to approximately @NCOMPLETE@ of 3,200 generated participants. Observed depression prevalence is @PREV@ and the mean treated fraction is @TREAT@ across architectures and reference strata. Pressure remains in an interpretable illustrative range while treatment, diagnostic error and missingness change its relationship with the observed depression indicator. Comorbidity also occurs in the no-genetic-effect architecture, consistent with the explicitly simulated residual dependence. Its presence alone is therefore insufficient evidence of a shared genetic cause.

@PHENOTYPE_TABLE@

Table 3. Mean descriptive phenotype summaries across equally weighted design cells. These are synthetic distributions, not estimated clinical prevalences.

![Synthetic phenotype distributions and observation quality](v3/figures/phenotype_checks.png)

Figure 2. Distributions of replicate-averaged cell summaries. Boxes describe variation across the 48 architecture-by-population-by-region cells in each regime, rather than confidence intervals for a patient population.

### Association and sharing under ideal phenotypes

Joint association detects many trait-specific effects because its alternative includes an effect on either phenotype. Requiring both marginal associations changes the target, but distinct variants in LD can still produce a conjunction at the same marker. Under concordant shared causes, conventional colocalisation declares sharing in @IDEAL_SHARED@ of datasets. In the distinct high-LD architecture it declares sharing in @IDEAL_DISTINCT@. Composite calibration addresses that ambiguity by using distinct-cause mechanisms as explicit no-sharing controls, with the corresponding power cost shown below.

@IDEAL_TABLE@

Table 4. Ideal-regime discovery proportions. Each entry uses 1,200 evaluation datasets. Joint association and tagging conjunction have different inferential targets from colocalisation. Their rejection rates must not be interpreted as competing estimates of true pleiotropy.

![Discovery proportions under ideal phenotypes](v3/figures/discovery_ideal.png)

Figure 3. Ideal-regime rates with calibration-aware bootstrap intervals. Higher rejection under a trait-specific architecture is expected for an omnibus association test; it is not evidence for a shared cause.

### Sharing under realistic phenotypes

Under realistic concordant shared effects, PP4 > 0.8 declares sharing in @SHARED_RAW@ of datasets, compared with @SHARED_CAL@ after composite calibration. For distinct high-LD causes, these rates are @DISTINCT_RAW@ and @DISTINCT_CAL@. The calibrated rule therefore needs to be interpreted together with its reduced sensitivity and with the irreducible ambiguity in the third region. The shared-plus-distinct architecture is a stress test of the single-causal model, not an additional setting in which that assumption is satisfied.

@REAL_TABLE@

Table 5. Realistic-regime discovery proportions, with the same 1,200-dataset denominators as Table 4.

![Discovery proportions under realistic phenotypes](v3/figures/discovery_realistic.png)

Figure 4. Realistic-regime rates with paired bootstrap uncertainty. The figure separates association discovery from shared-cause declarations and retains all simulated alternatives.

@CI_TABLE@

Table 6. Aggregate colocalisation rates with percentile intervals from 400 stratified bootstrap draws for interior rates. Boundary rates show a 95% finite-evaluation bound conditional on the cutoffs, rather than a degenerate bootstrap interval. Only the bootstrap intervals include calibration-bank resampling; none are clinical effect intervals.

### Calibration under distinct causal variants

The strongest difficulty is the distinct high-LD architecture. In contrast to a null with no associations, it contains genuine signals for both traits and can give high PP4 to a shared explanation. Cell-level results show where external-reference and study LD differ and where the shared and distinct configurations cannot be distinguished. The calibrated cutoff is allowed to be stringent in these cells; reducing it to recover apparent power would invalidate the stated no-sharing control.

![Cell level shared versus distinct variant resolution](v3/figures/identifiability.png)

Figure 5. Shared-cause declaration rates under a genuinely shared concordant variant and under distinct high-LD variants. Region 3 has perfectly correlated causal candidates in both study populations. Differences there arise from finite calibration and evaluation samples rather than identifiable causal information.

![Held out no sharing distributions](v3/figures/composite_null.png)

Figure 6. Held-out PP4 distributions under five mechanisms with no shared direct causal variant. The vertical line marks the conventional 0.8 rule. Strong PP4 in the high-LD distinct-cause architecture illustrates the difference between a model-based posterior and error control against realistic competing explanations.

### Prior choice and overlapping samples

Prior sensitivity changes the strength of sharing declarations, as expected from the explicit p12 multiplier in the H4 weight. The overlapping whole-sample diagnostic uses more information per phenotype than either split-sample arm and also shares estimation errors between traits. Its comparison with the primary split analysis therefore cannot isolate overlap bias. We retain it to show the consequence of a commonly used analysis choice, without presenting it as an independent confirmation of sharing.

![Prior and overlap sensitivity](v3/figures/prior_overlap_sensitivity.png)

Figure 7. PP4 > 0.8 rates under three sharing priors and the overlapping whole-sample diagnostic. The latter changes both effective sample size and sampling-error dependence.

### Copula diagnostics and numerical validation

The preserved varying-dependence copula test evaluates association at the fixed anchor. Its rejection under a cardiometabolic-only effect or a distinct high-LD cause is compatible with its omnibus target. The smaller diagnostic panel does not establish a ranking of region-level power and does not turn copula dependence into genetic correlation.

@COP_TABLE@

Table 7. Nominal anchor tests in the realistic-regime diagnostic panel. Each entry pools 20 replicates across six population-region cells. Both procedures test association with at least one phenotype; neither is a shared-causal-variant test.

The main audit recorded @MAIN_FAIL@ failed datasets out of 31,200 attempts and @COP_FAIL@ failures among 1,920 planned copula diagnostics. All original scientific files were preserved byte-for-byte. Twelve independently selected evaluation datasets were regenerated from their seeds and matched their saved statistics to numerical tolerance. R validation independently checked full-panel regression calculations and the established coloc posterior implementation. Detailed errors, package versions, raw results, cutoffs and audit records are available in the repository.

## Discussion

The experiment addresses a more specific question than whether modelling two outcomes can improve an association test. The relevant question for comorbidity genetics is whether the evidence distinguishes one causal variant affecting both traits from alternative explanations. Trait-specific effects, correlated but distinct causal variants and shared residual factors all need to be represented in that assessment.

Conjunction testing supplies explicit evidence that a marker is associated with both measured traits. Its interpretation stops at marginal association: LD can transmit two different causal signals to the same tag. Colocalisation adds a model of competing causal configurations, but a single-causal posterior depends on the available variants, prior weights and the possibility of resolving their genetic patterns. Perfect LD is an identifiability limit rather than a numerical problem. More complex likelihoods cannot recover distinctions absent from the observed data.

Composite-null calibration makes that limitation visible. A high fixed posterior cutoff can appear persuasive when judged only against no genetic effects, yet perform poorly against genuinely associated distinct causes. Calibrating against those causes can improve specificity for the tested settings, at the cost of sensitivity. The present thresholds are conditional on this simulation design and should not be transferred to a real cohort or a genome-wide scan.

African genomic diversity matters here through empirical allele frequencies and LD differences between YRI and LWK and between their small reference and study pools. These observations do not justify a generic claim that African ancestry guarantees greater fine-mapping resolution. The public populations, finite donor pools and three selected windows cover a narrow part of African variation. An applied study would need cohort-specific ancestry and relatedness models, appropriate reference panels, comprehensive variant coverage, phenotype harmonisation and independent replication.

The realistic phenotype regime tests several practical departures from a simple Gaussian/probit model, but it is not a fully calibrated representation of depression, hypertension or multimorbidity in an African cohort. Treatment, diagnostic accuracy and missingness were chosen for controlled sensitivity analysis. Their consequences should motivate measurement checks in a real study rather than clinical conclusions. The binary-liability and continuous-proxy setup also leaves two-binary-disease models, repeated measurements and survival outcomes for future work.

The shared-plus-distinct architecture identifies a further model boundary. Standard ABF coloc assumes at most one cause per trait in the candidate region. Multi-signal fine-mapping and methods such as coloc with SuSiE are appropriate next comparisons [6]; they were not executed in this version. Such methods would still face indistinguishability at perfectly correlated variants. Similarly, colocalisation cannot distinguish direct horizontal sharing from a shared signal caused through one phenotype influencing another without additional assumptions and evidence.

No real shared genetic cause of cardiometabolic and mental illness has been established here. Instead, the experiment evaluates how often specified procedures make correct or incorrect sharing declarations when the generating causes are known. This constitutes methodological evidence that could inform a prospective analysis plan for suitable real data.

## Conclusion

Public African haplotypes and realistic phenotype perturbations provide a reproducible test of shared-variant inference under known simulated mechanisms. Association with two traits is insufficient to establish a shared causal variant, and high LD can make shared and distinct causal explanations observationally equivalent. The resulting calibration and uncertainty analyses define the conditions under which sharing claims are supported and where causal resolution remains unavailable.

## Reproducibility

The frozen design, source identity, public-reference provenance, raw simulation chunks, all evaluation tables and numerical audits are versioned. The original study and the 92,000-dataset oracle experiment retain their own results and Word reports. The new report shares this manuscript's text, equations, figures and tables and names Vusumuzi Mabasa as author. It is generated with Pandoc and python-docx; it is not a Quarto execution.

Use Python 3.12 with the pinned root requirements, R 4.5.1, jsonlite 2.0.0 and coloc 5.2.3. Set BLAS thread counts to one when using process workers. Always run R with --vanilla. See [the complete workflow](workflow.txt) for Toshiba and Linux commands, reference retrieval, checkpoint handling, report generation and verification. The small smoke run is separate from the main results. The main analysis can take substantial time; existing chunks are reused only after source and output hashes pass validation. Public reference extraction is reproducible from the original checksum-verified files, without an absolute Toshiba path in the scientific code.

```bash
python -m pip install -r requirements.txt
python v3/src/freeze.py
python v3/src/test_model.py
python v3/src/run.py --smoke --workers 2
python v3/src/run.py --workers 2
python v3/src/analyse.py
python v3/src/audit.py --summaries
python v3/src/export_r.py
Rscript --vanilla v3/src/crosscheck.R v3
python v3/src/manuscript.py
python v3/src/build_word.py --pandoc /path/to/pandoc
```

Primary outputs include [per-cell rates](v3/results/rates.csv), [aggregate rates](v3/results/aggregate_rates.csv), [paired comparisons](v3/results/paired_contrasts.csv), [frozen cutoffs](v3/results/thresholds.json), [phenotype checks](v3/results/phenotypes.csv), [copula diagnostics](v3/results/copula_diagnostic.csv), [reference QC](v3/reference/qc.csv) and [the complete audit](v3/results/audit.json). Uncertainty draws and independent R validations are also retained. Main simulations and bootstrap resampling are implemented in Python; independent regressions and established-method validation are performed in R. There is no Genomic SEM component.

## References

1. Wits University. Uncovering links between depression and hypertension in African populations. 28 August 2025. [Programme announcement](https://www.wits.ac.za/news/latest-news/research-news/2025/2025-08/uncovering-links-between-depression-and-hypertension-in-african-populations.html).
2. 1000 Genomes Project Consortium. A global reference for human genetic variation. Nature 526, 68–74 (2015). [doi:10.1038/nature15393](https://doi.org/10.1038/nature15393).
3. International Genome Sample Resource. [Public data and reuse information](https://www.internationalgenome.org/data/) and [terms of use](https://www.internationalgenome.org/IGSR_disclaimer/). Exact Phase 3 source files and SHA256 values are recorded in [the provenance file](v3/reference/provenance.json).
4. Giambartolomei C et al. Bayesian test for colocalisation between pairs of genetic association studies using summary statistics. PLoS Genetics 10, e1004383 (2014). [doi:10.1371/journal.pgen.1004383](https://doi.org/10.1371/journal.pgen.1004383).
5. Wallace C et al. [coloc single-causal-variant methods](https://chr1swallace.github.io/coloc/articles/a03_enumeration.html) and [prior sensitivity](https://chr1swallace.github.io/coloc/articles/a04_sensitivity.html). Executed established comparator: R coloc 5.2.3.
6. Wallace C. Eliciting priors and relaxing the single causal variant assumption in colocalisation analyses. PLoS Genetics 16, e1008720 (2020). [doi:10.1371/journal.pgen.1008720](https://doi.org/10.1371/journal.pgen.1008720).
'''
    replacements={
        '@SHARED_RAW@':rate('realistic','shared_same',METHODS[3]),'@SHARED_CAL@':rate('realistic','shared_same',METHODS[4]),'@DISTINCT_RAW@':rate('realistic','distinct_high',METHODS[3]),'@DISTINCT_CAL@':rate('realistic','distinct_high',METHODS[4]),'@IDEAL_SHARED@':rate('ideal','shared_same',METHODS[3]),'@IDEAL_DISTINCT@':rate('ideal','distinct_high',METHODS[3]),
        '@CAUSAL_TABLE@':causal_table,'@TRUTH_TABLE@':truth_table,'@PHENOTYPE_TABLE@':phenotype_table,'@IDEAL_TABLE@':result_table('ideal'),'@REAL_TABLE@':result_table('realistic'),'@CI_TABLE@':ci_table,'@COP_TABLE@':cop_table,
        '@RMIN@':f"{min(float(r['study_r2_high']) for r in pairs):.3f}",'@NCOMPLETE@':f"{meanphen('realistic','n_complete'):.0f}",'@PREV@':pct(meanphen('realistic','prevalence')),'@TREAT@':pct(meanphen('realistic','treated_fraction')),'@MAIN_FAIL@':str(audit['main_failures']),'@COP_FAIL@':str(audit['copula_failures'])}
    for a,b in replacements.items():text=text.replace(a,b)
    import re
    assert not re.search(r'@[A-Z_]+@',text)
    (ROOT/'README.md').write_text(text,encoding='utf-8',newline='\n')
    dump(V/'results/manuscript_sources.json',dict(source='v3/src/manuscript.py',rates_sha256=sha(V/'results/aggregate_rates.csv'),audit_sha256=sha(V/'results/audit.json'),readme_sha256=sha(ROOT/'README.md')))
    print('Research manuscript generated')
if __name__=='__main__':main()
