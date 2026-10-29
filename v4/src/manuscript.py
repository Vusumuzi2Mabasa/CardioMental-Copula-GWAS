"""Integrate all three levels of simulated shared genetics in one research manuscript."""
from pathlib import Path
import sys,csv,json,subprocess,re
ROOT=Path(__file__).resolve().parents[2];sys.path.append(str(ROOT/'v3/src'));V=ROOT/'v4'
from prepare import dump,sha
def pct(x):return f'{100*float(x):.1f}%'
def mdtable(headers,rows):return '\n'.join(['| '+' | '.join(headers)+' |','| '+' | '.join(['---']*len(headers))+' |']+['| '+' | '.join(map(str,r))+' |' for r in rows])

def read(path):
    with path.open() as f:return list(csv.DictReader(f))
def main():
    # Importing v3 helpers by explicit path avoids the identical module basename.
    import importlib.util
    spec=importlib.util.spec_from_file_location('v3_manuscript',ROOT/'v3/src/manuscript.py');base_module=importlib.util.module_from_spec(spec);spec.loader.exec_module(base_module)
    base_module.main();text=(ROOT/'README.md').read_text(encoding='utf-8');base_sha=sha(ROOT/'README.md')
    audit=json.loads((V/'results/audit.json').read_text());assert audit['passed']
    fm=read(V/'results/finemap_aggregate.csv');comp=read(V/'results/finemap_matched_comparison.csv');pa=read(V/'results/polygenic_validation/polygenic_aggregate.csv');pc=read(V/'results/polygenic_validation/polygenic_rates.csv')
    def poly(n,g):return next(r for r in pa if int(r['n'])==n and float(r['genetic_truth'])==g and r['mode']=='observed_binary')
    def fine(s):return next(r for r in comp if r['regime']=='realistic' and r['LD']=='matched_reference' and r['scenario']==s)
    def fixed(x):return f'{float(x):.3f}'
    from analyse import LABELS
    fine_table=mdtable(['Realistic architecture','Single causal rule','SuSiE pair rule','Single causal calibrated','SuSiE calibrated'],[[LABELS[s]]+[pct(fine(s)[k]) for k in ['abf_conventional','susie_conventional','abf_calibrated','susie_calibrated']] for s in ['null','cardio_only','mental_only','shared_same','shared_opposing','distinct_low','distinct_high','shared_plus_distinct']])
    recovery_table=mdtable(['LD source','Cardio causal recovery','Mental causal recovery','Cardio largest set','Mental largest set'],[[r['LD'].replace('_',' '),pct(r['c_credible_coverage']),pct(r['b_credible_coverage']),f"{float(r['c_mean_max_set_size']):.1f}",f"{float(r['b_mean_max_set_size']):.1f}"] for r in fm if r['regime']=='realistic' and r['scenario']=='shared_plus_distinct'])
    poly_table=mdtable(['n','True genetic r','Mean estimate','RMSE','Interval coverage','Calibrated rejection'],[[n,f'{g:.1f}',fixed(poly(n,g)['mean_rg']),fixed(poly(n,g)['rmse_rg']),pct(poly(n,g)['coverage_attempted']),pct(poly(n,g)['calibrated_rate'])] for n in [600,1200] for g in [-.5,0.,.5]])
    nullrows=[]
    for n in [600,1200]:
        for e in [0.,.4]:
            rows=[r for r in pc if int(r['n'])==n and float(r['genetic_truth'])==0 and float(r['residual_truth'])==e and r['mode']=='observed_binary']
            avg=lambda k:sum(float(r[k]) for r in rows)/len(rows)
            nullrows.append([n,f'{e:.1f}',fixed(avg('mean_phenotypic_correlation')),fixed(avg('mean_rg')),pct(avg('calibrated_rate'))])
    null_table=mdtable(['n','Residual liability r','Observed phenotype r','Estimated genetic r','Calibrated null rate'],nullrows)
    title='Shared genetic associations causal variants and polygenic architecture in simulated cardiometabolic and mental health traits'
    text=text.replace(text.splitlines()[0],'# '+title,1)
    text=text.replace('report/CardioMental_Copula_GWAS_Report_v3.docx','report/CardioMental_Shared_Genetics_Report_v4.docx')
    text=text.replace('[Complete results](v3/results/aggregate_rates.csv)','[Local results](v3/results/aggregate_rates.csv) · [Fine mapping](v4/results/finemap_matched_comparison.csv) · [Polygenic validation](v4/results/polygenic_validation/polygenic_aggregate.csv)')
    abstract=f'''## Abstract

Shared genetic association, one shared causal variant and genome-scale genetic correlation are distinct inferential targets. We evaluated all three in synthetic cardiometabolic and mental-health data using public African haplotypes. The local experiment generated 31,200 datasets using Yoruba and Luhya LD, trait-specific and shared effects, distinct variants in LD and realistic phenotype perturbations. A fixed 3,120-dataset subset supported 6,240 ancestry-aware multi-signal fine-mapping and colocalisation analyses. A separate polygenic model used 22 independent synthetic LD blocks and bivariate relationship-kernel REML. Its initial 4,000-cohort experiment was retained, followed by 8,000 new calibration cohorts and 4,800 fresh evaluation cohorts. Across both experiments, 48,000 synthetic cohorts were generated; outcome modes and reused fine-mapping subsets are not additional independent cohorts. In fresh binary-outcome validation at n=1,200, mean estimated genetic correlation was {fixed(poly(1200,.5)['mean_rg'])} for truth 0.5 and {fixed(poly(1200,-.5)['mean_rg'])} for truth -0.5; calibrated rejection under zero genetic covariance was {pct(poly(1200,0)['calibrated_rate'])}. Local LD included configurations in which shared and distinct causal variants were observationally indistinguishable. The results provide methodological evidence for separating the three targets and for testing genetic correlation while allowing residual dependence. All clinical phenotypes and genetic effects are simulated; no shared genetic cause of real cardiometabolic and mental disorders is established.

'''
    text=re.sub(r'## Abstract\n.*?(?=## Introduction)',lambda _:abstract,text,flags=re.S)
    methods=r'''### Multi signal fine mapping for Level 2

The additional Level 2 protocol was frozen in commit 85b0bc4, after the local experiment had begun and before its fine-mapping analyses. We reused the first 20 calibration and evaluation replicates in every local stratum and architecture. The subset contains 1,200 calibration and 1,920 evaluation datasets; two LD conditions produce 6,240 analysis records. These datasets are included in the local experiment's denominator and are not new independent replications. Conventional ABF and SuSiE comparisons use the same subset and the same 20-replicate calibration banks in the secondary matched-subset analysis.

We ran established R SuSiE and coloc.susie, with up to five single-effect components and 95% credible sets [6,7]. SuSiE represents the effect vector as a sum of sparse components, each assigning one effect to a candidate variant:

$$
\beta=\sum_{\ell=1}^{5}b_\ell,\qquad z\ \mathrel{\dot\sim}\ N(\sqrt{N}R\beta,R).
$$

The RSS normal approximation conditions on aligned signed LD and the standardised effect vector. In this expression beta denotes standardised effects, not the unscaled generating probit coefficients. It is an approximation for the binary logistic GWAS. Residual variance estimation was disabled, and retained fits had to converge within 1,000 iterations. Prior effect variance was estimated by the package. Credible sets required the package's default purity checks; absence of a retained signal is a valid no-discovery outcome. Variant posterior inclusion probabilities and all retained credible sets are published, including known simulated causal indices.

The primary external LD matrix uses the 40 ancestry-matched reference donors, paired as their observed diploid genotypes. The comparison uses LD from the actual simulated study halves, residualised for the continuous GWAS covariates and weighted using the binary null-model information for the binary GWAS. Both LD matrices use the fixed regularisation

$$
R_{used}=0.95R+0.05I.
$$

This provides a stable positive-definite approximation and does not eliminate reference-panel sampling error. All columns retain ALT-allele orientation and the same variant order. Allele/LD sign-flip controls check this alignment. The external reference is ancestry-matched rather than interchangeable with in-sample LD; both arms are reported.

coloc.susie evaluates pairs of retained signals. The maximum pairwise PP4 is the sharing score. A maximum above 0.8 is a conventional declaration, not the overall posterior probability that any pair shares a cause. The calibrated score uses the largest observed maximum across five no-sharing banks per stratum. With 20 replicates per bank, the strict rank rule gives an averaged bound no greater than 1/21 for each included mechanism under exchangeability. This is not a conditional guarantee for one realised bank. The small banks imply coarse cutoffs and considerable uncertainty. Numerical failures are retained as non-discoveries and a calibration failure disables that calibrated stratum; they are never excluded to improve the reported rate. Exact causal recovery is the proportion of generating causal variants contained in the union of retained credible sets, including datasets without retained sets in its denominator. It is a recovery measure, not a claim that every empirical credible set has nominal frequentist coverage.

### Polygenic architecture for Level 3

Level 3 uses a separate genome-scale simulation, rather than extrapolating the three-region covariance decomposition. We extracted 32 common SNPs from each of 22 non-overlapping chromosome 22 reference windows spanning 20.0 to 41.1 Mb. Independent whole-haplotype sampling across windows defines 22 synthetic genomic blocks and 704 markers. All source positions remain identified as chromosome 22. This is a synthetic polygenic genome informed by African local LD, not a real chromosome 1–22 African reference panel. The small marker backbone and independent-block assumption limit extrapolation to actual genome-wide analyses.

Within YRI and LWK, we generated genotype kernels at n=600 and n=1,200. Marker dosages are standardised using generating-pool allele frequencies. Rows are normalised so the relationship kernel has unit diagonal. If Z is this genotype representation with m=704 columns,

$$
K=ZZ^\mathsf{T}/m,\qquad
\operatorname{Cov}\!\begin{pmatrix}g_C\\g_M\end{pmatrix}=\Sigma_g\otimes K.
$$

The two traits follow Y = X B + g + e, with independent genetic and residual random effects. Each genetic variance is 0.4, each residual variance is 0.6, genetic correlation is -0.5, 0 or 0.5, and residual correlation is 0 or 0.4. Genotypes and the kernel are fixed within population and sample size; random effects and errors vary independently across replicates. The target is the random-effects covariance parameter conditional on that kernel, rather than the realised correlation of two finite vectors of SNP effects in one replicate. Genetic correlation is

$$
r_g=\frac{\sigma_{g,CM}}{\sqrt{\sigma^2_{g,C}\sigma^2_{g,M}}}.
$$

An intercept and an independently generated nuisance covariate are projected out. Their generating coefficients are zero. For an orthonormal residual basis Q and eigendecomposition of the projected kernel, the restricted likelihood decomposes into independent bivariate modes:

$$
Q^\mathsf{T}KQ=U\Lambda U^\mathsf{T},\qquad
v_\ell=(U^\mathsf{T}Q^\mathsf{T}Y)_\ell,\quad
V_\ell=\lambda_\ell\Sigma_g+\Sigma_e.
$$

$$
\ell_R=-\frac{1}{2}\sum_\ell(\log|V_\ell|+v_\ell^\mathsf{T}V_\ell^{-1}v_\ell)+C.
$$

The implementation estimates both full covariance matrices using log standard deviations and Fisher-transformed correlations, with analytic gradients. The likelihood-ratio test compares unrestricted genetic covariance with zero genetic covariance while allowing residual covariance to remain free. Information-based intervals use the numerical Hessian of the analytic gradient. Boundary estimates, invalid information matrices and failed fits remain recorded. This is an independently implemented bivariate REML model of the type used in genomic relationship analyses [8,9]; the GCTA software itself was not run. Independent R optimisation verifies the restricted likelihood and estimates.

The exact benchmark observes continuous mental-health liability. The applied working analysis observes a binary depression proxy with population prevalence 0.15, sensitivity 0.90 and specificity 0.98. The same cohort supplies both modes, so they are paired and not independent replications. The Gaussian REML likelihood is exact for the continuous benchmark and a working approximation for the continuous/binary pair. This genome-scale experiment isolates covariance estimation and does not reuse the local experiment's treatment or missingness mechanisms.

With threshold t = Phi inverse of 0.85 and diagnostic attenuation q = 0.88, an approximate liability transformation uses kappa = q phi(t):

$$
\sigma_{g,CM}^{observed}=\kappa\sigma_{g,CM}^{liability},\qquad
\sigma_{g,M}^{2,observed}\ \mathrel{\dot\approx}\ \kappa^2\sigma_{g,M}^{2,liability}.
$$

The cross-covariance relation follows from the joint Gaussian liability model; the binary variance relation is a first-order approximation in between-person genetic relatedness. Genetic correlation approximately cancels the scale factors. We report transformed covariance separately and evaluate bias, RMSE, interval coverage and error rates rather than assuming this approximation is exact. There is no case-control ascertainment in this experiment.

### Fresh calibration and validation of Level 3

The initial protocol generated 1,600 calibration and 2,400 evaluation cohorts, retaining both outcome modes. Initial 200-replicate calibration banks showed visible cell-level uncertainty, including an observed 10% null rejection rate in one 100-cohort binary-outcome cell. Rather than selectively adjusting that cell, we froze a uniform extension in commit dd84837 before generating its new datasets. All generating parameters and estimators remained unchanged. The extension used 1,000 independent calibration cohorts in every population, sample-size and residual-correlation stratum, and 200 fresh evaluation cohorts for each genetic-correlation alternative. It adds 8,000 calibration and 4,800 evaluation cohorts. The initial experiment is retained and linked, and the fresh validation supplies the primary Level 3 results below.

Calibration and evaluation use disjoint replicate-index offsets of 1,000,000 and 2,000,000 in the frozen seed namespace. For 1,000 null statistics the strict cutoff is the 951st ascending likelihood-ratio statistic, giving average rejection probability at most 50/1001 under exchangeability and no numerical failures. The cutoff is computed separately for each outcome mode and stratum before its fresh evaluation begins. It is not selected using the new evaluation outcomes. As in the local experiment, this is averaged calibration-bank control; it does not imply that every realised cell rate must be below 5%.

Each aggregate Level 3 entry contains 800 fresh cohorts across the two populations and two residual-correlation conditions. Four hundred paired bootstrap draws resample calibration and evaluation cohorts and recompute thresholds. Outcome modes remain paired. All cell-level error rates, boundary estimates and interval failures are retained. The exact likelihood, analytic gradients, sign-flip invariance, seeded replays and independent R fits are checked before interpretation. No Genomic SEM model is used.

'''
    text=text.replace('## Results\n',methods+'## Results\n',1)
    results=f'''### Multi signal colocalisation and causal resolution

The matched-subset comparison separates conventional evidence thresholds from composite-null calibrated declarations. Under realistic concordant sharing, ancestry-matched SuSiE declares a shared signal in {pct(fine('shared_same')['susie_conventional'])} of the 120 evaluation datasets, compared with {pct(fine('shared_same')['abf_conventional'])} for single-causal ABFs. Under distinct high-LD causes, those rates are {pct(fine('distinct_high')['susie_conventional'])} and {pct(fine('distinct_high')['abf_conventional'])}. These rates must be interpreted with the small subset size, finite reference sampling and irreducible LD ambiguity. A multi-signal method can represent additional causes but cannot identify which perfectly correlated allele generated a signal.

{fine_table}

Table 8. Matched 20-replicate-per-cell comparison under realistic phenotypes and ancestry-matched LD. Each entry pools 120 evaluation datasets. Both calibrated methods use the same 20-replicate no-sharing banks in this secondary comparison. A SuSiE declaration refers to any retained signal pair, not a global posterior probability.

![Multi signal and single causal colocalisation](v4/figures/fine_mapping_comparison.png)

Figure 8. Matched-subset sharing declarations under ideal and realistic phenotypes. The small subset and coarse calibration support a methodological comparison rather than a claim of universal superiority. Per-cell rates and descriptive conditional intervals are published.

Credible sets describe variant resolution separately from sharing declarations. In the shared-plus-distinct architecture, recovery requires retaining the shared cardiometabolic cause and both mental-health causes. The quantities below include datasets with no retained credible set. The largest-set summaries assign size zero when no set is retained; they are descriptive and do not imply a precisely identified causal allele.

{recovery_table}

Table 9. Exact causal-variant recovery and mean largest retained credible-set size in the realistic shared-plus-distinct architecture. Each LD arm uses the same 120 datasets. Recovery is evaluated against known simulated causal variants.

![Variant posterior inclusion probabilities](v4/figures/fine_mapping_pip.png)

Figure 9. Variant inclusion probabilities for the prespecified first YRI replicate in region 1 under realistic shared-plus-distinct effects. Gold stars identify generating causal variants; they are simulation truth, not external annotations. The two LD sources use the same summary effects and allele orientation. This fixed example was not selected for a favourable result.

### Genome scale genetic correlation in fresh validation

The fresh observed-binary analysis estimates both positive and negative genetic correlations. At n=1,200, mean estimates are {fixed(poly(1200,.5)['mean_rg'])} for truth 0.5 and {fixed(poly(1200,-.5)['mean_rg'])} for truth -0.5. For zero genetic covariance the mean is {fixed(poly(1200,0)['mean_rg'])}. Information-based coverage and calibrated rejection are assessed against the known generating parameters, with all attempted cohorts retained. The exact liability benchmark remains available to distinguish estimation behaviour from information lost through thresholding and diagnostic error.

{poly_table}

Table 10. Fresh validation of the observed binary working model. Each entry contains 800 cohorts. Interval coverage counts an unavailable interval as non-coverage; parameter summaries use converged fits, with failure counts separately reported. At zero genetic correlation the rejection rate is a false positive rate; elsewhere it is power for the specified effect size.

![Genetic correlation estimation across repeated simulations](v4/figures/genetic_correlation.png)

Figure 10. Distribution of genetic-correlation estimates in the fresh cohorts. Boxes show medians and interquartile ranges; whiskers extend to 1.5 interquartile ranges, with outlying points omitted from the display and retained in all calculations. The exact continuous-liability benchmark and observed binary working analysis use the same simulated cohorts. Black crosses mark generating correlations. A zero correlation parameter does not mean that no SNP affects both traits: infinitesimal effects can contribute opposing cross-trait products that cancel in aggregate.

### Separating genetic from residual dependence

The zero-genetic-correlation scenarios allow positive residual correlation, which induces observed phenotypic association. Estimating residual covariance freely avoids assigning that dependence to the genetic component by construction of the fitted model. The table reports the actual finite-simulation behaviour rather than assuming perfect separation. Its genetic correlations target a covariance component of the synthetic genome and are distinct from either the phenotypic correlation or the local copula parameter.

{null_table}

Table 11. Fresh zero-genetic-correlation results, pooling the two populations within each sample-size and residual-correlation condition. Each row contains 400 independent cohorts. The residual correlation is on the generating liability scale; the observed phenotype correlation uses the binary diagnosis proxy.

![Environmental and genetic dependence](v4/figures/environmental_separation.png)

Figure 11. Observed phenotypic correlation and estimated genetic correlation when generating genetic covariance is zero. Error bars show 1.96 Monte Carlo standard errors for the stratified replicate means, conditional on the fixed kernels. Positive environmental dependence can produce comorbidity without a genetic-correlation parameter different from zero.

The larger independent calibration banks reduce threshold sampling uncertainty without guaranteeing a rejection proportion below 5% in every finite cell. In the fresh binary analysis the aggregate null rate is {pct(poly(600,0)['calibrated_rate'])} at n=600 and {pct(poly(1200,0)['calibrated_rate'])} at n=1,200. Calibration-aware bootstrap intervals and all individual cell rates are provided. These thresholds apply to the specified synthetic model and do not validate genome-wide GWAS significance tails.

![Fresh genetic covariance calibration](v4/figures/genetic_correlation_calibration.png)

Figure 12. Fresh held-out null rejection rates for nominal chi-square tests and empirical calibration. Error bars include calibration and evaluation resampling uncertainty for the calibrated rules; the dashed line marks 0.05. The original smaller-bank results remain available for comparison.

The final extension audit includes 6,240 fine-mapping analyses and 33,600 polygenic outcome-fit records across 16,800 unique polygenic cohorts. The two outcome modes are paired observations of each cohort. Independent R checks cover both exact-liability and observed-binary likelihoods, and dedicated established-method controls verify sharing with opposing directions and signed-LD harmonisation. Fit failures and boundary counts are reported in the audit files.

'''
    text=text.replace('## Discussion\n',results+'## Discussion\n',1)
    text=text.replace('Multi-signal fine-mapping and methods such as coloc with SuSiE are appropriate next comparisons [6]; they were not executed in this version. Such methods would still face indistinguishability at perfectly correlated variants.','The additional SuSiE comparison permits multiple signals and evaluates this limitation directly [6,7]. Its finite reference panels, LD regularisation and small comparison subset remain constraints. Multi-signal methods still face indistinguishability at perfectly correlated variants.')
    more_disc='''The polygenic experiment adds a third inferential level. Genetic covariance and correlation summarise effects across the synthetic genome; they do not identify a causal SNP, prove a disease-to-disease effect or count the fraction of variants shared between diseases. Conversely, local shared causal effects can coexist with near-zero genome-scale correlation when positive and negative contributions cancel. These targets require separate interpretation even when they are evaluated in one research programme.

The genome-scale validation is deliberately favourable in several respects: it uses common variants, 704 markers, 22 independent blocks, fixed ancestry-specific kernels and a generating covariance matched to the fitted relationship model. The continuous benchmark observes a liability that would be unavailable clinically. The binary model is explicitly a working approximation, with its error and coverage evaluated rather than assumed. These results are not evidence that the same sample sizes would give comparable precision in a real African GWAS with millions of variants, complex ascertainment or heterogeneous genetic architecture. A real application would require cohort-scale genomic QC, appropriate relationship/LD models, phenotype validation and independent replication.

The validation extension was designed after examining the first experiment's calibration uncertainty. Its new calibration and evaluation cohorts are independent, and the original findings are retained. This provides a transparent prospective check of the unchanged procedure, rather than treating a post hoc revision as an entirely preregistered study.

'''
    text=text.replace('No real shared genetic cause of cardiometabolic and mental illness has been established here.',more_disc+'No real shared genetic cause of cardiometabolic and mental illness has been established here.',1)
    text=text.replace('The resulting calibration and uncertainty analyses define the conditions under which sharing claims are supported and where causal resolution remains unavailable.','The additional multi-signal and polygenic analyses implement all three levels with distinct tests, explicit uncertainty and known synthetic truth. They define conditions under which each inference is supported and where causal resolution remains unavailable.')
    text=text.replace('The frozen design, source identity, public-reference provenance, raw simulation chunks, all evaluation tables and numerical audits are versioned.','The frozen designs, source identities, public-reference provenance, raw simulation chunks, all evaluation tables and numerical audits are versioned. Level 2 and the initial Level 3 experiment were frozen in commit 85b0bc4; the larger-bank fresh validation was frozen separately in commit dd84837.')
    repro='''The integrated Level 2 and Level 3 commands follow the local experiment. The full workflow explains the sequential calibration gates, R package versions, independent validation, checkpoint preservation and Word rendering. Fine-mapping is executed in R. Polygenic REML and repeated simulations are implemented in Python and independently checked in R. The original smaller-bank polygenic results remain at [the initial aggregate table](v4/results/polygenic_aggregate.csv); primary fresh validation is at [the confirmatory aggregate table](v4/results/polygenic_validation/polygenic_aggregate.csv).

```bash
python v4/src/freeze.py
python v4/src/export_finemap.py --phase cal
Rscript --vanilla v4/src/finemap.R v4 cal
python v4/src/finemap_summary.py --freeze
python v4/src/export_finemap.py --phase eval
Rscript --vanilla v4/src/finemap.R v4 eval
python v4/src/finemap_summary.py
python v4/src/polygenic.py
python v4/src/polygenic_summary.py
python v4/src/polygenic_validation.py
python v4/src/polygenic_summary.py --validation
Rscript --vanilla v4/src/check_reml.R v4
Rscript --vanilla v4/src/check_finemap.R v4
python v4/src/audit.py
python v4/src/figures.py
python v4/src/manuscript.py
python v4/src/build_word.py --pandoc /path/to/pandoc
```

These commands reuse already verified outputs in a clone. Fresh reference/kernel rebuilding and a completely fresh simulation should use a separate checkout and archived output directories, as detailed in workflow.txt. Each calibration phase must precede its evaluation. Published sources are frozen; scientific modifications belong in a new version.

'''
    text=text.replace('## References\n',repro+'## References\n',1)
    text=text.replace('python v3/src/manuscript.py\npython v3/src/build_word.py --pandoc /path/to/pandoc\n','')
    text+='\n7. Wallace C. A more accurate method for colocalisation analysis allowing for multiple causal variants. PLoS Genetics (2021). [Method and implementation](https://pmc.ncbi.nlm.nih.gov/articles/PMC8504726/).\n8. Lee SH et al. Estimation of pleiotropy between complex diseases using single-nucleotide polymorphism-derived genomic relationships and restricted maximum likelihood. Bioinformatics (2012). [Primary publication](https://pmc.ncbi.nlm.nih.gov/articles/PMC3463125/).\n9. Yang Lab. [GCTA bivariate GREML documentation](https://yanglab.westlake.edu.cn/software/gcta). Cited for methodological context; this study uses an independent implementation and R verification.\n'
    assert len(re.findall(r'^!\[',text,re.M))==12
    (ROOT/'README.md').write_text(text,encoding='utf-8',newline='\n')
    prior=ROOT/'v3/results/manuscript_sources.json';state=json.loads(prior.read_text());state['stage']='Intermediate v3-only manuscript before Level 2/3 integration; regenerate with v3/src/manuscript.py';dump(prior,state)
    dump(V/'results/manuscript_sources.json',dict(source='v4/src/manuscript.py',intermediate_v3_readme_sha256=base_sha,readme_sha256=sha(ROOT/'README.md'),local_audit_sha256=sha(ROOT/'v3/results/audit.json'),extension_audit_sha256=sha(V/'results/audit.json'),fine_comparison_sha256=sha(V/'results/finemap_matched_comparison.csv'),primary_polygenic_rates_sha256=sha(V/'results/polygenic_validation/polygenic_aggregate.csv')))
    print('Integrated three-level research manuscript generated')
if __name__=='__main__':main()
