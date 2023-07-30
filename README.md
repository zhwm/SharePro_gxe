# SharePro for joint fine-mapping and GxE analysis

SharePro is a computational method to account for genetic effect heterogeneity in fine-mapping and improve power to detect gene-environment interactions (GxE). For analysis conducted in [the paper](https://doi.org/10.1101/2023.07.27.550862), please refer to [SharePro_gxe_analysis](https://github.com/zhwm/SharePro_gxe_analysis).

## Overview 

Characterizing genetic effect heterogeneity across subpopulations with different environmental exposures is useful in understanding disease heterogeneity within a population and further pinpointing modifiable risk factors for disease prevention and management. 
Classical GxE analysis can be used to detect genetic effect heterogeneity. However, it can have a high multiple testing burden in the context of genome-wide association studies (GWAS) and requires a large sample size to achieve sufficient power.

We developed SharePro for GxE analysis to reduce multiple testing burden.
Below we showcase three scenarios of potential effect heterogeneity in a population with three environmental exposure statuses. 
In setting 1, there is no effect heterogeneity and through the combined analysis, the causal variant can be identified while a stratified GWAS is under-powered. 
Through joint analysis of exposure-stratified GWAS summary statistics, we can recover this signal with SharePro.
In setting 2, the causal variant has different effect sizes across exposure categories and combined approaches will be disadvantaged. 
Classical GxE analysis is also under-powered due to a high multiple testing burden. 
With a joint approach, SharePro can accurately identify the causal signal and detect effect heterogeneity. Setting 3 is similar to setting 2 where the joint approach is more favorable than the combined approach.

<p align="center">
  <img src="doc/SharePro_gxe_overview.png" alt="example image">
  <br>
  <em>Figure 1: SharePro for GxE analysis overview.</em>
</p>

## Installation

SharePro was developed under Python 3.9.7 environment but should be compatible with older versions of Python 3. The following Python modules are required:

* [numpy](http://www.numpy.org/)
* [scipy](http://www.scipy.org/)
* [pandas](https://pandas.pydata.org/getpandas.html)

To install SharePro for joint fine-mapping and GxE analysis:

```
git clone https://github.com/zhwm/SharePro_gxe.git
cd SharePro_gxe
pip install -r requirements.txt 
``` 

To test the installation and display basic usage:
```
python sharepro_gxe.py -h
```

## Input files

Example input files from simulation studies are included in the [dat](dat/) directory.

SharePro takes in a summary file with path to z-score files and LD files, exposure stratified z-scores files, LD files as inputs.

1. the **summary file** contains two mandatory columns: names of z-score file and ld files. Multiple files are allowed and should be separated by comma. An example can be found at [dat/CL.zld.txt](dat/CL.zld.txt).

2. the **exposure stratified zscore files** that contain two mandatory columns: variant IDs and z-scores. Examples are available at [dat/C21.z](dat/C21.z) and [dat/L21.z](dat/L21.z). Multiple exposure categories are allowed.

3. the **LD files** that contain correlation coefficient matrix. **Please make sure the REF/ALT alleles used in calculating LD are the same as the GWAS study!!** An example can be found at [dat/Locus1.ld](dat/Locus1.ld) and a working script for matching raw GWAS summary statistics and PLINK bim file is provided [here](match_bim_ss.py).

## Usage examples

We use `--zld` to indicate path to the summary file and `--zdir` to indicate path to zscore files.
Additionally, we specify the sample sizes for exposure categories with `--N`.
We use `--save` to specify path to save result and `--prefix` to specify prefix of output files. We set the max number of causal signals as 5 with `--K`.

```
python sharepro_gxe.py \
--zld dat/CL.zld.txt \
--zdir dat \
--N 25000 25000 \
--save res \
--prefix CL \
--verbose \
--K 5
```

## Output interpretation

In this simulated example, we have one causal variant: rs112819506 with causal effect sizes of 0.05 and 0.02 in the exposed and unexposed group.

From the output we obtained below, we have successfully identified one effect group consisting of two variants rs112819506 and rs138116565 with weights of 0.7345 and 0.2536. Since those two variants have a Pearson correlation of 0.99, they are nearly statistically indistinguishable.
Based on the estimated effect size of 0.0538 and 0.0207, a GxE p-value of 2.15e-04 was derived for this effect group.

```
**********************************************************************
* SharePro for joint fine-mapping and GxE analysis                   *
* Version 1.0.0                                                      *
* (C) Wenmin Zhang (wenmin.zhang@mail.mcgill.ca)                     *
**********************************************************************
Using locus fine-mapping mode with --zld
LD list with 1 LD blocks loaded

processing C21.z,L21.z
**********************************************************************
Iteration-->0 . Likelihood: 42.3 . KL_b: -4.9 . KL_c: -12.3 . KL_s: 31.2 . ELBO: 56.4
**********************************************************************
Iteration-->1 . Likelihood: 42.2 . KL_b: -4.8 . KL_c: -12.1 . KL_s: 31.3 . ELBO: 56.5
**********************************************************************
Iteration-->2 . Likelihood: 42.2 . KL_b: -4.8 . KL_c: -12.1 . KL_s: 31.3 . ELBO: 56.6
**********************************************************************
Iteration-->3 . Likelihood: 42.2 . KL_b: -4.8 . KL_c: -12.1 . KL_s: 31.3 . ELBO: 56.6
Attainable coverage for effect groups: [1.   0.25 0.02 0.01 0.75]
The 0-th effect contains effective variants:
causal variants: ['rs112819506', 'rs138116565']
variant probabilities for this effect group: [0.7345, 0.2536]
causal effect sizes for traits: [0.0538, 0.0207]
GxE p-value: 2.15e-04
```

We can additionally visualize both the raw GWAS summary statistics (A) and the GxE analysis results (B,D) in this locus.
In this example, the association test in the exposed and the unexposed categories exhibited clear differences (A) and the GxE test in GEM accurately detected the simulated effect heterogeneity (B). However, the power for this test diminished after adjusting for multiple testing (B). In contrast, SharePro identified a candidate effect group through joint fine-mapping (C), thus avoiding testing on every variant, resulting in a successful GxE detection (D).

<p align="center">
  <img src="doc/SharePro_gxe_example.png" alt="example image">
  <br>
  <em>Figure 1: SharePro usage example.</em>
</p>

## Output files

1. the **effect group summary** (cs) file contains four columns: 
`cs` for variant representations in effect groups; 
`p_diff` for GxE p-values;
`beta` for effect size estimates;
`variantProb` for variant representation weight in effect groups.

```
$> cat res/C21.z_L21.z.cs 
cs	p_diff	beta	variantProb
rs112819506/rs138116565	2.15e-04	0.0538,0.0207	0.7345/0.2536
```

2. the **variant summary** (snp) file contains zscores and one additional column of posterior inclusion probabilities.

```
$> head -5 res/C21.z_L21.z.snp
SNP	C21.z	L21.z	vProb
rs111073422	0.6018	0.0611	4.16e-04
rs10414006	-1.9531	1.83	5.90e-04
rs188970225	1.1552	-0.4448	4.33e-04
rs814535	-1.8894	1.6999	5.47e-04
```

3. the **hyperparameters summary** (h2) file adds two additional columns to the summary file to record the heritability and effect size variance estimates.

```
$> cat res/CL.h2 
z	ld	h2	varb
C21.z,L21.z	Locus1.ld,Locus1.ld	4.41e-04	2.98e-03
```

## Citations

If you find SharePro for GxE analysis useful, please cite:

[Wenmin Zhang, Robert Sladek, Yue Li, Hamed Najafabadi, Josée Dupuis. Accounting for genetic effect heterogeneity in fine-mapping and improving power to detect gene-environment interactions with SharePro.](https://doi.org/10.1101/2023.07.27.550862)