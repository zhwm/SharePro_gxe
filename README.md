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

1. the **exposure stratified summary statistic files** that contain at least four columns: SNP/N/BETA/SE. Examples are available at [dat/C21.txt](dat/C21.txt) and [dat/L21.txt](dat/L21.txt). Multiple exposure categories are allowed.

3. the **LD files** that contain correlation coefficient matrix. **Please make sure the REF/ALT alleles used in calculating LD are the same as the GWAS study!!** An example can be found at [dat/Locus1.ld](dat/Locus1.ld).

## Usage examples

We use `--z` and  `--ld` to indicate path to the summary statistic files and the LD file.
We use `--save` to specify path to save result.

```
python sharepro_gxe.py \
--z dat/C21.txt dat/L21.txt \
--ld dat/Locus1.ld \
--save dat/CL21
```

## Output interpretation

In this simulated example, we have one causal variant: rs112819506 with causal effect sizes of 0.05 and 0.02 in the exposed and unexposed group.

From the output we obtained below, we have successfully identified one effect group consisting of two variants rs112819506 and rs138116565 with weights of 0.7345 and 0.2536. Since those two variants have a Pearson correlation of 0.99, they are nearly statistically indistinguishable.
Based on the estimated effect size of 0.0538 and 0.0207, a GxE p-value of 2.15e-04 was derived for this effect group.

```
2024-06-03 20:22:38 - z: ['dat/C21.txt', 'dat/L21.txt']
2024-06-03 20:22:38 - ld: dat/Locus1.ld
2024-06-03 20:22:38 - save: dat/CL21
2024-06-03 20:22:38 - sigma: 0.99
2024-06-03 20:22:38 - K: 10
2024-06-03 20:22:38 - maxite: 100
2024-06-03 20:22:38 - eps: 0.01
2024-06-03 20:22:38 - ubound: 10000000000.0
2024-06-03 20:22:38 - cthres: 0.95
2024-06-03 20:22:38 - pthres: 0.8
2024-06-03 20:22:38 - **********************************************************************
2024-06-03 20:22:38 - Iteration-->0 . Likelihood: 47.5 . KL_b: -25.9 . KL_c: -0.3 . KL_s: 65.8 . ELBO: 87.1
2024-06-03 20:22:38 - **********************************************************************
2024-06-03 20:22:38 - Iteration-->1 . Likelihood: 46.7 . KL_b: -25.8 . KL_c: -0.3 . KL_s: 66.8 . ELBO: 87.5
2024-06-03 20:22:39 - **********************************************************************
2024-06-03 20:22:39 - Iteration-->2 . Likelihood: 46.6 . KL_b: -25.7 . KL_c: -0.3 . KL_s: 67.0 . ELBO: 87.6
2024-06-03 20:22:39 - **********************************************************************
2024-06-03 20:22:39 - Iteration-->3 . Likelihood: 46.6 . KL_b: -25.7 . KL_c: -0.3 . KL_s: 67.1 . ELBO: 87.7
2024-06-03 20:22:39 - **********************************************************************
2024-06-03 20:22:39 - Iteration-->4 . Likelihood: 46.6 . KL_b: -25.7 . KL_c: -0.3 . KL_s: 67.1 . ELBO: 87.7
2024-06-03 20:22:39 - Attainable coverage for effect groups: [1.   0.89 0.   0.   0.11 0.   0.   0.   0.   0.  ]
2024-06-03 20:22:39 - The 0-th effect group contains effective variants:
2024-06-03 20:22:39 - causal variants: ['rs112819506', 'rs138116565']
2024-06-03 20:22:39 - variant probabilities for this effect group: [0.7293, 0.2573]
```

We can additionally visualize both the raw GWAS summary statistics (A) and the GxE analysis results (B,D) in this locus.
In this example, the association test in the exposed and the unexposed categories exhibited clear differences (A) and the GxE test in GEM accurately detected the simulated effect heterogeneity (B). However, the power for this test diminished after adjusting for multiple testing (B). In contrast, SharePro identified a candidate effect group through joint fine-mapping (C), thus avoiding testing on every variant, resulting in a successful GxE detection (D).

<p align="center">
  <img src="doc/SharePro_gxe_example.png" alt="example image">
  <br>
  <em>Figure 1: SharePro usage example.</em>
</p>

## Output files

1. the **log** file. The log file for the example above has been provided at [dat/CL21.sharepro.gxe.log](dat/CL21.sharepro.gxe.log)

2. the **summary** file contains the posterior inclusion probabilities and credible sets information. The summary file for the example above has been provided at [dat/CL21.sharepro.gxe.txt](dat/CL21.sharepro.gxe.txt)


## Citations

If you find SharePro for GxE analysis useful, please cite:

[Wenmin Zhang, Robert Sladek, Yue Li, Hamed Najafabadi, Josée Dupuis. Accounting for genetic effect heterogeneity in fine-mapping and improving power to detect gene-environment interactions with SharePro.](https://doi.org/10.1101/2023.07.27.550862)