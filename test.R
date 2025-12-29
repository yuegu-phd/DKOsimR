# Title: Double-CRISPR Knockout Simulation (DKOsim) - R package testing
# Author: Yue Gu, Luis Leon Novelo, John Paul Shen, Traver Hart
# install and load required packages
if(!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")
devtools::install_github("yuegu-phd/DKOsimR")
devtools::install(dependencies = TRUE)
### Major Modifications Before Each Run
### 1. Simulation Name: sample_name
### 2. Simulation Settings: initialized parameter

################################### Initialize Parameters ##########################################################
# library(DKOsimR)
# set.seed(888) # random seed
#
# # name of the simulation run
# sample_name = "test"
#
# ## initialized library parameters
# coverage = 100 # coverage: cell representation per guide
# n = 120 # number of unique single gene
# n_guide_g = 3 # number of guide per gene
# n_gene_pairs = n * (n-1) / 2 + n  # number of unique gene pairs (both SKO and DKO)
# n_construct = (n*n_guide_g) * ((n-1)*n_guide_g) / 2 + n*n_guide_g  # total number of constructs
# library_size = n_construct * coverage # number of total cells in the initialized gene-level library
#
# moi = 0.3 # moi
# moi_pois = dpois(1, moi) # get the number of viral particles delivered per cell during transfection from Poisson(moi) to calculate resampling size
# p_gi = 0.03 # % of genetic interaction presence
# sd_gi = 1.5 # std dev of re-sampled phenoytpes w/ gi presence
# sd_freq0 = 1/3.29 # dispersion of initial counts
# ## guide parameters
# p_high = 1 # % of high-efficacy guides
# mode = "CRISPRn-100%Eff" # CRISPR mode: use CRISPRn-100%Eff for need 100% efficient guides; use CRISPRn for high-efficient guides draw from distribution
#
# ## gene class parameters:
# ## gene class parameters: % of theoretical phenotype to each gene class - add up to 1
# pt_neg = 0.15 # % negative phenotype
# pt_pos = 0.05 # % positive phenotype
# pt_wt = 0.75 # % wt
# pt_ctrl = 0.05 # % non-targeting control
# ### mean and std dev of theoretical phenotypes
# mu_neg = -0.75 # mean: negative phenotype
# sd_neg = 0.1 # std dev: negative phenotype
#
# mu_pos = 0.75 # mean: positive phenotype
# sd_pos = 0.1 # std dev: positive phenotype
#
# sd_wt = 0.25 # std dev: wildtype phenotype
#
# ## bottleneck parameters
# bottleneck = 2 * library_size # bottleneck size
# n.bottlenecks = 1 # how many times do we encountering bottlenecks?
# n.iterations = 30 # assuming a maximum of 30 doubling cycles if we didnt encounter bottleneck
# resampling = round(moi_pois * bottleneck) # determine resampling size based on moi and bottleneck size

################################### Run Simulations ##########################################################
# test on all parameters included
library(DKOsimR)
dkosim(sample_name="test",
       coverage=100,
       n=60,
       n_guide_g=3,
       sd_freq0 = 1/3.29,
       moi = 0.3,
       p_gi=0.03,
       sd_gi=1.5,
       p_high=1,
       mode="CRISPRn-100%Eff",
       pt_neg=0.15,
       pt_pos=0.05,
       pt_wt=0.75,
       pt_ctrl=0.05,
       mu_neg=-0.75,
       sd_neg=0.1,
       mu_pos=0.75,
       sd_pos=0.1,
       sd_wt=0.25,
       size.bottleneck = 3,
       n.bottlenecks= 2,
       n.iterations = 30,
       rseed = 111,
       path = ".")

# test on parameters without default values: sample_name, n
dkosim(sample_name="test", n=40, mode="CRISPRn")

# test on the lab approximating mode with all parameters included
dkosim_lab(sample_name="test_lab",
           coverage=100,
           n=60,
           n_guide_g=3,
           sd_freq0 = 1/3.29,
           moi = 0.3,
           p_gi=0.03,
           sd_gi=1.5,
           p_high=1,
           mode="CRISPRn-100%Eff",
           pt_neg=0.15,
           pt_unknown=0.80,
           pt_ctrl=0.05,
           mu_neg=-0.75,
           sd_neg=0.1,
           sd_unknown=0.25,
           size.bottleneck = 3,
           n.bottlenecks= 2,
           n.iterations = 30,
           rseed = 111,
           path = ".")



# test on the lab approximating mode with default parameter values used to mimic Fong-2024
dkosim_lab(sample_name="test_lab")
