Supplementary Code S4: Introgression
================
Tauana Cunha \| <https://tauanacunha.com>
August 2021

<font size="4">Cunha et al. 2021. Investigating Sources of Conflict in
Deep Phylogenomics of Vetigastropod Snails. Systematic Biology.</font>

This document contains the R code related to the introgression
analysis. See the main article for more information.

------------------------------------------------------------------------

# Required libraries

``` r
.libPaths(new = "libR/")
library(ggplot2)
library(tidyverse)
library(data.table)
```

# Read the data

Import observed concordance factors. Legend to column names are at the
bottom of this file.

``` r
d_cat = read.delim("../results/*concordance_factors/iqtree-cat_31taxa/ConcFact.cf.stat",
               header = T, comment.char = '#')
d_part = read.delim("../results/*concordance_factors/iqtree-part_31taxa/ConcFact.cf.stat",
               header = T, comment.char = '#')

d_cat$tree = "Topology 1" # M1-IQTREEcat
d_part$tree = "Topology 2" # M1-IQTREEpart

d = rbind(d_cat, d_part)
```

Import concordance factors calculated from resampling gene trees (with
replacement) 2000 times:

``` r
read_plus <- function(flnm){ # function to save filename in a new column
    fread(flnm) %>% 
        mutate(filename = basename(flnm))}

resampled_cat = list.files(path = "../results/introgression/concordance_factors/iqtree-cat_31taxa/all_CF_resamples/",
                           pattern = "*.cf.stat", full.names = T) %>% 
    map_df(~read_plus(.)) # use new function instead of fread

resampled_part = list.files(path = "../results/introgression/concordance_factors/iqtree-part_31taxa/all_CF_resamples/",
                           pattern = "*.cf.stat", full.names = T) %>% 
    map_df(~read_plus(.))

resampled_cat$tree = "Topology 1" # M1-IQTREEcat
resampled_part$tree = "Topology 2" # M1-IQTREEpart

resampled = rbind(resampled_cat, resampled_part)
```

# Test statistic delta

Calculate test statistic delta from Vanderpool et al. 2020 on observed
CFs:

``` r
observed_delta = d %>%
  rowwise() %>%
  mutate(prop_discord_trees = sum(gDF1, gDF2),
         delta = if_else((gDF1_N + gDF2_N)==0, 0,
                         abs((gDF1_N - gDF2_N)/(gDF1_N + gDF2_N))))
```

# Nodes of interest

Select nodes with more than 5% discordant gene trees to investigate
introgression:

``` r
# A lot of discordant gene trees can indicate a number of things, ILS, introgression, error in the gene trees etc
most_discord = observed_delta %>%
  filter(prop_discord_trees > 5)
most_discord # 43 of the 59 nodes have over 5% of discordant trees
```

    ## # A tibble: 86 x 22
    ## # Rowwise: 
    ##       ID   gCF gCF_N  gDF1 gDF1_N  gDF2 gDF2_N  gDFP gDFP_N    gN   sCF sCF_N
    ##    <int> <dbl> <int> <dbl>  <int> <dbl>  <int> <dbl>  <int> <int> <dbl> <dbl>
    ##  1    63  32.8   251  6.27     48  8.62     66  52.4    401   766  49.6  850.
    ##  2    64  20.2   167 12.1     100  7.02     58  60.6    501   826  30.7  355.
    ##  3    65  21.5    29  5.19      7  1.48      2  71.8     97   135  48.1  126.
    ##  4    67  19.4   161  4.81     40  5.65     47  70.2    584   832  37.0  504.
    ##  5    69  12.6    91  9.85     71  4.85     35  72.7    524   721  29.2  358.
    ##  6    70  29.9   181  3.64     22 10.9      66  55.5    336   605  44.9  503.
    ##  7    71  59.0   411  1.87     13  3.16     22  36.0    251   697  54    481.
    ##  8    72  46.3   374  4.33     35  7.18     58  42.2    341   808  45.8  563.
    ##  9    73  18.1   150 13.1     108 12.4     103  56.4    466   827  33.9  448.
    ## 10    74  14.7   108 15.0     110 16.1     118  54.2    398   734  34.7  606.
    ## # … with 76 more rows, and 10 more variables: sDF1 <dbl>, sDF1_N <dbl>,
    ## #   sDF2 <dbl>, sDF2_N <dbl>, sN <dbl>, Label <int>, Length <dbl>, tree <chr>,
    ## #   prop_discord_trees <dbl>, delta <dbl>

Just for plotting, select the two conflicting nodes of interest:

``` r
key_nodes = observed_delta %>%
  ungroup() %>%
  filter(ID %in% c(69,87) & tree %in% "Topology 1" |
           ID %in% c(69,70) & tree %in% "Topology 2") %>%
  mutate(nodeID = c(2,1,4,3),# Node numbers used in Fig 3 (phylogeny)
         node_res = c("Node 2 or 4","Node 1 or 3","Node 2 or 4","Node 1 or 3"))
```

Plot proportion of discordant trees, highlighting the two conflicting
nodes:  
(continuous line represents nodes 1 and 3 of Figure 3 in the paper,
dashed line represents nodes 2 and 4 of Figure 3)

``` r
p_discord = ggplot(observed_delta) +
  geom_histogram(aes(x = prop_discord_trees),
                 bins = 30, alpha = 0.7, fill = "#1b9e77") +
  facet_wrap(~tree) +
  geom_vline(data = key_nodes,
             aes(xintercept = prop_discord_trees, linetype = node_res),
             col = "#B50033", size = 0.15, show.legend = F) +
  labs(x = "Proportion of discordant trees", y = "Number of branches") +
  theme_minimal()
p_discord
```

![](introgression_files/figure-gfm/plot%20discordant%20trees-1.png)<!-- -->


# Z-score, resampling

Calculate delta score on resampled concordance factors:

``` r
resampled_delta = resampled %>%
  rowwise() %>%
  mutate(delta = if_else((gDF1_N + gDF2_N)==0, 0,
                         abs((gDF1_N - gDF2_N)/(gDF1_N + gDF2_N))))
```

Calculate mean and sd of delta for each node based on the 2000 resamples
(a null distribution):

``` r
nulldist = resampled_delta %>%
  group_by(tree, ID) %>% # group by node/resample
  summarize(mean = mean(delta), sd = sd(delta)) # mean and sd for all resamples of that node
nulldist %>%
  group_split
```

    ## <list_of<
    ##   tbl_df<
    ##     tree: character
    ##     ID  : integer
    ##     mean: double
    ##     sd  : double
    ##   >
    ## >[2]>
    ## [[1]]
    ## # A tibble: 59 x 4
    ##    tree          ID  mean     sd
    ##    <chr>      <int> <dbl>  <dbl>
    ##  1 Topology 1    63 0.164 0.0854
    ##  2 Topology 1    64 0.268 0.0757
    ##  3 Topology 1    65 0.570 0.272 
    ##  4 Topology 1    66 0.441 0.364 
    ##  5 Topology 1    67 0.109 0.0789
    ##  6 Topology 1    68 0.749 0.167 
    ##  7 Topology 1    69 0.339 0.0915
    ##  8 Topology 1    70 0.503 0.0919
    ##  9 Topology 1    71 0.263 0.150 
    ## 10 Topology 1    72 0.251 0.101 
    ## # … with 49 more rows
    ## 
    ## [[2]]
    ## # A tibble: 59 x 4
    ##    tree          ID   mean     sd
    ##    <chr>      <int>  <dbl>  <dbl>
    ##  1 Topology 2    63 0.164  0.0854
    ##  2 Topology 2    64 0.268  0.0757
    ##  3 Topology 2    65 0.570  0.272 
    ##  4 Topology 2    66 0.441  0.364 
    ##  5 Topology 2    67 0.109  0.0789
    ##  6 Topology 2    68 0.183  0.139 
    ##  7 Topology 2    69 0.0751 0.0570
    ##  8 Topology 2    70 0.202  0.0648
    ##  9 Topology 2    71 0.434  0.100 
    ## 10 Topology 2    72 0.263  0.150 
    ## # … with 49 more rows

For each branch, calculate a standardized z-score of observed deltas
relative to the null distribution created with the resampling, then
check probability of such a score and calculate the pvalue. One-tailed
test looking for CDF\_P&gt;0.95 (pvalue&lt;0.05) as evidence of
introgression:

``` r
test = most_discord %>%
  left_join(nulldist) %>%
  mutate(zscore = (delta - mean)/sd,
         CDF_prob = pnorm(zscore),
         pvalue = 1-CDF_prob)

test %>%
  select(-starts_with(c("g","s","m","L"))) %>%
  group_by(tree) %>%
  arrange(desc(CDF_prob)) %>%
  group_split()
```

    ## <list_of<
    ##   tbl_df<
    ##     ID                : integer
    ##     tree              : character
    ##     prop_discord_trees: double
    ##     delta             : double
    ##     zscore            : double
    ##     CDF_prob          : double
    ##     pvalue            : double
    ##   >
    ## >[2]>
    ## [[1]]
    ## # A tibble: 43 x 7
    ##       ID tree       prop_discord_trees delta   zscore CDF_prob pvalue
    ##    <int> <chr>                   <dbl> <dbl>    <dbl>    <dbl>  <dbl>
    ##  1    92 Topology 1              24.0  0.294  0.0569     0.523  0.477
    ##  2   101 Topology 1               9.98 0.613  0.0375     0.515  0.485
    ##  3    82 Topology 1              42.8  0.240  0.0159     0.506  0.494
    ##  4    69 Topology 1              14.7  0.340  0.00327    0.501  0.499
    ##  5    75 Topology 1              13.4  0.168  0.00130    0.501  0.499
    ##  6    90 Topology 1              10.2  0.7   -0.0144     0.494  0.506
    ##  7    81 Topology 1              41.9  0.446 -0.0153     0.494  0.506
    ##  8    98 Topology 1              12.6  0.3   -0.0238     0.491  0.509
    ##  9    64 Topology 1              19.1  0.266 -0.0240     0.490  0.510
    ## 10   100 Topology 1               7.7  0.571 -0.0290     0.488  0.512
    ## # … with 33 more rows
    ## 
    ## [[2]]
    ## # A tibble: 43 x 7
    ##       ID tree       prop_discord_trees delta   zscore CDF_prob pvalue
    ##    <int> <chr>                   <dbl> <dbl>    <dbl>    <dbl>  <dbl>
    ##  1    92 Topology 2              24.0  0.294  0.0569     0.523  0.477
    ##  2   101 Topology 2               9.98 0.613  0.0375     0.515  0.485
    ##  3    77 Topology 2              39.8  0.182  0.0163     0.506  0.494
    ##  4    74 Topology 2              29.5  0.228  0.00704    0.503  0.497
    ##  5    90 Topology 2              10.2  0.7   -0.0144     0.494  0.506
    ##  6    76 Topology 2              44.2  0.353 -0.0154     0.494  0.506
    ##  7    71 Topology 2              10.2  0.432 -0.0203     0.492  0.508
    ##  8    98 Topology 2              12.6  0.3   -0.0238     0.491  0.509
    ##  9    64 Topology 2              19.1  0.266 -0.0240     0.490  0.510
    ## 10    70 Topology 2              22    0.2   -0.0289     0.488  0.512
    ## # … with 33 more rows

# Extra resources

Column names/meaning for output of IQTREE concordance factors:

``` r
#   ID: Branch ID
#   gCF: Gene concordance factor (=gCF_N/gN %)
#   gCF_N: Number of trees concordant with the branch
#   gDF1: Gene discordance factor for NNI-1 branch (=gDF1_N/gN %)
#   gDF1_N: Number of trees concordant with NNI-1 branch
#   gDF2: Gene discordance factor for NNI-2 branch (=gDF2_N/gN %)
#   gDF2_N: Number of trees concordant with NNI-2 branch
#   gDFP: Gene discordance factor due to polyphyly (=gDFP_N/gN %)
#   gDFP_N: Number of trees decisive but discordant due to polyphyly
#   gN: Number of trees decisive for the branch
#   sCF: Site concordance factor averaged over 100 quartets (=sCF_N/sN %)
#   sCF_N: sCF in absolute number of sites
#   sDF1: Site discordance factor for alternative quartet 1 (=sDF1_N/sN %)
#   sDF1_N: sDF1 in absolute number of sites
#   sDF2: Site discordance factor for alternative quartet 2 (=sDF2_N/sN %)
#   sDF2_N: sDF2 in absolute number of sites
#   sN: Number of informative sites averaged over 100 quartets
#   Label: Existing branch label
#   Length: Branch length
```
