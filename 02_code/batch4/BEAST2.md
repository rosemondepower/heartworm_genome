# BEAST2 on mitochondria

BEAST version 2.7.7

Denver 2004 paper: https://doi.org/10.1038/nature02697
- mutation rate estimate of 2.1 × 10-8 mutations per site per generation (standard error of the mean, s.e.m. = ± 6.7 × 10-9)
-NUCLEAR rate

Denver 2009 paper: https://doi.org/10.1073/pnas.0904895106
- average base-substitution rate of 2.7 (±0.4) × 10−9 per site per generation
- NUCLEAR rate

Denver 2000 paper: https://doi.org/10.1126/science.289.5488.2342
- total mutation rate equal to 1.6 × 10−7 per site per generation (±3.1 × 10−8)
- MITO rate

More recent Konrad 2017 paper: https://doi.org/10.1093/molbev/msx051
- overall spontaneous mitochondrial mutation rate, μtotal, of 1.05 × 10−7 per site per generation (95% CI: 0.52 to 1.80 × 10−7).
- MITO rate

## BEAST RUN 1

- Sequences: Mito D. immitis & all outgroups (Dirofilaria & Oncho)
- Mutation rate: 1.6E-7 (upper 1.91E-7, lower 1.29E-7) (Denver 2000)
- Model: Birth-death, Yule
- Chain length 10000000, log every 1000

![](images/mito_outgroups_yule_tracer.PNG)

Some ESS values are below 200, need to run for more steps or alter parameters?

![](images/mito_outgroups_yule_BEASTTREE.png)


# BEAST RUN 2


- Sequences: Mito D. immitis & all outgroups (Dirofilaria & Oncho)
- Mutation rate: 1.6E-7 (upper 1.91E-7, lower 1.29E-7) (Denver 2000) -----> SLOW IT DOWN TO GENERATION = 1 YEAR
- C. elegans generation time is let's say ~ 1 week so D. immitis generation time of 1 year (52 weeks) is 52 times slower. So mutation rate would be 3.07E-9 (upper 3.67E-9, lower 2.48E-9)
- Model: Birth-death, Yule
- Chain length 10000000, log every 1000

![](images/mito_outgroups_yule_g1_tracer.PNG)

ESS values really aren't great. Run for longer.

![](images/mito_outgroups_yule_g1_BEASTTREE.png)

# BEAST RUN3 (RUN 2 FOR LONGER)

- Sequences: Mito D. immitis & all outgroups (Dirofilaria & Oncho)
- Mutation rate: 1.6E-7 (upper 1.91E-7, lower 1.29E-7) (Denver 2000) -----> SLOW IT DOWN TO GENERATION = 1 YEAR
- C. elegans generation time is let's say ~ 1 week so D. immitis generation time of 1 year (52 weeks) is 52 times slower. So mutation rate would be 3.07E-9 (upper 3.67E-9, lower 2.48E-9)
- Model: Birth-death, Yule
- Chain length 200000000, log every 20000