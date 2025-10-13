# Linkage disequilibrium

Get linkage disequilibrium statistics using vcftools.

## Obtain R2 stats per population

```bash
module load PaM/environment
module load bsub.py/0.42.1
module load vcftools/0.1.16-c4

cd /lustre/scratch127/pam/teams/team333/rp24/dirofilaria_immitis/R2/LD

ln -s ../vcf/dirofilaria_global.cohort.2025-06-18.nuclearSNPs..final.autosomal.recode.vcf.keep_samples.dimmitis-only.recode.vcf

# loop through each continent
# removed replicate samples
bsub.py 40 pop_ld "./pop_ld.sh"
#!/bin/bash
cd /lustre/scratch127/pam/teams/team333/rp24/dirofilaria_immitis/R2/LD
for i in ASIA AUS CENAM EUR USA; do
  vcftools --vcf dirofilaria_global.cohort.2025-06-18.nuclearSNPs..final.autosomal.recode.vcf.keep_samples.dimmitis-only.recode.vcf \
    --maf 0.02 \
    --max-missing 1 \
    --keep ${i}_samplelist.keep \
    --recode \
    --out dirofilaria_global.cohort.2025-06-18.nuclearSNPs..final.autosomal.recode.vcf.keep_samples.dimmitis-only.${i} ;
done


for i in ASIA AUS CENAM EUR USA; do
  vcftools --vcf dirofilaria_global.cohort.2025-06-18.nuclearSNPs..final.autosomal.recode.vcf.keep_samples.dimmitis-only.${i}.recode.vcf \
    --ld-window-bp 1000000 \
    --max-alleles 2 \
    --min-alleles 2 \
    --geno-r2 \
    --out ${i}_LD ;
done
```

--ld-window-bp: maximum number of physical bases between the SNPs being tested for LD
--thin: thin sites so that no two sites are within the specified distance from one another. --> EDIT: not appropriate to do this for LD
--min-r2 0.1 --> don't set a minimum so we capture the whole decay curve



## Plot

```bash
bsub.py 400 plot_ld "Rscript ./LD.R"
```


```R
version
# version 4.4.0

library(ggplot2)
library(dplyr)

# Function to process LD data
process_ld_data <- function(file_path) {
  # Read data
  ld_data <- read.table(file_path, header = TRUE)
  
  # Make distance column by calculating absolute distance between SNP1 and SNP2, converts this distance into kb
  ld_data <- ld_data %>% mutate(Distance = abs(POS2 - POS1) / 1000)
  
   # Separate into blocks
  ld_bin <- ld_data %>%
    mutate(Bin = cut(Distance, breaks = seq(0, max(Distance) + 0.1, by = 1), include.lowest = TRUE)) %>%
    group_by(Bin) %>%
    summarize(Mean_R.2 = mean(R.2, na.rm = TRUE),
              Mean_Distance = mean(Distance))
  # each bin is 100bp or 0.1 kb

  # Add cohort identifier
  ld_bin$Pop <- gsub("_LD.geno.ld", "", basename(file_path))
  
  return(ld_bin)
}

# File paths for the five cohorts
pop_files <- c("USA_LD.geno.ld", "EUR_LD.geno.ld", "CENAM_LD.geno.ld", "AUS_LD.geno.ld", "ASIA_LD.geno.ld")

# Process each cohort
pop_data <- lapply(pop_files, process_ld_data)

# Combine the results into a single data frame
all_data <- bind_rows(pop_data)

all_data <- all_data %>%
  mutate(Location = sub("_.*", "", Pop))

# Define custom colors for each cohort
pop_colours <- c("USA" = "tomato2", "EUR" = "forestgreen", "CENAM" = "purple", "AUS" = "cornflowerblue", "ASIA" = "hotpink")

# Plot
plot <- ggplot(all_data, aes(x = Mean_Distance, y = Mean_R.2, group = Pop, color = Pop)) +
  geom_line(size = 1) +
  guides(color = guide_legend(override.aes = list(size = 10, linewidth=3))) +
  labs(x = "Distance (kb)", y = expression(R^2), color = "Population") +
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        legend.key = element_rect(fill = "white"),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 18),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18)) +
  scale_y_continuous(limits = c(0, 0.8), breaks = seq(0, 1, by = 0.2)) +
  scale_color_manual(values = pop_colours)

# Save the plot
ggsave("plot_LD.pdf", plot, dpi = 300, height = 6, width = 10)
ggsave("plot_LD.tif", plot, dpi = 300, height = 6, width = 10)
ggsave("plot_LD.png", plot, dpi = 300, height = 6, width = 10)


# Plot with a smoother line
plot <- ggplot(all_data, aes(x = Mean_Distance, y = Mean_R.2, group = Pop, color = Pop)) +
  geom_point(alpha = 0.1, size = 1) +  # Adjust alpha for point visibility
  geom_smooth(method = "loess", se = FALSE) +  # Add smooth line for better visualization
  guides(color = guide_legend(override.aes = list(size = 10, linewidth=3))) +
  labs(x = "Distance (kb)", y = expression(R^2), title = "Linkage Disequilibrium") +
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        legend.key = element_rect(fill = "white"),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 18),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18)) +
  scale_y_continuous(limits = c(0, 0.8), breaks = seq(0, 1, by = 0.2)) +
  scale_color_manual(values = pop_colours)

# Save the plot
ggsave("plot_LD_smooth.pdf", plot, dpi = 300, height = 6, width = 8)
ggsave("plot_LD_smooth.tif", plot, dpi = 300, height = 6, width = 8)
ggsave("plot_LD_smooth.png", plot, dpi = 300, height = 6, width = 8)
```

LD for most populations remains quite high at long distances, which is as bit strange. It should decay towards zero. Run kinship analysis.


## Plot showing LD(1/2)

LD(1/2) is the distance between SNPs in which the LD is half the maximum value. Do for each population, it will allow us to be more quantitative about different rates of decay.

```R
version
# version 4.4.0

#install.packages("rlang")
#install.packages("dplyr", dependencies = TRUE)
#install.packages("ggplot2", dependencies = TRUE)

library(ggplot2)
library(dplyr)

# get points for each pop where LD is 1/2 the max value
LD_Half <- all_data %>%
group_by(Pop) %>%
summarise(LD_Background = mean(Mean_R.2[Mean_Distance > 800]), # get LD background for each pop
          LD_Max = max(Mean_R.2), #get max LD for each pop
          LD_Half = (0.5 * (LD_Max - LD_Background) + LD_Background), # Get half of max R2
          Distance_Half = min(Mean_Distance[Mean_R.2 <= LD_Half])) # get the ~distance where LD is half the maximum value

print(LD_Half)
'
# A tibble: 5 × 5
  Pop   LD_Background LD_Max LD_Half Distance_Half
  <chr>         <dbl>  <dbl>   <dbl>         <dbl>
1 ASIA          0.420  0.755   0.587          31.5
2 AUS           0.287  0.603   0.445          49.5
3 CENAM         0.359  0.643   0.501          42.5
4 EUR           0.315  0.586   0.451          37.5
5 USA           0.175  0.413   0.294          36.5
'

# Plot with vertical lines for each pop where the distance is half LD
plot <- ggplot(all_data, aes(x = Mean_Distance, y = Mean_R.2, group = Pop, color = Pop)) +
  geom_line(size = 1) +
  guides(color = guide_legend(override.aes = list(size = 10, linewidth=3))) +
  geom_segment(data=LD_Half, aes(x=Distance_Half, y=LD_Half, xend=Distance_Half, yend=0, color = Pop), size = 0.6) +
  labs(x = "Distance (kb)", y = expression(R^2), color = "Population") +
  theme_bw() +
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 2),
        legend.key = element_rect(fill = "white"),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 18),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18)) +
  scale_y_continuous(limits = c(0, 0.8), breaks = seq(0, 1, by = 0.2)) +
  scale_color_manual(values = pop_colours)

# Save the plot
ggsave("plot_LD_half_v2.pdf", plot, dpi = 300, height = 6, width = 7)
ggsave("plot_LD_half_v2.png", plot, dpi = 300, height = 6, width = 7)
ggsave("plot_LD_half_v2.tif", plot, dpi = 300, height = 6, width = 7)

```








### Kinship analysis

```bash
module load vcftools/0.1.16-c4

cd /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/05_ANALYSIS/KINSHIP

ln -s /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/04_VARIANTS/FILTER1/NO_OUTGROUPS/FINAL_SETS/nuclear_samples3x_missing0.9.recode.vcf

# calculate relatedness
bsub.py 20 vcf_kinship "vcftools --vcf nuclear_samples3x_missing0.9.recode.vcf \
--relatedness2 --max-missing 1"

# only want within population comparisons (so within each country)
for i in ASIA AUS EUR CENAM USA ; do
     awk -v name=$i '$9~name && $10~name {print $0}' OFS="\t" relatedness_with_pop.tsv >> allpops.relatedness2;
done
```
