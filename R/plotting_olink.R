# Load packages
source(here::here("R/package-loading.R"))

### Graph by Kleiner

# Load dataset
load(here::here("data/GALA.rda"))

GALA %>%
    filter(!is.na(as.numeric(kleinerFscore))) %>%
    filter(Assay == "IL8") %>%
    ggplot(aes(x = kleinerFscore, y = corr_NPX)) +
    geom_violin(draw_quantiles = c(0.5)) +
    geom_jitter() +
    ggtitle("IL8")

GALA %>%
    filter(Assay == "IL8") %>%
    ggplot(aes(x = kleiner, y = corr_NPX, color = cohort)) +
    geom_violin(draw_quantiles = c(0.5)) +
    geom_jitter(position = position_jitterdodge(), size = 0.5) +
    ggtitle("IL8")

# Graph by TE
GALA %>%
    filter(Assay == "IL8") %>%
    ggplot(aes(x = as.numeric(te), y = corr_NPX)) +
    geom_point(aes(color = cohort)) +
    geom_smooth(method = lm, color = "black") +
    scale_x_log10() +
    ggtitle("IL8")

# Graph by dichotomized fibrosis
GALA %>%
    filter(!is.na(te_fibrosis)) %>%
    filter(Assay == "IL8") %>%
    ggplot(aes(x = te_fibrosis, y = corr_NPX)) +
    geom_violin(draw_quantiles = c(0.5)) +
    geom_jitter() +
    ggtitle("IL8")

# Run linear mixed model
GALA$te <- as.numeric(GALA$te)

lm_GALA_te <- olink_lmer(GALA,
                         variable = "te",
                         outcome = "corr_NPX",
                         random = "cohort",
                         return.covariates = TRUE)

# Extract list of significant cytokines
sign_cytokines <- lm_GALA_te %>%
    filter(Adjusted_pval < 0.0001) %>%
    pull(Assay)

# Plotting ----------------------------------------------------------------

theme_set(theme_bw())

# histogram
ALCO %>%
    filter(Assay == "IL8") %>%
    ggplot(aes(x = corr_NPX)) +
    geom_density() +
    facet_grid(rows = vars(subgroup))

# boxplot
# time points grouped
### Loop over all cytokines
allcytokines_subgroup <- list()

for(i in cytokines[1:92]) {
    ploti <- ALCO %>%
        filter(Assay == i) %>%
        ggplot(aes(x = subgroup, y = corr_NPX)) +
        geom_boxplot() +
        ggtitle(paste(i))

    allcytokines_subgroup[[i]] <- ploti
    # print(plot1)
}

pdf(here::here("doc/images/ALCO_allcytokines.pdf"), height = 60, width = 40) # height = 5 for each row of 8 -> 60 for 12x8
do.call('grid.arrange', c(allcytokines_subgroup, ncol=8))
dev.off()



### plot significant cytokines
# run ANOVA
anova_ALCO_group <- olink_anova(ALCO, variable = "subgroup")

# Extract list of significant cytokines
sign_cytokines <- anova_ALCO_group %>%
    filter(Adjusted_pval < 0.0001) %>%
    pull(Assay)
# print(sign_cytokines) # 51 cytokines significant

signcytokines_subgroup <- list()

for(i in sign_cytokines[1:51]) {
    ploti <- ALCO %>%
        filter(Assay == i) %>%
        ggplot(aes(x = subgroup, y = corr_NPX)) +
        geom_violin() +
        ggtitle(paste(i))

    signcytokines_subgroup[[i]] <- ploti
    # print(plot1)
}

pdf(here::here("doc/images/ALCO_sign_cytokines_subgroup.pdf"), height = 35, width = 40) # height = 5 for each row of 8 -> 35 for 7x8
do.call('grid.arrange', c(signcytokines_subgroup, ncol=8))
dev.off()


### Individual graphs
# By subgroup
ALCO %>%
    filter(Assay == "IL8") %>%
    ggplot(aes(x = subgroup, y = corr_NPX)) +
    geom_boxplot(outlier.colour = "transparent") +
    geom_jitter() +
    ggtitle("IL8")

# time points separated
ALCO %>%
    filter(Assay == "IL8") %>%
    ggplot(aes(x = subgroup, y = corr_NPX, color = time_point)) +
    geom_boxplot() +
    ggtitle("IL8")

# sample types separated
ALCO %>%
    filter(Assay == "IL8") %>%
    ggplot(aes(x = subgroup, y = corr_NPX, color = sample_type)) +
    geom_boxplot() +
    ggtitle("IL8")

# only baseline
ALCO %>%
    filter(Assay == "IL8") %>%
    filter(time_point == 1) %>%
    ggplot(aes(x = subgroup, y = corr_NPX)) +
    geom_boxplot(outlier.colour = "transparent") +
    geom_jitter() +
    ggtitle("IL8_baseline")

### Top hits for GALA
GALA %>%
    # filter(!is.na(as.numeric(kleiner))) %>%
    filter(Assay == "HGF") %>%
    ggplot(aes(x = kleiner, y = corr_NPX)) +
    geom_violin(draw_quantiles = c(0.5)) +
    geom_jitter() +
    ggtitle("HGF")
ggsave(here::here("doc/images/HGF_kleiner.jpg"), height = 4, width = 5)

GALA %>%
    # filter(!is.na(as.numeric(kleiner))) %>%
    filter(Assay == "IL8") %>%
    ggplot(aes(x = kleiner, y = corr_NPX)) +
    geom_violin(draw_quantiles = c(0.5)) +
    geom_jitter() +
    ggtitle("IL8")
ggsave(here::here("doc/images/IL8_kleiner.jpg"), height = 4, width = 5)

GALA %>%
    filter(Assay == "HGF") %>%
    ggplot(aes(x = as.numeric(te), y = corr_NPX)) +
    geom_point(aes(color = cohort)) +
    geom_smooth(method = lm, color = "black") +
    scale_x_log10() +
    ggtitle("HGF")
ggsave(here::here("doc/images/HGF_te.jpg"), height = 4, width = 5)

GALA %>%
    filter(Assay == "IL8") %>%
    ggplot(aes(x = as.numeric(te), y = corr_NPX)) +
    geom_point(aes(color = cohort)) +
    geom_smooth(method = lm, color = "black") +
    scale_x_log10() +
    ggtitle("IL8")
ggsave(here::here("doc/images/IL8_te.jpg"), height = 4, width = 5)


### Plot using Olink package
# plot all cytokines
#ALCO_subgroup_1 <-
olink_boxplot(ALCO,
              variable = "subgroup",
              olinkid_list = cytokine_IDs[1:16],
              verbose = T,
              number_of_proteins_per_plot = 16)
# ggsave(here::here("doc/images/ALCO_subgroup_1.pdf"),
# ALCO_subgroup_1, width = 10, height = 6) # not working because it is not a ggplot

# plot significant cytokines
anova_ALCO_group <- olink_anova(ALCO, variable = "subgroup")
sign_cytokines <- anova_ALCO_group %>%
    filter(Adjusted_pval < 0.001) %>%
    pull(OlinkID)

olink_boxplot(ALCO,
              variable = "subgroup",
              olinkid_list = sign_cytokines,
              verbose = T,
              number_of_proteins_per_plot = 6) # only prints the last set of plots
