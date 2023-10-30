## Initialization ##
# Set seed
set.seed(2106725034)

# Set working directory
setwd("D:/OneDrive - UNIVERSITAS INDONESIA/Kuliah/Semester 5/[SCST603107] Sains Data Genom/UTS")

## Load the dataset ##
library(data.table)
dtgse <- fread("Breast_GSE45827.csv")

# Change to dataframe
dfgse <- as.data.frame(dtgse)
as.data.frame(table(dfgse$type))

# Randomly select 50% of the genes from the dataframe 
typesample <- subset(dfgse, select = c(samples, type))
rdgse <- subset(dfgse, select = -c(type, samples))

rdgse <- rdgse[, sample(ncol(rdgse), ncol(rdgse) * 0.5)]
rdgse <- rdgse[, order(names(rdgse))]

rdgse <- cbind(typesample, rdgse)

## Extract Expression (modify dataframe) ##
exprdgse <- rdgse
rownames(exprdgse) <- exprdgse$samples
types <- exprdgse$type

exprdgse <- subset(exprdgse, select = -c(type, samples))
exprdgse_t <- transpose(exprdgse)

rownames(exprdgse_t) <- colnames(exprdgse)
colnames(exprdgse_t) <- rownames(exprdgse)

exprdgse_mat <- as.matrix(exprdgse_t)

#### Nomor 1 ####
## Exploration ##
as.data.frame(table(rdgse$type))

library(RColorBrewer)
colour <- brewer.pal(5, "Set2")
barplot(table(rdgse$type),
        main = "Distribution of label types in breast cancer data",
        xlab = "Sample type",
        ylab = "Number of occurences",
        ylim = c(0,50),
        col = colour)

hist(exprdgse_mat, col = 'darkturquoise')

## Gene Filtering ##
library(Biobase)

# Preparing expression set
pData <- data.frame(type = rdgse$type)
rownames(pData) <- rdgse$samples

eset <- ExpressionSet(assayData = exprdgse_mat,
                      phenoData = AnnotatedDataFrame(pData),
                      annotation = 'hgu133plus2')

# Perform gene filtering
require(genefilter)
esetFilt <- nsFilter(eset)

# Filtering result
esetFilt

# Extract the Expression of the Filtered Dataset
exprdgseFilt_mat <- exprs(esetFilt$eset)

# Plot the original and filtered data
par(mfrow = c(1,2))
hist(exprdgse_mat, main = 'original')
hist(exprdgseFilt_mat, main = 'filtered')

# Convert to dataframe
exprdgseFilt <- as.data.frame(exprdgseFilt_mat)

rdgseFilt <- t(exprdgseFilt)
rdgseFilt <- cbind(typesample, rdgseFilt)

#### Nomor 2 ####
## LIMMA Multiple Analysis ##
library(limma)
library(ggplot2)
library(stringr)
library(tibble)
library(dplyr)

# Create dataframe with only cancer genes
crdgse <- rdgseFilt[!(rdgseFilt$type == "cell_line" | rdgseFilt$type == "normal"),]

cexprdgse <- crdgse
rownames(cexprdgse) <- cexprdgse$samples

ctypes <- cexprdgse$type

cexprdgse <- subset(cexprdgse, select = -c(type, samples))
cexprdgse_t <- transpose(cexprdgse)

rownames(cexprdgse_t) <- colnames(cexprdgse)
colnames(cexprdgse_t) <- rownames(cexprdgse)

des_mat <- model.matrix(~ ctypes + 0, data = cexprdgse_t)
colnames(des_mat) <- str_remove(colnames(des_mat), "ctypes")
head(des_mat)

# Apply linear model to data
fit <- lmFit(cexprdgse_t, design = des_mat)

# Apply empirical Bayes to smooth standard errors
fit <- eBayes(fit)

# Create contrast matrix
contrast_matrix <- makeContrasts(
  "basalvsOther" = basal - (HER + luminal_A + luminal_B) / 3,
  "HERvsOther" = HER - (basal + luminal_A + luminal_B) / 3,
  "luminal_AvsOther" = luminal_A - (basal + HER + luminal_B) / 3,
  "luminal_BvsOther" = luminal_B - (basal + HER + luminal_A) / 3,
  levels = des_mat
)

# Fit the model according to the contrasts matrix
contrasts_fit <- contrasts.fit(fit, contrast_matrix)

# Re-smooth the Bayes
contrasts_fit <- eBayes(contrasts_fit)

# Apply multiple testing correction and obtain stats
stats_df <- topTable(contrasts_fit, number = nrow(cexprdgse_t)) %>%
  rownames_to_column("Gene")
head(stats_df)

# Show top genes
top_gene_df <- dplyr::select(crdgse, samples, '216836_s_at', type)
head(top_gene_df)

ggplot(top_gene_df, aes(x = type, y = `216836_s_at`, color = type)) +
  geom_jitter(width = 0.2, height = 0) +
  theme_classic() +
  theme(legend.position = "none")

# Let's extract the contrast p values for each and transform them with -log10()
contrast_p_vals_df <- -log10(contrasts_fit$p.value) %>%
  # Make this into a data frame
  as.data.frame() %>%
  # Store genes as their own column
  tibble::rownames_to_column("Gene") %>%
  # Make this into long format
  tidyr::pivot_longer(dplyr::contains("vsOther"),
                      names_to = "contrast",
                      values_to = "neg_log10_p_val"
  )

# Let's extract the fold changes from `stats_df`
log_fc_df <- stats_df %>%
  # We only want to keep the `Gene` column as well
  dplyr::select("Gene", dplyr::contains("vsOther")) %>%
  # Make this a longer format
  tidyr::pivot_longer(dplyr::contains("vsOther"),
                      names_to = "contrast",
                      values_to = "logFoldChange"
  )

plot_df <- log_fc_df %>%
  dplyr::inner_join(contrast_p_vals_df,
                    by = c("Gene", "contrast"),
                    # This argument will add the given suffixes to the column names
                    # from the respective data frames, helping us keep track of which columns
                    # hold which types of values
                    suffix = c("_log_fc", "_p_val")
  )

# Print out what this looks like
head(plot_df)

# Convert p value cutoff to negative log 10 scale
p_val_cutoff <- -log10(0.05)

# Absolute value cutoff for fold changes
abs_fc_cutoff <- 0.5

plot_df <- plot_df %>%
  dplyr::mutate(
    signif_label = dplyr::case_when(
      abs(logFoldChange) > abs_fc_cutoff & neg_log10_p_val > p_val_cutoff ~ "p-val and FC",
      abs(logFoldChange) > abs_fc_cutoff ~ "FC",
      neg_log10_p_val > p_val_cutoff ~ "p-val",
      TRUE ~ "NS"
    )
  )

volcanoes_plot <- ggplot(
  plot_df,
  aes(
    x = logFoldChange, # Fold change as x value
    y = neg_log10_p_val, # -log10(p value) for the contrasts
    color = signif_label # Color code by significance cutoffs variable we made
  )
) +
  # Make a scatter plot with points that are 30% opaque using `alpha`
  geom_point(alpha = 0.3) +
  # Draw our `p_val_cutoff` for line here
  geom_hline(yintercept = p_val_cutoff, linetype = "dashed") +
  # Using our `abs_fc_cutoff` for our lines here
  geom_vline(xintercept = c(-abs_fc_cutoff, abs_fc_cutoff),
             linetype = "dashed") +
  # The default colors aren't great, we'll specify our own here
  scale_colour_manual(values = c("cadetblue3", "darkgray", "gray", "darkolivegreen3")) +
  # Let's be more specific about what this p value is in our y axis label
  ylab("Contrast -log10(p value)") +
  # This makes separate plots for each contrast!
  facet_wrap(~contrast) +
  # Just for making it prettier!
  theme_classic()

# Print out the plot!
volcanoes_plot

# Select only 50 top genes
topResult <- topTable(contrasts_fit, coef = 2, number = 50)

# Selected Genes
rownames(topResult)

# Extract selected gene names
cexprdgse_mat <- as.matrix(cexprdgse_t)
selected <- rownames(cexprdgse_mat) %in% rownames(topResult)

# Extract the expression of the selected genes
exprdgsesel <- cexprdgse_mat[selected, ]
colnames(exprdgsesel) <- ctypes

## Heatmap of the top genes ##
heatmap(exprdgsesel)

#### Nomor 3 ####
# Assign Group Code
group <- ifelse(types == "normal", 0, 1)

# Apply t-test
t.test(group, exprdgseFilt_mat[1,])$p.value

# Create box plot
boxplot(group, exprdgseFilt_mat[1,])

# Volcano plot
pval <- apply(exprdgseFilt_mat, 1, function(x) t.test(group, x)$p.value)
dift <- apply(exprdgseFilt_mat, 1,
              function(x) diff(t.test(x[1:4], x[5:8])$estimate))
plot(dift, -log10(pval))

library(RColorBrewer)
smoothScatter(dift, -log10(pval))

# Calculate significants
sum(pval < 0.05)

pvalBonf <- p.adjust(pval, method = "bonferroni" )
pvalHolm <- p.adjust(pval, method = "holm" )

sum(pvalBonf < 0.05)
sum(pvalHolm < 0.05)

## LIMMA Analysis ##
design <- model.matrix(~group)

# Apply linear model to data
fit <- lmFit(exprdgseFilt_mat, design)
fit

# Apply empirical Bayes to smooth standard errors
fitted.ebayes <- eBayes(fit)
toptab <- topTable(fitted.ebayes, coef=2, n = Inf)

# Visualize the volcano plot
toptab <- toptab %>%
  mutate(Significant = ifelse(
    abs(logFC) > 1 & adj.P.Val < 0.05,
    TRUE, FALSE))

ggplot(toptab, aes(x = logFC, 
                   y = -log10(P.Value),
                   colour = Significant)) +
  geom_point(size = 1, alpha = 0.6) +
  scale_color_manual(values = c("grey", "red")) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(x = "Log2 Fold Change", y = "-Log10 Adjusted P-value") +
  ggtitle("Volcano plot of differentially expressed genes") +
  theme(plot.title = element_text(hjust = 0.5))

# Show only top 50 genes
topResult <- topTable(fitted.ebayes, coef = 2, number = 50)

# Selected Genes
rownames(topResult)

# Extract selected gene names
selected <- rownames(exprdgseFilt_mat) %in% rownames(topResult)

# Extract the expression of the selected genes
exprdgsesel <- exprdgseFilt_mat[selected, ]
colnames(exprdgsesel) <- types

## Heatmap of the top genes ##
heatmap(exprdgsesel)

## Boxplot for the top 4 genes ##
group2 <- as.factor(group)
levels(group2) <- c('normal', 'cancer')

par(mfrow = c(2,2))
for (i in 1:4) plot(group2, exprdgsesel[i,],
                    main = rownames(exprdgsesel[i]))

## See gene name and description ##
library(annotate)
library(hgu133plus2.db)

GeneSelected <- select(hgu133plus2.db, rownames(topResult),
                       c("SYMBOL", "ENTREZID", "GENENAME"))
GeneSelected

ids <- rownames(topResult)
GeneSelected <- select(hgu133plus2.db, ids,
                       c("SYMBOL", "ENTREZID", "GENENAME", "GO"))

## Gene ontology for the top genes ##
library(GO.db)
GOselected <- select(GO.db, GeneSelected$GO, c("TERM", "GOID"))
head(GOselected)

# Combine the result
finalres <- cbind(GeneSelected, GOselected)

# Convert to csv
write.csv2(finalres, file = "GEOres_UTS Sains Data Genom.csv")
