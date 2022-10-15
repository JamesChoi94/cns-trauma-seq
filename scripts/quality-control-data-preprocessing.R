library(Seurat)
library(ggplot2)
library(dplyr)

dir.create(path = 'results/quality-control-data-processing/')
dir.create(path = 'results/quality-control-data-processing/emptyDrops/')

previous <- Read10X_h5(filename = '../MiamiProject/sci_scRNAseq/data/raw_feature_bc_matrix/raw_feature_bc_matrix_3dpi_sample1.h5')
current <- Read10X_h5(filename = 'data/raw-feature-barcode-matrices/3dpi_sample1-2022_raw_feature_bc_matrix.h5')


# emptyDrops QC -----------------------------------------------------------

# Since v6.0 of cellranger by 10X Genomics, GEM-barcodes with zero total UMIs are automatically removed from the count matrix before being written out to h5.

ranks_previous <- DropletUtils::barcodeRanks(
  m = previous,
  lower = 200,
  fit.bounds = c(500, 3e4)
)
drops_previous <- DropletUtils::emptyDrops(
  m = previous,
  retain = ranks_previous@metadata$knee,
  lower = ranks_previous@metadata$inflection,
  BPPARAM = BiocParallel::SnowParam(workers = 2, type = 'SOCK')
)
previous_results <- cbind(drops_previous, ranks_previous)
rm(ranks_previous, drops_previous)

ranks_current <- DropletUtils::barcodeRanks(
  m = current,
  lower = 200,
  fit.bounds = c(500, 3e4)
)
drops_current <- DropletUtils::emptyDrops(
  m = current,
  retain = ranks_current@metadata$knee,
  lower = ranks_current@metadata$inflection,
  BPPARAM = BiocParallel::SnowParam(workers = 2, type = 'SOCK')
)
current_results <- cbind(drops_current, ranks_current)
rm(ranks_current, drops_current)

all_results <- list(
  previous_results = previous_results,
  current_results = current_results
)
saveRDS(all_results, file = 'results/quality-control-data-processing/emptyDrops/emptyDrops-outs_barcodeRanks-outs.rds')

shuffle_rows <- function(x) return(x[sample(1:nrow(x), nrow(x)),])
all_results$previous_results$Sample <- 'previous'
all_results$current_results$Sample <- 'current'
tmp <- all_results %>% 
  Reduce(f = rbind) %>% 
  as.data.frame() %>% 
  group_by(Sample) %>% 
  filter(!duplicated(rank)) %>% 
  mutate(sig = ifelse(FDR < 0.05, 'cell', 'empty-droplet')) %>% 
  shuffle_rows() %>% 
  arrange(rank) %>% 
  ungroup()
p1 <- tmp %>% 
  ggplot(mapping = aes(x = rank, y = Total)) +
  geom_point(mapping = aes(color = Sample), size = 1) +
  scale_x_continuous(
    trans = 'log10', 
    labels = scales::label_number(scale_cut = scales::cut_short_scale()),
    breaks = 10^seq(0,10,1)
  ) +
  scale_y_continuous(
    trans = 'log10', 
    labels = scales::label_number(scale_cut = scales::cut_short_scale()),
    breaks = 10^seq(0,10,1)
  ) + 
  scale_color_manual(values = c('previous' = 'blue', 'current' = 'red')) +
  # scale_alpha_manual(values = c('cell' = 0.6, 'empty-droplet' = 0.05), 
  #                    na.value = 0.05) +
    theme_bw() +
    theme(axis.title = element_text(size = 12),
          axis.text.x = element_text(angle = 65, hjust = 1, size = 12),
          axis.text.y = element_text(size = 12),
          legend.text = element_text(size = 12))
p2 <- tmp %>% 
  ggplot(mapping = aes(x = rank, y = Total)) +
  geom_point(mapping = aes(color = sig), size = 1, alpha = 0.4) +
  facet_grid(. ~ Sample) +
  scale_x_continuous(
    trans = 'log10', 
    labels = scales::label_number(scale_cut = scales::cut_short_scale()),
    breaks = 10^seq(0,10,1)
  ) +
  scale_y_continuous(
    trans = 'log10', 
    labels = scales::label_number(scale_cut = scales::cut_short_scale()),
    breaks = 10^seq(0,10,1)
  ) + 
  scale_color_manual(values = c('cell' = 'red', 'empty-droplet' = 'black')) +
    theme_bw() +
    theme(axis.title = element_text(size = 12),
          axis.text.x = element_text(angle = 65, hjust = 1, size = 12),
          axis.text.y = element_text(size = 12),
          legend.text = element_text(size = 12))
p <- (p1 | p2) + patchwork::plot_layout(widths = c(0.5,1))
ggsave(filename = 'results/quality-control-data-processing/barcode-rank-plot_3dpi-sample1-comparison-joint.tiff', plot = p1, height = 3, width = 4.75)
ggsave(filename = 'results/quality-control-data-processing/barcode-rank-plot_3dpi-sample1-comparison-separate.tiff', plot = p2, height = 3, width = 6)
ggsave(filename = 'results/quality-control-data-processing/barcode-rank-plot_3dpi-sample1-comparison.tiff', plot = p, height = 3, width = 10.5)



# Gene-level QC -----------------------------------------------------------

# params
all_results <- readRDS(file = 'results/quality-control-data-processing/emptyDrops/emptyDrops-outs_barcodeRanks-outs.rds')
max_mt <- 25

# run
previous <- previous[, which(all_results$previous_results$FDR < 0.05)]
current <- current[, which(all_results$current_results$FDR < 0.05)]
colnames(previous) <- paste(colnames(previous), 'Milich2018', sep = '.')
colnames(current) <- paste(colnames(current), 'Choi2022', sep = '.')
dim(current)
dim(previous)

previous <- CreateSeuratObject(counts = previous, project = 'previous')
current <- CreateSeuratObject(counts = current, project = 'current')

previous <- PercentageFeatureSet(previous, pattern = '^mt-', col.name = 'percent_mt')
current <- PercentageFeatureSet(current, pattern = '^mt-', col.name = 'percent_mt')
previous$log.nCount_RNA <- log10(previous$nCount_RNA)
current$log.nCount_RNA <- log10(current$nCount_RNA)

p.meta.unfiltered <- rbind(previous@meta.data, current@meta.data) %>% 
  reshape2::melt(id.vars = c('orig.ident')) %>% 
  ggplot(mapping = aes(x = orig.ident, y = value)) + 
  geom_violin(scale = 'width') + 
  facet_wrap(. ~ variable, scales = 'free_y') +
  theme_bw() +
  labs(title = 'cellranger count v2.0 vs v6.1', subtitle = 'previous = Milich2019, current = Choi2022;\n(unfiltered)')
ggsave(filename = 'results/quality-control-data-processing/gene-level-qc-vln_3dpi-sample1-comparison_unfiltered.tiff', plot = p.meta.unfiltered, device = 'tiff', height = 4.5, width = 3.75, dpi = 320)

computeMAD <- function(x) {
  x.med <- median(x)
  x.mad <- median(abs(x - x.med))
  return(x.mad)
}
previous.umi <- computeMAD(previous$log.nCount_RNA) * 3
current.umi <- computeMAD(current$log.nCount_RNA) * 3

previous.keep <- which((abs(previous$log.nCount_RNA - median(previous$log.nCount_RNA)) < previous.umi) & (previous$percent_mt < max_mt))
current.keep <- which((abs(current$log.nCount_RNA - median(current$log.nCount_RNA)) < current.umi) & (current$percent_mt < max_mt))

previous <- previous[, previous.keep]
current <- current[, current.keep]

p.meta.filtered <- rbind(previous@meta.data, current@meta.data) %>% 
  reshape2::melt(id.vars = c('orig.ident')) %>% 
  ggplot(mapping = aes(x = orig.ident, y = value)) + 
  geom_violin(scale = 'width') + 
  facet_wrap(. ~ variable, scales = 'free_y') +
  theme_bw() +
  labs(title = 'cellranger count v2.0 vs v6.1', subtitle = 'previous = Milich2019, current = Choi2022\n(filtered)')
ggsave(filename = 'results/quality-control-data-processing/gene-level-qc-vln_3dpi-sample1-comparison_filtered.tiff', plot = p.meta.filtered, device = 'tiff', height = 4.5, width = 3.75, dpi = 320)

merged <- merge(current, previous) %>% 
  NormalizeData() %>% 
  FindVariableFeatures(nfeatures = 2500) %>% 
  ScaleData() %>% 
  RunPCA() %>% 
  FindNeighbors(dims = 1:10) %>% 
  RunUMAP(dims = 1:10)

DimPlot(merged)
merged <- list(current, previous)
merged <- lapply(merged, NormalizeData)
merged <- lapply(merged, FindVariableFeatures, nfeatures = 2500)
feats <- SelectIntegrationFeatures(merged, nfeatures = 2500)
anchors <- FindIntegrationAnchors(merged, anchor.features = feats)
