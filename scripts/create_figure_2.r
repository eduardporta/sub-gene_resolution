library (ggplot2)
library (pheatmap)
library (reshape2)
library (gridExtra)
library (gplots)
library (ggdendro)
library (grid)
library (ggrepel)

#WE CREATE THE PCA PLOTS WITH THE DIFFERENT METHODS, PART "A" OF THE FIGURE

#THIS IS A FUNCTION THAT EXTRACTS THE LEGEND FROM A GGPLOT
#TOOK IT FROM https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

#FIRST WE GET THE PCA DATA

#FOR BLCA
blca.data <- read.table ("../raw_data/matrix_pca_blca.txt", sep = "\t", check.names = FALSE, header = TRUE, row.names = 1)

#SUBSET THE COLUMNS WE WANT
blca.bk <- blca.data[,2:ncol(blca.data)]
blca.bk <- as.matrix (blca.bk)

#ADJUST THE LIMITS OF THE P-VALUE SO WE DON'T HAVE ANY 0s
blca.bk <- ifelse (blca.bk < 1e-20, 1e-20, blca.bk)

#CALC PCA AND CREATE DATA FRAME WITH THE RESULT
blca.prcomp <- prcomp (-log10(blca.bk), center = TRUE, scale = TRUE)
blca.pca <- cbind (blca.data, blca.prcomp$x)

#STORE THE PLOT
blca.plot <- ggplot (blca.pca, aes (x = PC1, y = PC2)) + geom_point (aes (color = Type)) + geom_text_repel (aes (label = row.names (blca.pca), color = Type), max.iter = 50000, size = 3, segment.alpha = 0.2)  + theme (axis.line.x = element_line (color = "black"), axis.line.y = element_line (color = "black"), axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "none", panel.background = element_blank()) + ggtitle ("BLCA") + theme (plot.title = element_text (face = "bold.italic", hjust = 0.5)) + xlim (min (blca.pca$PC1) - 1, max (blca.pca$PC1) + 10) + ylim (min (blca.pca$PC2) - 1, max (blca.pca$PC2) + 1) + scale_color_manual (values = c ("red", "orange", "#598234", "#4cb5f5", "black"))

#STORE THE LEGEND
blca.legend <- ggplot (blca.pca, aes (x = PC1, y = PC2)) + geom_point (aes (color = Type)) + theme (axis.line.x = element_line (color = "black"), axis.line.y = element_line (color = "black"), axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "bottom", panel.background = element_blank()) + ggtitle ("BLCA") + theme (plot.title = element_text (face = "bold.italic", hjust = 0.5)) + scale_color_manual (values = c ("red", "orange", "#598234", "#4cb5f5", "black"))
plot.legend <- g_legend (blca.legend)

#SAME FOR BRCA
brca.data <- read.table ("../raw_data/matrix_pca_brca.txt", sep = "\t", check.names = FALSE, header = TRUE, row.names = 1)
brca.bk <- brca.data[,2:ncol(brca.data)]
brca.bk <- as.matrix (brca.bk)
brca.bk <- ifelse (brca.bk < 1e-20, 1e-20, brca.bk)
brca.prcomp <- prcomp (-log10(brca.bk), center = TRUE, scale = TRUE)
brca.pca <- cbind (brca.data, brca.prcomp$x)
brca.plot <- ggplot (brca.pca, aes (x = PC1, y = PC2)) + geom_point (aes (color = Type)) + geom_text_repel (aes (label = row.names (brca.pca), color = Type), max.iter = 50000, size = 3, segment.alpha = 0.2)  + theme (axis.line.x = element_line (color = "black"), axis.line.y = element_line (color = "black"), axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "none", panel.background = element_blank()) + ggtitle ("BRCA") + theme (plot.title = element_text (face = "bold.italic", hjust = 0.5)) + xlim (min (brca.pca$PC1) - 1, max (brca.pca$PC1) + 10) + ylim (min (brca.pca$PC2) - 1, max (brca.pca$PC2) + 10) + scale_color_manual (values = c ("red", "orange", "#598234", "#4cb5f5", "black"))

#SAME FOR GBM
gbm.data <- read.table ("../raw_data/matrix_pca_gbm.txt", sep = "\t", check.names = FALSE, header = TRUE, row.names = 1)
gbm.bk <- gbm.data[,2:ncol(gbm.data)]
gbm.bk <- as.matrix (gbm.bk)
gbm.bk <- ifelse (gbm.bk < 1e-20, 1e-20, gbm.bk)
gbm.prcomp <- prcomp (-log10(gbm.bk), center = TRUE, scale = TRUE)
gbm.pca <- cbind (gbm.data, gbm.prcomp$x)
gbm.plot <- ggplot (gbm.pca, aes (x = PC1, y = PC2)) + geom_point (aes (color = Type)) + geom_text_repel (aes (label = row.names (gbm.pca), color = Type), max.iter = 50000, size = 3, segment.alpha = 0.2)  + theme (axis.line.x = element_line (color = "black"), axis.line.y = element_line (color = "black"), axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "none", panel.background = element_blank()) + ggtitle ("GBM") + theme (plot.title = element_text (face = "bold.italic", hjust = 0.5)) + xlim (min (gbm.pca$PC1) - 1, max (gbm.pca$PC1) + 1) + ylim (min (gbm.pca$PC2) - 1, max (gbm.pca$PC2) + 1) + scale_color_manual (values = c ("red", "orange", "#598234", "#4cb5f5", "black"))

#SAME FOR LUAD
luad.data <- read.table ("../raw_data/matrix_pca_luad.txt", sep = "\t", check.names = FALSE, header = TRUE, row.names = 1)
luad.bk <- luad.data[,2:ncol(luad.data)]
luad.bk <- as.matrix (luad.bk)
luad.bk <- ifelse (luad.bk < 1e-20, 1e-20, luad.bk)
luad.prcomp <- prcomp (-log10(luad.bk), center = TRUE, scale = TRUE)
luad.pca <- cbind (luad.data, luad.prcomp$x)
luad.plot <- ggplot (luad.pca, aes (x = PC1, y = PC2)) + geom_point (aes (color = Type)) + geom_text_repel (aes (label = row.names (luad.pca), color = Type), max.iter = 50000, size = 3, segment.alpha = 0.2)  + theme (axis.line.x = element_line (color = "black"), axis.line.y = element_line (color = "black"), axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "none", panel.background = element_blank()) + ggtitle ("LUAD") + theme (plot.title = element_text (face = "bold.italic", hjust = 0.5)) + xlim (min (luad.pca$PC1) - 1, max (luad.pca$PC1) + 20) + ylim (min (luad.pca$PC2) - 10, max (luad.pca$PC2) + 1) + scale_color_manual (values = c ("red", "orange", "#598234", "#4cb5f5", "black"))


#NOW WE GO FOR PART "B"

pipeline.res <- read.table ("../raw_data/matrix_gbm.txt", sep = "\t", header = TRUE, check.names = FALSE)

#THIS ADJUSTS THE ORDER OF THE "MODE OF ACTION" LEGEND
pipeline.res$Mode <- factor (pipeline.res$Mode, levels = c ("OG", "TSG", "OG/TSG", "Unk", "N/A"))

#WE SUBSET TO THE FRACTION OF GENES WE WANT
#I LIMITED THE PLOT TO GENES DETECTED BY FIVE OR MORE METHODS, OR KNOWN DRIVER GENES IN THAT TISSUE DETECTED AT LEAST ONCE
#FEEL FREE TO ADJUST TO ADD OTHER GENES
pipe.sig <- subset (pipeline.res, Detected > 3 | (Detected > 0 & KnownDriver != "No"))

#SUBSET ONLY THE PART WE WANT
mat.plot <- pipe.sig[,c ("Gene", "MutSigCV", "OncodriveFM", "OncodriveCLUST", "NMC", "MutSig-CL", "e-Driver", "e-Driver3D", "iPAC", "GraphPAC", "SpacePAC", "CLUMPS", "iSIMPRe", "ActiveDriver", "Hotspot", "LowMACA")]
mat.melt <- melt (mat.plot)

#BINARIZE THE PREDICTIONS INSTEAD OF P-VALS
mat.melt$value <- ifelse (mat.melt$value < 0.05, 1, 0)

#RE-CREATE THE MATRIX
mat.plot <- acast (mat.melt, Gene ~ variable, value.var = "value")

#CALC THE CLUSTERING FOR THE METHODS AND GET THE DENDROGRAM
model.meth <- hclust (dist (t(mat.plot)))
dhc.meth <- as.dendrogram (model.meth)
ddata.meth <- dendro_data (dhc.meth, type = "rectangle")
dendro.meth <- ggplot (segment (ddata.meth)) + geom_segment (aes (x = x, y = y, xend = xend, yend = yend)) + theme (panel.background = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())

#CALC THE CLUSTERING FOR THE GENES AND GET THE DENDROGRAM
model.gene <- hclust (dist (mat.plot))
dhc.gene <- as.dendrogram (model.gene)
ddata.gene <- dendro_data (dhc.gene, type = "rectangle")
dendro.gene <- ggplot (segment (ddata.gene)) + geom_segment (aes (x = -x, y = -y, xend = -xend, yend = -yend)) + theme (panel.background = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank()) + coord_flip ()

#RE-ORDER THE GENES
pipe.sig$Gene <- factor (pipe.sig$Gene, levels = rev(rownames (mat.plot)[model.gene$order]))
pipe.melt <- melt (pipe.sig)
pipe.npat <- pipe.melt

#BINARIZE THE PREDICTIONS INSTEAD OF P-VALS
pipe.melt$value <- ifelse (pipe.melt$value < 0.05, "Yes", "No")

#RE-ORDER THE METHODS
methods.sorted <- colnames (mat.plot)[model.meth$order]
methods.sorted[16] <- "NPat"
methods.sorted[17] <- "Detected"
pipe.melt$variable <- factor (pipe.melt$variable, levels = methods.sorted)

#PREPARE THE PLOTS
heat.plot <- ggplot (subset (pipe.melt, !(variable %in% c("NPat", "Detected"))), aes (x = variable, y = Gene)) + geom_tile (aes (fill = value), color = "black") + theme (axis.text = element_text (color = "black"), axis.text.x = element_text (angle = 45, hjust = 1), axis.ticks = element_blank(), axis.line = element_blank(), panel.background = element_blank(), axis.title = element_blank(), legend.position = "none") + scale_fill_manual (name = "Detected", values = c("white", "#FB6542"))
npat.plot <- ggplot (subset (pipe.npat, variable == "NPat"), aes (x = Gene, y = value)) + geom_bar (stat = "identity") + coord_flip() + theme (axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title = element_blank(), axis.text.x = element_text (color = "black"), axis.ticks.x = element_line (color = "black"), axis.line.x = element_line (color = "black"), panel.background = element_blank())
known.plot <- ggplot (pipe.melt, aes (x = "Known driver", y = Gene)) + geom_tile (aes (fill = Mode), color = "black") + theme (axis.text.y = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), axis.text.x = element_text (color = "black", angle = 45, hjust = 1), legend.position = "none", panel.background = element_blank()) + scale_fill_manual (values = c ("#FFBB00", "#375E97", "#3F681C", "#DDC5A2", "white"))
empty.plot <- rectGrob(gp=gpar(col="white"))

#CREATE THE SUMMARY PER CATEGORIES
pipeline.back <- pipeline.res
pipeline.back$Whole <- apply (pipeline.res[,c("OncodriveFM", "MutSigCV")], 1, min)
pipeline.back$TypeI <- apply (pipeline.res[,c("NMC", "OncodriveCLUST", "MutSig-CL", "Hotspot", "iSIMPRe")], 1, min)
pipeline.back$TypeII <- apply (pipeline.res[,c("GraphPAC", "iPAC", "SpacePAC", "CLUMPS")], 1, min)
pipeline.back$TypeIII <- apply (pipeline.res[,c("ActiveDriver", "e-Driver", "LowMACA")], 1, min)
pipeline.back$TypeIV <- pipeline.res$'e-Driver3D'

pipe.back.sig <- subset (pipeline.back, Detected > 3 | (Detected > 0 & KnownDriver != "No"))

mat.back.plot <- pipe.back.sig[,c ("Gene", "Whole", "TypeI", "TypeII", "TypeIII", "TypeIV")]
mat.back.melt <- melt (mat.back.plot)
mat.back.melt$value <- ifelse (mat.back.melt$value < 0.05, 1, 0)
mat.back.plot <- acast (mat.back.melt, Gene ~ variable, value.var = "value")

#CALC DENDROGRAM FROM CATEGORIES
model.type <- hclust (dist (t(mat.back.plot)))
dhc.type <- as.dendrogram (model.type)
ddata.type <- dendro_data (dhc.type, type = "rectangle")
dendro.type <- ggplot (segment (ddata.type)) + geom_segment (aes (x = x, y = y, xend = xend, yend = yend)) + theme (panel.background = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())

pipe.back.sig$Gene <- factor (pipe.back.sig$Gene, levels = rev(rownames (mat.plot)[model.gene$order]))
pipe.back.melt <- melt (pipe.back.sig[,c ("Gene", "Whole", "TypeI", "TypeII", "TypeIII", "TypeIV", "NPat", "Detected", "KnownDriver")])
pipe.back.melt$value <- ifelse (pipe.back.melt$value < 0.05, "Yes", "No")

#ORDER CATEGORIES
types.sorted <- colnames (mat.back.plot)[model.type$order]
types.sorted[6] <- "NPat"
types.sorted[7] <- "Detected"
pipe.back.melt$variable <- factor (pipe.back.melt$variable, levels = types.sorted)

#CREATE THE PLOT FOR THE CATEGORIES
type.plot <- ggplot (subset (pipe.back.melt, !(variable %in% c("NPat", "Detected"))), aes (x = variable, y = Gene)) + geom_tile (aes (fill = value), color = "black") + theme (axis.text = element_text (color = "black"), axis.text.x = element_text (angle = 45, hjust = 1), axis.ticks = element_blank(), axis.line = element_blank(), panel.background = element_blank(), axis.title = element_blank(), legend.position = "none") + scale_fill_manual (name = "Detected", values = c("white", "#FB6542")) + theme (axis.text.y = element_blank())

#ADJUST HEIGHTS OF THE METHODS
gA=ggplot_gtable(ggplot_build(heat.plot))
gB=ggplot_gtable(ggplot_build(type.plot))
gC=ggplot_gtable(ggplot_build(known.plot))
gD=ggplot_gtable(ggplot_build(npat.plot))
gE=ggplot_gtable(ggplot_build(dendro.gene))
maxHeight = grid::unit.pmax(gA$heights, gB$heights, gC$heights, gD$heights, gE$heights)
gA$heights <- as.list(maxHeight)
gB$heights <- as.list(maxHeight)
gC$heights <- as.list(maxHeight)
gD$heights <- as.list(maxHeight)
gE$heights <- as.list(maxHeight)

gF=ggplot_gtable(ggplot_build(dendro.meth))
gG=ggplot_gtable(ggplot_build(dendro.type))

maxWidth = grid::unit.pmax(gA$widths, gF$widths)
gA$widths <- as.list(maxWidth)
gF$widths <- as.list(maxWidth)

maxWidth = grid::unit.pmax(gB$widths, gG$widths)
gB$widths <- as.list(maxWidth)
gG$widths <- as.list(maxWidth)

#ARRANGE IT ALL
grid.arrange (arrangeGrob (arrangeGrob(blca.plot, brca.plot, gbm.plot, luad.plot, nrow = 1), plot.legend, ncol = 1, heights = c (10,1)), empty.plot, arrangeGrob(empty.plot, gF, gG, empty.plot, empty.plot, gE, gA, gB, gC, gD, nrow = 2, widths = c (1/10, 4/5, 3/10, 2/25, 2/5), heights = c (1/10, 9/10)), heights = c (6,1,10))