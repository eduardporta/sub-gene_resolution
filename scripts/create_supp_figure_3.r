library (ggplot2)
library (pheatmap)
library (reshape2)
library (gridExtra)
library (gplots)
library (ggdendro)
library (grid)

#THIS IS A FUNCTION THAT EXTRACTS THE LEGEND FROM A GGPLOT
#TOOK IT FROM https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}
  
#LOAD DATA
pipeline.res <- read.table ("../raw_data/matrix_brca.txt", sep = "\t", header = TRUE, check.names = FALSE)

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
pipe.melt$binarized <- ifelse (pipe.melt$value < 0.05, "q < 0.05", ifelse (pipe.melt$value < 0.1, "0.05 < q < 0.1", "q > 0.1"))
pipe.melt$binarized <- factor (pipe.melt$binarized, levels = c ("q < 0.05", "0.05 < q < 0.1", "q > 0.1"))

#RE-ORDER THE METHODS
methods.sorted <- colnames (mat.plot)[model.meth$order]
methods.sorted[16] <- "NPat"
methods.sorted[17] <- "Detected"
methods.sorted[18] <- "Coverage"
pipe.melt$variable <- factor (pipe.melt$variable, levels = methods.sorted)

#PREPARE THE PLOTS
heat.plot <- ggplot (subset (pipe.melt, !(variable %in% c("NPat", "Detected", "Coverage"))), aes (x = variable, y = Gene)) + geom_tile (aes (fill = binarized), color = "black") + theme (axis.text = element_text (color = "black"), axis.text.x = element_text (angle = 45, hjust = 1), axis.ticks = element_blank(), axis.line = element_blank(), panel.background = element_blank(), axis.title = element_blank(), legend.position = "none") + scale_fill_manual (name = "Detected", values = c("#FB6542", "#FFEC5C", "white"))
leg.plot.heat <- ggplot (subset (pipe.melt, !(variable %in% c("NPat", "Detected", "Coverage"))), aes (x = variable, y = Gene)) + geom_tile (aes (fill = binarized), color = "black") + theme (axis.text = element_text (color = "black"), axis.text.x = element_text (angle = 45, hjust = 1), axis.ticks = element_blank(), axis.line = element_blank(), panel.background = element_blank(), axis.title = element_blank(), legend.position = "right") + scale_fill_manual (guide = guide_legend (title = "Detected", title.position = "top", title.hjust = 0.5, title.theme = element_text (face = "italic", angle = 0, size = 10), ncol = 1), values = c("#FB6542", "#FFEC5C", "white"))
leg.heat <- g_legend (leg.plot.heat)
npat.plot <- ggplot (subset (pipe.npat, variable == "NPat"), aes (x = Gene, y = value)) + geom_bar (stat = "identity") + coord_flip() + theme (axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title = element_blank(), axis.text.x = element_text (color = "black"), axis.ticks.x = element_line (color = "black"), axis.line.x = element_line (color = "black"), panel.background = element_blank())
known.plot <- ggplot (pipe.melt, aes (x = "Known driver", y = Gene)) + geom_tile (aes (fill = Mode), color = "black") + theme (axis.text.y = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), axis.text.x = element_text (color = "black", angle = 45, hjust = 1), legend.position = "none", panel.background = element_blank()) + scale_fill_manual (values = c ("#FFBB00", "#375E97", "#3F681C", "#DDC5A2", "white"))
leg.plot.known <- ggplot (pipe.melt, aes (x = "Known driver", y = Gene)) + geom_tile (aes (fill = Mode), color = "black") + theme (axis.text.y = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), axis.text.x = element_text (color = "black", angle = 45, hjust = 1), legend.position = "right", panel.background = element_blank()) + scale_fill_manual (name = "Mode", values = c ("#FFBB00", "#375E97", "#3F681C", "#DDC5A2", "white"), guide = guide_legend (ncol = 2, title = "Mode of action", title.position = "top", title.hjust = 0.5, title.theme = element_text (face = "italic", angle = 0, size = 10)))
leg.known <- g_legend (leg.plot.known)
struc.plot <- ggplot (subset (pipe.melt, variable == "Coverage"), aes (x = "Structural coverage", y = Gene, fill = value)) + geom_tile (color = "black") + theme (axis.text.y = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), axis.text.x = element_text (color = "black", angle = 45, hjust = 1), legend.position = "none", panel.background = element_blank()) + scale_fill_continuous (low = "white", high = "#4CB5F5")
leg.plot.struc <- ggplot (subset (pipe.melt, variable == "Coverage"), aes (x = "Structural coverage", y = Gene, fill = value)) + geom_tile (color = "black") + theme (axis.text.y = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), axis.text.x = element_text (color = "black", angle = 45, hjust = 1), legend.position = "bottom", panel.background = element_blank()) + scale_fill_gradient (low = "white", high = "#4CB5F5", limits = c (0,1), labels = scales::percent, guide = guide_legend (title = "Structural coverage", title.position = "top", title.hjust = 0.5, label.position = "bottom", label.hjust = 0.5, title.theme = element_text (face = "italic", angle = 0, size = 10)))
leg.struc <- g_legend (leg.plot.struc)
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
pipe.back.melt$binarized <- ifelse (pipe.back.melt$value < 0.05, "q < 0.05", ifelse (pipe.back.melt$value < 0.1, "0.05 < q < 0.1", "q > 0.1"))
pipe.back.melt$binarized <- factor (pipe.back.melt$binarized, levels = c ("q < 0.05", "0.05 < q < 0.1", "q > 0.1"))

#ORDER CATEGORIES
types.sorted <- colnames (mat.back.plot)[model.type$order]
types.sorted[6] <- "NPat"
types.sorted[7] <- "Detected"
types.sorted[8] <- "Coverage"
pipe.back.melt$variable <- factor (pipe.back.melt$variable, levels = types.sorted)

#CREATE THE PLOT FOR THE CATEGORIES
type.plot <- ggplot (subset (pipe.back.melt, !(variable %in% c("NPat", "Detected", "Coverage"))), aes (x = variable, y = Gene)) + geom_tile (aes (fill = binarized), color = "black") + theme (axis.text = element_text (color = "black"), axis.text.x = element_text (angle = 45, hjust = 1), axis.ticks = element_blank(), axis.line = element_blank(), panel.background = element_blank(), axis.title = element_blank(), legend.position = "none") + scale_fill_manual (name = "Detected", values = c("#FB6542", "#FFEC5C", "white")) + theme (axis.text.y = element_blank())

#ADJUST HEIGHTS OF THE METHODS
gA=ggplot_gtable(ggplot_build(heat.plot))
gB=ggplot_gtable(ggplot_build(type.plot))
gC2=ggplot_gtable(ggplot_build(struc.plot))
gC=ggplot_gtable(ggplot_build(known.plot))
gD=ggplot_gtable(ggplot_build(npat.plot))
gE=ggplot_gtable(ggplot_build(dendro.gene))
maxHeight = grid::unit.pmax(gA$heights, gB$heights, gC$heights, gD$heights, gE$heights)
gA$heights <- as.list(maxHeight)
gB$heights <- as.list(maxHeight)
gC2$heights <- as.list(maxHeight)
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
grid.arrange (arrangeGrob(empty.plot, gF, gG, empty.plot, empty.plot, empty.plot, gE, gA, gB, gC2, gC, gD, nrow = 2, widths = c (1/10, 4/5, 3/10, 2/25, 2/25, 2/5), heights = c (1, 12)), empty.plot, arrangeGrob(empty.plot, leg.heat, leg.known, leg.struc, empty.plot, nrow = 1, widths = c (7/2,2,4,2,1/4)), ncol = 1, heights = c (14,1/4,2))