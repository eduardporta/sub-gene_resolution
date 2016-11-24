require (ggplot2)
require (reshape)
require (gridExtra)

egfr.data <- read.table ("../raw_data/matrix_EGFR_gbm.txt", sep = "\t", header = TRUE, check.names = FALSE)
egfr.data$Sample <- factor (egfr.data$Sample, levels = egfr.data[order(-egfr.data$Gene, -egfr.data$Hotspot, -egfr.data$'e-Driver', -egfr.data$'e-Driver3D', -egfr.data$OncodriveCLUST, -egfr.data$ActiveDriver, -egfr.data$iPAC),]$Sample)

methods <- c ("e-Driver3D", "e-Driver", "ActiveDriver", "CLUMPS", "SpacePAC", "GraphPAC", "SpacePAC", "NMC", "OncodriveCLUST", "Hotspot", "Gene")
egfr.melt <- melt (egfr.data)

sub.methods <- subset (egfr.melt, variable %in% methods)
sub.methods$variable <- factor (sub.methods$variable, levels = methods)
sub.methods$value <- factor (sub.methods$value, levels = c (2,1,0), labels = c("EGFR interface\nmutations", "Other EGFR\nmutations", "No EGFR\nmutations"))
plot.meth <- ggplot (sub.methods, aes (x = Sample, y = variable)) + geom_tile (aes (fill = as.factor(value)), color = "black") + theme (axis.text.x = element_blank(), axis.ticks.x =element_blank()) + scale_fill_manual (name = "", values = c ("#F62A00", "#4CB5F5", "white" )) + theme (legend.position = "right", axis.text.y = element_text (color = "black"), axis.line.y = element_line (color = "black"), axis.ticks.y = element_line (color = "black"), axis.title = element_blank())

mut.analysis <- c ("Sample", "SIFT", "Polyphen", "MutAssessor", "e-Driver3D")
sub.mut.analysis <- egfr.data[,mut.analysis]
sub.mut.analysis <- subset (sub.mut.analysis, SIFT != "NA")
sub.mut.analysis$'e-Driver3D' <- factor (sub.mut.analysis$'e-Driver3D', levels = c (2,1,0), labels = c("EGFR interface\nmutations", "Other EGFR\nmutations", "No EGFR\nmutations"))
sub.mut.melt <- melt (sub.mut.analysis)

prot.analysis <- c ("Sample", "e-Driver3D", "EGFR-R-C", "EGFR_pY1068-R-V", "EGFR_pY1173-R-C", "EGFR_pY992-R-V")
sub.prot.analysis <- egfr.data[,prot.analysis]
sub.prot.analysis <- subset (sub.prot.analysis, "EGFR-R-C" != "NA")
sub.prot.analysis$'e-Driver3D' <- factor (sub.prot.analysis$'e-Driver3D', levels = c (2,1,0), labels = c("EGFR interface\nmutations", "Other EGFR\nmutations", "No EGFR\nmutations"))
sub.prot.melt <- melt (sub.prot.analysis)
sub.prot.melt$variable <- factor (sub.prot.melt$variable, levels = c ("EGFR-R-C", "EGFR_pY992-R-V", "EGFR_pY1068-R-V", "EGFR_pY1173-R-C"), labels = c ("EGFR", "EGFR pY992", "EGFR pY1068", "EGFR pY1173"))


#RED-ORANGE
plot.meth <- ggplot (sub.methods, aes (x = Sample, y = variable)) + geom_tile (aes (fill = as.factor(value)), color = "white") + theme (axis.text.x = element_blank(), axis.ticks.x =element_blank()) + scale_fill_manual (name = "", values = c ("#FE0038", "#FFB74C", "lightgray")) + theme (legend.position = "none", axis.text.y = element_text (color = "black"), axis.line.y = element_line (color = "black"), axis.ticks.y = element_line (color = "black"), axis.title = element_blank())
plot.prot <- ggplot (sub.prot.melt, aes (x = sub.prot.melt$'e-Driver3D', y = value)) + geom_boxplot (outlier.shape = NA) + facet_wrap (~variable, scales = "free", nrow = 1) + geom_point (position = position_jitter (width = 0.2), aes (color = sub.prot.melt$'e-Driver3D')) + theme (axis.text = element_text (color = "black"), axis.ticks = element_line (color = "black"), axis.line.x = element_line (color = "black"), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.y = element_line (color = "black"), panel.background = element_blank(), axis.title.x = element_blank(), axis.title.y = element_text (face = "italic"), strip.background = element_blank(), strip.text = element_text (face = "italic")) + ylab ("Protein level (RPPA)\n") + scale_color_manual (name = "", values = c ("#FE0038", "#FFB74C", "lightgray")) + theme (legend.position = "bottom")
grid.arrange (plot.meth, plot.prot, ncol = 1, heights = c (2, 4))


