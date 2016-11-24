library (ggplot2)
library (reshape)
library (gridExtra)

#THIS IS PART "E" OF THE PLOT

#THIS IS A FUNCTION THAT EXTRACTS THE LEGEND FROM A GGPLOT
#TOOK IT FROM https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

#BLCA

pred.data <- read.table (file = "../raw_data/matrix_blca.txt", sep = "\t", header = TRUE, check.names = FALSE)
pred.melt <- melt (subset (pred.data, NPat > 4))
sub.analyze <- subset (pred.melt, !(variable %in% c ("NPat", "Group", "Detected")))

sub.analyze$AllDriver <- factor (sub.analyze$AllDriver, levels = c ("Missense somatic drivers in matching tissue",  "Missense somatic drivers in other tissues", "Other somatic drivers in other tissues", "Other somatic drivers in matching tissue", "Other germline drivers in other tissues", "No"))
sub.analyze$variable <- factor (sub.analyze$variable, levels = c ("MutSigCV", "OncodriveFM", "Hotspot", "OncodriveCLUST", "NMC", "MutSig-CL", "iSIMPRe", "iPAC", "GraphPAC", "SpacePAC", "CLUMPS", "e-Driver", "ActiveDriver", "LowMACA", "e-Driver3D"))

sub.blca.plot <- subset (sub.analyze, value < 0.05)
blca.plot <- ggplot (sub.blca.plot, aes (x = variable)) + geom_bar (position = "fill", aes (fill = AllDriver), color = "black") + theme (axis.text = element_text (color = "black"), axis.ticks.y = element_line (color = "black"), axis.ticks.x = element_blank(), axis.line.x = element_line (color = "black"), axis.line.y = element_line (color = "black"), panel.background = element_blank(), axis.text.x = element_blank(), legend.position = "none", axis.title = element_blank()) + scale_y_continuous (labels = scales::percent) + scale_fill_manual (values = c ("#FE0000", "orange","#2F496E", "#2988BC", "#5CC5EF", "white"), name = "Known driver role?") + ggtitle ("BLCA") + theme (plot.title = element_text (hjust = 0.5, size = 12, face = "italic"))

blca.legend <- ggplot (sub.blca.plot, aes (x = variable)) + geom_bar (position = "fill", aes (fill = AllDriver), color = "black") + theme (axis.text = element_text (color = "black"), axis.ticks = element_line (color = "black"), axis.line.x = element_line (color = "black"), axis.line.y = element_line (color = "black"), panel.background = element_blank(), axis.text.x = element_text (angle = 45, hjust = 1), legend.position = "right", axis.title = element_blank()) + scale_y_continuous (labels = scales::percent) + scale_fill_manual (values = c ("#FE0000", "orange","#2F496E", "#2988BC", "#5CC5EF", "white"), name = "Known driver role?") + ggtitle ("BLCA") + theme (plot.title = element_text (hjust = 0.5, size = 12, face = "italic"))
plot.legend <- g_legend (blca.legend)

#NOW BRCA

pred.data <- read.table (file = "../raw_data/matrix_brca.txt", sep = "\t", header = TRUE, check.names = FALSE)
pred.melt <- melt (subset (pred.data, NPat > 4))
sub.analyze <- subset (pred.melt, !(variable %in% c ("NPat", "Group", "Detected")))

sub.analyze$AllDriver <- factor (sub.analyze$AllDriver, levels = c ("Missense somatic drivers in matching tissue",  "Missense somatic drivers in other tissues", "Other somatic drivers in other tissues", "Other somatic drivers in matching tissue", "Other germline drivers in other tissues", "No"))
sub.analyze$variable <- factor (sub.analyze$variable, levels = c ("MutSigCV", "OncodriveFM", "Hotspot", "OncodriveCLUST", "NMC", "MutSig-CL", "iSIMPRe", "iPAC", "GraphPAC", "SpacePAC", "CLUMPS", "e-Driver", "ActiveDriver", "LowMACA", "e-Driver3D"))

sub.brca.plot <- subset (sub.analyze, value < 0.05)
brca.plot <- ggplot (sub.brca.plot, aes (x = variable)) + geom_bar (position = "fill", aes (fill = AllDriver), color = "black") + theme (axis.text = element_text (color = "black"), axis.ticks.y = element_line (color = "black"), axis.ticks.x = element_blank(), axis.line.x = element_line (color = "black"), axis.line.y = element_line (color = "black"), panel.background = element_blank(), axis.text.x = element_blank(), legend.position = "none", axis.title = element_blank()) + scale_y_continuous (labels = scales::percent) + scale_fill_manual (values = c ("#FE0000", "orange","#2F496E", "#2988BC", "#5CC5EF", "white"), name = "Known driver role?") + ggtitle ("BRCA") + theme (plot.title = element_text (hjust = 0.5, size = 12, face = "italic"))

#NOW GBM

pred.data <- read.table (file = "../raw_data/matrix_gbm.txt", sep = "\t", header = TRUE, check.names = FALSE)
pred.melt <- melt (subset (pred.data, NPat > 4))
sub.analyze <- subset (pred.melt, !(variable %in% c ("NPat", "Group", "Detected")))

sub.analyze$AllDriver <- factor (sub.analyze$AllDriver, levels = c ("Missense somatic drivers in matching tissue",  "Missense somatic drivers in other tissues", "Other somatic drivers in other tissues", "Other somatic drivers in matching tissue", "Other germline drivers in other tissues", "No"))
sub.analyze$variable <- factor (sub.analyze$variable, levels = c ("MutSigCV", "OncodriveFM", "Hotspot", "OncodriveCLUST", "NMC", "MutSig-CL", "iSIMPRe", "iPAC", "GraphPAC", "SpacePAC", "CLUMPS", "e-Driver", "ActiveDriver", "LowMACA", "e-Driver3D"))

sub.gbm.plot <- subset (sub.analyze, value < 0.05)
gbm.plot <- ggplot (sub.gbm.plot, aes (x = variable)) + geom_bar (position = "fill", aes (fill = AllDriver), color = "black") + theme (axis.text = element_text (color = "black"), axis.ticks.y = element_line (color = "black"), axis.ticks.x = element_blank(), axis.line.x = element_line (color = "black"), axis.line.y = element_line (color = "black"), panel.background = element_blank(), axis.text.x = element_blank(), legend.position = "none", axis.title = element_blank()) + scale_y_continuous (labels = scales::percent) + scale_fill_manual (values = c ("#FE0000","orange", "#2988BC", "white"), name = "Known driver role?")  + ggtitle ("GBM") + theme (plot.title = element_text (hjust = 0.5, size = 12, face = "italic"))

#NOW LUAD

pred.data <- read.table (file = "../raw_data/matrix_luad.txt", sep = "\t", header = TRUE, check.names = FALSE)
pred.melt <- melt (subset (pred.data, NPat > 4))
sub.analyze <- subset (pred.melt, !(variable %in% c ("NPat", "Group", "Detected")))

sub.analyze$AllDriver <- factor (sub.analyze$AllDriver, levels = c ("Missense somatic drivers in matching tissue",  "Missense somatic drivers in other tissues", "Other somatic drivers in other tissues", "Other somatic drivers in matching tissue", "Other germline drivers in other tissues", "No"))
sub.analyze$variable <- factor (sub.analyze$variable, levels = c ("MutSigCV", "OncodriveFM", "Hotspot", "OncodriveCLUST", "NMC", "MutSig-CL", "iSIMPRe", "iPAC", "GraphPAC", "SpacePAC", "CLUMPS", "e-Driver", "ActiveDriver", "LowMACA", "e-Driver3D"))

sub.luad.plot <- subset (sub.analyze, value < 0.05)
luad.plot <- ggplot (sub.luad.plot, aes (x = variable)) + geom_bar (position = "fill", aes (fill = AllDriver), color = "black") + theme (axis.text = element_text (color = "black"), axis.ticks = element_line (color = "black"), axis.line.x = element_line (color = "black"), axis.line.y = element_line (color = "black"), panel.background = element_blank(), axis.text.x = element_text (angle = 45, hjust = 1), legend.position = "none", axis.title = element_blank()) + scale_y_continuous (labels = scales::percent) + scale_fill_manual (values = c ("#FE0000", "orange", "#2F496E", "#5CC5EF", "white"), name = "Known driver role?") + ggtitle ("LUAD") + theme (plot.title = element_text (hjust = 0.5, size = 12, face = "italic"))

#PUT IT TOGETHER

grid.arrange (arrangeGrob (blca.plot, brca.plot, gbm.plot, luad.plot, ncol = 1, heights = c (3,3,3,5)), plot.legend, ncol = 2, widths = c (3,2))
