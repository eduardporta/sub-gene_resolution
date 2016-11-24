library (ggplot2)

cov.data <- read.table ("../raw_data/coverage_proteins.txt", sep = "\t", header = TRUE)

#Adjust order
cov.data$Variable <- factor (cov.data$Variable, levels = c ("Proteins", "Aminoacids"))

#Adjust labels
cov.data$Feature <- factor (cov.data$Feature, levels = c ("Linear regions", "Structures (>95%)", "Structures (e-value < 1e-9)"), labels = c ("Linear\nregions", "Structures\n(>95% identity)", "Structures\n(evalue < 1e-9)"))

#Plot
ggplot (cov.data, aes (x = Variable, y = Percent)) + geom_bar (stat = "identity", aes (fill = Covered), color = "black") + facet_wrap (~Feature) + theme (axis.text = element_text (color = "black"), axis.ticks = element_line (color = "black"), axis.line.x = element_line (color = "black"), axis.line.y = element_line (color = "black"), panel.background = element_blank(), axis.title = element_blank(), axis.text.x = element_text (angle = 45, hjust = 1), strip.background = element_blank(), strip.text = element_text (face = "italic", size = 10), legend.position = "bottom") + scale_y_continuous (labels = scales::percent) + scale_fill_manual (values = c ("orange", "gray"))
