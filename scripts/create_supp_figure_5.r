library (ggplot2)
library (reshape)
library (gridExtra)

tissues <- c ("blca", "brca", "gbm", "luad")

m.res <- matrix (nrow = 15*5, ncol = 6)
colnames (m.res) <- c ("Method", "Tissue", "nDetOG", "nDet", "nOG", "nKnown")
m.res <- as.data.frame (m.res)

for (t in 1:length (tissues)) {
	
	file.name <- paste ("../raw_data/matrix_", tissues[t], ".txt", sep = "")
	
	pred.data <- read.table (file = file.name, sep = "\t", header = TRUE, check.names = FALSE)
	pred.melt <- melt (subset (pred.data, NPat > 4))

	sub.analyze <- subset (pred.melt, !(variable %in% c ("NPat", "Group", "Detected")))
	methods <- unique (sub.analyze$variable)

	for (i in 1:length(methods)) {
	
		sub.m <- subset (sub.analyze, variable == methods[i])
	
		og.detected <- nrow (subset (sub.m, Mode %in% c ("OG", "OG/TSG") & value < 0.05 & KnownDriver == 'Yes'))
		total.detected <- nrow (subset (sub.m, value < 0.05 & KnownDriver == 'Yes'))
		total.og <- nrow (subset (sub.m, Mode %in% c ("OG", "OG/TSG") & KnownDriver == 'Yes'))
		total.known <- nrow (subset (sub.m, KnownDriver == "Yes"))
	
		m.res[i + (15 * (t-1)),1] <- as.character (methods[i])
		m.res[i + (15 * (t-1)),2] <- as.character (tissues[t])
		m.res[i + (15 * (t-1)),3] <- og.detected
		m.res[i + (15 * (t-1)),4] <- total.detected
		m.res[i + (15 * (t-1)),5] <- total.og
		m.res[i + (15 * (t-1)),6] <- total.known
	}
}

for (i in 1:length(methods)) {
	
	sub.m <- subset (m.res, Method == as.character (methods[i]))
	
	og.detected <- sum (sub.m$nDetOG)
	total.detected <- sum (sub.m$nDet)
	total.og <- sum (sub.m$nOG)
	total.known <- sum (sub.m$nKnown)
	
	m.res[i + (15 * (4)),1] <- as.character (methods[i])
	m.res[i + (15 * (4)),2] <- "Pancan"
	m.res[i + (15 * (4)),3] <- og.detected
	m.res[i + (15 * (4)),4] <- total.detected
	m.res[i + (15 * (4)),5] <- total.og
	m.res[i + (15 * (4)),6] <- total.known
}

calc.fish <- function (n11, n12, n21, n22) {
	
	m <- matrix (nrow = 2,, ncol = 2)
	
	m[1,1] <- n11
	m[1,2] <- n12 - n11
	m[2,1] <- n21 - n11
	m[2,2] <- n22 - n12 - n21 + 11
	
	result <- fisher.test (m, alternative = "greater")
 	
 	return (result$p.value)
}

calc.or <- function (n11, n12, n21, n22) {
	
	m <- matrix (nrow = 2,, ncol = 2)
	
	m[1,1] <- n11
	m[1,2] <- n12 - n11
	m[2,1] <- n21 - n11
	m[2,2] <- n22 - n12 - n21 + 11
	
	m <- m+1
	
	result <- fisher.test (m, alternative = "greater")
	
	if (abs (log10(result$conf.int[1])) < abs(log10(result$conf.int[2]))) {
		return (result$conf.int[1])
	} else {
		return (result$conf.int[2])
	}
}

m.res$p <- mapply (calc.fish, m.res$nDetOG, m.res$nDet, m.res$nOG, m.res$nKnown)
m.res$AdjP <- ifelse (m.res$p < 1e-5, 1e-5, m.res$p)
m.res$sig <- ifelse (m.res$p < 0.05, "Yes", "No")
m.res$OR <- mapply (calc.or, m.res$nDetOG, m.res$nDet, m.res$nOG, m.res$nKnown)

m.res$Method <- factor (m.res$Method, levels = c ("MutSigCV", "OncodriveFM", "Hotspot", "OncodriveCLUST", "NMC", "MutSig-CL", "iSIMPRe", "iPAC", "GraphPAC", "SpacePAC", "CLUMPS", "e-Driver", "ActiveDriver", "LowMACA", "e-Driver3D"))
m.res$Tissue <- factor (m.res$Tissue, levels = c ("blca", "brca", "gbm", "luad", "Pancan"), labels = c ("BLCA", "BRCA", "GBM", "LUAD", "Pancan"))

pred.data <- read.table ("../raw_data/matrix_pancan.txt", header = TRUE, check.names = FALSE, sep = "\t")
pred.melt <- melt (subset (pred.data, NPat > 4))

sub.analyze <- subset (pred.melt, !(variable %in% c ("NPat", "Group", "Detected")))
sub.analyze$variable <- factor (sub.analyze$variable, levels = c ("MutSigCV", "OncodriveFM", "Hotspot", "OncodriveCLUST", "NMC", "MutSig-CL", "iSIMPRe", "iPAC", "GraphPAC", "SpacePAC", "CLUMPS", "e-Driver", "ActiveDriver", "LowMACA", "e-Driver3D"))
sub.analyze$Mode <- factor (sub.analyze$Mode, levels = c ("OG", "TSG", "OG/TSG", "Unk", "N/A"), labels = c ("OG", "TSG", "OG/TSG", "Unknown", "N/A"))

sub.analyze$Tissue <- factor (sub.analyze$Tissue, levels = c ("blca", "brca", "gbm", "luad"), labels = c ("BLCA", "BRCA", "GBM", "LUAD"))

plot.1 <- ggplot (subset (sub.analyze, value < 0.05), aes (x = variable, fill = Mode)) + geom_bar (position = "fill", color = "black") + theme (axis.text = element_text (color = "black"), axis.ticks = element_line (color = "black"), axis.ticks.x = element_line (color = "black"), axis.line.x = element_line (color = "black"), axis.line.y = element_line (color = "black"), axis.title.y = element_text (face = "italic"), panel.background = element_blank(), legend.position = "bottom", axis.text.x = element_text (angle = 45, hjust = 1, vjust = 1), axis.title.x = element_blank()) + ylab ("Fraction of genes") + scale_y_continuous (labels = scales::percent) + scale_fill_manual (values = c ("#FFBB00", "#375E97", "#3F681C", "#DDC5A2", "white")) + facet_wrap (~Tissue, ncol = 2) + theme (strip.background = element_blank(), strip.text = element_text (face = "italic"))
plot.2 <- ggplot (subset (m.res, nDetOG > 0 & Tissue == "Pancan"), aes (x = Method, y = OR)) + geom_bar (stat = "identity", aes (alpha = sig, fill = sig), color = "black") + theme (axis.text.x = element_text (angle = 45, hjust = 1)) + geom_hline (yintercept = 1, linetype = 2, color = "black") + theme (panel.background = element_blank(), axis.text = element_text (color = "black"), axis.ticks = element_line (color = "black"), axis.line.x = element_line (color = "black"), axis.line.y = element_line (color = "black"), axis.title.x = element_blank(), axis.title.y = element_text (face = "italic"), legend.position = "right") + ylab ("Oncogene/dual role\n fold enrichment\n") + scale_alpha_manual (values = c (0.3, 1), name = "p < 0.05") + scale_fill_manual (values = c ("lightgray", "#FFBB00"), name = "p < 0.05")
grid.arrange (plot.1, plot.2, ncol = 1, heights = c (5,3))