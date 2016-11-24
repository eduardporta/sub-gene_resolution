library (ggplot2)
library (reshape)
library (reshape2)
library (pheatmap)

#THE SCRIPT STORES THE VALUES IN THIS MATRIX
m.res <- matrix (nrow = 15*4, ncol = 4)
colnames (m.res) <- c ("Method", "Precision", "Recall", "Tissue")
m.res <- as.data.frame (m.res)

#THE SCRIPT ITERATES THROUGH THE DIFFERENT MATRIXES AND CALCULATES THE PRECISION AND RECALL VALUES FOR EACH METHOD AND CATEGORY
tissues <- c ("blca", "brca", "gbm", "luad")

for (t in 1:length (tissues)) {
	
	file.name <- paste ("../raw_data/matrix_", tissues[t], ".txt", sep = "")
	
	pred.data <- read.table (file = file.name, sep = "\t", header = TRUE, check.names = FALSE)
	pred.melt <- melt (subset (pred.data, NPat > 4))

	sub.analyze <- subset (pred.melt, !(variable %in% c ("NPat", "Group", "Detected")))
	methods <- unique (sub.analyze$variable)

	for (i in 1:length(methods)) {
	
		sub.m <- subset (sub.analyze, variable == methods[i])
	
		true.pos <- nrow (subset (sub.m, KnownDriver == "Yes" & value < 0.05)) #NUMBER OF TRUE POSITIVES BY THE METHOD
		total.pos <- nrow (subset (sub.m, value < 0.05)) #NUMBER OF POSITIVES IDENTIFIED BY THE METHOD
		total.true <- nrow (subset (sub.m, KnownDriver == "Yes")) #NUMBER OF TOTAL TRUE GENES
	
		m.res[i + (15 * (t-1)),1] <- as.character (methods[i])
		m.res[i + (15 * (t-1)),2] <- true.pos/total.pos #PRECISION
		m.res[i + (15 * (t-1)),3] <- true.pos/total.true #RECALL
		m.res[i + (15 * (t-1)),4] <- as.character (tissues[t])
	}
}

#NOW THE SAME PER CATEGORIES

categories <- c ("Whole gene", "Type I", "Type II", "Type III", "Type IV")

for (t in 1:length (tissues)) {
	
	file.name <- paste ("../raw_data/matrix_", tissues[t], ".txt", sep = "")
	
	pred.data <- read.table (file = file.name, sep = "\t", header = TRUE, check.names = FALSE)
	pred.melt <- melt (subset (pred.data, NPat > 4))

	sub.analyze <- subset (pred.melt, !(variable %in% c ("NPat", "Group", "Detected")))
	
	for (c in 1:length (categories)) {
		
		sub.m <- sub.analyze
		
		if (categories[c] == "Whole gene") {
			sub.m <- subset (sub.analyze, variable %in% c("OncodriveFM", "MutSigCV"))
		}
		else if (categories[c] == "Type I") {
			sub.m <- subset (sub.analyze, variable %in% c("MutSig-CL", "OncodriveCLUST", "NMC", "Hotspot", "iSIMPRe"))
		}
		else if (categories[c] == "Type II") {
			sub.m <- subset (sub.analyze, variable %in% c("iPAC", "GraphPAC", "SpacePAC", "CLUMPS"))
		}
		else if (categories[c] == "Type III") {
			sub.m <- subset (sub.analyze, variable %in% c("e-Driver", "ActiveDriver", "LowMACA"))
		}
		else if (categories[c] == "Type IV") {
			sub.m <- subset (sub.analyze, variable %in% c("e-Driver3D"))
		}
		
		true.pos <- nrow (subset (sub.m, KnownDriver == "Yes" & value < 0.05))
		total.pos <- nrow (subset (sub.m, value < 0.05))
		total.true <- nrow (subset (sub.m, KnownDriver == "Yes"))
			
		m.res[15*4 + 5*(t-1) + c, 1] <- as.character (categories[c])
		m.res[15*4 + 5*(t-1) + c, 2] <- true.pos/total.pos
		m.res[15*4 + 5*(t-1) + c, 3] <- true.pos/total.true
		m.res[15*4 + 5*(t-1) + c, 4] <- as.character (tissues[t])
	}
}

m.res$Tissue <- factor (m.res$Tissue, levels = c ("blca", "brca", "gbm", "luad"), labels = c ("BLCA", "BRCA", "GBM", "LUAD"))

#HEATMAP METHODS

bk2 = unique(c(seq(0, 1, length=50)))
col1 <- colorRampPalette (c ("#4CB5F5", "white"))(15)
col2 <- colorRampPalette(c("white", "#F62A00"))(35)

#FIG 3A
recall.data <- acast (subset (m.res, Method %in% c ("Whole gene", "Type I", "Type II", "Type III", "Type IV")), Method ~ Tissue, value.var = "Recall")
pheatmap (recall.data, main = "Recall\n", breaks = bk2, color = c (col1, col2))

#FIG 3B
precision.data <- acast (subset (m.res, Method %in% c ("Whole gene", "Type I", "Type II", "Type III", "Type IV")), Method ~ Tissue, value.var = "Precision")
pheatmap (precision.data, main = "Precision\n", breaks = bk2, color = c (col1, col2))

#FIG 3C
recall.data <- acast (subset (m.res, !(Method %in% c ("Whole gene", "Type I", "Type II", "Type III", "Type IV"))), Method ~ Tissue, value.var = "Recall")
pheatmap (recall.data, main = "Recall\n", breaks = bk2, color = c (col1, col2))

#FIG 3D
precision.data <- acast (subset (m.res, !(Method %in% c ("Whole gene", "Type I", "Type II", "Type III", "Type IV"))), Method ~ Tissue, value.var = "Precision")
pheatmap (precision.data, main = "Precision\n", breaks = bk2, color = c (col1, col2))
