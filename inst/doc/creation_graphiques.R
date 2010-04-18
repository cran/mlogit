source("fgraphiques.R")

pdf("./graph/opt1.pdf",onefile=FALSE,width=8,height=6)
goptimisation(1)
dev.off()

pdf("./graph/opt2.pdf",onefile=FALSE,width=8,height=6)
goptimisation(2)
dev.off()

pdf("./graph/opt3.pdf",onefile=FALSE,width=8,height=6)
goptimisation(2)
dev.off()

pdf("./graph/halton1.pdf",onefile=FALSE,width=8,height=4)
ghalton(1)
dev.off()

pdf("./graph/halton2.pdf",onefile=FALSE,width=8,height=4)
ghalton(2)
dev.off()

pdf("./graph/halton3.pdf",onefile=FALSE,width=8,height=4)
ghalton(3)
dev.off()

pdf("./graph/draws.pdf",onefile=FALSE)
gunifgumbel(20)
dev.off()

pdf("./graph/correlation.pdf",onefile=FALSE,width=8,height=4)
gcorrelation(200)
dev.off()

pdf("./graph/haltonvsrandom.pdf",onefile=FALSE,width=8,height=4)
haltonvsrandom(100)
dev.off()
