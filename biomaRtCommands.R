library('affy')
library('gplots')
library('lattice')
library('biomaRt')

mitocarta = read.table("MitoCarta Mouse Ensembl List.txt")
ensemblHuman = useMart('ensembl', dataset='hsapiens_gene_ensembl')
ensemblMouse = useMart('ensembl', dataset='mmusculus_gene_ensembl')

mapmitocarta1 = getBM(attributes = c('ensembl_gene_id','hsapiens_homolog_ensembl_gene'), filters = c('ensembl_gene_id','with_homolog_hsap'), values = list( as.vector(mitocarta[,1]),TRUE ), mart = ensemblMouse)
mapmitocarta2 = getBM(attributes = c('ensembl_gene_id','entrezgene'), filters = c('ensembl_gene_id','with_entrezgene'), values = list( as.vector(mapmitocarta1[,2]),TRUE ), mart = ensemblHuman)

matchIdxs1 = match(intersect( as.vector(mapmitocarta1[,2]),as.vector(mapmitocarta2[,1]) ),mapmitocarta1[,2])
matchIdxs2 = match(intersect( as.vector(mapmitocarta1[,2]),as.vector(mapmitocarta2[,1]) ),mapmitocarta1[,2])

mapmitocarta3 = cbind(mapmitocarta1[matchIdxs1,1], mapmitocarta2[matchIdxs2,2])
write.table(mapmitocarta3, file="MitoCarta Mouse Ensembl To Human Entrez.csv", quote=FALSE, sep=",", row.names=FALSE)