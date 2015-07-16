library("RDAVIDWebService")
david = DAVIDWebService$new("yw595@cornell.edu")
setEmail(david,"yw595@cornell.edu")
connect(david)

data(demoList1)
result = addList(david, demoList1, idType = "AFFYMETRIX_3PRIME_IVT_ID", listName = "demoList1", listType = "Gene")

allEnsembl = read.table("AllEnsemblMouseIDs.txt")
lookSig6 = read.table("lookSigEnsemblIDs6.txt")
getIdTypes(david)
result = addList(david, as.vector(allEnsembl[,1]), idType = "ENSEMBL_GENE_ID", listName = "allEnsembl", listType = "Background")
result = addList(david, as.vector(lookSig6[,1]), idType = "ENSEMBL_GENE_ID", listName = "lookSig6", listType = "Gene")

termCluster = getClusterReport(david, type="Term")
annotationSummary = getAnnotationSummary(david)
geneList = getGeneListReport(david)
functionAnnotationChart = getFunctionalAnnotationChart(david)
functionAnnotationTable = getFunctionalAnnotationTable(david)