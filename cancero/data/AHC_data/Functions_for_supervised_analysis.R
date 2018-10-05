ExpCompare <- function(exp=NULL,
							refsamp=NULL,
							gpsamp=NULL,
							featureTab=NULL
							)
{
	res <- featureTab
	meanExp.ref <- apply(matrix(exp[,refsamp],ncol=length(refsamp)),1,mean,na.rm=T)
	meanExp.gp <- apply(matrix(exp[,gpsamp],ncol=length(gpsamp)),1,mean,na.rm=T)
	FC <- (2^meanExp.gp)/(2^meanExp.ref)
	log2FC <- log2(FC)
	pval <- sapply(1:nrow(exp),function(i){t.test(exp[i,refsamp],exp[i,gpsamp])$p.value})
	qval <- p.adjust(pval)
	res <- data.frame(res,meanExp.ref,meanExp.gp,FC,log2FC,pval,qval)
	res	
}

enrichmentTest <- function (gene.sets, mylist, possibleIds, sep = ";", silent = F) 
{
    possibleIds <- unique(possibleIds)
    mylist <- unique(mylist)
    gene.sets <- lapply(gene.sets, unique)
    nids <- length(possibleIds)
    gene.sets <- lapply(gene.sets, function(x) intersect(x, possibleIds))
    nref <- sapply(gene.sets, length)
    if (all(nref == 0)) stop("Error: no intersection between gene sets and possible IDs.")
    if (any(nref == 0)) print("Warning: some of the gene sets have no intersection with possibleIds")
    if (!all(mylist %in% possibleIds)) stop("Error: some genes in mylist are not in possibleIds")
    if (!silent) cat(paste("NB : enrichment tests are based on", nids, "distinct ids.\n"))
    gene.sets <- gene.sets[nref > 0]
    n <- length(mylist)
    fun <- function(x) {
        y <- intersect(x, mylist)
        nx <- length(x)
        ny <- length(y)
        pval <- phyper(ny - 1, nx, nids - nx, n, lower.tail = F)
        c(nx, ny, pval,paste(y, collapse = sep))
    }
    tmp <- as.data.frame(t(sapply(gene.sets, fun)))
    rownames(tmp) <- names(gene.sets)
    for (i in 1:3) tmp[,i] <- as.numeric(as.character(tmp[,i]))
    tmp <- data.frame(tmp[,1:3],p.adjust(tmp[,3],method="BH"),tmp[,4])
    names(tmp) <- c("Nb of genes in population","Nb of genes in list","Hypergeometric test p-value","BH adjusted p-value","corresponding elements in list")
    tmp
}
