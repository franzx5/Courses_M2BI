# layout.respect=FALSE for layout function, if TRUE, unit column-width is the same physical measurement on the device as a unit row-height (where high narrow heatmaps when nrow >> ncol
# aggl.method="complete" for hclustfun
# lh.plot: relative extra heigth of plotting area (usual is 4)
heatmapHW <- function (x, Rowv = NULL, Colv = if (symm) "Rowv" else NULL,distfun = dist, hclustfun = hclust, add.expr,
	 symm = FALSE, revC = identical(Colv, "Rowv"), scale = c("row", "column","none"), na.rm = TRUE, margins = c(5, 5), 
	 ColSideColors, RowSideColors, cexRow = 0.2 + 1/log10(nr), cexCol = 0.2 + 1/log10(nc), labRow = NULL, labCol = NULL, 
	 main = NULL,  xlab = NULL, ylab = NULL, lh.plot=0,keep.dendro = FALSE, verbose = getOption("verbose"),
	 layout.respect=FALSE,colors,aggl.method="complete",zlim=c(min(x,na.rm = TRUE),max(x,na.rm = TRUE)),    ...) {
    scale <- if (symm && missing(scale)) 
        "none"
    else match.arg(scale)
    if (length(di <- dim(x)) != 2 || !is.numeric(x))
        stop("`x' must be a numeric matrix")
    nr <- di[1]
    nc <- di[2]
    if (nr <= 1 || nc <= 1)
        stop("`x' must have at least 2 rows and 2 columns")
    if (!is.numeric(margins) || length(margins) != 2)
        stop("`margins' must be a numeric vector of length 2")
    doRdend <- !identical(Rowv, NA)
    doCdend <- !identical(Colv, NA)
    if (is.null(Rowv))
        Rowv <- rowMeans(x, na.rm = na.rm)
    if (is.null(Colv))
        Colv <- colMeans(x, na.rm = na.rm)
    if (doRdend) {
        if (inherits(Rowv, "dendrogram"))
            ddr <- Rowv
        else {
            hcr <- hclustfun(distfun(x),method=aggl.method)
            ddr <- as.dendrogram(hcr)
            #if (!is.logical(Rowv) || Rowv)
            #    ddr <- reorder(ddr, Rowv)
        }
        if (nr != length(rowInd <- order.dendrogram(ddr)))
            stop("row dendrogram ordering gave index of wrong length")
    }
    else rowInd <- 1:nr
    if (doCdend) {
        if (inherits(Colv, "dendrogram"))
            ddc <- Colv
        else if (identical(Colv, "Rowv")) {
            if (nr != nc)
                stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
            ddc <- ddr
        }
        else {
            hcc <- hclustfun(distfun(if (symm)
                x
            else t(x)))
            ddc <- as.dendrogram(hcc)
            if (!is.logical(Colv) || Colv)
                ddc <- reorder(ddc, Colv)
        }
        if (nc != length(colInd <- order.dendrogram(ddc)))
            stop("column dendrogram ordering gave index of wrong length")
    }
    else colInd <- 1:nc
    x <- x[rowInd, colInd]
    labRow <- if (is.null(labRow))
        if (is.null(rownames(x)))
            (1:nr)[rowInd]
        else rownames(x)
    else labRow[rowInd]
    labCol <- if (is.null(labCol))
        if (is.null(colnames(x)))
            (1:nc)[colInd]
        else colnames(x)
    else labCol[colInd]
    if (scale == "row") {
        x <- sweep(x, 1, rowMeans(x, na.rm = na.rm))
        sx <- apply(x, 1, sd, na.rm = na.rm)
        x <- sweep(x, 1, sx, "/")
    }
    else if (scale == "column") {
        x <- sweep(x, 2, colMeans(x, na.rm = na.rm))
        sx <- apply(x, 2, sd, na.rm = na.rm)
        x <- sweep(x, 2, sx, "/")
    }
    lmat <- rbind(c(NA, 3), 2:1)
    lwid <- c(if (doRdend) 1 else 0.05, 4)
    lhei <- c((if (doCdend) margins[1]/(5+lh.plot) else 0.05) + if (!is.null(main)) 0.2 else 0, 4)
    if (!missing(ColSideColors)) {
        if (!is.character(ColSideColors) || !(length(ColSideColors) ==  nc || nrow(ColSideColors) == nc))
            stop("'ColSideColors' must be a character vector of length ncol(x) or a matrix/dataframe with nrows ncol(x)")
        ColSideColors <- as.matrix(ColSideColors)
	lmat <- rbind(lmat[1, ] + ncol(ColSideColors), cbind(rep(NA,ncol(ColSideColors)), 1:ncol(ColSideColors)), lmat[2, ] + ncol(ColSideColors))
        lhei <- c(lhei[1], rep(0.2,ncol(ColSideColors)), lhei[2])
    }
    if (!missing(RowSideColors)) {
        if (!is.character(RowSideColors) || !(length(RowSideColors) == nr || nrow(RowSideColors) == nr))
            stop("'RowSideColors' must be a character vector of length nrow(x)")
	RowSideColors <- as.matrix(RowSideColors)
        lmat <- cbind(lmat[, 1] + ncol(RowSideColors), rbind(matrix(rep(NA, (nrow(lmat) - 1)*ncol(RowSideColors)),ncol=ncol(RowSideColors)),1:ncol(RowSideColors)), lmat[, 2] + ncol(RowSideColors))
        lwid <- c(lwid[1], rep(0.2,ncol(RowSideColors)), lwid[2])
    }
    lmat[is.na(lmat)] <- 0
    if (verbose) {
        cat("layout: widths = ", lwid, ", heights = ", lhei,
            "; lmat=\n")
        print(lmat)
    }
		
	lmat2 <- matrix(c(rep(0,3*ncol(lmat))),3,ncol(lmat),byrow=TRUE)
	lmat2[2,ncol(lmat)] <- max(lmat)+1
	lmat <- rbind(lmat,lmat2)
	lhei <- c(lhei,margins[1]/5,0.15,margins[1]/15)
	if(lh.plot){lhei[lhei==4]<-4+lh.plot}

    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    layout(lmat, widths = lwid, heights = lhei, respect = layout.respect)
    if (!missing(RowSideColors)) {
	for(r in 1:ncol(RowSideColors)){
		par(mar = c(margins[1], 0, 0, 0.5))
        	image(rbind(1:nr), col = RowSideColors[rowInd,r], axes = FALSE) }
    }
    if (!missing(ColSideColors)) {
        for(c in 1:ncol(ColSideColors)){
	  par(mar = c(0.5, 0, 0, margins[2]))
          image(cbind(1:nc), col = ColSideColors[colInd,c], axes = FALSE) }
    }
    par(mar = c(margins[1], 0, 0, margins[2]))
    if (!symm || scale != "none")
        x <- t(x)
    if (revC) {
        iy <- nr:1
        ddr <- rev(ddr)
        x <- x[, iy]
    }
    else iy <- 1:nr

    image(1:nc,1:nr,matrix(rep(1,nrow(x)*ncol(x)),nrow(x),ncol(x)), xlim = 0.5 + c(0, nc), ylim = 0.5 + c(0, nr), axes = FALSE, xlab = "", ylab = "",col="blue",...)
    image(1:nc,1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + c(0, nr), axes = FALSE, xlab = "", ylab = "",col=colors, add=TRUE,zlim=zlim,...)

    axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0,
        cex.axis = cexCol)
    if (!is.null(xlab))
        mtext(xlab, side = 1, line = margins[1] - 1.25)
    axis(4, iy ,labels = labRow, las = 2, line = 0, tick = 1, cex.axis = cexRow)
    if (!is.null(ylab))
        mtext(ylab,side = 4, line = margins[2] - 1.25)
    if (!missing(add.expr))
        eval(substitute(add.expr))
    par(mar = c(margins[1], 0, 0, 1))
    if (doRdend)
        plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
    else {
    	par(mar=c(0,0,0,0))
    	frame()    }
    par(mar = c(0, 0, if (!is.null(main)) 1 else 0, margins[2]))
    if (doCdend)
        plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
    else if (!is.null(main))
        frame()
    colorbar(seq(zlim[1],zlim[2],diff(zlim)/(length(colors)-1)),col=colors)
    if (!is.null(main))
	       title(main, outer=TRUE,cex.main = 2)
    invisible(list(rowInd = rowInd, colInd = colInd, Rowv = if (keep.dendro &&
        doRdend) ddr, Colv = if (keep.dendro && doCdend) ddc))
}
# Cordist function for calculation of distance for heatmap clustering (by Laurent)
cordist <- function(m, want.log = TRUE) {
  if (want.log) {
    m <- log(m)
  }
  d <- abs(1 - cor(t(m)))
  d <- as.dist(d)
  return(d)
}

colorbar <- function (x, col = heat.colors(50), scale = 1:length(x), k = 11, ...) {
    x <- x
    colmap <- col
    if (length(x) > k) 
        x.small <- seq(x[1], x[length(x)], length = k)
    else x.small <- x
    image(x, 1, matrix(x, length(x), 1), axes = FALSE, xlab = "", ylab = "", col = colmap, ...)
    axis(1, at = rev(x.small), labels = signif(rev(x.small),2), srt = 270)
    box()
}
