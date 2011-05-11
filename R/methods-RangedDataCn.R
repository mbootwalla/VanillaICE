setMethod("initialize", "RangedDataCn",
	  function(.Object,
		   st=integer(),
		   en=integer(),
		   nmarkers=integer(),
		   id=character(),
		   chrom=integer(),
		   state=integer(),
		   LLR=numeric(), ...){
	.Object <- callNextMethod()
	.Object$sampleId <- id
	.Object$chrom <- chrom
	.Object$numMarkers <- nmarkers
	.Object$state <- state
	.Object$LLR <- LLR
	.Object
})

setMethod("chromosome", "RangedDataCn", function(object) object$chrom)
setMethod("sampleNames", "RangedDataCn", function(object) object$sampleId)
setMethod("nMarkers", "RangedDataCn", function(object) object$numMarkers)
setMethod("LLR", "RangedDataCn", function(object) object$LLR)
setMethod("state", "RangedDataCn", function(object) object$state)

setMethod("rd2df", "RangedDataCn", function(object, palette, ...){
	isAmp <- state(object) > 3
	isHom <- state(object) == 1
	isHem <- state(object) == 2
	cols <- rep(NA, nrow(object))
	cols[isHom] <- palette[1]
	cols[isHem] <- palette[2]
	cols[isAmp] <- palette[3]
	h <- 0.75
	meanSegment <- apply(cbind(start(object), end(object)), 1, mean)
	data(chromosomeAnnotation)
	chr.size <- chromosomeAnnotation[1:22, "chromosomeSize"]
	chrom <- chromosome(object)
	chr.size <- chr.size[chrom]
	y <- split(sampleNames(object), chrom)
	y <- lapply(y, function(x){
		tmp <- as.numeric(as.factor(x))
		names(tmp) <- as.character(x)
		tmp
	})
	y <- unlist(y)
	nms2 <- paste(chromosome(object), sampleNames(object), sep=".")
	if(!identical(names(y), nms2)){
		y <- y[match(nms2, names(y))]
		stopifnot(identical(names(y), nms2))
	}
	dat <- data.frame(x0=start(object)/1e6,
			  x1=end(object)/1e6,
			  y0=y-h/2,
			  y1=y+h/2,
			  chr=chromosome(object),
			  coverage=nMarkers(object),
			  midpoint=meanSegment/1e6,
			  id=sampleNames(object),
			  chr.size=chr.size/1e6,
			  col=cols,
			  y=y,
			  stringsAsFactors=FALSE)
	dat$chr <- as.factor(dat$chr)
	dat <- as(dat, "data.frame.CN")
	return(dat)
})



setMethod("plot", signature(x="RangedDataCn"), function(x, y, palette, border="black",
			    labelAllSamples=TRUE,
			    show.coverage=TRUE,
			    sampleLabels.cex=0.5, ...){
	df <- rd2df(x, palette=palette)
	plot(df, palette=palette, border=border, labelAllSamples=labelAllSamples,
	     show.coverage=show.coverage,
	     sampleLabels.cex=sampleLabels.cex, ...)
})
setMethod("plot", signature(x="data.frame.CN"), function(x, y, palette, border="black",
			    labelAllSamples=TRUE,
			    show.coverage=TRUE,
			    sampleLabels.cex=0.5, ...){
	stopifnot(length(unique(df$chr))==1)
	mykey <- simpleKey(c("homo-del", "hemi-del", "amp")[palette %in% df$col], points=FALSE,
			   rectangles=TRUE, col=palette[palette %in% df$col], space="top")
	mykey$rectangles[["border"]] <- mykey$rectangles[["col"]] <- palette[palette %in% df$col]
	if(border=="black") border <- rep("black", nrow(df)) else border <- df$col
	if(labelAllSamples) {
		labels <- df$id
		ticks.at <- df$y
	} else {
		labels <- FALSE
		ticks.at <- pretty(df$y)
	}
	fig <- xyplot(y~midpoint, data=df,
		      panel=function(x, y, x0, x1, chr.size,
		      col, border, coverage, chr, show.coverage=TRUE, max.y,
		      ..., subscripts){
			      panel.grid(h=-1, v=10)
			      panel.xyplot(x, y, ..., subscripts)
			      h <- 0.75
			      lrect(xleft=x0[subscripts],
				    xright=x1[subscripts],
				    ybottom=y-h/2,
				    ytop=y+h/2,
				    border=border[subscripts],
				    col=col[subscripts], ...)
			      if(show.coverage)
				      ltext(x, y,labels=coverage[subscripts], cex=0.6)
			      ##plot centromere
			      chr <- unique(as.integer(as.character(df$chr)))
			      coords <- chromosomeAnnotation[chr, 1:2]/1e6
			      lrect(xleft=coords[1],
				    xright=coords[2],
				    ybottom=0,
				    ytop=max.y+h/2,
				    col="grey",
				    border="grey")
		      },
		      x0=df$x0,
		      x1=df$x1,
		      col=df$col,
		      border=border,
		      alpha=1,
		      chr.size=df$chr.size,
		      scales=list(y=list(labels=labels, at=ticks.at, cex=sampleLabels.cex)),
		      coverage=df$coverage,
		      xlab="Mb",
		      ylab="offspring index",
		      show.coverage=show.coverage,
		      key=mykey,
		      par.strip.text=list(lines=0.7, cex=0.6),
		      prepanel=prepanel.fxn,
		      max.y=max(df$y), ...)
	##axis=function(side, text.cex){
	##panel.axis(side, text.cex=text.cex)}, ...)
	return(fig)
})

prepanel.fxn <- function(x,y, chr, chr.size, ..., subscripts){
	list(xlim=c(0, unique(chr.size[subscripts])), ylim=range(as.integer(as.factor(y[subscripts]))))
}

