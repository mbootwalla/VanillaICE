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

setMethod("rd2df", "RangedDataCn", function(object, hmm.params, ...){
	require(SNPchip)
	data(chromosomeAnnotation)
	object <- object[object$state != normalIndex(hmm.params), ]
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
			  border=rep("black", nrow(object)),
			  col=as.integer(as.factor(state(object))),
			  state=state(object),
			  y=y,
			  stringsAsFactors=FALSE)
	dat$chr <- as.factor(dat$chr)
	dat <- new("DataFrameCN", dat)
	return(dat)
})


setMethod("plot", signature(object="RangedDataCn", hmm.params="HmmOptionList"),
	  function(object, hmm.params, ...){
		  df <- rd2df(object, hmm.params)
		  plot(object=df, ...)
	  })

setMethod("plot", signature(object="DataFrameCN", hmm.params="missing"),
	  function(object, ...){
		  data(chromosomeAnnotation)
		  df <- as.data.frame(object@.Data)
		  colnames(df) <- names(object)
		  df$x <- df$midpoint
		  fig <- xyplot(y~x, data=df,
				panel=my.xypanel,
				x0=df$x0,
				x1=df$x1,
				col=df$col,
				border=df$border,
				alpha=1,
				chr=df$chr,
				chr.size=df$chr.size,
				coverage=df$coverage,
				xlab="Mb",
				ylab="offspring index",
				key=getKey(df),
				par.strip.text=list(lines=0.7, cex=0.6),
				prepanel=prepanel.fxn,
				max.y=max(df$y),
				chromosomeAnnotation=chromosomeAnnotation,
				...)
		  return(fig)
	  })

getKey <- function(df){
	states <- unique(df$state)
	states <- states[order(states)]
	col <- df$col[match(states, df$state)]
	mykey <- simpleKey(c("homo-del", "hemi-del", "normal", "single-dup", "double-dup")[states], points=FALSE,
			   rectangles=TRUE, col=col, space="top")
	mykey$rectangles[["border"]] <- mykey$rectangles[["col"]] <- col
	mykey
}

my.xypanel <- function(x, y,
		       x0, x1, chr.size,
		       col, border, coverage,
		       chr, show.coverage=TRUE,
		       max.y,
		       chromosomeAnnotation,
		       addCentromere=TRUE,
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
	if(addCentromere){
		chr <- unique(as.integer(as.character(chr)))
		coords <- chromosomeAnnotation[chr, 1:2]/1e6
		lrect(xleft=coords[1],
		      xright=coords[2],
		      ybottom=0,
		      ytop=max.y+h/2,
		      col="grey",
		      border="grey")
	}
}

prepanel.fxn <- function(x,y, chr.size, ..., subscripts){
	list(xlim=c(0, unique(chr.size[subscripts])), ylim=range(as.integer(as.factor(y[subscripts]))))
}

