setwd("~/research/gelscore/src/")

file.out = F
output.type = 'svg'

fig.liquid.length = F
fig.liquid.structured.hydro = F
fig.liquid.hydro.corrs = F
domain.linear.models = T
fig.liquid.conservation = T
fig.liquid.slopes = F
fig.liquid.aadist = F
fig.pabp.sep.dist = F
fig.droplets = F
fig.compare.droplets = F
fig.layers = F
fig.pab1.pbp1.comp.covariation = F
fig.pab1.pbp1.hydro.covariation = F
fig.pab1.hydro.thermo.vs.meso = F
fig.pab1.sliding = F
fig.pabp.similarity = F
fig.coevolution = F
fig.pab1.sis1.coev = F
score.pc = F

all.aas = charlist('ACDEFGHIKLMNPQRSTVWY')

get.slope <- function(g) {
	res <- c(NA,NA)
	if (!is.null(g)) {
		s <- summary(g)$coefficients
		if (nrow(s)>1) {
			res <- s[2,1:2]
		}
	}
	res
}

load.species.data <- function(fname, id.fld, master.ids) {
	x <- read.table(fname, header=T, sep='\t', stringsAsFactor=F)
	x[match(master.ids, x[,id.fld]),]
	x
}

# Read overall species data
spec.data <- read.table("../data/species-info.txt", header=TRUE, quote='', sep='\t', stringsAsFactor=FALSE)
for (g in unique(spec.data$group1)) {
	spec.data[spec.data$group1==g,'group1.count'] <- nrow(spec.data[spec.data$group1==g,])
}

if (fig.liquid.length) {
	domain <- 'nsp1-fg'
	#domain <- 'pab1p'
	x <- load.species.data(paste("../data/", domain, "-composition.txt",sep=''), 'id', spec.data$short.name)
	x$col <- rep('black',nrow(x))
	uspec <- unique(spec.data$group1)
	cols <- myrainbow(length(uspec))
	for (xi in 1:length(uspec)) {
		g <- uspec[xi]
		col <- cols[xi]
		z <- na.omit(match(subset(spec.data,group1==g)$short.name, x$id))
		x[z,'col'] <- col
	}
	o <- order(x$total.length)
	names <- x$id[o]
	if (file.out) dev.out("fig.liquid.length", width=6, height=4, output.type=output.type)
	barplot(x$total.length[o], las=1, horiz=T, space=c(0,0), names.arg=NULL, col=x$col[o], xlim=c(0,max(x$total.length)), border=0, xlab='P domain length')
	#abline(h=min())
	legend(max(x$total.length)*2/3, 100, legend=uspec, col=cols, fill=cols, bty='n', border=0, cex=0.7)
	if (file.out) dev.off()
}

if (fig.liquid.structured.hydro) {
	liquid.domain <- 'nsp1-fg'
	structured.domain <- 'nsp1-struct'
	liquid.domain <- 'pab1p'
	structured.domain <- 'pab1rrm'
	#liquid.domain <- 'pbp1-pb'
	#structured.domain <- 'pbp1-lsmad'

	plot.groups <- T

	raw.structured.counts <- load.species.data(paste("../data/",structured.domain,"-composition.txt",sep=''), 'id', spec.data$short.name)
	raw.liquid.counts <- load.species.data(paste("../data/",liquid.domain,"-composition.txt",sep=''), 'id', spec.data$short.name)
	target.aas <- c('L','I','M','V') #,'G','P','N','Q')
	aacounts <- p.0(target.aas,'count')
	aafracs <- p.0('f',target.aas)

	summarize.proportions <- function(x, motifs, na.rm=T) {
		fracs <- x[,p.0(motifs,'count')]/x$total.length
		names(fracs) <- p.0('f',motifs)
		list(mean=sapply(fracs[,aafracs], mean, na.rm=na.rm), err=sapply(fracs[,aafracs], sderr, na.rm=na.rm))
	}

	# Get hydrophobicity information
	hyd <- read.table("~/research/base/data/chimera-hydrophobicities.txt", header=T)
	hyd.fld <- 'hh.hyd'
	z1 <- match(target.aas, hyd$aa)
	h <- hyd[z1,hyd.fld]
	o <- order(h)

	spec.groups <- NULL
	if (plot.groups) {
		# Calculate proportions for each species group
		# Only consider groups for which there are at least a given number of species
		min.species.count = 5
		spec.groups <- unique(subset(spec.data, group1.count>=min.species.count)$group1)
	}
	all.groups <- c('all',spec.groups)
	merge.proportions <- function(x) {
		res <- lapply(spec.groups, function(g) {
			y <- x[match(subset(spec.data, group1==g)$short.name, x$id),]
			summarize.proportions(y, target.aas)
		})
		names(res) <- spec.groups
		res$all <- summarize.proportions(x, target.aas)
		res
	}
	liquid.props <- merge.proportions(raw.liquid.counts)
	struct.props <- merge.proportions(raw.structured.counts)
	props <- list(struct.props, liquid.props)

	# Set up plotting information
	spec.cols <- c(myrainbow(length(spec.groups)),'black')
	names(spec.cols) <- c(spec.groups,'all')
	pchs <- c(14 + length(spec.groups):1)
	names(pchs) <- c(spec.groups)
	xlim <- c(min(h),max(h))
	prop <- props[[2]]
	max.prop <- max(sapply(props, function(m){sapply(all.groups, function(g){m[[g]]$mean})}) + sapply(props, function(m){sapply(all.groups, function(g){m[[g]]$err})}),na.rm=T)
	ylim <- c(0, 1.05*max.prop)
	#print(ylim)
	mars <- list(c(4,4,0,0),c(4,0,0,4))

	if (file.out) dev.out("fig.aa.usage.hydrophobicity", width=12, height=4, output.type=output.type)
	split.screen(c(1,3))
	for (xi in 1:2) {
		screen(xi)
		par(mar=mars[[xi]])
		ylab <- ''
		yaxt <- 'n'
		if (xi==1) {
			yaxt <- 's'
			ylab <- 'Mean proportion amino-acid usage'
		}
		plot(1, 1, type='n', xlim=xlim, ylim=ylim, yaxs='i', pch=16, xlab='Hydrophilicity', yaxt=yaxt, ylab=ylab, las=1)
		prop <- props[[xi]]
		for (spec in spec.groups) {
			y <- prop[[spec]]$mean
			err <- prop[[spec]]$err
			#lines(h[o], y[o], type='l', pch=16, col=tcol(spec.cols[[spec]],0.5),lwd=1)
			lines.err(h[o], y[o], y.lower=y[o]-err[o], y.upper=y[o]+err[o], type='l', pch=16, col=tcol(spec.cols[[spec]],0.5), lwd=1)
		}
		y <- prop[['all']]$mean
		err <- prop[['all']]$err
		lines.err(h[o], y[o], y.lower=y[o]-err[o], y.upper=y[o]+err[o], type='l', pch=16, col=spec.cols[['all']], lwd=2)
		text(h, y, label=target.aas, pos=3, cex=1.5)
		#if (xi==1)	legend(xlim[1], ylim[2], legend=all.groups, col=unlist(spec.cols[all.groups]), lwd=c(2,rep(1,length(all.groups)-1)), bty='n')
	}
	screen(3)
	par(mar=c(4,4,0,0))
	# Treat as Boltzmann: log(p/q) = -dGp.q
	p <- liquid.props[['all']]$mean
	prop.mat <- p %*% t(1/p)
	hv <- hyd[match(target.aas,hyd$aa),hyd.fld]
	n <- length(target.aas)
	dg.mat <- matrix(rep(hv,n),n,n) - t(matrix(rep(hv,n),n,n))

	d <- data.frame(dg=lt(dg.mat,diag=T), rel.prop=lt(prop.mat,diag=T))
	plot(d$dg, log(d$rel.prop), log='', pch=16, xlab='Hessa-von Heijne biological hydrophilicity\n\u0394\u0394G (kcal/mol)', ylab='Log relative proportion', las=1)
	pcor(ct <- cortest(d$dg, log(d$rel.prop), meth='p'))
	g <- lm(log(rel.prop)~dg+0, data=d)
	abline(g, lty='solid')
	abline(0,1,lty='dotted')

	close.screen(all=TRUE)
	if (file.out) dev.off()
}

if (fig.liquid.hydro.corrs) {
	y <- spec.data
	y$liquid.hcor <- apply(log(raw.liquid.counts[,aacounts][,o]+1), 1, cor, method='p', y=exp(h[o]))
	y$struct.hcor <- apply(log(raw.structured.counts[,aacounts][,o]+1), 1, cor, method='p', y=exp(h[o]))

	if (file.out) dev.out(p.0("fig",liquid.domain,"aa.usage.hydrophobicity.corrs"), width=4, height=4, output.type=output.type)
	par(mar=c(4,4,0,0))
	plot(1, 1, type='n', xlim=c(-0.4,1), ylim=c(-0.4,1), xaxs='r', yaxs='r', pch=16, xlab='Hydrophilicity--usage correlation, liquid domain', ylab='Hydrophilicity--usage correlation, structured domain(s)', las=1)
	points(y$liquid.hcor, y$struct.hcor, pch=unlist(pchs[y$group1]), col=unlist(spec.cols[y$group1]))
	abline(0,1,lty='dotted')
	legend(-0.4,1, legend=spec.groups,pch=unlist(pchs[spec.groups]), col=unlist(spec.cols[spec.groups]), bty='n')
	close.screen(all=TRUE)
	if (file.out) dev.off()
}

if (domain.linear.models) {
	domain <- 'pab1p'
	#domain <- 'nsp1-fg'
	#domain <- 'pbp1-pb'
	raw.counts <- read.table(paste("../data/",domain,"-composition.txt",sep=''), header=T)
	addl.mots <- all
	all.mots <- as.character(sapply(names(raw.counts), function(x){substring(x, 1, nchar(x)-6)})[3:ncol(raw.counts)-1])
	#compound.mots <- all.mots[sapply(all.mots, nchar)>1]
	compound.mots <- c(c("FY",'KR',"AGNPQ","ILMVW",'DE','ST'))
	motcount <- function(m, x=raw.counts) {
		res <- sapply(m, function(mot){
			rowSums(x[,p.0(sort(charlist(mot)),'count')])
			})
		colnames(res) <- p.0(m,'count')
		res
	}
	z <- match(p.0(compound.mots,'count'),colnames(raw.counts))
	if (any(is.na(z))) {
		raw.counts <- data.frame(raw.counts, motcount(compound.mots[is.na(z)]))
	}
	all.mots <- c(all.mots, compound.mots[is.na(z)])
	counts <- raw.counts
	#cat("# Only looking at FY.count>4\n")
	#counts <- subset(raw.counts, FY.count>4)
	fracs <- data.frame(counts[,p.0(all.mots,'count')]/counts$total.length, counts$total.length)
	names(fracs) <- c(p.0('f',all.mots),'total.length')
	dat <- counts

	log.fit <- T
	fn <- noop
	if (log.fit) {
		fn <- log.nz
	}
	lms <- lapply(all.mots, function(m){
		d <- data.frame(x=dat$total.length, y=dat[,p.0(m,'count')])
		#points(d$x, d$y, pch=16, col=tcol(cols[xi],0.5))
		#d <- subset(d, x>0 & y>0)
		g <- NULL
		if (nrow(d)>0) {
			g <- lm(fn(y)~fn(x), data=d)			
		}
		g
	})
	names(lms) <- all.mots
	dat <- fracs
	frac.lms <- lapply(all.mots, function(m){
		d <- data.frame(x=dat$total.length, y=dat[,p.0('f',m)])
		#points(d$x, d$y, pch=16, col=tcol(cols[xi],0.5))
		#d <- subset(d, x>0 & y>0)
		g <- NULL
		if (nrow(d)>0) {
			g <- lm(fn(y)~fn(x), data=d)			
		}
		g
	})
	names(frac.lms) <- all.mots
}

if (fig.liquid.conservation) {
	mots <- c('KR',"AGNPQ","ILMVW",'FY')

	cols <- myrainbow(length(mots))
	names(cols) <- mots
	more.cols <- myrainbow(length(all.mots)-length(mots))
	names(more.cols) <- all.mots[!(all.mots %in% mots)]
	all.cols <- c(cols,more.cols)

	dat <- counts

	mar <- c(4,4,1,1)

	if (file.out) dev.out(p.0("fig",domain,"aa.conservation"), width=8, height=8, output.type=output.type)
	split.screen(c(2,3))
	screen(1)
	par(mar=mar)
	xlim <- c(max(1,min(raw.counts$total.length)),max(raw.counts$total.length))
	xlim <- c(0,max(raw.counts$total.length))
	ylim <- c(1,max(raw.counts[,p.0(all.mots,'count')]))
	#xlim <- c(0,300)
	#ylim <- c(0,50)
	#xlim <- c(40,260)
	plot(1,1, type='n', log='', xlim=xlim, ylim=ylim, las=1, pch=16,
		xlab='Total domain length', ylab='Number of occurrences')
	for (xi in 1:length(mots)) {
		mot <- mots[xi]
		g <- lms[[mot]]
		# Background Poisson variability
		d <- data.frame(x=seq(1,max(xlim)*1.1,1))
		pred <- predict.lm(g, d)
		err <- 1.96*sqrt(abs(pred))
		poly.err(x=d$x, y=NULL, y.lower=pred-err, y.upper=pred+err, lines=F, col=tcol(cols[mot],0.2))

		# The points
		d <- data.frame(x=dat$total.length, y=dat[,p.0(mot,'count')])
		col <- tcol(cols[[mot]], 0.7)
		points(d$x, d$y, pch=16, col=col)

		# Lines
		d <- subset(d, x>0 & y>0)
		xbase <- seq(min(xlim), max(xlim), length.out=100)
		if (log.fit) {
			y <- exp(coef(g)[1])*xbase^(coef(g)[2])
			lines(xbase, y, col=col, lwd=2)
		} else {
			abline(g, col=col, lwd=2)
		}

		#lines(d$x, predict.lm(lms[[mot]], d)+sqrt(d$x), col=cols[mot])
	}
	legend(xlim[1],ylim[2], legend=mots, col=unlist(cols[mots]), pch=16, cex=0.8, bty='n')
	screen(2)
	par(mar=mar)
	aas <- charlist('ACDEFGHIKLMNPQRSTVWY')
	mots <- aas
	#rc <- corr.list(dat, 'total.length', p.0(mots,'count'))
	cts.med <- apply(dat[,p.0(mots,'count')]/dat$total.length, 2, median, na.rm=T)
	cts.sd <- apply(dat[,p.0(mots,'count')]/dat$total.length, 2, sd, na.rm=T)
	slopes <- sapply(lms[mots], function(g){s <- get.slope(g); s[1]})
	slopes.err <- sapply(lms[mots], function(g){s <- get.slope(g); s[2]})
	slopes.lower <- slopes-slopes.err
	slopes.upper <- slopes+slopes.err
	cts.lower <- (cts.med-cts.sd)
	cts.upper <- (cts.med+cts.sd)
	plot.err(cts.med, slopes, x.lower=cts.lower, x.upper=cts.upper, y.lower=slopes.lower, y.upper=slopes.upper, 
		yaxs='i', pch=16, col=unlist(all.cols[mots]), las=1, xlim=c(min(cts.lower),max(cts.upper)), ylim=c(min(slopes.lower),max(slopes.upper)),
		xlab='Median proportion of domain', ylab='Slope', log='')
	abline(h=1, lty='dotted', col='gray')
	abline(h=0, lty='dashed', col='gray')
	text(cts.med, slopes, label=mots, col=unlist(all.cols[mots]), pos=4, adj=c(0,-0.5))
	screen(3)
	par(mar=mar)
	mots <- compound.mots
	cts.med <- apply(dat[,p.0(mots,'count')]/dat$total.length, 2, median, na.rm=T)
	cts.sd <- apply(dat[,p.0(mots,'count')]/dat$total.length, 2, sd, na.rm=T)
	slopes <- sapply(lms[mots], function(g){s <- get.slope(g); s[1]})
	slopes.err <- sapply(lms[mots], function(g){s <- get.slope(g); s[2]})
	slopes.lower <- slopes-slopes.err
	slopes.upper <- slopes+slopes.err
	cts.lower <- (cts.med-cts.sd)
	cts.upper <- (cts.med+cts.sd)
	xlim <- c(min(cts.lower),max(cts.upper))
	plot.err(cts.med, slopes, x.lower=cts.lower, x.upper=cts.upper, y.lower=slopes.lower, y.upper=slopes.upper, 
		pch=16, col=unlist(all.cols[mots]), las=1, xlim=xlim, ylim=c(min(slopes.lower),max(slopes.upper)),
		xlab='Median proportion of domain', ylab='Slope', log='')
	abline(h=1, lty='dotted', col='gray')
	abline(h=0, lty='dashed', col='gray')
	text(cts.med, slopes, label=mots, col=unlist(all.cols[mots]), cex=1, pos=4, adj=c(0,-0.5))
	screen(4)
	par(mar=mar)
	bins <- seq(0,max(dat$FY.count)+1)
	h <- hist(dat$FY.count, breaks=bins-0.5, plot=F)
	barplot(h$count, las=1, col='gray', names.arg=bins[1:(length(bins)-1)], ylim=c(0,max(h$counts)), ylab='Number of FY occurrences')

	# Residuals. Fluctuations up in KR matched by fluctuations down in AGPNQ
	screen(5)
	par(mar=mar)
	mots <- compound.mots
	resid <- lapply(lms, residuals)
	fits <- lapply(lms, fitted)
	aa1 <- 'G'
	aa2 <- 'KR'
	nres.x <- resid[[aa1]]/fits[[aa1]]
	nres.y <- resid[[aa2]]/fits[[aa2]]
	plot(nres.x, nres.y, pch=16)
	abline(0,-1)
	#pairn <- combn(mots,2)
	#corrs <- apply(pairn,2,function(x){
	#	cortest(resid[[x[1]]]/fits[[x[1]]], resid[[x[2]]]/fits[[x[2]]], meth='p')$estimate
	#	})
	#names(corrs) <- apply(pairn, 2, paste, collapse='-')
	#pcor(cortest(nres.x,nres.y,meth='s'))

	# Selection strength: vs. Poisson
	screen(6)
	par(mar=mar)
	mots <- c("FY",'KR',"AGNPQ","ILMVW")
	dat <- fracs
	xlim <- c(max(1,min(raw.counts$total.length)),max(raw.counts$total.length))
	ylim <- c(0,1)
	#xlim <- c(0,300)
	#ylim <- c(0,50)
	#xlim <- c(40,260)
	plot(1,1, type='n', log='', xlim=xlim, ylim=ylim, las=1, pch=16,
		xlab='Total domain length', ylab='Fraction')
	for (xi in 1:length(mots)) {
		mot <- mots[xi]
		g <- frac.lms[[mot]]
		# Background Poisson variability
		#d <- data.frame(x=seq(1,max(xlim)*1.1,1))
		#pred <- predict.lm(g, d)
		#err <- 1.96*sqrt(abs(pred))
		#poly.err(x=d$x, y=NULL, y.lower=pred-err, y.upper=pred+err, lines=F, col=tcol(cols[mot],0.2))

		# The points
		d <- data.frame(x=dat$total.length, y=dat[,p.0('f',mot)])
		col <- tcol(cols[[mot]], 0.7)
		points(d$x, d$y, pch=16, col=col)

		# Lines
		d <- subset(d, x>0 & y>0)
		xbase <- seq(min(xlim), max(xlim), length.out=100)
		if (log.fit) {
			y <- exp(coef(g)[1])*xbase^(coef(g)[2])
			lines(xbase, y, col=col, lwd=2)
		} else {
			abline(g, col=col, lwd=2)
		}

		#lines(d$x, predict.lm(lms[[mot]], d)+sqrt(d$x), col=cols[mot])
	}
	legend(xlim[1],ylim[2], legend=mots, col=unlist(cols[mots]), pch=16, cex=0.8, bty='n')
	close.screen(all=TRUE)	
	if (file.out) dev.off()

	#d <- dat[,p.0(all.mots,'count')]/dat$total.length
	#rc <- rcormat(d)
}

if (fig.liquid.slopes) {
	resid <- lapply(lms, residuals)
	fits <- lapply(lms, fitted)
	mar <- c(4,4,1,1)
	mots <- c('AGNPQ','KR')
	#mots <- c('F','Y')
	#mots <- c('MI','VL')
	if (file.out) dev.out(p.0("fig",domain,"aa.slopes"), width=8, height=8, output.type=output.type)
	split.screen(c(1,2))
	screen(1)
	par(mar=mar)
	plot(1, 1, type='n', xlim=c(0,max(counts$total.length)), ylim=c(-1,1), xlab='Length', ylab='Fraction')
	for (xi in 1:length(mots)) {
		mot <- mots[xi]
		fmot <- resid[[mot]]/fits[[mot]] #counts[,p.0(mot,'count')]
		points(counts$total.length, fmot, col=cols[xi], pch=16)
	}
	screen(2)
	par(mar=mar)
	d <- sapply(mots, function(mot) {resid[[mot]]/fits[[mot]]})
	plot(d[,1], d[,2], pch=16)
	abline(0,-1)
	pcor(cortest(d[,1],d[,2],meth='s'))

	#multidens(list(counts$AGNPQ, counts$KR, counts$AGNPQ/counts$KR), log=T, legend.at=c(10,1))
	#for (mot in mots) {
	#	
	#}
	#plot(counts$total.length, fracs$f.ILMVW)
	close.screen(all=TRUE)	
	if (file.out) dev.off()
}

mot.freq <- function(dat, mot, inds) {
	y <- rowSums(x <- dat[[mot]][,1+inds])
	y
}

if (fig.liquid.aadist) {
	#liquid.domain <- 'nsp1-fgrepeat'
	#structured.domain <- 'nsp1-struct'
	liquid.domain <- 'pab1p'
	structured.domain <- 'pab1rrm'
	#liquid.domain <- 'nup57-fgrepeat'
	#structured.domain <- 'nup57-struct'

	all.mots <- c("FY",'KR','GKR','AGNPQ','ILMVW')
	dist.dat <- lapply(all.mots, function(m) {
		df <- read.table(paste("../data/",liquid.domain,"-",m,"-distances.txt",sep=''), header=T)
		df <- subset(df, distance>0)
		rownames(df) <- df$distance
		df
		})
	names(dist.dat) <- all.mots
	struct.dat <- lapply(all.mots, function(m) {
		df <- read.table(paste("../data/",structured.domain,"-",m,"-distances.txt",sep=''), header=T)
		df <- subset(df, distance>0)
		rownames(df) <- df$distance
		df
		})
	names(struct.dat) <- all.mots

	summarize.distance.proportions <- function(x) {
		list(means=apply(x, 1, mean, na.rm=T), errs=apply(x, 1, sderr, na.rm=T))
	}

	# Calculate proportions for each species group
	min.species.count = 5
	spec.groups <- unique(subset(spec.data, group1.count>=min.species.count)$group1)
	merge.distance.proportions <- function(x, motifs) {
		res <- lapply(motifs, function(m) {
			d <- x[[m]]
			# For each motif group...
			# ...summarize across taxonomic groups.
			rg <- lapply(spec.groups, function(g) {
				key <- p.0(subset(spec.data, group1==g)$short.name,m,'count')
				z <- match(key, colnames(d))
				y <- d[,na.omit(z),drop=FALSE]
				summarize.distance.proportions(y)
				})
			names(rg) <- spec.groups
			rg[['all']] <- summarize.distance.proportions(d)
			rg$distance <- d$distance
			rg
			})
		names(res) <- motifs
		res
	}
	liquid.dist <- merge.distance.proportions(dist.dat, all.mots)
	struct.dist <- merge.distance.proportions(struct.dat, all.mots)
	dats <- list(liquid.dist, struct.dist)

	cols <- c(myrainbow(length(spec.groups)),'black')
	names(cols) <- c(spec.groups,'all')

	xlim <- c(0.9,21)
	data.range <- seq(xlim[1],xlim[2],1)
	plot.species.groups <- F

	if (file.out) dev.out(p.0("fig",liquid.domain,"aa.distances"), width=8, height=6, output.type=output.type)
	split.screen(c(2,length(all.mots)))
	for (screen.row in 0:1) {
		cur.dat <- dats[[screen.row+1]]
		for (xi in 1:length(all.mots)) {
			mot <- all.mots[[xi]]
			x <- cur.dat[[mot]]
			sind <- screen.row*length(all.mots)+xi
			#print(paste(mot, xi, sind))
			screen(sind)
			par(mar=c(4,0,1,0.001))
			maxy <- 1
			props <- cur.dat[[mot]][['all']]
			means <- props$means[data.range]
			errs <- props$errs[data.range]
			maxm <- max(means)
			#print(paste(mot,min(means),xi))
			plot(1, 1, type='n', log='', xlim=xlim, xaxs='i', ylim=c(0.0,1.1*maxm), main=mot, yaxt='n', xlab='Distance in AA', ylab='')
			lines.err(x$distance[data.range], means, y.lower=(means - errs), y.upper=(means + errs), col='black', lwd=2)
			if (plot.species.groups) {
				for (g in spec.groups) {
					props <- x[[g]]
					means <- props$means[data.range]
					errs <- props$errs[data.range]
					maxm <- max(props$means)
					lines.err(x$distance[data.range], means/maxm, y.lower=(means - errs)/maxm, y.upper=(means + errs)/maxm, col=tcol(cols[[g]],0.3))
				}
			}
			#if (xi==1) legend(1,maxm, legend=c(spec.groups,'all'), lwd=c(rep(1,length(spec.groups)),2), col=c(unlist(cols[spec.groups]),'black'), cex=0.5, bty='n')
		}
	}
	close.screen(all=TRUE)
	if (file.out) dev.off()
}

if (fig.pabp.sep.dist) {
	post <- 'FY'
	pickets <- c('M','N','G','Q','P','PQ','PQK','PQKS','PQKNESDG','ACDEGHIKLMNPQRSTVW')
	sep.dat <- lapply(pickets, function(m) {
		df <- read.table(paste("../data/pab1p-",post,"_",m,"-separator-distributions.txt",sep=''), header=T)
		#df <- subset(df, distance>0)
		#rownames(df) <- df$distance
		df
		})
	names(sep.dat) <- pickets

	summarize.separators <- function(x) {
		list(means=apply(x, 2, mean), errs=apply(x, 2, sderr))
	}

	# Calculate proportions for each species group
	min.species.count = 5
	spec.groups <- unique(subset(spec.data, group1.count>=min.species.count)$group1)
	max.dist <- 20

	merge.separators <- function(x, motifs) {
		res <- lapply(motifs, function(m) {
			d <- x[[m]]
			#print(d)
			# For each motif group...
			# ...summarize across taxonomic groups.
			dist.flds <-  c(p.0('all',0:20),p.0(m,0:20))
			all.flds <- c(p.0(post,'count'), p.0(m,'count'), dist.flds)
			rg <- lapply(spec.groups, function(g) {
				z <- match(subset(spec.data, group1==g)$short.name, d$id)
				y <- d[z,all.flds]
				summarize.separators(y)
				})
			names(rg) <- spec.groups
			rg[['all']] <- summarize.separators(d[,all.flds])
			rg
			})
		names(res) <- motifs
		res
	}
	p.sep <- merge.separators(sep.dat, pickets)

	cols <- c(myrainbow(length(spec.groups)),'black')
	names(cols) <- c(spec.groups,'all')

	ds <- 0:max.dist
	n <- length(pickets)

	if (F) {
		if (file.out) dev.out("fig.pabp.separator.distances", width=8, height=6, output.type=output.type)
		split.screen(c(1,n))
		for (xi in 1:n) {
			screen(xi)
			par(mar=c(4,1,1,0.001))
			mot <- pickets[xi]
			plot(1,1, type='n', main=mot, xlim=c(0,max.dist), ylim=c(0,5), xlab=paste('Residues separating',post), ylab='Average count')
			for (g in spec.groups) {
				d <- p.sep[[mot]][[g]]$means
				lines(ds, d[p.0('all',ds)], col=tcol('red',0.3))
				lines(ds, d[p.0(mot,ds)], col=tcol('blue',0.3))
			}
			d <- p.sep[[mot]][['all']]$means
			lines(ds, d[p.0('all',ds)], col='red', lwd=2)
			lines(ds, d[p.0(mot,ds)], col='blue', lwd=2)
		}
		close.screen(all=TRUE)
		if (file.out) dev.off()	
	}

	if (file.out) dev.out("fig.pabp.separators", width=8, height=4, output.type=output.type)
	nspecs <- length(spec.groups)
	split.screen(c(1,nspecs))
	for (xi in 1:nspecs) {
		screen(xi)
		spec <- spec.groups[xi]
		wwo <- sapply(pickets, function(mot){
			d <- p.sep[[mot]][[spec]]$means
			baseline <- d['all.0']
			with.only.seps <- d[p.0(mot,0)]
			#d[c(p.0(mot,0),'all.0')]
			with.only.seps - baseline
			})
		names(wwo) <- pickets
		wwo <- rev(wwo)
		mids <- barplot(wwo, beside=T, main=spec, las=1, horiz=T, names.arg=NA, xlim=c(0,6), xlab='FY clusters remaining after adding')
		text(wwo, mids, names(wwo),pos=4)
	}
	close.screen(all=TRUE)
	if (file.out) dev.off()
}
# Norm
aa.comp.vector <- function(sequen, aas=all.aas, norm=TRUE) {
	res <- t(sapply(sequen, function(s) {
		tt <- (table(charlist(s))[aas])/nchar(s)
		#names(tt) <- all.aas
		tt[is.na(tt)] <- 0
		if (norm) tt <- normv(tt)
		tt}, USE.NAMES=FALSE))
	colnames(res) <- aas
	res
}

if (fig.pab1.sliding) {
	species <- c('scer','opar','hsap','dmel','xlae')
	#species <- c('dmel','xlae')
	profs <- lapply(species, function(s) {
		read.table(paste("~/research/projects/gelscore/data/pab1-",s,"-slidecomp.txt",sep=''), header=T)
		})
	names(profs) <- species
	print(names(profs))
	cols <- myrainbow(length(species))
	names(cols) <- species
	lwd <- 1
	if (file.out) dev.out("fig.scer.regions", width=5, height=5, output.type=output.type)
	plot(1,1, type='n', xlim=c(1,800), ylim=c(0.0,1), xaxs='i', las=1, xlab='Amino acid position', ylab='Motif similarity score')
	abline(h=0.75)
	for (s in species) {
		x <- profs[[s]]
		lines(x$pos, x$score, lwd=lwd, col=cols[[s]])
	}
	legend(10,1, legend=species, col=cols, lwd=lwd, cex=0.5)
	if (file.out) dev.off()
}

if (fig.pabp.similarity) {
	x <- read.table("../data/scer-pab1p-matches.txt", header=T)
	z <- match(yres$bg$orf, x$orf)
	d <- data.frame(yres$sL[,c('orf','gene','mrna','prot')], x[z,2:ncol(x)])
	xa <- subset(d, n.above>0)
	hist(xa$n.above, breaks=seq(1,25))
	
	if (file.out) dev.out("fig.p.similarity", width=5, height=5, output.type=output.type)
	plot(1,1, type='n', xlim=c(1,800), ylim=c(0.0,1), xaxs='i', las=1, xlab='Amino acid position', ylab='Motif similarity score')
	for (s in species) {
		x <- profs[[s]]
		lines(x$pos, x$score, lwd=lwd, col=cols[[s]])
	}
	if (file.out) dev.off()
}

if (fig.droplets) {
	#liquid.domain <- 'pab1p'
	liquid.domain <- 'nsp1-fg'
	x <- read.table(paste("../data/",liquid.domain,"-droplet-structure.txt",sep=''), header=T)
	species <- unique(x$id)
	n.spec <- length(species)
	aas <- unique(x$aas)
	n.aas <- length(aas)
	cols <- myrainbow(n.aas)
	mar <- c(4,4,1,1)
	side.dim <- ceiling(sqrt(n.spec+1))
	dim.per.droplet <- 6 # nm
	side.length <- side.dim * dim.per.droplet

	droplet.pos <- function(i, dim.per.droplet, side.dim) {
		# Return the center of droplet i
		xpos = ((i-1)%%side.dim + 0.5)*dim.per.droplet
		ypos = side.dim*dim.per.droplet - (floor((i-1)/side.dim)+0.5)*dim.per.droplet
		c(xpos,ypos)
	}

	if (file.out) dev.out(p.0("fig",liquid.domain,"droplets"), width=15, height=15, output.type=output.type)
	#split.screen(c(screen.dim,screen.dim))
	par(mar=mar)
	plot(1,1, type='n', xlim=c(0,side.length), ylim=c(0,side.length), xaxs='i', yaxs='i', las=1, xlab='Nanometers', ylab='Nanometers')
	cent <- droplet.pos(1, dim.per.droplet, side.dim)
	# Reorder from outside to inside
	sent <- subset(x, id=="P.falciparum")
	if (nrow(sent)==0) {
		sent <- subset(x, id==unique(x$id)[1])
	}
	sentinel.radii <- sent$outer.radius
	#ord <- order(sentinel.radii, decreasing=T)
	#cols <- cols[ord]
	#aas <- aas[ord]
	legend(0,side.length, legend=aas, col=cols, fill=cols, cex=0.7, bty='n')
	for (xi in 1:n.spec) {
		cent <- droplet.pos(xi+1, dim.per.droplet, side.dim)
		y <- subset(x, id==species[xi])
		radii <- y$outer.radius
		# Draw droplet layers
		for (ri in length(radii):1) {
			r <- radii[ri]
			draw.ellipse(cent[1],cent[2], r, r, col=cols[ri], border=NA)
		}
		text(cent[1],cent[2]+max(radii)+0.5, species[xi], cex=1, font=3)
	}
	close.screen(all=TRUE)
	if (file.out) dev.off()
}

if (fig.compare.droplets) {
	liquid.domains <- c('pab1p','pbp1-pb','nsp1-fg','nup57-fg')
	liquid.domains <- c('pab1p','nsp1-fg')
	#liquid.domains <- c('pab1p','nsp1-fg')
	domains <- lapply(liquid.domains, function(f) {
		read.table(paste("../data/",f,"-droplet-structure.txt",sep=''), header=T)
		})
	names(domains) <- liquid.domains

	# Find shared ids
	all.ids <- unique(unlist(sapply(domains, function(d){d$id})))
	all.aas <- unique(unlist(sapply(domains, function(d){d$aas})))
	id.data <- data.frame(id=all.ids, sapply(domains, function(d){d[match(all.ids,d$id),'id']}))
	species <- as.character(na.omit(id.data)$id) #as.character(na.omit(unique(x2[match(x1$id, x2$id),'id'])))
	n.spec <- length(species)
	aas <- unique(domains[[1]]$aas)
	n.aas <- length(aas)
	cols <- myrainbow(n.aas)
	# FY	670067
	# KR	5FB7BB
	# ILMV	228855
	# AGNPQ	CFC169
	# DE	0000CC
	# ST	CE7377
	# HCW	83C9E8
	col.list <- list(FY="#670067",KR="#5FB7BB",ILMV="#228855",AGNPQ="#CFC169",DE="#0000CC",ST="#CE7377",HCW="#83C9E8")
	cols <- unlist(col.list)

	mar <- c(4,4,1,1)
	side.dim <- ceiling(sqrt(n.spec+1))
	dim.per.droplet <- 6 # nm
	side.length <- side.dim * dim.per.droplet

	droplet.pos <- function(i, dim.per.droplet, side.dim) {
		# Return the center of droplet i
		xpos = ((i-1)%%side.dim + 0.5)*dim.per.droplet
		ypos = side.dim*dim.per.droplet - (floor((i-1)/side.dim)+0.5)*dim.per.droplet
		c(xpos,ypos)
	}

	keys <- as.character(sapply(species, p.0, all.aas))
	ordered.domains <- lapply(domains, function(d){
		d[match(keys, sapply(1:nrow(d),function(m){p.0(d[m,'id'],d[m,'aas'])})),]
		})
	names(ordered.domains) <- liquid.domains
	domain.lengths <- sapply(ordered.domains,function(d){subset(d,aas==all.aas[1])$length})
	rownames(domain.lengths) <- species

	dom.corrs <- sapply(all.aas, function(aa) {lt(rcormat(sapply(ordered.domains,function(d){subset(d,aas==aa)$num.aas}))$r)})
	dom.corrs.l <- lt(rcormat(domain.lengths)$r)
	dom.corrs.f <- sapply(all.aas, function(aa) {lt(rcormat(sapply(ordered.domains,function(d){subset(d,aas==aa)$num.aas/subset(d,aas==aa)$length}))$r)})

	base.dim <- 4
	if (file.out) dev.out(p.0("fig.droplets.comparison"), width=base.dim*(length(ordered.domains	)+1), height=base.dim, output.type=output.type)
	split.screen(c(1,length(ordered.domains)+1))
	screen(1)
	# Make legend
	par(mar=mar)
	plot(1,1, type='n', xlim=c(0,side.length), ylim=c(0,side.length), xaxs='i', yaxs='i', las=1, xlab='Nanometers', ylab='Nanometers')
	legend(0,side.length, legend=aas, col=cols, fill=cols, cex=0.7, bty='n')
	# Plot droplets
	for (di in 1:length(ordered.domains)) {
		screen(di+1)
		par(mar=mar)
		x <- ordered.domains[[di]]
		plot(1,1, type='n', xlim=c(0,side.length), ylim=c(0,side.length), xaxs='i', yaxs='i', las=1, xlab='Nanometers', ylab='Nanometers')
		cent <- droplet.pos(1, dim.per.droplet, side.dim)
		# Reorder from outside to inside
		sent <- subset(x, id=="S.cerevisiae")
		if (nrow(sent)==0) {
			sent <- subset(x, id==unique(x$id)[1])
		}
		sentinel.radii <- sent$outer.radius
		#ord <- order(sentinel.radii, decreasing=T)
		#cols <- cols[ord]
		#aas <- aas[ord]
		#legend(0,side.length, legend=aas, col=cols, fill=cols, cex=0.7, bty='n')
		for (xi in 1:n.spec) {
			cent <- droplet.pos(xi, dim.per.droplet, side.dim)
			y <- subset(x, id==species[xi])
			radii <- y$outer.radius
			# Draw droplet layers
			for (ri in length(radii):1) {
				r <- radii[ri]
				draw.ellipse(cent[1],cent[2], r, r, col=cols[ri], border=NA)
			}
			#text(cent[1],cent[2]+max(radii)+0.5, species[xi], cex=1, font=3)
		}
	}
	close.screen(all=TRUE)
	#grid.export("test.svg")
	if (file.out) dev.off()
}

if (fig.layers) {
	x <- read.table("../data/pab1p-droplet-structure.txt", header=T)
	species <- unique(x$id)
	n.spec <- length(species)
	aas <- unique(x$aas)
	n.aas <- length(aas)
	cols <- myrainbow(n.aas)
	mar <- c(4,4,1,1)

	if (file.out) dev.out("fig.pabp.layers", width=n.spec*3, height=3, output.type=output.type)
	par(mar=mar)
	fld <- 'num.aas'
	ylim <- c(1,max(x[,fld]))
	plot(1,1, type='n', log='y', xlim=c(0,300), ylim=ylim, xaxs='i', yaxs='r', las=1, xlab='Length', ylab='Thickness (nm)')
	for (xi in 1:n.aas) {
		y <- subset(x, aas==aas[xi])
		points(y$length, y[,fld], pch=16, col=cols[xi])
		abline(lm(log(num.aas)~length,data=y), col=cols[xi])
	}
	legend(10,max(ylim), legend=aas, pch=16, col=cols, bty='n')
	close.screen(all=TRUE)
	if (file.out) dev.off()

}

if (fig.pab1.pbp1.comp.covariation) {
	p1 <- 'pab1'
	p2 <- 'pbp1'
	# Putative interacting domains
	x1 <- read.table("../data/pab1p-composition.txt", header=T)
	x2 <- read.table("../data/pbp1-pb-composition.txt", header=T)
	# Controls
	x1c <- read.table("../data/pab1rrm-composition.txt", header=T)
	x2c <- read.table("../data/pbp1-lsmad-composition.txt", header=T)
	z <- match(x1$id, x2$id)
	colnames(x1) <- p.0(p1,colnames(x1))
	colnames(x2) <- p.0(p2,colnames(x2))
	zc <- match(x1c$id, x2c$id)
	colnames(x1c) <- p.0(p1,colnames(x1c))
	colnames(x2c) <- p.0(p2,colnames(x2c))

	all.mots <- as.character(sapply(names(x1), function(x){substring(x, 6, nchar(x)-6)})[3:ncol(x1)-1])
	d <- na.omit(data.frame(x1, x2[z,]))
	dc <- na.omit(data.frame(x1c, x2c[z,]))
	compound.mots <- all.mots[sapply(all.mots, nchar)>1]
	compound.mots <- c("FY",'KR',"AGNPQ","ILMVW")
	#compound.mots <- c(c("FY",'KR',"AGNPQ","ILMVW",'MI','VL','DE','ST','G'))
	# Absolute numbers
	# Fractions
	#frac1
	mots <- compound.mots
	#mots <- charlist('ACDEGHIKLMNPQRSTVW')
	flds1 <- p.0(p1,mots,'count')
	flds2 <- p.0(p2,mots,'count')
	rc.abs <- rcormat(d[,c(flds1,flds2)])$r[flds1,flds2]
	rc.abs.c <- rcormat(dc[,c(flds1,flds2)])$r[flds1,flds2]
	rc.frac <- rcormat(df <- data.frame(
		apply(d[,flds1],2,function(x){x/d[,p.0(p1,'total.length')]}),
		apply(d[,flds2],2,function(x){x/d[,p.0(p2,'total.length')]})))$r[flds1,flds2]
	rc.frac.c <- rcormat(dfc <- data.frame(
		apply(dc[,flds1],2,function(x){x/dc[,p.0(p1,'total.length')]}),
		apply(dc[,flds2],2,function(x){x/dc[,p.0(p2,'total.length')]})))$r[flds1,flds2]
	df <- data.frame(id=d[,p.0(p1,'id')], df)
	dfc <- data.frame(id=d[,p.0(p1,'id')], dfc)

	xlim <- ylim <- c(-1,1)
	mar <- c(4,4,1,1)

	if (file.out) dev.out("fig.pab1.pbp1.comp.covariation", width=6, height=6, output.type=output.type)
	split.screen(c(2,2))
	screen(1)
	par(mar=mar)
	#plot(rc.abs.c, rc.abs, xlim=xlim, ylim=ylim, pch=16, las=1, xlab='Composition correlation (non-interacting domains)', ylab='Composition correlation (interacting domains)')
	plot(l1 <- d[,p.0(p1,'total.length')], l2<- d[,p.0(p2,'total.length')], pch=16, las=1, xlab='Length of interacting domain 1', ylab='Length of interacting domain 2')
	abline(lm(l2~l1))
	screen(2)
	par(mar=mar)
	plot(l1 <- dc[,p.0(p1,'total.length')], l2<- dc[,p.0(p2,'total.length')], pch=16, las=1, xlab='Length of non-interacting domain 1', ylab='Length of non-interacting domain 2')
	abline(lm(l2~l1))
	screen(3)
	par(mar=mar)
	plot(rc.frac.c, rc.frac, xlim=xlim, ylim=ylim, pch=16, las=1, xlab='Composition correlation (non-interacting domains)', ylab='Composition correlation (interacting domains)')
	abline(0,1)
	screen(4)
	par(mar=mar)
	plot(rc.frac.c^2, rc.frac^2, xlim=c(0,1), ylim=c(0,1), pch=16, las=1, xlab='Composition correlation (non-interacting domains)', ylab='Composition correlation (interacting domains)')
	abline(0,1)
	close.screen(all=TRUE)
	if (file.out) dev.off()
}

if (fig.pab1.pbp1.hydro.covariation) {
	p1 <- 'pab1'
	p2 <- 'pbp1'
	# Putative interacting domains
	x1 <- read.table("../data/pab1p-composition.txt", header=T)
	x2 <- read.table("../data/pbp1-pb-composition.txt", header=T)
	# Controls
	x1c <- read.table("../data/pab1rrm-composition.txt", header=T)
	x2c <- read.table("../data/pbp1-lsmad-composition.txt", header=T)
	shared.specs <- as.vector(na.omit(x1[match(x2$id,x1$id),]$id))
	z1 <- match(shared.specs, x1$id)
	z2 <- match(shared.specs, x2$id)
	shared.specs.c <- as.vector(na.omit(x1c[match(x2c$id,x1c$id),]$id))
	z1c <- match(shared.specs, x1c$id)
	z2c <- match(shared.specs, x2c$id)


	aa.group <- charlist('ILVM')
	logo <- function(x, aas) {
		# Compute the log ratios of aas
		res <- apply(combn(aas,2),2,function(y){
			log.nz(x[,p.0(y[1],'count')]/x[,p.0(y[2],'count')])
			})
		colnames(res) <- p.0(pair.names(aas,sep='.'),'logr')
		res
	}

	lo1 <- logo(x1[z1,],aa.group)
	lo1.nm <- p.0(p1,colnames(lo1))
	colnames(lo1) <- lo1.nm
	lo2 <- logo(x2[z2,],aa.group)
	lo2.nm <- p.0(p2,colnames(lo2))
	colnames(lo2) <- lo2.nm
	d <- data.frame(id=x1[z1,]$id, lo1, lo2)

	# Control log odds ratios
	lo1c <- logo(x1c[z1c,],aa.group)
	colnames(lo1c) <- lo1.nm
	lo2c <- logo(x2c[z2c,],aa.group)
	colnames(lo2c) <- lo2.nm
	dc <- data.frame(id=x1c[z1c,]$id, lo1c, lo2c)

	hyd <- read.table("~/research/lib/data/chimera-hydrophobicities.txt", header=T)
	rownames(hyd) <- hyd$aa
	hyd.fld <- 'hh.hyd'
	hyd.diffs <- apply(combn(aa.group,2),2,function(y){hyd[y[1],hyd.fld]-hyd[y[2],hyd.fld]})


	flds1 <- lo1.nm
	flds2 <- lo2.nm
	rc.abs <- rcormat(d[,c(flds1,flds2)])$r[flds1,flds2]
	rc.abs.c <- rcormat(dc[,c(flds1,flds2)])$r[flds1,flds2]

	xlim <- ylim <- c(-1,1)
	mar <- c(4,4,1,1)

	thermo.list <- c('C.thermophilum','C.thermophila','T.terrestris','T.laguninosus','T.thermophilus','T.stellatus','A.fumigatus')
	meso.list <- c('S.cerevisiae','C.globosum','A.pullulans','P.stipitis')
	psychro.list <- c('S.kudriavzevii','G.destructans')

	# Coe
	l1 <- sapply(1:nrow(d), function(m){coef(lm(as.numeric(d[m,lo1.nm])~hyd.diffs))[2]})
	l2 <- sapply(1:nrow(d), function(m){coef(lm(as.numeric(d[m,lo2.nm])~hyd.diffs))[2]})
	l1c <- sapply(1:nrow(dc), function(m){coef(lm(as.numeric(dc[m,lo1.nm])~hyd.diffs))[2]})
	l2c <- sapply(1:nrow(dc), function(m){coef(lm(as.numeric(dc[m,lo2.nm])~hyd.diffs))[2]})
	# Correlations
	l1 <- sapply(1:nrow(d), function(m){cor(as.numeric(d[m,lo1.nm]),hyd.diffs,method='s')})
	l2 <- sapply(1:nrow(d), function(m){cor(as.numeric(d[m,lo2.nm]),hyd.diffs,method='s')})
	l1c <- sapply(1:nrow(d), function(m){cor(as.numeric(dc[m,lo1.nm]),hyd.diffs,method='s')})
	l2c <- sapply(1:nrow(d), function(m){cor(as.numeric(dc[m,lo2.nm]),hyd.diffs,method='s')})
	res <- data.frame(id=d$id, l1=l1,l2=l2,l1c=l1c,l2c=l2c)

	if (file.out) dev.out("fig.pab1.pbp1.comp.covariation", width=6, height=6, output.type=output.type)
	split.screen(c(2,2))
	screen(1)
	par(mar=mar)
	plot(c(l1,l1c),c(l2,l2c), type='n')
	points(l1,l2,pch=16,col='red')
	points(l1c,l2c,pch=16,col='blue')
	therm.sub <- res[match(thermo.list,res$id),]
	meso.sub <- res[match(meso.list,res$id),]
	psychro.sub <- res[match(psychro.list,res$id),]
	with(therm.sub, points(l1,l2,pch=16,col='darkred',cex=2))
	with(meso.sub, points(l1,l2,pch=16,col='red',cex=2))
	with(psychro.sub, points(l1,l2,pch=16,col='pink',cex=2))
	#points(l1c,l2c,pch=16,col='blue')
	#plot(rc.abs.c, rc.abs, xlim=xlim, ylim=ylim, pch=16, las=1, xlab='Composition correlation (non-interacting domains)', ylab='Composition correlation (interacting domains)')
	abline(0,1)
	screen(2)
	# Plot extremes.
	#plot()
	close.screen(all=TRUE)
	if (file.out) dev.off()
}

if (fig.pab1.hydro.thermo.vs.meso) {
	# Do thermophiles use more hydrophobic amino acids?

	x1 <- read.table("../data/pab1p-composition.txt", header=T)
	#x1 <- read.table("../data/pbp1-pb-composition.txt", header=T)
	thermo.list <- c('C.thermophilum','C.thermophila','T.terrestris','L.thermotolerans','T.lanuginosus','T.thermophilus','T.stellatus','O.parapolymorpha')#,'K.marxianus')
	meso.list <- c('S.cerevisiae','S.stipitis','C.glabrata')
	psychro.list <- c('S.kudriavzevii','A.pullulans')#,'G.destructans')
	temp.gps <- list(thermo=thermo.list, meso=meso.list, psychro=psychro.list)
	all.gps <- as.vector(unlist(temp.gps))


	aa.list <- charlist('ILVMA')#AGP')

	x.all <- x1[match(all.gps,x1$id),c('id',p.0(aa.list,'count'))]
	xa <- sapply(aa.list, function(m){
		r <- sapply(names(temp.gps), function(g){
			y <- x1[match(temp.gps[[g]],x1$id),p.0(m,'count')]
			mean(y,na.rm=T)
			})
		r
		})
	xad <- sapply(aa.list, function(m){
		r <- sapply(names(temp.gps), function(g){
			y <- x1[match(temp.gps[[g]],x1$id),p.0(m,'count')]
			sd(y,na.rm=T)/sqrt(length(y))
			})
		r
		})
	xf <- sapply(aa.list, function(m){
		r <- sapply(names(temp.gps), function(g){
			y <- x1[match(temp.gps[[g]],x1$id),]
			yr <- y[,p.0(m,'count')]/y[,'total.length']
			mean(yr,na.rm=T)
			})
		r
		})
	xfd <- sapply(aa.list, function(m){
		r <- sapply(names(temp.gps), function(g){
			y <- x1[match(temp.gps[[g]],x1$id),]
			yr <- y[,p.0(m,'count')]/y[,'total.length']
			sd(yr,na.rm=T)/sqrt(length(y))
			})
		r
		})

	hyd <- read.table("~/research/lib/data/chimera-hydrophobicities.txt", header=T)
	rownames(hyd) <- hyd$aa
	hyd.fld <- 'hh.hyd'

	cols <- myrainbow(length(aa.list))
	lwd <- 2

	if (file.out) dev.out("fig.pab1.pbp1.comp.covariation", width=6, height=6, output.type=output.type)
	split.screen(c(1,2))
	screen(1)
	par(mar=mar)
	x <- hyd[match(aa.list,hyd$aa),hyd.fld]
	#x <- 1:length(aa.list)
	plot.err(x, xa[1,], y.lower=xa[1,]-xad[1,], y.upper=xa[1,]+xad[1,], las=1, pch=16, ylim=c(0,max(xa+xad,na.rm=T)), xaxt='n', xlab='Hydrophilicity', ylab='Frequency')
	for (xi in 1:nrow(xf)) {
		lines.err(x, xa[xi,], y.lower=xa[xi,]-xad[xi,], y.upper=xa[xi,]+xad[xi,], las=1, pch=16, lwd=lwd, col=cols[xi])
	}
	mtext(aa.list, side=1, at=x)

	#matplot(hyd[match(aa.list,hyd$aa),hyd.fld], t(xf), type='l', lwd=lwd, lty='solid', col=cols, xaxt='n', ylab='Proportion of sequence', las=1)
	#matplot(t(xf), type='l', lwd=lwd, lty='solid', col=cols, xaxt='n', ylab='Proportion of sequence', las=1)
	#	matplot(hyd[match(aa.list,hyd$aa),hyd.fld], t(xa), type='l', lwd=lwd, lty='solid', col=cols, xaxt='n', ylab='Number of occurrences', las=1)
	#	mtext(aa.list, side=1, at=hyd[match(aa.list,hyd$aa),hyd.fld])
	#barplot(xa, horiz=T, beside=T)
	screen(2)
	par(mar=mar)
	x <- hyd[match(aa.list,hyd$aa),hyd.fld]
	#x <- 1:length(aa.list)
	plot.err(x, xf[1,], y.lower=xf[1,]-xfd[1,], y.upper=xf[1,]+xfd[1,], las=1, pch=16, ylim=c(0,max(xf+xfd,na.rm=T)), xaxt='n', xlab='Hydrophilicity', ylab='Proportion')
	for (xi in 1:nrow(xf)) {
		lines.err(x, xf[xi,], y.lower=xf[xi,]-xfd[xi,], y.upper=xf[xi,]+xfd[xi,], las=1, pch=16, lwd=lwd, col=cols[xi])
	}
	mtext(aa.list, side=1, at=x)

	#matplot(hyd[match(aa.list,hyd$aa),hyd.fld], t(xf), type='l', lwd=lwd, lty='solid', col=cols, xaxt='n', ylab='Proportion of sequence', las=1)
	#matplot(t(xf), type='l', lwd=lwd, lty='solid', col=cols, xaxt='n', ylab='Proportion of sequence', las=1)
	#plot.err(hyd[match(aa.list,hyd$aa),hyd.fld], xfd[])
	#barplot(xf, horiz=T, beside=T)
	close.screen(all=TRUE)
	if (file.out) dev.off()
}

if (fig.coevolution) {
	colloid.domains <- list("Pab1"='pab1p',"Pbp1"='pbp1-pb',"Nsp1"='nsp1-fg',"Nup57"='nup57-fg')
	#colloid.domains <- list("Pab1"='pab1rrm',"Pbp1"='pbp1-lsmad',"Nsp1"='nsp1-struct',"Nup57"='nup57-struct')
	#liquid.domains <- c('pab1p','nsp1-fg')
	domains <- lapply(colloid.domains, function(f) {
		read.table(paste("../data/",f,"-composition.txt",sep=''), header=T)
		})
	dom.names <- names(colloid.domains)
	#names(domains) <- colloid.domains

	# Find shared ids
	all.ids <- unique(unlist(sapply(domains, function(d){d$id})))
	shared.ids <- all.ids
	for (d in domains) {
		shared.ids <- intersect(shared.ids, d$id)
	}
	ordered.domains <- lapply(domains, function(d){d[match(shared.ids,d$id),]})
	# Compute correlations for each amino acid and class
	# For each pair of domains, generate all 
	combs <- combn(dom.names, 2)
	res <- apply(combs, 2, function(m){
		x1 <- ordered.domains[[m[1]]]
		x2 <- ordered.domains[[m[2]]]
		sapply(p.0(all.aas,'count'), function(aa.id){cortest(x1[,aa.id]/x1[,'total.length'], x2[,aa.id]/x2[,'total.length'], meth='s')$estimate})
		})
	colnames(res) <- pair.names(colloid.domains)

	motifs <- c('FY','AGNPQ','KR','ILMVW')
	mot.res <- apply(combs, 2, function(m){
		x1 <- ordered.domains[[m[1]]]
		x2 <- ordered.domains[[m[2]]]
		sapply(motifs, function(mot){
			cols <- p.0(charlist(mot),'count')
			cortest(rowSums(x1[,cols])/x1[,'total.length'], rowSums(x2[,cols])/x2[,'total.length'], meth='s')$estimate})
		})
	colnames(mot.res) <- pair.names(dom.names)

	# Relative probabilities
	rel.aas <- charlist('ILVMA')
	pseudocount <- 0.5
	dom.probs <- lapply(ordered.domains, function(d){
		probs <- (d[,p.0(rel.aas,'count')]+pseudocount)/d$total.length
		rel.probs <- apply(probs, 1, function(x){
			ut(x %*% t(1/x))
			})
		res <- t(rel.probs)
		dimnames(res) <- list(shared.ids, pair.names(rel.aas,sep='.'))
		res
		})

	# Correlation between 
	usage.corrs <- apply(combn(dom.names,2), 2, function(x){
		sapply(1:nrow(dom.probs[[x[1]]]), function(m){
			cortest(dom.probs[[x[1]]][m,],dom.probs[[x[2]]][m,],method='s')$estimate
			})
		})
	colnames(usage.corrs) <- pair.names(dom.names)

	# Lengths
	dom.lengths <- sapply(ordered.domains, function(d){d$total.length})

	# Get hydrophobicity information
	hyd <- read.table("~/research/lib/data/chimera-hydrophobicities.txt", header=T)
	hyd.fld <- 'hh.hyd'
	z1 <- match(rel.aas, hyd$aa)
	h <- hyd[z1,hyd.fld]
	hv <- hyd[match(rel.aas,hyd$aa),hyd.fld]
	n <- length(rel.aas)
	dg.mat <- matrix(rep(hv,n),n,n) - t(matrix(rep(hv,n),n,n))
	hyd.ddgs <- ut(dg.mat)
	ho <- order(hyd.ddgs)

	hyd.corrs <- lapply(dom.probs, function(d) {
		apply(d, 1, cor, y=hyd.ddgs, method='p')
		})

	hyd.lm <- lapply(dom.probs, function(d) {
		apply(d, 1, function(x){lm(log(x)~hyd.ddgs)})
		})
	hyd.slopes <- lapply(hyd.lm, function(x) {
		sapply(x, function(g){coef(summary(g))[2,c(1,4)]})
		})

	mar <- c(4,4,1,1)
	lwd <- 2
	if (file.out) dev.out("fig.conservation")
	split.screen(c(2,3))
	screen(1)
	# Example of usage--hydrophobicity correlation
	spec <- 'S.cerevisiae'
	plot(1, 1, type='n', las=1, xlim=c(-0.8,0), ylim=c(-2.5,0), xlab='Hydrophilicity', ylab='Log relative usage', log='')
	cols <- myrainbow(length(dom.names))
	#for (x in enum(dom.names)) {
	#	points(hyd.ddgs[ho], log(dom.probs[[x$x]][spec,][ho]), pch=16, col=cols[x$i])
	#	abline(hyd.lm[[x$x]][[spec]], col=cols[x$i])
	#}
	points(hyd.ddgs[ho], log(dom.probs[['Pab1']][spec,][ho]), pch=16) #, col=cols[x$i])
	abline(hyd.lm[['Pab1']][[spec]])#, col=cols[x$i])
	
	screen(2)
	par(mar=mar)
	multi.ecdf(usage.corrs, lwd=lwd, legend.at=c(-1,1), xlim=c(-1,1))
	screen(3)
	par(mar=mar)
	slopes <- lapply(hyd.slopes, function(x){x[1,]})
	slope.sds <- sapply(slopes, sderr)
	slope.means <- sapply(slopes, mean)
	#barplot.err(slope.means, lower=slope.means-slope.sds, upper=slope.means+slope.sds)
	multi.ecdf(lapply(hyd.slopes, function(x){x[1,]}), log=F, lwd=lwd, xlab="Slope of frequency\u2013usage association")
	screen(5)
	par(mar=mar)
	multi.ecdf(lapply(hyd.slopes, function(x){x[2,]}), log=T, lwd=lwd, xlab="P value for frequency\u2013usage association")
	screen(6)
	par(mar=mar)
	multi.ecdf(hyd.corrs, log=F, lwd=lwd, 
		xlab="Relative usage\u2013hydrophilicity correlation")
	abline(v=sapply(hyd.corrs, median), lwd=lwd, lty='dotted', col=myrainbow(length(hyd.corrs)))
	#multi.ecdf(lapply(hyd.lm, function(x){sapply(x, function(g){summary(g)$r.squared})}), log=F, lwd=lwd, 
	#	xlab="R2 for frequency\u2013usage association")
	close.screen()
	if (file.out) dev.off()

	#dom.corrs <- sapply(all.aas, function(aa) {lt(rcormat(sapply(ordered.domains,function(d){subset(d,aas==aa)$num.aas}))$r)})
	#dom.corrs.l <- lt(rcormat(domain.lengths)$r)
	#dom.corrs.f <- sapply(all.aas, function(aa) {lt(rcormat(sapply(ordered.domains,function(d){subset(d,aas==aa)$num.aas/subset(d,aas==aa)$length}))$r)})

	panel.smooth2 <- function (x, y, col = par("col"), bg = NA, pch = par("pch"), 
	    cex = 1, col.smooth = "red", span = 2/3, iter = 3, ...) 
	{
	    points(x, y, pch = pch, col = col, bg = bg, cex = cex)
	    ok <- is.finite(x) & is.finite(y)
	    if (any(ok)) 
	        lines(stats::lowess(x[ok], y[ok], f = span, iter = iter), 
	            col = col.smooth, lwd=2, ...)
	}

	panel.cor <- function(x, y, digits=2, prefix="", cex.cor) 
	{
	    usr <- par("usr"); on.exit(par(usr)) 
	    par(usr = c(0, 1, 0, 1)) 
	    r <- cor(x, y, method='s')
	    txt <- format(c(r, 0.123456789), digits=digits)[1] 
	    txt <- paste(prefix, txt, sep="") 
	    if(missing(cex.cor)) cex <- 0.8/strwidth(txt) 
	 
	    test <- cor.test(x,y) 
	    # borrowed from printCoefmat
	    Signif <- symnum(test$p.value, corr = FALSE, na = FALSE, 
	                  cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
	                  symbols = c("***", "**", "*", ".", " ")) 
	 
	    text(0.5, 0.5, txt, cex = cex * r) 
	    #text(.8, .8, Signif, cex=cex, col=2) 
	}
	par(las=1)
	par(pch=16)
	#par(lwd=2)

	pairs(dom.lengths, lower.panel=panel.smooth2, upper.panel=panel.cor)
	#dev.off()
}

# Normalize composition by dividing by total length
norm.comp <- function(x, L='total.length', omit=c('id','total.length')) {
	flds <- colnames(x)
	q.flds <- setdiff(flds,omit)
	qx <- x[,q.flds] # fields to quantify/normalize
	q.res <- sapply(1:length(q.flds), function(m){qx[,m]/x[,L]})
	colnames(q.res) <- q.flds
	res <- data.frame(x[,omit], q.res)
	# Return normalized version of x: same columns, in same order
	res[,colnames(x)]
}

if (fig.pab1.sis1.coev) {
	p1 <- 'pab1p'
	p2 <- 'sis1-lc'
	p1c <- 'pab1rrm'
	p2c <- 'sis1-struct'
	# Putative interacting domains
	x1 <- norm.comp(read.table(paste("../data/",p1,"-composition.txt",sep=''), header=T))
	x2 <- norm.comp(read.table(paste("../data/",p2,"-composition.txt",sep=''), header=T))
	# Controls
	x1c <- norm.comp(read.table(paste("../data/",p1c,"-composition.txt",sep=''), header=T))
	x2c <- norm.comp(read.table(paste("../data/",p2c,"-composition.txt",sep=''), header=T))

	species <- intersect(x1$id,x2$id)
	z1 <- match(species, x1$id)
	z2 <- match(species, x2$id)
	z1c <- match(species, x1c$id)
	z2c <- match(species, x2c$id)

	# All amino acids
	f <- p.0(charlist('ACDEFGHIKLMNPQRSTVWY'),'count')
	# These specific amino acids yield a stronger signal
	#f <- p.0(charlist('ILMVFYR'),'count')

	species.rand <- sample(1:length(species),length(species))

	meth='s'

	# Correlations between AA usage in protein 1 vs. protein 2, in first (test) region
	cr.tt <- sapply(1:length(species), function(s){
		cor(as.numeric(x1[z1,][s,f]), as.numeric(x2[z2,][s,f]),meth=meth)})
	# Correlations between AA usage in protein 1 vs. protein 2, in second (control) region
	cr.cc <- sapply(1:length(species), function(s){
		cor(as.numeric(x1c[z1c,][s,f]), as.numeric(x2c[z2c,][s,f]),meth=meth)})
	# Correlations between AA usage in protein 1 in test vs. control region
	cr.tc1 <- sapply(1:length(species), function(s){
		cor(as.numeric(x1[z1,][s,f]), as.numeric(x1c[z1c,][s,f]),meth=meth)})
	# Correlations between AA usage in protein 2 in test vs. control region
	cr.tc2 <- sapply(1:length(species), function(s){
		cor(as.numeric(x2[z2,][s,f]), as.numeric(x2c[z2c,][s,f]),meth=meth)})
	# Correlations between AA usage in test region of protein 1 vs. control region of protein 2
	cr.t1c2 <- sapply(1:length(species), function(s){
		cor(as.numeric(x1[z1,][s,f]), as.numeric(x2c[z2c,][s,f]),meth=meth)})
	# Correlations between AA usage in test region of protein 2 vs. control region of protein 1
	cr.t2c1 <- sapply(1:length(species), function(s){
		cor(as.numeric(x2[z2,][s,f]), as.numeric(x1c[z1c,][s,f]),meth=meth)})

	n.species.rand <- 5
	# Correlations between AA usage in protein 1 vs. protein 2, in first (test) region
	srand.cr.tt <- replicate(n.species.rand,
		sapply(1:length(species), function(s){cor(as.numeric(x1[z1,][species.rand,][s,f]), as.numeric(x2[z2,][s,f]),meth=meth)}))
	# Correlations between AA usage in protein 1 vs. protein 2, in second (control) region
	srand.cr.cc <- replicate(n.species.rand,
		sapply(1:length(species), function(s){cor(as.numeric(x1c[z1c,][species.rand,][s,f]), as.numeric(x2c[z2c,][s,f]),meth=meth)}))
	# Correlations between AA usage in protein 1 in test vs. control region
	#srand.cr.tc1 <- sapply(1:length(species), function(s){
	#	cor(as.numeric(x1[z1,][species.rand,][s,f]), as.numeric(x1c[z1c,][s,f]),meth=meth)})
	# Correlations between AA usage in protein 2 in test vs. control region
	#srand.cr.tc2 <- sapply(1:length(species), function(s){
	#	cor(as.numeric(x2[z2,][species.rand,][s,f]), as.numeric(x2c[z2c,][s,f]),meth=meth)})
	# Correlations between AA usage in test region of protein 1 vs. control region of protein 2
	#srand.cr.t1c2 <- sapply(1:length(species), function(s){
	#	cor(as.numeric(x1[z1,][species.rand,][s,f]), as.numeric(x2c[z2c,][s,f]),meth=meth)})
	# Correlations between AA usage in test region of protein 2 vs. control region of protein 1
	#srand.cr.t2c1 <- sapply(1:length(species), function(s){
	#	cor(as.numeric(x2[z2,][species.rand,][s,f]), as.numeric(x1c[z1c,][s,f]),meth=meth)})

	comparisons <- list(cr.tt, cr.cc, cr.tc1, cr.tc2, cr.t1c2, cr.t2c1)
	species.rand.comparisons <- list(
		sapply(1:n.species.rand, function(m){cr.tt-srand.cr.tt[,m]}), 
		sapply(1:n.species.rand, function(m){cr.cc-srand.cr.cc[,m]}))
	pd <- function(...){paste(...,sep='\u2013')}
	names(comparisons) <- c(pd(p1,p2), pd(p1c,p2c), pd(p1,p1c), pd(p2,p2c), pd(p1,p2c), pd(p2,p1c))
	names(species.rand.comparisons) <- c(pd(p1,p2), pd(p1c,p2c)) #, pd(p1c,p2c), pd(p1c,p2c,'[randomized]'))

	mar <- c(4,4,1,1)

	if (file.out) dev.out("fig_aa_usage_correlation", width=8, height=4, output.type=output.type)
	split.screen(c(1,2))
	screen(1)
	par(mar=mar)
	multi.ecdf(comparisons, xlim=c(-1,1), legend.at=c(-1,1), legend.cex=0.5, xlab='Spearman rank correlation of amino-acid usage')
	screen(2)
	par(mar=mar)
	multi.ecdf(species.rand.comparisons, col=myrainbow(6)[c(1,2)], lty='solid', 
		xlim=c(-1,1), legend.at=c(-1,1), legend.cex=0.5, xlab='Spearman rank correlation of amino-acid usage')
	abline(v=0)
	close.screen(all=TRUE)
	if (file.out) dev.off()

	t.res <- lapply(species.rand.comparisons,t.test, alternative='greater')
	print(t.res)
}