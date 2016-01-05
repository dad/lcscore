setwd("~/research/gelscore/src/")

file.out = T
output.type = 'svg'

fig.sim.size = F
fig.sim.nucleation.comparison = T
fig.sim.conservation = F

parse.num.in.assemblies <- function(s) {
	as.vector(sapply(strsplit(as.character(s), ','), as.numeric))
}

radii.from.counts <- function(x, monomer.volume=1.0, packing.density=0.74) {
	(0.75*x*monomer.volume/(pi*packing.density))^(1/3)
}

plot.assemblies <- function(radii, xlim, ylim, buffer=0, border.buffer=0, ...) {
	y <- random.circles(radii, xlim=xlim, ylim=ylim, buffer=buffer, border.buffer=border.buffer)
	plot.circles(y$circles, ...)
}

if (fig.sim.size) {
	x <- read.table("../data/sim-slow-nucleation.txt", header=T, stringsAsFactors=F)
	temp.jump <- subset(x, temperature != x[1,]$temperature)[1,]
	monomer.numbers <- parse.num.in.assemblies(x[nrow(x),'numbers.in.assemblies'])
	monomer.radius <- 1.0
	monomer.volume <- (4/3)*pi*monomer.radius^3
	radii <- radii.from.counts(monomer.numbers, monomer.volume, packing.density=0.74)

	set.seed(111)

	mar <- c(4,4,1,1)
	if (file.out) dev.out("fig.sim.results", width=12, height=8, output.type=output.type)
	split.screen(c(2,3))
	screen(1)
	par(mar=mar)
	plot(x$time, x$temperature, las=1, type='l', xlab='Time (sec)', ylab='Temperature (\u00B0C)')
	screen(2)
	par(mar=mar)
	ylim <- c(1,max(x$n.monomers))
	plot(x$time, x$n.monomers, ylim=ylim, log='y', type='l', las=1, xlab='Time (sec)', ylab='Number of monomers', yaxt='n')
	my.axis(2, ylim, log=T)
	abline(v=temp.jump$time, lty='dashed')
	screen(3)
	par(mar=mar)
	plot(x$time, x$n.assemblies, log='', type='l', las=1, xlab='Time (sec)', ylab='Number of assemblies')
	abline(v=temp.jump$time, lty='dashed')
	screen(4)
	par(mar=mar)
	plot(x$time, x$largest.radius, log='', type='l', lty='dotted', las=1, xlab='Time (sec)', ylab='Radius (nm)')
	lines(x$time, x$mean.radius, lty='solid')
	lines(x$time, x$intensity.weighted.radius, lty='dashed')
	abline(v=temp.jump$time, lty='dashed')
	screen(5)
	par(mar=mar)
	#multidens(radii, xlab='Radius (nm)', xlim=c(0,max(radii)*1.5), fill=T, col='gray', line.cols='black')
	h <- hist(radii, breaks=seq(0,ceiling(max(radii))), plot=T, col='gray', las=1, main='', ylab='Frequency', xlab='Radius (nm)', yaxs='i')
	#barplot(h$mids, h$count)
	screen(6)
	par(mar=mar)
	buffer <- 1
	border.buffer <- 2
	total.area <- sum(pi*(radii+buffer)^2)
	lim <- c(0,610) #c(0,2*sqrt(total.area)+2*border.buffer)
	plot(1,1,type='n', xlim=lim, ylim=lim, xaxs='i', yaxs='i', xlab='Nanometers', ylab='Nanometers', las=1)
	plot.assemblies(radii, xlim=lim, ylim=lim, col='gray', buffer=buffer, border.buffer=border.buffer, lwd=1.5)
	close.screen(all=TRUE)
	if (file.out) dev.off()
}

# Arrow from (x1,y1) to (x0,y0)
# 
down.arrow <- function(x0, y0=0, x1=x0, y1=NULL, prop, col='black', head.col=NULL, ...) {
	if (is.null(head.col)) {
		head.col <- col
	}
	segments(x0,y0,x1,y1, ...)
	arrow.length <- sqrt(sum((c(x0,y0)-c(x1,y1))^2))
	pu <- par('usr')
	xy.aspect.ratio <- abs((pu[2]-pu[1])/(pu[4]-pu[3]))
	polygon(c(x1,x1-prop*0.7*arrow.length*xy.aspect.ratio, x1+prop*0.7*arrow.length*xy.aspect.ratio),c(y1, y1+prop*arrow.length, y1+prop*arrow.length), col=head.col, ...)
}

if (fig.sim.nucleation.comparison) {
	rates <- c('1e-12','1e-14','1e-15')
	fnames <- paste("../data/sim-nucleation-",rates,".txt", sep='')
	runs <- lapply(fnames, function(f) { read.table(f, header=T, sep='\t')})

	set.seed(111)
	n.runs <- length(rates)
	t.max <- max(sapply(runs, function(x){max(subset(x, time<10000)$time)}))
	temp.jump.time <- 1000
	timelim <- c(0,t.max)
	spacelim <- c(0,1000)
	buffer <- 1
	border.buffer <- 2
	mar <- c(4,4,1,1)
	if (file.out) dev.out("fig.nucleation.results", width=12, height=8, output.type=output.type)
	split.screen(c(2,n.runs))
	for (xi in 1:n.runs){
		x <- runs[[xi]]
		monomer.numbers <- parse.num.in.assemblies(x[nrow(x),'numbers.in.assemblies'])
		monomer.radius <- 1.0
		monomer.volume <- (4/3)*pi*monomer.radius^3
		radii <- radii.from.counts(monomer.numbers, monomer.volume, packing.density=0.74)
		screen(xi)
		par(mar=mar)
		plot(x$time, x$prop.assembled, xlim=timelim, ylim=c(0,1), log='', lwd=2, type='l', las=1, xlab='Time (sec)', ylab='Proportion of monomers in assemblies')
		down.arrow(temp.jump.time, 0.25, temp.jump.time, 0, prop=0.15, lwd=2, col='red')
		time.95 <- min(subset(x, prop.assembled>0.95)$time)
		time.05 <- min(subset(x, prop.assembled>0.05)$time)
		abline(v=time.95, lty='dotted', lwd=1.5)
		abline(v=time.05, lty='dashed', lwd=1.5)
		text(5000, 0.5, paste("T 5% = ", format(time.05-temp.jump.time,digits=2), 's', sep=''))
		text(5000, 0.4, paste("T 95% = ", format(time.95-temp.jump.time,digits=2), 's', sep=''))

		screen(xi+n.runs)
		par(mar=mar)
		plot(1,1,type='n', xlim=spacelim, ylim=spacelim, xaxs='i', yaxs='i', xlab='Nanometers', ylab='Nanometers', las=1)
		plot.assemblies(radii, xlim=spacelim, ylim=spacelim, col='gray', buffer=buffer, border.buffer=border.buffer, lwd=1.5)
	}

	close.screen(all=TRUE)
	if (file.out) dev.off()
}

if (fig.sim.conservation) {
	x <- read.table("../data/sim-test.txt", header=T, stringsAsFactors=F)
	monomer.numbers <- parse.num.in.assemblies(x[nrow(x),'numbers.in.assemblies'])
	total.number <- function(xi) {
		x[xi,'n.monomers'] + sum(parse.num.in.assemblies(x[xi,'numbers.in.assemblies']))
	}
	plot(x$time, sapply(1:nrow(x),total.number))
}
