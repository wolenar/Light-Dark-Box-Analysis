# BarPlot for Sal vs EtOH Behavioral Data across RI panel
# #######################################################

# Requires Strain & Treatment columns

# Arguments
############
# pheno: columname that contains pheno data to be plotted
# factor: columnname that contains factor used to group strains
#			 (e.g. factor="Treatment")
# sort: Strain or factor
# x.rot: Logical, rotate x-axis annotations?
# y.lab: y-axis label
# x.lab: x-axis label
# percent: Logical, is y-axis in percentages
# sig.test: logical, compute t-tests and add sig markers?
# by: vector of grouping elements (column names)
# fill: Fill or no fill bars based on this binary factor.
#		Provide column index number or name
# bar.col: vector of colors to differentiate factor levels
#				(e.g., bar.cols=c("blue","red"))
#           (or, aaron.cols)
#ril.barplot(data=all.bxd, pheno="DTL", sort="pvalue")


# Version 1.2
# (DONE) Intelligently determine proper y-axis size
# (DONE) Intelligently identify factors instead of forcing sal/etoh labels
# (DONE) Remove requirement for gplots by manually drawing error bars

# Known issues
# currently can only handle 1 factor with 2 levels

ril.barplot<-function(data, pheno, factor, sort, text.size, x.rot, y.lab,  sig.test, title, bar.cols, fill, legend,...){
	
	###################
	# Variables #######
	###################
	if(missing(pheno)){
		stop(paste("You didn't specify a phenotype. How am I supposed to know",
		"what to plot? I'm a computer, not a mind reader"))
	}
	
	if(missing(factor)){
		stop("You didn't provide a factor. Try again.")
	}
	
	if(missing(sort)){
		sort<-"Strain"
	}
	
	if(missing(text.size)){
		text.size<-1
	}
	
	if(missing(x.rot)){
		x.rot<-FALSE
	}
	
	if(missing(sig.test)){
		sig.test<-FALSE
	}


	# Standard Error Function
	se<-function(x){
		sd(x)/sqrt(length(x)-1)
		}
	
	# Identify factor levels
	factor.levels<-unique(data[,factor])

	# Mean of phenotype for each strain and factor level
	pheno.means<-aggregate(data[,pheno], 
		by=list(Strain=data$Strain,
					data[,factor]), mean)
	colnames(pheno.means)=c("Strain",factor,"Mean")
		
	# Convert long format data frame into "wide" format,
	# where factor levels are split into different columns
	pheno.means<-reshape(pheno.means, timevar=factor, 
		idvar="Strain", direction="wide")

	# SE of phenotype for each strain and treatment	
	pheno.se<-aggregate(data[,pheno], 
		by=list(Strain=data$Strain,
				data[,factor]),se)
	colnames(pheno.se)=c("Strain",factor,"SE")

	pheno.se<-reshape(pheno.se, timevar=factor, 
		idvar="Strain", direction="wide")
		
	# Number of strains at each factor level
	strain.n<-aggregate(data[,pheno], 
		by=list(Strain=data$Strain,
					data[,factor]), length)
	colnames(strain.n)=c("Strain",factor,"n")

	strain.n<-reshape(strain.n, timevar=factor, 
		idvar="Strain", direction="wide")

	# Check that reshaped means and se data.frames have matching rownames
	if(sum(pheno.means$Strain!=pheno.se$Strain)>0) {
		stop("Rownames from means and se data.frames don't match!")
	}
	
	# Reorder means data.frame by user specified criteria
	pheno.means<-pheno.means[order(pheno.means[,sort]),]
	# Reorder se data.frame to mirror means data.frame
	pheno.se<-pheno.se[match(as.character(pheno.means$Strain), 
			as.character(pheno.se$Strain)),]
	strain.n<-strain.n[match(as.character(pheno.means$Strain), 
			as.character(strain.n$Strain)),]

	# Check that rownames match from means and se data.frames
	if(sum(pheno.means$Strain!=pheno.se$Strain)>0) {
		stop("Rownames from means and se data.frames don't match!")
	}

	plot.means<-t(data.frame(pheno.means[,-1], row.names=pheno.means$Strain))
	plot.se<-t(data.frame(pheno.se[,-1], row.names=pheno.se$Strain))
	plot.n<-t(data.frame(strain.n[,-1], row.names=strain.n$Strain))

	# Conduct treatment t-tests for each unique strain
	##################################################
	if(!missing(sig.test) & sig.test==TRUE){
		pvals<-as.numeric()
		for(i in 1:length(unique(data$Strain))){
			pvals[i]<-t.test(data[,pheno]~data[,factor], data=data,
			subset=Strain==unique(data$Strain)[i])$p.value
		}

		t.results<-data.frame(pvalue=pvals, row.names=unique(data$Strain))
		t.results<-t(t.results)
		t.results<-t.results[,match(colnames(plot.means), 
			colnames(t.results))]
		# Identify significant rows at various alpha levels
		#sig01<-t.results<.01
		sig05<-t.results<=.05
	}

	if(missing(bar.cols)){
		bar.cols<-c("red","blue")
	}
	
	# If fill factor is provided create vecetor of colors
	if(!missing(fill)){
		print("FILLLLLL")
		# Identify fill factor level for eachs train
		fill<-data[match(colnames(plot.means),
		data$Strain), fill]
		border.cols<-as.character()
		for(i in 1:length(fill)){
			if(fill[i]==unique(fill)[1]){
				border.cols<-c(border.cols,c(myblue,myred))	
			}
			if(fill[i]==unique(fill)[2]){
				border.cols<-c(border.cols,rep(rgb(0,0,0,0),2))
			}
		}
	} 

	# Calculate upper and lower error bar points
	ci.u<-plot.means+plot.se
	ci.l<-plot.means-plot.se
	
	# Calculate y limits for TLA and percentage behaviors
	if(pheno=="TLA"){
		ymax<-ceiling(max(plot.means)+(3*max(plot.se)))
	} else {
		ymax<-seq(0,100,10)[max(ci.u)<seq(0,100,10)][1]
		ymax<-ifelse((ymax-max(ci.u))<5,ymax+10,ymax)
	}
		
	# Empty plot first
	barplot(plot.means*0, beside=TRUE, axes=FALSE, ylim=c(0,ymax), 
		axisnames=FALSE,
		ylab=ifelse(!missing(y.lab),y.lab,""),
		font.lab=2, cex.lab=text.size, cex.names=text.size) 

		
	if(!missing(title)){
		title(main=title, cex.main=text.size*1.5, font=2, line=2)
		}

	# If ymax is >100 we're obviously not dealing w/ percentages
	# so adjust y-axis accordingly and don't use %
	if(ymax>100){
		axis(2, at=seq(0,ymax,300), labels=FALSE, 
			col="darkgrey",  tck=.98)
		
		axis(2, at=seq(0,ymax,300), labels=seq(0,ymax,300),
			las=2, cex.axis=text.size, tick=TRUE, lwd=text.size)
	}
	else {
	axis(2, at=seq(0,ymax,10), labels=FALSE, col="darkgrey",  tck=.98)
	
	axis(2, at=seq(0,ymax,10), labels=paste(seq(0,ymax,10),
		rep("%",length(seq(0,ymax,10))),sep=""), las=2, 
		cex.axis=text.size, tick=TRUE, lwd=text.size)
	}
	# If fill factor is provided make bar borders thicker
	if(!missing(fill)){
		par(lwd=2)
	}
	
	pheno.plot<-barplot(plot.means, beside=TRUE, axes=FALSE,
		axisnames=FALSE, col=bar.cols, 
		border=ifelse(exists("border.cols"),border.cols,NA),
		ylim=c(0,ymax),
		cex.names=text.size, add=TRUE)

	# Draw error bars
	######################
	
	# Distance between bars			
	bar.dist<-as.numeric(pheno.plot)[1]-as.numeric(pheno.plot)[2]

	segments(x0=as.numeric(pheno.plot)-(bar.dist*.25), y0=as.numeric(ci.u),
				x1=as.numeric(pheno.plot)+(bar.dist*.25), y1=as.numeric(ci.u),
				lwd=text.size)
	
	segments(x0=as.numeric(pheno.plot), y0=as.numeric(ci.l),
				x1=as.numeric(pheno.plot), y1=as.numeric(ci.u),
				lwd=text.size)
				
	segments(x0=as.numeric(pheno.plot)-(bar.dist*.25), y0=as.numeric(ci.l),
				x1=as.numeric(pheno.plot)+(bar.dist*.25), y1=as.numeric(ci.l),
				lwd=text.size)

	# Provides labels to bars if xlab isn't specified
	#axis(1, at=colMeans(pheno.plot), 
	#	labels=ifelse(!missing(x.lab),x.lab,FALSE),
	#	cex.axis=text.size, tick=FALSE,
	#	las=ifelse(!missing(x.rot) & x.rot==TRUE,2,1),
	#	line=text.size/1.5, font=2, lwd=text.size)
	
	# Add bar labels
	mtext(side = 1, at = colMeans(pheno.plot), line=1,
			text=colnames(plot.means), cex=text.size,
			las=ifelse(x.rot==TRUE,2,1))
	
	# Print number of strains per group beneath x-axis
	mtext(side = 1, at = as.numeric(pheno.plot), line = 0,
			text = paste("n=", as.numeric(plot.n),sep=""), 
			font=3, cex=text.size*.7)
	
	# Add significance points if t.tests were conducted
	if(exists("t.results")){
		# Print p-values below bar labels
		mtext(side = 1, at = colMeans(pheno.plot), 
				line = ifelse(x.rot==TRUE,4,2),
				text=formatC(t.results,digits=2,format="e"), 
				font=3, cex=text.size*.7)

		# Size of significance stars is a factor of text size
		star.size<-(text.size*1.2)
		# Pos is determined by number of rows in plot.means
		# because more columns equals smaller bars
		#star.offset<-ncol(plot.means)*.1/(ncol(plot.means)*star.size)
		star.offset<-ncol(plot.means)*.1/(star.size-.5)
		#print(plot.means)
		#print(colMeans(pheno.plot))
		#print(star.size)
		#print(star.offset)
		#if(sum(sig01)>0){
		#	points(x=(colMeans(pheno.plot)+(.5-star.offset))[sig01], 
		#		y=(ci.u[2,]+ymax*.07)[sig01], pch=8,
		#		cex=star.size, lwd=text.size)
		#	points(x=(colMeans(pheno.plot)+(.5+star.offset))[sig01], 
		#		y=(ci.u[2,]+ymax*.07)[sig01], pch=8,
		#		cex=star.size, lwd=text.size)
		#	}
			
		if(sum(sig05)>0){
			points(x=(colMeans(pheno.plot)+.5)[sig05], 
				y=(ci.u[2,]+ymax*.07)[sig05], pch=8,
				cex=star.size, lwd=text.size)
			}
		}
	# Add legend if specified
	if(!missing(legend)){
		if(legend){
			legend("top", legend=sub("Mean.","",rownames(plot.means)), 
			col=bar.cols, pch=15, pt.cex=2, bty="n", horiz=TRUE, cex=text.size)
		}
	}	
}

#out<-ril.barplot(ornl, pheno="DTL", sort="Strain", ymax=70, x.rot=TRUE, sig.test=TRUE, text.size=1.8,fill="Allele.Chr12")