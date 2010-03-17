ldbox.analysis<-function(data, plot, export.raw, sig.test=TRUE, ...) {
	
	# Capture input file name
	input.name<-strsplit(data,split="/")
	input.name<-input.name[[1]][length(input.name[[1]])]
	
	# Read in data
	data<-read.csv(file=data, stringsAsFactors=FALSE)
	
	# Check for blank cells or lines
	if(sum(is.na(data))>0){
		stop("So, there appears to be either a blank line or a blank cell in your data. Please delete the offending cell or line and try again.", call.=FALSE)
	}
	
	# If 'Include' column is provided exclude
	# mice with a FALSE value and get rid of 
	# the column
	if("Include" %in% colnames(data)){
		data<-data[data$Include==TRUE,-grep("Include",colnames(data))]
	}
	
	if(missing(export.raw)){
		export.raw<-FALSE
	}
	
	# Convert to class Date (doesn't work)
	#data$Date.run<-as.Date(as.character(data$Date.run), "%m%d%y")
	
	# Create folder to store results
	folder.name<-paste("BXD_LDbox_Analysis_", Sys.Date(), sep="")
	system(paste("mkdir ", folder.name, sep=""))

	# Change working directory to folder
	wd<-getwd()
	wd<-paste(wd,"/",folder.name,"/",sep="")

	# First identify columns containing time data
	time.cols<-grep("time", colnames(data), ignore.case=TRUE)
	# Identify columns with phenotypes
	phenotypes<-c("Weight", "Dist.trvl", "Time.amb", "Amb.counts", "Time.ster",
		"Ster.counts", "Time.rest", "Vert.counts", "Vert.time","Zone.entries",
		"Zone.time")
	pheno.cols<-match(phenotypes, colnames(data))
	
	# Convert MedAssociates time format to seconds for 
	# Time.amb, Time.ster, Time.rest, Vert.time and Zone.time
	##########################################################
	for(i in time.cols){
		raw.times<-data[,i]
		# Identify rows that need to be converted
		unconverted<-grep(":", raw.times)
		# Skip to next column if no rows need converting
		if(length(unconverted>0)){
			data[unconverted,i]<-sapply(strsplit(raw.times[unconverted],split=":"), 
				FUN=function(x) as.numeric(tail(x,n=3)[1])*3600 + 
				as.numeric(tail(x,n=3)[2])*60 + 
				as.numeric(tail(x,n=3)[3]))
			data[,i]<-as.numeric(data[,i])
			} else{
				next()
			}
		}

	# Export raw data with processed time columns
	if(export.raw==TRUE){	
		write.csv(data,
			file=paste(wd,sub(".csv","_convertedTime.csv",x=input.name),sep=""),
				row.names=FALSE)
	}

	# Identify all timepoints
	timepoint<-levels(factor(data$Time.point))

	# Data sanity check
	######################
	# Make sure zone totals add up to time points by calculating
	# variance within a timepoint. Any deviation from zero will
	# cause an error. 
	
	zone.totals<-aggregate(as.numeric(data$Zone.time), 
		by=list(Mouse=data$Mouse, Timepoint=data$Time.point), FUN=sum)
		
	for(i in 1:length(timepoint)){
		if(var(subset(zone.totals, 
				subset=zone.totals$Timepoint==timepoint[i], select="x"))>0) {
			stop("Zone totals differ within time points. Check your data!", 
				call.=FALSE)
			} else {
				writeLines(paste(timepoint[i],"data passes inspection..."))
			}
		}
	
	# Calculate individual mouse behaviors
	######################################
	# Iterate through all time points
	for(i in 1:length(timepoint)) {
	
		# Subset for for 5min or 10min timepoint
		data.timepoint<-subset(data, Time.point==timepoint[i])

		# Get light zone values for 'percentage behaviors'
		light<-subset(data.timepoint, subset=data.timepoint$Zone=="light",
				select=c("Dist.trvl","Time.rest","Zone.time"))

		# Get sum of 'percentage behaviors' across all zones
		all.zones<-aggregate(data.timepoint[,
		c("Dist.trvl","Time.rest","Zone.time")],
		by=list(data.timepoint$Mouse), function(x) sum(as.numeric(x)))
	
		# Remove first column created by aggregate function
		all.zones<-all.zones[,-1]
	
		# Calculate percentages
		behaviors<-(light/all.zones)*100

		# Add total locomotor activity to behaviors data frame
		behaviors<-cbind(behaviors, all.zones$Dist.trvl)

		# Change behavior names
		phenotype.names<-c("PDT", "PTR", "PTS", "TLA")
		colnames(behaviors)=phenotype.names

		# Add Cage, Strain, Weight and Treatment to data
		# by adding all columns preceding zone column
		behaviors<-cbind(subset(data.timepoint[,2:(grep("Zone"
			,colnames(data.timepoint))[1]-1)], 
			data.timepoint$Zone=="light"), behaviors)
	
		# Set mouse ID's as rownames
		row.names(behaviors)<-subset(data.timepoint[,"Mouse"],
			 data.timepoint$Zone=="light")
		
		# Store in a separate dataframe that will hold results for both 
		# timepoints, first add timepoint so data can be teased out later
		behaviors<-cbind(Time.point=rep(timepoint[i], nrow(behaviors)),behaviors)
		ifelse(exists("indv.behaviors"),
			indv.behaviors<-rbind(indv.behaviors,behaviors),
			indv.behaviors<-behaviors)	
		}
	
	# Export individual behavior results
	write.table(indv.behaviors,
		file=paste(wd,"Individual_Mouse_LDbox_Behaviors_",
		Sys.Date(),".csv",sep=""),row.names=FALSE,sep=",")
	
	# Calculate strain behaviors
	############################
	# Means
	strain.means<-aggregate(indv.behaviors[,phenotype.names], 
				by=list(Strain=indv.behaviors$Strain, 
				Treatment=indv.behaviors$Treatment,
				Time.point=indv.behaviors$Time.point),mean)
	
	# Standard Error Function
	se<-function(x){
		sd(x)/sqrt(length(x)-1)
		}
	
	# Standard Error
	strain.se<-aggregate(indv.behaviors[,phenotype.names], 
				by=list(Strain=indv.behaviors$Strain, 
				Treatment=indv.behaviors$Treatment,
				Time.point=indv.behaviors$Time.point),se)
				
	# Export BXD strain behaviors
	write.table(strain.means, file=paste(wd,"Mean_BXD_LDbox_Behaviors_",
				Sys.Date(),".csv",sep=""),row.names=FALSE,sep=",")

	# Export BXD strain behaviors
	write.table(strain.se, file=paste(wd,"SE_BXD_LDbox_Behaviors_",
				Sys.Date(),".csv",sep=""),row.names=FALSE,sep=",")
				
	# Plots
	########
	if(plot==TRUE) {
		writeLines("Generating plots just for you...")
		# Output plots for each behavioral measurment at each timepoint
		# Iterate through timepoints
		for(i in 1:length(timepoint)){
			plot.data<-indv.behaviors[indv.behaviors$Time.point==timepoint[i],]

			pdf(file=paste(wd,timepoint[i],"_Plots_BXD_LDbox_Behaviors_", 
			Sys.Date(),".pdf",sep=""), paper="USr")
			
			# Iterate through all analyzed phenotypes
			for(p in 1:length(phenotype.names)) {

				ril.barplot(data=plot.data, pheno=phenotype.names[p],
						factor="Treatment", sig.test=sig.test, 
						legend=TRUE, text.size=1.3,
						title=paste(timepoint[i],": ",phenotype.names[p]),...)
				}
			dev.off()
			}
		}
	return(indv.behaviors)
	}
