rm(list = ls());

library(dplyr);
library(plyr);
library(tidyr);
library(reshape);
library(ggh4x);
library(patchwork);
library(Polychrome);

#parameters
CORRECTED_CSV_FILE <- "SIL.csv_CorrectedJamie.csv"; #Exact IsoCorrectoR _Corrected.csv file
COMPOUND_FILE <- "compoundList.txt"; #Newline-delimited file of the compounds to plot (must match the IsoCorrectoR file)
METADATA_FILE <- "metadataJamie.txt"; #Tab-delimited file with:
	#first column = Name (exact match to the column name of the Compound File). Important: starts with a letter, no spaces, hyphens, plus signs
	#second column = Timepoint #any format, using a consistent unit (don't mix h and d)
	#third column = CellLine
	#fourth column = Treatment #could be BSA, Stearate, etc.
	#fifth column = RealSample #TRUE or FALSE -- exclude blanks, QCpools, etc.
	#the first row contains these headers names exactly as written
OUTFILE_NAME <- "barplots.pdf";
WIDTH <- 18;
HEIGHT <- 6;

#read in metadata on samples, compound list
groupdata <- read.delim(METADATA_FILE,stringsAsFactors=FALSE);
compoundsToPlot <- readLines(COMPOUND_FILE,warn=FALSE);

#read in IsoCorrectoR csv and rearrange for plotting
mydata <- read.csv(CORRECTED_CSV_FILE);
mydataFiltered <- mydata %>% plyr::rename(c("X"="CompoundName")) %>% separate(CompoundName, c("CompoundName","Isotopologue"), "_");
longData <- mydataFiltered %>% reshape::melt(., i=c("CompoundName", "Isotopologue")) %>%
plyr::rename( c("variable"="Name")) %>%
select(CompoundName, Name, Isotopologue, value);
longData$Isotopologue <- as.numeric(longData$Isotopologue);

#add additional info the the isotope dataframe
dataLabeledLong <- full_join(longData, groupdata, by="Name");

#filter out RealSample=FALSE samples
tracedDataLong <- dataLabeledLong[dataLabeledLong$RealSample,];

#get center and dispersion for plotting; can change median to mean and mad (median absolute deviation) to sd
summarizedDataLong <- tracedDataLong %>% ddply(., c("CompoundName","CellLine","Timepoint","Treatment","Isotopologue"), summarize, my_avg=median(value), my_dispersion=mad(value));

#make timepoint a factor with orders
myTimepoints <- groupdata[groupdata$RealSample,"Timepoint"];
orderedTimepoints <- unique(myTimepoints[order(as.numeric(gsub("\\D", "", myTimepoints)))]);
summarizedDataLong$Timepoint <- factor(summarizedDataLong$Timepoint, orderedTimepoints);

#make the colors then organize the colors
maxIsotopologues <- max(summarizedDataLong$Isotopologue)+1; #+1 for the M+0 isotopologue
IsotopologuePalette <- as.vector(createPalette(maxIsotopologues,c("#ff0000","#00ff00","#0000ff")));
pie(rep(1,length(IsotopologuePalette)),col=IsotopologuePalette); #check colors chosen

#Helper function
#Calculate error bar center as cumsum. This assumes data are by increasing heaviness of the Isotopologue, which is what IsoCorrectoR does
cumulative_center <- function(this_data) {
	this_pos <- 0;
	this_data$my_pos <- rep(0,nrow(this_data));
	for(j in 1:nrow(this_data)) {
		if(this_data[j,"Isotopologue"]==0) {
			this_pos <- this_data[j,"my_avg"];
		} else {
			this_pos <- this_pos + this_data[j,"my_avg"];
		}
		this_data[j,"my_pos"] <- this_pos;
	}
	return(this_data);
}

#Helper function
#Calculate the fraction of total for each abundance metric. This assumes data are by increasing heaviness of the Isotopologue, which is what IsoCorrectoR does
fraction_center <- function(this_data) {
	counter <- 0;
	mysum <- 0;
	this_data$my_total <- rep(0,nrow(this_data));
	this_data$my_fraction_avg <- rep(0,nrow(this_data));
	this_data$my_fraction_dispersion <- rep(0,nrow(this_data));
	this_data$my_fraction_pos <- rep(0,nrow(this_data));
	for(j in 1:nrow(this_data)) {
		if(j==nrow(this_data) | this_data[j+1,"Isotopologue"] < this_data[j,"Isotopologue"]) { #this is the last row of a set
			mysum <- mysum + this_data[j,"my_avg"];
			for(k in seq(j-counter,j)) {
				this_data[k,"my_total"] <- mysum;
				this_data[k,"my_fraction_avg"] <- this_data[k,"my_avg"]/mysum;
				this_data[k,"my_fraction_dispersion"] <- this_data[k,"my_dispersion"]/mysum;
				this_data[k,"my_fraction_pos"] <- this_data[k,"my_pos"]/mysum;
			}
			#reset
			counter <- 0;
			mysum <- 0;
		} else {
			counter <- counter+1;
			mysum <- mysum + this_data[j,"my_avg"];
		}
	}
	return(this_data);
}

#Plotting function -- defaults should work with the wrangling script above, just set the plottype="absolute" or plottype="percetnage"
#Makes barplots of absolute peak areas, loop through compounds
barplot_function <- function(plottype,mydata=summarizedDataLong,outname=OUTFILE_NAME,compounds=compoundsToPlot,myPalette=IsotopologuePalette) {
	pdf(paste(plottype,outname,sep=""),width=WIDTH,height=HEIGHT,onefile=TRUE);
	for(i in compounds) {
		plot_data <- subset(mydata, CompoundName == i);
		#Brute force calculate absolute error bar center
		plot_data <- cumulative_center(plot_data);
		#Brute force calculate percentage bars and error bars
		plot_data <- fraction_center(plot_data);
		plot_data$Isotopologue <- factor(plot_data$Isotopologue); #factor rather than numeric required for fill
		#plot
		if(plottype!="absolute" & plottype!="percentage") {
			stop("Wrong parameter for plottype in barplot_function!");
		}
		barplot <- 
			(if(plottype=="absolute") { #main
				ggplot(plot_data, aes (x=Timepoint, y=my_avg, fill=Isotopologue))
			} else { #percentage
				ggplot(plot_data, aes (x=Timepoint, y=my_fraction_avg, fill=Isotopologue))
			}
			)+
			geom_bar(position = position_stack(reverse=TRUE), stat="identity", color="black", linewidth=0.5)+ #bars
			(if(plottype=="absolute") { #errorbars
				geom_errorbar(aes(ymin= my_pos-my_dispersion/2,ymax=my_pos+my_dispersion/2),color="black",linewidth=0.5,width=0.2)
			} else { #percentage
				geom_errorbar(aes(ymin= my_fraction_pos-my_fraction_dispersion/2,ymax=my_fraction_pos+my_fraction_dispersion/2),color="black",linewidth=0.5,width=0.2)
			})+
			guides(fill=guide_legend(title=NULL))+
			scale_fill_manual(values = myPalette)+
			(if(plottype=="absolute") {
				labs(x=NULL, y="Corrected abundance")
			} else { #percentage
				labs(x=NULL, y="Corrected percentage")
			})+
			theme_classic()+
			theme(text = element_text(size=14),
			axis.text.x = element_text(angle=45, hjust=1))+
			facet_nested_wrap(vars(Treatment, CellLine), nrow=1, nest_line=element_line(linetype=1), scales="free_y", axes="all")+
			theme(strip.background=element_blank(), strip.text.x=element_text(size=14,color="black"))+
			(if(plottype=="absolute") {
				scale_y_continuous(labels = scales::scientific)
			} else { #percentage
				scale_y_continuous(labels = scales::percent_format())
			})
		Titleplot <- barplot + plot_annotation(title = paste(i, ""),
			theme = theme(plot.caption = element_text(size=12), plot.title=element_text(size=14,face="bold",hjust=0.5)))
		print(Titleplot)
	}
	dev.off(); #can now open the pdfs in Illustrator to copy-and-paste as EMF/SVG into Powerpoint
}

barplot_function(plottype="absolute");
barplot_function(plottype="percentage");
