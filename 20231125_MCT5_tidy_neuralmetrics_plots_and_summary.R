library(tidyverse)


##Instructions for dealing with NeuralMetrics Update:
#Delete three lines above line 26 in the neural metrics file
# make the "well averages" row line 27
# note also that you will get an error if the summary sheet output is left in the neuralmetrics folder

#Set wd, use ctrl shift H for shortcut
setwd("C:\\Users\\nmloj\\Desktop\\Desktop_input")
#makes a list of the file names in the working directory
###begin neural metrics tidying

#change to insert plate name!!
plate_num<-"MCT5_100-9636"


files = list.files()
#lapply seems to apply read.csv across the list above
#creates a large list from files in the wd
data2=lapply(files,read.csv, sep= ",",skip = 29,nrows = 57, check.names=FALSE, fill=TRUE, header= TRUE)
#ERROR in above line?: be sure that you are skipping the correct number of columns! "Well Averages" should be first row.
#for loop column bind files vector as a column for identification of rows
#also transpose each csv so that rows become columns
for (i in 1:length(data2)){data2[[i]]<-cbind(t(data2[[i]]),files[i])}
#runs rbind to concatenate CSVs across the list, returns one big csv
data2_f <- do.call("rbind", data2) 

#data2_f <-rownames_to_column(data2_f, "well_names")
##addition of header row
#extract column names
col_names_ve<-data.frame(data2_f[1,]) #extract first row from dataset

my_header<-as.character(col_names_ve[,1])#this is our header vector


#creates one large dataframe
data2<-lapply(X= data2, FUN=as.data.frame)
data2<-lapply(X= data2, FUN=rownames_to_column, "well_name")

data_headl<-lapply( X= data2, FUN=head,1)
data2<-lapply( X= data2, FUN=tail,-1)#removes header row for all csv in list, also need to remove header from every dataframe

data_3.1 <- do.call("rbind.data.frame", data2)# list into 
#add header to large csv
my_header<-append(my_header,"well_name",0)# need to add missing well name label
colnames(data_3.1)<-my_header
colnames(data_3.1)[59] <- "file_name" #replace column header for file names


#select only the columns that you want. You def want column 112 for naming, others are optional
#currently: 32 = Net burst freq, 51= AUNCC, 24/25 = Burst freq avg/ std
data_3.2<-data_3.1%>%
  select(1,2,59,4:12,24, 25, 32, 51)
#Enter manually number of wells wanted
wells_want<-c("A5", "B5","C5", "D5","E5","F5",
              "A6","B6","C6", "D6","E6","F6",
              "A7","B7","C7", "D7","E7","F7",
              "A8","B8","C8", "D8","E8","F8")
#gets the total number of wells
num_wells<-length(wells_want)
#
data_3.2<-data_3.2 %>%
  separate(file_name,c("plate", "junk","DIV", "junk1"), "_" )%>%#extract plate name
  mutate(DIV= str_sub(DIV,1,2))%>%
  mutate(well_name=as.factor(well_name))%>%
  mutate(`    Mean Firing Rate (Hz)`=as.numeric(`    Mean Firing Rate (Hz)`))%>%
  mutate(`    Weighted Mean Firing Rate (Hz)`=as.numeric(`    Weighted Mean Firing Rate (Hz)`))%>%
  mutate(`    Number of Active Electrodes`=as.numeric(`    Number of Active Electrodes`))%>%
  mutate(`    Network Burst Frequency`=as.numeric(`    Network Burst Frequency`))%>%
  filter(well_name!="")%>%
  filter(well_name %in% wells_want)

#summarise data, group by the DIV, then take mean and/or sd
MEA_summary_sheet<-data_3.2%>% 
  group_by(DIV) %>% 
  summarise( MFR=mean(`    Mean Firing Rate (Hz)`, na.rm = TRUE),
             MFR_sd = sd(`    Mean Firing Rate (Hz)`, na.rm = TRUE),
             WMFR= mean(`    Weighted Mean Firing Rate (Hz)`, na.rm = TRUE),
             WMFR_sd= sd(`    Weighted Mean Firing Rate (Hz)`, na.rm = TRUE),
             AEY = ((sum(`    Number of Active Electrodes`, na.rm = TRUE)/(n=n()*16))*100),
             Num_active_electrodes = ((sum(`    Number of Active Electrodes`, na.rm = TRUE))),
            Net_bust_freq=mean(`    Network Burst Frequency`, na.rm = TRUE), 
            Net_bust_freq_sd=sd(`    Network Burst Frequency`, na.rm = TRUE) )%>%
  mutate(DIV =as.numeric(DIV))%>%
  mutate(total_electrodes=(num_wells*16))

#for saving summary sheet
wd<-getwd()
filename_csv<-paste(plate_num,"_summary",".csv", sep = "") #concatenate file name + "_summary" ".csv"
filename_csv<-file.path(wd,filename_csv) #make file path
write.csv(MEA_summary_sheet,filename_csv)#save to WD as .csv


###Make plots
#x axis limit
x_day<-27
#AEY plot
AEY_plot<- ggplot(MEA_summary_sheet,aes(x=DIV, y=AEY))+
  geom_point(aes(),color='darkblue', show.legend = FALSE)+
  geom_line()+
  ylim(0,100)+
  xlim(0,x_day)+
  ggtitle("AEY", subtitle = plate_num)+
  theme_classic()

AEY_plot

#MFR plot
MFR_plot<- ggplot(MEA_summary_sheet,aes(x=DIV, y=MFR))+
  geom_point(aes(),color='red', show.legend = FALSE)+
  geom_line(color= "red")+
  geom_errorbar(aes(ymin=MFR -MFR_sd, ymax=MFR +MFR_sd), width=.2, position=position_dodge(0.2), color= "red")+
  ylim(0,17)+
  xlim(0,x_day)+
  ggtitle("MFR", subtitle = plate_num)+
  theme_classic()

MFR_plot


#WMFR plot
WMFR_plot<- ggplot(MEA_summary_sheet,aes(x=DIV, y=WMFR))+
  geom_point(aes(),color='orange', show.legend = FALSE)+
  geom_line(color= "orange")+
  geom_errorbar(aes(ymin=WMFR -WMFR_sd, ymax=WMFR +WMFR_sd), width=.2, position=position_dodge(0.2), color= "orange")+
  ylim(0,18)+
  xlim(0,x_day)+
  ggtitle("WMFR", subtitle = plate_num)+
  theme_classic()

WMFR_plot

NBF_plot<- ggplot(MEA_summary_sheet,aes(x=DIV, y=Net_bust_freq))+
  geom_point(aes(),color='olivedrab', show.legend = FALSE)+
  geom_line(color= "olivedrab")+
  geom_errorbar(aes(ymin=Net_bust_freq -Net_bust_freq_sd, ymax=Net_bust_freq +Net_bust_freq_sd), width=.2, position=position_dodge(0.2), color= "olivedrab")+
  ylim(-.2,1)+
  xlim(0,x_day)+
  ggtitle("Network Burst Frequency", subtitle = plate_num)+
  theme_classic()

NBF_plot
###########################################################################
#Statistical testing for stability on active electrodes
last_3_recordings_AEY<- MEA_summary_sheet%>%
  mutate(the_rank  = rank(-DIV, ties.method = "random"))%>%
  filter(the_rank <= 3) 


#two sample proportion test
prop.test(x = last_3_recordings_AEY$Num_active_electrodes, n = last_3_recordings_AEY$total_electrodes)
prop.test(x = c(316,325), n = c(384,384)) #enter manually from summary sheet
prop.test(x = c(322,325), n = c(384,384)) #enter manually from summary sheet


library(rstatix)

#means comparison
last_3_recordings_MFR<- data_3.2%>%
  #mutate(DIV=as.numeric(DIV))%>%
  slice_max(DIV,n=(num_wells*3), with_ties = TRUE)%>%
  mutate(`    Number of Spikes`)


pwc <- last_3_recordings_MFR %>% 
  dunn_test(`    Weighted Mean Firing Rate (Hz)` ~ DIV, p.adjust.method = "bonferroni") 
pwc
