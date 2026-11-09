library(tidyverse)


##Instructions for dealing with NeuralMetrics Update:
#Delete three lines above line 26 in the neural metrics file
# make the "well averages" row line 27
# note also that you will get an error if the summary sheet output is left in the neuralmetrics folder

#Set wd, use ctrl shift H for shortcut
#setwd("C:\\Users\\nmloj\\Desktop\\Desktop_input")
#setwd("C:/Users/Neal/OneDrive - UMass Lowell/Pain Research Lab/Grants/NeuroRad_GMHA (2024)/Data/MCT8/neuralmetrics")
setwd("C:/Users/Neal/OneDrive - UMass Lowell/Pain Research Lab/Grants/NeuroRad_GMHA (2024)/Data/MCT8/neuralmatrics_adjust_for_burst/")
#makes a list of the file names in the working directory
###begin neural metrics tidying

#change to insert plate name!!
plate_num<-"MCT8_108-3520"


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
wells_want<-c("A1", "B1","C1", "D1","E1","F1",
              "A2","B2","C2", "D2","E2","F2",
              "A3","B3","C3", "D3","E3","F3",
              "A4","B4","C4", "D4","E4","F4")
#gets the total number of wells
num_wells<-length(wells_want)
#
#####label wells condition
ami_wells<- c("A1", "B1","C1", "D1","E1","F1")

ami_wells<-c("A1"= "amifostine","A3" = "amifostine",
             "B1" = "amifostine","B3" = "amifostine",
             "C1" = "amifostine", "C3" = "amifostine",
             "D1" = "amifostine", "D3" = "amifostine",
             "E1" = "amifostine","E3" = "amifostine",
             "F1" = "amifostine","F3" = "amifostine")

vehicle_wells<-c("A2" = "vehicle","A4"= "vehicle",
                 "B2"= "vehicle", "B4"= "vehicle",
                 "C2"= "vehicle","C4"= "vehicle",
                 "C2"= "vehicle","C4"= "vehicle",
                 "D2"= "vehicle","D4"= "vehicle",
                 "E2"= "vehicle","E4"= "vehicle",
                 "F2"= "vehicle","F4"= "vehicle" )
#####Radiation status
rad_wells<-c("A1"= "10_Gy_eq","A2" = "10_Gy_eq",
             "B1" = "10_Gy_eq","B2" = "10_Gy_eq",
             "C1" = "10_Gy_eq", "C2" = "10_Gy_eq",
             "D1" = "10_Gy_eq", "D2" = "10_Gy_eq",
             "E1" = "10_Gy_eq","E2" = "10_Gy_eq",
             "F1" = "10_Gy_eq","F2" = "10_Gy_eq")

sham_wells<-c("A3" = "sham","A4"= "sham",
                 "B3"= "sham", "B4"= "sham",
                 "C3"= "sham","C4"= "sham",
                 "C3"= "sham","C4"= "sham",
                 "D3"= "sham","D4"= "sham",
                  "E3"= "sham","E4"= "sham",
                 "F3"= "sham","F4"= "sham" )
#

#In this version we only want to filter for the relevant time points
#we also must add conditional Identifiers for each timepoint
#must adjust filter for desired DIV 
timepoints<-data_3.2 %>%
  separate(file_name,c("plate", "junk","DIV", "junk1"), "_", remove = FALSE )%>%#extract plate name
  mutate(DIV= str_sub(DIV,1,2))%>%
  separate(file_name,c("plate", "junk","DIV2", "junk1"), "_", remove = TRUE )%>%#extract plate name
  mutate(time= str_sub(DIV2,4,4))%>%
  mutate(well_name=as.factor(well_name))%>%
  filter(DIV>=27)%>%

  mutate(`    Mean Firing Rate (Hz)`=as.numeric(`    Mean Firing Rate (Hz)`))%>%
  mutate(`    Weighted Mean Firing Rate (Hz)`=as.numeric(`    Weighted Mean Firing Rate (Hz)`))%>%
  mutate(`    Number of Active Electrodes`=as.numeric(`    Number of Active Electrodes`))%>%
  mutate(`    Network Burst Frequency`=as.numeric(`    Network Burst Frequency`))%>%
  mutate(`    Area Under Normalized Cross-Correlation`=as.numeric(`    Area Under Normalized Cross-Correlation`))%>%
  mutate(`    ISI Coefficient of Variation - Avg`=as.numeric(`    ISI Coefficient of Variation - Avg`))%>%
  mutate(`    Number of Bursts`=as.numeric(`    Number of Bursts`))%>%
  filter(well_name!="")%>%
  filter(`    Number of Spikes`!="")%>% #eliminate non-recored wells in timepoint 2 and beyond
  filter(well_name %in% wells_want) #filter for only the well vector specified above
#change names of time variable
timepoints<-timepoints%>%
  mutate(RP_status=str_replace_all(well_name,ami_wells,"amifostine"))%>%
  mutate(RP_status=str_replace_all(RP_status,vehicle_wells,"vehicle"))%>% #note that input column name must be changed here!
  filter(RP_status == c("amifostine" , "vehicle"))%>%
  
  mutate(rad_stat=str_replace_all(well_name,rad_wells,"10_Gy_eq"))%>%
  mutate(rad_stat=str_replace_all(rad_stat,sham_wells,"Sham"))%>%
  #filter(RP_status == c("10_Gy_eq" , "Sham"))%>%
  
  mutate(time=str_replace_all(time,"0","Pre-RP treatment"))%>%
  mutate(time=str_replace_all(time,"1$","1 hr")) %>% 
  mutate(time=str_replace_all(time,"2$","1 day"))%>% #remove for <40min timepoint
  mutate(time=str_replace_all(time,"3$","2 day"))%>%
  mutate(time=str_replace_all(time,"4$","3 day"))%>%
  mutate(time=str_replace_all(time,"5$","4 day"))%>%
  mutate(time=str_replace_all(time,"6$","5 day"))%>%
  mutate(time=str_replace_all(time,"7$","6 day"))%>%
  mutate(time=str_replace_all(time,"8$","7 day"))%>%
  
  
  
  unite("all_stat", rad_stat,RP_status, sep = "_", remove = FALSE)  %>% #create the all_stat column which combines RP and rad
  mutate(all_stat= as.factor(all_stat))%>%
  
  mutate(across(`    Network Burst Frequency`, ~ replace(., is.na(.), 0))) #change wells with na NBF to 0


  #########
#filter for unwanted wells
timepoints<-timepoints%>%
  filter(well_name != c("D4"))%>% #this well did not meet the minimum electrode criteria at pre-treatment timepoint
  filter(well_name != c("D3")) #right  now I am eliminating this well because it had major outlier values
  



######################################################
#create DoS 
#DoS = ((post-pre)/(post+pre))
all_Dos<-timepoints%>%
  group_by(well_name)%>%
  mutate(AE_DoS=((`    Number of Active Electrodes`-`    Number of Active Electrodes`[time=="Pre-RP treatment"])/(`    Number of Active Electrodes`+`    Number of Active Electrodes`[time=="Pre-RP treatment"])))%>%
  mutate(WMFR_DoS=((`    Weighted Mean Firing Rate (Hz)`-`    Weighted Mean Firing Rate (Hz)`[time=="Pre-RP treatment"])/(`    Weighted Mean Firing Rate (Hz)`+`    Weighted Mean Firing Rate (Hz)`[time=="Pre-RP treatment"])))%>%
  mutate(MFR_DoS=((`    Mean Firing Rate (Hz)`-`    Mean Firing Rate (Hz)`[time=="Pre-RP treatment"])/(`    Mean Firing Rate (Hz)`+`    Mean Firing Rate (Hz)`[time=="Pre-RP treatment"])))%>%
  mutate(AUNCC_DoS=((`    Area Under Normalized Cross-Correlation`-`    Area Under Normalized Cross-Correlation`[time=="Pre-RP treatment"])/(`    Area Under Normalized Cross-Correlation`+`    Area Under Normalized Cross-Correlation`[time=="Pre-RP treatment"])))%>%
  mutate(ISI_CoV_DoS=((`    ISI Coefficient of Variation - Avg`-`    ISI Coefficient of Variation - Avg`[time=="Pre-RP treatment"])/(`    ISI Coefficient of Variation - Avg`+`    ISI Coefficient of Variation - Avg`[time=="Pre-RP treatment"])))%>%
  mutate(NBF_DoS=((`    Network Burst Frequency`-`    Network Burst Frequency`[time=="Pre-RP treatment"])/(`    Network Burst Frequency`+`    Network Burst Frequency`[time=="Pre-RP treatment"])))%>%
  mutate(NumB_DoS=((`    Number of Bursts`-`    Number of Bursts`[time=="Pre-RP treatment"])/(`    Number of Bursts`+`    Number of Bursts`[time=="Pre-RP treatment"])))%>%
  

  
  filter(time != "Pre-RP treatment")%>% #no longer need the pre treamtent rows, remove here
  mutate(across(time, factor, levels=c("1 hr", "1 day", "2 day", "3 day", "4 day", "5 day", "6 day", "7 day")))

####################################################
#normality testing
library(rstatix)
library(FSA)

#shapiro_test() doesnt seem to like the long column names with spaces, make new df
stats_rename_all<-all_Dos%>%
  rename(AUNCC = `    Area Under Normalized Cross-Correlation`,
         MFR= `    Mean Firing Rate (Hz)`, 
         WMFR= `    Weighted Mean Firing Rate (Hz)`,
         AE= `    Number of Active Electrodes`,
         ISI_CoV= `    ISI Coefficient of Variation - Avg`,
         NBF= `    Network Burst Frequency`)


all_test_norm<-stats_rename_all%>%
  group_by(time, RP_status, rad_stat) %>%
  shapiro_test(AUNCC, 
               MFR, 
               WMFR,
               AE, 
               ISI_CoV,)
all_test_norm<-all_test_norm%>%
  arrange(variable)
#norm testing for DoS
all_test_norm_DoS<-stats_rename_all%>%
  group_by(time, RP_status, rad_stat) %>%
  shapiro_test(AUNCC_DoS, 
               MFR_DoS, 
               WMFR_DoS,
               AE_DoS, 
               ISI_CoV_DoS,)
all_test_norm_DoS<-all_test_norm%>%
  arrange(variable)
#RESULT: same stats as the Non- DoS


# R doesnt like NBF because of NAs for shapiro_test
#######################################
#Stats
stats_rename_all$id<-rep(1:46)

#normal dist
aov.MFR<-stats_rename_all %>%
  group_by(time ) %>%
  anova_test(MFR_DoS ~ all_stat )
get_anova_table(aov.MFR)
#NS
aov.AUNCC<-stats_rename_all %>%
  group_by(time ) %>%
  anova_test(AUNCC_DoS ~ all_stat )
get_anova_table(aov.AUNCC)
#NS

aov.ISI_CoV<-stats_rename_all %>%
  group_by(time ) %>%
  anova_test(ISI_CoV_DoS ~ all_stat )
get_anova_table(aov.ISI_CoV)
#sig

#non-normal dist

kruskal.test(WMFR_DoS ~ all_stat, data = stats_rename_all, subset = time == "1 hr")
kruskal.test(WMFR_DoS ~ all_stat, data = stats_rename_all, subset = time == "24 hr")
#NS
kruskal.test(AE_DoS ~ all_stat, data = stats_rename_all, subset = time == "1 hr")
kruskal.test(AE_DoS ~ all_stat, data = stats_rename_all, subset = time == "24 hr")
#NS
kruskal.test(NBF_DoS ~ all_stat, data = stats_rename_all, subset = time == "1 hr")
kruskal.test(NBF_DoS ~ all_stat, data = stats_rename_all, subset = time == "24 hr")
#NS
kruskal.test(AUNCC_DoS ~ all_stat, data = stats_rename_all, subset = time == "1 hr")
kruskal.test(AUNCC_DoS ~ all_stat, data = stats_rename_all, subset = time == "1 day")
kruskal.test(AUNCC_DoS ~ all_stat, data = stats_rename_all, subset = time == "2 day")
kruskal.test(AUNCC_DoS ~ all_stat, data = stats_rename_all, subset = time == "3 day")
kruskal.test(AUNCC_DoS ~ all_stat, data = stats_rename_all, subset = time == "4 day")
kruskal.test(AUNCC_DoS ~ all_stat, data = stats_rename_all, subset = time == "5 day")
dunnTest(data=stats_rename_all, AUNCC_DoS ~ all_stat)


dunnTest(data=MCT3_nuc, Mean ~ rad_RP)

##combine Dunn's testing
#Dunn's test for each variable
MFR_tuk<-stats_rename_all %>%
  group_by(time) %>%
  tukey_hsd(MFR_DoS ~ all_stat )
MFR_tuk<-MFR_tuk%>% mutate(MEA_metric = "MFR_tuk")

WMFR_tuk<-stats_rename_all %>%
  group_by(time) %>%
  tukey_hsd(WMFR_DoS ~ all_stat )
WMFR_tuk<-WMFR_tuk%>% mutate(MEA_metric = "WMFR_tuk")


ISI_CoV_tuk<-stats_rename_all %>%
  group_by(time) %>%
  tukey_hsd(ISI_CoV_DoS ~ all_stat )
ISI_CoV_tuk<-ISI_CoV_tuk%>% mutate(MEA_metric = "ISI_CoV_tuk")

AUNCC_tuk<-stats_rename_all %>%
  group_by(time) %>%
  tukey_hsd(AUNCC_DoS ~ all_stat )
AUNCC_tuk<-AUNCC_tuk%>% mutate(MEA_metric = "AUNCC_tuk")

NBF_tuk<-stats_rename_all %>%
  group_by(time) %>%
  tukey_hsd(NBF_DoS ~ all_stat )
NBF_tuk<-NBF_tuk%>% mutate(MEA_metric = "NBF_tuk")

NumB_tuk<-stats_rename_all %>%
  group_by(time) %>%
  tukey_hsd(NumB_DoS ~ all_stat )
NumB_tuk<-NumB_tuk%>% mutate(MEA_metric = "NumB_tuk")


compile_tuk<-rbind(MFR_tuk,WMFR_tuk,ISI_CoV_tuk,AUNCC_tuk, NumB_tuk, NBF_tuk)
 compile_tuk_p.05<-compile_tuk%>%
   filter(p.adj<= 0.05)
 #combine all variables into a single df and filter for only the significant ones
 
 ####
 ###to make a dummy 0 hr timepoint so as to bump facets one slot and leave the first facet open so that key can be added
   
dummy_0hr<-all_Dos%>%
  filter(time== "1 hr")%>%
  mutate(time=str_replace_all(time,"1 hr" , "0.1 hr"))
dummy_0hr2<-rbind(all_Dos, dummy_0hr)
 
  ###################################################
#plots
timepoints<-timepoints%>% arrange(desc(RP_status))
levelx_order <- c("sham_vehicle", "sham_amifostine", "10_Gy_eq_vehicle", "10_Gy_eq_amifostine")
group_order<- c("vehicle", "amifostine" )
stat_comp<-list(c("sham_vehicle", "sham_amifostine" ), c("sham_vehicle","10_Gy_eq_vehicle"), c("sham_amifostine", "10_Gy_eq_amifostine" ), c("10_Gy_eq_vehicle", "10_Gy_eq_amifostine")  )
################################################################

#DoS plots
library(ggpubr)

AE_DoS_plot_dummy<- ggplot(dummy_0hr2,aes(x=all_stat, y=AE_DoS, fill=rad_stat, color = RP_status ))+
  geom_boxplot(aes(), outlier.shape = NA, show.legend = FALSE)+
  geom_jitter(width = 0.1, size= 0.75, show.legend = FALSE)+
  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
        axis.line = element_line(colour = "black"))+
  theme(axis.text.x = element_text(angle = -45, vjust = 0.5, hjust=0.0, size = 12, colour = "black"))+
  theme(axis.text.x =  element_blank())+
  theme(axis.text.y = element_text(angle = 00, vjust = 0.5, hjust=0, size = 12, colour = "black"))+
  theme(plot.title = element_text(hjust = 0.5))+              # title posisiton
  theme(plot.subtitle = element_text(hjust = 0.5)) +
  theme(axis.title.y = element_text(size = 14))+
  
  scale_x_discrete(limits = levelx_order)+
  scale_color_manual(values= c( "red", "black"))+
  scale_fill_manual(values= c("#999999", "#ecececec"))+
  ylim( c(-1.01, 1.01))+
  facet_wrap(~time)+
  geom_rect(aes(xmin = 2.5, xmax = Inf, ymin = -Inf, ymax = Inf),
            alpha = .01,
            fill = "grey", color= 0)+
  
  #stat_pvalue_manual(MFR_dunn, label = "p.adj.signif",  y.position = c(29, 35))+
  #The above line is not working right now and I suspect that this is because fill/ color are different metrics from all_stat
  #this could be really useful though!
  
  #stat_compare_means(method = "wilcox.test",paired = FALSE
  #                  ,comparisons = stat_comp,
  #                 label= "p.signif", hide.ns = TRUE,
  #                symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")))+
  
  labs(title="", subtitle = "", x="", y= "Active Electrode DoS", caption = "")

AE_DoS_plot_dummy
ggsave("AE_DoS_plot_dummy0.pdf",width = 1400,  height = 1400, units = "px")



#DoS plots
AE_DoS_plot<- ggplot(all_Dos,aes(x=all_stat, y=AE_DoS, fill=rad_stat, color = RP_status ))+
  geom_boxplot(aes(), outlier.shape = NA, show.legend = FALSE)+
  geom_jitter(width = 0.1, size= 0.75, show.legend = FALSE)+
  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
        axis.line = element_line(colour = "black"))+
  theme(axis.text.x = element_text(angle = -45, vjust = 0.5, hjust=0.0, size = 12, colour = "black"))+
  theme(axis.text.x =  element_blank())+
  theme(axis.text.y = element_text(angle = 00, vjust = 0.5, hjust=0, size = 12, colour = "black"))+
  theme(plot.title = element_text(hjust = 0.5))+              # title posisiton
  theme(plot.subtitle = element_text(hjust = 0.5)) +
  theme(axis.title.y = element_text(size = 14))+
  
  scale_x_discrete(limits = levelx_order)+
  scale_color_manual(values= c( "red", "black"))+
  scale_fill_manual(values= c("#999999", "#ecececec"))+
  ylim( c(-1.01, 1.01))+
  facet_wrap(~time)+
  geom_rect(aes(xmin = 2.5, xmax = Inf, ymin = -Inf, ymax = Inf),
            alpha = .01,
            fill = "grey", color= 0)+
  
  #stat_pvalue_manual(MFR_dunn, label = "p.adj.signif",  y.position = c(29, 35))+
#The above line is not working right now and I suspect that this is because fill/ color are different metrics from all_stat
  #this could be really useful though!
  
    #stat_compare_means(method = "wilcox.test",paired = FALSE
   #                  ,comparisons = stat_comp,
    #                 label= "p.signif", hide.ns = TRUE,
     #                symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")))+
  
  labs(title="", subtitle = "", x="", y= "Active Electrode DoS", caption = "")

AE_DoS_plot
ggsave("AE_DoS_plot.pdf",width = 1400,  height = 1400, units = "px")
#*

MFR_DoS_plot<- ggplot(all_Dos,aes(x=all_stat, y=MFR_DoS, fill=rad_stat, color = RP_status ))+
  geom_boxplot(aes(),outlier.shape = NA, show.legend = FALSE)+
  geom_jitter(width = 0.1, size= 0.75, show.legend = FALSE)+
  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
        axis.line = element_line(colour = "black"))+
  theme(axis.text.x = element_text(angle = -45, vjust = 0.5, hjust=0.0, size = 12, colour = "black"))+
  theme(axis.text.x =  element_blank())+
  theme(axis.text.y = element_text(angle = 00, vjust = 0.5, hjust=0, size = 12, colour = "black"))+
  theme(plot.title = element_text(hjust = 0.5))+              # title posisiton
  theme(plot.subtitle = element_text(hjust = 0.5)) +
  theme(axis.title.y = element_text(size = 14))+
  
  scale_x_discrete(limits = levelx_order)+
  scale_color_manual(values= c( "red", "black"))+
  scale_fill_manual(values= c("#999999", "#ecececec"))+
  ylim( c(-1.01, 1.01))+
  facet_wrap(~time)+
  geom_rect(aes(xmin = 2.5, xmax = Inf, ymin = -Inf, ymax = Inf),
            alpha = .01,
            fill = "grey", color= 0)+
  # stat_compare_means(method = "wilcox.test",paired = FALSE
  #                    ,comparisons = stat_comp,
  #                    label= "p.signif", hide.ns = TRUE,
  #                    symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")))+
  labs(title="", subtitle = "", x="", y= "MFR DoS", caption = "")

  
MFR_DoS_plot
ggsave("MFR_DoS_plot.pdf",width = 1400,  height = 1400, units = "px")

#*



WMFR_DoS_plot<- ggplot(all_Dos,aes(x=all_stat, y=WMFR_DoS, fill=rad_stat, color = RP_status ))+
  geom_boxplot(aes(),outlier.shape = NA, show.legend = FALSE)+
  geom_jitter(width = 0.1, size= 0.75, show.legend = FALSE)+
  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
        axis.line = element_line(colour = "black"))+
  theme(axis.text.x = element_text(angle = -45, vjust = 0.5, hjust=0.0, size = 12, colour = "black"))+
  theme(axis.text.x =  element_blank())+
  theme(axis.text.y = element_text(angle = 00, vjust = 0.5, hjust=0, size = 12, colour = "black"))+
  theme(plot.title = element_text(hjust = 0.5))+              # title posisiton
  theme(plot.subtitle = element_text(hjust = 0.5)) +
  theme(axis.title.y = element_text(size = 14))+
  
  scale_x_discrete(limits = levelx_order)+
  scale_color_manual(values= c( "red", "black"))+
  scale_fill_manual(values= c("#999999", "#ecececec"))+
  ylim( c(-1.01, 1.01))+
  facet_wrap(~time)+
  geom_rect(aes(xmin = 2.5, xmax = Inf, ymin = -Inf, ymax = Inf),
            alpha = .01,
            fill = "grey", color= 0)+
  # stat_compare_means(method = "wilcox.test",paired = FALSE
  #                   ,comparisons = stat_comp,
  #                  label= "p.signif", hide.ns = TRUE,
  #                 symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")))+
  labs(title="", subtitle = "", x="", y= "WMFR DoS", caption = "")

WMFR_DoS_plot
ggsave("WMFR_DoS_plot.pdf",width = 1400,  height = 1400, units = "px")

#*

AUNCC_DoS_plot<- ggplot(all_Dos,aes(x=all_stat, y=AUNCC_DoS, fill=rad_stat, color = RP_status ))+
  geom_boxplot(aes(),outlier.shape = NA, show.legend = FALSE)+
  geom_jitter(width = 0.1, size= 0.75, show.legend = FALSE)+
  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
        axis.line = element_line(colour = "black"))+
  theme(axis.text.x = element_text(angle = -45, vjust = 0.5, hjust=0.0, size = 12, colour = "black"))+
  theme(axis.text.x =  element_blank())+
  theme(axis.text.y = element_text(angle = 00, vjust = 0.5, hjust=0, size = 12, colour = "black"))+
  theme(plot.title = element_text(hjust = 0.5))+              # title posisiton
  theme(plot.subtitle = element_text(hjust = 0.5)) +
  theme(axis.title.y = element_text(size = 14))+
  
  scale_x_discrete(limits = levelx_order)+
  scale_color_manual(values= c( "red", "black"))+
  scale_fill_manual(values= c("#999999", "#ecececec"))+
  ylim( c(-1.01, 1.51))+
  facet_wrap(~time)+
  geom_rect(aes(xmin = 2.5, xmax = Inf, ymin = -Inf, ymax = Inf),
            alpha = .01,
            fill = "grey", color= 0)+
  # stat_compare_means(method = "wilcox.test",paired = FALSE
  #                    ,comparisons = stat_comp,
  #                    label= "p.signif", hide.ns = TRUE,
  #                    symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")))+
  labs(title="", subtitle = "", x="", y= "AUNCC DoS", caption = "")

AUNCC_DoS_plot
ggsave("AUNCC_DoS_plot.pdf",width = 1400,  height = 1400, units = "px")


#*
ISI_CoV_plot<- ggplot(all_Dos,aes(x=all_stat, y=ISI_CoV_DoS, fill=rad_stat, color = RP_status ))+
  geom_boxplot(aes(),outlier.shape = NA, show.legend = FALSE)+
  geom_jitter(width = 0.1, size= 0.75, show.legend = FALSE )+
  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
        axis.line = element_line(colour = "black"))+
  theme(axis.text.x = element_text(angle = -45, vjust = 0.5, hjust=0.0, size = 12, colour = "black"))+
  theme(axis.text.x =  element_blank())+
  theme(axis.text.y = element_text(angle = 00, vjust = 0.5, hjust=0, size = 12, colour = "black"))+
  theme(plot.title = element_text(hjust = 0.5))+              # title posisiton
  theme(plot.subtitle = element_text(hjust = 0.5)) +
  theme(axis.title.y = element_text(size = 14))+
  
  scale_x_discrete(limits = levelx_order)+
  scale_color_manual(values= c( "red", "black"))+
  scale_fill_manual(values= c("#999999", "#ecececec"))+
  ylim( c(-1.01, 1.01))+
  facet_wrap(~time)+
  geom_rect(aes(xmin = 2.5, xmax = Inf, ymin = -Inf, ymax = Inf),
            alpha = .01,
            fill = "grey", color= 0)+
  # stat_compare_means(method = "wilcox.test",paired = FALSE
  #                    ,comparisons = stat_comp,
  #                    label= "p.signif", hide.ns = TRUE,
  #                    symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")))+
  labs(title="", subtitle = "", x="", y= "ISI CoV DoS", caption = "")

ISI_CoV_plot
ggsave("ISI_CoV_plot.pdf",width = 1400,  height = 1400, units = "px")

#*
NumB_plot<- ggplot(all_Dos,aes(x=all_stat, y=ISI_CoV_DoS, fill=rad_stat, color = RP_status ))+
  geom_boxplot(aes(),outlier.shape = NA, show.legend = FALSE)+
  geom_jitter(width = 0.1, size= 0.75, show.legend = FALSE )+
  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
        axis.line = element_line(colour = "black"))+
  theme(axis.text.x = element_text(angle = -45, vjust = 0.5, hjust=0.0, size = 12, colour = "black"))+
  theme(axis.text.x =  element_blank())+
  theme(axis.text.y = element_text(angle = 00, vjust = 0.5, hjust=0, size = 12, colour = "black"))+
  theme(plot.title = element_text(hjust = 0.5))+              # title posisiton
  theme(plot.subtitle = element_text(hjust = 0.5)) +
  theme(axis.title.y = element_text(size = 14))+
  
  scale_x_discrete(limits = levelx_order)+
  scale_color_manual(values= c( "red", "black"))+
  scale_fill_manual(values= c("#999999", "#ecececec"))+
  ylim( c(-1.01, 1.01))+
  facet_wrap(~time)+
  geom_rect(aes(xmin = 2.5, xmax = Inf, ymin = -Inf, ymax = Inf),
            alpha = .01,
            fill = "grey", color= 0)+
  # stat_compare_means(method = "wilcox.test",paired = FALSE
  #                    ,comparisons = stat_comp,
  #                    label= "p.signif", hide.ns = TRUE,
  #                    symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")))+
  labs(title="", subtitle = "", x="", y= "Num. Burst DoS", caption = "")

NumB_plot
ggsave("NumB_plot.pdf",width = 1400,  height = 1400, units = "px")
#*
NBF_plot<- ggplot(all_Dos,aes(x=all_stat, y=NBF_DoS, fill=rad_stat, color = RP_status ))+
  geom_boxplot(aes(),outlier.shape = NA, show.legend = FALSE)+
  geom_jitter(width = 0.1,  size= 0.75, show.legend = FALSE)+
  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
        axis.line = element_line(colour = "black"))+
  theme(axis.text.x = element_text(angle = -45, vjust = 0.5, hjust=0.0, size = 12, colour = "black"))+
  theme(axis.text.x =  element_blank())+
  theme(axis.text.y = element_text(angle = 00, vjust = 0.5, hjust=0, size = 12, colour = "black"))+
  theme(plot.title = element_text(hjust = 0.5))+              # title posisiton
  theme(plot.subtitle = element_text(hjust = 0.5)) +
  theme(axis.title.y = element_text(size = 14))+
  
  scale_x_discrete(limits = levelx_order)+
  scale_color_manual(values= c( "red", "black"))+
  scale_fill_manual(values= c("#999999", "#ecececec"))+
  ylim( c(-1.2, 1.01))+
  facet_wrap(~time)+
  geom_rect(aes(xmin = 2.5, xmax = Inf, ymin = -Inf, ymax = Inf),
            alpha = .01,
            fill = "grey", color= 0)+
  # stat_compare_means(method = "wilcox.test",paired = FALSE
  #                    ,comparisons = stat_comp,
  #                    label= "p.signif", hide.ns = TRUE,
  #                    symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")))+
  labs(title="", subtitle = "", x="", y= "NBF DoS", caption = "")

NBF_plot
ggsave("NBF_plot_nas_to_0.pdf",width = 1400,  height = 1400, units = "px")


###

#results summary matrix

MCT8_MEA_results_summary<- all_Dos%>%
  group_by(RP_status,rad_stat, time)%>%
  summarise(mean_AEY_summary= mean(AE_DoS, na.rm = TRUE),
            mean_AEY_summary_sd= sd(AE_DoS, na.rm = TRUE),
            mean_WMFR_summary= mean(WMFR_DoS, na.rm = TRUE),
            mean_WMFR_summary_sd= sd(AE_DoS, na.rm = TRUE),
            mean_AUNCC_summary= mean(AUNCC_DoS, na.rm = TRUE),
            mean_AUNCC_summary_sd= sd(AUNCC_DoS, na.rm = TRUE),
            mean_ISI_CoV_summary= mean(ISI_CoV_DoS, na.rm = TRUE),
            mean_ISI_CoV_summary_sd= sd(ISI_CoV_DoS, na.rm = TRUE),
            mean_NumB_summary= mean(NumB_DoS, na.rm = TRUE),
            mean_NumB_summary_sd= sd(NumB_DoS, na.rm = TRUE),
            mean_NBF_summary= mean(NBF_DoS, na.rm = TRUE),
            mean_NBF_summary_sd= sd(NBF_DoS, na.rm = TRUE))
#to save as CSV
wd<-getwd()
project<- "MCT8_MEA_results"
filename_csv<-paste(project,"_summary",".csv", sep = "") #concatenate file name + "_summary" ".csv"
filename_csv<-file.path(wd,filename_csv) #make file path
write.csv(MCT8_MEA_results_summary,filename_csv)#save to WD as .csv
#*
#compile_stats
filename_csv<-paste(project,"_summary",".csv", sep = "") #concatenate file name + "_summary" ".csv"
filename_csv<-file.path(wd,filename_csv) #make file path
write.csv(MCT8_MEA_results_summary,filename_csv)#save to WD as .csv
#*

########################################################################
#############################################End
########################################################################

