library(tidyverse)


setwd("C:/Users/Neal/OneDrive - UMass Lowell/Pain Research Lab/Grants/NeuroRad_GMHA (2023)/Data/MCT5_20231116/20240820_53BP1_gH2AX/Z_profile_intensity/gH2AX_intensity_by_depth_result")


files = list.files()
#lapply seems to apply read.table to the large list "files"
#creates a large list from files in the wd
data2=lapply(files, read.table, header=TRUE, sep=",")
for (i in 1:length(data2)){data2[[i]]<-cbind(data2[[i]],files[i])}
data_rbind <- do.call("rbind", data2) 
#colnames(data_rbind)[c(1,2,3,4)]<-c("name", "sex", "amount","year")

mmm<-data_rbind

## 
#Extract slice
newcols<-c("A", "B", "C", "D", "E")
newcols2<- c("GG", "G")
newcols2.1<- c("slice", "junk")
newcols3<- c("date", "project", "objective", "sample", "image")
MCT5_Z<- mmm %>%
  separate(Label, into= newcols, sep = ":", remove=TRUE  )%>%
  separate(E, into= newcols2, sep = "-", remove=TRUE  )%>%
  separate(GG, into= newcols2.1, sep = "/", remove=TRUE  )%>%
  separate(G, into= newcols3, sep = "_", remove=TRUE  )

MCT5_Z<-MCT5_Z%>%
  mutate(slice=as.numeric(slice))%>%
  mutate(depth = (slice*4))%>%
  unite("replicate",sample,image, sep = "_", remove = FALSE)

setwd("C:/Users/Neal/OneDrive - UMass Lowell/Pain Research Lab/Grants/NeuroRad_GMHA (2023)/Data/MCT5_20231116/20240820_53BP1_gH2AX/Z_profile_intensity/gH2AX_background_results")
files = list.files()
#lapply seems to apply read.table to the large list "files"
#creates a large list from files in the wd
data2=lapply(files, read.table, header=TRUE, sep=",")
for (i in 1:length(data2)){data2[[i]]<-cbind(data2[[i]],files[i])}
data_rbind <- do.call("rbind", data2) 
#colnames(data_rbind)[c(1,2,3,4)]<-c("name", "sex", "amount","year")

bbb<-data_rbind

## 
#Extract bg slice
MCT5_Z_bg<- bbb %>%
  separate(Label, into= newcols, sep = ":", remove=TRUE  )%>%
  separate(E, into= newcols2, sep = "-", remove=TRUE  )%>%
  separate(GG, into= newcols2.1, sep = "/", remove=TRUE  )%>%
  separate(G, into= newcols3, sep = "_", remove=TRUE  )
MCT5_Z_bg<-MCT5_Z_bg%>%
  mutate(Mean_bg = Mean)%>%
  mutate(slice=as.numeric(slice))
MCT5_Z_bg<-MCT5_Z_bg%>%
  select(6,11,12,18)%>%
  filter( sample!= "ctrl")
#combine by index nuclear measurments with background
MCT5_Z_comb<-left_join(MCT5_Z, MCT5_Z_bg)
#subtract background
MCT5_Z_comb<-MCT5_Z_comb%>%
  mutate(Mean_sub_bg= (Mean-Mean_bg))
#adjust depth to zero the min
MCT5_Z_comb<-MCT5_Z_comb%>%
  group_by(sample, image)%>%
  mutate(depth_min= min(depth))
#find zeroed de
MCT5_Z_comb<-MCT5_Z_comb%>%
  mutate(zeroed_depth = (depth-depth_min))


#Changes names of treatment
MCT5_Z_comb<-MCT5_Z_comb%>%
  mutate(sample=as.character(sample))%>% #need to make treatment a factor so that we can order x axis correctly
  mutate(treatment = str_replace_all(sample, c( "A"= "10 Gy", "B"="10 Gy", "C"= "Sham", "D"="10 Gy","E"="Sham","F"= "Sham") ))%>%
  filter(sample != "ctrl") #remove control wells from plot

#second separation to acquire plane
#results are product of mean - background
summary_MCT5_Z<-MCT5_Z_comb%>%
  group_by(treatment, sample,  zeroed_depth)%>%
  summarise(mean_gH2AX_by_slice= mean(Mean_sub_bg),
            sd_gH2AX_by_slice= sd(Mean_sub_bg),
            num_cells = n())
  #the below command was removed as was the "image" column in the group by as we are taking the mean of the two ROIs for each sample
  #unite("replicate",sample,image, sep = "_", remove = FALSE)
  
#all mean of mean. background has been already subtracted
all_mean<-summary_MCT5_Z%>%
  group_by(treatment,zeroed_depth)%>%
  summarise(treatment_mean = mean(mean_gH2AX_by_slice),
            treatment_sd = sd(mean_gH2AX_by_slice))
######
#########import for second axis which lots Bragg peak
TOP_Tissue_pLET <- read_csv("C:/Users/Neal/OneDrive - UMass Lowell/Proton Manuscript June 2024/NR_Proton_figure2_TOPASnBIO/From_Issac_meyer_updated/SampleTissue_pLET_converted.csv")
TOP_Tissue_pLET<-TOP_Tissue_pLET%>%
  mutate(treatment = ("10 Gy"))%>%
  filter(x<= 180)
#Reverse Bragg peak
TOP_Tissue_pLET <- TOP_Tissue_pLET%>%
  mutate(x = rev(x))

#plot

d3_plot<- ggplot()+
  geom_point(data = summary_MCT5_Z, aes(x=zeroed_depth, y=mean_gH2AX_by_slice, color = sample),size = 0.7,  show.legend = FALSE)+
  geom_line(data = summary_MCT5_Z, aes(x=zeroed_depth, y=mean_gH2AX_by_slice,  color = sample), linewidth = 1.0 ,  show.legend = FALSE )+
  geom_line(data =TOP_Tissue_pLET , aes(x = x, y = (value/3) ), linewidth = 1.0, linetype= 5)+
  #geom_line(data = all_mean, aes(x = depth, y = treatment_mean ),  linewidth = 1.0, linetype= 5 )+
  # geom_ribbon(aes(y = mean_gH2AX_by_slice, ymin = mean_gH2AX_by_slice - sd_gH2AX_by_slice, ymax = mean_gH2AX_by_slice + sd_gH2AX_by_slice,
  #                 fill = replicate), alpha = .2 , show.legend = FALSE) +
  # geom_errorbar(aes(ymin=mean_heal -sd_heal, ymax=mean_heal +sd_heal, color = cond),linewidth = 0.75,
  #               width= 1, position=position_dodge(0.2),  show.legend = FALSE)+
  # #geom_quantile(aes(color = cond), linetype = 2, quantiles = 0.5)+ #quantile for trajectory
  
  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
        axis.line = element_line(colour = "black"))+
  theme(axis.text.x = element_text(angle = -0, vjust = 0.5, hjust=0.5, size = 18, colour = "black"))+
  theme(axis.text.y = element_text(angle = 00, vjust = 0.5, hjust=0, size = 18, colour = "black"))+
  theme(plot.title = element_text(hjust = 0.5))+              # title posisiton
  theme(plot.subtitle = element_text(hjust = 0.5)) +
  theme(axis.title.y = element_text(size = 18))+
  theme(axis.title.x = element_text(size = 18))+
  facet_wrap(~treatment,ncol= 1)+
  #ylim(-1,30)+
  
  scale_color_manual(values= c(  "green","green", "grey", "green","grey","grey"))+
  #second axis
  scale_y_continuous( 
    sec.axis = sec_axis( ~. *3 , name = "Proton dose averaged LET (keV/μm) "))+ 
  
  labs(title="", subtitle = "", x=" Depth (μm)", y= "Mean nuclear γH2AX intensity", caption = "")


d3_plot



#end of code
###############
###################
########################

#slice keep interval
slice_want<-seq(from = 1, to = 34, by = 3)

MCT8_NeuN_by_3_slice<-MCT8_NeuN_count%>%
  filter(slice %in% slice_want)

NeuN_mean_threshold<-87 #above mean NeuN value in nuclei for this threshold number in "NeuN +" 
##this is the threshold by which we determine whether a cell is NeuN positive of not
#this was calculated by finding mean of posistve cells and negative cells
#determine threshold  
#used slices from sample 1 img 1. index nuclei to their number in column X from dataframe

MCT8_NeuN_by_3_slice_NeuN_positive<-MCT8_NeuN_by_3_slice%>%
  filter(Mean>NeuN_mean_threshold)

#filter for NeuN positive cells, consolidate replicates
all_positive_nuclei_summary<-MCT8_NeuN_by_3_slice_NeuN_positive%>%
  group_by(treatment,sample)%>%
  summarise(NeuN_nuclei = n())
#all cells, consolidate replicates
all_nuclei_summary<-MCT8_NeuN_by_3_slice%>%
  group_by(treatment,sample)%>%
  summarise(total_nuclei = n())
#join the NeuN cell and the all-cell datasets 
percent_NeuN<-  left_join(all_nuclei_summary, all_positive_nuclei_summary, by = c("treatment", "sample"))
#divide NeuN cells by total cells
percent_NeuN<- percent_NeuN%>%
  mutate(percent_NeuN_positive_cells = (NeuN_nuclei/ total_nuclei)*100)

#final export numbers by treatment group summary
summary_percent_NeuN<- percent_NeuN%>%
  group_by(treatment)%>%
  summarise(mean_percent_NeuN_cells= mean(percent_NeuN_positive_cells),
            sd_percent_NeuN_cells= sd(percent_NeuN_positive_cells),
            sum_total_cells= sum(total_nuclei))



#Statistics
library(rstatix)
test_norm<-percent_NeuN%>%
  group_by(treatment) %>%
  shapiro_test(percent_NeuN_positive_cells)
#can't reject normality

t.test(percent_NeuN_positive_cells~ treatment, data = percent_NeuN)
#p = 0.5911
