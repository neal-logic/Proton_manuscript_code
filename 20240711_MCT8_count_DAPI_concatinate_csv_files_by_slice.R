library(tidyverse)
#this files concatinates the CSV output from the GFAP dataset run in FIJI
#it takes the DAPI count from every 3rd slice and summs them for each image.
#this is also done for the TUJ1 imaging dataset
#the two datasets are the combine and the mead and sd are taken for each well. 

#######################First chunk##################
#GFAP dataset
setwd("C:/Users/Neal/OneDrive - UMass Lowell/Pain Research Lab/Grants/NeuroRad_GMHA (2024)/Data/MCT8/MCT8_secreted_markers_luminex/MCT8_count_results_GFAP_1/")


#makes a list of the file names in the working directory
files = list.files()
#lapply seems to apply read.table to the large list "files"
#creates a large list from files in the wd
data2=lapply(files, read.table, header=TRUE, sep=",")
for (i in 1:length(data2)){data2[[i]]<-cbind(data2[[i]],files[i])}
data_rbind <- do.call("rbind", data2) 
#colnames(data_rbind)[c(1,2,3,4)]<-c("name", "sex", "amount","year")

#shorten the file name
data_rbind<-data_rbind %>%
  mutate(`files[i]`= str_sub(`files[i]`,0,8), remove = FALSE)
proj<-"MCT8"
MCT8_GFAP_count_<-data_rbind

## 
#Extract slice
newcols<-c("Project", "Objective", "junk", "rad_stat", "sample", "image")
MCT8_GFAP_count_<-MCT8_GFAP_count_ %>%
  separate(Slice, into= newcols, sep = "_", remove=FALSE  )
#second separation to acquire plane
newcols2<-c("Junk", "plane_chunk", "Project2")
MCT8_GFAP_count_<-MCT8_GFAP_count_ %>%
  separate(Project, into= newcols2, sep = "/", remove=FALSE  )
#final separation to acquire plane
newcols3<-c("Junk", "plane")
MCT8_GFAP_count_<-MCT8_GFAP_count_ %>%
  separate(plane_chunk, into= newcols3, sep = ":", remove=TRUE  )


#mutate for correct well name
MCT8_GFAP_count_<-MCT8_GFAP_count_ %>%
  unite("well_ID",rad_stat,sample,sep = "_", remove = FALSE)%>%
  mutate(well_ID = str_replace_all(well_ID, "10Gyeq_1", "A5"))%>%
  mutate(well_ID = str_replace_all(well_ID, "10Gyeq_2", "B5"))%>%
  mutate(well_ID = str_replace_all(well_ID, "10Gyeq_3", "C5"))%>%
  mutate(well_ID = str_replace_all(well_ID, "sham_1", "D5"))%>%
  mutate(well_ID = str_replace_all(well_ID, "sham_2", "E5"))%>%
  mutate(well_ID = str_replace_all(well_ID, "sham_3", "F5"))


##
#make condition identifier
MCT8_GFAP_count_<-MCT8_GFAP_count_ %>%
  unite("well",sample,image, image,sep = "_", remove = FALSE)

##
#image identifier
MCT8_GFAP_count_<-MCT8_GFAP_count_ %>%
  mutate(Image_ID= str_sub(Slice,-25,-1))
##
#slice keep interval
slice_want<-seq(from = 1, to = 34, by = 3)

MCT8_GFAP_count_by_3_slice<-MCT8_GFAP_count_%>%
  filter(plane %in% slice_want)


#every 3rd slice count summary
MCT8_GFAP_count_by_3_summary<-MCT8_GFAP_count_by_3_slice%>%
  group_by(well_ID, image)%>%
  summarise(three_slice_nuclei_sum = sum(Count), rad_stat= first(rad_stat))%>% #note that the rad_stat column must be kept as an identifier
  mutate(ICC_session = "GFAP")


MCT8_GFAP_count_by_3_summary_all<-MCT8_GFAP_count_by_3_summary%>%
  group_by(well_ID)%>%
  summarise(three_slice_nuclei_mean = mean(three_slice_nuclei_sum),
            three_slice_nuclei_sd = sd(three_slice_nuclei_sum))


MCT8_GFAP_count_by_3_summary_all<-MCT8_GFAP_count_by_3_summary_all%>%
  mutate(total_nuclei_per_construct = ((three_slice_nuclei_mean/0.0648)*5.65))  #.581×.581×.192 =0.0648 mm3 for slice vol, total construct vol = 5.65 mm^3

########second chunk############################


#TUJ1 dataset
setwd("C:/Users/Neal/OneDrive - UMass Lowell/Pain Research Lab/Grants/NeuroRad_GMHA (2024)/Data/MCT8/MCT8_secreted_markers_luminex/MCT8_count_results_TUJ1_1/")



#makes a list of the file names in the working directory
files = list.files()
#lapply seems to apply read.table to the large list "files"
#creates a large list from files in the wd
data2=lapply(files, read.table, header=TRUE, sep=",")
for (i in 1:length(data2)){data2[[i]]<-cbind(data2[[i]],files[i])}
data_rbind <- do.call("rbind", data2) 
#colnames(data_rbind)[c(1,2,3,4)]<-c("name", "sex", "amount","year")

#shorten the file name
data_rbind<-data_rbind %>%
  mutate(`files[i]`= str_sub(`files[i]`,0,8), remove = FALSE)
proj<-"MCT8_TUJ1"
MCT8_TUJ1_count_<-data_rbind

## 
#Extract slice
newcols<-c("Project", "Objective", "junk", "rad_stat", "sample", "image")
MCT8_TUJ1_count_<-MCT8_TUJ1_count_ %>%
  separate(Slice, into= newcols, sep = "_", remove=FALSE  )
#second separation to acquire plane
newcols2<-c("Junk", "plane_chunk", "Project2")
MCT8_TUJ1_count_<-MCT8_TUJ1_count_ %>%
  separate(Project, into= newcols2, sep = "/", remove=FALSE  )
#final separation to acquire plane
newcols3<-c("Junk", "plane")
MCT8_TUJ1_count_<-MCT8_TUJ1_count_ %>%
  separate(plane_chunk, into= newcols3, sep = ":", remove=TRUE  )


#mutate for correct well name
MCT8_TUJ1_count_<-MCT8_TUJ1_count_ %>%
  unite("well_ID",rad_stat,sample,sep = "_", remove = FALSE)%>%
  mutate(well_ID = str_replace_all(well_ID, "10Gyeq_1", "A5"))%>%
  mutate(well_ID = str_replace_all(well_ID, "10Gyeq_2", "B5"))%>%
  mutate(well_ID = str_replace_all(well_ID, "10Gyeq_3", "C5"))%>%
  mutate(well_ID = str_replace_all(well_ID, "sham_1", "D5"))%>%
  mutate(well_ID = str_replace_all(well_ID, "sham_2", "E5"))%>%
  mutate(well_ID = str_replace_all(well_ID, "sham_3", "F5"))


##
#make condition identifier
MCT8_TUJ1_count_<-MCT8_TUJ1_count_ %>%
  unite("well",sample,image, image,sep = "_", remove = FALSE)

##
#image identifier
MCT8_TUJ1_count_<-MCT8_TUJ1_count_ %>%
  mutate(Image_ID= str_sub(Slice,-25,-1))
##
#slice keep interval
slice_want<-seq(from = 1, to = 34, by = 3)

MCT8_TUJ1_count_by_3_slice<-MCT8_TUJ1_count_%>%
  filter(plane %in% slice_want)


#every 3rd slice count
MCT8_TUJ1_count_by_3_summary<-MCT8_TUJ1_count_by_3_slice%>%
  group_by(well_ID, image)%>%
  summarise(three_slice_nuclei_sum = sum(Count), rad_stat= first(rad_stat))%>% #note that the rad_stat column must be kept as an identifier
  mutate(ICC_session = "TUJ1")

MCT8_TUJ1_count_by_3_summary_all<-MCT8_TUJ1_count_by_3_summary%>%
  group_by(well_ID)%>%
  summarise(three_slice_nuclei_mean = mean(three_slice_nuclei_sum),
            three_slice_nuclei_sd = sd(three_slice_nuclei_sum))


MCT8_TUJ1_count_by_3_summary_all<-MCT8_TUJ1_count_by_3_summary_all%>%
  mutate(total_nuclei_per_construct = ((three_slice_nuclei_mean/0.0648)*5.65))  #.581×.581×.192 =0.0648 mm3 for slice vol, total construct vol = 5.65 mm^3


#need to figure out how to alphabetize the rows to match count from z-stack count
#############third chunk#########################
#combine datasets

combine_both_MCT8_imaging_datasets_nuclei_count<- rbind(MCT8_GFAP_count_by_3_summary, MCT8_TUJ1_count_by_3_summary )


combine_both_MCT8_imaging_datasets_nuclei_count_summary<-combine_both_MCT8_imaging_datasets_nuclei_count%>%
  group_by(well_ID)%>%
  summarise(three_slice_nuclei_mean = mean(three_slice_nuclei_sum),
            three_slice_nuclei_sd = sd(three_slice_nuclei_sum))

combine_both_MCT8_imaging_datasets_nuclei_count_summary<-combine_both_MCT8_imaging_datasets_nuclei_count_summary%>%
  mutate(total_nuclei_per_construct = ((three_slice_nuclei_mean/0.0648)*5.65))  #.581×.581×.192 =0.0648 mm3 for slice vol, total construct vol = 5.65 mm^3
  

wd<-getwd()
filename_csv<-paste(proj,"_concatinated_DAPI_counts_by_slice_1",".csv", sep = "") #concatenate file name + "_summary" ".csv"
filename_csv<-file.path(wd,filename_csv) #make file path
write.csv(combine_both_MCT8_imaging_datasets_nuclei_count_summary,filename_csv)#save to WD as .csv
###########chunk 4 ############################
#extract cell counts without taking mean of ROIs

all_ROI_est_cell_num<-combine_both_MCT8_imaging_datasets_nuclei_count%>%
  mutate(total_nuclei_per_construct = ((three_slice_nuclei_sum/0.0648)*5.65))%>%
  mutate(rad_stat = as.factor(rad_stat))%>%
  unite("well_ID_by_MCT", well_ID, ICC_session, sep = "_", remove = FALSE)


#plots
#plot for total cell number in 3D CTM samples
library(ggpubr)
levelx_order4 <- c("sham", "10Gyeq")

MCT8_cell_num_3D<- ggplot(all_ROI_est_cell_num,aes(x=rad_stat, y=total_nuclei_per_construct, fill= rad_stat ))+
  geom_boxplot(aes(), outlier.shape = NA, show.legend = FALSE)+
  geom_jitter(aes(color = well_ID_by_MCT), width = 0.1, size= 1.5, show.legend = FALSE)+
  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
        axis.line = element_line(colour = "black"))+
  theme(axis.text.x = element_text(angle = -0, vjust = 0.5, hjust=0.5, size = 16, colour = "black"))+
  #theme(axis.text.x =  element_blank())+
  theme(axis.text.y = element_text(angle = 00, vjust = 0.5, hjust=0, size = 16, colour = "black"))+
  theme(plot.title = element_text(hjust = 0.5))+              # title posisiton
  theme(plot.subtitle = element_text(hjust = 0.5)) +
  theme(axis.title.y = element_text(size = 16))+
  
  scale_x_discrete(limits = levelx_order4)+
  scale_fill_manual(values= c("#999999", "#ecececec"))+
  ylim( c(0, 35000))+
  geom_rect(aes(xmin = 1.5, xmax = Inf, ymin = -Inf, ymax = Inf),
            alpha = .01,
            fill = "grey", color= 0)+


   labs(title="", subtitle = "", x="", y= "Nuclei per construct", caption = "MCT8 3D luminex, nuclei per construct")

MCT8_cell_num_3D

#stats
library(rstatix)
all_test_norm_3D<-all_ROI_est_cell_num%>%
  group_by(rad_stat) %>%
  shapiro_test(total_nuclei_per_construct)
#both are normal

all_ROI_est_cell_num<-all_ROI_est_cell_num%>%
  mutate(rad_stat = desc(rad_stat))
#ttest

t.test(total_nuclei_per_construct ~ rad_stat, data = all_ROI_est_cell_num)
#no statistical difference, p = 0.06676