library(tidyverse)
#for both the second ELISA
#note that plots for cell number for sup. figures can be found in the R file: 20240711_MCT8_count_DAPI_concatinate_csv_files_by_slice


MCT8 <- read_csv("C:/Users/Neal/OneDrive - UMass Lowell/Pain Research Lab/Grants/NeuroRad_GMHA (2024)/Data/MCT8/MCT8_secreted_markers_luminex/20240703_processed_Luminex_data_MCT8.csv")
MCT8_20240710 <- read_csv("C:/Users/Neal/OneDrive - UMass Lowell/Pain Research Lab/Grants/NeuroRad_GMHA (2024)/Data/MCT8/MCT8_secreted_markers_luminex/20240710_luminex_vol_corrected/20240710_MCT8_luminex_3D.csv")
#reading in 2D cell count data which is integrated into a single CSV elsewhere. 

#set WD
setwd("C:/Users/Neal/OneDrive - UMass Lowell/Pain Research Lab/Grants/NeuroRad_GMHA (2024)/Data/MCT8/MCT8_2D_DAPI_count/20240719_MCT8_DAPI_cell_count_result")
#import 2D cell count dataset as multiple CSVs
files = list.files()
#lapply seems to apply read.table to the large list "files"
#creates a large list from files in the wd
data2=lapply(files, read.table, header=TRUE, sep=",")
for (i in 1:length(data2)){data2[[i]]<-cbind(data2[[i]],files[i])}
data_rbind <- do.call("rbind", data2) 
newcols<-c("folder", "file")
newcols1<-c("Project", "Objective", "junk", "rad_stat", "sample", "image")
MCT8_2D_count<-data_rbind %>% #split naming column into individual columns
  separate(Slice, into= newcols, sep = "-", remove=FALSE  )%>%
  separate(file, into= newcols1, sep = "_", remove=FALSE  )
#take only the columns we need from 2D count dataset
MCT8_2D_count_narrow<- MCT8_2D_count%>%
  select(sample, Count)
#plot histogram of cell distribution
hist(MCT8_2D_count_narrow$Count,xlim = c(800, 2800), breaks = 20, main= "", xlab = "DAPI nuclei per image")

MCT8.1<-MCT8%>%
  mutate(IL1B=as.numeric(IL1B))%>%
  mutate(IL6=as.numeric(IL6))%>%
  mutate(KC=as.numeric(KC))%>%
  mutate(MCP1=as.numeric(MCP1))%>%
  mutate(MIP2=as.numeric(MIP2))%>%
  mutate(TNFA=as.numeric(TNFA))

#binds together both luminex datasets
MCT8.1<-rbind(MCT8.1,MCT8_20240710)
#joins the 2D nuclei count with the luminex data
MCT8.1<-left_join(MCT8.1, MCT8_2D_count_narrow, by = c("Exp_well" = "sample" ))

#adjust count for total well cell number. call this estimated cell count
#1.163^2 = 1.35 mm^2, this is 19.1% of total estimated seeded area
#total seeded area 3 mm = diameter, (1.5^2)*3.14 = 7.065 mm^2
#multiply counted cells by 5.22
MCT8.1<-MCT8.1%>%
  mutate(est_count = (Count*5.22))

MCT8_rel_wide<-MCT8.1%>%
  filter(sample_type == "sample",
         time != "media_2")%>%
  pivot_wider(names_from =  time, values_from = c(IL1B,IL6,KC,MCP1,MIP2,TNFA))

#create med_1 (pre-treatment) df
med_1<-MCT8.1%>%
  filter(time == "media_1",
         sample_type =="sample")%>%
  pivot_wider(names_from =  time, values_from = c(IL1B,IL6,KC,MCP1,MIP2,TNFA))

#create 48hr df, remove empty well
hr48<-MCT8.1%>%
  filter(time == "48_hr",
         sample_type =="sample", #remove standards
         )%>%
  pivot_wider(names_from =  time, values_from = c(IL1B,IL6,KC,MCP1,MIP2,TNFA))

#remove wells in 48 hr timepoint that will be duplicate in med_2 timepoint
hr48<-hr48%>%
  select(-(1:5),-(7:12))
#change name of column 1 to avoid creating duplicates when merging
colnames(hr48)[1] <- "Exp.well"

#NOTE: for other scripts great application of SAPPLY across a df
hr48[,2:7] <- sapply(hr48[,2:7],as.numeric)#need to make columns numeric
med_1[,13:18] <- sapply(med_1[,13:18],as.numeric)#need to make columns numeric


#merge med_2 and 48hr dfs
alltg<-cbind(hr48, med_1) #combine the two timepoint dfs
#remove first column of combine dataset to avoid redundancy 
alltg<-alltg%>%
  select(-1)

#NEXT: need to divide all cytokine values (med_1 and 48 hr) by adjusted cell numbers before calculate relative difference
#multiply by 1000 to convert from pg to fg
pg_cell<-alltg%>%
  mutate(IL6_48_hr = (IL6_48_hr*0.2/ est_count*1000),
           KC_48_hr= (KC_48_hr*0.2/ est_count*1000),
           MCP1_48_hr = (MCP1_48_hr*0.2/ est_count*1000),
           MIP2_48_hr = (MIP2_48_hr*0.2/ est_count*1000),
           TNFA_48_hr = (TNFA_48_hr*0.2/ est_count*1000),
           IL6_media_1 = (IL6_media_1*0.2/ est_count*1000),
           KC_media_1 = (KC_media_1*0.2/ est_count*1000),
           MCP1_media_1 =(MCP1_media_1*0.2/ est_count*1000),
           MIP2_media_1 =(MIP2_media_1*0.2/ est_count*1000),
           TNFA_media_1 = (TNFA_media_1*0.2/ est_count*1000) )
#PG/cell absolute med_1
pg_cell_48hr<-pg_cell%>%
  select(1:18)

hr48_abs_lng<-pg_cell_48hr%>%
  unite("all_stat" , rad_stat, amifostine_stat, sep = "_", remove = FALSE)%>%
  pivot_longer( cols = c(IL6_48_hr, KC_48_hr, MCP1_48_hr,MIP2_48_hr ), 
                names_to = "cytokine", values_to = "cell_adj_val")

pg_cell_med_1<-pg_cell%>%
  select(7:24)

med_1_abs_lng<-pg_cell_med_1%>%
  unite("all_stat" , rad_stat, amifostine_stat, sep = "_", remove = FALSE)%>%
  pivot_longer( cols = c(IL6_media_1, KC_media_1, MCP1_media_1,MIP2_media_1 ), 
                names_to = "cytokine", values_to = "cell_adj_val")
###
#extract mean and SD by condition for cell count for export
sum_of_cell_count<-alltg%>%
  group_by(rad_stat, amifostine_stat)%>%
  summarise(mean_group  = mean(Count),
            sd_group  = sd(Count))
## export data for supplemental table 1
proj<-"MCT8"
wd<-getwd()
filename_csv<-paste(proj,"_Luminex_supplemental_2D",".csv", sep = "") #concatenate file name + "_summary" ".csv"
filename_csv<-file.path(wd,filename_csv) #make file path
write.csv(pg_cell_48hr,filename_csv)#save to WD as .csv

#plots###
#######
#order variables for plotting
levelx_order <- c("Sham_No_RP", "Sham_amifostine", "10_Gy_No_RP", "10_Gy_amifostine")
group_order<- c("vehicle", "amifostine" )
stat_comp<-list(c("sham_vehicle", "sham_amifostine" ), c("sham_vehicle","10_Gy_eq_vehicle"), c("sham_amifostine", "10_Gy_eq_amifostine" ), c("10_Gy_eq_vehicle", "10_Gy_eq_amifostine")  )

#plot cell count by condition
hr48_abs_lng_IL6 <-hr48_abs_lng%>%
  filter(cytokine == "IL6_48_hr")
cell_count_plot_2D<- ggplot(hr48_abs_lng_IL6,aes(x=all_stat, y=est_count, fill= rad_stat, color = amifostine_stat ))+
  geom_boxplot(aes(), outlier.shape = NA, show.legend = FALSE)+
  geom_jitter(width = 0.1, size= 0.75, show.legend = FALSE)+
  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
        axis.line = element_line(colour = "black"))+
  theme(axis.text.x = element_text(angle = -45, vjust = 0.5, hjust=0.0, size = 16, colour = "black"))+
  theme(axis.text.x =  element_blank())+
  theme(axis.text.y = element_text(angle = 00, vjust = 0.5, hjust=0, size = 16, colour = "black"))+
  theme(plot.title = element_text(hjust = 0.5))+              # title posisiton
  theme(plot.subtitle = element_text(hjust = 0.5)) +
  theme(axis.title.y = element_text(size = 16))+
  
  scale_x_discrete(limits = levelx_order)+
  scale_color_manual(values= c( "red", "black"))+
  scale_fill_manual(values= c("#999999", "#ecececec"))+
  ylim( c(0, 15000))+
  geom_rect(aes(xmin = 2.5, xmax = Inf, ymin = -Inf, ymax = Inf),
            alpha = .01,
            fill = "grey", color= 0)+
  
  
  labs(title="", subtitle = "", x="", y= "Nuclei per well", caption = "MCT8 2D cell count by ICC for DAPI, 7 days post radiation")
cell_count_plot_2D


###
#

MCT8_hr48_abs<- ggplot(hr48_abs_lng,aes(x=all_stat, y=cell_adj_val, fill= rad_stat, color = amifostine_stat ))+
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
  #ylim( c(-1.01, 1.01))+
  facet_wrap(~cytokine, scales = "free")+
  geom_rect(aes(xmin = 2.5, xmax = Inf, ymin = -Inf, ymax = Inf),
            alpha = .01,
            fill = "grey", color= 0)+
  

  labs(title="", subtitle = "", x="", y= "Cytokine per cell (fg)", caption = "MCT8 2D luminex, 48hr normalized to cell number")

MCT8_hr48_abs
##################3

MCT8_med_1_abs<- ggplot(med_1_abs_lng,aes(x=all_stat, y=cell_adj_val, fill= rad_stat, color = amifostine_stat ))+
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
  #ylim( c(-1.01, 1.01))+
  facet_wrap(~cytokine, scales = "free")+
  geom_rect(aes(xmin = 2.5, xmax = Inf, ymin = -Inf, ymax = Inf),
            alpha = .01,
            fill = "grey", color= 0)+
  

  labs(title="", subtitle = "", x="", y= "Cytokine per cell (fg)", caption = "MCT8 2D luminex, Pre-treatment media normalized to cell number")

MCT8_med_1_abs





#relative cytokine
#create DoS ((post-pre)/(post+pre))
alltg_DoS<-pg_cell%>%
  mutate(IL6 = ((IL6_48_hr -IL6_media_1) / (IL6_48_hr +IL6_media_1)),
         KC = ((KC_48_hr -  KC_media_1)/ (KC_48_hr +  KC_media_1)),
         MCP1 = ( (MCP1_48_hr- MCP1_media_1)/ (MCP1_48_hr+ MCP1_media_1)),
         MIP2 = ((MIP2_48_hr - MIP2_media_1)/ (MIP2_48_hr + MIP2_media_1)),
         TNFA = ((TNFA_48_hr- TNFA_media_1)/ (TNFA_48_hr+TNFA_media_1)))%>%
  unite("all_stat" , rad_stat, amifostine_stat, sep = "_", remove = FALSE) #we DO need to make this variable again 

#pivot the dataset longer for plotting
alltg_DoS_lng<-alltg_DoS%>%
  #unite("stat" , rad_stat, amifostine_stat, sep = "_", remove = FALSE)%>%
  pivot_longer( cols = c(IL6, KC, MCP1,MIP2 ), 
                names_to = "cytokine", values_to = "DoS_values")

###
#plots
levelx_order <- c("Sham_No_RP", "Sham_amifostine", "10_Gy_No_RP", "10_Gy_amifostine")
group_order<- c("vehicle", "amifostine" )
stat_comp<-list(c("sham_vehicle", "sham_amifostine" ), c("sham_vehicle","10_Gy_eq_vehicle"), c("sham_amifostine", "10_Gy_eq_amifostine" ), c("10_Gy_eq_vehicle", "10_Gy_eq_amifostine")  )
#

MCT8_lum_DoS<- ggplot(alltg_DoS_lng,aes(x=all_stat, y=DoS_values, fill= rad_stat, color = amifostine_stat ))+
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
  facet_wrap(~cytokine)+
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
  
  labs(title="", subtitle = "", x="", y= "Cytokine per cell DoS ", caption = "MCT8 2D luminex, 48hr normalized to pre amifostine/radiation, per cell")

MCT8_lum_DoS





#STATS for DoS
library(rstatix)
library(FSA)
#stats for cell number by condition
count_norm_test<-hr48_abs_lng_IL6%>%
  group_by(amifostine_stat, rad_stat) %>%
  shapiro_test(est_count)
#normal dist
aov.count<- hr48_abs_lng_IL6 %>% #ANOVA
  anova_test(est_count ~ all_stat )
get_anova_table(aov.count)

#norm test on DoS
all_test_norm<-alltg_DoS%>%
  group_by(amifostine_stat, rad_stat) %>%
  shapiro_test(IL6, 
               KC, 
               MCP1,
               MIP2 
  )
all_test_norm<-all_test_norm%>%
  arrange(variable)
#MIP2 normal, all others not normal   Note that TNFA has too few values to test




#normal dist

aov.MIP2<- alltg_DoS %>% #ANOVA
  anova_test(MIP2 ~ all_stat )
get_anova_table(aov.MIP2)
#NS

MIP2_tuk<-alltg_DoS %>% #pairwise comparison
  group_by() %>%
  tukey_hsd(MIP2 ~ all_stat )
MIP2_tuk<-MIP2_tuk%>% mutate(cytokine = "MIP2")



###
#non-normal dist
kruskal.test(IL6 ~ all_stat, data = alltg_DoS) 
#NS
#Pairwise comparison
IL6_pwc<-dunnTest(data=alltg_DoS, IL6 ~ all_stat)
get_anova_table(IL6_pwc)

kruskal.test(KC ~ all_stat, data = alltg_DoS) 
#NS
#Pairwise comparison
KC_pwc<-dunnTest(data=alltg_DoS, KC ~ all_stat)
get_anova_table(KC_pwc)


kruskal.test(MCP1 ~ all_stat, data = alltg_DoS)
#significant
#Pairwise comparison
MCP1_pwc<-dunnTest(data=alltg_DoS, MCP1 ~ all_stat)
get_anova_table(MCP1_pwc)

kruskal.test(TNFA ~ all_stat, data = alltg_DoS)
#NS
#Pairwise comparison
MCP1_pwc<-dunnTest(data=alltg_DoS, TNFA ~ all_stat)
get_anova_table(MCP1_pwc)

#STATS for 2D absolute 48 hrs
pg_cell_48hr<-pg_cell_48hr%>%
 unite("all_stat" , rad_stat, amifostine_stat, sep = "_", remove = FALSE)

#normality test
h48_test_norm<-pg_cell_48hr%>%
  group_by(amifostine_stat, rad_stat) %>%
  shapiro_test(IL6_48_hr, 
               KC_48_hr, 
               MCP1_48_hr,
               MIP2_48_hr 
  )
h48_test_norm<-h48_test_norm%>%
  arrange(variable)
#IL6 and MCP not normal

aov.MIP2_48_hr<- pg_cell_48hr %>% #ANOVA
  anova_test(MIP2_48_hr ~ all_stat )
get_anova_table(aov.MIP2_48_hr)
#NS
aov.KC_48_hr<- pg_cell_48hr %>% #ANOVA
  anova_test(KC_48_hr ~ all_stat )
get_anova_table(aov.KC_48_hr)
#NS
#non normal means compare
kruskal.test(IL6_48_hr ~ all_stat, data = pg_cell_48hr)
#NS
kruskal.test(MIP2_48_hr ~ all_stat, data = pg_cell_48hr)
#NS
#
#
#
#
##
###for 3D samples

#create med_2 (pre-treatment) df
med_2_3D<-MCT8.1%>%
  filter(time == "media_2",
         sample_type =="3D.sample")%>%
  pivot_wider(names_from =  time, values_from = c(IL1B,IL6,KC,MCP1,MIP2,TNFA))

#create 48hr df, remove empty well
hr48_3D<-MCT8.1%>%
  filter(time == "48_hr",
         sample_type =="3D.sample", #remove standards
  )%>%
  pivot_wider(names_from =  time, values_from = c(IL1B,IL6,KC,MCP1,MIP2,TNFA))

#remove wells in 48 hr timepoint that will be duplicate in med_2 timepoint
hr48_3D<-hr48_3D%>%
  select(-(1:5),-(7:12))
#change name of column 1 to avoid creating duplicates when merging
colnames(hr48_3D)[1] <- "Exp.well"

#merge med_2 and 48hr dfs
alltg_3D<-cbind(hr48_3D, med_2_3D) #combine the two timepoint dfs
#remove first column of combine dataset to avoid redundancy 
alltg_3D<-alltg_3D%>%
  select(-1)

MCT8_3D_DAPI_count <- read_csv("C:/Users/Neal/OneDrive - UMass Lowell/Pain Research Lab/Grants/NeuroRad_GMHA (2024)/Data/MCT8/MCT8_secreted_markers_luminex/MCT8_concatinated_DAPI_counts_by_slice_1.csv")


alltg_3D<-left_join(alltg_3D, MCT8_3D_DAPI_count, by = c("Exp_well" = "well_ID" ))
alltg_3D[,1:6] <- sapply(alltg_3D[,1:6],as.numeric)#need to make columns numeric
alltg_3D[,19:24] <- sapply(alltg_3D[,19:24],as.numeric)#need to make columns numeric


# multiply times 1 because we are looking at content per cell, not per group of cells
#multiply by 0.2 to convert from pg/ mL to pg per well
#multiply by 1000 to convert to femto grams
alltg_relative_3D<-alltg_3D%>%
  mutate(IL6 = (IL6_48_hr* 0.2 /total_nuclei_per_construct*1000),
         KC = (KC_48_hr *0.2 /  total_nuclei_per_construct*1000),
         MCP1 = ( MCP1_48_hr *0.2/ total_nuclei_per_construct*1000),
         MIP2 = (MIP2_48_hr *0.2 / total_nuclei_per_construct*1000),
         TNFA = (TNFA_48_hr *0.2/ total_nuclei_per_construct*1000))%>%
  unite("all_stat" , rad_stat, amifostine_stat, sep = "_", remove = FALSE)%>%
  mutate(rad_stat = as.factor(rad_stat))


#pivot the dataset longer for plotting
alltg_FC_3D_lng<-alltg_relative_3D%>%
  #unite("stat" , rad_stat, amifostine_stat, sep = "_", remove = FALSE)%>%
  pivot_longer( cols = c(IL6, KC, MCP1,MIP2 ), 
                names_to = "cytokine", values_to = "rel_values")


#
###
####
###3D plots
levelx_order4 <- c("sham", "10_Gy")


MCT8_lum_FC_3D<- ggplot(alltg_FC_3D_lng,aes(x=rad_stat, y=rel_values, fill= rad_stat ))+
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
  
  scale_x_discrete(limits = levelx_order4)+
  #scale_color_manual(values= c( "red", "black"))+
  scale_fill_manual(values= c("#999999", "#ecececec"))+
  #ylim( c(-1.01, 1.01))+
  facet_wrap(~cytokine, scales =  "free")+
  geom_rect(aes(xmin = 1.5, xmax = Inf, ymin = -Inf, ymax = Inf),
            alpha = .05,
            fill = "grey", color= 0)+
  
  #stat_pvalue_manual(MFR_dunn, label = "p.adj.signif",  y.position = c(29, 35))+
  #The above line is not working right now and I suspect that this is because fill/ color are different metrics from all_stat
  #this could be really useful though!
  
  #stat_compare_means(method = "wilcox.test",paired = FALSE
  #                  ,comparisons = stat_comp,
  #                 label= "p.signif", hide.ns = TRUE,
  #                symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")))+
  
  labs(title="", subtitle = "", x="", y= "Cytokine per cell (fg)", caption = "MCT8 3D luminex, 48hr normalized estimated cell number")

MCT8_lum_FC_3D

#for export select out the media_2 data that we do not want to include
alltg_relative_3D<-alltg_relative_3D%>%
  select(-(1:11), -(17:24))

#export data for supplemental table

proj<-"MCT8"
wd<-getwd()
filename_csv<-paste(proj,"_Luminex_supplemental_3D",".csv", sep = "") #concatenate file name + "_summary" ".csv"
filename_csv<-file.path(wd,filename_csv) #make file path
write.csv(alltg_relative_3D,filename_csv)#save to WD as .csv




#stats

all_test_norm_3D<-alltg_relative_3D%>%
  group_by(rad_stat) %>%
  shapiro_test(IL6, 
               KC, 
               MCP1,
               MIP2 
  )
all_test_norm_3D<-all_test_norm_3D%>%
  arrange(variable)
#all are normal

#ttest
t.IL6<- alltg_relative_3D %>% 
    t_test(IL6 ~ rad_stat)
t.KC<- alltg_relative_3D %>%  
    t_test(KC ~ rad_stat)
t.MCP1<- alltg_relative_3D %>%
    t_test(MCP1 ~ rad_stat)
t.MIP2<- alltg_relative_3D %>%  
    t_test(MIP2 ~ rad_stat)
t.TNFA<- alltg_relative_3D %>%  
    t_test(TNFA ~ rad_stat)

compile_ttest_3D<-rbind(t.IL6, t.KC, t.MCP1, t.MIP2, t.TNFA)
#NS for all

