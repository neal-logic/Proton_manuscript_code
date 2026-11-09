library(tidyverse)

MCT8 <- read_csv("C:/Users/Neal/OneDrive - UMass Lowell/Pain Research Lab/Grants/NeuroRad_GMHA (2024)/Data/MCT8/MCT8_secreted_markers_luminex/20240703_processed_Luminex_data_MCT8.csv")
MCT8_20240710 <- read_csv("C:/Users/Neal/OneDrive - UMass Lowell/Pain Research Lab/Grants/NeuroRad_GMHA (2024)/Data/MCT8/MCT8_secreted_markers_luminex/20240710_luminex_vol_corrected/20240710_MCT8_luminex_3D.csv")

MCT8<-rbind(MCT8,MCT8_20240710)

MCT8.1<-MCT8%>%
  mutate(IL1B=as.numeric(IL1B))%>%
  mutate(IL6=as.numeric(IL6))%>%
  mutate(KC=as.numeric(KC))%>%
  mutate(MCP1=as.numeric(MCP1))%>%
  mutate(MIP2=as.numeric(MIP2))%>%
  mutate(TNFA=as.numeric(TNFA))



MCT8_rel_wide<-MCT8.1%>%
  filter(sample_type == "sample",
  time != "media_1")%>%
  pivot_wider(names_from =  time, values_from = c(IL1B,IL6,KC,MCP1,MIP2,TNFA))

#create med_2 (pre-treatment) df
med_2<-MCT8.1%>%
  filter(time == "media_2",
         sample_type =="sample")%>%
  pivot_wider(names_from =  time, values_from = c(IL1B,IL6,KC,MCP1,MIP2,TNFA))

#create 48hr df, remove empty well
hr48<-MCT8.1%>%
  filter(time == "48_hr",
         sample_type =="sample", #remove standards
         Exp_well != "F4")%>%
  pivot_wider(names_from =  time, values_from = c(IL1B,IL6,KC,MCP1,MIP2,TNFA))

#remove wells in 48 hr timepoint that will be duplicate in med_2 timepoint
hr48<-hr48%>%
  select(-(1:5),-(7:10))
#change name of column 1 to avoid creating duplicates when merging
colnames(hr48)[1] <- "Exp.well"

#merge med_2 and 48hr dfs
alltg<-cbind(hr48, med_2) #combine the two timepoint dfs
#remove first column of combine dataset to avoid redundancy 
alltg<-alltg%>%
  select(-1)

#create DoS ((post-pre)/(post+pre))
alltg_DoS<-alltg%>%
  mutate(IL6 = (IL6_48_hr- IL6_media_2) /(IL6_48_hr+ IL6_media_2),
         KC = (KC_48_hr- KC_media_2) / (KC_48_hr+ KC_media_2),
         MCP1 = ( MCP1_48_hr- MCP1_media_2)/ (MCP1_48_hr+ MCP1_media_2),
         MIP2 = (MIP2_48_hr- MIP2_media_2) / (MIP2_48_hr+ MIP2_media_2),
         TNFA = (TNFA_48_hr- TNFA_media_2 )/ (TNFA_48_hr+ TNFA_media_2))%>%
  unite("all_stat" , rad_stat, amifostine_stat, sep = "_", remove = FALSE)
  
#pivot the dataset longer for plotting
alltg_DoS_lng<-alltg_DoS%>%
  #unite("stat" , rad_stat, amifostine_stat, sep = "_", remove = FALSE)%>%
  pivot_longer( cols = c(IL6, KC, MCP1,MIP2, TNFA ), 
                names_to = "cytokine", values_to = "FC_values")

###
#plots
levelx_order <- c("Sham_No_RP", "Sham_amifostine", "10_Gy_No_RP", "10_Gy_amifostine")
group_order<- c("vehicle", "amifostine" )
stat_comp<-list(c("sham_vehicle", "sham_amifostine" ), c("sham_vehicle","10_Gy_eq_vehicle"), c("sham_amifostine", "10_Gy_eq_amifostine" ), c("10_Gy_eq_vehicle", "10_Gy_eq_amifostine")  )
#
  
MCT8_lum_DoS<- ggplot(alltg_DoS_lng,aes(x=all_stat, y=FC_values, fill= rad_stat, color = amifostine_stat ))+
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
  
  labs(title="", subtitle = "", x="", y= "Cytokine DoS ", caption = "MCT8 2D luminex, 48hr normalized to post amifostine pre radiation")

MCT8_lum_DoS

#stats
library(rstatix)
library(FSA)

all_test_norm<-alltg_DoS%>%
  group_by(amifostine_stat, rad_stat) %>%
  shapiro_test(IL6, 
               KC, 
               MCP1,
               MIP2 
            )
all_test_norm<-all_test_norm%>%
  arrange(variable)
#IL6 and MCP1 not-normal,   Note that TNFA has too few values to test



#normal dist
aov.KC<- alltg_DoS %>% #ANOVA
  anova_test(KC ~ all_stat )
get_anova_table(aov.KC)
#NS

KC_tuk<-alltg_DoS %>% #pairwise comparison
  group_by() %>%
  tukey_hsd(KC ~ all_stat )
KC_tuk<-KC_tuk%>% mutate(cytokine = "KC_tuk")

aov.MIP2<- alltg_DoS %>% #ANOVA
  anova_test(MIP2 ~ all_stat )
get_anova_table(aov.MIP2)
#NS

MIP2_tuk<-alltg_DoS %>% #pairwise comparison
  group_by() %>%
  tukey_hsd(MIP2 ~ all_stat )
MIP2_tuk<-MIP2_tuk%>% mutate(cytokine = "MIP2")

aov.TNFA<- alltg_DoS %>% #ANOVA
  anova_test(TNFA ~ all_stat )
get_anova_table(aov.TNFA)
#NS

TNFA_tuk<-alltg_DoS %>% #pairwise comparison
  group_by() %>%
  tukey_hsd(TNFA ~ all_stat )
TNFA_tuk<-TNFA_tuk%>% mutate(cytokine = "TNFA")

#compiled Tukey tests
compile_tukey_tests<-rbind(TNFA_tuk, MIP2_tuk, KC_tuk)

###
#non-normal dist
kruskal.test(IL6 ~ all_stat, data = alltg_DoS) 
#NS
#Pairwise comparison
IL6_pwc<-dunnTest(data=alltg_DoS, IL6 ~ all_stat)
get_anova_table(IL6_pwc)


kruskal.test(MCP1 ~ all_stat, data = alltg_DoS)
#NS
#Pairwise comparison
MCP1_pwc<-dunnTest(data=alltg_DoS, MCP1 ~ all_stat)
get_anova_table(MCP1_pwc)

#compile dunn test
compile_dunn_tests<-rbind(MCP1_pwc, IL6_pwc)

############
##################
#compare raw values

#pivot the dataset for time 1 vs time 2 raw value plotting
levelx_order3 <- c("Sham_No_RP", "Sham_amifostine", "10_Gy_No_RP", "10_Gy_amifostine")

alltg_raw_lng<-alltg%>%
  pivot_longer( cols = c(IL6_48_hr, KC_48_hr, MCP1_48_hr,MIP2_48_hr, TNFA_48_hr,
                         IL6_media_2, KC_media_2, MCP1_media_2, MIP2_media_2, TNFA_media_2), 
                names_to = "cytokine", values_to = "values")%>%
  separate(cytokine,c("cytokine", "media","time_pt"), "_", remove = TRUE )%>%#extract cytokine and timepoint
  mutate(time_pt=str_replace(time_pt,"hr","48hr_post_treatment"))%>%
  mutate(time_pt=str_replace(time_pt,"2","pre_radiation"))%>%
  unite("rad_ami" , rad_stat, amifostine_stat,  sep = "_", remove = FALSE)

#make pretimepoint df
alltg_raw_lng_pre<- alltg_raw_lng%>%
  filter(time_pt == "pre_radiation" )
#pre timpoint 
MCT8_lum_raw_pre<- ggplot(alltg_raw_lng_pre,aes(x=rad_ami, y=values, fill= rad_stat, color = amifostine_stat ))+
  geom_boxplot(aes(), outlier.shape = NA, show.legend = TRUE)+
  geom_jitter(width = 0.1, size= 0.75, show.legend = FALSE)+
  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
        axis.line = element_line(colour = "black"))+
  theme(axis.text.x = element_text(angle = -45, vjust = 0.5, hjust=0.0, size = 12, colour = "black"))+
  theme(axis.text.x =  element_blank())+
  theme(axis.text.y = element_text(angle = 00, vjust = 0.5, hjust=0, size = 12, colour = "black"))+
  theme(plot.title = element_text(hjust = 0.5))+              # title posisiton
  theme(plot.subtitle = element_text(hjust = 0.5)) +
  theme(axis.title.y = element_text(size = 14))+
  
  scale_x_discrete(limits = levelx_order3)+
  scale_color_manual(values= c( "red", "black"))+
  scale_fill_manual(values= c("#999999", "#ecececec"))+
  #ylim( c(-1.01, 1.01))+
  facet_wrap(~cytokine, scales = "free")+
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
  
  labs(title="", subtitle = "", x="", y= "Concentration (pg/mL)", caption = "MCT8 2D luminex, pre radiation, post 1hr in amifostine")

MCT8_lum_raw_pre
#############################

alltg_raw_lng_post<- alltg_raw_lng%>%
  filter(time_pt == "48hr_post_treatment" )
#pre timpoint 
MCT8_lum_raw_post<- ggplot(alltg_raw_lng_post,aes(x=rad_ami, y=values, fill= rad_stat, color = amifostine_stat ))+
  geom_boxplot(aes(), outlier.shape = NA, show.legend = TRUE)+
  geom_jitter(width = 0.1, size= 0.75, show.legend = FALSE)+
  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
        axis.line = element_line(colour = "black"))+
  theme(axis.text.x = element_text(angle = -45, vjust = 0.5, hjust=0.0, size = 12, colour = "black"))+
  theme(axis.text.x =  element_blank())+
  theme(axis.text.y = element_text(angle = 00, vjust = 0.5, hjust=0, size = 12, colour = "black"))+
  theme(plot.title = element_text(hjust = 0.5))+              # title posisiton
  theme(plot.subtitle = element_text(hjust = 0.5)) +
  theme(axis.title.y = element_text(size = 14))+
  
  scale_x_discrete(limits = levelx_order3)+
  scale_color_manual(values= c( "red", "black"))+
  scale_fill_manual(values= c("#999999", "#ecececec"))+
  #ylim( c(-1.01, 1.01))+
  facet_wrap(~cytokine, scales = "free")+
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
  
  labs(title="", subtitle = "", x="", y= "Concentration (pg/mL)", caption = "MCT8 2D luminex, post radiation 48hr")

MCT8_lum_raw_post




############
#################
#########################
#compare between med_1 and med_2 timepoint

MCT8_rel_wide_1v2<-MCT8.1%>%
  filter(sample_type == "sample",
         time != "48hr")%>%
  pivot_wider(names_from =  time, values_from = c(IL1B,IL6,KC,MCP1,MIP2,TNFA))

#create med_2 (pre-treatment) df
med_2<-MCT8.1%>%
  filter(time == "media_2",
         sample_type =="sample",
         amifostine_stat == "amifostine")%>%
  pivot_wider(names_from =  time, values_from = c(IL1B,IL6,KC,MCP1,MIP2,TNFA))

#create 48hr df, remove empty well
med_1<-MCT8.1%>%
  filter(time == "media_1",
         sample_type =="sample" )%>% #remove standards
  pivot_wider(names_from =  time, values_from = c(IL1B,IL6,KC,MCP1,MIP2,TNFA))

med_1<-med_1%>%
  select(-(1:5),-(7:10))
#change name of column 1 to avoid creating duplicates when merging
colnames(med_1)[1] <- "Exp.well"

#merge med_2 and 48hr dfs
alltg_1v2<-cbind(med_1, med_2) #combine the two timepoint dfs
#remove first column of combine dataset to avoid redundancy 
alltg_1v2<-alltg_1v2%>%
  select(-1)

#pivot the dataset for time 1 vs time 2 raw value plotting
alltg_1v2_lng<-alltg_1v2%>%
  pivot_longer( cols = c(IL6_media_1, KC_media_1, MCP1_media_1,MIP2_media_1, TNFA_media_1,
                         IL6_media_2, KC_media_2, MCP1_media_2, MIP2_media_2, TNFA_media_2), 
                names_to = "cytokine", values_to = "values")%>%
  separate(cytokine,c("cytokine", "media","time_pt"), "_", remove = TRUE )%>%#extract cytokine and timepoint
  mutate(time_pt=str_replace(time_pt,"1","pre_amifostine"))%>%
  mutate(time_pt=str_replace(time_pt,"2","1hr_post_amifostine"))


#time 1 vs time 2 plots
levelx_order1 <- c("pre_amifostine", "1hr_post_amifostine")


MCT8_lum_1v2<- ggplot(alltg_1v2_lng,aes(x=time_pt, y=values, color = time_pt ))+
  geom_boxplot(aes(), outlier.shape = NA)+
  geom_jitter(width = 0.1, size= 0.75, show.legend = FALSE)+
  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
        axis.line = element_line(colour = "black"))+
  theme(axis.text.x = element_text(angle = -45, vjust = 0.5, hjust=0.0, size = 12, colour = "black"))+
  theme(axis.text.x =  element_blank())+
  theme(axis.text.y = element_text(angle = 00, vjust = 0.5, hjust=0, size = 12, colour = "black"))+
  theme(plot.title = element_text(hjust = 0.5))+              # title posisiton
  theme(plot.subtitle = element_text(hjust = 0.5)) +
  theme(axis.title.y = element_text(size = 14))+
  
  scale_x_discrete(limits = levelx_order1)+
  scale_color_manual(values= c( "red", "purple"))+
  #scale_fill_manual(values= c("#999999", "#ecececec"))+
  #ylim( c(-1.01, 1.01))+
  facet_wrap(~cytokine, scales = "free")+
  
  #stat_pvalue_manual(MFR_dunn, label = "p.adj.signif",  y.position = c(29, 35))+
  #The above line is not working right now and I suspect that this is because fill/ color are different metrics from all_stat
  #this could be really useful though!
  
  #stat_compare_means(method = "wilcox.test",paired = FALSE
  #                  ,comparisons = stat_comp,
  #                 label= "p.signif", hide.ns = TRUE,
  #                symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")))+
  
  labs(title="", subtitle = "", x="", y= "Concentration (pg/mL)", caption = "MCT8 2D luminex, pre vs. post amifostine, prior to radiation")

MCT8_lum_1v2







MCT8_med_2_to_48hr_fc<-MCT8_rel_wide%>%
  mutate(fc_IL6= MCT_rel_wide$IL1B_48_hr/ MCT_rel_wide$IL6_media_2)

MCT8_tidy<-MCT8.1%>%
    pivot_longer( names_to = "cytokine", values_to = "cyto_conc")
  
MCT8_tidy_stan<- MCT8.1%>%
  filter(sample_type == "standard")%>%
  pivot_longer( cols = c(IL6, KC, MCP1,MIP2, TNFA ), 
                names_to = "cytokine", values_to = "values")
  
###################
###########################
######################################

#plot standards
stan_plot<- ggplot(MCT8_tidy_stan,aes(x=`Sample ID`, y= values, color =cytokine))+
  geom_point(aes(), size= 1.5, show.legend = FALSE)+

  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
        axis.line = element_line(colour = "black"))+
  theme(axis.text.x = element_text(angle = -90, vjust = 0.0, hjust=0.5, size = 12, colour = "black"))+
  theme(axis.text.x =  element_blank())+
  theme(axis.text.y = element_text(angle = 00, vjust = 0.5, hjust=0, size = 12, colour = "black"))+
  theme(plot.title = element_text(hjust = 0.5))+              # title posisiton
  theme(plot.subtitle = element_text(hjust = 0.5)) +
  theme(axis.title.y = element_text(size = 14))+
  
  facet_wrap(~cytokine  )+
  scale_y_log10()+
  #scale_colour_manual(values = c("#a6a6a6ff", "#d9d9d9ff"))+
 
  ylab( "Concentration (pg/mL)")+
  xlab( "")

#theme_classic()

stan_plot

#################
#########################
###################################
#compile results 

#DoS results
MCT8_2D_luminex_results_summary_DoS<- alltg_DoS_lng%>%
  group_by(amifostine_stat,rad_stat, cytokine)%>%
  summarise(mean_DoS= mean(FC_values, na.rm = TRUE),
            DoS_sd= sd(FC_values, na.rm = TRUE))
           
#to save as CSV, remember to set the WD
wd<-getwd()
project<- "MCT8_2D_luminex_results"
filename_csv<-paste(project,"_DoS_supplemental",".csv", sep = "") #concatenate file name + "_summary" ".csv"
filename_csv<-file.path(wd,filename_csv) #make file path
write.csv(MCT8_2D_luminex_results_summary_DoS,filename_csv)#save to WD as .csv

MCT8_2D_luminex_results_summary_raw<- alltg_raw_lng%>%
  group_by(amifostine_stat,rad_stat,time_pt, cytokine)%>%
  summarise(mean_DoS= mean(values, na.rm = TRUE),
            DoS_sd= sd(values, na.rm = TRUE))
#write csv
filename_csv<-paste(project,"_absolute_supplemental",".csv", sep = "") #concatenate file name + "_summary" ".csv"
filename_csv<-file.path(wd,filename_csv) #make file path
write.csv(MCT8_2D_luminex_results_summary_raw,filename_csv)#save to WD as .csv


