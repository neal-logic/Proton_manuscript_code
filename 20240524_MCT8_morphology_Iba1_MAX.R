library(tidyverse)
#MCT8
#20240524
#Iba1 AVG projections

setwd("C:/Users/Neal/OneDrive - UMass Lowell/Pain Research Lab/Grants/NeuroRad_GMHA (2024)/Data/MCT8/MCT8_MCT5_ICC_images/MCT8_Iba1_measure/IBA1_results_AVG/")
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
  mutate(`files[i]`= str_sub(`files[i]`,42,81))

newcols<-c("Project", "objective", "ab_stat", "treatment", "sample", "image")
data_rbind<-data_rbind %>%
  separate(`files[i]`, into= newcols, sep = "_", remove=FALSE  )

MCT8_Iba1<-data_rbind%>%
  mutate(treatment = str_replace(treatment,  "sham", "Sham" ))%>%
  mutate(treatment = str_replace(treatment,  "10Gyeq", "10 Gy-min" ))%>%
  mutate(treatment = as.factor(treatment))%>%
  mutate(Mean = as.numeric(Mean))%>%
  filter(Area<=15000) #must filter the single giant object bc it is an outlier in size



#MCT8 Iba1 Mean
##################
#########################
library(ggpubr)
levelx_order <- c( "Sham", "10 Gy-min")
stat_comp<-list(c("Sham", "10 Gy-min" ))

MCT8_Iba1_Mean<- ggplot(MCT8_Iba1,aes(x=treatment, y=Mean))+
  geom_violin(aes(fill= treatment), show.legend = FALSE)+
  geom_boxplot(aes( ), outlier.shape = NA,width = 0.2, show.legend = FALSE, notch = FALSE)+
  #geom_jitter(aes(),width = 0.075, alpha = 0.3, show.legend = FALSE)+
  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
        axis.line = element_line(colour = "black"))+
  theme(axis.text.x = element_text(angle = -0, vjust = 0.5, hjust=0.5, size = 18, colour = "black"))+
  #theme(axis.text.x =  element_blank())+
  theme(axis.text.y = element_text(angle = 00, vjust = 0.5, hjust=0, size = 18, colour = "black"))+
  theme(plot.title = element_text(hjust = 0.5))+              # title posisiton
  theme(plot.subtitle = element_text(hjust = 0.5)) +
  theme(axis.title.y = element_text(size = 16))+
  
  scale_x_discrete(limits = levelx_order)+
  #scale_color_manual(values= c(  "#999999", "#ecececec"))+
  scale_fill_manual(values= c("#999999", "#ecececec"))+
 # ylim( c(-1.01, 1.01))
  labs(title="", subtitle = "", x="", y= "Iba1 mean intensity per cell object", caption = "MCT8 7 day post IPR treatment")+
  
  stat_compare_means(method = "wilcox.test",paired = FALSE
                    ,comparisons = stat_comp,
                   label= "p.signif", hide.ns = FALSE,
                  symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")))
  
MCT8_Iba1_Mean
ggsave("MCT8_Iba1_mean_intensity_by_AVG_proj.png",width = 1500, height = 1500, units = "px" )
#######
#stats
test_norm_MCT8_Mean_Sham<-MCT8_Iba1%>%
    filter(treatment=="Sham")#
  shapiro.test(test_norm_MCT8_Mean_Sham$Mean) # not normal
  
  test_norm_MCT8_Mean_10Gyeq<-MCT8_Iba1%>%
    filter(treatment=="10 Gy-min")
  shapiro.test(test_norm_MCT8_Mean_10Gyeq$Mean) # not normal
  
testofw<-wilcox.test(data = MCT8_Iba1, Mean~treatment)
testofw #not significant
#
#
##results
MCT8_Iba1_summary<-MCT8_Iba1%>%
  group_by(treatment) %>% 
  summarise(mean_intensity= mean(Mean, na.rm = TRUE),
            mean_intensity_sd= sd(Mean, na.rm = TRUE),
            mean_area= mean(Area, na.rm = TRUE),
            mean_area_sd= sd(Area, na.rm = TRUE),
            mean_circ = mean(Circ., na.rm = TRUE),
            mean_circ_sd = sd(Circ., na.rm = TRUE))

#########################
MCT8_Iba1_Area<- ggplot(MCT8_Iba1,aes(x=treatment, y=Area))+
  geom_violin(aes(fill= treatment), show.legend = FALSE)+
  geom_boxplot(aes( ), outlier.shape = NA,width = 0.2, show.legend = FALSE, notch = FALSE)+
  #geom_jitter(aes(),width = 0.075, alpha = 0.3, show.legend = FALSE)+
  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
        axis.line = element_line(colour = "black"))+
  theme(axis.text.x = element_text(angle = -0, vjust = 0.5, hjust=0.5, size = 18, colour = "black"))+
  #theme(axis.text.x =  element_blank())+
  theme(axis.text.y = element_text(angle = 00, vjust = 0.5, hjust=0, size = 18, colour = "black"))+
  theme(plot.title = element_text(hjust = 0.5))+              # title posisiton
  theme(plot.subtitle = element_text(hjust = 0.5)) +
  theme(axis.title.y = element_text(size = 16))+
  
  scale_x_discrete(limits = levelx_order)+
  #scale_color_manual(values= c(  "#999999", "#ecececec"))+
  scale_fill_manual(values= c("#999999", "#ecececec"))+
  # ylim( c(-1.01, 1.01))
  labs(title="", subtitle = "", x="", y= "Iba1 area (μm)", caption = "MCT8 7 day post IPR treatment")+
  
  stat_compare_means(method = "wilcox.test",paired = FALSE
                     ,comparisons = stat_comp,
                     label= "p.signif", hide.ns = FALSE,
                     symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")))

MCT8_Iba1_Area
ggsave("MCT8_Iba1_Area.png",width = 1500, height = 1500, units = "px" )
#######
#stats
test_norm_MCT8_Area_Sham<-MCT8_Iba1%>%
  filter(treatment=="Sham")#
shapiro.test(test_norm_MCT8_Area_Sham$Area) # not normal

test_norm_MCT8_Area_10Gyeq<-MCT8_Iba1%>%
  filter(treatment=="10 Gy equiv.")# 
shapiro.test(test_norm_MCT8_Area_10Gyeq$Area) # not normal

testofw_a<-wilcox.test(data = MCT8_Iba1, Area~treatment)
testofw_a #not significant
#
#########################
#Circularity
#########################
MCT8_Iba1_Circ.<- ggplot(MCT8_Iba1,aes(x=treatment, y=Circ.))+
  geom_violin(aes(fill= treatment), show.legend = FALSE)+
  geom_boxplot(aes( ), outlier.shape = NA,width = 0.2, show.legend = FALSE, notch = FALSE)+
  #geom_jitter(aes(),width = 0.075, alpha = 0.3, show.legend = FALSE)+
  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
        axis.line = element_line(colour = "black"))+
  theme(axis.text.x = element_text(angle = -0, vjust = 0.5, hjust=0.5, size = 18, colour = "black"))+
  #theme(axis.text.x =  element_blank())+
  theme(axis.text.y = element_text(angle = 00, vjust = 0.5, hjust=0, size = 18, colour = "black"))+
  theme(plot.title = element_text(hjust = 0.5))+              # title posisiton
  theme(plot.subtitle = element_text(hjust = 0.5)) +
  theme(axis.title.y = element_text(size = 16))+
  
  scale_x_discrete(limits = levelx_order)+
  #scale_color_manual(values= c(  "#999999", "#ecececec"))+
  scale_fill_manual(values= c("#999999", "#ecececec"))+
  # ylim( c(-1.01, 1.01))
  labs(title="", subtitle = "", x="", y= "Iba1 circularity", caption = "MCT8 7 day post IPR treatment")+
  
  stat_compare_means(method = "wilcox.test",paired = FALSE
                     ,comparisons = stat_comp,
                     label= "p.signif", hide.ns = FALSE,
                     symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")))

MCT8_Iba1_Circ.
ggsave("MCT8_Iba1_Circularity.png",width = 1500, height = 1500, units = "px" )
#######
#stats
test_norm_MCT8_Circ._Sham<-MCT8_Iba1%>%
  filter(treatment=="Sham")# 
shapiro.test(test_norm_MCT8_Circ._Sham$Circ.) # not normal

test_norm_MCT8_Circ._10Gyeq<-MCT8_Iba1%>%
  filter(treatment=="10 Gy equiv.")#
shapiro.test(test_norm_MCT8_Circ._10Gyeq$Circ.) # not normal

testofw_cir<-wilcox.test(data = MCT8_Iba1, Circ.~treatment)
testofw_cir #not significant
#
#########################

count_MCT8_Iba1<-MCT8_Iba1%>%
  group_by(treatment, sample, image)%>%
  summarise(cell_count = n())#this counts the number of cells per ROI by condition

############### count
MCT8_Iba1_count_plot<- ggplot(count_MCT8_Iba1,aes(x=treatment, y=cell_count))+
  geom_boxplot(aes( fill= treatment), outlier.shape = NA,width = 0.2, show.legend = FALSE, notch = FALSE)+
  geom_jitter(aes(),width = 0.075, alpha = 0.3, show.legend = FALSE)+
  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
        axis.line = element_line(colour = "black"))+
  theme(axis.text.x = element_text(angle = -0, vjust = 0.5, hjust=0.5, size = 18, colour = "black"))+
  #theme(axis.text.x =  element_blank())+
  theme(axis.text.y = element_text(angle = 00, vjust = 0.5, hjust=0, size = 18, colour = "black"))+
  theme(plot.title = element_text(hjust = 0.5))+              # title posisiton
  theme(plot.subtitle = element_text(hjust = 0.5)) +
  theme(axis.title.y = element_text(size = 16))+
  
  scale_x_discrete(limits = levelx_order)+
  #scale_color_manual(values= c(  "#999999", "#ecececec"))+
  scale_fill_manual(values= c("#999999", "#ecececec"))+
  # ylim( c(-1.01, 1.01))
  labs(title="", subtitle = "", x="", y= "Iba1 cell object per ROI count", caption = "MCT8 7 day post IPR treatment")+
  
  stat_compare_means(method = "wilcox.test",paired = FALSE
                     ,comparisons = stat_comp,
                     label= "p.signif", hide.ns = FALSE,
                     symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")))

MCT8_Iba1_count_plot
ggsave("MCT8_Iba1_count_plot.png",width = 1500, height = 1500, units = "px" )

#######
#stats
test_norm_MCT8_count_Sham<-count_MCT8_Iba1%>%
  filter(treatment=="Sham")# 
shapiro.test(test_norm_MCT8_Circ._Sham$Circ.) # not normal

test_norm_MCT8_count_10Gyeq<-count_MCT8_Iba1%>%
  filter(treatment=="10 Gy equiv.")#
shapiro.test(test_norm_MCT8_Circ._10Gyeq$Circ.) # not normal

testofw_count<-wilcox.test(data = count_MCT8_Iba1, cell_count~treatment)
testofw_count #not significant

