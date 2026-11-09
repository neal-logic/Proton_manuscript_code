library(tidyverse)
#MCT5
#Iba1
#20240501


setwd("C:/Users/Neal/OneDrive - UMass Lowell/Pain Research Lab/Grants/NeuroRad_GMHA (2024)/Data/MCT8/MCT8_MCT5_ICC_images/MCT5_iba1_measure/AVG_results/")
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
  mutate(`files[i]`= str_sub(`files[i]`,55,83))

newcols<-c("Project", "objective", "ab_stat", "treatment", "sample", "image")
newcols1<-c("junk", "image")
data_rbind<-data_rbind %>%
  separate(`files[i]`, into= newcols, sep = "_", remove=FALSE  )%>%
  separate(image, into= newcols1, sep = "g", remove=TRUE  )%>%
  mutate(image= str_sub(image,1,2))%>%
  mutate(image=as.numeric(image))%>%
  unite("identifier", treatment,sample, image, sep = "_", remove = FALSE)
  
  
MCT5_Iba1<-data_rbind%>%
  mutate(treatment = str_replace(treatment,  "sham", "Sham" ))%>%
  mutate(treatment = str_replace(treatment,  "10Gyeq", "10 Gy-min" ))%>%
  mutate(treatment= as.factor(treatment))%>%
  filter(Area>=100) #filter the small objects



#MCT5 Iba1 Mean
##################
#########################
library(ggpubr)
levelx_order <- c( "Sham", "10 Gy-min")
stat_comp<-list(c("Sham", "10 Gy-min" ))

MCT5_Iba1_Mean<- ggplot(MCT5_Iba1,aes(x=treatment, y=Mean))+
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
  labs(title="", subtitle = "", x="", y= "Iba1 mean intensity per cell object", caption = "MCT5 48hr post")+
  
  stat_compare_means(method = "wilcox.test",paired = FALSE
                    ,comparisons = stat_comp,
                   label= "p.signif", hide.ns = FALSE,
                  symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")))
  
MCT5_Iba1_Mean

#######
#stats
test_norm_MCT5_Mean_Sham<-MCT5_Iba1%>%
    filter(treatment=="Sham")#
  shapiro.test(test_norm_MCT5_Mean_Sham$Mean) # normal
  
  test_norm_MCT5_Mean_10Gyeq<-MCT5_Iba1%>%
    filter(treatment=="10 Gy-min")# same number of nuclei as there are for 53BP1
  shapiro.test(test_norm_MCT5_Mean_10Gyeq$Mean) # not normal
  
testofw<-wilcox.test(data = MCT5_Iba1, Mean~treatment)
testofw #not significant
#
##
#results
MCT5_Iba1_summary<-MCT5_Iba1%>%
  group_by(treatment) %>% 
  summarise(mean_intensity= mean(Mean, na.rm = TRUE),
            mean_intensity_sd= sd(Mean, na.rm = TRUE),
            mean_area= mean(Area, na.rm = TRUE),
            mean_area_sd= sd(Area, na.rm = TRUE),
            mean_circ = mean(Circ., na.rm = TRUE),
            mean_circ_sd = sd(Circ., na.rm = TRUE))


#########################
MCT5_Iba1_Area<- ggplot(MCT5_Iba1,aes(x=treatment, y=Area))+
  geom_violin(aes(fill= treatment), show.legend = FALSE)+
  geom_boxplot(aes(), outlier.shape = NA,width = 0.2, show.legend = FALSE, notch = FALSE)+
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
  labs(title="", subtitle = "", x="", y= "Iba1 object area (μm)", caption = "MCT5 48hr post")+
  
  stat_compare_means(method = "wilcox.test",paired = FALSE
                     ,comparisons = stat_comp,
                     label= "p.signif", hide.ns = FALSE,
                     symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")))

MCT5_Iba1_Area

#######
#stats
test_norm_MCT5_Area_Sham<-MCT5_Iba1%>%
  filter(treatment=="Sham")# same number of nuclei as there are for 53BP1
shapiro.test(test_norm_MCT5_Area_Sham$Area) # not normal

test_norm_MCT5_Area_10Gyeq<-MCT5_Iba1%>%
  filter(treatment=="10 Gy-min")# same number of nuclei as there are for 53BP1
shapiro.test(test_norm_MCT5_Area_10Gyeq$Area) # not normal

testofw_a<-wilcox.test(data = MCT5_Iba1, Area~treatment)
testofw_a #not significant
#
#########################
#Circularity
#########################
MCT5_Iba1_Circ.<- ggplot(MCT5_Iba1,aes(x=treatment, y=Circ.))+
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
  labs(title="", subtitle = "", x="", y= "Iba1 object circularity", caption = "MCT5 48hr post")+
  
  stat_compare_means(method = "wilcox.test",paired = FALSE
                     ,comparisons = stat_comp,
                     label= "p.signif", hide.ns = FALSE,
                     symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")))

MCT5_Iba1_Circ.

#######
#stats
test_norm_MCT5_Circ._Sham<-MCT5_Iba1%>%
  filter(treatment=="Sham")# same number of nuclei as there are for 53BP1
shapiro.test(test_norm_MCT5_Circ._Sham$Circ.) # not normal

test_norm_MCT5_Circ._10Gyeq<-MCT5_Iba1%>%
  filter(treatment=="10 Gy-min")# same number of nuclei as there are for 53BP1
shapiro.test(test_norm_MCT5_Circ._10Gyeq$Circ.) # not normal

testofw_cir<-wilcox.test(data = MCT5_Iba1, Circ.~treatment)
testofw_cir #not significant
#
#########################
########
#############
count_MCT5_Iba1<-MCT5_Iba1%>%
  group_by(treatment, sample, image)%>%
  summarise(cell_count = n())#this counts the number of cells per ROI by condition

############### count
MCT5_Iba1_count_plot<- ggplot(count_MCT5_Iba1,aes(x=treatment, y=cell_count))+
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
  labs(title="", subtitle = "", x="", y= "Iba1 cell object per ROI count", caption = "MCT5 48hr post")+
  
  stat_compare_means(method = "wilcox.test",paired = FALSE
                     ,comparisons = stat_comp,
                     label= "p.signif", hide.ns = FALSE,
                     symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")))

MCT5_Iba1_count_plot
ggsave("MCT5_Iba1_count_plot.png",width = 400, height = 400, units = "px" )
#######
#stats
test_norm_MCT5_count_Sham<-count_MCT5_Iba1%>%
  filter(treatment=="Sham")# same number of nuclei as there are for 53BP1
shapiro.test(test_norm_MCT5_Circ._Sham$Circ.) # not normal

test_norm_MCT5_count_10Gyeq<-count_MCT5_Iba1%>%
  filter(treatment=="10 Gy equiv.")# same number of nuclei as there are for 53BP1
shapiro.test(test_norm_MCT5_Circ._10Gyeq$Circ.) # not normal

testofw_count<-wilcox.test(data = count_MCT5_Iba1, cell_count~treatment)
testofw_count #not significant

