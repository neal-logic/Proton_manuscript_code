library(tidyverse)
#MCT5
#20240425
MCT5_TUJ1_SNT_traces <- read_csv("C:/Users/Neal/OneDrive - UMass Lowell/Pain Research Lab/Grants/NeuroRad_GMHA (2024)/Data/MCT8/MCT8_MCT5_ICC_images/MCT5_TUJ1_measure/MCT5_TUJ1_SNT_traces.csv")


MCT5_TUJ1_SNT_trace_sep<-MCT5_TUJ1_SNT_traces %>%
  separate(image,c("project", "stain1","target2", "stain2", "target3", "stain3", "IR", "junk", "junk1", "junk2", "objective", "primary_stat", "IR_status", "sample", "image"), "_", remove = FALSE )%>%#extract plate name
  mutate(IR_status = str_replace(IR_status,  "10Gyeq", "10 Gy-min" ))%>%
  mutate(IR_status = str_replace(IR_status,  "sham", "Sham" ))%>%
  unite("image_name", sample,image,IR_status, sep = "_", remove = FALSE)  %>% #create the all_stat column which combines RP and rad
  mutate(image_name= as.factor(image_name))


image_summary<- MCT5_TUJ1_SNT_trace_sep%>%
  group_by(image_name)%>%
  summarise(neurite_length_sum_um= sum(neurite_length),
            image_area_sum=sum(image_area, na.rm = TRUE),
            num_cell_sum = sum(TUJ1_cell_count, na.rm = TRUE))

neurite_per_var<-image_summary%>%
  mutate(image_area_sum_mm= (image_area_sum * 1e-6))%>%
  mutate(neurite_sum_um_divide_area_mm2= (neurite_length_sum_um/image_area_sum_mm) )%>%
  separate(image_name, c("sample", "image", "IR_stat"), sep = "_", remove = FALSE)%>%
  mutate(neurite_sum_um_per_cell = ((neurite_length_sum_um/ num_cell_sum) ))



##################
#########################
library(ggpubr)
levelx_order <- c( "Sham", "10 Gy-min")
stat_comp<-list(c("Sham", "10 Gy-min" ))

neurite_plot<- ggplot(neurite_per_var,aes(x=IR_stat, y=neurite_sum_um_divide_area_mm2))+
  geom_boxplot(aes( fill= IR_stat), outlier.shape = NA,width = 0.2, show.legend = FALSE, notch = FALSE)+
  geom_jitter(aes(),width = 0.05, show.legend = FALSE)+
  
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
  labs(title="", subtitle = "", x="", y= "Summed neurite length (μm) per area (mm^2)", caption = "MCT5 48hr post")+
  
  stat_compare_means(method = "t.test",paired = FALSE
                    ,comparisons = stat_comp,
                   label= "p.signif", hide.ns = FALSE,
                  symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")))
  


  neurite_plot
#######
#stats
test_norm_MCT5_Sham<-neurite_sum_um_per_area_mm%>%
    filter(IR_stat=="Sham")# same number of nuclei as there are for 53BP1
  shapiro.test(test_norm_MCT5_Sham$neurite_sum_um_divide_area_mm2) #normal
  
  test_norm_MCT5_10Gyeq<-neurite_sum_um_per_area_mm%>%
    filter(IR_stat=="10 Gy-min")# same number of nuclei as there are for 53BP1
  shapiro.test(test_norm_MCT5_10Gyeq$neurite_sum_um_divide_area_mm2) #normal
  
testoft<-t.test(neurite_sum_um_per_area_mm$neurite_sum_um_divide_area_mm2~neurite_sum_um_per_area_mm$IR_stat)
#
#results
MCT5_TUJ1_n_summary<-neurite_per_var%>%
  group_by(IR_stat) %>% 
  summarise(mean_len= mean(neurite_length_sum_um, na.rm = TRUE),
            mean_len_sd= sd(neurite_length_sum_um, na.rm = TRUE))

#
#
#
#
#########################
### cell count
neurite_per_cell_plot<- ggplot(neurite_per_var,aes(x=IR_stat, y=neurite_sum_um_per_cell))+
  geom_boxplot(aes( fill= IR_stat), outlier.shape = NA,width = 0.2, show.legend = FALSE, notch = FALSE)+
  geom_jitter(aes(),width = 0.05, show.legend = FALSE)+
  
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
  labs(title="", subtitle = "", x="", y= "Summed neurite length (μm) per cell object", caption = "MCT5 48hr post")+
  
  stat_compare_means(method = "t.test",paired = FALSE
                     ,comparisons = stat_comp,
                     label= "p.signif", hide.ns = FALSE,
                     symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")))

neurite_per_cell_plot

#stats
test_norm_MCT5_Sham<-neurite_per_var%>%
  filter(IR_stat=="Sham")# same number of nuclei as there are for 53BP1
shapiro.test(test_norm_MCT5_Sham$neurite_sum_um_per_cell) # not normal

test_norm_MCT5_10Gyeq<-neurite_per_var%>%
  filter(IR_stat=="10 Gy-min")# same number of nuclei as there are for 53BP1
shapiro.test(test_norm_MCT5_10Gyeq$neurite_sum_um_per_cell) #normal

testoft<-wilcox.test(neurite_per_var$neurite_sum_um_per_cell~neurite_per_var$IR_stat)
testoft
