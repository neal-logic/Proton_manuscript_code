library(tidyverse)
#MCT8
#20240524
MCT8_TUJ1_SNT_traces <- read_csv("C:/Users/Neal/OneDrive - UMass Lowell/Pain Research Lab/Grants/NeuroRad_GMHA (2024)/Data/MCT8/MCT8_MCT5_ICC_images/MCT8_7day_TUJ1_measure/MCT8_Traces_Data/MCT8_Traces_Data.csv")


MCT8_TUJ1_SNT_trace_sep<-MCT8_TUJ1_SNT_traces %>%
  separate(image,c("project", "timept1","DAPI", "target1", "stain1", "target2", "junka", "objective",  "primary_stat", "IR_status", "sample", "image"), "_", remove = FALSE )%>%#extract plate name
  mutate(IR_status = str_replace(IR_status,  "1", "Sham" ))%>%
  mutate(IR_status = str_replace(IR_status,  "0", "10 Gy-min" ))%>%
  unite("image_name", sample,image,IR_status, sep = "_", remove = FALSE)  %>% #create the all_stat column which combines RP and rad
  mutate(image_name= as.factor(image_name))  

image_summary<- MCT8_TUJ1_SNT_trace_sep%>%
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
  labs(title="", subtitle = "", x="", y= "Summed neurite length (μm) per cell object", caption = "MCT8 48hr post")+
  
  stat_compare_means(method = "t.test",paired = FALSE
                     ,comparisons = stat_comp,
                     label= "p.signif", hide.ns = FALSE,
                     symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")))

neurite_per_cell_plot 
#note that the warning here is the two images that have 0 cells
setwd("C:/Users/Neal/OneDrive - UMass Lowell/Pain Research Lab/Grants/NeuroRad_GMHA (2024)/Data/MCT8/MCT8_MCT5_ICC_images/MCT8_7day_TUJ1_measure")
ggsave("MCT8_TUJ1_length_per_cell_object.png",width = 1500, height = 1500, units = "px" )
#stats
test_norm_MCT8_Sham<-neurite_per_var%>%
  filter(IR_stat=="Sham")# 
shapiro.test(test_norm_MCT8_Sham$neurite_sum_um_per_cell) # normal

test_norm_MCT8_10Gyeq<-neurite_per_var%>%
  filter(IR_stat=="10 Gy-min")# 
shapiro.test(test_norm_MCT8_10Gyeq$neurite_sum_um_per_cell) #normal

testoft<-t.test(neurite_per_var$neurite_sum_um_per_cell~neurite_per_var$IR_stat)
testoft
#results

MCT8_TUJ1_n_summary<-neurite_per_var%>%
  group_by(IR_stat) %>% 
  summarise(mean_len= mean(neurite_length_sum_um, na.rm = TRUE),
            mean_len_sd= sd(neurite_length_sum_um, na.rm = TRUE))
