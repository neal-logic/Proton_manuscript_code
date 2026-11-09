library(tidyverse)
library(ggpubr)

setwd("C:/Users/Neal/OneDrive - UMass Lowell/Pain Research Lab/Grants/NeuroRad_GMHA (2024)/Data/MCT8/MCT8_DAPI_NeuN_count/results")


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
newcols2<- c("slice", "G")
newcols3<- c("project", "objective", "junk", "treatment", "sample", "image")
MCT8_NeuN_count<- mmm %>%
  separate(Label, into= newcols, sep = ":", remove=TRUE  )%>%
  separate(E, into= newcols2, sep = "/", remove=TRUE  )%>%
  separate(G, into= newcols3, sep = "_", remove=TRUE  )
#second separation to acquire plane


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


###Plot
percent_NeuN<-percent_NeuN%>%
  mutate(treatment =str_replace_all(treatment, c("sham" = "Sham", 
                         "10Gyeq" = "10 Gy")))


x_order<- c( "Sham", "10 Gy")
stat_comp<- list(c( "Sham", "10 Gy"))

NeuN_plot<- ggplot(percent_NeuN, aes(x= treatment , y= percent_NeuN_positive_cells))+
  geom_boxplot(aes(fill = treatment), show.legend = FALSE)+
  #geom_jitter(width = 0.1)+

theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
      axis.line = element_line(colour = "black"))+
  theme(axis.text.x = element_text(angle = 0, vjust = 0.0, hjust=0.5, size = 16, colour = "black"))+
  theme(axis.text.y = element_text(angle = 00, vjust = 0.5, hjust=0, size = 18, colour = "black"))+
  theme(plot.title = element_text(hjust = 0.5))+              # title posisiton
  theme(plot.subtitle = element_text(hjust = 0.5)) +
  theme(axis.title.y = element_text(size = 14))+
  #theme( axis.text.x=element_blank())+
  scale_x_discrete(limits = x_order)+
  ylim(0, 100)+
  scale_fill_manual(values = c(   "#999999", "#cecece"))+
  
  stat_compare_means(method = "t.test",paired = FALSE,
                     comparisons = stat_comp,
                    label.y =  75,
                     label= "p.signif", hide.ns = FALSE,
                     symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")))+

  labs(title="", subtitle = "", x="", y= " Percent of NeuN positive nuclei", caption = "MCT12 NeuN nuclei count, per biological replicate")
#theme_bw()



NeuN_plot

#Statistics
library(rstatix)
test_norm<-percent_NeuN%>%
  group_by(treatment) %>%
  shapiro_test(percent_NeuN_positive_cells)
#can't reject normality

percent_NeuN<-percent_NeuN%>%
  mutate(treatment= as.factor(treatment))

NeuN_stats <-t_test(percent_NeuN_positive_cells~ treatment, data = percent_NeuN)
#p = 0.5911
NeuN_stats<-NeuN_stats%>%
  add_y_position()
