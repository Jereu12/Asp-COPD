#COPD Asp figures

library(ggtext)
library(ggplot2) 
library(permute)
library(lattice)
library(vegan)
library(RColorBrewer)
library(magrittr)
library(ggpubr)
library(FSA)
library(reshape)
library (plyr)
library(ComplexUpset)
library(UpSetR)
library(ggupset)
library(tidyverse, warn.conflicts = FALSE)
library(pheatmap)
library(table1)
library(flextable)
library(MASS)
library(jtools)
library(colorspace)



#Load data

## Fungi species relative abundance 
Master <- read.csv("~/Desktop/copd_fungi_ra.csv")

# Top 30 fungi species
biome.Species <- read.csv("~/Desktop/air/top/CSV/biome-Species.csv", row.names=1)

##Environment data
env_data <- read.csv("~/Desktop/env_data.csv", row.names=1)

##Clinical data
Clinical_data <- read.csv("~/Desktop/Clinical_data.csv", row.names=1)

## Aspergillus allergen data
Asp_allergen <- read.csv("~/Desktop/Asp_allergen.csv", row.names=1)

## Aspergillus abundance data
Asp_ab <- read.csv("~/Desktop/Asp_ab.csv")


##Forest plot data Exacerbation preceding year
asp_exc_b <- read.csv("~/Desktop/asp_exc_b.csv")

##Forest plot data Exacerbation 1 year follow up
asp_exc1yr <- read.csv("~/Desktop/asp_exc1yr.csv")

##Blank data
blank <- read.csv("~/Desktop/blank.csv")

## Sputum ITS data
sputum_ITS <- read.csv("~/Desktop/sputum_ITS.csv", row.names=1)

#Asp Prevalence summarized data
asp.prev<- read.csv("~/Desktop/asp.prev.csv")

######
###########

#Figure 1 upset plot
sg<-subset(Master, Master$Study=="Singapore")#only singapore data
#create prevalence true and false for Aspergillus spp
aa<-sg[,(4:17)]
at<-aa>0
ind1<-colSums(at)!=0
z<-at[,ind1]
aa.p<-cbind(sg[,c(1:3)],z)

#CONVERT TRUE TO NAME OF FUNGI NA TO EMPTY
linelist_sym_1 <- aa.p %>% 
  mutate(A.fumigatus = ifelse( A.fumigatus== "TRUE", "A.fumigatus", NA), 
         A.steynii = ifelse(A.steynii == "TRUE", "A.steynii", NA),
         A.glaucus = ifelse(A.glaucus == "TRUE", "A.glaucus", NA),
         A.nidulans = ifelse(A.nidulans == "TRUE", "A.nidulans", NA),
         A.niger = ifelse(A.niger == "TRUE", "A.niger", NA),
         A.heteromorphus = ifelse(A.heteromorphus == "TRUE", "A.heteromorphus", NA),
         A.clavatus = ifelse(A.clavatus == "TRUE", "A.clavatus", NA),
         A.aculeatinus = ifelse(A.aculeatinus == "TRUE", "A.aculeatinus", NA),
         A.pseudoglaucus = ifelse(A.pseudoglaucus == "TRUE", "A.pseudoglaucus", NA),
         A.awamori = ifelse(A.awamori == "TRUE", "A.awamori", NA),
         A.welwitschiae = ifelse(A.welwitschiae == "TRUE", "A.welwitschiae", NA))

#create list column 
linelist_sym_1 <- linelist_sym_1 %>% 
  unite(col = "all_asp",
        c(A.fumigatus,A.steynii,A.glaucus, A.nidulans, A.niger, A.heteromorphus,A.clavatus,A.aculeatinus,A.pseudoglaucus,A.awamori,A.welwitschiae), 
        sep = "; ",
        remove = TRUE,
        na.rm = TRUE) %>% 
  mutate(
    # make a copy of all_symptoms column, but of class "list" (which is required to use ggupset() in next step)
    all_asp_list = as.list(strsplit(all_asp, "; "))
  )

tiff("figure1.upset.tiff", units="in", width=8, height=6, res=400)

f1<-ggplot(
  data = linelist_sym_1,
  mapping = aes(x = all_asp_list)) +
  geom_bar(aes(fill=aa.p$SampleType),color="black") +
  scale_x_upset(order_by = "degree",
                reverse = FALSE,
                n_intersections = 3,#number of bar
                sets = c("A.fumigatus","A.steynii","A.glaucus", "A.nidulans", "A.niger", "A.heteromorphus","A.clavatus","A.aculeatinus","A.pseudoglaucus","A.awamori","A.welwitschiae"))+
  labs(
    x = "",
    y = "Frequency of homes \nwith detectable Aspergillus spp.")+scale_fill_manual(values = c("purple","orange","green3","green3","green3"))+labs(fill="")+theme(panel.background = element_blank()) +theme(axis.title.x = element_text(face="italic"))+theme_combmatrix(
      combmatrix.panel.point.color.empty = "darkgrey",
      combmatrix.panel.striped_background.color.one   = "white",
      combmatrix.panel.striped_background.color.two   = "grey90",combmatrix.label.text = element_text(face = "italic")
    )
f1+labs(y = ~ atop(paste("Frequency of homes with detectable"),paste('',italic("Aspergillus spp."))))
dev.off()

####Figure 2 

f2<-cbind(env_data,Master)
copd_bf<-subset(f2,f2$Study=="Singapore")

p2<-ggscatter(copd_bf, y = "A.fumigatus", x = "pm25", color = "SampleType",
              add = "reg.line", conf.int = TRUE, 
              cor.coef = TRUE, cor.method = "spearman",xlab = "PM 2.5 (µg/m3) ",palette = c("purple", "orange", "green3"),
              cor.coeff.args = list(method = "spearman"),add.params = list(color = "black", fill = "lightgray", linetype=2))+labs(color="")+scale_y_log10()

p2<-p2+labs(y = ~ atop(paste('',italic("A.fumigatus")),
                       paste("relative abundance (log",
                             scriptstyle("10"),")")))

p1<-ggscatter(copd_bf, y = "A.fumigatus", x = "pm10", color = "SampleType",
              add = "reg.line", conf.int = TRUE, 
              cor.coef = TRUE, cor.method = "spearman",xlab = "PM 10 (µg/m3) ",palette = c("purple", "orange", "green3"),
              cor.coeff.args = list(method = "spearman"), add.params = list(color = "black", fill = "lightgray", linetype=2))+yscale("log10", .format = TRUE)+labs(color="")+scale_y_log10()

p1<-p1+labs(y = ~ atop(paste('',italic("A.fumigatus")),
                       paste("relative abundance (log",
                             scriptstyle("10"),")")))


pti<-ggscatter(copd_bf, y = "A.fumigatus", x = "Temp", color = "SampleType",
               add = "reg.line", conf.int = TRUE,  
               cor.coef = TRUE, cor.method = "spearman",xlab = "Temperature (°C)",palette = c("purple", "orange", "green3"),
               cor.coeff.args = list(method = "spearman"),add.params = list(color = "black", fill = "lightgray", linetype=2))+yscale("log10", .format = TRUE)+labs(color="")+scale_y_log10()

pti<-pti+labs(y = ~ atop(paste('',italic("A.fumigatus")),
                         paste("relative abundance (log",
                               scriptstyle("10"),")")))

pri<-ggscatter(copd_bf, y = "A.fumigatus", x = "RH", color = "SampleType",
               add = "reg.line", conf.int = TRUE,  
               cor.coef = TRUE, cor.method = "spearman",xlab = "Relative Humidity (%)",palette = c("purple", "orange", "green3"),
               cor.coeff.args = list(method = "spearman"),add.params = list(color = "black", fill = "lightgray", linetype=2))+yscale("log10", .format = TRUE)+labs(color="")+scale_y_log10()

pri<-pri+labs(y = ~ atop(paste('',italic("A.fumigatus")),
                         paste("relative abundance (log",
                               scriptstyle("10"),")")))


tiff("fig2.tiff", units="in", width=10, height=8, res=400)
ggarrange(pti,p2,p1,pri,labels = NA,common.legend = T,legend="right", ncol=2,nrow=2)
dev.off()

#Figure 3 
#F3A-B

e1<-ggscatter(copd_bf, y = "A.fumigatus", x = "totalexac", color = "SampleType",
              add = "reg.line", conf.int = TRUE, 
              cor.coef = TRUE, cor.method = "spearman",ylab="A. fumigatus\n relative abundance (log 10)",xlab = "No. of Exacerbations\n in the preceding year",palette = c("purple", "orange", "green3"),
              cor.coeff.args = list(method = "spearman"),add.params = list(color = "black", fill = "lightgray", linetype=2))+ scale_x_continuous(breaks = seq(1, 12, 3))+labs(color="")+scale_y_log10()

e1<-e1+labs(y = ~ atop(paste('',italic("A.fumigatus")),
                       paste("relative abundance (log",
                             scriptstyle("10"),")")))

e2<-ggscatter(copd_bf, y = "A.fumigatus", x = "Texac1yr", color = "SampleType",
              add = "reg.line", conf.int = TRUE, 
              cor.coef = TRUE, cor.method = "spearman",ylab="A. fumigatus\n relative abundance (log 10)",xlab = "No. of Exacerbations\n during 1 year follow-up",palette = c("purple", "orange", "green3"),
              cor.coeff.args = list(method = "spearman"),add.params = list(color = "black", fill = "lightgray", linetype=2))+labs(color="")+scale_y_log10()

e2<-e2+labs(y = ~ atop(paste('',italic("A.fumigatus")),
                       paste("relative abundance (log",
                             scriptstyle("10"),")")))

#OR
my_y_title <- expression(paste("", italic("A.fumigatus "), " RA (log 10)"))
e1<-e1+labs(y =my_y_title)

tiff("fig3a.b.tiff", units="in", width=9, height=4, res=400)
ggarrange(e1,e2,labels = NA,common.legend = T,legend="right", ncol=2,nrow=1)
dev.off()

#F3C-D
in_copd<-subset(Clinical_data, Clinical_data$Study=="Singapore") #subset singapore data

kruskal.test(in_copd$Total_Exacerbation~in_copd$aspf_env)#0.02255
kruskal.test(in_copd$Texac1yr~in_copd$aspf_env)#0.002975

b<-ggplot(in_copd, aes(x=factor(aspf_env), y=Total_Exacerbation))+geom_boxplot(outlier.color = NA, color="black")+geom_point(aes(fill=factor(aspf_env)),position = position_jitter(0.1),size=2.5, shape=21)+theme_pubr()+scale_fill_manual(values = c("grey", "palevioletred1","blue"))+scale_x_discrete(labels=c("Moderate","High"))+xlab("A. fumigatus")+ylab("No. of Exacerbations\n in the preceding year")+ylim(0,6)+labs(fill="")+theme(axis.title.x = element_text(face = "italic"))
ex1<-ggplot(in_copd, aes(x=factor(aspf_env), y=Texac1yr))+geom_boxplot(outlier.color = NA, color="black")+geom_point(aes(fill=factor(aspf_env)),position = position_jitter(0.1),size=2.5, shape=21)+theme_pubr()+scale_fill_manual(values = c("grey", "palevioletred1","blue"))+scale_x_discrete(labels=c("Moderate","High"))+xlab("A. fumigatus")+ylab("No. of Exacerbations\nduring 1 year follow-up")+ylim(0,6)+labs(fill="")+theme(axis.title.x = element_text(face = "italic"))

tiff("fig3c.d.tiff", units="in", width=9, height=4, res=400)
ggarrange(b,ex1,labels = NA,common.legend = T,legend="none", ncol=2,nrow=1)
dev.off()


#F3E-F
#calculate IRR- exacerbation preceding year
x<-glm.nb(Total_Exacerbation~factor(aspf_env)+Age+BMI+Smoking.Status+FEV1_percent_predicted,data=in_copd)
summary(x)
exp(coef(x))
#calculate IRR- exacerbation at one year follow-up
x<-glm.nb(Texac1yr~factor(aspf_env)+Age+BMI+Smoking.Status+FEV1_percent_predicted+factor(FE),data=in_copd)
summary(x)

#Forest plot for exacerbation at one year
asp_exc1yr$yAxis<-factor(asp_exc1yr$yAxis, levels = c("Age","BMI","Ex Smoker","FEV1 %predicted","Non-FE","Asp f"))

p1 <- ggplot(asp_exc1yr  , aes(x = boxOdds, y = yAxis))
p1<-p1 + geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed") +
  geom_errorbarh(aes(xmax = boxCIHigh, xmin = bocCILow), size = .5, height = .2, color = "gray50") +
  geom_point(aes(fill=p),size = 2.0, shape=21) +
  theme_pubr() +
  theme(panel.grid.minor = element_blank()) +      
  ylab("") +
  xlab("Risk Ratio (RR) with 95% CI\n Exacerbations during 1 year follow-up") +scale_fill_manual(values = c("grey","red"))+theme(legend.position =" none",axis.title.x = element_text(size=9), axis.text.x = element_text(size = 9), axis.text.y = element_text(size=9))+scale_x_log10()

my.labels=c("Age","BMI",
            "Ex-smoker",expression(paste("FEV", scriptstyle("1")," % predicted")), "Non-FE",expression(paste("",italic("A.fumigatus"))))
p1<-p1+scale_y_discrete(labels=my.labels)

#Forest plot exacerbation at preceding year
asp_exc_b$yAxis<-factor(asp_exc_b$yAxis, levels = c("Age","BMI","Ex Smoker","FEV1 %predicted","A.fumigatus"))

p2 <- ggplot(asp_exc_b  , aes(x = boxOdds, y = yAxis))
p2<-p2 + geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed") +
  geom_errorbarh(aes(xmax = boxCIHigh, xmin = bocCILow), size = .5, height = .2, color = "gray50") +
  geom_point(aes(fill=p),size = 2.0, shape=21) +
  theme_pubr() +
  theme(panel.grid.minor = element_blank()) +      
  ylab("") +
  xlab("Risk Ratio (RR) with 95% CI\n Exacerbations in the preceding year") +scale_fill_manual(values = c("grey","red"))+theme(legend.position =" none",axis.title.x = element_text(size=9), axis.text.x = element_text(size = 9), axis.text.y = element_text(size=9))+scale_x_log10()

my.labels=c("Age","BMI",
            "Ex-smoker",expression(paste("FEV", scriptstyle("1")," % predicted")),expression(paste("",italic("A.fumigatus"))))
p2<-p2+scale_y_discrete(labels=my.labels)

tiff("fig3e_f.tiff", units="in", width=8, height=3, res=400)
ggarrange(p2,p1,ncol=2, nrow=1)
dev.off()

#Figure 4 A
sga<-subset(Asp_allergen,Asp_allergen$Study=="Singapore") #subset singapore data
f4<-sga[,c(20:42)]
f4<-t(f4)

tiff("fig4a.tiff", units="in", width=10, height=4, res=400)
pheatmap(f4,cluster_cols = F, cluster_rows = F, fontsize = 9,show_colnames=F, color = 
           c("aliceblue","aliceblue","white","lightskyblue1","dodgerblue","royalblue4"),scale = "row",gaps_col = c(42,42,86,86))
dev.off()

#Figure 4 B-D

# to exclude row with missing data 
data_subset <- Asp_allergen[ , c("af3_pos2","af5_pos2","af6_pos2")] 
df <- Asp_allergen[complete.cases(data_subset), ]

kruskal.test(Asp_allergen$Asp.f.3~Asp_allergen$af3_pos2)#0.043
kruskal.test(Asp_allergen$Asp.f.5~Asp_allergen$af5_pos2)#0.8486
kruskal.test(Asp_allergen$Asp.f.6~Asp_allergen$af6_pos2)#0.4669

#value 
nx<-subset(Asp_allergen,Asp_allergen$af3_pos2=="0")
sx<-subset(Asp_allergen,Asp_allergen$af3_pos2=="1")
summary(sx$Asp.f.3)
summary(nx$Asp.f.3)

ggplot(df, aes(x=factor(af1_pos2), y=Asp.f.1))+geom_boxplot(outlier.color = NA, color="black")+geom_point(aes(fill=factor(SampleType)),position = position_jitter(0.1),size=1.5, shape=21)+theme_pubr()+scale_fill_manual(values = c("purple", "orange2","green3"))+scale_y_log10()+scale_y_log10()+xlab("rAsp f3")+ylab("Asp f 3\nRelative Abundance (%)")+scale_x_discrete(labels=c("Non-Sensitized","Sensiized"))+labs(fill="")
f3<-ggplot(df, aes(x=factor(af3_pos2), y=Asp.f.3.RA))+geom_boxplot(outlier.color = NA, color="black")+geom_point(aes(fill=factor(SampleType)),position = position_jitter(0.1),size=1.5, shape=21)+theme_pubr()+scale_fill_manual(values = c("purple", "orange2","green3"))+scale_y_log10()+scale_y_log10()+xlab("Asp f 3")+ylab("Asp f 3\nRelative Abundance (%)")+scale_x_discrete(labels=c("Non-Sensitized","Sensitized"))+labs(fill="")
f5<-ggplot(df, aes(x=factor(af5_pos2), y=Asp.f.5.RA))+geom_boxplot(outlier.color = NA, color="black")+geom_point(aes(fill=factor(SampleType)),position = position_jitter(0.1),size=1.5, shape=21)+theme_pubr()+scale_fill_manual(values = c("purple", "orange2","green3"))+scale_y_log10()+scale_y_log10()+xlab("Asp f 5")+ylab("Asp f 5\nRelative Abundance (%)")+scale_x_discrete(labels=c("Non-Sensitized","Sensitized"))+labs(fill="")
f6<-ggplot(df, aes(x=factor(af6_pos2), y=Asp.f.6.RA))+geom_boxplot(outlier.color = NA, color="black")+geom_point(aes(fill=factor(SampleType)),position = position_jitter(0.1),size=1.5, shape=21)+theme_pubr()+scale_fill_manual(values = c("purple", "orange2","green3"))+scale_y_log10()+scale_y_log10()+xlab("Asp f 6")+ylab("Asp f 6\nRelative Abundance (%)")+scale_x_discrete(labels=c("Non-Sensitized","Sensitized"))+labs(fill="")

tiff("fig4b_d.tiff", units="in", width=12, height=4, res=400)
ggarrange(f3,f5,f6,labels = NA,common.legend = T,legend="right", ncol=3,nrow=1)
dev.off()

#Fig 4 E-H
in_copd<-subset(Clinical_data, Clinical_data$Study=="Singapore")

f3te<-ggplot(data=subset(in_copd,!is.na(f3_env)), aes(x=factor(f3_env), y=Total_Exacerbation))+geom_boxplot(aes(fill=factor(f3_env), alpha=0.3),outlier.color = NA, color="black")+geom_point(aes(fill=factor(f3_env)),position = position_jitter(0.1),size=1.5, shape=21)+theme_pubr()+scale_fill_manual(values = c("steelblue", "hotpink","green3"))+xlab("Asp f 3\n")+ylab("No. of Exacerbation\nin the preceding year")+scale_y_continuous(breaks = c(0,2,4,6,8,10))+labs(fill="")+scale_x_discrete(label=c("E", "E+S"))
f3t1<-ggplot(data=subset(in_copd,!is.na(f3_env)), aes(x=factor(f3_env), y=Texac1yr))+geom_boxplot(aes(fill=factor(f3_env), alpha=0.3),outlier.color = NA, color="black")+geom_point(aes(fill=factor(f3_env)),position = position_jitter(0.1),size=1.5, shape=21)+theme_pubr()+scale_fill_manual(values = c("steelblue", "hotpink","green3"))+xlab("Asp f 3\n")+ylab("No. of Exacerbation\nduring 1 year follow-up")+labs(fill="")+scale_x_discrete(label=c("E", "E+S"))
f3f<-ggplot(data=subset(in_copd,!is.na(f3_env)), aes(x=factor(f3_env), y=FEV1_percent_predicted))+geom_boxplot(aes(fill=factor(f3_env), alpha=0.3),outlier.color = NA, color="black")+geom_point(aes(fill=factor(f3_env)),position = position_jitter(0.1),size=1.5, shape=21)+theme_pubr()+scale_fill_manual(values = c("steelblue", "hotpink","green3"))+xlab("Asp f 3\n")+ylab("FEV1 % predicted")+labs(fill="")+scale_x_discrete(label=c("E", "E+S"))
f3f<-f3f+labs(y=expression(paste("FEV", scriptstyle("1")," % predicted")))
f3c<-ggplot(data=subset(in_copd,!is.na(f3_env)), aes(x=factor(f3_env), y=CAT_score))+geom_boxplot(aes(fill=factor(f3_env), alpha=0.3),outlier.color = NA, color="black")+geom_point(aes(fill=factor(f3_env)),position = position_jitter(0.1),size=1.5, shape=21)+theme_pubr()+scale_fill_manual(values = c("steelblue", "hotpink","green3"))+xlab("Asp f 3\n")+ylab("CAT score")+labs(fill="")+scale_x_discrete(label=c("E", "E+S"))

tiff("fig4e_h.tiff", units="in", width=8, height=6, res=400)
ggarrange(f3f,f3t1,f3c,f3te,labels = NA,common.legend = F,legend=F, ncol=2,nrow=2)
dev.off()

kruskal.test(in_copd$Total_Exacerbation~in_copd$f3_env)#0.2963
kruskal.test(in_copd$Texac1yr~in_copd$f3_env)#0.0443
kruskal.test(in_copd$FEV1_percent_predicted~in_copd$f3_env)#0.0277
kruskal.test(in_copd$CAT_score~in_copd$f3_env)#0.3998

#value
ee<-subset(in_copd,in_copd$f3_env=="E")
ees<-subset(in_copd,in_copd$f3_env=="ES")
summary(ees$FEV1_percent_predicted)
summary(ee$FEV1_percent_predicted)
summary(ees$Texac1yr)
summary(ee$Texac1yr)

#Fig 5 A

idr<-subset(Asp_ab, SampleType%in% c('Indoor air' , 'Surface'))#subset indoor data including air and surface

mAirData <- melt(idr)
mAirData$SampleType<-factor(mAirData$SampleType, levels = c("Indoor air", "Surface"))
mAirData$variable<-factor(mAirData$variable, levels = c("Aspergillus.terreus", "Aspergillus.mulundensis", "Aspergillus.novofumigatus","Aspergillus.heteromorphus","Aspergillus.niger","Aspergillus.clavatus","Aspergillus.glaucus", "Aspergillus.steynii","Aspergillus.pseudoglaucus","Aspergillus.awamori", "Aspergillus.welwitschiae","Aspergillus.aculeatinus", "Aspergillus.nidulans","Aspergillus.fumigatus" ))
mAirData$variable<-factor(mAirData$variable, labels = c("A.terreus", "A.mulundensis", "A.novofumigatus","A.heteromorphus","A.niger","A.clavatus","A.glaucus", "A.steynii","A.pseudoglaucus","A.awamori", "A.welwitschiae","A.aculeatinus", "A.nidulans","A.fumigatus" ))
indn<-c("Vancouver\nIndoor Air","Vancouver\nSurfaces","Singapore\nIndoor Air","Singapore\nSurfaces")
names(indn)<-c("CIHC","CSHC","SIHC","SSHC")

tiff("bubble_asp_sg_can.copd.tiff", units="in", width=10, height=6, res=400)
ggplot(data=mAirData,aes(x=SampleID, y=variable))+
  geom_point(aes(size=value,fill=variable), alpha=0.85, shape=21) +
  scale_fill_manual(values = c("indianred1","dodgerblue1","thistle","lightskyblue", "darkolivegreen1","cyan3","bisque", "royalblue", "khaki2", "pink", 
                                           "yellow", "limegreen", "coral2", "violetred", "indianred2", "darkblue", 
                                           "darkorange1", "cyan1", "royalblue4", "maroon2", "darkblue", "aquamarine", 
                                           "deepskyblue2", "slateblue", "dodgerblue3", "darkolivegreen1", "darkgoldenrod1", 
                                           "violetred", "grey","slategrey","black", "pink3", 
                                           "bisque", "lightblue", "darkblue", "cadetblue", "indianred1", "turquoise2", 
                                           "cyan1", "cyan3", "thistle", "salmon2", "blue2" )) +
                                             xlab("")+ylab("")+theme_classic()+labs(fill="") +scale_size_continuous(limits = c(0,1000), range = c(1,30), breaks = c(0,50,100,500,1000),name="Metagenome reads")+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())+facet_wrap(~grp, scales = "free_x",ncol=4,labeller = labeller(grp=indn))+  guides(fill="none")+theme(axis.text.y = element_text(face="italic"))
dev.off()

#Figure 5 B-D


f2<-cbind(env_data,Master)
canc<-subset(f2,f2$Study=="Vancouver")#subset vancouver data

gec<-ggscatter(canc, y = "A.fumigatus", x = "totalexac", color = "SampleType",
               add = "reg.line", conf.int = TRUE, 
               cor.coef = TRUE, cor.method = "spearman",xlab = "No. of Exacerbations\nin the preceding year",palette = c("purple", "orange", "green3"),
               cor.coeff.args = list(method = "spearman"),add.params = list(color = "black", fill = "lightgray", linetype=2))+labs(color="")+scale_y_continuous(limits = c(0,0.75))
gec<-gec+labs(y = ~ atop(paste('',italic("A.fumigatus")),
                         paste("relative abundance")))

gcc<-ggscatter(canc, y = "A.fumigatus", x = "Cat.score", color = "SampleType",
               add = "reg.line", conf.int = TRUE, 
               cor.coef = TRUE, cor.method = "spearman",xlab = "CAT score",palette = c("purple", "orange", "green3"),
               cor.coeff.args = list(method = "spearman"),add.params = list(color = "black", fill = "lightgray", linetype=2))+labs(color="")+scale_y_continuous(limits = c(0,0.75))
gcc<-gcc+labs(y = ~ atop(paste('',italic("A.fumigatus")),
                         paste("relative abundance")))

gfc<-ggscatter(canc, y = "A.fumigatus", x = "FEV1_percent_predicted", color = "SampleType",
               add = "reg.line", conf.int = TRUE, 
               cor.coef = TRUE, cor.method = "spearman",xlab = "FEV1 % predicted",palette = c("purple", "orange", "green3"),
               cor.coeff.args = list(method = "spearman"),add.params = list(color = "black", fill = "lightgray", linetype=2))+labs(color="")+scale_y_continuous(limits = c(0,0.75))
gfc<-gfc+labs(y = ~ atop(paste('',italic("A.fumigatus")),
                         paste("relative abundance")))
gfc<-gfc+labs(x=expression(paste("FEV", scriptstyle("1")," % predicted")))

tiff("fig5b_d.tiff", units="in", width=12, height=4, res=400)
ggarrange(gec,gcc,gfc,labels = NA,common.legend = T,legend="right", ncol=3,nrow=1)
dev.off()

#Figure 6
Asp_allergen$cat10<-Asp_allergen$CAT_score>9
aspc<-subset(Asp_allergen, SampleType%in% c('Indoor' , 'Outdoor')) #subset indoor and outdoor samples
aspc$Study<-factor(aspc$Study, levels = c("Vancouver","Singapore"))

ffe<-ggplot(aspc, aes(x=factor(FE), y=asp_al))+geom_boxplot(outlier.color = NA, color="black")+geom_point(aes(fill=factor(Study)),position = position_jitter(0.1),size=1.5, shape=21)+theme_pubr()+ylab("No.of Asp f allergens ")+xlab("Frequent Exacerbator\n")+scale_x_discrete(labels=c("No","Yes"))+scale_fill_manual(values = c("red", "blue","blue"), labels = c("Vancouver","Singapore"))+labs(fill="")
ce<-ggplot(aspc, aes(x=factor(Study), y=asp_al))+geom_boxplot(outlier.color = NA, color="black")+geom_point(aes(fill=factor(Study)),position = position_jitter(0.1),size=1.5, shape=21)+theme_pubr()+ylab("No.of Asp f allergens ")+xlab("Location\n")+scale_x_discrete(labels=c("Vancouver","Singapore"))+scale_fill_manual(values = c("red", "blue","blue"), labels = c("Vancouver","Singapore"))+labs(fill="")
ccat<-ggplot(aspc, aes(x=factor(cat10), y=asp_al))+geom_boxplot(outlier.color = NA, color="black")+geom_point(aes(fill=factor(Study)),position = position_jitter(0.1),size=1.5, shape=21)+theme_pubr()+ylab("No.of Asp f allergens ")+xlab("CAT score")+scale_x_discrete(labels=c("<10","≥10"))+scale_fill_manual(values = c("red", "blue","blue"), labels = c("Vancouver","Singapore"))+labs(fill="")
cfev<-ggplot(aspc, aes(x=factor(gold_stage), y=asp_al))+geom_boxplot(outlier.color = NA, color="black")+geom_point(aes(fill=factor(Study)),position = position_jitter(0.1),size=1.5, shape=21)+theme_pubr()+ylab("No.of Asp f allergens ")+xlab("GOLD (FEV1) group")+scale_fill_manual(values = c("red", "blue","blue"), labels = c("Vancouver","Singapore"))+labs(fill="")
cfev<-cfev+labs(x=expression(paste("GOLD ","(FEV", scriptstyle("1")," % predicted)"," group")))

tiff("Fig6.tiff", units="in", width=8, height=6, res=400)
ggarrange(ce,ffe,ccat,cfev,labels = NA,common.legend = T,legend="right", ncol=2,nrow=2)
dev.off()

#value
ffe<-subset(aspc,aspc$FE=="1")
nffe<-subset(aspc,aspc$FE=="0")
summary(nffe$asp_al)
summary(ffe$asp_al)

kruskal.test(asp_al~FE, data=aspc)#0.028
kruskal.test(asp_al~Study, data=aspc)#0.00047
kruskal.test(asp_al~cat10, data=aspc)#0.2531
kruskal.test(asp_al~gold_stage, data=aspc)#0.4647


####Supplementary Figure 1-blank

tiff("e-fig1a.tiff", units="in", width=8, height=5, res=400)
blank$Specimen.type<-factor(blank$Specimen.type, levels = c("Sequencing Blank","Filter Blank","Reagent Blank", "Sample"))
ggplot(blank, aes(x=ID, y=Read.count))+geom_point(aes(shape=Specimen.type,fill=Specimen.type),position = position_jitter(0.1),size=2, color="black")+theme_pubr()+scale_shape_manual(values = c(21,22,23,1))+scale_fill_manual(values=c("red","blue","green","white"))+xlab("Specimen")+theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())+ylab("Read counts")+theme(legend.title = element_text(size=10), legend.position = "right")
dev.off()

#e-fig1b
cc<-blank[,c(2,5,7:45)]#others as 0 
mAirData<-melt(cc)

tiff("e-fig1b.tiff", units="in", width=8, height=5, res=400)
ggplot(data=mAirData,aes(x=Type, y=value, fill=variable))+
  geom_bar(stat="identity", position="fill") +
  scale_fill_manual(values = c("slateblue", "maroon2","yellow","darkblue","cyan1", "darkorange3", "grey55", 
                                          "yellow3", "steelblue", "lightblue", "royalblue", "lavender","green3",
                                          "darksalmon", "cyan3", "red","violet", "seagreen4", 
                                          "pink1", "gold2", "dodgerblue3", "darkgoldenrod1", "red3", 
                                          "blue2","grey", "darkgreen","salmon2", "lightgreen", "lightskyblue", "darkorchid4",
                                          "bisque", "lightsteelblue2", "pink", "cadetblue", "indianred3", "turquoise3", 
                                          "lavenderblush2", "blueviolet", "honeydew2", "salmon2", "blue2" )) +
                                            xlab("")+ylab("Relative abundance (%)")+theme_classic()+labs(fill="")+scale_y_continuous(labels = scales::percent)+scale_x_discrete(labels=c("Sequencing blanks\n(n=6)","Filter Blanks\n(n=9)", "Reagent Blanks\n(n=6)","Samples\n(n=157)"))+theme(legend.text = element_text(face="italic"))
dev.off()



#Supplementary Figure 2
c1<-ggscatter(copd_bf, y = "A.fumigatus", x = "Cat.score", color = "SampleType",
              add = "reg.line", conf.int = TRUE, 
              cor.coef = TRUE, cor.method = "spearman",ylab="A. fumigatus\n relative abundance (log 10)",xlab = "CAT score",palette = c("purple", "orange", "green3"),
              cor.coeff.args = list(method = "spearman"),add.params = list(color = "black", fill = "lightgray", linetype=2))+labs(color="")+ scale_y_log10(limit=c(0.01,10))

c1<-c1+labs(y = ~ atop(paste('',italic("A.fumigatus")),
                       paste("relative abundance (log",
                             scriptstyle("10"),")")))

f1<-ggscatter(copd_bf, y = "A.fumigatus", x = "FEV1_percent_predicted", color = "SampleType",
              add = "reg.line", conf.int = TRUE, 
              cor.coef = TRUE, cor.method = "spearman",ylab="A. fumigatus\n relative abundance (log 10)",xlab = "FEV1 % predicted",palette = c("purple", "orange", "green3"),
              cor.coeff.args = list(method = "spearman"),add.params = list(color = "black", fill = "lightgray", linetype=2))+labs(color="")+ scale_y_log10(limit=c(0.01,10))

f1<-f1+labs(y = ~ atop(paste('',italic("A.fumigatus")),
                       paste("relative abundance (log",
                             scriptstyle("10"),")")))
f1<-f1+labs(x=expression(paste("FEV", scriptstyle("1")," % predicted")))

tiff("efig2.tiff", units="in", width=10, height=3, res=400)
ggarrange(c1,f1,labels = NA,common.legend = T,legend="right", ncol=2,nrow=1)
dev.off()

#Supplementary Figure 3

tt<-ggscatter(f3, y = "Temp", x = "Total_Exacerbation", 
              add = "reg.line", conf.int = TRUE, 
              cor.coef = TRUE, cor.method = "spearman",xlab = "No. of Exacerbations\nin the preceding year",
              cor.coeff.args = list(method = "spearman", label.x.npc = "left", label.y = 38),add.params = list(color = "blue", fill = "lightgray", linetype=2))+ scale_x_continuous(breaks = seq(1, 12, 3))+labs(color="")+ylab("Temperature (°C)")

tr<-ggscatter(f3, y = "RH", x = "Total_Exacerbation", 
              add = "reg.line", conf.int = TRUE, 
              cor.coef = TRUE, cor.method = "spearman",xlab = "No. of Exacerbations\nin the preceding year",
              cor.coeff.args = list(method = "spearman", label.x.npc = "left", label.y = 90),add.params = list(color = "blue", fill = "lightgray", linetype=2))+ scale_x_continuous(breaks = seq(1, 12, 3))+labs(color="")+ylab("Relative humidity (%)")

tp<-ggscatter(f3, y = "pm10", x = "Total_Exacerbation",
              add = "reg.line", conf.int = TRUE, 
              cor.coef = TRUE, cor.method = "spearman",xlab = "No. of Exacerbations\nin the preceding year",
              cor.coeff.args = list(method = "spearman", label.x.npc = "left", label.y = 70),add.params = list(color = "blue", fill = "lightgray", linetype=2))+ scale_x_continuous(breaks = seq(1, 12, 3))+labs(color="")+ylab("PM 10 (µg/m3)")

tp2<-ggscatter(f3, y = "pm25", x = "Total_Exacerbation",
               add = "reg.line", conf.int = TRUE, 
               cor.coef = TRUE, cor.method = "spearman",xlab = "No. of Exacerbations\nin the preceding year",
               cor.coeff.args = list(method = "spearman", label.x.npc = "left", label.y = 150),add.params = list(color = "blue", fill = "lightgray", linetype=2))+ scale_x_continuous(breaks = seq(1, 12, 3))+labs(color="")+ylab("PM 2.5 (µg/m3)")

tiff("e-fig3.tiff", units="in", width=10, height=3, res=400)
ggarrange(tt,tr,tp2,tp,labels = NA,common.legend = T,legend="bottom", ncol=4,nrow=1)
dev.off()

##Supplementary Figure 4

library(ggbreak)
mAirData <- melt(sputum_ITS )

tiff("e-fig4.tiff", units="in", width=6, height=6, res=400)

p1<-ggplot(data=mAirData,aes(x=Sputum, y=value, fill=variable))+
  geom_bar(stat="identity", position="fill") +
  scale_fill_manual(values = c("royalblue2","limegreen","red2","sandybrown","hotpink1","skyblue2","orangered","mediumspringgreen","indianred1","darkviolet","deepskyblue2","darkorchid4","darkcyan","orange","blue2","rosybrown1",
                                           "darkseagreen2","cyan2","maroon3","slateblue2","cadetblue3","black","chocolate2", "lavender","hotpink3", "seashell1","yellow2","grey","burlywood1","black")) +
                                             xlab("")+ylab("Relative abundance (%)")+theme_classic()+labs(fill="")+scale_y_continuous(labels = scales::percent)+theme(axis.text.x = element_blank())+theme(axis.ticks.x = element_blank())+theme(legend.text = element_text(face="italic"))

p1+ scale_y_break(c(0.4, 0.99))+scale_y_continuous(breaks = c(0,0.1,0.2,0.3,0.4,1.0),labels = scales::percent)
dev.off()

f2<-cbind(env_data,Master)
sp<-subset(f2, sputum_asp%in% c('Yes' , 'No'))

tiff("e-fig4b.tiff", units="in", width=5, height=4, res=400)
p2<-ggplot(sp, aes(x=factor(sputum_asp), y=A.fumigatus))+geom_boxplot(outlier.color = NA, color="black")+geom_point(aes(fill=factor(SampleType)),position = position_jitter(0.1),size=2.5, shape=21)+theme_pubr()+stat_compare_means( method="kruskal.test", aes(label = ..p.signif..), size=4,label.x = 1.5)+scale_fill_manual(values = c("purple", "orange","green2","grey"))+xlab("Sputum Aspergillus")+ylab("Environmental Aspergillus\n Relative abundance (%)")+labs(fill="")+ theme(legend.position="right")
p2+labs(y = ~ atop(paste('Environmental ',italic("Aspergillus")),
                   paste("Relative abundance (%)") ))+labs(x = ~ atop(paste('Sputum ',italic("Aspergillus")),paste('(by 18S ITS Sequencing)') ))
dev.off()

#sSupplementary Fig 5A
f2<-cbind(Master,env_data)

a<-f2[,c(4:89)]
pcoa<- vegdist(a, "bray")
pcoa<-as.matrix(pcoa)
BrayCurtMbiome=cmdscale(pcoa) 
ordiplot (BrayCurtMbiome, display = 'species', type = 'text') 
BCords<-scores(BrayCurtMbiome)
plot(scores(pcoa))
text(pcoa, row.names(pcoa), cex=0.6, pos=4, col="red") 

Master=a[order(row.names(a)),]
Master<-cbind(Master,BCords[order(row.names(BCords)),])
Master<-cbind(Master,labels = f2[order(row.names(f2)),])

gg <- data.frame(cluster=factor(Master$labels.Type), x=Master$Dim1, y=Master$Dim2)
centroids <- aggregate(cbind(x,y)~cluster,data=gg,mean)
gg <- merge(gg,centroids,by="cluster",suffixes=c("",".centroid"))



tiff("e-fig5a.tiff", units="in", width=6, height=4, res=400)
ggplot(gg) +
  scale_linetype_identity() + 
  geom_segment(aes(x=x.centroid, y=y.centroid, xend=x, yend=y, colour = cluster), alpha=0.4)+
  geom_point(data=centroids, aes(x=x, y=y,  fill=cluster, shape=cluster),shape=21, size=6) +geom_point(data=centroids, aes(x=x, y=y), shape=7, color="black", size=6)+
  geom_point(aes(x=x,y=y, colour = cluster ), size=2, alpha=0.4) +
  scale_colour_manual(values = c("blue","dodgerblue","skyblue2","red","red3","darkred","darkgreen","skyblue2"))+scale_fill_manual(values = c("blue","dodgerblue","skyblue","red","red3","darkred","darkgreen","skyblue2"))+
  labs(colour="",  
       y = "PC2", x = "PC1 ")+ labs(color='')+ 
  theme(legend.position = "bottom")+theme_classic()+guides(fill="none")
dev.off()

#Supplementary Fig 5B
a<-biome.Species[,c(4,6,50:79)]

mAirData <- melt(a)

tiff("species_country_um2.tiff", units="in", width=10, height=6, res=400)

ggplot(data=mAirData,aes(x=SampleType, y=value, fill=variable))+
  geom_bar(stat="identity", position="fill") +
  scale_fill_manual(values = c( "hotpink2","lightblue","darkorchid", "coral1", "darkcyan", "aquamarine3", 
                                          "cyan1", "yellow3", "yellow", "lightblue", "slateblue", 
                                          "darksalmon", "darkorange1", "royalblue1", "lavender","violet", "seagreen3", 
                                          "pink1", "gold2", "dodgerblue3","forestgreen", "darkgoldenrod1", 
                                          "red","red3", "salmon2", "lightskyblue","greenyellow", "grey44", 
                                          "bisque", "grey", "darkblue", "cadetblue", "indianred1", "turquoise2", 
                                          "cyan1", "cyan3", "gold2", "salmon2", "blue2" )) +
                                            xlab("")+ylab("Relative abundance (%)")+theme_classic()+labs(fill="")+scale_y_continuous(labels = scales::percent)+facet_grid(~Group)+theme(legend.text = element_text(face="italic"))

dev.off()


#Supplementary Figure 6
asp.country1$Species<-factor(asp.country1$Species, levels = c("A.terreus", "A.mulundensis", "A.novofumigatus","A.heteromorphus","A.niger","A.clavatus","A.glaucus", "A.steynii","A.pseudoglaucus","A.awamori", "A.welwitschiae","A.aculeatinus", "A.nidulans","A.fumigatus" ))

tiff("asp_barchart_sg.tiff", units="in", width=8, height=8, res=400)

sf5<-ggplot(data=asp.country1,aes(x=Species, y=Value,fill=Country))+
  geom_bar(stat="identity", position="dodge") +
  scale_fill_manual(values = c("blue", "red","green2" )) + xlab("")+ylab("")+theme_classic()+labs(fill="")+ylab("Percentage (%) of home with \ndetectable Aspergillus spp.")+coord_flip()+theme(axis.text.y = element_text(face ="italic", size = 12, color="red", family="Calibri"))+theme(axis.text.x = element_text(size = 12))+ theme(legend.position="top")+theme(legend.text = element_text(size=12))+theme(axis.title.x = element_text(size=12))

sf5+labs(y = ~ atop(paste('Percentage (%) of home with detectable ',italic("Aspergillus spp.")) ))

dev.off()


#Supplementary Figure 7

asp<-subset(Asp_allergen, SampleType%in% c('Indoor' , 'Outdoor'))#exclude surfaces

f4c<-asp[,c(20:34,36:42)]#exclude Asp.f17 as is zero in all samples
f4c<-t(f4c)

tiff("allergen_io_ab_pheatmap_all.tiff", units="in", width=10, height=4, res=400)
pheatmap(f4c,cluster_cols = F, cluster_rows = F, fontsize = 9, show_colnames=F,color = 
           c("aliceblue","aliceblue","white","lightskyblue1","dodgerblue","royalblue4"),scale = "row",gaps_col = c(43,43,86,86,98,98))
dev.off()


#Table 1
library(table1)
library(flextable)

#p-value
pvalue <- function(x, ...) {
  # Construct vectors of data y, and groups (strata) g
  y <- unlist(x)
  g <- factor(rep(1:length(x), times=sapply(x, length)))
  if (is.numeric(y)) {
    # For numeric variables, perform a standard 2-sample t-test
    p <- kruskal.test(y ~ g)$p.value
  } else {
    # For categorical variables, perform a chi-squared test of independence
    p <- chisq.test(table(y, g))$p.value
  }
  # Format the p-value, using an HTML entity for the less-than sign.
  # The initial empty string places the output on the line below the variable label.
  c("", sub("<", "&lt;", format.pval(p, digits=3, eps=0.001)))
}


table1(~ Age+Sex.at.birth+BMI+Smoking.Status+packyears+CAT_score+FE+FEV1_FVC_ratio_percent_predicted+FEV1_percent_predicted+Total_Exacerbation+RH+Temp+treatment+type_of_house| Study, render.continuous = c(.="Median [Q1, Q3]"),data=Clinical_data,extra.col=list(`P-value`=pvalue))

###############
