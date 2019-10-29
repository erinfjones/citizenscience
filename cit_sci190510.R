###### Citizen Science Spatial statistics ######
###### Adapted from BWAbbott and RJFrei
###### Compiled by EFJones (erinfjones3@gmail.com) 5/02/19
###### Last updated 10/22/19
library(ggthemes)
library(reshape2) #load package
require(dplyr)
library(ggplot2)
library(changepoint)
require(cowplot)
require(data.table)

setwd("~/Box Sync/BYU/CitSci/CitSciData/Rdata")
setwd("C:/Users/Erin/Box Sync/AbbottLab/Projects/Utah Lake Citizen Science/CitSciData/Rdata")
theme_set(theme_few())
# DOC NO3 PO4 DIC
#create vector with desired elements for the figure 
main = c("DOC", "PO4.P","DIN.N","TN","Sulfate","Chloride") 
ext = c("Fluoride", "Bromide", "Lithium", "Sodium", "Potassium")
scale_this <- function(x){(x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)}


############## CALCULATE SPATIAL STABILITY ###################
#import data

stoich = read.csv("ULRC_Rchem_new.csv")
#stoich[stoich == 0] <- NA #set all 0's to NA (assuming you want them to be NA)
stoich$Area_SS = as.numeric(stoich$Area_SS)
###step 1: melt and cast data by site and sampling date
ssmelt <- melt(stoich, id.vars = c("SiteID", "Event", "Category", "Area_SS"), measure.vars = c(3:13), na.rm = TRUE) #melt data
ssmelt$value <- as.numeric(ssmelt$value)
sscast <- dcast(ssmelt, SiteID + variable +Category~Event) #cast data to sort into 3 sampling date columns

#sscast$Fall <- as.numeric(sscast$Fall) # set columns to numeric
#sscast$Spring <- as.numeric(sscast$Spring)
#sscast$Summer <- as.numeric(sscast$Summer)
  

###step 2: do the Spearman rank correlation
#group data by variable and calculate spearman rank correlation for all variables between sampling dates (3 comparisons)
stb = sscast %>%
  subset(Category!="Utah Lake")%>%
  group_by(Category,variable)  %>%
  summarize(SumFllSp = cor(Summer, Fall, method="spearman", use="na.or.complete"), 
            SprFllSp = cor(Spring, Fall, method="spearman", use="na.or.complete"), 
            SprSumSp = cor(Spring, Summer, method="spearman", use="na.or.complete"))
stb # view the results to make sure it worked

##NOTE: there is controversy with the use="" part of cor(). There is a nice article online to read if you're worried about it in the future

#calculate the mean of the 3 spearman rank correlations
df.stab <- cbind(stb, SpearMean = rowMeans(stb[3:5])) #create new stability dataframe "df.stab"
df.stab  

stabstats=reshape2::melt(df.stab%>%rename(Solute=variable),id.vars=c("Category", "Solute"), measure.vars=c("SumFllSp", "SprFllSp", "SprSumSp","SpearMean") )
stabanova=aov(value~Solute*Category+ Category*variable +Solute*Category, data=stabstats)
summary(stabanova)
plot(TukeyHSD(stabanova))

###Step 3: plot the rankings as geom_pointrange() 
#plot the spearman rankings
#to use geom_pointrange need to get the max and min values
stabsum = df.stab %>%
  rowwise()%>%
  mutate(Min = min(SumFllSp, SprFllSp, SprSumSp), Max = max(SumFllSp, SprFllSp, SprSumSp))
dodge <- position_dodge(width=0.5)
#now do the plot
p.stab <- ggplot(stabsum[stabsum$variable %in% main,], aes(variable, SpearMean, color=Category)) +  # plot by mean spearman rank and variable (the solutes you're interested in)
  geom_hline(yintercept = .7, color='grey', size=1) + # add the stability line (0.7 is ~the sqrt(0.5))
  geom_errorbar(aes(ymin = Min, ymax = Max), position=dodge, width=0.1)+   # this adds the point range
  geom_pointrange(aes(ymin = Min, ymax = Max), position=dodge) + 
  ylab("Spearman rank correlation coefficient")
  #  geom_text(aes(label = Category), hjust = "inward", vjust = "inward") + # label your points (so you know which solute is which)
 # theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +  # remove x-axis labels
 #ylim(0,1)  # set y-axis range
p.stab

ggsave("figures/staball.png", width =6, height = 4, dpi = 300)

### calculate summary stats ####
SiteStats = ssmelt %>%
  group_by(SiteID, variable) %>%
  summarize(Mean = mean(value), SD = sd(value), CVpercent = sd(value)/mean(value)*100)
write.csv(SiteStats, "output/sitestats.csv")

categorystats= ssmelt[ssmelt$variable %in% main,] %>%
 # subset(Category!="Utah Lake")%>%
  group_by(Category, variable) %>%
  summarize(Mean = mean(value), SD = sd(value))
write.csv(categorystats, "output/categorystats.csv")

eventstats= ssmelt[ssmelt$variable %in% main,] %>%
  # subset(Category!="Utah Lake")%>%
  group_by(Event, variable) %>%
  summarize(Mean = mean(value), SD = sd(value))
write.csv(eventstats, "output/eventstats.csv")


########### SCALED CONCENTRATION #############
scalestats = ssmelt[ssmelt$variable %in% main,] %>%
  subset(Category!="Utah Lake")%>%
  group_by(variable) %>%
  mutate_each(funs(scale_this), value)


anova=aov(value~Event*variable*Category , data=scalestats)
summary(anova)
tukey=TukeyHSD(anova)
plot(tukey)
tukey
capture.output(print(tukey),file="output/tukey.csv")

boxplot= ggplot(scalestats, aes(Event, value, color=Category)) +
  geom_boxplot( notch=TRUE)+facet_wrap(.~variable) +ylim(-1,5) +ylab("Scaled Concentration")
boxplot
ggsave("Figures/scaledconc.png", plot=boxplot, width =8, height = 5, dpi = 300)

#group the data by catchment (individual scaling for each catchment)
streamsN = ssmelt %>%
  subset(Category!="Utah Lake")%>%
  group_by(Category,variable) %>%
  mutate_each(funs(scale_this), value)


plot1=ggplot(streamsN[streamsN$variable %in% main,], aes(x=Area_SS, y=value, color = Category, shape = Event)) + 
  geom_hline(yintercept=0)+
  geom_point(alpha=0.4)+
#  scale_color_manual( values= c("#000043", "#0033FF", "#01CCA4", "#BAFF12", "#FFCC00", "#FF3300", "#BAFF12", "#FFCC00", "#FF3300")) +
  facet_grid(variable~., scales="free")+
  scale_x_continuous(expand=c(.02,2))+
  xlab(expression("Catchment Area (km2)")) +
  ylab(expression("Concentration (scaled)"))
plot1
ggsave("Figures/variancecollapse_main.png", plot = plot1, width = 6, height = 8, dpi=300)

plot2=ggplot(streamsN[streamsN$variable %in% ext,], aes(x=Area_SS, y=value, color = Category, shape = Event)) + 
  geom_hline(yintercept=0)+
  geom_point(alpha=0.4)+
#  scale_color_manual( values= c("#000043", "#0033FF", "#01CCA4", "#BAFF12", "#FFCC00", "#FF3300", "#BAFF12", "#FFCC00", "#FF3300")) +
  facet_grid(variable~., scales="free")+
  scale_x_continuous(expand=c(.02,2))+
  xlab(expression("Catchment Area (km2)")) +
  ylab(expression("Concentration (scaled)"))
plot2
ggsave("Figures/variancecollapse_ext.png", plot = plot2, width = 6, height = 8, dpi=300)



############ CALCULATE LEVERAGE #############
#calculate leverage
lev= ssmelt %>% 
  subset(Category!="Utah Lake")%>%
  group_by( Category, variable) %>%
  mutate(leverage = ((value - value[which.max(Area_SS)])*100/ value[which.max(Area_SS)])* Area_SS/max(Area_SS, na.rm=TRUE))

#save the lev object as a csv
write.csv(lev, "output/leverage.csv")

#import data
lev.new<-read.csv('output/leverage.csv')
str(lev.new)

#create data frame with threshold data
#bpsmain= bps[bps$variable %in% main,]

#make plot of catchment area vs leverage of main parameters

lev.p=ggplot(lev[lev$variable %in% main,], aes(x=Area_SS, y=leverage, color = Category, shape = Event)) + 
  geom_hline(yintercept=0)+
  geom_point(alpha=0.5)+
#  scale_color_manual( values= c("#000043", "#0033FF", "#01CCA4", "#BAFF12", "#FFCC00", "#FF3300","#FFCC00", "#FF3300")) +
  facet_grid(variable~., scales="free")+
  scale_x_continuous(expand=c(.02,2))+
  xlab(expression("Catchment Area (km2)"))+
  ylab(expression("Leverage (%)"))
lev.p

ggsave("figures/leverageinspace_main.png", plot = lev.p, width = 6, height = 8, dpi=300)

lev.ext=ggplot(lev[lev$variable %in% ext,], aes(x=Area_SS, y=leverage, color = Category, shape = Event)) + 
  geom_hline(yintercept=0)+
  geom_point(alpha=0.5)+
#  scale_color_manual( values= c("#000043", "#0033FF", "#01CCA4", "#BAFF12", "#FFCC00", "#FF3300", "#FF3300")) +
  facet_grid(variable~., scales="free")+
  scale_x_continuous(expand=c(.02,2))+
#  coord_cartesian(ylim = c(0, 600))+
  xlab(expression("Catchment Area (km2)"))+
  ylab(expression("Leverage (%)"))
lev.ext
ggsave("figures/leverageinspace_ext.png", plot = lev.ext, width = 6, height = 8, dpi=300)


# Zoomed in plot
levz.p=ggplot(lev[lev$variable %in% main,], aes(x=Area_SS, y=leverage, color = Category, shape = Event)) + 
  geom_hline(yintercept=0)+
  geom_point(alpha=0.5)+
#  scale_color_manual( values= c("#000043", "#0033FF", "#01CCA4", "#BAFF12", "#FFCC00", "#FF3300","#000043")) +
  facet_grid(variable~., scales="free")+
  scale_x_continuous(expand=c(.02,2))+
  xlab(expression("Catchment Area (km2)"))+
  coord_cartesian(ylim = c(0, 600))+
  ylab(expression("Leverage (%)"))
levz.p

#define quantile functions
q95 <- function(x) {quantile(x,probs=0.95)}
q5 <- function(x) {quantile(x,probs=0.05)}

#make boxplot with full scale
leverbox =ggplot(lev[lev$variable %in% main,], aes(x=variable, y=leverage, fill = Category, color=Category)) + 
  geom_hline(yintercept=0)+
  geom_boxplot( alpha= 0.8, position = position_dodge(width=0.75))+
  stat_summary(fun.y=mean, color = "black", geom="crossbar", size=0.3, ymin=0, ymax=0, position = position_dodge(width=0.75))+
  stat_summary(fun.y=q95, color = "gray47", geom="crossbar", size=0.3, ymin=0, ymax=0, position = position_dodge(width=0.75))+
  stat_summary(fun.y=q5, color = "gray47", geom="crossbar", size=0.3, ymin=0, ymax=0, position = position_dodge(width=0.75))+
#  scale_fill_manual(values= c("#000043", "#0033FF", "#01CCA4", "#BAFF12", "#FFCC00", "#FF3300","#000043")) +
#  scale_color_manual(values= c("#000043", "#0033FF", "#01CCA4", "#BAFF12", "#FFCC00", "#FF3300","#000043")) +
 coord_cartesian(ylim = c(-100, 300))+
  xlab(expression(""))+
  ylab(expression("Leverage (%)"))
leverbox
ggsave("figures/leveragebox_main_zoom.png", plot = leverbox, width = 8, height = 6, dpi=300)

leverbox.ext =ggplot(lev[lev$variable %in% ext,], aes(x=variable, y=leverage, fill = Category, color=Category)) + 
  geom_hline(yintercept=0)+
  geom_boxplot( alpha= 0.8, position = position_dodge(width=0.75))+
  stat_summary(fun.y=mean, color = "black", geom="crossbar", size=0.3, ymin=0, ymax=0, position = position_dodge(width=0.75))+
  stat_summary(fun.y=q95, color = "gray47", geom="crossbar", size=0.3, ymin=0, ymax=0, position = position_dodge(width=0.75))+
  stat_summary(fun.y=q5, color = "gray47", geom="crossbar", size=0.3, ymin=0, ymax=0, position = position_dodge(width=0.75))+
 # scale_fill_manual(values= c("#000043", "#0033FF", "#01CCA4", "#BAFF12", "#FFCC00", "#FF3300","#000043")) +
#  scale_color_manual(values= c("#000043", "#0033FF", "#01CCA4", "#BAFF12", "#FFCC00", "#FF3300","#000043")) +
  coord_cartesian(ylim = c(-100, 300))+
  xlab(expression(""))+
  ylab(expression("Leverage (%)"))
leverbox.ext
ggsave("figures/leveragebox_ext_zoom.png", plot = leverbox.ext, width = 6, height = 6, dpi=300)


##### linear models #####
PO4.Plmodel=select(subset(stoich,Category!="Utah Lake"), Event, Category, PO4.P, X.DEV,X.IMP, Forest, X.Herb)
PO4.Plmodel=PO4.Plmodel%>%
  drop_na()
PO4.Pmodel=gls(PO4.P~X.DEV+X.IMP+Forest+X.Herb+Event ,data=PO4.Plmodel)
model=  dredge(PO4.Pmodel)
print(model)
PO4.P=get.models(model,  1)[[1]] ####best model
summary(PO4.P)
print(cor(fitted(PO4.P), PO4.Plmodel$PO4.P)^2)
hist(PO4.P$residuals)

DIN.Nlmodel=select(subset(stoich,Category!="Utah Lake"), Event, Category, DIN.N, X.DEV,X.IMP, Forest, X.Herb)
DIN.Nlmodel=DIN.Nlmodel%>%
  drop_na()
DIN.Nmodel=gls(DIN.N~X.DEV+X.IMP+Forest+X.Herb+Event ,data=DIN.Nlmodel)
model=  dredge(DIN.Nmodel)
print(model)
DIN.N=get.models(model,  1)[[1]] ####best model
summary(DIN.N)
print(cor(fitted(DIN.N), DIN.Nlmodel$DIN.N)^2)
hist(DIN.N$residuals)

TNlmodel=select(subset(stoich,Category!="Utah Lake"), Event, Category, TN, X.DEV,X.IMP, Forest, X.Herb)
TNlmodel=TNlmodel%>%
  drop_na()
TNmodel=gls(TN~X.DEV+X.IMP+Forest+X.Herb+Event ,data=TNlmodel)
model=  dredge(TNmodel)
print(model)
TN=get.models(model,  1)[[1]] ####best model
summary(TN)
print(cor(fitted(TN), TNlmodel$TN)^2)
hist(TN$residuals)

DOClmodel=select(subset(stoich,Category!="Utah Lake"), Event, Category, DOC, X.DEV,X.IMP, Forest, X.Herb)
DOClmodel=DOClmodel%>%
  drop_na()
DOCmodel=gls(DOC~X.DEV+X.IMP+Forest+X.Herb+Event ,data=DOClmodel)
model=  dredge(DOCmodel)
print(model)
DOC=get.models(model,  1)[[1]] ####best model
summary(DOC)
print(cor(fitted(DOC), DOClmodel$DOC)^2)
hist(DOC$residuals)

Chloridelmodel=select(subset(stoich,Category!="Utah Lake"), Event, Category, Chloride, X.DEV,X.IMP, Forest, X.Herb)
Chloridelmodel=Chloridelmodel%>%
  drop_na()
Chloridemodel=gls(Chloride~X.IMP+X.Herb ,data=Chloridelmodel)
model=  dredge(Chloridemodel)
print(model)
Chloride=get.models(model,  1)[[1]] ####best model
summary(Chloride)
print(cor(fitted(Chloride), Chloridelmodel$Chloride)^2)
hist(Chloride$residuals)

Sulfatelmodel=select(subset(stoich,Category!="Utah Lake"), Event, Category, Sulfate, X.DEV,X.IMP, Forest, X.Herb)
Sulfatelmodel=Sulfatelmodel%>%
  drop_na()
Sulfatemodel=gls(Sulfate~X.IMP+X.Herb+Event ,data=Sulfatelmodel)
model=  dredge(Sulfatemodel)
print(model)
Sulfate=get.models(model,  1)[[1]] ####best model
summary(Sulfate)
print(cor(fitted(Sulfate), Sulfatelmodel$Sulfate)^2)
hist(Sulfate$residuals)
############ TEMPORAL VARIANCE ########

co.var<-function(x) 100*(sd(x,na.rm=TRUE)/mean(x,na.rm=TRUE))# create function

db.cv = dcast(ssmelt, Category + Area_SS + SiteID ~ variable, co.var) 

#melt and stack again to calculate the mean for each parameter
db.cv.m=melt(db.cv, id.vars=c("Area_SS", "Category", "SiteID"), measure.vars=4:14, na.rm=FALSE)

db.cvmean = dcast(db.cv.m, Category ~ variable, mean) 

write.csv(db.cv, file="output/dbOverallCV.csv")
write.csv(db.cvmean, file="output/meanCV.csv")

#### Tukey label graph ####
# I need to group the treatments that are not different each other together.
generate_label_df <- function(TUKEY, variable){
  
  # Extract labels and factor levels from Tukey post-hoc 
  Tukey.levels <- TUKEY[[variable]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  
  #I need to put the labels in the same order as in the boxplot :
  Tukey.labels$treatment=rownames(Tukey.labels)
  Tukey.labels=Tukey.labels[order(Tukey.labels$treatment) , ]
  return(Tukey.labels)
}

# Apply the function on my dataset
LABELS <- generate_label_df(TUKEY , "data$treatment")


# A panel of colors to draw each group with the same color :
my_colors <- c( 
  rgb(143,199,74,maxColorValue = 255),
  rgb(242,104,34,maxColorValue = 255), 
  rgb(111,145,202,maxColorValue = 255)
)

# Draw the basic boxplot
a <- boxplot(data$value ~ data$treatment , ylim=c(min(data$value) , 1.1*max(data$value)) , col=my_colors[as.numeric(LABELS[,1])] , ylab="value" , main="")

# I want to write the letter over each box. Over is how high I want to write it.
over <- 0.1*max( a$stats[nrow(a$stats),] )

#Add the labels
text( c(1:nlevels(data$treatment)) , a$stats[nrow(a$stats),]+over , LABELS[,1]  , col=my_colors[as.numeric(LABELS[,1])] )
#########Broken? #########
dbt <- as.data.table(ssmelt)
dbtWM <- dbt[,lapply(.SD, weighted.mean, w=q),by=list(Category, SiteID)]
#add column called qquant with 0 in it
dbtWM$qquant = 0
str(dbt)
write.csv(dbtWM, file="output/dbOverallWMean.csv")

WM<-read.csv('output/dbOverallCV.csv')
str(WM)
colnames(WM)


#create catchment and quantile dataframes
pr=subset(WM, Watershed == "Provo River")
pr=droplevels(pr)
sf=subset(WM, Watershed == "Spanish Fork")
sf=droplevels(sf)
hc=subset(WM, Watershed == "Hobble Creek")
hc=droplevels(hc)
af=subset(WM, Watershed == "American Fork")
af=droplevels(af)
bs=subset(WM, Watershed == "Benjamin Slough")
bs=droplevels(bs)
o=subset(WM, Watershed == "Other")
o=droplevels(o)

#calculate breakpoints in WM for hr
WMav = apply(WM[, 5:17], 2, function(x) cpts(cpt.var(x, method="PELT")))

#calculate breakpoints in WM for hc
hcav = apply(hc[, 5:17], 2, function(x) cpts(cpt.var(x, method="PELT")))

#calculate breakpoints in quantiles
#importing quantile dataset
QM<-read.csv('dbQuantileMean.csv')
str(QM)
colnames(QM)

#create catchment and quantile dataframes
hr1=subset(QM, catchment == "HR" & qquant == "1")
hr2=subset(QM, catchment == "HR" & qquant == "2")
hr3=subset(QM, catchment == "HR" & qquant == "3")
hr4=subset(QM, catchment == "HR" & qquant == "4")
hc1=subset(QM, catchment == "HC" & qquant == "1")
hc2=subset(QM, catchment == "HC" & qquant == "2")
hc3=subset(QM, catchment == "HC" & qquant == "3")
hc4=subset(QM, catchment == "HC" & qquant == "4")

#determine the breakpoints
hrq1 = apply(hr1[, 6:16], 2, function(x) cpts(cpt.var(x, method="PELT")))
hrq2 = apply(hr2[, 6:16], 2, function(x) cpts(cpt.var(x, method="PELT")))
hrq3 = apply(hr3[, 6:16], 2, function(x) cpts(cpt.var(x, method="PELT")))
hrq4 = apply(hr4[, 6:16], 2, function(x) cpts(cpt.var(x, method="PELT")))

hcq1 = apply(hc1[, 6:16], 2, function(x) cpts(cpt.var(x, method="PELT")))
hcq2 = apply(hc2[, 6:16], 2, function(x) cpts(cpt.var(x, method="PELT")))
hcq3 = apply(hc3[, 6:16], 2, function(x) cpts(cpt.var(x, method="PELT")))
hcq4 = apply(hc4[, 6:16], 2, function(x) cpts(cpt.var(x, method="PELT")))
#compile the vectors
breakpoints = rbind(hrav, hcav, hrq1, hrq2, hrq3, hrq4, hcq1, hcq2, hcq3, hcq4)
write.csv(breakpoints, file="breakpoints.csv")



