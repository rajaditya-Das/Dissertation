setwd("C:/Dissertation/Data/Traits/")

library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)
library(LeafArea)
library(ggpmisc)
library(patchwork)
library(ggpubr)

install.packages("gtsummary")
library(gtsummary)

#-------------------------------
#FUNCTIONS FOR PLOTTING
#----------------------------------
ggplotRegression <- function (fit) {
  
  require(ggplot2)
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point(col = "gray7") +
    stat_smooth(method = "lm", col = "turquoise4") +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 2),
                       " P =",signif(summary(fit)$coef[2,4], 2)))
}

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#-----------------------------------------------------------------------------------
#function cannot process all leaf scans at once, needs images to be sub divided into folders
# 
#LAI.df1 <- run.ij(set.directory = "H:/Consolidated/F1", path.imagej = 'H:/ImageJ')
#LAI.df2 <- run.ij(set.directory = "H:/Consolidated/F2", path.imagej = 'H:/ImageJ')
#LAI.df3 <- run.ij(set.directory = "H:/Consolidated/F3", path.imagej = 'H:/ImageJ')
#LAI.df4 <- run.ij(set.directory = "H:/Consolidated/F4", path.imagej = 'H:/ImageJ')
#LAI.df5 <- run.ij(set.directory = "H:/Consolidated/F5", path.imagej = 'H:/ImageJ')
#LAI.df6 <- run.ij(set.directory = "H:/Consolidated/F6", path.imagej = 'H:/ImageJ')
#LAI.df7 <- run.ij(set.directory = "H:/Consolidated/F7", path.imagej = 'H:/ImageJ')
#LAI.df8 <- run.ij(set.directory = "H:/Consolidated/F8", path.imagej = 'H:/ImageJ')
#LAI.df9 <- run.ij(set.directory = "H:/Consolidated/F9", path.imagej = 'H:/ImageJ')
# 
# LAI.df <- Reduce(rbind, list(LAI.df1,LAI.df2,LAI.df3,LAI.df4,LAI.df5,LAI.df6,LAI.df7,LAI.df8, LAI.df9))
# write.csv(LAI.df, file = "LAI.csv")
# 'LAI.csv' file cleaned and saved as "Leaf Area Index.csv" file
#--------------------------------------------------------------------------------

dat <- read.csv("ms_data with dry.csv") ##considering only for the species which are above 5 individuals and has dry weights
LAI <- read.csv("Leaf Area Index.csv")
envt_data <- read.csv("PlotXwater.csv")
light_data <- read.csv("Light_data.csv")
dat1 <- left_join(dat, light_data)
MS_data <- left_join(left_join(dat1, envt_data), LAI)
str(MS_data)
MS_data <- distinct(MS_data)

##----------------------------------------------------------------------------------
#coerce to numeric; new colums for derived traits - RMF, LMF, SMF, SSL, SRL, SLA, LDMC

MS_data$Leaf.wet.mass <- as.numeric(as.character(MS_data$Leaf.wet.mass))
MS_data$Above.ground.mass <- as.numeric(as.character(MS_data$Above.ground.mass))
MS_data$Root.dry.mass <- as.numeric(as.character(MS_data$Root.dry.mass))
MS_data$Stem.dry.mass <-  as.numeric(as.character(MS_data$Stem.dry.mass))

MS_data$Stem.length <- as.numeric(as.character(MS_data$Stem.length))
MS_data$Root.length <- as.numeric(as.character(MS_data$Root.length))
MS_data$Leaf1 <- as.numeric(as.character(MS_data$Leaf1))
MS_data$Leaf2 <- as.numeric(as.character(MS_data$Leaf2))

MS_data$Leaf.other <- as.numeric(as.character(MS_data$Leaf.other))

#-------
MS_data$total.mass <- MS_data$Above.ground.mass + MS_data$Root.mass

MS_data <- MS_data %>% 
  rowwise() %>% 
  mutate(Leaf.dry.mass = sum(Leaf1,Leaf2, Leaf.other, na.rm = TRUE))

MS_data <- MS_data %>% 
  rowwise() %>% 
  mutate(total.dry.mass = sum(Leaf.dry.mass, Root.dry.mass, Stem.dry.mass))


MS_data$Leaf.mature <- sum(MS_data$Leaf1,MS_data$Leaf2, na.rm=T)
MS_data$Leaf.dry.matter.content = MS_data$Leaf.dry.mass / MS_data$Leaf.wet.mass

MS_data$SLA = MS_data$total.leaf.area / MS_data$Leaf.mature
MS_data$SLA1 = MS_data$total.leaf.area / MS_data$Leaf.wet.mass
MS_data[sapply(MS_data, is.infinite)] <- NA



MS_data <- mutate(MS_data, Leaf.mass.fraction1 = ((Above.ground.mass- Stem.mass) /total.mass))
MS_data <- mutate(MS_data, Leaf.mass.fraction = (Leaf.dry.mass)/total.mass)

MS_data <-mutate(MS_data, Root.mass.fraction1 = Root.mass/total.mass)
MS_data <-mutate(MS_data, Root.mass.fraction = Root.dry.mass/total.dry.mass)

MS_data <- mutate(MS_data, Stem.mass.fraction1 = Stem.mass/total.mass)
MS_data <- mutate(MS_data, Stem.mass.fraction = Stem.dry.mass/total.dry.mass)

MS_data <- mutate(MS_data, Root.specific.length = Root.length/Root.dry.mass)
MS_data <- mutate(MS_data, Stem.specific.length = Stem.length/Stem.dry.mass)
MS_data <- mutate(MS_data, Root.specific.length1 = Root.length/Root.mass)
MS_data <- mutate(MS_data, Stem.specific.length1 = Stem.length/Stem.mass)

MS_data$VWC = (MS_data$VWC1+MS_data$VWC2+MS_data$VWC3)/3 #Mean VWC at each plot


#-------------------------
#shapiro.test((MS_data$Stem.mass.fraction))
#ggdensity(MS_data$Stem.mass.fraction)


#------------------------------------------------------------------------
#Bootstrapping : sampling 999 times to get mean CV 
#-------------------------------------------------------------------------

Sampled.SD <- MS_data %>%
  group_by(Species.ID)

x <- group_split(Sampled.SD)
group_keys(Sampled.SD)

#---------------------------------------
F= list()
S1 <- list()

for(k in 1:14){                              # for each species
  lmf <- x[[k]] %>% pull(Leaf.mass.fraction1)
  smf <- x[[k]] %>% pull(Stem.mass.fraction)
  rmf <- x[[k]] %>% pull(Root.mass.fraction)
  ssl <- x[[k]] %>% pull(Stem.specific.length)
  rsl <- x[[k]] %>% pull(Root.specific.length)
  sla <- x[[k]] %>% pull(SLA)
  lai <- x[[k]] %>% pull(total.leaf.area)
  ldmc <- x[[k]] %>% pull(Leaf.dry.matter.content)
  
  for(i in 1:9999){                          # 9999 times 
    y <-sample(lmf, 4, replace = T)
    z <-sample(smf, 4, replace = T)
    w <- sample(rmf, 4, replace = T)
    v <- sample(ssl, 4, replace = T)
    u <- sample(rsl, 4, replace = T)
    t <- sample(sla, 4, replace = T)
    s <- sample(lai, 4, replace = T)
    q <- sample(ldmc, 4, replace = T)
    
    S1$ISTV1[i]= sd(y, na.rm= T)/mean(y, na.rm= T) # sd/mean = CV of 9999 samples
    S1$ISTV2[i]= sd(z, na.rm= T)/mean(z, na.rm= T)
    S1$ISTV3[i]= sd(w, na.rm= T)/mean(w, na.rm= T)
    S1$ISTV4[i]= sd(v, na.rm= T)/mean(v, na.rm= T)
    S1$ISTV5[i]= sd(u, na.rm= T)/mean(u, na.rm= T)
    S1$ISTV6[i]= sd(t, na.rm= T)/mean(t, na.rm= T)
    S1$ISTV7[i]= sd(s, na.rm= T)/mean(s, na.rm= T)
    S1$ISTV8[i]= sd(q, na.rm= T)/mean(q, na.rm= T)
  }
  F[k] = list(S1)
}

IST <- data.frame(matrix(ncol= 8, nrow = 14), "Species.ID" = c('Cin','ClaDen', 'Comp', 'GarMor', 'Litsea', 'morph1', 'morph3', 'morph4', 'morph2', 'PsyMac', 'PsyNig', 'SymRac', 'Syz', 'wavy'))
row <- list(c('Cin','ClaDen', 'Comp', 'GarMor', 'Litsea', 'morph1', 'morph3', 'morph4', 'morph2', 'PsyMac', 'PsyNig', 'SymRac', 'Syz', 'wavy'))
col <- c('ist1', 'ist2', 'ist3', 'ist4' , 'ist5', 'ist6', 'ist7', 'ist8', "Species.ID")
colnames(IST) <- col
for(i in 1:14){
  IST[i,1] <- mean(unlist(F[[i]][1]), na.rm= T)
  IST[i,2] <- mean(unlist(F[[i]][2]), na.rm= T)
  IST[i,3] <- mean(unlist(F[[i]][3]), na.rm= T)
  IST[i,4] <- mean(unlist(F[[i]][4]), na.rm= T)
  IST[i,5] <- mean(unlist(F[[i]][5]), na.rm= T)
  IST[i,6] <- mean(unlist(F[[i]][6]), na.rm= T)
  IST[i,7] <- mean(unlist(F[[i]][7]), na.rm= T)
  IST[i,8] <- mean(unlist(F[[i]][8]), na.rm= T)
}

#-------------------------------
#Functional position values 

G <-  list()
H <- list()
R1 <- list()
T1 <- list()
cwm <- list()
for(k in 1:14){                              # for each species
  lmf <- x[[k]] %>% pull(Leaf.mass.fraction)
  smf <- x[[k]] %>% pull(Stem.mass.fraction)
  rmf <- x[[k]] %>% pull(Root.mass.fraction)
  ssl <- x[[k]] %>% pull(Stem.specific.length)
  rsl <- x[[k]] %>% pull(Root.specific.length)
  sla <- x[[k]] %>% pull(SLA)
  lai <- x[[k]] %>% pull(total.leaf.area)
  ldmc <- x[[k]] %>% pull(Leaf.dry.matter.content)
  
  cwm[1] <- (mean(MS_data$Leaf.mass.fraction, na.rm = T))
  cwm[2] <- (mean(MS_data$Stem.mass.fraction, na.rm = T))
  cwm[3] <- (mean(MS_data$Root.mass.fraction, na.rm = T))
  cwm[4] <- (mean(MS_data$Stem.specific.length, na.rm = T))
  cwm[5] <- (mean(MS_data$Root.specific.length, na.rm = T))
  cwm[6] <- (mean(MS_data$SLA, na.rm = T))
  cwm[7] <- (mean(MS_data$total.leaf.area, na.rm = T))
  cwm[8] <- (mean(MS_data$Leaf.dry.matter.content, na.rm = T))
  
  for(i in 1:9999){                          # 9999 times 
    y <-sample(lmf, 4, replace = T)
    z <-sample(smf, 4, replace = T)
    w <- sample(rmf, 4, replace = T)
    v <- sample(ssl, 4, replace = T)
    u <- sample(rsl, 4, replace = T)
    t <- sample(sla, 4, replace = T)
    s <- sample(lai, 4, replace = T)
    q <- sample(ldmc, 4, replace = T)
    
    R1$fposx1[i]= (mean(y, na.rm= T) - cwm[[1]]) 
    R1$fposx2[i]= (mean(z, na.rm= T) - cwm[[2]])
    R1$fposx3[i]= (mean(w, na.rm= T) - cwm[[3]])
    R1$fposx4[i]= (mean(v, na.rm= T) - cwm[[4]])
    R1$fposx5[i]= (mean(u, na.rm= T) - cwm[[5]])
    R1$fposx6[i]= (mean(t, na.rm= T) - cwm[[6]])
    R1$fposx7[i]= (mean(s, na.rm= T) - cwm[[7]])
    R1$fposx8[i]= (mean(q, na.rm= T) - cwm[[8]])
    
    T1$fpos1[i]= abs(mean(y, na.rm= T) - cwm[[1]]) 
    T1$fpos2[i]= abs(mean(z, na.rm= T) - cwm[[2]])
    T1$fpos3[i]= abs(mean(w, na.rm= T) - cwm[[3]])
    T1$fpos4[i]= abs(mean(v, na.rm= T) - cwm[[4]])
    T1$fpos5[i]= abs(mean(u, na.rm= T) - cwm[[5]])
    T1$fpos6[i]= abs(mean(t, na.rm= T) - cwm[[6]])
    T1$fpos7[i]= abs(mean(s, na.rm= T) - cwm[[7]])
    T1$fpos8[i]= abs(mean(q, na.rm= T) - cwm[[8]])
    
  }
  G[k] = list(R1)
  H[k] = list(T1)
}

FPOSX <- data.frame(matrix(ncol= 8, nrow = 14), "Species.ID" =c('Cin','ClaDen', 'Comp', 'GarMor', 'Litsea', 'morph1', 'morph3', 'morph4', 'morph2', 'PsyMac', 'PsyNig', 'SymRac', 'Syz', 'wavy'))
row <- list(c('Cin','ClaDen', 'Comp', 'GarMor', 'Litsea', 'morph1', 'morph3', 'morph4', 'morph2', 'PsyMac', 'PsyNig', 'SymRac', 'Syz', 'wavy'))
col <- c('fposx1', 'fposx2', 'fposx3', 'fposx4' , 'fposx5', 'fposx6', 'fposx7', 'fposx8', "Species.ID")
colnames(FPOSX) <- col

for(i in 1:14){
  FPOSX[i,1] <- mean(unlist(G[[i]][1]), na.rm= T)
  FPOSX[i,2] <- mean(unlist(G[[i]][2]), na.rm= T)
  FPOSX[i,3] <- mean(unlist(G[[i]][3]), na.rm= T)
  FPOSX[i,4] <- mean(unlist(G[[i]][4]), na.rm= T)
  FPOSX[i,5] <- mean(unlist(G[[i]][5]), na.rm= T)
  FPOSX[i,6] <- mean(unlist(G[[i]][6]), na.rm= T)
  FPOSX[i,7] <- mean(unlist(G[[i]][7]), na.rm= T)
  FPOSX[i,8] <- mean(unlist(G[[i]][8]), na.rm= T)
}

FPOS <- data.frame(matrix(ncol= 8, nrow = 14), "Species.ID" =c('Cin','ClaDen', 'Comp', 'GarMor', 'Litsea', 'morph1', 'morph3', 'morph4', 'morph2', 'PsyMac', 'PsyNig', 'SymRac', 'Syz', 'wavy'))
row <- list(c('Cin','ClaDen', 'Comp', 'GarMor', 'Litsea', 'morph1', 'morph3', 'morph4', 'morph2', 'PsyMac', 'PsyNig', 'SymRac', 'Syz', 'wavy'))
col <- c('fpos1', 'fpos2', 'fpos3', 'fpos4' , 'fpos5', 'fpos6', 'fpos7', 'fpos8', "Species.ID")
colnames(FPOS) <- col

for(i in 1:14){
  FPOS[i,1] <- mean(unlist(H[[i]][1]), na.rm= T)
  FPOS[i,2] <- mean(unlist(H[[i]][2]), na.rm= T)
  FPOS[i,3] <- mean(unlist(H[[i]][3]), na.rm= T)
  FPOS[i,4] <- mean(unlist(H[[i]][4]), na.rm= T)
  FPOS[i,5] <- mean(unlist(H[[i]][5]), na.rm= T)
  FPOS[i,6] <- mean(unlist(H[[i]][6]), na.rm= T)
  FPOS[i,7] <- mean(unlist(H[[i]][7]), na.rm= T)
  FPOS[i,8] <- mean(unlist(H[[i]][8]), na.rm= T)
}

#[3] Abundance = number of entries in Species(i)

abundance <- MS_data %>%
  group_by(Species.ID) %>%
  summarise(abundance = length(Species.ID))

#[4] Niche width(i) = standard deviation of VWC at plots with Species(i)

total.niche <- max( MS_data$VWC, na.rm= T)- min( MS_data$VWC, na.rm= T)

niche.width.water <- MS_data %>%
  group_by(Species.ID) %>%
  summarise(niche.width.water = max(VWC, na.rm= T)- min(VWC, na.rm= T))

niche.width.ligth <- MS_data %>%
  group_by(Species.ID) %>%
  summarise(niche.width.light = max(Canopy, na.rm= T)- min(Canopy, na.rm= T))

niche.width <- left_join(niche.width.ligth, niche.width.water)

#[5] Niche Position =  mean VWC (plots of Species'i')- mean VWC(all sampled plots)

mean.VWC <- mean(MS_data$VWC, na.rm = T)
mean.canopy <- mean(MS_data$Canopy, na.rm = T)

niche.pos.water <- MS_data %>%
  group_by(Species.ID) %>%
  summarise(niche.pos.water = abs(mean(VWC, na.rm= T) - mean.VWC))

niche.pos.light <- MS_data %>%
  group_by(Species.ID) %>%
  summarise(niche.pos.light = abs(mean(Canopy, na.rm= T) - mean.canopy))

#-----------------------------------
#alternate defn of niche position without taking absolute value

niche.pos.water1 <- MS_data %>%
  group_by(Species.ID) %>%
  summarise(niche.pos.water1 = (mean(VWC, na.rm= T) - mean.VWC))

niche.pos.light1 <- MS_data %>%
  group_by(Species.ID) %>%
  summarise(niche.pos.light1 = (mean(Canopy, na.rm= T) - mean.canopy))

#------

niche.pos <- left_join(niche.pos.light, niche.pos.water)

Met.dat <- data.frame(Reduce(left_join, list(IST, FPOS, FPOSX, abundance, niche.width, niche.pos)))


###-----------#--------------------------------------------------
#     A       # Abundance X ISTV
##------------#-------------------------------------------------

A1 <- lm(abundance ~ ist1, data = Met.dat)
(summary(A1))
PA1 <- ggplotRegression(lm(abundance ~ ist1, data = Met.dat))+
  xlab(c('ISTV(Leaf mass fraction)'))+ ylab(c('Abundance'))

A2 <- lm(abundance ~ ist2, data = Met.dat)
summary(A2)
PA2 <- ggplotRegression(lm(abundance ~ ist2, data = Met.dat))+
  xlab(c('ISTV(Stem mass fraction)'))+ ylab(c('Abundance'))

A3 <-  lm(abundance ~ ist3, data = Met.dat)
summary(A3)
PA3 <- ggplotRegression(lm(abundance ~ ist3, data = Met.dat))+
  xlab(c('ISTV(Root mass fraction)'))+ ylab(c('Abundance'))

A4 <-  lm(abundance ~ ist4, data = Met.dat)
summary(A4)
PA4 <- ggplotRegression(lm(abundance ~ ist4, data = Met.dat))+
  xlab(c('ISTV(Specific stem length)'))+ ylab(c('Abundance'))

A5 <-  lm(abundance ~ ist5, data = Met.dat)
summary(A5)
PA5 <- ggplotRegression(lm(abundance ~ ist5, data = Met.dat))+
  xlab(c('ISTV(Specific root length)'))+ ylab(c('Abundance'))

A6 <-  lm(abundance ~ ist6, data = Met.dat)
summary(A6)
PA6 <- ggplotRegression(lm(abundance ~ ist6, data = Met.dat))+
  xlab(c('ISTV(SLA)'))+ ylab(c('Abundance'))

A7 <-  lm(abundance ~ ist7, data = Met.dat)
summary(A7)
PA7 <- ggplotRegression(lm(abundance ~ ist7, data = Met.dat))+
  xlab(c('ISTV(Leaf area index)'))+ ylab(c('Abundance'))


A8 <- lm(abundance ~ ist8, data = Met.dat)
summary(A8)
PA8 <- ggplotRegression(lm(abundance ~ ist8, data = Met.dat))+
  xlab(c('ISTV(Leaf dry matter content)'))+ ylab(c('Abundance'))

multiplot(PA1, PA2, PA3, PA4, PA5, PA6, PA7, PA8, cols = 2)

# aa <- glance(summary(A1))
# bb <- glance(summary(A2))
# cc <- glance(summary(A3))
# dd <- glance(summary(A4))
# ee <- glance(summary(A5))
# ff <- glance(summary(A6))
# gg <- glance(summary(A7))
# hh <- glance(summary(A8))
# 
# A <- rbind(aa,bb,cc,dd,ee,ff,gg,hh)[,c(1:5)]

#-----------#---------------------------------------------------
#     B     #   Niche.width X ISTV : WATER
#-----------#----------------------------------------------------

B1 <- lm(niche.width.water ~ ist1, dat = Met.dat)
summary(B1)
PB1 <- ggplotRegression(lm(niche.width.water ~ ist1, data = Met.dat))+
  xlab(c('ISTV(Leaf mass fraction)'))

B2 <- lm(niche.width.water ~ ist2, data = Met.dat)
summary(B2)
PB2 <- ggplotRegression(lm(niche.width.water ~ ist2, data = Met.dat))+
  xlab(c('ISTV(Stem mass fraction)'))

B3 <- lm(niche.width.water ~ ist3, data = Met.dat)
summary(B3)
PB3 <- ggplotRegression(lm(niche.width.water ~ ist3, data = Met.dat))+
  xlab(c('ISTV(Root mass fraction)'))


B4 <- lm(niche.width.water ~ ist4, data = Met.dat)
summary(B4)
PB4 <- ggplotRegression(lm(niche.width.water ~ ist4, data = Met.dat))+
  xlab(c('ISTV(Specific stem length)'))


B5 <- lm(niche.width.water ~ ist5, data = Met.dat)
summary(B5)
PB5 <- ggplotRegression(lm(niche.width.water ~ ist5, data = Met.dat))+
  xlab(c('ISTV(Specific root length)'))


B6 <- lm(niche.width.water ~ ist6, data = Met.dat)
summary(B6)
PB6 <- ggplotRegression(lm(niche.width.water ~ ist6, data = Met.dat))+
  xlab(c('ISTV(SLA)'))


B7 <- lm(niche.width.water ~ ist7, data = Met.dat)
summary(B7)
PB7 <- ggplotRegression(lm(niche.width.water ~ ist7, data = Met.dat))+
  xlab(c('ISTV(Leaf area index)'))

B8 <- lm(niche.width.water ~ ist8, data = Met.dat)
summary(B8)
PB8 <- ggplotRegression(lm(niche.width.water ~ ist8, data = Met.dat))+
  xlab(c('ISTV(Leaf dry matter content)'))

multiplot(PB1, PB2, PB3, PB4,PB5, PB6, PB7, PB8, cols = 2)

#-----------#---------------------------------------------------
#     B- light    #   Niche.width X ISTV : LIGHT
#-----------#----------------------------------------------------

BL1 <- lm(niche.width.light ~ ist1, dat = Met.dat)
summary(BL1)
PB1 <- ggplotRegression(lm(niche.width.light ~ ist1, data = Met.dat))+
  xlab(c('ISTV(Leaf mass fraction)'))

BL2 <- lm(niche.width.light ~ ist2, data = Met.dat)
summary(BL2)
PB2 <- ggplotRegression(lm(niche.width.light ~ ist2, data = Met.dat))+
  xlab(c('ISTV(Stem mass fraction)'))

BL3 <- lm(niche.width.light ~ ist3, data = Met.dat)
summary(BL3)
PB3 <- ggplotRegression(lm(niche.width.light ~ ist3, data = Met.dat))+
  xlab(c('ISTV(Root mass fraction)'))


BL4 <- lm(niche.width.light~ ist4, data = Met.dat)
summary(BL4)
PB4 <- ggplotRegression(lm(niche.width.light ~ ist4, data = Met.dat))+
  xlab(c('ISTV(Specific stem length)'))


BL5 <- lm(niche.width.light ~ ist5, data = Met.dat)
summary(BL5)
PB5 <- ggplotRegression(lm(niche.width.light ~ ist5, data = Met.dat))+
  xlab(c('ISTV(Specific root length)'))


BL6 <- lm(niche.width.light ~ ist6, data = Met.dat)
summary(BL6)
PB6 <- ggplotRegression(lm(niche.width.light ~ ist6, data = Met.dat))+
  xlab(c('ISTV(SLA)'))


BL7 <- lm(niche.width.light ~ ist7, data = Met.dat)
summary(BL7)
PB7 <- ggplotRegression(lm(niche.width.light~ ist7, data = Met.dat))+
  xlab(c('ISTV(Leaf area index)'))

B8 <- lm(niche.width.light ~ ist8, data = Met.dat)
summary(B8)
PB8 <- ggplotRegression(lm(niche.width.light ~ ist8, data = Met.dat))+
  xlab(c('ISTV(Leaf dry matter content)'))

multiplot(PB1, PB2, PB3, PB4,PB5, PB6, PB7, PB8, cols = 2)

#-----------#---------------------------------------------------
#     c     #     functional position x istv | Functional position is the absolute value of the difference from the mean community functional trait
#-----------#---------------------------------------------------

C1 <- lm(fpos1 ~ ist1, data = Met.dat)
summary(C1)
PC1 <- ggplotRegression(lm(fpos1 ~ ist1, data = Met.dat))+
  xlab(c('ISTV(Leaf mass fraction)'))

C2 <- lm(fpos2 ~ ist2, data = Met.dat)
summary(C2)
PC2 <- ggplotRegression(lm(fpos2 ~ ist2, data = Met.dat))+
  xlab(c('ISTV(Stem mass fraction)')) +ylab(c('FPOS(Stem mass fraction)'))

C3 <- lm(fpos3 ~ ist3, data = Met.dat)
summary(C3)
PC3 <- ggplotRegression(lm(fpos3 ~ ist3, data = Met.dat))+
  xlab(c('ISTV(Root mass fraction)'))


C4 <- lm(fpos4 ~ ist4, data = Met.dat)
summary(C4)
PC4 <- ggplotRegression(lm(fpos4 ~ ist4, data = Met.dat))+
  xlab(c('ISTV(Specific stem length)'))+ylab(c('FPOS(Specific stem length)'))


C5 <- lm(fpos5 ~ ist5, data = Met.dat)
summary(C5)
PC5 <- ggplotRegression(lm(fpos5 ~ ist5, data = Met.dat))+
  xlab(c('ISTV(Specific root length)')) +ylab(c('FPOS(Specific root length)'))


C6 <- lm(fpos6 ~ ist6, data = Met.dat)
summary(C6)
PC6 <- ggplotRegression(lm(fpos6 ~ ist6, data = Met.dat))+
  xlab(c('ISTV(Specific leaf area)'))+ ylab(c('FPOS(Specific leaf area)'))


C7 <- lm(fpos7 ~ ist7, data = Met.dat)
summary(C7)
PC7 <- ggplotRegression(lm(fpos7 ~ ist7, data = Met.dat))+
  xlab(c('ISTV(Leaf area index)')) +ylab(c('FPOS(Leaf area index)'))

C8 <- lm(fpos8 ~ ist8, data = Met.dat)
summary(C8)
PC8 <- ggplotRegression(lm(fpos8 ~ ist8, data = Met.dat))+
  xlab(c('ISTV(Leaf dry matter content)'))


multiplot(PC1, PC2, PC3, PC4,PC5, PC6, PC7, PC8, cols = 2)

# aa <- glance(summary(C1))
# bb <- glance(summary(C2))
# cc <- glance(summary(C3))
# dd <- glance(summary(C4))
# ee <- glance(summary(C5))
# ff <- glance(summary(C6))
# gg <- glance(summary(C7))
# hh <- glance(summary(C8))
# 
# D<- rbind(aa,bb,cc,dd,ee,ff,gg,hh)[,c(1:5)]


#-----------#---------------------------------------------------
#     c[1]     #     functional position x istv | alternate defn of functional position used here. Absolute value not taken
#-----------#---------------------------------------------------

C1 <- lm(fposx1 ~ ist1, data = Met.dat)
summary(C1)
PC1 <- ggplotRegression(lm(fposx1 ~ ist1, data = Met.dat))+
  xlab(c('ISTV(Leaf mass fraction)'))

C2 <- lm(fposx2 ~ ist2, data = Met.dat)
summary(C2)
PC2 <- ggplotRegression(lm(fposx2 ~ ist2, data = Met.dat))+
  xlab(c('ISTV(Stem mass fraction)'))

C3 <- lm(fposx3 ~ ist3, data = Met.dat)
summary(C3)
PC3 <- ggplotRegression(lm(fposx3 ~ ist3, data = Met.dat))+
  xlab(c('ISTV(Root mass fraction)'))


C4 <- lm(fposx4 ~ ist4, data = Met.dat)
summary(C4)
PC4 <- ggplotRegression(lm(fposx4 ~ ist4, data = Met.dat))+
  xlab(c('ISTV(Specific stem length)'))


C5 <- lm(fposx5 ~ ist5, data = Met.dat)
summary(C5)
PC5 <- ggplotRegression(lm(fposx5 ~ ist5, data = Met.dat))+
  xlab(c('ISTV(Specific root length)'))


C6 <- lm(fposx6 ~ ist6, data = Met.dat)
summary(C6)
PC6 <- ggplotRegression(lm(fposx6 ~ ist6, data = Met.dat))+
  xlab(c('ISTV(SLA)'))


C7 <- lm(fposx7 ~ ist7, data = Met.dat)
summary(C7)
PC7 <- ggplotRegression(lm(fposx7 ~ ist7, data = Met.dat))+
  xlab(c('ISTV(Leaf area index)'))


C8 <- lm(fposx8 ~ ist8, data = Met.dat)
summary(C8)
PC8 <- ggplotRegression(lm(fposx8 ~ ist8, data = Met.dat))+
  xlab(c('ISTV(LDMC)'))

multiplot(PC1, PC2, PC3, PC4,PC5, PC6, PC7, PC8, cols = 2)



#-------#---------------------------------------------------------------
#   D   #       Niche position x ISTV : Light
#-------#---------------------------------------------------------------

D1 <- lm(niche.pos.light ~ ist1, data = Met.dat)
summary(D1)
PD1 <- ggplotRegression(lm(niche.pos.light ~ ist1, data = Met.dat))+
  xlab(c('ISTV(Leaf mass fraction)'))

D2 <- lm(niche.pos.light ~ ist2, data = Met.dat)
summary(D2)
PD2 <- ggplotRegression(lm(niche.pos.light ~ ist2, data = Met.dat))+
  xlab(c('ISTV(Stem mass fraction)'))

D3 <- lm(niche.pos.light ~ ist3, data = Met.dat)
summary(D3)
PD3 <- ggplotRegression(lm(niche.pos.light ~ ist3, data = Met.dat))+
  xlab(c('ISTV(Root mass fraction)'))

D4 <- lm(niche.pos.light ~ ist4, data = Met.dat)
summary(D4)
PD4 <- ggplotRegression(lm(niche.pos.light ~ ist4, data = Met.dat))+
  xlab(c('ISTV(Specific stem length)'))


D5 <- lm(niche.pos.light ~ ist5, data = Met.dat)
summary(D5)
PD5 <- ggplotRegression(lm(niche.pos.light ~ ist5, data = Met.dat))+
  xlab(c('ISTV(Specific root length)'))

D6 <- lm(niche.pos.light ~ ist6, data = Met.dat)
summary(D6)
PD6 <- ggplotRegression(lm(niche.pos.light ~ ist6, data = Met.dat))+
  xlab(c('ISTV(SLA)'))

D7 <- lm(niche.pos.light ~ ist7, data = Met.dat)
summary(D7)
PD7 <- ggplotRegression(lm(niche.pos.light ~ ist7, data = Met.dat))+
  xlab(c('ISTV(Leaf area index)'))

D8 <- lm(niche.pos.light ~ ist8, data = Met.dat)
summary(D8)
PD8 <- ggplotRegression(lm(niche.pos.light ~ ist8, data = Met.dat))+
  xlab(c('ISTV(Leaf dry matter content)'))

multiplot(PD1, PD2, PD3, PD4, PD5, PD6, PD7, PD8, cols = 2)

# aa <- glance(summary(D1))
# bb <- glance(summary(D2))
# cc <- glance(summary(D3))
# dd <- glance(summary(D4))
# ee <- glance(summary(D5))
# ff <- glance(summary(D6))
# gg <- glance(summary(D7))
# hh <- glance(summary(D8))
# 
# E <- rbind(aa,bb,cc,dd,ee,ff,gg,hh)[,c(1:5)]

#-------#---------------------------------------------------------------
#   E   #       Niche position x ISTV : Water
#-------#---------------------------------------------------------------

E1 <- lm(niche.pos.water ~ ist1, data = Met.dat)
summary(E1)
PE1 <- ggplotRegression(lm(niche.pos.water ~ ist1, data = Met.dat))+
  xlab(c('ISTV(Leaf mass fraction'))

E2 <- lm(niche.pos.water ~ ist2, data = Met.dat)
summary(E2)
PE2 <- ggplotRegression(lm(niche.pos.water ~ ist2, data = Met.dat))+
  xlab(c('ISTV(Stem mass fraction)'))

E3 <- lm(niche.pos.water ~ ist3, data = Met.dat)
summary(E3)
PE3 <- ggplotRegression(lm(niche.pos.water ~ ist3, data = Met.dat))+
  xlab(c('ISTV(Root mass fraction)'))

E4 <- lm(niche.pos.water ~ ist4, data = Met.dat)
summary(E4)
PE4 <- ggplotRegression(lm(niche.pos.water ~ ist4, data = Met.dat))+
  xlab(c('ISTV(Specific stem length)'))

E5 <- lm(niche.pos.water ~ ist5, data = Met.dat)
summary(E5)
PE5 <- ggplotRegression(lm(niche.pos.water ~ ist5, data = Met.dat))+
  xlab(c('ISTV(Specific root length)'))

E6 <- lm(niche.pos.water ~ ist6, data = Met.dat)
summary(E6)
PE6 <- ggplotRegression(lm(niche.pos.water ~ ist6, data = Met.dat))+
  xlab(c('ISTV(SLA)'))

E7 <- lm(niche.pos.water ~ ist7, data = Met.dat)
summary(E7)
PE7 <- ggplotRegression(lm(niche.pos.water ~ ist7, data = Met.dat))+
  xlab(c('ISTV(Leaf area index)'))

E8 <- lm(niche.pos.water ~ ist8, data = Met.dat)
summary(E8)
PE8 <- ggplotRegression(lm(niche.pos.water ~ ist8, data = Met.dat))+
  xlab(c('ISTV(Leaf area index)'))

multiplot(PE1, PE2, PE3, PE4, PE5, PE6, PE7, PE8, cols = 2)


# TO output tables
# aa <- glance(summary(D1))
# bb <- glance(summary(D2))
# cc <- glance(summary(D3))
# dd <- glance(summary(D4))
# ee <- glance(summary(D5))
# ff <- glance(summary(D6))
# gg <- glance(summary(D7))
# hh <- glance(summary(D8))
# 
# FFf <- rbind(aa,bb,cc,dd,ee,ff,gg,hh)[,c(1:5)]


#-------------------------------------------------------------------
# Niche position(light and water) x Functional Position
#-----------------------------------------------------------

D1 <- lm(niche.pos.light ~ fpos1, data = Met.dat)
summary(D1)
PD1 <- ggplotRegression(lm(niche.pos.light ~ fpos1, data = Met.dat))+
  xlab(c('ISTV(Leaf mass fraction'))

D2 <- lm(niche.pos.light ~ fpos2, data = Met.dat)
summary(D2)
PD2 <- ggplotRegression(lm(niche.pos.light ~ fpos2, data = Met.dat))+
  ylab(c('NPOS(Light)')) + xlab(c('FPOS(Stem mass fraction)'))

D3 <- lm(niche.pos.light ~ fpos3, data = Met.dat)
summary(D3)
PD3 <- ggplotRegression(lm(niche.pos.light ~ fpos3, data = Met.dat))+
  xlab(c('ISTV(Root mass fraction)'))

D4 <- lm(niche.pos.light ~ fpos4, data = Met.dat)
summary(D4)
PD4 <- ggplotRegression(lm(niche.pos.light ~ fpos4, data = Met.dat))+
  xlab(c('ISTV(Specific stem length)'))

D5 <- lm(niche.pos.light ~ fpos5, data = Met.dat)
summary(D5)
PD5 <- ggplotRegression(lm(niche.pos.light ~ fpos5, data = Met.dat))+
  xlab(c('ISTV(Specific root length)'))

D6 <- lm(niche.pos.light ~ fpos6, data = Met.dat)
summary(D6)
PD6 <- ggplotRegression(lm(niche.pos.light ~ fpos6, data = Met.dat))+
  xlab(c('ISTV(SLA)'))

D7 <- lm(niche.pos.light ~ fpos7, data = Met.dat)
summary(D7)
PD7 <- ggplotRegression(lm(niche.pos.light ~ fpos7, data = Met.dat))+
  xlab(c('ISTV(Leaf area index)'))

D8 <- lm(niche.pos.light ~ fpos8, data = Met.dat)
summary(D8)
PD8 <- ggplotRegression(lm(niche.pos.light ~ fpos8, data = Met.dat))+
  xlab(c('ISTV(Ldmc)'))


multiplot(PD1, PD2, PD3, PD4, PD5, PD6, PD7, cols = 2)

#-----------------------------------------------------------

D1 <- lm(niche.pos.water ~ fpos1, data = Met.dat)
summary(D1)
PD1 <- ggplotRegression(lm(niche.pos.water ~ fpos1, data = Met.dat))+
  xlab(c('ISTV(Leaf mass fraction'))

D2 <- lm(niche.pos.water~ fpos2, data = Met.dat)
summary(D2)
PD2 <- ggplotRegression(lm(niche.pos.water ~ fpos2, data = Met.dat))+
  xlab(c('ISTV(Stem mass fraction)'))

D3 <- lm(niche.pos.water ~ fpos3, data = Met.dat)
summary(D3)
PD3 <- ggplotRegression(lm(niche.pos.water ~ fpos3, data = Met.dat))+
  xlab(c('ISTV(Root mass fraction)'))

D4 <- lm(niche.pos.water ~ fpos4, data = Met.dat)
summary(D4)
PD4 <- ggplotRegression(lm(niche.pos.water ~ fpos4, data = Met.dat))+
  xlab(c('ISTV(Specific stem length)'))

D5 <- lm(niche.pos.water ~ fpos5, data = Met.dat)
summary(D5)
PD5 <- ggplotRegression(lm(niche.pos.water ~ fpos5, data = Met.dat))+
  xlab(c('ISTV(Specific root length)'))

D6 <- lm(niche.pos.water ~ fpos6, data = Met.dat)
summary(D6)
PD6 <- ggplotRegression(lm(niche.pos.water ~ fpos6, data = Met.dat))+
  xlab(c('ISTV(SLA)'))

D7 <- lm(niche.pos.water ~ fpos7, data = Met.dat)
summary(D7)
PD7 <- ggplotRegression(lm(niche.pos.water ~ fpos7, data = Met.dat))+
  xlab(c('ISTV(Leaf area index)'))

multiplot(PD1, PD2, PD3, PD4, PD5, PD6, PD7, cols = 2)

#--------------------------------------------------------------------
# Niche Position X niche width
#------------------------------------------------

D1 <- lm(niche.pos.water ~ niche.width.water, data = Met.dat)
summary(D1)
PD1 <- ggplotRegression(lm(niche.pos.water1 ~ niche.width.water, data = Met.dat))

D2 <- lm(niche.pos.light ~ niche.width.light, data = Met.dat)
summary(D2)
PD2 <- ggplotRegression(lm(niche.pos.light1 ~ niche.width.light, data = Met.dat))

D3 <- lm(niche.pos.water ~ niche.width.light, data = Met.dat)
summary(D3)
PD3 <- ggplotRegression(lm(niche.pos.water1 ~ niche.width.light, data = Met.dat))

D4 <- lm(niche.pos.light ~ niche.width.water, data = Met.dat)
summary(D4)
PD4 <- ggplotRegression(lm(niche.pos.light1 ~ niche.width.water, data = Met.dat))


multiplot(PD1, PD2, PD3, PD4, cols = 2)

#----------------------------------
#   Fucntional position x Abundance
#------------------------------------

C1 <- lm(fpos1 ~ abundance, data = Met.dat)
summary(C1)
PC1 <- ggplotRegression(lm(fpos1 ~ abundance, data = Met.dat))+
  xlab(c('ISTV(Leaf mass fraction)'))

C2 <- lm(fpos2 ~ abundance, data = Met.dat)
summary(C2)
PC2 <- ggplotRegression(lm(fpos2 ~ abundance, data = Met.dat))+
  xlab(c('ISTV(Stem mass fraction)'))

C3 <- lm(fpos3 ~ abundance, data = Met.dat)
summary(C3)
PC3 <- ggplotRegression(lm(fpos3 ~ abundance, data = Met.dat))+
  xlab(c('ISTV(Root mass fraction)'))


C4 <- lm(fpos4 ~ abundance, data = Met.dat)
summary(C4)
PC4 <- ggplotRegression(lm(fpos4 ~ abundance, data = Met.dat))+
  xlab(c('ISTV(Specific stem length)'))


C5 <- lm(fpos5 ~ abundance, data = Met.dat)
summary(C5)
PC5 <- ggplotRegression(lm(fpos5 ~ abundance, data = Met.dat))+
  xlab(c('ISTV(Specific root length)'))


C6 <- lm(fpos6 ~ abundance, data = Met.dat)
summary(C6)
PC6 <- ggplotRegression(lm(fpos6 ~ abundance, data = Met.dat))+
  xlab(c('ISTV(SLA)'))


C7 <- lm(fpos7 ~ abundance, data = Met.dat)
summary(C7)
PC7 <- ggplotRegression(lm(fpos7 ~ abundance, data = Met.dat))+
  xlab(c('ISTV(Leaf area index)'))

C8 <- lm(fpos8 ~ abundance, data = Met.dat)
summary(C8)
PC8 <- ggplotRegression(lm(fpos8 ~ abundance, data = Met.dat))+
  xlab(c('ISTV(LDMC)'))


multiplot(PC1, PC2, PC3, PC4,PC5, PC6, PC7, PC8, cols = 2)


#---------

#work in progress
#-----------------------------------------------------------
#non parametric - pearson's coefficient

library("ggpubr")
cor.test(Met.dat$ist1, Met.dat$abundance)

ggscatter(Met.dat, x = "abundance", y = "ist1", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "ISTV", ylab = "ABUNDANCE")
ggscatter(Met.dat, x = "ist2", y = "abundance", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "ISTV", ylab = "abundance")
ggscatter(Met.dat, x = "ist3", y = "abundance", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "istv3", ylab = "ABUNDANCE")
ggscatter(Met.dat, x = "ist4", y = "abundance", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "istv4", ylab = "ABUNDANCE")
ggscatter(Met.dat, x = "ist5", y = "abundance", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "istv5", ylab = "ABUNDANCE")
ggscatter(Met.dat, x = "ist6", y = "abundance", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "istv6", ylab = "ABUNDANCE")
ggscatter(Met.dat, x = "ist7", y = "abundance", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "istv7", ylab = "ABUNDANCE")

#-------------------------------------------------------------------------
# PCA
#-----------------------

library(factoextra)
library("FactoMineR")


#-----------------------------------


