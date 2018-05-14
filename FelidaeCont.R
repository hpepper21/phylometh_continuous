getwd()
setwd("/Users/elipepper/Documents/GitHub")
install.packages("yearn")
yearn::yearn(ape)
yearn::yearn(geiger)
yearn::yearn(OUwie)
library(rotl)
cat.id<-tnrs_match_names("Felidae")$ott_id
cat.tree <- tol_subtree(ott_id=cat.id)
library(devtools)
plot.phylo(cat.tree, cex=0.3)
print(paste("The Felidae tree has ", Ntip(cat.tree), " terminals and ", 
            Nnode(cat.tree), " internal nodes out of ",Ntip(cat.tree)-2,
            " possible, which means it is ", 
            round(100*(Nnode(cat.tree)-1)/(Ntip(cat.tree)-3), 2),
            "% resolved", sep=""))
cat.studies <- studies_find_studies(property="ot:focalCladeOTTTaxonName",
                                    value="Felidae")
cat.studies.ids <- unlist(cat.studies$study_ids)
cat.study1.tree1 <- get_study(cat.studies.ids[[1]])
plot.phylo(cat.study1.tree1, cex=0.3)
print(paste("The Felidae tree has ", Ntip(cat.study1.tree1), " terminals and ", 
            Nnode(cat.study1.tree1), " internal nodes out of ",Ntip(cat.study1.tree1)-2,
            " possible, which means it is ", 
            round(100*(Nnode(cat.study1.tree1)-1)/(Ntip(cat.study1.tree1)-3), 2),
            "% resolved", sep=""))
cat.study2.tree2 <- get_study(cat.studies.ids[[2]])
plot.phylo(cat.study2.tree2, cex = 0.3)
print(paste("The Felidae tree has ", Ntip(cat.study2.tree2), " terminals and ", 
            Nnode(cat.study2.tree2), " internal nodes out of ",Ntip(cat.study2.tree2)-2,
            " possible, which means it is ", 
            round(100*(Nnode(cat.study2.tree2)-1)/(Ntip(cat.study2.tree2)-3), 2),
            "% resolved", sep=""))
cat.study3.tree3 <- get_study(cat.studies.ids[[3]])
plot.phylo(cat.study3.tree3, cex = 0.3)
print(paste("The Felidae tree has ", Ntip(cat.study3.tree3), " terminals and ", 
            Nnode(cat.study3.tree3), " internal nodes out of ",Ntip(cat.study3.tree3)-2,
            " possible, which means it is ", 
            round(100*(Nnode(cat.study3.tree3)-1)/(Ntip(cat.study3.tree3)-3), 2),
            "% resolved", sep=""))
cat.study4.tree4 <- get_study(cat.studies.ids[[4]])
plot.phylo(cat.study4.tree4, cex = 0.3)
print(paste("The Felidae tree has ", Ntip(cat.study4.tree4), " terminals and ", 
            Nnode(cat.study4.tree4), " internal nodes out of ",Ntip(cat.study4.tree4)-2,
            " possible, which means it is ", 
            round(100*(Nnode(cat.study4.tree4)-1)/(Ntip(cat.study4.tree4)-3), 2),
            "% resolved", sep=""))
cat.study5.tree5 <- get_study(cat.studies.ids[[5]])
plot.phylo(cat.study5.tree5, cex = 0.3)
print(paste("The Felidae tree has ", Ntip(cat.study5.tree5), " terminals and ", 
            Nnode(cat.study5.tree5), " internal nodes out of ",Ntip(cat.study5.tree5)-2,
            " possible, which means it is ", 
            round(100*(Nnode(cat.study5.tree5)-1)/(Ntip(cat.study5.tree5)-3), 2),
            "% resolved", sep=""))
cat.study6.tree6 <- get_study(cat.studies.ids[[6]])
plot.phylo(cat.study6.tree6, cex = 0.3)
print(paste("The Felidae tree has ", Ntip(cat.study6.tree6), " terminals and ", 
            Nnode(cat.study6.tree6), " internal nodes out of ",Ntip(cat.study6.tree6)-2,
            " possible, which means it is ", 
            round(100*(Nnode(cat.study6.tree6)-1)/(Ntip(cat.study6.tree6)-3), 2),
            "% resolved", sep=""))
cat.study7.tree7 <- get_study(cat.studies.ids[[7]])
plot.phylo(cat.study7.tree7, cex = 0.3)
print(paste("The Felidae tree has ", Ntip(cat.study7.tree7), " terminals and ", 
            Nnode(cat.study7.tree7), " internal nodes out of ",Ntip(cat.study7.tree7)-2,
            " possible, which means it is ", 
            round(100*(Nnode(cat.study7.tree7)-1)/(Ntip(cat.study7.tree7)-3), 2),
            "% resolved", sep=""))

mamm.dat<-read.delim("Mammal_lifehistories_v2.txt")
View(mamm.dat)
cat.sub<-subset(mamm.dat,mamm.dat$family=="Felidae")
cat.cont.sub<-cat.sub[,c("family","Genus","species","mass.g.")]
View(cat.cont.sub)
library(devtools)
install_github("phylotastic/datelife")
cat.cont.sub$name=paste(cat.cont.sub$Genus,cat.sub$species)
library(datelife)
chrono<-datelife_search(cat.cont.sub$name,summary_format = "phylo_sdm")
plot(chrono)
axisPhylo()

not.matching.positions <- c()
for (i in sequence(nrow(cat.cont.sub))) {
  if(cat.cont.sub$name[i] %in% chrono$tip.label) {
    
  } else {
    not.matching.positions <- c(not.matching.positions, i)
  }
}
cat.sub.clean <- cat.cont.sub[-not.matching.positions,]

rownames(cat.cont.sub)=cat.cont.sub$name
chrono.dat<-treedata(chrono,cat.cont.sub)

VizualizeData<-function(phy,data){
  plot(phy)
  print(data)
}
VizualizeData(phy = chrono.dat$phy,data = chrono.dat$data)

cat.mass<-chrono.dat$data[,"mass.g."]
cat.mass<-as.numeric(cat.mass)
log.cat.mass<-log(cat.mass)
names(log.cat.mass) <-rownames(chrono.dat$data)
View(log.cat.mass)

BM1 <- geiger::fitContinuous(chrono.dat$phy, log.cat.mass, model="BM")
BM1
print(paste("The rate of evolution is", BM1$opt$sigsq, "in units of", "log.cat.mass^2/my"))

OU1 <- fitContinuous(chrono.dat$phy, log.cat.mass, model="OU")
head(OU1)
print(paste("The rate of evolution is", OU1$opt$sigsq, "in units of", "log.cat.mass^2/my"))
par(mfcol=(c(1,1)))
plot(chrono.dat$phy, show.tip.label=FALSE)
ou.tree <- rescale(chrono.dat$phy, model="OU", 0.05408678)
plot(ou.tree)

#compare the trees
library(MuMIn)
AIC.BM1 <- BM1$opt$aicc
AIC.OU1 <- OU1$opt$aicc
delta.AIC.BM1 <- AIC.BM1-AIC.BM1
delta.AIC.OU1 <- AIC.OU1-AIC.BM1
#the comparison seems to indicate that the brownian model is the best fit

one.discrete.char <- c(0,0,1,0,1,1,0,0,1,1,1,0,1,0,0,1,0,1,0,0,1,0,1,0,0,0,0,0,0)
reconstruction.info <- ace(one.discrete.char, chrono.dat$phy, type = "discrete", method = "ML", CI = T)
best.states <- colnames(reconstruction.info$lik.anc)[apply(reconstruction.info$lik.anc, 1, which.max)]
#adding these labels to the tree

nodes <- nodelabels(text = best.states, node = chrono.dat$phy$Nnode, interactive = F, shape = "ellipse", cex = 0.1)

nodeBased.OUMV <- OUwie(chrono.dat$phy, log.cat.mass, model = "OUMV", simmap.tree = FALSE, diagn = FALSE)

#running the OUwie mods:
models <- c("BM1", "BMS", "OU1", "OUM", "OUMV", "OUMA", "OUMVA")
results <- lapply(models, OU1, phy = tree, data = trait1)

AICc.values <- sapply(results, "[[", "AICc")
names(AICc.values)<-models
AICc.values<-AICc.values-min(AICc.values)

print(AICc.values) #The best model is the one with smallest AICc score

best<-results[[which.min(AICc.values)]] #store for later

print(best) #prints info on best model



