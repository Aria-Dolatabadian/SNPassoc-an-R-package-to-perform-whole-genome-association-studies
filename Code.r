#https://cran.r-project.org/web/packages/SNPassoc/vignettes/SNPassoc.html

library(SNPassoc)

#Load data from the package

data(asthma, package = "SNPassoc")
str(asthma, list.len=9)

#Or read from WD
#asthma <- read.csv("asthma_data.csv")

asthma[1:5, 1:8]

asthma.s <- setupSNP(data=asthma, colSNPs=7:ncol(asthma), sep="")

idx <- grep("^rs", colnames(asthma))
asthma.s <- setupSNP(data=asthma, colSNPs=idx, sep="")

head(asthma.s$rs1422993)

class(asthma.s$rs1422993)

summary(asthma.s$rs1422993)

plot(asthma.s$rs1422993)

plot(asthma.s$rs1422993, type=pie)

summary(asthma.s, print=FALSE)

plotMissing(asthma.s, print.labels.SNPs = FALSE)

hwe <- tableHWE(asthma.s)
head(hwe)

hwe2 <- tableHWE(asthma.s, casecontrol)

#SNPs is HWE in the whole sample but not controls
snpNHWE <- hwe2[,1]>0.05 & hwe2[,2]<0.05
rownames(hwe2)[snpNHWE]


hwe2[snpNHWE,]



snps.ok <- rownames(hwe2)[hwe2[,2]>=0.001]
pos <- which(colnames(asthma)%in%snps.ok, useNames = FALSE)
asthma.s <- setupSNP(asthma, pos, sep="")



association(casecontrol ~ rs1422993, data = asthma.s)



maxstat(asthma.s$casecontrol, asthma.s$rs1422993)

association(casecontrol ~ rs1422993, asthma.s, model="dominant")


association(casecontrol ~ rs1422993 + country + smoke, asthma.s)


association(casecontrol ~ rs1422993 + survival::strata(gender), asthma.s)


association(casecontrol ~ rs1422993, asthma.s, 
                subset=country=="Spain")



association(bmi ~ rs1422993, asthma.s) 


ans <- WGassociation(casecontrol, data=asthma.s)
head(ans)



ans.adj <- WGassociation(casecontrol ~ country + smoke, asthma.s)
head(ans.adj)

ans.fast <- scanWGassociation(casecontrol, asthma.s)


devtools::install_github("isglobal-brge/SNPassoc")

plot(ans)

ans.max <- maxstat(asthma.s, casecontrol)
ans.max

#minimum P-value across SNPs
min(ans.max["Pr(>z)",])

infoTable <- WGstats(ans)

infoTable$rs1422993

library(xtable)
out <- getNiceTable(ans[c("rs1422993", "rs184448")])

nlines <- attr(out, "nlines")
hlines <- c(-1, -1, 0, cumsum(nlines+1), nrow(out), nrow(out))

print(xtable(out, caption='Genetic association using
                different genetic models from asthma 
                data example of rs1422993 and rs184448 
                SNPs obtained with SNPassoc.',
             label = 'tab-2SNPs'),
      tabular.enviroment="longtable", file="tableSNPs",
      floating=FALSE,  include.rownames = FALSE, 
      hline.after= hlines, sanitize.text.function=identity)


association(casecontrol ~ dominant(rs1422993)*factor(smoke), 
            data=asthma.s)

association(casecontrol ~ rs1422993*factor(rs184448), 
            data=asthma.s, model.interaction = "dominant" )



ans <- WGassociation(casecontrol, data=asthma.s)
mask <- apply(ans, 1, function(x) min(x, na.rm=TRUE)<0.1)
sig.snps <- names(mask[mask])
sig.snps

idx <- which(colnames(asthma)%in%sig.snps)
asthma.s2 <- setupSNP(asthma, colSNPs = idx, sep="")
ans.int <- interactionPval(casecontrol ~ 1, data=asthma.s2)
ans.int

plot(ans.int)

library(haplo.stats)
snpsH <- c("rs714588", "rs1023555",  "rs898070")
genoH <- make.geno(asthma.s, snpsH)
em <- haplo.em(genoH, locus.label = snpsH, miss.val = c(0, NA))
em


trait <- asthma.s$casecontrol
mod <- haplo.glm(trait ~ genoH,           
                 family="binomial", 
                 locus.label=snpsH,
                 allele.lev=attributes(genoH)$unique.alleles,
                 control = haplo.glm.control(haplo.freq.min=0.05))   
intervals(mod)



snpsH2 <- labels(asthma.s)[6:15]
genoH2 <- make.geno(asthma.s, snpsH2)
haplo.score <- list()
for (i in 4:7) {
 trait <- asthma.s$casecontrol
 haplo.score[[i-3]] <- haplo.score.slide(trait, genoH2, 
                          trait.type="binomial",
                          n.slide=i,
                          simulate=TRUE,
                          sim.control=score.sim.control(min.sim=100,
                                       max.sim=200)) 
 }


par(mfrow=c(2,2))
for (i in 4:7) {
    plot(haplo.score[[i-3]])
    title(paste("Sliding Window=", i, sep=""))
 }


snpsH3 <- snpsH2[4:7]
genoH3 <- make.geno(asthma.s, snpsH3)
mod <- haplo.glm(trait~genoH3,           
                 family="binomial", 
                 locus.label=snpsH3,
                 allele.lev=attributes(genoH3)$unique.alleles,
                 control = haplo.glm.control(haplo.freq.min=0.05))      
intervals(mod)

lrt <- mod$lrt
pchisq(lrt$lrt, lrt$df, lower=FALSE)


smoke <- asthma.s$smoke
mod.adj.ref <- glm(trait ~ smoke, family="binomial")
mod.adj <- haplo.glm(trait ~ genoH3 + smoke ,           
                 family="binomial", 
                 locus.label=snpsH3,
                 allele.lev=attributes(genoH3)$unique.alleles,
                 control = haplo.glm.control(haplo.freq.min=0.05))

lrt.adj <- mod.adj.ref$deviance - mod.adj$deviance
pchisq(lrt.adj, mod.adj$lrt$df, lower=FALSE)

sessionInfo()






