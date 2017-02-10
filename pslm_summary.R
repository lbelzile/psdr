# File:     pslm_summary.R
# Author:   Olli Saarela (olli.saarela@utoronto.ca)
# Date:     2016-02-06
# Summary:  Produce Table 1 based on the simulation study output
#  Saarela, O., Belzile, L. R. and D. A. Stephens. A Bayesian view of doubly robust causal inference,
#  Biometrika (2016), 103 (3): 667-681, doi:10.1093/biomet/asw025


rm(list=ls())
outpath <- ""
results <- matrix(unlist(t(sapply(list.files(outpath, pattern = "results"), function(file){read.table(file)}))),ncol=90)[,-1]

# results<-read.table('~/Dropbox/work/ps_olli/data/results_pslm_new.txt')

getpointest <- function(estimator, values, truevalue) {
    return(data.frame(Estimator=estimator,
                      Mean=mean(values),
                      RB=100.0 * (mean(values) - truevalue)/abs(truevalue),
                      SD=sd(values),
                      RMSE=sqrt(var(values) + (mean(values) - truevalue)^2)))
}

getvarest <- function(estimator, pevalues, vevalues, truevalue) {
    return(data.frame(Estimator=estimator,
                      Mean=mean(sqrt(vevalues)),
                      RB=100.0 * (mean(vevalues) - var(pevalues))/var(pevalues),
                      SDx100=100*sqrt(var(vevalues)),
                      RMSEx100=100*sqrt(var(vevalues) + (mean(vevalues) - var(pevalues))^2),
                      CProb=100.0 * mean((truevalue > (pevalues + qnorm(0.025) * sqrt(vevalues))) &
                                         (truevalue < (pevalues + qnorm(0.975) * sqrt(vevalues))))))
}

getciprob <- function(estimator, cilvalues, ciuvalues, truevalue) {
    return(data.frame(Estimator=estimator,
                      CProb=100.0 * mean((truevalue > cilvalues) &
                                         (truevalue < ciuvalues))))
}

#Output table

tab <- NULL
tab <- rbind(tab, cbind(getpointest('Oracle', results[,65], mean(results[,1])),
    getvarest('Oracle', results[,65], results[,66], mean(results[,1]))[c(2,3,6)]))
tab <- rbind(tab, cbind(getpointest('IPTW', results[,2], mean(results[,1])),
    getvarest('IPTW', results[,2], results[,68], mean(results[,1]))[c(2,3,6)]))
tab <- rbind(tab, cbind(getpointest('Doubly robust', results[,67], mean(results[,1])),
	    getvarest('IPTW-DR', results[,67], results[,69], mean(results[,1]))[c(2,3,6)]))
tab <- rbind(tab, cbind(getpointest('PS regression/unadjusted', results[,9], mean(results[,1])),
    getvarest('DRPS unadjusted', results[,17], results[,18], mean(results[,1]))[c(2,3,6)]))
tab <- rbind(tab, cbind(getpointest('PS regression/robust', results[,9], mean(results[,1])),
    getvarest('DRPS robust', results[,17], results[,19], mean(results[,1]))[c(2,3,6)]))
tab <- rbind(tab, cbind(getpointest('PS regression/adjusted', results[,9], mean(results[,1])),
    getvarest('DRPS adjusted', results[,17], results[,20], mean(results[,1]))[c(2,3,6)]))
interm<-cbind(getpointest('Clever covariate', results[,11], mean(results[,1])),t(c(NA,NA,NA)))
colnames(interm)<-colnames(tab)
tab<-rbind(tab,interm)

tab <- rbind(tab, cbind(getpointest('Joint', results[,21], mean(results[,1])),
    getvarest('Joint', results[,21], results[,22], mean(results[,1]))[c(2,3,6)]))
tab <- rbind(tab, cbind(getpointest('Forward sampling', results[,23], mean(results[,1])),
    getvarest('MI-PS', results[,23], results[,25], mean(results[,1]))[c(2,3,6)]))
tab <- rbind(tab, cbind(getpointest('Forward sampling (DR)', results[,24], mean(results[,1])),
    getvarest('MI-DRPS', results[,24], results[,26], mean(results[,1]))[c(2,3,6)]))
tab <- rbind(tab, cbind(getpointest('Forward sampling-PS0', results[,35], mean(results[,1])),
    getvarest('MI-PS', results[,35], results[,37], mean(results[,1]))[c(2,3,6)]))
tab <- rbind(tab, cbind(getpointest('Forward sampling (DR)-PS0', results[,36], mean(results[,1])),
    getvarest('MI-DRPS', results[,36], results[,38], mean(results[,1]))[c(2,3,6)]))
interm<-cbind(getpointest('Forward sampling CC-PS0', results[,39], mean(results[,1])),t(c(NA,NA,NA)))
colnames(interm)<-colnames(tab)
tab<-rbind(tab,interm)

interm<-cbind(getpointest('Forward sampling CC (DR)-PS0', results[,40], mean(results[,1])),t(c(NA,NA,NA)))
colnames(interm)<-colnames(tab)
tab<-rbind(tab,interm)
tab <- rbind(tab, cbind(getpointest('Pseudo-posterior', results[,27], mean(results[,1])),
    getvarest('PP-PS', results[,27], results[,29], mean(results[,1]))[c(2,3,6)]))
tab <- rbind(tab, cbind(getpointest('Pseudo-posterior (DR)', results[,28], mean(results[,1])),
    getvarest('PP-DRPS', results[,28], results[,30], mean(results[,1]))[c(2,3,6)]))
tab <- rbind(tab, cbind(getpointest('Bootstrap', results[,35], mean(results[,1])),
    getvarest('BS-PS', results[,35], results[,37], mean(results[,1]))[c(2,3,6)]))
tab <- rbind(tab, cbind(getpointest('Bootstrap (DR)', results[,36], mean(results[,1])),
    getvarest('BS-DRPS', results[,36], results[,38], mean(results[,1]))[c(2,3,6)]))
tab <- rbind(tab, cbind(getpointest('IS-Dirichlet', results[,49], mean(results[,1])),
    getvarest('IS-Dirichlet', results[,49], results[,51], mean(results[,1]))[c(2,3,6)]))
tab <- rbind(tab, cbind(getpointest('IS-Dirichlet-adj', results[,50], mean(results[,1])),
    getvarest('IS-Dirichlet-adj', results[,50], results[,52], mean(results[,1]))[c(2,3,6)]))
tab <- rbind(tab, cbind(getpointest('IS-Multinomial', results[,57], mean(results[,1])),
    getvarest('IS-Multinomial', results[,57], results[,59], mean(results[,1]))[c(2,3,6)]))
tab <- rbind(tab, cbind(getpointest('IS-Multinomial-adj', results[,58], mean(results[,1])),
    getvarest('IS-Multinomial-adj', results[,58], results[,60], mean(results[,1]))[c(2,3,6)]))
library(xtable)
rownames(tab)=sapply(1:nrow(tab),function(x){paste0("[",tolower(as.roman(x)),"] ",tab[x,1])})
tab<-tab[,-1]
colnames(tab)=c("Mean", "RB (\\%)", colnames(tab[3:4]),"SE", "RB$_{\\mathrm{SE}}$ (\\%)", "CP")

mat<-xtable(tab,digits=c(2,2,1,rep(2,3),1,1),align=c("l","c","r",rep("c",3),"r","c"))
            print(mat, sanitize.text.function = function(x){x},
                  floating=FALSE,
                  hline.after=NULL,
                  add.to.row=list(pos=list(-1,0, nrow(mat)),
                  command=c('\\toprule ',
                            '\\midrule ',
                            '\\bottomrule ')))
#Empirical coverage
cover <- NULL
cover <- rbind(cover,getciprob('Pseudo-posterior', results[,31], results[,32], mean(results[,1])))
cover <- rbind(cover,getciprob('Pseudo-posterior (DR)', results[,33], results[,34], mean(results[,1])))
cover <- rbind(cover,getciprob('Bootstrap', results[,39], results[,40], mean(results[,1])))
cover <- rbind(cover,getciprob('Bootstrap (DR)', results[,41], results[,42], mean(results[,1])))
cover <- rbind(cover,getciprob('IS-Dirichlet', results[,53], results[,54], mean(results[,1])))
cover <- rbind(cover,getciprob('IS-Dirichlet-adj', results[,55], results[,56], mean(results[,1])))
cover <- rbind(cover,getciprob('IS-Multinomial', results[,61], results[,62], mean(results[,1])))
cover <- rbind(cover,getciprob('IS-Multinomial-adj', results[,63], results[,64], mean(results[,1])))

mat2<-xtable(cover,digits=rep(3,3),align=c("l",rep("r",2)))
            print(mat2, sanitize.text.function = function(x){x},
                  floating=FALSE,
                  hline.after=NULL,
                  add.to.row=list(pos=list(-1,0, nrow(mat2)),
                  command=c('\\toprule ',
                            '\\midrule ',
                            '\\bottomrule ')))
