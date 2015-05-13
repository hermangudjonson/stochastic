
# Parsing and plots for stochastic gillespie runs

setwd('~/Documents/Dinner/fly eye/stochastic/gillespie')

load_gillespie = function(filename,ntrials) {
  con = file(filename,'r') # open connection 
  # using scan to retrieve header info line by line
  times = head(scan(con,nlines=1,sep='\t'),-1)
  species = head(scan(con,nlines=1,sep='\t',what=character()),-1)
  reactions = head(scan(con,nlines=1,sep='\t',what=character()),-1)
  parameters = head(scan(con,nlines=1,sep='\t',what=character()),-1) # could be parsed out

  concentrations.mat = matrix(nrow=ntrials*length(times),ncol=length(species)+2)
  colnames(concentrations.mat) = c(species,'times','trial')
  rxn_counts.mat = matrix(nrow=ntrials,ncol=length(reactions))
  colnames(rxn_counts.mat) = reactions
  for (i in 1:ntrials) {
    # concentration data in times x species format
    rdata = matrix(scan(con,nlines=length(times),sep='\t'),
                   nrow=length(times),ncol=length(species)+1,byrow=T)
    concentrations.mat[(1+(i-1)*length(times)):(i*length(times)),] = cbind(rdata[,1:length(species)],times,i)
    
    # reaction counts in one line 
    rxndata = head(scan(con,nlines=1,sep='\t'),-1)
    rxn_counts.mat[i,] = rxndata
  }
  
  close(con)
  
  ## return extracted data ##
  results = list()
  results$concentrations.mat = concentrations.mat
  results$rxn_counts.mat = rxn_counts.mat
  results$parameters = parameters
  return(results)
}

## load relevant data ##
tst = load_gillespie('results/yan_network.txt',3)
pre_full = load_gillespie('results/yan_egfr_pre_full.txt',100)
diff_full = load_gillespie('results/yan_egfr_diff_full.txt',100)
pre_half = load_gillespie('results/yan_egfr_pre_half.txt',100)
diff_half = load_gillespie('results/yan_egfr_diff_half.txt',100)
pre_increase = load_gillespie('results/yan_egfr_pre_increase.txt',100)
diff_increase = load_gillespie('results/yan_egfr_diff_increase.txt',100)

## plotting library things ##
require(ggplot2)
require(reshape2)
require(plyr)
require(ggthemes)

## basic plot ## 
bonus_concs.df = data.frame(concentrations.mat)
bonus_concs.df$YTotal = with(bonus_concs.df,Y + YP + M_Y + M_YP + 2*Y_Y)
conc.mlt = melt(bonus_concs.df,id.vars=c('times','trial'))
sub_cond = conc.mlt$variable %in% c('YTotal','P1','P2P','M','mR7','M_Y') & conc.mlt$trial %in% 1:10
# sub_cond = conc.mlt$variable %in% c('P1','P2P','mR7') & conc.mlt$trial %in% 1:10
# sub_cond = conc.mlt$variable %in% unique(conc.mlt$variable)
g = ggplot(conc.mlt[sub_cond,],aes(times/3600,value,group=interaction(trial,variable),col=variable))
g = g + geom_line()
g

# ggsave('plots/yan_egfr_diff_unstable.png')

## coefficient of variation plot ##
YT.df = bonus_concs.df[,c('times','YTotal')]
YT.summary = ddply(YT.df,c('times'),summarise,mean=mean(YTotal),sd=sd(YTotal),
                   cv=sd(YTotal)/mean(YTotal))

g = ggplot(YT.summary,aes(times,cv))
g = g + geom_line(aes(times/3600,cv))
# g = g + scale_y_continuous(limits=c(0,0.5))
g

# ggsave('plots/cv_yan_egfr_pre_half.png')

## convert concentrations to data frame with Ytotal ##
convert_concs = function(concentrations.mat) {
  # get total Yan column
  bonus_concs.df = data.frame(concentrations.mat)
  bonus_concs.df$YTotal = with(bonus_concs.df,Y + YP + M_Y + M_YP + 2*Y_Y)
  return(bonus_concs.df)
}

### nicer plot functions ###
plot_mean_course = function(concentrations.mat,outfile,plot_legend=F,trial=1) {
  bonus_concs.df = convert_concs(concentrations.mat)
  species_toplot = c('YTotal','M_Y','M','mR7','P2P','P1')
  species_names = c('Total Yan','Mae-Yan','Mae','miR-7','pPntP2','PntP1')
  sub_conc.df = bonus_concs.df[,c(species_toplot,'times','trial')]
  sub_conc.mlt = melt(sub_conc.df,id.vars=c('times','trial'))
  sub_conc.summary = ddply(sub_conc.mlt[,c('times','variable','value')],c('times','variable'),
                           summarise,mean=mean(value))
  
  # overlay single stochastic run
  trial_sub_conc.mlt = sub_conc.mlt[sub_conc.mlt$trial %in% trial,]
  
  # plot of these values
  g = ggplot(sub_conc.summary,aes(times/3600,mean,col=variable))
  g = g + geom_line(size=1.5)
  # plot stochastic trial overlay
  g = g + geom_line(data=trial_sub_conc.mlt,aes(times/3600,value,col=variable),size=0.7,alpha=0.5)
  ## axis ##
  g = g + scale_x_continuous(breaks=seq(0,60,10),limits=c(0,60))
  g = g + scale_y_continuous(breaks=seq(0,600,100),limits=c(0,600))
  ## labels ##
  if (plot_legend) {
    g = g + scale_colour_discrete(breaks=species_toplot,labels=species_names)
  }
  else {
    g = g + scale_colour_discrete(guide=F,breaks=species_toplot,labels=species_names)
  }
  ## theme things ##
  g = g + theme_minimal()
  g = g + theme(panel.border=element_rect(fill=NA,size=1))
  g = g + theme(legend.title=element_blank()) # no legend title
  g = g + theme(axis.title=element_blank()) # no value labels
  g = g + theme(axis.text=element_text(size=16,face='bold',family='Helvetica'))
  g = g + theme(legend.text=element_text(size=14,family='Helvetica'))
  g = g + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
  # g = g + theme(panel.grid.major = element_line(size = 0.25, color = "grey"))
  g
  
  ggsave(outfile)

}

# plot mean courses
plot_mean_course(pre_full$concentrations.mat,'plots/mean_yan_egfr_pre_full.pdf',trial=4)
plot_mean_course(diff_full$concentrations.mat,'plots/mean_yan_egfr_diff_full.pdf',trial=4)

plot_mean_course(pre_full$concentrations.mat,'plots/mean_yan_egfr_pre_full_legend.pdf',trial=4,plot_legend=T)

## plot mutant yan courses ##
consolidate_mutants = function(mut_data,mut_names) {
  ytotal.df = data.frame(times=numeric(0),trial=numeric(0),YTotal=numeric(0),type=character(0))
  for (i in 1:length(mut_data)) {
    conc.df = convert_concs(mut_data[[i]]$concentrations.mat)
    conc.df$type = mut_names[i]
    ytotal.df = rbind(ytotal.df,conc.df[,c('times','trial','YTotal','type')])
  }
  return(ytotal.df)
}

pre_ytotal.df = consolidate_mutants(list(pre_full,pre_half,pre_increase),c('WT','DECREASE','INCREASE'))
diff_ytotal.df = consolidate_mutants(list(diff_full,diff_half,diff_increase),c('WT','DECREASE','INCREASE'))

plot_yquantiles = function(ytotal.df,outfile,plot_legend=F) {
  ytotal.summary = ddply(ytotal.df,c('times','type'),summarise,
                         median=median(YTotal),low=quantile(YTotal,0.25),high=quantile(YTotal,0.75))
  ytotal.summary$type = factor(ytotal.summary$type,c('WT','DECREASE','INCREASE'))
  types = c('WT','DECREASE','INCREASE')
  type_labels = c('1.0*dpERK','0.3*dpERK','2.0*dpERK')
  type_colors = c('blue','orange','red1')
  
  g = ggplot(ytotal.summary,aes(times/3600,median,col=type))
  g = g + geom_ribbon(aes(ymin=low,ymax=high,fill=type),alpha=0.4,size=0)
  g = g + geom_line(size=1.5)
  ## labels ##
  if(plot_legend) {
    g = g + scale_colour_manual(breaks=types,labels=type_labels,values=type_colors)
    g = g + scale_fill_manual(breaks=types,labels=type_labels,values=type_colors)
  }
  else {
    g = g + scale_colour_manual(guide=F,breaks=types,labels=type_labels,values=type_colors)
    g = g + scale_fill_manual(guide=F,breaks=types,labels=type_labels,values=type_colors)
  }
  ## axis ##
  g = g + scale_x_continuous(breaks=seq(0,60,10),limits=c(0,60))
  g = g + scale_y_continuous(breaks=seq(0,600,100),limits=c(0,600))
  ## theme things ##
  g = g + theme_minimal()
  g = g + theme(panel.border=element_rect(fill=NA,size=1))
  g = g + theme(legend.title=element_blank()) # no legend title
  g = g + theme(axis.title=element_blank()) # no value labels
  g = g + theme(axis.text=element_text(size=16,face='bold',family='Helvetica'))
  g = g + theme(legend.text=element_text(size=14,family='Helvetica'))
  g = g + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
  # g = g + theme(panel.grid.major = element_line(size = 0.25, color = "grey"))
  g
  
  ggsave(outfile)
}

plot_yquantiles(pre_ytotal.df,'plots/yquantiles_pre.pdf')
plot_yquantiles(diff_ytotal.df,'plots/yquantiles_diff.pdf')
# plot_yquantiles(pre_ytotal.df,'plots/yquantiles_pre_legend.pdf',plot_legend=T)

plot_ycov = function(ytotal.df,outfile,plot_legend=F) {
  ytotal.summary = ddply(ytotal.df,c('times','type'),summarise,
                         cv=sd(YTotal)/mean(YTotal))
  ytotal.summary$type = factor(ytotal.summary$type,c('WT','DECREASE','INCREASE'))
  types = c('WT','DECREASE','INCREASE')
  type_labels = c('1.0*dpERK','0.3*dpERK','2.0*dpERK')
  type_colors = c('blue','orange','red1')
  
  g = ggplot(ytotal.summary,aes(times/3600,cv,col=type))
  g = g + geom_line(size=1.5)
  g
  ## labels ##
  if(plot_legend) {
    g = g + scale_colour_manual(breaks=types,labels=type_labels,values=type_colors)
    g = g + scale_fill_manual(breaks=types,labels=type_labels,values=type_colors)
  }
  else {
    g = g + scale_colour_manual(guide=F,breaks=types,labels=type_labels,values=type_colors)
    g = g + scale_fill_manual(guide=F,breaks=types,labels=type_labels,values=type_colors)
  }
  ## axis ##
  g = g + scale_x_continuous(breaks=seq(0,60,10),limits=c(0,60))
  g = g + scale_y_continuous(limits=c(0,0.4))
  ## theme things ##
  g = g + theme_minimal()
  g = g + theme(panel.border=element_rect(fill=NA,size=1))
  g = g + theme(legend.title=element_blank()) # no legend title
  g = g + theme(axis.title=element_blank()) # no value labels
  g = g + theme(axis.text=element_text(size=16,face='bold',family='Helvetica'))
  g = g + theme(legend.text=element_text(size=14,family='Helvetica'))
  g = g + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
  # g = g + theme(panel.grid.major = element_line(size = 0.25, color = "grey"))
  g
  
  ggsave(outfile)
}

plot_ycov(pre_ytotal.df,'plots/ycov_pre.pdf')
plot_ycov(diff_ytotal.df,'plots/ycov_diff.pdf')
# plot_ycov(pre_ytotal.df,'plots/ycov_pre_legend.pdf',plot_legend=T)

### nicer plot functions end ###

## binding site proportion calculation ##
to_conc = function(c,nsize=78e-15){
  return(c/(nsize * 6.022e23))
}

YTotal = 500
PMax = 500

NStates = 4
KbT = 0.593
dGE = 10.0 # 10 original
dGNS = 5.8
dGSAM = 7
dGP = 12.4 # 9.7 original -> 12.4 max
Kd = 7e-6
KdCopy = Kd*6.022e23*78e-15 # Kd * NA * V
YCopy = (1/4)*(sqrt(KdCopy^2 + 8 * KdCopy * YTotal) - KdCopy)
PCopies = seq(0,PMax,length.out=1000)
fractions = matrix(NA,nrow=length(PCopies),ncol=NStates+1)
colnames(fractions) = c('unbound','single_yan','dimer_yan','pnt','pnt_copy')
for (i in 1:length(PCopies)) {
  bweights = c(1,
               to_conc(YCopy) * (exp(dGE/KbT) + exp(dGNS/KbT)),
               to_conc(YCopy)^2 * (exp((dGE+dGNS+dGSAM)/KbT) + exp((dGE+dGNS)/KbT)/Kd),
               to_conc(PCopies[i]) * exp(dGP/KbT)
               )
  fractions[i,] = c(bweights/sum(bweights),PCopies[i])
}

fractions.mlt = melt(data.frame(fractions),id.vars='pnt_copy')
g = ggplot(fractions.mlt,aes(x=pnt_copy,y=value,fill=variable))
g = g + geom_bar(stat='identity')
g = g + theme(axis.text.x = element_text(size=16),axis.text.y = element_text(size=16),
              text = element_text(size=14))
g

# ggsave('plots/ypsite_fractions_y500_1nM.png')

### calculations for yan promoter ###
YTotal = 500
PMax = 250

NStates = 4
KbT = 0.593
dGE = 5.8 # 10
dGNS = 5.8
dGSAM = 7
dGP = 14.0 # 9.7 original
Kd = 7e-6
KdCopy = Kd*6.022e23*78e-15 # Kd * NA * V
YCopy = (1/4)*(sqrt(KdCopy^2 + 8 * KdCopy * YTotal) - KdCopy)
PCopies = seq(0,PMax,length.out=1000)
fractions = matrix(NA,nrow=length(PCopies),ncol=NStates+1)
colnames(fractions) = c('unbound','single_yan','dimer_yan','pnt','pnt_copy')
for (i in 1:length(PCopies)) {
  bweights = c(1,
               to_conc(YCopy) * (exp(dGE/KbT) + exp(dGNS/KbT)),
               to_conc(YCopy)^2 * (exp((dGE+dGNS+dGSAM)/KbT) + exp((dGE+dGNS)/KbT)/Kd),
               to_conc(PCopies[i]) * exp(dGP/KbT)
  )
  fractions[i,] = c(bweights/sum(bweights),PCopies[i])
}

fractions.mlt = melt(data.frame(fractions),id.vars='pnt_copy')
g = ggplot(fractions.mlt,aes(x=pnt_copy,y=value,fill=variable))
g = g + geom_bar(stat='identity')
g = g + theme(axis.text.x = element_text(size=16),axis.text.y = element_text(size=16),
              text = element_text(size=14))
g
