setwd("/Users/pengfoen/OneDrive - University of Connecticut/poolseq_LoF")

library(data.table)


lof_filelist = list.files(path = "./Data/Variants", pattern=".unfiltered_lof.bcf", full.names = T)
files = sapply(lof_filelist, function(x) strsplit(strsplit(x,split = "/")[[1]][4],split="_")[[1]][1])
lof_all <- do.call(rbind,Map(cbind, lapply(lof_filelist, data.table::fread, sep=" "), lakes = files))
colnames(lof_all) <- c("LG", "Pos", "REF", "ALT", "DP","DP4","QUAL","MQ","IDV","gene","conseq","lake")

# change indel read depth to numeric
lof_all[IDV==".",IDV:=NA]
lof_all[,IDV:=as.numeric(IDV)]
lof_all[, LGn := as.numeric(as.roman(substr(LG, 4, nchar(LG))))]

# calculate the number of alterative allele read depth, DP4: ref forwad reads, ref reverse reads, alt forward, alt reverse
alf_calculation<-function(DP4){
  temp<-lapply(DP4, function(x) as.numeric(unlist(strsplit(x,split=","))))
  dp_alt<-unlist(lapply(temp, function(x) x[4]+x[3]))
  return(dp_alt)
}
lof_all[, DP_ALT := alf_calculation(DP4) ]
lof_all[,ALF:=DP_ALT/DP]

# filter low quality records, cannot use IDV==1,as it drops NA as well
#lof_all[QUAL<10|DP<10|MQ<30|DP<10|IDV<3 ,.N]
lof_all_filtered<-lof_all[!(QUAL<10|DP<11| IDV %in% 1) ,]

# for the ones that are within a couple of bps, only keep one
setorderv(lof_all_filtered,c("lake","LGn","Pos"))
lof_all_filtered[,Grouping:=Pos]
lof_all_filtered[,Pos_lead := shift(Pos,1L,type="lead")]
lof_all_filtered[,Pos_lag := shift(Pos,1L,type="lag")]
lof_all_filtered[ Pos-Pos_lag<5 | Pos_lead-Pos<5, Grouping:= 0L]
lof_all_filtered<-lof_all_filtered[lof_all_filtered[,.I[DP_ALT==max(DP_ALT)],by=rleid(Grouping)]$V1]
lof_all_filtered[,c("Pos_lead","Pos_lag","LG"):=NULL]

# create a by locus data.table with allele frequency of each lake as columns
lof_all_filtered[,paste0("ALF_", files):=lapply(files, function(i) ALF*(as.integer(lake==i)))]
lof_all_filtered[,(16:23):=lapply(.SD, max),.SDcols=16:23, by=list(LGn, Pos)]
lof_all_filtered_bylocus<-unique(lof_all_filtered, by=c("LGn","Pos"))
lof_all_filtered_bylocus<-lof_all_filtered_bylocus[,c("LGn","Pos","gene","REF","ALT","conseq","DP","IDV","DP_ALT",
                                     "ALF_boot","ALF_echo","ALF_fred","ALF_gos","ALF_law","ALF_pach","ALF_rob","ALF_say")]


# create a by gene data.table, the allele frequency is filled with max value from any LOF alleles in that gene
lof_all_filtered_bylocus[,(10:17):=lapply(.SD, max),.SDcols=10:17, by=list(LGn, gene)]
lof_all_filtered_bygene<-unique(lof_all_filtered_bylocus, by=c("LGn","gene"))
