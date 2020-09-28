setwd("/Users/pengfoen/OneDrive - University of Connecticut/poolseq_LoF")

library(data.table)
{
gene_info<-fread("protein_coding_gene.csv")
# gene_info = fread("../Analysis_Expression/Data files RAW/protein_coding_gene.csv")
# gene_info[,V1 := as.numeric(as.roman(substr(V1, 6, nchar(V1))))]
# gene_info[,c("V2","V7"):=NULL]
# gene_info[V8=="ensembl",V8:=NA]
# setnames(gene_info,c("V1","V3","V4","V5","V6","V8"),c("LGn","start","stop","strand","GeneID","gene.name"))
# fwrite(gene_info,"protein_coding_gene.csv")

CDS_info<-fread("protein_coding_CDS.csv")
# CDS_info = fread("../Analysis_Expression/Data files RAW/protein_coding_CDS.csv")
# CDS_info[,V1 := as.numeric(as.roman(substr(V1, 6, nchar(V1))))]
# CDS_info[,V2:=NULL]
# setnames(CDS_info,c("V1","V3","V4","V5","V6","V7","V8"),c("LGn","start","stop","strand","gene_ID","transcript_ID","exon_number"))
# fwrite(CDS_info,"protein_coding_CDS.csv")

lof_filelist = list.files(path = "./Data/Variants", pattern=".unfiltered_lof.bcf", full.names = T)
files = sapply(lof_filelist, function(x) strsplit(strsplit(x,split = "/")[[1]][4],split="_")[[1]][1])
lof_all <- do.call(rbind,Map(cbind, lapply(lof_filelist, data.table::fread, sep=" "), lakes = files))
colnames(lof_all) <- c("LG", "Pos", "REF", "ALT", "DP","DP4","QUAL","MQ","IDV","gene","conseq","lake")
}

#### remove duplicated record at the same position, keep the row with the largest IDV (indel supporting reads)
lof_all_indel<-lof_all[IDV!=".",]
lof_all_snp<-lof_all[IDV==".",]
lof_all_indel<-lof_all_indel[lof_all_indel[,.I[which.max(IDV)],by=rleid(lake,LG,gene,Pos,conseq)]$V1]
lof_all<-merge(lof_all_indel,lof_all_snp,all=T)
rm(lof_all_indel,lof_all_snp)
#duplicated<-lof_all_noduplication[duplicated(lof_all_noduplication,by=c("lake","LG","gene","Pos","conseq")),]

# parse the consequence
lof_all[,conseq_parsed:=lapply(conseq,strsplit,split="&")]
unique(rapply(lof_all[,conseq_parsed], function(x) unlist(x))) # list all types of consequence

lof_all[,frameshift:=unlist(lapply(conseq_parsed, function(x) "frameshift" %in% x | "frameshift" %in% unlist(x)))]
lof_all[,stop_gained:=unlist(lapply(conseq_parsed, function(x) "stop_gained" %in% x | "stop_gained" %in% unlist(x)))]
lof_all[,stop_lost:=unlist(lapply(conseq_parsed, function(x) "stop_lost" %in% x | "stop_lost" %in% unlist(x)))]
lof_all[,start_lost:=unlist(lapply(conseq_parsed, function(x) "start_lost" %in% x | "start_lost" %in% unlist(x)))]
lof_all[,splice:=unlist(lapply(conseq_parsed, function(x) any(c("splice_acceptor","splice_donor") %in% x) | any(c("splice_acceptor","splice_donor") %in% unlist(x)) ))]

lof_all[,conseq_parsed:=NULL]


# change indel read depth to numeric
lof_all[IDV==".",IDV:=NA]
lof_all[,IDV:=as.numeric(IDV)]
lof_all[, LGn := as.numeric(as.roman(substr(LG, 4, nchar(LG))))]

# calculate the number of alternative allele read depth, DP4: ref forwad reads, ref reverse reads, alt forward, alt reverse
DP4_calculation<-function(DP4){
  temp<-lapply(DP4, function(x) as.numeric(unlist(strsplit(x,split=","))))
  dp_alt_f<-unlist(lapply(temp, function(x) x[3]))
  dp_alt_r<-unlist(lapply(temp, function(x) x[4]))
  dp4_sum<-unlist(lapply(temp, function(x) x[1]+x[2]+x[3]+x[4]))
  return(list(dp_alt_f,dp_alt_r,dp4_sum))
}
lof_all[, c("DP_ALT_f","DP_ALT_r","DP4_sum") := DP4_calculation(DP4) ]
lof_all[,ALF:=(DP_ALT_f+DP_ALT_r)/DP4_sum]

# filter low quality records, cannot use IDV==1,as it drops NA as well
lof_all_filtered<-lof_all[!(QUAL<20 | DP<20 | IDV %in% c(1:3) | MQ<30 | DP_ALT_f<1 | DP_ALT_r<1),]

# find gene ID for all genes, first based on gene name, then based on position
lof_all_filtered[,gene_ID:= ifelse(grepl("^ENS",gene), gene,NA )]
lof_all_filtered[is.na(gene_ID),gene_name:= gene]
lof_all_filtered[is.na(gene_ID),gene_ID:= gene_info[match(gene_name,gene_info[,gene.name]),GeneID]]
temp<-gene_info[lof_all_filtered[is.na(gene_ID)],on=c("LGn", "start<=Pos", "stop>=Pos"),mult="first"]
lof_all_filtered[is.na(gene_ID),gene_ID:= temp[match(gene_name,temp[,gene]),GeneID]]
rm(temp)

### remove lof alleles that map to alternative spliced exons
CDS_info[,exon_length:=stop-start]
CDS_info[,cum_length:=cumsum(exon_length),by=transcript_ID]
CDS_info[,full_length:=sum(exon_length),by=transcript_ID]
CDS_info[,total_exon:=max(exon_number),by=transcript_ID]

mult_transcript<-CDS_info[,unique(transcript_ID),by=gene_ID][,if(.N>1) .SD,,by=gene_ID]
unique_transcript<-CDS_info[!mult_transcript,on=c("transcript_ID"="V1")][,unique(transcript_ID),by=gene_ID]
mult_CDS<-CDS_info[mult_transcript,on=c("transcript_ID"="V1")]
uniq_CDS<-CDS_info[!mult_CDS,on=c("transcript_ID")]
mult_CDS[,core:=.N>1,by=list(gene_ID,start,stop)]
mult_CDS_uniq<-unique(mult_CDS[core==T,],by=c("gene_ID","start","stop"))
core_CDS<-rbind(uniq_CDS,mult_CDS_uniq,fill=T)
core_CDS[,c("start_join","stop_join"):=.(start-2,stop+2)]
lof_all_filtered[,Pos_join:=Pos]
lof_all_filtered[,Pos_end_join:=Pos+nchar(REF)-1]
setkeyv(core_CDS,c("LGn", "gene_ID","start_join","stop_join"))
lof_all_filtered_core_CDS<-foverlaps(lof_all_filtered,core_CDS,by.x=c("LGn", "gene_ID","Pos_join","Pos_end_join"),by.y=c("LGn", "gene_ID","start_join","stop_join"),type="any",mult="all",nomatch=0L)
#temp_duplicated<-temp[duplicated(temp,by=c("gene_ID","start","stop","Pos","lake","LGn")),]
#lof_excluded<-lof_all_filtered_core_CDS[!temp,on=c("lake","LGn","Pos")]
rm(mult_transcript,unique_transcript,mult_CDS,mult_CDS_uniq,uniq_CDS,core_CDS)

# calculate the relative position of the mutation in the gene
#lof_all_filtered_core_CDS<-CDS_info[lof_all_filtered,on=c("LGn", "start_join<=Pos_join", "end_join>=Pos_join")]
lof_all_filtered_core_CDS[Pos >= start & Pos <= stop,Pos_approx:=Pos]
lof_all_filtered_core_CDS[!(Pos >= start & Pos <= stop), Pos_approx:= ifelse(abs(Pos-start)<abs(Pos-stop), start, stop)]
lof_all_filtered_core_CDS[,Pos_gene:=cum_length -(stop-Pos_approx)]
lof_all_filtered_core_CDS[strand=="-",Pos_gene:=full_length -Pos_gene]
lof_all_filtered_core_CDS[,Pos_gene_perc:=Pos_gene/full_length]
# png("./Results/Histogram of distribution of the LoF alleles in a gene.png",width= 3000, height= 1500, res=300)
# hist(lof_all_filtered_core_CDS[,Pos_gene_perc],breaks=50)
# dev.off()
lof_all_filtered_core_CDS<-lof_all_filtered_core_CDS[Pos_gene_perc<=0.95,]
lof_all_filtered_core_CDS[,c("i.gene_ID","core","start_join","stop_join","LG"):=NULL]


# for the ones that are within 20  bps, only keep one with the largest alternative allele depth
setorderv(lof_all_filtered_core_CDS,c("lake","LGn","Pos"))
lof_all_filtered_core_CDS[,Grouping:=Pos]
lof_all_filtered_core_CDS[,Pos_lead := shift(Pos,1L,type="lead")]
lof_all_filtered_core_CDS[,Pos_lag := shift(Pos,1L,type="lag")]
lof_all_filtered_core_CDS[ (Pos-Pos_lag>0 & Pos-Pos_lag<20) | (Pos_lead-Pos_lag>0 &Pos_lead-Pos<20), Grouping:= 0L]
lof_all_filtered_core_CDS<-lof_all_filtered_core_CDS[lof_all_filtered_core_CDS[,.I[which.max(DP_ALT_f+DP_ALT_r)],by=rleid(Grouping)]$V1]
lof_all_filtered_core_CDS[,c("Pos_lead","Pos_lag"):=NULL]


# find the list of the genes which frameshif mutation has been corrected by frame restoring indels
temp<-lof_all_filtered_core_CDS[frameshift==T,.(No_frameshift=sum(frameshift==T),Pos=Pos,Pos_gene=Pos_gene,conseq=conseq,ALF=ALF,frameshift=frameshift,REF=REF,ALT=ALT,Pos_gene),by=c("lake","gene")][No_frameshift>1]
temp[,ref_bp:=nchar(REF)]
temp[,alt_bp:=nchar(ALT)]
temp[grepl(",",ALT),alt_bp:=NA]
frame_restored_loci[,gap_length:=max(Pos_gene)-min(Pos_gene),by=.(lake,gene)]
frame_restored_loci<-temp[,diff:=if(!anyNA(alt_bp)) abs(sum(ref_bp)-sum(alt_bp)) ,by=c("lake","gene")][diff%%3==0,]
lof_all_filtered_core_CDS<-lof_all_filtered_core_CDS[!frame_restored_loci,on=.(lake,gene,Pos)]
rm(temp)

# create a by locus data.table with allele frequency of each lake as columns
lof_all_filtered_core_CDS[,paste0("ALF_", files):=lapply(files, function(i) ALF*(as.integer(lake==i)))]
lof_all_filtered_core_CDS[,(39:46):=lapply(.SD, max),.SDcols=39:46, by=list(LGn, Pos)]
lof_all_filtered_bylocus<-unique(lof_all_filtered_core_CDS, by=c("LGn","Pos"))
lof_all_filtered_bylocus<-lof_all_filtered_bylocus[,c("LGn" ,          "gene_ID",   "gene_name",    "start",         "stop",          "strand",       
                                                      "transcript_ID", "full_length",   "total_exon",    "Pos",  "Pos_gene_perc",        "REF",           
                                                      "ALT",          "conseq",            "IDV",                    
                                                      "DP_ALT_f",      "DP_ALT_r" ,     "DP4_sum" ,"ALF",      
                                                      "ALF_boot",      "ALF_echo", "ALF_fred", "ALF_gos", "ALF_law", "ALF_pach","ALF_rob","ALF_say"  )] 

# create a by gene data.table, the allele frequency is filled with max value from any LOF alleles in that gene
temp<-lof_all_filtered_bylocus
temp[,(20:27):=lapply(.SD, max),.SDcols=20:27, by=list(LGn, gene_ID)]
lof_all_filtered_bylocus[,.N,by=c("LGn","gene_ID")]
lof_all_filtered_bygene<-unique(temp, by=c("LGn","gene_ID"))
lof_all_filtered_bygene<-lof_all_filtered_bygene[lof_all_filtered_bylocus[,.N,by=c("LGn","gene_ID")],on=c("LGn","gene_ID")]
lof_all_filtered_bygene[,min_ALF:=pmin(ALF_boot,ALF_echo,ALF_fred,ALF_gos,ALF_law,ALF_pach,ALF_rob,ALF_say)]
lof_all_filtered_bygene[,max_ALF:=pmax(ALF_boot,ALF_echo,ALF_fred,ALF_gos,ALF_law,ALF_pach,ALF_rob,ALF_say)]
lof_all_filtered_bygene[, col_min:= colnames(.SD)[max.col(-.SD, ties.method = "first")],.SDcols=20:27]
lof_all_filtered_bygene[, col_max:= colnames(.SD)[max.col(.SD, ties.method = "first")],.SDcols=20:27]
lof_all_filtered_bygene[, No_pres_pop:=rowSums((sapply(.SD, function(x) x>0))),.SDcols=20:27]
rm(temp)


# rob has the most pop specific lof alleles
png("./Results/Number of pop specific LoF genes per population.png",width=3000,height=1500,res=300)
barplot(table(lof_all_filtered_bygene[No_pres_pop==1, col_max])) 
dev.off()

lof_all_filtered_bygene[No_pres_pop==8,.N] # 989 lof genes are shared by all pops
lof_all_filtered_bygene[No_pres_pop==7 & col_min=="ALF_say",gene_ID] # 72 lof genes are only present in freshwater lakes
lof_all_filtered_bygene[No_pres_pop==1 & col_max=="ALF_say",gene_ID] # 111 lof genes are only present in say



rob0.6_gos0_lof<-lof_all_filtered_bygene[ALF_gos==0 & ALF_rob>0.8,]
gos0.6_rob0_lof<-lof_all_filtered_bygene[ALF_rob==0 & ALF_gos>0.8,]

fwrite(rob0.6_gos0_lof,"./Results/rob0.6_gos0_lof.csv")
fwrite(gos0.6_rob0_lof,"./Results/gos0.6_rob0_lof.csv")

### validation with Brian's Rob and Gos expression data
expression_rob_gos_population<-fread("./Data/Lohman_expression/Love_Population.csv")

## analysis of unnormalized raw counts
# counts<-fread("./Data/Lohman_expression/98_counts.csv")
# design<-fread("./Data/Lohman_expression/98_rnaDesign.csv")
# select<-design[V5 == "Control", V1]
# select_rob<-select[unlist(lapply(select, grepl,pattern="^RR"))]
# select_gos<-select[unlist(lapply(select, grepl,pattern="^GG"))]
# select_counts<-counts[,c(..select,"V1")]
# select_counts[,RR:=rowMeans(.SD),.SDcols=select_rob]
# select_counts[,GG:=rowMeans(.SD),.SDcols=select_gos]
# select_counts<-select_counts[,.(V1,RR,GG)]
# select_counts[,log2fold:=log2(RR/GG)]

rob0.6_gos0_lof_expression<-expression_rob_gos_population[rob0.6_gos0_lof,on=c("V1"="gene_ID")][!is.na(baseMean),]
gos0.6_rob0_lof_expression<-expression_rob_gos_population[gos0.6_rob0_lof,on=c("V1"="gene_ID")][!is.na(baseMean),]
    
t.test(rob0.6_gos0_lof_expression[,stat],gos0.6_rob0_lof_expression[,stat])
t.test(rob0.6_gos0_lof_expression[,log2FoldChange],gos0.6_rob0_lof_expression[,log2FoldChange])
boxplot(rob0.6_gos0_lof_expression[pvalue<0.05,log2FoldChange],gos0.6_rob0_lof_expression[pvalue<0.05,log2FoldChange])
boxplot(rob0.6_gos0_lof_expression[,log2FoldChange],gos0.6_rob0_lof_expression[,log2FoldChange])

### run a simulation to draw genes with no lof mutations
no_lof_gene<-gene_info[!lof_all_filtered_bygene,on=c("GeneID"="gene_ID")]
no_lof_gene_expressed<-expression_rob_gos_population[no_lof_gene,on=c("V1"="GeneID"),nomatch=0L]
length_t<-dim(no_lof_gene_expressed)[1]
length_rob<-dim(rob0.6_gos0_lof_expression)[1]
length_gos<-dim(gos0.6_rob0_lof_expression)[1]
t_storage<-vector(mode = "list", length = 5000)
for(i in 1:5000 ){
  rob_control<-no_lof_gene_expressed[sample(1:length_t,length_rob,replace=F)]
  gos_control<-no_lof_gene_expressed[sample(1:length_t,length_rob,replace=F)]
  if(dim(rob_control[gos_control,on="V1",nomatch=0L])[1]>0) next # make sure there is no overlap between the two control groups
  t_storage[i]<-t.test(rob_control[,stat],gos_control[,stat])$statistic[["t"]]
}
hist(unlist(t_storage),breaks=100)
quantile(unlist(t_storage),probs = c(0.05,0.95))

### overlapping with population genetics analysis
PBS_slide_20kb<-fread("../Analyses_Poolseq/Results/Foen/PBS_slide_20kb.csv")
PBS_slide_20kb[,GRFst_ranking:=rank(-GRFst, na.last=T,ties.method = "random" )]
PBS_slide_500snp<-fread("../Analyses_Poolseq/Results/Foen/PBS_slide_500snp.csv")
PBS_slide_500snp[,GRFst_ranking:=rank(-GRFst, na.last=T,ties.method = "random" )]

lof_20kb<-PBS_slide_20kb[GRFst.focal==T,][lof_all_filtered_bygene,on=c("LGn","window_start<=Pos","window_end>=Pos"),nomatch=0L]
lof_20kb<-lof_20kb[lof_20kb[,.I[which.max(GRFst)],by=.(LGn, gene_ID)]$V1]
lof_20kb<-lof_20kb[,RG_AFD:=ALF_rob-ALF_gos][abs(RG_AFD)>0.3 & (ALF_rob==0 | ALF_gos==0)]
setorder(lof_20kb,"GRFst_ranking")

lof_500snp<-PBS_slide_500snp[GRFst.focal==T,][lof_all_filtered_bygene,on=c("LGn","Pos_start<=Pos","Pos_end>=Pos"),nomatch=0L]
lof_500snp<-lof_500snp[lof_500snp[,.I[which.max(GRFst)],by=.(LGn, gene_ID)]$V1]
lof_500snp<-lof_500snp[,RG_AFD:=ALF_rob-ALF_gos][abs(RG_AFD)>0.3 & (ALF_rob==0 | ALF_gos==0)]
setorder(lof_500snp,"GRFst_ranking")

fwrite(lof_500snp,"./Results/lof_500snpGRFst0.99_overlap.csv")
fwrite(lof_20kb,"./Results/lof_20kbGRFst0.99_overlap.csv")

expression_rob_gos_population[lof_20kb,on=c("V1"="gene_ID"),nomatch=0L][,keyby=GRFst]
expression_rob_gos_population[lof_500snp,on=c("V1"="gene_ID"),nomatch=0L][,keyby=GRFst]

# rob0.6_gos0_lof<-lof_all_filtered_bygene[ALF_gos==0 & ALF_rob>0.6,]
# gos0.6_rob0_lof<-lof_all_filtered_bygene[ALF_rob==0 & ALF_gos>0.6,]
# 
# rob0.6_gos0_lof_PBS<-PBS_slide_20kb[rob0.6_gos0_lof,on=c("LGn","window_start<=Pos","window_end>=Pos")]
# rob0.6_gos0_lof_PBS<-rob0.6_gos0_lof_PBS[rob0.6_gos0_lof_PBS[,.I[which.max(GRFst)],by=.(LGn, window_start)]$V1]
# gos0.6_rob0_lof_PBS<-PBS_slide_20kb[gos0.6_rob0_lof,on=c("LGn","window_start<=Pos","window_end>=Pos")]
# gos0.6_rob0_lof_PBS<-gos0.6_rob0_lof_PBS[gos0.6_rob0_lof_PBS[,.I[which.max(GRFst)],by=.(LGn, window_start)]$V1]


