Sys.setlocale("LC_ALL","English_United States.utf8")
args<-commandArgs(trailingOnly=TRUE);worker<-if(length(args)>0)as.integer(args[1]) else 0L;workers<-if(length(args)>1)as.integer(args[2]) else 1L
root<-normalizePath(".",winslash="/");.libPaths(c(file.path(root,"Rlib"),.libPaths()))
Sys.setenv(OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1")
suppressPackageStartupMessages({library(affy);library(frma);library(hgu133plus2cdf);library(hgu133plus2.db);library(AnnotationDbi)})
lock<-jsonlite::read_json("config/expression_preprocessing_lock.json",simplifyVector=TRUE)
members<-read.csv("metadata/new_NAM_membership_and_split.csv",colClasses="character",fileEncoding="UTF-8-BOM",check.names=FALSE)
members<-members[(seq_len(nrow(members))-1L)%%workers==worker,,drop=FALSE]
selection<-read.csv("metadata/PHH_download_plan.csv",colClasses="character",fileEncoding="UTF-8-BOM",check.names=FALSE)
outdir<-file.path(root,"derived/expression_by_parent");dir.create(outdir,showWarnings=FALSE)
probes<-ls(hgu133plus2cdf)
ann<-AnnotationDbi::select(hgu133plus2.db,keys=intersect(probes,keys(hgu133plus2.db,keytype="PROBEID")),columns=c("ENTREZID","SYMBOL"),keytype="PROBEID")
count<-tapply(ann$ENTREZID,ann$PROBEID,function(v)length(unique(v[!is.na(v)])))
ann$unique_gene_mapping<-unname(count[ann$PROBEID])==1
ann$mapping_reason<-ifelse(ann$unique_gene_mapping,"one_unique_Entrez_gene",ifelse(unname(count[ann$PROBEID])==0,"unmapped","ambiguous_multi_gene_excluded"))
if(worker==0L)write.csv(ann,"metadata/PHH_probe_gene_mapping.csv",row.names=FALSE)
map<-unique(ann[ann$unique_gene_mapping&!is.na(ann$ENTREZID),c("PROBEID","ENTREZID")]);stopifnot(!anyDuplicated(map$PROBEID))
genes<-sort(unique(map$ENTREZID));geneindices<-split(match(map$PROBEID,probes),map$ENTREZID)
h<-read.csv("metadata/Hallmark_members_original.csv",colClasses="character",check.names=FALSE)
sets<-split(h$db_ncbi_gene,h$gs_name);sets<-lapply(sets,unique)
support<-data.frame(geneset=names(sets),original_unique_genes=vapply(sets,length,integer(1)),mapped_genes=vapply(sets,function(x)length(intersect(x,genes)),integer(1)))
support$fraction<-support$mapped_genes/support$original_unique_genes;support$eligible<-support$mapped_genes>=10&support$fraction>=.2
if(worker==0L)write.csv(support,"metadata/Hallmark_platform_support.csv",row.names=FALSE)
sets<-lapply(sets[support$eligible],function(x)intersect(x,genes));if(worker==0L)write.csv(data.frame(geneset=rep(names(sets),lengths(sets)),Entrez_gene=unlist(sets,use.names=FALSE)),"metadata/Hallmark_fixed_mapped_members.csv",row.names=FALSE)
if(worker==0L)jsonlite::write_json(list(saved_before_model_fitting=TRUE,method=lock$method,hard_QC=c("complete source selected treatment/control CEL", "chip/CDF equals U133Plus2", "all raw intensities finite and nonnegative with nonzero range", "finite normalized probes and correct IDs", "unique gene mapping only"),repeat_correlation="report_only_not_outcome_based_exclusion",viability="metadata_only_not_exclusion",geneset_rule="platform unique mapping >=10 genes and >=20 percent coverage",genesets=names(sets)),"config/expression_QC_and_sets_lock.json",pretty=TRUE,auto_unbox=TRUE)
writegz<-function(mat,idname,path){co<-gzfile(path,"wt");on.exit(close(co));d<-data.frame(id=rownames(mat),mat,check.names=FALSE);names(d)[1]<-idname;write.csv(d,co,row.names=FALSE)}
process_one<-function(pid,d) {
 t0<-proc.time()[3];expected<-selection[selection$parent_id==pid,]
 if(nrow(d)!=nrow(expected)||!setequal(d$sample_id,expected$sample_id))stop("Incomplete acquisition row mapping")
 if(any(d$acquisition_status!="available"))stop("Selected highest/control CEL unavailable; no lower concentration replacement")
 ab<-ReadAffy(filenames=d$CEL_path,verbose=FALSE)
 if(annotation(ab)!="hgu133plus2"||cdfName(ab)!="HG-U133_Plus_2")stop("Incompatible actual chip or CDF")
 raw<-intensity(ab);if(any(!is.finite(raw))||any(raw<0)||any(apply(raw,2,function(x)diff(range(x)))==0))stop("Hard raw intensity QC failure")
 raw_qc<-data.frame(BARCODE=d$BARCODE,role=d$role,raw_min=apply(raw,2,min),raw_median=apply(raw,2,median),raw_max=apply(raw,2,max));rm(raw)
 if(lock$method=="compatible_fixed_reference_fRMA") {
   pe<-sapply(seq_along(sampleNames(ab)),function(j)exprs(frma(ab[,j,drop=FALSE]))[geneNames(ab),1]);rownames(pe)<-geneNames(ab)
 } else pe<-exprs(rma(ab,verbose=FALSE))
 if(!setequal(rownames(pe),probes)||any(!is.finite(pe)))stop("Invalid normalized output IDs/values")
 pe<-pe[probes,,drop=FALSE];colnames(pe)<-d$BARCODE
 ge<-t(vapply(geneindices,function(ix)apply(pe[ix,,drop=FALSE],2,median),numeric(ncol(pe))))
 colnames(ge)<-colnames(pe)
 ti<-which(d$role=="treated");ci<-which(d$role=="control");stopifnot(length(ti)>0,length(ci)>0)
 response<-rowMeans(ge[,ti,drop=FALSE])-rowMeans(ge[,ci,drop=FALSE]);z<-vapply(sets,function(ix)mean(response[ix]),numeric(1));stopifnot(all(is.finite(z)))
 prefix<-file.path(outdir,pid)
 writegz(pe,"PROBEID",paste0(prefix,"_probe_log2.csv.gz"));writegz(ge,"Entrez_gene",paste0(prefix,"_gene_log2.csv.gz"))
 write.csv(data.frame(Entrez_gene=names(response),response=as.numeric(response)),paste0(prefix,"_gene_response.csv"),row.names=FALSE)
 write.csv(data.frame(geneset=names(z),response=as.numeric(z)),paste0(prefix,"_geneset_response.csv"),row.names=FALSE)
 write.csv(raw_qc,paste0(prefix,"_raw_QC.csv"),row.names=FALSE)
 cors<-cor(pe);write.csv(cors,paste0(prefix,"_probe_repeat_correlations.csv"))
 list(parent_id=pid,status="success",source_EXP_ID=d$EXP_ID[1],arrays=nrow(d),probe_sets=nrow(pe),unique_genes=nrow(ge),gene_sets=length(z),n_treated=length(ti),n_controls=length(ci),method=lock$method,elapsed_seconds=unname(proc.time()[3]-t0),completed_utc=format(Sys.time(),tz="UTC",usetz=TRUE),prediction_generated=FALSE)
}
done<-character();start<-Sys.time()
repeat {
 for(pid in members$parent_id[!members$parent_id%in%done]) {
   statusfile<-file.path(outdir,paste0(pid,"_status.json"))
   if(file.exists(statusfile)){done<-c(done,pid);next}
   file<-file.path(root,"metadata/archive_directories",paste0(pid,"_acquired.csv"))
   if(!file.exists(file))next
   d<-tryCatch(read.csv(file,colClasses="character",fileEncoding="UTF-8-BOM",check.names=FALSE),error=function(e)NULL)
   if(is.null(d)||nrow(d)!=sum(selection$parent_id==pid))next
   cat("PROCESS",pid,"\n");flush.console()
   result<-tryCatch(process_one(pid,d),error=function(e)list(parent_id=pid,status="failed",error=conditionMessage(e),prediction_generated=FALSE,completed_utc=format(Sys.time(),tz="UTC",usetz=TRUE)))
   jsonlite::write_json(result,statusfile,pretty=TRUE,auto_unbox=TRUE);done<-c(done,pid);cat("DONE",pid,result$status,"\n");flush.console();gc()
 }
 if(length(done)==nrow(members))break
 if(as.numeric(difftime(Sys.time(),start,units="hours"))>3)stop("Bounded expression wait expired; actual partial outputs retained")
 Sys.sleep(5)
}
writeLines(capture.output(sessionInfo()),paste0("logs/R_expression_session_worker",worker,".txt"))
cat("Finished actual expression processing",length(done),"parents\n")
