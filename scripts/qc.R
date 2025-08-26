#!/usr/bin/env Rscript

library(dplyr)

##working on collect all QC in one single file
##replace the rRNA in qc53.txt
options(stringsAsFactors=FALSE)


checknames<-function(x,title) {
  if(sum(rownames(x)!=samples)==0) return(1)
  #to fix the problem of multqc that the order of S1 is later than S10
  if(sum(sort(rownames(x))!=sort(samples))==0) return(2)
  stop(paste0(title," is wrong"))
}


myTryCatch <- function(expr) {
  warn <- err <- NULL
  value <- withCallingHandlers(
    tryCatch(expr, error=function(e) {
      err <<- e
      NULL
    }), warning=function(w) {
      warn <<- w
      invokeRestart("muffleWarning")
    })
  list(value=value, warning=warn, error=err)
}


collect.qc <- function(snakemake.params){
  
  print(paste("read star qc info"))
  file_path<-snakemake.params[["star_align_merge_all"]]
  if (!file.exists(file_path)) {
    print (paste("The specified file does not exist - ",file_path))
    #return(NULL)
  }else{
    #read star qc info the multiqc misses chimeric%, go with manully collected data from Log.final.out 
    star <- t(read.delim(file_path,sep="\t",row.names=1,strip.white=TRUE,check.names=FALSE))
    star<-star[,c(5,6,8:16,24,26,28,29,30,33),drop=FALSE]
    colnames(star)<-c("reads","avg_input_read_length","uniquely_mapped","%uniquely_mapped","avg_mapped_read_length",
                      "num_splices","num_annotated_splices","num_GTAG_splices","num_GCAG_splices","num_ATAC_splices",
                      "num_noncanonical_splices","%multimapped","%multimapped_toomany","%unmapped_mismatches","%unmapped_tooshort",
                      "%unmapped_other","%chimeric")
    star<-sub("%$","",star)
    rownames(star)<-sub(" \\|$","",sub("^star_align \\| ","",rownames(star)))
    
    star<-star[sort(rownames(star)),,drop=FALSE]
    print(paste("star df dims: ","nrow - ",nrow(star), "ncol - ",ncol(star)))
    # assertthat::are_equal(rownames(star), samples)
  }
  
  #----------------------------------------------------------------------------------------------
  
  
  
  print(paste("read fastqc info for fastqc_raw"))
  #read fastqc info for fastqc_raw
  file_path<-snakemake.params[["multiqc_fastqc_raw_txt"]]
  if (!file.exists(file_path)) {
    print (paste("The specified file does not exist - ",file_path))
    #return(NULL)
  }else{
    
    fastqc_raw<-read.delim(file_path,sep="\t",head=TRUE,row.names=1,check.names=FALSE,strip.white=TRUE)
    
    fastqc_raw<- fastqc_raw[,c("Total Sequences","%GC","total_deduplicated_percentage")]
    All_single<-TRUE
    if(length(grep("_R2$",rownames(fastqc_raw)))>0){
      #some are paired
      All_single<-FALSE
      fastqc_raw$name<- row.names(fastqc_raw)
      fastqc_raw$name <- sub("_R2$","_R1",fastqc_raw$name)
      fastqc_raw<-aggregate(cbind(`Total Sequences`,`%GC`,`total_deduplicated_percentage`) ~ name, data = fastqc_raw, FUN = mean, na.rm = TRUE)
      row.names(fastqc_raw)<-fastqc_raw$name
      fastqc_raw<-fastqc_raw[,-grep("name",names(fastqc_raw))]
    }
    rownames(fastqc_raw)<-sub(".* fastqc_raw \\| (.*)_R1$","\\1",rownames(fastqc_raw))
    
    samples <- snakemake.params[["samples"]]
    fastqc_raw <- fastqc_raw [rownames(fastqc_raw)%in%samples,] 
    
    print(paste("fastqc_raw df dims: ","nrow - ",nrow(fastqc_raw), "ncol - ",ncol(fastqc_raw)))
    rm(All_single)
  }
  #----------------------------------------------------------------------------------------------
  
  
  
  print(paste("read clean-up fastqc"))  
  #read clean-up fastqc
  file_path<-snakemake.params[["multiqc_fastqc_trim_txt"]]
  if (!file.exists(file_path)) {
    print (paste("The specified file does not exist - ",file_path))
    #return(NULL)
  }else{
    
    fastqc<-read.delim(file_path,sep="\t",head=TRUE,row.names=1,check.names=FALSE,strip.white=TRUE)
    fastqc<- fastqc[,c("Total Sequences","%GC","total_deduplicated_percentage")]
    if(length(grep("_R2$",rownames(fastqc)))>0){
      #some are paired
      All_single<-FALSE
      fastqc$name<- row.names(fastqc)
      fastqc$name <- sub("_R2$","_R1",fastqc$name)
      fastqc<-aggregate(cbind(`Total Sequences`,`%GC`,`total_deduplicated_percentage`) ~ name, data = fastqc, FUN = mean, na.rm = TRUE)
      row.names(fastqc)<-fastqc$name
      fastqc<-fastqc[,-grep("name",names(fastqc))]
    }
    rownames(fastqc)<-sub(".*fastqc_trim \\| (.*)_R1","\\1",rownames(fastqc))
    samples <- snakemake.params[["samples"]]
    fastqc <- fastqc [rownames(fastqc)%in%samples,] 
    print(paste("multiqc_fastqc_trim_txt df dims: ","nrow - ",nrow(fastqc), "ncol - ",ncol(fastqc)))
  }
  #----------------------------------------------------------------------------------------------
  
  
  
  print(paste("read the trim info"))    
  file_path<-snakemake.params[["multiqc_fastqc_cutadapt_txt"]]
  if (!file.exists(file_path)) {
    print (paste("The specified file does not exist - ",file_path))
    #return(NULL)
  }else{
    #read the trim info
    trim <-read.delim(file_path,sep="\t",head=TRUE,row.names=1,check.names=FALSE,strip.white=TRUE)
    R1R2_present<-snakemake.params[["samples_R1R2_present"]]
    # rownames(trim)<-sub(".*tool_logs \\| (.*)_R1$","\\1",rownames(trim))  # orignal code 
    if (R1R2_present){
      rownames(trim)<-sub(".*tool_logs \\| (.*?)_R2", "\\1", rownames(trim)) # SR: searches for the last appearance of the "tool_logs |" that followed by R2 and replaced by the 1st captured group
    }else{
      rownames(trim)<-sub(".*tool_logs \\| (.*?)_R1", "\\1", rownames(trim)) # SR: searches for the last appearance of the "tool_logs |" that followed by R1 and replaced by the 1st captured group
    }
    samples <- snakemake.params[["samples"]]
    trim <- trim[rownames(trim)%in%samples,] 
    
    print(paste("trim df dims: ","nrow - ",nrow(trim), "ncol - ",ncol(trim)))
  }
  #----------------------------------------------------------------------------------------------
  
  
  
  print(paste("read qc53"))    
  ##read qc53,note that the multiqc is not working well for collectRNAmetrics
  file_path<-snakemake.params[["post_align_qc53_txt"]]
  if (!file.exists(file_path)) {
    print (paste("The specified file does not exist - ",file_path))
    #return(NULL)
  }else{
    
    qc53<-t(read.delim(file_path,sep="\t",head=TRUE,row.names=1,check.names=FALSE))
    qc53<-qc53[,-match(c("RIBOSOMAL_BASES","PCT_RIBOSOMAL_BASES","SAMPLE","LIBRARY","READ_GROUP"),
                       colnames(qc53)),drop=FALSE]
    Nqc53<-ncol(qc53)
    nqc53<-nrow(qc53)
    
    samples <- snakemake.params[["samples"]]
    NS<-length(samples)
    if(nqc53==NS){
      qc53<-qc53[samples,,drop=FALSE]
    }else{
      stop("star_qc and qc53 do not match")
    }
    
    
    id<-c("CODING","UTR","INTRONIC","INTERGENIC","MRNA")
    loc<-match(c(paste0("PCT_",id,"_BASES"),"MEDIAN_5PRIME_TO_3PRIME_BIAS"),
               colnames(qc53)[1:Nqc53])
    qc53<-qc53[,loc,drop=FALSE]
    id2<-tolower(id)
    colnames(qc53)[1:6]<-c(paste0("%",id2),"median_5_3_bias")
    
    id<-grep("^%",colnames(qc53))
    qc53[,id]<-qc53[,id]*100
    qc53<-round(qc53,dig=3)
    rm(id,loc,id2,NS,Nqc53,nqc53)
    print(paste("qc53 df dims: ","nrow - ",nrow(qc53), "ncol - ",ncol(qc53)))
  }
  #----------------------------------------------------------------------------------------------
  
  
  
  print(paste("collect data for rRNA"))   
  ##collect data for rRNA, globin and phix and duplicates when present
  UMI <- snakemake.params[["samples_R1I1_present"]]
  samples <- snakemake.params[["samples"]]
  NS<-length(samples)
  misc<-matrix(NA,NS,6)
  colnames(misc)<-c("globin","rRNA","phix","picard_dup","UMI_dup","adapter_detected")
  for(i in 1:NS){
    SID<-samples[i]
    
    file_path<-snakemake.params[["globin"]]
    if (!file.exists(file_path)) {
      print (paste("The specified file does not exist - ",file_path))
      #return(NULL)
    }else{
      
      file <- file.path(file_path,paste0(SID,".txt"))
      v<-tail(readLines(file))
      v <-gsub("(.*)\\%.*","\\1",v[6])
      misc[i,"globin"]<-v
    }
    
    file_path<-snakemake.params[["rRNA"]]
    if (!file.exists(file_path)) {
      print (paste("The specified file does not exist - ",file_path))
      #return(NULL)
    }else{
      
      file <- file.path(file_path,paste0(SID,".txt"))
      v<-tail(readLines(file))
      v <-gsub("(.*)\\%.*","\\1",v[6])
      misc[i,"rRNA"]<-v
    }
    
    file_path<-snakemake.params[["phix"]]
    if (!file.exists(file_path)) {
      print (paste("The specified file does not exist - ",file_path))
      #return(NULL)
    }else{
      
      file <- file.path(file_path,paste0(SID,".txt"))
      v<-tail(readLines(file))
      v <-gsub("(.*)\\%.*","\\1",v[6])
      misc[i,"phix"]<-v
    }
    
    
    
    file_path<-snakemake.params[["mark_dup_metrics"]]
    if (!file.exists(file_path)) {
      print (paste("The specified file does not exist - ",file_path))
      #return(NULL)
    }else{
      
      misc[i,"picard_dup"]<-read.delim(file.path(file_path,paste0(SID,".dup_metrics")),skip=6,head=TRUE,row=1)[1,"PERCENT_DUPLICATION"]*100
    }
    
    
    
    if(UMI){
      file_path<-snakemake.params[["umi_dup"]]
      if (!file.exists(file_path)) {
        print (paste("The specified file does not exist - ",file_path))
        #return(NULL)
      }else{
        
        file <- file.path(file_path,SID,"dup_log.txt")
        v<-tail(readLines(file))
        v <-gsub(".*\\t(.*)","\\1",v[2])
        misc[i,"UMI_dup"]<-as.numeric(v)*100
      }
    }
    
    file_path<-snakemake.params[["trim_log"]]
    if (!file.exists(file_path)) {
      print (paste("The specified file does not exist - ",file_path))
      #return(NULL)
    }else{
      
      file <- file.path(file_path,paste0("log.",SID))
      v<-readLines(file)
      v<-grep("with adapter",v,value = T)
      r1<-gsub(".*\\((.*)%\\)$","\\1",v[1])
      r2<-gsub(".*\\((.*)%\\)$","\\1",v[2])
      ##future plan is to consider the paired info for each indivudal files
      #if(length(zval)!= PAIRED+1){
      #    stop("the fastq_trim log for the contained adapter% is not consistent with the pairedness of the fastq data for sample", SID)
      #}
      misc[i,"adapter_detected"]<-mean(c(as.numeric(r1),as.numeric(r2)))
    }
  }
  rm(NS,samples,v,UMI,i,SID,file)
  
  # misc <- as.data.frame(apply(misc, 2,function(x) as.numeric(as.character(x)))) #re-implement
  # misc<-round(misc,dig=3)
  colnames(misc)<-paste0("%",colnames(misc))
  print(paste("misc df dims: ","nrow - ",nrow(misc), "ncol - ",ncol(misc)))
  
  #----------------------------------------------------------------------------------------------
  
  
  
  
  print(paste("Read the chr_info.txt"))    
  ##Read the chr_info.txt
  samples <- snakemake.params[["samples"]]
  NS<-length(samples)
  chr_info<-matrix(NA,NS,5)
  colnames(chr_info)<-c("chrX","chrY","chrM","chrAuto","contig")
  
  readchr<-function(sid){
    file_path<-snakemake.params[["chr_info"]]
    if (!file.exists(file_path)) {
      print (paste("The specified file does not exist - ",file_path))
      #return(NULL)
    }else{
      
      x<-read.delim(paste0(file_path,"/",sid,"/chr_info.txt"),
                    sep="\t",row.names=1,header=FALSE)
      xt<-sum(x[,2])
      y<-x[c("chrX","chrY","chrM"),2]/xt*100
      y<-c(y,sum(x[grep("^chr[1-9]",rownames(x)),2])/xt*100)
      c(y,100-sum(y))
    }
  }
  
  for(i in 1:NS){
    chr_info[i,1:5]<-readchr(samples[i])
  }
  colnames(chr_info)<-paste0("%",colnames(chr_info))
  ##keep chry with four digits
  chr_info[,2]<-round(chr_info[,2],dig=5)
  chr_info[,-2]<-round(chr_info[,-2],dig=3)
  rm(NS,samples,readchr,i,file_path)
  print(paste("chr_info df dims: ","nrow - ",nrow(chr_info), "ncol - ",ncol(chr_info)))
  #----------------------------------------------------------------------------------------------
  
  
  
  
  print(paste("Now putting all together"))    
  #Now putting all together
  
  print(paste("First convert all individual qc tables ro DFs"))    
  
  fastqc_raw <- as.data.frame(fastqc_raw)
  fastqc <- as.data.frame(fastqc)
  chr_info <- as.data.frame(chr_info)
  misc <- as.data.frame(misc)
  qc53 <- as.data.frame(qc53)
  star <- as.data.frame(star)
  trim <- as.data.frame(trim)
  fastqc_raw <- as.data.frame(fastqc_raw)
  
  print(paste("Done. Now combind them"))    
  
  qc<-NULL
  # qc<-cbind("reads_raw"=fastqc_raw[,1,drop=FALSE],"%adapter_detected"=misc[,6,drop=FALSE],
  #           "%trimmed"=round(100-as.numeric(star[,1,drop=FALSE])/fastqc_raw[,"Total Sequences",drop=FALSE]*100,dig=3),
  #           "%trimmed_bases"=round(trim[,"percent_trimmed",drop=FALSE],dig=3),"reads"=star[,1,drop=FALSE],"%GC"=round(fastqc[,"%GC",drop=FALSE],dig=3),
  #           "%dup_sequence"=round(100-fastqc[,"total_deduplicated_percentage",drop=FALSE],dig=3),misc[,1:4,drop=FALSE])
  
  qc<-bind_cols(
    "reads_raw"=fastqc_raw[,1,drop=FALSE],
    
    "%adapter_detected"=misc[,6,drop=FALSE],
    
    "%trimmed"=round(100- data.frame(Col1 = as.numeric(star[[1]]))/fastqc_raw[,"Total Sequences",drop=FALSE] *100,dig=3),
    
    "%trimmed_bases"=round(trim[,"percent_trimmed",drop=FALSE],dig=3),"reads"=star[,1,drop=FALSE],"%GC"=round(fastqc[,"%GC",drop=FALSE],dig=3),
    
    "%dup_sequence"=round(100-fastqc[,"total_deduplicated_percentage",drop=FALSE],dig=3),misc[,1:4,drop=FALSE]
  )
  
  if(snakemake.params[["samples_R1I1_present"]]){
    # qc<-cbind(qc,"%umi_dup"=misc[,5,drop=FALSE])
    qc<-bind_cols(qc,"%umi_dup"=misc[,5,drop=FALSE])
  }
  # qc<-cbind(qc,star[,-1,drop=FALSE],chr_info,qc53)
  qc<-bind_cols(qc,star[,-1,drop=FALSE],chr_info,qc53)
  
  print(paste("Done. Finalizing and returning the combine QC table"))    
  
  id<-grep("_percent$",colnames(qc))
  colnames(qc)[id]<-paste0("%",sub("_percent$","",colnames(qc)[id]))
  colnames(qc)<-sub("^%","pct_",colnames(qc))
  
  return(qc)
}







if(exists("snakemake")){
  print(paste("the 'snakemake' object is available in the ENV","script invoked by snakemake"))
}else{
  print(paste("the 'snakemake' object is NOT evailable in the ENV","script invoked manually"))
  dir<-"/sc/arion/projects/sealfs01a/sealfonlab_seq_projects/processed_data/Public/rnaseq/Sepsis_GSE185263"
  print(paste("setting run folder to: ",dir))
  setwd(dir)
  # setwd("/sc/arion/projects/sealfs01a/sealfonlab_seq_projects/processed_data/ECHO/rnaseq/ECRE-Radiation-bulkRNAseq_198_202407241745_SR1")
  snakemake <- readRDS(file = "snakemake.RData")
}


log_file <- snakemake@log[[1]]
sink(log_file)

print(paste("log_file",log_file))
print(paste("initiating","snakemake.params"))

snakemake.params <- list()
snakemake.params[["multiqc_fastqc_raw_txt"]] <- snakemake@input[["multiqc_fastqc_raw_txt"]]
snakemake.params[["multiqc_fastqc_trim_txt"]] <- snakemake@input[["multiqc_fastqc_trim_txt"]]
snakemake.params[["multiqc_fastqc_cutadapt_txt"]] <- snakemake@input[["multiqc_fastqc_cutadapt_txt"]]
snakemake.params[["samples"]] <- snakemake@params[["samples"]]
snakemake.params[["pipeline_warning_file_path"]] <- snakemake@params[["pipeline_warning_file_path"]]
snakemake.params[["samples_R1I1_present"]] <- snakemake@params[["samples_R1I1_present"]]
snakemake.params[["samples_R1R2_present"]] <- snakemake@params[["samples_R1R2_present"]]
snakemake.params[["phix"]] <- snakemake@params[["phix"]]
snakemake.params[["chr_info"]] <- snakemake@params[["chr_info"]]
snakemake.params[["trim_log"]] <- snakemake@params[["trim_log"]]

snakemake.params[["bismark_summary_report"]] <- snakemake@input[["bismark_summary_report"]]
snakemake.params[["bismark_lambda_summary_report"]] <- snakemake@input[["bismark_lambda_summary_report"]]
snakemake.params[["bismark_4strand"]] <- snakemake@input[["bismark_4strand"]]
snakemake.params[["bismark_lambda_4strand"]] <- snakemake@input[["bismark_lambda_4strand"]]

snakemake.params[["debug"]] <- snakemake@params[["debug"]]


print(paste("done initiating","snakemake.params"))
print(snakemake.params)


# if(tolower(snakemake@params[["debug"]])=="true"){
if(snakemake@params[["debug"]]){
  rds_file <- snakemake@params[["rds_file"]]  # snakemake@output[["rds_file"]]
  print(paste("saving the snakemake object to",rds_file))
  saveRDS(snakemake,file = rds_file)
  print(paste("done"))
}

stop("Test interruption")  # added to exit script

#- run the main function
res <- myTryCatch(collect.qc(snakemake.params))


lines<-""
if(is.null(res$error)){ # no errors
  qc <- res$value
  # Extract row names
  row_names <- rownames(qc)
  
  # Convert all values to numeric
  qc <- as.data.frame(apply(qc, 2, function(x) as.numeric(as.character(x))))
  
  # Restore row names
  rownames(qc) <- row_names
  
  
  
  qc <- round(qc,dig=3)
  #We could have write a date and the project folder here
  out_file <- snakemake@output[["qc_info_csv"]]
  write.table(qc,out_file,row.names=TRUE,col.names=NA,quote=FALSE,sep=",")
}else{# there was an error
  lines <-paste0(lines,"\nqc.R errors:\n")
  lines <- paste0(lines,res$error$message)
}
if(!is.null(res$warning)){
  lines <-paste0(lines,"\nqc.R warnings:\n")
  lines <- paste0(lines,res$warning$message)
}

if(lines!=""){ # there are either errors or warnings
  file <- snakemake@params[["pipeline_warning_file_path"]]
  # file <- "pipeline_warningGN.txt"
  write(lines,file=file,append=TRUE)
  print(lines)
}

# Close the connection to the log file
sink()

# Optionally, return to normal console output
sink(stdout())
