get_biosample_id<-function(biosample_name){
  biosample_id<-system(glue('bs get biosample -n {biosample_name} --terse'),intern = T)
  return(biosample_id)
}

get_all_biosamples<-function(){
  biosamples_text<-system('bs biosample list -f csv',intern=T)
  return(read.table(text=biosamples_text,sep=',',header = T))
}

get_run_biosamples<-function(run_id=NA,run_experiment_name=NA){
  all_runs<-get_all_runs()
  if (!is.na(run_experiment_name)){
    run_id<-all_runs%>%filter(ExperimentName==run_experiment_name)%>%pull(Id)
  }
  all_biosamples<-get_all_biosamples()
  run_biosamples_ids<-system(glue('bs run property get --id {run_id} --property-name Input.BioSamples --field Id -f csv'),intern=T)
  run_biosamples<-all_biosamples%>%filter(Id %in% run_biosamples_ids)
  return(run_biosamples)
}

get_biosample_fastqs<-function(biosample_id){
  message(glue('retrieving fastq file names for biosample: {biosample_id}'))
  biosample_content<-system(glue('bs biosample content --id {biosample_id} -f csv'),intern=T)
  biosample_content<-read.table(text=biosample_content,sep=',',header = T)
  biosample_fastqs<-biosample_content%>%filter(grepl('fastq.gz',FilePath))
  if (nrow(biosample_fastqs)==0){message(glue('Biosample {biosample_id} has no fastqs (perhaps it was transfered to a different workgroup)'))}
  return(biosample_fastqs)
}

get_biosamples_id_by_name<-function(biosample_names,stop_if_not_found=T){
  message(glue('getting the biosample ids for the given biosmaples ({paste0(biosample_names,collapse=", ")})..'))
  all_biosamples<-get_all_biosamples()
  # check that all names are in the list
  biosamples_not_found<-biosample_names[!(biosample_names %in% (all_biosamples%>%pull(BioSampleName)))]
  if (stop_if_not_found){
    if (length(biosamples_not_found)>0){stop(glue('could not find biosamples: {paste0(biosamples_not_found,collapse=",")}'))}
  }else{
    message(glue('could not find biosamples: {paste0(biosamples_not_found,collapse=",")}'))
  }
  biosamples_ids<-all_biosamples%>%filter(BioSampleName%in%biosample_names)%>%pull(Id)
  message(glue('found: {paste0(biosamples_ids,collapse=", ")}'))
  return(biosamples_ids)
}

get_biosamples_fastqs_from_attributes<-function(run_attributes){
  biosamples<-run_attributes%>%filter(attribute=='input_list.sample-id')%>%pull(value)%>%str_split(pattern = ',',simplify = F)%>%unlist()
  bs_fqs_df<-list()
  for (bs in biosamples){
    bs_fqs<-get_biosample_fastqs(bs)
    bs_fqs_df[[bs]]<-bs_fqs
  }
  return(bs_fqs_df)
}

validate_biosample_ids<-function(){
  bs_fqs_list<-get_biosamples_fastqs_from_attributes(run_attributes)
  # get the missing biosamples names
  all_biosamples<-get_all_biosamples()
  biosamples_without_fqs<-data.frame(Id=as.integer(names(bs_fqs_list)[which(lapply(bs_fqs_list,nrow)==0)]))
  biosamples_without_fqs<-biosamples_without_fqs%>%left_join(all_biosamples)
}

get_biosample_attributes<-function(biosample_id){
  attributes_text<-system(glue('bs biosample get --id {biosample_id} -f csv'),intern=T)
  return(read.table(text=attributes_text,sep=',',header = T))
}

get_project_ids_from_biosamples_sheet<-function(biosamples_sheet){
  biosamples_df<-readr::read_delim(biosamples_sheet,delim = '\t')
  projects<-biosamples_df%>%pull(project_id)%>%unique()
}
