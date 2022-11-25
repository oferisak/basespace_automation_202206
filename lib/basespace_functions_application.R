remove_biosample_id_from_run_attributes<-function(biosample_ids_to_remove,run_attributes){
  biosamples_in_run_attributes<-run_attributes%>%filter(attribute=='input_list.sample-id')%>%pull(value)%>%str_split(pattern = ',',simplify = F)%>%unlist()
  biosamples_fixed<-setdiff(biosamples_in_run_attributes,biosample_ids_to_remove)
  run_attributes_fixed<-run_attributes
  run_attributes_fixed[run_attributes_fixed$attribute=='input_list.sample-id','value']<-paste0(biosamples_fixed,collapse=',')
  return(run_attributes_fixed)
}

parse_run_sheet<-function(run_sheet_file,attributes_list=list()){
  run_sheet<-readr::read_delim(run_sheet_file,col_names = c('attribute','value') )
  samples<-run_sheet%>%dplyr::filter(attribute=='input_list.sample-id')%>%pull(value)
  # convert biosample names into sample ids
  #validate_biosamples(sample_ids)
  attributes<-run_sheet%>%dplyr::filter(attribute!='input_list.sample-id')
  if (length(samples)>0){
    sample_ids<-get_biosamples_id_by_name(samples,stop_if_not_found = T)
    attributes<-attributes%>%bind_rows(data.frame(attribute='input_list.sample-id',value=paste0(sample_ids,collapse=',')))
  }
  for (missing_attribute in attributes%>%filter(is.na(value))%>%pull(attribute)){
    if (missing_attribute %in% names(attributes_list)){
      attributes[attributes$attribute==missing_attribute,'value']<-attributes_list[[missing_attribute]]
    }else{
      usr_val<-readline(prompt=glue("Enter {missing_attribute}: "))
      attributes[attributes$attribute==missing_attribute,'value']<-usr_val
    }
  }
  verify_project_name(attributes)
  return(attributes)
}

launch_application<-function(run_attributes,app_name,app_version,command_line_args=NA){
  cmd_attributes<-paste0('-o ',run_attributes%>%pull(attribute),":",run_attributes%>%pull(value),collapse=' ')
  run_command<-paste0(glue('bs launch application -n "{app_name}" --app-version {app_version} '),
                      cmd_attributes,
                      collapse=' ')
  if (!is.na(command_line_args)){
    run_command<-paste0(c(run_command,command_line_args),collapse=' ')
  }
  message(glue('Running {app_name} with command: {run_command}'))
  system(run_command)
}

get_application_options<-function(app_id,argument_name=NA){
  bs_command<-as.character(glue('bs launch application --id {app_id} --list -f csv'))
  message(glue('running: {bs_command}'))
  application_options<-system(bs_command,intern=T)
  all_options<-read.table(text=application_options,sep=',',header = T)
  if (!is.na(argument_name)){
    return(all_options%>%filter(Option==argument_name))
  }else{return(all_options)}
}

get_all_applications<-function(){
  applications_text<-system('bs application list -f csv',intern=T)
  return(read.table(text=applications_text,sep=',',header = T))
}

#### DRAGEN Enrichment ####
generate_dragen_enrichment_runsheet <- function(run_sheet_name_prefix,
                                                app_session_name,
                                                project_id,
                                                ht_ref = 'hg38_altmaskedv2-cnv-hla-graph-anchored.v8',
                                                fixed_bed = 'custom',
                                                target_bed_id,
                                                vc_target_bed_padding = 50,
                                                high_coverage = 1,
                                                lowpass = 1,
                                                input_list_sample_ids = c(),
                                                validate_input_sample_ids = T) {
  
  app_session_name<-glue('app-session-name\t{app_session_name}')
  project_id<-glue('project-id\t{project_id}')
  ht_ref<-glue('ht-ref\t{ht_ref}')
  fixed_bed<-glue('fixed-bed\t{fixed_bed}')
  target_bed<-glue('target_bed_id\t{target_bed_id}')
  vc_target_bed_padding<-glue('vc-target-bed-padding\t{vc_target_bed_padding}')
  high_coverage<-glue('high-coverage\t{high_coverage}')
  lowpass<-glue('lowpass\t{lowpass}')
  input_list_sample_ids<-paste0('input_list.sample-id\t',input_list_sample_ids,collapse='\n')
  run_sheet<-paste0(c(app_session_name,
                      project_id,
                      ht_ref,
                      fixed_bed,
                      target_bed,
                      vc_target_bed_padding,
                      high_coverage,
                      lowpass,
                      input_list_sample_ids),
                    collapse='\n')
  write(run_sheet,file = glue('./data/run_sheets/{run_sheet_name_prefix}.dragen_enrichment.csv'))
}

#### DRAGEN Germline ####
generate_dragen_germline_runsheet <- function(run_sheet_name_prefix,
                                                app_session_name,
                                                project_id,
                                                ht_ref = 'hg19_altmaskedv2-cnv-hla-graph-anchored.v8',
                                                cnv_checkbox =1,
                                                cnv_segmentation_mode='slm',
                                                sv_checkbox = 1,
                                                eh_checkbox = 1,
                                                input_list_sample_ids = c(),
                                                validate_input_sample_ids = T) {
  
  app_session_name<-glue('app-session-name\t{app_session_name}')
  project_id<-glue('project-id\t{project_id}')
  ht_ref<-glue('ht-ref\t{ht_ref}')
  cnv_checkbox<-glue('cnv_checkbox\t{cnv_checkbox}')
  cnv_segmentation_mode<-glue('cnv_segmentation_mode\t{cnv_segmentation_mode}')
  sv_checkbox<-glue('sv_checkbox\t{sv_checkbox}')
  eh_checkbox<-glue('eh_checkbox\t{eh_checkbox}')
  input_list_sample_ids<-paste0('input_list.sample-id\t',input_list_sample_ids,collapse='\n')
  run_sheet<-paste0(c(app_session_name,
                      project_id,
                      ht_ref,
                      cnv_checkbox,
                      cnv_segmentation_mode,
                      sv_checkbox,
                      eh_checkbox,
                      input_list_sample_ids),
                    collapse='\n')
  write(run_sheet,file = glue('./data/run_sheets/{run_sheet_name_prefix}.dragen_germline.csv'))
}


#### DRAGEN PopGen ####
generate_dragen_popgen_runsheet <- function(run_sheet_name_prefix,
                                           app_session_name,
                                           project_id,
                                           pipeline_mode = 'datasets',
                                           source_projects,
                                           ht_ref = 'hg38',
                                           output_file_prefix='dragen-popgen') {
  
  app_session_name<-glue('app-session-name\t{app_session_name}')
  project_id<-glue('project-id\t{project_id}')
  ht_ref<-glue('ht-ref\t{ht_ref}')
  source_projects<-glue('source-projects\t{source_projects}')
  pipeline_mode<-glue('pipeline-mode\t{pipeline_mode}')
  output_file_prefix<-glue('output-file-prefix\t{output_file_prefix}')
  basespace_disclaimer<-glue('basespace-labs-disclaimer\tAccepted')
  #gg_remove_nonref<-'"arbitraty:--gg-remove-nonref\ttrue'
  commandline_disclaimer<-'commandline-disclaimer\ttrue'

  run_sheet<-paste0(c(app_session_name,
                      project_id,
                      ht_ref,
                      pipeline_mode,
                      source_projects,
                      #source_datasets,
                      output_file_prefix,
                      #gg_remove_nonref,
                      commandline_disclaimer,
                      basespace_disclaimer),
                    collapse='\n')
  output_file<-glue('./data/run_sheets/{run_sheet_name_prefix}.dragen_popgen.csv')
  message(glue('saving run sheet: {output_file}'))
  write(run_sheet,file = output_file)
}
