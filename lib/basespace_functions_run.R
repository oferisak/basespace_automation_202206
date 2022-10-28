get_all_runs<-function(){
  runs_text<-system('bs run list -f csv',intern=T)
  return(read.table(text=runs_text,sep=',',header = T))
}

get_run_id_by_name<-function(experiment_name){
  all_runs<-get_all_runs()
  run_id<-all_runs%>%filter(ExperimentName==experiment_name)%>%pull(Id)
  if (length(run_id)==0){stop(glue('There is no run id associated with experiment (run name) {experiment_name}'))}
  return(run_id)
}

get_run_statistics<-function(run_id){
  run_stats<-system(glue('bs run get --id {run_id} -f csv'),intern=T)
  return(read.table(text=run_stats,sep=',',header = T))
}
