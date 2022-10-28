get_whoami<-function(){
  whoami_txt<-system('bs whoami -f csv',intern = T)
  return(read.table(text=whoami_txt,sep=',',header = T))
}

verify_workgroup<-function(workgroup){
  whoami<-get_whoami()
  cur_workgroup<-whoami%>%pull(Name)
  if (workgroup!=cur_workgroup){stop(glue('workgroup is not {workgroup} | current workgroup is {cur_workgroup}'))}
  else{message(glue('Current workgroup is {cur_workgroup}'))}
}

get_all_datasets<-function(){
  datasets_text<-system('bs datasets list -f csv',intern=T)
  return(read.table(text=datasets_text,sep=',',header = T))
}
