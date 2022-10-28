get_all_projects<-function(){
  projects_text<-system('bs projects list -f csv',intern=T)
  return(read.table(text=projects_text,sep=',',header = T))
}

create_project<-function(project_name){
  project_details<-system(glue('bs project create -n {project_name} -f csv'),intern=T)
  return(read.table(text=project_details,sep=',',header = T))
}

verify_project_name<-function(run_attributes){
  project_id<-run_attributes%>%filter(attribute=='project-id')%>%pull(value)
  all_projects<-get_all_projects()
  if (!(project_id%in%all_projects$Id)){stop(glue('Project id {project_id} is not found in the workspace'))}
  project_name<-all_projects%>%filter(Id==project_id)%>%pull(Name)
  message(glue('The run will be saved in project: {project_name} ({project_id})'))
}

# if create_project_if_not_available=T create the project and return the ID
get_project_id_by_name<-function(project_name,create_project_if_not_available=F){
  all_projects<-get_all_projects()
  if (project_name %in% all_projects$Name){
    project_details<-all_projects%>%filter(Name==project_name)
    project_id<-project_details$Id
    message(glue('Project {project_name} exists. Its ID is {project_id}'))
    return(project_id)
  }else{
    if (create_project_if_not_available){
      project_details<-create_project(project_name)
      project_id<-project_details$Id
      message(glue('Project {project_name} did not exist. Created it with ID {project_id}'))
      return(project_id)
    }else{
      stop(glue('There is no project named {project_name} in the workgroup'))
    }
  }
}
