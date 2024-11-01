readProteins<-function(wd){
  wd<-wd
  h<-list.files(wd,pattern="Proteins.txt")|>as.list()
  if(length(h)>0){
    if(any(stringr::str_detect(h[[1]],".txt"))){
      h<-list.files(wd,pattern="Proteins.txt")|>as.list()

      Proteins<-purrr::map(h,function(x) read.delim(paste0(wd,"/",x))) |> dplyr::bind_rows()
      Proteins<-Proteins[,!stringr::str_detect(names(Proteins),"Found")]
    }
  }else{
    h<-list.files(wd,pattern="Proteins.xlsx")|>as.list()
    Proteins<-purrr::map(h,function(x) readxl::read_xlsx(path=paste0(wd,"/",x),.name_repair = "universal")) |> dplyr::bind_rows()
    Proteins<-Proteins[,!stringr::str_detect(names(Proteins),"Found")]
  }

  return(Proteins)
}
