
readPSMs<-function(wd){
  wd<-wd
  h<-list.files(wd,pattern="PSMs.txt")|>as.list()
  if(length(h)>0){
    if(any(stringr::str_detect(h[[1]],"txt"))){
      h<-list.files(wd,pattern="PSMs.txt")|>as.list()
      setwd(wd)
      PSMs<-purrr::map(h,function(x) read.delim(x))|>dplyr::bind_rows()
    }
  }else{
    h<-list.files(wd,pattern="PSMs.xlsx")|>as.list()
    PSMs<-purrr::map(h,function(x) readxl::read_xlsx(path=paste0(wd,"/",x),.name_repair = "universal")) |> dplyr::bind_rows()
  }
  return(PSMs)
}
