#' Read in protein data from proteome discoverer
#'
#' Read in Excel file and apply minimal pre-processing
#'
#'
#' @param  f.  Path of Excel file output from Proteome Discoverer
#' @return a dataframe containing extracted information
#'
#' @importFrom readxl read_excel
#' @import dplyr
#' @importFrom tidyr gather
#' @importFrom stringr str_extract
#' @export
#'
#'
read_cetsa_concentrations <- function(protein_path,peptide_path,
                       Frac=TRUE,
                       concentrations=set_concentration(16,concentrations=unique(CETSA_OnePot_annotation$Condition),replicates=3),
                       solvent="DMSO",
                       CARRIER=TRUE){
  file.list<-protein_path
  i=1
  peptide_path<-as.character(peptide_path)
  protein_path<-as.character(protein_path)

  find<-c('[:digit:][:digit:][:digit:][N|C]|126|131|131N')
  df<-list()
  df1<-list()
  df3<-data.frame()
  df2<-list()
  df.raw<-list(data.frame())

  #read_PSMs and proteins

  Proteins<-readProteins(protein_path)
  PSMs<-readPSMs(peptide_path)

  if(isTRUE(Frac)){
    Proteins <- Proteins |>
      tidylog::pivot_longer(cols=colnames(Proteins)[stringr::str_detect(colnames(Proteins),"[:digit:][:digit:][:digit:][N|C]|[:digit:][:digit:][:digit:]")],
                            names_to = "id",
                            values_to ="value") |>
      dplyr::mutate(File.ID = stringr::str_extract(id,"[:upper:][[:digit:]]+"),
                    Channel = stringr::str_extract(id,'[:digit:][:digit:][:digit:][N|C]|[:digit:][:digit:][:digit:]'),
                    Abundance = as.numeric(value),
                    sample_id = stringr::str_extract(File.ID,"[:upper:][[:digit:]]+"))
    if(any(stringr::str_detect(stringr::str_to_lower(Proteins$id),"control"))){
      Proteins$treatment<-ifelse(stringr::str_detect(stringr::str_to_lower(Proteins$id),"control"),"vehicle","treated")
    }

    Mapping_PSMs<-PSMs|>
      tidylog::pivot_longer(cols=colnames(PSMs)[stringr::str_detect(colnames(PSMs),"[:digit:][:digit:][:digit:][N|C]|[:digit:][:digit:][:digit:]")],
                            names_to = "id",
                            values_to ="value") |>
      dplyr::select(id,File.ID,Spectrum.File)|>dplyr::distinct()|>
      dplyr::mutate(Fraction=stringr::str_remove(File.ID,"[:upper:][[:digit:]]+."),
                    File.ID=stringr::str_extract(File.ID,"[:upper:][[:digit:]]+"),
                    Channel = stringr::str_extract(id,'[:digit:][:digit:][:digit:][N|C]|[:digit:][:digit:][:digit:]'))|>
      dplyr::select(File.ID,Spectrum.File,Channel)|>
      dplyr::distinct()

    Proteins<-Proteins |>
      dplyr::filter(!is.na(File.ID))|>
      dplyr::inner_join(Mapping_PSMs,relationship="many-to-many")

  }else{
    Mapping_PSMs<-PSMs|>
      tidylog::pivot_longer(cols=colnames(PSMs)[stringr::str_detect(colnames(PSMs),"[:digit:][:digit:][:digit:][N|C]|[:digit:][:digit:][:digit:]")],
                            names_to = "id",
                            values_to ="value") |>
      dplyr::select(id,File.ID,Spectrum.File)|>dplyr::distinct()|>
      dplyr::mutate(File.ID=stringr::str_extract(File.ID,"[:upper:][[:digit:]]+"),
                    Channel = stringr::str_extract(id,'[:digit:][:digit:][:digit:][N|C]|[:digit:][:digit:][:digit:]'))|>
      dplyr::select(File.ID,Spectrum.File,Channel)|>
      dplyr::distinct()
    Proteins <- Proteins %>%
      tidylog::pivot_longer(cols=colnames(Proteins)[stringr::str_detect(colnames(Proteins),"[:digit:][:digit:][:digit:][N|C]|[:digit:][:digit:][:digit:]")],
                            names_to = "id",
                            values_to ="value") %>%
      dplyr::mutate(File.ID = stringr::str_extract(id,"[:upper:][[:digit:]]+"),
                    Channel = stringr::str_extract(id,'[:digit:][:digit:][:digit:][N|C]|[:digit:][:digit:][:digit:]'),
                    Abundance = as.numeric(value),
                    Fraction = stringr::str_remove(File.ID,"[:upper:][[:digit:]]+."),
                    sample_id = stringr::str_extract(File.ID,"[:upper:][[:digit:]]+"))
    Proteins<-Proteins|> dplyr::filter(!is.na(File.ID)) |>
      dplyr::inner_join(Mapping_PSMs,relationship="many-to-many")
  }

  #if this is a protein file
  if(any(stringr::str_detect(names(PSMs),"Accessions"))){
    master<-names(PSMs)[stringr::str_detect(names(PSMs),"Master")]
  }else{
    master<-names(PSMs)[stringr::str_detect(names(PSMs),"Accessions")]
  }

  #Select protein columns
  if(any(names(Proteins)=="Biological.Process")&any(names(Proteins[[1]]=="MW.kDa."))){
    Proteins<-dplyr::bind_rows(Proteins) |>
      dplyr::rename("Cell_Component"="Cellular.Component","MW_kDa"="MW.kDa.","Bio_Process"="Biological.Process","Molecular_Function"="Molecular.Function")
  }
  #Add additional information from the experiment

  Proteins$CC <-ifelse(stringr::str_detect(Proteins$id,solvent),0,1)
  Proteins$Run<-Proteins$Spectrum.File
  Proteins<-dplyr::bind_rows(Proteins) |> dplyr::distinct()
  return(Proteins)

}

