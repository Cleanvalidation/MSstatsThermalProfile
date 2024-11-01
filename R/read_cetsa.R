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
#' @importFrom stringr stringr::str_extract
#' @export
#'
read_cetsa <- function(protein_path,peptide_path,Frac=TRUE,temps=set_temps(11,c(37,41, 44,47,50,53,56,59,63,67,68)),solvent="DMSO",CARRIER=TRUE){
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
  #shorten spectrum.file column
  if(any(stringr::str_detect(names(PSMs),"Spectrum"))){
    PSMs$Spectrum.File<-stringr::str_remove_all(PSMs$Spectrum.File,"[[:digit:]]+.raw")
    PSMs$Spectrum.File<-stringr::str_remove(PSMs$Spectrum.File,"Fx[[:digit:]]+")
    PSMs$Spectrum.File<-stringr::str_remove(PSMs$Spectrum.File,"Fx[:digit:]")

    PSMs$File.ID<-stringr::str_extract(PSMs$File.ID,"F[[:digit:]]")
    if(any(names(PSMs)=="Master.Protein.Accessions")){
      PSMs<-PSMs|>dplyr::select(Master.Protein.Accessions,File.ID,Spectrum.File)|>dplyr::distinct()
    }else if(any(names(PSMs)=="Protein.Accessions")){
      PSMs<-PSMs|>dplyr::select(Protein.Accessions,File.ID,Spectrum.File)|>dplyr::distinct()
    }
  }
  if(isTRUE(Frac)){

    Proteins<-Proteins |>
      dplyr::select(Accession, dplyr::starts_with("Abundance")) |>
      tidylog::pivot_longer(cols = colnames(Proteins)[stringr::str_detect(colnames(Proteins), "[:digit:][:digit:][:digit:][N|C]|[:digit:][:digit:][:digit:]")],
                            names_to = "id",
                            values_to = "value") |>
      dplyr::mutate(File.ID = stringr::str_extract(id, "[:upper:][[:digit:]]+"),
                    Channel = stringr::str_extract(id, '[:digit:][:digit:][:digit:][N|C]|[:digit:][:digit:][:digit:]'),
                    Abundance = as.numeric(value),
                    sample_id = stringr::str_extract(File.ID, "[:upper:][[:digit:]]+"))|>
      dplyr::distinct()

    checkProteins<-Proteins|>dplyr::group_by(Accession)|>dplyr::mutate(missing_channels=sum(is.na(Abundance)))|>dplyr::ungroup()
    if(any(checkProteins$missing_channels>=44)){
      Proteins<-Proteins|>
        dplyr::arrange(Abundance)|>
        dplyr::group_by(Accession,File.ID)|>
        dplyr::group_split()
      Proteins<-lapply(Proteins,function(x) x|>
                         head(11)|>
        dplyr::mutate(missing_channels=sum(is.na(Abundance)))|>
        dplyr::ungroup())|>dplyr::bind_rows()
    }
    if(any(stringr::str_detect(stringr::str_to_lower(Proteins$id),"control"))){
      Proteins$treatment<-ifelse(stringr::str_detect(stringr::str_to_lower(Proteins$id),"control"),"vehicle","treated")
    }

    Mapping_PSMs<-PSMs|>
      dplyr::select(Spectrum.File,File.ID)|>
      dplyr::distinct()|>
      dplyr::mutate(File.ID=stringr::str_extract(File.ID,"[:upper:][[:digit:]]+"),
                                                 treatment=ifelse(
                                                   stringr::str_detect(Spectrum.File,solvent),
                                                                  "vehicle","treated"))|>
      dplyr::group_by(treatment)|>
      dplyr::mutate(TechRepMixture=seq(1,dplyr::n()),
                    Experiment = paste0(treatment,"_",TechRepMixture))|>
      dplyr::distinct()|>
      dplyr::ungroup()
    Proteins<-Proteins |>
      dplyr::filter(!is.na(File.ID))|>
      dplyr::inner_join(temps, by ="Channel")|>
      dplyr::inner_join(Mapping_PSMs)

  }else{
    Proteins <- Proteins |>
      dplyr::select(Accession,dplyr::starts_with("Abundance"))|>
      dplyr::distinct()|>
      tidylog::pivot_longer(cols=colnames(Proteins)[stringr::str_detect(colnames(Proteins),"[:digit:][:digit:][:digit:][N|C]|[:digit:][:digit:][:digit:]")],
                            names_to = "id",
                            values_to ="value")|>
      dplyr::mutate(File.ID = stringr::str_extract(id,"[:upper:][[:digit:]]+"),
                    Channel = stringr::str_extract(id,'[:digit:][:digit:][:digit:][N|C]|[:digit:][:digit:][:digit:]'),
                    Abundance = as.numeric(value),
                    sample_id = stringr::str_extract(File.ID,"[:upper:][[:digit:]]+"))|>
      dplyr::distinct()
    #check proteins for missing values
      checkProteins<-Proteins|>dplyr::group_by(Accession)|>dplyr::mutate(missing_channels=sum(is.na(Abundance)))|>dplyr::ungroup()
      if(any(checkProteins$missing_channels>=44)){
        Proteins<-Proteins|>
         # dplyr::filter(!is.na(Abundance))|>
          dplyr::group_by(Accession)|>
          dplyr::mutate(missing_channels=sum(is.na(Abundance)))|>
          dplyr::ungroup()
      }
      if(any(stringr::str_detect(stringr::str_to_lower(Proteins$id),"control"))){
        Proteins$treatment<-ifelse(stringr::str_detect(stringr::str_to_lower(Proteins$id),"control"),"vehicle","treated")

        Proteins<-Proteins |>
          dplyr::select(-id) |>
          dplyr::distinct()
      }


    Mapping_PSMs<-PSMs|>
      dplyr::select(Spectrum.File,File.ID)|>
      dplyr::distinct()|>
      dplyr::mutate(File.ID=stringr::str_extract(File.ID,"[:upper:][[:digit:]]+"),
                    treatment=ifelse(stringr::str_detect(Spectrum.File,solvent),"vehicle","treated"))|>
      dplyr::group_by(treatment)|>
      dplyr::mutate(TechRepMixture=seq(1,dplyr::n()),
                    Experiment = paste0(treatment,"_",TechRepMixture))
    Proteins<-Proteins |>
      dplyr::distinct()|>
      dplyr::filter(!is.na(File.ID))|>
      dplyr::inner_join(temps, by ="Channel")|>
      dplyr::inner_join(Mapping_PSMs)|>
      dplyr::select(Accession,File.ID,Channel,Abundance,sample_id,treatment,temperature,Spectrum.File)|>
      dplyr::distinct()

  }

  #if this is a protein file
  if(any(stringr::str_detect(names(PSMs),"Accessions"))){
    master<-names(PSMs)[stringr::str_detect(names(PSMs),"Master")]
  }else{
    master<-names(PSMs)[stringr::str_detect(names(PSMs),"Accessions")]
  }

  #Proteins<-Proteins[!is.na(Proteins$File.ID),]|>dplyr::group_by(Accession,File.ID)|>dplyr::filter(!is.na(value))|>dplyr::ungroup()

  #Select protein columns
  if(any(names(Proteins)=="Biological.Process")&any(names(Proteins[[1]]=="MW.kDa."))){
    Proteins<-dplyr::bind_rows(Proteins) |>
      dplyr::rename("Cell_Component"="Cellular.Component","MW_kDa"="MW.kDa.","Bio_Process"="Biological.Process","Molecular_Function"="Molecular.Function")
  }
  #Add additional information from the experiment

  Proteins$CC <-ifelse(stringr::str_detect(Proteins$Spectrum.File,solvent),0,1)
  if(any(stringr::str_detect(Proteins$Spectrum.File,solvent))){
    Proteins$treatment<-as.factor(ifelse(stringr::str_detect(Proteins$Spectrum.File,solvent),"vehicle","treated"))
    Proteins$Condition<-paste0(Proteins$Channel,"_",Proteins$treatment)
  }
  Proteins$Subject<-Proteins$File.ID
  if(any(names(Proteins)=="BioReplicate")){
    Proteins$Mixture<-Proteins$BioReplicate
  }


  Proteins<-Proteins|>
    dplyr::left_join(temps)|>
    dplyr::group_by(Accession,File.ID)|>
    #check missing channels in plex
    dplyr::mutate(missing_channels=sum(is.na(Abundance)))|>
    dplyr::ungroup()|>
    dplyr::distinct()

  return(Proteins)

}

