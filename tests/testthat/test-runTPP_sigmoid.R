#@ import future
backports::import(future)
test_that("Normalization and sigmoid works in TPP", {

  data(hdacTR_smallExample, package = "TPP")
  Data<-purrr::map(hdacTR_data,function(x)x %>% dplyr::filter(gene_name %in% c("HDAC1","HDAC2","HDAC3","HDAC4")))
  Data<-list(TPPdata=Data,TPPconfig=hdacTR_config)
  testNorm<-expect_warning(runTPP_sigmoid(Data,NORM=TRUE))


})
