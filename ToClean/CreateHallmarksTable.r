H_featuers


na.pad <- function(x,len){
  x[1:len]
}

makePaddedDataFrame <- function(l,...){
  maxlen <- max(sapply(l,length))
  data.frame(lapply(l,na.pad,len=maxlen),...)
}

hallmark_names


tbl <- makePaddedDataFrame(list("Sustaining Proliferative Signaling"=H_featuers[[1]],
                         "Evading Growth Suppressors"=H_featuers[[2]],
                         "Avoiding Immune Destruction" =H_featuers[[3]],
                         "Enabling Replicative Immortality"=H_featuers[[4]],
                         "Tumour-Promoting Inflammation"=H_featuers[[5]],
                         "Activating Invasion and Metastasis"=H_featuers[[6]],
                         "Inducing Angiogenesis"=H_featuers[[7]],
                         "Genome Instability and Mutation"=H_featuers[[8]],
                         "Resisting Cell Death"=H_featuers[[9]],
                         "Deregulating Cellular Energetics"=H_featuers[[10]],
                         "Senescent cells"=H_featuers[[11]],
                         "Nonmutational Epigenetic reprogramming"=H_featuers[[12]],
                         "Unlocking Phenotypic Plasticity"=H_featuers[[13]]))

View(tbl)
write.table(tbl, "Desktop/IJC/datasets/IGTP/figuresPaper/tables/Supp_hallmarks.xls", sep = "\t", row.names = F)
