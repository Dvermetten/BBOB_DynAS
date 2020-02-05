# Used regex: cat file_list.txt | grep -oP http[^\[]*tgz > t4

library(stringi)
setwd("~/IOH/BBOB_Data")
dt2 <- fread("t4", header = F)
x <- dt2[[1]]
# substr(x[[1]], start = stri_locate_last(x[[1]], fixed = '/')[[1]] + 1, stop = nchar(x[[1]]) - 4)
for (i in seq(length(x))) {
  # print(x[[i]])
  if (!file.exists(paste0("Data_tgz/", substr(x[[i]],
                                              start = stri_locate_last(x[[i]],
                                                                       fixed = '/')[[1]] + 1,
                                              stop = nchar(x[[i]]))))) {
    download.file(x[[i]], destfile = paste0("Data_tgz/", substr(x[[i]],
                                                                start = stri_locate_last(x[[i]],
                                                                                         fixed = '/')[[1]] + 1,
                                                                stop = nchar(x[[i]]))))
  }
  else
    print(paste0(x[[i]], " Already exists"))
}
#
for (f in list.files("Data_tgz", full.names = T)) {
  untar(f, exdir = "Data")
}
#

failed <- c()
dsl_coco <- DataSetList()
for (d in list.dirs("Data", recursive = F)) {
  tryCatch({
    dsl_temp <- DataSetList(d,
                            maximization = F,
                            format = 'COCO',
                            subsampling = T)
    dsl_coco <- c(dsl_coco, dsl_temp)},
    error = function(e) {
      print(d)
      print(e)
      failed <- c(failed, basename(d))
    }
  )
}


dsl_coco <- DataSetList("Data", maximization = F, format = COCO, subsampling = T)
saveRDS(dsl_coco, file = "dsl_coco.rds")


#Split targets from 1e2 to 1e-8
splits <- rev(seq_FV(c(100,1e-8), length.out = 51, scale = 'log'))

calc_best_split <- function(x){
  #Ignore algorithms which have less then 15 runs
  algs <- get_overview(x)[get_overview(x)$runs >= 15]$algId
  x <- subset(x, algId %in% algs)

  #Get data.table containing the ERTs for all targets
  full_dt <- get_RT_summary(x, ftarget = splits)[,c('algId', 'target', 'ERT')]
  #Reduce to just ERT of the final target
  final_erts <- get_RT_summary(x, 1e-8)[,c('algId','ERT')] %>% set_colnames(c('algId', 'ERT_final'))

  #Get info of best algorithm without split
  best_ERT <- min(final_erts$ERT_final)
  best_algid <- final_erts[which.min(final_erts$ERT)]$algId

  #Data.table manipulations to get ERTs until splitpoint and between final target and splitpoint
  full_dt2 <- merge(full_dt, final_erts, by = 'algId')
  full_dt2 <- transform(full_dt2, ERT_p2 = ERT_final - ERT)
  a3 <- full_dt2[, lapply(.SD, min, na.rm = T), by = c('target'), .SDcols = c('ERT', 'ERT_p2')]
  a4 <- transform(a3, total_ERT = ERT + ERT_p2)

  #Get best split-ERT
  split_row <- a4[, .SD[which.min(total_ERT)], ]

  #Extract corresponding p1 and p2 algorithm names
  dt2_usefull <- full_dt2[which(full_dt2$target == split_row$target)]
  p1_alg <- dt2_usefull[which(dt2_usefull$ERT == split_row$ERT)]$algId
  p2_alg <- dt2_usefull[which(dt2_usefull$ERT_p2 == split_row$ERT_p2)]$algId

  #Only use first algorithm in case of ties
  split_row <- split_row[, `:=`(p1_alg = p1_alg[1], p2_alg = p2_alg[1], best_alg = best_algid[1],
                                best_ert = best_ERT, funcId = get_funcId(x), DIM = get_dim(x))]
  #Print to show progress
  print(split_row)
  return(split_row)
}

#All downloaded data from coco in datasetlist-format (from IOHanalyzer)
# loadRDS(file = "dsl_coco.rds")

funcs <- get_funcId(dsl_coco)
dims <- get_dim(dsl_coco)
dt_all_funcs_dims <- rbindlist(apply(expand.grid(funcs, dims), 1, function(x) {calc_best_split(subset(dsl_coco, funcId == x[[1]], DIM == x[[2]]))}))
dt_all_funcs_dims[, impr := 1 - (total_ERT / best_ert) ]
fwrite(dt_all_funcs_dims, file = "Split_COCO_V3.csv")
