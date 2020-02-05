dt <- fread("Split_COCO_V3.csv")

#Basic algorithm number statistics
bests <- unique(dt$best_alg)
best_p1s <- unique(dt$p1_alg)
best_p2s <- unique(dt$p2_alg)
print(paste0("length bests: ", length(bests), "; length best_p1s: ", length(best_p1s), "; length best_p2s: ", length(best_p2s)))
print(paste0("length total used: ", length(unique(c(bests, best_p1s, best_p2s)))))
print(paste0("length total available: ", length(unique(get_algId(dsl_coco)))))
o <- get_FV_overview(dsl_coco, target = 1)
algs <- o[o$runs >= 15]$algId %>% unique
print(paste0("length total usable algs (>= 15 runs): ", length(algs)))

#Basic info about improvement
print(paste0("Mean improvement: ", format(mean(dt$impr), digits = 3), " (median: ", format(median(dt$impr), digits = 3), ")"))
dt[]
#Some initial plots (uses the function from the file test_aggregated_relative_ert)
p <- plot_general_data(dt, x_attr = 'funcname', y_attr = "impr", legend_attr = "funcname")
save_plotly(p, "Impr_distr_per_func.pdf", 1000, 600)

x <- subset(dsl_coco, DIM == 40 && funcId == 24)
algs <- get_overview(x)[get_overview(x)$runs >= 15]$algId
x2 <- subset(x, algId %in% algs)
rts <- c()
for (alg in get_algId(x2)) {
  x3 <- subset(x2, algId == alg)
  rts_temp <- get_runtimes(x3)
  rts <- c(rts, max(rts_temp))
}
rts

### Added 17/12/19

splits <- rev(seq_FV(c(100,1e-8), length.out = 51, scale = 'log'))
library(magrittr)

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

calc_split_v2 <- function(x1, x2){
  dt_a1 <- get_RT_summary(x1, ftarget = splits)[,c('target', 'ERT')]
  dt_a2 <- get_RT_summary(x2, ftarget = splits)[,c('target', 'ERT')]
  final_ert = dt_a2[51,'ERT']
  dt_a2 = dt_a2[,ERT2 := final_ert$ERT - ERT][,c('target', 'ERT2')]
  dt = merge(dt_a1, dt_a2, by = 'target')[, total_ert := ERT + ERT2]
  # print(dt)
  res <- c(which.min(dt$total_ert), min(dt$total_ert))
  # print(res)
  return(res)
}

#All downloaded data from coco in datasetlist-format (from IOHanalyzer)
load(file = "dsl_coco.rda")

funcs <- get_funcId(dsl_coco)
dims <- get_dim(dsl_coco)
dt_all_funcs_dims <- rbindlist(apply(expand.grid(funcs, dims), 1, function(x) {calc_best_split(subset(dsl_coco, funcId == x[[1]], DIM == x[[2]]))}))
dt_all_funcs_dims[, impr := 1 - (total_ERT / best_ert) ]
fwrite(dt_all_funcs_dims, "../BBOB_Data/Split_COCO_V3.csv")


for (fid in get_funcId(dsl_coco)) {
  for (dim in get_dim(dsl_coco)) {
    dsl_sub <- subset(dsl_coco, funcId == fid && DIM == dim)
    algs_ <- get_overview(dsl_sub)[get_overview(dsl_sub)$runs >= 15]$algId

    temp <- get_RT_summary(dsl_sub, ftarget = 1e-8)[,c('algId', 'ERT')]
    algs <- temp[!is.infinite(temp$ERT)]$algId

    erts <- array(NaN, dim = c(length(algs), length(algs)))
    splitpoints <- array(-1, dim = c(length(algs), length(algs)))

    best_ert <- min(temp$ERT, na.rm = T)

    for (i in seq_along(algs)) {
      for (j in seq_along(algs)) {
        tryCatch({
          temp = calc_split_v2(subset(dsl_sub, algId == algs[[i]]), subset(dsl_sub, algId == algs[[j]]))
          # temp[,  impr := 1 - (total_ERT / best_ert) ]
          if (!is.na(temp[[1]]))
            erts[i,j] = best_ert/temp[[2]]
          splitpoints[i,j] = temp[[1]]
        }, error = function(e) {})
      }
      print(i)
    }
  fwrite(erts, file = paste0("ERTS/ERTS_F", fid, "_D", dim), row.names = F, col.names = F)
  fwrite(splitpoints, file = paste0("ERTS/Splits_F", fid, "_D", dim), row.names = F, col.names = F)
  }
}

fid <- 10
dim <- 10

dsl_sub <- subset(dsl_coco, funcId == fid && DIM == dim)
algs_ <- get_overview(dsl_sub)[get_overview(dsl_sub)$runs >= 15]$algId

temp <- get_RT_summary(dsl_sub, ftarget = 1e-8)[,c('algId', 'ERT')]
algs <- temp[!is.infinite(temp$ERT)]$algId

erts <- array(NaN, dim = c(length(algs), length(algs)))
splitpoints <- array(-1, dim = c(length(algs), length(algs)))

best_ert <- min(temp$ERT, na.rm = T)

for (i in seq_along(algs)) {
  for (j in seq_along(algs)) {
    tryCatch({
    temp = calc_split_v2(subset(dsl_sub, algId == algs[[i]]), subset(dsl_sub, algId == algs[[j]]))
    # temp[,  impr := 1 - (total_ERT / best_ert) ]
    if (!is.na(temp[[1]]))
      erts[i,j] = best_ert/temp[[2]]
    splitpoints[i,j] = temp[[1]]
    }, error = function(e) {})
  }
  print(i)
}


temp = calc_best_split(c(subset(x1, funcId == fid && DIM == dim) ,subset(x2, funcId == fid && DIM == dim)))[,c('total_ERT', 'best_ert')]
temp[,  impr := 1 - (total_ERT / best_ert) ]

get_all_split_imprs <- function(x) { #x is DataSetList of 1 func, 1 dim, filetered alg to hit target
  algnames <- get_algId(x)
  dt <- get_RT_summary(x, splits)[,c('algId','target', 'ERT')]
  dt <- dt %>% group_by(algId) %>% mutate( ERT2 = max(ERT) - ERT )
  dt <- dt %>% group_by(target) %>% merge(dt, by = 'target') %>% as.data.table
  dt <- dt[, total_ERT := `ERT.x` + `ERT2.y`] %>% group_by(algId.x, algId.y) %>% slice(which.min(total_ERT)) %>% as.data.table
  return(dt[, c('target', 'algId.x', 'algId.y', 'total_ERT')])
}

for (fid in get_funcId(dsl_coco)) {
  for (dim in get_dim(dsl_coco)) {
    dsl_sub <- subset(dsl_coco, funcId == fid && DIM == dim)
    algs_ <- get_overview(dsl_sub)[get_overview(dsl_sub)$runs >= 15]$algId

    temp <- get_RT_summary(subset(dsl_sub, algId %in% algs_), ftarget = 1e-8)[,c('algId', 'ERT')]
    algs <- temp[!is.infinite(temp$ERT)]$algId

    dt <- get_all_split_imprs(subset(dsl_sub, algId %in% algs))
    fwrite(dt, file = paste0("ERTS/DataTable_F", fid, "_D", dim), row.names = F, col.names = T)
  }
}

cutoff_val <- -0.01




dt <- fread("ERTS/DataTable_F1_D20")
erts <- dt %>% acast(algId.x ~ algId.y, value.var = 'total_ERT')
rel_imprs <- 1 - erts/min(diag(erts))

rel_imprs[rel_imprs < cutoff_val] <- cutoff_val
plot_ly(x = dt$algId.x, y = dt$algId.y, z = rel_imprs, type = 'heatmap')
rowSums(rel_imprs) %>% which.max %>% names #Best A1
colSums(rel_imprs) %>% which.max %>% names #Best A2


get_aggr_imprs <- function(fid){
  algs <- c()
  for (dim in get_dim(dsl_coco)) {
    dt <- fread(paste0("ERTS/DataTable_F", fid, "_D", dim))

    algs <- c(algs, unique(dt$algId.x))
  }

  algs <- Filter(function(elem) length(which(algs == elem)) == 6, algs) %>% unique

  temp_impr <- array(0, dim = c(length(algs), length(algs)))
  for (dim in get_dim(dsl_coco)) {
    dt <- fread(paste0("ERTS/DataTable_F", fid, "_D", dim))
    dt <- dt[dt[, algId.x %in% algs]]
    dt <- dt[dt[, algId.y %in% algs]]

    erts <- dt %>% acast(algId.x ~ algId.y, value.var = 'total_ERT')
    rel_imprs <- 1 - erts/min(diag(erts))

    rel_imprs[rel_imprs < cutoff_val] <- cutoff_val
    temp_impr <- temp_impr + rel_imprs
  }
  temp_impr <- temp_impr / length(get_dim(dsl_coco))
  rowSums(temp_impr) %>% which.max %>% names %>% print #Best A1?
  colSums(temp_impr) %>% which.max %>% names %>% print #Best A2?
  diag(temp_impr) %>% which.max %>% names %>% print #Best static?
  p <- plot_ly(x = algs, y = algs, z = temp_impr, type = 'heatmap')
  save_plotly(p, paste0("Heatmaps/Aggregated_F", fid, "_v2.pdf"), 1000, 1000)
}


# fid <- 1
# dim <- 5
# dt <- fread(paste0("ERTS/DataTable_F", fid, "_D", dim))
# algs <- unique(dt$algId.x)
# alg <- algs[[1]]
# means <- lapply(algs, function(alg) {
# dt_temp <- dt[algId.x == alg]
# rels <- lapply(algs, function(alg2) {
#   curr_ert <- dt_temp[algId.y == alg2]$total_ERT
#   rel <- curr_ert / rmean[[alg2]]
#   rel
# }) %>% unlist %>% mean
# })
# erts <- dt %>% acast(algId.x ~ algId.y, value.var = 'total_ERT')
# cmean <- colMeans(erts)
#
# dt[, cmean := cmean[algId.y]]
# dt[, impr := 1 - total_ERT / cmean]
# dt[, ]

### Value as A1 -> interpretation 1
get_impr_A1_mean <- function(fid, dim, which = 'a1') {
  dt <- fread(paste0("ERTS/DataTable_F", fid, "_D", dim))
  erts <- dt %>% acast(algId.x ~ algId.y, value.var = 'total_ERT')
  # algs <- unique(dt$algId.x)
  # p <- plot_ly(x = algs, y = algs, z = erts, type = 'heatmap')

  if (which == 'a1') {
    cmean <- colMeans(erts)
    dt[, cmean := cmean[algId.y]]
  }
  else {
    cmean <- rowMeans(erts)
    dt[, cmean := cmean[algId.x]]
  }
  dt[, impr := 1 - total_ERT / cmean]
  if (which == 'a1')
    dt <- dt[, .(impr = mean(impr)), by = algId.x]
  else
    dt <- dt[, .(impr = mean(impr)), by = algId.y]
  dt[, funcId := fid]
  dt[, dim := dim]
  dt[order(-impr)]
}

### Value as A1 -> interpretation 2
get_impr_A1_rel_diag <- function(fid, dim, which = 'a1') {
  dt <- fread(paste0("ERTS/DataTable_F", fid, "_D", dim))
  erts <- dt %>% acast(algId.x ~ algId.y, value.var = 'total_ERT')
  if (which == 'a2') {
    T1 <- colMins(erts, T)
    names(T1) <- colnames(erts)
  }
  else {
    T1 <- rowMins(erts, T)
    names(T1) <- rownames(erts)
  }
  T1 <- 1 - (T1 / min(diag(erts)))
  dt <- data.table(T1, names(T1), fid, dim)
  colnames(dt) <- c("impr", "algId", "funcId", "DIM")
  dt
}

### Value as A1 -> interpretation 3
get_impr_A1 <- function(fid, dim, which = 'a1') {
  dt <- fread(paste0("ERTS/DataTable_F", fid, "_D", dim))
  erts <- dt %>% acast(algId.x ~ algId.y, value.var = 'total_ERT')
  if (which == 'a2') {
    T1 <- colMins(erts, T)
    names(T1) <- colnames(erts)
  }
  else {
    T1 <- rowMins(erts, T)
    names(T1) <- rownames(erts)
  }
  T1 <- T1 / min(erts)
  dt <- data.table(T1, names(T1), fid, dim)
  colnames(dt) <- c("impr", "algId", "funcId", "DIM")
  dt
}

plot_heatmap_fid_dim <- function(fid, dim) {
  dt <- fread(paste0("ERTS/DataTable_F", fid, "_D", dim))
  erts <- dt %>% acast(algId.x ~ algId.y, value.var = 'total_ERT')
  algs <- unique(dt$algId.x)
  p <- plot_ly(x = algs, y = algs, z = log(erts), type = 'heatmap')
  p
}

get_impr_heatmap_per_dim <- function(dim, which = 'a1', subgroup = T, smaller = T) {
  dt2 <- lapply(seq_len(24), function(x) get_impr_A1(x, dim, which)) %>% rbindlist
  if (subgroup) {
    dt3 <- dt2[, .SD[which.min(impr)], by = funcId]
    algs <- unique(dt3$algId)
    # print(length(unique(algs)))

    arr <- dt2 %>% acast(algId ~ funcId, value.var = 'impr')
    arr[arr > 3] <- 3
    arr[is.na(arr)] <- 3 #arbitrary
    dt3 <- sort(rowSums(arr), F)
    idx <- 1
    while (length(algs) < 25) {
      algs <- unique(c(algs, names(dt3)[[idx]]))
      # print(names(dt3)[[idx]])
      idx <- idx + 1
      if (idx > length(dt3)) {
        break
      }
    }
    print(algs)
    if (smaller) {
      algs <- names(dt3)[1:15]
    }

    dt2 <- dt2[algId %in% algs]
    # print(dt2)
  }
  arr <- dt2 %>% acast(algId ~ funcId, value.var = 'impr')
  arr[arr < -1] <- -1
  # colscales <- ifelse(which == "a1", 'RdBu', 'Viridis')
  colscales <- 'Viridis'
  print(min(arr, na.rm = T))
  p <- plot_ly(x = colnames(arr), y = rownames(arr), z = arr, type = 'heatmap',
               colorscale = colscales, zmin = 1, zmax = 3, zauto = F)
  p %<>% layout(xaxis = list(type = 'category', title = "Function ID"), yaxis = list(tickfont = list(size = 10)))
  save_plotly(p, paste0("Heatmaps/Heatmap_Mean_Impr_dim", dim, "_as", which, ifelse(subgroup, "_subgroup_", "_"), "v4.pdf"), 750, 500)
  # arr
  # p
}

### Verification of dsl_coco (14/01/2020)
dsl_coco <- DataSetList()
for (f in list.files("~/repository/bbob/", full.names = T)) {
  temp <- readRDS(f)
  dsl_coco <- c(dsl_coco, temp)
}

### Figure generation for ERT distribution of static algorithms (14/01/2020)
dt <- get_RT_summary(dsl_coco, 1e-8)
dt <- dt[target == 1e-8]
dt <- dt[, funcId := paste0('F', sprintf("%02d", funcId))]

for (dim in get_dim(dsl_coco)) {
  # TODO: Turn into boxplot
  p <- IOH_plot_ly_default(x.title = "Function ID", y.title = "ERT")
  p %<>% add_boxplot(data = dt[DIM == dim], x = ~funcId, y = ~ERT, showlegend = F, color = ~funcId, colors = get_color_scheme(seq_len(24))) %>%
    layout(yaxis = list(type = "log", ticklen = 3))
  dt3 <- rbindlist(lapply(seq_len(24), function(fid) {
    dt2 <- fread(paste0("ERTS/DataTable_F", fid, "_D", dim))
    erts <- dt2 %>% acast(algId.x ~ algId.y, value.var = 'total_ERT')
    vbss <- min(diag(erts))
    vbsdyn <- min(erts)
    funcId <- paste0('F', sprintf("%02d", fid))
    data.table(vbss, vbsdyn, funcId)
  }))
  p %<>% add_trace(data = dt3, x = ~funcId, y = ~vbss, mode = 'lines', name = 'VBS_Static', showlegend = T, color = 'black', linetypes = 'dash', linetype = 'dash', line = list(width = 3), type = 'scatter')
  p %<>% add_trace(data = dt3, x = ~funcId, y = ~vbsdyn, mode = 'lines', name = 'VBS_Dynamic', showlegend = T, color = 'black', linetypes = 'dot', linetype = 'dot', line = list(width = 3), type = 'scatter')
  # p <- plot_general_data(dt[DIM == dim], 'funcId', 'ERT', 'violin', 'funcId', scale.ylog = T, points = F)
  p %<>% layout(legend = list(x = 0, y = 1, orientation = 'v', font = list(size = getOption("IOHanalyzer.legend_fontsize", default = 18),
                family = 'Old Standard TT, serif')))
  save_plotly(p, paste0("../BBOB_Data/Figures/Boxplot_ERTs_dim", dim, ".pdf"), 700, 500)
}



### Number of 'valid' algs per func / dim
dt <- get_RT_summary(dsl_coco, c(Inf,1e-8))
x <- rbindlist(lapply(get_dim(dsl_coco), function(dim) {
  dt1 <- dt[target == Inf]
  dt2 <- dt[target == 1e-8]
  dt2 <- dt2[DIM == dim]
  dt1 <- dt1[DIM == dim]
  dt2[, x := (ps > 0)]
  dt2 <- dt2[dt1[['runs']] >= 15]
  dt2 <- dt2[, .(succ = sum(x)), by = funcId]
  dt1[runs >= 15]
  dt2[, DIM := paste0(sprintf("%02d", dim),'D')]
  dt2[, funcId := paste0('F', sprintf("%02d", funcId))]
}))
p <- plot_general_data(x, 'funcId', 'succ', legend_attr = 'funcId')
save_plotly(p, "../BBOB_Data/FiguresViolins_nr_finished_algs.pdf", 1000, 700)

p <- IOH_plot_ly_default(x.title = "Function ID", y.title = "Number of algorithms")
p %<>% add_trace(data = x, x = ~funcId, y = ~succ, color = ~DIM, marker = list(size = 10), type = 'scatter')
p %<>% layout(legend = list(y = 1.1, orientation = 'h',
                            font = list(size = getOption("IOHanalyzer.legend_fontsize", default = 18),
                                        family = 'Old Standard TT, serif')))
save_plotly(p, "../BBOB_Data/Figures/Scatter_nr_finished_algs_v2.pdf", 600, 400)

### Heatmap_imprs_dim_funcid generation
dt_split <- fread("../BBOB_Data/Split_COCO_V3.csv")
dt_split[, rel := best_ert / total_ERT]

#Relative impr
arr <- dt_split %>% acast(funcId ~ DIM, value.var = "impr")

#Relative ERT
arr <- dt_split %>% acast(funcId ~ DIM, value.var = "rel")

p <- plot_ly(y = paste0('F',rownames(arr)), x = paste0(colnames(arr),'D'), z = arr,
             type = 'heatmap', zauto = F, zmin = 1, zmax = 3, colorscale = 'Viridis')
p %<>% layout(yaxis = list(autorange = "reversed"))

# arr_percent <- arr * 100
for (fval in seq_len(24)) {
  yspot <- fval - 1
  # text <- arr_percent[fval, ] %>% format(digits = 3) %>% paste0(., '%')
  text <- arr[fval, ]
  i <- 1
  for (dval in paste0(colnames(arr),'D')) {
    tval <- text[[i]]
    i <- i + 1
    fcol <- if (tval < 2) 'white' else 'black'
    print(c(tval, fcol, tval < 2))
    p %<>% add_annotations(x = dval,
                          y = yspot,
                          text = as.character(format(tval, digits = 2)),
                          showarrow = FALSE,
                          font = list(color = fcol))
  }
}
p

#Relative impr
save_plotly(p, "../BBOB_Data/Heatmaps/Impr_per_funcid_dim_perc.pdf", 600, 600)

#Relative ERT
save_plotly(p, "../BBOB_Data/Heatmaps/Rel_ERT_per_funcid_dim_times_v6.pdf", 600, 600)

### Added 20/01

get_a1_a2_medians <- function(dim, use_method = "Median") {
  mds <- lapply(c('a1', 'a2'), function(which) {

    dt2 <- lapply(seq_len(24), function(x) {get_impr_A1(x, dim, which)}) %>% rbindlist

    arr <- dt2 %>% acast(algId ~ funcId, value.var = 'impr')

    nr_unfinished <- is.na(arr) %>% rowSums()
    algs <- names(nr_unfinished[nr_unfinished < 5]) # quantile(nr_finished, 0.2)])
    arr[arr > 3] <- 3
    arr[is.na(arr)] <- 3 #arbitrary
    if (use_method == "Median")
      mds <- rowMedians(arr)
    else
      mds <- rowMeans(arr)
    mds <- mds[rownames(arr) %in% algs]
    names(mds) <- algs
    mds
  })
  x <- data.table('algId' = names(mds[[1]]), 'I_1' = mds[[1]], 'I_2' = mds[[2]])

  x[, med := pmin(I_1, I_2) ]
  x <- x[order(med)][1:min(20, nrow(x))]
  x[, med := NULL]

  x <- melt(x, value.name = c("algId"))
  p <- IOH_plot_ly_default(x.title = "Algorithm", y.title = paste0(use_method, " I-value"))
  p %<>% add_trace(data = x, x = ~algId, y = ~`algId.1`, color = ~variable, marker = list(size = 10))
  p %<>% layout(yaxis = list(range = c(1,2)))
  save_plotly(p, paste0("../BBOB_Data/Figures/", use_method, "_I1_I2_", dim, "_dim.pdf"), 800, 500)
  # p
}

plot_heatmap_fid_dim <- function(fid, dim, subset = T) {
  dt <- fread(paste0("ERTS/DataTable_F", fid, "_D", dim))
  erts <- dt %>% acast(algId.x ~ algId.y, value.var = 'total_ERT')

  if (subset) {
    colrowmed <- colMedians(erts) + rowMedians(erts)
    names(colrowmed) <- rownames(erts)
    algs <- sort(colrowmed)[1:20] %>% names

    dt <- dt[algId.x %in% algs][algId.y %in% algs]

    erts <- dt %>% acast(algId.x ~ algId.y, value.var = 'total_ERT')
  }

  erts <- min(diag(erts)) / erts
  algs <- unique(dt$algId.x)
  p <- plot_ly(x = algs, y = algs, z = erts, type = 'heatmap', zmin = 0, zmax = 2,
               zauto = F, colorscale = 'RdBu')
  p %<>% layout(xaxis = list(title = "A2"), yaxis = list(title = "A1"))
  # p
  save_plotly(p, paste0("Heatmaps/All/F", fid, "_", dim, "D_v2.pdf"), 700, 700)
  # TRUE
}

### added 22/01

get_matrix_barplots <- function(dim){
  algs <- c("NELDERDOERR", "DE-AUTO",
            "BIPOP-aCMA-STEP", "HMLSL", "PSO-BFGS")

  dt_full <- rbindlist(lapply(seq_len(24), function(fid) {
    dt <- fread(paste0("ERTS/DataTable_F", fid, "_D", dim))
    dt_temp <- dt[algId.x %in% algs][algId.y %in% algs]
    min_ert <- dt_temp %>% acast(algId.x ~ algId.y, value.var = 'total_ERT') %>% diag %>% min
    dt_temp[, impr := min_ert / total_ERT]
    dt_temp[, fid := fid]
    dt_temp
  }))

  plist <- apply(expand.grid(algs, rev(algs)), 1, function(x) {
    dt_temp <- dt_full[algId.x == x[[1]]][algId.y == x[[2]]]
    dt_temp[, i2 := pmax(impr, 0)]
    p_temp <- plot_ly(data = dt_temp, x = ~fid, y = ~i2, type = 'bar')
    p_temp %>% layout(showlegend = FALSE, yaxis = list(linecolor = toRGB("black"),
                                                       linewidth = 1, mirror = "allticks",
                                                       range = c(1, 2), title = x[[2]], showline = T),
                      xaxis = list(linecolor = toRGB("black"),
                                   linewidth = 1, mirror = "allticks", title = x[[1]],
                                   showline = T, gridcolor = '#D1CDCC'))
  })

  p <- subplot(plist, nrows = length(algs), shareX = T, shareY = T, titleX = T, titleY = T, margin = 0)

  save_plotly(p, paste0("Figures/Matrix_of_bars_", dim, "D_v3.pdf"), 1000, 1000)
  # p
}

for (dim in get_dim(dsl_coco)) { get_matrix_barplots(dim)}


### Calculate SBS per dim
dt <- get_RT_summary(dsl_coco, c(1e-8, Inf))

sbss <- lapply(get_dim(dsl_coco), function(dim) {
  dt1 <- dt[target == Inf]
  dt2 <- dt[target == 1e-8]
  dt2 <- dt2[DIM == dim]
  dt1 <- dt1[DIM == dim]
  dt2[, x := (ps > 0)]
  dt2 <- dt2[dt1[['runs']] >= 15]
  dt2[, rank := frank(ERT, na.last = T), by = .(funcId)]
  dt3 <- dt2[, .(succ = sum(x), rsum = sum(rank)), by = c('algId')]
  dt3 <- dt3[succ == 24]
  sbs <- dt3[which.min(dt3[['rsum']]), ][['algId']]
  sbs
})

y <- rbindlist(lapply(c(2,3,5,10,20), function(dim) {
  dt_temp <- dt_split[DIM == dim]
  vbs_erts <- dt_temp[, c('best_ert', 'funcId')]
  sbs_erts <- dt[target == 1e-8][DIM == dim][algId == sbss[[paste0(dim, 'D')]]][, c('ERT', 'funcId')]
  dt4 <- merge(vbs_erts, sbs_erts)
  dt4[, DIM := dim]
  dt4[, rel := ERT / best_ert]
}))

arr <- y %>% acast(funcId ~ DIM, value.var = "rel")

p <- plot_ly(y = paste0('F',rownames(arr)), x = paste0(colnames(arr),'D'), z = arr,
             type = 'heatmap', zauto = F, zmin = 1, zmax = 3, colorscale = 'Viridis')
p %<>% layout(yaxis = list(autorange = "reversed"))

for (fval in seq_len(24)) {
  yspot <- fval - 1
  text <- arr[fval, ]
  i <- 1
  for (dval in paste0(colnames(arr),'D')) {
    tval <- text[[i]]
    i <- i + 1
    fcol <- if (tval < 2) 'white' else 'black'
    print(c(tval, fcol, tval < 2))
    p %<>% add_annotations(x = dval,
                           y = yspot,
                           text = as.character(format(tval, digits = 2)),
                           showarrow = FALSE,
                           font = list(color = fcol))
  }
}
p

save_plotly(p, "../BBOB_Data/Heatmaps/Rel_ERT_SBS_over_VBS.pdf", 600, 600)


get_i1_i2_box <- function(dim) {
  algs <- lapply(c('a1', 'a2'), function(which) {

    dt5 <- lapply(seq_len(24), function(x) {get_impr_A1(x, dim, which)}) %>% rbindlist

    arr <- dt5 %>% acast(algId ~ funcId, value.var = 'impr')

    nr_unfinished <- is.na(arr) %>% rowSums()
    algs <- names(nr_unfinished[nr_unfinished < 5]) # quantile(nr_finished, 0.2)])

    arr[arr > 3] <- 3
    arr[is.na(arr)] <- 3 #arbitrary
    mds <- rowMeans(arr)
    mds <- mds[algs]
    algs <- mds[order(mds)][1:10] %>% names
    algs
  }) %>% unlist %>% unique

  dt51 <- lapply(seq_len(24), function(x) {get_impr_A1(x, dim, 'a1')}) %>% rbindlist
  dt52 <- lapply(seq_len(24), function(x) {get_impr_A1(x, dim, 'a2')}) %>% rbindlist

  dt51[, type := 'I_1']
  dt51 <- dt51[algId %in% algs]
  dt52[, type := 'I_2']
  dt52 <- dt52[algId %in% algs]

  dt6 <- rbindlist(list(dt51, dt52))

  p <- IOH_plot_ly_default(x.title = "Algorithm", y.title = " I-value")
  p %<>% add_boxplot(data = dt6, x = ~algId, y = ~impr, color = ~type) %>% layout(boxmode = 'group')
  p %<>% layout(yaxis = list(range = c(1,3)))
  save_plotly(p, paste0("../BBOB_Data/Figures/boxplot_I1_I2_", dim, "_dim.pdf"), 800, 500)
  # p
}

### Plot of Fixed-target behavior
algs <- c("NELDERDOERR", "DE-AUTO", "BIPOP-aCMA-STEP", "HMLSL", "PSO-BFGS")

dsl_temp <- subset(dsl_coco, algId %in% algs)

dim <- 3
fid <- 12
dsl_sub <- subset(dsl_temp, funcId == fid && DIM == dim)

splits <- rev(seq_FV(c(100,1e-8), length.out = 51, scale = 'log'))
dt_ert <- get_RT_summary(dsl_sub, splits)

p <- plot_general_data(dt_ert, 'target', 'ERT', 'line',
                       scale.xlog = T, scale.ylog = T, scale.reverse = T,
                       x_title = "Target", show.legend = T)
p
temp <- fread("ERTS/DataTable_F12_D3")
temp <- temp[algId.x %in% algs][algId.y %in% algs]

x <- list()

x <- rbindlist(lapply(seq_len(25), function(i) {
  te <- temp[i,]
  y <- dt_ert[algId == te$algId.x][round(target, 9) == round(te$target, 9)][['ERT']]
  te$y <- y
  te
}))

symbols <- c('star', 'triangle-up', 'x', 'square', 'pentagon')
names(symbols) <- algs
for (alg in algs) {
  symbol <- symbols[[alg]]
  x1 <- x[algId.y == alg]
  p %<>% add_trace(data = x1, x = ~target, y = ~y, color = ~algId.x,
                 marker = list(size = 13, symbol = symbol), type = 'scatter', showlegend = F)
}
p %>% save_plotly("Figures/fixed_target_curve_with_switchmarks_F12_3D.pdf", 700, 500)


### Table for DIM 5
dtl <- dt_split[DIM == 5]
dtl <- dtl[, list(FID = funcId, VBS_static = best_alg, VBS_static_ERT = best_ert, A1 = p1_alg, A2 = p2_alg, tau = log10(target), total_ERT, speedup_factor = rel)]
dtl_x <- xtable(dtl, "Caption", "Dim5_overview_table", NULL, 2)
digits(dtl_x) <- c(0, 0, 0, 1, 1, 1, 1, 1, 2)
print(dtl_x, include.rownames = F)

