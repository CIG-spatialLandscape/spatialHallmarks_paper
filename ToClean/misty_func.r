

################################################################################

#General pipeline for misty model
run_misty_seurat <- function(visium.slide,
                             # Seurat object with spatial transcriptomics data.
                             view.assays,
                             # Named list of assays for each view.
                             view.features = NULL,
                             # Named list of features/markers to use.
                             # Use all by default.
                             view.types,
                             # Named list of the type of view to construct
                             # from the assay.
                             view.params,
                             # Named list with parameters (NULL or value)
                             # for each view.
                             spot.ids = NULL,
                             # spot IDs to use. Use all by default.
                             out.alias = "results"
                             # folder name for output
) {
  
  # Extracting geometry
  geometry <- GetTissueCoordinates(visium.slide,
                                   cols = c("row", "col"), scale = NULL
  )
  
  # Extracting data
  view.data <- map(view.assays,
                   extract_seurat_data,
                   geometry = geometry,
                   visium.slide = visium.slide
  )
  
  # Constructing and running a workflow
  build_misty_pipeline(
    view.data = view.data,
    view.features = view.features,
    view.types = view.types,
    view.params = view.params,
    geometry = geometry,
    spot.ids = spot.ids,
    out.alias = out.alias
  )
}


################################################################################

# Extracts data from an specific assay from a Seurat object
# and aligns the IDs to the geometry
extract_seurat_data <- function(visium.slide,
                                assay,
                                geometry) {
  data <- GetAssayData(visium.slide, assay = assay) %>%
    t() %>%
    as_tibble(rownames = NA)
  
  return(data %>% slice(match(rownames(.), rownames(geometry))))
}

################################################################################

# Filters data to contain only features of interest
filter_data_features <- function(data,
                                 features) {
  if (is.null(features)) features <- colnames(data)
  
  return(data %>% rownames_to_column() %>%
           select(rowname, all_of(features)) %>% rename_with(make.names) %>%
           column_to_rownames())
}

################################################################################

# Builds views depending on the paramaters defined
create_default_views <- function(data,
                                 view.type,
                                 view.param,
                                 view.name,
                                 spot.ids,
                                 geometry) {
  view.data.init <- create_initial_view(data)
  
  if (!(view.type %in% c("intra", "para", "juxta"))) {
    view.type <- "intra"
  }
  
  if (view.type == "intra") {
    data.red <- view.data.tmp$data %>%
      rownames_to_column() %>%
      filter(rowname %in% spot.ids) %>%
      select(-rowname)
  } else if (view.type == "para") {
    view.data.tmp <- view.data.init %>%
      add_paraview(geometry, l = view.param)
    
    data.ix <- paste0("paraview.", view.param)
    data.red <- view.data.tmp[[data.ix]]$data %>%
      mutate(rowname = rownames(data)) %>%
      filter(rowname %in% spot.ids) %>%
      select(-rowname)
  } else if (view.type == "juxta") {
    view.data.tmp <- view.data.init %>%
      add_juxtaview(
        positions = geometry,
        neighbor.thr = view.param
      )
    
    data.ix <- paste0("juxtaview.", view.param)
    data.red <- view.data.tmp[[data.ix]]$data %>%
      mutate(rowname = rownames(data)) %>%
      filter(rowname %in% spot.ids) %>%
      select(-rowname)
  }
  
  if (is.null(view.param) == TRUE) {
    misty.view <- create_view(
      paste0(view.name),
      data.red
    )
  } else {
    misty.view <- create_view(
      paste0(view.name, "_", view.param),
      data.red
    )
  }
  
  return(misty.view)
}

################################################################################

#misty pipline for just a pathway as intra
build_misty_pipeline <- function(view.data,
                                 view.features,
                                 view.types,
                                 view.params,
                                 geometry,
                                 spot.ids = NULL,
                                 out.alias = "default") {
  
  # Adding all spots ids in case they are not defined
  if (is.null(spot.ids)) {
    spot.ids <- rownames(view.data[[1]])
  }
  
  # First filter the features from the data
  view.data.filt <- map2(view.data, view.features, filter_data_features)
  
  # Create initial view
  views.main <- create_initial_view(view.data.filt[[1]] %>%
                                      rownames_to_column() %>%
                                      filter(rowname %in% spot.ids) %>%
                                      select(-rowname))
  
  # Create other views
  view.names <- names(view.data.filt)
  
  all.views <- pmap(list(
    view.data.filt[-1],
    view.types[-1],
    view.params[-1],
    view.names[-1]
  ),
  create_default_views,
  spot.ids = spot.ids,
  geometry = geometry
  )
  
  pline.views <- add_views(
    views.main,
    unlist(all.views, recursive = FALSE)
  )
  
  
  # Run MISTy
  run_misty(pline.views, out.alias)
}

################################################################################

#misty piple for seurat object
para_ppln_seurat = function(visium_slide,
                            intra_assay, 
                            intra_features = NULL,
                            para_assay,
                            para_features = NULL,
                            l,
                            spot_ids = NULL,
                            out_alias = "default"){
  
  # Getting data ready to create views
  geometry = visium_slide@images$prostate@coordinates[,c(2,3)]
  
  # Intrinsic data
  intra_df = as.matrix(visium_slide@assays[[intra_assay]]@data)
  intra_df = intra_df %>%
    t %>% data.frame(check.names = F)
  
  intra_df = intra_df[rownames(geometry),]
  
  # Para data
  para_df = as.matrix(visium_slide@assays[[para_assay]]@data)
  para_df = para_df %>%
    t %>% data.frame(check.names = F)
  para_df = para_df[rownames(geometry),]
  
  
  para_pipeline(intra_df = intra_df,
                intra_features = intra_features,
                para_df = para_df,
                para_features = para_features,
                geometry = geometry,
                l = l,
                spot_ids = spot_ids,
                out_alias = out_alias)
  
}

################################################################################

para_pipeline = function(intra_df, 
                         intra_features = NULL,
                         para_df,
                         para_features = NULL,
                         geometry,
                         l,
                         spot_ids = NULL,
                         out_alias = "default"){
  
  plan(multiprocess, workers = 4)
  
  clear_cache()
  
  if(is.null(intra_features)){
    intra_features = colnames(intra_df)
  }
  
  if(is.null(para_features)){
    para_features = colnames(para_df)
  }
  
  if(is.null(spot_ids)){
    spot_ids = rownames(intra_df)
  }
  
  # Defining useful data intra
  intra_df = intra_df[spot_ids,intra_features]
  colnames(intra_df) = gsub("-","_", colnames(intra_df))
  
  views_main = create_initial_view(intra_df, 
                                   unique.id = "intra")
  
  
  # Defining useful data para
  para_df = para_df[rownames(geometry),para_features]
  
  colnames(para_df) = gsub("-","_", colnames(para_df))
  
  views_para = create_initial_view(para_df, 
                                   unique.id = paste0("para_",l)) %>% 
    add_paraview(geometry, l)
  
  # Fetching actual para info to be used
  
  # Spot specific view comes from the view above
  data_red = views_para[[3]]$data
  rownames(data_red) = rownames(para_df) #we named rows just for easy access
  data_red = data_red[rownames(intra_df),]
  
  # Define frankenstein of views
  views = views_main %>% 
    add_views(create_view(paste0("para_",l),
                          data_red))
  
  MISTy_run = run_misty(views,paste0(out_alias,"_",l))
  
}

################################################################################

MISTy_aggregator = function(results_folder,
                            p.cutoff = 0.05){
  
  images = results_folder
  
  # Improvement 
  impr = images %>% map_dfc(function(image) {
    performance <- read_delim(paste0(image, .Platform$file.sep,
                                     "performance.txt"),
                              delim = " ", col_types = cols()
    ) %>% distinct()
    
    targets <<- unique(performance$target)
    performance %>%
      arrange(target) %>%
      transmute(RMSE = (intra.RMSE - multi.RMSE) / intra.RMSE, R2
                = (multi.R2 - intra.R2))
  })
  
  avg <- ((images %>% map(function(image) {
    coefficients <- read_delim(paste0(image, .Platform$file.sep, "coefficients.txt"),
                               delim = " ", col_types = cols()
    )
    
    targets <<- coefficients %>%
      pull(target) %>%
      sort
    
    coefficients %>%
      distinct() %>%
      arrange(target) %>%
      dplyr::select(-target, -contains("intercept")) %>%
      mutate_at(vars(starts_with("p.")), ~ as.numeric(. <= p.cutoff)) %>%
      mutate_at(vars(-starts_with("p.")), abs)
  }) %>% 
    purrr::reduce(`+`)) / length(images)) %>% 
    mutate(target = targets)
  
  ctotals <- avg %>% 
    dplyr::select(-starts_with("p."), -"target") %>%
    rowSums
  
  coefs <- avg %>% 
    dplyr::select(-starts_with("p.")) %>%
    mutate_if(is.numeric, ~./ctotals) %>%
    tidyr::pivot_longer(-target, names_to = "view")
  
  maps = images %>% map(function(image) {
    coefficients <- read_delim(paste0(image, .Platform$file.sep, "coefficients.txt"),
                               delim = " ", col_types = cols()
    ) %>% distinct()
    
    targets <- unique(coefficients$target)
    views <<- (coefficients %>% dplyr::select(-target, -starts_with("p."), -intercept) %>% colnames())
    
    # one heatmap per view
    maps <- views %>% map(function(view) {
      all.importances <- targets %>% map(~ read_csv(paste0(
        image, .Platform$file.sep, "importances_",
        .x, "_", view, ".txt"
      ),
      col_types = cols()
      ) %>%
        distinct() %>%
        filter(!grepl("_2$", target)))
      
      features <- unique(all.importances %>% map(~ .x$target) %>% unlist())
      
      pview <- paste0("p.", view)
      ps <- coefficients %>%
        dplyr::select(target, !!pview) %>%
        mutate(!!pview := (1 - !!sym(pview)))
      
      
      # importances are standardized for each target an multiplied by 1-pval(view)
      result <- all.importances %>%
        imap_dfc(~
                   tibble(target = features, zero.imp = 0) %>%
                   left_join(.x, by = "target") %>%
                   transmute(feature = target, importance = (zero.imp + scale(imp)[, 1]) *
                               (ps %>% filter(target == targets[.y]) %>% pull(pview))) %>%
                   dplyr::select(importance)) %>%
        `colnames<-`(targets) %>%
        mutate(Predictor = features)
      
      # in order for aggregation
      result %>%
        arrange(Predictor) %>%
        dplyr::select((order(colnames(result))))
      #dplyr::select(noquote(order(colnames(result))))
    })
  })
  
  aggregated = maps %>% purrr::reduce(function(acc, l) {
    map2(acc, l, ~ (((.x %>% dplyr::select(-Predictor)) + (.y %>% dplyr::select(-Predictor))) %>%
                      mutate(Predictor = .x %>% pull(Predictor))))
  })
  
  performance_df = read_delim(paste0(results_folder,"/performance.txt"), delim = " ")
  
  result_list = list("impr" = impr,
                     "targets" = targets,
                     "coefs" = coefs,
                     "importance" = aggregated,
                     "performance" = performance_df)
  
  return(result_list)
  
}

################################################################################



plot_misty_performance = function(MISTy_out,
                                  predicted_features = NULL){
  #Generating figures
  #Overall performance
  performance_df = MISTy_out$performance %>% tidyr::pivot_longer(cols = -target) %>%
    dplyr::filter(grepl("R2",name))
  
  R2_impr = MISTy_out$impr %>% dplyr::select(contains("R2")) %>%
    mutate(target = targets) %>%
    pivot_longer(cols = -target, 
                 names_to = "name", 
                 values_to = "value") %>% 
    arrange(desc(value))
  
  #Look at importances of views
  contribution_df = MISTy_out$coefs
  
  if(length(predicted_features)>0){
    
    performance_df = performance_df %>%
      dplyr::filter(target %in% predicted_features)
    
    R2_impr = R2_impr %>%
      dplyr::filter(target %in% predicted_features)
    
    contribution_df = contribution_df %>%
      dplyr::filter(target %in% predicted_features)
    
  }
  
  performance_plt = ggplot(performance_df,
                           aes(fill = name, y = target, x = value)) +
    geom_bar(stat = "identity",position="dodge") + 
    theme_minimal()
  
  impr_plot = ggplot(R2_impr) +
    geom_point(aes(x = target, y = value * 100)) +
    theme_classic() +
    ylab("Change in variance explained") +
    xlab("Target") +
    theme(axis.title = element_text(size=11),
          axis.text = element_text(size=10),
          axis.text.x = element_text(angle = 90, hjust = 1))
  
  coefs_plot = ggplot(contribution_df) + 
    geom_col(aes(x=target, 
                 y=value, group=view, fill=view)) +
    xlab("Target") +
    ylab("Contribution") +
    theme_classic() +
    theme(axis.text = element_text(size=13),
          axis.title = element_text(size=13),
          legend.text = element_text(size=11)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1,
                                     vjust = 0.5)) +
    scale_fill_brewer(palette="Dark2",
                      labels = c("Intrinsic","Para pathway")) +
    labs(fill = "View")
  
  plot(performance_plt)
  plot(impr_plot)
  plot(coefs_plot)
}

cluster_importance = function(importance_obj,
                              predictors_vars,
                              predicted_vars){
  
  importance_intra_mat = importance_obj %>% 
    dplyr::filter(Predictor %in% predictors_vars) %>%
    as.data.frame()
  
  rownames(importance_intra_mat) = importance_intra_mat$Predictor
  importance_intra_mat = as.matrix(importance_intra_mat[,-which(colnames(importance_intra_mat) == "Predictor")])
  importance_intra_mat = importance_intra_mat[,predicted_vars]
  
  predictor_clust =  hclust(d = t(dist(importance_intra_mat)))
  predicted_clust =  hclust(d = t(dist(t(importance_intra_mat))))
  
  return(list("Predictor_order" = predictor_clust$labels[predictor_clust$order],
              "Predicted_order" =  predicted_clust$labels[predicted_clust$order]))
  
  
}

################################################################################

plot_misty_importance = function(MISTy_out,
                                 predicted_features = NULL,
                                 predictors_features = NULL,
                                 importance_cut = 1,
                                 make_clust = T,
                                 width_pdf = 15,
                                 height_pdf = 13){
  
  # We have a list of importance heatmaps in the MISTy_out file
  n_views = length(MISTy_out$importance)
  importance_out = list()
  
  for(i in 1:n_views){
    
    importance_df = tidyr::gather(MISTy_out$importance[[i]], 
                                  "Predicted",
                                  "Importance", -Predictor)
    
    predicted_features_v = predicted_features
    
    predictors_features_v = predictors_features[[i]]
    
    if(length(predicted_features_v)>0){
      importance_df = importance_df %>% 
        dplyr::filter(Predicted %in% predicted_features_v)
    }
    
    if(length(predictors_features_v)>0){
      importance_df = importance_df %>% 
        dplyr::filter(Predictor %in% predictors_features_v)
    }
    
    #Here we extract predictors that have a predictor importance of at least importance_cut
    
    importance_summ = importance_df %>% 
      dplyr::mutate(importance_bool = Importance >= importance_cut) %>%
      group_by(Predictor) %>% 
      summarize(predictor_summ = sum(importance_bool,na.rm = T)) %>%
      dplyr::filter(predictor_summ >= 1)
    
    importance_predictors = importance_summ %>% select(Predictor) %>% pull()
    
    importance_df = importance_df %>% 
      dplyr::filter(Predictor %in% importance_predictors) %>%
      dplyr::arrange(Predicted, -Importance)
    
    # Clustering the plots
    
    if(make_clust){
      clust_order = cluster_importance(importance_obj = MISTy_out$importance[[i]],
                                       predictors_vars = unique(importance_predictors),
                                       predicted_vars = unique(importance_df$Predicted))
      
      imp_plot = importance_df %>% 
        ggplot(aes(x = factor(Predictor,
                              levels = clust_order$Predictor_order), 
                   y = factor(Predicted,
                              levels = clust_order$Predicted_order) ,
                   fill = Importance)) + geom_tile() + 
        theme(panel.grid = element_blank(),
              axis.title = element_text(size=11),
              axis.text.x = element_text(angle = 90,
                                         hjust = 1,
                                         vjust = 0.5),
              axis.text = element_text(size=10),
              legend.key.size = unit(.6, "cm"),
              legend.title = element_text(size = 10),
              legend.text = element_text(size = 10),
              legend.position = "bottom") +
        scale_fill_gradient2(low = "white", 
                             mid = "white", 
                             high = scales::muted("blue"),
                             midpoint = 0.9) +
        xlab("Predictor features") + ylab("Predicted markers")
    }else{
      
      imp_plot = importance_df %>% 
        ggplot(aes(x = Predictor, 
                   y = Predicted,
                   fill = Importance)) + geom_tile() + 
        theme(panel.grid = element_blank(),
              axis.title = element_text(size=11),
              axis.text.x = element_text(angle = 90,
                                         hjust = 1,
                                         vjust = 0.5),
              axis.text = element_text(size=10),
              legend.key.size = unit(.6, "cm"),
              legend.title = element_text(size = 10),
              legend.text = element_text(size = 10),
              legend.position = "bottom") +
        scale_fill_gradient2(low = "white", 
                             mid = "white", 
                             high = scales::muted("blue"),
                             midpoint = 0.9) +
        xlab("Predictor features") + ylab("Predicted markers")
      
    }
    
    plot(imp_plot)
    
    importance_out[[i]] = importance_df
    
  }
  return(importance_out)
}

################################################################################

get_optimal = function(out_dir_name,
                       ls){
  
  system(paste0("mkdir ", out_dir_name,"_optim"))
  
  l = ls
  
  perf = (l) %>% map_dfr(function(p){
    performance <- read_delim(paste0(out_dir_name, "_",p, "/performance.txt"),
                              delim = " ", col_types = cols()
    ) %>% distinct()
    
    #selection criterion maximum R2 performance improvement per marker
    performance %>% arrange(target) %>%
      mutate(impr = multi.R2 - intra.R2) %>% 
      dplyr::select(target,impr) %>% 
      column_to_rownames("target") %>% t %>% as.data.frame
  })
  
  # For each target get the maximum of improvement
  # distinct ls
  optimal.l = colnames(perf) %>% enframe(name = NULL, value = "target") %>% 
    mutate(l = apply(perf, 2, which.max) %>% (l)[.])
  
  # Copy the relevant documents to the new directory, deleting the l parameter info
  optimal.l %>% pull(target) %>% walk2(optimal.l %>% pull(l), function(.x, .y){
    files <- list.files(paste0(out_dir_name,"_", .y, "/"), 
                        paste0("importances_", .x, '*'), 
                        full.names = TRUE)
    
    files %>% walk(~file.copy(., paste0(out_dir_name,"_optim/", 
                                        str_replace(last(str_split(., "/")[[1]]), 
                                                    "a([\\._][0-9]+)+", "a"))))
  })
  
  #very suboptimal
  optimal.l %>% pull(target) %>% walk2(optimal.l %>% pull(l), function(.x, .y){
    performance <- read_delim(paste0(out_dir_name,"_", .y, "/performance.txt"),
                              delim = " ", col_types = cols()
    ) %>% distinct()
    
    coeff <- read_delim(paste0(out_dir_name,"_", .y, "/coefficients.txt"),
                        delim = " ", col_types = cols()
    ) %>% distinct()
    
    
    
    if(!file.exists(paste0(out_dir_name,"_optim/performance.txt"))){
      write(colnames(performance) %>% paste0(collapse=" "), 
            paste0(out_dir_name,"_optim/performance.txt"))
    }
    
    if(!file.exists(paste0(out_dir_name,"_optim/coefficients.txt"))){
      write(str_remove_all(colnames(coeff) %>% paste0(collapse=" "), "([\\._][0-9]+)+"), 
            paste0(out_dir_name,"_optim/coefficients.txt"))
    }
    
    write(performance %>% filter(target == .x) %>% unlist %>% unname %>% paste0(collapse=" "), 
          paste0(out_dir_name,"_optim/performance.txt"), append = T)
    
    write(coeff %>% filter(target == .x) %>% unlist %>% unname %>% paste0(collapse=" "), 
          paste0(out_dir_name,"_optim/coefficients.txt"), append = T)
    
  })
  
}

################################################################################


# data collection functions
get_niche_views <- function(visium.slide,
                            view.assays,
                            view.features = NULL,
                            view.types,
                            view.params,
                            spot.ids = NULL) {
  
  
  # Extracting geometry
  geometry <- GetTissueCoordinates(visium.slide,
                                   cols = c("row", "col"), scale = NULL
  )
  
  # Extracting data
  view.data <- map(view.assays,
                   extract_seurat_data,
                   geometry = geometry,
                   visium.slide = visium.slide
  )
  
  get_misty_views(
    view.data = view.data,
    view.features = view.features,
    view.types = view.types,
    view.params = view.params,
    geometry = geometry,
    spot.ids = spot.ids
  )
  
}

get_blocked_misty_views <- function(slide_path,
                                    assay_name,
                                    spot_group,
                                    group_id,
                                    # Define spatial context of each view -----
                                    view_types = list("main" = "intra", 
                                                      "juxta" = "juxta",
                                                      "para" = "para"),
                                    # Define additional parameters (l in case of paraview,
                                    # n of neighbors in case of juxta) --------
                                    view_params = list("main" = NULL, 
                                                       "juxta" = 5,
                                                       "para" = 15),
                                    feature_exception = NULL) {
  
  print(slide_path)
  # Reading data
  visium_slide <- readRDS(slide_path)
  
  # This identifies the specific spots to be used
  Idents(visium_slide) <- spot_group
  spots <- WhichCells(visium_slide, idents = group_id)
  
  # Get useful features by default all the features in the assay are included
  cell_feats <- rownames(visium_slide@assays[[assay_name]])
  if(!is.null(feature_exception)) {
    cell_feats <- cell_feats[!cell_feats %in% feature_exception]
  }
  
  # Define assay of each view ---------------
  view_assays <- list("main" = assay_name,
                      "juxta" = assay_name,
                      "para" = assay_name)
  
  # Define features of each view ------------
  view_features <- list("main" = cell_feats, 
                        "juxta" = cell_feats,
                        "para" = cell_feats)
  
  get_niche_views(visium.slide = visium_slide,
                  view.assays = view_assays,
                  view.features = view_features,
                  view.types = view_types,
                  view.params = view_params,
                  spot.ids = spots)
  
}

niche_misty_pplne <- function(niche_selection, 
                              view_names, 
                              out_alias = "./visium_results_manuscript/niche_misty/") {
  
  print(niche_selection)
  
  # Get view data from each slide for the selected niche
  niche_dat <- visium_df %>%
    group_by(niche) %>%
    nest() %>%
    dplyr::filter(niche == niche_selection) %>%
    mutate(data = map2(niche, data, function(niche_id, df) {
      
      slide_views <- df %>%
        mutate(niche_data = map(prostate_if, get_blocked_misty_views,
                                assay_name = "c2l_major_props",
                                spot_group = "opt_clust_integrated",
                                group_id = niche_id,
                                # Define spatial context of each view -----
                                view_types = list("main" = "intra", 
                                                  "juxta" = "juxta",
                                                  "para" = "para"),
                                # Define additional parameters (l in case of paraview,
                                # n of neighbors in case of juxta) --------
                                view_params = list("main" = NULL, 
                                                   "juxta" = 5,
                                                   "para" = 15),
                                feature_exception = NULL))
      
    }
    ))
  
  # Second step is to create alternative joint views:
  # First per view name you extract data
  
  niche_dat <- niche_dat %>%
    mutate(data = map(data, function(slide_dat) {
      
      view_dat <- map(view_names, function(vn) { 
        
        print(vn)
        
        map(slide_dat$niche_data, function(misty_views) {
          
          first_level <- misty_views[[vn]]
          second_level_names <- names(first_level)
          second_level_ix <- second_level_names[grepl(vn, second_level_names)]
          second_level <- first_level[[second_level_ix]]
          second_level[["data"]]
          
        })
        
      })
      
    })) %>% 
    unnest_wider("data")
  
  # Third step is to format the joint views for MISTy
  main = bind_rows(niche_dat$main)
  juxta = bind_rows(niche_dat$juxta)
  para = bind_rows(niche_dat$para)
  
  views_main <- create_initial_view(main, unique.id = NULL)
  
  other_views <- list("juxta" = create_view(name = "juxta", 
                                            data = juxta,
                                            abbrev = "juxta"),
                      "para" = create_view(name = "para", 
                                           data = para, 
                                           abbrev = "para"))
  
  views_main <- add_views(
    views_main,
    unlist(other_views, recursive = FALSE)
  )
  
  # Define path
  misty_out <- paste0(out_alias, niche_selection)
  
  run_misty(views_main, results.folder = misty_out)
  
  
}
