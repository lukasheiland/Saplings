
# combineData -------------------------------------------------------------
combineData <- function(B, S, E) {
  
  regclass_smaller <- c("h[20,50)", "h[50,130)")
  S %<>% 
    dplyr::filter(regclass %in% regclass_smaller) %>%
    dplyr::filter(obsid == "DE_BWI_2012") %>%
    mutate(count_ha = count/countarea) %>%
    group_by(plotid, taxid, tax, obsid, plotobsid, methodid, time) %>%
    summarize(count_ha = sum(count_ha), .groups = "drop")
  
  # taxa_small <- unique(S$tax)
  
  B %<>%
    dplyr::filter(obsid == "DE_BWI_2012") %>%
    # mutate(istaxoninsmall = tax %in% taxa_small) %>%
    group_by(plotid, tax) %>%
    mutate(ba = pi * (dbh/2)^2 * 1e-6) %>% # mm^2 to m^2)
    mutate(count_ha = count/countarea) %>% ## just verbose, actually all (B$countarea == 1)
    mutate(ba_ha = ba * count_ha) %>% ## m^2 ha^-1 per tree
    dplyr::summarize(ba_ha_species = sum(ba_ha, na.rm = T), .groups = "drop") %>%
    
    group_by(plotid) %>%
    mutate(ba_ha_plot = sum(ba_ha_species, na.rm = T))
  
  D <- S %>% 
    dplyr::left_join(B, by = c("plotid", "tax")) ## joining with these two ensures species identity for ba_ha_species

  D %<>%
    dplyr::left_join(E, by = "plotid") %>%
    st_as_sf() %>%
    bind_cols(st_coordinates(.)) %>%
    st_drop_geometry() %>%
    droplevels()
  
  return(D)
}



# subsetData -------------------------------------------------------------
subsetData <- function(species, D, predictor, n_plots = NULL) {
  
  D %<>% filter(tax == species) %>%
    dplyr::select(all_of(c(predictor, "count_ha", "plotid", "X", "Y"))) %>% ## selection for drop_na
    tidyr::drop_na()
  
  plots <- unique(D$plotid) 
  
  if(!is.null(n_plots) & n_plots < length(plots)) {
    plots %<>% sample(n_plots)
    D %<>% dplyr::filter(plotid %in% plots)
  }
  
  Scalings <- data.frame(predictor = predictor, mean = colMeans(D[predictor]), sd = apply(D[predictor], 2, FUN = sd))
  
  D %<>% mutate_at(predictor, function(x) c(scale(x))) ## !!! to have vectors as data.frame columns
  
  attr(D, "Scalings") <- Scalings
  attr(D, "species") <- species
  
  return(D)
}



# fitModel -------------------------------------------------------------
fitModel <- function(Data, predictor_fixed, predictor_randomfield) {
  
  mesh <- make_mesh(as.data.frame(Data), xy_cols = c("X", "Y"), n_knots = 200)
  
  predictor_all <- c(predictor_fixed, predictor_randomfield)
  form_randomfield <- as.formula(paste("~ 0 + ", paste(predictor_randomfield, collapse = "+")))
  form_fixed <- as.formula(paste("count_ha ~ 1 +", paste(predictor_all, collapse = "+")))
  
  fit <- sdmTMB(formula = form_fixed, spatial_varying = form_randomfield,
                mesh = mesh, data = Data, family = nbinom2(), spatial = TRUE)
  
  attr(fit, "Scalings") <- attr(Data, "Scalings")
  attr(fit, "species") <- attr(Data, "species")
  
  return(fit)
}



# predictFit -------------------------------------------------------------
predictFit <- function(fit, Newdata) {
  
  
  Scalings <- attr(fit, "Scalings")
  
  predictor_fit <- Scalings$predictor
  predictor_stack <- names(Stack)
  predictor_diff <- setdiff(predictor_fit, predictor_stack)
  predictor_common <- intersect(predictor_fit, predictor_stack)
  Scalings_sub <- Scalings[na.exclude(match(Scalings$predictor, predictor_common)),]

  ## For method with RasterStack
  # Stack <- raster::subset(Stack, predictor_common)
  # Stack_scaled <- raster::scale(Stack, center = Scalings_sub$mean, scale = Scalings_sub$sd)
  
  Scaled <- scale(Newdata[,predictor_common], center = Scalings_sub$mean, scale = Scalings_sub$sd)
  
  Newdata %<>%
    dplyr::select(-all_of(predictor_common)) %>%
    dplyr::bind_cols(Scaled) %>%
    dplyr::select(all_of(c(predictor_common, "X", "Y"))) %>%
    dplyr::mutate(ba_ha_plot = 0, ba_ha_species = 0) %>%
    tidyr::drop_na()
  
  P <- predict(fit, newdata = Newdata, type = "response")

  attr(P, "species") <- attr(fit, "species")
  return(P)
}



# plotMap -------------------------------------------------------------
plotMap <- function(Predictions, var = est) {
  
  plot <- ggplot(Predictions, aes(X, Y, fill = {{ var }})) +
    geom_raster() +
    coord_fixed() +
    scale_fill_viridis_c(name = "Response") +
    theme_minimal() +
    ggtitle(attr(Predictions, "species"))
  
  return(plot)
}

