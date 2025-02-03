
#' Prepare for running an SDM
#'
#' The background sampling includes code based on a
#' [Geographic Information Systems stack exchange](https://gis.stackexchange.com/)
#' [answer](https://gis.stackexchange.com/a/224347)
#' by user [Spacedman](https://gis.stackexchange.com/users/865/spacedman).
#'
#' Options for managing memory are `terra_options`, `max_cells_in_memory` and
#' `do_gc`.
#'
#'
#' @param this_taxa Character. Name of taxa. Only used to print some messages.
#' Ignored if NULL
#' @param out_dir FALSE or character. If FALSE the result of prep_sdm will be
#' saved to a temporary folder. If character, a file 'prep.rds' will be created
#' at the path defined by out_dir.
#' @param return_val Character: "object" or "path". Both return a named list. In
#' the case of "path" the named list is simply list(prep = out_dir). Will be set
#' to "object" if `out_dir` is FALSE.
#' @param presence Dataframe of presences with columns `pres_x` and `pres_y`.
#' @param pres_crs Anything that will return a legitimate crs when passed to the
#' crs attribute of `sf::st_transform()` or `sf::st_as_sf()`.
#' @param pres_x,pres_y Character. Name of the columns in `presence` that have
#' the x and y coordinates
#' @param pred_limit Limit the background points and predictions?
#' Can be `TRUE` (use `presence` to generate a minimum convex polygon to use as
#' a limit. Not recommended as the points in `presence` have usually been
#' filtered to very accurate spatial reliability and thus may be missing a large
#' number of legitimate records); `FALSE` (the full extent of the predictors
#'  will be used); path to existing .parquet to use; or sf object.
#' @param limit_buffer Numeric. Apply this buffer to `pred_limit`. Only used if
#' `pred_limit` is `TRUE`. Passed to the `dist` argument of `sf::st_buffer()`.
#' @param pred_clip sf. Optional sf to clip the pred_limit back to (e.g. to
#' prevent prediction into ocean).
#' @param predictors Character. Vector of paths to predictor `.tif` files.
#' @param is_env_pred Logical. Does the naming of the directory and files in
#' `predictors` follow the pattern required by `envRaster::parse_env_tif()`?
#' @param cat_preds Character. Vector of predictor names that are character.
#' @param max_cells_in_memory Numeric passed to `exactextractr::exactextract()`
#' argument with the same name. `prep_sdm()` will be quicker with larger values,
#' but if running on many cores, memory issues are likely. The default of
#' 30000000 is the default for `exactextractr::exactextract()` at the time of
#' writing.
#' @param terra_options Passed to `terra::terraOptions()`. e.g. list(memfrac =
#' 0.6)
#' @param num_bg Numeric. How many background points?
#' @param prop_abs Character. Is `num_bg` a proportion (`prop`) of the number of
#' records in `presence` or an absolute (`abs`) number?
#' @param many_p_prop Numeric. Ensure the number of background points is at
#' least `many_p_prop * number of presences`. e.g. If there are more than 5000
#' presences and num_bg is set at `10000` and `many_p_prop` is `2`, then num_bg
#' will be increased to `many_p_prop * nrow(presences)`
#' @param folds Numeric. How many folds to use in cross validation? Will be
#' adjusted downwards if number of presences do not support `folds * min_fold_n`
#' @param spatial_folds Logical. Use spatial folds? Even if `TRUE`, can resort
#' to non-spatial cv if presences per fold do not meet `min_fold_n` or there are
#' not enough presences to support more than one fold.
#' @param hold_prop Numeric. Proportion of data to keep back for
#' validating the final model.
#' @param min_fold_n Numeric. Sets the minimum number of presences, and
#' min_fold_n * max_hold_prop sets the minimum number of presences required
#' for a model.
#' @param stretch_value Numeric. Stretch the density raster to this value.
#' @param dens_res `NULL` or numeric. Resolution (in metres) of density raster.
#' Set to `NULL` to use the same resolution as the predictors.
#' @param block_iterations Numeric. Number of attempts to achieve `min_fold_n`
#' presences in each of `folds` folds before lowering `folds` by one and trying
#' again. If attempts fail down to folds == 2, non-spatial cross validation is
#' used.
#' @param reduce_env_thresh Numeric. Threshold used to flag highly correlated
#' and low importance variables. Set to 0 to skip this step. If > 0, highly
#' correlated and low importance variables will be removed. In the case of
#' highly correlated pairs of variables, only one is removed.
#' @param do_gc Logical. Run `base::rm(list = ls)` and `base::gc()` at end of
#' function? Useful when running SDMs for many, many taxa, especially if done in
#' parallel.
#' @param force_new Logical. If outputs already exist, should they be remade?
#'
#' @return If `return_val` is "object" a named list. If `return_val` is "path"
#' a named list `list(prep = out_dir)`. If `out_dir` is a valid path, the 'full
#' result' (irrespective of `return_val`) is also saved to
#' `fs::path(out_dir, "prep.rds")`. The 'full result' is a named list with
#' elements:
#' * log:
#'     + a log of (rough) timings and other information from the process
#' * abandoned:
#'     + Logical indicating if the sdm was abandoned. If abandoned is TRUE, some
#'     list elements may not be present
#' * presence_ras:
#'     + tibble with two columns ('x' and 'y') representing unique cell
#'     centroids on the predictors at presences supplied in argument `presence`
#' * predict_boundary:
#'     + sf used to limit the background points and used by `predict_sdm()` to
#'     generate the 'mask'ed output
#' * bg_points:
#'    + sf of cell centroids representing unique cell centroids for background
#'    points
#' * blocks
#'     + data.frame with columns:
#'       + `pa`: presence (1) or absence/background (0)
#'       + `x` and `y`: cell centroids for each presence and absence
#'       + `block`: the spatial block to which the row belongs
#'       +  a column with values for each of `predictors` at `x` and `y`
#' * spatial_folds_used:
#'     + logical indicating if spatial folds were used. This may differ from
#'     the `spatial_folds` argument provided to `prep_sdm()` if an attempt to
#'     use spatial folds failed to meet desired `folds` and `min_fold_n`
#' * correlated:
#'     + list with elements as per `envModel::reduce_env()`, or, if
#'     `reduce_env` is `FALSE`, a list with elements `remove_env` which is
#'     empty, and `env_var` and `keep`, which both contain the names of all
#'     predictors.
#'
#' @export
#'
#' @example inst/examples/prep_sdm_ex.R
  prep_sdm <- function(this_taxa = NULL
                       , out_dir = FALSE
                       , return_val = "path"
                       , presence
                       , pres_crs = 4326
                       , pres_x = "long"
                       , pres_y = "lat"
                       , pred_limit = TRUE
                       , limit_buffer = 0
                       , pred_clip = NULL
                       , predictors
                       , is_env_pred = TRUE
                       , max_cells_in_memory = 3e+07
                       , terra_options = NULL
                       , cat_preds = NULL
                       , num_bg = 10000
                       , prop_abs = "abs"
                       , many_p_prop = 2
                       , folds = 5
                       , spatial_folds = TRUE
                       , hold_prop = 0.2
                       , min_fold_n = 8
                       , stretch_value = 10
                       , dens_res = 1000
                       , block_iterations = 30
                       , reduce_env_thresh = 0.9
                       , do_gc = FALSE
                       , force_new = FALSE
                       ) {

    # setup -------
    return_val <- if(any(isFALSE(out_dir), return_val == "object")) "prep" else "prep_file"

    if(isFALSE(out_dir)) {

      out_dir <- tempfile()

      delete_out <- TRUE

    } else {

      delete_out <- FALSE

    }

    log_file <- fs::path(out_dir, "prep.log")

    if(is.character(out_dir)) {

      fs::dir_create(out_dir)

      if(dir.exists(out_dir)) {

        prep_file <- fs::path(out_dir
                              , "prep.rds"
                              )

        if(file.exists(prep_file)) {

          prep <- rio::import(prep_file)

        }

      } else stop("can't create out_dir")

    }

    if(!exists("prep", inherits = FALSE)) prep <- list(abandoned = FALSE, finished = FALSE)

    if(nrow(presence) == 0) prep$abandoned <- TRUE

    run <- if(any(prep$abandoned, prep$finished)) force_new else TRUE

    if(run) {

      if(is.null(this_taxa)) this_taxa <- basename(out_dir)

      prep$this_taxa <- this_taxa

      message(paste0("prep for "
                     , this_taxa
                     , "\nout_dir is "
                     , out_dir
                     )
              )

      ## start timer -----
      start_time <- Sys.time()

      readr::write_lines(paste0("\n"
                                , this_taxa
                                , "\nprep start at "
                                , start_time
                                , ".\n"
                                , nrow(presence), " original record"
                                , if(nrow(presence) > 1) "s"
                                )
                         , log_file
                         , append = TRUE
                         )

      ## predictors -----
      if(is_env_pred) {

        pred_names <- envRaster::name_env_tif(tibble::tibble(path = predictors), parse = TRUE) %>%
          dplyr::pull(name)

        prep_preds <- terra::rast(predictors)
        names(prep_preds) <- pred_names

      } else {

        prep_preds <- terra::rast(predictors)

      }

      readr::write_lines(paste0(length(names(prep_preds))
                                , " predictors"
                                )
                         , file = log_file
                         , append = TRUE
                         )

      # presence --------
      prep$original <- presence

      ## raster presence ------

      prep$presence_ras <- terra::cellFromXY(prep_preds
                                             , presence %>%
                                               sf::st_as_sf(coords = c(pres_x, pres_y)
                                                            , crs = pres_crs
                                                            ) %>%
                                               sf::st_transform(crs = sf::st_crs(prep_preds)) %>%
                                               sf::st_coordinates()
                                             ) %>%
        terra::xyFromCell(prep_preds, .) %>%
        stats::na.omit() %>%
        tibble::as_tibble()


      # predict_boundary -------

      ## TRUE ------
      # create new mcp from points
      if(isTRUE(pred_limit)) {

        prep$predict_boundary <- presence %>%
          dplyr::distinct(!!rlang::ensym(pres_x), !!rlang::ensym(pres_y)) %>%
          sf::st_as_sf(coords = c(pres_x, pres_y)
                       , crs = pres_crs
                       ) %>%
          sf::st_transform(crs = sf::st_crs(prep_preds[[1]])) %>%
          sf::st_union() %>%
          sf::st_convex_hull() %>%
          sf::st_sf() %>%
          sf::st_buffer(limit_buffer) %>%
          sf::st_make_valid()

      }

      ## file path -------
      # use existing MCP polygon file for the predict boundary
      if(is.character(pred_limit)) {

        if(file.exists(pred_limit)) {

          prep$predict_boundary <- sfarrow::st_read_parquet(pred_limit) %>%
            sf::st_geometry() %>%
            sf::st_sf() %>%
            sf::st_transform(crs = sf::st_crs(prep_preds[[1]])) %>%
            sf::st_make_valid()

        }

      }

      ## sf --------
      # existing sf
      if("sf" %in% class(pred_limit)) {

        prep$predict_boundary <- pred_limit %>%
          sf::st_transform(crs = sf::st_crs(prep_preds[[1]])) %>%
          sf::st_make_valid() %>%
          sf::st_intersection(sf::st_bbox(prep_preds) %>%
                                sf::st_as_sfc() %>%
                                sf::st_sf()
                              ) %>%
          sf::st_make_valid()

      }

      # recast pred_limit from here on as T/F
      if(!isFALSE(pred_limit)) pred_limit <- TRUE

      ## FALSE ------
      # use the full extent of the predictors for the predict boundary
      # No point buffering here
      # captures anything not already captured above
      if(isFALSE(pred_limit)) {

        prep$predict_boundary <- sf::st_bbox(prep_preds) %>%
          sf::st_as_sfc() %>%
          sf::st_sf()  %>%
          sf::st_make_valid()

      }

      # clip predict_boundary? -------
      if(!is.null(pred_clip)) {

        prep$predict_boundary <- prep$predict_boundary %>%
          sf::st_intersection(pred_clip %>%
                                sf::st_transform(crs = sf::st_crs(prep$predict_boundary))
                              ) %>%
          sf::st_make_valid()

      }

      # predict_boundary only on env rasters
      prep$predict_boundary <- prep$predict_boundary %>%
        sf::st_intersection(sf::st_bbox(prep_preds) %>%
                              sf::st_as_sfc() %>%
                              sf::st_sf()
                            )


      # prep$presence_ras clip to predict_boundary ---------
      # catch some cases where there are presence records on the predictors
      # but outside the predict boundary.
      prep$presence_ras <- prep$presence_ras %>%
        sf::st_as_sf(coords = c("x", "y")
                     , crs = sf::st_crs(prep_preds)
                     ) %>%
        sf::st_filter(prep$predict_boundary %>%
                        sf::st_transform(crs = sf::st_crs(prep_preds))
                      ) %>%
        sf::st_coordinates() %>%
        tibble::as_tibble() %>%
        purrr::set_names(c("x", "y"))

      n_p <- nrow(prep$presence_ras)

      readr::write_lines(paste0(n_p
                                , " presences in predict_boundary and on env rasters"
                                )
                         , file = log_file
                         , append = TRUE
                         )

      # Adjust args ---------
      ## folds adj -------

      k_folds <- folds

      while(n_p < min_fold_n * k_folds) {

        k_folds <- k_folds - 1

      }

      if(k_folds < folds) {

        readr::write_lines(paste0("warning: too few records to support original folds ("
                                  , folds
                                  , "). Folds reduced to "
                                  , k_folds
                                  )
                           , file = log_file
                           , append = TRUE
                           )

      }

      if(k_folds < 3) {

        readr::write_lines(paste0("WARNING: there are only "
                                  , n_p
                                  , " presences. This is not enough to support the bare minimum"
                                  , " of three folds (one as 'test' set and two as 'training' sets"
                                  , " with "
                                  , min_fold_n
                                  , " in each fold. SDM abandoned"
                                  )
                           , file = log_file
                           , append = TRUE
                           )

        prep$abandoned <- TRUE

      }

      if(!prep$abandoned) {

        # num_bg adj------

        if(prop_abs == "prop") {

          num_bg <- num_bg * n_p

        }

        if(num_bg < many_p_prop * n_p) num_bg <- many_p_prop * n_p

        if(n_p > min_fold_n) {

          # subset predictors? ------
          if(pred_limit) {

            start_subset <- Sys.time()

            if(!is.null(terra_options)) {

              do.call(terra::terraOptions
                      , args = terra_options
                      )

            }

            terra::window(prep_preds) <- terra::ext(terra::vect(prep$predict_boundary))

            readr::write_lines(paste0("subsetting predictors done in "
                                      , round(difftime(Sys.time(), start_subset, units = "mins"), 2)
                                      , " minutes"
                                      )
                               , file = log_file
                               , append = TRUE
                               )

          }

          # density raster ------

          density_file <- fs::path(out_dir
                                   , "density.tif"
                                   )

          run <- if(file.exists(density_file)) force_new else TRUE

          if(run) {

            start_dens_ras <- Sys.time()

            if(all(!is.null(dens_res), !terra::is.lonlat(prep_preds[[1]]))) {

              # resolution of density raster < pred raster

              use_res <- if(terra::res(prep_preds)[[1]] < dens_res) dens_res else terra::res(prep_preds)[[1]]

              temp_ras <- terra::rast(resolution = use_res
                                      , crs = terra::crs(prep_preds[[1]])
                                      , extent = terra::ext(prep$predict_boundary)
                                      , vals = 1
                                      )

              bw <- MASS::kde2d(as.matrix(prep$presence_ras[,1])
                                , as.matrix(prep$presence_ras[,2])
                                , n = c(nrow(temp_ras), ncol(temp_ras))
                                , lims = terra::ext(temp_ras) %>% as.vector()
                                )

              target_density <- raster::raster(bw) %>%
                terra::rast() %>%
                terra::resample(temp_ras) %>%
                terra::focal(3
                             , mean
                             , na.policy = "only"
                             , na.rm = TRUE
                             )

            } else {

              # density raster = predict raster
              # STU METHOD
              # weighted bandwidth
              ## same method as spatialEco::sf.kde

              temp_ras <- terra::rast(resolution = terra::res(prep_preds)
                                      , crs = terra::crs(prep_preds[[1]])
                                      , extent = terra::ext(prep$predict_boundary)
                                      , vals = 1
                                      )

              pres_ras <- terra::rasterize(prep$presence_ras %>%
                                             sf::st_as_sf(coords = c("x", "y")
                                                          , crs = sf::st_crs(prep_preds[[1]])
                                                          )
                                           , y = temp_ras
                                           , fun = length
                                           , touches = TRUE
                                           )

              pres <- terra::as.data.frame(pres_ras, xy = TRUE) %>%
                tibble::as_tibble()

              colnames(pres)[3] <- "COUNT"

              pres <- pres %>%
                dplyr::mutate(scaled = COUNT * (length(COUNT) / sum(COUNT)))

              if(all(pres$scaled == 1)) pres$scaled <- sample(c(0.999, 1.001), length(pres$scaled), replace = TRUE) # nw

              bw <- ks::Hpi.diag(x = as.matrix(pres[, c("x", "y", "scaled")])
                                 , pilot = "dscalar"
                                 )

              target_density <- terra::setValues(pres_ras
                                                 , matrix(ks::kde(as.matrix(pres[, c("x", "y")])
                                                                  , eval.points = terra::xyFromCell(pres_ras, 1:terra::ncell(pres_ras))
                                                                  , gridsize = c(nrow(pres_ras), ncol(pres_ras))
                                                                  # use bw or NULL if error
                                                                  ## NULL will use default method - maybe unweighted?
                                                                  , h = ifelse(inherits(bw, "error"), NULL, bw)
                                                                  , w = pres[["scaled"]]
                                                                  , density = TRUE
                                                                  )$estimate
                                                          , nrow = nrow(pres_ras)
                                                          , ncol = ncol(pres_ras)
                                                          , byrow = TRUE
                                                          )
                                                 ) %>%
                terra::focal(3
                             , mean
                             , na.policy = "only"
                             , na.rm = TRUE
                             )

            }

            if(pred_limit) {

              target_density <- target_density %>%
                terra::mask(terra::vect(prep$predict_boundary)
                              , touches = TRUE
                            )

            }

            ## stretch-------
            target_density <- target_density %>%
              terra::stretch(1, stretch_value) %>%
              terra::subst(NA,0)

            if(FALSE) {

              terra::plot(prep_preds[[1]])
              terra::plot(target_density, add = TRUE)

              ps <- prep$presence_ras %>%
                tibble::as_tibble() %>%
                sf::st_as_sf(coords = c("x", "y"), crs = sf::st_crs(prep_preds)) %>%
                terra::vect()

              terra::plot(ps, add = TRUE)

              terra::plot(prep$predict_boundary %>% terra::vect()
                          , add = TRUE
                          )

            }

            terra::writeRaster(target_density
                               , density_file
                               , overwrite = TRUE
                               )

            readr::write_lines(paste0("density raster done in "
                                      , round(difftime(Sys.time(), start_dens_ras, units = "mins"), 2)
                                      , " minutes"
                                      )
                               , file = log_file
                               , append = TRUE
                               )

          }


          # background points -------
          run <- if(exists("bg_points", where = prep)) force_new else TRUE

          if(run) {

            start_bg <- Sys.time()

            if(!exists("target_density")) target_density <- terra::rast(density_file)

            ptscell <- sample(1:terra::ncell(target_density)
                             , num_bg * 2 # over generate here as terra::window only crops - it does not mask
                             , prob = target_density[]
                             , replace = TRUE
                             )

            centres <- terra::xyFromCell(target_density, ptscell)

            hs <- terra::res(target_density) / 2

            prep$bg_points <- cbind(runif(nrow(centres), centres[, 1] - hs[1], centres[, 1] + hs[1])
                                    , runif(nrow(centres), centres[, 2] - hs[2], centres[, 2] + hs[2])
                                    ) %>%
              as.matrix() %>%
              terra::cellFromXY(prep_preds[[1]], .) %>%
              terra::xyFromCell(prep_preds[[1]], .) %>%
              tibble::as_tibble() %>%
              na.omit() %>%
              dplyr::distinct() %>%
              dplyr::anti_join(prep$presence_ras) %>% # background not on presences
              sf::st_as_sf(coords = c("x", "y")
                           , crs = sf::st_crs(target_density)
                           ) %>%
              # filter back to only those within predict_boundary
              sf::st_filter(prep$predict_boundary)

            if(nrow(prep$bg_points) > num_bg) {

              prep$bg_points <- prep$bg_points %>%
                dplyr::sample_n(num_bg)

            }

            readr::write_lines(paste0(num_bg
                                      , " background points attempted. "
                                      , nrow(prep$bg_points)
                                      , " achieved in: "
                                      , round(difftime(Sys.time(), start_bg, units = "mins"), 2)
                                      , " minutes"
                                      )
                               , file = log_file
                               , append = TRUE
                               )


          } else {

          target_density <- terra::rast(density_file)

        }


        # env--------

        run <- if(exists("env", prep)) force_new else TRUE

        if(run) {

          start_env <- Sys.time()

          spp_pa <- dplyr::bind_rows(prep$presence_ras %>%
                                       tibble::as_tibble() %>%
                                       sf::st_as_sf(coords = c("x", "y")
                                                    , crs = sf::st_crs(prep_preds[[1]])
                                                    ) %>%
                                       dplyr::mutate(pa = "presence")
                                     , prep$bg_points %>%
                                       dplyr::mutate(pa = "absence") %>%
                                       dplyr::select(pa)
                                     ) %>%
            sf::st_buffer(terra::res(prep_preds)[[1]] / 100)

          prep$env <- exactextractr::exact_extract(prep_preds
                                                     , y = spp_pa
                                                     , include_cols = "pa"
                                                     , include_cell = TRUE
                                                     , max_cells_in_memory = max_cells_in_memory
                                                     ) %>%
            dplyr::bind_rows() %>%
            stats::na.omit() %>%
            dplyr::select(-coverage_fraction) %>%
            dplyr::bind_cols(terra::xyFromCell(prep_preds
                                               , .$cell
                                               )
                             ) %>%
            dplyr::mutate(pa = forcats::fct_relevel(pa, "presence"))

          if(!is.null(cat_preds)) {

            prep$env <- prep$env %>%
              dplyr::mutate(dplyr::across(tidyselect::any_of(cat_preds), as.factor))

          }

          readr::write_lines(paste0("env data extracted in: "
                                    , round(difftime(Sys.time(), start_bg, units = "mins"), 2)
                                    , " minutes"
                                    )
                             , file = log_file
                             , append = TRUE
                             )

        }

        # Abandon if too few presences with env data
        if(nrow(prep$env[prep$env$pa == "presence",]) < min_fold_n) {

          readr::write_lines(paste0("warning: too few presences ("
                                    , nrow(prep$env[prep$env$pa == 1,])
                                    , ") with environmental variables. SDM abandoned"
                                    )
                             , file = log_file
                             , append = TRUE
                             )

          prep$abandoned <- TRUE

        }


        # reduce env -------

        if(all(reduce_env_thresh, !prep$abandoned)) {

          run <- if(exists("reduce_env", prep)) force_new else TRUE

          if(run) {

            start_reduce_env <- Sys.time()

            prep$reduce_env <- envModel::reduce_env(env_df = prep$env
                                                    , env_cols = names(prep_preds)
                                                    , y_col = "pa"
                                                    , imp_col = "presence"
                                                    , thresh = reduce_env_thresh
                                                    , remove_always = c(pres_x, pres_y, "x", "y", "pa", "block", "cell")
                                                    )

            readr::write_lines(paste0("reduce_env completed in: "
                                      , round(difftime(Sys.time(), start_reduce_env, units = "mins"), 2)
                                      , " minutes"
                                      )
                               , file = log_file
                               , append = TRUE
                               )

          }

        } else {

          if(!prep$abandoned) {

            prep$reduce_env$remove <- ""
            prep$reduce_env$keep <- names(prep_preds)
            prep$reduce_env$env_cols <- names(prep_preds)

          }

        }

          # folds------

        if(!prep$abandoned) {

          run <- if(exists("blocks", prep)) force_new else TRUE

          if(run) {

            start_blocks <- Sys.time()

            ## split -------
            prep$split_check <- list(presence_test = 0)

            while(prep$split_check$presence_test <= 0) {

              prep$split <- rsample::initial_split(prep$env %>%
                                                     dplyr::select(x, y, pa, prep$reduce_env$keep) %>%
                                                     dplyr::mutate(pa = factor(pa))
                                                   , strata = "pa"
                                                   , prop = hold_prop
                                                   )

              prep$split_check <- tidysdm::check_splits_balance(prep$split
                                                                , pa
                                                                )

            }

            ## folds -------
            prep$blocks_check <- list(presence_assessment = 0)

            counter <- 0

            block_folds <- k_folds

            while(all(min(prep$blocks_check$presence_assessment) <= min_fold_n, block_folds > 1)) {

              while(all(counter < block_iterations, min(prep$blocks_check$presence_assessment) <= min_fold_n)) {

                prep$blocks <- spatialsample::spatial_block_cv(data = rsample::testing(prep$split) %>%
                                                                 sf::st_as_sf(coords = c("x", "y")
                                                                              , crs = sf::st_crs(prep_preds)
                                                                              )
                                                               , v = block_folds
                                                               )

                prep$blocks_check <- tidysdm::check_splits_balance(prep$blocks
                                                                     , pa
                                                                     )

                counter <- counter + 1

              }

              if(min(prep$blocks_check$presence_assessment) <= min_fold_n) {

                prep$blocks_check <- list(presence_assessment = 0)

                block_folds <- block_folds - 1

                counter <- 0

              }

            }

            if(min(prep$blocks_check$presence_assessment) <= min_fold_n) {

              prep$spatial_folds <- FALSE

              prep$blocks <- rsample::vfold_cv(data = rsample::testing(prep$split)
                                               , v = k_folds
                                               , strata = pa
                                               )


            }

            readr::write_lines(paste0("validation split and blocks completed in: "
                                      , round(difftime(Sys.time(), start_blocks, units = "mins"), 2)
                                      , " minutes"
                                      )
                               , file = log_file
                               , append = TRUE
                               )

          }

        }

          # end timer ------
          readr::write_lines(paste0("prep completed. elapsed time: "
                                    , round(difftime(Sys.time(), start_time, units = "mins"), 2)
                                    , " minutes"
                                    )
                             , file = log_file
                             , append = TRUE
                             )

      }

    }

    # save / clean up-------
    # export before gc()
    prep$finished <- TRUE
    prep$log <- if(file.exists(log_file)) readr::read_lines(log_file) else NULL

    if(delete_out) {

      fs::dir_delete(out_dir)

    } else {

      rio::export(prep, prep_file)

    }

    # clean up -------
    if(do_gc) {

      stuff <- ls()

      stuff <- stuff[! stuff %in% c(return_val, "return_val")]

      rm(list = stuff)

      gc()

    }

    res <- if(return_val == "prep") get("prep") else list(prep_file = get("prep_file"))

    return(res)

    }

  }

