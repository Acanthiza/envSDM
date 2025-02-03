
#' Tune, and evaluate, species distribution models
#'
#' @param prep Character or named list. If character, the path to an existing
#' `prep.rds`. Otherwise, the result of a call to prep_sdm with return_val =
#' "object"
#' @param out_dir FALSE or character. If FALSE the result of tune_sdm will be
#' saved to a temporary folder. If character, a file 'tune.rds' will be created
#' at the path defined by out_dir.
#' @param return_val Character: "object" or "path". Both return a named list. In
#' the case of "path" the named list is simply list(tune = out_dir). Will be set
#' to "object" if `out_dir` is FALSE.
#' @param algo Character. Name of algorithm to use.
#' @param fc Character. Used to generate levels of `classes` argument to
#' `maxnet::maxnet()` that are tuned.
#' @param limit_p `TRUE`, `FALSE` or number of predictor variables above which
#' to limit the use of `p` in the classes argument used in `maxnet::maxnet()`.
#' Useful with many predictor variables when it becomes unwieldy to generate
#' interactions for all predictors.
#' @param rm Numeric. Used to generate levels of `regmult` argument to
#' `maxnet::maxnet()` that are tuned.
#' @param trees Used to generate the levels of `ntree` argument to
#' `randomForest::randomForest()` that are tuned. `TRUE` (tune with default
#' `trees`), `FALSE` (don't tune `trees`) or numeric (the `trees` values to tune
#'  with).
#' @param mtry Used to generate the levels of `mtry` argument to
#' `randomForest::randomForest()` that are tuned. `TRUE` (tune with sensible guesses for
#' `mtry`), `FALSE` (only use default `randomForest::randomForest()` `mtry`) or
#' numeric (the `mtry` values to tune with).
#' @param limit_spat_mtry Numeric. If `mtry` is `TRUE` and if using spatial
#' cross validation, the values of `mtry` to tune will be limited to less than
#' or equal to `limit_spat_mtry`.
#' @param nodesize Used to generate the levels of `nodesize` argument to
#' `randomForest::randomForest()` that are tuned. `TRUE` (tune with default
#' `nodesize`), `FALSE` (only use default `randomForest::randomForest()`
#' `nodesize`) or numeric (the `nodesize` values to tune with).
#' @param keep_model Logical. If `TRUE` the model results will be appended as a
#' list column in the returned tibble (as column `m`)
#' @param best_run Logical. If `TRUE` this alters the behaviour of the
#' `tune_sdm()` by, well, not tuning. :). Sets all blocks to the same value so
#' no cross-validation.
#' @param metrics_df Dataframe. Defines which metrics to use when deciding on
#' 'good' SDMs.
#' @param use_metrics Character. Vector of values in metrics_df$metric to use
#' when finding the 'best' model.
#' @param do_gc Logical. Run `base::rm(list = ls)` and `base::gc()` at end of
#' function? Useful when running SDMs for many, many taxa, especially if done in
#' parallel.
#' @param force_new Logical. If outputs already exist, should they be remade?
#' @param ... Passed to `evaluate_sdm()`. e.g. thresholds for use in
#' `predicts::pa_evaluate()` (as `tr` argument, although if used, the values of
#' the `thresholds` element of the `pa_ModelEvaluation` object returned by
#' `predicts::pa_evaluate()` will be limited to the values in `tr`).
#'
#' @return If `return_val` is "object" a named list. If `return_val` is "path"
#' a named list `list(prep = out_dir)`. If `out_dir` is a valid path, the 'full
#' result' (irrespective of `return_val`) is also saved to
#' `fs::path(out_dir, "prep.rds")`. The 'full result' is a named list with
#' elements:
#'
#' @export
#'
#' @example inst/examples/tune_sdm_ex.R
  tune_sdm <- function(prep
                       , out_dir = FALSE
                       , return_val = "path"
                       , algo = c("all", "maxnet", "bioclim", "envelope", "rf") # envelope = bioclim
                       , keep_model = FALSE
                       , metrics_df = envSDM::sdm_metrics
                       , use_metrics = c("auc_po", "CBI_rescale", "IMAE")
                       , do_gc = TRUE
                       , force_new = FALSE
                       , ...
                       ) {

    # setup -------
    ## return ------
    return_val <- if(any(isFALSE(out_dir), return_val == "object")) "tune" else "tune_file"

    if(isFALSE(out_dir)) out_dir <- tempfile()

    ## log file ------
    log_file <- fs::path(out_dir
                         , "tune.log"
                         )

    ## out_dir ------
    if(is.character(out_dir)) {

      fs::dir_create(out_dir)

      if(dir.exists(out_dir)) {

        tune_file <- fs::path(out_dir
                              , "tune.rds"
                              )

        if(file.exists(tune_file)) {

          tune <- rio::import(tune_file)

        }

      } else stop("can't create out_dir")

    }

    ## prep -------
    if(! "list" %in% class(prep)) prep <- rio::import(prep)

    ## tune ---------
    if(!exists("tune", inherits = FALSE)) tune <- list(finished = FALSE)

    # run?-----
    run <- all(!prep$abandoned
               , prep$finished
               , if(tune$finished) force_new else TRUE
               )

    if(run) {

      this_taxa <- prep$this_taxa

      message("tuning "
              , this_taxa
              , "\nout_dir is "
              , out_dir
              )

      # tune -----
      if(exists("blocks", where = prep)) {

        # start timer ------
        start_time <- Sys.time()

        readr::write_lines(paste0("\n\n"
                                  , this_taxa
                                  # , "\nbest_run = "
                                  # , best_run
                                  , "\ntune start "
                                  , start_time
                                  )
                           , file = log_file
                           , append = TRUE
                           )

        ## recipe ------
        tune$rec <- recipes::recipe(x = prep$env
                                    , vars = c("pa"
                                               , prep$reduce_env$keep
                                               )
                                    , roles = c("outcome"
                                                , rep("predictor", length(prep$reduce_env$keep))
                                                )
                                    )

        ## workflow --------
        tune$wf <-
          # create the workflow_set
          workflowsets::workflow_set(
            preproc = list(default = tune$rec)
            , models = list(
              # rf specs with tuning
              rf = tidysdm::sdm_spec_rf(mode = "classification"
                                        , trees = 3001
                                        ) %>%
                parsnip::set_engine(engine = "randomForest"
                                    , sampsize = rep(min(prep$blocks_check$presence_assessment)
                                                     , 2
                                                     )
                                    )
              # maxent specs with tuning
              , maxent = tidysdm::sdm_spec_maxent()
              )
            , # make all combinations of preproc and models
            cross = TRUE
            ) %>%
          # tweak controls to store information needed later to create the ensemble
          workflowsets::option_add(control = tidysdm::control_ensemble_grid()) %>%
          workflowsets::workflow_map("tune_grid"
                                     , resamples = prep$blocks
                                     , grid = 10
                                     , metrics = tidysdm::sdm_metric_set()
                                     , verbose = TRUE
                                     )

        tune$mod <- tidysdm::simple_ensemble() %>%
          tidysdm::add_member(tune$wf, metric = "boyce_cont")

        terra::window(prep_preds) <- terra::ext(prep$predict_boundary)

       tune$pred <- tidysdm::predict_raster(tune$mod, prep_preds) %>%
         terra::mask(mask = prep$predict_boundary)

       tune$eval <- predicts::pa_evaluate(p = predict(tune$mod
                                                      , new_data = prep$split %>% training() %>% dplyr::filter(pa == "presence")
                                                      )[,1]
                                          , a = predict(tune$mod
                                                        , new_data = prep$split %>% training() %>% dplyr::filter(pa == "absence")
                                                        )[,1]
                                          )

       tune$thresh <- tune$eval@thresholds$max_spec_sens

       tune$pred_thresh <- tune$pred > tune$thresh

    }

    # save -------
    # export before gc()
    tune$finished <- TRUE
    tune$log <- if(file.exists(log_file)) readr::read_lines(log_file) else NULL

    rio::export(tune, tune_file)

    # clean up --------

    if(do_gc) {

      stuff <- ls()

      delete_stuff <- stuff[! stuff %in% c(return_val, "return_val")]

      rm(list = delete_stuff)

      gc()

    }

    res <- if(return_val == "tune") get("tune") else list(tune_file = get("tune_file"))

    return(res)

  }

}
