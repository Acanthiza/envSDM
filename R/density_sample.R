
#' Create a spatially thickened background sample
#'
#' Same idea as `terra::spatSample(weights = weight_raster)`, but can use a
#' density raster with a different spatial resolution. Can be much quicker than
#' `terra::spatSample()` for large areas with small resolutions.
#'
#' @param x spatRaster. Sample will be cell XY values from this raster.
#' @param dens_rast spatRaster. Density raster. This is equivalent to the
#' `weights` argument of `terra::spatSample()` but does not have to be at the
#' same resolution as `x`. Usually made within `prep_sdm()`
#' @param samp_df Dataframe. Usually built by `prep_sdm()`
#' @param boundary sf. Boundary within which to sample cells from `x`
#' @param mult Numeric. In each iteration try sampling `mult` *
#' `sum(samp_df$target)` points
#' @param max_sample Numeric. How many iterations to try to meet the density
#' targets in `samp_df`.
#' @param verbose Logical. Prints how many iterations have been attempted.
#'
#' @return Dataframe
#' @export
#' @keywords internal
#'
  density_sample <- function(x
                             , dens_rast
                             , samp_df
                             , boundary
                             , mult = 2
                             , max_sample = 200
                             , verbose = FALSE
                             ) {

    cell_ids <- terra::cells(x, terra::vect(boundary))[,2]

    a_sample <- function(x, dens_rast, mult, samp_df) {

      sample(cell_ids, sum(samp_df$target) * mult) %>%
        terra::xyFromCell(x, .) %>%
        as.data.frame %>%
        cbind(terra::extract(dens_rast, as.matrix(.))) %>%
        stats::na.omit() %>%
        dplyr::rename(value = 3) %>%
        merge(samp_df[, c("value", "target")])

    }

    new <- a_sample(x, dens_rast, mult, samp_df)

    todo <- TRUE
    counter <- 1

    while(all(counter < max_sample, todo)) {

      new <- rbind(new[,1:4], a_sample(x, dens_rast, mult, samp_df)) %>%
        dplyr::add_count(value
                         , name = "n"
                         ) %>%
        dplyr::group_by(value) %>%
        dplyr::slice(1:unique(target)) %>%
        dplyr::ungroup() %>%
        dplyr::distinct() %>%
        dplyr::add_count(value, name = "n")

      todo <- new %>%
        dplyr::count(value, target) %>%
        dplyr::mutate(todo = n < target) %>%
        dplyr::pull(todo) %>%
        sum()

      counter <- counter + 1

      if(verbose) message("counter: ", counter)

    }

    stuff <- ls() %>%
      grep("res", value = TRUE, invert = TRUE)

    rm(list = stuff)

    gc()

    return(new)

  }
