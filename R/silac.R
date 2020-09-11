
#' Determine the SILAC incorporation rate given intensity values for light and heavy
#'
#' @description SILAC incorporation is estimated by H / (H + L), where H=Heavy
#' and L=Light labelled peptide intensity. If L is NA, incorporation==1, if H is
#' NA, incorproration==0
#'
#' @param light `numeric` Light peptide intensity
#' @param heavy `numeric` Heavy peptide intensity
#' @return `incorporation`
#' @export
get_incorporation <- function(light, heavy){

  if (is.na(light)|is.na(heavy)){

    if (is.na(light) & is.na(heavy)){
      incorporation = NA
    }
    else if (is.na(heavy)){
      incorporation = 0
    }
    else {
      incorporation = 1
    }
  }
  else{
    incorporation = heavy/(light+heavy)
  }
  return(incorporation)
}


#' Plot annotated histogram of incorporation values
#'
#' @description Incorporation is estimated from the set of observed features. If
#' the test sample contains cells grown in 'heavy' media only, incorporation is
#' simply the mean observed incorporation. If the test sample contains a mixture of
#' 'Heavy' and 'Light', the incorporation is best calculated from the median
#' incorporation using the corrected 'Light' intensity.
#' From a `data.frame`, with
#'
#' @param obj `data.frame` containing features with column 'incorporation'.
#' If mix>0, must also contain column 'corrected_incorporation' where incorporation
#' is calculated using the corrected light peptide intensity
#' @param level `string`. Name for feature level. e.g 'Peptide' or 'Protein'
#' @param mix `numeric`. Default is 0. Otherwise use a number greater than zero
#' (i.e. the 'Heavy' to 'Light' mix ratio).
#' @return `list` with `p`=`ggplot` plot and `incorporation_estimates`=`list` of
#' incorporation estimates
#' @export
plot_incorporation <- function(obj, level='Peptide', mix=0){

  p <- obj %>%
    ggplot() +
    geom_histogram(aes(100*.data$incorporation), fill='grey') +
    theme_camprot(base_size=10) +
    xlab('Observed incorporation (%)') +
    ylab('Count')

  total_count <- nrow(obj)

  if(mix>0){

    median_incorporation <- obj %>%
      pull(.data$incorporation) %>%
      stats::median(na.rm=TRUE) %>%
      "*"(100)

    median_corrected_incorporation <- obj %>%
      pull(.data$corrected_incorporation) %>%
      stats::median(na.rm=TRUE) %>%
      "*"(100)

    incorporation_estimates <- list(
      'Median incorporation' = median_incorporation,
      'Median incorporation (corrected)' = median_corrected_incorporation)
    text <- paste0('%ss: %s\n',
                   'Median incorporation: %#.2f %%\n',
                   'Corrected median incorporation:%#.2f %%')

    sprintf_values <- list(text, level, total_count,
                           median_incorporation,
                           median_corrected_incorporation)

  } else{
    median_incorporation <- obj %>%
      pull(.data$incorporation) %>%
      stats::median(na.rm=TRUE) %>%
      "*"(100)

    mean_incorporation <- obj %>%
      pull(.data$incorporation) %>%
      mean(na.rm=TRUE) %>%
      "*"(100)


    incorporation_estimates <- list(
      'Median incorporation' = median_incorporation,
      'Mean incorporation' = mean_incorporation)

    text <- paste0('%ss: %s\n',
                   'Median incorporation: %.2f %%\n',
                   'Mean incorporation:%.2f %%')
    sprintf_values <- list(text, level, total_count,
                           median_incorporation, mean_incorporation)
  }

  label <- do.call('sprintf', sprintf_values)
  p <- p + annotate(geom='text', x=0, y = Inf, label=label, vjust=1.5, hjust=0,
                    size=3, colour=get_cat_palette(1))

  invisible(list('p'=p, 'incorporation_estimates'=incorporation_estimates))
}

