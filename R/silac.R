#' Determine SILAC (2-plex) incorporation rate given intensity values for Light and Heavy
#'
#' @description SILAC incorporation is estimated by H / (H + L), where H = Heavy
#' and L = Light labelled peptide intensity. If L is `NA`, incorporation == 1
#' and if H is `NA`, incorporation == 0.
#'
#' @param light `numeric`. Light peptide intensity.
#' @param heavy `numeric`. Heavy peptide intensity.
#' @return Returns a `numeric` from 0 to 1 (or `NA` if `light` and `heavy` are both `NA`)
#' @export
get_incorporation <- function(light, heavy) {
  if (is.na(light) | is.na(heavy)) {
    if (is.na(light) & is.na(heavy)) {
      incorporation = NA
    }
    else if (is.na(heavy)) {
      incorporation = 0
    }
    else {
      incorporation = 1
    }
  }
  else {
    incorporation = heavy / (light + heavy)
  }
  return(incorporation)
}


#' Plot annotated histogram of incorporation values
#'
#' @description SILAC (2-plex) label incorporation is estimated from the set of
#' observed features (peptides or proteins). If the test sample contains cells
#' grown in 'Heavy' media only, incorporation is simply the mean observed
#' incorporation. If the test sample contains a mixture of
#' 'Heavy' and 'Light', the incorporation is best calculated from the median
#' incorporation using the corrected 'Light' intensity.
#'
#' @param data `data.frame` containing features with column 'incorporation'.
#' If mix>0, must also contain column 'corrected_incorporation' where incorporation
#' is calculated using the corrected light peptide intensity
#' @param incorporation_col `string`. Name of the column containing incorporation values
#' which can be calculated with \code{\link{get_incorporation}}.
#' @param corrected_col `string` Optional, but required if mix > 0. Name of the
#' column containing mix-corrected incorporation values which can be calculated
#' with \code{\link{get_incorporation}}.
#' @param level `string`. Name for feature level, either 'peptide' or 'protein'.
#' @param mix `numeric`. The 'Heavy' to 'Light' mix ratio. Default is 0, otherwise
#' use a number greater than 0.
#' @return Returns a `list` with `p` = `ggplot` plot and
#' `incorporation_estimates` = `list` of incorporation estimates.
#' @export
plot_incorporation <- function(
  data,
  incorporation_col = "Incorporation",
  corrected_col = "Incorporation.corrected",
  level = "peptide",
  mix = 0
) {
  # basic input checking
  match.arg(level, c("peptide", "protein"))
  stopifnot(mix >= 0)

  p <- data %>%
    ggplot() +
    geom_histogram(aes(x = 100 * !!sym(incorporation_col)), fill = "grey") +
    theme_csd(base_size = 10) +
    xlab('Observed incorporation (%)') +
    ylab('Count')

  total_count <- nrow(data)

  median_incorporation <- stats::median(data[, incorporation_col], na.rm = TRUE) * 100

  if (mix > 0) {
    median_corrected_incorporation <- stats::median(data[, corrected_col], na.rm = TRUE) * 100

    incorporation_estimates <- list(
      'Median incorporation' = median_incorporation,
      'Median incorporation (corrected)' = median_corrected_incorporation
    )

    plot_text <- paste0(
      '%ss: %s\n',
      'Median incorporation: %#.2f %%\n',
      'Corrected median incorporation:%#.2f %%'
    )

    sprintf_values <- list(
      plot_text, tools::toTitleCase(level), total_count,
      median_incorporation,
      median_corrected_incorporation
    )
  } else {
    mean_incorporation <- mean(data[, incorporation_col], na.rm = TRUE) * 100

    incorporation_estimates <- list(
      'Median incorporation' = median_incorporation,
      'Mean incorporation' = mean_incorporation
    )

    plot_text <- paste0(
      '%ss: %s\n',
      'Median incorporation: %.2f %%\n',
      'Mean incorporation:%.2f %%'
    )
    sprintf_values <- list(
      plot_text, tools::toTitleCase(level), total_count,
      median_incorporation, mean_incorporation
    )
  }

  plot_label <- do.call('sprintf', sprintf_values)
  p <- p + annotate(
    geom = 'text', x = 0, y = Inf, label = plot_label, vjust = 1.5, hjust = 0,
    size = 3, colour = get_cat_palette(1)
  )

  invisible(list('p' = p, 'incorporation_estimates' = incorporation_estimates))
}
