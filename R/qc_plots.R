#' Plot distributions for feature intensities per sample.
#'
#' @description Given an `MSnSet`, return a plot
#' of the feature intensities per sample.
#'
#' @param obj `MSnSet`.
#' @param method `string`. Plot style. Choice of box, density or histogram plot.
#' @param facet_by_sample `logical`. Facet the plot by sample.
#'
#' @return `ggplot` object.
#' @export
plot_quant <- function(obj,
                       method=c('box', 'density', 'histogram'),
                       facet_by_sample=FALSE){
  # check method argument
  method = match.arg(method)

  e_data <- data.frame(exprs(obj))

  e_data[e_data==""] <- NA
  e_data <- e_data %>% gather(key='sample', value='intensity') %>%
    mutate(sample = factor(remove_x(sample), levels = remove_x(colnames(e_data))))

  p <- ggplot(e_data) + theme_bw()

  if(method=='box'){
    p <- p +
      geom_boxplot(aes(sample, .data$intensity)) +
      theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1)) +
      ylab("Feature intensity") +
      xlab("")
  }
  else if(method=='density'){
    p <- p +
      geom_density(aes(.data$intensity, col=sample)) +
      xlab("Feature intensity") +
      ylab("Density")
  }
  else if(method=='histogram'){
    p <- p +
      geom_histogram(aes(.data$intensity), bins=100) +
      scale_y_continuous(expand = c(0, 0)) +
      xlab("Feature intensity") +
      ylab("Count")
  }

  if(facet_by_sample){
    p <- p + facet_wrap(~sample)
  }
  return(p)
}
