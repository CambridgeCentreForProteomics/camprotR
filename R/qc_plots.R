#' Plot distributions for feature intensities per sample.
#'
#' @description Given an `MSnSet`, return a plot
#' of the feature intensities per sample. 
#'
#' @param obj `MSnSet`.
#' @param method `string` Plot style. Choice of box, density or histogramplot.
#'
#' @return `ggplot` object.
#' @export
plot_label_quant <- function(obj, method=c('box', 'density', 'histogram')){
  e_data <- data.frame(exprs(obj))

  e_data[e_data==""] <- NA
  e_data <- e_data %>% gather(key='sample', value='intensity')

  p <- ggplot(e_data) + theme_bw()
    
  if(method=='box'){
    p <- p +
      geom_boxplot(aes(sample, intensity)) +
      theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1)) +
      ylab("Feature intensity") +
      xlab("")
  }
  else if(method=='density'){
    p <- p +
      geom_density(aes(intensity, col=sample)) +
      xlab("Feature intensity") +
      ylab("Density")
  }
  else if(method=='histogram'){
    p <- p +
      geom_histogram(aes(intensity), bins=50) +
      xlab("Feature intensity") +
      ylab("Density") +
      facet_wrap(~sample)
  }

  return(p)
}
