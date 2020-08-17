#' Plot distributions for feature intensities per sample.
#'
#' @description Given an `MSnSet`, return a plot
#' of the feature intensities per sample. 
#'
#' @param obj `MSnSet`.
#' @param method `string` Plot style. Choice of box, density or histogramplot.
#' @param facet_by_sample Facet the plot by sample
#' 
#' @return `ggplot` object.
#' @export
plot_quant <- function(obj,
                       method=c('box', 'density', 'histogram'),
                       facet_by_sample=FALSE){
  e_data <- data.frame(exprs(obj))

  e_data[e_data==""] <- NA
  e_data <- e_data %>% gather(key='sample', value='intensity') %>%
    mutate(sample=removeX(sample))

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
      geom_histogram(aes(intensity), bins=100) +
      scale_y_continuous(expand = c(0, 0)) +
      xlab("Feature intensity") +
      ylab("Count")
  }
  
  if(facet_by_sample){
    p <- p + facet_wrap(~sample)
  }
  return(p)
}

#' Summarise percentage of tag intensities below notch and missing
#'
#' @description Given an `MSnSet`, return a `data.frame` summarising the
#' percentage of tag intensities below notch and missing
#' and the median tag intensity
#'
#' @param obj `MSnSet`.
#' @param group_by_sample Group the summaries by sample
#' @param notch_lower `numeric` Lower boundary of notch
#' @param notch_upper `numeric` Upper boundary of notch
#'  
#' @return `data.frame` object.
#' @export
get_psm_metrics <- function(obj,
                            notch_lower=3.75,
                            notch_upper=5.75,
                            group_by_sample=FALSE){
  e_data <- data.frame(exprs(obj))
  
  e_data[e_data==""] <- NA
  
  e_data <- e_data %>% gather(key='sample', value='intensity')
  
  if(group_by_sample){
    e_data_grouped <- e_data %>% group_by(sample=removeX(sample))
    } else{
  e_data_grouped <- e_data
  }

  e_data_grouped %>%
    summarise(
      missing=sum(is.na(intensity), na.rm=TRUE),
      not_missing=sum(!is.na(intensity), na.rm=TRUE),
      below_notch=sum(intensity<notch_upper, na.rm=TRUE),
      above_notch=sum(intensity>=notch_upper, na.rm=TRUE),
      perc_below=round(100*below_notch/(below_notch+above_notch),1),
      perc_missing=round(100*missing/(missing+not_missing),1),
      median_intensity=median(intensity, na.rm=TRUE)) %>%
    ungroup
}

#' Plot histograms for TMT tag intensities per sample
#'
#' @description Given an `MSnSet`, return a histogram
#' of TMT tag intensities per sample, annotated with the proportion of
#' intensities below the upper boundary of the notch region
#'
#' @param obj `MSnSet`.
#' @param facet_by_sample Facet the plot by sample
#' @param notch_lower `numeric` Lower boundary of notch
#' @param notch_upper `numeric` Upper boundary of notch
#'  
#' @return `ggplot` object.
#' @export
plot_TMT_notch <- function(obj, notch_lower=3.75, notch_upper=5.75, facet_by_sample=FALSE){
  
  colours = get_cat_palette(2)
  
  p <- obj %>%
    log(base=2) %>%
    plot_quant(method='histogram', facet_by_sample=facet_by_sample) +
    scale_y_continuous(expand = c(0, 0)) +
    xlab('PSM intensity (log2)') +
    geom_vline(xintercept=log2(notch_lower), size=0.5, colour=colours[2]) +
    geom_vline(xintercept=log2(notch_upper), size=0.5, colour=colours[1])
  
  if(facet_by_sample){
    psm_metrics <- get_psm_metrics(obj, group_by_sample=TRUE) %>%
      mutate(label=sprintf('%s %%  ', perc_below))

    p <- p + geom_text(aes(label=label), data=psm_metrics,
                       family='serif',
                       x=log2(notch_lower), y=Inf,
                       vjust=2, hjust=1, size=4) +
      theme_camprot(border=FALSE, base_size=10)
  } else{
    
    psm_metrics <- get_psm_metrics(obj)
    
    p <- p + annotate(geom='text',
                      label=sprintf('%s %%  ', psm_metrics$perc_below),
                      family='serif',
                      x=log2(notch_lower), y=Inf,
                      vjust=2, hjust=1, size=8) +
      theme_camprot(border=FALSE)
  }
  
  invisible(p)
}
