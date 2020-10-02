#' Summarise percentage of tag intensities below a given threshold and missing
#'
#' @description Given an `MSnSet`, return a `data.frame` summarising the
#' percentage of tag intensities below a given threshold, missing values
#' and the median tag intensity
#'
#' @param obj `MSnSet`.
#' @param threshold `numeric`. Minimum intensity threshold
#' @param group_by_sample `logical`. Group the summaries by sample
#'
#' @return `data.frame` object.
#' @export
get_psm_metrics <- function(obj,
                            threshold=5.75,
                            group_by_sample=FALSE){
  e_data <- data.frame(exprs(obj))

  e_data[e_data==""] <- NA

  e_data <- e_data %>% gather(key='sample', value='intensity')

  if(group_by_sample){
    e_data_grouped <- e_data %>% group_by(sample=remove_x(sample))
  } else{
    e_data_grouped <- e_data
  }

  e_data_grouped %>%
    summarise(
      missing=sum(is.na(.data$intensity), na.rm=TRUE),
      not_missing=sum(!is.na(.data$intensity), na.rm=TRUE),
      below_thresh=sum(.data$intensity<threshold, na.rm=TRUE),
      above_thresh=sum(.data$intensity>=threshold, na.rm=TRUE),
      total=.data$missing+.data$not_missing,
      perc_below=round(100*.data$below_thresh/(.data$total),1),
      perc_missing=round(100*.data$missing/(.data$total),1),
      median_intensity=stats::median(.data$intensity, na.rm=TRUE)) %>%
    ungroup
}

#' Plot histograms for TMT tag intensities including notch annotation
#'
#' @description The distribution of TMT intensities obtained by Orbitrap has
#' been observed to have a 'notch' where few values are observed. This function
#' returns a histogram of TMT tag intensities per sample, annotated with the
#' proportion of intensities below the upper boundary of the notch region
#'
#' @param obj `MSnSet` containing PSM level TMT intensities
#' @param notch_lower `numeric`. Lower boundary of notch.
#' @param notch_upper `numeric`. Upper boundary of notch.
#' @param facet_by_sample `logical`. Facet the plot by sample.
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
    psm_metrics <- get_psm_metrics(obj, threshold=notch_upper, group_by_sample=TRUE) %>%
      mutate(label=sprintf('%s %%  ', .data$perc_below))

    p <- p + geom_text(aes(label=.data$label), data=psm_metrics,
                       family='serif',
                       x=log2(notch_lower), y=Inf,
                       vjust=2, hjust=1, size=4) +
      theme_camprot(border=FALSE, base_size=10)
  } else{
    psm_metrics <- get_psm_metrics(obj, threshold=notch_upper)

    p <- p + annotate(geom='text',
                      label=sprintf('%s %%  ', psm_metrics$perc_below),
                      family='serif',
                      x=log2(max(exprs(obj), na.rm=TRUE)), y=Inf,
                      vjust=2, hjust=1, size=5) +
      theme_camprot(border=FALSE)
  }

  invisible(p)
}

#' Summarise the number of sub-notch PSM intensities per protein
#'
#' @description The distribution of TMT intensities obtained by Orbitrap has
#' been observed to have a 'notch' where few values are observed. This function
#' summarises the number of sub-notch values observed per protein
#'
#' @param obj `MSnSet` containing PSM-level TMT intensities.
#' @param notch_upper `numeric`. Upper boundary of notch.
#' @param master_prot_col `character`. Column name for master proteins.
#'
#' @return `data.frame` Summarises # PSMs above/below notch for each protein
#' @export
get_notch_per_protein <- function(obj,
                                  master_prot_col='Master.Protein.Accessions',
                                  notch_upper=5.75){

  retain_prot <- obj %>%
    fData() %>%
    group_by_at(master_prot_col) %>%
    tally() %>%
    filter(n>1) %>%
    pull(!!sym(master_prot_col))

  notch_per_protein <- data.frame(exprs(obj)<log2(notch_upper)) %>%
    merge(fData(obj)[,master_prot_col, drop=FALSE], by='row.names') %>%
    gather(key='sample', value='below_notch', -c(.data$Row.names, !!sym(master_prot_col))) %>%
    filter(!!sym(master_prot_col) %in% retain_prot) %>%
    filter(!is.na(.data$below_notch)) %>%
    #tibble::rownames_to_column('sample') %>%
    group_by(!!sym(master_prot_col), sample) %>%
    summarise(n_psm=length(.data$below_notch),
              n_below=sum(.data$below_notch),
              fraction_below=.data$n_below/.data$n_psm) %>%
    arrange(desc(.data$fraction_below))

  return(notch_per_protein)
}


#' Plot the number of sub-notch PSM intensities per protein
#'
#' @description The distribution of TMT intensities obtained by Orbitrap has
#' been observed to have a 'notch' where few values are observed. This function
#' plots the number of sub-notch values observed per protein, as obtained using
#' `get_notch_per_protein`
#'
#' @param notch_per_protein `data.frame`. Notch summary as generated by
#' `get_notch_per_protein`
#'
#' @return `ggplot` stacked bar plot with # sub-notch PSM intensities per protein
#' @export
plot_below_notch_per_prot <- function(notch_per_protein){

  p_notch_per_protein <- notch_per_protein %>%
    group_by(.data$n_below, sample=remove_x(.data$sample)) %>%
    tally() %>%
    ggplot(aes(sample, n, fill=Hmisc::cut2(
      .data$n_below, cuts=c(0, 1, 3, 5, 8, 12, 20, 43)))) +
    geom_bar(stat='identity', position='fill') +
    scale_fill_manual(name='# PSMs below notch', values=c('grey', get_cat_palette(6))) +
    xlab('Tag') +
    scale_y_continuous(name='Fraction of proteins', expand = c(0, 0)) +
    theme_camprot(base_size=12, border=FALSE)

  return(p_notch_per_protein)

}


#' Plot the fraction of sub-notch PSM intensities per protein
#'
#' @description The distribution of TMT intensities obtained by Orbitrap has
#' been observed to have a 'notch' where few values are observed. This function
#' plots the fraction of sub-notch values observed per protein, as obtained using
#' `get_notch_per_protein`
#'
#' @param notch_per_protein `data.frame`. Notch summary as generated by
#' `get_notch_per_protein`
#'
#' @return `ggplot` stacked bar plot with # sub-notch PSM intensities per protein
#' @export
plot_fraction_below_notch_per_prot <- function(notch_per_protein){
  p_fraction_psm_below_notch <- notch_per_protein %>%
    mutate(sample=remove_x(sample)) %>%
    ggplot(aes(.data$fraction_below)) +
    geom_histogram(bins=10) +
    theme_camprot(base_size=10) +
    facet_wrap(~sample) +
    xlab('Fraction at/below notch PSMs') +
    ylab('Count')

  return(p_fraction_psm_below_notch)
}


#' Plot the missing values vs signal:noise
#'
#' @description Missing values are more frequent with low signal:noise (S:N).
#' This function visualises this relationship to aid selection of thresholds for
#' minimal S:N filtering
#'
#' @param obj `MSnSet` containing PSM level TMT intensities
#' @param sn_column `character` column name for Signal:noise values
#' @param bins `numeric` Number of bins to plot
#'
#' @return `ggplot` stacked bar plot to show S:N vs # missing values
#' @export
#' @importFrom grDevices colorRampPalette
plot_missing_SN <- function(obj,
                          sn_column="Average.Reporter.SN",
                          bins=20){

  pal <- colorRampPalette(get_cat_palette(2))

  n_missing <- obj %>% exprs() %>% is.na() %>% rowSums()

  p <- data.frame('n_missing'=n_missing,
                  'sn'=fData(obj)[[sn_column]]) %>%
    mutate(binned_sn=Hmisc::cut2(.data$sn, g=bins, digits=1)) %>%
    filter(is.finite(.data$sn)) %>%
    group_by(.data$n_missing, .data$binned_sn) %>%
    tally() %>%
    ggplot(aes(.data$binned_sn, n, fill=factor(n_missing))) +
    geom_bar(stat='identity', position='fill', colour='grey20', lwd=0.1) +
    scale_fill_manual(values=c('grey70', pal(max(n_missing))), name='Missing values') +
    guides(colour=guide_legend(override.aes = list(size = 1.5))) +
    xlab('Signal:Noise') +
    scale_y_continuous(name='Fraction', expand = c(0, 0)) +
    theme_camprot(base_size=12, border=FALSE) +
    theme(axis.text.x=element_text(size=10, angle=45, vjust=1, hjust=1))

  invisible(p)

}


#' Plot the missing values vs signal:noise for each sample
#'
#' @description Missing values are more frequent with low signal:noise (S:N).
#' This function visualises this relationship for each sample to aid selection
#' of thresholds for minimal S:N filtering.
#'
#' @param obj `MSnSet` containing PSM level TMT intensities
#' @param sn_column `character` column name for Signal:noise values
#' @param bins `numeric` Number of bins to plot
#'
#' @return `ggplot` tile plot to show S:N vs # missing values for each sample
#' @export
plot_missing_SN_per_sample <- function(obj,
                                       sn_column="Average.Reporter.SN",
                                       bins=20){

  sn_per_sample <- obj %>%
    exprs() %>%
    data.frame() %>%
    tibble::rownames_to_column('PSM_id') %>%
    gather(key='sample', value="value", -.data$PSM_id) %>%
    merge(fData(obj)[,sn_column, drop=FALSE], by.x='PSM_id', by.y='row.names') %>%
    filter(is.finite(.data$Average.Reporter.SN))

  sn_per_sample$binned_sn <- Hmisc::cut2(sn_per_sample[[sn_column]], g=bins, digits=1)

  p <- sn_per_sample %>%
    group_by(.data$sample, .data$binned_sn,
             missing=ifelse(is.na(.data$value), 'missing', 'present')) %>%
    tally() %>%
    spread(key=.data$missing, value=.data$n, fill=0) %>%
    mutate(percentage_missing=(100*.data$missing)/(.data$missing+.data$present),
           sample=remove_x(.data$sample)) %>%
    ggplot(aes(.data$binned_sn, sample, fill=.data$percentage_missing)) +
    geom_tile(colour='grey20', lwd=0.1) +
    scale_fill_gradient(low='grey97', high=get_cat_palette(1), name='Missing (%)',
                        limits=c(0,100)) +
    theme_camprot(base_size=10) +
    theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1)) +
    xlab('Signal:Noise') +
    ylab('Sample')

  return(p)
}





