require(ggplot2)

camprot_theme <- theme_bw() + theme(
  text=element_text(size=20,  family="serif"),
  panel.grid=element_blank(), aspect.ratio=1,
  panel.background = element_blank())

camprot_theme_no_box <- camprot_theme + theme(
  panel.border = element_blank(),
  axis.line = element_line(colour = "black"))

#' Generate a colour-blind friendly palette for categorical colour encoding
#'
#' @description For a given number of categories, identify a suitable
#' palette of colours which are colour-blind friendly. Palletes are derived from  
#' \url{http://mkweb.bcgsc.ca/colorblind/palettes.mhtml#page-container}
#'
#' @param n `numeric`. The number of colours required
#'
#' @return Hex codes for colour pallete
#' @export
get_cat_palette <- function(n){
  
  if((n>12 | n<1)) stop('n must be >=1 & <=12')
  
  # colour-blind friendly palettes from:
  # http://mkweb.bcgsc.ca/colorblind/palettes.mhtml#page-container
  
  cbPalette7 <- c("#2271B2", "#f0e442", "#359B73", "#d55e00",
                  "#3DB7E9", "#e69f00", "#F748A5")
  
  cbPalette12 <- c("#9F0162", "#009F81", "#FF5AAF", "#00FCCF",
                   "#8400CD", "#008DF9", "#00C2F9", "#FFB2FD",
                   "#A40122", "#E20134", "#FF6E3A", "#FFC33B")
  
  cbPalette12 <- c("#9F0162", "#009F81", "#FF5AAF", "#00FCCF",
                   "#8400CD", "#008DF9", "#00C2F9", "#FFB2FD",
                   "#A40122", "#E20134", "#FF6E3A", "#FFC33B")
  
  if(n<=7){
    return(cbPalette7[1:n])
  } else{
    return(cbPalette12[1:n])
  }
}

