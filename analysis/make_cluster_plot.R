require(intervals)
require(gridExtra)

cbbPalette <- c( "#000000", "#009E73","#F0E442","#D55E00", "#E69F00", "#56B4E9",  "#0072B2",  "#CC79A7")

#' Make cluster level plot
#'
#' @import ggplot2
#' @import gridExtra
#' @import intervals
#' @export 
make_cluster_plot <- function(
  exons_table = NULL, 
  meta = NULL, 
  counts = NULL, 
  intron_meta = NULL,
  snp_pos=NA,
  curv = 0.1,
  labelTextSize=3.5, # orignally set as 5
  curveExponent = 1,
  centreLineWidth = 3,   # horizontal white line to clean up the edges of the curves
  exon_height = 6,
  junction_colors = NULL,
  legend_title = NULL,
  include_legend = T,
  yOffset = 0,
  yFactor = 0.65,   # originally set as 0.65
  yConstant = 0, #-0.25 # originally set as 0.5
  geom_for_labels=geom_text,
  # gene_colors=NULL,
  length_transform = function(g){ log(g+1) }){

  meta$group=as.factor(meta$group)
  group_names=levels(meta$group)
  
  x <- meta$group

  # make sure intron_meta has "chr" in front of chromosome name so it plays nice with the exon table
  intron_meta$chr <- add_chr(as.character(intron_meta$chr))

  new_theme_empty <- theme_bw(base_size = 14 )
  new_theme_empty$panel.background = element_rect(fill="white", colour = "white")
  new_theme_empty$line <- element_blank()
  new_theme_empty$rect <- element_blank()
  new_theme_empty$strip.text <- element_blank()
  new_theme_empty$axis.text <- element_blank()
  
  group_names=sort(unique(x))
  
  exons_here=NULL
  if (!is.null(exons_table)) {
    # find the exons
    exons_here <- exons_table %>% 
      filter(chr==intron_meta$chr[1],
             ( end %in% intron_meta$start | 
                 start %in% intron_meta$end ))
  }
  
  intron_meta$id=as.factor(1:nrow(intron_meta)) # number each junction

  original_coords=c(intron_meta$start, intron_meta$end) 
  if (!is.na(snp_pos)) original_coords=c(original_coords, snp_pos)
  if (!is.null(exons_here)) original_coords=c(original_coords, exons_here$start, exons_here$end)
  original_coords=sort(unique(original_coords))
  d=original_coords[2:length(original_coords)]-original_coords[1:length(original_coords)-1] # get the difference between pairs of coordinates
  trans_d <- length_transform(d) # e.g. log(d+1), sqrt(d), d ^ .33 # apply a trasnformation function - doesn't work on negative numbers!
  coords <- c(0,cumsum(trans_d))
  names(coords)=original_coords
  
  snp_coord=coords[as.character(snp_pos)]
  
  total_length=sum(trans_d) # ==max(coords)
  
  my_xlim=c(-.05*total_length,1.05*total_length)
  
  invert_mapping=function(pos){
    if (pos %in% original_coords) coords[as.character(pos)] else
      if (pos < min(original_coords)) my_xlim[1] else
        if (pos > max(original_coords)) my_xlim[2] else {
          w=which( pos < original_coords[2:length(original_coords)] & pos > original_coords[1:(length(original_coords)-1)] )
          stopifnot(length(w)==1)
          coords[w] + (coords[w+1]-coords[w])*(pos - original_coords[w])/(original_coords[w+1]-original_coords[w])
        }
  }

  # sweep is dividing each entry in each row by the sum of all entries in that row and then apply is finding the mean value of each column
  summary_func=function(a) apply( sweep(a,1,rowSums(a),"/"),2, function(g) mean(g, na.rm=T) ) 

  last_group=group_names[length(group_names)]
  first_plot=T
  plots = foreach( plot_index = 1:length(group_names) ) %do% {

    intron_meta$prop=summary_func(counts[ group_names[plot_index]==x,,drop=F])

    group_sample_size=sum(group_names[plot_index]==x)
    #print(intron_meta$prop)

    # for each junction
    allEdges=foreach(start_i=1:2) %do% {
      foreach (i=seq(start_i,nrow(intron_meta),2), .combine = bind_rows) %do% {
        data.frame(startv = intron_meta$start[i],
          endv = intron_meta$end[i],
          start = coords[ as.character(intron_meta$start[i]) ], 
          end = coords[ as.character(intron_meta$end[i]) ],
          prop=intron_meta$prop[i],
          label=format(intron_meta$prop[i],digits=2, scientific=FALSE),
          clu= intron_meta$clu[i],
          Group = i,
          color = intron_meta$color[i], 
          stringsAsFactors = F) %>%
        mutate( l=end-start,
                xtext =start+l/2,
                ytext = ifelse(start_i==1,1,-1) * ( l^(yFactor) / 2  + yConstant)
                )
      }
    } %>% set_names(c("odd","even"))

    first_plot <- FALSE

    YLIMP <- 1.25 * max( allEdges$odd$ytext)
    YLIMN <- 1.25 * min( allEdges$even$ytext) 

    g <- ggplot() +
      geom_curve(data=allEdges$odd, aes(x = start, xend = xtext, y =  yOffset, yend = ytext, group = Group, colour = color, size = prop^curveExponent ),
                 angle=90, curvature=-curv,lineend="round") +
      geom_curve(data=allEdges$odd, aes(x = xtext, xend = end, y = ytext, yend =  yOffset, group = Group, colour = color, size = prop ^curveExponent ),
                 angle=90, curvature=-curv,lineend="round") +
      geom_curve(data=allEdges$even, aes(x = start, xend = xtext, y = -  yOffset, yend = ytext, group = Group, colour = color, size = prop ^curveExponent ),
                 angle=90,curvature=curv,lineend="round") +
      geom_curve(data=allEdges$even, aes(x = xtext, xend = end, y = ytext, yend = - yOffset, group = Group, colour = color, size = prop ^curveExponent ),
                 angle=90,curvature=curv,lineend="round") +
      
      new_theme_empty +
      ylab(NULL) +
      xlab(NULL) +
      xlim(my_xlim) +
      ggtitle(NULL, subtitle=paste0(group_names[plot_index]," (n=",group_sample_size,")" ) ) +

      # horizontal line - smooth out the ends of the curves
      geom_hline(yintercept=0, size = centreLineWidth, colour = "white") +
      geom_hline(yintercept=0,alpha=.9, size=1) +

      # label the junctions
      geom_for_labels(data=allEdges$even,aes(x=xtext,y=ytext,label=label), size = labelTextSize, colour = "black" ) +
      geom_for_labels(data=allEdges$odd,aes(x=xtext,y=ytext,label=label), size= labelTextSize,  colour = "black" ) +
      #
      ylim(YLIMN,YLIMP) +
      scale_size_continuous(limits=c(0,1),range=c(0.1,5),guide='none')
    # is this used for anything? color is currently set to clu which doesn't change for each junction

    if (!is.na(snp_coord)) {
      g <- g +geom_vline(xintercept=snp_coord) # + geom_segment(data=df,aes(x=x,y=y,xend=xend,yend=yend)) #
    }
    
    g 
  }

    # if any exons survive the cull
    if ( nrow(exons_here) > 0) {

      # fit exons within the cluster scale
      exons_here = exons_here %>% 
        mutate( x=sapply(start,invert_mapping), xend=sapply(end,invert_mapping)) %>%
        group_by(x,xend) %>% # get distinct
        sample_n(1) %>%
        ungroup()
        
      # add exons to plots
      for (i in 1:length(plots) ){ 
        plots[[i]] <- plots[[i]] +
          #geom_segment( data=exon_df, aes(x=x,y=y,xend=xend,yend=yend, colour = label), alpha=1, size=exon_height) # +
          geom_segment( data=exons_here, aes(x=x,y=0,xend=xend,yend=0), color="black", alpha=1, size=exon_height) # +
          #geom_segment( data = exon_df, aes(x = x, xend = x+0.01, y = y, yend = yend), colour = "white", size = 6, alpha = 1) +
          #geom_segment( data = exon_df, aes(x = xend-0.01, xend=xend, y = y, yend = yend), colour = "white", size = 6, alpha = 1)
      }
    }
    
    # add colour palette - different depending on whether exons are included or not - hide in top plot
  # add junction colours 
  if (is.null(junction_colors)) {
    junction_colors=cbbPalette[1+seq_along(levels(intron_meta$color))] # %>% set_names(levels(intron_meta$color))
  }

  print(levels(intron_meta$color))
  print(junction_colors)
  for(i in seq_len(length(plots)-1)) {
    plots[[i]] <- plots[[i]] + 
      scale_colour_manual("", breaks=levels(intron_meta$color), values = junction_colors ) + guides(colour=FALSE) # don't show colour legend in top plot
  }
  plots[[length(plots)]] <- plots[[length(plots)]] +
    scale_colour_manual(legend_title,  breaks=levels(intron_meta$color), values = junction_colors ) + if( include_legend ) theme(legend.position="bottom") else guides(colour=FALSE)
  
  list( plots=plots, exons=exons_here, coords=coords )
}
  

