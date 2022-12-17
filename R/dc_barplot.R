
dc_barplot <- function(dc_obj,
                       samp_notes,
                       annots = NULL,
                       proportional = TRUE,
                       binned_graph = TRUE) {

  if (length(annots) > 4){
    warning(cat('Too many provided annotations, will graph only first 4'))
    annots <- annots[1:4]
  }
  samp_notes$Sample_ID <- rownames(samp_notes)

  if (proportional == TRUE & binned_graph == TRUE) {
    keepers <- which(apply(dc_obj$prop_of_all, 2, function(x) any(!is.na(x))))
    prop_df <- dcdat_df <- as.data.frame(dc_obj$prop_of_all[, keepers])
    yaxLab <- 'Proportional Cell Abundance'
  } else if (proportional == FALSE & binned_graph == FALSE) {
    keepers <- which(apply(dc_obj$beta.granular, 2, function(x) any(!is.na(x))))
    dcdat_df <- as.data.frame(dc_obj$beta.granular[, keepers])
    prop_df <- sweep(dc_obj$beta.granular, 2, apply(dc_obj$beta.granular, 2, sum), '/')
    prop_df <- as.data.frame(prop_df[, keepers])
    yaxLab <- 'Estimated Cell Abundance'
  } else if (proportional == FALSE & binned_graph == TRUE) {
    keepers <- which(apply(dc_obj$beta, 2, function(x) any(!is.na(x))))
    dcdat_df <- as.data.frame(dc_obj$beta[, keepers])
    prop_df <- as.data.frame(dc_obj$prop_of_all[, keepers])
    yaxLab <- 'Estimated Cell Abundance'
  } else {
    #create granular proportional
    dcdat_df <- as.data.frame(dc_obj$beta.granular)
    dcdat_df <- sweep(dcdat_df, 2, apply(dcdat_df, 2, sum), '/')
    keepers <- which(apply(dcdat_df, 2, function(x) any(!is.na(x))))
    prop_df <- dcdat_df <- dcdat_df[, keepers]
    yaxLab <- 'Proportional Cell Abundance'
  }
  # Filter to just valid cell types
  valcells <- rownames(dcdat_df)[which(rowSums(dcdat_df) > 0)]
  dcdat_df <- dcdat_df[valcells, ]
  prop_df <- prop_df[valcells, ]
  samp_notes <- samp_notes[keepers, ]

  #set colors for cell types from available broad palette
  av_clrs <- c('#8DD3C7','#FFFFB3','#BEBADA','#FB8072','#80B1D3','#FDB462','#B3DE69',
               '#FCCDE5','#A6CEE3','#1F78B4','#B2DF8A','#33A02C','#FB9A99','#E31A1C',
               '#FDBF6F','#FF7F00','#1B9E77','#D95F02','#7570B3','#E7298A','#66A61E',
               '#E6AB02','#A6761D','#666666', sample(grDevices::colors(), 150))
  cell_clrs <- av_clrs[1:length(valcells)]
  names(cell_clrs) <- valcells

  cat('Producing Annotated Bargraph...\n')
  # prep the matrix for ggplot
  grph_df <- dcdat_df %>%
    dplyr::mutate(cells = rownames(dcdat_df)) %>%
    reshape2::melt(id="cells")

  # tack on columns for labeling, final renaming of columns
  sampmatch <- match(grph_df$variable, samp_notes$Sample_ID)
  grph_df[, annots] <- samp_notes[sampmatch, annots]
  colnames(grph_df) <- c('Cells','Sample_ID','Proportions',str_to_title(annots))

  # cull colors to just what is being graphed
  # cell types is already accurate, but gather the annotation colors
  ann_clrs <- annot_cols(samp_notes, annots)
  # combined color group
  comb_col <- c(cell_clrs, ann_clrs)

  #build out a dataframe for an annotation bar, compatible with the bar plot above
  #y-axis placement must be programmatically derived
  if (proportional == T) {
    brkdown <- c(0, 0.2, 0.4, 0.6, 0.8, 1)
    y_step <- 0.05
    yval <- -0.025
    laby_df <- data.frame(brk = brkdown,
                          laby = c('0%', '20%', '40%', '60%', '80%', '100%'))
  } else {
    brkdown <- as.numeric(pretty_breaks()(0:(max(apply(dcdat_df, 2, sum)))))
    y_step <- max(brkdown) * 0.05
    yval <- -y_step/2
    laby_df <- data.frame(brk = c(brkdown))
  }

  if(!is.null(annots)) {
    annot_df <- data.frame()
    for (i in annots) {
      annot_df <- rbind(annot_df, data.frame(grp = samp_notes[[i]],
                                             x = seq(1:nrow(samp_notes)),
                                             y = yval))
      laby_df <- rbind(laby_df, data.frame(brk = yval, laby = str_to_title(i)))
      yval <- yval - y_step
    }
  }
  laby_df <- laby_df %>% arrange(brk)

  # Basic BarPlot, no annotations
  deconplot <- ggplot(grph_df, aes_string(x = 'Sample_ID', y = 'Proportions',
                                          fill = 'Cells')) +
    geom_col(color = 'white', size = 0.05) +
    scale_fill_manual(values = comb_col) +
    theme_bw() +
    scale_x_discrete(labels = as.character(samp_notes[['Sample_ID']])) +
    labs(x = 'ROI #', y = yaxLab, title = 'Cell Deconvolution') +
    theme(axis.text.y = element_text(size = 10),
          axis.text.x = element_text(size = ifelse(nrow(samp_notes) < 65, 10, 8),
                                     angle = ifelse(nrow(samp_notes) < 65, 45, 90),
                                     vjust = 0.5, hjust = 0.5),
          axis.title = element_text(size = 12),
          plot.title = element_text(size = 16, face = 'bold', hjust = 0.5),
          legend.position="none") +
    scale_y_continuous(labels = laby_df$laby,
                       breaks = laby_df$brk,
                       expand = expansion(mult = 0)) #+
    # geom_text(x=0, y=-0.075, label='test ', hjust = 1)

  deconplot <- deconplot +
    geom_tile(data = annot_df, aes(x = x, y = y, height = y_step, fill = grp))

  basicbar <- ggplot(data = grph_df, aes_string(x = 'Sample_ID', y = 'Proportions', fill = 'Cells')) +
    geom_col() +
    scale_fill_manual(values = cell_clrs, name = "Cell Types")
  bb_leg <- get_legend(basicbar)

  annot_tile <- ggplot(data = annot_df, aes(x = x, y = y, fill = grp)) +
    geom_tile() +
    scale_fill_manual(values = ann_clrs, name = "Annotations")
  ann_leg <- get_legend(annot_tile)

  #determine grid ratios for plot

  legw <- max(c(unlist(bb_leg$widths)[5],unlist(ann_leg$widths)[5]))*0.15/5.5
  wds <- c(1-legw, legw)
  legh <- as.numeric(ann_leg$heights[[3]])/(as.numeric(ann_leg$heights[[3]]) +
                                              as.numeric(bb_leg$heights[[3]]))
  hts <- c(1-legh, legh)

  #create gridded graph
  final_gr <- grid.arrange(deconplot, bb_leg, ann_leg,
                           ncol=2, nrow = 2,
                           layout_matrix = rbind(c(1,2), c(1,3)),
                           widths = wds, heights = hts)

  return(final_gr)

}

annot_cols <- function(samp_notes,
                       annots) {
  require(RColorBrewer)
  # empty list, palette order
  col_vct <- c()
  pals <- c('Set3','Set1','Set2','Dark2')

  for (i in seq(1:length(annots))) {
    col_gr <- brewer.pal.info[pals[i], "maxcolors"]
    getPalette <- colorRampPalette(brewer.pal(col_gr, pals[i]))
    anns <- unique(samp_notes[[annots[i]]])
    anns <- anns[!is.na(anns)]
    colz <- getPalette(max(length(anns), col_gr))[1:length(anns)]
    names(colz) <- anns

    #append
    col_vct <- c(col_vct, colz)
  }

  return(col_vct)
}

