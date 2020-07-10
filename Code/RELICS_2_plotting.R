#' @title given a set of posteriors, ignore the training layer, color the remaining ones by layer and regulatory set
#' @param input.pp: matrix of posteriors for each layer
#' @param input.seg.info: seg info
#' @param rs.cutoff: regulatory set cutoff
#' @return gene / bin output file
#' @export pp_layer_color_generator()

pp_layer_color_generator <- function(input.pp, input.seg.info, rs.cutoff = 0.1){

  sum.of.pp <- colSums(input.pp[c(2:nrow(input.pp)),])

  sum.of.pp[sum.of.pp > 1] <- 1

  input.seg.info$genomeScore <- sum.of.pp

  input.seg.info$label <- 'chr'

  color.df <- data.frame(label = 'chr',
                         color = '#000000',
                         stringsAsFactors = F)

  colors.to.use <- brewer.pal(n = 8, name = "Dark2")

  if(nrow(input.pp) > 8){
    colors.to.use <- c(colors.to.use, brewer.pal(n = 8, name = "Set2"))
  }
  if(nrow(input.pp) > 16){
    colors.to.use <- c(colors.to.use, brewer.pal(n = 8, name = "Paired"))
  }
  if(nrow(input.pp) > 24){
    print('Warning, now enough colors available!')
    break()
    # extend using: https://www.datanovia.com/en/blog/top-r-color-palettes-to-know-for-great-data-visualization/
  }


  for(i in 2:nrow(input.pp)){

    temp.pos <- which(input.pp[i,] > rs.cutoff)

    input.seg.info$label[temp.pos] <- paste0('FS ', i - 1)

    color.df <- rbind(color.df, c(paste0('FS ', i - 1), colors.to.use[i - 1]))
  }


  out.list <- list(pp_df = input.seg.info,
                   pp_cols = color.df)

  return(out.list)

}

#' @title given a set of posteriors, ignore the training layer, color the remaining ones by layer and regulatory set
#' @param input.pp: matrix of posteriors for each layer
#' @param input.seg.info: seg info
#' @param rs.cutoff: regulatory set cutoff
#' @return gene / bin output file
#' @export pp_layer_color_generator()

pp_FS_color_generator_wFS0 <- function(input.pp, input.seg.info, rs.cutoff = 0.1){

  sum.of.pp <- colSums(input.pp[c(1:nrow(input.pp)),])

  sum.of.pp[sum.of.pp > 1] <- 1

  input.seg.info$genomeScore <- sum.of.pp

  input.seg.info$label <- 'chr'

  color.df <- data.frame(label = 'chr',
                         color = '#000000',
                         stringsAsFactors = F)

  colors.to.use <- brewer.pal(n = 8, name = "Dark2")

  if(nrow(input.pp) > 8){
    colors.to.use <- c(colors.to.use, brewer.pal(n = 8, name = "Set2"))
  }
  if(nrow(input.pp) > 16){
    colors.to.use <- c(colors.to.use, brewer.pal(n = 8, name = "Paired"))
  }
  if(nrow(input.pp) > 24){
    print('Warning, now enough colors available!')
    break()
    # extend using: https://www.datanovia.com/en/blog/top-r-color-palettes-to-know-for-great-data-visualization/
  }

  temp.pos <- which(input.pp[1,] > rs.cutoff)
  input.seg.info$label[temp.pos] <- 'FS 0'
  color.df <- rbind(color.df, c('FS 0', '#FF0000'))

  for(i in 2:nrow(input.pp)){

    temp.pos <- which(input.pp[i,] > rs.cutoff)

    input.seg.info$label[temp.pos] <- paste0('FS ', i - 1)

    color.df <- rbind(color.df, c(paste0('FS ', i - 1), colors.to.use[i - 1]))
  }


  out.list <- list(pp_df = input.seg.info,
                   pp_cols = color.df)

  return(out.list)

}


#' @title given a set of posteriors, ignore the training FS, color the remaining ones by FS and regulatory set
#' @param input.pp: matrix of posteriors for each layer
#' @param input.seg.info: seg info
#' @param rs.cutoff: regulatory set cutoff
#' @return gene / bin output file
#' @export pp_layer_color_generator()

pp_FS_color_generator <- function(input.pp, input.seg.info, rs.cutoff = 0.1){

  sum.of.pp <- colSums(input.pp[c(2:nrow(input.pp)),])

  sum.of.pp[sum.of.pp > 1] <- 1

  input.seg.info$genomeScore <- sum.of.pp

  input.seg.info$label <- 'chr'

  color.df <- data.frame(label = 'chr',
                         color = '#000000',
                         stringsAsFactors = F)

  colors.to.use <- brewer.pal(n = 8, name = "Dark2")

  if(nrow(input.pp) > 8){
    colors.to.use <- c(colors.to.use, brewer.pal(n = 8, name = "Set2"))
  }
  if(nrow(input.pp) > 16){
    colors.to.use <- c(colors.to.use, brewer.pal(n = 8, name = "Paired"))
  }
  if(nrow(input.pp) > 24){
    print('Warning, now enough colors available!')
    break()
    # extend using: https://www.datanovia.com/en/blog/top-r-color-palettes-to-know-for-great-data-visualization/
  }


  for(i in 2:nrow(input.pp)){

    temp.pos <- which(input.pp[i,] > rs.cutoff)

    input.seg.info$label[temp.pos] <- paste0('FS', i - 1)

    color.df <- rbind(color.df, c(paste0('FS', i - 1), colors.to.use[i - 1]))
  }


  out.list <- list(pp_df = input.seg.info,
                   pp_cols = color.df)

  return(out.list)

}


#' @title given a set of posteriors, ignore the training layer, color the remaining ones by layer and regulatory set
#' @param input.pp: matrix of posteriors for each layer
#' @param input.seg.info: seg info
#' @return data frame: min to max location of the promoter
#' @export promoter_location_extraction()

promoter_location_extraction <- function(input.pps, input.seg.info){

  promoter.segs <- input.seg.info[which(input.pps[1,] == 1),]

  out.promoter.df <- data.frame(chrom = promoter.segs$chrom[1],
                                start = min(promoter.segs$start),
                                end = max(promoter.segs$end),
                                label = 'promoter', stringsAsFactors = F)

  return(out.promoter.df)
}

#' @title given a set of posteriors, ignore the training layer, color the remaining ones by layer and regulatory set, return list of per-layer colored
#' @param input.pp: matrix of posteriors for each layer
#' @param input.seg.info: seg info
#' @param rs.cutoff: regulatory set cutoff
#' @return list(): each element is a posterior, colored cutoffs
#' @export per_layer_coloring()

per_layer_coloring <- function(input.pps, input.seg.info, rs.cutoff = 0.1, include.L0 = FALSE){

  out.list <- list()

  colors.to.use <- brewer.pal(n = 8, name = "Dark2")

  if(nrow(input.pps) > 8){
    colors.to.use <- c(colors.to.use, brewer.pal(n = 8, name = "Set2"))
  }
  if(nrow(input.pps) > 16){
    colors.to.use <- c(colors.to.use, brewer.pal(n = 8, name = "Paired"))
  }
  if(nrow(input.pps) > 24){
    print('Warning, now enough colors available!')
    break()
    # extend using: https://www.datanovia.com/en/blog/top-r-color-palettes-to-know-for-great-data-visualization/
  }

  for(i in 2:nrow(input.pps)){
    temp.scores <- input.seg.info
    temp.scores$genomeScore <- input.pps[i,]
    temp.scores$label <- 'chr'
    temp.scores$label[which(temp.scores$genomeScore > rs.cutoff)] <- paste0('Layer ', i - 1)

    layer.color.df <- data.frame(label = 'chr',
                                 color = '#000000',
                                 stringsAsFactors = F)

    layer.color.df <- rbind(layer.color.df, c(paste0('Layer ', i - 1), colors.to.use[i - 1]))

    out.list[[i - 1]] <- list(score_df = temp.scores, color_df = layer.color.df)
  }

  # if the layer0 should be included
  if(include.L0){
    temp.scores <- input.seg.info
    temp.scores$genomeScore <- input.pps[1,]
    temp.scores$label <- 'chr'
    temp.scores$label[which(temp.scores$genomeScore > rs.cutoff)] <- paste0('Layer ', 0)

    layer.color.df <- data.frame(label = 'chr',
                                 color = '#000000',
                                 stringsAsFactors = F)

    layer.color.df <- rbind(layer.color.df, c(paste0('Layer ', 0), '#ff0000'))

    out.list[[nrow(input.pps)]] <- list(score_df = temp.scores, color_df = layer.color.df)

  }

  return(out.list)
}


#' @title given a set of posteriors, determine the % overlap with the training data for each rs
#' @param input.pp: matrix of posteriors for each layer
#' @param rs.cutoff: regulatory set cutoff
#' @return list(): each element is a posterior, colored cutoffs
#' @export layer_train_overl()

layer_train_overl <- function(input.pps, rs.cutoff = 0.1){

  out.overl.prct <- vector('numeric', length = nrow(input.pps) - 1)

  training.idx <- which(input.pps[1,] == 1)

  for(i in 2:nrow(input.pps)){
    temp.rs.idx <- which(input.pps[i,] > rs.cutoff)
    temp.prct.overl <- length(which(temp.rs.idx %in% training.idx)) / length(temp.rs.idx)
    out.overl.prct[i - 1] <- temp.prct.overl
  }

  return(out.overl.prct)
}

#' @title given a set of posteriors, the max correlation that any layer has to another layer
#' @param input.pp: matrix of posteriors for each layer
#' @param rs.cutoff: regulatory set cutoff
#' @return list(): each element is a posterior, colored cutoffs
#' @export pps_corr_overl_only()

pps_corr_overl_only <- function(input.pp, rs.cutoff, input.corr.df, layer.nr){

  if(length(unique(input.pp[1,])) > 2){
    print('Layer matrix matrix is not standard format')
    return()
  }

  layer.cor.tests <- combinations(n=nrow(input.pp), r=2,v=c(1:nrow(input.pp)), repeats.allowed=F)
  layer.cors <- matrix(0, nrow = nrow(layer.cor.tests), ncol = ncol(layer.cor.tests) + 1)

  layers.to.test <- unique(layer.cor.tests[,1])

  for(i in layers.to.test){

    layer.tests <- layer.cor.tests[which(layer.cor.tests[,1] == i), , drop = F]

    max.layer.corr <- -3
    max.layer.corr.nr <- -1

    for(j in 1:nrow(layer.tests)){

      temp.corr <- cor(input.pp[layer.tests[j,1],], input.pp[layer.tests[j,2],])

      if(temp.corr > max.layer.corr){
        max.layer.corr <- temp.corr
        max.layer.corr.nr <- layer.tests[j,2]
      }

    }

    layer.cors[i,1] <- i - 1
    layer.cors[i,2] <- max.layer.corr.nr - 1
    layer.cors[i,3] <- max.layer.corr

  }

  layer.cors <- as.data.frame(layer.cors)
  colnames(layer.cors) <- c('layer_1', 'layer_2', 'corr')
  layer.cors$layer_comb <- paste(layer.cors$layer_1, layer.cors$layer_2, sep = '_')
  layer.cors$layer_Iteration <- rep(as.character(layer.nr), nrow(layer.cors))

  out.corr.df <- rbind(input.corr.df, layer.cors)

  out.list <- list(corr_df = out.corr.df)

  return(out.list)

}
