
print('Loading Libraries...')
library(ggplot2)
library(grid)
library(gridExtra)
library(GenomicRanges)


if("biomaRt" %in% rownames(installed.packages()) == TRUE){
  library(biomaRt)
}

#' @title Generic wrapper function for plotting
#' @param input.df: data frame containing the scores and the coordinates, $chrom, $start, $end, $label
#' @param score.name: name of the column containing the scores input.df[[input.analysis.type]]
#' @param label.colors: data frame containing the color to label key, $color and $label
#' @return NA
#' @export plot_bedgraphFormat()

plot_bedgraphFormat <- function(input.df.orig, score.name, label.colors = NULL){

  input.df <- input.df.orig
  if('nrSupportGuides' %in% names(input.df)){
    no.guides.index <- which(input.df$nrSupportGuides == 0)
    if(length(no.guides.index) > 0){
      input.df <- input.df[-no.guides.index,]
    }
  }

  chrom.df <- c()  # data frame containing only valid chromosomes
  na.df <- NULL
  if(length(which(is.na(input.df$start))) > 0){
    chrom.df <- input.df[which(! is.na(input.df$start)),] #filter(input.df, grepl('chr',chrom))
    #chrom.df$midpoint <- round((chrom.df$start + chrom.df$end) / 2)
    na.df <- input.df[which(is.na(input.df$start)),] #filter(input.df, grepl('NA',chrom))

    all.na.labels <- unique(na.df$label)
    temp.inter.chr.space <- round(1/8* (max(chrom.df$end) - min(chrom.df$start))) + 1
    for(lab in all.na.labels){
      na.lab.df <- na.df[na.df$label == lab,]
      na.lab.start <- max(chrom.df$end) + temp.inter.chr.space
      na.lab.stepSize <- median(sort(chrom.df$start)[2:length(chrom.df$start)] - sort(chrom.df$end)[1:(length(chrom.df$end) - 1)])  # adjust step size to median step size of targeting guides
      if(is.na(na.lab.stepSize)){
        na.lab.stepSize <- 1
      }
      na.lab.x.start <- na.lab.start + na.lab.stepSize*c(1:nrow(na.lab.df))
      na.lab.x.end <- na.lab.x.start + na.lab.stepSize - 1

      na.lab.df$start <- na.lab.x.start
      na.lab.df$end <- na.lab.x.end
      na.lab.df$chrom <- rep(chrom.df$chrom[1], nrow(na.lab.df))

      chrom.df <- rbind(chrom.df, na.lab.df)
      # temp.x.pos <- c(temp.x.pos, na.lab.x.pos)
      # temp.y.pos <- c(temp.y.pos, na.lab.df[[score.name]])
      # temp.labels <- c(temp.labels, na.lab.df$label)
    }

  } else {
    chrom.df <- input.df
    #chrom.df$midpoint <- round((chrom.df$start + chrom.df$end) / 2)
  }

  nr.chroms <- unique(chrom.df$chrom)

  if(length(nr.chroms) > 1){
    print('Warning! Multi-chromosome plotting not established yet')
    break()
  } else {
    # plot.values <- singlePlotSteup(chrom.df, na.df, score.name, label.colors)
    bedgraph_plotting(chrom.df, score.name, label.colors = label.colors)
  }
}

#' @title Generic wrapper function for plotting
#' @param input.df: data frame containing the scores and the coordinates, $chrom, $start, $end, $label
#' @param score.name: name of the column containing the scores input.df[[input.analysis.type]]
#' @param label.colors: data frame containing the color to label key, $color and $label
#' @return NA
#' @export bedgraph_plotting()

bedgraph_plotting <- function(input.df, score.name, label.colors = NULL){

  chrom.info <- c()
  if(! is.null(label.colors)){
    if(sum(! input.df$label %in% label.colors$label)){
      print('Warning, color key does not contain all labels')
    }

    chrom.info <- input.df

    chrom.info$color <- rep(label.colors$color[1], nrow(chrom.info))

    for(i in 2:nrow(label.colors)){
      chrom.info$color[which(chrom.info$label == label.colors$label[i])] <- label.colors$color[i]
    }

    chrom.info$label <- factor(chrom.info$label, levels = label.colors$label[1:nrow(label.colors)])
    chrom.info$color <- factor(chrom.info$color, levels = label.colors$color[1:nrow(label.colors)])

    out.chrom.plot <- ggplot(chrom.info)+
      geom_rect(aes(xmin = start, xmax = end, ymin = 0, ymax = genomeScore,
                    color = label, fill = label )) +
      scale_color_manual(guide = FALSE, values = as.character(label.colors$color),
                         breaks = as.character(label.colors$label))+
      scale_fill_manual(name = 'Labels',
                        values = adjustcolor(as.character(label.colors$color), alpha.f = 0.5)) +
      guides(fill = guide_legend(override.aes = list(alpha = 0.5))) +
      theme_bw()
    print(out.chrom.plot)

  } else {
    chrom.info <- input.df
    out.chrom.plot <- c()

    out.chrom.plot <- ggplot(chrom.info)+
      geom_rect(aes(xmin = start, xmax = end, ymin = 0, ymax = genomeScore,
                    color = label, fill = label))
    out.chrom.plot <- out.chrom.plot + theme_bw()

    print(out.chrom.plot)
  }
}

#' @title Generic wrapper function for plotting
#' @param input.df: data frame containing the scores and the coordinates, $chrom, $start, $end, $label
#' @param score.name: name of the column containing the scores input.df[[input.analysis.type]]
#' @param label.colors: data frame containing the color to label key, $color and $label
#' @return NA
#' @export
#' plot_scores()

plot_scores <- function(input.df.orig, score.name, label.colors = NULL){

  input.df <- input.df.orig
  if('nrSupportGuides' %in% names(input.df)){
    no.guides.index <- which(input.df$nrSupportGuides == 0)
    if(length(no.guides.index) > 0){
      input.df <- input.df[-no.guides.index,]
    }
  }

  chrom.df <- c()  # data frame containing only valid chromosomes
  na.df <- NULL
  if(length(which(is.na(input.df$start))) > 0){
    chrom.df <- input.df[which(! is.na(input.df$start)),] #filter(input.df, grepl('chr',chrom))
    chrom.df$midpoint <- round((chrom.df$start + chrom.df$end) / 2)
    na.df <- input.df[which(is.na(input.df$start)),] #filter(input.df, grepl('NA',chrom))
  } else {
    chrom.df <- input.df
    chrom.df$midpoint <- round((chrom.df$start + chrom.df$end) / 2)
  }

  nr.chroms <- unique(chrom.df$chrom)

  if(length(nr.chroms) > 1){
    print('Warning! Multi-chromosome plotting not established yet')
    break()
  } else {
    plot.values <- singlePlotSteup(chrom.df, na.df, score.name, label.colors)
  }
}

#' @title Wrapper for plotting the scores for a single chromosome
#' @param input.df: data frame containing the scores and the coordinates, $chrom, $start, $end
#' @param input.na.df: data frame containing the scores and the coordinates for non-targeting guides, $chrom, $start, $end
#' @param score.name: name of the column containing the scores input.df[[input.analysis.type]]
#' @param label.colors: data frame containing the color to label key, $color and $label
#' @return plot
#' @export NA
#' singlePlotSteup()

singlePlotSteup <- function(input.df, input.na.df, score.name, label.colors = NULL){

  chrom.df <- input.df
  na.df <- input.na.df

  temp.x.pos <- chrom.df$midpoint
  temp.y.pos <- chrom.df[[score.name]]
  temp.labels <- chrom.df$label
  temp.inter.chr.space <- round(1/8* (max(temp.x.pos) - min(temp.x.pos))) + 1

  if(! is.null(na.df) ){
    all.na.labels <- unique(na.df$label)
    for(lab in all.na.labels){
      na.lab.df <- na.df[na.df$label == lab,]
      na.lab.start <- max(temp.x.pos) + temp.inter.chr.space
      na.lab.stepSize <- median(sort(temp.x.pos)[2:length(temp.x.pos)] - sort(temp.x.pos)[1:(length(temp.x.pos) - 1)])  # adjust step size to median step size of targeting guides
      if(is.na(na.lab.stepSize)){
        na.lab.stepSize <- 1
      }
      na.lab.x.pos <- na.lab.start + na.lab.stepSize*c(1:nrow(na.lab.df))

      temp.x.pos <- c(temp.x.pos, na.lab.x.pos)
      temp.y.pos <- c(temp.y.pos, na.lab.df[[score.name]])
      temp.labels <- c(temp.labels, na.lab.df$label)
    }
  }

  create_chromosome_scorePlot(temp.x.pos, temp.y.pos, temp.labels, label.colors)

}

#' @title Plot the scores for a set of coordinates and scores
#' @param input.df: data frame containing the scores and the coordinates, $chrom, $start, $end
#' @param input.na.df: data frame containing the scores and the coordinates for non-targeting guides, $chrom, $start, $end
#' @param score.name: name of the column containing the scores input.df[[input.analysis.type]]
#' @param label.colors: data frame containing the color to label key, $color and $label
#' @return plot
#' @export create_chromosome_scorePlot()

create_chromosome_scorePlot <- function(input.x.values, input.y.values, input.labels, label.colors = NULL){

  chrom.info <- c()
  if(! is.null(label.colors)){
    if(sum(! input.labels %in% label.colors$label)){
      print('Warning, color key does not contain all labels')
    }

    chrom.info <- cbind.data.frame(coords = input.x.values, scores = input.y.values,
                                   labels = input.labels, stringsAsFactors = F)
    chrom.info$color <- rep(label.colors$color[1], nrow(chrom.info))

    for(i in 2:nrow(label.colors)){
      chrom.info$color[which(chrom.info$labels == label.colors$label[i])] <- label.colors$color[i]
    }

    chrom.info$labels <- factor(chrom.info$labels, levels = label.colors$label[1:nrow(label.colors)])
    chrom.info$color <- factor(chrom.info$color, levels = label.colors$color[1:nrow(label.colors)])
    #
    # label.colors$color <- factor(label.colors$color, levels = label.colors$color[1:nrow(label.colors)])
    # label.colors$label <- factor(label.colors$label, levels = label.colors$label[1:nrow(label.colors)])
    #
    # out.chrom.plot <- ggplot()+
    #   geom_point(data = chrom.info, aes(x = coords, y = scores, colour = labels), size = 1) +
    #   scale_color_manual(values = label.colors$color, breaks = label.colors$label) +
    #   theme_bw()
    # print(out.chrom.plot)
    #
    out.chrom.plot <- ggplot()+
      geom_point(data = chrom.info, aes(x = coords, y = scores, colour = labels), size = 1) +
      scale_color_manual(values = as.character(label.colors$color),
                         breaks = as.character(label.colors$label)) +
      theme_bw()
    print(out.chrom.plot)

  } else {
    chrom.info <- cbind.data.frame(coords = input.x.values, scores = input.y.values, labels = input.labels)
    out.chrom.plot <- c()

    out.chrom.plot <- ggplot()+
      geom_point(data = chrom.info, aes(x = coords, y = scores, colour = labels), size = 1)

    out.chrom.plot <- out.chrom.plot + theme_bw() + scale_colour_discrete(name = 'Scores')

    print(out.chrom.plot)
  }

}

#' @title Adjust the labels which overlap
#' @param input.query.df: data frame containing: $chrom, $start, $end, $formatScores
#' @param subject.df: data frame containing: $chrom, $start, $end, $formatScores
#' @param score.name: name of the column containing the scores input.df[[input.analysis.type]]
#' @return plot
#' @export label_overlaps()

label_overlaps <- function(input.query.df, input.subject.df, overl.label){

  query.df <- c()
  non.range.query.df <- c()
  subject.df <- c()
  non.range.subject.df <- c()

  if(length(which(is.na(input.query.df$start))) > 0){
    query.df <- input.query.df[-which(is.na(input.query.df$start)),]
    non.range.query.df <- input.query.df[which(is.na(input.query.df$start)),]
  } else {
    query.df <- input.query.df
  }
  if(length(which(is.na(input.subject.df$start))) > 0){
    subject.df <- input.subject.df[-which(is.na(input.subject.df$start)),]
    non.range.subject.df <- input.subject.df[which(is.na(input.subject.df$start)),]
  } else {
    subject.df <- input.subject.df
  }


  query.ranges <- GRanges(seqnames = query.df$chrom, ranges = IRanges(query.df$start, query.df$end))
  subject.ranges <- GRanges(seqnames = subject.df$chrom, ranges = IRanges(subject.df$start, subject.df$end))

  range.overlaps <- as.data.frame(findOverlaps(query.ranges, subject.ranges, type = 'any'))

  query.df$label[unique(range.overlaps$queryHits)] <- overl.label

  out.df <- query.df
  if(length(which(is.na(input.query.df$start))) > 0){
    out.df <- rbind(out.df, non.range.query.df)
  }
  return(out.df)
}

#' @title Adjust the labels which overlap
#' @param input.df: data frame containing: $chrom, $start, $end, and scores
#' @param score.name: column name which contains the scores
#' @param out.name: name of the output bedgraph
#' @param out.path: path for output
#' @param start.chrom: default: NULL, will obtain from data. if multiple chromosomes present, specify view of chromosome when in browser
#' @param split.score: value used to split scores into top and bottom tracks
#' @param color.id: color of the bedgraph track
#' @return two bedgraph files, split into top and bottom halves
#' @export create_splitBedGraph()

create_splitBedGraph <- function(input.df, score.name, out.name, out.path, start.chrom = NULL, split.score = 0, color.id = 'red'){

  filtered.df <- input.df[-which(is.na(input.df$start)),]

  top.df <- filtered.df[which(filtered.df[[score.name]] > split.score),]
  bottom.df <- filtered.df[which(filtered.df[[score.name]] <= split.score),]

  rgb.color <- col2rgb(color.id)
  bg.color <- paste(rgb.color[1,1], rgb.color[2,1], rgb.color[3,1], sep = ',')

  if(is.null(start.chrom)){
    start.chrom <- filtered.df$chrom[1]
  }

  # the assumption is that there are a lot more bottom hits.
  # will use the range of bottom df to guide plot range

  start.chrom.range <- bottom.df[which(bottom.df$chrom == start.chrom),]
  min.start <- min(start.chrom.range$start)
  max.end <- max(start.chrom.range$end)

  header1 <- paste0('browser position ', start.chrom, ':', min.start,
                         '-', max.end, '\n')
  top.header2 <- paste0("track type=bedGraph name='", out.name, "_top' description='",
                        out.name, "_top' visibility=full color=", bg.color, '\n')

  bottom.header2 <- paste0("track type=bedGraph name='", out.name, "_bottom' description='",
                        out.name, "_bottom' visibility=full color=", bg.color, '\n')

  top.bedgraph <- set_up_bedGraphTable(top.df, score.name)
  bottom.bedgraph <- set_up_bedGraphTable(bottom.df, score.name)

  top.file <- paste0(out.path, out.name, '_top', ".bedgraph")
  cat(header1, file = top.file)
  cat(top.header2, file = top.file, append = TRUE)
  write.table(top.bedgraph, file = top.file, append = TRUE, row.names = F, col.names = F, quote = F, sep = '\t')

  bottom.file <- paste0(out.path, out.name, '_botom', ".bedgraph")
  cat(header1, file = bottom.file)
  cat(bottom.header2, file = bottom.file, append = TRUE)
  write.table(bottom.bedgraph, file = bottom.file, append = TRUE, row.names = F, col.names = F, quote = F, sep = '\t')

}

#' @title Adjust the labels which overlap
#' @param input.df: data frame containing: $chrom, $start, $end, and scores
#' @param score.name: column name which contains the scores
#' @return date frame in bedgraph format
#' @export set_up_bedGraphTable()

set_up_bedGraphTable <- function(input.df, score.name){

  temp.score.gr <- GRanges(seqnames = input.df$chrom, ranges = IRanges(input.df$start, input.df$end),
                           score = input.df[[score.name]])
  temp.score.gr <- sortSeqlevels(temp.score.gr)
  temp.score.gr <- sort(temp.score.gr)
  df.temp.score <- as.data.frame(temp.score.gr)
  df.temp.score.final <- df.temp.score[,c(1,2,3,6)]

  return(df.temp.score.final)
}

#' @title Adjust the labels which overlap
#' @param data.list: list of data frames: data frame containing: $chrom, $start, $end, and formatScores
#' @param main.data.name: name of data frame who's peak are displayed at the top
#' @param data.to.correlate: vector of names of other data frames which to overlap with
#' @param input.colors: colors for each of the tracks
#' @param output.file: location and name where to place the pdf and the .csv table containing all the values
#' @param top.percent: top fraction of the data to plot (1 is 100% of data)
#' @param input.binarized: binary vector whether data to overlap with is just binary or not
#' @param score.name: y-axix label of main plot
#' @return date frame in bedgraph format
#' @export display_peak_max()

display_peak_max <- function(data.list, main.data.name, data.to.correlate, input.colors,
                             output.file, top.percent, input.binarized, score.name){

  #
  # library(ggplot2)
  # library(grid)
  # library(gtable)
  # library(dplyr)

  overlap.score.list <- list()
  main.data <- data.list[[main.data.name]]
  if(! 'formatScores' %in% colnames(main.data)){
    if('genomeScore' %in% colnames(main.data)){
      main.data$formatScores <- main.data$genomeScore
    } else if('guideScore' %in% colnames(main.data)){
      main.data$formatScores <- main.data$guideScore
    } else {
      print('formatScores column required!')
      break()
    }
  }
  ordered.main.data <- main.data[order(-main.data$formatScores),]
  ordered.main.data$ranks <- c(1:nrow(ordered.main.data))
  ordered.main.data <- ordered.main.data[c(1:round(nrow(ordered.main.data) * top.percent)),]

  ordered.main.ranges <- GRanges(seqnames = ordered.main.data$chrom,
                                 ranges = IRanges(ordered.main.data$start, ordered.main.data$end))
  save.overlap.scores <- ordered.main.data

  for(i in 1:length(data.to.correlate)){
    temp.data.name <- data.to.correlate[i]
    temp.data <- data.list[[temp.data.name]]
    if(! 'formatScores' %in% colnames(temp.data)){
      if('genomeScore' %in% colnames(temp.data)){
        temp.data$formatScores <- temp.data$genomeScore
      } else if('guideScore' %in% colnames(temp.data)){
        temp.data$formatScores <- temp.data$guideScore
      } else {
        print('formatScores column required!')
        break()
      }
    }
    if(length(which(is.na(temp.data$start))) > 0){
      temp.data <- temp.data[-which(is.na(temp.data$start)),]
    }
    temp.data.ranges <- GRanges(seqnames = temp.data$chrom,
                                ranges = IRanges(temp.data$start, temp.data$end))

    temp.overlaps <- as.data.frame(findOverlaps(ordered.main.ranges,
                                                temp.data.ranges, type = 'any'))

    # iterate through all the main peaks
    overl.main.peaks <- unique(temp.overlaps$queryHits)
    overl.temp.peaks <- unique(temp.overlaps$subjectHits)
    temp.out.scores <- vector('numeric', length = nrow(ordered.main.data))
    temp.out.start <- vector('numeric', length = nrow(ordered.main.data))
    temp.out.end <- vector('numeric', length = nrow(ordered.main.data))

    for(j in 1:nrow(ordered.main.data)){
      if(j %in% overl.main.peaks){
        temp.overl.rows <- which(temp.overlaps$queryHits == j)
        temp.data.rows <- temp.overlaps$subjectHits[temp.overl.rows]
        temp.overl.data <- temp.data[temp.data.rows,]
        temp.overl.data.max <- max(temp.overl.data$formatScores)
        temp.out.scores[j] <- temp.overl.data.max
      } else {
        temp.out.scores[j] <- 0
      }
    }

    save.overlap.scores[[temp.data.name]] <- temp.out.scores
    temp.data.scores <- data.frame(chrom = ordered.main.data$chrom,
                                   start = ordered.main.data$start, end = ordered.main.data$end, formatScores = temp.out.scores)
    temp.data.scores$ranks <- c(1:nrow(ordered.main.data))
    overlap.score.list[[temp.data.name]] <- temp.data.scores

  }

  # 3. and 4. Find overlaps and plot
  plot.list <- list()
  plot.color.gradient <- ggplot(ordered.main.data) +
    geom_segment(aes(x = ranks, xend = ranks, y = 0, yend = formatScores, colour = formatScores)) + scale_colour_gradient(low = 'green', high = 'purple') +
    theme(legend.position = 'none', axis.title.x=element_blank(), axis.text.x = element_blank(), axis.ticks.x=element_blank(),
          axis.title.y = element_text(angle = 0)) +  #theme(legend.position = c(0.9, 0.2))
    ylab(score.name)  # scale_x_continuous(limits = c(0, nrow(ordered.data.scores))) +
  plot.list$gradientScore <- plot.color.gradient
  to.plot0 <- plot.color.gradient
  gtab0 <- ggplotGrob(to.plot0)
  main.range <- 3
  gtab0$heights[7] <- unit(main.range, 'null')

  i <- 0
  if(length(overlap.score.list) > 1){
    for(i in 1:(length(overlap.score.list) - 1) ){
      if(input.binarized[i] == 0){  # if the data is not a binary (i.e. winthing TAD or not)
        eval(parse(text = paste0('to.plot', i, ' <- ggplot(overlap.score.list[[i]]) +
          geom_segment(aes(x = ranks, xend = ranks, y = 0, yend = formatScores), colour = input.colors[i] ) +
          ylab(names(overlap.score.list)[i]) + theme(axis.title.x=element_blank(), axis.text.x = element_blank(),
                                 axis.ticks.x=element_blank(), axis.title.y = element_text(angle = 0))' )))
        eval(parse(text = paste0('gtab', i, ' <- ggplotGrob(to.plot', i, ')' ) ))
        temp.range <- main.range
        eval(parse(text = paste0('gtab', i, "$heights[7] <- unit(temp.range, 'null')" ) ))
      } else {
        eval(parse(text = paste0('to.plot', i, ' <- ggplot(overlap.score.list[[i]]) +
          geom_segment(aes(x = ranks, xend = ranks, y = 0, yend = formatScores), colour = input.colors[i] ) +
          ylab(names(overlap.score.list)[i]) + theme(axis.title.x=element_blank(), axis.text.x = element_blank(), axis.ticks.x=element_blank(),
          axis.text.y = element_blank(), axis.ticks.y=element_blank(), axis.title.y = element_text(angle = 0))' )))
        eval(parse(text = paste0('gtab', i, ' <- ggplotGrob(to.plot', i, ')' ) ))
        temp.range <- 1
        eval(parse(text = paste0('gtab', i, "$heights[7] <- unit(temp.range, 'null')" ) ))
      }

    }
  }

  i <- i + 1
  if(input.binarized[i] == 0){  # if the data is not a binary (i.e. withing TAD or not)
    eval(parse(text = paste0('to.plot', i, ' <- ggplot(overlap.score.list[[i]]) +
      geom_segment(aes(x = ranks, xend = ranks, y = 0, yend = formatScores), colour = input.colors[i] ) +
      theme(axis.title.y = element_text(angle = 0)) +
      ylab(names(overlap.score.list)[i])' )))
    eval(parse(text = paste0('gtab', i, ' <- ggplotGrob(to.plot', i, ')' ) ))
    temp.range <- main.range
    eval(parse(text = paste0('gtab', i, "$heights[7] <- unit(temp.range, 'null')" ) ))
  } else {
    eval(parse(text = paste0('to.plot', i, ' <- ggplot(overlap.score.list[[i]]) +
      geom_segment(aes(x = ranks, xend = ranks, y = 0, yend = formatScores), colour = input.colors[i] )  +
      theme(axis.title.y = element_text(angle = 0)) +
      ylab(names(overlap.score.list)[i])' )))  #+ theme(axis.text.y = element_blank(), axis.ticks.y=element_blank())
    eval(parse(text = paste0('gtab', i, ' <- ggplotGrob(to.plot', i, ')' ) ))
    temp.range <- 1
    eval(parse(text = paste0('gtab', i, "$heights[7] <- unit(temp.range, 'null')" ) ))
  }

  eval(parse(text = paste0('gtab_out <- rbind(', paste0(paste0('gtab', c(0:length(overlap.score.list))), collapse = ','), ", size = 'first')")))
  eval(parse(text = paste0('gtab_out$widths <- unit.pmax(' ,paste0(paste0('gtab', c(0:length(overlap.score.list)), '$widths'), collapse = ','), ')')))

  pdf(paste0(output.file,'_corrPlots.pdf'), width = 14, height = 7)
  grid.newpage()
  grid.draw(gtab_out)
  dev.off()

  write.csv(save.overlap.scores, file = paste0(output.file, '_corrTable.csv'), row.names = F)

}

#' @title Generic wrapper function for plotting
#' @param input.df.orig.list: list of data frames containing the scores and the coordinates, $chrom, $start, $end, $label
#' @param score.names: vector of names of the column containing the scores input.df[[input.analysis.type]] or 'bedgraph' to specify y-axis height does not matter
#' @param label.color.list: data frame containing the color to label key, $color and $label
#' @param x.min: optional parameter: manually specify the minimum range of the plot
#' @param x.max: optional parameter: manually specify the maximum range of the plot
#' @param track.height: vector length of input data but only applied to bedgraph. Height of each bedgraph track: default 3
#' @param genome.to.use: default hg19
#' @param only.gene: if many genes in region, display name of only this gene
#' @return NA
#' @export plot_tracks()

plot_tracks <- function(input.df.orig.list, score.names, label.color.list,
                              x.min = NULL, x.max = NULL, track.height = NULL,
                              genome.to.use = NULL, only.gene = NULL ){

  input.df.list <- lapply(input.df.orig.list, function(x){
    temp.df <- x
    if('nrSupportGuides' %in% names(temp.df)){
      no.guides.index <- which(temp.df$nrSupportGuides == 0)
      if(length(no.guides.index) > 0){
        temp.df <- temp.df[-no.guides.index,]
      }
    }
    temp.df
  })

  chrom.df.list <- lapply(input.df.list, clean_plotting_df)

  browserStyle_plotting(chrom.df.list, score.names, label.color.list = label.color.list,
                        x.min, x.max, track.height, genome.to.use, only.gene)

}

#' @title Adjust data frames for non-targeting guides and add them to the end of the plot
#' @param input.df: data frame containing the scores and the coordinates, $chrom, $start, $end, $label
#' @return cleaned data frame
#' @export clean_plotting_df()

clean_plotting_df <- function(input.df){

  chrom.df <- c()  # data frame containing only valid chromosomes
  na.df <- NULL
  if(length(which(is.na(input.df$start))) > 0){
    chrom.df <- input.df[which(! is.na(input.df$start)),]
    na.df <- input.df[which(is.na(input.df$start)),]

    all.na.labels <- unique(na.df$label)
    temp.inter.chr.space <- round(1/8* (max(chrom.df$end) - min(chrom.df$start))) + 1
    for(lab in all.na.labels){
      na.lab.df <- na.df[na.df$label == lab,]
      na.lab.start <- max(chrom.df$end) + temp.inter.chr.space
      na.lab.stepSize <- median(sort(chrom.df$start)[2:length(chrom.df$start)] - sort(chrom.df$end)[1:(length(chrom.df$end) - 1)])  # adjust step size to median step size of targeting guides
      if(is.na(na.lab.stepSize)){
        na.lab.stepSize <- 1
      }
      na.lab.x.start <- na.lab.start + na.lab.stepSize*c(1:nrow(na.lab.df))
      na.lab.x.end <- na.lab.x.start + na.lab.stepSize - 1

      na.lab.df$start <- na.lab.x.start
      na.lab.df$end <- na.lab.x.end
      na.lab.df$chrom <- rep(chrom.df$chrom[1], nrow(na.lab.df))

      chrom.df <- rbind(chrom.df, na.lab.df)
    }

  } else {
    chrom.df <- input.df
  }

  nr.chroms <- unique(chrom.df$chrom)
  if(length(nr.chroms) > 1){
    print('Warning! Multi-chromosome plotting not established yet')
    break()
  }

  return(chrom.df)

}

#' @title plot the data frame as bed file, ignoring size.
#' @param input.df: list of data frames data frame containing the scores and the coordinates, $chrom, $start, $end, $label
#' @param plot.name: name of the track
#' @return grob object
#' @export plot_browser_bed()

plot_browser_bed <- function(input.df, plot.min, plot.max, plot.name){
  # set up the rectangles
  out.plot <- ggplot(input.df) + geom_rect(aes(xmin = start, xmax = end, ymin = 0.45, ymax = 0.55,
                                               color = as.factor(plot.name), fill = as.factor(plot.name) )) + #color = 'black', fill = 'black' )) +
    # set the range of the entire plot on the axis
    scale_x_continuous(limits=c(plot.min, plot.max)) +
    scale_y_continuous(limits=c(0.45, 0.55)) +
    scale_color_manual(guide = FALSE, values = as.factor(plot.name), #'black',
                       breaks = as.character(input.df$chrom[1]))+
    scale_fill_manual(name = plot.name,
                      values = as.factor(plot.name)) + #adjustcolor(as.factor(plot.name), alpha.f = 1)
    theme(axis.title.x=element_blank(), axis.text.x = element_blank(), axis.ticks.x=element_blank(),
          axis.text.y = element_blank(), axis.ticks.y=element_blank(),
          panel.background = element_rect(fill = "white", colour = "white", size = 0.5, linetype = "solid"),
          panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "white"),
          legend.title = element_blank() )

  out.grob <- ggplotGrob(out.plot)
  out.grob$heights[7] <- unit(0.2, 'null')
  return(out.grob)
}

#' @title plot the data frame as bedgraph
#' @param input.df: list of data frames data frame containing the scores and the coordinates, $chrom, $start, $end, $label
#' @param plot.min: minimum x axis range
#' @param plot.max: maximum x axis range
#' @param plot.name: name of the track
#' @param plot.colors: data frame containing the color to label key, $color and $label
#' @return grob object
#' @export
#' plot_browser_bedGraph()

plot_browser_bedGraph <- function(input.df, plot.min, plot.max, plot.name, plot.colors = NULL, track.height = NULL){
  # if no plot colors are specified, set up the data frame to contain simply the names of the labels
  label.colors <- plot.colors
  if(is.null(plot.colors)){
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    label.colors <- data.frame(color = col_vector[1:length(unique(input.df$label))], label = unique(input.df$label), stringsAsFactors = F)
  }

  out.plot <- ggplot(input.df)+
    # set up the rectangles of the plot
    geom_rect(aes(xmin = start, xmax = end, ymin = 0, ymax = genomeScore, color = label, fill = label )) +
    # fill the rectangels and their borders as specified
    scale_color_manual(guide = FALSE, values = as.character(label.colors$color),
                       breaks = as.character(label.colors$label))+
    scale_fill_manual(name = 'Labels',
                      values = adjustcolor(as.character(label.colors$color), alpha.f = 1)) +
    guides(fill = guide_legend(override.aes = list(alpha = 1))) + # not sure what this does
    # specifiy the size of the plot
    scale_x_continuous(limits=c(plot.min, plot.max)) +
    theme(axis.title.x=element_blank(), axis.text.x = element_blank(), axis.ticks.x=element_blank(),
          panel.background = element_rect(fill = "white", colour = "white", size = 0.5, linetype = "solid"),
          panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "white"),
          panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "white"),
          legend.title = element_blank(), axis.line.y = element_line(colour = 'black', size = 0.5, linetype = 'solid')) +
    ylab(plot.name) + geom_hline(yintercept = 0, size = 0.1)

  out.grob <- ggplotGrob(out.plot)
  if(is.null(track.height)){
    out.grob$heights[7] <- unit(3, 'null') # height arbitrarily set to 3, could be a variable..
  } else {
    out.grob$heights[7] <- unit(track.height, 'null') # height arbitrarily set to 3, could be a variable..
  }
  return(out.grob)
}

#' @title plot the data frame as bedgraph
#' @param input.df: list of data frames data frame containing the scores and the coordinates, $chrom, $start, $end, $label
#' @param plot.min: minimum x axis range
#' @param plot.max: maximum x axis range
#' @param plot.name: name of the track
#' @param plot.colors: data frame containing the color to label key, $color and $label
#' @return grob object
#' @export
#' plot_pointCloud()

plot_pointCloud <- function(input.df, plot.min, plot.max, plot.name, plot.colors = NULL, track.height = NULL){
  # if no plot colors are specified, set up the data frame to contain simply the names of the labels
  label.colors <- plot.colors
  if(is.null(plot.colors)){
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    label.colors <- data.frame(color = col_vector[1:length(unique(input.df$label))], label = unique(input.df$label), stringsAsFactors = F)
  }

  input.df$avgStart <- (input.df$start + input.df$end) / 2

  out.plot <- ggplot(input.df)+
    # set up the rectangles of the plot
    geom_point(aes(x = avgStart, y = genomeScore, color = label, fill = label )) +
    # fill the rectangels and their borders as specified
    scale_color_manual(guide = FALSE, values = as.character(label.colors$color),
                       breaks = as.character(label.colors$label))+
    scale_fill_manual(name = 'Labels',
                      values = adjustcolor(as.character(label.colors$color), alpha.f = 1)) +
    guides(fill = guide_legend(override.aes = list(alpha = 1))) + # not sure what this does
    # specifiy the size of the plot
    scale_x_continuous(limits=c(plot.min, plot.max)) +
    theme(axis.title.x=element_blank(), axis.text.x = element_blank(), axis.ticks.x=element_blank(),
          panel.background = element_rect(fill = "white", colour = "white", size = 0.5, linetype = "solid"),
          panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "white"),
          panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "white"),
          legend.title = element_blank(), axis.line.y = element_line(colour = 'black', size = 0.5, linetype = 'solid')) +
    ylab(plot.name) + geom_hline(yintercept = 0, size = 0.1)

  out.grob <- ggplotGrob(out.plot)
  if(is.null(track.height)){
    out.grob$heights[7] <- unit(3, 'null') # height arbitrarily set to 3, could be a variable..
  } else {
    out.grob$heights[7] <- unit(track.height, 'null') # height arbitrarily set to 3, could be a variable..
  }
  return(out.grob)
}


#' @title Generic wrapper function for plotting
#' @param input.exons: data frame containing info about exons: ensembl_transcript_id ensembl_gene_id exon_chrom_start exon_chrom_end   gene_biotype, external_gene_name
#' @return major isoform transcripts
#' @export extract_majorIso()

extract_majorIso <- function(input.exons, min.pos, max.pos){
  plot.range <- GRanges(seqnames = 'arbitrary', ranges = IRanges(min.pos, max.pos))
  unique.genes <- unique(input.exons$ensembl_gene_id)
  out.genes <- c()

  if(length(unique.genes) == 0){
    return(NULL)
  } else {
    for(i in 1:length(unique.genes)){
      temp.exons <- input.exons[which(input.exons$ensembl_gene_id == unique.genes[i]),]
      if(temp.exons$gene_biotype[1] == 'protein_coding'){
        major.id <- temp.exons$ensembl_transcript_id[1]
        major.exons <- temp.exons[which(temp.exons$ensembl_transcript_id == major.id),]
        major.ranges <- GRanges(seqnames = 'arbitrary', ranges = IRanges(major.exons$exon_chrom_start, major.exons$exon_chrom_end))
        # check if major isoform is actually within defined range:
        major.region.overlaps <- as.data.frame(findOverlaps(major.ranges, plot.range, type = 'any'))
        if(nrow(major.region.overlaps) > 0){
          major.exons.ordered <- major.exons[order(major.exons$exon_chrom_start),]

          out.genes <- rbind(out.genes, major.exons.ordered)
        }
      }

    }

    return(out.genes)
  }

}

#' @title Generic wrapper function for plotting
#' @param input.list: list of data frames data frame containing the scores and the coordinates, $chrom, $start, $end, $label
#' @param score.names: vector of names of the column containing the scores input.df[[input.analysis.type]] or 'bed' to specify y-axis height does not matter
#' @param label.color.list: list of data frames containing the color to label key, $color and $label
#' @param genome.to.use: default hg19
#' @param only.gene: if many genes in region, display name of only this gene
#' @return NA
#' @export browserStyle_plotting()

browserStyle_plotting <- function(input.list, score.names, label.color.list,
                                  x.min = NULL, x.max = NULL, track.height = NULL,
                                  genome.to.use = NULL, only.gene){

  # 1. need to extablish the max and min to plot
  all.start.min <- c()
  all.start.max <- c()
  if(is.null(x.min)){
    all.start.min <- min(unlist(lapply(input.list, function(x) min(x$start))))
  } else {
    all.start.min <- x.min
  }

  if(is.null(x.max)){
    all.start.max <- max(unlist(lapply(input.list, function(x) max(x$end))))
  } else {
    all.start.max <- x.max
  }

  # 2. Start adding the different plots
  list.to.plot <- list()


  # 2.1 first plot initial plot
  if(score.names[1] == 'bed'){
    list.to.plot[[1]] <- plot_browser_bed(input.list[[1]], all.start.min, all.start.max,
                                          names(input.list)[1])
  } else if(score.names[1] == 'pointCloud'){
    if(! 'genomeScore' %in% names(input.list[[1]])){
      print("cloud plotting requires 'genomeScore as score column")
      break()
    } else {
      list.to.plot[[1]] <- plot_pointCloud(input.list[[1]], all.start.min, all.start.max,
                      names(input.list)[1], label.color.list[[1]], track.height[1])
    }
  } else {
    input.list[[1]]$genomeScore <- input.list[[1]][[score.names[1]]]
    list.to.plot[[1]] <- plot_browser_bedGraph(input.list[[1]], all.start.min, all.start.max,
                                               names(input.list)[1], label.color.list[[1]], track.height[1])
  }

  # 2.2 if more ate to be plotted, add them
  if(length(input.list) > 1){
    for(i in 2:length(input.list)){
      if(score.names[i] == 'bed'){
        list.to.plot[[i]] <- plot_browser_bed(input.list[[i]], all.start.min, all.start.max,
                                              names(input.list)[i])
      } else if(score.names[i] == 'pointCloud'){
        if(! 'genomeScore' %in% names(input.list[[1]])){
          print("cloud plotting requires 'genomeScore as score column")
          break()
        } else {
          list.to.plot[[i]] <- plot_pointCloud(input.list[[i]], all.start.min, all.start.max,
                                               names(input.list)[i], label.color.list[[i]], track.height[i])
        }
      } else {
        input.list[[i]]$genomeScore <- input.list[[i]][[score.names[i]]]
        list.to.plot[[i]] <- plot_browser_bedGraph(input.list[[i]], all.start.min, all.start.max,
                                                   names(input.list)[i], label.color.list[[i]], track.height[i])
      }
    }
  }

  # 2.3 add different tracks to plot
  gtab_out <- list.to.plot[[1]]
  if(length(list.to.plot) > 1){
    for(i in 2:length(list.to.plot)){
      gtab_out <- gtable_rbind(gtab_out, list.to.plot[[i]])
    }
  }


  # 3. extract the major isoforms of all genes in this section
  # 3.1 identify chromosome and region of interest, extract all exons
  # http://bowtie-bio.sourceforge.net/recount/biomaRt.pdf
  # to add another mart:
  #   listMarts()
  #   ensembl = useMart("ENSEMBL_MART_ENSEMBL")
  #   listDatasets(ensembl)
  mart.to.use <- c()
  plot.exons <- c()
  plot.chrom <- strsplit(input.list[[1]]$chrom[1], 'chr')[[1]][2]
  if(is.null(genome.to.use)){
    mart.to.use <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org",
               path="/biomart/martservice",dataset="hsapiens_gene_ensembl")

    plot.exons <- getBM( attributes=c("ensembl_transcript_id", "ensembl_gene_id","exon_chrom_start",
                                      "exon_chrom_end", 'transcription_start_site','gene_biotype', 'external_gene_name'),
                         filters=c('chromosome_name','start','end'),
                         values=list(plot.chrom, all.start.min, all.start.max),
                         mart=mart.to.use)
  } else if(genome.to.use == 'mm9') {
    # pulling archived version: https://support.bioconductor.org/p/51418/
    mart.to.use <- useMart(host='may2012.archive.ensembl.org',
                                     biomart='ENSEMBL_MART_ENSEMBL', dataset = 'mmusculus_gene_ensembl')

    plot.exons <- getBM( attributes=c("ensembl_transcript_id", "ensembl_gene_id","exon_chrom_start",
                                      "exon_chrom_end", 'transcript_start', 'gene_biotype', 'external_gene_id'),
                         filters=c('chromosome_name','start','end'),
                         values=list(plot.chrom, all.start.min, all.start.max),
                         mart=mart.to.use)

    for(e in 1:(nrow(plot.exons)-1) ){
     if(plot.exons$exon_chrom_start[e] < plot.exons$exon_chrom_end[e + 1]){
       plot.exons$transcript_start[e] <- plot.exons$exon_chrom_start[e]
     } else {
       plot.exons$transcript_start[e] <- plot.exons$exon_chrom_end[e]
     }
    }

    #'X', 7119039, 7268372

    colnames(plot.exons) <- c("ensembl_transcript_id", "ensembl_gene_id","exon_chrom_start",
                              "exon_chrom_end", 'transcription_start_site','gene_biotype', 'external_gene_name')

  } else if(genome.to.use == 'mm10') {

    mart.to.use <- useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl")

    plot.exons <- getBM( attributes=c("ensembl_transcript_id", "ensembl_gene_id","exon_chrom_start",
                                      "exon_chrom_end", 'transcription_start_site','gene_biotype', 'external_gene_name'),
                         filters=c('chromosome_name','start','end'),
                         values=list(plot.chrom, all.start.min, all.start.max),
                         mart=mart.to.use)

    colnames(plot.exons) <- c("ensembl_transcript_id", "ensembl_gene_id","exon_chrom_start",
                              "exon_chrom_end", 'transcription_start_site','gene_biotype', 'external_gene_name')
  } else if(genome.to.use == 'hg38'){
    #hsapiens_gene_ensembl
    mart.to.use <- useMart(biomart='ensembl', dataset="hsapiens_gene_ensembl", host="uswest.ensembl.org")

    plot.exons <- getBM( attributes=c("ensembl_transcript_id", "ensembl_gene_id","exon_chrom_start",
                                      "exon_chrom_end", 'transcription_start_site','gene_biotype', 'external_gene_name'),
                         filters=c('chromosome_name','start','end'),
                         values=list(plot.chrom, all.start.min, all.start.max),
                         mart=mart.to.use)

  } else if(genome.to.use == 'none'){

    grid.newpage()
    grid.draw(gtab_out)
    return()

  } else {
    print('mart not implemented yet!')
    break()
  }

  # 3.2 identify major isoform (first one of each unique gene)
  plot.transcripts <- extract_majorIso(plot.exons, all.start.min, all.start.max)

  if(is.null(plot.transcripts)){
    grid.newpage()
    grid.draw(gtab_out)
  } else {
    # 4. draw space grid and layer overlapping sections
    exon.layer <- layer_exons(plot.transcripts, input.list[[1]]$chrom[1], all.start.min, all.start.max)
    exon.grobs <- exon_layer_toGrob(exon.layer, all.start.min, all.start.max, only.gene)
    # 5. add grobs to plot

    gtab_out <- gtable_rbind(gtab_out, exon.grobs)
    grid.newpage()
    grid.draw(gtab_out)
  }


}

#' @title layer exons such that they don't overlap
#' @param input.layer: list of lists, each sub-list is a layer and each element of a sub list is a GRanges object
#' @param all.min: minimum coordiante of plot
#' @param all.max: maximum coordinate of plot
#' @param only.gene: if many genes in region, display name of only this gene
#' @return grob object
#' @export exon_layer_toGrob()

exon_layer_toGrob <- function(input.layer, all.min, all.max, only.gene){

  gene.layer <- c()
  exon.lines <- c()
  name.layer <- c()
  arrow.layer <- c()
  line.layer <- c()

  plot.arrows <- FALSE # arrows will only be plotted if resolution is high enough

  layer.jump <- 1
  for(i in 1:(length(input.layer) - 1)){ # last layer is 'NULL'
    temp.layer <- input.layer[[i]]
    temp.exon.df <- c()
    temp.exon.line.df <- c()
    temp.gene.name.df <- c()
    temp.arrow.df <- c()
    temp.line.df <- c()
    for(l in 1:length(temp.layer)){
      temp.gene <- as.data.frame(temp.layer[[l]], stringsAsFactors = F)

      # add ending block:
      gene.end.df <- c()
      gene.reversed <- FALSE
      temp.complete.gene <- c()

      # not displaying untranslated regions
      temp.gene$minY <- rep(0.5 - i*layer.jump - 0.1, nrow(temp.gene))
      temp.gene$maxY <- rep(0.5 - i*layer.jump + 0.1, nrow(temp.gene))
      temp.complete.gene <- temp.gene

      if(temp.gene$end[1] == temp.gene$TSS[1]){
        gene.reversed <- TRUE
      } else if(temp.gene$end[nrow(temp.gene)] == temp.gene$TSS[1]){
        gene.reversed <- TRUE
      }

      # if(nrow(temp.gene) > 1){
      #   if(temp.gene$start[1] > temp.gene$start[2]){ # ending is on left side
      #     gene.end.df <- temp.gene[nrow(temp.gene),]
      #     temp.gene$start[nrow(temp.gene)] <- temp.gene$end[nrow(temp.gene)] - 10
      #     gene.reversed <- TRUE
      #   } else {# ending is on right side
      #     gene.end.df <- temp.gene[nrow(temp.gene),]
      #     temp.gene$end[nrow(temp.gene)] <- temp.gene$start[nrow(temp.gene)] + 10
      #   }
      #
      #   # add layer range to gene
      #   temp.gene$minY <- rep(0.5 - i*layer.jump - 0.1, nrow(temp.gene))
      #   temp.gene$maxY <- rep(0.5 - i*layer.jump + 0.1, nrow(temp.gene))
      #   gene.end.df$minY <- rep(0.5 - i*layer.jump - 0.05, 1)
      #   gene.end.df$maxY <- rep(0.5 - i*layer.jump + 0.05, 1)
      #   temp.complete.gene <- rbind(temp.gene, gene.end.df)
      # } else {
      #   temp.gene$minY <- rep(0.5 - i*layer.jump - 0.1, nrow(temp.gene))
      #   temp.gene$maxY <- rep(0.5 - i*layer.jump + 0.1, nrow(temp.gene))
      #   temp.complete.gene <- temp.gene
      # }

      # create exon connecting lines
      # all.min, all.max
      arrow.step.size <- (all.max - all.min) / 100
      temp.line <- c()
      temp.arrows <- c()

      temp.exon.lines <- c()
      if(nrow(temp.gene) > 1){
        if(gene.reversed){

          temp.exon.lines <- temp.gene[1:(nrow(temp.gene) - 1),]
          temp.exon.lines$start <- temp.gene$start[1:(nrow(temp.gene) - 1)]
          temp.exon.lines$end <- temp.gene$end[2:nrow(temp.gene)]
          temp.exon.lines$yCoord <- rep(0.5 - i*layer.jump, nrow(temp.exon.lines))

          temp.line <- data.frame(start = temp.exon.lines$start[1], end = temp.exon.lines$end[nrow(temp.exon.lines)],
                                yCoord = 0.5 - i*layer.jump)

          # arrow.starts <- seq(temp.exon.lines$start[1], temp.exon.lines$end[nrow(temp.exon.lines)], by = -arrow.step.size)
          arrow.starts <- seq(temp.exon.lines$end[nrow(temp.exon.lines)], temp.exon.lines$start[1], by = -arrow.step.size)
          if(arrow.step.size < abs(temp.line$end - temp.line$start) ){
            plot.arrows <- TRUE
            temp.arrows<- data.frame(start = arrow.starts[1:(length(arrow.starts) - 1)], end = arrow.starts[2:length(arrow.starts)],
                                   yCoord = rep(0.5 - i*layer.jump, length(arrow.starts) - 1) )
          }

        } else {
          temp.exon.lines <- temp.gene[1:(nrow(temp.gene) - 1),]
          temp.exon.lines$start <- temp.gene$end[1:(nrow(temp.gene) - 1)]
          temp.exon.lines$end <- temp.gene$start[2:nrow(temp.gene)]
          temp.exon.lines$yCoord <- rep(0.5 - i*layer.jump, nrow(temp.exon.lines))

          temp.line <- data.frame(start = temp.exon.lines$start[1], end = temp.exon.lines$end[nrow(temp.exon.lines)],
                                  yCoord = 0.5 - i*layer.jump)

          arrow.starts <- seq(temp.exon.lines$start[1], temp.exon.lines$end[nrow(temp.exon.lines)], by = arrow.step.size)
          if(arrow.step.size < abs(temp.line$end - temp.line$start)){
            plot.arrows <- TRUE
            temp.arrows<- data.frame(start = arrow.starts[1:(length(arrow.starts) - 1)], end = arrow.starts[2:length(arrow.starts)],
                                     yCoord = rep(0.5 - i*layer.jump, length(arrow.starts) - 1) )
          }
        }
      } else {
        temp.exon.lines <- temp.gene
        temp.exon.lines$yCoord <- rep(0.5 - i*layer.jump, nrow(temp.exon.lines))
      }

      # position the gene name
      temp.name <- c()
      if(is.null(only.gene)){
        temp.name <- data.frame(chrom = temp.gene$seqnames[1], xposition = (min(temp.gene$start) + max(temp.gene$end)) / 2,
                                yposition = 0.5 - i*layer.jump + 0.55*layer.jump, stringsAsFactors = F)
        temp.gene.name <- names(temp.layer)[l]
        if(gene.reversed){
          temp.gene.name <- paste0('<', temp.gene.name)
        } else {
          temp.gene.name <- paste0(temp.gene.name, '>')
        }
        temp.name$geneName <- temp.gene.name
      } else {
        if(names(temp.layer)[l] == only.gene){
          temp.name <- data.frame(chrom = temp.gene$seqnames[1], xposition = (min(temp.gene$start) + max(temp.gene$end)) / 2,
                                  yposition = 0.5 - i*layer.jump + 0.55*layer.jump, stringsAsFactors = F)
          temp.gene.name <- names(temp.layer)[l]
          if(gene.reversed){
            temp.gene.name <- paste0('<', temp.gene.name)
          } else {
            temp.gene.name <- paste0(temp.gene.name, '>')
          }
          temp.name$geneName <- temp.gene.name
        }
      }

      # check plot limits
      if(min( c(temp.complete.gene$start, temp.complete.gene$end) ) < all.min){

        if(gene.reversed){
          temp.complete.gene <- temp.complete.gene[-which(temp.complete.gene$end < all.min),]

          temp.line$start <- all.min + 10

          if(plot.arrows){
            temp.arrows <- temp.arrows[-which(temp.arrows$end < all.min),]
          }

        } else {
          temp.complete.gene <- temp.complete.gene[-which(temp.complete.gene$start < all.min),]

          temp.line$end <- all.min + 10

          if(plot.arrows){
            temp.arrows <- temp.arrows[-which(temp.arrows$start < all.min),]
          }
        }

      } else if(max( c(temp.complete.gene$start, temp.complete.gene$end) ) > all.max){

        if(gene.reversed){
          temp.complete.gene <- temp.complete.gene[-which(temp.complete.gene$start > all.max),]

          temp.line$start <- all.max - 10

          if(plot.arrows){
            temp.arrows <- temp.arrows[-which(temp.arrows$start > all.max),]
          }
        } else {
          temp.complete.gene <- temp.complete.gene[-which(temp.complete.gene$end > all.max),]

          temp.line$end <- all.max - 10

          if(plot.arrows){
            temp.arrows <- temp.arrows[-which(temp.arrows$end > all.max),]
          }
        }

      }
      # check plot limits
      # if(min( c(temp.complete.gene$start, temp.complete.gene$end) ) < all.min){
      #
      #   browser()
      #   temp.complete.gene <- temp.complete.gene[-which(temp.complete.gene$start < all.min),]
      #
      #   out.of.bounds.lines <- which(temp.exon.lines$start < all.min)
      #
      #   if(gene.reversed){
      #
      #     temp.exon.lines <- temp.exon.lines[-out.of.bounds.lines[2:length(out.of.bounds.lines)],]
      #     temp.exon.lines$start[nrow(temp.exon.lines)] <- all.min + 1
      #   } else {
      #     temp.exon.lines <- temp.exon.lines[-out.of.bounds.lines[1:(length(out.of.bounds.lines) - 1)],]
      #     temp.exon.lines$start[1] <- all.min + 1 # after trimming the first position of start is in bounds
      #   }
      #
      #   temp.name$xposition <- max(temp.complete.gene$start)
      #
      # } else if(max( c(temp.complete.gene$start, temp.complete.gene$end) ) > all.max){
      #   browser()
      #   if(gene.reversed){
      #     temp.complete.gene <- temp.complete.gene[-which(temp.complete.gene$start > all.max),]
      #     out.of.bounds.lines <- which(temp.exon.lines$start > all.max)
      #
      #     temp.exon.lines <- temp.exon.lines[-out.of.bounds.lines[1:(length(out.of.bounds.lines) - 1)],]
      #     temp.exon.lines$end[1] <- all.max - 1
      #   } else {
      #     temp.complete.gene <- temp.complete.gene[-which(temp.complete.gene$end > all.max),]
      #     out.of.bounds.lines <- which(temp.exon.lines$end > all.max)
      #
      #     temp.exon.lines <- temp.exon.lines[-out.of.bounds.lines[2:length(out.of.bounds.lines)],]
      #     temp.exon.lines$end[nrow(temp.exon.lines)] <- all.max - 1
      #   }
      #
      #   temp.name$xposition <- min(temp.complete.gene$end)
      #
      # }

      # add exons to overall exon df
      temp.exon.df <- rbind(temp.exon.df, temp.complete.gene)
      temp.exon.line.df <- rbind(temp.exon.line.df, temp.exon.lines)
      temp.gene.name.df <- rbind(temp.gene.name.df, temp.name)
      temp.arrow.df <- rbind(temp.arrow.df, temp.arrows)
      temp.line.df <- rbind(temp.line.df, temp.line)

    }
    gene.layer <- rbind(gene.layer, temp.exon.df)
    exon.lines <- rbind(exon.lines, temp.exon.line.df)
    name.layer <- rbind(name.layer, temp.gene.name.df)
    arrow.layer <- rbind(arrow.layer, temp.arrow.df)
    line.layer <- rbind(line.layer, temp.line.df)
  }

  # adding the exons
  out.plot <- ggplot(gene.layer) + geom_rect(aes(xmin = start, xmax = end,
                                                 ymin = minY, ymax = maxY,
                                                 color = as.factor('Genes'), fill = as.factor('Genes') )) + #color = 'black', fill = 'black'
    scale_x_continuous(limits=c(all.min, all.max)) +
    scale_y_continuous(limits=c(min(gene.layer$minY), max(gene.layer$maxY) + layer.jump)) +
    scale_color_manual(guide = FALSE, values = as.factor('Genes'), breaks = as.character(gene.layer$seqnames[1]))+
    scale_fill_manual(name = 'Genes', values = as.factor('Genes')) + # adjustcolor('black', alpha.f = 1)) +
    theme(axis.title.y=element_blank(), axis.text.y = element_blank(), axis.ticks.y=element_blank(), axis.title.x = element_blank(),
          panel.background = element_rect(fill = "white", colour = "white", size = 0.5, linetype = "solid"),
          panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "white"),
          panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "white"),
          legend.title = element_blank())

  # adding the exon lines
  if(plot.arrows){
    out.plot <- out.plot + geom_segment(data = arrow.layer,
                                        aes(x = start, xend = end, y = yCoord, yend = yCoord ), size = 0.05,
                                        arrow = arrow(angle = 75, length = unit(0.1, 'cm'))) # arrow(length = unit(0.5,"cm"), type = "closed")
  }

  out.plot <- out.plot + geom_segment(data = line.layer, aes(x = start, xend = end, y = yCoord, yend = yCoord ), size = 0.5)

  # out.plot <- out.plot + geom_segment(data = exon.lines,
  #                                     aes(x = start, xend = end, y = yCoord, yend = yCoord ), size = 0.5,
  #                                     arrow = arrow(angle = 75))

  # adding Gene names
  out.plot <- out.plot + annotate('text', x = name.layer$xposition , y = name.layer$yposition, label = name.layer$geneName)

  # conver to grob
  out.plot <- ggplotGrob(out.plot)
  out.plot$heights[7] <- unit(0.6, 'null')
  return(out.plot)
}

#' @title layer exons such that they don't overlap
#' @param input.exons: data frame containing the exon name info for each gene: ensembl_transcript_id, ensembl_gene_id, exon_chrom_start, exon_chrom_end, transcription_start_site, gene_biotype, external_gene_name
#' @param input.chrom: chromsome to plot
#' @param all.min: minimum coordiante of plot
#' @param all.max: maximum coordinate of plot
#' @return list of lists, each sub-list is a layer and each element of a sub list is a GRanges object
#' @export layer_exons()

layer_exons <- function(input.exons, input.chrom, all.min, all.max){

  unique.genes <- unique(input.exons$external_gene_name)

  # each list element contains a list of GRanges elements of genes:
  # ex: all.gene.layers[[1]] = $IL2RA, $RBM17
  #     all.gene.layers[[2]] = $PFKFB3
  all.gene.layers <- list('NULL')

  for(i in 1:length(unique.genes)){
    temp.gene <- input.exons[which(input.exons$external_gene_name == unique.genes[i]),]
    temp.gene.range <- GRanges(seqnames = rep(input.chrom, nrow(temp.gene)),
                               ranges = IRanges(temp.gene$exon_chrom_start, temp.gene$exon_chrom_end),
                               TSS = temp.gene$transcription_start_site)
    # for every layer, check if there is an overlap with existing genes:
    all.gene.layers <- check_gene_layer(all.gene.layers, temp.gene.range, unique.genes[i])
  }

  return(all.gene.layers)
}

#' @title layer exons such that they don't overlap
#' @param input.gene.list: list, each list element contains a list of GRanges elements of genes:
#' @param input.gene: gene to add to input.gene.list in GRanges format
#' @param input.gene.name: name of gene to add
#' @return GRanges object of the gene
#' @export check_gene_layer()

check_gene_layer <- function(input.gene.list, input.gene, input.gene.name){
  for(i in 1:length(input.gene.list)){
    # if the layer is empty: add Granges object to list
    if(input.gene.list[[i]][1] == 'NULL'){
      input.gene.list[[i]] <- list()
      input.gene.list[[i]][[input.gene.name]] <- input.gene
      input.gene.list[[i + 1]] <- 'NULL' # assignment of null for 2 or higher does not seem to work
      return(input.gene.list)
    } else { # check if there is overlap with existing genes

      temp.layer.genes <- input.gene.list[[i]]

      # vector of 0s and 1s. Latter indicated an overlap and gene has to be passed to next layer
      has.overlaps <- lapply(temp.layer.genes, function(x){
        temp.overlaps <- as.data.frame(findOverlaps(x, input.gene, type = 'any'))
        if(nrow(temp.overlaps) > 0){
          return(1)
        } else {
          return(0)
        }
      })

      if(sum(unlist(has.overlaps)) > 0){

      } else { # add gene to the layer
        input.gene.list[[i]][[input.gene.name]] <- input.gene
        return(input.gene.list)
      }
      #need to initiate next layer as null if new layer is started

    }
  }
  print('something went wrong with gene layer addition')
}

#' @title layer exons such that they don't overlap
#' @param input.score.list: list, each list element contains the name of the bedgraph to be plotted and a data frame with: $chrom, $start, $end, $formatScores
#' @param bg.name: name to add to bedgraph
#' @return bedgraph file
#' @export create_bedgraphs()

create_bedgraphs <- function(input.score.list, bg.name){
  score.names <- names(input.score.list)
  nr.bg <- length(score.names)  # create unique colors for each bedgraph to plot
  bedgraph.colors <- col2rgb(rainbow(nr.bg,s = 1, v = 1, start = 0, end =  max(1, nr.bg - 1)/nr.bg, alpha = 1))
  for(i in 1:length(score.names)){
    intit.temp.score.df <- input.score.list[[i]]
    temp.score.df <- intit.temp.score.df[which(! is.na(intit.temp.score.df$start)),]

    if(nrow(temp.score.df) == 0){
      print('bedraph plotting failed. No non-targeting guides to remove. Fix code to account for that')
      break()
    }

    if(! 'formatScores' %in% names(temp.score.df)){
      if('genomeScore' %in% names(temp.score.df)){
        temp.score.df$formatScores <- temp.score.df$genomeScore
      } else if('guideScore' %in% names(temp.score.df)){
        temp.score.df$formatScores <- temp.score.df$guideScore
      } else if('score' %in% names(temp.score.df)){
        temp.score.df$formatScores <- temp.score.df$score
      } else{
        print('failed to identify score column')
        break()
      }
    }

    temp.score.name <- score.names[i]
    print(paste0('Creating Bedgraph for: ', temp.score.name))
    temp.bg.color <- paste(bedgraph.colors[1,i], bedgraph.colors[2,i], bedgraph.colors[3,i], sep = ',')
    first.chrom <- temp.score.df[which(temp.score.df$chrom == temp.score.df$chrom[1]),]
    temp.header1 <- paste0('browser position ', first.chrom$chrom[1], ':', min(first.chrom$start, na.rm = T),
                           '-', max(first.chrom$end, na.rm = T), '\n')
    temp.header2 <- paste0("track type=bedGraph name='", bg.name, temp.score.name, "_bg' description='",
                           bg.name, temp.score.name, "_bg' visibility=full color=", temp.bg.color, '\n')

    temp.score.gr <- GRanges(seqnames = temp.score.df$chrom, ranges = IRanges(temp.score.df$start, temp.score.df$end),
                             score = temp.score.df$formatScores)
    temp.score.gr <- sortSeqlevels(temp.score.gr)
    temp.score.gr <- sort(temp.score.gr)
    df.temp.score <- as.data.frame(temp.score.gr)
    df.temp.score.final <- df.temp.score[,c(1,2,3,6)]

    tmp.file <- paste0(bg.name, '_',temp.score.name, ".bedgraph")
    cat(temp.header1, file = tmp.file)
    cat(temp.header2, file = tmp.file, append = TRUE)

    write.table(df.temp.score.final, file = tmp.file, append = TRUE, row.names = F, col.names = F, quote = F, sep = '\t')

    temp.method.name <- strsplit(score.names[i],'_')[[1]]

  }
}
