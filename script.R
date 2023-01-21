## setwd("absolute path of a directory containing the input data files")
options(warn=-1)
library(circlize)
library(RColorBrewer)
library(GenomicRanges)
library(data.table)
library(RLumShiny)
library(grDevices)

plotcircos <- function(x, color, height, plotTypes, units, rotation, gap.width, labeltextchr, poslabelschr, heightlabelschr, marginlabelschr, data.CN){
  circos.par("start.degree"=90-rotation, "gap.degree"=gap.width, cell.padding=c(0,0,0,0), track.margin=c(0,0))
  circos.genomicInitialize.new(x,plotType=plotTypes,unit=units)
  if(!is.null(data.CN) && ncol(data.CN)==4 && labeltextchr==1 && poslabelschr=="outer"){
    circos.genomicLabels(data.CN, labels.column=4, connection_height=heightlabelschr, track.margin=c(0.01,marginlabelschr), side="outside")
  }		
  circos.genomicTrackPlotRegion(ylim = c(0, 1),bg.col = color, bg.border = NA, track.height = height)	
  if(!is.null(data.CN) && ncol(data.CN)==4 && labeltextchr==1 && poslabelschr=="inner"){
    circos.genomicLabels(data.CN, labels.column=4, connection_height=heightlabelschr, track.margin=c(0.01,marginlabelschr), side="inside")
  }		
}

plotcircos.notrack <- function(x, plotTypes, units, rotation, gap.width, data.CN, labeltextchr, poslabelschr, heightlabelschr, marginlabelschr){
  circos.par("start.degree"=90-rotation, "gap.degree"=gap.width, cell.padding=c(0,0,0,0), track.margin=c(0,0))
  circos.genomicInitialize.new(x,plotType=plotTypes,unit=units)
  if(!is.null(data.CN) && ncol(data.CN)==4 && labeltextchr==1 && poslabelschr=="inner"){
    circos.genomicLabels(data.CN, labels.column = 4, connection_height = heightlabelschr, track.margin = c(0.01,marginlabelschr), side = "inside")
  }
}

plotcircos.font <- function(x, color, height, plotTypes, units, rotation, gap.width, cexLabel, labeltextchr, poslabelschr, heightlabelschr, marginlabelschr, data.CN){
  circos.par("start.degree"=90-rotation, "gap.degree"=gap.width, cell.padding=c(0,0,0,0), track.margin=c(0,0))
  circos.genomicInitialize.new.font(x, plotType=plotTypes, unit=units, cexlabel=cexLabel)
  if(!is.null(data.CN) && ncol(data.CN)==4 && labeltextchr==1 && poslabelschr=="outer"){
    circos.genomicLabels(data.CN, labels.column = 4, connection_height = heightlabelschr, track.margin = c(0.01,marginlabelschr), side = "outside")
  }	
  circos.genomicTrackPlotRegion(ylim = c(0, 1),bg.col = color, bg.border = NA, track.height = height)
  if(!is.null(data.CN) && ncol(data.CN)==4 && labeltextchr==1 && poslabelschr=="inner"){
    circos.genomicLabels(data.CN, labels.column = 4, connection_height = heightlabelschr, track.margin = c(0.01,marginlabelschr), side = "inside")
  }	
}

plotcircos.notrack.font <- function(x, plotTypes, units, rotation, gap.width, cexLabel, data.CN, labeltextchr, poslabelschr, heightlabelschr, marginlabelschr){  
  circos.par("start.degree"=90-rotation, "gap.degree"=gap.width, cell.padding=c(0,0,0,0), track.margin=c(0,0))
  circos.genomicInitialize.new.font(x, plotType=plotTypes, unit=units, cexlabel=cexLabel)
  if(!is.null(data.CN) && ncol(data.CN)==4 && labeltextchr==1 && poslabelschr=="inner"){
    circos.genomicLabels(data.CN, labels.column = 4, connection_height = heightlabelschr, track.margin = c(0.01,marginlabelschr), side = "inside")
  }	
}

plotcircos.cyto <- function(x, height, plotTypes, units, rotation, gap.width, labeltextchr, poslabelschr, heightlabelschr, marginlabelschr, data.CN){ 
  circos.par("start.degree"=90-rotation, "gap.degree"=gap.width, cell.padding=c(0,0,0,0), track.margin=c(0,0))
  circos.genomicInitialize.new(x, plotType = plotTypes, unit=units)
  if(!is.null(data.CN) && ncol(data.CN)==4 && labeltextchr==1 && poslabelschr=="outer"){
    circos.genomicLabels(data.CN, labels.column = 4, connection_height = heightlabelschr, track.margin = c(0.01,marginlabelschr), side = "outside")
  }	
  circos.genomicTrackPlotRegion(x, ylim = c(0, 1), bg.border = NA, 
                                track.height = height, panel.fun = function(region, value, ...){
                                  col = cytoband.col(value[[2]])
                                  circos.genomicRect(region, value, ybottom = 0, 
                                                     ytop = 1, col = col, border = NA, ...)
                                  xlim = get.cell.meta.data("xlim")
                                  circos.rect(xlim[1], 0, xlim[2], 1, border = "black")
                                }, cell.padding = c(0, 0, 0, 0))  
  if(!is.null(data.CN) && ncol(data.CN)==4 && labeltextchr==1 && poslabelschr=="inner"){
    circos.genomicLabels(data.CN, labels.column = 4, connection_height = heightlabelschr, track.margin = c(0.01,marginlabelschr), side = "inside")
  }	
}

plotcircos.cyto.font <- function(x, height, plotTypes, units, rotation, gap.width, cexLabel, labeltextchr, poslabelschr, heightlabelschr, marginlabelschr, data.CN){
  circos.par("start.degree"=90-rotation, "gap.degree"=gap.width, cell.padding=c(0,0,0,0), track.margin=c(0,0))
  circos.genomicInitialize.new.font(x, plotType=plotTypes, unit=units, cexlabel=cexLabel)  
  if(!is.null(data.CN) && ncol(data.CN)==4 && labeltextchr==1 && poslabelschr=="outer"){
    circos.genomicLabels(data.CN, labels.column = 4, connection_height = heightlabelschr, track.margin = c(0.01,marginlabelschr), side = "outside")
  }	
  circos.genomicTrackPlotRegion(x, ylim = c(0, 1), bg.border = NA, 
                                track.height = height, panel.fun = function(region, value, ...){
                                  col = cytoband.col(value[[2]])
                                  circos.genomicRect(region, value, ybottom = 0, 
                                                     ytop = 1, col = col, border = NA, ...)
                                  xlim = get.cell.meta.data("xlim")
                                  circos.rect(xlim[1], 0, xlim[2], 1, border = "black")
                                }, cell.padding = c(0, 0, 0, 0))
  if(!is.null(data.CN) && ncol(data.CN)==4 && labeltextchr==1 && poslabelschr=="inner"){
    circos.genomicLabels(data.CN, labels.column = 4, connection_height = heightlabelschr, track.margin = c(0.01,marginlabelschr), side = "inside")
  }	
}

circos.genomicInitialize.new <- 
  function (data, sector.names = NULL, major.by = NULL, unit = "", plotType, tickLabelsStartFromZero = TRUE, track.height = 0.05, 
            ...) 
  {
    if(is.factor(data[[1]])){
      fa = levels(data[[1]])
    }
    else {
      fa = unique(data[[1]])
    }
    if(!is.null(sector.names)){
      if(length(sector.names) != length(fa)){
        stop("length of `sector.names` and length of sectors differ.")
      }
    }
    else {
      sector.names = fa
    }
    names(sector.names) = fa
    x1 = tapply(data[[2]], data[[1]], min)[fa]
    x2 = tapply(data[[3]], data[[1]], max)[fa]
    op = circos.par("cell.padding")
    ow = circos.par("points.overflow.warning")
    circos.par(cell.padding = c(0, 0, 0, 0), points.overflow.warning = FALSE)
    circos.initialize(factor(fa, levels = fa), xlim = cbind(x1, 
                                                            x2), ...)
    if(any(plotType %in% c("axis", "labels"))){
      circos.genomicTrackPlotRegion(data, ylim = c(0, 1), bg.border = NA, 
                                    track.height = track.height, panel.fun = function(region, 
                                                                                      value, ...){
                                      sector.index = get.cell.meta.data("sector.index")
                                      xlim = get.cell.meta.data("xlim")
                                      if(tickLabelsStartFromZero){
                                        offset = xlim[1]
                                        if(is.null(major.by)){
                                          xlim = get.cell.meta.data("xlim")
                                          major.by = .default.major.by()
                                        }
                                        major.at = seq(xlim[1], xlim[2], by = major.by)
                                        major.at = c(major.at, major.at[length(major.at)] + 
                                                       major.by)
                                        if(major.by > 1e+06){
                                          major.tick.labels = paste((major.at - offset)/1e+06, 
                                                                    "MB", sep = "")
                                        }
                                        else if(major.by > 1000){
                                          major.tick.labels = paste((major.at - offset)/1000, 
                                                                    "KB", sep = "")
                                        }
                                        else {
                                          major.tick.labels = paste((major.at - offset), 
                                                                    "bp", sep = "")
                                        }
                                      }
                                      else {
                                        if(is.null(major.by)){
                                          xlim = get.cell.meta.data("xlim")
                                          major.by = .default.major.by()
                                        }
                                        major.at = seq(floor(xlim[1]/major.by) * major.by, 
                                                       xlim[2], by = major.by)
                                        major.at = c(major.at, major.at[length(major.at)] + 
                                                       major.by)
                                        if(major.by > 1e+06){
                                          major.tick.labels = paste(major.at/1e+06, 
                                                                    "MB", sep = "")
                                        }
                                        else if(major.by > 1000){
                                          major.tick.labels = paste(major.at/1000, 
                                                                    "KB", sep = "")
                                        }
                                        else {
                                          major.tick.labels = paste(major.at, "bp", 
                                                                    sep = "")
                                        }
                                      }
                                      
                                      if(unit==""){ major.tick.labels <- gsub("[mkbp]","",major.tick.labels,ignore.case = T)}
                                      
                                      if(all(c("axis", "labels") %in% plotType)){
                                        circos.axis(h = 0, major.at = major.at, labels = major.tick.labels, 
                                                    labels.cex = 0.49 * par("cex"), labels.facing = "clockwise", 
                                                    major.tick.percentage = 0.2)
                                        circos.text(mean(xlim), 1.2, labels = sector.names[sector.index], 
                                                    cex = par("cex")-0.1, adj = c(0.5, -0.1*par("cex")*6-(par("cex")-1)*3), niceFacing = TRUE)
                                      }
                                      else if("labels" %in% plotType){
                                        circos.text(mean(xlim), 0, labels = sector.names[sector.index], 
                                                    cex = par("cex")-0.1, adj = c(0.5, -0.1*par("cex")*6-(par("cex")-1)*3), niceFacing = TRUE)
                                      }
                                      else if("axis" %in% plotType){
                                        circos.axis(h = 0, major.at = major.at, labels = major.tick.labels, 
                                                    labels.cex = 0.49 * par("cex"), labels.facing = "clockwise", 
                                                    major.tick.percentage = 0.2)
                                      }
                                    })
    }
    circos.par(cell.padding = op, points.overflow.warning = ow)
    return(invisible(NULL))
  }
  
circos.genomicInitialize.new.font <- 
  function (data, sector.names = NULL, major.by = NULL, unit = "", plotType, tickLabelsStartFromZero = TRUE, track.height = 0.05, cexlabel, 
            ...) 
  {
    if(is.factor(data[[1]])){
      fa = levels(data[[1]])
    }
    else {
      fa = unique(data[[1]])
    }
    if(!is.null(sector.names)){
      if(length(sector.names) != length(fa)){
        stop("length of `sector.names` and length of sectors differ.")
      }
    }
    else {
      sector.names = fa
    }
    names(sector.names) = fa
    x1 = tapply(data[[2]], data[[1]], min)[fa]
    x2 = tapply(data[[3]], data[[1]], max)[fa]
    op = circos.par("cell.padding")
    ow = circos.par("points.overflow.warning")
    circos.par(cell.padding = c(0, 0, 0, 0), points.overflow.warning = FALSE)
    circos.initialize(factor(fa, levels = fa), xlim = cbind(x1, 
                                                            x2), ...)
    if(any(plotType %in% c("axis", "labels"))){
      circos.genomicTrackPlotRegion(data, ylim = c(0, 1), bg.border = NA, 
                                    track.height = track.height, panel.fun = function(region, 
                                                                                      value, ...){
                                      sector.index = get.cell.meta.data("sector.index")
                                      xlim = get.cell.meta.data("xlim")
                                      if(tickLabelsStartFromZero){
                                        offset = xlim[1]
                                        if(is.null(major.by)){
                                          xlim = get.cell.meta.data("xlim")
                                          major.by = .default.major.by()
                                        }
                                        major.at = seq(xlim[1], xlim[2], by = major.by)
                                        major.at = c(major.at, major.at[length(major.at)] + 
                                                       major.by)
                                        if(major.by > 1e+06){
                                          major.tick.labels = paste((major.at - offset)/1e+06, 
                                                                    "MB", sep = "")
                                        }
                                        else if(major.by > 1000){
                                          major.tick.labels = paste((major.at - offset)/1000, 
                                                                    "KB", sep = "")
                                        }
                                        else {
                                          major.tick.labels = paste((major.at - offset), 
                                                                    "bp", sep = "")
                                        }
                                      }
                                      else {
                                        if(is.null(major.by)){
                                          xlim = get.cell.meta.data("xlim")
                                          major.by = .default.major.by()
                                        }
                                        major.at = seq(floor(xlim[1]/major.by) * major.by, 
                                                       xlim[2], by = major.by)
                                        major.at = c(major.at, major.at[length(major.at)] + 
                                                       major.by)
                                        if(major.by > 1e+06){
                                          major.tick.labels = paste(major.at/1e+06, 
                                                                    "MB", sep = "")
                                        }
                                        else if(major.by > 1000){
                                          major.tick.labels = paste(major.at/1000, 
                                                                    "KB", sep = "")
                                        }
                                        else {
                                          major.tick.labels = paste(major.at, "bp", 
                                                                    sep = "")
                                        }
                                      }
                                      
                                      if(unit==""){ major.tick.labels <- gsub("[mkbp]","",major.tick.labels,ignore.case = T)}
									  
                                      if(all(c("axis", "labels") %in% plotType)){
                                        circos.axis(h = 0, major.at = major.at, labels = major.tick.labels, 
                                                    labels.cex = 0.49 * cexlabel, labels.facing = "clockwise", 
                                                    major.tick.percentage = 0.2)
                                        circos.text(mean(xlim), 1.2, labels = sector.names[sector.index], 
                                                    cex = cexlabel, adj = c(0.5, -0.1*cexlabel*6-(cexlabel-1)*3), niceFacing = TRUE)
                                      }
                                      else if("labels" %in% plotType){
                                        circos.text(mean(xlim), 0, labels = sector.names[sector.index], 
                                                    cex = cexlabel, adj = c(0.5, -0.1*cexlabel*6-(cexlabel-1)*3), niceFacing = TRUE)
                                      }
                                      else if("axis" %in% plotType){
                                        circos.axis(h = 0, major.at = major.at, labels = major.tick.labels, 
                                                    labels.cex = 0.49 * cexlabel, labels.facing = "clockwise", 
                                                    major.tick.percentage = 0.2)
                                      }
                                    })
    }
    circos.par(cell.padding = op, points.overflow.warning = ow)
    return(invisible(NULL))
  }
  
.default.major.by = function(sector.index = get.cell.meta.data("sector.index"),
	track.index = get.cell.meta.data("track.index")){
	d = circos.par("major.by.degree")
	cell.start.degre = get.cell.meta.data("cell.start.degree", sector.index, track.index)
	tm = reverse.circlize(c(cell.start.degre, cell.start.degre-d), rep(get.cell.meta.data("cell.bottom.radius", sector.index = sector.index, track.index = track.index), 2))
	major.by = abs(tm[1, 1] - tm[2, 1])
	digits = as.numeric(gsub("^.*e([+-]\\d+)$", "\\1", sprintf("%e", major.by)))
	major.by = round(major.by, digits = -1*digits)
	return(major.by)
}

get_most_inside_radius = function() {
	tracks = get.all.track.index()
	if(length(tracks) == 0) {
	   1
	}else{
	   n = length(tracks)
	   get.cell.meta.data("cell.bottom.radius", track.index = tracks[n]) - get.cell.meta.data("track.margin", track.index = tracks[n])[1] - circos.par("track.margin")[2]
	}
}

data.C.name <- "example_data_chromosome_general.csv"
data.C <- data.frame(fread(data.C.name),stringsAsFactors=F)
data.C[,2] <- as.numeric(data.C[,2])
data.C[,3] <- as.numeric(data.C[,3])
data.T.file <- c("example_data_point.csv","example_data_heatmap.csv","example_data_line.csv","example_data_links_gradual_color.csv")
data.T <- lapply(1:length(data.T.file),function(x){
		  if(!is.null(data.T.file[x])){
		  data.frame(fread(data.T.file[x]),stringsAsFactors=F)
		  }
		  })
data.CN.name <- "example_data_gene_label.csv"
data.CN <- data.frame(fread(data.CN.name),stringsAsFactors=F)
data.CN[,2] <- as.numeric(data.CN[,2])
  data.CN[,3] <- as.numeric(data.CN[,3])
data.N.file <- c("example_data_gene_label.csv","","","","","","","","","")
uploadtrack <- c(2,2,2,1,1,1,1,1,1,1)
data.N <- lapply(1:10,function(x){
			 if(uploadtrack[x] == 2 && nchar(data.N.file[x])>0){	  
		     data.frame(fread(data.N.file[x]),stringsAsFactors=F)
			 }
			 })
trackindx <- c(1,2,3)
data.N <- data.N[trackindx]
data.L <- data.frame(fread("example_data_links.csv"),stringsAsFactors=F)
data.L1 <- data.L[,1:3]
		     data.L2 <- data.L[,4:6]
		     data.L1[,2] <- as.numeric(data.L1[,2])
		     data.L1[,3] <- as.numeric(data.L1[,3])
		     data.L2[,2] <- as.numeric(data.L2[,2])
		     data.L2[,3] <- as.numeric(data.L2[,3])	  
		     data.L1$num <- 1:nrow(data.L1)
             data.L2$num <- 1:nrow(data.L2)
		     rownames(data.L1) <- data.L1$num
		     rownames(data.L2) <- data.L2$num
hltdataLinks <- "1,40000000,120000000,navy
6,40000000,80000000,red"
tmpL <- matrix(strsplit(hltdataLinks, "\n")[[1]])
		   colnamesL <- c("chr","start","end","color")
		   datL <- matrix(0, length(tmpL), length(colnamesL))
		   colnames(datL) <- colnamesL
		   for(l in 1:length(tmpL)){
		      rowL <- strsplit(tmpL[l], ",")[[1]]
                  if(length(rowL)==4){                                        
                    datL[l,] <- rowL
                  }
		   }
		   datL <- data.frame(datL,stringsAsFactors=F)
		   datL$start <- as.numeric(datL$start)
		   datL$end <- as.numeric(datL$end)
		   datL$color <- datL$color
		   queryL <- GRanges(seqnames = datL$chr,ranges=IRanges(start=datL$start,end=datL$end),seqinfo=NULL)
		   subj1 <- GRanges(seqnames = data.L1[,1],ranges=IRanges(start=data.L1[,2],end=data.L1[,3]),seqinfo=NULL)
		   subj2 <- GRanges(seqnames = data.L2[,1],ranges=IRanges(start=data.L2[,2],end=data.L2[,3]),seqinfo=NULL)
           indx1 <- findOverlaps(queryL,subj1)
           indx1 <- data.frame(indx1,stringsAsFactors=F)
           indx1$queryHits <- as.numeric(indx1$queryHits)
           indx1$subjectHits <- as.numeric(indx1$subjectHits)
           hltregion1 <- data.L1[indx1$subjectHits,]
           data.LL1 <- data.L1
           hltregion1$color <- datL$color[indx1[,1]]
           indx2 <- findOverlaps(queryL,subj2)
           indx2 <- data.frame(indx2,stringsAsFactors=F)
           indx2$queryHits <- as.numeric(indx2$queryHits)
           indx2$subjectHits <- as.numeric(indx2$subjectHits)
           hltregion2 <- data.L2[indx2$subjectHits,]
           data.LL2 <- data.L2
           hltregion2$color <- datL$color[indx2[,1]]
for(i in 1:length(data.T.file)){
  assign(paste("hltdata",i,sep=""),"")
}
hltdata1 <- "1,1,90000000,blue
6,2,90000000,black
9,1,90000000,cyan"
hltdata2 <- ""
hltdata3 <- ""
hltregion.List <- list()
if(!is.null(data.T)){
			for(k in 1:length(data.T)){
			data.TT <- data.T[[k]]
			hltregion.List[[k]] <- ""
if(nchar(get(paste("hltdata",k,sep="")))>0){
tmp <- matrix(strsplit(get(paste("hltdata",k,sep="")), "\n")[[1]])
            myColnames <- c("chr","start","end","color")
            data <- matrix(0, length(tmp), length(myColnames))
            colnames(data) <- myColnames
            for(p in 1:length(tmp)){
                 myRow <- strsplit(tmp[p], ",")[[1]]
                  if(length(myRow)==4){                                        
                    data[p,] <- myRow
                  }
               }
            data <- data.frame(data,stringsAsFactors=F)
            data$start <- as.numeric(data$start)
            data$end <- as.numeric(data$end)
			query <- GRanges(seqnames = data$chr,ranges=IRanges(start=data$start,end=data$end),seqinfo=NULL)
            subj <- GRanges(seqnames = data.TT[,1],ranges=IRanges(start=data.TT[,2],end=data.TT[,3]),seqinfo=NULL) 
            indx <- findOverlaps(query,subj)
            indx <- data.frame(indx,stringsAsFactors=F)
			indx$queryHits <- as.numeric(indx$queryHits)
			indx$subjectHits <- as.numeric(indx$subjectHits)
            hltregion <- data.TT[indx$subjectHits,]
			hltregion$color <- data$color[indx[,1]]
			hltregion$id <- paste(hltregion[,1],hltregion[,2],hltregion[,3],sep="")
			hltregion.List[[k]] <- hltregion
			}
			}
			}

pdf("shinyCircos.pdf", width=750/72, height=750/72)
## svg("shinyCircos.svg", width=750/72, height=750/72)
fontSize <- 1.1
par(oma=c(0,0,0,0), mar=c(9,0.5,1,9.5), xpd=TRUE, cex=fontSize-0.05)
trackChr <- "track"
plotTypes <- "labels"
plotTypes <- "axis"
unitChr <- "unit"
rotation <- 0.5
gap.width <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
labeltextchr <- 1
poslabelschr <- "outer"
heightlabelschr <- 0.06
marginlabelschr <- 0.02
colorChr <- c("gold","gold","gold","gold","gold","gold","gold","gold","gold","gold","gold","gold","gold","gold","gold","gold","gold","gold","gold","gold","gold","gold","gold","gold")
heightChr <- 0.05
plotcircos(data.C, height=heightChr, color=colorChr, plotTypes=plotTypes, units=unitChr, rotation=rotation, gap.width=gap.width, labeltextchr=labeltextchr, poslabelschr=poslabelschr, heightlabelschr=heightlabelschr, marginlabelschr=marginlabelschr, data.CN=data.CN)
takindx <- 1
takindx <- takindx+2
typeTrack <- c("point","heatmap","line")
i <- 1
data.TT <- data.T[[i]]
	tktype <- typeTrack[i]
	data.TT[,2] <- as.numeric(data.TT[,2])
	data.TT[,3] <- as.numeric(data.TT[,3])
	data.NN <- data.N[[i]]
	data.TT$num <- 1:nrow(data.TT)
data.TTC <- NULL
coltypeTrack <- 2
tkcolor <- c("orange","blue")
data.TT$num <- NULL
tkbgcol <- c("grey98","grey98","grey98","grey98","grey98","grey98","grey98","grey98","grey98","grey98","grey98","grey98","grey98","grey98","grey98","grey98","grey98","grey98","grey98","grey98","grey98","grey98","grey98","grey98")
tkmargin <- 0.01
tkheight <- 0.1
tklinecoord <- c(0.25,0.75)
tklinecolor <- c("grey","grey")
hmapcols <- c("blue","white","red")
tkborder <- ""
innergap <- 0.5
tkbordercol <- NA
tkbardir <- 1
tkrectcol <- 1
selrectcol <- 1
rectcols <- c("#EDEDFD","#6969F5","#00008B")
tktransparency <- 1
tkcolor <- c("#FFA500FF","#0000FFFF")
data.TTT <- data.T[[i]]
	data.TTT$id <- paste(data.TTT[,1],data.TTT[,2],data.TTT[,3],sep="")
	data.TTT$num <- 1:nrow(data.TTT)
transparencyHlt <- c(1,1,1)
lkmargin <- 0
tkborder <- NA
columns <- c(4)
data.TT[,ncol(data.TT)] <- as.numeric(data.TT[,ncol(data.TT)])
			circos.genomicTrackPlotRegion(data.TT, track.height = tkheight, track.margin = c(lkmargin,tkmargin),
                              bg.col = tkbgcol, bg.border = tkborder, panel.fun = function(region,value,...){
							  if(nchar(tklinecolor[1])!=0){               
					            xlim <- get.cell.meta.data("xlim")
                                ylim <- get.cell.meta.data("ylim")
					            for(k in 1:length(tklinecoord)){
                                y1 <- as.numeric(quantile(ylim,probs=tklinecoord[k]))
                                circos.lines(x=xlim,y=c(y1,y1), col=tklinecolor[k], lwd=0.1)
					            }
				              }
							  if(!("cex" %in% colnames(data.T[[i]])) && !("pch" %in% colnames(data.T[[i]])) && ((coltypeTrack==1 && !("color" %in% colnames(data.T[[i]]))) | coltypeTrack==2)){
							    if(length(columns)==1){
							       tkcolor <- tkcolor[1]
							    }else{
								   tkcolor <- c(tkcolor,rep("grey",length(columns)))
								   tkcolor <- tkcolor[1:length(columns)]
								}
							  circos.genomicPoints(region, value, numeric.column=columns-3, col=tkcolor, cex=0.6, pch=16, ...)
							  }
                             })
assign("hltregion",hltregion.List[[i]])
				  hlttransparency <- transparencyHlt[i]
				  hltregion$color <- adjustcolor(hltregion$color, alpha.f = hlttransparency)
			      hltregion$color <- gsub("0x","#", hltregion$color)
                  chrr <- unique(hltregion[,1])
                  lapply(chrr, function(x){
				  datt <- hltregion[hltregion[,1] %in% x,]
                  trackk <- data.TTT[data.TTT[,1] %in% x,]
                  trackk <- trackk[!trackk$id %in% datt$id,]
                  col <- unique(datt$color)
				  if(trackChr=="track"){
				      circos.updatePlotRegion(sector.index = x, track.index=takindx+2, bg.col = tkbgcol[which(unique(data.C[,1])==x)], bg.border = tkborder)
				  }else{
				      circos.updatePlotRegion(sector.index = x, track.index=takindx+1, bg.col = tkbgcol[which(unique(data.C[,1])==x)], bg.border = tkborder)
				  }
				  if(nchar(tklinecolor[1])!=0){               
					xlim <- get.cell.meta.data("xlim")
                    ylim <- get.cell.meta.data("ylim")
					for(k in 1:length(tklinecoord)){
                    y1 <- as.numeric(quantile(ylim,probs=tklinecoord[k]))
                    circos.lines(x=xlim,y=c(y1,y1), col=tklinecolor[k], lwd=0.1)
					}
				  }
                  if(!"pch" %in% colnames(data.TTT)){                  
                     dattt_pch <- 16
                     trackk_pch <- 16
                  }else{
                     trackk_pch <- trackk$pch
                  }
				  if(!"cex" %in% colnames(data.TTT)){                  
                     dattt_cex <- 0.6
					 trackk_cex <- 0.6
                  }else{
				     trackk_cex <- trackk$cex
				  }
                  lapply(col, function(m){
                       dattt <- datt[datt$color %in% m,]
					   if("pch" %in% colnames(data.TTT)){
                          dattt_pch <- dattt$pch
                       }
					   if("cex" %in% colnames(data.TTT)){
                          dattt_cex <- dattt$cex
                       }					   
                       circos.points((dattt[,2]+dattt[,3])/2,dattt[,4], col=m, cex=dattt_cex, pch=dattt_pch)
                  })
				  if("color" %in% colnames(trackk) && coltypeTrack!=2){
				     data.TTC$id <- paste(data.TTC[,1],data.TTC[,2],data.TTC[,3],sep="")			  
				     trackkk <- data.TTC[data.TTC$id %in% trackk$id,]
				     trackkk <- trackkk[trackkk[,1] %in% x,]
				     lapply(unique(trackkk$cols),function(f){
				        trackkkk <- trackkk[trackkk$cols %in% f,]
				        if("cex" %in% colnames(trackkkk)){
				           trackk_cex <- trackkkk$cex
				        }
				        if("pch" %in% colnames(trackkkk)){
				           trackk_pch <- trackkkk$pch
				        }
                        circos.points((trackkkk[,2]+trackkkk[,3])/2,trackkkk[,4], col=adjustcolor(f,alpha.f = tktransparency), cex=trackk_cex, pch=trackk_pch)				  
				     })
				  }else{
                     circos.points((trackk[,2]+trackk[,3])/2,trackk[,4], col=tkcolor[1], cex=trackk_cex, pch=trackk_pch)
				  }
                  })
poslabels <- c("outer","outer","outer")
if(poslabels[i]=="inner"){
			    takindx <- takindx+3
			}else{
			    takindx <- takindx+1
			}
i <- 2
data.TT <- data.T[[i]]
	tktype <- typeTrack[i]
	data.TT[,2] <- as.numeric(data.TT[,2])
	data.TT[,3] <- as.numeric(data.TT[,3])
	data.NN <- data.N[[i]]
	data.TT$num <- 1:nrow(data.TT)
tkbgcol <- c("grey95","grey95","grey95","grey95","grey95","grey95","grey95","grey95","grey95","grey95","grey95","grey95","grey95","grey95","grey95","grey95","grey95","grey95","grey95","grey95","grey95","grey95","grey95","grey95")
tkmargin <- 0.02
tkheight <- 0.1
tklinecoord <- c(0.25,0.75)
tklinecolor <- c("grey","grey")
hmapcols <- c("#0016DB","#FFFFFF","#FFFF00")
heightlines <- 0.04
marginlines <- 0.01
tkborder <- ""
innergap <- 0.5
tkbordercol <- NA
tkbardir <- 1
tkrectcol <- 1
selrectcol <- 1
rectcols <- c("#EDEDFD","#6969F5","#00008B")
tktransparency <- 1
data.TTT <- data.T[[i]]
	data.TTT$id <- paste(data.TTT[,1],data.TTT[,2],data.TTT[,3],sep="")
	data.TTT$num <- 1:nrow(data.TTT)
transparencyHlt <- c(1,1,1)
lkmargin <- 0
tkborder <- NA
columns <- c(4,5,6,7,8,9,10,11,12,13,14)
data.TT$num <- NULL
break1 <- min(as.numeric(as.matrix(data.TT[,-c(1:3)])))
			break2 <- max(as.numeric(as.matrix(data.TT[,-c(1:3)])))
			midpoint <- (break1+break2)/2
			f <- colorRamp2(breaks = c(break1, midpoint, break2), colors = hmapcols)
circos.genomicPosTransformLines(data.TT, posTransform = posTransform.default,
                                horizontalLine = "top", track.height = heightlines, track.margin = c(0,marginlines))
            circos.genomicTrackPlotRegion(data.TT, track.height = tkheight, track.margin = c(lkmargin,tkmargin), stack = TRUE,
                              panel.fun = function(region, value, ...){
                                i = getI(...)
								circos.genomicRect(region, value, col = f(value[[1]]), ybottom = i - innergap, ytop = i + innergap,
                                                   border = f(value[[1]]), posTransform = posTransform.default, ...)
                              }, bg.border = NA)
poslabels <- c("outer","outer","outer")
if(poslabels[i]=="inner"){
			    takindx <- takindx+3
			}else{
			    takindx <- takindx+1
			}
i <- 3
data.TT <- data.T[[i]]
	tktype <- typeTrack[i]
	data.TT[,2] <- as.numeric(data.TT[,2])
	data.TT[,3] <- as.numeric(data.TT[,3])
	data.NN <- data.N[[i]]
	data.TT$num <- 1:nrow(data.TT)
data.TTC <- NULL
coltypeTrack <- 2
tkcolor <- c("yellowgreen")
data.TT$num <- NULL
tkbgcol <- c("grey99","grey99","grey99","grey99","grey99","grey99","grey99","grey99","grey99","grey99","grey99","grey99","grey99","grey99","grey99","grey99","grey99","grey99","grey99","grey99","grey99","grey99","grey99","grey99")
tkmargin <- 0.01
tkheight <- 0.06
tklinecoord <- c(0.25,0.75)
tklinecolor <- c("grey","grey")
hmapcols <- c("blue","white","red")
tkborder <- ""
innergap <- 0.5
tkbordercol <- NA
tkbardir <- 1
tkrectcol <- 1
selrectcol <- 1
rectcols <- c("#EDEDFD","#6969F5","#00008B")
tktransparency <- 1
tkcolor <- c("#9ACD32FF")
data.TTT <- data.T[[i]]
	data.TTT$id <- paste(data.TTT[,1],data.TTT[,2],data.TTT[,3],sep="")
	data.TTT$num <- 1:nrow(data.TTT)
transparencyHlt <- c(1,1,1)
lkmargin <- 0.01
tkborder <- NA
columns <- c(4)
selreaTrack <- c(1,1,1)
fillareaTrack=""
fillareaTrack=""
fillareaTrack=""
area <- FALSE
borderset <- NA
lwdnum <- 1
data.TT[,ncol(data.TT)] <- as.numeric(data.TT[,ncol(data.TT)])			
			    circos.genomicTrackPlotRegion(data.TT, track.height = tkheight, track.margin = c(lkmargin,tkmargin),
                              bg.col = tkbgcol, bg.border = tkborder, panel.fun = function(region,value,...){
							  if(nchar(tklinecolor[1])!=0){               
							    xlim <- get.cell.meta.data("xlim")
                                ylim <- get.cell.meta.data("ylim")
								for(k in 1:length(tklinecoord)){
                                y1 <- as.numeric(quantile(ylim,probs=tklinecoord[k]))
                                circos.lines(x=xlim,y=c(y1,y1), col=tklinecolor[k], lwd=0.1)
					            }
				              }
							  if((coltypeTrack==1 && !("color" %in% colnames(data.T[[i]]))) | coltypeTrack==2){
							    if(length(columns)==1){
							       tkcolor <- tkcolor[1]
							    }else{
								   tkcolor <- c(tkcolor,rep("grey",length(columns)))
								   tkcolor <- tkcolor[1:length(columns)]
								}								
								if(selreaTrack[i]==1 | fillareaTrack[i]!="add"){
								   borderset <- adjustcolor(tkcolor,alpha.f = tktransparency)
								}
								circos.genomicLines(region, value, numeric.column=columns-3, col=borderset, area=area, border=tkcolor, lwd=lwdnum, lty=1, ...)
							  }
                              })
poslabels <- c("outer","outer","outer")
if(poslabels[i]=="inner"){
			    takindx <- takindx+3
			}else{
			    takindx <- takindx+1
			}
transparencyLinks <- 0.5
rou <- get_most_inside_radius()
			rou <- rou[1]
colLinks <- adjustcolor(selcolorLinks, alpha.f = transparencyLinks)
transparencyhltLinks <- 1
colL <- unique(datL[,4])
					colL <- adjustcolor(colL, alpha.f = transparencyhltLinks)
		            colL <- gsub("0x","#", colL)
					hltregion1$color <- adjustcolor(hltregion1$color, alpha.f = transparencyhltLinks)
		            hltregion1$color <- gsub("0x","#", hltregion1$color)
					hltregion2$color <- adjustcolor(hltregion2$color, alpha.f = transparencyhltLinks)
		            hltregion2$color <- gsub("0x","#", hltregion2$color)
                    linkk1 <- data.LL1[!data.LL1$num %in% c(hltregion1$num,hltregion2$num),][,c(1:3)]
                    linkk2 <- data.LL2[!data.LL2$num %in% c(hltregion1$num,hltregion2$num),][,c(1:3)]
					colindx <- (!data.LL1$num %in% c(hltregion1$num,hltregion2$num)) & (!data.LL2$num %in% c(hltregion1$num,hltregion2$num))
circos.genomicLink(linkk1, linkk2, rou = rou, col = colLinks, border = NA)
lapply(colL, function(x){
                        hltregion11 <- hltregion1[hltregion1$color %in% x,]
                        hltregion12 <- data.LL2[data.L2$num %in% hltregion11$num,]
                        hltregion12 <- hltregion12[,c(1:3)]
                        hltregion11 <- hltregion11[,c(1:3)]
                        hltregion22 <- hltregion2[hltregion2$color %in% x,]
                        hltregion21 <- data.LL1[data.L1$num %in% hltregion22$num,]
                        hltregion21 <- hltregion21[,c(1:3)]
                        hltregion22 <- hltregion22[,c(1:3)]
                        if(nrow(hltregion11)!=0){
                            circos.genomicLink(hltregion11, hltregion12, rou = rou, col = x, border = NA)
                        }
                        if(nrow(hltregion22)!=0){
                            circos.genomicLink(hltregion21, hltregion22, rou = rou, col = x, border = NA)
                        }
				    })
n <- 1
xleft <- 1.225
xright <- 1.245
ybottom <- -0.03286
ytop <- 0.02714
len <- 0.06
gap <- 0.3
legendtext <- c("1. Chromosome")
for(i in 1:n){
			   assign(paste("n",i,sep=""),legendtext[i])
			}
rect(xleft, ybottom, xright, ytop, col = "black")
            polygon(x=c(xleft-0.01,(xleft+xright)/2,xright+0.01), y=c(ybottom,ybottom-0.02,ybottom), col="black")
            text(x=xleft-0.08, y=ybottom, labels="inner", cex=0.95)
text(x=xleft-0.08, y=ytop-0.01, labels="outer", cex=0.95)
			    text(x=xright+0.025, y=ytop-0.04, labels=get("n1"), cex=1, adj=c(0,0))
dev.off()
circos.clear()
