## Create a OncoMatrix for MetaPlot
## Create MetaPlot
## 
## inspired by: https://github.com/BIMIB-DISCo/TRONCO/blob/master/R/visualization.R
## and: http://bioconductor.org/packages/devel/bioc/vignettes/maftools/inst/doc/oncoplots.html
## 
## 01/21/22
## chris-kreitzer


clean()


## Input:
## see script: COAD_KRASmut.R [mut table]
# write.table(mut, file = 'Data/binary_mutation_test.txt', sep = '\t', row.names = F)
mut = read.csv('Data/binary_mutation_test.txt', sep = '\t')


## Functions: 
#' Sort the gene matrix
.Exclusivity.sort = function(M) {
  geneOrder = sort(rowSums(M),
                   decreasing = TRUE,
                   index.return = TRUE)$ix
  scoreCol = function(x) {
    score = 0;
    for (i in 1:length(x)) {
      if(x[i]) {
        score = score + 2^(length(x)-i);
      }
    }
    return(score);
  }
  scores = apply(M[geneOrder, , drop = FALSE ], 2, scoreCol);
  sampleOrder = sort(scores, decreasing=TRUE, index.return=TRUE)$ix
  res = list()
  res$geneOrder = geneOrder
  res$sampleOrder = sampleOrder
  res$M = M[geneOrder, sampleOrder]
  
  return(res);
}


#' create OncoMatrix
OncoMatrix = function(M){
  unique_sites = as.character(unique(unlist(M)))
  unique_sites = unique_sites[unique_sites != 'na']
  #' matrix
  ma = matrix(nrow = length(unique_sites), 
              ncol = nrow(M), 
              dimnames = list(unique_sites, 
                              paste0('Patient', seq(1, nrow(M)))))
  
  #' fill matrix with values
  for(row in 1:nrow(M)){
    string = as.character(M[row, ])
    for(i in string){
      if(i %in% row.names(ma)){
        if(length(string[which(string == i)]) > 1){
          ma[which(row.names(ma) == i), row] = 2
        } else {
          ma[which(row.names(ma) == i), row] = 1
        }
      }
    }
  }
  
  #' post modifications
  x = which(apply(ma, 2, function(x) sum(x, na.rm = T)) == 0)
  ma = ma[,-x]
  ma[is.na(ma)] = 0
  
  rm(i, x, string)
  return(ma)
  
}


#' Visualization
MetaPlot = function(M, 
                    fontSize = 0.8,
                    legendFontSize = 1){
  mat_origin = M
  mat_origin[mat_origin == 0] = ""
  mat_origin[mat_origin == 1] = 'single'
  mat_origin[mat_origin == 2] = 'multi'
  
  nsamps = as.numeric(ncol(mat_origin))
  percent_alt = apply(M, 1, function(x) length(x[x != 0]))
  #percent_alt = paste0(rev(round(percent_alt * 100/nsamps)), "%")
  
  #' color code
  vc_col = c('red', 'black')
  names(vc_col) = c('single', 'multi')
  
  vc_codes = c("", 'single', 'multi')
  names(vc_codes) = c(0, 1, 2)
  
  #' Plot Layout
  gene_mar = 5
  mat_lo = matrix(data = c(1, 2), 
                  nrow = 2, 
                  ncol = 1, 
                  byrow = TRUE)
  lo = graphics::layout(mat = mat_lo, heights = c(12, 2))
  par(mar = c(0.5, gene_mar, 1, 2.5), xpd = TRUE)
  
  #' sorted matrix
  ma.sorted = .Exclusivity.sort(M = M)
  
  #' layout
  nm = t(apply(ma.sorted$M, 2, rev))
  nm[nm == 0] = NA
  image(x = 1:nrow(nm), 
        y = 1:ncol(nm), 
        z = nm, 
        axes = FALSE, 
        xaxt = "n", 
        yaxt = "n",
        xlab = "", 
        ylab = "", 
        col = "white")
  
  vc_codes_temp = vc_codes
  for(i in 2:length(names(vc_codes_temp))){
    vc_code = vc_codes_temp[i]
    col = vc_col[vc_code]
    nm = t(apply(ma.sorted$M, 2, rev))
    nm[nm != names(vc_code)] = NA
    image(x = 1:nrow(nm), 
          y = 1:ncol(nm), 
          z = nm, 
          axes = FALSE, 
          xaxt = "n", 
          yaxt = "n",
          xlab = "", 
          ylab = "", 
          col = col, 
          add = TRUE)
  }
  
  bgCol = 'grey75'
  borderCol = 'white'
  sepwd_genes = 0.5
  sepwd_samples = 0.20
  
  #' adding background
  nm = t(apply(ma.sorted$M, 2, rev))
  nm[nm != 0] = NA
  image(x = 1:nrow(nm), 
        y = 1:ncol(nm), 
        z = nm, 
        axes = FALSE, 
        xaxt = "n", 
        yaxt = "n", 
        xlab = "", 
        ylab = "", 
        col = bgCol, 
        add = TRUE)
  
  #' add grids
  abline(h = (1:ncol(nm)) + 0.5, col = borderCol, lwd = sepwd_genes)
  # abline(v = (1:nrow(nm)) + 0.5, col = borderCol, lwd = sepwd_samples)
  
  
  #' column annotations: MetsRank
  mtext(text = colnames(nm)[ma.sorted$geneOrder], 
        side = 2, 
        at = ma.sorted$geneOrder,
        font = 3, 
        line = 0.4, 
        cex = fontSize, 
        las = 2)
  
  #' add percentage at right side
  mtext(text = rev(paste0(round(percent_alt[ma.sorted$geneOrder] / nsamps * 100), '%')), 
        side = 4, 
        at = 1:length(ma.sorted$geneOrder),
        font = 3, 
        line = 0.4, 
        cex = fontSize, 
        las = 2)
  
  
  lep = legend("bottomright", 
               legend = names(vc_col[vc_codes[2:length(vc_codes)]]),
               col = vc_col[vc_codes[2:length(vc_codes)]], 
               border = NA, 
               bty = "n",
               pch = 15, 
               xpd = TRUE, 
               xjust = 0, yjust = 0, 
               cex = legendFontSize)
  
}

dev.off()
pdf(file = 'MetaPlot_KRASmut.pdf', width = 14, height = 8.5)
MetaPlot(M = x)
title(main = 'MetaPlot for KRAS mutant Colorectal Cancer patients', line = 2, outer = T,
      sub = 'n = 1660 patients')
dev.off()


## Hierchacial clustering
library(ade4)
#' Jaccard Similarity matrix
dist_matrix = as.matrix(dist.binary(df = x, method = 1))
hc = hclust(dist.binary(df = x, method = 1), method = 'ward.D')
plot(hc, hang = 0.02, frame.plot = T, xlab = '', ylab = '', yaxt = 'n', main = '')
rect.hclust(hc , k = 4, border = 2:6)
title(main = 'Ward.D Clustering on KRAS mut patients')


## Cooccur further methods:
#' cooccur: Probabilistic Species Co-Occurrence
#' Analysis in Rfile:///Users/chriskreitzer/Downloads/v69c02%20(2).pdf
#' https://cran.r-project.org/web/packages/cooccur/cooccur.pdf


#' out