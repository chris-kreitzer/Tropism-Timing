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
mut


#' Function: create oncoMatrix-like
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
        ma[which(row.names(ma) == i), row] = 1
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

u = OncoMatrix(M = mut)
View(u)


#' make the MetaPlot
mat_origin = ma
mat_origin[mat_origin == 1] = 'singleMet'
mat_origin[is.na(mat_origin)] = ""


#mat_origin = mat_origin[rownames(numMat), colnames(numMat), drop = FALSE]
nsamps = as.numeric(ncol(mat_origin))
percent_alt = apply(ma, 1, function(x) length(x[x != 0]))
percent_alt = paste0(rev(round(percent_alt * 100/nsamps)), "%")

vc_col = c('red', 'green')
names(vc_col) = names = c('singleMet', 'multiHit')

vc_codes = c("", 'singleMet', 'multiHit')
names(vc_codes) = c(0, 1, 2)

#' Plot layout
gene_mar = 5
mat_lo = matrix(data = c(1, 2), nrow = 2, ncol = 1, byrow = TRUE)
lo = graphics::layout(mat = mat_lo, heights = c(12, 2))
par(mar = c(0.5, gene_mar, 1, 2.5), xpd = TRUE)


exclusivity.sort = function(M) {
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

y = exclusivity.sort(M = ma)

nm = t(apply(y$M, 2, rev))
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
  nm = t(apply(y$M, 2, rev))
  nm[nm != names(vc_code)] = NA
  image(x = 1:nrow(nm), 
        y = 1:ncol(nm), 
        z = nm, 
        axes = FALSE, 
        xaxt="n", 
        yaxt="n",
        xlab="", 
        ylab="", 
        col = col, 
        add = TRUE)
}

# Add blanks
bgCol = 'grey85'
borderCol = 'black'
sepwd_genes = 0.5
sepwd_samples = 0.25
fontSize = 0.8

nm = t(apply(y$M, 2, rev))
nm[nm != 0] = NA
image(x = 1:nrow(nm), 
      y = 1:ncol(nm), 
      z = nm, 
      axes = FALSE, 
      xaxt="n", 
      yaxt="n", 
      xlab="", 
      ylab="", 
      col = bgCol, add = TRUE)

# Add grids
abline(h = (1:ncol(nm)) + 0.5, col = borderCol, lwd = sepwd_genes)
abline(v = (1:nrow(nm)) + 0.5, col = borderCol, lwd = sepwd_samples)

mtext(text = colnames(nm)[y$geneOrder], side = 2, at = y$geneOrder,
      font = 3, line = 0.4, cex = fontSize, las = 2)
mtext(text = percent_alt[y$geneOrder], side = 4, at = y$geneOrder,
      font = 3, line = 0.4, cex = fontSize, las = 2)

lep = legend("bottomright", legend = names(vc_col[vc_codes[2:length(vc_codes)]]),
             col = vc_col[vc_codes[2:length(vc_codes)]], border = NA, bty = "n",
             pch = 15, xpd = TRUE, xjust = 0, yjust = 0, cex = legendFontSize)

#' bottom annotation
mtext(text = colnames(ma), side = 1,
      font = 0.1, line = 0.4, cex = fontSize, las = 2, at = 1:ncol(ma))
