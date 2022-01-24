library(data.table)


Interactions = function(M, pvalue = c(0.05, 0.01), 
                        colPal = "BrBG", 
                        limitColorBreaks = T, 
                        fontSize = 0.8,
                        sigSymbolsSize = 2,
                        sigSymbolsFontSize = 0.9, 
                        pvSymbols = c(46,42),
                        colNC = 9, 
                        nShiftSymbols = 5){
  #' create binary Matrix and transpose;
  binaryMatrix = OncoMatrix(M = M)
  binaryMatrix = t(binaryMatrix)
  binaryMatrix = binaryMatrix[rowSums(binaryMatrix) > 0, ]
  
  #' convert any integer > 1 to 1
  binaryMatrix[binaryMatrix > 0] = 1
  
  #' pairwise fisher test
  interactions = sapply(1:ncol(binaryMatrix), 
                        function(i) sapply(1:ncol(binaryMatrix), 
                                           function(j){
                                             f = try(fisher.test(binaryMatrix[,i], 
                                                                 binaryMatrix[,j]), 
                                                     silent = TRUE); 
                                             if(class(f) == "try-error") NA 
                                             else ifelse(f$estimate > 1, -log10(f$p.val), log10(f$p.val))}))
  #' ODDS
  oddsRatio = oddsGenes = sapply(1:ncol(binaryMatrix), 
                                 function(i) sapply(1:ncol(binaryMatrix), 
                                                    function(j){
                                                      f = try(fisher.test(binaryMatrix[, i], 
                                                                          binaryMatrix[, j]), 
                                                              silent = TRUE); 
                                                      if(class(f) == "try-error") f = NA
                                                      else f$estimate}))
  #' naming
  rownames(interactions) = colnames(interactions) = rownames(oddsRatio) = colnames(oddsRatio) = colnames(binaryMatrix)
  
  sigPairs = which(x = 10^-abs(interactions) < 1, arr.ind = TRUE)
  sigPairs2 = which(x = 10^-abs(interactions) >= 1, arr.ind = TRUE)
  
  #' modify output
  sigPairs = rbind(sigPairs, sigPairs2)
  sigPairsTbl = data.table::rbindlist(
    lapply(X = seq_along(1:nrow(sigPairs)), function(i) {
      x = sigPairs[i, ]
      g1 = rownames(interactions[x[1], x[2], drop = FALSE])
      g2 = colnames(interactions[x[1], x[2], drop = FALSE])
      tbl = as.data.frame(table(apply(X = binaryMatrix[, c(g1, g2), drop = FALSE], 1, paste, collapse = "")))
      combn = data.frame(t(tbl$Freq))
      colnames(combn) = tbl$Var1
      pval = 10^-abs(interactions[x[1], x[2]])
      fest = oddsRatio[x[1], x[2]]
      d = data.table::data.table(Organ1 = g1,
                                 Organ2 = g2,
                                 pValue = pval, 
                                 oddsRatio = fest)
      d = cbind(d, combn)
      d
    }), 
    fill = TRUE)
  
  sigPairsTbl = sigPairsTbl[!Organ1 == Organ2] #Remove diagonal elements
  sigPairsTbl[is.na(sigPairsTbl)] = 0
  sigPairsTbl$Event = ifelse(test = sigPairsTbl$oddsRatio > 1, 
                             yes = "Co_Occurence", 
                             no = "Mutually_Exclusive")
  sigPairsTbl$pair = apply(X = sigPairsTbl[,.(Organ1, Organ2)], 
                           MARGIN = 1, 
                           FUN = function(x) paste(sort(unique(x)), collapse = ", "))
  sigPairsTbl[,event_ratio := `01`+`10`]
  sigPairsTbl[,event_ratio := paste0(`11`, '/', event_ratio)]
  sigPairsTblSig = sigPairsTbl[order(as.numeric(pValue))][!duplicated(pair)]
  
  
  #' Visualization:
  if(nrow(interactions) >= 5){
    diag(interactions) = 0
    m = nrow(interactions)
    n = ncol(interactions)
    
    
    col_pal = RColorBrewer::brewer.pal(9, colPal)
    col_pal = grDevices::colorRampPalette(colors = col_pal)
    col_pal = col_pal(m*n-1)
    
    
    interactions[lower.tri(x = interactions, diag = TRUE)] = NA
    
    par(bty = "n", 
        mar = c(1, 4, 4, 2) + .1, 
        las = 2, 
        fig = c(0, 1, 0, 1))
    
    # adjust breaks for colors according to predefined legend values
    breaks = NA
    if(limitColorBreaks){
      minLog10pval = 3
      breaks = seq(-minLog10pval, minLog10pval, length.out = m * n+1)
      interactions4plot  = interactions
      interactions4plot[interactions4plot < (-minLog10pval)] = -minLog10pval
      interactions4plot[interactions4plot > minLog10pval] = minLog10pval
      interactions = interactions4plot
    }
    
    image(x = 1:n, 
          y = 1:m, 
          z = interactions, 
          col = col_pal,
          xaxt = "n", 
          yaxt = "n",
          xlab = "",
          ylab = "", 
          xlim = c(0, n+1), 
          ylim = c(0, n+1),
          breaks = seq(-3, 3, length.out = (nrow(interactions) * ncol(interactions))))
    
    abline(h = 0:n + .5, col = "white", lwd = .5)
    abline(v = 0:n + .5, col = "white", lwd = .5)
    
    mtext(side = 2, at = 1:m, text = rownames(interactions), cex = fontSize, font = 3)
    mtext(side = 3, at = 1:n, text = rownames(interactions), cex = fontSize, font = 3)
  
    
   #' show significance symbols
    w = arrayInd(which(10^-abs(interactions) < min(pvalue)), rep(m, 2))
    points(w, pch = pvSymbols[2], col = "black", cex = sigSymbolsSize)
    w = arrayInd(which((10^-abs(interactions) < max(pvalue)) & (10^-abs(interactions) > min(pvalue))), rep(m,2))
    points(w, pch=pvSymbols[1], col = "black", cex = sigSymbolsSize)
    
    points(x = n - nShiftSymbols, 
             y = 0.7*n, 
             pch = pvSymbols[2], 
             cex = sigSymbolsSize) # "*"
    text(x = n - nShiftSymbols, 
           y = 0.7*n, 
           paste0(" P < ", min(pvalue)), 
           pos = 4, 
           cex = sigSymbolsFontSize, 
           adj = 0)
    points(x = n - nShiftSymbols, 
             y = 0.65*n,
             pch = pvSymbols[1],
             cex = sigSymbolsSize) # "."
    text(x = n - nShiftSymbols,
           y = 0.65*n, 
           paste0(" P < ", max(pvalue)), pos = 4, 
         cex = sigSymbolsFontSize)
    
    #' add legend
    par(fig = c(0.4, 0.7, 0, 0.4), new = TRUE)
    image(
      x = c(0.8, 1),
      y = seq(0, 1, length.out = 200),
      z = matrix(seq(0,1,length.out = 200), nrow = 1),
      col = col_pal, 
      xlim = c(0, 1), 
      ylim = c(0, 1), 
      axes = FALSE, 
      xlab = NA, 
      ylab = NA
    )
    
    atLims = seq(0, 1, length.out = 7)
    axis(side = 4, 
         at = atLims,  
         tcl = -.15, 
         labels = c("> 3 (Mutually exclusive)", 2, 1, 0, 1, 2, ">3 (Co-occurence)"), 
         lwd = .5, 
         cex.axis = sigSymbolsFontSize, 
         line = 0.2)
    text(x = 0.4, 
         y = 0.5, 
         labels = "-log10(P-value)", 
         srt = 90, 
         cex = sigSymbolsFontSize, xpd = TRUE)
  }
  
  return(sigPairsTblSig)
}
  
  
 



