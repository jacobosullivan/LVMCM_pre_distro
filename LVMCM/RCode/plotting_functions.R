require(data.table)
source("LVMCM_src/LVMCM/RCode/new_aes.R")
require(ggplot2)
require(reshape2)
require(ggnetwork)
require(network)
require(spatgraphs)
require(gridExtra)
require(colorRamps)

### Plot metacommunity assembly (regional diversity as a function of invasion event)
plotAssembly <- function(out) {
  plot(out$S.R[,1], type='l', xlab="Invasion", ylab=expression(gamma), ylim=c(0, max(1.1*out$S.R[,1])))
  if (max(out$S.R[,2])>0) {
    lines(out$S.R[,2], col="blue")
  }
  if (max(out$S.R[,2])>0) {
    legend("topleft", legend=c("prod", "cons"), lty=1, col=c("black", "blue"))
  }
}

### Image the competitive overlap and trophic interaction matrices

plotInteractionMat <- function(out) {
  p <- list()
  if (!is.null(out$C)) {
    p[[length(p)+1]] <- ggplot(data = melt(t(apply(out$C, 2, rev))), aes(x=Var1, y=Var2, fill=value)) + 
      geom_tile() +
      scale_fill_gradient2(high = "steelblue3", low = "white", 
                           midpoint = 0, limit = c(0,1), space = "Lab", 
                           name=expression(bold(C)[ij])) +
      theme_bw() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            axis.title=element_blank(),
            axis.text=element_blank(),
            axis.ticks=element_blank()
      ) +
      coord_fixed(ratio = 1)
  }
  
  if (!is.null(out$A)) {
    p[[length(p)+1]] <- ggplot(data = melt((apply(log10(out$A), 2, rev))), aes(x=Var1, y=Var2, fill=value)) + 
      # geom_tile(colour="grey70") +
      geom_tile() +
      scale_fill_gradient(low = "white", high = "firebrick3", na.value = "white",
                           limit = c(-5,max(log10(out$A))), space = "Lab", 
                           name=expression(paste(log[10],bold(A)[ik]))) +
      theme_bw() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            axis.title=element_blank(),
            axis.text=element_blank(),
            axis.ticks=element_blank()
      ) +
      coord_fixed(ratio = nrow(out$B_c)/nrow(out$B_c))
    
  } 
  
  if(length(p) > 1) {
    print(grid.arrange(p[[1]], p[[2]], ncol=2))
  } else {
    print(p[[1]])  
  }
  return(p)
}

### Plot spatially explicit growth rate distribution (producers)

plotGrowthDist <- function(out) {
  
  p <- list()
  
  p[[1]] <- ggplot(mapping = aes(x=Var1, y=Var2)) +
    geom_tile(aes(fill=value), data = melt(t(apply(out$R, 2, rev)))) +
    # scale_fill_gradient(low="steelblue3", high='firebrick3', name=expression(bold(R)[ix]), limit=c(min(out$R),max(out$R)), na.value = "white") +
    scale_fill_gradient(low="white", high='firebrick3', name=expression(bold(R)[ix]), limit=c(min(out$R),max(out$R)), na.value = "white") +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank()) +
    coord_fixed(ratio = ncol(out$R)/nrow(out$R)) +
    xlab("local community") +
    ylab("species")
  
  out$D[out$D != 0] <- 1
  diag(out$D) <- 0
  gra <- network.initialize(nrow(out$X), directed=F)
  network.adjacency(out$D, gra)
  n <- ggnetwork(gra, layout=out$X, arrow.gap=0)
  
  randspp <- sample(1:nrow(out$R), 1, replace=F)
  r_i <- data.frame(R=out$R[randspp,], X=n$x[1:nrow(out$X)], Y=n$y[1:nrow(out$X)], xend=rep(0, nrow(out$X)), yend=rep(0, nrow(out$X)))
  
  p[[2]] <- ggplot(n, aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_edges(color = "grey80") +
    geom_point(size=5, col="grey80") +
    geom_point(size=3.5, col="white") +
    geom_point(dat=r_i, aes(x=X, y=Y, color = R), size=3.5) +
    # scale_colour_gradient(low="steelblue3", high='firebrick3', name=expression(bold(R)[ix]), limit=c(min(out$R[randspp,]),max(out$R[randspp,])), na.value = "white") +
    scale_colour_gradient(low="white", high='firebrick3', name=expression(bold(R)[ix]), limit=c(min(out$R[randspp,]),max(out$R[randspp,])), na.value = "white") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    theme(panel.background = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())
  
  print(grid.arrange(p[[1]],p[[2]],ncol=2))
  
  return(p)
}

### Plot biomass distribution

plotBiomassDist <- function(out) {
  
  p <- list()
  
  ## Producer species
  thresh <- 1e-4
  out$B <- log10(out$B)
  out$B[out$B < -4] <- NaN
  src <- snk <- out$B
  src[out$B_src != 1] <- NaN
  snk[out$B_src != -1] <- NaN
  
  p[[1]] <- ggplot(mapping = aes(x=Var1, y=Var2)) + 
    geom_tile(aes(fill=value), data = melt(t(apply(src, 2, rev)))) +
    scale_fill_gradient(low = "white", high = "steelblue3",
                        limit = c(-4,1), space = "Lab",
                        name=expression(atop(paste(log[10], B[list(p,ix)]), "Source")),
                        na.value = "transparent") +
    new_scale_fill() +
    geom_tile(data = melt(t(apply(snk, 2, rev))), aes(fill=value)) +
    scale_fill_gradient(low = "white", high = "firebrick3",
                        limit = c(-4,1), space = "Lab",
                        name="Sink",
                        na.value = "transparent") +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank()) +
    coord_fixed(ratio = ncol(out$B)/nrow(out$B)) +
    xlab("local community") +
    ylab("species")
  
  out$D[out$D != 0] <- 1
  diag(out$D) <- 0
  gra <- network.initialize(nrow(out$X), directed=F)
  network.adjacency(out$D, gra)
  n <- ggnetwork(gra, layout=out$X, arrow.gap=0)
  
  randspp <- sample(1:nrow(out$B), 1, replace=F)
  src_i <- data.frame(B=src[randspp,], X=n$x[1:nrow(out$X)], Y=n$y[1:nrow(out$X)], xend=rep(0, nrow(out$X)), yend=rep(0, nrow(out$X)))
  snk_i <- data.frame(B=snk[randspp,], X=n$x[1:nrow(out$X)], Y=n$y[1:nrow(out$X)], xend=rep(0, nrow(out$X)), yend=rep(0, nrow(out$X)))
  
  p[[2]] <- ggplot(n, aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_edges(color = "grey80") +
    geom_point(size=5, col="grey80") +
    geom_point(size=3.5, col="white") +
    geom_point(dat=src_i, aes(x=X, y=Y, color = B), size=3.5) +
    scale_color_gradient(low = "white", 
                         high = "steelblue3", 
                         limits=c(log10(thresh),1),
                         na.value = "transparent",
                         name="Source") +
    new_scale_colour() + # see https://eliocamp.github.io/codigo-r/2018/09/multiple-color-and-fill-scales-with-ggplot2/
    geom_point(dat=snk_i, aes(x=X, y=Y, color = B), size=3.5) +
    scale_color_gradient(low = "white", 
                         high = "firebrick3", 
                         limits=c(log10(thresh),1),
                         na.value = "transparent",
                         name="Sink") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    theme(panel.background = element_blank(), 
          axis.title = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank())
  
  ### Consumer species
  if (!is.null(out$B_c)) {
    out$B_c <- log10(out$B_c)
    out$B_c[out$B_c < -4] <- NaN
    src <- snk <- out$B_c
    src[out$B_c_src != 1] <- NaN
    snk[out$B_c_src != -1] <- NaN
    
    p[[3]] <- ggplot(mapping = aes(x=Var1, y=Var2)) + 
      geom_tile(aes(fill=value), data = melt(t(apply(src, 2, rev)))) +
      scale_fill_gradient(low = "white", high = "steelblue3",
                          limit = c(-4,1), space = "Lab",
                          name=expression(atop(paste(log[10], B[list(c,ix)]), "Source")),
                          na.value = "transparent") +
      new_scale_fill() +
      geom_tile(data = melt(t(apply(snk, 2, rev))), aes(fill=value)) +
      scale_fill_gradient(low = "white", high = "firebrick3",
                          limit = c(-4,1), space = "Lab",
                          name="Sink",
                          na.value = "transparent") +
      theme_bw() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            axis.text=element_blank(),
            axis.ticks=element_blank()) +
      coord_fixed(ratio = ncol(out$B_c)/nrow(out$B_c)) +
      xlab("local community") +
      ylab("species")
    
    randspp <- sample(1:nrow(out$B_c), 1, replace=F)
    src_i <- data.frame(B=src[randspp,], X=n$x[1:nrow(out$X)], Y=n$y[1:nrow(out$X)], xend=rep(0, nrow(out$X)), yend=rep(0, nrow(out$X)))
    snk_i <- data.frame(B=snk[randspp,], X=n$x[1:nrow(out$X)], Y=n$y[1:nrow(out$X)], xend=rep(0, nrow(out$X)), yend=rep(0, nrow(out$X)))
    
    p[[4]] <- ggplot(n, aes(x = x, y = y, xend = xend, yend = yend)) +
      geom_edges(color = "grey80") +
      geom_point(size=5, col="grey80") +
      geom_point(size=3.5, col="white") +
      geom_point(dat=src_i, aes(x=X, y=Y, color = B), size=3.5) +
      scale_color_gradient(low = "white", 
                           high = "steelblue3", 
                           limits=c(log10(thresh),1),
                           na.value = "transparent",
                           name="Source") +
      new_scale_colour() + # see https://eliocamp.github.io/codigo-r/2018/09/multiple-color-and-fill-scales-with-ggplot2/
      geom_point(dat=snk_i, aes(x=X, y=Y, color = B), size=3.5) +
      scale_color_gradient(low = "white", 
                           high = "firebrick3", 
                           limits=c(log10(thresh),1),
                           na.value = "transparent",
                           name="Sink") +
      coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
      theme(panel.background = element_blank(), 
            axis.title = element_blank(), 
            axis.text = element_blank(), 
            axis.ticks = element_blank())
  }
  
  if (length(p) == 2) {
    print(grid.arrange(p[[1]],p[[2]],ncol=2))
  } else {
    print(grid.arrange(p[[1]],p[[2]],p[[3]],p[[4]],ncol=2))
  }
  return(p)
}
