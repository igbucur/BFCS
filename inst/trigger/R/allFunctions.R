####################
# Functions
####################

trigger.build <- function(exp = exp, exp.pos = exp.pos, marker = marker, marker.pos = marker.pos){
  if(nrow(exp)!=nrow(exp.pos)||nrow(marker)!=nrow(marker.pos)||ncol(exp)!=ncol(marker)){
    stop("Check matrix dimensions")
  }
  triggerobj <- new("trigger", exp = exp, exp.pos = exp.pos, marker = marker, marker.pos = marker.pos)
}


solvemat<-function(Y,x){
  coef = solve.qr(qr(x), t(Y))
  Yhat = crossprod(coef, t(x))
}

## compute likage for one expression trait with one locus (could be a sex chromosome locus)
link.stat.c <- function(exp, genotype, gender=NULL) {
  n <- as.integer(length(exp))
  if(is.null(gender)){ gender <- rep(1, n) }
  oo <- NULL
  for (g in unique(gender)){
    idx <- which(gender==g)
    geno <- genotype[idx]
    trait <- exp[idx]
    ni <- as.integer(length(trait))
    ng <- as.integer(length(table(geno)))
    geno2 <- rep(0, length(geno))
    for (i in 1:ng){
      geno2[geno==unique(as.vector(geno))[i]] <- i
    }
                                        #   depend on the proper type coming in!
    storage.mode(trait) <- "double"
    storage.mode(geno2) <- "integer"
    lik0 <- double(1)
    lik1 <- double(1)
                                        #.C("link_stat", ni, trait, geno2, ng, lik0, lik1, DUP=FALSE, PACKAGE="trigger")
    out <- .C("link_stat", ni, trait, geno2, ng, lik0, lik1, DUP=TRUE)
    lik0 <- out[[5]]
    lik1 <- out[[6]]
    oo <- cbind(oo, c(lik0, lik1))
  }
  foo <- n*log(rowSums(oo))
  return(foo[1]-foo[2])
}

# As link.stat.c, but with matrix arguments.  Much faster than apply()
link.stat.xx.c <- function(exp, genotype, gender=NULL) {
  n <- as.integer(ncol(exp))
  if(is.null(gender)){ gender <- rep(1, n) }
  nexp <- as.integer(nrow(exp))
  ngeno <- as.integer(nrow(genotype))
  
  oo <- NULL
  for (g in unique(gender)){
    idx <- which(gender==g)
    geno <- genotype[, idx]
    trait <- exp[, idx]
    ni <- as.integer(length(idx))
    ng <- as.integer(length(table(geno)))
    geno2 <- rep(0, nrow(geno)*ncol(geno))
    for (i in 1:ng){
      geno2[geno==unique(as.vector(geno))[i]] <- i
    }
    dim(geno2) <- c(nrow(geno), ncol(geno))
    storage.mode(trait) <- "double"
    storage.mode(geno2) <- "integer"
    lik0 <- double(nexp*ngeno)
    lik1 <- double(nexp*ngeno)
    
#.C("link_stat_xx", ni , nexp, trait, ngeno,geno2, ng, lik0, lik1, DUP=FALSE, PACKAGE="trigger")
    
    out <- .C("link_stat_xx", ni , nexp, trait, ngeno,
       geno2, ng, lik0, lik1, DUP=TRUE)
    lik0 <- out[[7]]
    lik1 <- out[[8]]
    oo <- cbind(oo, c(lik0, lik1))
    
  }
  foo <- n*log(rowSums(oo))
  foo <- matrix(foo,nrow=nexp*ngeno)
  stat <- foo[,1]-foo[,2]
  dim(stat)<-c(nexp, ngeno)
  return(stat)
}

sec.link.stat.c <- function(exp, cisexp, genotype, gender=NULL) {
  n <- as.integer(length(exp))
  if(is.null(gender)){ gender <- rep(1, n) }
  oo <- NULL
  for (g in unique(gender)){
    idx <- which(gender==g)
    geno <- genotype[idx]
    trait <- exp[idx]
    cistrait <- cisexp[idx]
    ni <- as.integer(length(trait))
    ng <- as.integer(length(table(geno)))
    geno2 <- rep(0, length(geno))
    for (i in 1:ng){
      geno2[geno==unique(as.vector(geno))[i]] <- i
    }
    lik0 <- double(1)
    lik1 <- double(1)
    storage.mode(trait) <- storage.mode(cistrait) <- "double"
    storage.mode(geno2) <- "integer"
		#.C("sec_link_stat", ni, trait, cistrait,geno2, ng, lik0, lik1, DUP=FALSE, PACKAGE="trigger")
    out <- .C("sec_link_stat", ni, trait, cistrait,
       geno2, ng, lik0, lik1, DUP=TRUE)
    lik0 <- out[[6]]
    lik1 <- out[[7]]
    oo <- cbind(oo, c(lik0, lik1))
  }
  foo <- n*log(rowSums(oo))
  return(foo[1]-foo[2])
}

# As sec.link.stat but with matrix for exp
sec.link.stat.x.c <- function(exp, cisexp, genotype, gender=NULL) {
  n <- as.integer(ncol(exp))
  if(is.null(gender)){ gender <- rep(1, n) }
  nexp <- as.integer(nrow(exp))
  
  oo <- NULL
  for (g in unique(gender)){
    idx <- which(gender==g)
    geno <- genotype[idx]
    trait <- exp[, idx]
    cistrait <- cisexp[idx]
    ni <- as.integer(length(idx))
    ng <- as.integer(length(table(geno)))
    geno2 <- rep(0, length(geno))
    for (i in 1:ng){
      geno2[geno==unique(as.vector(geno))[i]] <- i
    }
    storage.mode(trait) <- storage.mode(cistrait) <- "double"
    storage.mode(geno2) <- "integer"
    lik0 <- double(nexp)
    lik1 <- double(nexp)
    
#.C("sec_link_stat_x", ni, nexp,trait, cistrait, geno2, ng, lik0, lik1, DUP=FALSE, PACKAGE="trigger")
    out <- .C("sec_link_stat_x", ni, nexp, trait, cistrait, geno2, ng, lik0, lik1, DUP=TRUE)
    lik0 <- out[[7]]
    lik1 <- out[[8]]
    oo <- cbind(oo, c(lik0, lik1))
    
  }
  foo <- n*log(rowSums(oo))
  foo <- matrix(foo,nrow=nexp)
  stat <- foo[,1]-foo[,2]
  return(stat)
}


condi.indep.stat.c <- function(exp, cisexp, genotype, gender=NULL) {
  n <- as.integer(length(exp))
  if(is.null(gender)){ gender <- rep(1, n) }
  oo <- 0
  for (g in unique(gender)){
    idx <- which(gender==g)
    geno <- genotype[idx]
    trait <- exp[idx]
    cistrait <- cisexp[idx]
    ni <- as.integer(length(trait))
    ng <- as.integer(length(table(geno)))
    geno2 <- rep(0, length(geno))
    for (i in 1:ng){
      geno2[geno==unique(as.vector(geno))[i]] <- i
    }
    r <- double(1)
    storage.mode(trait) <- storage.mode(cistrait) <- "double"
    storage.mode(geno2) <- "integer"
#.C("condi_indep_stat", ni, trait, cistrait,geno2, ng, r, DUP=FALSE, PACKAGE="trigger")
    out <- .C("condi_indep_stat", ni, trait, cistrait,
       geno2, ng, r, DUP=TRUE)
    r <- out[[6]]
    oo <- oo + r
  }
  return(oo)
}

condi.indep.stat.x.c <- function(exp, cisexp, genotype, gender=NULL) {
  n <- as.integer(ncol(exp))
  if(is.null(gender)){ gender <- rep(1, n) }
  nexp <- as.integer(nrow(exp))
  
  oo <- rep(0, nexp)
  for (g in unique(gender)){
    idx <- which(gender==g)
    geno <- genotype[idx]
    trait <- exp[, idx]
    cistrait <- cisexp[idx]
    ni <- as.integer(length(idx))
    ng <- as.integer(length(table(geno)))
    geno2 <- rep(0, length(geno))
    for (i in 1:ng){
      geno2[geno==unique(as.vector(geno))[i]] <- i
    }
    storage.mode(trait) <- storage.mode(cistrait) <- "double"
    storage.mode(geno2) <- "integer"
    r <- double(nexp)		
#.C("condi_indep_stat_x", ni, nexp,trait, cistrait, geno2, ng, r, DUP=FALSE, PACKAGE="trigger")
    out <- .C("condi_indep_stat_x", ni, nexp,
       trait, cistrait, geno2, ng, r, DUP=TRUE)
    r <- out[[7]]
    oo <- oo+ r
  }
  return(oo)
}

condi.indep.stat.rx.c <- function(exp, cisexp, genotype, gender=NULL) {
  n <- as.integer(ncol(cisexp))
  if(is.null(gender)){ gender <- rep(1, n) }
  nexp <- as.integer(nrow(cisexp))
  
  oo <- rep(0, nexp)
  for (g in unique(gender)){
    idx <- which(gender==g)
    geno <- genotype[idx]
    trait <- exp[idx]
    cistrait <- cisexp[, idx]
    ni <- as.integer(length(idx))
    ng <- as.integer(length(table(geno)))
    geno2 <- rep(0, length(geno))
    for (i in 1:ng){
      geno2[geno==unique(as.vector(geno))[i]] <- i
    }
    storage.mode(trait) <- storage.mode(cistrait) <- "double"
    storage.mode(geno2) <- "integer"
    r <- double(nexp)		
    out <- .C("condi_indep_stat_rx", ni, nexp,
       trait, cistrait, geno2, ng, r, DUP=TRUE, PACKAGE = "trigger")
    r <- out[[7]]
    oo <- oo+ r
  }
  return(oo)
}



### MULTI-LOCUS LINKAGE FUNCTIONS  
filt <- function(marker, w=NULL) {
  m <- nrow(marker)
  k <- ncol(marker)/2
  u <- rep(TRUE,m)
  u[w] <- FALSE
  ww <- w-1
  while((ww > 0) && (sum(marker[w,] != marker[ww,]) < k)) {
    u[ww] <- FALSE
    ww <- ww-1 
  }
  ww <- w+1
  while((ww <= m) && (sum(marker[w,] != marker[ww,]) < k)) {
    u[ww] <- FALSE
    ww <- ww+1 
  }
  return(u)
}




########## Identify the eQTLs (linkage hotspots). 
### This function returns an eigen-R2 measure for each locus in the genome. The eigen-R2 measure of a locus quantifies the proportion 
### of genome-wide expression variation that can be explained by that locus.

get.v <- function(x, eigenG, SST, df0, weights, adjust,N) {
  model <- model.matrix(~1+x)
  df <- ncol(model)
  
  H <- model %*% solve(t(model) %*% model) %*% t(model) 
  res <- eigenG - t(H %*% t(eigenG))
  SSR <- apply(res, 1, function(x) {sum((x-mean(x))^2)})
  R2s <- 1- SSR/SST
  v <- sum(R2s*weights)
  
  if (adjust == TRUE) {
    v <- 1-(1-v)*(N-df0)/(N-df)        
  }
  return(v)
}

#### Plots
eqtlPlot <- function(object, ylim=NULL){
  marker.pos <- object@marker.pos
  exp.pos <- object@exp.pos
  marker <- object@marker
  eqtl.R2 <- object@eqtl.R2
  if(length(eqtl.R2)==0){
    stop("Run eqtl() function before plotting")
  }
  uchr.name <- unique(as.vector(marker.pos[, 1]))
  nchr <- length(uchr.name)
  idx <- which(uchr.name=="X")
  if (length(idx)>0) {
    chr.name <- c(sort(as.integer(uchr.name[-idx])), uchr.name[idx])
  } else {
    chr.name <- sort(as.integer(uchr.name))
  }
  
  markerp <- rep(0, nrow(marker.pos))
  ranget <- 0
  for (i in 1:nchr){
    idxm <- which(marker.pos[,1]==chr.name[i])
    ss <- min(marker.pos[idxm,2])
    ee <- max(marker.pos[idxm,2])
    markerp[idxm] <- ranget[i] + marker.pos[idxm,2]-ss
    ranget <- c(ranget, ranget[i] + (ee-ss))
  }
  
  totl <- max(ranget)
  if (is.null(ylim)) ylim <- c(0, max(eqtl.R2)*1.1)
  plot(markerp, eqtl.R2, xlim=c(1, totl), ylim=ylim, xlab="marker position", ylab="eqtl-R2 value",xaxt="n", "n", pch=16, cex=.5, main = "Eigen-R2")
  lines(markerp, eqtl.R2)
#abline(h=1/(ncol(marker)-1))
  axis(1,at= ranget, labels=F) 
  axis(1,at= c(ranget[-1] + ranget[-(1+nchr)])/2, tick=F, labels = chr.name, cex.axis=.7)
  
}

linkPlot<-function(x, cutoff = 0){
  marker.pos <- x@marker.pos
  exp.pos <- x@exp.pos
  if(length(x@pvalue)==0){
    stop("Run link() function before plotting")
  }
  p <- x@pvalue
  
  pt <- which(p<=cutoff)
  if (length(pt)<10) stop("Please choose a more liberal cutoff.")
  pt.y <- as.integer(pt/nrow(p))+1
  pt.x <- pt-(pt.y-1)*nrow(p)
  pt.y[pt.x==0] <- pt.y[pt.x==0] -1
  pt.x[pt.x==0] <- nrow(p)
  
  gene.mid <- (exp.pos[,2]+exp.pos[,3])/2
  
  uchr.m <- unique(as.vector(marker.pos[, 1]))
  nchr.m <- length(uchr.m)
  idx <- which(uchr.m=="X")
  if (length(idx)>0) {
    chr.m <- c(sort(as.integer(uchr.m[-idx])), uchr.m[idx])
  } else {
    chr.m <- sort(as.integer(uchr.m))
  }
  
  uchr.e <- unique(as.vector(exp.pos[, 1]))
  nchr.e <- length(uchr.e)
  idx <- which(uchr.e=="X")
  if (length(idx)>0) {
    chr.e <- c(sort(as.integer(uchr.e[-idx])), uchr.e[idx])
  } else {
    chr.e <- sort(as.integer(uchr.e))
  }
  markerp <- rep(0, nrow(marker.pos))
  genep <- rep(0, nrow(exp.pos))
  rangem <- rangee <- 0
  for (i in 1:nchr.m){
    idxm <- which(marker.pos[,1]==chr.m[i])
    ss <- min(marker.pos[idxm,2])
    ee <- max(marker.pos[idxm,2])
    markerp[idxm] <- rangem[i] + marker.pos[idxm,2]-ss+1
    rangem <- c(rangem, rangem[i] + (ee-ss)+1)
  }
  
  for (i in 1:nchr.e){
    idxe <- which(exp.pos[,1]==chr.e[i])
    ss <- min(exp.pos[idxe,2])
    ee <- max(exp.pos[idxe,3])
    genep[idxe] <- rangee[i] + gene.mid[idxe]-ss+1
    rangee <- c(rangee, rangee[i] + (ee-ss)+1)
  }
  
  
  totm <- max(rangem)
  tote <- max(rangee)
  plot(markerp[pt.y], genep[pt.x], xlim=c(1, totm), ylim=c(1, tote), xlab="marker position",
       ylab="expression-trait position",xaxt="n",yaxt="n",pch=16, cex=.5,
       main = paste("Cutoff =", cutoff, sep = ""))
  axis(2,at= rangee, labels=F)
  axis(2,at= c(rangee[-1] + rangee[-(1+nchr.e)])/2, tick=F, labels = chr.e, cex.axis=.7)
  axis(1,at= rangem, labels=F) 
  axis(1,at= c(rangem[-1] + rangem[-(1+nchr.m)])/2, tick=F, labels = chr.m, cex.axis=.7)
  
}

mlinkPlot <- function(object, qcut = qcut, bin.size=NULL) {
  marker.pos <- object@marker.pos
  exp.pos <- object@exp.pos
  mlink.obj <-object@mlink
  idx <- which(mlink.obj$qvalue <= qcut)
  if (length(idx)<10) stop("Please choose a more liberal q-value cutoff.")
  qtl <- (mlink.obj$qtl)[idx,]
  qtl <- t(apply(qtl, 1, sort, decreasing=T))
  
  uchr.m <- unique(as.vector(marker.pos[, 1]))
  nchr.m <- length(uchr.m)
  idx <- which(uchr.m=="X")
  if (length(idx)>0) {
    chr.m <- c(sort(as.integer(uchr.m[-idx])), uchr.m[idx])
  } else {
    chr.m <- sort(as.integer(uchr.m))
  }
  
  markerp <- rep(0, nrow(marker.pos))
  if (is.null(bin.size)) bin.size <-1
  binpos <- as.integer(marker.pos[,2]/bin.size)+1
  rangem <- 0
  for (i in 1:nchr.m){
    idxm <- which(marker.pos[,1]==chr.m[i])
    ss <- min(binpos[idxm])
    ee <- max(binpos[idxm])
    markerp[idxm] <- rangem[i] + binpos[idxm]-ss+1
    rangem <- c(rangem, rangem[i] + (ee-ss)+1)
  }
  
  totm <- max(rangem)
  plot(markerp[qtl[,1]], markerp[qtl[,2]], xlim=c(1, totm), ylim=c(1, totm), xlab="marker position", ylab="marker position",xaxt="n",yaxt="n","n", cex=.5)
  oov <- apply(cbind(markerp[qtl[,1]], markerp[qtl[,2]]), 1, paste, collapse=",")
  oon <- table(oov)
  noon <- names(oon)
  for (i in 1:length(names(oon))) {
    pos <- as.integer(strsplit(noon[i],",")[[1]])
    cnt <- as.character(oon[i])
    text(pos[1]-0.5, pos[2]-0.5, cnt, cex=0.5)
  }
  abline(0, 1, lty=2)
  axis(2,at= rangem, labels=F)
  axis(2,at= c(rangem[-1] + rangem[-(1+nchr.m)])/2, tick=F, labels = chr.m, cex.axis=.7)
  axis(1,at= rangem, labels=F) 
  axis(1,at= c(rangem[-1] + rangem[-(1+nchr.m)])/2, tick=F, labels = chr.m, cex.axis=.7)
  
}

### GET PERMUTATION indEX FUNCTION ###

get.loc.markers <- function(exp.index, exp.pos, marker, marker.pos, window.size = NULL) {

    marker.pos <- as.matrix(marker.pos)
    exp.pos <- as.matrix(exp.pos)
    exp.chr <- exp.pos[exp.index, 1]
    if (is.null(window.size)) {
      midx <- which(marker.pos[, 1] == exp.chr)
      return(midx)
    } else{
      exp.coord <- range(as.numeric(exp.pos[exp.index, 2:3]))
      exp.length <- exp.coord[2] - exp.coord[1]
      exp.range <- exp.coord + c(-1, 1) * (window.size - exp.length)/2
      midx <- which((marker.pos[, 1] == exp.chr) & (as.numeric(marker.pos[, 2]) >= exp.range[1]) & (as.numeric(marker.pos[, 2]) <= exp.range[2]))
      if (length(midx) > 1) {
        return(midx)
      } else {
        tidxs <- which(as.vector(marker.pos[, 1]) == exp.chr)
        mleft <- as.numeric(marker.pos[tidxs, 2]) - exp.range[1]
        mleft[mleft>0] <- min(mleft)
        mright <- as.numeric(marker.pos[tidxs, 2])- exp.range[2]
        mright[mright<0] <- max(mright)
        mm1 <- tidxs[order(mleft, decreasing=T)[1:2]]
        mm2 <- tidxs[order(mright, decreasing=F)[1:2]]
        
        mm <- unique(c(mm1, mm2)) 
        if (length(mm)>1) {
          return(mm)
        }else{
          stop("Please select a larger window.size.")
        }
      } }
  }



perm.idx.fun <- function(n, B) {
  Idx.matrix <- matrix(rep(1:n, B), byrow = T, nrow = B)
  Idx.matrix <- apply(Idx.matrix, 1, sample)
  return(Idx.matrix)
}

order.c <- function(x) {
  out <- integer(length(x))
  storage.mode(x) <- "double"
#.C("order_c",x,as.integer(length(x)),out,DUP=FALSE, PACKAGE="trigger")
  hold <- .C("order_c",x,as.integer(length(x)),out,DUP=TRUE)
  out <- hold[[3]]
  return(out)
}

mergeorder.c <- function(x1, x2) {
  storage.mode(x1) <- "double";
  storage.mode(x2) <- "double";
  n1<-as.integer(length(x1))
  n2<-as.integer(length(x2))
  out<-integer(n1+n2)
#.C("mergeorder",n1,x1,n2,x2,out,DUP=FALSE, PACKAGE="trigger")
  hold <- .C("mergeorder",n1,x1,n2,x2,out,DUP=TRUE)
  out <- hold[[5]]
  return(out)
}



get.pi0 <- function (p = NULL, lambda = seq(0, 0.9, 0.05), 
                     pi0.method = "smoother", fdr.level = NULL, robust = FALSE, 
                     gui = FALSE, smooth.df = 3, smooth.log.pi0 = FALSE)
{
  if (min(p) < 0 || max(p) > 1) {
    print("ERROR: p-values not in valid range.", quote=FALSE)
    return(0)
  }
  if (length(lambda) > 1 && length(lambda) < 4) {
    print(paste("ERROR: If length of lambda greater than 1, ",
		"you need at least 4 values."), quote=FALSE)
    return(0)
  }
  if (length(lambda) > 1 && (min(lambda) < 0 || max(lambda) >= 1)) {
    print("ERROR: Lambda must be within [0, 1).", quote=FALSE)
    return(0)
  }
  m <- length(p)
  if (length(lambda) == 1) {
    if (lambda < 0 || lambda >= 1) {
      print("ERROR: Lambda must be within [0, 1).", quote=FALSE)
      return(0)
    }
    pi0 <- mean(p >= lambda)/(1 - lambda)
    pi0 <- min(pi0, 1)
    pi0 <- max(pi0, 0)
  } else {
    pi0 <- rep(0, length(lambda))
    for (i in 1:length(lambda)) {
      pi0[i] <- mean(p >= lambda[i])/(1 - lambda[i])
    }
    if (pi0.method == "smoother") {
      if (smooth.log.pi0) pi0 <- log(pi0)
      spi0 <- smooth.spline(lambda, pi0, df = smooth.df)
      pi0 <- predict(spi0, x = max(lambda))$y
      if (smooth.log.pi0) pi0 <- exp(pi0)
      pi0 <- min(pi0, 1)
      pi0 <- max(pi0, 0)
    } else if (pi0.method == "bootstrap") {
      minpi0 <- min(pi0)
      mse <- rep(0, length(lambda))
      pi0.boot <- rep(0, length(lambda))
      for (i in 1:100) {
        p.boot <- sample(p, size = m, replace = TRUE)
        for (i in 1:length(lambda)) {
          pi0.boot[i] <- mean(p.boot > lambda[i])/(1 - lambda[i])
        }
        mse <- mse + (pi0.boot - minpi0)^2
      }
      pi0 <- min(pi0[mse == min(mse)])
      pi0 <- min(pi0, 1)
      pi0 <- max(pi0, 0)
    } else {
      print("ERROR: 'pi0.method' must be one of 'smoother' or 'bootstrap'.", quote=FALSE)
      return(0)
    }
  }
  return(pi0)
}


find.cis.genes<-function(out.em, cross, chr.s, marker, marker.pos, exp, exp.pos, thr = thr){
  
  out.em.s = summary(out.em)
  n.sig = length(unique(out.em.s$chr))
  pos.s = as.numeric(max(out.em, chr = chr.s)$pos)
  int.bayes = bayesint(out.em, chr.s, 0.95)
  int.lod = lodint(out.em, chr.s, 1.5)
  int.comb = rbind(int.bayes, int.lod)
  pos.left = min(int.comb$pos); pos.right = max(int.comb$pos)
    
  qtl.left = find.marker(cross, chr.s, pos.left)
  qtl.right = find.marker(cross, chr.s, pos.right)
  qtl.max = find.marker(cross, chr.s, pos.s)
  
  markernames = rownames(marker)
  pos.left = marker.pos[which(markernames == qtl.left), 2]
  pos.right = marker.pos[which(markernames == qtl.right), 2]
  
  marker.max = which(markernames == qtl.max)
  prox.genes = which(exp.pos[,1]==chr.s&exp.pos[,2]>=pos.left&exp.pos[,3]<=pos.right)
  
  result = list(marker.max = marker.max, prox.genes = prox.genes)
  
}

sva.fit<-function(exp, qtl, traitid = NULL, nsv = 0){
  m = nrow(exp)
  n = ncol(exp)
  
  exp.comb = exp
  mod = model.matrix(~1+as.factor(qtl))
  
  if(nsv == 0){
    H =	sva(exp.comb, mod)
  }
  if(nsv != 0){
    H = twostepsva.build(exp.comb, mod, n.sv = nsv)
  }
  
  
  n.sv <- H$n.sv
  resids <- exp.comb
  
  
  if (n.sv>0){
    use.var <- matrix(rep(FALSE, n.sv*m), ncol = n.sv)
    pp <- matrix(rep(FALSE, n.sv*m), ncol = n.sv)
    for(i in 1:n.sv) {
      pp[,i] <- sv.support(exp,H$sv[,i])
      use.var[,i] <- rank(pp[,i]) <= (m-m*get.pi0(pp[,i]))
    }
    for (i in 1:nrow(exp.comb)){
      if (sum(use.var[i,])>0) {
        resids[i, ] <- lsfit(H$sv[, use.var[i,]], exp.comb[i,])$residuals 
      }
    }
    resids <- t(apply(resids, 1, function(x) qnorm(rank(x)/(n+1))))
  } else {
    resids <- t(apply(exp, 1, function(x) qnorm(rank(x)/(n+1))))
  }
  
  result = list(t1.fit = resids[1:m,], t2.fit = resids[traitid,], n.sv = n.sv)
  return(result)
}

traitmap.fun <- function(exp, trait, cis.list, markermorph, B= 5, norm=TRUE, df0){
  
  m.id = as.numeric(unique(cis.list[,1]))
  n = ncol(exp)
  if (norm) exp <- t(apply(exp, 1, function(x) qnorm(rank(x)/(n+1))))
  if (norm) trait <- qnorm(rank(trait)/(n+1))
  qtl = markermorph[m.id, ]
  qtl.mod = model.matrix(~1 + qtl)
  
  ps.list <- list()
  
#for(i in 1:nrow(cis.list)){
#	g.id = as.numeric(cis.list[i,2])
  condi <- condi.indep.stat.rx.c(exp = trait, cisexp = exp,  genotype = qtl)
  condi <- condi[!is.na(condi)]
  
  stat0 <- NULL
  for (k in 1:B){
    
    trait.null = c()
    for (l in unique(qtl)){ trait.null[qtl==l] <- sample(trait[qtl==l])}
    resids0 <- qnorm(rank(trait.null)/(n+1))
    oo <- condi.indep.stat.rx.c(exp = trait.null, cisexp = exp,  genotype = qtl)
    stat0 = c(stat0, oo)
  }
  ps <- unlist(lapply(condi, function(x,stat0){ mean(stat0<=x)} , stat0=stat0),use.names=FALSE)[as.numeric(cis.list[,2])]
  names(ps) = cis.list[, 3]
#ps.list[[cis.list[i ,3]]] <- ps
  
  
  
  return(ps)
}


########################################################################################
#Title: Extraction and Analysis of Differential Gene Expression 
#Version: 2.0
#Author: John D. Storey <jstorey@princeton.edu>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
########################################################################################

edge.lfdr <- function(p, trunc=TRUE, monotone=TRUE, transf=c("probit", "logit"), adj=1.5, eps=10^-8, pi0 = NULL, ...) {
	
	if(is.null(pi0)) pi0 = edge.qvalue(p, ...)$pi0
	n = length(p)
	transf = match.arg(transf)
	
	if(transf=="probit") {
          p = pmax(p, eps)
          p = pmin(p, 1-eps)
          x = qnorm(p)
          myd = density(x, adjust=adj)
          mys = smooth.spline(x=myd$x, y=myd$y)
          y = predict(mys, x)$y
          lfdr = pi0*dnorm(x)/y
	}
	
	if(transf=="logit") {
          x = log((p+eps)/(1-p+eps))
          myd = density(x, adjust=adj)
          mys = smooth.spline(x=myd$x, y=myd$y)
          y = predict(mys, x)$y
          dx = exp(x)/(1+exp(x))^2
          lfdr = pi0*dx/y
	}
	
	if(trunc) {lfdr[lfdr > 1] = 1}
	if(monotone) {	
          lfdr = lfdr[order(p)]
          for(i in 2:n) {
            if(lfdr[i] < lfdr[i-1]) {lfdr[i] = lfdr[i-1]}
          }
          lfdr = lfdr[rank(p)]
	}
	
	return(lfdr)
      }

edge.qvalue <- function(p,lambda = seq(0, 0.9, 0.05), pi0.method = "smoother", fdr.level = NULL, robust = FALSE,smooth.df = 3, smooth.log.pi0 = FALSE, ...) {
  
  err.func <- "edge.qvalue"
  if (min(p) < 0 || max(p) > 1) {
    err.msg(err.func,"P-values not in valid range.")
    return(invisible(1))
  }
  if (length(lambda) > 1 && length(lambda) < 4) {
    err.msg(err.func,"If length of lambda greater than 1, you need at least 4 values.")
    return(invisible(1))
  }
  if (length(lambda) > 1 && (min(lambda) < 0 || max(lambda) >= 1)) {
    err.msg(err.func,"Lambda must be in [0,1).")
    return(invisible(1))
  }
  m <- length(p)
  if (length(lambda) == 1) {
    if (lambda < 0 || lambda >= 1) {
      err.msg(err.func,"Lambda must be in [0,1).")
      return(invisible(1))
    }
    pi0 <- mean(p >= lambda)/(1 - lambda)
    pi0 <- min(pi0, 1)
  } else {
    pi0 <- rep(0, length(lambda))
    for (i in 1:length(lambda)) {
      pi0[i] <- mean(p >= lambda[i])/(1 - lambda[i])
    }
    if (pi0.method == "smoother") {
      if (smooth.log.pi0){ 
        pi0 <- log(pi0)
        spi0 <- smooth.spline(lambda, pi0, df = smooth.df)
        pi0 <- predict(spi0, x = max(lambda))$y
      }
      if (smooth.log.pi0) {
        pi0 <- exp(pi0)
      }
      pi0 <- min(pi0, 1)
    }
    else if (pi0.method == "bootstrap") {
      minpi0 <- min(pi0)
      mse <- rep(0, length(lambda))
      pi0.boot <- rep(0, length(lambda))
      for (i in 1:100) {
        p.boot <- sample(p, size = m, replace = TRUE)
        for (i in 1:length(lambda)) {
          pi0.boot[i] <- mean(p.boot > lambda[i])/(1 - lambda[i])
        }
        mse <- mse + (pi0.boot - minpi0)^2
      }
      pi0 <- min(pi0[mse == min(mse)])
      pi0 <- min(pi0, 1)
    }
    else {
      err.msg(err.func,"'pi0.method' must be one of 'smoother' or 'bootstrap'")
      return(invisible(1))
    }
  }
  if (pi0 <= 0) {
    err.msg(err.func,"The estimated pi0 <= 0. Check that you have valid\np-values or use another lambda method.")
    return(invisible(1))
  }
  if (!is.null(fdr.level) && (fdr.level <= 0 || fdr.level > 1)) {
    err.msg(err.func,"'fdr.level' must be within (0,1].")
    return(invisible(1))
  }
  u <- order(p)
  qvalue.rank <- function(x) {
    idx <- sort.list(x)
    fc <- factor(x)
    nl <- length(levels(fc))
    bin <- as.integer(fc)
    tbl <- tabulate(bin)
    cs <- cumsum(tbl)
    tbl <- rep(cs, tbl)
    tbl[idx] <- tbl
    return(tbl)
  }
  v <- qvalue.rank(p)
  qvalue <- pi0 * m * p/v
  if (robust) {
    qvalue <- pi0 * m * p/(v * (1 - (1 - p)^m))
  }
  qvalue[u[m]] <- min(qvalue[u[m]], 1)
  for (i in (m - 1):1) {
    qvalue[u[i]] <- min(qvalue[u[i]], qvalue[u[i + 1]], 1)
  }
  if (!is.null(fdr.level)) {
    retval <- list(call = match.call(), pi0 = pi0, qvalues = qvalue, pvalues = p, fdr.level = fdr.level, significant = (qvalue <= fdr.level), lambda = lambda)
  }
  else {
    retval <- list(call = match.call(), pi0 = pi0, qvalues = qvalue, pvalues = p, lambda = lambda)
  }
  class(retval) <- "qvalue"
  return(retval)
}


err.msg <- function(err.func = "edge",msg) {
  cat('\n')
  cat('\t')
  cat('ERROR in the',err.func,'function: ','\n')
  cat('\t',msg,'\n\n')
}

edge.pvalue <- function(stat, stat0, pool=TRUE) {
  err.func <- "edge.pvalue"
  m <- length(stat)
  if(pool==TRUE) {
    if(is.matrix(stat0)) {stat0 <- as.vector(stat0)}
    m0 <- length(stat0) 
    v <- c(rep(T, m), rep(F, m0))
    v <- v[order(c(stat,stat0), decreasing = TRUE)]
    u <- 1:length(v)
    w <- 1:m
    p <- (u[v==TRUE]-w)/m0
    p <- p[rank(-stat)]
    p <- pmax(p,1/m0)
  } else {
    if(is.vector(stat0)) {
      err.msg(err.func,"stat0 must be a matrix.")
      return(invisible(1))
    }
    if(ncol(stat0)==m) {stat0 <- t(stat0)}
    if(nrow(stat0)!=m){
      err.msg(err.func,"Number of rows of stat0 must equal length of stat.")
      return(invisible(1))
    }
    stat0 <- (stat0 - matrix(rep(stat,ncol(stat0)),byrow=FALSE,nrow=m)) >= 0
    p <- apply(stat0,1,mean)
    p <- pmax(p,1/ncol(stat0))
  }
  return(p)
}

sv.support <- function(dat,y){
### Function by Jeff Leek and John Storey
  n <- dim(dat)[2]
  m <- dim(dat)[1]
  p <- rep(0,m)
  one <- rep(1,n)
  Id <- diag(one)
  X <- cbind(one,y)
  XXi <- solve(t(X)%*%X)
  resid <- t((Id- X%*%XXi%*%t(X))%*%t(dat))
  est <- (XXi%*%t(X)%*%t(dat))[2,]
  se <- sqrt(XXi[2,2]*resid^2%*%rep(1,n)/(n-2))
  p <-  2*(1-pt(abs(est/se),(n-2)))
  return(p)
}



