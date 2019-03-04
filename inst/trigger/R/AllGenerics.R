
setMethod("initialize", signature ="trigger",
function(.Object, ..., marker, exp, marker.pos, exp.pos){
	
	if(!is.null(rownames(marker))){marker.nam <- rownames(marker)}
	if(is.null(rownames(marker))){marker.nam <- paste("M", c(1:nrow(marker)), sep = "")}
	if(!is.null(rownames(exp))){gene.nam <- rownames(exp)}
	if(is.null(rownames(exp))){gene.nam <- paste("G", c(1:nrow(exp)), sep = "")}

	marker <- as.matrix(marker);
	exp <- as.matrix(exp); rownames(exp) <- gene.nam
	if (sum(c(is.na(marker), is.na(exp)))>0) stop("No missing values found in input data.")
	ng <- length(table(marker)) 
	if (length(ng)>3) stop("More than three genotype classes are found.")
	

	
	marker2 <- matrix(0, ncol=ncol(marker), nrow=nrow(marker))
	for (i in 1:ng){
	    marker2[marker==unique(as.vector(marker))[i]] <- i
	}
	
	chr <- as.vector(marker.pos[, 1])
	u <- (chr=="x" | chr=="X")
	chr[u] <- "X"
	pos <- as.numeric(as.vector(marker.pos[, 2]))
	marker.pos2 <- data.frame(chr, pos)
	rownames(marker2) <-rownames(marker.pos2)<- marker.nam
	
	
	chr <- as.vector(exp.pos[, 1])
	u <- (chr=="x" | chr=="X")
	chr[u] <- "X"
	pos1 <- as.numeric(as.vector(exp.pos[, 2]))
	pos2 <- as.numeric(as.vector(exp.pos[, 3]))
	start <- apply(cbind(pos1, pos2), 1, min)
	end <- apply(cbind(pos1, pos2), 1, max)
	
	exp.pos2 <- data.frame(chr, start, end)
	rownames(exp.pos2) <- gene.nam

	.Object <- callNextMethod(.Object, ...)
	.Object@marker<- marker2
	.Object@exp<- exp
	.Object@marker.pos<- marker.pos2
	.Object@exp.pos<- exp.pos2
	
	.Object
}
)

setMethod("show", "trigger", 
	function(object){
		cat("*** TRIGGER object *** \n" )
		cat("Marker matrix with ", nrow(object@marker), "rows and ", ncol(object@marker), "columns \n")
		cat("Expression matrix with ", nrow(object@exp), "rows and ", ncol(object@exp), "columns \n")
	}
)
setGeneric("plot")
setMethod("plot", signature = c("trigger", "missing"), 
 function(x, y, type = c("link", "mlink", "eqtl"), cutoff = 3.3e-4, qcut = 0.1, bin.size = NULL) {
    if(type ==  "eigenR2"){
		eqtlPlot(x)
	} else if(type ==  "link"){
		linkPlot(x, cutoff)
	} else if (type ==  "mlink"){
		mlinkPlot(x, qcut = qcut, bin.size = bin.size)
	}
	
	}
)

link.trigger <- function(triggerobj, gender = NULL, norm = TRUE){
		marker <- triggerobj@marker
		exp <- triggerobj@exp
		marker.pos <- triggerobj@marker.pos
		
		if (norm == TRUE) {
        	n <- ncol(exp)
        	exp <- t(apply(exp, 1, function(x) qnorm(rank(x)/(n+1)) ) ) 
			storage.mode(exp) <- "double"
    	}
		df <- length(table(marker))-1
		marker.chr <- as.vector(marker.pos[,1])
		if (sum(marker.chr == "X")>0) {
			if (is.null(gender)) stop("Please specifiy the gender of each sample.")
			idx.achr <- which(marker.chr!= "X")
			idx.schr <- which(marker.chr == "X")
			g1 <- gender == unique(gender)[1]
			g2 <- gender == unique(gender)[2]
			ng1 <- length(table(marker[idx.schr,g1]))
			ng2 <- length(table(marker[idx.schr,g2]))
			if (ng1!= ng2){
				stata <- link.stat.xx.c(exp, marker[idx.achr,])
				stats <- link.stat.xx.c(exp, marker[idx.schr,], gender = gender)
				df1 <- max(c(ng1,ng2))-1
				df2 <- ng1+ng2-2
				pvaluea <- 1-pchisq(abs(stata), df = df1)
		    	pvalues <- 1-pchisq(abs(stats), df = df2)
				stat <- matrix(0, nrow = nrow(exp), ncol = nrow(marker))
				pvalue <- stat
				stat[,idx.achr] <- stata
				stat[,idx.schr] <- stats
				pvalue[, idx.achr] <- pvaluea
		    	pvalue[, idx.schr] <- pvalues
			} 
			else{
			stat <- link.stat.xx.c(exp, marker)
			pvalue <- 1-pchisq(abs(stat), df = df)
			}
		} 
		else{
	    	stat <- link.stat.xx.c(exp, marker)
			pvalue <- 1-pchisq(abs(stat), df = df)
			}
			rownames(stat) <- rownames(pvalue) <- rownames(exp)
			colnames(stat) <- colnames(pvalue) <- rownames(marker)
			triggerobj@stat <-stat
			triggerobj@pvalue <-pvalue
			
			triggerobj
}
setGeneric("trigger.link",function(triggerobj, gender = NULL, norm = TRUE) {standardGeneric("trigger.link")})
setMethod("trigger.link", signature = "trigger", link.trigger)

mlink.trigger <- function(triggerobj, prob.cut = 0.9, gender = NULL,  idx = NULL, B = 5, seed = 123){
  set.seed(seed)
  marker <- triggerobj@marker
  exp <- triggerobj@exp
  marker.pos <- triggerobj@marker.pos
  loc.obj <- triggerobj@loc.obj

  if(length(loc.obj)==0){
	cat(paste("Running Local Linkage Scan", "\n"))
	triggerobj = trigger.loclink(triggerobj)
	loc.obj <- triggerobj@loc.obj
	}
	
  ne <- nrow(exp)
  m <- nrow(marker)
  nk <- length(table(marker))
  n <- ncol(exp)

  if(is.null(idx)) {idx <- 1:ne}
  if (length(idx)<100) stop("Please select at least 100 genes to compute multi-locus linkage for them")
  
  obs <- null <- mqtl <- NULL
  count <- 1
  print("Start to calculate multi-locus linkage statistics ..." ,quote = FALSE)

  gen.gender <- NULL
  idx.achr <- which(as.vector(marker.pos[,1])!= "X")
  idx.schr <- which(as.vector(marker.pos[,1]) == "X")
  if (length(idx.schr)>0) {
	if (is.null(gender)) stop("Please specifiy the gender of each sample.")
	g1 <- which(gender == unique(gender)[1])
	g2 <- which(gender == unique(gender)[2])
	ng1 <- length(table(marker[idx.schr, g1]))
	ng2 <- length(table(marker[idx.schr, g2]))
	if (ng1!= ng2) {
	    gen.gender <- gender
		df1 <- max(c(ng1, ng2))-1
		df2 <- ng1+ng2-2
	}
  }
  idxf <- NULL
  for(i in idx) {
	if(loc.obj$prob.loc[i] >= prob.cut){
	#cat(i)
	idxf = c(idxf, i)
	oov <- i
	if (floor(oov/ round(length(idx)/10)) ==  count & count <= 9){

		print(paste(count*10,"% completed", collapse = "", sep = ""),quote = FALSE)
		count <- count +1
		
	}
	if (oov ==  length(idx)){
		print(paste(100,"% completed", collapse = "", sep = ""),quote = FALSE)
	}
	
    
    stat <- rep(0, m)
	exp0 <- t(apply(matrix(rep(exp[i,], each = B), byrow = F, nrow = B), 1, sample))
	stat0 <- matrix(0, nrow = B, ncol = m)
	pstat <- rep(NA, m)
	pstat0 <- matrix(NA, nrow = B, ncol = m)
		
	if (!is.null(gen.gender)) {
		marka <- matrix(marker[idx.achr,],nrow = length(idx.achr))
		stata <- link.stat.xx.c(exp = matrix(exp[i,],nrow = 1), genotype = marka)
		stat[idx.achr] <- stata
		
		pstat[idx.achr] <- pchisq(abs(stata), df = df1, log.p = T)
		marks <- matrix(marker[idx.schr,],nrow = length(idx.schr))
	    stats <- apply(marks, 1, link.stat.c, exp = exp[i,], gender = gen.gender)
		stat[idx.schr] <- stats
		pstat[idx.schr] <- pchisq(abs(stats), df = df2, log.p = T)
		m1 <- which.max(pstat)
		
		stata0 <- link.stat.xx.c(exp = exp0, genotype = marka)
		stat0[, idx.achr] <- stata0
		pstat0[, idx.achr] <- pchisq(abs(stata0), df = df1, log.p = T)
		stats0 <-  link.stat.xx.c(exp = exp0, genotype = marks, gender = gen.gender)
		stat0[, idx.schr] <- stats0
		pstat0[, idx.schr] <- pchisq(abs(stats0), df = df2, log.p = T)
		stat0.i <- apply(cbind(stat0, pstat0), 1, function(x) {
						s <- x[1:m]; p <- x[m+1:m]
						return(s[which.max(p)])
						})
			
	}else {
	    #stat <- apply(marker, 1, link.stat.c, exp = exp[i,])
		stat <- link.stat.xx.c(exp = matrix(exp[i,], nrow = 1), genotype = marker)
		m1 <- which.max(stat)
		
		stat0 <- link.stat.xx.c(exp = exp0, genotype = marker)
		stat0.i <- apply(stat0, 1, function(x) x[which.max(x)])
	}
    

    stat.i <- stat[m1]
    u <- filt(marker,w = m1)
	
	gen.gender2 <- NULL
	idxu <- which(as.vector(marker.pos[u,1]) == "X")
    if (length(idxu)>0) {
	    ng1 <- length(table((marker[u,])[idxu, g1]))
	    ng2 <- length(table((marker[u,])[idxu, g2]))
	if (ng1!= ng2) {
	    gen.gender2 <- gender
		df1 <- max(c(ng1, ng2))-1
		df2 <- ng1+ng2-2
	}
	}

	m2 <- NULL
	for (j in unique(marker[m1,])){
        v <- marker[m1,] == j
		exp0v <- t(apply(matrix(rep(exp[i, v], each = B), byrow = F, nrow = B), 1, sample))
		mark <- matrix(marker[u,v],nrow = sum(u))
        stat <- stat*0
		pstat <- rep(NA, m)
		stat0 <- matrix(0, nrow = B, ncol = m)
		pstat0 <- matrix(NA, nrow = B, ncol = m)
		
	    if (!is.null(gen.gender2)) {
		    marka <- mark[-idxu, ] 
			stata <- link.stat.xx.c(exp = matrix(exp[i,v],nrow = 1), genotype = marka)
			idx1 <- which(u == TRUE)[-idxu]
			idx2 <- which(u == TRUE)[idxu]
	        stat[idx1] <- stata
			
			marks <- mark[idxu, ]
	        #stats <- apply(marks, 1, link.stat.c, exp = exp[i,v], gender = gen.gender2[v])
		    stats <- link.stat.xx.c(exp = matrix(exp[i,v], nrow = 1), gender = gen.gender2[v], genotype = marks )
			stat[idx2] <- stats
			
			pstat[idx1] <- pchisq(abs(stata), df = df1, log.p = T)
			pstat[idx2] <- pchisq(abs(stats), df = df2, log.p = T)
			m2j <- which.max(pstat)
			
			stata0 <- link.stat.xx.c(exp = exp0v, genotype = marka)
			stat0[, idx1] <- stata0
			stats0 <- link.stat.xx.c(exp = exp0v, genotype = marks, gender = gen.gender2[v])  
			stat0[, idx2] <- stats0
			pstat0[, idx1] <- pchisq(abs(stata0), df = df1, log.p = T)
			pstat0[, idx2] <- pchisq(abs(stats0), df = df2, log.p = T)

			stat2.0 <- apply(cbind(stat0, pstat0), 1, function(x) {
						s <- x[1:m]; p <- x[m+1:m]
						return(s[which.max(p)])})
		} else {
		    stat[u] <- link.stat.xx.c(exp = matrix(exp[i,v],nrow = 1), genotype = mark)
			m2j <- which.max(stat)
			
			stat0[, u] <- link.stat.xx.c(exp = exp0v, genotype = mark)
			stat2.0 <- apply(stat0[,u], 1, function(x) x[which.max(x)])
        }
		m2 <- c(m2, m2j)
        stat.i <- c(stat.i, stat[m2j])
		stat0.i <- c(stat0.i, stat2.0)
	}

    stat0.i <- matrix(stat0.i, byrow = F, nrow = B)
	obs <- rbind(obs, stat.i)
	mqtl <- rbind(mqtl, c(m1,m2)) 

	null <- rbind(null, stat0.i)
    }
  }
  
  
  rownames(obs) <- rownames(mqtl) <- paste("gene",idxf)
  colnames(obs) <- colnames(null) <- colnames(mqtl) <- c("first", paste("secondary", 1:nk))
  
  rownames(null) <- rep(paste("gene",idxf), each = B)
  
  probs <- matrix(rep(0, length(idxf)*2),ncol = 2)

  FF <- obs[, 1]
  FF0 <- as.vector(null[, 1])
  FF[FF > max(FF0)] <- max(FF0)
  probs[,1] <- 1 - edge.lfdr(edge.pvalue(FF,FF0), pi0 = 1)

  sec.idx <- apply(cbind(obs[,-1], mqtl[, -1]),1, function(x) (x[nk+(1:nk)])[which.max(x[1:nk])])
  mul.qtl <- cbind(mqtl[,1], sec.idx)
  
  oon <- rownames(exp)[idxf]
  if (is.null(oon)){ gnam <- paste("gene",idxf)
  } else { gnam <- oon}
  rownames(mul.qtl) <- gnam
  colnames(mul.qtl) <- c("Major Locus","Secondary Locus")
  
  FF <- apply(obs[,-1],1,function(x) x[which.max(x)])
  FF0 <- apply(null[,-1],1, function(x) x[which.max(x)])

  FF[FF > max(FF0)] <- max(FF0)
  #p0 <- get.pi0(getp.s(sort(FF),sort(FF0)))
  probs[,2] <- 1 - edge.lfdr(edge.pvalue(FF,FF0), pi0 = 1)

  jprobs <- probs[,1]*probs[,2]

  vv <- cumsum(sort(1-jprobs))
  jqvals <- vv[rank(1-jprobs)]/rank(1-jprobs)

  jprobs <- round(jprobs,4)
  probs1 <- round(probs[,1],4)
  probs2 <- round(probs[,2],4)
  jqvals <- round(jqvals,4)
  prob.list <- cbind(probs1, probs2, jprobs)
  colnames(prob.list) <- c("Major Locus", "Secondary Locus", "Joint Linkage")
  rownames(prob.list) <- names(jqvals) <- gnam
  
  out <- list(qtl = mul.qtl, prob = prob.list, qvalue = jqvals)
  triggerobj@mlink <-out
  return(triggerobj)
}
setGeneric("trigger.mlink",function(triggerobj, prob.cut = 0.9, gender = NULL,  idx = NULL, B = 5, seed = 123){standardGeneric("trigger.mlink")})
setMethod("trigger.mlink", "trigger", mlink.trigger)

eigenR2.trigger <- function(triggerobj, adjust = FALSE, meanR2 = FALSE){  
## Input
## == == == == == == == == == == == == == == == == == 
## inobj: An triggerobj from trigger function.
## adjust: If TRUE, the eigen R-square estimates will be adjusted for sample size effect. For sample size less than 100, we suggest to adjust.
## meanR2: By default, it is FALSE. If TRUE, we standardize each gene expression, and that will yield an average of the R2s for the locus on each gene.
## For more complicated options of eqtl function, please check out our R package eigenR2.

    if (meanR2 == TRUE){
	    exp <- t(apply(triggerobj@exp, 1, function(x) (x-mean(x))/sd(x)))
	} else {
        exp <- t(apply(triggerobj@exp, 1, function(x) x-mean(x)))
    }
	svd.t <- svd(exp) 
    eigenGenes <- t(svd.t$v)  ## each column represents an eigengene
    N <- nrow(eigenGenes)
    eigenG <- eigenGenes[-N, ]
    ds <- svd.t$d[-N]
    weights <- ds^2 / sum(ds^2)

    df0 <- 1
    res0 <- eigenG
    SST <- apply(res0, 1, function(x) sum((x-mean(x))^2))
    eqtl.R2 <- apply(triggerobj@marker, 1, get.v, eigenG = eigenG, SST = SST,df0 = df0, weights = weights, adjust = adjust,N = N)
	
	names(eqtl.R2) <- rownames(triggerobj@marker)
	triggerobj@eqtl.R2 <- eqtl.R2
	return(triggerobj)

}
setGeneric("trigger.eigenR2",function(triggerobj, adjust = FALSE, meanR2 = FALSE){standardGeneric("trigger.eigenR2")})
setMethod("trigger.eigenR2", signature = "trigger", eigenR2.trigger)

loclink.trigger <- function(triggerobj, gender = NULL, window.size = 30000){
    exp <- triggerobj@exp
	marker <- triggerobj@marker
	exp.pos <- triggerobj@exp.pos
	marker.pos <- triggerobj@marker.pos
	stat <- triggerobj@stat
	if(length(stat) ==0 ){
		print(paste("Computing linkage scan..."), quote = F)
		triggerobj <- trigger.link(triggerobj)
		stat <- triggerobj@stat
	}
	pval <- triggerobj@pvalue
	nr <- nrow(exp); nr.dec = round(nr/10)
	pval.loc = rep(0, nr); loc.idx = rep(0, nr)
	print(paste("Computing local-linkages with a window size of ", window.size/1000, " kb", sep = ""), quote = F)
	count <- 1
	for(i in 1:nr){
		if(floor(i/nr.dec) == count){
			print(paste(count*10, "% completed", sep = ""), quote = F)
			count = count + 1
		}
		loc.markers = get.loc.markers(i, exp.pos, marker, marker.pos, window.size = window.size)
		pval.loc[i] = min(pval[i, loc.markers ])
		loc.idx[i] = loc.markers[which.min(pval[i, loc.markers ])]
	}
	#hist(pval.loc)
	prob.loc <- 1- edge.lfdr(pval.loc, lambda = 0)
	#print(paste("Estimated proportion of non-null local linkages is",1 - edge.qvalue(pval.loc, lambda = 0)$pi0, sep = " " ), quote = F)
    loc.obj <- list(prob.loc = prob.loc, loc.idx = loc.idx)
	triggerobj@loc.obj <- loc.obj
    return(triggerobj)
	
}
setGeneric("trigger.loclink",function(triggerobj, gender = NULL, window.size = 30000){standardGeneric("trigger.loclink")})
setMethod("trigger.loclink", signature = "trigger", loclink.trigger)

net.trigger <- function(triggerobj, gender = NULL, idx = NULL, Bsec = 100, prob.cut = 0.7,  include.loc = TRUE, seed = 123, inputfile = NULL){
       ofile = "net_trigg_prob.txt" 
	if (length(inputfile > 0)){
		cat("Scanning pre-computed Trigger probabilities","\n")
		prob <- matrix(scan(inputfile),byrow = TRUE, ncol  = nrow(exp))
		return(prob)}
	
    exp <- triggerobj@exp
	marker <- triggerobj@marker
	exp.pos <- triggerobj@exp.pos
	marker.pos <- triggerobj@marker.pos
	loc.obj<- triggerobj@loc.obj
	if(length(loc.obj)==0){
		cat(paste("Running Local Linkage Scan", "\n"))
		triggerobj = trigger.loclink(triggerobj)
		loc.obj <- triggerobj@loc.obj
	}
	exp <- t(apply(exp, 1, function(x) qnorm(rank(x)/(length(x)+1)) ) ) 
	storage.mode(exp) <- "double"
	
	gnams <- names(loc.obj$prob.loc)[idx]
	
	if (include.loc == TRUE) {
        prob.loc <- loc.obj$prob.loc
	} else{
	    prob.loc <- rep(1, length(loc.obj$prob.loc))
	}
	loc.idx <- loc.obj$loc.idx
	
	set.seed(seed)
	nc <- ncol(exp)
    nr <- nrow(exp)    
	
	exp0.list <- vector("list", Bsec)
	for (j in 1:Bsec){
	   perm.idx.sec <- perm.idx.fun(nc, 1)	
       exp0.list[[j]] <- exp[, perm.idx.sec]
    }
		

    if (is.null(idx)) idx <- 1:nr
	ni <- length(idx); nr.dec = round(nr/10)

	print(paste("Computing network-Trigger regulatory probabilities ..."),quote = FALSE)
    count <- 1
	prob <- NULL
	for (i in idx){
		
	  jp <- rep(0, nr)
		
	  if (prob.loc[i] > prob.cut) {

		sec.stat0 <- ind.stat0 <- NULL
		
		locexp <- exp[i,]
		genotype <- marker[loc.idx[i],]
	    gen.chr <- as.vector(exp.pos[,1])[i]
		gen.gender <- NULL
		if (gen.chr == "X") {
		    if (is.null(gender)) {
				stop("Please specifiy the gender of each sample.") 
			} else{
				ng1 <- length(table(genotype[gender == unique(gender)[1]]))
				ng2 <- length(table(genotype[gender == unique(gender)[2]]))
				if (ng1!= ng2) gen.gender <- gender
			}
		}

		# sec.link.stat.x.c processes a matrix of expression, returns a vector.
		sec.stat <- sec.link.stat.x.c(exp[-i,], locexp, genotype, gender = gen.gender)
		# As with sec.link.stat, process a matrix, return a vector. 
		ind.stat <- condi.indep.stat.x.c(exp[-i,], locexp, genotype, gender = gen.gender)
		


		### CALCULATE NULL STATISTICS ###
		## first, sample index, to produce a index matrix of #.of.exp by B. 
		#  Each time to generate null statistics for a locexp, 
		## one can use the same sets of index.
        for (j in 1:Bsec){
		    secexp0 <- (exp0.list[[j]])[-i,]
		    locexp0 <- sample((exp0.list[[j]])[i,])
		    sec.stat0 <- c(sec.stat0, sec.link.stat.x.c(secexp0, locexp, genotype, gender = gen.gender))
		    ind.stat0 <- c(ind.stat0, condi.indep.stat.x.c(secexp0, locexp0, genotype, gender = gen.gender))
	    }
		ind.stat0.nosort <- ind.stat0
		ind.stat0 <- sort(ind.stat0)

		tstat <- sec.stat
		sec.pval <- edge.pvalue(tstat, sec.stat0)
		pi0.sec <- edge.qvalue(sec.pval)$pi0
		np.prob.sec <- rep(0,length(tstat))
		np.prob.sec <- 1 - edge.lfdr(sec.pval, pi0 = pi0.sec)
		if (pi0.sec < .98) {
				istat <- ind.stat
				isorder <- order.c(istat)
				pind <- edge.pvalue(istat, ind.stat0)
				ttcut <- sort(tstat, decreasing = T)[round(length(tstat)*(1-pi0.sec))]
				index <- which(tstat>= ttcut)
				pi0.ind <- edge.qvalue(pind[index], lambda = 0)$pi0
				
			} 
			else {
				pi0.ind <- 0 
			}
	  
			np.prob.ind <- rep(0, length(tstat))
			
			if (pi0.ind > 0) {
			    FF <- istat[index]
				FF0 <- sort(sample(ind.stat0, Bsec*length(index)))
				#FF[FF>max(FF0)] <- max(FF0)
				ind.pval = edge.pvalue(FF, FF0)
				
				np.prob.ind[index] <- edge.lfdr(ind.pval, pi0 = pi0.ind, lambda = 0)
			}
			
       
			jointp <- prob.loc[i]*np.prob.sec*np.prob.ind
			jointp[jointp<0] <- 0
			
			if (i>1) jp[1:(i-1)] <- jointp[1:(i-1)]
		    if (i < nr) jp[(i+1):nr] <- jointp[i:(nr-1)]  
		} 
	        cat(c(jp,"\n"),file = ofile,append = T) 
	        #prob <- rbind(prob, jp)
		
		# Keep track of our progress in the log file.
		oon <- i
		if (floor(oon/round(ni/10)) ==  count & count <= 9){

			print(paste(count*10,"% completed", collapse = "", sep = ""),quote = FALSE)
			count <- count +1
		
        }
		if (oon ==  ni){
			print(paste(100,"% completed", collapse = "", sep = ""),quote = FALSE)
		}

	}
	if(length(idx)>1) {
		prob <- matrix(scan(ofile),byrow = TRUE, ncol  = nrow(exp))
	} else {
		prob <- matrix(jp, nrow = 1)
	}
        #file.remove(ofile)
	rownames(prob) <- gnams
	colnames(prob) <- rownames(exp)
	prob <- round(prob, 4)
	return(prob)
}
setGeneric("trigger.net",function(triggerobj, gender = NULL, idx = NULL, Bsec = 100, prob.cut = 0.7,  include.loc = TRUE, seed = 123, inputfile = NULL){standardGeneric("trigger.net")})
setMethod("trigger.net", signature = "trigger", net.trigger)

netplot2ps.trigger <- function(triggerobj, trig.prob, filenam = NULL, pcut = 0.95, layout = c("radial", "energy-minimized", "circular","hierarchical"), node.color = NULL, edge.color = NULL, node.shape = NULL, nreg = 20){
	#trig.prob = triggerobj@trig.prob
    nam.r <- colnames(trig.prob)
	nam.c <- colnames(trig.prob)
	nam.r = gsub("-","", nam.r)
	nam.c = gsub("-","", nam.c)
	prob <- trig.prob*0
	prob[trig.prob>= pcut] <- 1
	nedge <- sum(prob == 1)
	if (nedge>= 1000) {
		shape <- "point"
	}else{
		if (is.null(node.shape))	{
			shape <- "box"
		} else {
			shape <- node.shape
		}
	}
	
	if (is.null(filenam)) {
		dotfilenam <- "temp"
	} else{
		dotfilenam <- filenam
	}
	rsprob <- rowSums(prob)
	reg.idx <- which(rsprob>0)
	od <- order(rsprob, decreasing = T)
	
	if(nreg>0) main.reg <- nam.r[od[1:min(nreg, length(nam.r))]]
	
	
    if (is.null(node.color)) node.color <- "green"
	if (is.null(edge.color)) edge.color <- "blue"
	dotfile <- paste(c(dotfilenam,".dot"),sep = "",collapse = "")
	cat(c("digraph g {", "\n"), file = dotfile, append = T, sep = "")
	cat(c("size = ", shQuote("7.5,11", type = "cmd"),";", "\n"), file = dotfile, append = T, sep = "")
	cat(c("node[shape = ", shape,", color = ", node.color, "];", "\n"), file = dotfile, append = T, sep = "")
	if (nreg>0) {cat(c(paste(main.reg, "[shape = ellipse, color = red];"), "\n"), file = dotfile, append = T, sep = "")}
	cat(c("edge[arrowsize = 0.5, minlen = 0.1,style = filled, color  = ", edge.color, "];", "\n"), file = dotfile, append = T, sep = "")
	cat(c("graph[center = true];", "\n"), file = dotfile, append = T, sep = "")
	cat(c("ratio = file;", "\n"), file = dotfile, append = T, sep = "")

	od2 <- order(rowSums(prob[reg.idx,]), decreasing = F)
	for (i in reg.idx[od2]){
	   regee.idx <- which(prob[i,] == 1)
	   cat(c(nam.r[i], " -> {"), file = dotfile, append = T, sep = "")
	   for (j in regee.idx){
	      cat(c(nam.c[j], ";"), file = dotfile, append = T, sep = "")
	   }
	   cat(c("}", "\n"), file = dotfile, append = T, sep = "")
	}
	cat(c("}", "\n"), file = dotfile, append = T, sep = "")
	layout <- substr(layout[1],1,1)
	if (layout == "r"|layout == "R") cmd.mode <- "twopi -Tps "
	if (layout == "e"|layout == "E") cmd.mode <- "neato -Tps -Gmaxiter = 1000 "
	if (layout == "c"|layout == "C") cmd.mode <- "circo -Tps "
	if (layout == "h"|layout == "H") cmd.mode <- "dot -Tps "
	cmd <- paste(c(cmd.mode, dotfile, " -o ", dotfilenam, ".ps &"), sep = "", collapse = "")
	
    system(cmd)
}
setGeneric("trigger.netPlot2ps",function(triggerobj, trig.prob, filenam = NULL, pcut = 0.95, layout = c("radial", "energy-minimized", "circular","hierarchical"), node.color = NULL, edge.color = NULL, node.shape = NULL, nreg = 20){standardGeneric("trigger.netPlot2ps")})
setMethod("trigger.netPlot2ps", signature = "trigger", netplot2ps.trigger)

exportdata.trigger <- function(triggerobj, plotarg = TRUE, verbose = TRUE, warning = FALSE){
	exp <- triggerobj@exp
	marker <- triggerobj@marker
	marker.pos <- triggerobj@marker.pos
	
	traitmat = cbind(rownames(exp), "","", exp)
	colnames(traitmat) = c()
	genomat = cbind(rownames(marker), marker.pos, marker)
	colnames(genomat) = c()
	
	if(verbose){cat("Writing genotype and phenotype data to file", "\n")}
	write.table(traitmat, file = "geno_trait_data.csv", quote = F,sep = ",", row.names = F, col.names = F)
	write.table(genomat, file = "geno_trait_data.csv", quote = F, row.names = F, sep = ",", append = T, col.names = F)
	if(verbose){cat("Importing genotype and phenotype data", "\n")}
	if (warning) crossfile = read.cross("csvr", ".", "geno_trait_data.csv", genotypes = c(1,2))
	if (!warning) crossfile = suppressWarnings(read.cross("csvr", ".", "geno_trait_data.csv", genotypes = c(1,2)))
	crossfile = jittermap(rescalemap(crossfile, scale = 5e-4))
	if(plotarg){plot.map(crossfile)}
	crossfile <- calc.genoprob(crossfile, step = 1)
	return(crossfile)
	
}
setGeneric("trigger.export2cross", function(triggerobj, plotarg = T, verbose = T, warning = FALSE) {standardGeneric("trigger.export2cross")})
setMethod("trigger.export2cross", signature = "trigger", exportdata.trigger)

trait.trigger <- function(triggerobj, trait, cross = NULL, thr = 3, n.sv = NULL, addplot = TRUE){
	exp <- triggerobj@exp
	exp.pos <- triggerobj@exp.pos
	marker <- triggerobj@marker
	marker.pos <- triggerobj@marker.pos
	genes = rownames(exp)
	n.ind = ncol(marker)
	out.em.s = chr.s = cisgenes = marker.max = proxx.gene = t1.fit = t2.fit = 	cis.list = result = search.marker = resids = cis.find = c()	
	
	if(is.null(cross)) cross = trigger.export2cross(triggerobj)
	
	if(is.character(trait)){
		pheno = which(genes ==  trait); traitentry = F
	} else if(is.vector(trait)){
		pheno.vec = matrix(trait, nrow = length(trait)/n.ind, byrow = T)
		if (ncol(pheno.vec)!=  nind(cross)) stop("Number of entries in trait matrix should be the same as that of genotype matrix")
		cross$pheno = data.frame(t(pheno.vec)); pheno = 1; traitentry = T
	} else cat(stop("Trait should either be a gene-name or a numeric vector"))
	
	cross <- calc.genoprob(cross, step = 1)
	out.em <- scanone(cross, pheno.col = pheno, method = "em")
	out.em.s = summary(out.em, thr = thr)
	n.sig = length(unique(out.em.s$chr))
	chr.sig = as.numeric(unique(out.em.s$chr))
	
	if(n.sig ==  0){
		stop("No QTL found for the trait")
	} else {
		if(addplot & !traitentry) plot(out.em, chr = chr.sig, main = paste(genes[pheno],":chr", chr.sig, sep = ""), bandcol = "gray70", alternate.chrid = T)
		if(addplot & traitentry) plot(out.em, chr = chr.sig, main = paste("chr", chr.sig, sep = ""), bandcol = "gray70", alternate.chrid = T)
		result = c()
		for(i in 1:n.sig){
			chr.s = chr.sig[i]
			cisgenes = find.cis.genes(out.em, cross, chr.s, marker, marker.pos, exp, exp.pos, thr = thr)
			marker.max = cisgenes$marker.max
			prox.genes = cisgenes$prox.genes
			traitmatch = which(prox.genes ==  pheno)
			if(length(traitmatch)>0){
				prox.genes = prox.genes[-traitmatch]
			}
			if (length(prox.genes) == 0) stop("No putative causal regulator found in proximity of the QTL")
			
			if(is.null(n.sv)) {
				nsv = 0
				} else{
					nsv = n.sv
				}
			#resids = tryCatch(sva.fit(exp = exp, qtl = marker[marker.max, ], traitid = pheno, nsv = nsv), error = function(err){1})
			resids = sva.fit(exp = exp, qtl = marker[marker.max, ], traitid = pheno, nsv = nsv)
			#if(length(resids) ==  1) stop("Caught SVD error", "\n")
			if(length(resids) == 3){
				t2.fit = resids$t2.fit
				t1.fit = resids$t1.fit
				n.sv = resids$n.sv
				df0 = ncol(t1.fit) - resids$n.sv - 1
	
				cis.list = cbind(marker.max, prox.genes, genes[prox.genes])
				cis.list = matrix(cis.list, nrow = length(prox.genes))
				cat(paste("Fitting", nrow(cis.list), "genes on chromosome", chr.sig[i], sep = " "), "\n")
	
				result = c(result, traitmap.fun(t1.fit, t2.fit, cis.list , marker, df0 = df0, B = 5))
					
			}	
				
		}
	}
	return(result)
}
setGeneric("trigger.trait", function(triggerobj, trait, cross = NULL, thr = 3, n.sv = NULL, addplot = TRUE) {standardGeneric("trigger.trait")})
setMethod("trigger.trait", signature = "trigger", trait.trigger)

