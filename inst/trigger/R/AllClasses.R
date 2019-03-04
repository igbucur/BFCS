setClass("trigger",
		representation=representation(
		marker = "matrix", 
		exp = "matrix",
		marker.pos = "data.frame",
		exp.pos = "data.frame",
		stat = "matrix",
		pvalue = "matrix",
		eqtl.R2 = "numeric",
		loc.obj = "list",
		mlink = "list"
	), 
	prototype = prototype(
		marker = cbind(matrix(0,10,5), matrix(1,10,5)),
		marker.pos = data.frame(chr = c(1:10), pos = sample(c(1:1e4),10)),
		exp = matrix(rnorm(100),10,10), 
		exp.pos = data.frame(chr = c(1:10), start = sample(c(1:1e4),10), end = sample(c(1:1e4),10))
		),
	
	validity=function(object){
			if(nrow(marker) != ncol(marker.pos)){
				"Check marker matrix dimensions"}
			else
			if(nrow(exp) != nrow(exp.g)){
				"Check expression matrix dimensions"}
			else
			if(ncol(exp) != ncol(exp.g)){
				"Check marker and expression matrix dimensions"}
			else
				{return(TRUE)}
					}
)

	
