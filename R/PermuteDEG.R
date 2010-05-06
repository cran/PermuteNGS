PermuteDEG <-  function(dataA,dataB,samplingNumber, alternative ) {
                           numRow = length(dataA)
                           #z_mine = length(dataA)
                           #pvalue = length(dataA)
                           if (alternative == "one.sided"){
                               alternative = 1 }
                           else if (alternative == "two.sided"){
                               alternative =2 }
                           Pvalue = .Call('PermuteDEG', 
                           PACKAGE='PermuteNGS',
	                   as.integer(dataA), 
		           as.integer(dataB), 
                           as.integer(numRow),
   	                   as.integer(samplingNumber),
                           as.integer(alternative)
        ) 

	Qvalue_Benjamini_1995 <-p.adjust( Pvalue ,'BH')
	
	getQvalue_2 <- function(Pvalue_v){
		Qvalue_2 <- rep(NA,length(Pvalue_v))
		tmp <- try(library(qvalue))
		if(class(tmp) !="try-error") {
		res <- qvalue(Pvalue_v)
		Qvalue_2<- res$qvalues
		}
		else{
		cat("Warning:fail to load library qvalue(Storey et al. 2003)!\n")
		packages_tmp <-try(install.packages("qvalue"))
		if(class(packages_tmp) != "try-error"){ 
		cat("qvalue library installed, try agrain  \n")
		tmp <- try(library(qvalue))
		if(class(tmp) !="try-error") {
		res <- qvalue(Pvalue_v)
		Qvalue_Storey_2003 <- res$qvalues
		}
	}
      }
      Qvalue_2
}

tmp2 <- try(Qvalue_Storey_2003 <- getQvalue_2(Pvalue))
if(class(tmp2) !="try-error") {
                Qvalue_Storey_2003 <- getQvalue_2(Pvalue)
return(data.frame(dataA,dataB,Pvalue, Qvalue_Benjamini_1995)) #, Qvalue_Storey_2003))
}
else{
return(data.frame(dataA,dataB,Pvalue, Qvalue_Benjamini_1995))}
}
