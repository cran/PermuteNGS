\name{PermuteDEG}
\alias{PermuteDEG}
\title{PremuteDEG in PermuteNGS package}
\description{DEG test with un-replicated NGS expression data }
\usage{PermuteDEG(dataA, dataB, samplingNumber, alternative ) }
\arguments{
  \item{dataA}{ One of Input data, file format is raw read count as gene}
  \item{dataB}{Other Input data, file format is raw read count as gene}
  \item{samplingNumber}{ sampling Number from pooled set}
  \item{alternative}{ one-tail analysis : alternative= "one.sided", two-tail analysis : alternative="two.sided" }
}
\details{If necessary more details than the description above, see the manual or manuscprit, please. }
\references{ 
  Taeheon Lee and BO-Young Lee, et al (2010), PermuteNGS: The significance testing of transcriptome profiling for RNA-sequencing data 

  Benjamini,Y. and Hochberg,Y. (1995) Controlling the false discovery rate: a practical and
  powerful approach to multiple testing. \emph{J. R. Stat. Soc. Ser. B} \bold{57}, 289-300.

  Storey,J.D. and Tibshirani,R. (2003) Statistical significance for genomewide studies. \emph{Proc. Natl. Acad. Sci.} \bold{100}, 9440-9445.
}
\author{Taeheon Lee, Bo-Young Lee, Hee-Seok Oh,Heebal Kim}
\examples{
library(PermuteNGS)
data(Sample250)
attach(Sample250)
result <- PermuteDEG(R1L1Kidney, R1L2Liver, 1000, alternative = "two.sided")
rownames(result)= rownames(data)
write.csv(result, "PermuteNGS.result.csv")
}
\keyword{htest}%
