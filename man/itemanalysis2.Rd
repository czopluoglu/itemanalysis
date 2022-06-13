\name{itemanalysis2}
\alias{itemanalysis2}
\title{
Classical Test Theory Item Analysis for Polytomous Items
}
\description{
Classical Test Theory Item Analysis for Polytomous Items
}
\usage{
itemanalysis2(data, options, ngroup = ncol(data) + 1, correction = TRUE,
span.par=.3,verbose=T)
}
\arguments{
  \item{data}{
a data frame with \emph{N} rows and \emph{m} columns, with \emph{N} denoting the number of subjects and \emph{m} denoting the number of items.
}
  \item{options}{
a vector of numerical code of the response categories available for the items such as c(0,1,2,3). The minumum score is assumed to be 0.
}
  \item{ngroup}{
number of score groups to be use for plotting the item trace lines
}
  \item{correction}{
TRUE or FALSE. If it is TRUE, then an adjustment is made for point-biserial correlation. 
}
  \item{span.par}{
 a smoothing parameter to pass to ggplots when creating empirical ICCs  
}
  \item{verbose}{
TRUE or FALSE. If it is TRUE, text output is printed.
}
}
\details{
to be added later
}

\value{
    \item{plots}{a list object storing the item trace line plots for each item}
    \item{item.stat}{a matrix of basic item statistics}
    \item{dist.sel}{a matrix of distractor selection proportion statistics}
    \item{dist.disc}{a matrix of corrected point-biserial statistics for distractors}
    \item{dist.disc}{a matrix of corrected biserial statistics for distractors}
}

\author{
Cengiz Zopluoglu
}

\seealso{
\code{\link{itemanalysis1}} for classical item analysis of multiple-choice test items
}
\examples{

 \dontrun{
    data(timss2011_usa)

    timss2011_usa$Q14B <- recode(var = timss2011_usa$Q14B,
                                 recodes = "c(0)=3;c(1)=2;c(2)=1;c(3)=0")

    timss2011_usa$Q14C <- recode(var = timss2011_usa$Q14C,
                                 recodes = "c(0)=3;c(1)=2;c(2)=1;c(3)=0")

    item.analysis <- itemanalysis2(data=timss2011_usa,
                                   options=c(0,1,2,3),
                                   ngroup=18,
                                   correction=FALSE)
                                   
    item.analysis$item.stat
      
    item.analysis$dist.sel
      
    item.analysis$dist.disc

    item.analysis$plots[[1]]   # Item Trace Line for the first item
    item.analysis$plots[[2]]   # Item Trace Line for the second item
    item.analysis$plots[[3]]   # Item Trace Line for the third item
    item.analysis$plots[[4]]   # Item Trace Line for the fourth item
    item.analysis$plots[[5]]   # Item Trace Line for the fifth item
    item.analysis$plots[[6]]   # Item Trace Line for the sixth item
  }
}
