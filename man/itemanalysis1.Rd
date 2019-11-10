\name{itemanalysis1}
\alias{itemanalysis1}
\title{
Classical Test Theory Item Analysis for Multiple-Choice Test Items
}
\description{
Classicial Test Theory Item Analysis for Multiple-Choice Test Items
}
\usage{
itemanalysis1(data, key, options, ngroup = ncol(data) + 1, correction = TRUE,span.par=.3)
}
\arguments{
  \item{data}{
a data frame with \emph{N} rows and \emph{m} columns, with \emph{N} denoting the number of subjects and \emph{m} denoting the number of items.
}
  \item{key}{
a vector of answer key with a length of \emph{m}. The order of elements in this vector must match with the order of columns in the input dataset
}
  \item{options}{
a vector of response options in the input dataset, e.g. c("A","B","C","D")
}
  \item{ngroup}{
number of score groups to be use for plotting the empirical item trace lines
}
  \item{correction}{
TRUE or FALSE. If it is TRUE, then an adjustment is made for point-biserial correlation.
}
  \item{span.par}{
Smoothing parameter to pass to geom_smooth() in ggplot2 for drawing lines in empirical item characteristic curves.
}

}
\details{

This function computes the traditional statistics for the classical test theory item analysis for item response data that can be scored as correct/incorrect. Missing data is treated as incorrect. The score groups are created based on equal width intervals depending on the total score range and number of requested score groups.

}
\value{

    \item{item.stat}{a matrix with overall item statistics}
    \item{dist.sel}{a matrix with distractor selection proportions}
    \item{dist.disc}{a matrix with distractor point-biserial correlations}
    \item{dist.disc2}{a matrix with distractor biserial correlations}
    \item{plots}{a list object storing the item trace line plots for each item}
}

\author{
Cengiz Zopluoglu
}

\seealso{
\code{\link{itemanalysis2}} for classical item analysis of polytomously scored items
}
\examples{

    data(dichotomous)
    head(dichotomous)
    str(dichotomous)

    # Key response vector

    key <- c("A","D","C","B","C","B","C","D","A","D","C","A","D","C","A",
            "B","D","B","A","C","A","A","C","B","C","B","D","A","A","A",
            "C","B","B","A","B","D","D","A","D","C","D","A","B","B","C",
            "D","B","C","C","B","D","A","C","B","A","D")

    # Use itemanalysis1 function to run the item analysis

      # In order to reduce running time for the example below,
      # I specify "data=dichotomous[,1:10]", so it only analyze the 
      # first 10 items.
      # You should specify "data=dichotomous" to analyze based on 56 items.

    item.analysis <- itemanalysis1(data=dichotomous[,1:10],
                           key=key[1:10],
                           options=c("A","B","C","D"),
                           ngroup=10,
                           correction=FALSE,
                           span.par = .3)
                           
    # item.analysis$item.stat    # Item statistics 
    
    # item.analysis$dist.sel     # Distractor Selection Proportions
    
    # item.analysis$dist.disc    # Distractor point-biserial correlations
    
    # item.analysis$dist.disc2   # Distractor biserial correlations

    # item.analysis$plots[[1]]   # Item Trace Line for the first item
    # item.analysis$plots[[2]]   # Item Trace Line for the second item
    # item.analysis$plots[[3]]   # Item Trace Line for the third item
    # item.analysis$plots[[4]]   # Item Trace Line for the fourth item
    # item.analysis$plots[[5]]   # Item Trace Line for the fifth item
    # item.analysis$plots[[6]]   # Item Trace Line for the sixth item
    # item.analysis$plots[[7]]   # Item Trace Line for the seventh item
    # item.analysis$plots[[8]]   # Item Trace Line for the eigth item
    # item.analysis$plots[[9]]   # Item Trace Line for the ninth item
    # item.analysis$plots[[10]]  # Item Trace Line for the tenth item

}


