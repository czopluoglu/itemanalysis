itemanalysis1 <- function (data, key, options,ngroup=ncol(data)+1,correction=TRUE,span.par=.3, verbose = T) {
  
  #########################################################################################
  # data, a data frame with N rows and n columns, where N denotes the number of subjects 
  #        and n denotes the number of items. All items should be scored using nominal 
  #        response categories. All variables (columns) must be "character". 
  #        Missing values ("NA") are allowed and scored as incorrect in item analysis
  
  # key,     a character vector of length n, where n denotes the number of items. 
  
  # options, number of possible nominal options for items (e.g., "A","B","C","D")
  #          make sure each item is consistent, and includes the same response options
  
  # ngroup, number of score groups
  #
  # correction, TRUE or FALSE, if TRUE item and distractor discrimination is corrected for
  # 		  spuriousnes by removing the item score from the total score
  #
  # span.par, this is a smoothing parameter to pass to ggplots when creating empirical ICCs
  #########################################################################################
  
  
  # Checks whether or not the columns are characters
  # Make it character if not.
  
  for (i in 1:ncol(data)) {
    if (is.character(data[, i]) != TRUE) {
      data[,i]=as.character(data[,i])
    }
  }
  
  # Compare each column to the key response for that column
  # the order of key responses should be aligned with the columns in the dataset
  
  scored.data <- as.data.frame(matrix(nrow = nrow(data), ncol = ncol(data)))
  
  for (i in 1:ncol(scored.data)) {
    
    scored.data[,i] <- ifelse(data[,i] == key[i],1,0)
    
    if(length(which(is.na(scored.data[,i])))!=0) {
      scored.data[which(is.na(scored.data[,i])==TRUE),i]=0
    }
  }
  
  total.score <- rowSums(scored.data)
  ybar <- mean(total.score)
  sdt <- sd(total.score)
  p <- colMeans(scored.data)
  
  pbis <- c()
  pbis.corrected <- c()
  bis  <- c()
  bis.corrected <- c()
  
  for(k in 1:ncol(data)) {  
    pbis[k]=cor(scored.data[,k],total.score,use="pairwise.complete.obs")
    pbis.corrected[k]=cor(scored.data[,k],
                          rowMeans(scored.data[,-k],na.rm=TRUE)*(ncol(scored.data)-1),
                          use="pairwise.complete.obs")
    bis[k]=polyserial(total.score,scored.data[,k])
    bis.corrected[k]=polyserial(rowMeans(scored.data[,-k],na.rm=TRUE)*(ncol(scored.data)-1),scored.data[,k])
  }
  
  
  item.stat <- matrix(nrow=ncol(data),ncol=4)
  colnames(item.stat) <- c("Item Difficulty","Item Threshold","Point-Biserial","Biserial")
  
  rownames(item.stat) <- colnames(data)
  
  item.stat[,1]=p
  item.stat[,2]=qnorm(1-p)
  if(correction==TRUE){ item.stat[,3]=pbis.corrected } else { item.stat[,3]=pbis }
  if(correction==TRUE){ item.stat[,4]=bis.corrected} else {item.stat[,4]=bis}
  
  
  #    eff.size           <- matrix(nrow=ncol(data),ncol=4)
  #    colnames(eff.size) <- c("Incorrect","Correct","p","d")
  #    for(i in 1:ncol(data)){
  #      gr1 = rowSums(scored.data[scored.data[,i]==0,])
  #      gr2 = rowSums(scored.data[scored.data[,i]==1,])
  #      t = t.test(gr1,gr2)
  #      eff.size[i,1] = mean(gr1)
  #      eff.size[i,2] = mean(gr2)
  #      eff.size[i,3] = t$p.value
  #      eff.size[i,4] =  abs(cohen.d(gr1,gr2)$estimate)
  #    }
  
  #    eff.size <- as.data.frame(eff.size)
  #    
  #    eff.size$p2 <- ifelse(eff.size$p<.001,"p<.001",
  #                          ifelse(eff.size$p>.001 & eff.size$p <.01,"p<.01",
  #                                 ifelse(eff.size$p > .01 & eff.size$p <.05,"p<.05","Not significant")))
  #    
  #    eff.size <- eff.size[,c(1,2,5,4)]
  #    eff.size[,1]=round(eff.size[,1],2)
  #    eff.size[,2]=round(eff.size[,2],2)
  #    eff.size[,4]=round(eff.size[,4],2)
  #    colnames(eff.size) <- c("Incorrect","Correct","P.value","Cohen.d")
  
  
  # Create the score groups with equal width based on the distribution of total score
  # and number of groups
  
  sgroups <- cut(total.score,breaks=ngroup)
  slevels <- levels(sgroups)
  
  sgnum <- rowMeans(cbind(lower = as.numeric( sub("\\((.+),.*", "\\1", slevels) ),
                          upper = as.numeric( sub("[^,]*,([^]]*)\\]","\\1",slevels))))
  
  
  
  SG <- vector("list",ngroup)
  
  for(j in 1:ngroup){
    SG[[j]]=which(sgroups==slevels[j])
  }
  
  # Compute the proportion of selecting each response option within 
  # each score group 
  
  prop <- vector("list",ncol(data))
  names(prop) <- colnames(data)
  
  for(i in 1:ncol(data)) {
    
    dist <- matrix(nrow=length(options),ncol=ngroup)
    colnames(dist) <- slevels
    rownames(dist) <- options 
    
    for(g in 1:ngroup){
      for(o in 1:length(options)){
        dist[o,g]=length(which(data[SG[[g]],i]==options[o]))/length(SG[[g]])
      }
    }
    
    prop[[i]]=dist
    
  }
  
  dist.sel <- matrix(nrow=ncol(data),ncol=length(options))  
  dist.disc <- matrix(nrow=ncol(data),ncol=length(options))
  dist.disc2 <- matrix(nrow=ncol(data),ncol=length(options))
  colnames(dist.disc) <- options
  rownames(dist.disc) <- colnames(data)
  colnames(dist.disc2) <- options
  rownames(dist.disc2) <- colnames(data)
  colnames(dist.sel) <- options
  rownames(dist.sel) <- colnames(data)
  
  for(i in 1:ncol(data)){
    for(o in 1:length(options)) {
      temp <- ifelse(data[,i]==options[o],1,0)
      temp[is.na(temp)]=0
      dist.sel[i,o]=mean(temp,na.rm=TRUE)
      if(correction==FALSE){
        dist.disc[i,o]=cor(temp,total.score,use="pairwise.complete.obs")
        dist.disc2[i,o]=polyserial(total.score,temp)
      } else {
        dist.disc[i,o]=cor(temp,rowMeans(scored.data[,-i],na.rm=TRUE)*(ncol(scored.data)-1),use="pairwise.complete.obs")
        dist.disc2[i,o]=polyserial(rowMeans(scored.data[,-i],na.rm=TRUE)*(ncol(scored.data)-1),temp)
      }
    }
  }
  dist.sel <- as.data.frame(dist.sel)
  dist.sel$Missing <- 1-rowSums(dist.sel)
  
  plots <- vector("list",ncol(data))
  
  for(i in 1:ncol(data)) {
    
    options.d <- c()
    for(u in 1:length(options)){ 
      options.d[u] <- paste(options[u],"( ",round(dist.disc2[i,u],2)," )",sep="")
    }
    
    d <- as.data.frame(cbind(sg=sgnum,p=prop[[i]][1,]))
    for(u in 2:length(options)){ d <- rbind(d,cbind(sg=sgnum,p=prop[[i]][u,]))}
    optt <- c()
    for(u in 1:length(options)){ optt <- c(optt,rep(options.d[u],ngroup))}
    d$opt <- optt
    
    
    pp <- ggplot(data=d,aes_string(x="sg",y="p",group="opt",shape="opt"))+
      geom_point(size=2)+
      geom_smooth(span=span.par,col="black",lwd=.5)+
      ggtitle(colnames(data)[i])+
      theme(panel.background = element_blank(),legend.title=element_blank(),legend.key = element_blank())+
      theme(legend.justification=c(0,1),legend.position=c(0,1),legend.text=element_text(size=12,face="bold"))+
      scale_x_continuous(limits = c(0,ncol(data)),breaks=seq(0,ncol(data),ceiling(ncol(data)/10)))+
      scale_y_continuous(limits = c(0,1))+xlab("Score Groups")+ylab("Proporion of Being Selected")
    
    plots[[i]] <- pp
  }
  
  ###############################################################
  if (verbose==T){
    
    
    cat("************************************************************************","\n")
    cat("itemanalysis: An R package for Classical Test Theory Item Analysis","\n")
    cat("","\n")
    cat("Cengiz Zopluoglu","\n")
    cat("","\n")
    cat("University of Oregon","\n")
    cat("College of Education","\n")
    cat("","\n")
    cat("cen.zop@gmail.com","\n")
    cat("","\n")
    cat("Please report any programming bug or problem you experience to improve the code.","\n")
    cat("*************************************************************************","\n")
    
    cat("Processing Date: ",date(),"\n\n")
    
    cat(sprintf("%50s","ITEM STATISTICS"),"\n")
    cat("","\n")
    print(round(item.stat,3))
    cat("","\n")
    cat("","\n")
    cat("","\n")
    
    cat(sprintf("%50s","DISTRACTOR SELECTION PROPORTIONS"),"\n")
    cat("","\n")
    print(round(dist.sel,3))
    cat("","\n")
    cat("","\n")
    cat("","\n")
    
    cat(sprintf("%50s","DISTRACTOR Point-Biserial"),"\n")
    cat("","\n")
    print(round(dist.disc,3))
    cat("","\n")
    cat("","\n")
    cat("","\n")
    
    cat(sprintf("%50s","DISTRACTOR Biserial"),"\n")
    cat("","\n")
    print(round(dist.disc2,3))	
    cat("","\n")
    cat("","\n")
    cat("","\n")
    
  } else {
    
  }
  
  return(list(item.stat=item.stat,
              dist.sel=dist.sel,
              dist.disc=dist.disc,
              dist.disc2=dist.disc2,
              plots=plots))
  
  
}