getClustProfiles <- function(expDataIn,clustIn,clusterCutoff){
  numClusters <- ncol(clustIn$membership)
  clustProfiles <- matrix(data=0,ncol=ncol(expDataIn),nrow=numClusters)
  for(i in 1:numClusters){
    SampleInClust <- c(1:nrow(expDataIn))[clustIn$membership[,i]>=clusterCutoff]
    if(length(SampleInClust)>0){
      if(length(SampleInClust)>1){
        clustProfiles[i,]<-as.numeric(apply(expDataIn[SampleInClust,],2,median))
      }else{
        clustProfiles[i,]<-as.numeric(expDataIn[SampleInClust,])
      }
    }
  }
  whichCluster <- c(1:numClusters)[apply(clustProfiles,1,var)>0]
  clustProfiles <- clustProfiles[apply(clustProfiles,1,var)>0,]
  return(list(profiles=clustProfiles,whichClust=whichCluster))
}

makeCollapsedProfiles <- function(expDataIn,clustProfsIn, clusterCutoff,grouping,origClust,clustIn){                                                                   
  collapsedProfiles <- matrix(data=0,ncol=ncol(expDataIn),nrow=max(grouping))
  for(i in 1:max(grouping)){
    inClust <- c(1:length(grouping))[grouping==i]
    if(length(inClust)==1){
      collapsedProfiles[i,]<-clustProfsIn[inClust,]
      SampleInClust <- c(1:nrow(expDataIn))[clustIn$membership[,origClust[inClust]]>=clusterCutoff]                                                                     
    }else{
      clustLookAt <- origClust[inClust]
      SampleInClust <- c(1:nrow(expDataIn))[apply(clustIn$membership[,clustLookAt],1,max)>= clusterCutoff]
      collapsedProfiles[i,]<-apply(expDataIn[SampleInClust,],2,median)
    }
  }
  return(collapsedProfiles)
}
correlateSampleToProfiles <- function(expDataIn,profilesIn){
  profileCor <- matrix(data=0,ncol=nrow(profilesIn),nrow=nrow(expDataIn))
  for(i in 1:nrow(profilesIn)){
    profileCor[,i] <- apply(expDataIn,1,cor,y=profilesIn[i,])
  }
  rownames(profileCor)<-rownames(expDataIn)
  colnames(profileCor)<-rownames(profilesIn)
  return(profileCor)
}

GetMemberList <- function(patternCor,minDist){
  tmpSample <- row.names(patternCor)
  tmpList <- list()[colnames(patternCor)]
  for(i in 1:nrow(patternCor)){
    tmpIndex = which.max(patternCor[i,])
    if(patternCor[i,tmpIndex] > minDist){
      tmpList[[tmpIndex]] <- c(tmpList[[tmpIndex]],tmpSample[i])
    }
  }
  tmpMat <- matrix(NA,nrow = length(unlist(tmpList)),ncol = 2)
  n = 1
  for(i in 1:length(tmpList)){
    tmpSamplet <- tmpList[[i]]
    for(j in 1:length(tmpSamplet)){
      tmpMat[n,1] <- tmpSamplet[j]
      tmpMat[n,2] <- paste("Pattern",i,sep = "_")
      n = n+1
    }
  }
  return(tmpMat)
}