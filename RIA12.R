#-------------------------------------------------------------------------------------------------------------------------
#
#																		      ANALYSE IMPLICATIVE DE RASCH
#
#																			Yvonnick Noël, Univ Rennes, 2015-2020
#      VERSION AVEC :
#       - crossprod parallélisé
#       - suppression minimales des colonnes/lignes à distances manquantes dans D
#       - uniquement du MDS numérique
#       - corrige un bug sur le rangement des variables par réussite (ne prenait pas en compte les valeurs manquantes)
#       - permet de choisir le type de modèles de classes ("EEV","VEE",...) en mélange de gaussiennes (package Mclust)
#
#-------------------------------------------------------------------------------------------------------------------------

library(R6)
library(mclust,quietly=TRUE)

RIA = R6Class("RIA",
  
  # Public members
  public = list(
  
    data = NULL,
    N = 0,
    p = 0,
    coord = NULL,
    true.coord = NULL,
    alpha = NULL,
    GoF = Inf,
    eig = NULL,
    D = NULL,
    weights = NULL,
    clustfit = NULL,
    mdsfit = NULL,
    Iota = NULL,
    scales = NULL,
    ndim = 2,
    ncores = NULL,
    cls = NULL,
    
    initialize = function(X=NULL,N=300,p=50,outliers=NULL,scale=1,theta=NULL,sort.theta=FALSE,slopes=FALSE,ndim=3,missing.prob=0.0) {

      self$N = N
      self$p = p
      
      # No data provided: Simulate a Rasch
      if(is.null(X)) {
        if(is.null(theta)) {
          theta = rnorm(N)
          if(sort.theta) theta = sort(theta)
        }
        else {
          N = length(theta)
          self$N = N
        }
        delta = seq(-2,2,len=p)
        self$true.coord = delta
        D = outer(theta,delta,"-")*scale
        
        if(slopes) {
          self$alpha = rlnorm(p)
          D = D %*% diag(self$alpha)
        }

        pi = exp(D)/(1+exp(D))
        self$data = (matrix(runif(N*p),N,p) < pi) + 0.0
        self$N = nrow(self$data)
        self$p = ncol(self$data)

        # Add missing values
        if(missing.prob > 0) {
          W = matrix(runif(prod(dim(self$data))),self$N,self$p) < missing.prob
          self$data[W] = NA
        }
        # Add outliers, i.e. uncorrelated items with given easiness, for instances c(.25,.5,.75)
        if(length(outliers) && is.numeric(outliers))  {
          for(j in 1:length(outliers)) {
            self$data = cbind(self$data,rbinom(N,1,outliers[j]))
            self$p = self$p + 1
          }
        }

        colnames(self$data) = paste("item",1:self$p,sep="")
      }

      else {
        self$data = X
        self$N = nrow(X)
        self$p = ncol(X)
      }

      self$ndim = ndim
      if(is.null(ndim)) self$ndim = self$p/2

    },
    setDim = function(n) {
      self$ndim = n
    },
    #---------------------------------------------- Run estimation
    run = function(clustModel="VVE") {
    
      # Order by increasing difficulty
      X1 = as.matrix(self$data)
  
      # Create new ordered names (a hack to reorder in the end)
      if(length(dimnames(X1)[[2]])) true.names = colnames(X1)
      nc = nchar(as.character(ncol(X1)))
      num.prefix = paste("I",formatC(1:ncol(X1),width=nc,format="d",flag="0"),sep="") # I01, I02 etc. if nc=2
      colnames(X1) = paste(num.prefix,colnames(X1),sep="-")
  
      # Reorder by successes
      k = order(colMeans(X1,na.rm=TRUE),decreasing=TRUE)
      X1 = X1[,k]

      # Remarque : crossprod fonctionne en parallèle sans package spécialisé
      notmissing = crossprod(!is.na(X1),!is.na(X1))
      X0 = 1-X1
      
      # Deal with missing values
      X1[is.na(X1)] = 0
      X0[is.na(X0)] = 0

      # Co-occurrences
      n11 = crossprod(X1,X1) ; n11[n11==0] = 0.5
      n10 = crossprod(X1,X0) ; n10[n10==0] = 0.5
      n00 = crossprod(X0,X0) ; n00[n00==0] = 0.5

      f11 = n11/notmissing
      f10 = n10/notmissing
      f00 = n00/notmissing
  
      D = log(f11)+log(f00)-2*log(f10)

      # Deal with missing distances (iteratively remove variables with large number of missing values)
      keep = 1:ncol(X1)
      while(any( (S <- colSums(is.na(D))) > 0 )) {
        sel = which.max(S)
        D = D[-sel,-sel]
        keep = keep[-sel]
      }
      X1 = X1[,keep]

      diag(D) = 0
      absD = abs(D)  

      # Replace Inf distances by their alternative estimates
      infinite = !is.finite(absD)
      absD[infinite] = absD[t(infinite)]

      # Symmetrize
      upper = upper.tri(absD)
      absD = absD*upper + t(absD*upper)
      
      self$D = absD

      # Implication coefficients
      Iota = tanh(D)

      # MDS
      # [TODO] Voir le package BMDS ou Rdimtools
      self$mdsfit = cmdscale(absD,k=min(self$ndim,ncol(X1)/2),eig=TRUE)
      self$coord = self$mdsfit$points
      self$eig = self$mdsfit$eig
      self$GoF = self$mdsfit$GOF

      # Coordinates (reorder coordinates in the original variable ordering)
      k = order(colnames(X1))
      X1 = X1[,k]
      colnames(X1) = gsub("^I[0-9]*-","",colnames(X1)) # Remove label numbering
      self$coord = self$coord[k,]
      rownames(self$coord) = colnames(X1)
      
			# Cosmetic: Flip coordinates on dim1 by increasing difficulty
			if(cor(self$coord[,1],colMeans(X1,na.rm=TRUE)) > 0) self$coord[,1] = -self$coord[,1]

      # Dimensionality reduction via automatic scree test
      # self$reduce()
      
      # Automatic scale detection
      self$cluster(model=clustModel)
      
      # Implication indicies
      self$Iota = Iota[k,k] ; rownames(self$Iota) = colnames(X1) ; colnames(self$Iota) = colnames(X1)

      # Print first positive eigenvalues
      cat("First eigenvalues:\n")
      e = self$eig[1:10]
      e = e[e>0]
      print(e)
    },
    # Automatic detection of scales by model-based clustering
    cluster = function(model="VVE") {
    
      if(self$ndim == 1) {
        self$scales = rep(1,self$p)
        return()
      }

      # Subscale detection by model based clustering
      self$clustfit = Mclust(self$coord[,1:3],modelNames=model)
      self$scales = self$clustfit$classification
      names(self$scales) = rownames(self$coord)

      cat(length(unique(self$scales)),"scales detected.\n")
    },
    plot = function(axes=NULL,type="confplot",plot.labels=TRUE,plot.grid=TRUE,pch=19,cex=.7,...) {

      if(type=="eig") {
        plot(self$eig,type="h",xlab="Rank",ylab="Eigenvalue",main="Eigenvalues")
        return()
      }

      colors = c("#4285F4","#DB4437","#F4B400","#0F9D58","#9CCBDF","black")
    
      if(is.null(axes)) axes = 1:2
      plot.Mclust(self$clustfit,what="uncertainty",dimens=axes,asp=1,xlim=range(self$coord[,1]),col=colors[self$scales])
      title(main="Scale detection",xlab=paste("Dim.",axes[1]),ylab=paste("Dim.",axes[2]),col="darkgrey")
      points(self$coord[,axes[1]],self$coord[,axes[2]],ylim=range(self$coord[,1]))

      if(plot.grid)   abline(h=seq(-5,5,by=.5),v=seq(-5,5,by=.5),lty=2,col="grey",lwd=1)
      if(plot.labels) {
        shift = diff(range(self$coord[,1])) * 0.015
        text(self$coord[,axes[1]],self$coord[,axes[2]]+shift,rownames(self$coord),col=colors[self$scales],cex=cex,...)
      }

    }

    # # Plot item model
    # plot = function(type="basic",axes=c(1,2),plot.labels=TRUE,pch=19,threshold=0.9,arrows=TRUE,from=NULL,to=NULL,success=TRUE,relsize=0.9,box.size=0.06,box.type="ellipse",curve=.4,shadow.size=0,box.prop=.4,...) {
    
    #   if(type=="basic") {
    #     plot(self$coord[,axes],type="n",asp=1,xlab=paste("Dimension",axes[1]),ylab=paste("Dimension",axes[2]),main="Implication plane",ylim=c(min(self$coord[,axes]),max(self$coord[,axes])))
    #     if(!plot.labels) { points(self$coord[,axes],pch=pch,cex=0.5) }
    #     else             { text(self$coord[,axes],rownames(self$coord),cex=.7) }
    #     abline(h=c(-.25,.25,.5,-.5),lty=2,col="grey",lwd=2)
    #   }
    #   else if(type=="graph") {
 
    #     require(diagram)

    #     items = colnames(self$Iota)
    #     success.rates = round(colMeans(self$data[,items],na.rm=TRUE),2)
    #     graph0 = round(self$Iota,2)
        
    #     # Don't plot below this threshold
    #     graph0[self$Iota < threshold] = -10
        
    #     # Don't plot recursive links
    #     diag(graph0) = -10
    #     if(!arrows) graph0 = matrix(-10,nrow(graph0),ncol(graph0))
                
    #     scaledpos = (self$coord - min(self$coord)) / (max(self$coord)-min(self$coord))
		# 		rownames(scaledpos) = rownames(self$coord)
    #     labels = paste(items,paste("(",success.rates,")",sep=""),sep="\n")
		# 		names(labels) = items
				
		# 		if(is.null(to) && is.null(from)) {
	  #       plotmat(t(graph0),pos=scaledpos[,axes],absent=-10,name=labels,box.size=box.size,box.type=box.type,curve=curve,shadow.size=shadow.size,box.prop=box.prop,relsize=relsize,...)
		# 		}
		# 		else {
		# 			if(length(from) == 1) sel.from = from
		# 			else {
		# 				sel.from = items[which(colSums(as.matrix(graph0[to,]) != -10)>0)]
		# 			}
		# 			if(length(to) == 1) sel.to = to
		# 			else {
		# 				sel.to = items[which(rowSums(as.matrix(graph0[,from]) != -10)>0)]
		# 			}
		# 			subset = unique(c(from,to,sel.to,sel.from))
		# 			graph1 = graph0[subset,subset,drop=FALSE] # Don't drop, to keep matrix format and indexing
		# 			print(graph1)
		# 			if(!is.null(from)) { 
		# 				exclude.from = subset[!(subset %in% from)]
		# 				graph1[,exclude.from] = -10
		# 			}
		# 			if(!is.null(to)) {
		# 				exclude.to = subset[!(subset %in% to)]
		# 				graph1[exclude.to,] = -10
		# 			}
		# 			print(graph1)
		# 			graph0[graph0 != -10] = -10
		# 		  plotmat(t(graph0),pos=scaledpos[,axes],lcol="lightgrey",txt.col="lightgrey",arr.col="lightgrey",absent=-10,name=labels,box.size=box.size,box.type=box.type,curve=curve,shadow.size=shadow.size,box.prop=box.prop,relsize=relsize,...)
		# 			plotmat(t(graph1),pos=scaledpos[subset,axes],add=TRUE,absent=-10,name=labels[subset],box.size=box.size,box.type=box.type,curve=curve,shadow.size=shadow.size,box.prop=box.prop,relsize=relsize,...)
		# 		}
    #   }
    #   else if(type=="cluster") {
    #     opar = par(mfrow=c(1,2))
    #     plot.Mclust(self$clustfit,what="classification",dimens=c(1,2),asp=1)
    #     plot.Mclust(self$clustfit,what="classification",dimens=c(3,2),asp=1)
    #     par(opar)
    #   }
    #   else if(type=="eig") {
    #     plot(self$eig,type="h",xlab="Rank",ylab="Eigenvalue",main="Eigenvalues")
    #   }
    #   else return()
    # }
  )
)