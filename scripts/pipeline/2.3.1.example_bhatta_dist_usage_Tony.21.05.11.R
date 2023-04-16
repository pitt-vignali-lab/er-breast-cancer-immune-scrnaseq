## Example usage of Bhattacharrya Distance 
## May 11 2021

library(Seurat)

## This example assumes one is comparing the distance between PBMC and TIL for certain cell types, using PC embeddings to measure the distance. 

# create list of pca embeddings for Bhattacharrya
pbmc.split <- SplitObject(pbmc,split.by="cell_types")

pbmc.embed <- lapply(covid.pbmc.split,function(x) Embeddings(x,reduction="pca"))

til.split <- SplitObject(til,split.by="new_cell_types")

til.embed <- lapply(til.split,function(x) Embeddings(x,reduction="pca"))


# Bhattcharrya distance function
# Example: use 50 PCs and sample 100 cells
Db <- function(pca.mat.x,pca.mat.y,pcs.use=1:50,num.sample=100,random.sample=FALSE) {

	if (random.sample==FALSE) {

		samp1 <- sample(1:nrow(pca.mat.x),num.sample)
		samp2 <- sample(1:nrow(pca.mat.y),num.sample)

		mat.x <- pca.mat.x[samp1,pcs.use]
		mat.y <- pca.mat.y[samp2,pcs.use]

		mean.mat1 <- apply(mat.x,2,mean)
		mean.mat2 <- apply(mat.y,2,mean)

		cov.mat1 <- cov(mat.x)
		cov.mat2 <- cov(mat.y)

		comb.cov <- (cov.mat1 + cov.mat2)/2

		(1/8)*t(mean.mat1-mean.mat2)%*%solve(comb.cov)%*%(mean.mat1-mean.mat2)+0.5*log(det(comb.cov)/sqrt(det(cov.mat1)*det(cov.mat2)))

	} else {

		comb.mat <- rbind(pca.mat.x,pca.mat.y)

		samp1 <- sample(1:nrow(comb.mat),num.sample)
		samp2 <- sample(1:nrow(comb.mat),num.sample)

		mat.x <- comb.mat[samp1,pcs.use]
		mat.y <- comb.mat[samp2,pcs.use]

		mean.mat1 <- apply(mat.x,2,mean)
		mean.mat2 <- apply(mat.y,2,mean)

		cov.mat1 <- cov(mat.x[,pcs.use])
		cov.mat2 <- cov(mat.y[,pcs.use])

		comb.cov <- (cov.mat1 + cov.mat2)/2

		(1/8)*t(mean.mat1-mean.mat2)%*%solve(comb.cov)%*%(mean.mat1-mean.mat2)+0.5*log(det(comb.cov)/sqrt(det(cov.mat1)*det(cov.mat2)))


	}

}

# Set up vectors for results

Db.samp <- Db.rand <- vector("logical",length=100)
Db.list <- vector("list",length=5)

# Subset datasets to only cell types we want to measure Bhattacharrya distance for 

pbmc.embed.test <- pbmc.embed[c("B cells","CD4+ T cells","CD8+ T cells","Classical monocytes","NK cells")]

til.embed.test <- til.embed[c("B cells","CD4+ T cells","CD8+ T cells","Classical monocytes","NK cells")]

names(Db.list) <- names(til.embed.test)

# Measure BD by looping through 5 cell types 100 times
# Random sample creates random distance

for (a in 1:5) {

for (i in 1:100) {

	Db.samp[i] <- Db(pca.mat.x=pbmc.embed.test[[a]],pca.mat.y=til.embed.test[[a]],pcs.use=1:50,num.sample=100,random.sample=FALSE)

	Db.rand[i] <- Db(pca.mat.x=pbmc.embed.test[[a]],pca.mat.y=til.embed.test[[a]],pcs.use=1:50,num.sample=100,random.sample=TRUE)

}

Db.list[[a]] <- list(Db.samp,Db.rand)

}

# Median distances
sapply(Db.list,function(x) sapply(x,function(y) median(y)))

# Create data.frame from distances for each cell type

sampled.distances <- sapply(Db.list,function(x) sapply(x[[1]],function(y) cbind(y)))
random.distances <- sapply(Db.list,function(x) sapply(x[[2]],function(y) cbind(y)))

colnames(sampled.distances) <- paste(colnames(sampled.distances),"_","sampled",sep="")

colnames(random.distances) <- paste(colnames(random.distances),"_","random",sep="")

Db.frame <- cbind(sampled.distances,random.distances)

Db.frame.melt <- reshape2::melt(Db.frame)
