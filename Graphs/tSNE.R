library(Rtsne)
train<-
Labels<-train$Gene.ID
train$label<-as.factor(train$Gene.ID)
colors = rainbow(length(unique(train$Gene.ID)))
names(colors) = unique(train$Gene.ID)
tsne <- Rtsne(train[,-1], dims = 2, perplexity=17, verbose=TRUE, max_iter = 500)
exeTimeTsne<- system.time(Rtsne(train[,-1], dims = 2, perplexity=17, verbose=TRUE, max_iter = 500))
## Plotting
plot(tsne$Y, t='n', main="tsne")
text(tsne$Y, labels=train$Gene.ID, col=colors[train$Gene.ID])
