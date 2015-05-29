library(ggplot2)
library(dplyr)

data <- read.csv("leaf_389_sample.csv", header=T, row.names=1)
data <- data + 1
data <- log2(data)

leaf <- tbl_df(data)
##============================================================##
replicate1 <- read.table("leaf_replicate1.txt", header=T)
rep1 <- as.character(replicate1$samp1)

rep1_df <- as.data.frame(0)
gene_data1 <- as.data.frame(c())
for (j in 1:length(rep1))
{
    
    gene_data1 <- as.data.frame(select(data, starts_with(rep1[j])) ) 
    
    rep1_df <- as.data.frame(cbind(rep1_df, gene_data1))
}
rep1_df <- rep1_df[-1]
rm(gene_data1)
##============================================================##

replicate2 <- read.table("leaf_replicate2.txt", header=T)
rep2 <- as.character(replicate2$samp2)

rep2_df <- as.data.frame(0)
gene_data2 <- as.data.frame(c())
for (i in 1:length(rep2))
{
    
    gene_data2 <- as.data.frame(select(data, starts_with(rep2[i])) ) 
    
    rep2_df <- as.data.frame(cbind(rep2_df, gene_data2))
}
rep2_df <- rep2_df[-1]
rm(gene_data2)

##==============================================================##
t_rep1 <- t(rep1_df)
t_rep2 <- t(rep2_df)
gene <- row.names(data)



corr_value <- matrix(nrow=nrow(data), ncol=2)
rownames(corr_value) <- row.names(data)

colnames(corr_value) <- c("cor_value", "extra")

t_rep1_df <- as.data.frame(0)
t_rep2_df <- as.data.frame(0)



for (k in 1:nrow(data))
{
    
    t_rep1_df <- as.data.frame(t_rep1[, gene[k]])
    t_rep2_df <- as.data.frame(t_rep2[, gene[k]])
    
    t_gene_df <- as.data.frame(cbind(t_rep1_df, t_rep2_df))
    
    t_gene_df$Avg <- apply(t_gene_df, 1, mean, na.rm =T)
    names(t_gene_df) <- c("trep1", "trep2")
    
    rep1_val <- t_gene_df$trep1
    rep2_val <- t_gene_df$trep2
    
    cor_value=cor(rep1_val,rep2_val)
    
    corr_value[k,1] <- cor_value 
    
    
    
    q <- ggplot(t_gene_df, aes(rep1_val, rep2_val)) + geom_point(shape=1)
    q
    q <- ggplot(t_gene_df, aes(rep1_val, rep2_val)) + geom_point(color="red", shape=1) + geom_smooth(method = "lm", se=TRUE)
    
    #q <- q+ p + theme(axis.title.y = element_text(size=17, colour = rgb(0,0,0)))
    #Make the y axis title black and size 17:
    q = q + theme(axis.title.y = element_text(size=17, colour = rgb(0,0,0)))
    
    #Make the x axis title black and size 17:
    q = q + theme(axis.title.x = element_text(size=17, colour=rgb(0,0,0)));
    
    q =q + labs(x="Replicate 1", y="Replicate 2" )
    
    #Add a title:
    q= q + ggtitle(quote("FPKM values per Gene in Biological Replicate"))
    q
    
    mypath <- file.path("leaf_gene_plot",paste(gene[k],"_plot.png", sep=''))
    
    ggsave(file=mypath, width=6, height=4, dpi=72)
    
    
}
write.csv(corr_value, file="leaf_gene_correlation_coefficent_value.csv", row.names=T)   

comm_leaf <- read.csv("Leaf_filt_2_fold_change_stat_Oct_10.csv", header=T, row.names=1)
leaf_stat <- comm_leaf
comm_leaf_FPKM <- read.csv("common_set_leaf_FPKM.csv", header=T, row.names=1)
comm_leaf_id <- row.names(comm_leaf_FPKM)
comm_leaf_stat <- leaf_stat[comm_leaf_id, ] 
comm_leaf_xorr <- corr_value[comm_leaf_id, ]
comm_leaf_stat_n_corr <- cbind(comm_leaf_stat,comm_leaf_xorr)
comm_leaf_stat_n_corr <- comm_leaf_stat_n_corr[,-10]
write.csv(comm_leaf_stat_n_corr,file="Leaf_common_set_stat_n_replicate_correlation.csv", row.names=T)

exclusive_new <- read.csv("exclusive_new_leaf_FPKM.csv", header=T, row.names=1)
exclusive_leaf_id <- row.names(exclusive_new)
exclusive_leaf_stat <- leaf_stat[exclusive_leaf_id, ]
exclusive_leaf_corr <- corr_value[exclusive_leaf_id, ]
exclusive_leaf_stat_n_corr <- cbind(exclusive_leaf_stat, exclusive_leaf_corr)
write.csv(exclusive_leaf_stat_n_corr,file="Leaf_exclusive_set_stat_n_replicate_correlation.csv", row.names=T)
