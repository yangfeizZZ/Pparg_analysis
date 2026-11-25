g <- pheatmap(all_normalize_DEG,cluster_rows = TRUE,show_rownames = FALSE,scale = "row",cluster_cols = FALSE,use_raster =FALSE)
row_cluster <- cutree(g$tree_row, k = 4)
cluster_genes <- split(names(row_cluster), row_cluster)
gene_order <- rownames(all_normalize_DEG)[g$tree_row$order]
all_normalize_DEG <- all_normalize_DEG[gene_order,]

#row_anno <- data.frame(Cluster = factor(row_cluster))
#rownames(row_anno) <- names(row_cluster)
#row_ha <- rowAnnotation(Cluster = row_cluster,col = list(Cluster = c("1" = "#E41A1C","2" = "#377EB8","3" = "#4DAF4A"),"4"="black"))
row_split <- row_cluster[gene_order]
Heatmap(t(scale(t(all_normalize_DEG))),cluster_rows = FALSE,cluster_columns = FALSE,
        show_row_names = FALSE,show_column_names = TRUE,col = colorRamp2(c(-2,0,2),c("blue","grey","red3")),
        use_raster = FALSE,row_split = row_split)


pheatmap(all_normalize_DEG,cluster_rows = FALSE,show_rownames = FALSE,
  scale = "row",cluster_cols = FALSE,use_raster = FALSE)

#-----------------------------------

#---------------diff bind
library(DiffBind)
D4KO_D4WT <- dba(sampleSheet = "E:/scRNA-seq/yzb/ATAC/D4_downsample/diff_bind/diff_bind.csv")
D4KO_D4WT <- dba.count(D4KO_D4WT,bRemoveDuplicates = TRUE)
D4KO_D4WT <- dba.contrast(D4KO_D4WT,minMembers = 2,categories = DBA_CONDITION)
D4KO_D4WT <- dba.analyze(D4KO_D4WT,method = DBA_DESEQ2)
D4KO_D4WT <- dba.report(D4KO_D4WT,contrast = 1)
D4KO_D4WT <- as.data.frame(D4KO_D4WT)
D4KO_D4WT$chang < "None"
D4KO_D4WT_gain <- D4KO_D4WT[D4KO_D4WT$change == "Up", c("seqnames","start","end")]
D4KO_D4WT_loss <- D4KO_D4WT[D4KO_D4WT$change == "Down", c("seqnames","start","end")]
write.table(D4KO_D4WT_gain,file = "E:/scRNA-seq/yzb/ATAC/D4_downsample/diff_bind/D4KO_D4WT_gain_peaks.bed",sep = "\t",quote = FALSE, row.names = FALSE,col.names = FALSE)
write.table(D4KO_D4WT_loss,file = "E:/scRNA-seq/yzb/ATAC/D4_downsample/diff_bind/D4KO_D4WT_loss_peaks.bed",sep = "\t",quote = FALSE,row.names = FALSE,col.names = FALSE)
rm(D4KO_D4WT_gain,D4KO_D4WT_loss)

D14KO_D14WT <- dba(sampleSheet = "E:/scRNA-seq/yzb/ATAC/D14_downsample/diff_bind/diff_bind.csv")
D14KO_D14WT <- dba.count(D14KO_D14WT,bRemoveDuplicates = TRUE)
D14KO_D14WT <- dba.contrast(D14KO_D14WT,minMembers = 2,categories = DBA_CONDITION)
D14KO_D14WT <- dba.analyze(D14KO_D14WT,method = DBA_DESEQ2)
D14KO_D14WT <- dba.report(D14KO_D14WT,contrast = 1)
D14KO_D14WT <- as.data.frame(D14KO_D14WT)
D14KO_D14WT$change < "None"
D14KO_D14WT_gain <- D14KO_D14WT[D14KO_D14WT$change == "Up", c("seqnames","start","end")]
D14KO_D14WT_loss <- D14KO_D14WT[D14KO_D14WT$change == "Down", c("seqnames","start","end")]
write.table(D14KO_D14WT_gain,file = "E:/scRNA-seq/yzb/ATAC/D14_downsample/diff_bind/D14KO_D14WT_gain_peaks.bed",sep = "\t",quote = FALSE, row.names = FALSE,col.names = FALSE)
write.table(D14KO_D14WT_loss,file = "E:/scRNA-seq/yzb/ATAC/D14_downsample/diff_bind/D14KO_D14WT_loss_peaks.bed",sep = "\t",quote = FALSE,row.names = FALSE,col.names = FALSE)
rm(D14KO_D14WT_gain,D14KO_D14WT_loss)


#------------annotation of chipseeker--------
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
D4_gain_chipseeker <- annotatePeak(D4_gain,tssRegion = c(-3000, 3000),  #启动子区域,
                         flankDistance = 5000,TxDb = txdb,annoDb = "org.Ss.eg.db")
D4_gain_chipseeker_df <- as.data.frame(D4_gain_chipseeker)
D4_gain_chipseeker_df_go <- enrichGO(D4_gain_chipseeker_df$symbol[grep("Promoter",D4_gain_chipseeker_df$annotation)],
                                     keyType = "SYMBOL",OrgDb = org.Mm.eg.db,ont = "BP",readable = FALSE) %>% as.data.frame()

D4_loss_chipseeker <- annotatePeak(D4_loss,tssRegion = c(-3000, 3000),  #启动子区域
                                   flankDistance = 5000,TxDb = txdb,annoDb = "org.Ss.eg.db")
D4_loss_chipseeker_df <- as.data.frame(D4_loss_chipseeker)
D4_loss_chipseeker_df_go <- enrichGO(D4_loss_chipseeker_df$symbol,
                                     keyType = "SYMBOL",OrgDb = org.Mm.eg.db,ont = "BP",readable = FALSE) %>% as.data.frame()


D14_gain_chipseeker <- annotatePeak(D14_gain,tssRegion = c(-3000, 3000),  #启动子区域,
                                   flankDistance = 5000,TxDb = txdb,annoDb = "org.Ss.eg.db")
D14_gain_chipseeker_df <- as.data.frame(D14_gain_chipseeker)
D14_gain_chipseeker_df_go <- enrichGO(D14_gain_chipseeker_df$symbol[grep("Promoter",D14_gain_chipseeker_df$annotation)],
                                     keyType = "SYMBOL",OrgDb = org.Mm.eg.db,ont = "BP",readable = FALSE) %>% as.data.frame()

D14_loss_chipseeker <- annotatePeak(D14_loss,tssRegion = c(-3000, 3000),  #启动子区域
                                   flankDistance = 5000,TxDb = txdb,annoDb = "org.Ss.eg.db")
D14_loss_chipseeker_df <- as.data.frame(D14_loss_chipseeker)
D14_loss_chipseeker_df_go <- enrichGO(D14_loss_chipseeker_df$symbol,
                                     keyType = "SYMBOL",OrgDb = org.Mm.eg.db,ont = "BP",readable = FALSE) %>% as.data.frame()


