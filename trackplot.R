source("~/YZB/trackplot.R")
Sys.setenv(PATH=paste("/public/home/fyang/YZB/bwtool/bin", Sys.getenv("PATH"), sep=":"))
Sys.setenv(PATH=paste("/public/home/fyang/miniconda3/bin/", Sys.getenv("PATH"), sep=":"))
ATAC_bigWigs = c("your bigwig file")
RNA_bigWigs = c("your bigwig file")
cuttag_bigWigs = c("your bigwig file")

ATAC_bigWigs = read_coldata(bws = ATAC_bigWigs, build = "mm10")
RNA_bigWigs = read_coldata(bws = RNA_bigWigs, build = "mm10")
cuttag_bigWigs = read_coldata(bws = cuttag_bigWigs, build = "mm10")

track_plot(summary_list = track_extract(colData = ATAC_bigWigs, 
                                        loci = "",ideoTblName = "cytoBandIdeo"),
           col = c("#d35400","#2980b9"),y_min = 0,y_max =3