x <- org.Hs.egGO
gene_names <- mappedkeys(x)
xx <- as.list(x[gene_names])
#生成与xx列表等长的annotation_list列表用于存放以后提取的数据，并对此list用基因名称命名
annotation_list=list(NULL)
length(annotation_list)<-length(xx)
names(annotation_list)<-gene_names
gene_names[18600]
