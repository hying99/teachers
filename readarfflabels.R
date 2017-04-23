#用于将CLUS生成的test.pred.arff文件读入并计算其结果的评价指标
#20170421
setwd(data.path)
#输入待处理的数据集 名称
file.prefix="cellcycle"
#设置arff文件的存储文件夹名
fold.name=""

labels.result=ReadPredictArff(data.path,fold.name,file.prefix,part.num=3)
true.labels=labels.result[[1]]
original.plabels=labels.result[[2]]
pruned.plabels=labels.result[[3]]
