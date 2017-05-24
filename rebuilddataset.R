#用于重新生成训练集 验证集 测试集 20170329
#加入了数据集选择函数，修改factor.col为factor.col.index
work.path="H://R//CODES"
data.path="H://R//DATA"
setwd(work.path)#设置工作路径

#引入库函数以及自定义的函数
source("headfile.R")
BP.univ.graph <- Build.universal.graph.ontology.down(ontology = "BP")#得到BP全图
annotation.final.BP=AnnotationFinal("BP")#提取BP结构下的所有基因及注释信息
go.general.list.BP=Get.GO.all.classes(annotation.final.BP)#根据基因的注释GO标签得到基因的全部GO标签
go.general.table.BP=Build.GO.class.labels(go.general.list.BP)#生成基因及注释信息数据表

file.savepath="H://R//DATA//rebuilddata"
write.data.enable=TRUE
row.names.enable=TRUE



#首先选择使用的数据集，1 cellcycle 2 derisi 3 eisen 4 gasch1 5 gasch2 6 church 7 spo 8 seq 9 struc 10 hom

dataset.result=DatasetSelect(dataset.index = 1)
file.prefix=dataset.result[[1]]#得到数据集前缀名称
factor.col.index=dataset.result[[2]]#得到内容为类别信息的列号
factor.col.num=dataset.result[[3]]#得到factor列转化为数值后各列的总数
factor.levels=dataset.result[[4]]

DatasetRebuild(file.prefix=file.prefix,factor.col.index = factor.col.index,factor.levels = factor.levels,data.path,file.savepath, write.data.enable,row.names.enable)

# #生成eisen数据集
# file.prefix="eisen"
# factor.col=c(0)
# DatasetRebuild(file.prefix=file.prefix,factor.col.index = factor.col.index,factor.levels = factor.levels,data.path,file.savepath, write.data.enable,row.names.enable)
# #生成cellcycle数据集
# file.prefix="cellcycle"
# factor.col=c(0)
# DatasetRebuild(file.prefix=file.prefix,factor.col.index = factor.col.index,factor.levels = factor.levels,data.path,file.savepath, write.data.enable,row.names.enable)
# #生成gasch2数据集
# file.prefix="gasch2"
# factor.col=c(0)
# DatasetRebuild(file.prefix=file.prefix,factor.col.index = factor.col.index,factor.levels = factor.levels,data.path,file.savepath, write.data.enable,row.names.enable)
# 
