#用于生成MIML中所需的训练及测试文件
work.path="H://R//CODES"
data.path="H://R//DATA"
setwd(work.path)#设置工作路径

#引入库函数以及自定义的函数
source("headfile.R")
BP.univ.graph <- Build.universal.graph.ontology.down(ontology = "BP")#得到BP全图
annotation.final.BP=AnnotationFinal("BP")#提取BP结构下的所有基因及注释信息
go.general.list.BP=Get.GO.all.classes(annotation.final.BP)#根据基因的注释GO标签得到基因的全部GO标签
go.general.table.BP=Build.GO.class.labels(go.general.list.BP)#生成基因及注释信息数据表

read.original=FALSE
file.prefix="derisi"
factor.col=c(0)
work.path="H://R//DATA"
file.savepath="H://R//DATA//matlabfile"
delete.outlier=FALSE
replace.outlier = FALSE
no.del.replace = FALSE
NAreplace=TRUE
Zrescale=TRUE
write.data.enable=TRUE
write.class.enable=TRUE
BuildMatlabDataset(read.original=read.original,file.prefix=file.prefix, factor.col=factor.col, data.path=data.path,file.savepath=file.savepath,
                   delete.outlier=delete.outlier,replace.outlier=replace.outlier, no.del.replace =no.del.replace,NAreplace=NAreplace,
                   Zrescale=Zrescale, write.data.enable=write.data.enable,write.class.enable=write.class.enable)

#read.original=TRUE
file.prefix="eisen"
factor.col=c(0)
#file.savepath="H://R//DATA//matlabfile"
# delete.outlier=FALSE
# replace.outlier = FALSE
# no.del.replace = FALSE
# NAreplace=TRUE
# Zrescale=TRUE
# write.data.enable=TRUE
# write.class.enable=TRUE
BuildMatlabDataset(read.original=read.original,file.prefix=file.prefix, factor.col=factor.col, data.path=data.path,file.savepath=file.savepath,
                   delete.outlier=delete.outlier,replace.outlier=replace.outlier,no.del.replace =no.del.replace, NAreplace=NAreplace,
                   Zrescale=Zrescale, write.data.enable=write.data.enable,write.class.enable=write.class.enable)
read.original=TRUE
file.prefix="cellcycle"
factor.col=c(0)
# file.savepath="H://R//DATA//matlabfile"
# delete.outlier=FALSE
# replace.outlier = FALSE
# no.del.replace = TRUE
# NAreplace=TRUE
# Zrescale=TRUE
# write.data.enable=TRUE
# write.class.enable=TRUE
BuildMatlabDataset(read.original=read.original,file.prefix=file.prefix, factor.col=factor.col, data.path=data.path,file.savepath=file.savepath,
                   delete.outlier=delete.outlier,replace.outlier=replace.outlier,no.del.replace =no.del.replace, NAreplace=NAreplace,
                   Zrescale=Zrescale, write.data.enable=write.data.enable,write.class.enable=write.class.enable)
#read.original=TRUE
file.prefix="gasch2"
factor.col=c(0)
# file.savepath="H://R//DATA//matlabfile"
# delete.outlier=FALSE
# replace.outlier = FALSE
# no.del.replace = TRUE
# NAreplace=TRUE
# Zrescale=TRUE
# write.data.enable=TRUE
# write.class.enable=TRUE
BuildMatlabDataset(read.original=read.original,file.prefix=file.prefix, factor.col=factor.col,data.path=data.path,file.savepath=file.savepath,
                   delete.outlier=delete.outlier,replace.outlier=replace.outlier,no.del.replace =no.del.replace, NAreplace=NAreplace,
                   Zrescale=Zrescale, write.data.enable=write.data.enable,write.class.enable=write.class.enable)
