#用于生产arff文件的程序源文件，更改后不需运行rebuilddataprocess.R文件，全整合为一个文件
#读入的arff文件用于提取arff属性设置表头
#可以生成用于CLUS的数据文件,但是仅限于不存在factor类型数据的数据集，20170720
#节点中包含根结点，故验证集、测试集中使用select.node 进行设置，而非except.root.labels
#20170424进行第一次修改，按照rebuilddataprocess.R文件中内容修改、规范变量名称
#20170614进行第二次修改，根据修改后的DatasetSelect函数增加factor.col.index，factor.col.num，factor.levels等参数

#第一次使用，需先运行第三步，得到全局变量
####第一步设置代码存放路径 以及数据存放路径
work.path="H://R//CODES"
data.path="H://R//DATA"
setwd(work.path)#首先将工作路径设置为代码的路径
####第二步 引入库函数以及自定义的函数
source("headfile.R")

####第三步 得到BO结构信息及注释信息
# ontology.sel="BP"
# univ.graph <- Build.universal.graph.ontology.down(ontology = ontology.sel)#得到BP全图
# annotation.final=AnnotationFinal(ontology.sel)#提取BP结构下的所有基因及注释信息
# go.general.list=Get.GO.all.classes(annotation.final)#根据基因的注释GO标签得到基因的全部GO标签
# go.general.table=Build.GO.class.labels(go.general.list)#生成基因及注释信息数据表

####第四步 确定待处理的数据集
#首先选择使用的数据集，1 cellcycle 2 derisi 3 eisen 4 gasch1 5 gasch2 6 church 7 spo 8 seq 9 struc 10 hom
dataset.result=DatasetSelect(dataset.index = 6)
file.prefix=dataset.result[[1]]#得到数据集前缀名称
factor.col.index=dataset.result[[2]]#得到内容为类别信息的列号
factor.col.num=dataset.result[[3]]#得到factor列转化为数值后各列的总数
factor.levels=dataset.result[[4]]
####第五步 读取训练集数据
#此处读取原始数据，替换NA值 进行归一化
select.attributes.en=FALSE
delete.outlier=FALSE
replace.outlier = FALSE
no.del.replace = TRUE
NAreplace=FALSE
Zrescale=FALSE
setwd(data.path)
read.original=TRUE#此处选择读取原数据还是读取经过重采样的csv数据
if(read.original)#选择读取原始的训练文件
{
  train.original=ReadData(paste("originaldata//",file.prefix,"0.train",sep = ""),factor.col.index = factor.col.index,factor.levels = factor.levels)
  train.data=train.original[[1]]#基因的数据信息
  
}else#选择读取重采样后的训练文件
{
  file.name=paste("rebuilddata//",file.prefix,"re_traindataset.csv",sep = "")
  train.data=read.csv(file.name,header = FALSE,row.names = 1)#第一列是基因的名称
}
trainscale.result=TraindataScale(train.data,factor.col.num,delete.outlier=delete.outlier,replace.outlier=replace.outlier,
                                 no.del.replace =no.del.replace,NAreplace=NAreplace,Zrescale=Zrescale)
remain.data=trainscale.result[[1]]
sp=trainscale.result[[2]]

####第六步 得到目前有注释信息的基因名称，并将其转换为list及table形式

#得到训练集共有的基因的名称列表
common.genes <- Get.all.common.genes(go.general.table, remain.data)
remain.select.data=remain.data[common.genes,]#此句不可删去，原始数据时仍需使用
#得到common genes中每个基因对应的全部GO标签列表
match.go.general=go.general.list[common.genes]

#将标签列表转换为TABLE形式，行名称为基因名称，列名称为GO标签
match.go.table=Build.GO.class.labels(match.go.general)

####第七步 对基因注释信息进行处理，根据GO标签注释的样本数选择适当的GO标签

#得到目前的所有基因共有多少个不重复的GO标签
all.go.labels=Get.classes(match.go.general)

#得到DAG图中包含select.num=100个样本以上的节点列表及对应的基因名称
go.label.list=DataCleaning(match.go.general,match.go.table,select.num = 100)
all.go.labels.sel=Get.classes(go.label.list)
#except.root.labels=setdiff(all.go.labels.50,"GO:0008150")#除去根结点后所剩节点

####第八步 根据选择的GO节点生成DAG图，并得到DAG图中的层级信息

#得到修改后的DAG图中各节点的层级
graph.general.sel <- subGraph(all.go.labels.sel, univ.graph)
graph.level.sel=GraphLevel(graph.general.sel)
go.level.statistics=LevelStatistics(graph.level.sel)
go.for.each.level=go.level.statistics[[1]]
total.levels=length(go.for.each.level)


####第九步 按照层级选择一定数量的节点，生成这些节点的子图

#选择根结点至total.levels的所有go节点，select.node为向量形式
select.node=NodeSelectByLevel(go.level.statistics,total.levels,add.root.node = TRUE)
except.root.labels=setdiff(select.node,"GO:0008150")#去掉所选节点中的根结点，剩余节点的集合为需求的标签
sub.graph <- subGraph(select.node, univ.graph)#绘制子图


####第十步 读取验证集文件

if(read.original)#选择读取原始的验证集文件
{
  #读入valid基因特征属性
  valid.original=ReadData(paste("originaldata//",file.prefix,"0.valid",sep = ""),factor.col.index = factor.col.index,factor.levels = factor.levels)
  valid.data=valid.original[[1]]#基因的数据信息
  write.data.fname=paste(file.prefix,"0_validdataset.csv",sep = "")
  write.class.fname=paste(file.prefix,"0_validclass.csv",sep = "")
  
}else#选择读取重采样后的验证集文件
{
  file.name=paste("rebuilddata//",file.prefix,"re_validdataset.csv",sep = "")
  valid.data=read.csv(file.name,header = FALSE,row.names = 1)
  write.data.fname=paste(file.prefix,"1_validdataset.csv",sep = "")
  write.class.fname=paste(file.prefix,"1_validclass.csv",sep = "")
}




valid.scaled.data=ValiddataScale(valid.data,factor.col.num,sp,replace.outlier=replace.outlier,
                                 no.del.replace = no.del.replace,NAreplace=NAreplace,Zrescale=Zrescale)
#此处except.root.labels变量设置为select.node，包含根结点
valid.data.total=BuildValidset(valid.scaled.data,go.general.table,go.general.list,select.node,
                               write.data.enable=TRUE,write.class.enable=TRUE,write.data.fname=write.data.fname,
                               write.class.fname=write.class.fname,select.attributes.en=select.attributes.en,select.attributes)
valid.select.data=valid.data.total[[1]]
valid.select.table=valid.data.total[[2]]
valid.go.label.list=valid.data.total[[3]]
#直接将所有的标签传给每个样本作为标签，而不是像原文件一样只选择最具体的标签
valid.leaf.label.list=valid.go.label.list

####第十一步 读取测试集文件
setwd(data.path)
if(read.original)#选择读取原始的测试集文件
{
  #读入test基因特征属性
  test.original=ReadData(paste("originaldata//",file.prefix,"0.propertest",sep = ""),factor.col.index = factor.col.index,factor.levels = factor.levels)
  test.data=test.original[[1]]#基因的数据信息
  write.data.fname=paste(file.prefix,"0_testdataset.csv",sep = "")
  write.class.fname=paste(file.prefix,"0_testclass.csv",sep = "")
  
}else#选择读取重采样后的测试集文件
{
  file.name=paste("rebuilddata//",file.prefix,"re_testdataset.csv",sep = "")
  test.data=read.csv(file.name,header = FALSE,row.names = 1)
  write.data.fname=paste(file.prefix,"1_testdataset.csv",sep = "")
  write.class.fname=paste(file.prefix,"1_testclass.csv",sep = "")
}

test.scaled.data=ValiddataScale(test.data,factor.col.num,sp,replace.outlier=replace.outlier,
                                no.del.replace = no.del.replace,NAreplace=NAreplace,Zrescale=Zrescale)
test.data.total=BuildValidset(test.scaled.data,go.general.table,go.general.list,select.node,
                              write.data.enable=TRUE,write.class.enable=TRUE,write.data.fname=write.data.fname,
                              write.class.fname=write.class.fname,select.attributes.en=select.attributes.en,select.attributes)
test.select.data=test.data.total[[1]]
test.select.table=test.data.total[[2]]
test.go.label.list=test.data.total[[3]]
test.leaf.label.list=test.go.label.list


####第十二步 生成新的arff文件
setwd(paste(data.path,"//arff",sep = ""))

arff.result=WriteArff(original.filename=paste(file.prefix,"_GO.train.arff",sep = ""),write.filename= paste(file.prefix,"train.arff",sep = ""),
                      modified.data=remain.select.data,common.genes=common.genes,ontology=ontology.sel,
                      class.label=go.label.list,class.graph=sub.graph,genename.exist=FALSE) 

arff.result=WriteArff(original.filename=paste(file.prefix,"_GO.valid.arff",sep = ""),write.filename= paste(file.prefix,"valid.arff",sep = ""),
                      modified.data=valid.select.data,common.genes=NULL,ontology=ontology.sel,
                      class.label=valid.leaf.label.list,class.graph=sub.graph,genename.exist=FALSE)

arff.result=WriteArff(original.filename=paste(file.prefix,"_GO.test.arff",sep = ""),write.filename= paste(file.prefix,"test.arff",sep = ""),
                      modified.data=test.select.data,common.genes=NULL,ontology=ontology.sel,
                      class.label=test.leaf.label.list,class.graph=sub.graph,genename.exist=FALSE) 
