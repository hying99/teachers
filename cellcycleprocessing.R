#用于cellcycle数据的处理，训练前期数据处理，用于生成合适的训练集以及测试集
work.path="H://R//CODES"
data.path="H://R//DATA"
setwd(work.path)#设置工作路径


#引入库函数以及自定义的函数
source("headfile.R")
BP.univ.graph <- Build.universal.graph.ontology.down(ontology = "BP")#得到BP全图
annotation.final.BP=AnnotationFinal("BP")#提取BP结构下的所有基因及注释信息
go.general.list.BP=Get.GO.all.classes(annotation.final.BP)#根据基因的注释GO标签得到基因的全部GO标签
go.general.table.BP=Build.GO.class.labels(go.general.list.BP)#生成基因及注释信息数据表

file.prefix="cellcycle"
#matrix.cellcycle=ReadData("originaldata//cellcycle0.train")#读入基因特征属性
#读取数据时，01调整Zrescale和标准化ZNormal不能同时使用
factor.col=c(0)
NAreplace = FALSE
Zrescale = FALSE
ZNormal = FALSE
setwd(data.path)
matrix.cellcycle=ReadData(paste("originaldata//",file.prefix,"0.train",sep = ""),factor.col = factor.col)
matrix.cellcycle.data=matrix.cellcycle[[1]]#基因的数据信息

training.cellcycle.data=matrix.cellcycle.data
#training.cellcycle.name=rownames(training.cellcycle.data)#基因名称列表
#对输入数据剔除异常值，并进行归一化处理
trainscale.result=TraindataScale(training.cellcycle.data,factor.col,delete.outlier=FALSE,replace.outlier = FALSE,NAreplace=TRUE,Zrescale=TRUE)
remain.data=trainscale.result[[1]]
sp=trainscale.result[[2]]

#得到几个数据集共有的基因的名称列表
common.genes <- Get.all.common.genes(go.general.table.BP, remain.data)
#得到common genes中每个基因对应的全部GO标签列表
match.go.general=go.general.list.BP[common.genes]

#将标签列表转换为TABLE形式，行名称为基因名称，列名称为GO标签
match.go.table=Build.GO.class.labels(match.go.general)


#得到目前的所有基因共有多少个不重复的GO标签
all.go.labels=Get.classes(match.go.general)



#得到DAG图中包含50个样本以上的节点列表及对应的基因名称
go.label.list=DataCleaning(match.go.general,match.go.table,select.num = 100)
all.go.labels.50=Get.classes(go.label.list)
#except.root.labels=setdiff(all.go.labels.50,"GO:0008150")#除去根结点后所剩节点


#得到DAG选择节点后各基因所包含的底层节点
#go.each.gene=Setfinalannotation(univ.graph = BP.univ.graph,go.label.list = go.label.list)

#得到修改后的DAG图中各节点的层级
graph.BP.general.50 <- subGraph(all.go.labels.50, BP.univ.graph)
graph.BP.level.50=GraphLevel(graph.BP.general.50)
go.level.statistics=LevelStatistics(graph.BP.level.50)
go.for.each.level=go.level.statistics[[1]]
each.level.nodes.num=go.level.statistics[[2]]
total.levels=length(go.for.each.level)


select.node.3=NodeSelectByLevel(go.level.statistics,total.levels,add.root.node = TRUE)#选择根结点至total.levels的所有go节点
except.root.labels.3=setdiff(select.node.3,"GO:0008150")
graph.select.node.3 <- subGraph(select.node.3, BP.univ.graph)
PlotLabelGraph(except.root.labels.3,BP.univ.graph,num.only=TRUE,plot.en = TRUE,output.en = FALSE,write.pic.name = "go.graph.level.ps")

each.go.level.num=graph.BP.level.50[except.root.labels.3]#表示每个节点所属于的层级
each.go.weight=unname(each.go.level.num)

for (i in 1:length(each.go.level.num))
{
  each.go.weight[i]=(total.levels+1-each.go.weight[i])/(total.levels+1)
}

#go.leaf.nodes.3=GetLeafNode1(graph.select.node.3)

for (i in 1:length(go.label.list))
{
  go.label.list[[i]]=intersect(go.label.list[[i]],select.node.3)
}

root.table.3=Build.GO.class.labels(go.label.list)
setwd("H://R//DATA//traindata")

select.attributes.en=FALSE

data.total=BuildTrainDataset(root.table.3, except.root.labels.3, data.matrix=remain.data,
                             ontology = "BP", adjust.ratio=0.2,ratio.negative = 0, common.genes = common.genes,
                             seed = 1,select.attributes.en=select.attributes.en,write.en=TRUE)



setwd("H://R//DATA")
original.valid.file=paste("originaldata//",file.prefix,"0.valid",sep = "")
select.attributes.en=FALSE


if(select.attributes.en==FALSE)
{
  write.data.fname=paste(file.prefix,"0_validdataset.csv",sep = "")
  write.class.fname=paste(file.prefix,"0_validclass.csv",sep = "")
  valid.cellcycle=ReadData(original.valid.file,factor.col = factor.col)#读入valid基因特征属性
  valid.cellcycle.data=valid.cellcycle[[1]]#valid基因的数据信息
  valid.scaled.data=ValiddataScale(valid.cellcycle.data,factor.col,sp,replace.outlier=FALSE,NAreplace=TRUE,Zrescale=TRUE)
  
  valid.data.total=BuildValidset(valid.scaled.data,go.general.table.BP,go.general.list.BP,except.root.labels.3,
                                 write.data.enable=TRUE,write.class.enable=TRUE,write.data.fname=write.data.fname,
                                 write.class.fname=write.class.fname,select.attributes.en=select.attributes.en,select.attributes)
  valid.select.data=valid.data.total[[1]]
  valid.select.table=valid.data.total[[2]]
}

#original.test.file="originaldata//cellcycle0.propertest"
original.test.file=paste("originaldata//",file.prefix,"0.propertest",sep = "")
select.attributes.en=FALSE
if(select.attributes.en==FALSE)
{
  
  write.data.fname=paste(file.prefix,"0_testdataset.csv",sep = "")
  write.class.fname=paste(file.prefix,"0_testclass.csv",sep = "")
  test.cellcycle=ReadData(original.test.file,factor.col = factor.col)#读入test基因特征属性
  test.cellcycle.data=test.cellcycle[[1]]#test基因的数据信息
  test.scaled.data=ValiddataScale(test.cellcycle.data,factor.col,sp,replace.outlier=FALSE,NAreplace=TRUE,Zrescale=TRUE)
  test.data.total=BuildValidset(test.scaled.data,go.general.table.BP,go.general.list.BP,except.root.labels.3,
                                write.data.enable=TRUE,write.class.enable=TRUE,write.data.fname=write.data.fname,
                                write.class.fname=write.class.fname,select.attributes.en=select.attributes.en,select.attributes)
  test.select.data=test.data.total[[1]]
  test.select.table=test.data.total[[2]]
}




