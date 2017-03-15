#用于deris数据的处理，训练前期数据处理，用于生成合适的训练集以及测试集

####第一步设置代码存放路径 以及数据存放路径
work.path="H://R//CODES"
data.path="H://R//DATA"
setwd(work.path)#首先将工作路径设置为代码的路径


####第二步 引入库函数以及自定义的函数
source("headfile.R")

####第三步 得到BO结构信息及注释信息
BP.univ.graph <- Build.universal.graph.ontology.down(ontology = "BP")#得到BP全图
annotation.final.BP=AnnotationFinal("BP")#提取BP结构下的所有基因及注释信息
go.general.list.BP=Get.GO.all.classes(annotation.final.BP)#根据基因的注释GO标签得到基因的全部GO标签
go.general.table.BP=Build.GO.class.labels(go.general.list.BP)#生成基因及注释信息数据表

####第四步 读入待处理的训练数据
file.prefix="derisi"
#matrix.cellcycle=ReadData("originaldata//cellcycle0.train")#读入基因特征属性
#读取数据时，01调整Zrescale和标准化ZNormal不能同时使用
factor.col=c(0)

setwd(data.path)#将工作路径设置为存储训练训练数据的路径，以便读取训练数据
matrix.cellcycle=ReadData(paste("originaldata//",file.prefix,"0.train",sep = ""),factor.col = factor.col)
matrix.cellcycle.data=matrix.cellcycle[[1]]#基因的数据信息

training.cellcycle.data=matrix.cellcycle.data

####第五步 对输入数据剔除异常值，并进行归一化处理
trainscale.result=TraindataScale(training.cellcycle.data,factor.col,delete.outlier=FALSE,replace.outlier = FALSE,NAreplace=TRUE,Zrescale=TRUE)
remain.data=trainscale.result[[1]]
sp=trainscale.result[[2]]

####第六步 得到目前有注释信息的基因名称，并将其转换为list及table形式

#得到几个数据集共有的基因的名称列表
common.genes <- Get.all.common.genes(go.general.table.BP, remain.data)
#得到common genes中每个基因对应的全部GO标签列表
match.go.general=go.general.list.BP[common.genes]

#将标签列表转换为TABLE形式，行名称为基因名称，列名称为GO标签
match.go.table=Build.GO.class.labels(match.go.general)

####第七步 对基因注释信息进行处理，根据GO标签注释的样本数选择适当的GO标签

#得到目前的所有基因共有多少个不重复的GO标签
all.go.labels=Get.classes(match.go.general)

#得到DAG图中包含100个样本以上的节点列表及对应的基因名称
go.label.list=DataCleaning(match.go.general,match.go.table,select.num = 100)
all.go.labels.50=Get.classes(go.label.list)
#except.root.labels=setdiff(all.go.labels.50,"GO:0008150")#除去根结点后所剩节点

####第八步 根据选择的GO节点生成DAG图，并得到DAG图中的层级信息

#得到修改后的DAG图中各节点的层级
graph.BP.general.50 <- subGraph(all.go.labels.50, BP.univ.graph)
graph.BP.level.50=GraphLevel(graph.BP.general.50)
go.level.statistics=LevelStatistics(graph.BP.level.50)
go.for.each.level=go.level.statistics[[1]]
each.level.nodes.num=go.level.statistics[[2]]
total.levels=length(go.for.each.level)

####第九步 按照层级选择一定数量的节点，生成并绘制这些节点的子图，给出这些节点所属的层级

#选择根结点至total.levels的所有go节点，select.node.3为向量形式
select.node.3=NodeSelectByLevel(go.level.statistics,total.levels,add.root.node = TRUE)
except.root.labels.3=setdiff(select.node.3,"GO:0008150")#去掉所选节点中的根结点，剩余节点的集合为需求的标签
graph.select.node.3 <- subGraph(select.node.3, BP.univ.graph)#绘制子图

PlotLabelGraph(except.root.labels.3,BP.univ.graph,num.only=FALSE,plot.en = TRUE,output.en = FALSE,write.pic.name = "go.graph.level.ps")
each.go.level.num=graph.BP.level.50[except.root.labels.3]#给出每个节点所属于的层级

####第十步 得到所选节点的层级信息列表，叶子节点集合，父节点、子节点、祖先节点以及子孙节点列表
####       此步产生的信息用于分类结果的后处理

go.leaf.nodes.3=GetLeafNode1(graph.select.node.3)#得到子图中的叶子节点
#按照层级选择一定数量的节点，go.for.level.3所选层级与select.node.3相同，内容也相同，但存储为list形式
go.for.level.3=go.for.each.level[1:total.levels]
#用于生产各节点的编号，以及节点与子节点的编号映射列表
total.index=MakeIndex(except.root.labels.3)
nodes.to.index=total.index[[1]]
nodes.to.children=total.index[[2]]
nodes.to.ancestors=total.index[[3]]
nodes.to.parents=total.index[[4]]
nodes.to.descendants=total.index[[5]]

#这段代码是产生TPR算法中权值信息所需，目前暂不使用 20170311
# each.go.weight=unname(each.go.level.num)
# for (i in 1:length(each.go.level.num))
# {
#   each.go.weight[i]=(total.levels+1-each.go.weight[i])/(total.levels+1)
# }

####第十一步 生成训练集，若想用不同方法生成训练集，则在此处替换BuildTrainDataset函数

# for (i in 1:length(go.label.list))
# {
#   go.label.list[[i]]=intersect(go.label.list[[i]],select.node.3)
# }
# root.table.3=Build.GO.class.labels(go.label.list)
#产生构建训练集所需的GO标签
root.table.3=match.go.table[,except.root.labels.3]
#将工作路径改为将要存放生成的训练集csv文件的文件夹
#setwd("H://R//DATA//traindata")
setwd(paste(data.path,"//traindata",sep = ""))
select.attributes.en=FALSE#不进行属性选择

data.total=BuildTrainDataset(root.table.3, except.root.labels.3, data.matrix=remain.data,
                             ontology = "BP", adjust.ratio=0.2,ratio.negative = 0, common.genes = common.genes,
                             seed = 1,select.attributes.en=select.attributes.en,write.en=TRUE)


####第十二步 生成验证集和测试集

setwd(data.path)
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




