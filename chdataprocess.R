##在rebuilddataprocess基础上增加了生成生成train.hf部分
##该文件用于支持chloss2文件以及DAGlabel文件
##文件的运行流程为rebuilddataprocess+chloss2+DAGlabel
##20190308创建

##
####第一步设置代码存放路径 以及数据存放路径
work.path="D://R//CODES"
data.path="D://R//DATA"
setwd(work.path)#首先将工作路径设置为代码的路径


####第二步 引入库函数以及自定义的函数
source("headfile.R")

####第三步 得到GO结构信息及注释信息
ontology.sel="BP"
univ.graph <- Build.universal.graph.ontology.down(ontology = ontology.sel)#得到BP全图
annotation.final=AnnotationFinal(ontology.sel)#提取BP结构下的所有基因及注释信息
go.general.list=Get.GO.all.classes(annotation.final)#根据基因的注释GO标签得到基因的全部GO标签
go.general.table=Build.GO.class.labels(go.general.list)#生成基因及注释信息数据表

####第四步 读入待处理的训练数据
#首先选择使用的数据集，1 cellcycle 2 derisi 3 eisen 4 gasch1 5 gasch2 6 church 7 spo 8 seq 9 struc 10 hom

dataset.result=DatasetSelect(dataset.index = 3)
file.prefix=dataset.result[[1]]#得到数据集前缀名称
factor.col.index=dataset.result[[2]]#得到内容为类别信息的列号
factor.col.num=dataset.result[[3]]#得到factor列转化为数值后各列的总数
factor.levels=dataset.result[[4]]
#读取数据时，01调整Zrescale和标准化ZNormal不能同时使用
setwd(data.path)#将工作路径设置为存储训练训练数据的路径，以便读取训练数据
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

####第五步 对输入数据剔除异常值，并进行归一化处理
#delete.outlier 与replace.outlier no.del.replace只能有一个为TRUE，如果均为FALSE，则进行异常值上下界替换
#选择为no.del.replace时，不进行任何替换操作
#设置scale函数的操作选项，在此处保证对训练集 验证集 测试集的设置内容均相同
delete.outlier=FALSE
replace.outlier = FALSE
no.del.replace = FALSE
NAreplace=TRUE
Zrescale=TRUE
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

#得到DAG图中包含select.num=100个样本以上的节点列表及对应的GO名称
go.label.list=DataCleaning(match.go.general,match.go.table,select.num = 100)
all.go.labels.sel=Get.classes(go.label.list)
#except.root.labels=setdiff(all.go.labels.50,"GO:0008150")#除去根结点后所剩节点

####第八步 根据选择的GO节点生成DAG图，并得到DAG图中的层级信息

#得到修改后的DAG图中各节点的层级
graph.general.sel <- subGraph(all.go.labels.sel, univ.graph)
graph.level.sel=GraphLevel(graph.general.sel)
go.level.statistics=LevelStatistics(graph.level.sel)
go.for.each.level=go.level.statistics[[1]]
each.level.nodes.num=go.level.statistics[[2]]
total.levels=length(go.for.each.level)

####第九步 按照层级选择一定数量的节点，生成并绘制这些节点的子图，给出这些节点所属的层级

#选择根结点至total.levels的所有go节点，select.node为向量形式
select.node=NodeSelectByLevel(go.level.statistics,total.levels,add.root.node = TRUE)
except.root.labels=setdiff(select.node,"GO:0008150")#去掉所选节点中的根结点，剩余节点的集合为需求的标签
sub.graph <- subGraph(select.node, univ.graph)#绘制子图

graph.new=PlotLabelGraph(except.root.labels,univ.graph,num.only=TRUE,plot.en = TRUE,output.en = FALSE,write.pic.name = "go.graph.level.ps")
each.go.level.num=graph.level.sel[except.root.labels]#给出每个节点所属于的层级

####第十步 得到所选节点的层级信息列表，叶子节点集合，父节点、子节点、祖先节点以及子孙节点列表
####       此步产生的信息用于分类结果的后处理

go.leaf.nodes=GetLeafNode1(sub.graph)#得到子图中的叶子节点
#按照层级选择一定数量的节点，go.for.level所选层级与select.node相同，内容也相同，但存储为list形式
go.for.level=go.for.each.level[1:total.levels]
#用于生产各节点的编号，以及节点与子节点的编号映射列表
total.index=MakeIndex(except.root.labels)
nodes.to.index=total.index[[1]]
nodes.to.children=total.index[[2]]
nodes.to.ancestors=total.index[[3]]
nodes.to.parents=total.index[[4]]
nodes.to.descendants=total.index[[5]]
 
#用于生成train.hf

with.root.labels=c(except.root.labels,"GO:0008150")
withroot.total.index=MakeIndex(with.root.labels)
withroot.nodes.to.children=withroot.total.index[[2]]
withroot.nodes.to.parents=withroot.total.index[[4]]

#这段代码是产生TPR算法中权值信息所需，目前暂不使用 20170311
# each.go.weight=unname(each.go.level.num)
# for (i in 1:length(each.go.level.num))
# {
#   each.go.weight[i]=(total.levels+1-each.go.weight[i])/(total.levels+1)
# }



# for (i in 1:length(go.label.list))
# {
#   go.label.list[[i]]=intersect(go.label.list[[i]],select.node.3)
# }
# root.table.3=Build.GO.class.labels(go.label.list)
####第十一步 生成构建训练集所需的GO标签，并且计算平均每个样本所含有的标签数量
#产生构建训练集所需的GO标签
except.root.table=match.go.table[,except.root.labels]
root.table=match.go.table[,select.node]
#训练集样本数量
sample.num=nrow(root.table)
#平均每个样本所含有的标签数量
average.label=sum(root.table)/sample.num

####第十二步 生成训练集，若想用不同方法生成训练集，则在此处替换BuildTrainDataset函数
#将工作路径改为将要存放生成的训练集csv文件的文件夹
#setwd("H://R//DATA//traindata")
setwd(paste(data.path,"//traindata",sep = ""))
select.attributes.en=FALSE#不进行属性选择

data.total=BuildTrainDataset(root.table, except.root.labels, data.matrix=remain.select.data,
                             ontology = ontology.sel, adjust.ratio=0.2,ratio.negative = 0, common.genes = common.genes,
                             seed = 1,select.attributes.en=select.attributes.en,write.en=TRUE)


####第十三步 生成验证集

setwd(data.path)
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


#归一化
  
valid.scaled.data=ValiddataScale(valid.data,factor.col.num,sp,replace.outlier=replace.outlier,
                                   no.del.replace = no.del.replace,NAreplace=NAreplace,Zrescale=Zrescale)
  
valid.data.total=BuildValidset(valid.scaled.data,go.general.table,go.general.list,except.root.labels,
                                 write.data.enable=TRUE,write.class.enable=TRUE,write.data.fname=write.data.fname,
                                 write.class.fname=write.class.fname,select.attributes.en=select.attributes.en,select.attributes)
valid.select.data=valid.data.total[[1]]
valid.select.table=valid.data.total[[2]]

####第十四步 生成测试集
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
test.data.total=BuildValidset(test.scaled.data,go.general.table,go.general.list,except.root.labels,
                                write.data.enable=TRUE,write.class.enable=TRUE,write.data.fname=write.data.fname,
                                write.class.fname=write.class.fname,select.attributes.en=select.attributes.en,select.attributes)
test.select.data=test.data.total[[1]]
test.select.table=test.data.total[[2]]





