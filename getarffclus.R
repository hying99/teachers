#用于将数据集转化为arff文件的程序源文件，20170720编写
#其输出的arff文件仅供Clus软件进行预测使用
#程序内部进行了很多简化，去除了不必要的语句，故不能他用
#改程序可以处理带有factor因子的数据集，也可以处理一般数据集


#第一次使用，需先运行第二步，得到全局变量
####第一步设置代码存放路径 以及数据存放路径，引入库函数以及自定义的函数
work.path="H://R//CODES"
data.path="H://R//DATA"
setwd(work.path)#首先将工作路径设置为代码的路径
source("headfile.R")

####第二步 得到BO结构信息及注释信息
# ontology.sel="BP"
# univ.graph <- Build.universal.graph.ontology.down(ontology = ontology.sel)#得到BP全图
# annotation.final=AnnotationFinal(ontology.sel)#提取BP结构下的所有基因及注释信息
# go.general.list=Get.GO.all.classes(annotation.final)#根据基因的注释GO标签得到基因的全部GO标签
# go.general.table=Build.GO.class.labels(go.general.list)#生成基因及注释信息数据表


####第三步 确定待处理的数据集
#首先选择使用的数据集，1 cellcycle 2 derisi 3 eisen 4 gasch1 5 gasch2 6 church 7 spo 8 seq 9 struc 10 hom
dataset.result=DatasetSelect(dataset.index = 6)
file.prefix=dataset.result[[1]]#得到数据集前缀名称


####第四步 读取训练集数据
setwd(data.path)
train.filename=paste("originaldata//",file.prefix,"0.train",sep = "")
train.mydata <- scan(train.filename, what = list("", "", ""))#读入待处理数据,用于得到基因名称列表
#write.table(mydata[[1]],file="newdata",quote = FALSE,row.names = FALSE,col.names = FALSE)
train.filedata=read.table(train.filename,sep = ",",quote = "")#读入待处理数据,用于得到基因数据
train.matrix.data=as.matrix(train.filedata[,-dim(train.filedata)[2]])#去除样本标签这一列
rownames(train.matrix.data)=toupper(substring(train.mydata[[3]],3))#删除基因名开头的yt字母，并将基因名转为大写

####第五步 得到有注释信息的基因名称及其注释信息
#得到训练集共有的基因的名称列表
common.genes <- Get.all.common.genes(go.general.table, train.matrix.data)

remain.select.data=train.matrix.data[common.genes,]#此处即为新生成的训练集数据，此句不可删去，原始数据时仍需使用
#得到common genes中每个基因对应的全部GO标签列表
match.go.general=go.general.list[common.genes]

#将标签列表转换为TABLE形式，行名称为基因名称，列名称为GO标签
match.go.table=Build.GO.class.labels(match.go.general)

####第六步 对基因注释信息进行处理，根据GO标签注释的样本数选择适当的GO标签

#得到目前的所有基因共有多少个不重复的GO标签
all.go.labels=Get.classes(match.go.general)

#得到DAG图中包含select.num=100个样本以上的节点列表及对应的基因名称
go.label.list=DataCleaning(match.go.general,match.go.table,select.num = 100)
all.go.labels.sel=Get.classes(go.label.list)
#except.root.labels=setdiff(all.go.labels.50,"GO:0008150")#除去根结点后所剩节点

####第七步 根据选择的GO节点生成DAG图，并得到DAG图中的层级信息

#得到修改后的DAG图中各节点的层级
graph.general.sel <- subGraph(all.go.labels.sel, univ.graph)
graph.level.sel=GraphLevel(graph.general.sel)
go.level.statistics=LevelStatistics(graph.level.sel)
go.for.each.level=go.level.statistics[[1]]
total.levels=length(go.for.each.level)


####第八步 按照层级选择一定数量的节点，生成这些节点的子图,用于后面的WriteArff函数

#选择根结点至total.levels的所有go节点，select.node为向量形式
select.node=NodeSelectByLevel(go.level.statistics,total.levels,add.root.node = TRUE)
except.root.labels=setdiff(select.node,"GO:0008150")#去掉所选节点中的根结点，剩余节点的集合为需求的标签
sub.graph <- subGraph(select.node, univ.graph)#绘制子图

####第九步 读入验证集数据
valid.filename=paste("originaldata//",file.prefix,"0.valid",sep = "")
valid.mydata <- scan(valid.filename, what = list("", "", ""))#读入待处理数据,用于得到基因名称列表
#write.table(mydata[[1]],file="newdata",quote = FALSE,row.names = FALSE,col.names = FALSE)
valid.filedata=read.table(valid.filename,sep = ",",quote = "")#读入待处理数据,用于得到基因数据
valid.matrix.data=as.matrix(valid.filedata[,-dim(valid.filedata)[2]])#去除样本标签这一列
rownames(valid.matrix.data)=toupper(substring(valid.mydata[[3]],3))#删除基因名开头的yt字母，并将基因名转为大写


####第十步 对验证集数据进行处理
#此处except.root.labels变量设置为select.node，包含根结点
write.data.fname=paste(file.prefix,"0_validdataset.csv",sep = "")
write.class.fname=paste(file.prefix,"0_validclass.csv",sep = "")
valid.data.total=BuildValidset(valid.matrix.data,go.general.table,go.general.list,select.node,
                               write.data.enable=TRUE,write.class.enable=TRUE,write.data.fname=write.data.fname,
                               write.class.fname=write.class.fname,select.attributes.en=FALSE,select.attributes)
valid.select.data=valid.data.total[[1]]
valid.select.table=valid.data.total[[2]]
valid.go.label.list=valid.data.total[[3]]
#直接将所有的标签传给每个样本作为标签，而不是像原文件一样只选择最具体的标签
valid.leaf.label.list=valid.go.label.list


####第十一步 读入测试集数据
test.filename=paste("originaldata//",file.prefix,"0.propertest",sep = "")
test.mydata <- scan(test.filename, what = list("", "", ""))#读入待处理数据,用于得到基因名称列表
#write.table(mydata[[1]],file="newdata",quote = FALSE,row.names = FALSE,col.names = FALSE)
test.filedata=read.table(test.filename,sep = ",",quote = "")#读入待处理数据,用于得到基因数据
test.matrix.data=as.matrix(test.filedata[,-dim(test.filedata)[2]])#去除样本标签这一列
rownames(test.matrix.data)=toupper(substring(test.mydata[[3]],3))#删除基因名开头的yt字母，并将基因名转为大写

####第十二步 对测试集数据进行处理
write.data.fname=paste(file.prefix,"0_testdataset.csv",sep = "")
write.class.fname=paste(file.prefix,"0_testclass.csv",sep = "")
test.data.total=BuildValidset(test.matrix.data,go.general.table,go.general.list,select.node,
                              write.data.enable=TRUE,write.class.enable=TRUE,write.data.fname=write.data.fname,
                              write.class.fname=write.class.fname,select.attributes.en=FALSE,select.attributes)
test.select.data=test.data.total[[1]]
test.select.table=test.data.total[[2]]
test.go.label.list=test.data.total[[3]]
test.leaf.label.list=test.go.label.list


####第十三步 生成新的arff文件
setwd(paste(data.path,"//arff",sep = ""))

arff.result=WriteArff(original.filename=paste(file.prefix,"_GO.train.arff",sep = ""),write.filename= paste(file.prefix,"train_clus.arff",sep = ""),
                      modified.data=remain.select.data,common.genes=common.genes,ontology=ontology.sel,
                      class.label=go.label.list,class.graph=sub.graph,genename.exist=FALSE) 

arff.result=WriteArff(original.filename=paste(file.prefix,"_GO.valid.arff",sep = ""),write.filename= paste(file.prefix,"valid_clus.arff",sep = ""),
                      modified.data=valid.select.data,common.genes=NULL,ontology=ontology.sel,
                      class.label=valid.leaf.label.list,class.graph=sub.graph,genename.exist=FALSE)

arff.result=WriteArff(original.filename=paste(file.prefix,"_GO.test.arff",sep = ""),write.filename= paste(file.prefix,"test_clus.arff",sep = ""),
                      modified.data=test.select.data,common.genes=NULL,ontology=ontology.sel,
                      class.label=test.leaf.label.list,class.graph=sub.graph,genename.exist=FALSE) 
