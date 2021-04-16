####采用neural netwrok作为基础分类器，调用了rebuilddataprocess.R和neuralnetwork.R文件
####采用topdown训练方式，每个节点的训练数据都增加了其祖先节点的预测值作为属性
####采用一个文件，统一得到所有数据集的结果，并存储为txt文件
####20171116
work.path="D://R//CODES"
data.path="D://R//DATA"
setwd(work.path)#首先将工作路径设置为代码的路径
source("headfile.R")

####第三步 得到BO结构信息及注释信息
# ontology.sel="BP"
# univ.graph <- Build.universal.graph.ontology.down(ontology = ontology.sel)#得到BP全图
# annotation.final=AnnotationFinal(ontology.sel)#提取BP结构下的所有基因及注释信息
# go.general.list=Get.GO.all.classes(annotation.final)#根据基因的注释GO标签得到基因的全部GO标签
# go.general.table=Build.GO.class.labels(go.general.list)#生成基因及注释信息数据表
#首先选择使用的数据集，1 cellcycle 2 derisi 3 eisen 4 gasch1 5 gasch2 6 church 7 spo 8 seq 9 struc 10 hom
dataset.result=DatasetSelect(dataset.index = 1)
file.prefix=dataset.result[[1]]#得到数据集前缀名称
cat('The processing dataset is ' ,file.prefix,"\n")
setwd(work.path)#首先将工作路径设置为代码的路径
source("rebuilddataprocess.R")
setwd(work.path)#首先将工作路径设置为代码的路径
source("tdancestornn.R")

dataset.result=DatasetSelect(dataset.index = 2)
file.prefix=dataset.result[[1]]#得到数据集前缀名称
cat('The processing dataset is ' ,file.prefix,"\n")
setwd(work.path)#首先将工作路径设置为代码的路径
source("rebuilddataprocess.R")
setwd(work.path)#首先将工作路径设置为代码的路径
source("tdancestornn.R")

dataset.result=DatasetSelect(dataset.index = 3)
file.prefix=dataset.result[[1]]#得到数据集前缀名称
cat('The processing dataset is ' ,file.prefix,"\n")
setwd(work.path)#首先将工作路径设置为代码的路径
source("rebuilddataprocess.R")
setwd(work.path)#首先将工作路径设置为代码的路径
source("tdancestornn.R")

dataset.result=DatasetSelect(dataset.index = 4)
file.prefix=dataset.result[[1]]#得到数据集前缀名称
cat('The processing dataset is ' ,file.prefix,"\n")
setwd(work.path)#首先将工作路径设置为代码的路径
source("rebuilddataprocess.R")
setwd(work.path)#首先将工作路径设置为代码的路径
source("tdancestornn.R")

dataset.result=DatasetSelect(dataset.index = 5)
file.prefix=dataset.result[[1]]#得到数据集前缀名称
cat('The processing dataset is ' ,file.prefix,"\n")
setwd(work.path)#首先将工作路径设置为代码的路径
source("rebuilddataprocess.R")
setwd(work.path)#首先将工作路径设置为代码的路径
source("tdancestornn.R")

dataset.result=DatasetSelect(dataset.index = 6)
file.prefix=dataset.result[[1]]#得到数据集前缀名称
cat('The processing dataset is ' ,file.prefix,"\n")
setwd(work.path)#首先将工作路径设置为代码的路径
source("rebuilddataprocess.R")
setwd(work.path)#首先将工作路径设置为代码的路径
source("tdancestornn.R")

dataset.result=DatasetSelect(dataset.index = 7)
file.prefix=dataset.result[[1]]#得到数据集前缀名称
cat('The processing dataset is ' ,file.prefix,"\n")
setwd(work.path)#首先将工作路径设置为代码的路径
source("rebuilddataprocess.R")
setwd(work.path)#首先将工作路径设置为代码的路径
source("tdancestornn.R")

dataset.result=DatasetSelect(dataset.index = 8)
file.prefix=dataset.result[[1]]#得到数据集前缀名称
cat('The processing dataset is ' ,file.prefix,"\n")
setwd(work.path)#首先将工作路径设置为代码的路径
source("rebuilddataprocess.R")
setwd(work.path)#首先将工作路径设置为代码的路径
source("tdancestornn.R")