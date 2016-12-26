#用于生产arff文件的程序源文件，运行前需先运行dataprocessing等文件，获得相应变量内容
#读入的arff文件用于提取arff属性设置表头

#original.valid.file="originaldata//cellcycle0.valid"
original.valid.file=paste("originaldata//",file.prefix,"0.valid",sep = "")
select.attributes.en=FALSE
write.data.fname=paste(file.prefix,"0_validdataset.csv",sep = "")
write.class.fname=paste(file.prefix,"0_validclass.csv",sep = "")
valid.cellcycle=ReadData(original.valid.file,factor.col = factor.col,NAreplace = TRUE,Zrescale = TRUE,ZNormal = FALSE)#读入valid基因特征属性
valid.cellcycle.data=valid.cellcycle[[1]]#valid基因的数据信息
#此处except.root.labels变量设置为select.node.3，包含根结点
valid.data.total=BuildValidset(valid.cellcycle.data,go.general.table.BP,go.general.list.BP,
                               select.node.3,write.data.enable=FALSE,write.class.enable=FALSE,
                               write.data.fname=write.data.fname, write.class.fname=write.class.fname,
                               select.attributes.en=select.attributes.en,select.attributes)
valid.select.data=valid.data.total[[1]]
valid.select.table=valid.data.total[[2]]
valid.go.label.list=valid.data.total[[3]]
valid.leaf.label.list=valid.go.label.list
# valid.leaf.label.list=list()
# for (i in 1:length(valid.go.label.list))
# {
#   valid.leaf.label.list[[i]]=intersect(valid.go.label.list[[i]],go.leaf.nodes.3)
#   
# }

#original.test.file="originaldata//cellcycle0.propertest"
original.test.file=paste("originaldata//",file.prefix,"0.propertest",sep = "")
write.data.fname=paste(file.prefix,"0_testdataset.csv",sep = "")
write.class.fname=paste(file.prefix,"0_testclass.csv",sep = "")
test.cellcycle=ReadData(original.test.file,factor.col = factor.col,NAreplace = NAreplace,Zrescale = Zrescale,ZNormal = ZNormal)#读入test基因特征属性
test.cellcycle.data=test.cellcycle[[1]]#test基因的数据信息
#此处except.root.labels变量设置为select.node.3，包含根结点
test.data.total=BuildValidset(test.cellcycle.data,go.general.table.BP,go.general.list.BP,
                               select.node.3,write.data.enable=FALSE,write.class.enable=FALSE,
                               write.data.fname=write.data.fname, write.class.fname=write.class.fname,
                               select.attributes.en=select.attributes.en,select.attributes)
test.select.data=test.data.total[[1]]
test.select.table=test.data.total[[2]]
test.go.label.list=test.data.total[[3]]
test.leaf.label.list=test.go.label.list
# test.leaf.label.list=list()
# for (i in 1:length(test.go.label.list))
# {
#   test.leaf.label.list[[i]]=intersect(test.go.label.list[[i]],go.leaf.nodes.3)
#   
# }

setwd("H://R//DATA//arff")

arff.result=WriteArff(original.filename=paste(file.prefix,"_GO.train.arff",sep = ""),write.filename= paste(file.prefix,"train.arff",sep = ""),
                      modified.data=training.cellcycle.data,common.genes=common.genes,ontology="BP",
                      class.label=go.label.list,class.graph=graph.select.node.3,genename.exist=FALSE) 

arff.result=WriteArff(original.filename=paste(file.prefix,"_GO.valid.arff",sep = ""),write.filename= paste(file.prefix,"valid.arff",sep = ""),
                      modified.data=valid.select.data,common.genes=NULL,ontology="BP",
                      class.label=valid.leaf.label.list,class.graph=graph.select.node.3,genename.exist=FALSE)

arff.result=WriteArff(original.filename=paste(file.prefix,"_GO.test.arff",sep = ""),write.filename= paste(file.prefix,"test.arff",sep = ""),
                      modified.data=test.select.data,common.genes=NULL,ontology="BP",
                      class.label=test.leaf.label.list,class.graph=graph.select.node.3,genename.exist=FALSE) 
