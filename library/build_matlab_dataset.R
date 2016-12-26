#用于生成MIML中所需训练及测试文件的函数
BuildMatlabDataset<-function (file.prefix, factor.col, work.path,file.savepath,delete.outlier,
                              replace.outlier, NAreplace,Zrescale, write.data.enable,
                              write.class.enable)
{
  setwd(work.path)
  matrix.cellcycle=ReadData(paste("originaldata//",file.prefix,"0.train",sep = ""),factor.col = factor.col)
  matrix.cellcycle.data=matrix.cellcycle[[1]]#基因的数据信息
  
  training.cellcycle.data=matrix.cellcycle.data
  #training.cellcycle.name=rownames(training.cellcycle.data)#基因名称列表
  #对输入数据剔除异常值，并进行归一化处理
  
  trainscale.result=TraindataScale(training.cellcycle.data,factor.col,delete.outlier=delete.outlier,replace.outlier = replace.outlier,NAreplace=NAreplace,Zrescale=Zrescale)
  remain.data=trainscale.result[[1]]
  sp=trainscale.result[[2]]
  
  #得到几个数据集共有的基因的名称列表
  common.genes <- Get.all.common.genes(go.general.table.BP, remain.data)
  
  remain.select.data=remain.data[common.genes,]
  #得到common genes中每个基因对应的全部GO标签列表
  match.go.general=go.general.list.BP[common.genes]
  
  #将标签列表转换为TABLE形式，行名称为基因名称，列名称为GO标签
  match.go.table=Build.GO.class.labels(match.go.general)
  
  
  #得到目前的所有基因共有多少个不重复的GO标签
  all.go.labels=Get.classes(match.go.general)
  
  
  #得到DAG图中包含50个样本以上的节点列表及对应的基因名称
  go.label.list=DataCleaning(match.go.general,match.go.table,select.num = 100)
  all.go.labels.50=Get.classes(go.label.list)
  
  
  #得到修改后的DAG图中各节点的层级
  graph.BP.general.50 <- subGraph(all.go.labels.50, BP.univ.graph)
  graph.BP.level.50=GraphLevel(graph.BP.general.50)
  go.level.statistics=LevelStatistics(graph.BP.level.50)
  go.for.each.level=go.level.statistics[[1]]
  each.level.nodes.num=go.level.statistics[[2]]
  total.levels=length(go.for.each.level)
  
  
  select.node.3=NodeSelectByLevel(go.level.statistics,total.levels,add.root.node = TRUE)#选择根结点至第三层的所有go节点
  except.root.labels.3=setdiff(select.node.3,"GO:0008150")
  #graph.select.node.3 <- subGraph(select.node.3, BP.univ.graph)
  #PlotLabelGraph(except.root.labels.3,BP.univ.graph,num.only=FALSE,plot.en = TRUE,output.en = FALSE,write.pic.name = "go.graph.level.ps")
  #each.go.level.num=graph.BP.level.50[except.root.labels.3]
  
  
  for (i in 1:length(go.label.list))
  {
    go.label.list[[i]]=intersect(go.label.list[[i]],except.root.labels.3)
  }
  
  except.root.table=Build.GO.class.labels(go.label.list)
  
  if(write.data.enable)
  {
    setwd(file.savepath)
    write.data.fname=paste(file.prefix,"0_traindataset.csv",sep = "")
    
    write.table(remain.select.data,file=write.data.fname,sep = ",",eol="\n",quote=FALSE,na="0",row.names = FALSE,col.names = FALSE)
    
  }
  if(write.class.enable)
  {
    setwd(file.savepath)
    write.class.fname=paste(file.prefix,"0_trainclass.csv",sep = "")
    write.table(except.root.table,file=write.class.fname,sep = ",",eol="\n",quote=FALSE,na="0",row.names = FALSE,col.names = FALSE)
  }
  
  
  setwd(work.path)
  original.valid.file=paste("originaldata//",file.prefix,"0.valid",sep = "")
  select.attributes.en=FALSE
  
  
  if(select.attributes.en==FALSE)
  {
    
    write.data.fname=paste(file.prefix,"0_validdataset.csv",sep = "")
    write.class.fname=paste(file.prefix,"0_validclass.csv",sep = "")
    valid.cellcycle=ReadData(original.valid.file,factor.col = factor.col)#读入valid基因特征属性
    valid.cellcycle.data=valid.cellcycle[[1]]#valid基因的数据信息
    valid.scaled.data=ValiddataScale(valid.cellcycle.data,factor.col,sp,replace.outlier=replace.outlier,NAreplace=NAreplace,Zrescale=Zrescale)
    setwd(file.savepath)
    valid.data.total=BuildValidset(valid.scaled.data,go.general.table.BP,go.general.list.BP,except.root.labels.3,
                                   write.data.enable=TRUE,write.class.enable=TRUE,write.data.fname=write.data.fname,
                                   write.class.fname=write.class.fname,select.attributes.en=select.attributes.en,select.attributes)
    valid.select.data=valid.data.total[[1]]
    valid.select.table=valid.data.total[[2]]
  }
  
  setwd(work.path)
  original.test.file=paste("originaldata//",file.prefix,"0.propertest",sep = "")
  select.attributes.en=FALSE
  if(select.attributes.en==FALSE)
  {
    
    write.data.fname=paste(file.prefix,"0_testdataset.csv",sep = "")
    write.class.fname=paste(file.prefix,"0_testclass.csv",sep = "")
    test.cellcycle=ReadData(original.test.file,factor.col = factor.col)#读入test基因特征属性
    test.cellcycle.data=test.cellcycle[[1]]#test基因的数据信息
    test.scaled.data=ValiddataScale(test.cellcycle.data,factor.col,sp,replace.outlier=replace.outlier,NAreplace=NAreplace,Zrescale=Zrescale)
    setwd(file.savepath)
    test.data.total=BuildValidset(test.scaled.data,go.general.table.BP,go.general.list.BP,except.root.labels.3,
                                  write.data.enable=TRUE,write.class.enable=TRUE,write.data.fname=write.data.fname,
                                  write.class.fname=write.class.fname,select.attributes.en=select.attributes.en,select.attributes)
    test.select.data=test.data.total[[1]]
    test.select.table=test.data.total[[2]]
  }
  
}



