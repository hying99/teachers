setwd(data.path)#设置工作路径
#将SVM的概率结果读入
file.type=""
file.middle="0"
#setwd("H://R//DATA//matfile")#设置mat文件存储路径
mat.file=paste(file.prefix,file.middle,file.type,"_decision.mat",sep = "")
probability.data=readMat(mat.file,fixNames = FALSE)
#probability.data=readMat("cellcycle0replaceprocess_decision.mat",fixNames = FALSE)
#prob.for.genes=probability.data$decision
prob.for.genes=probability.data$d_decision_test

#fm.data=readMat("fm_t.mat",fixNames = FALSE)
#fm.values=fm.data$FM_t
# downtop.prob=prob.for.genes
# topdown.prob=prob.for.genes
#setwd("H://R//DATA")#设置工作路径
go.for.level.3=go.for.each.level[1:total.levels]#选取前三层节点进行处理
go.leaf.nodes.3=GetLeafNode1(graph.select.node.3)
#用于生产各节点的编号，以及节点与子节点的编号映射列表
total.index=MakeIndex(except.root.labels.3)
nodes.to.index=total.index[[1]]
nodes.to.children=total.index[[2]]

violate.detect.result=ViolateDetect(go.for.level.3, go.leaf.nodes.3,nodes.to.index,nodes.to.children,prob.for.genes)
#TPR 两步计算公式
downtop.w.prob=DownTopParent(go.for.level.3,go.leaf.nodes.3,nodes.to.index,nodes.to.children,prob.for.genes,each.go.weight)
topdown.w.prob=TopDownStep(go.for.level.3,go.leaf.nodes.3,nodes.to.index,nodes.to.children,downtop.w.prob)

downtop.prob=DownTopStep(go.for.level.3,go.leaf.nodes.3,nodes.to.index,nodes.to.children,prob.for.genes)
topdown.prob=TopDownStep(go.for.level.3,go.leaf.nodes.3,nodes.to.index,nodes.to.children,downtop.prob)
violate.detect.down=ViolateDetect(go.for.level.3, go.leaf.nodes.3,nodes.to.index,nodes.to.children,downtop.prob)
violate.detect.top=ViolateDetect(go.for.level.3, go.leaf.nodes.3,nodes.to.index,nodes.to.children,topdown.prob)
valid.en=FALSE
if(valid.en==TRUE)
{
  predict.labels=matrix(0,nrow(valid.select.table),ncol(valid.select.table))
  predict.scores=matrix(0,nrow(valid.select.table),ncol(valid.select.table))
  rownames(predict.labels)=rownames(valid.select.table)
  colnames(predict.labels)=colnames(valid.select.table)
  rownames(predict.scores)=rownames(valid.select.table)
  colnames(predict.scores)=colnames(valid.select.table)
  for(i in 1:nrow(topdown.prob))
  {
    for(j in 1:length(except.root.labels.3))
    {
      predict.scores[i,j]=topdown.prob[i,(2*j-1)]
      if(topdown.prob[i,(2*j-1)]>0.5)
      {
        predict.labels[i,except.root.labels.3[j]]=1
      }
    }
  }
  
  valid.label.index=nodes.to.index[colnames(valid.select.table)]
  measure.result=MHevaluate(predict.labels,valid.select.table)
  F.each.class=F.measure.single.over.classes(valid.select.table, predict.labels)
  prauc_result=PRAUCCalculate(predict.scores,valid.select.table)
  
}
test.en=TRUE
if(test.en==TRUE)
{
  predict.labels=matrix(0,nrow(test.select.table),ncol(test.select.table))
  predict.scores=matrix(0,nrow(test.select.table),ncol(test.select.table))
  rownames(predict.labels)=rownames(test.select.table)
  colnames(predict.labels)=colnames(test.select.table)
  rownames(predict.scores)=rownames(test.select.table)
  colnames(predict.scores)=colnames(test.select.table)
  for(i in 1:nrow(topdown.prob))
  {
    for(j in 1:length(except.root.labels.3))
    {
      predict.scores[i,j]=topdown.prob[i,(2*j-1)]
      if(topdown.prob[i,(2*j-1)]>0.5)
      {
        predict.labels[i,except.root.labels.3[j]]=1
      }
    }
  }
  
  test.label.index=nodes.to.index[colnames(test.select.table)]
  measure.result=MHevaluate(predict.labels,test.select.table)
  F.each.class=F.measure.single.over.classes(test.select.table, predict.labels)
  prauc_result=PRAUCCalculate(predict.scores,test.select.table)
  
}

if(test.en==TRUE)
{
  predict.w.labels=matrix(0,nrow(test.select.table),ncol(test.select.table))
  predict.w.scores=matrix(0,nrow(test.select.table),ncol(test.select.table))
  rownames(predict.w.labels)=rownames(test.select.table)
  colnames(predict.w.labels)=colnames(test.select.table)
  rownames(predict.w.scores)=rownames(test.select.table)
  colnames(predict.w.scores)=colnames(test.select.table)
  for(i in 1:nrow(topdown.w.prob))
  {
    for(j in 1:length(except.root.labels.3))
    {
      predict.w.scores[i,j]=topdown.w.prob[i,(2*j-1)]
      if(topdown.w.prob[i,(2*j-1)]>0.5)
      {
        predict.w.labels[i,except.root.labels.3[j]]=1
      }
    }
  }
  
  test.label.index.w=nodes.to.index[colnames(test.select.table)]
  measure.result.w=MHevaluate(predict.w.labels,test.select.table)
  F.each.class.w=F.measure.single.over.classes(test.select.table, predict.w.labels)
  prauc_result.w=PRAUCCalculate(predict.w.scores,test.select.table)
  
}
result.output.en=TRUE
result.savepath="H://R//DATA//result"
setwd(result.savepath)
today <-Sys.Date()
run.num=1
output.fname=paste(file.prefix,file.middle,file.type,"_result",".txt",sep = "")
if(result.output.en==TRUE)
{
  
  write.table(today,file=output.fname,sep = " , ",eol="\n",quote=FALSE,row.names = FALSE,col.names = FALSE,append = FALSE)
  write.table(measure.result,file=output.fname,sep = " , ",eol="\n",quote=FALSE,row.names = FALSE,col.names = TRUE,append = TRUE)
  write.table(measure.result.w,file=output.fname,sep = " ,",eol="\n",quote=FALSE,row.names = FALSE,col.names = TRUE,append = TRUE)
  write.table(prauc_result,file=output.fname,sep = " , ",eol="\n",quote=FALSE,row.names = FALSE,col.names = TRUE,append = TRUE)
  write.table(prauc_result.w,file=output.fname,sep = " , ",eol="\n",quote=FALSE,row.names = FALSE,col.names = TRUE,append = TRUE)
  
}

# setwd("H://R//DATA//fm")
# fm.data=readMat("fm_result.mat",fixNames = FALSE)
# fm.values=fm.data$FM
# rec.data=readMat("fm_result.mat",fixNames = FALSE)
# rec.values=fm.data$FM
# pre.data=readMat("fm_result.mat",fixNames = FALSE)
# pre.values=fm.data$FM
# setwd("H://R//DATA")
# 
# each.node.pre=F.each.class$per.class[,1]
# each.node.rec=F.each.class$per.class[,2]
# each.node.fm=F.each.class$per.class[,4]




