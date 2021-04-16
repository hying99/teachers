##用于计算直接topdown方法和downtop方法的结果
##20170823

####第一步 将SVM的概率结果读入
file.middle="0"
file.type=""

#设置mat文件存储路径
setwd(paste(data.path,"//matfile",sep = ""))
mat.file=paste(file.prefix,file.middle,file.type,"_decision.mat",sep = "")
probability.data=readMat(mat.file,fixNames = FALSE)
prob.for.genes=probability.data$decision_test

####第二步 进行TPR规则处理

#用于检测结果是否符合TPR规则
#violate.detect.result=ViolateDetectprob(go.for.level, go.leaf.nodes,nodes.to.index,nodes.to.children,prob.for.genes)
#output.prob=DownTopStep(go.for.level,go.leaf.nodes,nodes.to.index,nodes.to.children,prob.for.genes)
#topdown.prob=TopDownStep(go.for.level,go.leaf.nodes,nodes.to.index,nodes.to.children,prob.for.genes)
output.prob=NaiveDownTop(go.for.level,go.leaf.nodes,nodes.to.index,nodes.to.parents,prob.for.genes)
#output.prob=prob.for.genes
violate.detect.top=ViolateDetectprob(go.for.level, go.leaf.nodes,nodes.to.index,nodes.to.children,output.prob)




####第三步 计算处理结果的评价指标

test.en=TRUE
if(test.en==TRUE)
{
  predict.labels=matrix(0,nrow(test.select.table),ncol(test.select.table))
  predict.scores=matrix(0,nrow(test.select.table),ncol(test.select.table))
  rownames(predict.labels)=rownames(test.select.table)
  colnames(predict.labels)=colnames(test.select.table)
  rownames(predict.scores)=rownames(test.select.table)
  colnames(predict.scores)=colnames(test.select.table)
  for(i in 1:nrow(output.prob))
  {
    for(j in 1:length(except.root.labels))
    {
      predict.scores[i,j]=output.prob[i,(2*j-1)]
      if(output.prob[i,(2*j-1)]>0.5)
      {
        predict.labels[i,except.root.labels[j]]=1
      }
    }
  }
  
  test.label.index=nodes.to.index[colnames(test.select.table)]
  measure.result=MHevaluate(predict.labels,test.select.table)
  F.each.class=F.measure.single.over.classes(test.select.table, predict.labels)
  prauc_result=PRAUCCalculate(predict.scores,test.select.table)
  
}