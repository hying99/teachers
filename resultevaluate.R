##采用TPR方法对结果进行后处理的文件 20170311

####第一步 将SVM的概率结果读入

file.type=""
file.middle="0"
#设置mat文件存储路径
setwd(paste(data.path,"//matfile",sep = ""))
mat.file=paste(file.prefix,file.middle,file.type,"_decision.mat",sep = "")
probability.data=readMat(mat.file,fixNames = FALSE)
#probability.data=readMat("cellcycle0replaceprocess_decision.mat",fixNames = FALSE)
#prob.for.genes=probability.data$decision
prob.for.genes=probability.data$decision_test

####第二步 进行TPR规则处理

#用于检测结果是否符合TPR规则
violate.detect.result=ViolateDetectprob(go.for.level, go.leaf.nodes,nodes.to.index,nodes.to.children,prob.for.genes)
#TPR 两步计算公式
# downtop.w.prob=DownTopParent(go.for.level.3,go.leaf.nodes.3,nodes.to.index,nodes.to.children,prob.for.genes,each.go.weight)
# topdown.w.prob=TopDownStep(go.for.level.3,go.leaf.nodes.3,nodes.to.index,nodes.to.children,downtop.w.prob)


downtop.prob=DownTopStep(go.for.level,go.leaf.nodes,nodes.to.index,nodes.to.children,prob.for.genes)
topdown.prob=TopDownStep(go.for.level,go.leaf.nodes,nodes.to.index,nodes.to.children,downtop.prob)
violate.detect.down=ViolateDetectprob(go.for.level, go.leaf.nodes,nodes.to.index,nodes.to.children,downtop.prob)
violate.detect.top=ViolateDetectprob(go.for.level, go.leaf.nodes,nodes.to.index,nodes.to.children,topdown.prob)


# valid.en=FALSE #用于对验证集进行结果评价
# if(valid.en==TRUE)
# {
#   predict.labels=matrix(0,nrow(valid.select.table),ncol(valid.select.table))
#   predict.scores=matrix(0,nrow(valid.select.table),ncol(valid.select.table))
#   rownames(predict.labels)=rownames(valid.select.table)
#   colnames(predict.labels)=colnames(valid.select.table)
#   rownames(predict.scores)=rownames(valid.select.table)
#   colnames(predict.scores)=colnames(valid.select.table)
#   for(i in 1:nrow(topdown.prob))
#   {
#     for(j in 1:length(except.root.labels.3))
#     {
#       predict.scores[i,j]=topdown.prob[i,(2*j-1)]
#       if(topdown.prob[i,(2*j-1)]>0.5)
#       {
#         predict.labels[i,except.root.labels.3[j]]=1
#       }
#     }
#   }
#   
#   valid.label.index=nodes.to.index[colnames(valid.select.table)]
#   measure.result=MHevaluate(predict.labels,valid.select.table)
#   F.each.class=F.measure.single.over.classes(valid.select.table, predict.labels)
#   prauc_result=PRAUCCalculate(predict.scores,valid.select.table)
#   
# }

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
  for(i in 1:nrow(topdown.prob))
  {
    for(j in 1:length(except.root.labels))
    {
      predict.scores[i,j]=topdown.prob[i,(2*j-1)]
      if(topdown.prob[i,(2*j-1)]>0.5)
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

#使用newpathrule方法处理的结果
final.result=NewPathrule(prob.for.genes,except.root.labels,go.for.level,go.leaf.nodes,test.select.table)
final.predict.labels=final.result[[1]]
final.predict.scores=final.result[[2]]
prauc_result.final=PRAUCCalculate(final.predict.scores,test.select.table)
measure.result.final=MHevaluate(final.predict.labels,test.select.table)

# if(test.en==TRUE)#此部分代码用于实现权值TPR算法
# {
#   predict.w.labels=matrix(0,nrow(test.select.table),ncol(test.select.table))
#   predict.w.scores=matrix(0,nrow(test.select.table),ncol(test.select.table))
#   rownames(predict.w.labels)=rownames(test.select.table)
#   colnames(predict.w.labels)=colnames(test.select.table)
#   rownames(predict.w.scores)=rownames(test.select.table)
#   colnames(predict.w.scores)=colnames(test.select.table)
#   for(i in 1:nrow(topdown.w.prob))
#   {
#     for(j in 1:length(except.root.labels.3))
#     {
#       predict.w.scores[i,j]=topdown.w.prob[i,(2*j-1)]
#       if(topdown.w.prob[i,(2*j-1)]>0.5)
#       {
#         predict.w.labels[i,except.root.labels.3[j]]=1
#       }
#     }
#   }
#   
#   test.label.index.w=nodes.to.index[colnames(test.select.table)]
#   measure.result.w=MHevaluate(predict.w.labels,test.select.table)
#   F.each.class.w=F.measure.single.over.classes(test.select.table, predict.w.labels)
#   prauc_result.w=PRAUCCalculate(predict.w.scores,test.select.table)
#   
# }

####第四步 用于将输出结果进行存储

result.output.en=TRUE
result.savepath=paste(data.path,"//result",sep = "")
setwd(result.savepath)
today <-Sys.Date()
run.num=1
output.fname=paste(file.prefix,file.middle,file.type,"_result",".txt",sep = "")
if(result.output.en==TRUE)
{
  
  write.table(today,file=output.fname,sep = " , ",eol="\n",quote=FALSE,row.names = FALSE,col.names = FALSE,append = FALSE)
  write.table(measure.result,file=output.fname,sep = " , ",eol="\n",quote=FALSE,row.names = FALSE,col.names = TRUE,append = TRUE)
  #write.table(measure.result.w,file=output.fname,sep = " ,",eol="\n",quote=FALSE,row.names = FALSE,col.names = TRUE,append = TRUE)
  write.table(prauc_result,file=output.fname,sep = " , ",eol="\n",quote=FALSE,row.names = FALSE,col.names = TRUE,append = TRUE)
  #write.table(prauc_result.w,file=output.fname,sep = " , ",eol="\n",quote=FALSE,row.names = FALSE,col.names = TRUE,append = TRUE)
  
}






