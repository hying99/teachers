#用于实现后期概率值的处理 同时也给出概率转化为标签的评价指标 20170315

####第一步 将SVM的概率结果读入

file.type="change"
file.middle="0"
#设置mat文件存储路径
setwd(paste(data.path,"//matfile",sep = ""))
mat.file=paste(file.prefix,file.middle,file.type,"_decision.mat",sep = "")
probability.data=readMat(mat.file,fixNames = FALSE)
#probability.data=readMat("cellcycle0replaceprocess_decision.mat",fixNames = FALSE)
#prob.for.genes=probability.data$decision
prob.for.genes=probability.data$decision_test



first.prob=prob.for.genes
second.prob=first.prob
####第二步 将概率值转换为标签 1为正 0为负
first.predict.labels=matrix(0,nrow(test.select.table),ncol(test.select.table))
second.predict.labels=matrix(0,nrow(test.select.table),ncol(test.select.table))
final.predict.labels=matrix(0,nrow(test.select.table),ncol(test.select.table))
#存储各节点预测为1的概率值
final.predict.scores=matrix(0,nrow(test.select.table),ncol(test.select.table))
for(i in 1:nrow(first.prob))
{
  for(j in 1:length(except.root.labels.3))
  {
    #将初步预测的概率值转化为0 1标签
    if(first.prob[i,(2*j-1)]>=0.5)
    {
      first.predict.labels[i,j]=1
    }
  }
}
#计算初始标签值的评价指标
first.measure.result=MHevaluate(first.predict.labels,test.select.table)
#violate.detect.first=ViolateDetectlabel(go.for.level.3, go.leaf.nodes.3,nodes.to.index,nodes.to.children,first.predict.labels)

second.predict.labels=first.predict.labels
for(k in 1:nrow(prob.for.genes))#对于每一个样本k
{
  for (i in 1:total.levels)#按层自上而下遍历节点
  {
    for(j in 1:length(go.for.level.3[[i]]))
    {
      gene.name=(go.for.level.3[[i]])[[j]]#得到(go.for.level[[i]])[[j]]节点的GO标签
      gene.index=nodes.to.index[[(go.for.level.3[[i]])[j]]]#得到该节点的索引号
      ancestor.index=nodes.to.ancestors[[gene.index]]#得到其祖先节点的索引号
      ancestor.labels=first.predict.labels[k,ancestor.index]#得到其祖先节点的预测值
      children.index=nodes.to.children[[gene.index]]#得到此节点的所有子节点的序号
      children.labels=first.predict.labels[k,children.index]#得到其子节点的预测值，此时子节点可能有多个
      parent.index=nodes.to.parents[[gene.index]]#得到其父节点的索引号
      parent.labels=first.predict.labels[k,parent.index]#得到其父节点的预测值，此时父节点可能有多个
      is.leafnode=go.for.level.3[[i]][[j]] %in% go.leaf.nodes.3#判断该节点是否则为叶子节点
      
      if(i==1)#对于第一层节点，此层节点假定均属于根结点
      {
        if(first.predict.labels[k,gene.index]==0)#如果该节点的预测值为负值
        {
          if(is.leafnode!=TRUE)#如果该节点不是叶子节点
          {
            if(1 %in% children.labels)#如果该节点的子节点有为1的情况
            {
              pos.nums=1+sum(children.labels)#求取父节点中预测为1的个数，及子节点中为1的个数之和
              neg.nums=1
              #可以看出，当节点属于第一层且不为叶子节点时，
              #对预测值0是否更改仅取决于其子节点是否有为1的情况
              if(pos.nums>neg.nums)
              {
                child.pos.prob=c()
                for(m in 1:length(children.index))
                {
                  if(first.prob[k,(2*children.index[m]-1)]>=0.5)#若此子节点的概率大于等于0.5
                  {
                    #将此子节点的概率存储在child.prob中
                    child.pos.prob=c(child.pos.prob,first.prob[k,(2*(children.index[m])-1)])
                  }
                }
                second.predict.labels[k,gene.index]=1
                #将子节点中的最大概率赋值给该节点当做初始概率
                second.prob[k,(2*gene.index-1)]=max(child.pos.prob)
                second.prob[k,(2*gene.index)]=1-max(child.pos.prob)
              } 
            }
          }
        }
      } else#对于属于其它层的节点
      {
        if(first.predict.labels[k,gene.index]==0)#如果该节点的预测值为负值
        {
          if(is.leafnode!=TRUE)#如果这个节点不是叶子节点
          {
            if(1 %in% children.labels)#并且这个节点存在预测值为1的子节点
            {
              pos.nums=sum(parent.labels)+sum(children.labels)#求取父节点中预测为1的个数，及子节点中为1的个数
              neg.nums=length(parent.index)-sum(parent.labels)#求取父节点中预测为0的个数
              
              if(pos.nums>neg.nums)
              {
                child.pos.prob=c()
                for(m in 1:length(children.index))
                {
                  if(first.prob[k,(2*children.index[m]-1)]>=0.5)#若此子节点的概率大于等于0.5
                  {
                    #将此子节点的概率存储在child.prob中
                    child.pos.prob=c(child.pos.prob,first.prob[k,(2*(children.index[m])-1)])
                  }
                }
                second.predict.labels[k,gene.index]=1
                second.prob[k,(2*gene.index-1)]=max(child.pos.prob)
                second.prob[k,(2*gene.index)]=1-max(child.pos.prob)
              } 
            } 
          }
        } else #如果该节点的预测值为正值
        {
          if(is.leafnode!=TRUE)#如果这个节点不是叶子节点
          {
            if(0 %in% parent.labels)#并且这个节点存在预测值为0的父节点
            {
              pos.nums=sum(parent.labels)+sum(children.labels)#求取父节点中预测为1的个数，及子节点中为1的个数
              neg.nums=length(parent.index)-sum(parent.labels)#求取父节点中预测为0的个数
              if(neg.nums>=pos.nums)
              {
                parent.pos.prob=c()
                parent.neg.prob=c()
                for(n in 1:length(parent.index))
                {
                  
                  inter.prob=prob.for.genes[k,(2*parent.index[n]-1)]
                  if(inter.prob>=0.5)#若此子节点的概率大于等于0.5
                  {
                    #将此子节点的概率存储在child.prob中
                    parent.pos.prob=c(parent.pos.prob,first.prob[k,(2*(parent.index[n])-1)])
                  } else
                  {
                    parent.neg.prob=c(parent.neg.prob,first.prob[k,(2*(parent.index[n]))])
                  }
                }
                second.predict.labels[k,gene.index]=0
                #将父节点的最大负概率赋值给该节点作为初始负概率
                second.prob[k,(2*gene.index-1)]=1-max(parent.neg.prob)
                second.prob[k,(2*gene.index)]=max(parent.neg.prob)
              } 
               
            } 
          }else#如果这个节点是叶子节点,叶子节点没有子节点
          {
            
            pos.nums=sum(parent.labels)#求取父节点中预测为1的个数
            neg.nums=length(parent.index)-sum(parent.labels)#求取父节点中预测为0的个数
            if(neg.nums>=pos.nums)
            {
              parent.pos.prob=c()
              parent.neg.prob=c()
              for(n in 1:length(parent.index))
              {
                inter.prob=prob.for.genes[k,(2*parent.index[n]-1)]
                if(inter.prob>=0.5)#若此子节点的概率大于等于0.5
                {
                  #将此子节点的概率存储在child.prob中
                  parent.pos.prob=c(parent.pos.prob,first.prob[k,(2*(parent.index[n])-1)])
                } else
                {
                  parent.neg.prob=c(parent.neg.prob,first.prob[k,(2*(parent.index[n]))])
                }
              }
              second.predict.labels[k,gene.index]=0
              second.prob[k,(2*gene.index-1)]=1-max(parent.neg.prob)
              second.prob[k,(2*gene.index)]=max(parent.neg.prob)
            } 
          }
        }
      }
    }
    
  }
}

#计算二次标签的评价指标
second.measure.result=MHevaluate(second.predict.labels,test.select.table)


##自顶而下遍历样本，使得底层节点为正的概率小于等于上层节点为正的概率############
final.prob=second.prob

# downtop.prob=DownTopStep(go.for.level.3,go.leaf.nodes.3,nodes.to.index,nodes.to.children,second.prob)
# final.prob=TopDownStep(go.for.level.3,go.leaf.nodes.3,nodes.to.index,nodes.to.children,downtop.prob)

for(k in 1:nrow(final.prob))
{
  for (i in 1:(length(go.for.level.3)-1))
  {
    for(j in 1:length(go.for.level.3[[i]]))
    {
      gene.name=go.for.level.3[[i]][[j]]
      is.leafnode=gene.name %in% go.leaf.nodes.3
      gene.index=nodes.to.index[[gene.name]]
      pos.index=2*gene.index-1
      neg.index=2*gene.index
      if(is.leafnode!=TRUE)
      {
        inter.neg.prob=final.prob[k,neg.index]
        inter.pos.prob=final.prob[k,pos.index]
        children.index=nodes.to.children[[gene.name]]
        for(m in 1:length(children.index))
        {
            #若此时子节点children.index[m]为正的概率比父节点为正的概率大
            if(final.prob[k,(2*(children.index[m])-1)]>inter.pos.prob)
            {
              final.prob[k,(2*(children.index[[m]])-1)]=inter.pos.prob
              final.prob[k,2*(children.index[[m]])]=1-inter.pos.prob
            }
          
        }
      }
    }
  }
}

#将计算的最终概率转化为标签值
for(i in 1:nrow(final.prob))
{
  for(j in 1:length(except.root.labels.3))
  {
    final.predict.scores[i,j]=final.prob[i,(2*j-1)]
    if(final.prob[i,(2*j-1)]>=0.5)
    {
      final.predict.labels[i,j]=1
    }
  }
}

#### 计算结果的评价指标


#F.each.class=F.measure.single.over.classes(test.select.table, final.predict.labels)
prauc_result.final=PRAUCCalculate(final.predict.scores,test.select.table)
measure.result.final=MHevaluate(final.predict.labels,test.select.table)

# compare.predict.labels=matrix(0,nrow(test.select.table),4)
# for(i in 1:nrow(test.select.table))
# {
#   compare.predict.labels[i,1]=sum(test.select.table[i,])
#   compare.predict.labels[i,2]=sum(first.predict.labels[i,])
#   compare.predict.labels[i,3]=sum(second.predict.labels[i,])
#   compare.predict.labels[i,4]=sum(final.predict.labels[i,])
#   
# }


















