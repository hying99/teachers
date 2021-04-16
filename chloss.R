####实现chloss的程序文件 其中pi按文章中给出的算法实现
####20181218
options(digits=21)
w1=1
w3=1
w2=2
w4=2

#输入带有根结点的节点列表，根结点的序号为0
total.index.ch=MakeIndex(select.node,include.root=TRUE)
nodes.to.index.ch=total.index.ch[[1]]
nodes.to.children.ch=total.index.ch[[2]]
nodes.to.ancestors.ch=total.index.ch[[3]]
nodes.to.parents.ch=total.index.ch[[4]]
nodes.to.descendants.ch=total.index.ch[[5]]


#### 将SVM的概率结果读入
file.middle="0"
file.type="change"
#file.type=""

#设置mat文件存储路径
setwd(paste(data.path,"//matfile",sep = ""))
mat.file=paste(file.prefix,file.middle,file.type,"_decision.mat",sep = "")
probability.data=readMat(mat.file,fixNames = FALSE)
#prob.for.genes  存储时，对每个节点来说，第一位是为label1的概率，而后是为0的概率
prob.for.genes=probability.data$decision_test


nodes.total.num=length(select.node)

c.matrix=matrix(0,nodes.total.num,1)
c.matrix[1,1]=1
c.matrix=matrix(1,nodes.total.num,1)
#计算Ci,该值仅与节点有关，对每个样本均相同
for(i in 1:(nodes.total.num-1))
{
  if(each.go.level.num[i]==1)
  {
    sibling.num=each.level.nodes.num[[1]]
    c.matrix[(i+1),1]=1/sibling.num
  }else
  {
    
    parentnode.num=length(nodes.to.parents[[i]])
    current.c=0
    for(j in 1:parentnode.num)
    {
      parentnode.index=nodes.to.parents[[i]][j]
      sibling.num=length(nodes.to.children[[parentnode.index]])
      current.c=current.c+c.matrix[parentnode.index,1]/sibling.num
    }
    c.matrix[(i+1),1]=current.c/parentnode.num
  }
}
c.matrix=matrix(1,nodes.total.num,1)
#提取各节点为1的概率矩阵
prob.is.one=matrix(0,nrow(prob.for.genes),(nodes.total.num-1))
for(i in 1:nrow(prob.for.genes))
{
  for(j in 1:(nodes.total.num-1))
  {
    prob.is.one[i,j]=prob.for.genes[i,(2*j-1)]
  }
}
#计算pi
p.matrix=matrix(0,nrow(prob.for.genes),nodes.total.num)

for(k in 1:nrow(prob.for.genes))
{
  #根结点的pi为1
  p.matrix[k,1]=1
  for(i in 1:(nodes.total.num-1))
  {
    if(each.go.level.num[i]==1)
    {
      p.matrix[k,(i+1)]=prob.is.one[k,i]
    } else
    {
      parents.index=nodes.to.parents[[i]]
      parents.num=length(nodes.to.parents[[i]])
      
      par.prob=1
      for(j in 1:parents.num)
      {
        par.prob=par.prob*p.matrix[k,(parents.index[j]+1)]
      }
      p.matrix[k,(i+1)]=prob.is.one[k,i]*par.prob
    }
  }
}
# for(k in 1:nrow(prob.for.genes))
# {
#   #根结点的pi为1
#   p.matrix[k,1]=1
#   for(i in 1:(nodes.total.num-1))
#   {
#     if(each.go.level.num[i]==1)
#     {
#       p.matrix[k,(i+1)]=prob.is.one[k,i]
#     } else
#     {
#       #当前节点的祖先样本数量
#       ancestor.index=nodes.to.ancestors[[i]]
#       #当前节点的祖先样本个数
#       ancestor.num=length(ancestor.index)
#       anc.prob=1
#       for(j in 1:ancestor.num)
#       {
#         anc.prob=anc.prob*prob.is.one[k,ancestor.index[j]]
#       }
#       p.matrix[k,(i+1)]=prob.is.one[k,i]*anc.prob
#     }
#   }
# }

#sigma=matrix(0,nrow(prob.for.genes),nodes.total.num)
sigma.one=matrix(0,nrow(prob.for.genes),nodes.total.num)
sigma.two=matrix(0,nrow(prob.for.genes),nodes.total.num)
#求解sigma(1)
for(k in 1:nrow(prob.for.genes))
{
  for(i in 1:(nodes.total.num))
  {
    cur.node.index=i-1
    if(cur.node.index==0)
    {
      children.index=c(1:each.level.nodes.num[[1]])
      children.num=length(children.index)
      inter.value=0
      for(j in 1:children.num)
      {
        inter.value=inter.value+ c.matrix[(children.index[j]+1),1]*p.matrix[k,(children.index[j]+1)]
      }
      sigma.one[k,i]=(w2-w1)*inter.value
      
    }else
    {
      cur.node=except.root.labels[cur.node.index]
      #判断是否为叶子节点
      is.leafnode=cur.node %in% go.leaf.nodes
      if(is.leafnode==TRUE)
      {
        sigma.one[k,i]=0
      } else
      {
        children.index=nodes.to.children[[cur.node.index]]
        children.num=length(children.index)
        inter.value=0
        for(j in 1:children.num)
        {
          inter.value=inter.value+ c.matrix[(children.index[j]+1),1]*p.matrix[k,(children.index[j]+1)]
        }
        sigma.one[k,i]=(w2-w1)*inter.value
      }
    }
  }
}
#求解sigma(2)
for(k in 1:nrow(prob.for.genes))
{
  for(i in 1:(nodes.total.num))
  {
    cur.node.index=i-1
    if(cur.node.index==0)
    {
      sigma.two[k,i]=0
      
    }else
    {
      if(each.go.level.num[cur.node.index]==1)
      {
        sigma.two[k,i]=w1*c.matrix[i,1]*p.matrix[k,i]-w3*c.matrix[i,1]*(1-p.matrix[k,i])
      } else
      {
        parents.index=nodes.to.parents[[cur.node.index]]
        parents.num=length(parents.index)
        inter.value=1
        for(j in 1:parents.num)
        {
          inter.value=inter.value*p.matrix[k,(parents.index[j]+1)]
        }
        value.one=w1*c.matrix[i,1]*p.matrix[k,i]*parents.num
        value.two=w3*c.matrix[i,1]*parents.num*(inter.value-p.matrix[k,i])
        value.three=w4*c.matrix[i,1]*parents.num*(1-inter.value)
        sigma.two[k,i]=value.one-value.two-value.three
      }
    }
  }
}

sigma=sigma.one+sigma.two

