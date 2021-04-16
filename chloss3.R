####实现chloss的程序文件 其中pi直接取为svm的输出结果
####sigma计算过程为论文中修改过的公式
####20190424

options(digits=21)
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

w1=1
w3=1
w2=2
w4=2
nodes.total.num=length(select.node)

c.matrix=matrix(0,nodes.total.num,1)
c.matrix[1,1]=1

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


#计算pi，直接使用概率
p.matrix=cbind(matrix(1,nrow(prob.for.genes),1),prob.is.one)


#计算pi，使用条件概率
# p.matrix=matrix(0,nrow(prob.for.genes),nodes.total.num)
# 
# for(k in 1:nrow(prob.for.genes))
# {
#   #根结点的pi为1
#   p.matrix[k,1]=1
#   for(i in 1:(nodes.total.num-1))
#   {
#     parents.index=nodes.to.parents.ch[[i+1]]
#     parents.num=length(nodes.to.parents.ch[[i+1]])
#     par.prob=1
#     for(j in 1:parents.num)
#     {
#       par.prob=par.prob*p.matrix[k,(parents.index[j]+1)]
#     }
#     p.matrix[k,(i+1)]=prob.is.one[k,i]*par.prob
#   }
# }
#定义sigma矩阵
sigma.one=matrix(0,nrow(prob.for.genes),nodes.total.num)
sigma.two=matrix(0,nrow(prob.for.genes),nodes.total.num)
#求解sigma(1)
for(k in 1:nrow(prob.for.genes))
{
  for(i in 1:nodes.total.num)
  {
    cur.node=select.node[i]
    is.leafnode=cur.node %in% go.leaf.nodes
    if(is.leafnode==TRUE)
    {
      sigma.one[k,i]=0
    } else
    {
      
      children.index=nodes.to.children.ch[[i]]
      children.num=length(children.index)
      sigma1.inter.1=1
      sigma1.inter.2=0
      for(j in 1:children.num)
      {
        sigma1.inter.2=sigma1.inter.2+ c.matrix[(children.index[j]+1),1]*p.matrix[k,(children.index[j]+1)]
        sigma1.inter.1=sigma1.inter.1* c.matrix[(children.index[j]+1),1]*p.matrix[k,(children.index[j]+1)]
      }
      
      sigma.one[k,i]=w2*sigma1.inter.2-w1*sigma1.inter.1
    }
  }
}

#求解sigma(2)
for(k in 1:nrow(prob.for.genes))
{
  for(i in 1:(nodes.total.num))
  {
    cur.node.index=i-1
    if(nodes.to.index.ch[i]==0)
    {
      sigma.two[k,i]=0
      
    }else
    {
      parents.index=nodes.to.parents.ch[[i]]
      parents.num=length(parents.index)
      sigma2.inter.1=1
      sigma2.inter.2=0
      for(j in 1:parents.num)
      {
        sigma2.inter.1=sigma2.inter.1*p.matrix[k,(parents.index[j]+1)]
        sigma2.inter.2=sigma2.inter.2+(1-p.matrix[k,(parents.index[j]+1)])
      }
      value.one=w1*c.matrix[i,1]*p.matrix[k,i]
      value.two=w3*c.matrix[i,1]*(sigma2.inter.1-p.matrix[k,i])
      value.three=w4*c.matrix[i,1]*sigma2.inter.2
      sigma.two[k,i]=value.one-value.two-value.three
    }
  }
}

sigma=sigma.one+sigma.two
