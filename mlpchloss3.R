####采用neural netwrok作为基础分类器，运行前需先运行chddataprocess或其他文件，得到数据及类别信息
####后处理方法采用CHLOSS
####20190308
cat('The neural network starts' ,"\n")
label.length=length(except.root.labels)
#label.length=15
model.list=list()
nn.predict.labels=matrix(0,nrow(test.select.table),ncol(test.select.table))
nn.predict.scores=matrix(0,nrow(test.select.table),ncol(test.select.table))
rownames(nn.predict.labels)=rownames(test.select.table)
colnames(nn.predict.labels)=colnames(test.select.table)
rownames(nn.predict.scores)=rownames(test.select.table)
colnames(nn.predict.scores)=colnames(test.select.table)
nn.results.evaluation=matrix(0,ncol(test.select.table),6)
rownames(nn.results.evaluation)=colnames(test.select.table)
colnames(nn.results.evaluation)=c('Prec','Rec','Spe','F','Acc','Npos')
batch.size=32
epochs=100
patience=5
for (i in 1:label.length)
  #for (i in 15:15)
{
  gene.index=i#得到处理节点的序号
  cat('The processing class is ' ,gene.index,"\n")
  input.data =data.total[[gene.index]]$X#训练数据
  input.label=t(t(as.numeric(data.total[[gene.index]]$labels))) #训练数据的标签
  input.label[input.label==2]=0
  valid.input.data=valid.select.data
  valid.input.label=t(t(valid.select.table[,gene.index]))
  col.num=ncol(input.data)#训练数据属性的个数，即维度
  first.hneuron.num=round((col.num+1)/2)#隐藏层节点的数量
  second.hneuron.num=first.hneuron.num
  #second.hneuron.num=round((first.hneuron.num+1)/2)#隐藏层节点的数量
  model <- keras_model_sequential()
  
  # add layers and compile the model
  model %>%
    layer_dense(units = first.hneuron.num, activation = 'relu', input_shape = c(col.num)) %>%
    layer_dropout(rate = 0.5) %>%
    layer_dense(units = second.hneuron.num, activation = 'relu') %>%
    layer_dropout(rate = 0.5) %>%
    layer_dense(units = 1, activation = 'sigmoid')
  model %>% compile(
    #loss = loss_mean_squared_error,
    loss = loss_binary_crossentropy,
    optimizer = optimizer_adagrad(),
    #optimizer = optimizer_rmsprop(),
    metrics = c(metric_binary_accuracy)
  )
  
  
  history <- model %>% fit(
    input.data, input.label,
    epochs = epochs, batch_size = batch.size,
    #validation_split = 0.2
    validation_data = list(valid.input.data,valid.input.label),
    callbacks=list(callback_early_stopping(monitor="val_loss",min_delta=0,patience=patience,verbose=1,mode=c("auto")))
  )
  
  nn.predict.labels[,gene.index]=predict_classes(model,test.select.data,batch_size = batch.size)
  nn.predict.scores[,gene.index]=predict_proba(model,test.select.data,batch_size = batch.size)
  model.list[[i]]=model
  test.true.label=t(t(test.select.table[,i]))
  nn.results.evaluation[i,]=F.measure.single(nn.predict.labels[,i],test.true.label)
}
names(model.list)=except.root.labels

prob.for.genes=matrix(0,nrow(test.select.table),(ncol(test.select.table)*2))
for(i in 1:nrow(test.select.table))
{
  for(j in 1:label.length)
  {
    prob.for.genes[i,(2*j-1)]=nn.predict.scores[i,j]
    prob.for.genes[i,(2*j)]=1-nn.predict.scores[i,j]
  }
}


####实现chloss的程序文件 其中pi直接取为MLP的输出结果
####20190308

options(digits=21)
#输入带有根结点的节点列表，根结点的序号为0
total.index.ch=MakeIndex(select.node,include.root=TRUE)
nodes.to.index.ch=total.index.ch[[1]]
nodes.to.children.ch=total.index.ch[[2]]
nodes.to.ancestors.ch=total.index.ch[[3]]
nodes.to.parents.ch=total.index.ch[[4]]
nodes.to.descendants.ch=total.index.ch[[5]]

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


#运行DAGlable
setwd(work.path)
source("DAGlabel.R")
