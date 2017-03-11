#生成每个节点的祖先节点列表，各节点的祖先节点列表均用序号表示
nodes.to.ancestors=list()
for(i in 1:length(except.root.labels.3))
{
  inter.ids=AnnotationDbi::get(except.root.labels.3[i], GOBPANCESTOR)
  ancestor.ids=intersect(except.root.labels.3,inter.ids)
  if(length(ancestor.ids)>0)
  {
    inter.ancestor.vec=vector()
    for(j in 1:length(ancestor.ids))
    {
      inter.ancestor.vec=c(inter.ancestor.vec,nodes.to.index[[ancestor.ids[j]]])
    }
    nodes.to.ancestors[[i]]=inter.ancestor.vec
  }    else
  {
    nodes.to.ancestors[[i]]=NA
  }
}  
names(nodes.to.ancestors)=except.root.labels.3#用节点名称给list中各元素命名

#生成每个节点的父节点列表，各节点的父节点列表均用序号表示
nodes.to.parents=list()
for(i in 1:length(except.root.labels.3))
{
  inter.ids=AnnotationDbi::get(except.root.labels.3[i], GOBPPARENTS)
  parent.ids=intersect(except.root.labels.3,inter.ids)
  if(length(parent.ids)>0)
  {
    inter.parent.vec=vector()
    for(j in 1:length(parent.ids))
    {
      inter.parent.vec=c(inter.parent.vec,nodes.to.index[[parent.ids[j]]])
    }
    nodes.to.parents[[i]]=inter.parent.vec
  }    else
  {
    nodes.to.parents[[i]]=NA
  }
}  
names(nodes.to.parents)=except.root.labels.3

#生成每个节点的子孙节点列表，各节点的父节点列表均用序号表示
nodes.to.descendants=list()
for(i in 1:length(except.root.labels.3))
{
  inter.ids=AnnotationDbi::get(except.root.labels.3[i], GOBPOFFSPRING)
  descendant.ids=intersect(except.root.labels.3,inter.ids)
  if(length(descendant.ids)>0)
  {
    inter.descendant.vec=vector()
    for(j in 1:length(descendant.ids))
    {
      inter.descendant.vec=c(inter.descendant.vec,nodes.to.index[[descendant.ids[j]]])
    }
    nodes.to.descendants[[i]]=inter.descendant.vec
  }    else
  {
    nodes.to.descendants[[i]]=NA
  }
}  
names(nodes.to.descendants)=except.root.labels.3
# nodes.to.parents=list()
# for(i in 1:length(except.root.labels.3))
# {
#   inter.ids=AnnotationDbi::get(except.root.labels.3[i], GOBPPARENTS)
#   nodes.to.parents[[i]]=intersect(except.root.labels.3,inter.ids)
#   
# }  


total.levels=length(go.for.level.3)#节点形成的层数

first.predict.labels=matrix(0,nrow(test.select.table),ncol(test.select.table))
for(i in 1:nrow(prob.for.genes))
{
  for(j in 1:length(except.root.labels.3))
  {
    
    if(prob.for.genes[i,(2*j-1)]>=0.5)
    {
      first.predict.labels[i,j]=1
    }
  }
}

first.measure.result=MHevaluate(first.predict.labels,test.select.table)
violate.detect.first=ViolateDetectlabel(go.for.level.3, go.leaf.nodes.3,nodes.to.index,nodes.to.children,first.predict.labels)

second.predict.labels=first.predict.labels#预测的最终结果将存于final.predict.labels中

for(k in 1:nrow(prob.for.genes))#对于每一个样本k
{
  for (i in (total.levels):2)#按层自下而上遍历节点,并且所选层次不为最顶层时
  {      
    
      for(j in 1:length(go.for.level.3[[i]]))
      {
          gene.name=(go.for.level.3[[i]])[[j]]#得到(go.for.level[[i]])[[j]]节点的GO标签
        
          gene.index=nodes.to.index[[(go.for.level.3[[i]])[j]]]#得到该节点的索引号
          if(first.predict.labels[k,gene.index]==1)#如果该节点的预测值为正值
          {
            parent.index=nodes.to.parents[[gene.index]]#得到其父节点的索引号
            if(length(parent.ids)>0)
            {
              parent.labels=first.predict.labels[k,parent.index]#得到其父节点的预测值，此时父节点可能有多个
              
              if(0 %in% parent.labels)#如果有一个父节点预测值为0，即与该节点的预测值相冲突
              {
                ancestor.index=nodes.to.ancestors[[gene.index]]#得到其祖先节点的索引号
                ancestor.labels=first.predict.labels[k,ancestor.index]#得到其祖先节点的预测值
                
                pos.nums=sum(ancestor.labels)#求取祖先节点中预测为1的个数
                neg.nums=length(ancestor.index)-pos.nums#求取祖先节点中预测为0的个数
                #如果预测为0个节点个数少有预测为1的节点个数
                if(neg.nums<(pos.nums+1))#+1表示将该节点的预测结果计算在内
                {
                  second.predict.labels[k,ancestor.index]=1#则其所有祖先节点均预测为1
                } else
                {
                  second.predict.labels[k,gene.index]=0
                }
              }
            }
          }
          
      }
      
  }
}
second.measure.result=MHevaluate(final.predict.labels,test.select.table)
final.predict.labels=second.predict.labels

for(k in 1:nrow(final.predict.labels))#遍历每一个样本
{
  for (i in 1:(length(go.for.level.3)-1))#按层自上而下遍历所有节点
  {
    for(j in 1:length(go.for.level.3[[i]]))#分别遍历每层的节点
    {
      is.leafnode=go.for.level.3[[i]][[j]] %in% go.leaf.nodes.3#判断该节点是否则为叶子节点
      gene.name=go.for.level.3[[i]][[j]]#得到此时处理节点的GO标签名称
      gene.index=nodes.to.index[[(go.for.level.3[[i]])[j]]]#得到该节点在节点列表中的序号
      if(is.leafnode!=TRUE)
      {
        if(final.predict.labels[k,gene.index]==0)#如果该节点的预测值为负值
        {
          descendant.index=nodes.to.descendants[[gene.index]]#得到其子孙节点的索引
          ancestor.index=nodes.to.ancestors[[gene.index]]#得到其祖先节点的索引号号
          children.index=nodes.to.children[[gene.index]]#得到此节点的所有子节点的序号
          children.labels=final.predict.labels[k,children.index]#得到其子节点的预测值，此时子节点可能有多个
          if(1 %in% children.labels)#如果子节点有预测值为正值的情况
          {
            if(i==1)
            {
              pos.nums=sum(final.predict.labels[k,children.index])+1
              neg.nums=0
              if(pos.nums>(neg.nums+1))
              {
                final.predict.labels[k,gene.index]=1
              } else
              {
                final.predict.labels[k,descendant.index]=0
              }
            }else
            {
              pos.nums=sum(final.predict.labels[k,children.index])+sum(final.predict.labels[k,ancestor.index])
              neg.nums=length(ancestor.index)-sum(final.predict.labels[k,ancestor.index])
              if(pos.nums>(neg.nums+1))
              {
                final.predict.labels[k,gene.index]=1
              } else
              {
                final.predict.labels[k,descendant.index]=0
              }
            }
          }
        }
      }
    }
    
  }
}



for(k in 1:nrow(final.predict.labels))#遍历每一个样本
{
  for (i in 1:(length(go.for.level.3)-1))#按层自上而下遍历所有节点
  {
    for(j in 1:length(go.for.level.3[[i]]))#分别遍历每层的节点
    {
      is.leafnode=go.for.level.3[[i]][[j]] %in% go.leaf.nodes.3#判断该节点是否则为叶子节点
      gene.name=go.for.level.3[[i]][[j]]#得到此时处理节点的GO标签名称
      gene.index=nodes.to.index[[(go.for.level.3[[i]])[j]]]#得到该节点在节点列表中的序号
      if(is.leafnode!=TRUE)
      {
        if(final.predict.labels[k,gene.index]==0)#如果该节点的预测值为负值
        {
          descendant.index=nodes.to.descendants[[gene.index]]#得到其子孙节点的索引
          ancestor.index=nodes.to.ancestors[[gene.index]]#得到其祖先节点的索引号号
          children.index=nodes.to.children[[gene.index]]#得到此节点的所有子节点的序号
          children.labels=final.predict.labels[k,children.index]#得到其子节点的预测值，此时子节点可能有多个
          if(1 %in% children.labels)#如果子节点有预测值为正值的情况
          {
            #final.predict.labels[k,gene.index]=1
            final.predict.labels[k,descendant.index]=0
          }
        }
      }
    }
  }
}
final.measure.result=MHevaluate(final.predict.labels,test.select.table)
violate.detect.final=ViolateDetectlabel(go.for.level.3, go.leaf.nodes.3,nodes.to.index,nodes.to.children,final.predict.labels)

# violate.list=list()#输出的结果
# sample.index=0#用于计算存在冲突结果的样本数量
# sample.name=c()#用于存储存在冲突结果的样本名
# 
# for(k in 1:nrow(final.predict.labels))#遍历每一个样本
# {
#   
#   inter.list=list()#用于存储单一样本的冲突结果
#   temp.index=0#用于计算单一样本中有多少个冲突结果
#   inter.list.name=c()#用于存储单一样本中哪些节点存在冲突的子节点
#   for (i in 1:(length(go.for.level.3)-1))#按层自上而下遍历所有节点
#   {
#     for(j in 1:length(go.for.level.3[[i]]))#分别遍历每层的节点
#     {
#       is.leafnode=go.for.level.3[[i]][[j]] %in% go.leaf.nodes.3#判断该节点是否则为叶子节点
#       gene.name=go.for.level.3[[i]][[j]]#得到此时处理节点的GO标签名称
#       gene.index=nodes.to.index[[(go.for.level.3[[i]])[j]]]#得到该节点在节点列表中的序号
#       
#       temp.result=c()#用于计算针对某个节点，哪几个子节点与其结果相冲突
#       
#       if(is.leafnode!=TRUE)
#       {
#         children.index=nodes.to.children[[gene.name]]#得到此节点的所有子节点的序号
#         for(m in 1:length(children.index))#遍历所有子节点
#         {
#           if(final.predict.labels[k,gene.index]==0)#当对于父节点，此样本为负时
#           {
#             if(final.predict.labels[k,children.index[m]]==1)#若此时子节点children.index[m]为正
#             {
#               temp.result=c(temp.result,children.index[m])#出现结果冲突，记录该子节点的序号
#             }
#           } 
#         }
#         if(!is.null(temp.result))#若确实存在冲突的父子节点
#         {
#           temp.index=temp.index+1#有冲突的节点数+1
#           inter.list[[temp.index]]=temp.result#将与此节点冲突的所有子节点存入
#           inter.list.name=c(inter.list.name,gene.index)#将此节点的序号存入
#         }
#       }
#     }
#     
#   }
#   if(length(inter.list)>0)
#   {
#     names(inter.list)=inter.list.name
#     sample.index=sample.index+1
#     violate.list[[sample.index]]=inter.list#将此样本的所有节点存入
#     sample.name=c(sample.name,k)
#   }
# }
# 
# names(violate.list)=sample.name#按照样本序号给结果命名
