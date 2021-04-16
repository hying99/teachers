options(digits=21)
pri.result=matrix(0,nrow(prob.for.genes),(nodes.total.num-1))
for(i in 1:nrow(prob.for.genes))
{
  
  for(j in 1:(nodes.total.num-1))
  {
    if(sigma[i,(j+1)]>0)
    {
      pri.result[i,j]=1
    } 
  }
}

#violate.list.pri=ViolateDetectlabel(go.for.level,go.leaf.nodes,nodes.to.index,nodes.to.children,pri.result)
final.result=pri.result
#final.result=matrix(0,nrow(prob.for.genes),(nodes.total.num-1))

for(k in 1:nrow(prob.for.genes))
{
  for (i in 1:(length(go.for.level)-1))
  {
    for(j in 1:length(go.for.level[[i]]))
    {
      gene.name=go.for.level[[i]][[j]]
      is.leafnode=gene.name %in% go.leaf.nodes
      gene.index=nodes.to.index[[gene.name]]
      if(is.leafnode!=TRUE)
      {
        children.index=nodes.to.children[[gene.name]]
        for(m in 1:length(children.index))
        {
          if(final.result[k,gene.index]==0)
          {
            final.result[k,children.index[m]]=0
          }
        }
      }
    }
  }
}


# for(i in 1:nrow(prob.for.genes))
# {
#   for (j in 1:(nodes.total.num-1))
#   {
#     if(pri.result[i,j]==0)
#     {
#       
#       cur.node=except.root.labels[j]
#       #判断是否为叶子节点
#       is.leafnode=cur.node %in% go.leaf.nodes
#       if(is.leafnode==FALSE)
#       {
#         descants.index=nodes.to.descendants[[j]]
#         descants.num=length(descants.index)
#         for(K in 1:children.num)
#         {
#           final.result[i,(descants.index[k])]=0
#         }
#       }
#     }
#   }
# }
#violate.list.final=ViolateDetectlabel(go.for.level,go.leaf.nodes,nodes.to.index,nodes.to.children,final.result)
measure.result.pri=MHevaluate(pri.result,test.select.table)
measure.result.fin=MHevaluate(final.result,test.select.table)
