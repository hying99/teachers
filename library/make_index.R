#用于生成个GO标签的索引序列号，以及每个节点与其子节点的序列号对应关系
MakeIndex<-function (except.root.labels)
{
  nodes.to.index=list()
  for(i in 1:length(except.root.labels))
  {
    nodes.to.index[i]=i
  }
  
  names(nodes.to.index)=except.root.labels
  no.root.nodes.num=length(except.root.labels)
  #names(except.root.labels)=c(1:no.root.nodes.num)
  
  nodes.to.children=list()
  for(i in 1:length(except.root.labels))
  {
    inter.ids=AnnotationDbi::get(except.root.labels[[i]], GOBPCHILDREN)
    child.ids=intersect(except.root.labels,inter.ids)
    if(length(child.ids)>0)
    {
      inter.child.vec=vector()
      for(j in 1:length(child.ids))
      {
        inter.child.vec=c(inter.child.vec,nodes.to.index[[child.ids[j]]])
      }
      nodes.to.children[[i]]=inter.child.vec
    }
    else
    {
      nodes.to.children[[i]]=NA
    }
  }
  names(nodes.to.children)=except.root.labels
  total.index=list(nodes.to.index,nodes.to.children)
  return (total.index)
  
}




