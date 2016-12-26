PlotLabelGraph<-function (except.root.nodes,BP.univ.graph,plot.en=TRUE,num.only=FALSE,output.en=FALSE,write.pic.name)
{
select.nodes=c(except.root.nodes,"GO:0008150")
graph.select.node <- subGraph(select.nodes, BP.univ.graph)
if(num.only==TRUE)
{
  graph.new.nodes=graph.select.node@nodes
  graph.new.edgeL=graph.select.node@edgeL
  for(i in 1:length(graph.new.nodes))
  {
    if(graph.new.nodes[i]=="GO:0008150")
    {
      graph.new.nodes[i]=0
    }    else
    {
      graph.new.nodes[i]=which(except.root.nodes==graph.new.nodes[i])
    }
  }
  
  names(graph.new.edgeL)=graph.new.nodes
} else
{
  graph.new.nodes=graph.select.node@nodes
  graph.new.edgeL=graph.select.node@edgeL
  for(i in 1:length(graph.new.nodes))
  {
    if(graph.new.nodes[i]=="GO:0008150")
    {
      graph.new.nodes[i]=paste(0,graph.new.nodes[i],sep = ":")
    }    else
    {
      graph.new.nodes[i]=paste(which(except.root.nodes==graph.new.nodes[i]),graph.new.nodes[i],sep = ":")
    }
  }
  names(graph.new.edgeL)=graph.new.nodes
  
}
graph.new=graphNEL(graph.new.nodes, graph.new.edgeL, edgemode = "directed")
if(plot.en==TRUE)
{
  x11()
  if(num.only==TRUE)
  {
    Pretty.plot.graph(graph.new,fontsize=15,fillcolor="transparent",height=6,
                      width=8,color="transparent", fontcolor="black")
  }  else
  {
    Pretty.plot.graph(graph.new,fontsize=11,fillcolor="transparent",height=6,
                      width=8,color="transparent", fontcolor="black")
  }
}

if(output.en==TRUE)
{
  postscript(write.pic.name,onefile=FALSE);
  Pretty.plot.graph(graph.new,fontsize=11,fillcolor="transparent",height=6,width=8,color="transparent", fontcolor="black");
  dev.off();
}

}
