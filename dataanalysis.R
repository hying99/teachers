##用于统计各数据集样本及标签的基本信息
##使用时需承接rebuilddataprocess.R文件
##已将数据集基本信息写入图片中
##该文件 20180423创建

#用于选择进行分析的是1 训练集、2 验证集还是 3 测试集
dataset.select=2
if(dataset.select==1)
{
  #用于分析训练集数据
  input.data=remain.select.data
  except.root.label.matrix=except.root.table
  except.root.label.set=except.root.labels
  root.label.set=select.node
}else if(dataset.select==2)
{
  #用于分析验证集数据
  input.data=valid.select.data
  except.root.label.matrix=valid.select.table
  except.root.label.set=except.root.labels
  root.label.set=select.node
}else if(dataset.select==3)
{
  #用于分析测试集数据
  input.data=test.select.data
  except.root.label.matrix=test.select.table
  except.root.label.set=except.root.labels
  root.label.set=select.node
}

#####数据集样本及标签基本信息统计

#训练集样本数量
sample.num=nrow(input.data)
#平均每个样本所含有的标签数量
average.label=(sum(except.root.label.matrix)+sample.num)/sample.num
total.label.num=length(except.root.label.set)
each.label.sample.num=matrix(0,1,total.label.num)
colnames(each.label.sample.num)=except.root.label.set

for (i in 1:total.label.num)
{
  each.label.sample.num[1,i]=sum(except.root.label.matrix[,i])
  
}

#绘制带有标签和样本信息的DAG图
num.only=TRUE
plot.en=TRUE
graph.select.node <- subGraph(root.label.set, univ.graph)
if(num.only==TRUE)
{
  graph.new.nodes=graph.select.node@nodes
  graph.new.edgeL=graph.select.node@edgeL
  for(i in 1:length(graph.new.nodes))
  {
    if(graph.new.nodes[i]=="GO:0008150")
    {
      graph.new.nodes[i]=paste(0,sample.num,sep = ":")
    }    else
    {
      graph.new.nodes[i]=paste(which(except.root.label.set==graph.new.nodes[i]),each.label.sample.num[1,which(except.root.label.set==graph.new.nodes[i])],sep = ":")
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
      graph.new.nodes[i]=paste(graph.new.nodes[i],sample.num,sep = ":")
    }    else
    {
      graph.new.nodes[i]=paste(graph.new.nodes[i],each.label.sample.num[1,which(except.root.label.set==graph.new.nodes[i])],sep = ":")
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
  
  legend(13000,10000,legend="sample number" , bty ="n", pch=NA)
  legend(17000,10000,legend=sample.num, bty ="n", pch=NA)
  legend(13000,9500,legend="label number" , bty ="n", pch=NA)
  legend(17000,9500,legend=total.label.num, bty ="n", pch=NA) 
  legend(13000,9000,legend="feature number" , bty ="n", pch=NA)
  legend(17000,9000,legend=ncol(input.data), bty ="n", pch=NA) 
  legend(13000,8500,legend="level number" , bty ="n", pch=NA)
  legend(17000,8500,legend=total.levels, bty ="n", pch=NA) 
  legend(13000,8000,legend="average labels per sample" , bty ="n", pch=NA)
  legend(17000,8000,legend=average.label, bty ="n", pch=NA)
  legend(13000,7000,legend="level index" , bty ="n", pch=NA)
  legend(17000,7000,legend="each level label number" , bty ="n", pch=NA)
  for(i in 1:total.levels)
  {
    legend(13000,(7000-500*i),legend=i , bty ="n", pch=NA)
    legend(17000,(7000-500*i),legend=each.level.nodes.num[[i]] , bty ="n", pch=NA)
  }
}

