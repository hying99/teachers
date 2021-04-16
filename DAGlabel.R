#DAGlabel的实现，在rebuilddataprocess.R，chloss2.R之后运行
#2019年1月9日最后一次修改
#主要内容有寻找sigma最大节点，寻找其待结合的父节点，节点合成和修改层级关系几部分。
work.path="D://R//CODES"
data.path="D://R//DATA"
setwd(work.path)
source("headfile.R")


y=matrix(data = 0, nrow = nrow(sigma), ncol = (ncol(sigma)+100)) #y矩阵 列有扩充，用于存放合成的节点
y[,ncol(sigma)]=1   #根节点为1

for(k in 1:nrow(sigma)) #K为样本数
{
  #每次循环时，初始化
  nodes_to_children = withroot.nodes.to.children  #子节点list
  nodes_to_parents= withroot.nodes.to.parents     #父节点list
  
  #每一个样本的sigma矩阵（无根节点）,寻找最大值时，会有改动
  si=matrix(data = 0, nrow =ncol(sigma)-1, ncol = 2)  
  si[,1]=1:(ncol(sigma)-1)        
  si[,2]=sigma[k,2:ncol(sigma)]   #提取sigma矩阵的第k行
  #每一个样本的sigma矩阵（有根节点，在最后），无改动
  total_si=matrix(data = 0, nrow =ncol(sigma), ncol = 2)
  total_si[,1]=1:ncol(sigma)
  total_si[ncol(sigma),1]=ncol(sigma)
  total_si[ncol(sigma),2]=50      #假设给根节点一个较大的sigma值，值无影响，较大即可
  total_si[(1:ncol(sigma)-1),]=si #提取sigma矩阵的第k行
  
  bignodes=matrix(data = 0, nrow =1, ncol = 2) #合成节点矩阵
  relation=matrix(data = 1, nrow =1, ncol = 2) #存放合成前的节点与合成后的节点
  rei=matrix(data = 1, nrow =2, ncol = 2) #存放合成前的节点与合成后的节点
  
  find.max=which(si[,2]==max(si[,2]))  #寻找sigma最大值的下角标
  find.max=sort(find.max,decreasing = F) #可能有多个值最大且相等，优先取下角标最小的
  find.max=find.max[1]  #sigma最大值的下角标
  
  p=1 #中间变量的初始化 用于合成节点的标号
  list_parentnodes=list() #用于存放待结合的父节点的下角标 
  
  while((si[find.max,2]>=0)&&(sum(y[k,1:ncol(sigma)])<ncol(sigma))) #sigma最大值小于0或者所有数据都已标记后，退出
  {
    num_parentnodes=length(nodes_to_parents[[si[find.max,1]]]) #sigma最大的节点的父节点个数
    labeled_count=0    #已标记父节点个数
    for(i in 1:num_parentnodes)
    {
      if(y[k,nodes_to_parents[[si[find.max,1]]][i]]==1)
      {
        labeled_count=labeled_count+1
      }
    }
    if(labeled_count== num_parentnodes)  #所有父节点已经标记
    {
      y[k,si[find.max,1]]=1   #该节点标记为1
      si=si[-find.max,]       #移出si矩阵
      
      
    }else   #找到sigma小于0的父节点，一个或一个以上
    {
      #计算父节点中sigma值小于0的个数
      num_combinenodes=0 #待合成父节点的个数
      list_parentnodes=list() #存放父节点的下角标
      for(m in 1:num_parentnodes)
      {
        row_parentnode=which(total_si[,1]==nodes_to_parents[[si[find.max,1]]][m])
        if(total_si[row_parentnode,2]<0)
        {
          num_combinenodes=num_combinenodes+1
          list_parentnodes=append(list_parentnodes,row_parentnode) #pa里存放的是待合成的父节点下角标，
          
        }
      }
      
      if(num_combinenodes<=1) #找最小sigma，要返回一个1
      {
        list_parentnodes=list() #存放父节点的下角标
        min_panodes_sigma=10  #数值是多少无意义，只要比真实sigma值大就可以
        #找sigma值最小的父节点
        for(i in 1:num_parentnodes) 
        {
          row_parentnode=which(total_si[,1]==nodes_to_parents[[si[find.max,1]]][i]) #父节点的下角标 nodes_to_parents已经更改
          sigama_of_parent=total_si[row_parentnode,2]  #父节点的sigma值
          
          if(min_panodes_sigma>=sigama_of_parent)
          {
            min_panodes_sigma = sigama_of_parent   #父节点的sigma的最小值
            node_parent=total_si[row_parentnode,1]    #sigma取最小值时的父节点
          }
        }
        
        row_parentnode=which(total_si[,1]==node_parent)    #父节点下角标
        list_parentnodes=append(list_parentnodes,row_parentnode)   #pa里存放的是待合成的父节点下角标
        num_combinenodes=1
      }
      
      
      #合成节点
      rei[1,1]=si[find.max,1]
      value_bignodes=si[find.max,2] #子节点的sigma值
      m=1 #循环变量
      while(m<=num_combinenodes)
      {
        #构建合成节点
        value_bignodes=(value_bignodes*m+total_si[list_parentnodes[[m]],2])/(m+1)
        bignodes[1,1]=ncol(sigma)+p   #合成节点
        bignodes[1,2]=value_bignodes
        
        si=rbind(si,bignodes)   #合成节点拼接到si后面
        total_si=rbind( total_si,bignodes) 
        # 小节点与合成节点对应关系
        
        rei[1,2]=ncol(sigma)+p
        rei[2,1]=total_si[list_parentnodes[[m]],1]
        rei[2,2]=ncol(sigma)+p
        relation=rbind(relation,rei) #合成节点矩阵
        childnode=rei[1,1]
        parentnode=total_si[list_parentnodes[[m]],1]
        
        #correct hf 合成节点后，修改父子节点关系
        #合成节点的父节点，就是上述父节点的父节点加上上述子节点的父节点（去掉上述父节点）
        list_nodes_to_parents=list()
        list_nodes_to_parents[[1]]=nodes_to_parents[[parentnode]]
        nodes_to_parents[[childnode]]=setdiff( nodes_to_parents[[childnode]],parentnode)
        if(length( nodes_to_parents[[childnode]])!=0)
        {
          for(l in 1:length( nodes_to_parents[[childnode]]))
          {
            if(!is.element(nodes_to_parents[[childnode]][l],list_nodes_to_parents[[1]]))
            {
              list_nodes_to_parents[[1]][length(list_nodes_to_parents[[1]])+1]=nodes_to_parents[[childnode]][l]
            }
          }
          
        }
        nodes_to_parents=c(nodes_to_parents,list_nodes_to_parents) 
        
        #合成节点的子节点，就是上述父节点的子节点（去掉上述子节点）加上上述子节点的子节点
        nodes_to_children[[parentnode]]=setdiff(nodes_to_children[[parentnode]],childnode)
        list_nodes_to_children=list()
        if( (length(nodes_to_children[[parentnode]])==0)&&is.na(nodes_to_children[[childnode]][1]))
        {
          list_nodes_to_children[[1]]=nodes_to_children[[childnode]]
          nodes_to_children=c(nodes_to_children,list_nodes_to_children)
        }
        
        if( (length(nodes_to_children[[parentnode]])==0)&& !(is.na(nodes_to_children[[childnode]][1])))
        {
          list_nodes_to_children[[1]]=nodes_to_children[[childnode]]
          nodes_to_children=c(nodes_to_children,list_nodes_to_children)
        }
        
        if( (length(nodes_to_children[[parentnode]])!=0)&& (is.na(nodes_to_children[[childnode]][1])))
        {
          list_nodes_to_children[[1]]=nodes_to_children[[parentnode]]
          nodes_to_children=c(nodes_to_children,list_nodes_to_children)
        }
        
        if( (length(nodes_to_children[[parentnode]])!=0)&& !(is.na(nodes_to_children[[childnode]][1])))
        {
          list_nodes_to_children[[1]]=nodes_to_children[[childnode]]
          for(l in 1:length( nodes_to_children[[parentnode]]))
          {
            if(!is.element(nodes_to_children[[parentnode]][l],list_nodes_to_children[[1]]))
            {
              list_nodes_to_children[[1]][length(list_nodes_to_children[[1]])+1]=nodes_to_children[[parentnode]][l]
            }
          }
          nodes_to_children=c(nodes_to_children,list_nodes_to_children)  
        }
        
        
        # 修改上述子节点的子节点的父为合成节点  
        if(!is.na(nodes_to_children[[childnode]][1]))
        {
          for(j in 1:length(nodes_to_children[[childnode]]))
          {
            if(is.element(parentnode,nodes_to_parents[[nodes_to_children[[childnode]][j]]]))
            {
              nodes_to_parents[[nodes_to_children[[childnode]][j]]]=setdiff(nodes_to_parents[[nodes_to_children[[childnode]][j]]],parentnode)
            }
            
            aa=nodes_to_parents[[nodes_to_children[[childnode]][j]]]
            if(childnode %in% aa)
            {
              aa[which(aa==childnode)]=ncol(sigma)+p
              nodes_to_parents[[nodes_to_children[[childnode]][j]]]=aa
            }
          }
        }
        
        # 修改上述父节点的子节点
        if(length(nodes_to_children[[parentnode]])!=0)
        {
          for(j in 1:length(nodes_to_children[[parentnode]]))
          {
            aa=nodes_to_parents[[nodes_to_children[[parentnode]][j]]]
            if(parentnode %in% aa)
            {
              aa[which(aa==parentnode)]=ncol(sigma)+p
              nodes_to_parents[[nodes_to_children[[parentnode]][j]]]=aa
            }
          }  
        }
        
        # 修改上述父节点的父节点 
        if(length(nodes_to_parents[[parentnode]])>0)
        {
          for(j in 1:length(nodes_to_parents[[parentnode]]))
          {
            if(length(nodes_to_children[[nodes_to_parents[[parentnode]][j]]])!=0)
            {
              
              nodes_to_children[[nodes_to_parents[[parentnode]][j]]]=setdiff(nodes_to_children[[nodes_to_parents[[parentnode]][j]]],childnode)
              
              
              
              aa=nodes_to_children[[nodes_to_parents[[parentnode]][j]]]
              if(parentnode %in% aa)
              {
                aa[which(aa==parentnode)]=ncol(sigma)+p
                nodes_to_children[[nodes_to_parents[[parentnode]][j]]]=aa
              }
            }
            
          }
        }
        
        # 修改上述子节点的父节点的子节点
        
        if(length(nodes_to_parents[[childnode]])!=0)
        {
          
          for(j in 1:length(nodes_to_parents[[childnode]]))
          {
            aa=nodes_to_children[[nodes_to_parents[[childnode]][j]]]
            if(childnode %in% aa)
            {
              aa[which(aa==childnode)]=ncol(sigma)+p
              nodes_to_children[[nodes_to_parents[[childnode]][j]]]=aa
            }
          }
        }
        #考虑变成环形的特殊情况
        if(length( nodes_to_parents[[childnode]])!=0 && length(nodes_to_children[[parentnode]])!=0)
        {
          for(r in 1:length(nodes_to_children[[parentnode]]))
          {
            for(q in 1:length(nodes_to_parents[[childnode]]))
            {
              
              
              if(nodes_to_children[[parentnode]][r]==nodes_to_parents[[childnode]][q])
              {
                cyclenode=nodes_to_parents[[childnode]][q]
                
                #去环结构，合成节点作子节点
                nodes_to_children[[ncol(sigma)+p]]=setdiff( nodes_to_children[[ncol(sigma)+p]],cyclenode)
                nodes_to_parents[[cyclenode]]=setdiff( nodes_to_parents[[cyclenode]],ncol(sigma)+p)
                
                
                list_parentnodes=append(list_parentnodes,which(total_si[,1]==cyclenode))
                num_combinenodes=num_combinenodes+1
              }
              
            }
          }
        }
        
        rei[1,1]=ncol(sigma)+p  
        p=p+1
        si=si[-which(si[,1]==childnode),] #在sigma矩阵里，去掉子节点
        
        for(dd in 1:nrow(si))            #在sigma矩阵里，去掉父节点，父节点可能已被去除，加上判断
        {
          if(si[dd,1]==parentnode)
          {
            si=si[-which(si[,1]==parentnode),]
            break
          }
          
        }
        m=m+1
      }
      
    }
    find.max=which(si[,2]==max(si[,2])) #再次寻找最大值
    find.max=sort(find.max,decreasing = F) # #可能有多个值最大且相等，优先取下角标最小的
    find.max=find.max[1]  #sigma最大值的下角标
    
  }
  #针对每个k值 将合成节点矩阵的信息反馈回y矩阵
  for(z in nrow(relation):2)
  {
    if(y[k,relation[z,2]]==1)
    {
      y[k,relation[z,1]]=1
    }
  }
}
y=y[,1:(ncol(sigma)-1)]
#冲突函数
violate.list.final=ViolateDetectlabel(go.for.level,go.leaf.nodes,nodes.to.index,nodes.to.children,y) 
measure.result.fin=MHevaluate(y,test.select.table)


