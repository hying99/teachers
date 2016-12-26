
#用于数据训练之前的分析
matrix.data=matrix.cellcycle.data
training.original.data=matrix.cellcycle.data
col.num=ncol(matrix.data)
row.num=nrow(matrix.data)
factor.col=c()
numeric.col=setdiff(c(1:col.num),factor.col)
#绘制散点图和箱形图
sp=list()
sp_name=c()
sp.eachnum=matrix(0,nrow = col.num,ncol = 1)
for(j in numeric.col)
{
  #dotchart(matrix.cellcycle.data[1:10,1])
  sp[[j]]=boxplot(matrix.data[,j],horizontal=TRUE,plot=TRUE)
  sp_name=c(sp_name,names(sp[[j]]$out))
  sp.eachnum[j,1]=length(names(sp[[j]]$out))
}

for(j in numeric.col)
{
  each.col.min=sp[[j]]$stats[1]
  each.col.max=sp[[j]]$stats[5]
  for(i in 1:row.num)
  {
    if(!is.na(training.original.data[i,j]))
    {
      if(training.original.data[i,j]>each.col.max)
      {
        matrix.data[i,j]=NA
      } else 
      {
        if(training.original.data[i,j]<each.col.min)
        {
          matrix.data[i,j]=NA
        }
      } 
    }
    
  }
}

each.col.status=matrix(0,nrow = col.num,ncol = 3)
#计算每列的最小值、最大值、平均值
for(j in numeric.col)
{
  each.col.status[j,1]=min(matrix.data[,j],na.rm = TRUE)
  each.col.status[j,2]=max(matrix.data[,j],na.rm = TRUE)
  each.col.status[j,3]=mean(matrix.data[,j],na.rm = TRUE)
}

#连续属性，将NA值用对应列的均值替换
for(i in 1:row.num)
{
  for(j in numeric.col)
  {
    if(is.na(matrix.data[i,j]))
    {
      matrix.data[i,j]=each.col.status[j,3]
    }
  }
}
#剔除具有异常值的样本
sp_name=unique(sp_name)
matrix.name=rownames(matrix.data)
select.name=setdiff(matrix.name,sp_name)
remain.data=matrix.data[select.name,]

j=4
boxplot(matrix.data[,1],horizontal=TRUE,plot=TRUE)





#判断各属性的类型，并记录相应位置
for (j in 1:ncol(matrix.data))
{
  col.data=matrix.data[,j]
  if(is.factor(col.data))
  {
    factor.col=c(factor.col,j)
  }
  else
  {
    numeric.col=c(numeric.col,j)
  }
}

#将factor类型数据转化为数字型
for (j in factor.col)
{
  col.data=matrix.data[,j]
  matrix.data[,j]=as.numeric(col.data)#将factor类型数据转化为数字型
  cat("factor\n")
}




#计算每列的最小值、最大值、平均值
remain.col.status=matrix(0,nrow = col.num,ncol = 3)
for(j in numeric.col)
{
  remain.col.status[j,1]=min(remain.data[,j],na.rm = TRUE)
  remain.col.status[j,2]=max(remain.data[,j],na.rm = TRUE)
  remain.col.status[j,3]=mean(remain.data[,j],na.rm = TRUE)
}

valid.col.status=matrix(0,nrow = col.num,ncol = 3)
for(j in numeric.col)
{
  valid.col.status[j,1]=min(valid.cellcycle.data[,j],na.rm = TRUE)
  valid.col.status[j,2]=max(valid.cellcycle.data[,j],na.rm = TRUE)
  valid.col.status[j,3]=mean(valid.cellcycle.data[,j],na.rm = TRUE)
}




sub=list()
for(i in 1:col.num)
{
  sub[[i]]=which(is.na(matrix.cellcycle.data[1,i]))
  cat(i,length(sub[[i]]))
  cat("\n")
  
}
total_index=sub[[1]]
for(i in 2:col.num)
{
  total_index=union(total_index,sub[[i]])
  
}



