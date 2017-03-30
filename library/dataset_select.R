####用于选择使用哪一个数据集的函数 20170330
#### 每个数字对应一个数据集
#### 1 cellcycle 2 derisi 3 eisen 4 gasch1 5 gasch2 6 church 7 spo 8 seq 9 struc 10 hom
#### 目前只实现了1-3 5


DatasetSelect<-function (dataset.index)
{
  if(dataset.index==1)
  {
    file.prefix="cellcycle"
    factor.col=c(0)
  }else if(dataset.index==2)
  {
    file.prefix="derisi"
    factor.col=c(0) 
  }else if(dataset.index==3)
  {
    file.prefix="eisen"
    factor.col=c(0)
  }else if(dataset.index==4)
  {
    
  }else if(dataset.index==5)
  {
    file.prefix="gasch2"
    factor.col=c(0)
  }else if(dataset.index==6)
  {
    
  }else if(dataset.index==7)
  {
    
  }else if(dataset.index==8)
  {
    
  }else if(dataset.index==9)
  {
    
  }else if(dataset.index==10)
  {
    
  }
  result=list(file.prefix=file.prefix,factor.col=factor.col)
  return(result)
}