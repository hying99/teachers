#将原始基因数据读入,目前NA替换，数据归一化、标准化工作已经改在scale函数里实现
ReadData<-function (filename,factor.col,Nametoupper=TRUE,Tonum=TRUE)
{
  mydata <- scan(filename, what = list("", "", ""))#读入待处理数据,用于得到基因名称列表
  #write.table(mydata[[1]],file="newdata",quote = FALSE,row.names = FALSE,col.names = FALSE)
  filedata=read.table(filename,sep = ",",quote = "")#读入待处理数据,用于得到基因数据
  matrix.data=as.matrix(filedata[,-dim(filedata)[2]])#去除样本标签这一列
  col.num=ncol(matrix.data)
  row.num=nrow(matrix.data)
  each.col.status=matrix(0,nrow = col.num,ncol = 1)
  numeric.col=setdiff(c(1:col.num),factor.col)#根据输入的factor属性列号，求numeric属性的列号
  
  
  #将factor类型数据转化为数字型
  for (j in factor.col)
  {
    col.data=filedata[,j]
    matrix.data[,j]=as.numeric(col.data)#将factor类型数据转化为数字型
    #cat("factor\n")
  }
  
  if(Tonum==TRUE)
  {
    matrix.data.final=matrix(as.numeric(matrix.data),nrow=nrow(matrix.data))#将字符矩阵转化为数字矩阵
  }  else
  {
    matrix.data.final=matrix.data 
  }
  
  
  
  if(Nametoupper==TRUE)
  {
    rownames(matrix.data.final)=toupper(substring(mydata[[3]],3))#删除基因名开头的yt字母，并将基因名转为大写
  }  else
  {
    rownames(matrix.data.final)=substring(mydata[[3]],3)#删除基因名开头的yt字母
  }
  matrix.name=rownames(matrix.data.final)
  readdata.result=list(matrix.data.final,matrix.name)
}