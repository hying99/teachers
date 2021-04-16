####采用neural netwrok作为基础分类器，运行前需先运行rebuilddataprocess或其他文件，得到数据及类别信息
####采用top down训练方式，每个节点的数据添加其父节点的预测值作为属性
####20171115

model.list=list()
nn.predict.labels=matrix(0,nrow(test.select.table),ncol(test.select.table))
nn.predict.scores=matrix(0,nrow(test.select.table),ncol(test.select.table))
rownames(nn.predict.labels)=rownames(test.select.table)
colnames(nn.predict.labels)=colnames(test.select.table)
rownames(nn.predict.scores)=rownames(test.select.table)
colnames(nn.predict.scores)=colnames(test.select.table)
batch.size=32
epochs=100
patience=5
#每个节点训练一个神经网络，下层神经网络输入数据包含上层神经网络节点的预测值
for(i in 1:(length(go.for.level)))
{
  if(i==1)
  {
    
    for(j in 1:length(go.for.level[[i]]))
    {
      gene.name=go.for.level[[i]][[j]]
      gene.index=nodes.to.index[[gene.name]]
      cat('The processing level is ' ,i,"\n")
      cat('The processing class is ' ,gene.index,"\n")
      train.input.data =data.total[[gene.index]]$X#训练数据
      train.input.label=t(t(as.numeric(data.total[[gene.index]]$labels))) #训练数据的标签
      train.input.label[train.input.label==2]=0
      valid.input.data=valid.select.data
      valid.input.label=t(t(valid.select.table[,gene.index]))
      
      feature.num=ncol(train.input.data)#训练数据属性的个数，即维度
      hneuron.num=round((feature.num+1)/2)#隐藏层节点的数量
      
      model <- keras_model_sequential()
      
      # add layers and compile the model
      model %>%
        layer_dense(units = hneuron.num, activation = 'relu', input_shape = c(feature.num)) %>%
        layer_dropout(rate = 0.5) %>%
        layer_dense(units = hneuron.num, activation = 'relu') %>%
        layer_dropout(rate = 0.5) %>%
        layer_dense(units = 1, activation = 'sigmoid')
      model %>% compile(
        #loss = loss_mean_squared_error,
        loss = loss_binary_crossentropy,
        optimizer = optimizer_adagrad(),
        metrics = c(metric_binary_accuracy)
      )
      history <- model %>% fit(
        train.input.data, train.input.label,
        epochs = epochs, batch_size = batch.size,
        verbose=0,
        #validation_split = 0.2
        validation_data = list(valid.input.data,valid.input.label),
        callbacks=list(callback_early_stopping(monitor="val_loss",min_delta=0,patience=patience,verbose=0,mode=c("auto")))
      )
      
      model.list[[gene.index]]=model
    }
    
  } else
  {
    for(j in 1:length(go.for.level[[i]]))
    {
      gene.name=go.for.level[[i]][[j]]
      gene.index=nodes.to.index[[gene.name]]
      cat('The processing level is ' ,i,"\n")
      cat('The processing class is ' ,gene.index,"\n")
      original.train.data=data.total[[gene.index]]$X#原始训练数据
      train.input.data=NnGetDataset(gene.name,original.train.data,except.root.labels, model.list,go.for.level,each.go.level.num,
                              nodes.to.index,nodes.to.ancestors,nodes.to.parents,batch.size)
      original.valid.data=valid.select.data#原始验证数据
      valid.input.data=NnGetDataset(gene.name,original.valid.data,except.root.labels, model.list,go.for.level,each.go.level.num,
                                    nodes.to.index,nodes.to.ancestors,nodes.to.parents,batch.size)
      train.input.label=t(t(as.numeric(data.total[[gene.index]]$labels))) #训练数据的标签
      train.input.label[train.input.label==2]=0
      valid.input.label=t(t(valid.select.table[,gene.index]))#验证集数据的标签
      
      feature.num=ncol(train.input.data)#训练数据属性的个数，即维度
      hneuron.num=round((feature.num+1)/2)#隐藏层节点的数量
      
      model <- keras_model_sequential()
      
      # add layers and compile the model
      model %>%
        layer_dense(units = hneuron.num, activation = 'relu', input_shape = c(feature.num)) %>%
        layer_dropout(rate = 0.5) %>%
        layer_dense(units = hneuron.num, activation = 'relu') %>%
        layer_dropout(rate = 0.5) %>%
        layer_dense(units = 1, activation = 'sigmoid')
      model %>% compile(
         #loss = loss_mean_squared_error,
       loss = loss_binary_crossentropy,
        optimizer = optimizer_adagrad(),
        metrics = c(metric_binary_accuracy)
      )
      
      history <- model %>% fit(
        train.input.data, train.input.label,
        epochs = epochs, batch_size = batch.size,
        verbose=0,
        #validation_split = 0.2
        validation_data = list(valid.input.data,valid.input.label),
        callbacks=list(callback_early_stopping(monitor="val_loss",min_delta=0,patience=patience,verbose=0,mode=c("auto")))
      )
      
      model.list[[gene.index]]=model
    }
  }
}
names(model.list)=except.root.labels
##按层自上而下进行预测
for(i in 1:(length(go.for.level)))
{
  if(i==1)
  {
    for(j in 1:length(go.for.level[[i]]))
    {
      gene.name=go.for.level[[i]][[j]]
      gene.index=nodes.to.index[[gene.name]]
      cat('The prediction processing level is ' ,i,"\n")
      cat('The prediction processing class is ' ,gene.index,"\n")
      nn.predict.labels[,gene.index]=predict_classes(model.list[[gene.index]],test.select.data,batch_size = batch.size)
      nn.predict.scores[,gene.index]=predict_proba(model.list[[gene.index]],test.select.data,batch_size = batch.size)
    }
  }else
  {
    for(j in 1:length(go.for.level[[i]]))
    {
      gene.name=go.for.level[[i]][[j]]
      gene.index=nodes.to.index[[gene.name]]
      cat('The prediction processing level is ' ,i,"\n")
      cat('The prediction processing class is ' ,gene.index,"\n")
      parent.index=nodes.to.parents[[gene.name]]
      parent.num=length(parent.index)
      test.input.data=test.select.data
      for(k in 1:parent.num)
      {
        test.input.data=cbind(test.input.data,nn.predict.scores[,parent.index[k]])
      }
      nn.predict.labels[,gene.index]=predict_classes(model.list[[gene.index]],test.input.data,batch_size = batch.size)
      nn.predict.scores[,gene.index]=predict_proba(model.list[[gene.index]],test.input.data,batch_size = batch.size)
    }
  }
}


prob.for.genes=matrix(0,nrow(test.select.table),(ncol(test.select.table)*2))
label.length=length(except.root.labels)
for(i in 1:nrow(test.select.table))
{
  for(j in 1:label.length)
  {
    prob.for.genes[i,(2*j-1)]=nn.predict.scores[i,j]
    prob.for.genes[i,(2*j)]=1-nn.predict.scores[i,j]
  }
}



measure.result.nn=MHevaluate(nn.predict.labels,test.select.table)
prauc.result.nn=PRAUCCalculate(nn.predict.scores,test.select.table)
bn.result=BNcompute(prob.for.genes,except.root.labels,go.for.level,go.leaf.nodes,test.select.table)
bn.first.labels=bn.result[[1]]
bn.first.scores=bn.result[[2]]
bn.predict.labels=bn.result[[3]]
bn.predict.scores=bn.result[[4]]
prauc.first.bn=PRAUCCalculate(bn.first.scores,test.select.table)
measure.first.bn=MHevaluate(bn.first.labels,test.select.table)
prauc.result.bn=PRAUCCalculate(bn.predict.scores,test.select.table)
measure.result.bn=MHevaluate(bn.predict.labels,test.select.table)

result.output.en=TRUE
result.savepath=paste(data.path,"//result",sep = "")
setwd(result.savepath)
today <-Sys.Date()

output.fname=paste(file.prefix,"parent_nn_result",".txt",sep = "")
if(result.output.en==TRUE)
{
  
  write.table(today,file=output.fname,sep = " , ",eol="\n",quote=FALSE,row.names = FALSE,col.names = FALSE,append = FALSE)
  write("\n The results given by NN\n",file = output.fname,append = TRUE)
  write.table(measure.result.nn,file=output.fname,sep = " , ",eol="\n",quote=FALSE,row.names = FALSE,col.names = TRUE,append = TRUE)
  write.table(prauc.result.nn,file=output.fname,sep = " , ",eol="\n",quote=FALSE,row.names = FALSE,col.names = TRUE,append = TRUE)
  
  write("\n The results given by BN in the first step\n",file = output.fname,append = TRUE)
  write.table(measure.first.bn,file=output.fname,sep = " , ",eol="\n",quote=FALSE,row.names = FALSE,col.names = TRUE,append = TRUE)
  write.table(prauc.first.bn,file=output.fname,sep = " , ",eol="\n",quote=FALSE,row.names = FALSE,col.names = TRUE,append = TRUE)
  
  write("\n The final results given by BN \n",file = output.fname,append = TRUE)
  write.table(measure.result.bn,file=output.fname,sep = " , ",eol="\n",quote=FALSE,row.names = FALSE,col.names = TRUE,append = TRUE)
  write.table(prauc.result.bn,file=output.fname,sep = " , ",eol="\n",quote=FALSE,row.names = FALSE,col.names = TRUE,append = TRUE)
}

