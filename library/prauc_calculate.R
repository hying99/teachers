PRAUCCalculate <- function(predict,target)
{
  #求取AU（PRC）的值
  pr.result=PrecisionRecallCalculate(predict,target,matrix.en=TRUE)
  #precision.recall.curves.plot(list(pr.result),plot.precision=TRUE)
  AU.prc=AUPRC(list(pr.result), comp.precision=TRUE)[[1]]
  #求取各类AUPRC的平均值和加权值
  label.num=ncol(target)
  each.class.auc=numeric(length = label.num)
  each.class.freq=numeric(length = label.num)
  each.class.weight=numeric(length = label.num)
  for(i in 1:label.num)
  {
    single.predict.scores=predict[,i]
    single.target.label=target[,i]
    single.pr.result=PrecisionRecallCalculate(single.predict.scores,single.target.label,matrix.en=FALSE)
    each.class.auc[i]=AUPRC(list(single.pr.result), comp.precision=TRUE)
    each.class.freq[i]=sum(single.target.label)
  }
  average.auprc=sum(each.class.auc)/label.num
  
  average.auprc.w=0
  for(i in 1:label.num)
  {
    each.class.weight[i]=each.class.freq[i]/sum(each.class.freq)
    average.auprc.w=average.auprc.w+each.class.auc[i]*each.class.weight[i]
  }
  result=list(AU.prc=AU.prc,average.auprc=average.auprc,average.auprc.w=average.auprc.w)
  return(result)
}