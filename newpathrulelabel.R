#用于实现后期标签值的处理 20170310


setwd(data.path)
# mat.result=readMat("test_predict_labels.mat",fixNames = FALSE)
# test.predict.labels=mat.result$predict_labels
# first.predict.labels=test.predict.labels
####第一步 将SVM的概率结果读入

file.type="change"
file.middle="0"
#设置mat文件存储路径
setwd(paste(data.path,"//matfile",sep = ""))
mat.file=paste(file.prefix,file.middle,file.type,"_decision.mat",sep = "")
probability.data=readMat(mat.file,fixNames = FALSE)
#probability.data=readMat("cellcycle0replaceprocess_decision.mat",fixNames = FALSE)
#prob.for.genes=probability.data$decision
prob.for.genes=probability.data$decision_test

####第二步 将概率值转换为标签 1为正 0为负
first.predict.labels=matrix(0,nrow(test.select.table),ncol(test.select.table))
for(i in 1:nrow(prob.for.genes))
{
  for(j in 1:length(except.root.labels.3))
  {
    
    if(prob.for.genes[i,(2*j-1)]>=0.5)
    {
      first.predict.labels[i,j]=1
    }
  }
}

first.measure.result=MHevaluate(first.predict.labels,test.select.table)
violate.detect.first=ViolateDetectlabel(go.for.level.3, go.leaf.nodes.3,nodes.to.index,nodes.to.children,first.predict.labels)

second.predict.labels=first.predict.labels#预测的最终结果将存于final.predict.labels中

##父子节点投票法进行概率修改
# for(k in 1:nrow(prob.for.genes))#对于每一个样本k
# {
#   for (i in (total.levels-1):1)#按层自下而上遍历节点,并且所选层次不为最顶层时
#   { 
#     for(j in 1:length(go.for.level.3[[i]]))
#     {
#       gene.name=(go.for.level.3[[i]])[[j]]#得到(go.for.level[[i]])[[j]]节点的GO标签
#       gene.index=nodes.to.index[[(go.for.level.3[[i]])[j]]]#得到该节点的索引号
#       is.leafnode=go.for.level.3[[i]][[j]] %in% go.leaf.nodes.3#判断该节点是否则为叶子节点
#       if(is.leafnode!=TRUE)#如果这个节点不是叶子节点
#       {
#         if(second.predict.labels[k,gene.index]==0)#如果该节点的预测值为负值
#         {
#           children.index=nodes.to.children[[gene.index]]#得到此节点的所有子节点的序号
#           children.labels=second.predict.labels[k,children.index]#得到其子节点的预测值，此时子节点可能有多个
#           if(1 %in% children.labels)#并且这个节点存在预测值为1的子节点
#           {
#             pos.child.num=sum(children.labels)
#             sup.parent.per=(length(children.index)-pos.child.num)/(length(children.index))
#             pos.child.index=which(children.labels==1)
#             sup.child.per=rep(0,time=pos.child.num)
#             for(m in 1:pos.child.num)
#             {
#               if(i==1)
#               {
#                 sup.child.per[m]=0.5 
#               }else
#               {
#                 parent.index=nodes.to.parents[[(children.index[pos.child.index[m]])]]
#                 parent.labels=second.predict.labels[k,parent.index]
#                 if(length(parent.index)==1)
#                 {
#                   grandpa.index=nodes.to.parents[[parent.index]]
#                   grandpa.labels=second.predict.labels[k,grandpa.index]
#                   sup.child.per[m]=(sum(grandpa.labels))/(length(grandpa.index)+1)
#                 } else
#                 {
#                   sup.child.per[m]=(sum(parent.labels))/(length(parent.index))
#                 }
#               }
#             }
#             for(m in 1:pos.child.num)
#             {
#               if(sup.parent.per>=sup.child.per[m])
#               {
#                 second.predict.labels[k,(children.index[pos.child.index[m]])]=0
#               } else
#               {
#                 second.predict.labels[k,gene.index]=1
#               }
#             }
#           }
#         }
#       }
#     }
#   }
# }
# second.measure.result=MHevaluate(second.predict.labels,test.select.table)
# violate.detect.second=ViolateDetectlabel(go.for.level.3, go.leaf.nodes.3,nodes.to.index,nodes.to.children,second.predict.labels)

#####
# for(k in 1:nrow(prob.for.genes))#对于每一个样本k
# {
#   for (i in 1:total.levels)#按层自上而下遍历节点
#   {
#     for(j in 1:length(go.for.level.3[[i]]))
#     {
#       gene.name=(go.for.level.3[[i]])[[j]]#得到(go.for.level[[i]])[[j]]节点的GO标签
#       gene.index=nodes.to.index[[(go.for.level.3[[i]])[j]]]#得到该节点的索引号
#       ancestor.index=nodes.to.ancestors[[gene.index]]#得到其祖先节点的索引号
#       ancestor.labels=first.predict.labels[k,ancestor.index]#得到其祖先节点的预测值
#       children.index=nodes.to.children[[gene.index]]#得到此节点的所有子节点的序号
#       children.labels=first.predict.labels[k,children.index]#得到其子节点的预测值，此时子节点可能有多个
#       parent.index=nodes.to.parents[[gene.index]]#得到其父节点的索引号
#       parent.labels=first.predict.labels[k,parent.index]#得到其父节点的预测值，此时父节点可能有多个
#       is.leafnode=go.for.level.3[[i]][[j]] %in% go.leaf.nodes.3#判断该节点是否则为叶子节点
#       
#       if(i==1)#对于第一层节点，此层节点假定均属于根结点
#       {
#         if(first.predict.labels[k,gene.index]==0)#如果该节点的预测值为负值
#         {
#           if(is.leafnode!=TRUE)#如果该节点不是叶子节点
#           {
#             if(1 %in% children.labels)#如果该节点的子节点有为1的情况
#             {
#               pos.nums=1+sum(children.labels)#求取父节点中预测为1的个数，及子节点中为1的个数之和
#               neg.nums=length(children.index)-sum(children.labels)
#               #可以看出，当节点属于第一层且不为叶子节点时，
#               #对预测值0是否更改仅取决于其子节点是否有为1的情况
#               if(pos.nums>neg.nums)
#               {
#                 second.predict.labels[k,gene.index]=1
#               } 
#             }
#           }
#           #else{}#如果该节点同时是叶子节点，则保持原值不变
#         }
#       } else#对于属于其它层的节点
#       {
#         if(first.predict.labels[k,gene.index]==0)#如果该节点的预测值为负值
#         {
#           if(is.leafnode!=TRUE)#如果这个节点不是叶子节点
#           {
#             if(1 %in% children.labels)#并且这个节点存在预测值为1的子节点
#             {
#               pos.nums=sum(parent.labels)+sum(children.labels)#求取父节点中预测为1的个数，及子节点中为1的个数
#               #求取父节点中预测为0的个数以及子节点中为0的个数
#               neg.nums=length(parent.index)-sum(parent.labels)+length(children.index)-sum(children.labels)
#               if(sum(children.labels)>0)
#               {
#                 if(pos.nums>neg.nums)
#                 {
#                   second.predict.labels[k,gene.index]=1
#                 } 
#               } else
#               {
#                 if(pos.nums>neg.nums)
#                 {
#                   second.predict.labels[k,gene.index]=1
#                 } 
#               }
#             } 
#           }
#           #else如果这个节点预测值为0，同时是叶子节点，则值保持不变
#           
#         } else #如果该节点的预测值为正值
#         {
#           if(is.leafnode!=TRUE)#如果这个节点不是叶子节点
#           {
#             if(0 %in% parent.labels)#并且这个节点存在预测值为0的父节点
#             {
#               pos.nums=sum(parent.labels)+sum(children.labels)#求取父节点中预测为1的个数，及子节点中为1的个数
#               #求取父节点中预测为0的个数
#               neg.nums=length(parent.index)-sum(parent.labels)+length(children.index)-sum(children.labels)
#               if(sum(children.labels)>0)
#               {
#                 if(neg.nums>pos.nums)
#                 {
#                   second.predict.labels[k,gene.index]=0
#                 } 
#               } else
#               {
#                 if(neg.nums>pos.nums)
#                 {
#                   second.predict.labels[k,gene.index]=0
#                 } 
#               }
#             } 
#           }else#如果这个节点是叶子节点,叶子节点没有子节点
#           {
#             pos.nums=sum(parent.labels)#求取父节点中预测为1的个数
#             neg.nums=length(parent.index)-sum(parent.labels)#求取父节点中预测为0的个数
#             if(neg.nums>pos.nums)
#             {
#               second.predict.labels[k,gene.index]=0
#             } 
#           }
#         }
#       }
#     }
#     
#   }
# }
# 
# 
# second.measure.result=MHevaluate(second.predict.labels,test.select.table)

##自顶而下遍历样本，确定一个节点预测为0或者1时是否需要修改############
for(k in 1:nrow(first.predict.labels))#对于每一个样本k
{
  for (i in 1:total.levels)#按层自上而下遍历节点
  {
    for(j in 1:length(go.for.level.3[[i]]))
    {
      gene.name=(go.for.level.3[[i]])[[j]]#得到(go.for.level[[i]])[[j]]节点的GO标签
      gene.index=nodes.to.index[[(go.for.level.3[[i]])[j]]]#得到该节点的索引号
      ancestor.index=nodes.to.ancestors[[gene.index]]#得到其祖先节点的索引号
      ancestor.labels=first.predict.labels[k,ancestor.index]#得到其祖先节点的预测值
      children.index=nodes.to.children[[gene.index]]#得到此节点的所有子节点的序号
      children.labels=first.predict.labels[k,children.index]#得到其子节点的预测值，此时子节点可能有多个
      parent.index=nodes.to.parents[[gene.index]]#得到其父节点的索引号
      parent.labels=first.predict.labels[k,parent.index]#得到其父节点的预测值，此时父节点可能有多个
      is.leafnode=go.for.level.3[[i]][[j]] %in% go.leaf.nodes.3#判断该节点是否则为叶子节点
      
      if(i==1)#对于第一层节点，此层节点假定均属于根结点
      {
        if(first.predict.labels[k,gene.index]==0)#如果该节点的预测值为负值
        {
          if(is.leafnode!=TRUE)#如果该节点不是叶子节点
          {
            if(1 %in% children.labels)#如果该节点的子节点有为1的情况
            {
              pos.nums=1+sum(children.labels)#求取父节点中预测为1的个数，及子节点中为1的个数之和
              neg.nums=1
              #可以看出，当节点属于第一层且不为叶子节点时，
              #对预测值0是否更改仅取决于其子节点是否有为1的情况
              if(pos.nums>neg.nums)
              {
                second.predict.labels[k,gene.index]=1
              } 
            }
          }
          #else{}#如果该节点同时是叶子节点，则保持原值不变
        }
      } else#对于属于其它层的节点
      {
        if(first.predict.labels[k,gene.index]==0)#如果该节点的预测值为负值
        {
          if(is.leafnode!=TRUE)#如果这个节点不是叶子节点
          {
            if(1 %in% children.labels)#并且这个节点存在预测值为1的子节点
            {
              pos.nums=sum(parent.labels)+sum(children.labels)#求取父节点中预测为1的个数，及子节点中为1的个数
              neg.nums=length(parent.index)-sum(parent.labels)#求取父节点中预测为0的个数
              if(sum(children.labels)>0)
              {
                if(pos.nums>neg.nums)
                {
                  second.predict.labels[k,gene.index]=1
                } 
              } else
              {
                if(pos.nums>neg.nums)
                {
                  second.predict.labels[k,gene.index]=1
                } 
              }
            } 
          }
          #else如果这个节点预测值为0，同时是叶子节点，则值保持不变
          
        } else #如果该节点的预测值为正值
        {
          if(is.leafnode!=TRUE)#如果这个节点不是叶子节点
          {
            if(0 %in% parent.labels)#并且这个节点存在预测值为0的父节点
            {
              pos.nums=sum(parent.labels)+sum(children.labels)#求取父节点中预测为1的个数，及子节点中为1的个数
              neg.nums=length(parent.index)-sum(parent.labels)#求取父节点中预测为0的个数
              if(sum(children.labels)>0)
              {
                if(neg.nums>=pos.nums)
                {
                  second.predict.labels[k,gene.index]=0
                } 
              } else
              {
                if(neg.nums>=pos.nums)
                {
                  second.predict.labels[k,gene.index]=0
                } 
              }
            } 
          }else#如果这个节点是叶子节点,叶子节点没有子节点
          {
            pos.nums=sum(parent.labels)#求取父节点中预测为1的个数
            neg.nums=length(parent.index)-sum(parent.labels)#求取父节点中预测为0的个数
            if(neg.nums>=pos.nums)
            {
              second.predict.labels[k,gene.index]=0
            } 
          }
        }
      }
    }
    
  }
}


second.measure.result=MHevaluate(second.predict.labels,test.select.table)

final.predict.labels=second.predict.labels

# # ##自底而上遍历样本，如果为0的节点存在预测为1的子节点，则判断此节点的预测值是否需要修改############
# for(k in 1:nrow(prob.for.genes))#对于每一个样本k
# {
#   for (i in (total.levels):1)#按层自下而上遍历节点,并且所选层次不为最顶层时
#   { 
#     for(j in 1:length(go.for.level.3[[i]]))
#     {
#       gene.name=(go.for.level.3[[i]])[[j]]#得到(go.for.level[[i]])[[j]]节点的GO标签
#       gene.index=nodes.to.index[[(go.for.level.3[[i]])[j]]]#得到该节点的索引号
#       children.index=nodes.to.children[[gene.index]]#得到此节点的所有子节点的序号
#       children.labels=final.predict.labels[k,children.index]#得到其子节点的预测值，此时子节点可能有多个
#       parent.index=nodes.to.parents[[gene.index]]#得到其父节点的索引号
#       parent.labels=final.predict.labels[k,parent.index]#得到其父节点的预测值，此时父节点可能有多个
#       descendant.index=nodes.to.descendants[[gene.index]]#得到其子孙节点的索引
#       is.leafnode=go.for.level.3[[i]][[j]] %in% go.leaf.nodes.3#判断该节点是否则为叶子节点
#       if(is.leafnode!=TRUE)#如果该节点的预测值为负值
#       {
#         if(final.predict.labels[k,gene.index]==0)#如果这个节点不是叶子节点
#         {
#           if(1 %in% children.labels)#并且这个节点存在预测值为1的子节点
#           {
#             if(i==1)
#             {
#               pos.nums=sum(children.labels)+1
#               neg.nums=0
#               if(pos.nums>(neg.nums+1))
#               {
#                 final.predict.labels[k,gene.index]=1
#               } else
#               {
#                 final.predict.labels[k,descendant.index]=0
#               }
#             }else
#             {
#               pos.nums=sum(parent.labels)+sum(children.labels)#求取父节点中预测为1的个数，及子节点中为1的个数
#               neg.nums=length(parent.index)-sum(parent.labels)#求取父节点中预测为0的个数
#               if(pos.nums>=(neg.nums+1))
#               {
#                 final.predict.labels[k,gene.index]=1
#               }
#               else
#               {
#                 final.predict.labels[k,descendant.index]=0
#               }
#             }
#           }
#         }
#       }
#     }
#   }
# }

##自顶而下遍历样本，将为0的节点的所有子孙节点均置为0############
for(k in 1:nrow(final.predict.labels))#遍历每一个样本
{
  for (i in 1:(length(go.for.level.3)-1))#按层自上而下遍历所有节点
  {
    for(j in 1:length(go.for.level.3[[i]]))#分别遍历每层的节点
    {
      is.leafnode=go.for.level.3[[i]][[j]] %in% go.leaf.nodes.3#判断该节点是否则为叶子节点
      gene.name=go.for.level.3[[i]][[j]]#得到此时处理节点的GO标签名称
      gene.index=nodes.to.index[[(go.for.level.3[[i]])[j]]]#得到该节点在节点列表中的序号
      if(is.leafnode!=TRUE)
      {
        if(final.predict.labels[k,gene.index]==0)#如果该节点的预测值为负值
        {
          #ancestor.index=nodes.to.ancestors[[gene.index]]#得到其祖先节点的索引号号
          children.index=nodes.to.children[[gene.index]]#得到此节点的所有子节点的序号
          children.labels=final.predict.labels[k,children.index]#得到其子节点的预测值，此时子节点可能有多个
          if(1 %in% children.labels)#如果子节点有预测值为正值的情况
          {
            descendant.index=nodes.to.descendants[[gene.index]]#得到其子孙节点的索引
            #final.predict.labels[k,gene.index]=1
            final.predict.labels[k,descendant.index]=0
          }
        }
      }
    }
  }
}

#### 计算结果的评价指标

final.measure.result=MHevaluate(final.predict.labels,test.select.table)

compare.predict.labels=matrix(0,nrow(test.select.table),4)
for(i in 1:nrow(test.select.table))
{
  compare.predict.labels[i,1]=sum(test.select.table[i,])
  compare.predict.labels[i,2]=sum(first.predict.labels[i,])
  compare.predict.labels[i,3]=sum(second.predict.labels[i,])
  compare.predict.labels[i,4]=sum(final.predict.labels[i,])
  
}




















