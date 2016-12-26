work.path="H://R//CODES"
data.path="H://R//DATA"
setwd(work.path)#设置工作路径

#引入库函数以及自定义的函数
source("headfile.R")
BP.univ.graph <- Build.universal.graph.ontology.down(ontology = "BP")#得到BP全图
annotation.final.BP=AnnotationFinal("BP")#提取BP结构下的所有基因及注释信息
go.general.list.BP=Get.GO.all.classes(annotation.final.BP)#根据基因的注释GO标签得到基因的全部GO标签
go.general.table.BP=Build.GO.class.labels(go.general.list.BP)#生成基因及注释信息数据表


file.prefix="derisi"
factor.col=c(0)
file.savepath="H://R//DATA//processalldata"
delete.outlier=FALSE
replace.outlier = FALSE
NAreplace=TRUE
Zrescale=TRUE
write.enable=TRUE
result.process=TRUE
matfile.path="H://R//DATA//matfile"
file.type="change"
file.middle="0"
result.savepath="H://R//RESULT"
  
final.result=DataProcess(file.prefix=file.prefix, factor.col=factor.col, work.path=work.path,data.path=data.path,file.savepath=file.savepath,
                   delete.outlier=delete.outlier,replace.outlier=replace.outlier, NAreplace=NAreplace,
                   Zrescale=Zrescale, write.enable=write.enable,result.process=result.process,
            matfile.path=matfile.path,file.type=file.type, file.middle=file.middle,result.savepath=result.savepath)

