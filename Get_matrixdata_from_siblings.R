`Get.matrix.data.from.siblings.only` <-
  function(Table.gene.class, classid, data.matrix, ontology="BP", ratio.negative=0, common.genes=NULL, seed=1){
    set.seed(seed);
    
    if (ontology=="BP")
    {
      parent.id <-  AnnotationDbi::get(classid, GOBPPARENTS)
      children.id=AnnotationDbi::get(parent.id, GOBPCHILDREN)
      siblings.id=setdiff(children.id,classid)
    }
     
    else if (ontology=="MF")
    {
      parent.id <-  AnnotationDbi::get(classid, GOMFPARENTS)
      siblings.id=AnnotationDbi::get(parent.id, GOMFCHILDREN)
      siblings.id=setdiff(siblings.id,classid)
    }
      
    else if (ontology=="CC")
    {
      parent.id <-  AnnotationDbi::get(classid, GOCCPARENTS)
      siblings.id=AnnotationDbi::get(parent.id, GOCCCHILDREN)
      siblings.id=setdiff(siblings.id,classid)
    }
      
    else
      stop("Get.matrix.data.from.siblings.only: ontology not valid");
    
    cl.genes <- Extract.class(Table.gene.class,classid);
    positive.gene.names <- as.character(cl.genes$gene.names);
    positive.negative.gene.names <- character(0);
    all.gene.names <- row.names(Table.gene.class);
    
    siblings.id <- intersect(siblings.id,colnames(Table.gene.class));
      	
    
    
    for (i in 1:length(siblings.id)) {
      buffer <- all.gene.names[which(Table.gene.class[,siblings.id[i]] == 1)];
      positive.negative.gene.names <- c(positive.negative.gene.names, buffer);
    }
    positive.negative.gene.names <- unique(positive.negative.gene.names);
    negative.gene.names <- setdiff(positive.negative.gene.names, positive.gene.names);
    
    if (is.null(common.genes)) {
      #gene.names <- c(positive.gene.names, negative.gene.names); 
      gene.names.data<-rownames(data.matrix);
      
      positive.available.genes <- intersect(positive.gene.names,gene.names.data);
      negative.available.genes <- intersect(negative.gene.names,gene.names.data);
    } else {
      positive.available.genes <- intersect(positive.gene.names,common.genes);
      negative.available.genes <- intersect(negative.gene.names,common.genes);	
    }
    
    if (ratio.negative!=0) 
    {
      
      np <- length(positive.available.genes);
      tot.n <- length(negative.available.genes);
      if(tot.n>np)
      {
        nn <- round(np*ratio.negative);
        #if ((np==0) || (nn==0))
          
         # stop("Get.matrix.data.for.classid: there are no examples for at least one of the classes");
        if (nn < tot.n) 
        {
          ind <- sample(tot.n,nn);
          negative.available.genes <- negative.available.genes[ind];	
        }
      }
     
    }	
    
    np <- length(positive.available.genes);
    nn <- length(negative.available.genes);
    
    if ((np==0) || (nn==0))
    {
      cat(classid,np,nn,"\n")
    }
      # stop("Get.matrix.data.from.siblings.only: number of positive or negative examples equal to 0.");
    
    
    exprs.values <- data.matrix[c(positive.available.genes,negative.available.genes),];
    exprs.values[is.na(exprs.values)] <- 0;
    
    y <- as.factor(c(rep(1,np),rep(2,nn)));
    return (list(X=exprs.values,labels=y, n.pos=np, n.neg=nn));
  }

