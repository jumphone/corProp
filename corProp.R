

.simple_combine <- function(exp_sc_mat1, exp_sc_mat2, FILL=FALSE){    
    FILL=FILL
    exp_sc_mat=exp_sc_mat1
    exp_ref_mat=exp_sc_mat2
    ##############################################
    if(FILL==TRUE){ 
        gene1=rownames(exp_sc_mat)
        gene2=rownames(exp_ref_mat)
        gene12=gene2[which(!gene2 %in% gene1)]
        gene21=gene1[which(!gene1 %in% gene2)]
        exp_sc_mat_add=matrix(0,ncol=ncol(exp_sc_mat),nrow=length(gene12))
        rownames(exp_sc_mat_add)=gene12
        colnames(exp_sc_mat_add)=colnames(exp_sc_mat)
        exp_ref_mat_add=matrix(0,ncol=ncol(exp_ref_mat),nrow=length(gene21))
        rownames(exp_ref_mat_add)=gene21
        colnames(exp_ref_mat_add)=colnames(exp_ref_mat)
        exp_sc_mat=rbind(exp_sc_mat, exp_sc_mat_add)
        exp_ref_mat=rbind(exp_ref_mat, exp_ref_mat_add)
    }
    ############################################ 
    exp_sc_mat=exp_sc_mat[order(rownames(exp_sc_mat)),]
    exp_ref_mat=exp_ref_mat[order(rownames(exp_ref_mat)),]
    gene_sc=rownames(exp_sc_mat)
    gene_ref=rownames(exp_ref_mat)
    gene_over= gene_sc[which(gene_sc %in% gene_ref)]
    exp_sc_mat=exp_sc_mat[which(gene_sc %in% gene_over),]
    exp_ref_mat=exp_ref_mat[which(gene_ref %in% gene_over),]
    colname_sc=colnames(exp_sc_mat)
    colname_ref=colnames(exp_ref_mat)
    OUT=list()
    OUT$exp_sc_mat1=exp_sc_mat
    OUT$exp_sc_mat2=exp_ref_mat
    OUT$combine=cbind(exp_sc_mat,exp_ref_mat)
    return(OUT)
    }
    

corProp <- function( REF, BULK, S1=1, S2=10, N=1000, SEED=123){
   REF=REF
   BULK=BULK
   SHAPE1=S1
   SHAPE2=S2
   N=N
   SEED=SEED
   ##############################
   print('Beta shape1:')
   print(SHAPE1)
   print('Beta shape2:')
   print(SHAPE2)
   print('Number of random bdg:')
   print(N)
   print('Seed:')
   print(SEED)
   ##############################
   print('building random background...')
   set.seed(SEED)
   BKG=matrix(0,ncol=N,nrow=nrow(REF))
   rownames(BKG)=rownames(REF)
   colnames(BKG)=paste0('BKG',c(1:ncol(BKG)))
   i=1
   while(i<=ncol(BKG)){
       BKG[,i]= as.matrix(REF) %*% runif(ncol(REF))
        i=i+1}
   ##############################
   print('calculating correlation...')
   COM=.simple_combine(BULK, REF)
   COR=cor(COM$exp_sc_mat1, COM$exp_sc_mat2)
   COM.BKG=.simple_combine(BULK, BKG)
   COR.BKG=cor(COM.BKG$exp_sc_mat1, COM.BKG$exp_sc_mat2)
   ##############################
   print('dim(REF):')
   print(dim(REF))
   print('dim(BULK):')
   print(dim(REF))
   print('dim(COM):')
   print(dim(COM$combine))
   ##############################
   print('estimating proportion...')
   PCOR=COR
   NCOR=COR
   i=1
   while(i<=nrow(COR)){
       PCOR[i,]=1-pnorm(COR[i,],mean=mean(COR.BKG[i,]),sd=sd(COR.BKG[i,]))
       NCOR[i,]=qbeta(1-PCOR[i,],shape1=SHAPE1, shape2=SHAPE2)
       i=i+1}
   ##############################
   PNCOR=NCOR/rowSums(NCOR) *100
   ##############################
   print('finished!')
   return(PNCOR)
   }    
    
    
    
