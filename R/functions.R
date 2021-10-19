# library("Rtsne")
# library("ADPclust")
# library("mclust")
# library("ggplot2")
# library("pheatmap")
# library("pROC")
# library("umap")
# library(parallel)



#' Matrix splitting based on mutation types
#' @title data_prepare
#'
#' @description Split the mitochondrial genotype matrix into 12 submatrices according to mutation types.
#'
#' @param data A dataframe containing mitochondrial variants frequency matrix and two other column "altAllele" and "refAllele". 
#' "altAllele" column represents the mutant allele; "refAllele" column represents the reference allele. Required.
#'
#' @examples
#' data("TF1_clones")
#' data=TF1_clones$data
#' d=data_prepare(data)
#'
#' @return A list containing 12 submatrices with different mutation types. ("altAllele" column and "refAllele" column are removed)
#'
#' @export
#'
data_prepare=function(data){
    df_a=data[grep("A",data$altAllele),]
    df_t=data[grep("T",data$altAllele),]
    df_c=data[grep("C",data$altAllele),]
    df_g=data[grep("G",data$altAllele),]
    df_a_rc=df_a[grep("C",df_a$refAllele),]
    df_a_rt=df_a[grep("T",df_a$refAllele),]
    df_a_rg=df_a[grep("G",df_a$refAllele),]
    df_c_ra=df_c[grep("A",df_c$refAllele),]
    df_c_rt=df_c[grep("T",df_c$refAllele),]
    df_c_rg=df_c[grep("G",df_c$refAllele),]
    df_t_ra=df_t[grep("A",df_t$refAllele),]
    df_t_rc=df_t[grep("C",df_t$refAllele),]
    df_t_rg=df_t[grep("G",df_t$refAllele),]
    df_g_ra=df_g[grep("A",df_g$refAllele),]
    df_g_rc=df_g[grep("C",df_g$refAllele),]
    df_g_rt=df_g[grep("T",df_g$refAllele),]
    d=list()
    d[[1]]=subset(df_t_rg,select=-c(altAllele, refAllele))
    d[[2]]=subset(df_t_rc,select=-c(altAllele, refAllele))
    d[[3]]=subset(df_t_ra,select=-c(altAllele, refAllele))
    d[[4]]=subset(df_a_rg,select=-c(altAllele, refAllele))
    d[[5]]=subset(df_a_rt,select=-c(altAllele, refAllele))
    d[[6]]=subset(df_a_rc,select=-c(altAllele, refAllele))
    d[[7]]=subset(df_c_rt,select=-c(altAllele, refAllele))
    d[[8]]=subset(df_c_ra,select=-c(altAllele, refAllele))
    d[[9]]=subset(df_c_rg,select=-c(altAllele, refAllele))
    d[[10]]=subset(df_g_rc,select=-c(altAllele, refAllele))
    d[[11]]=subset(df_g_rt,select=-c(altAllele, refAllele))
    d[[12]]=subset(df_g_ra,select=-c(altAllele, refAllele))
    return(d)
}



#' Mitochondrial genotype matrix modification
#' @title correct
#'
#' @description Using mitochondrial variants frequency matrix as input, transform the zeroes into ones when the median of the non-zero frequencies in the same row >= 0.6.
#'
#' @param x A dataframe containing mitochondrial variants frequency matrix where a column represented a single cell and a row represented variants frequency of a specific mitochondrial genotype. Required.
#'
#' @examples
#' data("TF1_clones")
#' data=TF1_clones$data
#' d=data_prepare(data)
#' d2=data.frame()
#' for(i in 1:12){
#'     di=d[[i]]
#'     di2=as.data.frame(t(apply(di,1,correct)))
#'     d2=rbind(d2,di2)
#' }
#'
#' @return A dataframe containing a modified mitochondrial genotype matrix.
#'
#' @importFrom stats median
#'
#' @export
#'
correct=function(x){
    x1=x[x!=0]
    if(!is.na(median(x1)) & median(x1)>=0.6){
        x[x==0]=1
    }
    return(t(x))
}



#' Highly variable sites identification
#' @title top_sites
#'
#' @description Using mitochondrial variants frequency matrix as input, the 20/50 highest variable sites were called by identifying the rows with highest standard deviation (sd) across columns(cells).
#'
#' @param d A dataframe containing mitochondrial variants frequency matrix where a column represented a single cell and a row represented variants frequency of a specific mitochondrial genotype. Required.
#'
#' @examples
#' data("TF1_clones")
#' data=TF1_clones$data
#' d=data_prepare(data)
#' tops=data.frame()
#' tops2=data.frame()
#' for(i in 1:12){
#'     di=d[[i]]
#'     tops=rbind(tops,top_sites(di)[[1]])
#'     tops2=rbind(tops2,top_sites(di)[[2]])
#' }
#'
#' @return A dataframe containing 2 modified mitochondrial genotype matrix with 20/50 highest variable sites.
#'
#' @importFrom stats sd
#'
#' @export
#'
top_sites=function(d){
    d_1=as.data.frame(t(apply(d,1,correct)))
    names(d_1)=names(d)
    a1=apply(d_1,1,sd)
    sites=row.names(d)
    sites=sites[order(a1,decreasing=TRUE)]
    d1=d[sites[1:20],]
    d2=d[sites[1:50],]
    return(list(d1,d2))
}



#' Distance calculation in subspaces
#' @title learn_distance
#'
#' @description By integrating distance information and consensus information, learn_distance generate a more integrative distance matrix.
#'
#' @param c Integer. 1:(numer of selected subspaces). Required.
#' @param pclusters A list containing clustering results in selected subspaces with kmeans. Required.
#' @param diss A list containing Euclidean distance matrixes of selected subspaces. Required.
#'
#' @examples
#' data("pcluster")
#' data("diss")
#' data("cs")
#' c=1:length(cs)
#' dis=learn_distance(c,pcluster,diss)
#'
#' @return A dataframe containing a more integrative distance matrix.
#'
#' @export
#'
learn_distance=function(c,pclusters,diss){
    length=length(pclusters[[1]])
    w=matrix(0,length,length)
    dis=matrix(0,length,length)
    for(i in c){
        cluster=pclusters[[i]]
        pclusters[[i]]=NA
        j1=1
        while (j1 <length){
            j2=j1+1
            while(j2 <= length){
                if(cluster[j1]==cluster[j2]){
                    w[j1,j2]=w[j1,j2]+1
                    w[j2,j1]=w[j2,j1]+1
                }
                j2=j2+1
            }
            j1=j1+1
        }
        dis=dis+as.matrix(diss[[i]])
    }
    dis=as.matrix(dis)*(1/(3*w+1))
    return(dis)
}



#' Marker variants identification
#' @title find_markers
#'
#' @description find_markers performs evaluation of variant as marker for the consensus clustering result by calculating the receiver operating characteristics (AUC) score and Pearson's correlation coefficient between the binary cluster labels and the frequency distribution.
#'
#' @param f A list containing the frequency distribution of a variant. Required.
#' @param cluster A list containing the consensus clustering result. Required.
#'
#' @examples
#' data("TF1_clones")
#' data("scc")
#' data=TF1_clones$data
#' d=data_prepare(data)
#' tops2=data.frame()
#' for(i in 1:12){
#'     di=d[[i]]
#'     tops2=rbind(tops2,top_sites(di)[[2]])
#' }
#' find_markers_2=function(x,cluster=scc){
#'     a=find_markers(as.numeric(as.character(x)),cluster=as.numeric(as.character(scc)))
#'     a=a[a$p==min(a$p),]
#'     return(a)
#' }
#' markers=apply(tops2,1,find_markers_2)
#'
#' @return A dataframe containing AUC score and p-value of the given variant.
#'
#' @importFrom pROC roc
#' @importFrom stats cor.test
#'
#' @export
#'
find_markers=function(f,cluster){
    clu=max(cluster)
    clus=1:clu
    marker=data.frame()
    for(i in clus){
        ncluster=cluster
        nncluster=ncluster
        ncluster[nncluster!=as.numeric(i)]=0
        ncluster[nncluster==as.numeric(i)]=1
        rm(nncluster)
        comp=data.frame(f,ncluster)
        names(comp)=c("f","ncluster")
        roc1=roc(comp$ncluster,comp$f,quiet=TRUE)
        auc1=as.numeric(roc1$auc)
        p=cor.test(comp$ncluster,comp$f)$p.value
        result=c(auc1,p)
        marker=rbind(marker,result)
        }
        names(marker)=c("cor","p")
        return(marker)
}



#' Sscore and Dscore calculation for optimization
#' @title sum_score
#'
#' @description Calculate Sscore and Dscore for iterative optimization.
#'
#' @param score A dataframe containing adjusted rand indexes between the refined clustering results (resulted from t-SNE or UMAP-Kmeans procedures) and the clustering results in the selected subspaces. Required.
#' @param cs Indexes of selected subspaces. Required.
#' @param aris A dataframe containing ARI values among effective subspaces. Required.
#'
#' @examples
#' data("score")
#' data("cs")
#' data("aris")
#' score_sug=sum_score(score, cs, aris)
#'
#' @return A list containing Sscore, suggested dimensional reduction method, number of effective subspaces and average ARI values used in Dscore calculation.
#'
#' @export
#'
sum_score=function(score,cs,aris){
    temp=cs
    temp2=aris
    s=apply(score,2,max)
    cs=cs[s>0.1]
    aris=aris[cs]
    if(length(cs)==1){
        saris=0.5
    }else{
        aris[aris<0]=0.5
        saris=mean(aris)
    }
    score1=score[1,]
    sum_score1=length(cs)+2*(1-saris)+max(score1[score1>0.1])
    score2=score[2,]
    sum_score2=length(cs)+2*(1-saris)+max(score2[score2>0.1])
    if(is.na(sum_score1)){
        sum_score1=0
    }
    if(is.na(sum_score2)){
        sum_score2=0
    }
    if(sum_score1>sum_score2){
        sum_score=sum_score1
        suggest="tsne"
    }else{
        sum_score=sum_score2
        suggest="umap"
    }
    result=list(sum_score=sum_score,suggest=suggest,lcs=length(cs),saris=saris)
    return(result)
}



#' Main function
#' @title main
#'
#' @description main is the main part of LINEAGE for lineage analysis.
#'
#' @param d A list containing 12 submatrices with different mutation types. Output of data_prepare(). Required.
#' @param centers Integer. The number of clusters used in Kmeans procedure. Default: 3.
#' @param nmarker Integer. The number of markers showed in final result. Default: 16.
#'
#' @examples
#' data("TF1_clones")
#' data=TF1_clones$data
#' d=data_prepare(data)
#' a=dim(d[[1]])[2]
#' if(a>100){
#'     nmarker=c(15,20)
#'     centers=3
#' }else{
#'     nmarker=c(10,15)
#'     centers=2
#' }
#' result=main(d, centers, nmarker)
#'
#' @return A list containing lineage analysis result.
#'
#' @importFrom Rtsne Rtsne
#' @importFrom umap umap
#' @importFrom ADPclust adpclust
#' @importFrom mclust adjustedRandIndex
#' @importFrom ggplot2 ggplot aes facet_wrap geom_point
#' @importFrom stats cutree
#' @importFrom stats dist hclust kmeans as.dist cor
#'
#' @export
#'
main=function(d=d,centers=3,nmarker=16){
    # Highly variable sites identification
    tops=data.frame()
    tops2=data.frame()
    for(i in 1:12){
        di=d[[i]]
        tops=rbind(tops,top_sites(di)[[1]])
        tops2=rbind(tops2,top_sites(di)[[2]])
    }

    # Subspace separation based on the frequency dynamic patterns of highly variable sites
    corr=cor(t(tops))
    diss=1-corr
    hc=hclust(as.dist(diss))
    dds=list()
    for(i in 1:20){
        dd=tops[cutree(hc,20)==i,]
        dds[[i]]=dd
    }

    # Clustering in subspaces with kmeans
    clusters=list()
    dis1=list()
    tds1=list()
    perplexity=round((dim(dds[[1]])[2]-1)/3,0)-1
    if(perplexity>30){
        perplexity=30
    }else if(perplexity>20){
        perplexity=20
    }
    for(i in 1:20){
        dis1[[i]]=dist(t(dds[[i]]))
        tds1[[i]]=as.data.frame(Rtsne(dis1[[i]],perplexity=perplexity,is_distance=TRUE,theta=0)$Y)
        temp=kmeans(tds1[[i]],centers,nstart=20)
        clusters[[i]]=temp$cluster
    }

    # Entropy evaluation
    aris=rep(0,20*20)
    dim(aris)=c(20,20)
    for(i in 1:19){
        for(j in (i+1):20){
            ari=adjustedRandIndex(clusters[[i]],clusters[[j]])
            aris[i,j]=ari
            aris[j,i]=ari
        }
    }
    ari=sort(apply(aris,1,max),decreasing=TRUE)
    if(ari[6]>=0.05){
        cs=order(apply(aris,1,max),decreasing=TRUE)[1:6]
    }else{
        cs=cs[ari>0.05]
    }

    # Plot the clustering in selected subspaces
    tds2=tds1[cs]
    tsnes=do.call(rbind,tds2)
    tsnes$group=rep(cs,each=dim(tds1[[1]])[1])
    names(tsnes)=c("tsne_1","tsne_2","group")
    p1=ggplot(data=tsnes,aes(x=tsne_1,y=tsne_2))+facet_wrap(.~group,ncol=6)+geom_point(size=0.5)

    #sub_tsnes=do.call(rbind,tds1)
    #sub_tsnes$group=rep(1:20,each=dim(tds1[[1]])[1])

    # Consensus clustering
    pcluster=list()
    diss=dis1[cs]
    nu=1
    for(i in cs){
        data=tds1[[i]][,1:2]
        sc =kmeans(data,centers,nstart=20)$cluster
        pcluster[[nu]]=sc
        nu=nu+1
    }
    c=1:length(cs)
    dis=learn_distance(c,pcluster,diss)

    # Refinement based on marker variants
    perplexity=round((dim(dis)[1]-1)/3,0)-1
    if(perplexity>30){
        perplexity=30
    }else if(perplexity>20){
        perplexity=20
    }
    td=Rtsne(dis,perplexity=perplexity,is_distance=TRUE,theta=0)$Y
    td=as.data.frame(td)
    center=max(adpclust(td[,1:2], centroids = "auto",htype ="AMISE", nclust = 3:8)$clusters)
    scc=kmeans(td[,1:2],center,nstart=20)$cluster
    td$group=as.character(scc)
    names(td)=c("tsne1","tsne2","group")
    p2=ggplot(td,aes(x=tsne1,y=tsne2,col=group))+geom_point()

    # Group marker identification
    find_markers_2=function(x,cluster=scc){
        a=find_markers(as.numeric(as.character(x)),cluster=as.numeric(as.character(scc)))
        a=a[a$p==min(a$p),]
        return(a)
    }
    markers=apply(tops2,1,find_markers_2)
    markers2=do.call(rbind,markers)
    markers3=markers2[markers2$p<0.05,]
    markers3=markers3[order(markers3$cor,decreasing=TRUE)[1:nmarker],]
    markers1=tops2[row.names(markers3),]
    use=apply(markers1,2,sum)!=0
    markers1=markers1[,use]

    # Dimensionality reduction for visualization
    #diss=dist(t(markers1))
    diss=1-cor(markers1)
    perplexity=round((dim(diss)[1]-1)/3,0)-1
    if(perplexity>30){
        perplexity=30
    }else if(perplexity>20){
        perplexity=20
    }
    td2=Rtsne(diss,is_distance=TRUE,perplexity=perplexity,theta=0)$Y
    td2=as.data.frame(td2)
    names(td2)=c("tsne1","tsne2")

    td3=umap(t(markers1))$layout
    td3=as.data.frame(td3)
    names(td3)=c("umap1","umap2")


    ######
    clu=list()
    clu[[1]]=kmeans(td2[,1:2],centers,nstart=20)$cluster
    clu[[2]]=kmeans(td3[,1:2],centers,nstart=20)$cluster
    length=length(cs)
    score=rep(0,length*2)
    dim(score)=c(2,length)
    for(i in 1:2){
        for(j in 1:length){
            ari=adjustedRandIndex(clu[[i]],clusters[[cs[j]]][use])
            score[i,j]=ari
        }
    }
    #print(score)
    #print(cs)
    score_sug=sum_score(score,cs,aris)
    score=score_sug$sum_score+2.5*sum(markers3$cor[1:10])/10
    score1=score_sug$sum_score
    score2=sum(markers3$cor[1:10])
    suggest=score_sug$suggest
    #print(score)
    result=list(tsnes=td2,umaps=td3,markers=markers1,use=use,suggest=suggest,score=score,fig1=p1,fig2=p2,centers=centers, nmarker=nmarker,score1=score1,score2=score2,nsubspace=score_sug$lcs,saris=score_sug$saris)
    return(result)
}



#' Non-parallel iterative optimization
#' @title best_iter
#'
#' @description best_iter performs a non-parallel iteration process of main to guarantee a more stable and reliable cell clustering/clone tracing result.
#'
#' @param d A list containing 12 submatrices with different mutation types. Output of data_prepare(). Required.
#' @param centers Integer. The number of clusters used in Kmeans procedure. Default: c(2,3).
#' @param nmarker Integer. The number of markers showed in final result. Default: c(10,15,20).
#' @param repeats Integer. The number of iterations. Default: 30.
#'
#' @examples
#' data("TF1_clones")
#' data=TF1_clones$data
#' d=data_prepare(data)
#' a=dim(d[[1]])[2]
#' if(a>100){
#'     nmarker=c(15,20)
#'     centers=3
#' }else{
#'     nmarker=c(10,15)
#'     centers=2
#' }
#' all_results=best_iter(d, centers, nmarker, repeats=30)
#'
#' @return A list containing results of all repeats, including the best result.
#'
#' @importFrom ADPclust adpclust
#' @importFrom stats kmeans
#'
#' @export
#'
best_iter=function(d,centers=c(2,3),nmarker=c(10,15,20),repeats=30){
    tresult=list()
    scores=c()
    nu=1
    print("start to choose best parameters for c(centers, nmarker)")
    for(i in centers){
        for(j in nmarker){
            for(k in 1:repeats){
                print(sprintf("run the %d th iteration",nu))
                fit=try(main(d=d,centers=i,nmarker=j))
                trys=1
                while(("try-error" %in% class(fit)) & trys<=10){
                    fit=try(main(d=d,centers=i,nmarker=j))
                    trys=trys+1
                }
                if(trys>10){
                    tresult[[nu]]=NA
                    scores[nu]=0
                }else{
                    tresult[[nu]]=fit
                    scores[nu]=tresult[[nu]]$score
                }
                nu=nu+1
            }
        }
    }
    best=tresult[[order(scores,decreasing=TRUE)[1]]]
    i=best$centers
    j=best$nmarker
    print(sprintf("the best center is %d",i))
    print(sprintf("the best num of markers id %d",j))

    if(best$suggest=="tsne"){
        td=best$tsnes
        fcenter=max(adpclust(td[,1:2], centroids = "auto",htype ="AMISE", nclust = 3:8)$clusters)
        label=kmeans(td[,1:2],fcenter,nstart=20)$cluster
    }else{
        td=best$umaps
        fcenter=max(adpclust(td[,1:2], centroids = "auto",htype ="AMISE", nclust = 3:8)$clusters)
        label=kmeans(td[,1:2],fcenter,nstart=20)$cluster
    }
    print(">score results of all repeats with different 'centers' and 'num of markers'")
    print("centers num_of_markers Sscore(score) Dscore(score1) averageAUC(score2) num_of_selected_subspaces")
    for(i in 1:(nu-1)){
        print(c(tresult[[i]]$centers,tresult[[i]]$nmarker,as.numeric(tresult[[i]]$score),as.numeric(tresult[[i]]$score1),as.numeric(tresult[[i]]$score2),as.numeric(tresult[[i]]$nsubspace)))
    }
    print(">max Sscore of all repeats with different 'centers' and 'num of markers'")
    print(scores[order(scores,decreasing=TRUE)[1]])
    #print("*********")
    label=as.data.frame(label)
    label$Sample=names(d[[1]])[best$use]
    names(label)=c("Label","Sample")
    result=list(best=best,label=label,suggest=best$suggest,results=tresult,scores=scores)
    return(result)
}



#' Parallel iterative optimization
#' @title multi_best_iter
#'
#' @description multi_best_iter performs a parallel iteration process of main to guarantee a more stable and reliable cell clustering/clone tracing result.
#'
#' @param d A list containing 12 submatrices with different mutation types. Output of data_prepare(). Required.
#' @param centers Integer. The number of clusters used in Kmeans procedure. Default: c(2,3).
#' @param nmarker Integer. The number of markers showed in final result. Default: c(10,15,20).
#' @param repeats Integer. The number of iterations. Default: 30.
#' @param thread Integer. Integer. The number of threads to run multi_best_iter. Default: 10.
#'
#' @examples
#' data("TF1_clones")
#' data=TF1_clones$data
#' d=data_prepare(data)
#' a=dim(d[[1]])[2]
#' if(a>100){
#'     nmarker=c(15,20)
#'     centers=3
#' }else{
#'     nmarker=c(10,15)
#'     centers=2
#' }
#' # all_results=multi_best_iter(d, centers, nmarker, repeats=30, thread=10)
#'
#' @return A list containing results of all repeats, including the best result.
#'
#' @importFrom parallel makeCluster clusterExport parLapply stopCluster
#' @importFrom stats kmeans
#'
#' @export
#'

multi_best_iter=function(d,centers=c(2,3),nmarker=c(10,15,20),repeats=30,thread=10){
    tresult=list()
    scores=c()
    nu=1
    print("start to choose best parameters for c(centers, nmarker)")
    main_try=function(k){
        fit=try(main(d=d,centers=i,nmarker=j))
        trys=1
        while(("try-error" %in% class(fit)) & trys<=10){
            fit=try(main(d=d,centers=i,nmarker=j))
            trys=trys+1
        }
        if(trys>10){
            fit=list(result="NA",score=0)
        } 
        return(fit)
    }

    cl <- makeCluster(getOption("cl.cores", thread),type="FORK")
    for(i in centers){
        for(j in nmarker){
            clusterExport(cl,"d",envir = environment())
            clusterExport(cl,"i",envir = environment())
            clusterExport(cl,"j",envir = environment())
            c=parLapply(cl,1:repeats,main_try)
            for(z in 1:repeats){
                tresult[[nu]]=c[[z]]
                scores[nu]=tresult[[nu]]$score
                nu=nu+1
            }
        }
    }
    stopCluster(cl)
    best=tresult[[order(scores,decreasing=TRUE)[1]]]
    i=best$centers
    j=best$nmarker
    print(sprintf("the best center is %d",i))
    print(sprintf("the best num of markers id %d",j))
    if(best$suggest=="tsne"){
        td=best$tsnes
        fcenter=max(adpclust(td[,1:2], centroids = "auto",htype ="AMISE", nclust = 3:8)$clusters)
        label=kmeans(td[,1:2],fcenter,nstart=20)$cluster
    }else{
        td=best$umaps
        fcenter=max(adpclust(td[,1:2], centroids = "auto",htype ="AMISE", nclust = 3:8)$clusters)
        label=kmeans(td[,1:2],fcenter,nstart=20)$cluster
    }
    print(">score results of all repeats with different 'centers' and 'num of markers'")
    print("centers num_of_markers Sscore(score) Dscore(score1) averageAUC(score2) num_of_selected_subspaces")
    for(i in 1:(nu-1)){
        print(c(tresult[[i]]$centers,tresult[[i]]$nmarker,as.numeric(tresult[[i]]$score),as.numeric(tresult[[i]]$score1),as.numeric(tresult[[i]]$score2),as.numeric(tresult[[i]]$nsubspace)))
    }
    label=as.data.frame(label)
    label$Sample=names(d[[1]])[best$use]
    names(label)=c("Label","Sample")
    print(">max Sscore of all repeats with different 'centers' and 'num of markers'")
    print(scores[order(scores,decreasing=TRUE)[1]])
    #print("*********")
    result=list(best=best,label=label,suggest=best$suggest,results=tresult,scores=scores)
    return(result)
}

