# library("ggplot2")
# library("pheatmap")

#' Visualization of the result
#' @title traceplot
#'
#' @description Using the result of best_iter or multi_best_iter as input, traceplot will export a heatmap and a suggested 2D visualization plot for the lineage analysis result.
#' @param result The result of best_iter or multi_best_iter. Required.
#' @param label A dataframe containing lineage label with at least 2 columns named "Sample" and "Label". "Sample" column contains sample ids; "label" column is the lineage label of each sample. Required.
#'
#' @examples
#' data("TF1_clones")
#' data=TF1_clones$data
#' rlabel=TF1_clones$rlabel
#' ## performs a non-parallel iteration process 
#' result=lineage(data, repeats=30, thread=NULL)
#' ## performs a parallel iteration process
#' # result=lineage(data, repeats=30, thread=10)
#' ## plots with inferred clone labels
#' plots=traceplot(result = result, label= result$label)
#' ## plots with reference clone labels
#' plots=traceplot(result = result, label= rlabel)
#'
#' @return A list containing a heatmap and a 2D visualization plot for the result of LINEAGE.
#'
#' @importFrom ggplot2 ggplot aes geom_point
#' @importFrom pheatmap pheatmap
#' @importFrom grDevices dev.off pdf
#'
#' @export
#'
traceplot=function(result,label){
    best=result$best
    tsnes=best$tsnes
    umaps=best$umaps
    rownames(label)=label$Sample
    use=as.data.frame(best$use)
    names(use)="use"
    label=label[rownames(use),][best$use,]
    pdf(file="2d.pdf")
    if(best$suggest=="tsne"){
        td=tsnes
        td$label=as.character(label$Label)
        names(td)=c("tsne1","tsne2","label")
        p1=ggplot(td,aes(x=tsne1,y=tsne2,col=label))+geom_point()
    }else{
        td=umaps
        td$label=as.character(label$Label)
        names(td)=c("umap1","umap2","label")
        p1=ggplot(td,aes(x=umap1,y=umap2,col=label))+geom_point()
    }
    print(p1)
    dev.off()
    markers=best$markers
    annotation_col=as.data.frame(as.character(label$Label))
    row.names(annotation_col)=label$Sample
    names(annotation_col)=c('label')
    pdf(file="heatmap.pdf")
    if(best$suggest=="tsne"){
        p2=pheatmap(markers,annotation_col=annotation_col,scale="row",clustering_distance_cols="correlation",clustering_distance_rows="correlation")
    }else{
        p2=pheatmap(markers,annotation_col=annotation_col,scale="row",clustering_distance_cols="euclidean")
    }
    print(p2)
    dev.off()
    plots=list(d2d=p1,heatmap=p2)
    return(plots)
}
