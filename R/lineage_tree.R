#' trace the lineage tree of the lineage analysis result
#' @title lineage_tree
#'
#' @description Using the result of best_iter or multi_best_iter as input, lineage_tree will export  lineage tree of the lineage analysis result.
#' @param result The result of best_iter or multi_best_iter. Required.
#'
#' @examples
#' data("TF1_clones")
#' data=TF1_clones$data
#' ## performs a non-parallel iteration process 
#' result=lineage(data, repeats=30, thread=NULL)
#' ## performs a parallel iteration process
#' # result=lineage(data, repeats=30, thread=10)
#' hc= lineage_tree(result)
#'
#' @return An object of class hclust which describes the lineage tree produced by LINEAGE.
#'
#' @importFrom stats dist hclust as.dist cor
#'
#' @export
#'
lineage_tree=function(result){
    best=result$best
    markers=best$markers
    if(best$suggest=="tsne"){
        diss=1-cor(markers)
    }else{
        diss=dist(t(markers))
    }
    hc=hclust(as.dist(diss))
    return(hc)
}
