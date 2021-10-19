#' @title lineage
#'
#' @description 'lineage' can perform label-free identification of endogenous informative single cell mitochondrial RNA mutation for lineage analysis and clonal evolution. 
#' Using mitochondrial genotype matrix as input, select the informative variants of mitochondrial RNA for clones and perform lineage analysis with a consensus clustering method.
#'
#' @param data A dataframe containing mitochondrial variants frequency matrix, where a column represents a single cell and a row represents variants frequency of a specific mitochondrial genotype, and two other columns "altAllele" and "refAllele". 
#' "altAllele" column represents the mutant allele; "refAllele" column represents the reference allele. Required.
#' @param repeats Number of iterations. Default:30.
#' @param thread Number of threads to use. NULL for non-parallel iterative optimization. Default:10.
#'
#' @examples
#' data("TF1_clones")
#' data=TF1_clones$data
#' ## performs a non-parallel iteration process
#' result=lineage(data, repeats=30, thread=NULL)
#' ## performs a parallel iteration process 
#' # result=lineage(data, repeats=30, thread=10)
#'
#' @return A Rdata with the following elements.
#' \item{best}{List, containing the best result in all repeats with max Sscore}
#' \item{label}{List, containing the name/id and inferredlineage label of each sample}
#' \item{suggest}{Character, recommended dimensional reduction method, tsne or umap}
#' \item{results}{List, containing results of all repeats with different 'centers' and 'num of markers'}
#' \item{scores}{Numeric, Sscores of all repeats}
#'
#' @export
#'

lineage=function(data,repeats=30,thread=10){
    # Split the variants frequency matrix according to mutation types
    d=data_prepare(data)

    # Provide alternative values for 'centers' and 'nmarker' according to the input data
    a=dim(d[[1]])[2]
    if(a>100){
        nmarker=c(15,20)
        centers=3
    }else{
        nmarker=c(10,15)
        centers=2
    }

    # perform lineage tracing analysis
    if (is.null(thread)){
        print("non-parallel processing")
        all_results=best_iter(d=d,repeats=repeats,nmarker=nmarker,centers=centers)
    }else{
        print("parallel processing")
        all_results=multi_best_iter(d=d,repeats=repeats,nmarker=nmarker,centers=centers,thread=thread)
    }
    save(all_results,file="all_results.Rdata")
    return(all_results)
}
