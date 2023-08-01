library(ape)



#given the taxon tree this returns a list with : 
# [[1]] data frame - pthr taxon parent group, child group, 5 letter name of species if existing, taxon id if existing
# [[2]] vector - number of children I believe
recursiveTreeParse = function(tree){
    res = data.frame()
    res2 = c()
    recurseInner = function(tree){
        #get the name of current node
        n = tree$name
        n2 = n
        if(!is.null(tree$short_name)) n = str_c(n,'---',tree$short_name)
        else n = str_c(n,'---')
        if(!is.null(tree$taxon_id)) n = str_c(n,'@@@',tree$taxon_id)
        else n = str_c(n,'@@@')
        #get the list of direct childrens names
        childList = c()
        if('children' %in% names(tree)) {
            for(i in 1:length(tree$children$node)) {
                childtree = tree$children$node[[i]]
                childList = c(childList,recurseInner(childtree))
            }
        }
        #add each combo to res
        res2 <<- c(res2,length(childList))
        if(length(childList) > 0) {
            for(i in 1:length(childList)){
                res <<- rbind(res,c(n2,childList[i]))
            }
        }
        return(n)
    }
    recurseInner(tree)
    res3 = cbind(res[,1],str_split(res[,2],'---|@@@',simplify=T))
    return(list(res3,res2))
} 

 

#convert panther tree files to phylo objects
read_panther <- function(x, tree.reader = ape::read.tree, ...) {
  # Reading the data-in
  x  <- readLines(x)
  x[1] = gsub('SAR/','SAR-',x[1])
  x[1] = gsub('Ev=UNK','Ev=0>0',x[1]) #Neither of these two are recognized I think
  # Obtaining extra info and processing internal nodes labels
  rgxp <- "(?:[:])?([0-9.]+)?\\[\\&\\&NHX:Ev=([0-9><]{3})(?::S=([a-zA-Z_.-]*))?:ID=([a-zA-Z0-9]+)\\]"
  
  dat <- gregexpr(rgxp, x[1], perl = TRUE)
  dat <- regmatches(x[1], dat)
  dat <- do.call(rbind, regmatches(dat[[1]], regexec(rgxp, dat[[1]], perl = TRUE)))
  dat[dat == ""] <- NA_character_
  
  # Rewriting the file so that labels of inner nodes can be read in
  for (i in 1:nrow(dat))
    x[1] <- gsub(
      x           = x[1],
      pattern     = dat[i,1],
      fixed       = TRUE,
      replacement = ifelse(is.na(dat[i,2]), dat[i,5], paste(dat[i,5],dat[i,2], sep=":"))
      )
  
  # Getting the labels
  labs <- data.frame(
    id    = unlist(regmatches(x[-1], regexec(pattern = "^.+(?=\\:)", x[-1], perl = TRUE))),
    label = unlist(regmatches(x[-1], regexec(pattern = "(?<=\\:).+(?=\\;$)", x[-1], perl = TRUE))),
    stringsAsFactors = FALSE
  )
  
  # Reading the tree
  # tree <- ape::read.tree(text=x[1], ...)
  tmptree <- tempfile()
  write(x[1], tmptree)
  tree <- tree.reader(tmptree, ...)
  file.remove(tmptree)
  
  if (Nnode(tree) != nrow(dat))
    stop(
      "The number of nodes read does not coincide with that reported by the ",
      "tree.reader. This could be an updated version of PantherDB.",
      call. = FALSE
      )
  
  tree$tip.label <- paste(
    tree$tip.label,
    labs$label[match(tree$tip.label, labs$id)],
    sep=":"
  )
  
  # Creating a nice data-frame
  ans <- data.frame(
    branch_length    = as.numeric(dat[,2]),
    type             = ifelse(dat[,3] == "0>1", "S",
                              ifelse(dat[,3] == "1>0", "D", "T")),
    ancestor         = dat[,4],
    row.names        = dat[,5], 
    stringsAsFactors = FALSE
  )
  
  # Which ones are duplication nodes
  ans$duplication <- ifelse(ans$type %in% c("D", "T"), TRUE, FALSE)
  
  # Sorting and returning
  list(
    tree = tree,
    internal_nodes_annotations  = ans[order(as.integer(gsub("[a-zA-Z]+","",rownames(ans)))),]
  )
}

#' @rdname panther-tree
#' @export
read.panther <- read_panther

