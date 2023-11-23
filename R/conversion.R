#' Conversion between Representations of Responses or States
#' 
#' @description
#' Converts between binary matrix and pattern representations of response
#' patterns or knowledge states.
#'
#' @usage
#' as.pattern(R, freq = FALSE, as.letters = FALSE, as.set = FALSE)
#' as.binmat(N.R, uniq = TRUE, col.names = NULL)
#' 
#' @param R an indicator matrix of response patterns or knowledge states.
#' @param N.R either a (named) vector of absolute frequencies of response
#'            patterns; or a character vector of response patterns or 
#'            knowledge states;or a \code{set} of sets representing the 
#'            knowledge structure.
#' @param freq logical, should the frequencies of response patterns be
#'             reported?
#' @param uniq logical, if \code{TRUE}, only the unique response patterns are
#'             returned.
#' @param as.letters logical, return response patterns as combinations of letters.
#' @param as.set logical, return response patterns as set of sets.
#' @param col.names column names for the state or response matrix.
#'
#' @details 
#'   Functions \code{as.pattern} and \code{as.binmat} were adapted from the 
#'   package \code{pks} to deal with missing data.
#'   
#'   Missings in the binmat have to be coded as \code{NA} and are  represented in the
#'   pattern as \code{M}. 
#'   \code{as.letters} or \code{as.set} ignore NA values.
#' 
#' @return
#'   \code{as.pattern} returns a vector of integers named by the response
#'   patterns if \code{freq} is \code{TRUE}, else a character vector. If
#'   \code{as.set} is \code{TRUE}, the return value is of class \code{set}.
#'   
#'   \code{as.binmat} returns an indicator matrix.
#'
#' @aliases 
#'  as.pattern
#'  as.binmat
#'  
#' @name Conversion
#' @examples 
#' data(DoignonFalmagne7)
#' as.pattern(DoignonFalmagne7$K)
#' as.pattern(DoignonFalmagne7$K, freq = TRUE)
#' as.pattern(DoignonFalmagne7$K, as.letters = TRUE)
#' as.pattern(DoignonFalmagne7$K, as.set = TRUE)
#' 
#' dim(as.binmat(DoignonFalmagne7$N.R))
#' dim(as.binmat(DoignonFalmagne7$N.R, uniq = FALSE))
#' 
#' ## Knowledge structure as binary matrix
#' as.binmat(c("000", "100", "101", "111"))
#' as.binmat(set(set(), set("a"), set("a", "c"), set("a", "b", "c")))
#' 
#' ## deal with missings
NULL



# Convert binary matrix including missing values (NA) to vector of response patterns 
# adapted from pks package

#' @export
as.pattern <- function(R, freq = FALSE, as.letters = FALSE, as.set = FALSE){
  
  #handle  missings
  missings <- any(is.na(R)) # are there any missings?
  R[is.na(R)] <- "M"
  
  if(freq){
    N.R <- table(apply(R, 1, paste, collapse=""))
    setNames(as.integer(N.R), names(N.R))          # convert to named int
  }else
    if(as.letters | as.set){
      if(missings){
        warning("Caveat: as.letters or as.set ignore NA values!")
      }
      nitems <- ncol(R)
      item.names <- 
        make.unique(c("a", letters[(seq_len(nitems) %% 26) + 1])[-(nitems + 1)],
                    sep="")
      lett <- apply(R, 1, function(r) paste(item.names[which(r == 1)],
                                            collapse=""))
      lett[lett == ""] <- "0"
      
      if(as.set){
        # Separate elements in lett by "_", remove leading "_",
        # then strsplit along "_" (trailing "_" are ignored by strsplit)
        setfam <- as.set(lapply(strsplit(
          gsub("^_(.+)", "\\1", gsub("([0-9]*)", "\\1_", unname(lett))),
          "_"), as.set))
        if (set_contains_element(setfam, set("0")))
          setfam[[set("0")]] <- set()  # proper empty set
        setfam  # return family of sets, class set
      }else
        lett    # return letters, class character
    }else
      unname(apply(R, 1, paste, collapse=""))
}


# Convert vector of response patterns including missing values (M) to named 
# binary matrix adapted from pks package

#' @export
as.binmat <- function(N.R, uniq = TRUE, col.names = NULL){
  
  if (is.set(N.R)) {
    states <- lapply(N.R, as.character)
    items <- sort(unique(unlist(states)))
    R <- matrix(0, length(N.R), length(items),
                dimnames=list(NULL,
                              if(is.null(col.names)) items else col.names))
    for (i in seq_len(nrow(R))) R[i, states[[i]]] <- 1
  } else {
    pat <- if(is.null(names(N.R))) N.R else names(N.R)
    R   <- if(uniq) strsplit(pat, "") else strsplit(rep(pat, N.R), "")
    R   <- do.call(rbind, R)
    
    colnames(R) <- 
      if(is.null(col.names)){
        nitems <- ncol(R)
        make.unique(c("a", letters[(seq_len(nitems) %% 26) + 1])[-(nitems + 1)],
                    sep="")
      }else
        col.names
  }
  R[R == "M"] <- NA 
  storage.mode(R) <- "integer"
  R
}
