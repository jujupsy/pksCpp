



data(DoignonFalmagne7)
context("as.pattern")
test_that("as.pattern DoignonFalmagne7", {                                
  expect_equal(
    as.pattern(DoignonFalmagne7$K), 
    c("00000", "10000", "01000", "11000", "11100", "11010", "11110", 
      "11101", "11111")
  )   
  expect_equal(
    as.pattern(DoignonFalmagne7$K, freq = TRUE),  
    c("00000" = 1, "01000" = 1, "10000" = 1, "11000" = 1, 
      "11010" = 1, "11100" = 1, "11101" = 1, "11110" = 1, "11111" = 1)
  )
  expect_equal(
    as.pattern(DoignonFalmagne7$K, as.letters = TRUE),  
    c("00000" = "0",    "10000" = "a",   "01000" = "b",   "11000" = "ab", 
      "11100" = "abc",  "11010" = "abd", "11110" = "abcd", 
      "11101" = "abce", "11111" = "abcde")
  )   
  expect_equal(
    as.pattern(DoignonFalmagne7$K, as.set = TRUE),  
    set(set(), set("a"), set("b"), set("a", "b"), set("a", "b", "c"), 
        set("a", "b", "d"), set("a", "b", "c", "d"), set("a", "b", "c", "e"), 
        set("a", "b", "c", "d", "e"))
  )   
})


data(endm)
test_that("as.pattern endm", {                                
  expect_equal(
    as.pattern(endm$K), 
    c("0000", "0110", "0101", "1110", "1101", "1011", "1111")
  )   
  expect_equal(
    as.pattern(endm$K, freq = TRUE),  
    c("0000" = 1L, "0101" = 1L, "0110" = 1L, "1011" = 1L, "1101" = 1L, 
      "1110" = 1L, "1111" = 1L)
  )
  expect_equal(
    as.pattern(endm$K, as.letters = TRUE),  
    c("0000" = "0", "0110" = "bc", "0101" = "bd", "1110" = "abc", 
      "1101" = "abd", "1011" = "acd", "1111" = "abcd")
  )   
  expect_equal(
    as.pattern(endm$K, as.set = TRUE),
    structure(list(structure(list(), class = c("set", "gset", "cset"
    )), structure(list("b", "c"), class = c("set", "gset", "cset"
    )), structure(list("b", "d"), class = c("set", "gset", "cset"
    )), structure(list("a", "b", "c"), class = c("set", "gset", "cset"
    )), structure(list("a", "b", "d"), class = c("set", "gset", "cset"
    )), structure(list("a", "c", "d"), class = c("set", "gset", "cset"
    )), structure(list("a", "b", "c", "d"), class = c("set", "gset", 
                           "cset"))), class = c("set", "gset", "cset"))
 )   
})


data(DoignonFalmagne7)
dim(as.binmat(DoignonFalmagne7$N.R))
dim(as.binmat(DoignonFalmagne7$N.R, uniq = FALSE))

## Knowledge structure as binary matrix
as.binmat(c("000", "100", "101", "111"))
as.binmat(set(set(), set("a"), set("a", "c"), set("a", "b", "c")))


data(DoignonFalmagne7)
context("as.binmat")
test_that("as.binmat DoignonFalmagne7", {                                
  expect_equal(
    as.binmat(c("000", "100", "101", "111")), 
    structure(c(0L, 1L, 1L, 1L, 0L, 0L, 0L, 1L, 0L, 0L, 1L, 1L), .Dim = 4:3, .Dimnames = list(
      NULL, c("a", "b", "c")))
  )   
  expect_equal(
    as.binmat(set(set(), set("a"), set("a", "c"), set("a", "b", "c"))),  
    structure(c(0L, 1L, 1L, 1L, 0L, 0L, 0L, 1L, 0L, 0L, 1L, 1L), .Dim = 4:3, .Dimnames = list(
      NULL, c("a", "b", "c")))
  )
  expect_equal(
    as.binmat(set(set(), set("a"), set("a", "c"), set("a", "b", "c"))),  
    as.binmat(c("000", "100", "101", "111"))
  )
  expect_equal(
    as.binmat(c("000", "100", "101", "111"), col.names = c("r", "q", "jk")),  
    structure(c(0L, 1L, 1L, 1L, 0L, 0L, 0L, 1L, 0L, 0L, 1L, 1L), .Dim = 4:3, .Dimnames = list(
      NULL, c("r", "q", "jk")))
  )   
  expect_equal(
    as.binmat(c("000" = 10, "100" = 5, "101" = 3, "111" = 1), col.names = c("r", "q", "jk"), uniq = TRUE),  
    structure(c(0L, 1L, 1L, 1L, 0L, 0L, 0L, 1L, 0L, 0L, 1L, 1L), .Dim = 4:3, .Dimnames = list(
      NULL, c("r", "q", "jk")))
  )
  expect_equal(
    as.binmat(c("000" = 10, "100" = 5, "101" = 3, "111" = 1), col.names = c("r", "q", "jk"), uniq = FALSE),  
    structure(c(0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 1L, 1L, 
                1L, 1L, 1L, 1L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
                0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
                0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L), .Dim = c(19L, 
                  3L), .Dimnames = list(NULL, c("r", "q", "jk")))
  ) 
})


