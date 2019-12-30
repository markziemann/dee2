library("getDEE2")
library("testthat")

test_that("multiplication works", {
    expect_equal(2 * 2, 4)
})


# E. coli 
x<-getDEE2("ecoli",c("SRR1613487","SRR1613488"))

test_that("eco works", {
    expect_equal( sum(x$GeneCounts) , 20624168 )
})


# C. elegans
x<-getDEE2("celegans",c("SRR051935","SRR051934"))

test_that("cel works", {
    expect_equal( sum(x$TxCounts) , 29067 )
})


# M. musculus
x<-getDEE2("mmusculus",c("SRR1022283","SRR1022284"))
x<-Tx2Gene(x)

test_that("mmu works", {
    expect_equal( sum(x$GeneCounts) , 14378076 )
    expect_equal( sum(x$TxCounts) , 2872578 )
    expect_equal( sum(x$Tx2Gene) , 2872578 )
    expect_equal( nrow(x$Tx2Gene) , 34839 )
})



