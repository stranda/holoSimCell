

##
## minor allele freq calcs
##
test_that("mafreq is working",
{
    expect_equal(length(mafreq(out)), length(names(out@data))-2) #ensures the right data are being produced
    expect_equal(unname(mafreq(out)),unname(apply(out@data[,-1:-2],2,function(x) {tb=table(x);sort(tb/sum(tb))[1]}))) #compares two ways to calc mafreq
    expect_equal_to_reference(mafreq(out),"mafreq-out.rds")# compares to previous run on same data
})

##
####strataNames  tests
###
test_that("population strata names are coming out of strataNames correctly",
{
    expect_match(strataNames(out), "^[a-zA-Z]") #make sure the first character of each strata is char not num
    expect_equal(as.character(sort(unique(out@data[[2]]))),as.character(strataNames(out))) #make sure it is extracting from the second column of out@data
})

test_that("pairwise Fst calculation is working",
{
###set up for the call to pwise.fst.glob 
    allMAF <- mafreq(out)
    split_out <- strataSplit(out) #list of strataG objects for each pop
    locMAF <- do.call(cbind,lapply(split_out,function(o){mafreq(o)}))
    colnames(locMAF) <- names(split_out)
    locHe <- colMeans(2*locMAF*(1-locMAF))
    varlocHe <- apply(2*locMAF*(1-locMAF),2,var)
    locN <- sapply(split_out,function(o){length(o@data$ids)})
    names(locN) <- strataNames(out)
###ok data for testing pwise.fst.glob are ready (locMAF, allMAF, locN)

    pw.out = pwise.fst.glob(locMAF,allMAF,locN)
    diag(pw.out)=NA
    expect_equal_to_reference(pw.out,"pwise-fst-glob-out.rds")# compares to previous run on same data

###now compare the output of this function to another calculator
    pairwiseTest(out,stats="gst",nrep=NULL,quietly=T) -> pw.strataG
    ###check that the correlation between strataG gst calc and our version is high (>0.96)
    expect_gt(cor(c(pw.out)[!is.na(c(pw.out))],c(pw.strataG$pair.mat[[1]])[!is.na(c(pw.strataG$pair.mat[[1]]))]) , 0.96)
    
})
