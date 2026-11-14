# GPL >= 2. Original BioGeoBEARS coo2crs validation; Rust compression kernel.
coo2crs <- function (ia, ja, a, n = NA, transpose_needed = FALSE) 
{
    defaults = "\n\t# Qmat from Psychotria (16x16)\n\tvals = c(0,0.01,0.01,0.01,0.01,0,0,0,0,0,0,0,0,0,0,0,0,-0.04,0,0,0,0.01,0.01,0.01,0,0,0,0,0,0,0,0,0,0,-0.04,0,0,0.01,0,0,0.01,0.01,0,0,0,0,0,0,0,0,0,-0.04,0,0,0.01,0,0.01,0,0.01,0,0,0,0,0,0,0,0,0,-0.04,0,0,0.01,0,0.01,0.01,0,0,0,0,0,0,0.01,0.01,0,0,-0.06,0,0,0,0,0,0.01,0.01,0,0,0,0,0.01,0,0.01,0,0,-0.06,0,0,0,0,0.01,0,0.01,0,0,0,0.01,0,0,0.01,0,0,-0.06,0,0,0,0,0.01,0.01,0,0,0,0,0.01,0.01,0,0,0,0,-0.06,0,0,0.01,0,0,0.01,0,0,0,0.01,0,0.01,0,0,0,0,-0.06,0,0,0.01,0,0.01,0,0,0,0,0.01,0.01,0,0,0,0,0,-0.06,0,0,0.01,0.01,0,0,0,0,0,0,0.02,0.02,0,0.02,0,0,-0.06,0,0,0,0.01,0,0,0,0,0,0.02,0,0.02,0,0.02,0,0,-0.06,0,0,0.01,0,0,0,0,0,0,0.02,0.02,0,0,0.02,0,0,-0.06,0,0.01,0,0,0,0,0,0,0,0,0.02,0.02,0.02,0,0,0,-0.06,0.01,0,0,0,0,0,0,0,0,0,0,0,0.03,0.03,0.03,0.03,-0.04)\n\t\n\tQmat = matrix(vals, nrow=16, ncol=16, byrow=FALSE)\n\t#library(SparseM)\n\t#library(kexpmv)\n\tcrsQmat = SparseM::as.matrix.csr(Qmat)\n\tcrsQmat\n\t\n\t#library(BioGeoBEARS)\n\tcooQmat = rexpokit::mat2coo(Qmat)\n\tcooQmat\n\t\n\t# Test coo2crs\n\tia=cooQmat[,\"ia\"]\n\tja=cooQmat[,\"ja\"]\n\ta=cooQmat[,\"a\"]\n\ttranspose_needed = FALSE\n\tn=16\n\tcrsR_Qmat = BioGeoBEARS::coo2crs(ia=ia, ja=ja, a=a, n=16, transpose_needed=FALSE)\n\tcrsR_Qmat\n\t\n\t# Check:\n\tcrsQmat@ia\n\tcrsR_Qmat$ia\n\tall.equal(crsR_Qmat$ia, crsQmat@ia)\n\t"
    if (any(length(ia) == 0, length(ja) == 0, length(a) == 0) == TRUE) {
        txt = paste0("WARNING from coo2crs: you have input a COO-formated matrix with 0 rows. This suggests you have a 1x1 transition matrix. If so, you really don't need to be doing sparse matrix exponentiation!  In BioGeoBEARS, change force_sparse to FALSE.")
        cat("\n")
        cat(txt)
        warning(txt)
    }
    if (is.na(n)) {
        n = max(max(ia), max(ja))
        txt = paste0("WARNING in coo2crs(). The input 'n', the assumed dimension of the transition matrix, was guessed from max(max(ia), max(ja)). This guess was n=", n, ". If this is wrong, input n manually.")
        cat("\n\n")
        cat(txt)
        cat("\n\n")
        warning(txt)
    }
    if (all.equal(length(ia), length(ja), length(a)) == FALSE) {
        txt = paste0("STOP ERROR in coo2crs(). The inputs 'ia', 'ja', and 'a' must all have the same length. Instead, you have length(ia)=", length(ia), ", length(ja)=", length(ja), ", length(a)=", length(a), ".")
        cat("\n\n")
        cat(txt)
        cat("\n\n")
        stop(txt)
    }
    {
        r <- .bgb_kernel(14, list(ia, ja, a, n, transpose_needed))
        names(r) <- c("ia", "ja", "a", "dimension")
        if (is.integer(if (transpose_needed) ia else ja)) 
            r$ja <- as.integer(r$ja)
        if (is.integer(a)) 
            r$a <- as.integer(r$a)
        r$dimension <- c(n, n)
        return(r)
    }
}

# DEmat's historical COO convention stores Q^T (ia=destination, ja=source).
# The Rust action accepts a column vector. Forward probability propagation
# therefore consumes that COO matrix directly. The old kexpmv adapter's extra
# swap cannot be carried over to a conventional column-vector action.
calc_prob_forward_onebranch_sparse <- function(relprobs_branch_bottom, branch_length, tmpQmat_in_REXPOKIT_coo_fmt, coo_n, anorm, check_for_0_rows=TRUE, TRANSPOSE_because_forward=TRUE, printlevel=1) {
    coo<-tmpQmat_in_REXPOKIT_coo_fmt
    if(printlevel>=1)cat('U')
    .bgb_kernel(18,list(coo_n,coo[,1],coo[,2],coo[,3],branch_length,relprobs_branch_bottom,!TRANSPOSE_because_forward))[[1]]
}
