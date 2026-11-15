# Real noninteractive subprocesses must render output without opening a viewer.
test_that("BSM diagnostics render PDFs without launching viewers in Rscript", {
    exercise <- function() {
        folder <- tempfile("bgb-viewer-check-"); dir.create(folder)
        bin <- file.path(folder, "bin"); dir.create(bin)
        marker <- file.path(folder, "viewer-invoked.txt")
        old <- Sys.getenv(c("PATH", "RUSTBEARS_VIEWER_MARKER", "R_TESTS"), unset=NA_character_)
        on.exit({
            for (key in names(old)) {
                if (is.na(old[[key]])) Sys.unsetenv(key)
                else do.call(Sys.setenv, setNames(list(old[[key]]),key))
            }
            unlink(folder,recursive=TRUE)
        }, add=TRUE)
        # A native stub works with Windows CreateProcess as well as POSIX shells.
        # It records every invocation and returns an unmistakable failure code.
        cfile <- file.path(folder,"viewer.c")
        writeLines(c('#include <stdio.h>', '#include <stdlib.h>',
            'int main(void) { const char *p = getenv("RUSTBEARS_VIEWER_MARKER");',
            'FILE *f = p ? fopen(p, "ab") : NULL;',
            'if (f) { fputs("viewer invoked\\n", f); fclose(f); } return 99; }'), cfile)
        compiler <- Sys.which(if (.Platform$OS.type == "windows") "gcc" else "cc")
        expect_true(nzchar(compiler),info="a source-package check requires its C compiler")
        if (!nzchar(compiler)) stop("C compiler unavailable for the viewer process regression")
        executable <- file.path(bin, if (.Platform$OS.type == "windows") "open.exe" else "open")
        build_log <- file.path(folder,"build.log")
        status <- system2(compiler,c(shQuote(cfile),"-o",shQuote(executable)),
                          stdout=build_log,stderr=build_log,timeout=30)
        expect_identical(status,0L,info=paste(readLines(build_log,warn=FALSE),collapse="\n"))
        if (status != 0L) stop("viewer stub compilation failed")
        Sys.setenv(PATH=paste(bin,old[["PATH"]],sep=.Platform$path.sep),
                   RUSTBEARS_VIEWER_MARKER=marker,R_TESTS="")
        # Positive control proves that the actual 'open' command is intercepted.
        status <- suppressWarnings(system("open viewer-positive-control",ignore.stdout=TRUE,
                                          ignore.stderr=TRUE,timeout=10))
        expect_identical(status,99L)
        expect_true(file.exists(marker))
        unlink(marker)
        script <- file.path(folder,"render.R")
        writeLines(c(
            paste0(".libPaths(",paste(deparse(.libPaths()),collapse=""),")"),
            "suppressPackageStartupMessages(library(RustBEARS))",
            paste0("setwd(",deparse(folder),")"),
            "stopifnot(!interactive())",
            "counts <- list(founder_totals_list=c(0,1,0,1),ana_totals_list=c(1,2,1,3),",
            " subsetSymp_totals_list=c(1,0,1,0),vicariance_totals_list=c(0,1,0,1),",
            " sympatry_totals_list=c(1,1,1,1),all_totals_list=c(3,5,3,6))",
            "hist_event_counts(counts,pdffn='histogram.pdf')",
            "stopifnot(!file.exists(Sys.getenv('RUSTBEARS_VIEWER_MARKER')))",
            "tr <- ape::read.tree(text='((a:1,b:1):1,c:2);')",
            "res <- list(inputs=list(stratified=FALSE),",
            " ML_marginal_prob_each_state_at_branch_top_AT_node=matrix(c(.7,.3,.2,.8,.4,.6,.6,.4,.3,.7),ncol=2,byrow=TRUE))",
            "maps <- lapply(1:4,function(i) data.frame(node=1:5,node.type=c(rep('tip',3),'root','internal'),",
            " sampled_states_AT_nodes=c(1,2,1,c(1,1,2,1)[i],c(2,2,1,2)[i]),",
            " sampled_states_AT_brbots=c(1,2,1,1,2),clado_event_txt=rep('sympatry',5)))",
            "check_ML_vs_BSM(res,maps,model_name='comparison',tr=tr,MultinomialCI=FALSE)",
            "stopifnot(!file.exists(Sys.getenv('RUSTBEARS_VIEWER_MARKER')))",
            "for (p in c('histogram.pdf','comparison_ML_vs_BSM.pdf')) {",
            " bytes <- readBin(p,'raw',n=file.info(p)$size)",
            " stopifnot(length(bytes)>100,rawToChar(bytes[1:4])=='%PDF',",
            "           grepl('%%EOF',rawToChar(tail(bytes,20)),fixed=TRUE))",
            "}",
            "cat('NONINTERACTIVE_VIEWER_TEST_PASSED\\n')"),script)
        log <- file.path(folder,"render.log")
        status <- system2(file.path(R.home("bin"),"Rscript"),c("--vanilla",shQuote(script)),
                          stdout=log,stderr=log,timeout=30)
        text <- readLines(log,warn=FALSE)
        expect_identical(status,0L,info=paste(text,collapse="\n"))
        expect_true(any(grepl("NONINTERACTIVE_VIEWER_TEST_PASSED",text,fixed=TRUE)))
        expect_false(file.exists(marker))
        expect_true(all(file.info(file.path(folder,c("histogram.pdf","comparison_ML_vs_BSM.pdf")))$size>100))
    }
    exercise()
})
