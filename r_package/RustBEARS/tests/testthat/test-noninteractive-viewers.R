# Render real PDFs in Rscript while observing only the external-viewer boundary.
test_that("BSM diagnostics render PDFs without launching viewers in Rscript", {
    exercise <- function() {
        folder <- tempfile("bgb-viewer-check-"); dir.create(folder)
        marker <- file.path(folder,"viewer-invoked.txt")
        old_tests <- Sys.getenv("R_TESTS",unset=NA_character_)
        on.exit({
            if (is.na(old_tests)) Sys.unsetenv("R_TESTS") else Sys.setenv(R_TESTS=old_tests)
            unlink(folder,recursive=TRUE)
        },add=TRUE)
        Sys.setenv(R_TESTS="")
        # Absolute Rscript invocation avoids Windows SearchPath resolving a
        # different 'open' executable before a test's PATH entry.
        writeLines(c("args <- commandArgs(trailingOnly=TRUE)",
                     "writeLines('viewer invoked',args[[1]])",
                     "quit(save='no',status=99L,runLast=FALSE)"),
                   file.path(folder,"sentinel.R"))
        script_body <- quote({
            suppressPackageStartupMessages(library(RustBEARS))
            stopifnot(!base::interactive())
            marker <- file.path(getwd(),"viewer-invoked.txt")
            sentinel <- file.path(getwd(),"sentinel.R")
            rscript <- file.path(R.home("bin"),"Rscript")
            allow_viewer <- FALSE
            viewer_boundary <- function(command, ...) {
                stopifnot(grepl("^open[[:space:]]",command))
                status <- base::system2(rscript,c("--vanilla",shQuote(sentinel),shQuote(marker)),
                                        stdout="sentinel.log",stderr="sentinel.log",timeout=10)
                if (!allow_viewer) stop("unexpected external viewer invocation",call.=FALSE)
                status
            }
            clone_plot <- function(name) {
                original <- getExportedValue("RustBEARS",name)
                cloned <- original
                environment(cloned) <- list2env(list(system=viewer_boundary,
                    interactive=base::interactive),parent=environment(original))
                stopifnot(identical(body(cloned),body(original)),
                          identical(formals(cloned),formals(original)),
                          !identical(environment(cloned),environment(original)))
                list(fn=cloned,original=original,namespace_env=environment(original),name=name)
            }
            hist <- clone_plot("hist_event_counts")
            comparison <- clone_plot("check_ML_vs_BSM")
            counts <- list(founder_totals_list=c(0,1,0,1),ana_totals_list=c(1,2,1,3),
                subsetSymp_totals_list=c(1,0,1,0),vicariance_totals_list=c(0,1,0,1),
                sympatry_totals_list=c(1,1,1,1),all_totals_list=c(3,5,3,6))
            tr <- ape::read.tree(text="((a:1,b:1):1,c:2);")
            res <- list(inputs=list(stratified=FALSE),
                ML_marginal_prob_each_state_at_branch_top_AT_node=
                    matrix(c(.7,.3,.2,.8,.4,.6,.6,.4,.3,.7),ncol=2,byrow=TRUE))
            maps <- lapply(1:4,function(i) data.frame(node=1:5,
                node.type=c(rep("tip",3),"root","internal"),
                sampled_states_AT_nodes=c(1,2,1,c(1,1,2,1)[i],c(2,2,1,2)[i]),
                sampled_states_AT_brbots=c(1,2,1,1,2),clado_event_txt=rep("sympatry",5)))
            plots <- list(histogram=list(clone=hist,draw=function()
                hist$fn(counts,pdffn="histogram.pdf")),
                comparison=list(clone=comparison,draw=function()
                comparison$fn(res,maps,model_name="comparison",tr=tr,MultinomialCI=FALSE)))
            positive <- list()
            # Each unchanged production function must reach the real sentinel
            # process under a controlled interactive branch, proving interception.
            allow_viewer <- TRUE
            for (name in names(plots)) {
                plot <- plots[[name]]
                environment(plot$clone$fn)$interactive <- function() TRUE
                status <- plot$draw()
                positive[[name]] <- list(status=status,marker=file.exists(marker))
                stopifnot(identical(status,99L),file.exists(marker))
                unlink(marker)
                environment(plot$clone$fn)$interactive <- base::interactive
            }
            saveRDS(positive,"positive-controls.rds")
            pdfs <- c("histogram.pdf","comparison_ML_vs_BSM.pdf")
            unlink(pdfs)
            stopifnot(!any(file.exists(pdfs)))
            # The actual Rscript condition must create fresh PDFs with those
            # same function bodies; positive-control output cannot satisfy it.
            allow_viewer <- FALSE
            for (plot in plots) {
                stopifnot(!environment(plot$clone$fn)$interactive())
                plot$draw()
                stopifnot(!file.exists(marker))
                installed <- getExportedValue("RustBEARS",plot$clone$name)
                stopifnot(identical(body(installed),body(plot$clone$original)),
                          identical(formals(installed),formals(plot$clone$original)),
                          identical(environment(installed),plot$clone$namespace_env))
            }
            for (p in pdfs) {
                bytes <- readBin(p,"raw",n=file.info(p)$size)
                stopifnot(length(bytes)>100,rawToChar(bytes[1:4])=="%PDF",
                          grepl("%%EOF",rawToChar(tail(bytes,20)),fixed=TRUE))
            }
            cat("NONINTERACTIVE_VIEWER_TEST_PASSED\n")
        })
        script <- file.path(folder,"render.R")
        writeLines(c(paste0(".libPaths(",paste(deparse(.libPaths()),collapse=""),")"),
                     paste0("setwd(",deparse(folder),")"),deparse(script_body)),script)
        log <- file.path(folder,"render.log")
        status <- system2(file.path(R.home("bin"),"Rscript"),c("--vanilla",shQuote(script)),
                          stdout=log,stderr=log,timeout=30)
        text <- readLines(log,warn=FALSE)
        expect_identical(status,0L,info=paste(text,collapse="\n"))
        positive <- readRDS(file.path(folder,"positive-controls.rds"))
        for (name in c("histogram","comparison")) {
            expect_identical(positive[[name]]$status,99L,info=name)
            expect_true(positive[[name]]$marker,info=name)
        }
        expect_true(any(grepl("NONINTERACTIVE_VIEWER_TEST_PASSED",text,fixed=TRUE)))
        expect_false(file.exists(marker))
        expect_true(all(file.info(file.path(folder,c("histogram.pdf","comparison_ML_vs_BSM.pdf")))$size>100))
    }
    exercise()
})
