reportTGIF <-
function (result = NA, Comp.Matrix = NA, nDim = 10, label = NA, 
    out.dir = NA, html.open = FALSE, title = "Result of scTGIF", 
    author = "") 
{
    if (is.na(out.dir)) {
        stop("Please specify the output directory for saving your analysis result")
    }
    else {
        if (!file.exists(out.dir)) {
            dir.create(out.dir)
        }
    }
    Comp.Matrix <- data.frame(Comp.Matrix)
    if (class(result) != "scTGIF") {
        stop("Please specify the output object of calcTGIF!\n")
    }
    if (ncol(Comp.Matrix) <= 1) {
        stop("Larger number of dimention is specified in Comp.Matrix!\n")
    }
    if (result$NoCell != nrow(Comp.Matrix)) {
        stop("No. of Cell used in calcTGIF and that of Comp.Matrix are different!\n")
    }
    if (ncol(Comp.Matrix) == 2) {
        colnames(Comp.Matrix) <- paste0("Dim", 1:2)
        us <- lapply(1:nDim, function(x) {
            eval(parse(text = paste0("u", x, " <<- plot_ly(Comp.Matrix, x=~Dim1, y=~Dim2, text=result$NameCell, type='scatter', mode='markers', marker=list(color=smoothPalette(result$u[,", 
                x, "], palfunc=colorRampPalette(c('blue', 'gray', 'red')))), opacity=0.5)")))
        })
    }
    if (ncol(Comp.Matrix) == 3) {
        colnames(Comp.Matrix) <- paste0("Dim", 1:3)
        us <- lapply(1:nDim, function(x) {
            eval(parse(text = paste0("u", x, " <<- plot_ly(Comp.Matrix[,1:3], x=~Dim1, y=~Dim2, z=~Dim3, text=result$NameCell, type='scatter3d', mode='markers', marker=list(color=smoothPalette(result$u[,", 
                x, "], palfunc=colorRampPalette(c('blue', 'gray', 'red')))), opacity=0.5)")))
        })
    }
    if (ncol(Comp.Matrix) >= 4) {
        us <- lapply(1:nDim, function(x) {
            eval(parse(text = paste0("q <- sapply(seq(0,1,0.01), function(xx){quantile(result$u[, ", 
                x, "], xx)});group.CM <- sapply(result$u[, ", 
                x, "], function(xx){max(which(xx >= q))});col.CM <- smoothPalette(q, palfunc=colorRampPalette(c('blue', 'gray', 'red')));u", 
                x, " <- pairsD3(Comp.Matrix, width=150*nDim, big=TRUE, group=group.CM, col=col.CM)")))
        })
    }
    vs <- lapply(1:nDim, function(x) {
        eval(parse(text = paste0("v", x, " <<- plot_ly(x=result$NameTerm, y=result$v[,", 
            x, "], marker = list(color = smoothPalette(result$d[", 
            x, "] * result$v[,", x, "], palfunc=colorRampPalette(c('blue', 'gray', 'red')))), type='bar', opacity=0.5) %>% layout(yaxis=list(range = c(-1, 1)))")))
    })
    temp <- tempdir()
    sapply(1:nDim, function(x) {
        saveWidget(widget = us[[x]], file = paste0(temp, "/u", 
            x, ".html"), selfcontained = TRUE, background = "gray")
        saveWidget(widget = vs[[x]], file = paste0(temp, "/v", 
            x, ".html"), selfcontained = TRUE, background = "gray")
    })
    sapply(1:nDim, function(x) {
        file.copy(from = paste0(temp, "/u", x, ".html"), to = paste0(out.dir, 
            "/"), recursive = TRUE, overwrite = TRUE)
        file.copy(from = paste0(temp, "/v", x, ".html"), to = paste0(out.dir, 
            "/"), recursive = TRUE, overwrite = TRUE)
        file.copy(from = paste0(temp, "/u", x, "_files"), to = paste0(out.dir, 
            "/"), recursive = TRUE, overwrite = TRUE)
        file.copy(from = paste0(temp, "/v", x, "_files"), to = paste0(out.dir, 
            "/"), recursive = TRUE, overwrite = TRUE)
    })
    save(Comp.Matrix, file = paste0(temp, "/Comp.Matrix.RData"))
    save(label, file = paste0(temp, "/Group.RData"))
    file.copy(from = paste0(temp, "/Comp.Matrix.RData"), to = out.dir, 
        overwrite = TRUE)
    file.copy(from = paste0(temp, "/Group.RData"), to = out.dir, 
        overwrite = TRUE)
    png(filename = paste0(temp, "/EigenValue.png"))
    barplot(result$d^2, main = "Eigen Value Distribution")
    dev.off()
    file.copy(from = paste0(temp, "/EigenValue.png"), to = out.dir, 
        overwrite = TRUE)
    png(filename = paste0(temp, "/CumEigenValue.png"))
    barplot(cumsum(result$d^2)/sum(result$d^2), main = "Cumulative Eigen Value Distribution")
    dev.off()
    file.copy(from = paste0(temp, "/CumEigenValue.png"), to = out.dir, 
        overwrite = TRUE)
    HEADER <- "---\ntitle: XXXXX\nauthor: YYYYY\ndate: \"`r Sys.time()`\"\noutput: BiocStyle::html_document\nvignette: >\n %\\VignetteIndexEntry{Vignette Title}\n %\\VignetteEngine{knitr::rmarkdown}\n %\\VignetteEncoding{UTF-8}\n---\n"
    HEADER[grep("XXXXX", HEADER)] <- sub("XXXXX", title, HEADER[grep("XXXXX", 
        HEADER)])
    HEADER[grep("YYYYY", HEADER)] <- sub("YYYYY", author, HEADER[grep("YYYYY", 
        HEADER)])
    BODY1 <- "# Dataset details<U+FF08>X<U+FF09>\n<div style=\"width:100px\">\n![](Matrix1.tiff)\n</div>\n**AAAAA** Gene <U+00D7> **BBBBB** Cell Matrix\n\n# Gene set details<U+FF08>Y<U+FF09>\n\n<div style=\"width:100px\">\n![](Matrix2.tiff)\n</div>\n**CCCCC** Gene <U+00D7> **DDDDD** Term Matrix\n\n# About the Algorithm of scTGIF\n\n![](Algorithm.tiff)\n\nTransposed matrix $X$ are multiplied by matrix $Y$, and then Singular Value Decomposition (SVD) is performed for compressing the matrix to specified rank ($K$=**EEEEE**). Matrix $U$ contains some eigen vector $u_i$ ($1 \\leq i \\leq K$), which represent cellular pattern. Matrix $V$ contains some eigen vector $v_i$ ($1 \\leq i \\leq K$), which represent fuctional pattern. This is a kind of Partial Least Square (PLS). $\\Sigma^{1/2}$ is the diagonal matrix, which contains the squared eigen value in each diagonal element. Distribution of eigen value is bellow.\n\n# Global statistics and plots\n<div style=\"width:300px\">\n![](EigenValue.png)![](CumEigenValue.png)\n</div>\n"
    BODY1[grep("AAAAA", BODY1)] <- sub("AAAAA", result$NoGene, 
        BODY1[grep("AAAAA", BODY1)])
    BODY1[grep("BBBBB", BODY1)] <- sub("BBBBB", result$NoCell, 
        BODY1[grep("BBBBB", BODY1)])
    BODY1[grep("CCCCC", BODY1)] <- sub("CCCCC", result$NoGene, 
        BODY1[grep("CCCCC", BODY1)])
    BODY1[grep("DDDDD", BODY1)] <- sub("DDDDD", result$NoTerm, 
        BODY1[grep("DDDDD", BODY1)])
    BODY1[grep("EEEEE", BODY1)] <- sub("EEEEE", nDim, BODY1[grep("EEEEE", 
        BODY1)])
    BODY2 <- 1:nDim
    for (i in 1:nDim) {
        BODY2[i] <- paste0("## Spectral", i, " (", round((result$d^2)[i]/sum(result$d^2) * 
            100, 3), " % variance)\n\nCellular Pattern : [u", 
            i, "](u", i, ".html)\n\nFunctional Pattern : [v", 
            i, "](v", i, ".html)\n")
    }
    BODY2 <- paste(BODY2, collapse = "\n")
    BODY3 <- "# Interactive Cell Browsing\nRun the following R script. This script needs [pairsD3](https://cran.r-project.org/web/packages/pairsD3/index.html) package.\n\n```{r eval=FALSE}\nlibrary(\"pairsD3\")\n\nload(\"/PATH/TO/DIR/YOU/SPECIFIED/Comp.Matrix.RData\")\nload(\"/PATH/TO/DIR/YOU/SPECIFIED/Group.RData\")\nsp <- shinypairs(Comp.Matrix, group=label, labels=rownames(Comp.Matrix))\nsp\n```\n![](1.tiff)\n\nIf you want to get the list of cells surrounded by the cursor in shinypairs window, select **Table of data?** as **Yes** and **Include all variables in table?** as **Yes**.\n![](2.tiff)\n\nThis interactive interface can help you list the ID of cells you are interested in.\n![](3.tiff)\n\n\n\n# License\nCopyright (c) 2017 Koki Tsuyuzaki and RIKEN Bioinformatics Research Unit Released under the [MIT license](https://opensource.org/licenses/mit-license.php)\n"
    sink(file = paste0(temp, "/index.Rmd"))
    cat(paste0(HEADER), "\n")
    cat(paste0(BODY1), "\n")
    cat(paste0(BODY2), "\n")
    cat(paste0(BODY3), "\n")
    sink()
    file.copy(from = paste0(temp, "/index.Rmd"), to = paste0(out.dir, 
        "/index.Rmd"), overwrite = TRUE)
    file.copy(from = system.file("extdata", "Matrix2.tiff", package = "scTGIF"), 
        to = paste0(out.dir, "/Matrix2.tiff"), overwrite = TRUE)
    file.copy(from = system.file("extdata", "Matrix1.tiff", package = "scTGIF"), 
        to = paste0(out.dir, "/Matrix1.tiff"), overwrite = TRUE)
    file.copy(from = system.file("extdata", "Algorithm.tiff", 
        package = "scTGIF"), to = paste0(out.dir, "/Algorithm.tiff"), 
        overwrite = TRUE)
    file.copy(from = system.file("extdata", "3.tiff", package = "scTGIF"), 
        to = paste0(out.dir, "/3.tiff"), overwrite = TRUE)
    file.copy(from = system.file("extdata", "2.tiff", package = "scTGIF"), 
        to = paste0(out.dir, "/2.tiff"), overwrite = TRUE)
    file.copy(from = system.file("extdata", "1.tiff", package = "scTGIF"), 
        to = paste0(out.dir, "/1.tiff"), overwrite = TRUE)
    e <- try(rmarkdown::render(paste0(out.dir, "/index.Rmd")))
    if (class(e) == "try-error") {
        e <- try(rmarkdown::render(paste0(out.dir, "/index.Rmd")))
        if (class(e) == "try-error") {
            e <- try(rmarkdown::render(paste0(out.dir, "/index.Rmd")))
            if (class(e) == "try-error") {
                e <- try(rmarkdown::render(paste0(out.dir, "/index.Rmd")))
            }
        }
    }
    if (html.open) {
        browseURL(paste0(out.dir, "/index.html"))
    }
}
