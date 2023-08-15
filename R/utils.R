writehtmlmod <- function (x, links, tiny.pic, tiny.pic.size = 100, title = "",
    file, digits = 3, col.align='center', cellpadding=10, sorttable, dragtable)
{
    stopifnot(is(x, 'data.frame'))
    if (missing(links))
        links <- vector("list", ncol(x))
    if (missing(tiny.pic))
        tiny.pic <- vector("list", ncol(x))
    stopifnot(is(links, 'list'))
    stopifnot(is(tiny.pic, 'list'))
    stopifnot(length(links) == ncol(x))
    stopifnot(length(tiny.pic) == ncol(x))
    stopifnot(!missing(file))
    column.class <- unlist(lapply(x, class))
    for (j in seq_len(ncol(x))) {
        if (column.class[j] == "factor")
            x[, j] <- as.character(x[, j])
        if (column.class[j] == "numeric")
            x[, j] <- round(x[, j], digits = digits)
    }
    cat(paste0("<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01//EN\"",
            "\"http://www.w3.org/TR/html4/strict.dtd\">\n"),
        sep = "", file = file)
    cat("<html>\n", file = file, append = TRUE)
    cat("<body>\n", file = file, append = TRUE)
    cat(paste("<CAPTION ALIGN=\"top\"><center><B>",
            title, "</B></center></CAPTION><BR>\n"), sep = "", file = file,
        append = TRUE)
    cat(paste("<TABLE border=1 cellpadding=", cellpadding,">\n", sep=''),
        file = file, append = TRUE)
    cat("<TR>\n", file = file, append = TRUE)
    for (j in seq_len(ncol(x))) {
        cat("<TH>", file = file, append = TRUE)
        cat(colnames(x)[j], file = file, append = TRUE)
        cat("</TH>\n", file = file, append = TRUE)
    }
    cat("</TR>\n", file = file, append = TRUE)
    for (i in seq_len(nrow(x))) {
        cat("<TR>\n", file = file, append = TRUE)
        for (j in seq_len(ncol(x))) {
            cat(paste("<TD align=", col.align, ">", sep=''), file = file,
                append = TRUE)
            if (is.null(links[[j]]) & is.null(tiny.pic[[j]])) {
                cat(x[i, j], file = file, append = TRUE)
            }
            else if (is.null(links[[j]]) & !is.null(tiny.pic[[j]])) {
                cat(paste("<A HREF=\"", links[[j]][[i]], "\"><img src=\"",
                    tiny.pic[[j]][[i]], "\" height=\"", tiny.pic.size,
                    "\" width=\"", tiny.pic.size, "\" /></A>",
                    sep = ""), file = file, append = TRUE)
            }
            else if (!is.null(links[[j]]) & is.null(tiny.pic[[j]])) {
                cat(paste("<A HREF=\"", links[[j]][[i]], "\">",
                    x[i, j], "</A>", sep = ""), file = file, append = TRUE)
            }
            else if (!is.null(links[[j]]) & !is.null(tiny.pic[[j]])) {
                cat(paste("<A HREF=\"", links[[j]][[i]], "\"><img src=\"",
                    tiny.pic[[j]][[i]], "\" height=\"", tiny.pic.size,
                    "\" width=\"", tiny.pic.size, "\" /></A>",
                    sep = ""), file = file, append = TRUE)
            }
            cat("</TD>\n", file = file, append = TRUE)
        }
        cat("</TR>\n", file = file, append = TRUE)
    }
    cat("</TABLE>\n", file = file, append = TRUE)
    cat("</body>\n", file = file, append = TRUE)
    cat("</html>\n", file = file, append = TRUE)
    sortDragHtmlTableint(filename = file, sorttable, dragtable)
}

sortDragHtmlTableint <- function (filename, sorttable, dragtable)
{
    lastSlashPos <-
        gregexpr(.Platform$file.sep, filename)[[1]][length(gregexpr(
            .Platform$file.sep,filename)[[1]])]
    outputdir <- ifelse(lastSlashPos == -1, getwd(), substr(filename,
        0, lastSlashPos))
    fileCode <- ""

    writeLines(sorttable, con = paste(outputdir, "sorttable.js",
        sep = .Platform$file.sep))

    writeLines(dragtable, con = paste(outputdir, "dragtable.js",
        sep = .Platform$file.sep))
    tmpTxt <- readLines(filename)
    tmpTxt[1] <- paste("<script src=\"sorttable.js\"></script>\n",
        tmpTxt[1])
    tmpTxt[1] <- paste("<script src=\"dragtable.js\"></script>\n",
        tmpTxt[1])
    tmpTxt <- sub("TABLE", "TABLE class=\"draggable sortable\"",
        tmpTxt)
    writeLines(tmpTxt, con = filename)
}
