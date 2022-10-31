#' @importFrom magrittr %>%
NULL

#' Save DE results as JSON tables for viewing in browser
#'
#' @param de.raw List of DE results from e.g. cacoa, conos
#' @param sample.groups Sample groups as named list, each element containing a vector of samples. Can be retrieved from e.g. package cacoa (default=NULL)
#' @param saveprefix Prefix for created files (default=NULL)
#' @param dir.name Name for directory with results. If it doesn't exist, it will be created. To disable, set as NULL (default="JSON")
#' @param gene.metadata (default=NULL) # Needs explanation
#' @param verbose Show progress (default=TRUE)
#' @return JSON files, table of content, and viewer files for viewing DE results in browser
#' @examples
#' \dontrun{
#' saveDeAsJson(de.raw, sample.groups)
#' }
#' ## The results can be viewed in a webbrowser by opening toc.html
#'
#' @export
saveDeAsJson <- function(de.raw, sample.groups=NULL, saveprefix=NULL, dir.name="JSON", gene.metadata=NULL,verbose=TRUE) {

  checkPackageInstalled("jsonlite", cran=TRUE)
  if (!is.null(dir.name)) {
    if (!dir.exists(dir.name)) dir.create(dir.name)
  } else {
    dir.name <- "."
  }

  if (is.null(gene.metadata)) {
    gene.metadata <- lapply(de.raw, function(x) {
      if (!is.null(x)) {
        rownames(as.data.frame(x$res))
      } else {
        NULL
      }
    }) %>%
      unlist() %>%
      unique() %>%
      {data.frame(geneid = .)}
  } else {
    if(is.null(gene.metadata$gene.id)) stop("gene.metadata must contain $gene.id field")
  }

  lapply(sn(de.raw %>% names()), function(ncc) {
    if(verbose) print(ncc)
    res.celltype <- de.raw[[ncc]]
    res.table <- res.celltype$res %>% as.data.frame()
    res.table$gene <- rownames(res.table)
    res.table$significant <- res.table$padj < 0.05
    res.table$log2FoldChange[is.na(res.table$log2FoldChange)] <- 0
    res.table$rowid <- 1:nrow(res.table)

    all.genes <- rownames(res.table)
    cm <- res.celltype$cm

    ilev <- lapply(sample.groups, function(sg) {
      sg <- sg[sg %in% colnames(cm)]
      cm.tmp <- cm[,sg]
      cm.tmp <- as.matrix(cm.tmp)
      rownames(cm.tmp) <- rownames(cm)

      ## calculate cpm
      cpm <- sweep(cm.tmp, 2, apply(cm.tmp,2, sum), FUN='/')
      cpm <- log10(cpm * 1e6 + 1)
      snames1 <- colnames(cpm)

      ## Put genes in order
      cpm <- cpm[all.genes,]
      colnames(cpm) <- NULL;
      rownames(cpm) <- NULL;

      list(snames=snames1, val=as.matrix(cpm))
    })

    snames <- names(res.celltype$sample.groups)

    ## convert to json
    tojson <- list(
      res = res.table,
      genes = all.genes,
      ilev = ilev,
      snames = snames)
    y <- jsonlite::toJSON(tojson)
    file <- paste0(dir.name, "/", saveprefix, make.names(ncc), ".json")
    write(y, file)
    NULL
  })

  toc.file <- paste0(dir.name, "/toc.html")
  s <- c(list('<html><head><style>
    table {
    font-family: arial, sans-serif;
    border-collapse: collapse;
    width: 100%;
    }

    td, th {
    border: 1px solid #dddddd;
    text-align: left;
    padding: 8px;
    }

    tr:nth-child(even) {
    background-color: #dddddd;
    }

    </style></head><body><table>'),
         lapply(names(de.raw), function(n)
           paste0('<tr><td><a href="deview.2.html?d=', saveprefix, make.names(n),'.json">', n, '</a></td></tr>')),
         list('</table></body></html>')
  ) %>% paste(collapse='\n')

  write(s, file=toc.file)

  # Copy viewer files
  file.copy(from = system.file("extdata", c("deview.2.html", "deview.2.js"), package="sccore"),
            to = paste0(dir.name),
            overwrite = TRUE)

}
