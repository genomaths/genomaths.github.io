# ================= To generate HTML manual ================================#


dir <- "~/MethylIT_new/MethylIT/man/"
rds <- list.files(path = dir, pattern = ".Rd")
out <- "/media/sf_D/MethylIT_HTML_Manual/"
css <- "/media/sf_D/MethylIT_HTML_Manual/R.css"

# ------------------------ To write the HTML files ---------------------------#
for(k in 1:length(rds)) {
  tools::Rd2HTML(Rd = paste0(dir, rds[k]), out = paste0(out, sub(".Rd", ".html", rds[k])),
                 package = "MethylIT", stylesheet = css)
}

# -------------------- To write the links used in rmarkdown -------------------#

url <- "https://genomaths.github.io/MethylIT_HTML_Manual/"
x <- paste0("# ", "[",gsub(".Rd", "", rds),"]", "(", url, list.files(out),")")

# ---------------- Add link to the HTML file --------------------------------- #
files <- sapply(1:length(rds), function(k) paste0(out, sub(".Rd", ".html", rds[k])))

repls <- paste0("<a href=", "'https://genomaths.github.io/MethylIT_HTML_Manual/MethylIT_Manual.html'",
                ">", "MethylIT</a>")
for(k in 1:length(files)) {
  tx  <- readLines(files[k])
  tx2  <- gsub(pattern = "MethylIT", replacement = repls, x = tx)
  writeLines(tx2, con = files[k])
}






