library(knitr)
library(kableExtra)

for( fname in list.files(path = "./tex_shit/csvs", full.names = F)){
  df <- read.csv(paste("./tex_shit/csvs/", fname, sep=""))
  new_names <- gsub("h2\\.\\.\\.", "h2=", names(df))
  new_names[1] <- ""
  names(df) <- new_names
  # print(fname)
  # print(df)
  table <- kable(df, "latex", booktabs = TRUE, digits = 3) %>%
    kable_styling(latex_options = "hold_position")
  table <- gsub("\\\\begin\\{tabular\\}", "\\\\hspace*{-2.3cm}\\\\begin\\{tabular\\}\n", table)
  writeLines(table, paste("./tex_shit/tables/", substr(fname ,1, nchar(fname)-4), ".tex" ,sep=""))
}