library(knitr)
library(kableExtra)

df <- read.csv('./tex_shit/df.csv')[ , 1:9]
new_names <- gsub("h2\\.\\.\\.", "h2=", names(df))
new_names[1] <- ""
names(df) <- new_names


table <- kable(df, "latex", booktabs = TRUE, digits = 3) %>%
  kable_styling(latex_options = "hold_position")


correct_document <- paste(
"\\documentclass{article}",
"\\usepackage[utf8]{inputenc}",
"\\usepackage[russian]{babel}",
"\\usepackage{booktabs}",
"\\usepackage{amsmath}",
"\\usepackage{graphicx}",
"\\usepackage{float}",
"\\begin{document}",
"",
table,
"",
"\\end{document}",
sep = "\n"
)


writeLines(correct_document, "./tex_shit/table.tex")
