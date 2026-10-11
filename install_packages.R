# Install every package used by the pipeline (run once)
pk <- c("readxl", "dplyr", "tidyr", "purrr", "stringr", "forcats", "tibble", "ggplot2", "patchwork", "scales",
        "ragg", "systemfonts", "colorspace", "ggridges", "ggdist", "ggalluvial", "ggrepel", "ggtext", "FactoMineR",
        "mclust", "mgcv", "relaimpo", "quantreg", "MASS", "splines", "qgraph", "igraph", "ggraph", "tidygraph",
        "lavaan", "partykit", "pROC", "survival", "flextable", "officer", "car", "AER", "mice", "janitor", "broom",
        "Hmisc")
new <- setdiff(pk, rownames(installed.packages()))
if (length(new)) install.packages(new, repos = "https://cloud.r-project.org")
