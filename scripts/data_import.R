# Data origin
# https://doi.org/10.7717/peerj.14142
# 

df = read.csv(here::here("data",
                         "dataset.csv")) |>
  janitor::clean_names() |>
  dplyr::select(code, group, x1rm_pre, x1rm_post) 
