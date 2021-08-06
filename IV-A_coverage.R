library(tidyverse)

selected_cols <- cols_only(position = col_integer(), coverage = col_integer(),
                           segment = col_integer(), molis = col_integer())

data_path <- '/Volumes/analyses/Diagnostics/experiments/tables-190626_IV-A'
files <- dir(path = data_path, pattern = "*.csv")

# create a data frame holding the file names
data <- tibble(filename = files) %>%
  # a new data column
  mutate(file_contents = map(filename, ~ read_csv(file.path(data_path, .),
                                                  col_types = selected_cols))) %>%
  separate(filename, c('molis', 'segment', 'ext'))

#coverage_avg <- data %>%
#  unnest() %>%
#  group_by(molis, segment) %>%
#  summarise(avg = round(mean(coverage), 1)) %>%
#  write_csv(paste0(data_path, '/cov_average.csv'))

for (molis_nr in unique(data$molis)){
  print(molis_nr)
 p <- data %>%
    filter(molis == molis_nr) %>%
    unnest() %>%
    ggplot(aes(x = position, y = coverage)) +
    geom_bar(stat = 'identity') +
    facet_wrap(~ segment, scales = "free") +
    labs(caption = molis_nr)
 ggsave(paste0('covplot-', molis_nr, '.pdf'), path = data_path)
}
