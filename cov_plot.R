args <- commandArgs(TRUE)

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(cowplot)

# set path to data (terminal slash is important!)
#path <- "/Volumes/data/Diagnostics/experiments/170815/"
path <- paste0(args[1],"/")

files = list.files(path, pattern = "depth")
data = data.frame()

for (i in files) {

        if (file.info(paste0(path, i))$size == 0) {
                depth_i = c(0) #if depth file is empty
                names(depth_i) = c(1)
            } else {
                depth_i = read_delim(paste0(path, i), "\t", col_names = FALSE, trim_ws = TRUE, col_types = "cii") %>% pull(X3)
                names(depth_i) = read_delim(paste0(path, i), "\t", col_names = FALSE, trim_ws = TRUE) %>% pull(X2)
            }
        data_i = data.frame(pos = 1:max(as.numeric(names(depth_i)))) %>%
                mutate(cov = ifelse(is.na(depth_i[as.character(pos)]), 0, depth_i[as.character(pos)])) %>%
                mutate(file = str_sub(i, 1, -9)) %>%
                mutate(iteration = str_sub(i, -7, -7))
        data = rbind(data, data_i)
}

plot = data %>%
         ggplot(aes(x = pos, y = cov, color = iteration)) +
                geom_line(size = .2) +
                xlab("genome position") +
                ylab("coverage (reads)") +
                scale_y_log10(breaks = c(1,10,100,1000,10000,100000)) +
                facet_wrap( ~ file) +
                panel_border() +
                background_grid(major = "xy") +
                theme(axis.title = element_text(size = 7.5)) +
                theme(axis.text  = element_text(size = 7.5)) +
                theme(strip.text = element_text(size = 7.5)) +
                theme(legend.title = element_text(size = 7.5)) +
                theme(legend.text = element_text(size = 7.5)) +
                theme(legend.position = "bottom")

ggsave(filename = paste0(path, "coverage.pdf"), plot, width = 30/2.54, height = 21/2.54)
