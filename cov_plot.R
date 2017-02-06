path = "/Volumes/huber.michael/HCV/experiments/170202/"

library(tidyverse)
library(stringr)
library(cowplot)

files = list.files(path, pattern = "cov.list")
data = data.frame()

for (i in files) {
        data_i = read.delim(paste0(path, i), header=FALSE)
        data_i = data_i %>%
                mutate(pos = 1:nrow(data_i)) %>%
                mutate(file = str_sub(i, 1, -12)) %>%
                mutate(iteration = str_sub(i, -10, -10)) %>%
                rename(cov = V1)
        data = rbind(data, data_i)
}

plot = data %>%
        ggplot(aes(x=pos, y=cov, color=iteration)) +
                geom_line(size=.2) +
                xlab('genome position') +
                ylab('coverage (reads)') +
                scale_y_log10(breaks = c(1,10,100,1000,10000,100000)) +
                facet_wrap( ~ file) +
                panel_border() +
                background_grid(major = "xy") +
                theme(axis.title = element_text(size = 7.5)) +
                theme(axis.text  = element_text(size = 7.5)) +
                theme(strip.text = element_text(size = 7.5)) +
                theme(legend.title = element_text(size = 7.5)) +
                theme(legend.text = element_text(size = 7.5)) +
                theme(legend.position="bottom")
print(plot)

ggsave(filename=paste0(path, "coverage.pdf"), plot, width = 30/2.54, height = 21/2.54)