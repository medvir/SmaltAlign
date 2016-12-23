library(tidyverse)

path = "/Volumes/huber.michael/HCV/experiments/161212/"
files = list.files(path, pattern = "3_cov.list")

data = data.frame()
for (i in files) {
        data_i = read.delim(paste0(path, i), header=FALSE)
        data_i$V2 = seq.int(nrow(data_i))
        data_i$V3 = sub("_3_cov.list", "", i)
        data = rbind(data, data_i)
}

plot = data %>%
        ggplot(aes(x=V2, y=V1, color=V3)) +
                geom_line(size=.2) +
                xlab('genome position') +
                ylab('coverage (reads)') +
                scale_y_log10() + 
                facet_wrap(~ V3) +
                theme(legend.position="none")
print(plot)
ggsave(filename=paste0(path, "coverage.pdf"), plot, width = 30/2.54, height = 21/2.54)