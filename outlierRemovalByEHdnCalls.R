# This script is to generate the qq-plot in extended figure 3.
# It generally plot the distribution of number of EHdn calls per sample.
# Then try to see which cutoff is more appropriate to exclude outliers.

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(ggplot2)
library(cowplot)
library(qqplotr)
rm(list=ls())

al.stats <- read.delim("MSSNG_alignment_stats.csv", stringsAsFactors = F)
tmp.count <- read.delim("call_counts.txt", stringsAsFactors = F)
sample.info <- read.delim("Data_Merged_MSSNG.SSC.1KG.tsv", stringsAsFactors = F)
sample.info <- sample.info[sample.info$Library.type == "PCR-free" & sample.info$Dataset != "1000G", ]
tmp.count <- tmp.count[tmp.count$Sample %in% sample.info$Sample.ID, ]
tmp.count <- merge(tmp.count, al.stats[, c("samle", "AVG_READ_DEPTH")], by.x = "Sample", by.y = "samle", all.x = T)
readdepth <- read.delim("EHDN_read_depths.tsv", stringsAsFactors = F)

tmp.count <- merge(tmp.count, readdepth, by = "Sample", all.x = T)

mean <- mean(tmp.count$EHDN.call.count)
sd <- sd(tmp.count$EHDN.call.count)
tmp.count$outlier <- F
tmp.count$outlier[tmp.count$EHDN.call.count < mean-3*sd |
                    tmp.count$EHDN.call.count > mean+3*sd] <- T

p1 <- ggplot(tmp.count, aes(x = EHDN.call.count, y = ..density..)) + geom_histogram(fill = "grey", bins = 100) + 
  theme_bw() + ylab("Density") + 
  geom_vline(xintercept = c(mean, mean+sd, mean+2*sd, mean+3*sd, mean-sd, mean-2*sd, mean-3*sd), 
             lty = 2, lwd = 0.1, color = "black") + xlab("EHdn detected tandem repeats") +
  annotate("text", x = c(mean - 200, mean+sd - 200, mean+2*sd - 200, mean+3*sd - 200,
                         mean-sd - 200, mean-2*sd - 200, mean-3*sd - 200), 
           y =  0.0004, label = 
             c("mean", "+1SD", "+2SD", "+3SD", "-1SD", "-2SD", "-3SD"), angle = 90, cex = 3)
p2 <- ggplot(tmp.count, aes(sample = EHDN.call.count)) + stat_qq_band() + stat_qq_point(color = "steelblue") + stat_qq_line() +
  scale_color_manual(values = c("steelblue", "red")) + ylab("Sample Quantiles") + xlab("Theoretical Quantiles") +
  geom_vline(xintercept = c(mean, mean+sd, mean+2*sd, mean+3*sd, mean-sd, mean-2*sd, mean-3*sd), lty = 2, lwd = 0.1, color = "black") + 
  annotate("text", x = c(mean - 200, mean+sd - 200, mean+2*sd - 200, mean+3*sd - 200,
                         mean-sd - 200, mean-2*sd - 200, mean-3*sd - 200), y = 10000, label = 
             c("mean", "+1SD", "+2SD", "+3SD", "-1SD", "-2SD", "-3SD"), angle = 90, cex = 3) 
plot_grid(p1, p2, nrow = 1, rel_widths = c(1.5, 1))
ggsave("ext.fig2a.pdf", width = 10, height = 4)

### 2SD
tmp.clean <- tmp.count[tmp.count$EHDN.call.count > mean-2*sd &
                         tmp.count$EHDN.call.count < mean + 2*sd, ]


p1 <- ggplot(tmp.clean, aes(x = EHDN.call.count, y = ..density..)) + geom_histogram(fill = "grey", bins = 100) +  xlab("EHdn detected tandem repeats") +
  theme_bw() + ylab("Density")
p2 <- ggplot(tmp.clean, aes(sample = EHDN.call.count)) + stat_qq_band() +
  stat_qq_point(color = "steelblue") + stat_qq_line() +
  ylab("Sample Quantiles") + xlab("Theoretical Quantiles")
plot_grid(p1, p2, nrow = 1, rel_widths = c(1.5, 1))
ggsave("ext.fig2b.pdf", width = 10, height = 4)

### 3SD
tmp.clean <- tmp.count[tmp.count$EHDN.call.count > mean-3*sd &
                         tmp.count$EHDN.call.count < mean + 3*sd, ]


p1 <- ggplot(tmp.clean, aes(x = EHDN.call.count, y = ..density..)) + geom_histogram(fill = "grey", bins = 100) +  xlab("EHdn detected tandem repeats") +
  theme_bw() + ylab("Density")
p2 <- ggplot(tmp.clean, aes(sample = EHDN.call.count)) + stat_qq_band() + 
  stat_qq_point(color = "steelblue") + stat_qq_line() +
  ylab("Sample Quantiles") + xlab("Theoretical Quantiles")
plot_grid(p1, p2, nrow = 1, rel_widths = c(1.5, 1))
ggsave("ext.fig2c.pdf", width = 10, height = 4)

tmp.clean <- merge(tmp.clean, sample.info, by.x = "Sample", by.y = "Sample.ID", all.x =T)
t.test(tmp.clean$EHDN.call.count[tmp.clean$Dataset == "SSC"], 
            tmp.clean$EHDN.call.count[tmp.clean$Dataset == "MSSNG"])


tmp.aff <- tmp.clean[tmp.clean$Relation %in% c("affectedsibling", "proband", "unaffectedsibling"), ]
ggplot(tmp.aff, aes(sample = EHDN.call.count, color = Dataset)) + stat_qq_band() + stat_qq_point() + stat_qq_line()
d1 <- ggplot(tmp.aff, aes(x = EHDN.call.count, y = ..density..,fill = Dataset)) + geom_density(alpha = .5) +
  theme_bw()
d2 <- ggplot(tmp.aff[tmp.aff$Dataset == "SSC", ], aes(x = EHDN.call.count, y = ..density..,fill = Dataset)) + 
  geom_histogram(bins = 100, fill = "turquoise") +xlim(min(tmp.aff$EHDN.call.count), max(tmp.aff$EHDN.call.count)) +
  theme_bw()
d3 <- ggplot(tmp.aff[tmp.aff$Dataset == "MSSNG", ], aes(x = EHDN.call.count, y = ..density..,fill = Dataset)) + 
  geom_histogram(bins = 100, fill = "brown2") +xlim(min(tmp.aff$EHDN.call.count), max(tmp.aff$EHDN.call.count)) +
  theme_bw()
  
plot_grid(d1, d2, d3, nrow = 3)
ggsave("3SD.affected.png", width = 5, height = 7)
