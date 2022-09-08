library(ggplot2)
library(ggtext)
library(gridExtra)
library(scales)

# quast and busco analysis; CA
figure_data = data.frame()


# extract data from quast
setwd("~/Desktop/quast_stats_dir")

for (report in dir()){
  file_data = read.fwf(file = report, widths = c(28, 50), skip = 2, comment.char = "")
  figure_data = rbind(figure_data, c(file_data[1, 2], as.numeric(file_data[14, 2]), as.numeric(file_data[18, 2])))
}
colnames(figure_data) = c("name", "contigs", "N50")
figure_data$contigs = as.numeric(figure_data$contigs)
figure_data$N50 = as.numeric(figure_data$N50)



# extract data from busco
readShortSummary = function(df, directory, dataset){
  setwd(directory)
  df$data = 0
  for (summary in dir()){
    file_data = read.csv(file = summary, skip = 9, sep = "\t", header = FALSE)[, 2:3]
    completeness = signif((as.numeric(file_data[1, 1]) / as.numeric(file_data[6, 1])) * 100, digits = 4) ## check this
    species_name = strsplit(summary, split = ".", fixed = TRUE)[[1]][4]
    df$data[df$name == species_name] = completeness
    ## note to make sure that the names match up; for the short read assembly the names didn't match when the output came from quast and this so had to change name
  }
  names(df)[names(df) == "data"] = dataset
  return(df)
}

figure_data = readShortSummary(figure_data, directory = "~/Desktop/eukaryota_busco_short_summaries", dataset = "eukaryotaBUSCO")
figure_data = readShortSummary(figure_data, directory = "~/Desktop/metazoa_busco_short_summaries", dataset = "metazoaBUSCO")
figure_data[9, 5] = 93.3
figure_data[10, 5] = 91.2
figure_data[11, 5] = 97.7

# color coding for the plot

figure_data$group = "tunicate"
figure_data$group[figure_data$name == "masurca006_assembly"] = "our_hybrid"
figure_data$group[figure_data$name == "flyelongread_assembly"] = "our_long"
figure_data$group[figure_data$name == "bviol_SR_01072021_scaffolds"] = "our_short"
figure_data$group[figure_data$name == "B_leachii" | figure_data$name == "B_schlosseri"] = "botryllid"
figure_data$group[figure_data$name == "urchin_genome" | 
                    figure_data$name == "aplysia" | 
                    figure_data$name == "helobdella"] = "other"
figure_data[9, 6] = "other" # not sure why it would not change like the others with the above line of code

# contigs plot

speciesToPlot = figure_data[c(1, 3, 4, 6, 7, 9:12, 14, 16, 18, 20),]
speciesToPlot$logs = log(speciesToPlot$contigs)
speciesToPlot$group = factor(speciesToPlot$group, levels = c("our_hybrid", "our_long", "our_short", "botryllid", "tunicate", "other"))


## metazoa

ggplot(speciesToPlot) +
  
  geom_point(aes(x = metazoaBUSCO, y = contigs, color = group, shape = group), size = 5) +
  
  scale_color_manual(name = "Genomes", 
                     labels = c(expression(paste(italic("B. violaceus"), " hybrid assembly")), expression(paste(italic("B. violaceus"), " long-read assembly")), expression(paste(italic("B. violaceus"), " short-read assembly")), "Other botryllids", "Selected other tunicates", "Selected other invertebrates"), 
                     values = c("blue", "#0096FF", "#88ACE0", "#FFC000", "#380474", "#228B22")) +
  
  scale_shape_manual(name = "Genomes", 
                     labels = c(expression(paste(italic("B. violaceus"), " hybrid assembly")), expression(paste(italic("B. violaceus"), " long-read assembly")), expression(paste(italic("B. violaceus"), " short-read assembly")), "Other botryllids", "Selected other tunicates", "Selected other invertebrates"), values = c(17, 17, 17, 19, 19, 19)) +
  
  scale_y_log10(limits = c(10, 1e+05), labels = trans_format(trans = "log10", format = math_format(10^.x))) +
  
  scale_x_continuous(breaks = seq(from = 50, to = 100, by = 10), labels = paste(seq(from = 50, to = 100, by = 10), "%", sep = ""), limits = c(50, 100)) +
  
  labs(x = "BUSCO Score", y = "Scaffolds and Contigs") +
  
  theme_classic() +
  
  theme(legend.box.background = element_rect(color = "black", size = 1), 
        legend.text.align = 0,
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 20),
        axis.text.x = element_text(size = 18, color = "black"), 
        axis.text.y = element_text(size = 18, color = "black"), 
        axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 20))

## make fonts and overall plot larger

## completeness plot
ggplot(figure_data) + geom_bar(aes(x = name, y = BUSCO), stat = "identity") + coord_flip()

setwd("~/Desktop/")
pdf("quast_busco_figure.pdf")
dev.off()