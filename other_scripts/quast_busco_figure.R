library(ggplot2);
library(ggtext);
library(gridExtra);
library(scales);
library(stringr);

# quast and busco analysis; CA
figure_data = data.frame();

# extract data from quast
setwd("quast_stats_dir");

for (report in dir()){
  file_data = read.fwf(file = report, widths = c(28, 50), skip = 2, comment.char = "");
  figure_data = rbind(figure_data, c(str_trim(file_data[1, 2]), file_data[2, 2], file_data[8, 2], file_data[14, 2], file_data[16, 2], file_data[18, 2]));
}
colnames(figure_data) = c("name", "contigs_all", "length_all", "contigs_500", "length_500", "N50");
## note that N50 calculated with contigs > 500 bp (QUAST --min-contig parameter)


# extract data from busco
readShortSummary = function(df, directory, dataset){
  setwd(directory);
  df$data = 0;
  for (summary in dir()){
    file_data = read.csv(file = summary, skip = 9, sep = "\t", header = FALSE)[1:6, 2:3];
    completeness = round((as.numeric(file_data[1, 1]) / as.numeric(file_data[6, 1])) * 100, digits = 1); 
    species_name = strsplit(summary, split = ".", fixed = TRUE)[[1]][4];
    species_name <- gsub(pattern = '_busco_output', replacement = '', x = species_name); # had to do this for exceptions: aplysia, helobdella, or urchin_genome
    df$data[df$name == species_name] = completeness;
    ## note to make sure that the names match up
  }
  names(df)[names(df) == "data"] = dataset;
  return(df);
}

## note do not have the eukaryota short summaries for aplysia, helobdella, or urchin_genome
figure_data = readShortSummary(figure_data, directory = "../eukaryota_busco_short_summaries", dataset = "eukaryotaBUSCO");
figure_data = readShortSummary(figure_data, directory = "../metazoa_busco_short_summaries", dataset = "metazoaBUSCO");

# color coding for the plot

figure_data$group = "tunicate";
figure_data$group[figure_data$name == "masurca006_assembly"] = "our_hybrid";
figure_data$group[figure_data$name == "flyelongread_assembly"] = "our_long";
figure_data$group[figure_data$name == "bviol_SR_01072021_scaffolds"] = "our_short";
figure_data$group[figure_data$name == "B_leachii" | figure_data$name == "B_schlosseri"] = "botryllid";
figure_data$group[figure_data$name == "urchin_genome" | 
                    figure_data$name == "aplysia" | 
                    figure_data$name == "helobdella"] = "other";

# contigs plot
speciesToPlot = figure_data[c(1, 3, 4, 6, 7, 9:12, 14, 16, 18, 20),];
speciesToPlot$log10 = log10(speciesToPlot$contigs_all);
speciesToPlot$group = factor(speciesToPlot$group, levels = c("our_hybrid", "our_long", "our_short", "botryllid", "tunicate", "other"));
cols.to.numeric = c('contigs_all', 'length_all', 'contigs_500', 'length_500', 'N50', 'eukaryotaBUSCO', 'metazoaBUSCO');
speciesToPlot[cols.to.numeric] <- sapply(speciesToPlot[cols.to.numeric], as.numeric);

## metazoa

## contigs_500, using the --min-contig parameter
ggplot(speciesToPlot) +
  
  geom_point(aes(x = metazoaBUSCO, y = contigs_500, color = group, shape = group), size = 5) +
  
  scale_color_manual(name = expression(underline("Genomes")), 
                     labels = c(expression(paste(italic("B. violaceus"), " hybrid assembly")), expression(paste(italic("B. violaceus"), " long-read assembly")), expression(paste(italic("B. violaceus"), " short-read assembly")), "Other botryllids", "Selected other tunicates", "Selected other invertebrates"), 
                     values = c("blue", "#0096FF", "#88ACE0", "#FFC000", "#380474", "#228B22")) +
  
  scale_shape_manual(name = expression(underline("Genomes")), 
                    labels = c(expression(paste(italic("B. violaceus"), " hybrid assembly")), expression(paste(italic("B. violaceus"), " long-read assembly")), expression(paste(italic("B. violaceus"), " short-read assembly")), "Other botryllids", "Selected other tunicates", "Selected other invertebrates"),
                    values = c(17, 17, 17, 19, 19, 19)) +
  
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
        axis.title.y = element_text(size = 20));


## contigs_all, using all contigs regardless of their size (affects only the short and long read assembly)
contigs_all_plot <- ggplot(speciesToPlot) +
  
  geom_point(aes(x = metazoaBUSCO, y = contigs_all, color = group, shape = group), size = 5) +
  
  scale_color_manual(name = expression(underline("Genomes")), 
                     labels = c(expression(paste(italic("B. violaceus"), " hybrid assembly")), expression(paste(italic("B. violaceus"), " long-read assembly")), expression(paste(italic("B. violaceus"), " short-read assembly")), "Other botryllids", "Selected other tunicates", "Selected other invertebrates"), 
                     values = c("blue", "#0096FF", "#88ACE0", "#FFC000", "#380474", "#228B22")) +
  
  scale_shape_manual(name = expression(underline("Genomes")), 
                    labels = c(expression(paste(italic("B. violaceus"), " hybrid assembly")), expression(paste(italic("B. violaceus"), " long-read assembly")), expression(paste(italic("B. violaceus"), " short-read assembly")), "Other botryllids", "Selected other tunicates", "Selected other invertebrates"),
                    values = c(17, 17, 17, 19, 19, 19)) +
  
  scale_y_log10(limits = c(10, 1e+06), breaks = c(10^1, 10^2, 10^3, 10^4, 10^5, 10^6), labels = trans_format(trans = "log10", format = math_format(10^.x))) +
  
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
        axis.title.y = element_text(size = 20));


ggsave(filename = 'bviol_BUSCO_contigs_plot.png',
    plot = contigs_all_plot,
    width = 12,
    height = 8
    );

sessionInfo();