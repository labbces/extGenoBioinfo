library(ggplot2)

args<-commandArgs(trailingOnly = TRUE)
#print(args)

data<-read.delim(args[1])
data$repeatNumber<-data$forward_repeat_number + data$reverse_repeat_number
#head(data)

ggplot(data, aes(x=window,y=repeatNumber)) + 
  geom_point(size=0.5)+
  facet_grid(id ~ .) +
  theme(
    strip.text.y = element_text(size=3, face="bold"),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text = element_text(size = 5),
  )+
  ylab("Number of repeats") +
  xlab("Contig/Chromosome position")
