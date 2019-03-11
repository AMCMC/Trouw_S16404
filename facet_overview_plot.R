#### init ####
library(dplyr)
library(ggplot2)
library(ggvegan)
library(phyloseq)

ps <- readRDS("ps.S16404.rds")

#### adonis bray ####

d = phyloseq::distance(ps, "bray")
d = phyloseq::distance(ps, "jaccard", binary=T)
df = as(sample_data(ps), "data.frame")
pvals=rep("NA",length(colnames(ps@sam_data)))
names(pvals) <- colnames(ps@sam_data)
Fstat <- pvals
R2s <- pvals
samples_used <- pvals
for (i in 1:length(colnames(ps@sam_data))){
  tryCatch({
    samples_used[colnames(ps@sam_data)[i]] <- sum(!is.na(df[,i]))
    test.out <- adonis(as.matrix(d)[!is.na(df[,i]),!is.na(df[,i])] ~ get(colnames(ps@sam_data)[i])[!is.na(df[,i])], df)
    pvals[colnames(ps@sam_data)[i]] <- test.out$aov.tab$`Pr(>F)`[1]
    Fstat[colnames(ps@sam_data)[i]] <- test.out$aov.tab$F.Model[1]
    R2s[colnames(ps@sam_data)[i]] <- test.out$aov.tab$R2[1]
  }, error=function(e){})
}

df.out <- cbind(samples_used[pvals!="NA"],pvals[pvals!="NA"],Fstat[pvals!="NA"],R2s[pvals!="NA"])[order(as.numeric(Fstat[pvals!="NA"]), decreasing = T),]
df.out <- setNames(data.frame(df.out), c("samples_used","p.value","F.model","R2"))
df.out$F.model <- round(as.numeric(as.character(df.out$F.model)), digits = 4)
df.out$p.value <- formatC(as.numeric(as.character(df.out$p.value)), format = "e", digits = 3)
df.out$R2 <- round(as.numeric(as.character(df.out$R2)), digits = 4)

#### DESEQ ####

library("DESeq2")
packageVersion("DESeq2")

ps.deseq <- ps
diagdds = phyloseq_to_deseq2(ps.deseq, ~ AX*CELL)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
DESeq2::resultsNames(diagdds)
res1 = results(diagdds, cooksCutoff = FALSE, name = "AX_lc.AXOS_vs_Control")
res2 = results(diagdds, cooksCutoff = FALSE, name = "CELL_CELL_vs_Control")
res3 = results(diagdds, cooksCutoff = FALSE, name = "AXlc.AXOS.CELLCELL")
res4 = results(diagdds, cooksCutoff = FALSE, name = "Intercept")


resc <- cbind(res1, res2, res3)
alpha = 0.05
resc.red <- resc[(resc[,6]< alpha & !is.na(resc[,6])) | (resc[,12]< alpha & !is.na(resc[,12])) | (resc[,18]< alpha & !is.na(resc[,18])),]

sigtab = resc.red
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps)[rownames(sigtab), ], "matrix"))
head(sigtab)
dim(sigtab)

#write.csv(sigtab, file = "DESEQ_fatorial_design.csv")


#### taxaplot ####

ps.taxaplot <- ps

# relabel NA clades
for (rank in 2:dim(ps.taxaplot@tax_table)[2]){
    if (sum(is.na(ps.taxaplot@tax_table[, rank])) > 0) {
      ps.taxaplot@tax_table[is.na(ps.taxaplot@tax_table[, rank]), rank] <-
        paste(ps.taxaplot@tax_table[is.na(ps.taxaplot@tax_table[, rank]), rank-1], "_Unclassified", sep = "")
    }
  }
  # Remove repetitive Unclassified.
  ps.taxaplot@tax_table <-
    gsub("_Unclassified_.*Unclassified",
         "_Unclassified",
         ps.taxaplot@tax_table)

# color vector to use
colours <-
  c(
    "#F0A3FF","#0075DC","#993F00","#4C005C","#2BCE48","#FFCC99","#808080","#94FFB5",
    "#8F7C00","#9DCC00","#C20088","#003380","#FFA405","#FFA8BB","#426600","#FF0010",
    "#5EF1F2","#00998F","#740AFF","#990000","#FFFF00"
  )

# How many clades to plot
N <- 21

# make relative & prune taxa
ps.taxaplot <- transform_sample_counts(ps.taxaplot, function(x) x/sum(x))
ps.taxaplot <- prune_taxa(taxa_names(ps.taxaplot) %in% names(sort(taxa_sums(ps.taxaplot), decreasing = TRUE)[1:N]), ps.taxaplot)

# generate dataframe
df <- psmelt(ps.taxaplot)
df$Taxa <- df$Genus

# order Taxa by global abundance
x <- aggregate(df$Abundance, by=list(df$Taxa), FUN=sum)
df$Taxa <- factor(df$Taxa, levels=x$Group.1[order(x$x, decreasing = F)])

# generate plot
p <- ggplot(df, aes_string("piglet2", "Abundance", fill = "Taxa"))
p <- p + geom_bar(stat = "identity")
p <- p +  scale_fill_manual(values = c(rev(colours[0:(N)])))


facet1 = "AX"
facet2 = "CELL"

p <- p + facet_grid(reformulate(facet1,facet2), scales="free_x", space="free")
p <- p + theme_bw() + ylab("Proportions")
p <- p + scale_y_continuous(expand = c(0, 0)) +
  theme(strip.background = element_rect(fill = "gray85")) +
  theme(panel.spacing = unit(0.3, "lines"))
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
p <- p + guides(fill = guide_legend(ncol = 1, reverse = F))
p <- p + theme(panel.spacing = unit(1, "lines"))
p <- p + theme(text = element_text(size = 20),
               axis.text.x = element_text(angle = 90, hjust = 1))
p.taxplot <- p.taxplot + theme(axis.title.x=element_blank(),
               axis.text.x=element_blank(),
               axis.ticks.x=element_blank())

#### CCA ggplot ####

#select variables
sam_data.nf <- data.frame(ps@sam_data[,!unlist(lapply(data.frame(ps@sam_data),class))=="factor"])
select <- c("BiW","BW.d23.24.kg","reladgtot_inv","AA","pH","VA",
            "PA","BA","LA","dm","IBU","IVA","TSCFA","NH3Cae")
sam_data.nf.sc <- sam_data.nf[,colnames(sam_data.nf) %in% select]

#rename variable for plotting purpose
colnames(sam_data.nf.sc)[colnames(sam_data.nf.sc)=="BiW"] <- "BWb"
colnames(sam_data.nf.sc)[colnames(sam_data.nf.sc)=="BW.d23.24.kg"] <- "BWa"
colnames(sam_data.nf.sc)[colnames(sam_data.nf.sc)=="reladgtot_inv"] <- "WGr"

#perform cca
x <- ps@otu_table
vare.cca <- cca(t(x) ~ . , sam_data.nf.sc)
df <- fortify(vare.cca)

#plot cca
p <- ggvegan:::autoplot.cca(df) + 
  theme_bw() + 
  NULL
p <- p + geom_point(aes(x=df$CCA1[df$Score=="sites"], 
                        y=df$CCA2[df$Score=="sites"], 
                        color=droplevels(ps@sam_data[df[df$Score=="sites",]$Label,]$Treat)), size=4)
p <- ggedit::remove_geom(p = p, "point", idx = 1) 
p <- p + labs(colour="Treatment") + scale_colour_viridis_d()
p <- p + labs(x=paste0("CCA1: ",round(as.vector(vare.cca$CCA$eig/sum(vare.cca$CCA$eig)*100)[1]),"% Variance explained ")) + 
  labs(y=paste0("CCA2: ",round(as.vector(vare.cca$CCA$eig/sum(vare.cca$CCA$eig)*100)[2]),"% Variance explained "))

p.cca <- p

#### adiv plot + betadisp ####

mod <- betadisper(vegdist(unclass(t(ps@otu_table))), ps@sam_data$Treat)
ps@sam_data$Beta_dispersion <- mod$distances[sample_names(ps)]

df <- cbind(ps@sam_data,estimate_richness(ps, measures = c("Observed","Shannon")))
df <- df[df$Sample.Type=="Mid-Colon",]

df.l <- reshape2::melt(df)
df.l <- df.l[df.l$variable %in% c("Observed","Shannon","Beta_dispersion"),]
df.l$Treatment <- df.l$Treat
df.l$variable <- factor(df.l$variable, levels = c("Observed","Shannon","Beta_dispersion"))

p1 <- ggplot(df.l, aes(x=Treatment, y=value, color=Treatment))+ 
  geom_boxplot(alpha = 0.8, outlier.colour = NA, size=1) +
  ggbeeswarm::geom_beeswarm(size=3, cex = 3) + 
  facet_wrap(~variable, scales="free_y") +
  theme(legend.position="none") + 
  scale_fill_viridis_d() + 
  scale_color_viridis_d() + 
  theme_bw() + 
  labs(x="") + 
  NULL

p.adiv <- p1

#### AX*CELL boxplot (Deseq top 6) ####

# read DESEQ results

sigtab <- read.csv("DESEQ_fatorial_design.csv", sep=";", header = T, skip=1)

# taxglom and select top 6 genera
ps.genus <- tax_glom(ps, taxrank = "Genus")
ps.genus <- transform_sample_counts(ps.genus, function(x) x/100000)

ps.genus.sub <- prune_taxa(ps.genus@tax_table[,"Genus"] %in% sigtab$Genus[1:6], ps.genus)
ps.genus.sub.melt <- psmelt(ps.genus.sub)

# abundance boxplots
p.boxplots <- ggplot(ps.genus.sub.melt, aes(x=Treat, y=Abundance, color=Treat)) + 
  geom_boxplot() + 
  ggbeeswarm::geom_beeswarm(size=3, cex = 3) +
  facet_wrap(~Genus, scales="free_y") + 
  scale_colour_viridis_d() + 
  labs(y="Relative Abundance", x="") + 
  labs(colour="Treatment") + 
  theme_bw() + 
  scale_y_log10() + 
  theme(strip.text.x = element_text(size = 8)) + 
  NULL

p.boxplots

#### merge plots ####


p1 <- p.taxplot
p2 <- p.cca
p3 <- p.adiv
p4 <- p.boxplots

p1 <- p1 + 
  theme(text = element_text(size = 15)) + 
  labs(title="A) Genus composition") + 
  theme(legend.text=element_text(size=8)) + 
  theme(legend.key.size=unit(4,'mm')) +
  theme(panel.spacing = unit(0.1, "lines"))

p2 <- p2 + 
  theme(text = element_text(size = 15)) + 
  labs(title="B) Beta-Diversity (CCA)")

p3 <- p3 + 
  theme(text = element_text(size = 15)) + 
  labs(title="C) Alpha-Diversity", y="") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  labs(x="")

p4 <- p4 + theme(text = element_text(size = 15)) + labs(title="D) Responding Genera") + 
  theme(legend.text=element_text(size=8)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  theme(legend.key.size=unit(4,'mm'))

lay <- rbind(c(1,1,1,1,1,2,2,2),
             c(3,3,3,3,4,4,4,4))

merged <- gridExtra::grid.arrange(p1 + theme(plot.margin = unit(c(1,1,0.8,0.8), "cm"))+ 
                          theme(legend.text=element_text(size=8)) +
                          guides(fill = guide_legend(keyheight=0.5, keywidth=0.6, ncol = 1)), 
                        p2 + theme(plot.margin = unit(c(1,0,0,0), "cm")), 
                        p3 + theme(legend.position="none") + 
                          theme(plot.margin = unit(c(0.2,1,0,1), "cm")), 
                        p4 + theme(legend.position="none") + 
                          theme(plot.margin = unit(c(0.2,1,0,0), "cm")),
                        #  grobs = gl,
                        layout_matrix = lay)

pdf("Overview_microbiota_4P_13x8.6.pdf", width = 13, height = 8.6)
plot(merged)
dev.off()