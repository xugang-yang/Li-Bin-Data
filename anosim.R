# Bacterial stack diagram
setwd('~/Desktop/1')
otu<-read.csv('0.datas.txt',sep = '\t',header = T)
otu$ID <- factor(otu$ID, levels = otu$ID, ordered = TRUE)
otu1<-melt(otu)
head(otu1)
p<-ggplot(otu1, aes(variable,  value, fill =ID)) +geom_col(position = 'stack', width = .9) +
  scale_fill_manual(values =  rev(c( "#CACCC8", "lightskyblue","#92A1D3","cyan1", "khaki2", "firebrick","darkseagreen", "darkorchid","#ddc3c3","darkgreen","#11A187B2", "#F39B7FB2","#91D1C2B2","#8491B4B2", "#BCEE99",  "#7E6148B2", "yellow","pink","purple", 
                                     "royalblue4", "darksalmon","royalblue1","darkorange","cyan3", "#00861B","darkorange4", "#5F8E93","#6A4C72","#898456","#008B8B","#DC0000B2")))
p+labs(x = '', y = 'Relative abundance', fill = "Top 20 genara") +
  theme(panel.grid = element_blank(),
        panel.background = element_rect(color = '#636363',  fill = 'transparent'),
        strip.text = element_text(size = 20),
        strip.text.x = element_text(size = 20, colour = "black", face = "bold"),
        strip.background = element_rect(color="grey95",  fill="grey95",size=1.5,linetype="solid"),
        legend.title = element_text(size = 15,face = "bold", vjust = 1, hjust = 0),
        legend.text = element_text(size = 11), legend.key.size=unit(.4,'cm'),
        axis.title.x = element_text(size = 16, vjust = 0.5,  hjust = 0.5),
        axis.title.y = element_text(size = 15, vjust = 0.5, hjust = 0.5),
        axis.text.x = element_text(angle = 90, size = 10,vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(size = 15,vjust = 0.5,   hjust = 0.5),
        legend.position="bottom") + guides(fill=guide_legend(nrow=31,reverse=TRUE))  

Alpha diversity  
library("vegan")
setwd('~/Desktop/1')
otu1<-read.table('0.datas.txt',sep = '\t',header = T,row.names = 1)
phe <- read.delim('0.tra.txt', sep = '\t',row.names = 1,header = T)
phe$Genus_Shannon=diversity(t(otu1),"shannon")  
phe$Genus_Simpson=diversity(t(otu1),"simpson") 
phe$Genus_InSimpson=diversity(t(otu1),"inv") 
-----
PCoA 
#Pcoa
library('ape')
library('ggplot2')
library("ggpubr")
dis.bray<-vegdist(t(otu1), "bray",na.rm =T)#PCoA
dis.bray<-vegdist(t(otu1), "euclidean",na.rm =T)#euclidean
PCOA <- pcoa(dis.bray, correction="none", rn=NULL) 
result <-PCOA$values[,"Relative_eig"]
pro1 = as.numeric(sprintf("%.3f",result[1]))*100
pro2 = as.numeric(sprintf("%.3f",result[2]))*100
x = PCOA$vectors
sample_names = rownames(x)
pc = as.data.frame(PCOA$vectors)
pc$names = sample_names
legend_title = ""
group = phe$Group
pc$group = group
pchs = phe$Sample
pc$pch = pchs
xlab=paste("PCOA1(",pro1,"%)",sep="") 
ylab=paste("PCOA2(",pro2,"%)",sep="")
p<-ggplot(pc,aes(Axis.1,Axis.2)) + geom_point(size=3,aes(color=group)) + 
  labs(x=xlab,y=ylab,title="PCoA of bray distance in genus level",color=legend_title,shape=legend_title) + 
  geom_hline(yintercept=0,linetype=4,color="grey") + 
  geom_vline(xintercept=0,linetype=4,color="grey") + theme_bw()
P
----
boxplot  
library("ggplot2")
library("reshape2") 
library('pheatmap')
library("ggpubr")
library(rstatix)
otu<-read.csv('0.datas.txt',sep = '\t',header = T,row.names = 1)
otu<-na.omit(otu)
otu1<-melt(otu)
ggboxplot(otu1, x = "Group", y = "value",color = "Group", add = 'jitter', alpha=0.1 ) + facet_wrap(~variable,scales="free",nrow=2)+stat_compare_means(comparisons = my_comparisons1, label = "p.format",method = 't.test')
4.heatmap:
library(pheatmap)
library(psych)
library(reshape2)
setwd('~/Desktop/1')
met <- read.table(file = "0.datas.txt", sep = "\t", header = T,row.names=1)
phy <-read.table(file = "0.datas.txt", sep = "\t", header = T,row.names=1)
cor <-corr.test(phy,met,method = "spearman",adjust="none")
cmt <-na.omit(cor$r)
pmt <- na.omit(cor$p)
write.table(cmt,'MHO_Species_c.txt',sep='\t')
write.table(pmt,'MHO_Species_p.txt',sep='\t')
if (!is.null(pmt)){
  ssmt <- pmt< 0.01
  pmt[ssmt] <-'+'
  smt <- pmt >0.01& pmt <0.05
  pmt[smt] <- '*'
  pmt[!ssmt&!smt]<- ''
} else {
  pmt <- F
}
pheatmap(cmt,cluster_row = T, cluster_col = T,color = colorRampPalette(c("purple","white", "firebrick3"))(50),display_numbers = pmt)
# note: cofig order must the same pr
args <- commandArgs(T)ames(pr)
if (length(args) < 4) {rm", group[1:len], sep = "_"))
? stop("Rscript *.R [profile] [config] [alpha] [out]\n ? ? ? profile: profile.\n ? ? ? config: config.\n ? ? ? alpha: dunn alpha.\n ? ? ? out: out.")
}n <- c(cn, paste0("p", "_", group[row(dtp)[lower.tri(dtp)]], "_", group[col(dtp)[lower.tri(dtp)]]))
:
pr <- read.table(args[1],sep='\t',header = T,row.names = 1)
config <- read.table(args[2],sep='\t',header = T,row.names = 1)
alpha <- as.numeric(args[3])? # 0.05
out.dir <- args[4]
pr <- pr[, pmatch(rownames(config), colnames(pr))]
pr <- as.matrix(pr)
num <- nrow(pr)
fr <- as.factor(config[, 1])
group <- levels(fr)
len <- length(group)
print(fr)
out <- matrix(NA, num, (len + 1) * (len + 4)/2)

#ktp <- apply(pr, 1, function(x) {
#? kruskal.test(x ~ fr)$p.value
#})
library(PMCMRplus)
library(multcompView)
for (i in 1:num) {
? rk <- rank(pr[i, ])
? #
? ktp<-kruskal.test(pr[i,] ~ fr)$p.value?
? res <- c(ktp, tapply(rk, fr, mean), tapply(pr[i, ] > 0, fr, mean))
??
? # T
? if (!is.nan(ktp) & ktp < 2) {
? ? #
? ? dtp <- kwAllPairsDunnTest(pr[i, ], fr, p.adjust.method = "none")$p.value
?# ? dtp <- posthoc.kruskal.dunn.test(pr[i, ], fr)$p.value
? ? dtp[is.nan(dtp)] <- 1
? ? dtp <- cbind(dtp, rep(NA, len - 1))
? ? dtp <- rbind(rep(NA, len), dtp)
? ? dtp[upper.tri(dtp)] <- t(dtp)[upper.tri(dtp)]
? ? rownames(dtp)[1] <- colnames(dtp)[1]
? ? colnames(dtp)[len] <- rownames(dtp)[len]
?? ?
? ? #?
? ? res <- c(res, dtp[lower.tri(dtp)])
?? ?
? ? #?
? ? or <- order(res[1:len + 1])
? ? op <- rep(1, len - 1)
? ? for (j in 1:(len - 1)) {
? ? ? op[j] <- dtp[or[j], or[j + 1]]
? ? }
? ? symbol <- rep("<=", len - 1)
? ? symbol[!is.na(op) & op <= alpha] <- "<"
? ? symbol[is.na(op) | op == 1] <- "<=>"
? ? conclude <- rep(0, 2 * len - 1)
? ? conclude[(1:len) * 2 - 1] <- group[or]
? ? conclude[(1:(len - 1)) * 2] <- symbol
? ? res <- c(res, paste(conclude, collapse = " "))
?? ?
? ? #?
? ? res <- c(res, multcompLetters(dtp[or, or], threshold = alpha)$Letters[order(or)])
? ? # make sure 'a' < 'b'
? ? out[i, ] <- res
? } else {
? ? out[i, 1:(2 * len + 1)] <- res
? }
}
rownames(out) <- row.names(pr)
cn <- c("kw.p", paste("rm", group[1:len], sep = "_"))
cn <- c(cn, paste("or", group[1:len], sep = "_"))
cn <- c(cn, paste0("p", "_", group[row(dtp)[lower.tri(dtp)]], "_", group[col(dtp)[lower.tri(dtp)]]))
cn <- c(cn, "nearby", paste("sym", group[1:len], sep = "_"))
colnames(out) <- cn
write.table(out, out.dir, col.names = NA, sep = "\t", quote = F)?

