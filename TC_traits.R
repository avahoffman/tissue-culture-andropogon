###########################################################################################
##
## R source code to accompany Hoffman and Smith (2018) Thinking inside the box: tissue culture 
## for plant propagation in a key ecological species, Andropogon gerardii, last updated 19 July 2018.
## Please contact Ava Hoffman (avamariehoffman@gmail.com) with questions.
##
## If you found this code useful, please use the citation below:
## 
## Hoffman, A.M., Smith, M.D., 2018. Thinking inside the Box: Tissue Culture for Plant 
## Propagation in a Key Ecological Species, Andropogon gerardii. Am. J. Plant Sci. 09, 1987–2003. 
## https://doi.org/10.4236/ajps.2018.910144
##
###########################################################################################

wd <- "/Users/avahoffman/Dropbox/Research/Tissue\ Culture/tissue-culture-andropogon"
setwd(wd)

###########################################################################################

library(ggplot2)
#library(gridExtra)
#  require(gridExtra)
library(cowplot)
#library(multcomp)
library(viridis)
library(MASS)

###########################################################################################


TCG=read.table("TC_traits.txt",header=T) 
	#make sure no spaces in column titles
#RM$trt=as.factor(RM$trt)
	#check your variables
head(TCG);str(TCG)
TCG <- TCG[(TCG$day == 6),]
TCG$geno <- as.factor(TCG$geno)

#for the weights, putting it in grams because mg was too big of units to look good on the figure
TCG$wt=TCG$wt*0.001

#code for determining standard error for graphs.
min.mean.sd.max <- function(x) {
  r <- c(min(x), mean(x) - sd(x), mean(x), mean(x) + sd(x), max(x))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

plot1=ggplot(aes(y = wt, x = factor(geno), color=geno), data = TCG) + 
  #xlab("Genotype") + 
  ylab("Mass (g)") +
  theme_classic() +
  theme(legend.position="none") +
  scale_color_viridis(discrete=T, end=0.9, name="Genotype") +
  stat_summary(fun.data = min.mean.sd.max, geom = "boxplot") +
  geom_jitter(size=1.0) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank())
plot2=ggplot(aes(y = ht, x = factor(geno), color=geno), data = TCG) + 
  #xlab("Genotype") + 
  ylab("Height (cm)") +
  theme_classic() +
  theme(legend.position="none") +
  scale_color_viridis(discrete=T, end=0.9, name="Genotype") +
  stat_summary(fun.data = min.mean.sd.max, geom = "boxplot") +
  geom_jitter(size=1.0) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank())
plot3=ggplot(aes(y = til, x = factor(geno), color=geno), data = TCG) + 
  xlab("Genotype") + 
  ylab("Tiller count") +
  theme_classic() +
  theme(legend.position="none") +
  scale_color_viridis(discrete=T, end=0.9, name="Genotype") +
  stat_summary(fun.data = min.mean.sd.max, geom = "boxplot") +
  geom_jitter(size=1.0)
   
#pdf(file="fig3.pdf",height=7,width=3)
grid=plot_grid( plot1,plot2,plot3,
                   align = 'v',
                   labels = c("(a)", "(b)", "(c)"),
                   label_size=12,
                   hjust = 0,
                   nrow = 3
)
plot_grid(grid)
ggsave(file="fig3.jpg",height=7,width=3)
#dev.off()

#########################
#
#
# analyses looking for differences among TC genotypes
# 
#
#
#########################

# Tiller differences
mod <- rlm(til~geno, data=TCG)
### Note, Tukey means the type of contrast matrix.  See ?contrMat
glht.mod <- glht(mod, mcp(geno = "Tukey"))
###summaryize information
###applying the false discovery rate adjustment
summary(glht.mod, test=adjusted("fdr"))

# height differences
mod1 <- rlm(ht~geno, data=TCG)
glht.mod1 <- glht(mod1, mcp(geno = "Tukey"))
summary(glht.mod1, test=adjusted("fdr"))

# weight differences
mod2 <- rlm(wt~geno, data=TCG)
glht.mod2 <- glht(mod2, mcp(geno = "Tukey"))
summary(glht.mod2, test=adjusted("fdr"))



data2=TCG[(TCG$geno=="2"),]
data5=TCG[(TCG$geno=="5"),]
data11=TCG[(TCG$geno=="11"),]

#differences in variance
v1=var.test(data2$wt,data5$wt)$p.value
v2=var.test(data2$wt,data11$wt)$p.value
v3=var.test(data5$wt,data11$wt)$p.value
vwt=c(v1,v2,v3)
p.adjust(vwt,method="bonferroni")

v4=var.test(data2$ht,data5$ht)$p.value
v5=var.test(data2$ht,data11$ht)$p.value
v6=var.test(data5$ht,data11$ht)$p.value
vht=c(v4,v5,v6)
p.adjust(vht,method="bonferroni")

v7=var.test(data2$til,data5$til)$p.value
v8=var.test(data2$til,data11$til)$p.value
v9=var.test(data5$til,data11$til)$p.value
vtil=c(v7,v8,v9)
p.adjust(vtil,method="bonferroni")
