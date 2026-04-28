library(cowplot)
library(ggplot2)
library(ggpubr)


contextfile = "/Users/sharminm2/OneDrive - National Institutes of Health/context/data/results/hPPIN_DE2/share_with_gulden/HR_summary_context_panel.csv"
defile = "/Users/sharminm2/OneDrive - National Institutes of Health/context/data/results/hPPIN_DE2/share_with_gulden/HR_summary_de_panel.csv"
hpafile = "/Users/sharminm2/OneDrive - National Institutes of Health/context/data/results/hPPIN_DE2/share_with_gulden/HR_summary_hpa_panel.csv"
disfile = "/Users/sharminm2/OneDrive - National Institutes of Health/context/data/results/hPPIN_DE2/share_with_gulden/HR_summary_discarded_genes_panel.csv"


HR_summary_context_panel.PFI <- read.csv(hpafile, row.names=1)

dat = data.frame()
count = 1
for(i in 1:23){
  for(j in 1:14){
    dat[count, 1] = rownames(HR_summary_context_panel.PFI)[i]
    dat[count, 2] = colnames(HR_summary_context_panel.PFI)[j]
    dat[count, 3] = HR_summary_context_panel.PFI[i,j]
    count = count + 1
  }
}
dat = na.omit(dat)
ggplot(dat, aes(x = V1, y = V2, fill = V3)) +geom_tile(color = "white",linetype = 1) +coord_fixed()+ theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 16), axis.text.y = element_text(size =16), legend.position = "top")+theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+scale_fill_gradient2(low = "#000677",mid = "#ffffff", high = "#0b6d03", limits = c(-0.7, 0.7))
ggsave2(glue('{PROJECT_location}/plots/{CUR_DST}/Fig_5SB_survival.png'))


ggboxplot(dat,y = 'V3', x = 'V2',fill = 'V2', alpha = 0.5, ggtheme = theme_minimal(),  ylab = 'HR positive')+font("xy.text")+font('ylab', color = 'black')+font('xlab', color = 'black')+rremove('xlab')+rremove('legend')+rotate_x_text(angle = 90)+theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 14), axis.text.y = element_text(size =14))+xlab('')+rotate()+theme(axis.text.y=element_blank())+ geom_hline(yintercept=0, linetype="dashed", color = "red", lwd = 1)+ylim(-0.1, 0.1)

ggboxplot(dat,y = 'V3', x = 'V1',fill = 'V1', alpha = 0.5, ggtheme = theme_minimal(), ylab = 'HR positive')+font("xy.text")+font('xlab', color = 'black', size = 14)+font('ylab', color = 'black', size = 14)+rremove('legend')+rotate_x_text(angle = 90)+theme(axis.text.x=element_blank())+theme(axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=1, size = 14))+xlab('')+ geom_hline(yintercept=0, linetype="dashed", color = "red", lwd = 1)

