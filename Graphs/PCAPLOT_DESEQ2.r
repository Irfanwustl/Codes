plotPCA(rld, intgroup = "groups")+ geom_point(size = 8)+ scale_color_manual(values=c("#0547ff", "#ff0505", "#ff9000","#c4c4c4"))+theme_bw() +
+     theme(panel.grid.major = element_blank(),
+           panel.grid.minor = element_blank(),
+           strip.background = element_blank(),
+           panel.border = element_rect(colour = "black"))
