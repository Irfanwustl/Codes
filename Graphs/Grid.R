GSE44770$case[GSE44770$type == "1"] = "Case"
GSE44771$case[GSE44771$type == "1"] = "Case"
GSE4835017$case[GSE4835017$type == "1"] = "Case"
GSE4835018$case[GSE4835018$type == "1"] = "Case"
GSE4835019$case[GSE4835019$type == "1"] = "Case"
GSE4835020$case[GSE4835020$type == "1"] = "Case"
GSE44770$case[GSE44770$type == "0"] = "Control"
GSE44771$case[GSE44771$type == "0"] = "Control"
GSE4835017$case[GSE4835017$type == "0"] = "Control"
GSE4835018$case[GSE4835018$type == "0"] = "Control"
GSE4835019$case[GSE4835019$type == "0"] = "Control"
GSE4835020$case[GSE4835020$type == "0"] = "Control"


G1 <- ggplot(GSE4835017, aes(x = case, y = IL6))+geom_point(size=0.5)+labs(x="",y="GE_Scaled",title = "IL6",subtitle = "P-value = 0.3545")
G2 <- ggplot(GSE4835017, aes(x = case, y = HAMP))+geom_point(size=0.5)+labs(x="",y="GE_Scaled",title = "HAMP",subtitle ="P-value = 0.4078")
G3 <- ggplot(GSE4835017, aes(x = case, y = SLC39A14))+geom_point(size=0.5)+labs(x="",y="GE_Scaled",title = "SLC39A14",subtitle ="P-value = 0.07916")
G4 <- ggplot(GSE4835017, aes(x = case, y = TF))+geom_point(size=0.5)+labs(x="",y="GE_Scaled",title = "TF",subtitle ="P-value = 0.1374")
G5 <- ggplot(GSE4835017, aes(x = case, y = TFRC))+geom_point(size=0.5)+labs(x="",y="GE_Scaled",title = "TFRC",subtitle ="P-value = 0.008064")
G6 <- ggplot(GSE4835017, aes(x = case, y = SLC11A2))+geom_point(size=0.5)+labs(x="",y="GE_Scaled",title = "SLC11A2",subtitle ="P-value = 0.4476")
G7 <- ggplot(GSE4835017, aes(x = case, y = SLC40A1))+geom_point(size=0.5)+labs(x="",y="GE_Scaled",title = "SLC40A1",subtitle ="P-value = 0.5071")
G8 <- ggplot(GSE4835017, aes(x = case, y = SLC39A8))+geom_point(size=0.5)+labs(x="",y="GE_Scaled",title = "SLC39A8",subtitle ="P-value = 0.5252")
G9 <- ggplot(GSE4835017, aes(x = case, y = FTL))+geom_point(size=0.5)+labs(x="",y="GE_Scaled",title = "FTL",subtitle ="P-value = 0.7856")
grid.arrange(G1,G2,G3,G4,G5,G6,G7,G8,G9, top = "GSE4835017")


#pushViewport(viewport(layout = grid.layout(4,3)))
#grid.text("GSE4835017", vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
#print(G1 + labs(title = "IL6",subtitle = "P-value = 4.389e-10"),
#  vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
#print(G2 + labs(title = "HAMP",subtitle ="P-value < 2.2e-16"),
#  vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
#print(G3 + labs(title = "SLC39A14",subtitle ="P-value 2.506e-08"),
#  vp = viewport(layout.pos.row = 2, layout.pos.col = 3))
#print(G4 + labs(title = "TF",subtitle ="P-value = 0.02687"),
#  vp = viewport(layout.pos.row = 3, layout.pos.col = 1))
#print(G5 + labs(title = "TFRC",subtitle ="P-value = 0.2372"),
#  vp = viewport(layout.pos.row = 3, layout.pos.col = 2))
#print(G6 + labs(title = "SLC11A2",subtitle ="P-value = 0.001357"),
#  vp = viewport(layout.pos.row = 3, layout.pos.col = 3))
#print(G7 + labs(title = "SLC40A1",subtitle ="P-value = 0.4289"),
#  vp = viewport(layout.pos.row = 4, layout.pos.col = 1))
#print(G8 + labs(title = "SLC39A8",subtitle ="P-value = 0.007142"),
#  vp = viewport(layout.pos.row = 4, layout.pos.col = 2))
#print(G9 + labs(title = "FTL",subtitle ="P-value = 0.01823"),
#  vp = viewport(layout.pos.row = 4, layout.pos.col = 3))
