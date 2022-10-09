library(ggplot2)
slope_isolates<-read_excel("/Users/arijitmukherjee/Documents/soil_D_stat/slope_chart_input.xlsx",sheet = "isolates",col_names = T,skip = 0)
slope_MAGs<-read_excel("/Users/arijitmukherjee/Documents/soil_D_stat/slope_chart_input.xlsx",sheet = "MAGs",col_names = T,skip = 0)
p<-newggslopegraph(slope_isolates, Dataset,D_stat,BGC,
                   Title = "",SubTitle = "",Caption = "",
                   YTextSize = 4,DataTextSize = 3.5,LineThickness = 1.5)
p

q<-newggslopegraph(slope_MAGs, Dataset,D_stat,BGC,Title = "",SubTitle = "",Caption = "",
                   YTextSize = 4,DataTextSize = 3.5,LineThickness = 1.5)
q
ggsave(
  "slope_isolates/MAGs.tiff",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 8,
  height = 6,
  units = "in",
  dpi = 400,
)
dev.off()

