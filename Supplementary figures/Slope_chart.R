slope<-readRDS("/Users/arijitmukherjee/Documents/Phytobiome/BGC_isolates/New_final_data/MAGs/Earth_MAGs/latest_ANI_analysis/Fig.2/Slope_chart_final.rds")
write.table(slope$data,"slope_data.tsv",sep = "\t")

slope_isolates<-read_excel("/Users/arijitmukherjee/Downloads/slope_chart_input.xlsx",sheet = "isolates",col_names = T,skip = 0)
slope_MAGs<-read_excel("/Users/arijitmukherjee/Downloads/slope_chart_input.xlsx",sheet = "MAGs",col_names = T,skip = 0)
p<-newggslopegraph(slope_isolates, Dataset,D_stat,BGC,Title = "",SubTitle = "",Caption = "")
p

q<-newggslopegraph(slope_MAGs, Dataset,D_stat,BGC,Title = "",SubTitle = "",Caption = "")
q
ggsave(
  "slope_MAGs.tiff",
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

