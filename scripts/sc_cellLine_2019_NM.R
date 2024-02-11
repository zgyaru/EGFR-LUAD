library(Seurat)
library(cowplot)

meta = read.csv('./data/external/cell_lines/2019_NM/github/sc_mixology-data-csv/sc_10x.metadata.csv.gz')
unique(meta$cell_line_demuxlet)
counts = read.csv('./data/external/cell_lines/2019_NM/github/sc_mixology-data-csv/sc_10x.count.csv.gz')
egfr = 'ENSG00000146648'
t.test(counts[egfr,which(meta$cell_line_demuxlet == 'HCC827')],
       counts[egfr,which(meta$cell_line_demuxlet == 'H1975')]
       )

meta = read.csv('./data/external/cell_lines/2019_NM/github/sc_mixology-data-csv/sc_10x_5cl.metadata.csv.gz')
unique(meta$cell_line_demuxlet)
counts = read.csv('./data/external/cell_lines/2019_NM/github/sc_mixology-data-csv/sc_10x_5cl.count.csv.gz')
egfr = 'ENSG00000146648'
egfr = 'EGFR'
t.test(counts[egfr,which(meta$cell_line_demuxlet == 'HCC827')],
       counts[egfr,which(meta$cell_line_demuxlet == 'H1975')]
)
data = CreateSeuratObject(counts)
data@meta.data = meta
data = NormalizeData(data)
data = FindVariableFeatures(data)
data = ScaleData(data)
data = RunPCA(data)
data = RunUMAP(data,dims=1:13)
DimPlot(data, group.by = 'cell_line_demuxlet')
FeaturePlot(data, features = 'EGFR',cols = c('gray90','brown3'))
FeaturePlot(data, features = egfr,cols = c('gray90','brown3'))
data=readRDS('./data/external/cell_lines/2019_NM_yes/github/sc_mixology-data-csv/sc_10x_5cl.rds')
gene = 'DMBT1'
gene = 'GNAS'
gene = 'ERBB2'
t.test(data@assays$RNA@data[gene,which(data$cell_line_demuxlet %in% c('HCC827'))],
       data@assays$RNA@data[gene,which(data$cell_line_demuxlet %in% c('H1975'))]
)
p1 = DimPlot(data, group.by = 'cell_line_demuxlet',label = T,
             cols = c('#A49375','#73A1B0','#F9D46C','#85496F','#F1956E'))+NoLegend()
p2 = FeaturePlot(data, features = gene,cols = c('gray90','brown3'))
plot_grid(p1,p2)

VlnPlot(data[,which(data$cell_line_demuxlet %in% c('HCC827','H1975'))],features = gene,
        pt.size = 0,group.by = 'cell_line_demuxlet')+
  scale_fill_manual(values = c('#73A1B0','#F1956E'))


VlnPlot(data[,which(data$cell_line_demuxlet %in% c('HCC827','H1975'))],features = 'EGFR',
        pt.size = 0,group.by = 'cell_line_demuxlet')+
  scale_fill_manual(values = c('#73A1B0','#F1956E'))

data = readRDS('./data/external/cell_lines/2019_NM_yes/github/sc_mixology-data-csv/sc_10x_5cl.rds')
gene = 'ENSG00000187908'
gene = 'FN1'
gene = 'GNAS'
t.test(data@assays$RNA@data[gene,which(data$cell_line_demuxlet %in% c('HCC827'))],
       data@assays$RNA@data[gene,which(data$cell_line_demuxlet %in% c('H1975'))]
       )
median(data@assays$RNA@data['EGFR',which(data$cell_line_demuxlet %in% c('HCC827'))])
median(data@assays$RNA@data['EGFR',which(data$cell_line_demuxlet %in% c('H1975'))])

data=readRDS('./data/external/cell_lines/2019_NM_yes/github/sc_mixology-data-csv/sc_10x_3cl.rds')
t.test(data@assays$RNA@data[gene,which(data$cell_line_demuxlet %in% c('HCC827'))],
       data@assays$RNA@data[gene,which(data$cell_line_demuxlet %in% c('H1975'))]
)

gene = 'ENSG00000087460' ## GNAS
print(sum(data$cell_line_demuxlet == 'HCC827'))
print(sum(data$cell_line_demuxlet == 'H2228'))
print(sum(data$cell_line_demuxlet == 'H1975'))
p1 = DimPlot(data, group.by = 'cell_line_demuxlet',label = T,
             cols = c('#73A1B0','#F9D46C','#F1956E'))+NoLegend()
p2 = FeaturePlot(data, features = egfr,cols = c('gray90','brown3'))
plot_grid(p1,p2)

VlnPlot(data[,which(data$cell_line_demuxlet %in% c('HCC827','H1975'))],features = gene,
        pt.size = 0,group.by = 'cell_line_demuxlet')+
  scale_fill_manual(values = c('#73A1B0','#F1956E'))


print(sum(data$cell_line_demuxlet == 'HCC827'))
print(sum(data$cell_line_demuxlet == 'H2228'))
print(sum(data$cell_line_demuxlet == 'H1975'))

meta = read.csv('./data/external/cell_lines/2019_NM/github/sc_mixology-data-csv/sc_celseq2.metadata.csv.gz')
unique(meta$cell_line_demuxlet)
counts = read.csv('./data/external/cell_lines/2019_NM/github/sc_mixology-data-csv/sc_celseq2.count.csv.gz')
data=readRDS('./data/external/cell_lines/2019_NM_yes/github/sc_mixology-data-csv/sc_celseq2_3cl.rds')
ensg = 'ENSG00000187908'
t.test(data@assays$RNA@data[ensg,which(data$cell_line_demuxlet %in% c('HCC827'))],
       data@assays$RNA@data[ensg,which(data$cell_line_demuxlet %in% c('H1975'))]
)
print(sum(data$cell_line_demuxlet == 'HCC827'))
print(sum(data$cell_line_demuxlet == 'H2228'))
print(sum(data$cell_line_demuxlet == 'H1975'))
p1 = DimPlot(data, group.by = 'cell_line_demuxlet',label = T,
             cols = c('#73A1B0','#F9D46C','#F1956E'))+NoLegend()
p2 = FeaturePlot(data, features = egfr,cols = c('gray90','brown3'))
plot_grid(p1,p2)

VlnPlot(data[,which(data$cell_line_demuxlet %in% c('HCC827','H1975'))],features = gene,
        pt.size = 0,group.by = 'cell_line_demuxlet')+
  scale_fill_manual(values = c('#73A1B0','#F1956E'))


meta = read.csv('./data/external/cell_lines/2019_NM/github/sc_mixology-data-csv/sc_celseq2_5cl_p1.metadata.csv.gz')
unique(meta$cell_line_demuxlet)
counts = read.csv('./data/external/cell_lines/2019_NM/github/sc_mixology-data-csv/sc_celseq2_5cl_p1.count.csv.gz')
p1 = DimPlot(data, group.by = 'cell_line_demuxlet',label = T,
             cols = c('#F1956E','#F9D46C','#85496F'))+NoLegend()
p2 = FeaturePlot(data, features = egfr,cols = c('gray90','brown3'))
plot_grid(p1,p2)


meta = read.csv('./data/external/cell_lines/2019_NM/github/sc_mixology-data-csv/sc_celseq2_5cl_p1.metadata.csv.gz')
unique(meta$cell_line_demuxlet)
counts = read.csv('./data/external/cell_lines/2019_NM/github/sc_mixology-data-csv/sc_celseq2_5cl_p1.count.csv.gz')
data=readRDS('./data/external/cell_lines/2019_NM_yes/github/sc_mixology-data-csv/sc_celseq2_5cl_p1.rds')
for(c in unique(data$cell_line_demuxlet)){
  print(c, str(sum(data$cell_line_demuxlet == c)))
}
t.test(data@assays$RNA@data[egfr,which(data$cell_line_demuxlet %in% c('HCC827'))],
       data@assays$RNA@data[egfr,which(data$cell_line_demuxlet %in% c('H1975'))]
)
p1 = DimPlot(data, group.by = 'cell_line_demuxlet',label = T,
             cols = c('#A49375','#73A1B0','#F9D46C','#85496F','#F1956E'))+NoLegend()
p2 = FeaturePlot(data, features = egfr,cols = c('gray90','brown3'))
plot_grid(p1,p2)

VlnPlot(data[,which(data$cell_line_demuxlet %in% c('HCC827','H1975'))],features = gene,
        pt.size = 0,group.by = 'cell_line_demuxlet')+
  scale_fill_manual(values = c('#73A1B0','#F1956E'))




meta = read.csv('./data/external/cell_lines/2019_NM/github/sc_mixology-data-csv/sc_celseq2_5cl_p2.metadata.csv.gz')
unique(meta$cell_line_demuxlet)
counts = read.csv('./data/external/cell_lines/2019_NM/github/sc_mixology-data-csv/sc_celseq2_5cl_p2.count.csv.gz')
data=readRDS('./data/external/cell_lines/2019_NM_yes/github/sc_mixology-data-csv/sc_celseq2_5cl_p2.rds')
for(c in unique(data$cell_line_demuxlet)){
  print(c, str(sum(data$cell_line_demuxlet == c)))
}
t.test(data@assays$RNA@data[egfr,which(data$cell_line_demuxlet %in% c('HCC827'))],
       data@assays$RNA@data[egfr,which(data$cell_line_demuxlet %in% c('H1975'))]
)
p1 = DimPlot(data, group.by = 'cell_line_demuxlet',label = T,
             cols = c('#A49375','#73A1B0','#F9D46C','#85496F','#F1956E'))+NoLegend()
p2 = FeaturePlot(data, features = egfr,cols = c('gray90','brown3'))
plot_grid(p1,p2)

VlnPlot(data[,which(data$cell_line_demuxlet %in% c('HCC827','H1975'))],features = gene,
        pt.size = 0,group.by = 'cell_line_demuxlet')+
  scale_fill_manual(values = c('#73A1B0','#F1956E'))


meta = read.csv('./data/external/cell_lines/2019_NM/github/sc_mixology-data-csv/sc_celseq2_5cl_p3.metadata.csv.gz')
unique(meta$cell_line_demuxlet)
counts = read.csv('./data/external/cell_lines/2019_NM/github/sc_mixology-data-csv/sc_celseq2_5cl_p3.count.csv.gz')
data=readRDS('./data/external/cell_lines/2019_NM_yes/github/sc_mixology-data-csv/sc_celseq2_5cl_p3.rds')
for(c in unique(data$cell_line_demuxlet)){
  print(c, str(sum(data$cell_line_demuxlet == c)))
}
t.test(data@assays$RNA@data[egfr,which(data$cell_line_demuxlet %in% c('HCC827'))],
       data@assays$RNA@data[egfr,which(data$cell_line_demuxlet %in% c('H1975'))]
)
p1 = DimPlot(data, group.by = 'cell_line_demuxlet',label = T,
             cols = c('#A49375','#73A1B0','#F9D46C','#85496F','#F1956E'))+NoLegend()
p2 = FeaturePlot(data, features = egfr,cols = c('gray90','brown3'))
plot_grid(p1,p2)

VlnPlot(data[,which(data$cell_line_demuxlet %in% c('HCC827','H1975'))],features = gene,
        pt.size = 0,group.by = 'cell_line_demuxlet')+
  scale_fill_manual(values = c('#73A1B0','#F1956E'))



meta = read.csv('./data/external/cell_lines/2019_NM/github/sc_mixology-data-csv/sc_dropseq.metadata.csv.gz')
unique(meta$cell_line_demuxlet)
counts = read.csv('./data/external/cell_lines/2019_NM/github/sc_mixology-data-csv/sc_dropseq.count.csv.gz')
data=readRDS('./data/external/cell_lines/2019_NM_yes/github/sc_mixology-data-csv/sc_dropseq_3cl.rds')
for(c in unique(data$cell_line_demuxlet)){
  print(c, str(sum(data$cell_line_demuxlet == c)))
}
t.test(data@assays$RNA@data[egfr,which(data$cell_line_demuxlet %in% c('HCC827'))],
       data@assays$RNA@data[egfr,which(data$cell_line_demuxlet %in% c('H1975'))]
)
p1 = DimPlot(data, group.by = 'cell_line_demuxlet',label = T,
             cols = c('#73A1B0','#F9D46C','#F1956E'))+NoLegend()
p2 = FeaturePlot(data, features = egfr,cols = c('gray90','brown3'))

plot_grid(p1,p2)



