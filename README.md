# corProp

Use the correlation between expression profiles to estimate the composition of mixed expression profile

    source('https://raw.githubusercontent.com/jumphone/corProp/main/corProp.R')

    OUT=corProp(REF,BULK)
     
    library('ComplexHeatmap')
    library('circlize')
    col_fun =colorRamp2(c(0, 2.5, 5, 10, 50, 90 ), c('royalblue','grey90','grey90','lightsalmon','indianred1','brown1'))
    
    o.mat=as.matrix(OUT)
    ht=Heatmap(o.mat,row_title='',name="%",cluster_rows=F,cluster_columns=F,
        show_column_dend = F, show_row_dend = F,
        show_column_names=T, show_row_names=T,
        col=col_fun, border = TRUE,
        row_names_side='right',
        row_names_gp = gpar(fontsize = 10, lineheight=NULL),
        column_names_gp = gpar(fontsize = 10, lineheight=NULL),
        row_names_max_width = max_text_width(rownames(o.mat),gp = gpar(fontsize = 15)),
        column_names_max_height = max_text_width(colnames(o.mat),gp = gpar(fontsize = 15)),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.1f", o.mat[i, j]), x, y, gp = gpar(fontsize = 10))
          }
        )
        
    print(ht)

    
    
