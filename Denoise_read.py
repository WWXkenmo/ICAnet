import randomly
import pandas as pd
def denoise_read(readPath,index_col=0,sep=" ",min_tp=0,min_genes_per_cell=0,min_cells_per_gene=0):
  fileIN = pd.read_table(readPath, index_col=index_col, sep=sep)
  model = randomly.Rm()
  model.preprocess(fileIN, min_tp=min_tp, 
                                    min_genes_per_cell=min_genes_per_cell, 
                                    min_cells_per_gene=min_cells_per_gene,
                                    refined=True)
  model.refining(min_trans_per_gene=min_trans_per_gene)
  model.fit()
  fileOut = model.return_cleaned(fdr=fdr)
  return fileOut
