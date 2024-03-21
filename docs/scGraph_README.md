# scGraph
ScGraph is a GNN-based automatic cell identification algorithm leveraging gene interaction relationships to enhance the performance of the cell type identification.

## Requirements
- python = 3.6.7
- pytorch = 1.1.0
- pytorch-geometric = 1.3.1
- sklearn
  

## Instructions
 
### Data Preprocessing
#### Uncompress data
    cd ./data/origin
    unzip expr.zip -d ./scGraph/data
    unzip origiSTRINGDB.graph.csv.zip -d ./scGraph/data
    # get 3 files: expr.label.subset.csv,  expr.mat.subset.csv,  STRINGDB.graph.csv

#### View unzipped files

        expr_mat_file: scRNA-seq expression matrix with genes as rows and cells as columns (csv format)
            e.g.    ————————————————————————————————————————————————————
                    geneID, barocode1,barocode2,barocode3,barocode3 ...
                    ————————————————————————————————————————————————————
                    5685, 1,0,0,0 ... 
                    5692, 0,0,0,0 ...
                    6193, 0,0,0,1 ...
        
        network_backbone_file: gene interactin network backbone (csv format)
            e.g.    STRING database,1~3 cloumns indicate gene1, gene2 and combined_score respectively. Genes are in Entrez ID format.
                    —————————————————————————————
                    gene1, gene2, combined_score
                    —————————————————————————————
                    23521,6193,999
                    5692,5685,999
                    5591,2547,999
                    6222,25873,999

        expr_label_file: gene interactin network backbone (csv format)
            e.g.    ————————————————
                    Barcodes, label
                    ————————————————
                    TTH233-160909-ZL  CD4_C04-TCF7
                    NTY112-20170825   CD4_C04-TCF7
                    NTH96-0909-ZL     CD4_C04-TCF7
  
         
#### Preprocessing data for model training
    python gen_data.py  -expr ./data/expr.mat.subset.csv \
                        -net ./data/STRINGDB.graph.csv \
                        -label ./data/expr.label.subset.csv \
                        -out ./data/dataset.npz \

    Arguments:
        - expr: scRNA-seq expression matrix with genes as rows and cells as columns (csv format)
        - net: gene interactin network backbone (csv format)
        - label: gene interactin network backbone (csv format)
        - out: preprocessed data for model training (npz format)

    Optional Arguments:
        - q: <float> the top q quantile of network edges are used (default: 0.99 for STRING database)


### Run scGraph model 
    python scGraph.py   -in ./data/dataset.npz \
                        -out-dir ./results

    Arguments:
        - in: preprocessed data for model training (npz format)  
        - out-dir: the folder in which prediction results are saved 

    Optional Arguments:
        - batch_size : batch size for model training
