from __future__ import division
from itertools import izip
from functools import partial

import numpy as np
from scipy.stats import hypergeom
import pandas as pd


mm9GOFile = "/nas3/lovci/projects/GO/mm9.ENSG_to_GO.txt.gz"
hg19GOFile = "/nas3/lovci/projects/GO/hg19.ENSG_to_GO.txt.gz"
ce10GOFile = "/nas3/lovci/projects/GO/ce10.ENSG_to_GO.txt.gz"


def hypergeometric(row, lenAllGenes, lenTheseGenes):
        if (row['inBoth'] <= 3) or (row['expressedGOGenes'] < 5):
            return np.nan
        else:
            return hypergeom.sf(row['inBoth'], lenAllGenes, row['expressedGOGenes'], lenTheseGenes)


class GO(object):
    
    def __init__(self, GOFile):
        self.GO_to_ENSG = pd.read_table(GOFile, compression="gzip").dropna()

        GO, allGenes = self._generateOntology()
        self.GO = GO
        self.allGenes = allGenes
        self.gene_id_to_name = dict(izip(self.GO_to_ENSG['Ensembl Gene ID'], self.GO_to_ENSG['Associated Gene Name']))

    def enrichment(self, geneList, background=None):
        if background is None:
            background = self.allGenes
        return self.GO_enrichment(geneList, expressedGenes= background)

    def _generateOntology(self):
        """
        :return: returns dict of ontologeis and all genes in all go ontologeis as a background
        """
        allGenesInOntologies = set(self.GO_to_ENSG['Ensembl Gene ID'])
        ontology = self.GO_to_ENSG.groupby("GO Term Accession")
        ontology = ontology.aggregate(lambda x: set(x))
        ontology['nGenes'] = ontology['Ensembl Gene ID'].apply(len)

        return ontology, allGenesInOntologies



    def GO_enrichment(self, geneList, expressedGenes = None):
        geneList  = set(list(geneList))
        expressedGenes = set(list(expressedGenes))

        lenAllGenes  = len(expressedGenes)
        lenTheseGenes = len(geneList)

        df = self.GO.copy()

        if lenTheseGenes > lenAllGenes:
            raise ValueError("Length of genes examined should not be larger than the total number of genes in organism")

        df['inBoth'] = df['Ensembl Gene ID'].apply(lambda x: len(x & geneList))
        df['expressedGOGenes'] = df['Ensembl Gene ID'].apply(lambda x: len(x & expressedGenes))

        hypergeometric_partial = partial(hypergeometric,
                                         lenAllGenes=lenAllGenes,
                                         lenTheseGenes=lenTheseGenes)

        df['Hypergeometric p-Value'] = df.apply(hypergeometric_partial, axis=1)
        num_tests = len(df['Hypergeometric p-Value'].dropna())
        df['Bonferroni-corrected Hypergeometric p-Value'] = df['Hypergeometric p-Value'].apply(lambda x: min(x * num_tests, 1)).dropna()

        #Compute various value for backwards compatabality
        df['Ensembl Gene IDs in List'] = df['Ensembl Gene ID'].apply(lambda x: x & geneList)
        df['Gene symbols in List'] = df['Ensembl Gene IDs in List'].apply(lambda x: {self.gene_id_to_name[gene_id] for gene_id in x})

        df['Ensembl Gene IDs in List'] = df['Ensembl Gene IDs in List'].apply(",".join)
        df['Gene symbols in List'] = df['Gene symbols in List'].apply(",".join)
        df['GO Term Description'] = df['GO Term Name'].apply(",".join)
        df['GO domain'] = df['GO domain'].apply(",".join)
        df['N Genes in GO category'] = df['Ensembl Gene ID'].apply(len)

        #Rename stuff for backwards compatabality
        df['N Expressed Genes in GO Category'] = df['expressedGOGenes']
        df['N Genes in List and GO Category'] = df['inBoth']

        #Sort
        df = df.sort('Bonferroni-corrected Hypergeometric p-Value')

        #Reorder for presentation purposes
        df = df[[
            'GO Term Description',
            'Bonferroni-corrected Hypergeometric p-Value',
            'N Genes in List and GO Category',
            'N Expressed Genes in GO Category',
            'N Genes in GO category',
            'Ensembl Gene IDs in List',
            'Gene symbols in List',
            'GO domain',
            'Hypergeometric p-Value',
           ]]
        return df


class hg19GO(GO):
    def __init__(self, *args, **kwargs):
        super(hg19GO, self).__init__(hg19GOFile, *args, **kwargs)

class mm9GO(GO):
    def __init__(self, *args, **kwargs):
        super(mm9GO, self).__init__(mm9GOFile, *args, **kwargs)

class ce10GO(GO):
    def __init__(self, *args, **kwargs):
        super(ce10GO, self).__init__(ce10GOFile, *args, **kwargs)

#from yan(gene symbols) -> mouse gene id
neuro_genes_mouse = """ENSMUSG00000020932
ENSMUSG00000030310
ENSMUSG00000057182
ENSMUSG00000021609        
ENSMUSG00000020838
ENSMUSG00000027168
ENSMUSG00000004872
ENSMUSG00000041607
ENSMUSG00000031144
ENSMUSG00000005360
ENSMUSG00000024406
ENSMUSG00000041309
ENSMUSG00000028736
ENSMUSG00000018411
ENSMUSG00000020886
ENSMUSG00000039830
ENSMUSG00000024935
ENSMUSG00000062380
ENSMUSG00000007946
ENSMUSG00000033006
ENSMUSG00000038331
ENSMUSG00000005917
ENSMUSG00000032126
ENSMUSG00000031217
ENSMUSG00000021848
ENSMUSG00000059003
ENSMUSG00000032446
ENSMUSG00000035033
ENSMUSG00000020052
ENSMUSG00000048450
ENSMUSG00000028280
ENSMUSG00000020950
ENSMUSG00000001566
ENSMUSG00000057880
ENSMUSG00000042453
ENSMUSG00000042589
ENSMUSG00000064329
ENSMUSG00000023328
ENSMUSG00000022705
ENSMUSG00000079994
ENSMUSG00000032259
ENSMUSG00000000214
ENSMUSG00000074637
ENSMUSG00000047976
ENSMUSG00000024304
ENSMUSG00000037771
ENSMUSG00000048251
ENSMUSG00000000247
ENSMUSG00000004151
ENSMUSG00000029595
ENSMUSG00000022952
ENSMUSG00000031285
ENSMUSG00000019230
ENSMUSG00000027273
ENSMUSG00000022055
ENSMUSG00000070880
ENSMUSG00000070691
ENSMUSG00000029580
ENSMUSG00000004891
ENSMUSG00000031425
ENSMUSG00000020262
ENSMUSG00000025037
ENSMUSG00000026959
ENSMUSG00000073640
ENSMUSG00000063316
ENSMUSG00000027967
ENSMUSG00000001018
ENSMUSG00000042258
ENSMUSG00000030067
ENSMUSG00000035187
ENSMUSG00000029563
ENSMUSG00000070570
ENSMUSG00000045994
ENSMUSG00000032318
ENSMUSG00000030516
ENSMUSG00000030500
ENSMUSG00000033208
ENSMUSG00000033981
ENSMUSG00000028936
ENSMUSG00000020524
ENSMUSG00000052915
ENSMUSG00000023945
ENSMUSG00000019874
ENSMUSG00000062070
ENSMUSG00000027951""".split("\n")

neuro_genes_human = """ENSG00000131095
ENSG00000157103
ENSG00000153253
ENSG00000142319
ENSG00000108576
ENSG00000007372
ENSG00000135903
ENSG00000197971
ENSG00000102003
ENSG00000079215
ENSG00000204531
ENSG00000206349
ENSG00000206454
ENSG00000148826
ENSG00000009709
ENSG00000186868
ENSG00000132535
ENSG00000205927
ENSG00000106688
ENSG00000198211
ENSG00000165462
ENSG00000100146
ENSG00000119042
ENSG00000115507
ENSG00000149397
ENSG00000090776
ENSG00000165588
ENSG00000136531
ENSG00000183454
ENSG00000163508
ENSG00000136535
ENSG00000139352
ENSG00000163132
ENSG00000146276
ENSG00000176165
ENSG00000130675
ENSG00000183044
ENSG00000189056
ENSG00000104267
ENSG00000111249
ENSG00000184845
ENSG00000144285
ENSG00000087085
ENSG00000151577
ENSG00000182968
ENSG00000149295
ENSG00000180176
ENSG00000181449
ENSG00000111262
ENSG00000170558
ENSG00000101438
ENSG00000127152
ENSG00000106689
ENSG00000006468
ENSG00000089116
ENSG00000159216
ENSG00000077279
ENSG00000143355
ENSG00000132639
ENSG00000104725
ENSG00000128683
ENSG00000020633
ENSG00000075624
ENSG00000132688
ENSG00000123560
ENSG00000197381
ENSG00000189221
ENSG00000176884
ENSG00000131469
ENSG00000178403
ENSG00000143553
ENSG00000016082
ENSG00000114861
ENSG00000163623
ENSG00000128573
ENSG00000104888
ENSG00000109956
ENSG00000159556
ENSG00000104067
ENSG00000091664
ENSG00000160307
ENSG00000120251
ENSG00000116251
ENSG00000155511
ENSG00000188895
ENSG00000078018
ENSG00000187714
ENSG00000115665
ENSG00000164434
ENSG00000102144
ENSG00000160710""".split("\n")
