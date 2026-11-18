# Data

The raw data are not included in this repository. They are publicly available from the NCI Genomic Data Commons (GDC).

## How to download

1. Go to the [GDC Data Portal](https://portal.gdc.cancer.gov/) and open the **Repository**.
2. Filter by:
   - Project: `TCGA-COAD`
   - Data category: `transcriptome profiling`
   - Data type: `miRNA Expression Quantification`
   - Workflow: `BCGSC miRNA Profiling`
3. Add all files to the cart, then download:
   - the **files** (manifest + `gdc-client`, or direct download), and
   - the **sample sheet**.
4. Place them as follows:

```
data/
├── raw/                       # one folder per File ID, each with a
│   ├── <file-id>/             #   *.mirbase21.mirnas.quantification.txt file
│   └── ...
└── gdc_sample_sheet.tsv       # GDC sample sheet
```

The thesis analysis used 458 samples (450 primary tumour, 8 adjacent normal) and 1,881 miRNAs.
