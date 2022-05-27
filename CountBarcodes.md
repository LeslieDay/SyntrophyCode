This script is for finding primer indexes and barcodes within a single read 
  - finds sequences that have both a index and barcode with 29-32 nucleotides between
  - checks for minimum PHRED quality based on input cutoff (recommended between 20-30)
    - PHRED quality must be met for both the index sequence and the barcode sequence to be counted
  - Counts for each index/barcode combo are written to OUTPUT FILE ID _BarcodeCounts.csv 
    - OUTPUT FILE ID argument provided with input to identify sample analyzed 

Default parameters use for Barcode ratio analysis 
```bash
python CountBarcodes.py LD012_S80_R1_001.fastq 10 index.txt barcodes.txt LD012
```

Script saved as CountBarcodes.py

```python
import sys
from Bio import SeqIO
import regex
import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("[FASTQ FILE]", help='R1_file.fastq- Fastq file of forward reads')
parser.add_argument("[PHRED QUALITY]", help='minimum PHRED quality for barcode sequence')
parser.add_argument("[INDEX FILE]",
                    help='two column tab delimited file with primer index name and barcode sequence')
parser.add_argument("[BARCODE FILE]",
                    help='two column tab delimited file with barcode name and barcode sequence')
parser.add_argument("[OUTPUT FILE ID]",
                    help='identifier added to each output file to associate with input fastq')
args = parser.parse_args()

outputID = sys.argv[5]
indexFile = pd.read_csv(sys.argv[3], sep='\t', names=['indexName', 'indexSequence'])  # save index info as a data frame
print("Indexes to be found:", "\n", indexFile)  # visual check imported correctly
dfBarcode = pd.read_csv(sys.argv[4], sep='\t', names=['barcodeName', 'barcodeSequence'])
print("Barcodes to be found:", "\n", dfBarcode)  # visual check imported correctly

countList = []
for index, row in indexFile.iterrows():
    count_idx: int = 0
    qualityIndex: int = 0
    IndexName = row['indexName']
    nameFile = IndexName + outputID + '_qualityIndex.fastq'
    for barcode, line in dfBarcode.iterrows():
        qualityBarcode: int = 0
        searchFor = (row['indexSequence']) + ".{29,32}" + line['barcodeSequence']
        print("Finding...", "\n",
              "Index: ", row['indexName'], "\n",
              "Index Sequence: ", row['indexSequence'], "\n",
              "Barcode: ", line['barcodeName'], "\n",
              "Barcode Sequence: ", line['barcodeSequence'])
        with open(sys.argv[1], "r") as fastqFile:
            for record in SeqIO.parse(fastqFile, "fastq"):  # parse the provided fastq file and search for each index
                idx = regex.search(searchFor, str(record.seq))  # index must match 100%
                if idx:
                    count_idx += 1
                    indexLength: int = len(row['indexSequence'])
                    IndexQualityData = record[idx.start(0):idx.start(0) + indexLength]
                    if min(IndexQualityData.letter_annotations["phred_quality"]) >= int(sys.argv[2]):
                        qualityIndex += 1  # only count/keep reads with PHRED >= provided phred score
                        barcodeLength: int = len(line['barcodeSequence'])
                        BarcodeQualityData = record[idx.end(0) - barcodeLength:idx.end(0)]
                        if min(BarcodeQualityData.letter_annotations["phred_quality"]) >= int(sys.argv[2]):
                            qualityBarcode += 1
        countList.append([row['indexName'], line['barcodeName'], qualityBarcode])
        print("Counted: index name, barcode name, counts:", countList)

dfCounts = pd.DataFrame(countList, columns=['index', 'barcode', 'counts'])
dfSorted = dfCounts.sort_values(['index', 'barcode'], ignore_index=True)
CountFileName = outputID + "_BarcodeCounts.csv"
dfSorted.to_csv(CountFileName, index=False)
print("Wrote quality index and barcode counts to ", CountFileName, " BarcodeCounts.csv")
print(dfSorted)
```
