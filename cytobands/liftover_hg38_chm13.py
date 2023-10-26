from liftover import get_lifter

converter = get_lifter("hg38", "hs1")
chrom = "1"
pos = 103786442
converter[chrom][pos]

# other synonyms for the lift call
converter.convert_coordinate(chrom, pos)


# Abandoned - ran the following
# mamba create -n liftover -c conda-forge -c bioconda ucsc-liftover
# curl -LO https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHs1.over.chain.gz
# liftOver cytoBand_HG38.tsv hg38ToHs1.over.chain cytoband_chm13.bed cytoband.unmapped
