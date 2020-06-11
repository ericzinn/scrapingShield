# scrapingShield - version 0.1.1
A basic python tool to scrape the Shield database at HMS for FACS-Sorted RNASeq data.  Includes some basic analysis and visualization capabilities

## Included Files
### shieldScraper.py
This is the executable version of the script.  By default, it will scrape charts and FACS-Sorted RNA-Seq raw-data associated with all of the genes that the user specifies (in a single-column UTF-8 encoded csv file).  It then massages the data a little bit and outputs the data as csv files.  Additionally, it renders three different kinds of heatmaps for each experimental condition in the database: a basic heatmap (unsorted), a clustered heatmap, and a log-10 transformed clustered heatmap.  The user can comment out sections of the pipeline that they're not interested in.

### Demo - Scraper Notebook.ipynb
This is jupyter (ipython) notebook which breaks down every block of code and every function in plain-english.  Available as a static html file in github in case people are interested in reading up on the code but don't have jupyter notebooks set up on their own machines.

### shortList.csv
This is a demonstration list of 10 genes which can be used as a demo of the script.  The low number makes execution REALLY fast.

### Deafness_gene_list.csv
This is a list of deafness-related genes compiled by Cole Peters.

## Usage
```bash
python shieldScraper.py -i {path/to/gene_list.csv}
```

## Contact info
For questions, please contact Eric Zinn...who isn't going to put his contact information on a public-facing website for fear of invoking the wrath of spammers...but chances are you know how to get in touch if you're reading this.
