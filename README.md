# Topic 2 : The role of tissue-specific antigens in different cancer entities

**Supervisor:**

* Dr. Maria Dinkelacker (m.dinkelacker@dkfz.de)

**Introduction:**

Tissue restricted antigens (TRAs) are genes, which are highly expressed in a few tissues, but not in others. These genes are up-regulated in the medullary part of the thymus (mTECs) in order to teach developing T-cells in a process called negative selection not to react to them. These T cells are out selected by apoptosis. Only 5% of all T cells survive these selection procedures (positive selection and negative selection). While they are tested in the positive selection on self-MHC they are tested in the negative selection on self-antigens-self MHC.
If these selection procedures don’t work well, human suffer of multiple autoimmune diseases. [Kyewski et al. 2004]. For some reason these genes are often up-regulated in cancer, although this atopic gene expression should not take place and can due to the negative selection of T cells in the thymus not react to these genes. The tumor down regulates the MHC molecules and escapes the immune system. 

For this reason tissue-restricted antigens (TRAs) are potentially good drug targets both in cancer therapy as also in cancer immunotherapy. Among these genes, there are several gene groups, which have been already in focus of cancer immunotherapy, such as testis-specific antigens (CTAs), but also oocyte-specific genes, others are still very unknown, such as kallikrein genes, casein genes, skin-specific genes, pancreas-specific genes, and the role of TRAs in developmental processes in general.

**Literature:** 

Kyewski et al. 2004, Self-representation in the thymus: an extended view, Nat. Rev. Immunol. 2004 Sep, 4 (9) 699-698. https://pubmed.ncbi.nlm.nih.gov/15343368/

**Objectives:**

Please get ***TRA data*** from Dr. Dinkelacker, 2007, 2019 (6 different datasets) Su et al. 2002, 2004 (mouse Novartis data), Roth et al. 2008 (human Roth data), Lattin et al. 2006 (mouse Lattin data), human GTEX data 2015 (RNAseq data), and protein atlas data (protein data).

Try to import these tables (.csv) files into R and get an overview of tissue-distribution, and the role of the gene group of your interest in the datasets. (Follow here the basic steps, denoted in the R course from Carl Hermann, lecture 1, 2, 3).

-	**Klk genes (Project 1) – any cancer**
-	**Csn genes (Project 2) – any cancer**
-	**Skin spec. genes (Project 3) – skin cancer**
-	**Pancreas spec. genes (Project 4) – pancreas CA**
-	**Thyroid spec. genes (Project 5) – thyroid cancer**

 Piechart of distribution of tissue-type (all TRAs)
 
 Piechart of distribution of tissue-type (your gene group of interest, only Klk, Csn genes)
 
 Barplot of your genes of interest per chromosome (distribution of the genes of interest per chromosome, please be aware, that the X chromosome is a ***“character”***, while the others are numbers, and R does not like both datatypes in one object. Also note that the order is correct 1,2,3, not 1,10,11,2, …, use an order function or sorting vector for this)

**Data retrieval and data analysis:**

Please get the R course on ***“Tumor microenvironment”*** from Dr. Dinkelacker (Dropbox folder) and follow the instructions there, how to analyse Microarray data.

For this please select for each project a GEO dataset (different cancer entities) with Affymetrix hgu133plus2 chips and download the rawdata. If the dataset is too small, please limitate your search to the first 10 patients (20 patients) depending on your computer calculation power.

Start with the breast cancer dataset, denoted in the R course (“Tumor microenvironment”).

Do all steps, which are important for the quality control and data analysis and document your work well.

Try to plot instead of the IL genes (in the course) the gene list of your interest.

Download Ensembl biomart, annotation table, including the Affymetrix ID of your chip, with Chromosome, Gene Symbol, Ensembl Gene ID, Ensembl Transcript ID, Startsite (document the Version no.)

Install packages (vsn, affy, including dependencies, Bioconductor, Bioconductor and R packages, how to install these is described in the R course ***“Dinkelacker”***)

Boxplot of gene expression data before and after vsnrma normalization

Scatterplot, of single chip quality control

RNAdeg plot

Meansd plot (Quality control, explained in the R course material)

extract the gene group of interest out of the normalized data matrix, by extracting the rownames(normalized.data) and change them into gene symbols and re-apply this to the data matrix, check with head(normalized.data). 

Now extract from this vector via an index (rownumbers) the genes, you want to select and extract a matrix only including the genes you want to study normalized.data[ind,] this is the data matrix you want to work with.

transform this data matrix with t(data.matrix) and plot these genes with 
boxplot(transformed.data), sort these genes alphabetically and save the plot/the plots, give a title, and also title the axes (gene expression log2 – y axis), gene symbols – x axis, turn these by 90 degrees (all these steps are explained in the R course ***“Dinkelacker”***).

- Heatmap of these genes *
- Cluster analysis of your genes *
- PCA*
- k-means*
- hierarchical clustering*
- statistical test (limma analysis, t-test, F-test)*
- Venn diagram*

*If this is too difficult on the resulting data (cancer data), please apply all these methods to the TRA dataset in the beginning and compare here your genes of interest to all others, apply the methods you find suited and work.

*try to find out in the following datasets, comparisons between woman, men, ill, healthy, mutations, non mutations (depending on the dataset you used)

Then try to go and find in GEO a second dataset, different cancer entity (Affymetrix hgu133plus2 chips) and analyze this the same way.

- **For Project 1 – klk genes (any cancer type is good)**
- **For Project 2 – Csn genes (any cancer type is good)**
- **For Project 3 – Skin spec. genes (please select skin specific datasets only, melanoma, any other skin cancer)**
- **For Project 4 - Pankreas spec. genes (Pancreas CA)**
- **For Project 5 – Thyroid spec. genes (Thyroid spec. cancers)**

Proceed with this dataset of your choice (you might recheck with me, if you are in doupt), the same way as the breast cancer dataset and do the exact same analysis, use the same R script, which you saved before as session.R, also save your session.rda for your work documentation, this is very important. If the data is too big for your computer, please only use 10 chips (20 chips). The microarrays have to be downloaded from GEO (rawdata, .CEL files), and extracted twice (gunzip, .gz and untar .tar) put only the chips, you want to evaluated into an extra folder called rawdata, and read these chips in with the read.affy() function (all these steps are in the course “Dinkelacker”, task 6).

Try to document your work well in a written format as well as in a talk, in order to present the topic in the end.

If you need any help, please feel free to contact me or your tutor any time, and we will help you. Contact us early enough to help.

***Remark:*** In the R script of the course material there are some mistakes, don’t use the Brainarray packages due to R version conflicts and don’t change the variable for these packages. Please contact me, if you have problems with this, if the vsnrma normalization doesn’t work, this is usually the mistake.

**Literature review:**

For your research, please also do a literature research on the gene group of interest, the cancer type of your subtopic and look, if any differential gene expression of the gene group of interest has already been observed in cancer in general or the specific cancer in detail.

-	Please treat the TRA data confidential - 

Dr. Maria Dinkelacker,  m.dinkelacker@dkfz.de

