**Objective**: To use different screening techniques to distinguish patients with high risks from those with low risks for colon cancer recurrence and predicting their prognosis.

**Methods**: Colon cancer samples were collected from the Gene Expression Omnibus database along with The Cancer Genome Atlas. Data from GSE17537 were analyzed using the LIMMA method to identify the differentially expressed genes (DEGs). A PPI network was created and support vector machine (SVM) analysis was performed on the DEGs for screen for the feature genes. 

**Results**: A total of 1207(109 in our analysis) genes were identified as DEGs between recurrence and no-recurrence samples. Using SVM analysis and five gene expression profile data confirmation, a 15-gene signature were identified as a predictor of recurrence risk and prognosis for colon cancer patients. Our results differed from the results of the original research. We question the repeatability of the study. 

**Conclusion**: Xu et al., identified a 15-gene signature that may be useful to classify colon cancer patients with different prognosis and some genes in the their signature may represent new therapeutic targets. The SVM model we generated has an error rate of 48%. Based on the use case in trying to differentiate recurrence for colon cancer patients this may be a useful tool to use among others.

**folder**

| **folder** | **desc**|
|------|----------|
| code | major code |
| data | macroarray data |
| document | the paper and related material |
