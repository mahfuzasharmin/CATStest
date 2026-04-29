USERNAME = "sharminm2"


# REPLACE THESE LOCATIONS WITH YOUR DIRECTORIES
PROJECT_location = f"CATSEvalData"


SIGNATURES = ["Angiogenesis", "Apoptosis", "Cell_Cycle", "Differentiation",
               "DNA_damage", "DNA_repair", "EMT", "Hypoxia", "Inflammation",
               "Invasion", "Metastasis", "Proliferation", "Quiescence", "Stemness"] 

CANCERTYPES_SUBTYPES = ['ACC', 'BLCA', 'BRCA', 'BRCA_Basal', 'BRCA_LumA', 'BRCA_LumB', 'BRCA_Her2', 'CESC', 'COAD', 
               'ESCA', 'GBM', 'KICH', 'KIRC', 'KIRP', 'LGG', 'LIHC', 'LUAD', 'LUSC', 'OV', 'PAAD', 'PRAD', 
               'READ', 'THCA', 'SKCM', 'STAD', 'UCEC', 'UCS']

CANCER_TYPES = ['ACC', 'BLCA', 'BRCA', 'CESC', 'COAD', 'ESCA', 'GBM', 'KICH', 'KIRC', 'KIRP', 
                'LGG', 'LIHC', 'LUAD', 'LUSC', 'OV', 'PAAD', 'PRAD', 'READ', 'THCA', 'SKCM', 
                'STAD', 'UCEC', 'UCS']

BREAST_CANCER_SUBTYPES = ['BRCA_Basal', 'BRCA_LumA', 'BRCA_LumB', 'BRCA_Her2']
BRAIN_CANCER_SUBTYPES = ['GBM', 'LGG']

CANCERTYPE_TO_TISSUE = {'ACC': 'Adrenal', 'BLCA': 'Bladder', 'BRCA': 'breast', 
                        'BRCA_Basal': 'BRCA_Basal', 'BRCA_LumA': 'BRCA_LumA', 
                        'BRCA_LumB': 'BRCA_LumB', 'BRCA_Her2': 'BRCA_Her2',
                       'CESC':'cervix', 'COAD':'colon', 'ESCA':'Esophagus', 'GBM':'brain', 'Brain':'brain_all', 
                       'KICH':'kidney', 'KIRC':'kidney', 'KIRP':'kidney', 'LGG':'brain', 'LIHC':'liver', 
                       'LUAD':'lung', 'LUSC':'lung', 'OV':'ovary', 'PAAD':'pancreas', 'PRAD':'Prostate', 
                       'READ':'colon', 'SKCM':'skin' , 'STAD':'Stomach', 'THCA':'thyroid', 'UCEC':'ovary', 'UCS':'ovary'}

#'adrenal gland', 'breast', 'cervix', 'colon', 'esophagus', 'kidney', 
            # 'liver', 'lung', 'ovary', 'pancreas', 'prostate', 'rectum', 'skin', 
            # 'stomach', 'thyroid gland'
CANCERTYPE_TO_HPA_TISSUE = {'ACC': 'adrenal gland', 'BRCA': 'breast', 
                       'CESC':'cervix', 'COAD':'colon', 'ESCA':'esophagus', 
                       'KICH':'kidney', 'KIRC':'kidney', 'KIRP':'kidney',  'LIHC':'liver', 
                       'LUAD':'lung', 'LUSC':'lung', 'OV':'ovary', 'PAAD':'pancreas', 'PRAD':'prostate', 
                       'READ':'rectum', 'SKCM':'skin' , 'STAD':'stomach', 'THCA':'thyroid gland'}


TCGA_TO_DepMap_CANCERTYPE = {
        #'ACC': {'Adrenal Cancer': ['Carcinoma']},
        'BLCA': {'Bladder Cancer': ['Carcinoma']},
        'BRCA': {'Breast Cancer': ['Breast Ductal Carcinoma']}, # Breast invasive carcinoma
    'BRCA_Basal': {'Breast Cancer': ['Breast Ductal Carcinoma']}, 
    'BRCA_LumA': {'Breast Cancer': ['Breast Ductal Carcinoma']}, 
    'BRCA_LumB': {'Breast Cancer': ['Breast Ductal Carcinoma']},
    'BRCA_Her2': {'Breast Cancer': ['Breast Ductal Carcinoma']}, 
    
        'CESC': {'Cervical Cancer': ['Squamous Cell Carcinoma',  'Cervical Adenocarcinoma', 'Endocervical Carcinoma']},
        
        #'COAD': {'Colon/Colorectal Cancer': ['Colorectal Carcinoma', 'Colon Carcinoma', 'Colon Adenocarcinoma', 'Caecum Adenocarcinoma']},
        'ESCA':  {'Esophageal Cancer': ['Squamous Cell Carcinoma']},
        'GBM': {'Brain Cancer': ['Glioblastoma', 'Glioma']},
        'Brain': {'Brain Cancer': ['Glioblastoma', 'Glioma', 'meningioma', 'medulloblastoma']},
        'KICH': {'Kidney Cancer': ['Renal Cell Carcinoma']}, # Kidney Chromophobe
        'KIRC': {'Kidney Cancer': ['Renal Carcinoma, clear cell']}, # Kidney renal clear cell carcinoma
        'KIRP': {'Kidney Cancer': ['Renal Cell Carcinoma']}, # Kidney renal papillary cell carcinoma
        'LGG': {'Brain Cancer': ['Glioma']}, # Brain Lower Grade Glioma
        'LIHC': {'Liver Cancer': ['Hepatocellular Carcinoma']}, # Liver hepatocellular carcinoma
        'LUSC': {'Lung Cancer': ['Non-Small Cell Lung Cancer (NSCLC), Squamous Cell Carcinoma']}, # Lung squamous cell carcinoma
        'OV':  {'Ovarian Cancer': ['Small Cell Carcinoma of the Ovary Hypercalcemic Type (SCCOHT)', 
                                       'Cystadenocarcinoma', 'Cystadenocarcinoma, high grade serous',
                                       'Cystadenocarcinoma, clear cell', 
                                       'Cystadenocarcinoma, mucinous', 'Cystadenocarcinoma, endometrioid',
                                       'Carcinoma, brenner', 'Krukenberg Tumor']},
        #'PAAD': {'Pancreatic Cancer': ['Ductal Adenocarcinoma']},
        #'PRAD':  {'Prostate Cancer': ['Prostate Hyperplasia']},
        #'READ': {'Rectum': ['']},
        
        'SKCM': {'Skin Cancer': ['Melanoma']}, # Skin Cutaneous Melanoma
        
        #'STAD': {'Gastric Cancer': ['']},
        #'UCEC':  {'Uterine Corpus Endometrial Ca': ['']},
        #'UCS': {'Uterine Carcinosarcoma': ['']}
    }

CANCERTYPE_to_GWAS = {'PRAD': 'Prostate cancer', 'BRCA': 'Breast cancer', 'COAD': 'Colorectal cancer',
                          'LGG': 'Glioma', 'GBM': 'Glioma', 'LUSC': 'Lung cancer', 'LUAD': 'Lung cancer', 
                          'OV': 'Ovarian cancer', 'SKCM': 'Melanoma', 'PAAD': 'Pancreatic cancer',
                          'ESCA': 'Esophageal cancer',  'BLCA': 'Bladder cancer', 'THCA': 'Thyroid cancer'}