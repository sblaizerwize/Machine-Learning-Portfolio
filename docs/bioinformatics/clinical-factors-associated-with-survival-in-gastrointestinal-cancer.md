# Clinical Factors Associated with Survival in Gastrointestinal Cancer

by [Sefora Conti](https://es.linkedin.com/in/sefora-conti-6b313047) and [Salomon Marquez](https://www.linkedin.com/in/sblaizer/)

15/06/2025 

## 1 Objective 

Gastrointestinal cancer encompasses a group of malignant and heterogeneous neoplasms that affect different organs of the digestive system, including the biliary tract, the stomach, the pancreas, the colon, and the rectum. These diseases represent a significant proportion of the global cancer burden and often have unfavorable prognoses (1, 2). Identifying clinico-pathological factors that influence patients’ overall survival is essential for improving risk stratification, therapeutic decision-making, and the design of personalized follow-up strategies.

In recent years, several approaches have been developed to improve this prediction, including in vitro models (such as patient-derived organoids)(3–5), in vivo models (for example, xenografts in animal models)(6), and in silico strategies(7–9). The latter rely on different types of data, such as histopathological images(10), magnetic resonance imaging(11), or clinical information, to generate tools that help personalize cancer treatment.

In the present work, we analyze a publicly available clinical dataset from the [TCGA](https://portal.gdc.cancer.gov/) (The Cancer Genome Atlas Program) repository using R, with the aim of identifying potential clinical factors associated with the survival of patients with gastrointestinal tumors. This analysis will explore possible correlations between pathological and demographic variables and overall survival, with the aim of providing evidence that improves patient stratification.

The data used are open access, anonymized, and their use for academic purposes is permitted under the policies of the [GDC Data Portal](https://gdc.cancer.gov/about-gdc/gdc-policies).

Visit the :simple-github: [repository of the project](https://github.com/sblaizerwize/My-Professional-Portfolio/tree/main/docs/bioinformatics/clinical-factors-associated-with-survival-in-gastrointestinal-cancer) to check out:  
- :fontawesome-solid-file-csv: TCGA data  
- :octicons-code-16: R code  
- :material-application-outline: Shiny application  

---

## 2 Data Preparation

### 2.1 Initial Dataset Exploration

The following table shows a general description of the TCGA dataset. The following sections show how this information was obtained.

| Key Points | Description |
|----|----|
| **File Type and Data** |The dataset was imported from a semicolon-delimited .csv file, containing clinical and pathological information of cancer patients. It includes variables such as age, cancer type, presence of venous or perineural invasion, microsatellite instability, and follow-up variables such as days to last contact or death. Due to its combined origin (different clinical sources within TCGA), the dataset requires considerable preprocessing to make it suitable for analysis.|
| **Available Variables** |Categorical variables (e.g., cancer type, sex) and continuous variables (age, weight, height, etc.)|
| **Dataset Dimensions** |The dataset includes 154 clinical variables from 1,308 patients, of which we have selected 34 variables considered most relevant for the study objective.|
| **Detection of Missing Values** | There are 135,291 missing values in the original dataset and 21,007 in the new dataset (with the 34 selected variables).|
| **Inconsistencies** |Typical inconsistencies in clinical data were detected, such as mismatches between vital status and follow-up or death dates. To calculate overall survival, it was necessary to create a new survival time variable (`os_time`) combining `days_to_death` and `days_to_last_followup`, as well as a binary event variable (`os_event`).|

``` r
# Cargar el conjunto de datos (csv)
df_TCGA <- read.csv("Combined_TCGA_Clinical_300525.csv",header=TRUE, sep = ";")
# Visualizar registros
head(df_TCGA[,1:5], n=3) #mostramos los primeros 3 registros
```

    ##   bcr_patient_barcode additional_studies tumor_tissue_site
    ## 1        TCGA-3L-AA1B               <NA>             Colon
    ## 2        TCGA-4N-A93T               <NA>             Colon
    ## 3        TCGA-4T-AA8H               <NA>             Colon
    ##               histological_type other_dx
    ## 1          Colon Adenocarcinoma       No
    ## 2          Colon Adenocarcinoma       No
    ## 3 Colon Mucinous Adenocarcinoma       No

``` r
tail(df_TCGA[,1:5], n=3) #mostramos los ultimos 3 registros
```

    ##      bcr_patient_barcode additional_studies tumor_tissue_site
    ## 1306        TCGA-Z5-AAPL               <NA>          Pancreas
    ## 1307                <NA>               <NA>              <NA>
    ## 1308                                                         
    ##                        histological_type other_dx
    ## 1306 Pancreas-Adenocarcinoma Ductal Type       No
    ## 1307                                <NA>     <NA>
    ## 1308

``` r
# Obtener el número de variables y sus nombres
dim(df_TCGA) #dimensiones de df
```

    ## [1] 1308  154

``` r
names(df_TCGA) #nombres de las variables
```

    ##   [1] "bcr_patient_barcode"                                       
    ##   [2] "additional_studies"                                        
    ##   [3] "tumor_tissue_site"                                         
    ##   [4] "histological_type"                                         
    ##   [5] "other_dx"                                                  
    ##   [6] "gender"                                                    
    ##   [7] "vital_status"                                              
    ##   [8] "days_to_birth"                                             
    ##   [9] "days_to_last_known_alive"                                  
    ##  [10] "days_to_death"                                             
    ##  [11] "days_to_last_followup"                                     
    ##  [12] "race_list"                                                 
    ##  [13] "tissue_source_site"                                        
    ##  [14] "patient_id"                                                
    ##  [15] "bcr_patient_uuid"                                          
    ##  [16] "history_of_neoadjuvant_treatment"                          
    ##  [17] "informed_consent_verified"                                 
    ##  [18] "icd_o_3_site"                                              
    ##  [19] "icd_o_3_histology"                                         
    ##  [20] "icd_10"                                                    
    ##  [21] "tissue_prospective_collection_indicator"                   
    ##  [22] "tissue_retrospective_collection_indicator"                 
    ##  [23] "days_to_initial_pathologic_diagnosis"                      
    ##  [24] "age_at_initial_pathologic_diagnosis"                       
    ##  [25] "year_of_initial_pathologic_diagnosis"                      
    ##  [26] "person_neoplasm_cancer_status"                             
    ##  [27] "ethnicity"                                                 
    ##  [28] "weight"                                                    
    ##  [29] "height"                                                    
    ##  [30] "day_of_form_completion"                                    
    ##  [31] "month_of_form_completion"                                  
    ##  [32] "year_of_form_completion"                                   
    ##  [33] "residual_tumor"                                            
    ##  [34] "anatomic_neoplasm_subdivision"                             
    ##  [35] "primary_lymph_node_presentation_assessment"                
    ##  [36] "lymph_node_examined_count"                                 
    ##  [37] "number_of_lymphnodes_positive_by_he"                       
    ##  [38] "number_of_lymphnodes_positive_by_ihc"                      
    ##  [39] "preoperative_pretreatment_cea_level"                       
    ##  [40] "non_nodal_tumor_deposits"                                  
    ##  [41] "circumferential_resection_margin"                          
    ##  [42] "venous_invasion"                                           
    ##  [43] "lymphatic_invasion"                                        
    ##  [44] "perineural_invasion_present"                               
    ##  [45] "microsatellite_instability"                                
    ##  [46] "number_of_loci_tested"                                     
    ##  [47] "number_of_abnormal_loci"                                   
    ##  [48] "kras_gene_analysis_performed"                              
    ##  [49] "kras_mutation_found"                                       
    ##  [50] "kras_mutation_codon"                                       
    ##  [51] "braf_gene_analysis_performed"                              
    ##  [52] "braf_gene_analysis_result"                                 
    ##  [53] "synchronous_colon_cancer_present"                          
    ##  [54] "history_of_colon_polyps"                                   
    ##  [55] "colon_polyps_present"                                      
    ##  [56] "loss_expression_of_mismatch_repair_proteins_by_ihc"        
    ##  [57] "loss_expression_of_mismatch_repair_proteins_by_ihc_results"
    ##  [58] "number_of_first_degree_relatives_with_cancer_diagnosis"    
    ##  [59] "radiation_therapy"                                         
    ##  [60] "postoperative_rx_tx"                                       
    ##  [61] "primary_therapy_outcome_success"                           
    ##  [62] "has_new_tumor_events_information"                          
    ##  [63] "has_drugs_information"                                     
    ##  [64] "has_radiations_information"                                
    ##  [65] "has_follow_ups_information"                                
    ##  [66] "project"                                                   
    ##  [67] "stage_event_system_version"                                
    ##  [68] "stage_event_clinical_stage"                                
    ##  [69] "stage_event_pathologic_stage"                              
    ##  [70] "stage_event_tnm_categories"                                
    ##  [71] "stage_event_psa"                                           
    ##  [72] "stage_event_gleason_grading"                               
    ##  [73] "stage_event_ann_arbor"                                     
    ##  [74] "stage_event_serum_markers"                                 
    ##  [75] "stage_event_igcccg_stage"                                  
    ##  [76] "stage_event_masaoka_stage"                                 
    ##  [77] "cancer_type"                                               
    ##  [78] "patient_death_reason"                                      
    ##  [79] "anatomic_neoplasm_subdivision_other"                       
    ##  [80] "neoplasm_histologic_grade"                                 
    ##  [81] "country_of_procurement"                                    
    ##  [82] "city_of_procurement"                                       
    ##  [83] "reflux_history"                                            
    ##  [84] "antireflux_treatment"                                      
    ##  [85] "antireflux_treatment_types"                                
    ##  [86] "barretts_esophagus"                                        
    ##  [87] "h_pylori_infection"                                        
    ##  [88] "family_history_of_stomach_cancer"                          
    ##  [89] "number_of_relatives_with_stomach_cancer"                   
    ##  [90] "relative_family_cancer_history"                            
    ##  [91] "cancer_first_degree_relative"                              
    ##  [92] "blood_relative_cancer_history_list"                        
    ##  [93] "history_hepato_carcinoma_risk_factors"                     
    ##  [94] "post_op_ablation_embolization_tx"                          
    ##  [95] "eastern_cancer_oncology_group"                             
    ##  [96] "primary_pathology_tumor_tissue_site"                       
    ##  [97] "primary_pathology_histological_type"                       
    ##  [98] "primary_pathology_specimen_collection_method_name"         
    ##  [99] "primary_pathology_history_prior_surgery_type_other"        
    ## [100] "primary_pathology_days_to_initial_pathologic_diagnosis"    
    ## [101] "primary_pathology_age_at_initial_pathologic_diagnosis"     
    ## [102] "primary_pathology_year_of_initial_pathologic_diagnosis"    
    ## [103] "primary_pathology_neoplasm_histologic_grade"               
    ## [104] "primary_pathology_residual_tumor"                          
    ## [105] "primary_pathology_vascular_tumor_cell_type"                
    ## [106] "primary_pathology_perineural_invasion_present"             
    ## [107] "primary_pathology_child_pugh_classification_grade"         
    ## [108] "primary_pathology_ca_19_9_level"                           
    ## [109] "primary_pathology_ca_19_9_level_lower"                     
    ## [110] "primary_pathology_ca_19_9_level_upper"                     
    ## [111] "primary_pathology_fetoprotein_outcome_value"               
    ## [112] "primary_pathology_fetoprotein_outcome_lower_limit"         
    ## [113] "primary_pathology_fetoprotein_outcome_upper_limit"         
    ## [114] "primary_pathology_platelet_result_count"                   
    ## [115] "primary_pathology_platelet_result_lower_limit"             
    ## [116] "primary_pathology_platelet_result_upper_limit"             
    ## [117] "primary_pathology_prothrombin_time_result_value"           
    ## [118] "primary_pathology_inter_norm_ratio_lower_limit"            
    ## [119] "primary_pathology_intern_norm_ratio_upper_limit"           
    ## [120] "primary_pathology_albumin_result_specified_value"          
    ## [121] "primary_pathology_albumin_result_lower_limit"              
    ## [122] "primary_pathology_albumin_result_upper_limit"              
    ## [123] "primary_pathology_bilirubin_upper_limit"                   
    ## [124] "primary_pathology_bilirubin_lower_limit"                   
    ## [125] "primary_pathology_total_bilirubin_upper_limit"             
    ## [126] "primary_pathology_creatinine_value_in_mg_dl"               
    ## [127] "primary_pathology_creatinine_lower_level"                  
    ## [128] "primary_pathology_creatinine_upper_limit"                  
    ## [129] "primary_pathology_fibrosis_ishak_score"                    
    ## [130] "primary_pathology_cholangitis_tissue_evidence"             
    ## [131] "adenocarcinoma_invasion"                                   
    ## [132] "histological_type_other"                                   
    ## [133] "tumor_type"                                                
    ## [134] "initial_pathologic_diagnosis_method"                       
    ## [135] "init_pathology_dx_method_other"                            
    ## [136] "surgery_performed_type"                                    
    ## [137] "histologic_grading_tier_category"                          
    ## [138] "maximum_tumor_dimension"                                   
    ## [139] "source_of_patient_death_reason"                            
    ## [140] "tobacco_smoking_history"                                   
    ## [141] "year_of_tobacco_smoking_onset"                             
    ## [142] "stopped_smoking_year"                                      
    ## [143] "number_pack_years_smoked"                                  
    ## [144] "alcohol_history_documented"                                
    ## [145] "alcoholic_exposure_category"                               
    ## [146] "frequency_of_alcohol_consumption"                          
    ## [147] "amount_of_alcohol_consumption_per_day"                     
    ## [148] "history_of_diabetes"                                       
    ## [149] "days_to_diabetes_onset"                                    
    ## [150] "history_of_chronic_pancreatitis"                           
    ## [151] "days_to_pancreatitis_onset"                                
    ## [152] "family_history_of_cancer"                                  
    ## [153] "relative_cancer_types"                                     
    ## [154] "history_prior_surgery_type_other"

``` r
# Visualizar estructura del conjunto de datos y un resumen estadístico
str(df_TCGA)
```

    ## 'data.frame':    1308 obs. of  154 variables:
    ##  $ bcr_patient_barcode                                       : chr  "TCGA-3L-AA1B" "TCGA-4N-A93T" "TCGA-4T-AA8H" "TCGA-5M-AAT4" ...
    ##  $ additional_studies                                        : chr  NA NA NA NA ...
    ##  $ tumor_tissue_site                                         : chr  "Colon" "Colon" "Colon" "Colon" ...
    ##  $ histological_type                                         : chr  "Colon Adenocarcinoma" "Colon Adenocarcinoma" "Colon Mucinous Adenocarcinoma" "Colon Adenocarcinoma" ...
    ##  $ other_dx                                                  : chr  "No" "No" "No" "No" ...
    ##  $ gender                                                    : chr  "FEMALE" "MALE" "FEMALE" "MALE" ...
    ##  $ vital_status                                              : chr  NA NA NA "Dead" ...
    ##  $ days_to_birth                                             : int  -22379 -24523 -15494 -27095 -14852 -27870 -16512 -31329 -30237 -26292 ...
    ##  $ days_to_last_known_alive                                  : int  NA NA NA NA NA NA 424 NA NA NA ...
    ##  $ days_to_death                                             : int  NA NA NA 49 290 NA NA 1126 NA NA ...
    ##  $ days_to_last_followup                                     : num  475 146 385 -Inf -Inf ...
    ##  $ race_list                                                 : chr  "BLACK OR AFRICAN AMERICAN" "BLACK OR AFRICAN AMERICAN" "BLACK OR AFRICAN AMERICAN" "BLACK OR AFRICAN AMERICAN" ...
    ##  $ tissue_source_site                                        : chr  "3L" "4N" "4T" "5M" ...
    ##  $ patient_id                                                : chr  "AA1B" "A93T" "AA8H" "AAT4" ...
    ##  $ bcr_patient_uuid                                          : chr  "A94E1279-A975-480A-93E9-7B1FF05CBCBF" "92554413-9EBC-4354-8E1B-9682F3A031D9" "A5E14ADD-1552-4606-9FFE-3A03BCF76640" "1136DD50-242A-4659-AAD4-C53F9E759BB3" ...
    ##  $ history_of_neoadjuvant_treatment                          : chr  "No" "No" "No" "No" ...
    ##  $ informed_consent_verified                                 : chr  "YES" "YES" "YES" "YES" ...
    ##  $ icd_o_3_site                                              : chr  "C18.0" "C18.2" "C18.6" "C18.2" ...
    ##  $ icd_o_3_histology                                         : chr  "01/03/8140" "01/03/8140" "01/03/8480" "01/03/8140" ...
    ##  $ icd_10                                                    : chr  "C18.0" "C18.2" "C18.6" "C18.2" ...
    ##  $ tissue_prospective_collection_indicator                   : chr  "YES" "YES" "NO" "NO" ...
    ##  $ tissue_retrospective_collection_indicator                 : chr  "NO" "NO" "YES" "YES" ...
    ##  $ days_to_initial_pathologic_diagnosis                      : int  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ age_at_initial_pathologic_diagnosis                       : int  61 67 42 74 40 76 45 85 82 71 ...
    ##  $ year_of_initial_pathologic_diagnosis                      : int  2013 2013 2013 2009 2009 2011 2009 2009 2009 2009 ...
    ##  $ person_neoplasm_cancer_status                             : chr  "TUMOR FREE" "WITH TUMOR" "TUMOR FREE" "WITH TUMOR" ...
    ##  $ ethnicity                                                 : chr  "NOT HISPANIC OR LATINO" "NOT HISPANIC OR LATINO" "NOT HISPANIC OR LATINO" "HISPANIC OR LATINO" ...
    ##  $ weight                                                    : num  63.3 134 108 NA 99.1 ...
    ##  $ height                                                    : num  173 168 168 NA 162 ...
    ##  $ day_of_form_completion                                    : int  22 1 5 27 27 27 14 28 4 15 ...
    ##  $ month_of_form_completion                                  : int  4 10 6 1 1 1 10 1 10 10 ...
    ##  $ year_of_form_completion                                   : int  2014 2014 2014 2015 2015 2015 2010 2011 2010 2010 ...
    ##  $ residual_tumor                                            : chr  "R0" "R0" "R0" "R0" ...
    ##  $ anatomic_neoplasm_subdivision                             : chr  "Cecum" "Ascending Colon" "Descending Colon" "Ascending Colon" ...
    ##  $ primary_lymph_node_presentation_assessment                : chr  "YES" "YES" "YES" "YES" ...
    ##  $ lymph_node_examined_count                                 : int  28 25 24 3 11 15 22 27 29 20 ...
    ##  $ number_of_lymphnodes_positive_by_he                       : int  0 NA 0 0 10 0 NA 3 1 7 ...
    ##  $ number_of_lymphnodes_positive_by_ihc                      : int  0 2 NA 0 0 0 NA NA NA 0 ...
    ##  $ preoperative_pretreatment_cea_level                       : num  NA 2 NA 550 2.61 2.91 NA 17.4 3.4 32.8 ...
    ##  $ non_nodal_tumor_deposits                                  : chr  "NO" "YES" "NO" "NO" ...
    ##  $ circumferential_resection_margin                          : num  NA 30 20 NA NA NA NA NA NA NA ...
    ##  $ venous_invasion                                           : chr  "NO" "NO" "NO" "YES" ...
    ##  $ lymphatic_invasion                                        : chr  "NO" "NO" "NO" NA ...
    ##  $ perineural_invasion_present                               : chr  "NO" "NO" "NO" NA ...
    ##  $ microsatellite_instability                                : chr  "NO" NA "NO" NA ...
    ##  $ number_of_loci_tested                                     : int  NA NA NA NA NA NA NA NA NA NA ...
    ##  $ number_of_abnormal_loci                                   : int  NA NA NA NA NA NA NA NA NA NA ...
    ##  $ kras_gene_analysis_performed                              : chr  "NO" "NO" "NO" "NO" ...
    ##  $ kras_mutation_found                                       : chr  NA NA NA NA ...
    ##  $ kras_mutation_codon                                       : int  NA NA NA NA NA NA NA NA NA NA ...
    ##  $ braf_gene_analysis_performed                              : chr  "NO" "NO" "NO" "NO" ...
    ##  $ braf_gene_analysis_result                                 : chr  NA NA NA NA ...
    ##  $ synchronous_colon_cancer_present                          : chr  "NO" "YES" "NO" "NO" ...
    ##  $ history_of_colon_polyps                                   : chr  "YES" "NO" "NO" "NO" ...
    ##  $ colon_polyps_present                                      : chr  "YES" "YES" "NO" "YES" ...
    ##  $ loss_expression_of_mismatch_repair_proteins_by_ihc        : chr  "YES" "YES" "YES" "NO" ...
    ##  $ loss_expression_of_mismatch_repair_proteins_by_ihc_results: chr  "MLH1-ExpressedMSH2-ExpressedPMS2-ExpressedMSH6-Expressed" "MLH1-ExpressedMSH2-ExpressedPMS2-ExpressedMSH6-Expressed" "MLH1-ExpressedMSH2-ExpressedPMS2-ExpressedMSH6-Expressed" NA ...
    ##  $ number_of_first_degree_relatives_with_cancer_diagnosis    : int  0 0 0 0 NA 0 0 0 0 0 ...
    ##  $ radiation_therapy                                         : chr  "NO" "NO" "NO" "NO" ...
    ##  $ postoperative_rx_tx                                       : chr  "NO" "YES" "NO" "NO" ...
    ##  $ primary_therapy_outcome_success                           : chr  "Complete Remission/Response" "Stable Disease" "Complete Remission/Response" "Progressive Disease" ...
    ##  $ has_new_tumor_events_information                          : chr  "NO" "NO" "NO" "NO" ...
    ##  $ has_drugs_information                                     : chr  "NO" "YES" "NO" "NO" ...
    ##  $ has_radiations_information                                : chr  "NO" "NO" "NO" "NO" ...
    ##  $ has_follow_ups_information                                : chr  "YES" "YES" "YES" "YES" ...
    ##  $ project                                                   : chr  "TCGA-COAD" "TCGA-COAD" "TCGA-COAD" "TCGA-COAD" ...
    ##  $ stage_event_system_version                                : chr  "7th" "7th" "7th" "6th" ...
    ##  $ stage_event_clinical_stage                                : logi  NA NA NA NA NA NA ...
    ##  $ stage_event_pathologic_stage                              : chr  "Stage I" "Stage IIIB" "Stage IIA" "Stage IV" ...
    ##  $ stage_event_tnm_categories                                : chr  "T2N0M0" "T4aN1bM0" "T3N0MX" "T3N0M1b" ...
    ##  $ stage_event_psa                                           : logi  NA NA NA NA NA NA ...
    ##  $ stage_event_gleason_grading                               : logi  NA NA NA NA NA NA ...
    ##  $ stage_event_ann_arbor                                     : logi  NA NA NA NA NA NA ...
    ##  $ stage_event_serum_markers                                 : logi  NA NA NA NA NA NA ...
    ##  $ stage_event_igcccg_stage                                  : logi  NA NA NA NA NA NA ...
    ##  $ stage_event_masaoka_stage                                 : logi  NA NA NA NA NA NA ...
    ##  $ cancer_type                                               : chr  "TCGA-COAD" "TCGA-COAD" "TCGA-COAD" "TCGA-COAD" ...
    ##  $ patient_death_reason                                      : chr  NA NA NA NA ...
    ##  $ anatomic_neoplasm_subdivision_other                       : chr  NA NA NA NA ...
    ##  $ neoplasm_histologic_grade                                 : chr  NA NA NA NA ...
    ##  $ country_of_procurement                                    : chr  NA NA NA NA ...
    ##  $ city_of_procurement                                       : chr  NA NA NA NA ...
    ##  $ reflux_history                                            : chr  NA NA NA NA ...
    ##  $ antireflux_treatment                                      : chr  NA NA NA NA ...
    ##  $ antireflux_treatment_types                                : chr  NA NA NA NA ...
    ##  $ barretts_esophagus                                        : chr  NA NA NA NA ...
    ##  $ h_pylori_infection                                        : chr  NA NA NA NA ...
    ##  $ family_history_of_stomach_cancer                          : chr  NA NA NA NA ...
    ##  $ number_of_relatives_with_stomach_cancer                   : int  NA NA NA NA NA NA NA NA NA NA ...
    ##  $ relative_family_cancer_history                            : chr  NA NA NA NA ...
    ##  $ cancer_first_degree_relative                              : int  NA NA NA NA NA NA NA NA NA NA ...
    ##  $ blood_relative_cancer_history_list                        : chr  NA NA NA NA ...
    ##  $ history_hepato_carcinoma_risk_factors                     : chr  NA NA NA NA ...
    ##  $ post_op_ablation_embolization_tx                          : chr  NA NA NA NA ...
    ##  $ eastern_cancer_oncology_group                             : int  NA NA NA NA NA NA NA NA NA NA ...
    ##  $ primary_pathology_tumor_tissue_site                       : chr  NA NA NA NA ...
    ##  $ primary_pathology_histological_type                       : chr  NA NA NA NA ...
    ##  $ primary_pathology_specimen_collection_method_name         : chr  NA NA NA NA ...
    ##  $ primary_pathology_history_prior_surgery_type_other        : chr  NA NA NA NA ...
    ##   [list output truncated]

``` r
summary(df_TCGA)
```

    ##  bcr_patient_barcode additional_studies tumor_tissue_site  histological_type 
    ##  Length:1308         Length:1308        Length:1308        Length:1308       
    ##  Class :character    Class :character   Class :character   Class :character  
    ##  Mode  :character    Mode  :character   Mode  :character   Mode  :character  
    ##                                                                              
    ##                                                                              
    ##                                                                              
    ##                                                                              
    ##    other_dx            gender          vital_status       days_to_birth   
    ##  Length:1308        Length:1308        Length:1308        Min.   :-32873  
    ##  Class :character   Class :character   Class :character   1st Qu.:-27321  
    ##  Mode  :character   Mode  :character   Mode  :character   Median :-24655  
    ##                                                           Mean   :-24204  
    ##                                                           3rd Qu.:-21319  
    ##                                                           Max.   :-10659  
    ##                                                           NA's   :17      
    ##  days_to_last_known_alive days_to_death    days_to_last_followup
    ##  Min.   :   0.0           Min.   :   0.0   Min.   :-Inf         
    ##  1st Qu.: 343.5           1st Qu.: 145.0   1st Qu.:   0         
    ##  Median : 992.0           Median : 366.0   Median : 485         
    ##  Mean   :1635.1           Mean   : 503.5   Mean   :-Inf         
    ##  3rd Qu.:3223.0           3rd Qu.: 641.0   3rd Qu.: 942         
    ##  Max.   :3920.0           Max.   :3042.0   Max.   :4502         
    ##  NA's   :1289             NA's   :1071     NA's   :205          
    ##   race_list         tissue_source_site  patient_id        bcr_patient_uuid  
    ##  Length:1308        Length:1308        Length:1308        Length:1308       
    ##  Class :character   Class :character   Class :character   Class :character  
    ##  Mode  :character   Mode  :character   Mode  :character   Mode  :character  
    ##                                                                             
    ##                                                                             
    ##                                                                             
    ##                                                                             
    ##  history_of_neoadjuvant_treatment informed_consent_verified icd_o_3_site      
    ##  Length:1308                      Length:1308               Length:1308       
    ##  Class :character                 Class :character          Class :character  
    ##  Mode  :character                 Mode  :character          Mode  :character  
    ##                                                                               
    ##                                                                               
    ##                                                                               
    ##                                                                               
    ##  icd_o_3_histology     icd_10          tissue_prospective_collection_indicator
    ##  Length:1308        Length:1308        Length:1308                            
    ##  Class :character   Class :character   Class :character                       
    ##  Mode  :character   Mode  :character   Mode  :character                       
    ##                                                                               
    ##                                                                               
    ##                                                                               
    ##                                                                               
    ##  tissue_retrospective_collection_indicator days_to_initial_pathologic_diagnosis
    ##  Length:1308                               Min.   :0                           
    ##  Class :character                          1st Qu.:0                           
    ##  Mode  :character                          Median :0                           
    ##                                            Mean   :0                           
    ##                                            3rd Qu.:0                           
    ##                                            Max.   :0                           
    ##                                            NA's   :58                          
    ##  age_at_initial_pathologic_diagnosis year_of_initial_pathologic_diagnosis
    ##  Min.   :30.00                       Min.   :1996                        
    ##  1st Qu.:58.00                       1st Qu.:2008                        
    ##  Median :67.00                       Median :2010                        
    ##  Mean   :65.85                       Mean   :2009                        
    ##  3rd Qu.:74.00                       3rd Qu.:2011                        
    ##  Max.   :90.00                       Max.   :2013                        
    ##  NA's   :56                          NA's   :53                          
    ##  person_neoplasm_cancer_status  ethnicity             weight      
    ##  Length:1308                   Length:1308        Min.   : 34.00  
    ##  Class :character              Class :character   1st Qu.: 65.00  
    ##  Mode  :character              Mode  :character   Median : 78.05  
    ##                                                   Mean   : 80.44  
    ##                                                   3rd Qu.: 91.85  
    ##                                                   Max.   :175.30  
    ##                                                   NA's   :936     
    ##      height      day_of_form_completion month_of_form_completion
    ##  Min.   : 80.3   Min.   : 1.00          Min.   : 1.000          
    ##  1st Qu.:162.0   1st Qu.:11.00          1st Qu.: 4.000          
    ##  Median :170.0   Median :18.00          Median : 6.000          
    ##  Mean   :168.9   Mean   :16.62          Mean   : 5.881          
    ##  3rd Qu.:176.0   3rd Qu.:22.00          3rd Qu.: 7.000          
    ##  Max.   :193.0   Max.   :31.00          Max.   :12.000          
    ##  NA's   :957     NA's   :41             NA's   :5               
    ##  year_of_form_completion residual_tumor     anatomic_neoplasm_subdivision
    ##  Min.   :2010            Length:1308        Length:1308                  
    ##  1st Qu.:2011            Class :character   Class :character             
    ##  Median :2011            Mode  :character   Mode  :character             
    ##  Mean   :2012                                                            
    ##  3rd Qu.:2013                                                            
    ##  Max.   :2015                                                            
    ##  NA's   :5                                                               
    ##  primary_lymph_node_presentation_assessment lymph_node_examined_count
    ##  Length:1308                                Min.   :  0.00           
    ##  Class :character                           1st Qu.: 12.00           
    ##  Mode  :character                           Median : 18.00           
    ##                                             Mean   : 21.42           
    ##                                             3rd Qu.: 27.00           
    ##                                             Max.   :109.00           
    ##                                             NA's   :138              
    ##  number_of_lymphnodes_positive_by_he number_of_lymphnodes_positive_by_ihc
    ##  Min.   : 0.000                      Min.   : 0.0000                     
    ##  1st Qu.: 0.000                      1st Qu.: 0.0000                     
    ##  Median : 1.000                      Median : 0.0000                     
    ##  Mean   : 3.493                      Mean   : 0.2833                     
    ##  3rd Qu.: 4.000                      3rd Qu.: 0.0000                     
    ##  Max.   :57.000                      Max.   :12.0000                     
    ##  NA's   :143                         NA's   :1188                        
    ##  preoperative_pretreatment_cea_level non_nodal_tumor_deposits
    ##  Min.   :   0.000                    Length:1308             
    ##  1st Qu.:   1.700                    Class :character        
    ##  Median :   3.160                    Mode  :character        
    ##  Mean   :  65.074                                            
    ##  3rd Qu.:   8.982                                            
    ##  Max.   :7868.000                                            
    ##  NA's   :906                                                 
    ##  circumferential_resection_margin venous_invasion    lymphatic_invasion
    ##  Min.   :  0.00                   Length:1308        Length:1308       
    ##  1st Qu.:  2.50                   Class :character   Class :character  
    ##  Median : 13.00                   Mode  :character   Mode  :character  
    ##  Mean   : 22.96                                                        
    ##  3rd Qu.: 30.00                                                        
    ##  Max.   :165.00                                                        
    ##  NA's   :1185                                                          
    ##  perineural_invasion_present microsatellite_instability number_of_loci_tested
    ##  Length:1308                 Length:1308                Min.   : 0.000       
    ##  Class :character            Class :character           1st Qu.: 5.000       
    ##  Mode  :character            Mode  :character           Median : 5.000       
    ##                                                         Mean   : 4.944       
    ##                                                         3rd Qu.: 5.000       
    ##                                                         Max.   :10.000       
    ##                                                         NA's   :1236         
    ##  number_of_abnormal_loci kras_gene_analysis_performed kras_mutation_found
    ##  Min.   :0.0000          Length:1308                  Length:1308        
    ##  1st Qu.:0.0000          Class :character             Class :character   
    ##  Median :0.0000          Mode  :character             Mode  :character   
    ##  Mean   :0.5915                                                          
    ##  3rd Qu.:0.0000                                                          
    ##  Max.   :9.0000                                                          
    ##  NA's   :1237                                                            
    ##  kras_mutation_codon braf_gene_analysis_performed braf_gene_analysis_result
    ##  Min.   :12.00       Length:1308                  Length:1308              
    ##  1st Qu.:12.00       Class :character             Class :character         
    ##  Median :12.00       Mode  :character             Mode  :character         
    ##  Mean   :13.83                                                             
    ##  3rd Qu.:12.00                                                             
    ##  Max.   :61.00                                                             
    ##  NA's   :1278                                                              
    ##  synchronous_colon_cancer_present history_of_colon_polyps colon_polyps_present
    ##  Length:1308                      Length:1308             Length:1308         
    ##  Class :character                 Class :character        Class :character    
    ##  Mode  :character                 Mode  :character        Mode  :character    
    ##                                                                               
    ##                                                                               
    ##                                                                               
    ##                                                                               
    ##  loss_expression_of_mismatch_repair_proteins_by_ihc
    ##  Length:1308                                       
    ##  Class :character                                  
    ##  Mode  :character                                  
    ##                                                    
    ##                                                    
    ##                                                    
    ##                                                    
    ##  loss_expression_of_mismatch_repair_proteins_by_ihc_results
    ##  Length:1308                                               
    ##  Class :character                                          
    ##  Mode  :character                                          
    ##                                                            
    ##                                                            
    ##                                                            
    ##                                                            
    ##  number_of_first_degree_relatives_with_cancer_diagnosis radiation_therapy 
    ##  Min.   :0.0000                                         Length:1308       
    ##  1st Qu.:0.0000                                         Class :character  
    ##  Median :0.0000                                         Mode  :character  
    ##  Mean   :0.1654                                                           
    ##  3rd Qu.:0.0000                                                           
    ##  Max.   :3.0000                                                           
    ##  NA's   :770                                                              
    ##  postoperative_rx_tx primary_therapy_outcome_success
    ##  Length:1308         Length:1308                    
    ##  Class :character    Class :character               
    ##  Mode  :character    Mode  :character               
    ##                                                     
    ##                                                     
    ##                                                     
    ##                                                     
    ##  has_new_tumor_events_information has_drugs_information
    ##  Length:1308                      Length:1308          
    ##  Class :character                 Class :character     
    ##  Mode  :character                 Mode  :character     
    ##                                                        
    ##                                                        
    ##                                                        
    ##                                                        
    ##  has_radiations_information has_follow_ups_information   project         
    ##  Length:1308                Length:1308                Length:1308       
    ##  Class :character           Class :character           Class :character  
    ##  Mode  :character           Mode  :character           Mode  :character  
    ##                                                                          
    ##                                                                          
    ##                                                                          
    ##                                                                          
    ##  stage_event_system_version stage_event_clinical_stage
    ##  Length:1308                Mode:logical              
    ##  Class :character           NA's:1308                 
    ##  Mode  :character                                     
    ##                                                       
    ##                                                       
    ##                                                       
    ##                                                       
    ##  stage_event_pathologic_stage stage_event_tnm_categories stage_event_psa
    ##  Length:1308                  Length:1308                Mode:logical   
    ##  Class :character             Class :character           NA's:1308      
    ##  Mode  :character             Mode  :character                          
    ##                                                                         
    ##                                                                         
    ##                                                                         
    ##                                                                         
    ##  stage_event_gleason_grading stage_event_ann_arbor stage_event_serum_markers
    ##  Mode:logical                Mode:logical          Mode:logical             
    ##  NA's:1308                   NA's:1308             NA's:1308                
    ##                                                                             
    ##                                                                             
    ##                                                                             
    ##                                                                             
    ##                                                                             
    ##  stage_event_igcccg_stage stage_event_masaoka_stage cancer_type       
    ##  Mode:logical             Mode:logical              Length:1308       
    ##  NA's:1308                NA's:1308                 Class :character  
    ##                                                     Mode  :character  
    ##                                                                       
    ##                                                                       
    ##                                                                       
    ##                                                                       
    ##  patient_death_reason anatomic_neoplasm_subdivision_other
    ##  Length:1308          Length:1308                        
    ##  Class :character     Class :character                   
    ##  Mode  :character     Mode  :character                   
    ##                                                          
    ##                                                          
    ##                                                          
    ##                                                          
    ##  neoplasm_histologic_grade country_of_procurement city_of_procurement
    ##  Length:1308               Length:1308            Length:1308        
    ##  Class :character          Class :character       Class :character   
    ##  Mode  :character          Mode  :character       Mode  :character   
    ##                                                                      
    ##                                                                      
    ##                                                                      
    ##                                                                      
    ##  reflux_history     antireflux_treatment antireflux_treatment_types
    ##  Length:1308        Length:1308          Length:1308               
    ##  Class :character   Class :character     Class :character          
    ##  Mode  :character   Mode  :character     Mode  :character          
    ##                                                                    
    ##                                                                    
    ##                                                                    
    ##                                                                    
    ##  barretts_esophagus h_pylori_infection family_history_of_stomach_cancer
    ##  Length:1308        Length:1308        Length:1308                     
    ##  Class :character   Class :character   Class :character                
    ##  Mode  :character   Mode  :character   Mode  :character                
    ##                                                                        
    ##                                                                        
    ##                                                                        
    ##                                                                        
    ##  number_of_relatives_with_stomach_cancer relative_family_cancer_history
    ##  Min.   :0.0000                          Length:1308                   
    ##  1st Qu.:0.0000                          Class :character              
    ##  Median :0.0000                          Mode  :character              
    ##  Mean   :0.2632                                                        
    ##  3rd Qu.:0.0000                                                        
    ##  Max.   :2.0000                                                        
    ##  NA's   :1232                                                          
    ##  cancer_first_degree_relative blood_relative_cancer_history_list
    ##  Min.   :1.000                Length:1308                       
    ##  1st Qu.:1.000                Class :character                  
    ##  Median :1.000                Mode  :character                  
    ##  Mean   :1.696                                                  
    ##  3rd Qu.:2.000                                                  
    ##  Max.   :4.000                                                  
    ##  NA's   :1285                                                   
    ##  history_hepato_carcinoma_risk_factors post_op_ablation_embolization_tx
    ##  Length:1308                           Length:1308                     
    ##  Class :character                      Class :character                
    ##  Mode  :character                      Mode  :character                
    ##                                                                        
    ##                                                                        
    ##                                                                        
    ##                                                                        
    ##  eastern_cancer_oncology_group primary_pathology_tumor_tissue_site
    ##  Min.   :0.0000                Length:1308                        
    ##  1st Qu.:0.0000                Class :character                   
    ##  Median :0.0000                Mode  :character                   
    ##  Mean   :0.3514                                                   
    ##  3rd Qu.:1.0000                                                   
    ##  Max.   :3.0000                                                   
    ##  NA's   :1271                                                     
    ##  primary_pathology_histological_type
    ##  Length:1308                        
    ##  Class :character                   
    ##  Mode  :character                   
    ##                                     
    ##                                     
    ##                                     
    ##                                     
    ##  primary_pathology_specimen_collection_method_name
    ##  Length:1308                                      
    ##  Class :character                                 
    ##  Mode  :character                                 
    ##                                                   
    ##                                                   
    ##                                                   
    ##                                                   
    ##  primary_pathology_history_prior_surgery_type_other
    ##  Length:1308                                       
    ##  Class :character                                  
    ##  Mode  :character                                  
    ##                                                    
    ##                                                    
    ##                                                    
    ##                                                    
    ##  primary_pathology_days_to_initial_pathologic_diagnosis
    ##  Min.   :0                                             
    ##  1st Qu.:0                                             
    ##  Median :0                                             
    ##  Mean   :0                                             
    ##  3rd Qu.:0                                             
    ##  Max.   :0                                             
    ##  NA's   :1260                                          
    ##  primary_pathology_age_at_initial_pathologic_diagnosis
    ##  Min.   :29.00                                        
    ##  1st Qu.:58.00                                        
    ##  Median :66.00                                        
    ##  Mean   :63.64                                        
    ##  3rd Qu.:73.00                                        
    ##  Max.   :82.00                                        
    ##  NA's   :1263                                         
    ##  primary_pathology_year_of_initial_pathologic_diagnosis
    ##  Min.   :2005                                          
    ##  1st Qu.:2010                                          
    ##  Median :2011                                          
    ##  Mean   :2011                                          
    ##  3rd Qu.:2012                                          
    ##  Max.   :2013                                          
    ##  NA's   :1260                                          
    ##  primary_pathology_neoplasm_histologic_grade primary_pathology_residual_tumor
    ##  Length:1308                                 Length:1308                     
    ##  Class :character                            Class :character                
    ##  Mode  :character                            Mode  :character                
    ##                                                                              
    ##                                                                              
    ##                                                                              
    ##                                                                              
    ##  primary_pathology_vascular_tumor_cell_type
    ##  Length:1308                               
    ##  Class :character                          
    ##  Mode  :character                          
    ##                                            
    ##                                            
    ##                                            
    ##                                            
    ##  primary_pathology_perineural_invasion_present
    ##  Length:1308                                  
    ##  Class :character                             
    ##  Mode  :character                             
    ##                                               
    ##                                               
    ##                                               
    ##                                               
    ##  primary_pathology_child_pugh_classification_grade
    ##  Length:1308                                      
    ##  Class :character                                 
    ##  Mode  :character                                 
    ##                                                   
    ##                                                   
    ##                                                   
    ##                                                   
    ##  primary_pathology_ca_19_9_level primary_pathology_ca_19_9_level_lower
    ##  Min.   :   1.0                  Min.   :0.0000                       
    ##  1st Qu.:  26.5                  1st Qu.:0.0000                       
    ##  Median :  54.9                  Median :0.0000                       
    ##  Mean   : 345.8                  Mean   :0.1351                       
    ##  3rd Qu.: 200.0                  3rd Qu.:0.0000                       
    ##  Max.   :6910.0                  Max.   :5.0000                       
    ##  NA's   :1268                    NA's   :1271                         
    ##  primary_pathology_ca_19_9_level_upper
    ##  Min.   : 7.00                        
    ##  1st Qu.:35.00                        
    ##  Median :55.00                        
    ##  Mean   :45.56                        
    ##  3rd Qu.:55.00                        
    ##  Max.   :55.00                        
    ##  NA's   :1267                         
    ##  primary_pathology_fetoprotein_outcome_value
    ##  Min.   : 1.380                             
    ##  1st Qu.: 2.350                             
    ##  Median : 3.100                             
    ##  Mean   : 3.803                             
    ##  3rd Qu.: 4.325                             
    ##  Max.   :14.000                             
    ##  NA's   :1280                               
    ##  primary_pathology_fetoprotein_outcome_lower_limit
    ##  Min.   :0                                        
    ##  1st Qu.:0                                        
    ##  Median :0                                        
    ##  Mean   :0                                        
    ##  3rd Qu.:0                                        
    ##  Max.   :0                                        
    ##  NA's   :1275                                     
    ##  primary_pathology_fetoprotein_outcome_upper_limit
    ##  Min.   : 6.000                                   
    ##  1st Qu.: 6.000                                   
    ##  Median : 6.000                                   
    ##  Mean   : 7.579                                   
    ##  3rd Qu.: 9.000                                   
    ##  Max.   :15.000                                   
    ##  NA's   :1275                                     
    ##  primary_pathology_platelet_result_count
    ##  Min.   :   134.0                       
    ##  1st Qu.:   211.0                       
    ##  Median :   269.5                       
    ##  Mean   : 68013.6                       
    ##  3rd Qu.:179500.0                       
    ##  Max.   :354000.0                       
    ##  NA's   :1264                           
    ##  primary_pathology_platelet_result_lower_limit
    ##  Min.   :   140                               
    ##  1st Qu.:   150                               
    ##  Median :   150                               
    ##  Mean   : 32560                               
    ##  3rd Qu.: 15000                               
    ##  Max.   :150000                               
    ##  NA's   :1263                                 
    ##  primary_pathology_platelet_result_upper_limit
    ##  Min.   :   400                               
    ##  1st Qu.:   450                               
    ##  Median :   450                               
    ##  Mean   :126091                               
    ##  3rd Qu.:400000                               
    ##  Max.   :450000                               
    ##  NA's   :1263                                 
    ##  primary_pathology_prothrombin_time_result_value
    ##  Min.   : 0.800                                 
    ##  1st Qu.: 1.000                                 
    ##  Median : 1.100                                 
    ##  Mean   : 2.012                                 
    ##  3rd Qu.: 1.100                                 
    ##  Max.   :12.200                                 
    ##  NA's   :1266                                   
    ##  primary_pathology_inter_norm_ratio_lower_limit
    ##  Min.   : 0.000                                
    ##  1st Qu.: 0.800                                
    ##  Median : 0.900                                
    ##  Mean   : 1.865                                
    ##  3rd Qu.: 0.900                                
    ##  Max.   :10.400                                
    ##  NA's   :1274                                  
    ##  primary_pathology_intern_norm_ratio_upper_limit
    ##  Min.   : 1.100                                 
    ##  1st Qu.: 1.200                                 
    ##  Median : 1.200                                 
    ##  Mean   : 2.553                                 
    ##  3rd Qu.: 1.200                                 
    ##  Max.   :13.100                                 
    ##  NA's   :1274                                   
    ##  primary_pathology_albumin_result_specified_value
    ##  Min.   :2.400                                   
    ##  1st Qu.:3.775                                   
    ##  Median :4.150                                   
    ##  Mean   :3.998                                   
    ##  3rd Qu.:4.325                                   
    ##  Max.   :4.800                                   
    ##  NA's   :1268                                    
    ##  primary_pathology_albumin_result_lower_limit
    ##  Min.   :3.300                               
    ##  1st Qu.:3.500                               
    ##  Median :3.500                               
    ##  Mean   :3.467                               
    ##  3rd Qu.:3.500                               
    ##  Max.   :3.500                               
    ##  NA's   :1268                                
    ##  primary_pathology_albumin_result_upper_limit
    ##  Min.   :4.50                                
    ##  1st Qu.:4.80                                
    ##  Median :5.00                                
    ##  Mean   :4.91                                
    ##  3rd Qu.:5.00                                
    ##  Max.   :5.20                                
    ##  NA's   :1268                                
    ##  primary_pathology_bilirubin_upper_limit
    ##  Min.   : 0.200                         
    ##  1st Qu.: 0.400                         
    ##  Median : 0.700                         
    ##  Mean   : 2.859                         
    ##  3rd Qu.: 1.025                         
    ##  Max.   :84.000                         
    ##  NA's   :1264                           
    ##  primary_pathology_bilirubin_lower_limit
    ##  Min.   : 0.000                         
    ##  1st Qu.: 0.100                         
    ##  Median : 0.100                         
    ##  Mean   : 2.221                         
    ##  3rd Qu.: 0.150                         
    ##  Max.   :78.000                         
    ##  NA's   :1265                           
    ##  primary_pathology_total_bilirubin_upper_limit
    ##  Min.   : 0.300                               
    ##  1st Qu.: 1.000                               
    ##  Median : 1.000                               
    ##  Mean   : 5.505                               
    ##  3rd Qu.: 1.200                               
    ##  Max.   :96.000                               
    ##  NA's   :1265                                 
    ##  primary_pathology_creatinine_value_in_mg_dl
    ##  Min.   :0.5000                             
    ##  1st Qu.:0.8000                             
    ##  Median :0.9000                             
    ##  Mean   :0.8651                             
    ##  3rd Qu.:1.0000                             
    ##  Max.   :1.4000                             
    ##  NA's   :1265                               
    ##  primary_pathology_creatinine_lower_level
    ##  Min.   :0.4000                          
    ##  1st Qu.:0.6000                          
    ##  Median :0.6000                          
    ##  Mean   :0.6349                          
    ##  3rd Qu.:0.7000                          
    ##  Max.   :0.9000                          
    ##  NA's   :1265                            
    ##  primary_pathology_creatinine_upper_limit
    ##  Min.   :1.000                           
    ##  1st Qu.:1.100                           
    ##  Median :1.100                           
    ##  Mean   :1.198                           
    ##  3rd Qu.:1.300                           
    ##  Max.   :1.400                           
    ##  NA's   :1265                            
    ##  primary_pathology_fibrosis_ishak_score
    ##  Length:1308                           
    ##  Class :character                      
    ##  Mode  :character                      
    ##                                        
    ##                                        
    ##                                        
    ##                                        
    ##  primary_pathology_cholangitis_tissue_evidence adenocarcinoma_invasion
    ##  Length:1308                                   Length:1308            
    ##  Class :character                              Class :character       
    ##  Mode  :character                              Mode  :character       
    ##                                                                       
    ##                                                                       
    ##                                                                       
    ##                                                                       
    ##  histological_type_other  tumor_type        initial_pathologic_diagnosis_method
    ##  Length:1308             Length:1308        Length:1308                        
    ##  Class :character        Class :character   Class :character                   
    ##  Mode  :character        Mode  :character   Mode  :character                   
    ##                                                                                
    ##                                                                                
    ##                                                                                
    ##                                                                                
    ##  init_pathology_dx_method_other surgery_performed_type
    ##  Length:1308                    Length:1308           
    ##  Class :character               Class :character      
    ##  Mode  :character               Mode  :character      
    ##                                                       
    ##                                                       
    ##                                                       
    ##                                                       
    ##  histologic_grading_tier_category maximum_tumor_dimension
    ##  Length:1308                      Min.   : 0.300         
    ##  Class :character                 1st Qu.: 2.925         
    ##  Mode  :character                 Median : 3.500         
    ##                                   Mean   : 3.840         
    ##                                   3rd Qu.: 4.500         
    ##                                   Max.   :14.000         
    ##                                   NA's   :1138           
    ##  source_of_patient_death_reason tobacco_smoking_history
    ##  Length:1308                    Min.   :1.000          
    ##  Class :character               1st Qu.:1.000          
    ##  Mode  :character               Median :2.000          
    ##                                 Mean   :2.201          
    ##                                 3rd Qu.:3.000          
    ##                                 Max.   :5.000          
    ##                                 NA's   :1159           
    ##  year_of_tobacco_smoking_onset stopped_smoking_year number_pack_years_smoked
    ##  Min.   :1948                  Min.   :1952         Min.   : 0.30           
    ##  1st Qu.:1960                  1st Qu.:1980         1st Qu.:15.00           
    ##  Median :1971                  Median :1988         Median :25.00           
    ##  Mean   :1971                  Mean   :1990         Mean   :26.84           
    ##  3rd Qu.:1982                  3rd Qu.:2004         3rd Qu.:40.00           
    ##  Max.   :1993                  Max.   :2013         Max.   :75.00           
    ##  NA's   :1261                  NA's   :1258         NA's   :1251            
    ##  alcohol_history_documented alcoholic_exposure_category
    ##  Length:1308                Length:1308                
    ##  Class :character           Class :character           
    ##  Mode  :character           Mode  :character           
    ##                                                        
    ##                                                        
    ##                                                        
    ##                                                        
    ##  frequency_of_alcohol_consumption amount_of_alcohol_consumption_per_day
    ##  Min.   :0.500                    Min.   :0.500                        
    ##  1st Qu.:3.000                    1st Qu.:1.000                        
    ##  Median :6.500                    Median :1.000                        
    ##  Mean   :4.812                    Mean   :1.581                        
    ##  3rd Qu.:7.000                    3rd Qu.:2.000                        
    ##  Max.   :7.000                    Max.   :4.000                        
    ##  NA's   :1276                     NA's   :1277                         
    ##  history_of_diabetes days_to_diabetes_onset history_of_chronic_pancreatitis
    ##  Length:1308         Min.   :-9070.00       Length:1308                    
    ##  Class :character    1st Qu.: -163.00       Class :character               
    ##  Mode  :character    Median :  -52.00       Mode  :character               
    ##                      Mean   : -706.79                                      
    ##                      3rd Qu.:  -17.75                                      
    ##                      Max.   :  504.00                                      
    ##                      NA's   :1294                                          
    ##  days_to_pancreatitis_onset family_history_of_cancer relative_cancer_types
    ##  Min.   :-18029.0           Length:1308              Length:1308          
    ##  1st Qu.:  -231.5           Class :character         Class :character     
    ##  Median :   -71.0           Mode  :character         Mode  :character     
    ##  Mean   : -1744.2                                                         
    ##  3rd Qu.:   -38.0                                                         
    ##  Max.   :     1.0                                                         
    ##  NA's   :1297                                                             
    ##  history_prior_surgery_type_other
    ##  Length:1308                     
    ##  Class :character                
    ##  Mode  :character                
    ##                                  
    ##                                  
    ##                                  
    ## 

``` r
# Verificar valores nulos
print(table(is.na(df_TCGA)))
```

    ## 
    ##  FALSE   TRUE 
    ##  66141 135291

``` r
# Verificar valores nulos por columnas
# print(colSums(is.na(df_TCGA)))

# Verificar valores nulos por filas
# print(rowSums(is.na(df_TCGA)))
```

Let's identify the column that could serve as the main record ID of the dataset. 

``` r
# Visualizar los primeros 5 elementos de las siguientes columnas tipo ID
head(df_TCGA[,c("bcr_patient_barcode", "bcr_patient_uuid", "patient_id")], 5)
```

    ##   bcr_patient_barcode                     bcr_patient_uuid patient_id
    ## 1        TCGA-3L-AA1B A94E1279-A975-480A-93E9-7B1FF05CBCBF       AA1B
    ## 2        TCGA-4N-A93T 92554413-9EBC-4354-8E1B-9682F3A031D9       A93T
    ## 3        TCGA-4T-AA8H A5E14ADD-1552-4606-9FFE-3A03BCF76640       AA8H
    ## 4        TCGA-5M-AAT4 1136DD50-242A-4659-AAD4-C53F9E759BB3       AAT4
    ## 5        TCGA-5M-AAT6 CE00896A-F7D2-4123-BB95-24CB6E53FC32       AAT6

We observed three different identifiers: the full TCGA barcode, UUID (Universal Unique Identifier), and internal patient code, respectively. In this project, we will use the first one, as it is the standard patient-level identifier in TCGA. We also observed that it contains the internal patient code.

``` r
# Crear un subconjunto del dataset TCGA con variables relevantes desde una perspectiva oncologica
df_gastro <-subset(df_TCGA, select = c(bcr_patient_barcode, tumor_tissue_site, histological_type,
                               gender, vital_status, days_to_death, days_to_last_followup,
                               tissue_source_site,
                               age_at_initial_pathologic_diagnosis, person_neoplasm_cancer_status,
                               weight, height, residual_tumor, anatomic_neoplasm_subdivision,
                               number_of_lymphnodes_positive_by_he, number_of_lymphnodes_positive_by_ihc,
                               preoperative_pretreatment_cea_level, non_nodal_tumor_deposits,
                               circumferential_resection_margin, venous_invasion, lymphatic_invasion,
                               perineural_invasion_present, microsatellite_instability,
                               number_of_loci_tested, number_of_abnormal_loci, kras_gene_analysis_performed,
                               kras_mutation_found, kras_mutation_codon, braf_gene_analysis_performed,
                               braf_gene_analysis_result, synchronous_colon_cancer_present,
                               stage_event_pathologic_stage, stage_event_tnm_categories, cancer_type))
```

We then visualized in more detail the customized TCGA dataset with the 34 variables of interest.

``` r
# Usar función skimr()
skim(df_gastro)
```

    ## Warning: There was 1 warning in `dplyr::summarize()`.
    ## ℹ In argument: `dplyr::across(tidyselect::any_of(variable_names),
    ##   mangled_skimmers$funs)`.
    ## ℹ In group 0: .
    ## Caused by warning:
    ## ! There was 1 warning in `dplyr::summarize()`.
    ## ℹ In argument: `dplyr::across(tidyselect::any_of(variable_names),
    ##   mangled_skimmers$funs)`.
    ## Caused by warning in `inline_hist()`:
    ## ! Variable contains Inf or -Inf value(s) that were converted to NA.

|                                                  |           |
|:-------------------------------------------------|:----------|
| Name                                             | df_gastro |
| Number of rows                                   | 1308      |
| Number of columns                                | 34        |
| \_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_   |           |
| Column type frequency:                           |           |
| character                                        | 22        |
| numeric                                          | 12        |
| \_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_ |           |
| Group variables                                  | None      |

Data summary

**Variable type: character**

| skim_variable | n_missing | complete_rate | min | max | empty | n_unique | whitespace |
|:---|---:|---:|---:|---:|---:|---:|---:|
| bcr_patient_barcode | 1 | 1.00 | 0 | 12 | 1 | 1307 | 0 |
| tumor_tissue_site | 6 | 1.00 | 0 | 9 | 1 | 6 | 0 |
| histological_type | 17 | 0.99 | 0 | 65 | 1 | 19 | 0 |
| gender | 2 | 1.00 | 0 | 6 | 1 | 3 | 0 |
| vital_status | 636 | 0.51 | 4 | 5 | 0 | 4 | 0 |
| tissue_source_site | 1 | 1.00 | 0 | 2 | 1 | 99 | 0 |
| person_neoplasm_cancer_status | 194 | 0.85 | 0 | 10 | 1 | 3 | 0 |
| residual_tumor | 142 | 0.89 | 0 | 2 | 1 | 5 | 0 |
| anatomic_neoplasm_subdivision | 81 | 0.94 | 0 | 25 | 1 | 19 | 0 |
| non_nodal_tumor_deposits | 1010 | 0.23 | 0 | 3 | 1 | 3 | 0 |
| venous_invasion | 762 | 0.42 | 0 | 3 | 1 | 3 | 0 |
| lymphatic_invasion | 740 | 0.43 | 0 | 3 | 1 | 3 | 0 |
| perineural_invasion_present | 1030 | 0.21 | 0 | 3 | 1 | 3 | 0 |
| microsatellite_instability | 1190 | 0.09 | 0 | 3 | 1 | 3 | 0 |
| kras_gene_analysis_performed | 733 | 0.44 | 0 | 3 | 1 | 3 | 0 |
| kras_mutation_found | 1245 | 0.05 | 0 | 3 | 1 | 3 | 0 |
| braf_gene_analysis_performed | 747 | 0.43 | 0 | 3 | 1 | 3 | 0 |
| braf_gene_analysis_result | 1272 | 0.03 | 0 | 8 | 1 | 3 | 0 |
| synchronous_colon_cancer_present | 744 | 0.43 | 0 | 3 | 1 | 3 | 0 |
| stage_event_pathologic_stage | 52 | 0.96 | 0 | 10 | 1 | 15 | 0 |
| stage_event_tnm_categories | 3 | 1.00 | 0 | 9 | 1 | 141 | 0 |
| cancer_type | 1 | 1.00 | 0 | 9 | 1 | 6 | 0 |

**Variable type: numeric**

| skim_variable | n_missing | complete_rate | mean | sd | p0 | p25 | p50 | p75 | p100 | hist |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|:---|
| days_to_death | 1071 | 0.18 | 503.46 | 524.96 | 0.0 | 145.0 | 366.00 | 641.00 | 3042.0 | ▇▂▁▁▁ |
| days_to_last_followup | 205 | 0.84 | -Inf | NaN | -Inf | 0.0 | 485.00 | 942.00 | 4502.0 | ▇▃▁▁▁ |
| age_at_initial_pathologic_diagnosis | 56 | 0.96 | 65.85 | 11.88 | 30.0 | 58.0 | 67.00 | 74.00 | 90.0 | ▁▃▆▇▃ |
| weight | 936 | 0.28 | 80.44 | 20.88 | 34.0 | 65.0 | 78.05 | 91.85 | 175.3 | ▃▇▃▁▁ |
| height | 957 | 0.27 | 168.87 | 11.64 | 80.3 | 162.0 | 170.00 | 176.00 | 193.0 | ▁▁▁▇▆ |
| number_of_lymphnodes_positive_by_he | 143 | 0.89 | 3.49 | 6.26 | 0.0 | 0.0 | 1.00 | 4.00 | 57.0 | ▇▁▁▁▁ |
| number_of_lymphnodes_positive_by_ihc | 1188 | 0.09 | 0.28 | 1.64 | 0.0 | 0.0 | 0.00 | 0.00 | 12.0 | ▇▁▁▁▁ |
| preoperative_pretreatment_cea_level | 906 | 0.31 | 65.07 | 468.84 | 0.0 | 1.7 | 3.16 | 8.98 | 7868.0 | ▇▁▁▁▁ |
| circumferential_resection_margin | 1185 | 0.09 | 22.96 | 28.31 | 0.0 | 2.5 | 13.00 | 30.00 | 165.0 | ▇▂▁▁▁ |
| number_of_loci_tested | 1236 | 0.06 | 4.94 | 2.14 | 0.0 | 5.0 | 5.00 | 5.00 | 10.0 | ▁▁▇▁▁ |
| number_of_abnormal_loci | 1237 | 0.05 | 0.59 | 1.70 | 0.0 | 0.0 | 0.00 | 0.00 | 9.0 | ▇▁▁▁▁ |
| kras_mutation_codon | 1278 | 0.02 | 13.83 | 8.92 | 12.0 | 12.0 | 12.00 | 12.00 | 61.0 | ▇▁▁▁▁ |

``` r
# Obtener el número de variables y sus nombres
dim(df_gastro) 
```

    ## [1] 1308   34

``` r
# Verificar valores nulos
print(table(is.na(df_gastro)))
```

    ## 
    ## FALSE  TRUE 
    ## 23465 21007

**Data Transformation**

The first transformation to the dataframe `df_gastro` consisted in creating a new patient survival variable `os_time` to identify those who are still alive and those who are not. 

``` r
# Visualizar algunas variables relacionadas a la supervivencia de pacientes
df_gastro[1:10, c("vital_status", "days_to_last_followup", "days_to_death")]
```

    ##    vital_status days_to_last_followup days_to_death
    ## 1          <NA>                   475            NA
    ## 2          <NA>                   146            NA
    ## 3          <NA>                   385            NA
    ## 4          Dead                  -Inf            49
    ## 5          Dead                  -Inf           290
    ## 6          <NA>                  1200            NA
    ## 7          <NA>                   775            NA
    ## 8          Dead                  1126          1126
    ## 9          <NA>                  1419            NA
    ## 10         <NA>                  1331            NA

where,

- **vital_status**: patient’s vital status at the time of follow-up closure
- **days_to_last_followup**: number of days from diagnosis to the last clinical contact with the patient
- **days_to_death**: number of days from diagnosis to death

The `os_time` variable is a combination of `days_to_last_followup` and `days_to_death` where NA values of `days_to_death` are imputed by those in `days_to_last_followup`.

``` r
# Crear la variable de supervivencia de pacientes (los que siguen vivos y los que no)
df_gastro$os_time <- ifelse(
  !is.na(df_gastro$days_to_death), #si hay un valor que no es nulo en la variable "days_to_death"
  df_gastro$days_to_death,
  df_gastro$days_to_last_followup #en caso contrario, el valor es lo que aparece en la columna "days_to_last_followup"
)

# Transformar los valores infinitos de os_time en NA
df_gastro$os_time[is.infinite(df_gastro$os_time)] <- NA
```

In TCGA clinical data, `vital_status` is sometimes missing or incorrectly analyzed, but if the patient follow-up data are present, it can be assumed that the patient was alive at their last follow-up. This assumption was used in this project, although a more conservative approach would be to exclude those patients.

``` r
# Asignar "Alive" a los NAs de vital_status si la variable "days_to_last_followup" tiene un valor valido.
df_gastro$vital_status[is.na(df_gastro$vital_status) 
                    & !is.na(df_gastro$days_to_last_followup)] <- "Alive"  

# Definir os_event, si el paciente está vivo con 0, caso contrario con 1. 
df_gastro$os_event <- ifelse(df_gastro$vital_status == "Dead", 1, 0) 

# Eliminar la especificación TCGA del tipo de cáncer para simplificar.
df_gastro$cancer_type <-gsub("TCGA-","", df_gastro$cancer_type) #ref. The R book page 124

# Eliminar la especificación Stage para simplificar.
df_gastro$stage_event_pathologic_stage <-gsub("Stage ","", df_gastro$stage_event_pathologic_stage) #ref. The R book page 124

# Estructura del conjunto de datos y resumen estadístico
str(df_gastro)
```

    ## 'data.frame':    1308 obs. of  36 variables:
    ##  $ bcr_patient_barcode                 : chr  "TCGA-3L-AA1B" "TCGA-4N-A93T" "TCGA-4T-AA8H" "TCGA-5M-AAT4" ...
    ##  $ tumor_tissue_site                   : chr  "Colon" "Colon" "Colon" "Colon" ...
    ##  $ histological_type                   : chr  "Colon Adenocarcinoma" "Colon Adenocarcinoma" "Colon Mucinous Adenocarcinoma" "Colon Adenocarcinoma" ...
    ##  $ gender                              : chr  "FEMALE" "MALE" "FEMALE" "MALE" ...
    ##  $ vital_status                        : chr  "Alive" "Alive" "Alive" "Dead" ...
    ##  $ days_to_death                       : int  NA NA NA 49 290 NA NA 1126 NA NA ...
    ##  $ days_to_last_followup               : num  475 146 385 -Inf -Inf ...
    ##  $ tissue_source_site                  : chr  "3L" "4N" "4T" "5M" ...
    ##  $ age_at_initial_pathologic_diagnosis : int  61 67 42 74 40 76 45 85 82 71 ...
    ##  $ person_neoplasm_cancer_status       : chr  "TUMOR FREE" "WITH TUMOR" "TUMOR FREE" "WITH TUMOR" ...
    ##  $ weight                              : num  63.3 134 108 NA 99.1 ...
    ##  $ height                              : num  173 168 168 NA 162 ...
    ##  $ residual_tumor                      : chr  "R0" "R0" "R0" "R0" ...
    ##  $ anatomic_neoplasm_subdivision       : chr  "Cecum" "Ascending Colon" "Descending Colon" "Ascending Colon" ...
    ##  $ number_of_lymphnodes_positive_by_he : int  0 NA 0 0 10 0 NA 3 1 7 ...
    ##  $ number_of_lymphnodes_positive_by_ihc: int  0 2 NA 0 0 0 NA NA NA 0 ...
    ##  $ preoperative_pretreatment_cea_level : num  NA 2 NA 550 2.61 2.91 NA 17.4 3.4 32.8 ...
    ##  $ non_nodal_tumor_deposits            : chr  "NO" "YES" "NO" "NO" ...
    ##  $ circumferential_resection_margin    : num  NA 30 20 NA NA NA NA NA NA NA ...
    ##  $ venous_invasion                     : chr  "NO" "NO" "NO" "YES" ...
    ##  $ lymphatic_invasion                  : chr  "NO" "NO" "NO" NA ...
    ##  $ perineural_invasion_present         : chr  "NO" "NO" "NO" NA ...
    ##  $ microsatellite_instability          : chr  "NO" NA "NO" NA ...
    ##  $ number_of_loci_tested               : int  NA NA NA NA NA NA NA NA NA NA ...
    ##  $ number_of_abnormal_loci             : int  NA NA NA NA NA NA NA NA NA NA ...
    ##  $ kras_gene_analysis_performed        : chr  "NO" "NO" "NO" "NO" ...
    ##  $ kras_mutation_found                 : chr  NA NA NA NA ...
    ##  $ kras_mutation_codon                 : int  NA NA NA NA NA NA NA NA NA NA ...
    ##  $ braf_gene_analysis_performed        : chr  "NO" "NO" "NO" "NO" ...
    ##  $ braf_gene_analysis_result           : chr  NA NA NA NA ...
    ##  $ synchronous_colon_cancer_present    : chr  "NO" "YES" "NO" "NO" ...
    ##  $ stage_event_pathologic_stage        : chr  "I" "IIIB" "IIA" "IV" ...
    ##  $ stage_event_tnm_categories          : chr  "T2N0M0" "T4aN1bM0" "T3N0MX" "T3N0M1b" ...
    ##  $ cancer_type                         : chr  "COAD" "COAD" "COAD" "COAD" ...
    ##  $ os_time                             : num  475 146 385 49 290 ...
    ##  $ os_event                            : num  0 0 0 1 1 0 0 1 0 0 ...

### 2.2 Target Questions

1.  What is the correlation between the pathological stage of the tumor at diagnosis and survival?

2.  What survival differences are observed by patient gender and age?

3.  Are there significant differences in survival according to the type of gastrointestinal cancer?

4.  What is the impact of venous, lymphatic, or perineural invasion on survival?

---

## 3 Exploratory Data Analysis

### 3.1 Descriptive Statistics and Visualization

In this section we implement an exploratory and visualization analysis to address the target questions from section 2.2.

``` r
# Convertir a factor las variables tipo `chr`
df_gastro$gender<- factor(df_gastro$gender)
df_gastro$residual_tumor<- factor(df_gastro$residual_tumor)
df_gastro$venous_invasion<- factor(df_gastro$venous_invasion)
df_gastro$lymphatic_invasion<- factor(df_gastro$lymphatic_invasion)
df_gastro$perineural_invasion_present<- factor(df_gastro$perineural_invasion_present)
df_gastro$microsatellite_instability<- factor(df_gastro$microsatellite_instability)
df_gastro$kras_mutation_found<- factor(df_gastro$kras_mutation_found)
df_gastro$stage_event_pathologic_stage<- factor(df_gastro$stage_event_pathologic_stage)
df_gastro$cancer_type<- factor(df_gastro$cancer_type)
```

``` r
# Visualizar la distribucion de edades en los distintos tipos de cancer
ggplot(data=subset(df_gastro, !is.na(age_at_initial_pathologic_diagnosis)), aes(x = cancer_type, y = age_at_initial_pathologic_diagnosis, fill = cancer_type))+  
geom_boxplot(col = 'black')+
  labs(title = "Edad Paciente al Momento del Diagnóstico", x ="Tipo de Cáncer", y = "Edad (años)") +
  theme_minimal()+
  scale_fill_brewer(palette="PRGn")+
  theme(plot.title = element_text(size=16, color='Darkblue', face='bold', hjust = 0.5))
```

<figure>
<img
src="../images/unnamed-chunk-10-1.png"
style="width:100.0%" height="500"/>
</figure>

#### Question 1: What is the correlation between the pathological stage of the tumor at diagnosis and survival?

``` r
# Correlacion entre el estadio tumoral y la supervivencia de los pacientes
ggplot(data=subset(df_gastro, !is.na(stage_event_pathologic_stage)), aes(x = stage_event_pathologic_stage, y = os_time, fill = stage_event_pathologic_stage))+  
geom_boxplot(col = 'black')+
  labs(title = "Supervivencia en función del Estadio Tumoral", x ="Estadio Tumoral", y = "Supervivencia (días)") +
  theme_minimal()+
  #scale_fill_brewer(palette="PRGn")+
  theme(plot.title = element_text(size=16, color='Darkblue', face='bold', hjust = 0.5))
```

    ## Warning: Removed 225 rows containing non-finite outside the scale range
    ## (`stat_boxplot()`).

<figure>
<img
src="../images/unnamed-chunk-11-1.png"
style="width:100.0%" height="500"/>
</figure>

The survival plot as a function of tumor stage shows that survival progressively decreases as tumor stage increases, as indicated by Singh et al. (https://doi.org/10.5114/pg.2024.141834). Hence the critical importance of early diagnosis in cancer patients. To provide stronger statistical support for this observation, a more rigorous analysis using a Cox proportional-hazards model would be needed, which is beyond the scope of this project.

#### Question 2: Are there significant differences in survival according to the type of gastrointestinal cancer?

``` r
# Supervivencia en función del tipo de cáncer 
ggplot(data=subset(df_gastro, !is.na(os_time)), aes(x = cancer_type, y = os_time, fill = cancer_type))+ 
geom_boxplot(col = 'black')+
  labs(title = "Supervivencia en función del tipo de cáncer ", x ="Tipo Cáncer", y = "Superviviencia (días)") +
  theme_minimal()+
  #scale_fill_brewer(palette="RdBu")+
  theme(plot.title = element_text(size=16, color='Darkblue', face='bold', hjust = 0.5))
```

<figure>
<img
src="../images/unnamed-chunk-12-1.png"
style="width:100.0%" height="500"/>
</figure>

The survival plot as a function of cancer type shows that the mean expected survival in PAAD patients is lower than that observed in other cancer types, consistent with the literature on pancreatic ductal adenocarcinoma (PDAC). Singh et al. (https://doi.org/10.5114/pg.2024.141834) point out that PDAC patients have the worst prognosis, with a five-year survival rate barely reaching 12–13%. In contrast, other gastrointestinal cancers such as colorectal cancer have a more favorable prognosis, partly thanks to early detection programs and the development of more specific and effective therapies.

To identify statistically significant differences between cancer types and patient survival time, we implemented the following tests: ANOVA test, validation of ANOVA residual normality assumptions, Kruskal-Wallis test, and Kaplan-Meier survival analysis.

##### ANOVA Test

An ANOVA test was implemented assuming normality, then residual normality was checked to analyze whether there are statistically significant differences between cancer types.

``` r
# Eliminar NAs en variables de interés y crear un nuevo dataframe 
df_anova <- df_gastro[!is.na(df_gastro$os_time) & !is.na(df_gastro$cancer_type), ]

# Implementar test ANOVA variables os_time y cancer_type 
modelo_anova <- aov(os_time ~ cancer_type, data = df_anova)
summary(modelo_anova)
```

    ##               Df    Sum Sq Mean Sq F value   Pr(>F)    
    ## cancer_type    4  17419280 4354820   10.36 3.12e-08 ***
    ## Residuals   1052 442183622  420327                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

From the results, the p-value is 3.12e-08 ($p < 0.05$), indicating that there are significant differences between at least two groups. Furthermore, the F-value is 10.36 (> 1), and the higher this value, the stronger the evidence of significant differences between groups. Since ANOVA is significant, we next identified which groups differ from each other using TukeyHSD post-hoc test.

``` r
# Aplicar test TukeyHSD
TukeyHSD(modelo_anova)
```

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = os_time ~ cancer_type, data = df_anova)
    ## 
    ## $cancer_type
    ##                 diff       lwr        upr     p adj
    ## COAD-CHOL  139.08982 -161.7348  439.91445 0.7137666
    ## PAAD-CHOL -147.08580 -469.7097  175.53807 0.7244215
    ## READ-CHOL   89.86958 -235.1946  414.93381 0.9430886
    ## STAD-CHOL -127.06130 -430.0807  175.95811 0.7820316
    ## PAAD-COAD -286.17562 -457.6440 -114.70723 0.0000559
    ## READ-COAD  -49.22024 -225.2373  126.79680 0.9407864
    ## STAD-COAD -266.15112 -397.0557 -135.24654 0.0000003
    ## READ-PAAD  236.95538   25.8329  448.07786 0.0188148
    ## STAD-PAAD   20.02450 -155.2659  195.31487 0.9979407
    ## STAD-READ -216.93088 -396.6732  -37.18855 0.0089008

From the `TukeyHSD post-hoc test` results, the adjusted p-values 0.0000003, 0.0000559, 0.0089008, and 0.0188148 (all < 0.05) show that the greatest significant differences are between groups STAD-COAD, PAAD-COAD, STAD-READ, and READ-PAAD, respectively.

##### Validation of Model Residual Normality Assumptions

To validate ANOVA results, we checked the following assumptions: normality of residuals, homoscedasticity, and independence.

``` r
# Validar test ANOVA
plot(modelo_anova)
```

<figure>
<img
src="../images/unnamed-chunk-15-1.png"
style="width:100.0%" height="500"/>
<img
src="../images/unnamed-chunk-15-2.png"
style="width:100.0%" height="500"/>
<img
src="../images/unnamed-chunk-15-3.png"
style="width:100.0%" height="500"/>
<img
src="../images/unnamed-chunk-15-4.png"
style="width:100.0%" height="500"/>
</figure>

**Normality**

Observing the Q-Q residuals plot (Quantile-Quantile), we can compare the theoretical quantiles of a normal distribution with the actual residual quantiles. In our case, points tend to align close to the diagonal line, especially in the middle part. However, there is evidence of asymmetry in the right tail where an S-shaped curve appears. This may affect the validity of p-values and post-hoc comparisons.

**Homoscedasticity**

The homoscedasticity assumption requires that residual variance be constant for all fitted values of the model for each group. From the Scale-Location plot, there is no evident relationship between residuals and fitted values (the mean of each group). The band is almost horizontal and flat. Thus, we can assume homogeneity of variances.

Some ideas for answering this question were taken from the post **ANOVA in R** section [Another method to test normality and homogeneity](https://statsandr.com/blog/anova-in-r/#another-method-to-test-normality-and-homogeneity).

**Independence**

Since the data are not temporal, we can assume a priori there is no autocorrelation. Independence is assumed by design.

##### Kruskal-Wallis Test

If we assume that the normality assumption is not met, we can apply the following non-parametric tests: `Kruskal-Wallis` and `Dunn’s post-hoc test` to validate our results.

``` r
# Aplicar test Kruskal-Wallis
kruskal.test(os_time ~ cancer_type, data = df_anova)
```

    ## 
    ##  Kruskal-Wallis rank sum test
    ## 
    ## data:  os_time by cancer_type
    ## Kruskal-Wallis chi-squared = 50.438, df = 4, p-value = 2.925e-10

``` r
# Aplicar post-hoc test Dunn
dunn.test(df_anova$os_time, df_anova$cancer_type, method = "bonferroni")
```

    ##   Kruskal-Wallis rank sum test
    ## 
    ## data: x and group
    ## Kruskal-Wallis chi-squared = 50.4378, df = 4, p-value = 0
    ## 
    ## 
    ##                            Comparison of x by group                            
    ##                                  (Bonferroni)                                  
    ## Col Mean-|
    ## Row Mean |       CHOL       COAD       PAAD       READ
    ## ---------+--------------------------------------------
    ##     COAD |  -0.670817
    ##          |     1.0000
    ##          |
    ##     PAAD |   1.809028   4.580635
    ##          |     0.3522    0.0000*
    ##          |
    ##     READ |  -0.718359  -0.180178  -3.870495
    ##          |     1.0000     1.0000    0.0005*
    ##          |
    ##     STAD |   1.892153   5.921548  -0.058624   4.489049
    ##          |     0.2924    0.0000*     1.0000    0.0000*
    ## 
    ## alpha = 0.05
    ## Reject Ho if p <= alpha/2

The non-parametric tests confirmed that the following groups have statistically significant differences: PAAD-COAD, READ-PAAD, STAD-COAD, and STAD-READ.

##### Kaplan-Meier Survival Studies

``` r
# Implementar Kaplan-Meier por tipo de cancer
# Crear objeto de supervivencia
survival_obj <- Surv(time = df_gastro$os_time, event=df_gastro$os_event)

# Ajustar el model Kaplan-Meier por género
ajuste_cancer_type <- survfit(survival_obj ~ cancer_type, type="kaplan-meier", data=df_gastro)

# Para evitar nombres arbitrarios para el tipo de cáncer necesitamos utilizar los nombres de grupo que figuran en los datos
nombres_tipos <- names(ajuste_cancer_type$strata)

# Visualizar gráfico
plot(ajuste_cancer_type, ylab = "Probabilidad", xlab='Tiempo (días)',
     mark.time = TRUE, col = hue_pal ()(length(nombres_tipos)), main= "Supervivencia en Función del Tipo de Cancer")

legend("bottomleft", legend = nombres_tipos,
       fill = hue_pal ()(length(nombres_tipos)), bty = "n")
```

<figure>
<img
src="../images/unnamed-chunk-17-1.png"
style="width:100.0%" height="500"/>
</figure>

``` r
# Comparamos la supervivencia entre categorías de tipo de cancer
comparar_cancer_type<- survdiff(survival_obj ~ cancer_type, data=df_gastro)
comparar_cancer_type
```

    ## Call:
    ## survdiff(formula = survival_obj ~ cancer_type, data = df_gastro)
    ## 
    ## n=1057, 251 observations deleted due to missingness.
    ## 
    ##                    N Observed Expected (O-E)^2/E (O-E)^2/V
    ## cancer_type=CHOL  38       20     8.76      14.4     15.01
    ## cancer_type=COAD 397       57   101.75      19.7     34.75
    ## cancer_type=PAAD 146       66    28.17      50.8     57.88
    ## cancer_type=READ 136        9    33.79      18.2     21.23
    ## cancer_type=STAD 340       87    66.53       6.3      8.79
    ## 
    ##  Chisq= 111  on 4 degrees of freedom, p= <2e-16

Kaplan-Meier survival test results also show significant differences in survival according to cancer type, with a Chi-square value of 111, 4 degrees of freedom, and a p-value below 2e-16.

As previously reported in the literature, patients with pancreatic adenocarcinoma (PAAC) and cholangiocarcinoma have the worst prognoses. In our analysis, we confirmed that patients with pancreatic and biliary tract cancers have a number of deaths clearly higher than expected, reflecting a significantly lower overall survival compared to other tumor types such as colon or rectal cancer.

#### Question 3: What survival differences are observed by patient gender and age?

Before implementing survival analysis by age, it is recommended to categorize the variable `age_at_initial_pathologic_diagnosis` ensuring a balanced distribution between groups.

``` r
# Categorizar la variable age_at_initial_pathologic_diagnosis
df_gastro$age_category <- cut(df_gastro$age_at_initial_pathologic_diagnosis,
                       breaks = c(-Inf,59,69,79,Inf),
                       labels = c("Personas_60-","Personas_60","Personas_70","Personas_80+"),
                       include.lowest = FALSE)

# Verificar distribución de categorias
table(df_gastro$age_category)
```

    ## 
    ## Personas_60-  Personas_60  Personas_70 Personas_80+ 
    ##          368          363          368          153

``` r
# Distribucion de edades de pacientes por género
ggplot(data=subset(df_gastro, !is.na(age_category)), aes(x= age_category, fill = gender))+   
  geom_bar(position = "dodge", alpha = 0.7)+
  labs(title= "Distribución de edades de pacientes por género", 
       x= "Edad (años)",
       y= "Frecuencia")+
  theme(plot.title = element_text(size=16, color='Darkblue', face='bold', hjust = 0.5))
```

<figure>
<img
src="../images/unnamed-chunk-20-1.png"
style="width:100.0%" height="500"/>
</figure>

##### Kaplan-Meier Survival Studies

Kaplan-Meier survival studies by patient gender and age are shown below. For implementation, the following resources were consulted: [Survival
analysis with low-dimensional input
data](https://ocbe-uio.github.io/survomics/survomics.html#survival-analysis-with-low-dimensional-input-data),
[Analysis of Cancer Genome Atlas in
R](https://www.costalab.org/wp-content/uploads/2020/11/R_class_D3.html#5_Survival_Analysis)
and [The R Book, 3rd Edition by Elinor Jones, Simon Harden, Michael J.
Crawley, cap. 15 Survial
Analysis](https://learning.oreilly.com/library/view/the-r-book/9781119634324/c15.xhtml#head-2-134).

``` r
# Implementar Kaplan-Meier por género 
# Crear objeto de supervivencia
survival_obj <- Surv(time = df_gastro$os_time, event=df_gastro$os_event)

# Ajustar el model Kaplan-Meier por género
ajuste_genero <- survfit(survival_obj ~ gender, type="kaplan-meier", data=df_gastro)

# Visualizar gráfico
plot(ajuste_genero, ylab = "Probabilidad", xlab='Tiempo (días)',
     mark.time = TRUE, col = hue_pal ()(2)[1:2], main= "Supervivencia en Función del Género")

legend("bottomleft", legend = c("Female", "Male"),
       fill = hue_pal ()(2)[1:2], bty = "n")
```

<figure>
<img
src="../images/unnamed-chunk-21-1.png"
style="width:100.0%" height="500"/>
</figure>

``` r
# Comparamos la supervivencia entre géneros
comparar_generos <- survdiff(survival_obj ~ gender, data=df_gastro)
comparar_generos
```

    ## Call:
    ## survdiff(formula = survival_obj ~ gender, data = df_gastro)
    ## 
    ## n=1057, 251 observations deleted due to missingness.
    ## 
    ##                 N Observed Expected (O-E)^2/E (O-E)^2/V
    ## gender=FEMALE 433       86      102      2.54      4.44
    ## gender=MALE   624      153      137      1.89      4.44
    ## 
    ##  Chisq= 4.4  on 1 degrees of freedom, p= 0.04

Regarding survival and gender, in the Kaplan-Meier plot, the curve for women reaches a stable value of approximately 0.65 at 1,800 days, indicating that about 65% survive beyond that time. In contrast, the men’s curve shows a longer decline, stabilizing at different points, then dropping to about 0.35 at 3,000 days. This implies that only 35% of men reach that survival time threshold. The log-rank test (Chi² = 4.4, p = 0.04) supports this visual inspection, where the observed differences are statistically significant. Women have fewer events (deaths) than expected, suggesting better clinical outcomes compared to men.

``` r
# Implementar Kaplan-Meier por edad
# Ajustar el model Kaplan-Meier por edad
ajuste_edad <- survfit(survival_obj ~ age_category, type="kaplan-meier", data=df_gastro)

# Visualizar gráfico
plot(ajuste_edad, ylab = "Probabilidad", xlab='Tiempo (días)',
     mark.time = TRUE, col = hue_pal ()(4)[1:4], main= "Superviviencia en Función de la Edad")

legend("bottomleft", legend = c("Personas_60-","Personas_60","Personas_70","Personas_80+"),
       fill = hue_pal ()(4)[1:4], bty = "n")
```

<figure>
<img
src="../images/unnamed-chunk-23-1.png"
style="width:100.0%" height="500"/>
</figure>

``` r
# Comparamos la supervivencia entre categorías de edades
comparar_edades <- survdiff(survival_obj ~ age_category, data=df_gastro)
comparar_edades
```

    ## Call:
    ## survdiff(formula = survival_obj ~ age_category, data = df_gastro)
    ## 
    ## n=1016, 292 observations deleted due to missingness.
    ## 
    ##                             N Observed Expected (O-E)^2/E (O-E)^2/V
    ## age_category=Personas_60- 311       51     70.2     5.268     7.789
    ## age_category=Personas_60  298       58     64.7     0.695     0.989
    ## age_category=Personas_70  298       78     65.3     2.482     3.540
    ## age_category=Personas_80+ 109       32     18.8     9.290    10.197
    ## 
    ##  Chisq= 17.8  on 3 degrees of freedom, p= 5e-04

Regarding survival and age, in the Kaplan-Meier plot, younger groups maintain higher and more sustained survival curves over time, while curves for those over 70 and 80 years show a sharper decline. The log-rank test (Chi² = 17.8, df = 3, p = 0.0005) confirms that there are statistically significant differences between age groups. In particular, patients over 80 years have an observed mortality much higher than expected (32 observed vs 18.8 expected). By contrast, patients under 60 have a lower observed mortality than expected (51 vs 70.2), indicating better relative survival.

#### Question 4: What is the impact of venous, lymphatic, or perineural invasion on survival?

``` r
# Reorganizar los datos 
df_long <- df_gastro %>%
  select(os_time, perineural_invasion_present, venous_invasion, lymphatic_invasion) %>%
  pivot_longer(
    cols = c(perineural_invasion_present, venous_invasion, lymphatic_invasion),
    names_to = "tipo_invasion",
    values_to = "presencia"
  ) %>%
  filter(!is.na(presencia))  

# Crear el gráfico combinado
ggplot(df_long, aes(x = presencia, y = os_time, fill = presencia)) +
  geom_boxplot(color = "black") +
  facet_wrap(~ tipo_invasion, scales = "free_x") +
  labs(
    title = "Supervivencia según tipo de invasión tumoral",
    x = "Presencia de invasión",
    y = "Supervivencia (días)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, color = 'Darkblue', face = 'bold', hjust = 0.5)
  )
```

    ## Warning: Removed 195 rows containing non-finite outside the scale range
    ## (`stat_boxplot()`).

<figure>
<img
src="../images/unnamed-chunk-25-1.png"
style="width:100.0%" height="500"/>
</figure>

##### Kaplan-Meier Survival Studies

``` r
# Implementar Kaplan-Meier invásion venosa
# Crear objeto de supervivencia
survival_obj <- Surv(time = df_gastro$os_time, event=df_gastro$os_event)

# Ajustar el model Kaplan-Meier por género
ajuste_invasion_venosa<- survfit(survival_obj ~ venous_invasion, type="kaplan-meier", data=df_gastro)

# Visualizar gráfico
plot(ajuste_invasion_venosa, ylab = "Probabilidad", xlab='Tiempo (días)',
     mark.time = TRUE, col = hue_pal ()(2)[1:2], main= "Supervivencia en función de la invásion venosa")

legend("bottomleft", legend = c("Present", "Absent"),
       fill = hue_pal ()(2)[1:2], bty = "n")
```

<figure>
<img
src="../images/unnamed-chunk-26-1.png"
style="width:100.0%" height="500"/>
</figure>

``` r
# Comparamos la supervivencia entre pacientes con o sin invásion
comparar_invasion_venosa <- survdiff(survival_obj ~ venous_invasion, data=df_gastro)
comparar_invasion_venosa
```

    ## Call:
    ## survdiff(formula = survival_obj ~ venous_invasion, data = df_gastro)
    ## 
    ## n=467, 841 observations deleted due to missingness.
    ## 
    ##                       N Observed Expected (O-E)^2/E (O-E)^2/V
    ## venous_invasion=NO  361       29     42.1      4.08        20
    ## venous_invasion=YES 106       24     10.9     15.77        20
    ## 
    ##  Chisq= 20  on 1 degrees of freedom, p= 8e-06

``` r
# Implementar Kaplan-Meier invasion perineural
# Crear objeto de supervivencia
survival_obj <- Surv(time = df_gastro$os_time, event=df_gastro$os_event)

# Ajustar el model Kaplan-Meier por género
ajuste_invasion_perineural<- survfit(survival_obj ~ perineural_invasion_present, type="kaplan-meier", data=df_gastro)

# Visualizar gráfico
plot(ajuste_invasion_perineural, ylab = "Probabilidad", xlab='Tiempo (días)',
     mark.time = TRUE, col = hue_pal ()(2)[1:2], main= "Supervivencia en Función de la Invasión perineural")

legend("bottomleft", legend = c("Present", "Absent"),
       fill = hue_pal ()(2)[1:2], bty = "n")
```

<figure>
<img
src="../images/unnamed-chunk-28-1.png"
style="width:100.0%" height="500"/>
</figure>

``` r
# Comparamos la supervivencia entre los pacientes con y sin invasión perineural
comparar_invasion_perineural <- survdiff(survival_obj ~ perineural_invasion_present, data=df_gastro)
comparar_invasion_perineural
```

    ## Call:
    ## survdiff(formula = survival_obj ~ perineural_invasion_present, 
    ##     data = df_gastro)
    ## 
    ## n=247, 1061 observations deleted due to missingness.
    ## 
    ##                                   N Observed Expected (O-E)^2/E (O-E)^2/V
    ## perineural_invasion_present=NO  182       25    31.74      1.43      6.43
    ## perineural_invasion_present=YES  65       16     9.26      4.91      6.43
    ## 
    ##  Chisq= 6.4  on 1 degrees of freedom, p= 0.01

``` r
# Implementar Kaplan-Meier Invasión linfática
# Crear objeto de supervivencia
survival_obj <- Surv(time = df_gastro$os_time, event=df_gastro$os_event)

# Ajustar el model Kaplan-Meier por género
ajuste_invasion_linfatica<- survfit(survival_obj ~ lymphatic_invasion, type="kaplan-meier", data=df_gastro)

# Visualizar gráfico
plot(ajuste_invasion_linfatica, ylab = "Probabilidad", xlab='Tiempo (días)',
     mark.time = TRUE, col = hue_pal ()(2)[1:2], main= "Supervivencia en Función de la Invasión Linfática")

legend("bottomleft", legend = c("Present", "Absent"),
       fill = hue_pal ()(2)[1:2], bty = "n")
```

<figure>
<img
src="../images/unnamed-chunk-30-1.png"
style="width:100.0%" height="500"/>
</figure>

``` r
# Comparamos la supervivencia entre los pacientes con y sin invasión perineural
comparar_invasion_linfatica <- survdiff(survival_obj ~ lymphatic_invasion, data=df_gastro)
comparar_invasion_linfatica
```

    ## Call:
    ## survdiff(formula = survival_obj ~ lymphatic_invasion, data = df_gastro)
    ## 
    ## n=483, 825 observations deleted due to missingness.
    ## 
    ##                          N Observed Expected (O-E)^2/E (O-E)^2/V
    ## lymphatic_invasion=NO  305       25       36      3.37      10.3
    ## lymphatic_invasion=YES 178       29       18      6.75      10.3
    ## 
    ##  Chisq= 10.3  on 1 degrees of freedom, p= 0.001

Survival analysis by histopathological variables revealed significant differences associated with venous, lymphatic, and perineural invasion. These three features are known to reflect greater tumor aggressiveness, and our results statistically confirmed this.

First, venous invasion showed a highly significant association with survival (Chi-square = 20, p = 8e-06), with a number of deaths much higher than expected in patients with positive invasion. Similarly, lymphatic invasion was associated with worse prognosis (Chi-square = 10.3, p = 0.001), with higher mortality in the invasion-positive group.

Finally, perineural invasion was significantly associated with reduced survival (Chi-square = 6.4, p = 0.01), although with a smaller cohort size.

These results support the clinical value of these variables as potential negative prognostic markers in gastrointestinal cancer patients. They also highlight the importance of including them in risk stratification and therapeutic decision-making.

#### Question 5: What is the BMI distribution of patients by gender? 

BMI (Body Mass Index) is a widely used indicator in oncology to assess metabolic risk, prognosis, and treatment response, according to the IARC Working Group report [Body Fatness and Cancer — Viewpoint of the IARC Working Group](https://www.nejm.org/doi/10.1056/NEJMsr1606602?url_ver=Z39.88-2003&rfr_id=ori:rid:crossref.org&rfr_dat=cr_pub%20%200www.ncbi.nlm.nih.gov). BMI at diagnosis may reflect comorbidities (diabetes, hypertension) that can influence cancer progression. For this reason, we included BMI as a key variable for patient stratification, referencing NIH-proposed categories. For more details, consult [Obesity and Cancer](https://www.cancer.gov/espanol/cancer/causas-prevencion/riesgo/obesidad/hoja-informativa-obesidad).

``` r
# Función para calcular el IMC de un paciente 
calcular_imc <- function(peso, altura){
  imc <- peso / (altura/100)^2
  return (imc)
}

# calcular_imc <- function(peso, altura){
#   tryCatch(
#     {
#       # Convertir en valores numéricos por default
#       #peso      <- as.numeric(peso)
#       #altura    <- as.numeric(altura)
#       
#       # Validar rango de valores de entrada. 
#       if(altura <= 0.0){
#         stop("Error: La variable altura no debe ser negativa o cero")
#       }
#       
#       imc <- peso / (altura^2)
#       return(imc)
#     },
#     error = function(e){
#         message("Error en el cálculo: ", e$message)
#         return(NA)
#     }
#   )
# }

# Crear la variable imc_at_initial_pathologic_diagnosis
df_gastro$imc <- calcular_imc(df_gastro$weight, df_gastro$height)

# Categorizar la variable imc según las categorías propuestas por NIH
df_gastro$imc_category <- cut(df_gastro$imc,
                       breaks = c(-Inf,18.5,25.0,30.0,Inf),
                       labels = c("Bajo_Peso","Normal","Sobrepeso","Obesidad"),
                       include.lowest = FALSE)

# Verificar distribución de categorias
table(df_gastro$imc_category)
```

    ## 
    ## Bajo_Peso    Normal Sobrepeso  Obesidad 
    ##         6       108       138        99

``` r
# Distribucion de IMC de pacientes por género
ggplot(data=subset(df_gastro, !is.na(imc_category)), aes(x= imc_category, fill = gender))+   
  geom_bar(position = "dodge", alpha = 0.7)+
  labs(title= "Distribución de IMC de pacientes por género", 
       x= "Edad (años)",
       y= "Frecuencia")+
  theme(plot.title = element_text(size=16, color='Darkblue', face='bold', hjust = 0.5))
```

<figure>
<img
src="../images/unnamed-chunk-33-1.png"
style="width:100.0%" height="500"/>
</figure>

### 3.2 Probability Study

This section proposes three probability problems related to the TCGA dataset.

#### Exercise A

In the customized TCGA dataset (`df_gastro`), approximately 14.15% of patients have pancreatic cancer, according to the `cancer_type` variable. If we take a random sample of 100 patients:

- What is the probability that exactly 5 patients have pancreatic cancer?

- What is the probability that at most 10 patients have pancreatic cancer?

``` r
# Verificar distribución de categorias por tipo de cáncer
frecuencias <- table(df_gastro$cancer_type)
frecuencias
```

    ## 
    ##      CHOL COAD PAAD READ STAD 
    ##    1   48  459  185  171  443

``` r
# Verificar porcentaje de distribución de categorías por tipo de cáncer
# porcentajes <- prop.table(frecuencias) * 100
porcentajes <- round(prop.table(frecuencias) * 100, 2)
porcentajes
```

    ## 
    ##        CHOL  COAD  PAAD  READ  STAD 
    ##  0.08  3.67 35.12 14.15 13.08 33.89

``` r
# Describir datos para calcular P(X = 5)
n <- 100
p <- 0.1415
x <- 5

# Aplicar distribución binomial para P(X = 5)
paciente_5 <- dbinom(x,n,p)

print(paste("La probabilidad de que exactamente 5 pacientes tengan cáncer de páncreas es:",
            round(paciente_5,4)*100,"%"))
```

    ## [1] "La probabilidad de que exactamente 5 pacientes tengan cáncer de páncreas es: 0.22 %"

``` r
# Describir datos para calcular P(X <= 10)
n <- 100
p <- 0.1415
x <- 10

# Aplicar distribución binomial para P(X <= 10)
paciente_10 <- pbinom(x, n, p)

print(paste("La probabilidad de que un máximo de 10 pacientes tengan cáncer de páncreas es:",
            round(paciente_10, 4)*100,"%"))
```

    ## [1] "La probabilidad de que un máximo de 10 pacientes tengan cáncer de páncreas es: 14.62 %"

#### Exercise B

Of the patients with perineural invasion, what is the probability that they have an advanced pathological stage, equal to or greater than type III?

This is a conditional probability problem where we are asked: P(Stage >= III | Perineural Invasion = YES)

``` r
# Visualizar variables categóricas de interés
print("Distribución de categorías para la variable perineural_invasion_present")
```

    ## [1] "Distribución de categorías para la variable perineural_invasion_present"

``` r
table(df_gastro$perineural_invasion_present)
```

    ## 
    ##      NO YES 
    ##   1 206  71

``` r
print("Distribución de categorías para la variable stage_event_pathologic_stage")
```

    ## [1] "Distribución de categorías para la variable stage_event_pathologic_stage"

``` r
table(df_gastro$stage_event_pathologic_stage)
```

    ## 
    ##         I   IA   IB   II  IIA  IIB  IIC  III IIIA IIIB IIIC   IV  IVA  IVB 
    ##    1  131   22   56   82  248  193    2   36   94  149   93  115   27    7

``` r
# Eliminar NAs en variables de interés y crear un nuevo dataframe 
df_estate_peri <- df_gastro[!is.na(df_gastro$perineural_invasion_present) & 
                      !is.na(df_gastro$stage_event_pathologic_stage), ]
```

To study the relationship between categorical variable groups, the next step is to build a 14 x 2 contingency table, showing counts for each invasion_perineural–stage combination, resulting in 28 observed frequencies.

``` r
# Crear tabla de contingencia 
tabla_conti <- table(df_estate_peri$perineural_invasion_present,
               df_estate_peri$stage_event_pathologic_stage)
print("Tabla de contingencia de variables categóricas")
```

    ## [1] "Tabla de contingencia de variables categóricas"

``` r
tabla_conti
```

    ##      
    ##           I IA IB II IIA IIB IIC III IIIA IIIB IIIC IV IVA IVB
    ##        1  0  0  0  0   0   0   0   0    0    0    0  0   0   0
    ##   NO   0 53  1  0  9  52   4   1   7    9   31   14  6  10   4
    ##   YES  0  8  0  0  4  13   5   0   1    0   13    4 10   6   3

``` r
# Calcular proporciones para cada grupo de pacientes con invasión perineural 
prop_conti <- prop.table(tabla_conti, margin = 1)  # margin = 1 obtiene proporciones por fila
print("Tabla de proporciones")
```

    ## [1] "Tabla de proporciones"

``` r
prop_conti
```

    ##      
    ##                             I          IA          IB          II         IIA
    ##       1.000000000 0.000000000 0.000000000 0.000000000 0.000000000 0.000000000
    ##   NO  0.000000000 0.263681592 0.004975124 0.000000000 0.044776119 0.258706468
    ##   YES 0.000000000 0.119402985 0.000000000 0.000000000 0.059701493 0.194029851
    ##      
    ##               IIB         IIC         III        IIIA        IIIB        IIIC
    ##       0.000000000 0.000000000 0.000000000 0.000000000 0.000000000 0.000000000
    ##   NO  0.019900498 0.004975124 0.034825871 0.044776119 0.154228856 0.069651741
    ##   YES 0.074626866 0.000000000 0.014925373 0.000000000 0.194029851 0.059701493
    ##      
    ##                IV         IVA         IVB
    ##       0.000000000 0.000000000 0.000000000
    ##   NO  0.029850746 0.049751244 0.019900498
    ##   YES 0.149253731 0.089552239 0.044776119

``` r
# Obtener la proporción solo para el grupo con invasión perineural = 'YES' y estadios > III
print("Tabla de proporciones para pacientes con invasión perineural = 'YES' y estadios > III")
```

    ## [1] "Tabla de proporciones para pacientes con invasión perineural = 'YES' y estadios > III"

``` r
prop_conti["YES", c("III", "IIIA", "IIIB", "IIIC", "IV", "IVA", "IVB")]
```

    ##        III       IIIA       IIIB       IIIC         IV        IVA        IVB 
    ## 0.01492537 0.00000000 0.19402985 0.05970149 0.14925373 0.08955224 0.04477612

``` r
# Obtener la probabilidad total conjunta
prob_conjunta <- sum(prop_conti["YES", c("III", "IIIA", "IIIB", "IIIC", "IV", "IVA", "IVB")])
paste("Probabilidad estadio igual o superior a tipo III:", round(prob_conjunta, 4)*100,"%")
```

    ## [1] "Probabilidad estadio igual o superior a tipo III: 55.22 %"

#### Exercise C

The mean age at diagnosis is 65.85 years with a standard deviation of 11.88 years. Assuming a normal distribution, what is the probability that a patient is between 60 and 70 years old?

In this problem we will use the normal distribution to obtain P(60 < X < 70)

``` r
# Describir datos para obtener P(60 < X < 70)
media  <- mean(df_gastro$age_at_initial_pathologic_diagnosis, na.rm = TRUE)
sigma  <- sd(df_gastro$age_at_initial_pathologic_diagnosis, na.rm = TRUE)
x1     <- 60
x2     <- 70

# Aplicar distribución normal para obtener la probabilidad P(60 < X < 70)
menor_a_70   <- pnorm(x2, media, sigma)
menor_a_60   <- pnorm(x1, media, sigma)
entre_70_60 <- menor_a_70 - menor_a_60

print(paste("La probabilidad de que un paciente tenga entre 60 y 70 años es:",
            round(entre_70_60, 4)*100,"%"))
```

    ## [1] "La probabilidad de que un paciente tenga entre 60 y 70 años es: 32.55 %"

---

## 4 Machine Learning Models

Since our goal is to identify factors associated with cancer patient survival, we opted to apply supervised learning models, as we have a well-defined dependent variable: overall survival (`os_event`).

We trained a supervised learning model using Support Vector Machine (SVM). For this, we selected a subset of clinical variables including age at diagnosis, gender, residual tumor type, pathological stage, number of positive nodes, cancer type, and presence of venous, lymphatic, or perineural invasion.

### 4.1 Supervised Model (SVM)

``` r
# Filtramos filas sin valores NA en variables clave
model_df_supervised <- na.omit(df_gastro[, c("age_at_initial_pathologic_diagnosis",
                                 "gender", "residual_tumor", "stage_event_pathologic_stage",
                                "number_of_lymphnodes_positive_by_he",
                                "cancer_type", 
                                "venous_invasion", "lymphatic_invasion",
                                "perineural_invasion_present",
                                "os_event")])

# Convertimos os_event en factor (para clasificación)
model_df_supervised$os_event <- as.factor(model_df_supervised$os_event)

# Dividimos el conjunto de datos en entrenamiento y prueba
set.seed(123)  # Para reproducibilidad
indices_entrenamiento <- sample(1:nrow(model_df_supervised), 0.5 * nrow(model_df_supervised))
conjunto_entrenamiento <- model_df_supervised[indices_entrenamiento, ]
conjunto_prueba <- model_df_supervised[-indices_entrenamiento, ]

# Cargamos el paquete "kernlab" para SVM
library(kernlab)
```

    ## Warning: package 'kernlab' was built under R version 4.4.1

    ## 
    ## Attaching package: 'kernlab'

    ## The following object is masked from 'package:ggplot2':
    ## 
    ##     alpha

    ## The following object is masked from 'package:scales':
    ## 
    ##     alpha

``` r
# Entrenamos un modelo SVM
modelo_svm <- ksvm(os_event ~ ., data = conjunto_entrenamiento, kernel = "rbfdot", type = "C-svc") #indicamos que el objetivo es clasificacion

# Realizar predicciones en el conjunto de prueba
predicciones <- predict(modelo_svm, newdata = conjunto_prueba)

# Calcular la matriz de confusión
confusion_matriz <- table(Real = conjunto_prueba$os_event, Predicción =
predicciones)

# Calcular la precisión
precision <- sum(diag(confusion_matriz)) / sum(confusion_matriz)
cat("Precisión del modelo SVM:", precision, "\n")
```

    ## Precisión del modelo SVM: 0.9354839

``` r
# Mostrar la matriz de confusión
confusion_matriz
```

    ##     Predicción
    ## Real  0  1
    ##    0 87  0
    ##    1  6  0

**Dataset Preprocessing**

Rows with NA values in the selected variables were filtered out. Additionally, `os_event` was converted into a factor to enable classification. The resulting dataset was randomly split into two subsets: training (50%) and testing (50%).

**Model Training**

The `ksvm()` function from the `kernlab` package was used to train an SVM model with an `rbfdot` radial kernel. The model was fitted using the training data and then evaluated on the test set.

**Results**

The model achieved an accuracy of 93.5% on the test set. However, the confusion matrix reveals that the model failed to correctly classify any deceased patients: all deceased patients were classified as alive. This indicates that the model is biased toward the majority class (alive patients), likely due to a strong imbalance in the alive/deceased distribution. In other words, while the model appears accurate overall, its ability to detect patients with poor prognosis is zero. To correct this imbalance, we would need to add more patients with poor prognosis.

### 4.2 PCA - Unsupervised Model 

In this part of the work we used an unsupervised learning model to explore whether patients’ clinical data show any natural pattern or grouping without telling the model what to predict. For this, we used Principal Component Analysis (PCA).

``` r
# # Filtramos filas sin valores NA en variables clave
model_df_unsupervised <- na.omit(df_gastro[, c("age_at_initial_pathologic_diagnosis",
                                 "gender", "residual_tumor", "stage_event_pathologic_stage",
                                "number_of_lymphnodes_positive_by_he",
                                "cancer_type", 
                                "venous_invasion", "lymphatic_invasion",
                                "perineural_invasion_present")])
 
# Transformamos las variables categoricas en numericas 
# Ref. "The model.matrix function" The R Book page 255
df_numeric <- model.matrix(~ gender + age_at_initial_pathologic_diagnosis + 
                              residual_tumor + 
                              number_of_lymphnodes_positive_by_he + 
                              venous_invasion + 
                              lymphatic_invasion + 
                              perineural_invasion_present + 
                              stage_event_pathologic_stage - 1, data = model_df_unsupervised)
 
# Aplicamos PCA
pca_resultado <- prcomp(df_numeric, center = TRUE)

# Resumen de la varianza explicada por cada variable
summary(pca_resultado)
```

    ## Importance of components:
    ##                           PC1    PC2     PC3    PC4     PC5     PC6     PC7
    ## Standard deviation     12.855 4.6032 0.81452 0.7006 0.52639 0.48330 0.43598
    ## Proportion of Variance  0.875 0.1122 0.00351 0.0026 0.00147 0.00124 0.00101
    ## Cumulative Proportion   0.875 0.9872 0.99070 0.9933 0.99477 0.99600 0.99701
    ##                            PC8     PC9    PC10   PC11    PC12   PC13    PC14
    ## Standard deviation     0.41054 0.29120 0.28441 0.2380 0.22723 0.1941 0.15327
    ## Proportion of Variance 0.00089 0.00045 0.00043 0.0003 0.00027 0.0002 0.00012
    ## Cumulative Proportion  0.99790 0.99835 0.99878 0.9991 0.99935 0.9996 0.99968
    ##                           PC15    PC16    PC17    PC18    PC19    PC20
    ## Standard deviation     0.14545 0.11274 0.10065 0.07723 0.07598 0.07337
    ## Proportion of Variance 0.00011 0.00007 0.00005 0.00003 0.00003 0.00003
    ## Cumulative Proportion  0.99979 0.99986 0.99991 0.99994 0.99997 1.00000
    ##                             PC21      PC22      PC23      PC24      PC25
    ## Standard deviation     9.744e-16 8.971e-16 8.971e-16 8.971e-16 8.971e-16
    ## Proportion of Variance 0.000e+00 0.000e+00 0.000e+00 0.000e+00 0.000e+00
    ## Cumulative Proportion  1.000e+00 1.000e+00 1.000e+00 1.000e+00 1.000e+00
    ##                             PC26      PC27      PC28      PC29
    ## Standard deviation     8.971e-16 8.971e-16 8.971e-16 8.112e-16
    ## Proportion of Variance 0.000e+00 0.000e+00 0.000e+00 0.000e+00
    ## Cumulative Proportion  1.000e+00 1.000e+00 1.000e+00 1.000e+00

``` r
# Representamos la varianza explicada por cada componente
plot(pca_resultado)
```

<figure>
<img
src="../images/unnamed-chunk-42-1.png"
style="width:100.0%" height="500"/>
</figure>

``` r
cargos <- pca_resultado$rotation
#print(cargos)
```

**PCA Results**

The first component (PC1) explains 87.5% of the total variance of the dataset, while the second component (PC2) adds an additional 11.2%. Together, PC1 and PC2 explain approximately 98.7% of the total variance, indicating that the main information in the data can be effectively summarized in two dimensions.

Loadings analysis revealed that:

- PC1 is strongly influenced by age at diagnosis, suggesting that this variable is the main driver of variability among patients.

- PC2 is dominated by the number of positive lymph nodes, followed by a smaller contribution from variables related to tumor invasion (lymphatic, venous, and perineural).

``` r
# Realizamos agrupamiento jerárquico aglomerativo
dist_matrix <- dist(df_numeric)
# Aplicamos agrupamiento jerárquico aglomerativo con el método de Ward
hc_aglomerativo <- hclust(dist_matrix, method = "ward.D2") 
# Reprsentamos el dendrograma
plot(hc_aglomerativo, main = "Dendrograma de Agrupamiento Jerárquico
Aglomerativo", xlab = "Pacientes")
```

<figure>
<img
src="../images/unnamed-chunk-44-1.png"
style="width:100.0%" height="500"/>
</figure>

---

## 5 Visualization with Shiny

The Shiny application was developed to enable interactive analysis of clinical data from gastrointestinal cancer patients, using a CSV file uploaded by the user and an exploratory and visual approach based on the `.Rmd file`. The `gastro_app.R` file includes the application source code (attached in the PEC4 submission), and some screenshots are shown below to illustrate its operation.

Overall, the application incorporates key dataset transformations, including calculating survival time (`os_time`), death event (`os_event`), and BMI, from variables such as weight and height. In addition, new categorical variables are generated for age (`age_category`) and BMI (`imc_category`), which facilitate comparative analysis among different patient groups.

<figure>
<img
src="../images/PEC4_Shiny1.png"
style="width:75.0%" height="300"
alt="Vista general de la aplicación después de cargar el archivo de datos CSV" />
<figcaption aria-hidden="true">General view of the application after loading the CSV data file</figcaption>
</figure>

The application contains a sidebar with selectors that allow the user to choose variables for three main sections:

- **Survival by Clinical Variables**: Shows a boxplot representing the relationship between a categorical variable (such as cancer type or tumor stage) and a quantitative variable (for example, survival time or age).

<figure>
<img
src="../images/PEC4_Shiny3.png"
style="width:75.0%" height="300"
alt="Supervivencia según variables clínicas" />
<figcaption aria-hidden="true">Survival by Clinical Variables</figcaption>
</figure>

- **Age and BMI Distribution**: Uses a bar chart to analyze frequency distribution of variables such as age and BMI, stratified by gender.  

<figure>
<img
src="../images/PEC4_Shiny4.png"
style="width:75.0%" height="300" alt="Distribución de Edades e IMC" />
<figcaption aria-hidden="true">Age and BMI Distribution</figcaption>
</figure>

- **Kaplan-Meier Survival Analysis**: Using the survival library, a Kaplan-Meier model is fitted to assess the probability of survival over time, based on selected variables such as gender or age group.  

<figure>
<img
src="../images/PEC4_Shiny5.png"
style="width:75.0%" height="300"
alt="Análisis de Supervivencia Kaplan-Meier" />
<figcaption aria-hidden="true">Kaplan-Meier Survival Analysis</figcaption>
</figure>

---

## 6 Conclusions

The project’s objective was to analyze clinical data from gastrointestinal cancer patients obtained from the public TCGA repository (The Cancer Genome Atlas), to study possible clinical factors influencing patient survival. From the dataset of more than 1,300 cases: 1) We selected relevant variables, cleaned and transformed the data, 2) Created exploratory plots and performed statistical analyses, 3) Implemented Kaplan-Meier survival analysis, 4) Trained both a supervised model (SVM) to predict survival and an unsupervised model (PCA and hierarchical clustering) to explore natural groupings among patients, and 5) Implemented a Shiny application to facilitate interactive visualization of these data.

Below we describe the main findings of the TCGA dataset analysis. We also discuss limitations and areas of opportunity, and close with a comment on the team experience developing this project.

**Most Interesting Findings of Our Analysis**

- Results from parametric `ANOVA` and `TukeyHSD` tests, as well as from `non-parametric Kruskal-Wallis` and `Dunn tests`, confirmed that there are statistically significant differences between gastrointestinal cancer types and survival days for the following groups: Pancreas-Colon, Pancreas-Rectum, Stomach-Colon, and Stomach-Rectum.

- Kaplan-Meier plots showed that women have greater survival than men, with about 65% surviving beyond 1,800 days, compared to 35% of men at 3,000 days. The log-rank test (p = 0.04) confirms that these differences are statistically significant, suggesting a more favorable clinical course for women.

- Regarding age, younger patients have higher and more sustained survival curves, while those over 70 and 80 show faster decline. The log-rank test (p = 0.0005) shows evidence of significant differences between age groups. This is consistent with the fact that older age can be associated with higher mortality and lower treatment tolerance.

- Finally, conditional probability analysis showed that 55.22% of patients with perineural invasion had an advanced pathological stage (type III or higher) at diagnosis. This suggests a strong association between perineural invasion and greater tumor progression.

**Limitations of the Data and Possible Improvements or Extensions**

One of the main limitations of the analysis was the high proportion of missing data, especially in categorical variables such as `venous_invasion`, `lymphatic_invasion`, and `perineural_invasion_present`, as well as some continuous variables such as `days_to_death`, `weight`, and `height`. This reduced the effective sample size. Moreover, since this is a retrospective observational study, i.e., data come from an existing repository rather than a controlled experimental design, associations could be identified but not causal relationships. As an improvement, we could apply imputation techniques to handle missing data and validate the results with external cohorts. Additionally, some studies of this type also incorporate genomic data to enrich the analysis and provide a more comprehensive view of survival factors.

**On Teamwork and Project Evaluation**

Salomon’s Comment: Teamwork with Sefora was key to selecting the dataset, thanks to her experience in cancer research. This helped us define key questions about how certain clinical variables could affect survival in gastrointestinal cancer patients. I really enjoyed this project because it allowed me to put theoretical concepts into practice with a real dataset. I think the collaboration was effective as we divided tasks, maintained good communication, and reviewed each other’s work.

Sefora’s Comment: The project allowed us to apply multiple data analysis tools in R and deepen our study of survival in gastrointestinal cancer, which is the focus of my research at Vall d’Hebron Institute of Oncology. I hope to apply the same tools to analyze data from patients at my center and advance understanding of these devastating diseases. It was a pleasure working with Salomon, as we have complementary approaches and perspectives that allowed us to mutually improve our work. We had good communication and divided tasks equally.

## 7 Bibliography

1. M. F. Bijlsma, A. Sadanandam, P. Tan, L. Vermeulen, Molecular subtypes in cancers of the gastrointestinal tract. Nat Rev Gastroenterol Hepatol 14, 333–342 (2017).

2. R. L. Siegel, K. D. Miller, H. E. Fuchs, A. Jemal, Cancer statistics, 2022. CA: A Cancer Journal for Clinicians 72, 7–33 (2022).

3. J. Drost, H. Clevers, Organoids in cancer research. Nat Rev Cancer 18, 407–418 (2018).

4. P. W. Nagle, J. Th. M. Plukker, C. T. Muijs, P. van Luijk, R. P. Coppes, Patient-derived tumor organoids for prediction of cancer treatment response. Seminars in Cancer Biology 53, 258–264 (2018).

5. C. A. Pasch, P. F. Favreau, A. E. Yueh, C. P. Babiarz, A. A. Gillette, J. T. Sharick, M. R. Karim, K. P. Nickel, A. K. DeZeeuw, C. M. Sprackling, P. B. Emmerich, R. A. DeStefanis, R. T. Pitera, S. N. Payne, D. P. Korkos, L. Clipson, C. M. Walsh, D. Miller, E. H. Carchman, M. E. Burkard, K. K. Lemmon, K. A. Matkowskyj, M. A. Newton, I. M. Ong, M. F. Bassetti, R. J. Kimple, M. C. Skala, D. A. Deming, Patient-Derived Cancer Organoid Cultures to Predict Sensitivity to Chemotherapy and Radiation. Clinical Cancer Research 25, 5376–5387 (2019).

6. E. R. Zanella, E. Grassi, L. Trusolino, Towards precision oncology with patient-derived xenografts. Nat Rev Clin Oncol 19, 719–732 (2022).

7. X. Sun, B. Hu, Mathematical modeling and computational prediction of cancer drug resistance. Briefings in Bioinformatics 19, 1382–1399 (2018).

8. E. J. Mucaki, J. Z. L. Zhao, D. J. Lizotte, P. K. Rogan, Predicting responses to platin chemotherapy agents with biochemically-inspired machine learning. Sig Transduct Target Ther 4, 1–12 (2019).

9. C. Yang, X. Huang, Y. Li, J. Chen, Y. Lv, S. Dai, Prognosis and personalized treatment prediction in TP53-mutant hepatocellular carcinoma: an in silico strategy towards precision oncology. Briefings in Bioinformatics 22, bbaa164 (2021).

10. F. Li, Y. Yang, Y. Wei, P. He, J. Chen, Z. Zheng, H. Bu, Deep learning-based predictive biomarker of pathological complete response to neoadjuvant chemotherapy from histological images in breast cancer. J Transl Med 19, 348 (2021). 11. Z. Liu, Z. Li, J. Qu, R. Zhang, X. Zhou, L. Li, K. Sun, Z. Tang, H. Jiang, H. Li, Q. Xiong, Y. Ding, X. Zhao, K. Wang, Z. Liu, J. Tian, Radiomics of Multiparametric MRI for Pretreatment Prediction of Pathologic Complete Response to Neoadjuvant Chemotherapy in Breast Cancer: A Multicenter Study. Clinical Cancer Research 25, 3538–3547 (2019).