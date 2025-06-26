gtex_tissues = [
    'Adipose_Subcutaneous',
    'Adipose_Visceral_Omentum',
    'Adrenal_Gland',
    'Artery_Aorta',
    'Artery_Coronary',
    'Artery_Tibial',
    'Brain_Amygdala',
    'Brain_Anterior_cingulate_cortex_BA24',
    'Brain_Caudate_basal_ganglia',
    'Brain_Cerebellar_Hemisphere',
    'Brain_Cerebellum',
    'Brain_Cortex',
    'Brain_Frontal_Cortex_BA9',
    'Brain_Hippocampus',
    'Brain_Hypothalamus',
    'Brain_Nucleus_accumbens_basal_ganglia',
    'Brain_Putamen_basal_ganglia',
    'Brain_Spinal_cord_cervical_c_1',
    'Brain_Substantia_nigra',
    'Breast_Mammary_Tissue',
    'Cells_Cultured_fibroblasts',
    'Cells_EBV_transformed_lymphocytes',
    'Colon_Sigmoid',
    'Colon_Transverse',
    'Esophagus_Gastroesophageal_Junction',
    'Esophagus_Mucosa',
    'Esophagus_Muscularis',
    'Heart_Atrial_Appendage',
    'Heart_Left_Ventricle',
    'Kidney_Cortex',
    'Liver',
    'Lung',
    'Minor_Salivary_Gland',
    'Muscle_Skeletal',
    'Nerve_Tibial',
    'Ovary',
    'Pancreas',
    'Pituitary',
    'Prostate',
    'Skin_Not_Sun_Exposed_Suprapubic',
    'Skin_Sun_Exposed_Lower_leg',
    'Small_Intestine_Terminal_Ileum',
    'Spleen',
    'Stomach',
    'Testis',
    'Thyroid',
    'Uterus',
    'Vagina',
    'Whole_Blood',
]

brain_tissues = [x for x in gtex_tissues if 'Brain' in x and 'Cerebell' not in x]
cerebellum_tissues = [x for x in gtex_tissues if 'Cerebell' in x]
heart_tissues = [x for x in gtex_tissues if 'Heart' in x]
kidney_tissues = [x for x in gtex_tissues if 'Kidney' in x]
liver_tissues = [x for x in gtex_tissues if 'Liver' in x]
ovary_tissues = [x for x in gtex_tissues if 'Ovary' in x]
testis_tissues = [x for x in gtex_tissues if 'Testis' in x]


kaesmann_gtex_map = {
    'brain': brain_tissues,
    'cerebellum': cerebellum_tissues,
    'heart': heart_tissues,
    'kidney': kidney_tissues,
    'liver': liver_tissues,
    'ovary': ovary_tissues,
    'testis': testis_tissues,
    'other_gtex_tissues': [x for x in gtex_tissues if x not in sorted(set([
        *brain_tissues, 
        *cerebellum_tissues,
        *heart_tissues,
        *kidney_tissues,
        *liver_tissues,
        *ovary_tissues,
        *testis_tissues,
        ]))]
}