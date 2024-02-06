import pandas as pd
file='/data/jzr5814/sourmash_dnds_estimation/tests/test/create_sequence_using_NG_assumption/0.001/fmh_dnds/negative_selection/approximating_AAI_from_sourmash_translate.xlsx'
df = pd.read_excel(file).head()
print(df)
df = df[['AAI','Protein Cfrac (translate)','ksize','Approx AAI']]
print(df)
df.to_latex("/data/jzr5814/sourmash_dnds_estimation/tests/test/create_sequence_using_NG_assumption/0.001/fmh_dnds/negative_selection/approx_AAI_from_SM_translate_for_overleaf.tex")
