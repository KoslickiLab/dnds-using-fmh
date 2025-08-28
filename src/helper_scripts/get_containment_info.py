from fmh_omega import helperfuncs

wd='/data/jzr5814/sourmash_dnds_estimation/tests/test/create_sequence_using_NG_assumption/0.1/fmh_dnds_sketch_protein'

positive_protein_df = helperfuncs.extract_containment_matrix(f'{wd}/positive_selection_redo_sketch_protein_using_faa/compare.protein.7.csv')
positive_protein_df = positive_protein_df.rename(columns={'containment':'AA_Cfrac'})
positive_protein_df.to_csv(f'{wd}/positive_selection_redo_sketch_protein_using_faa/containment.protein.7.csv')
negative_protein_df = helperfuncs.extract_containment_matrix(f'{wd}/negative_selection_redo_sketch_protein_using_faa/compare.protein.7.csv')
negative_protein_df = negative_protein_df.rename(columns={'containment':'AA_Cfrac'})
positive_protein_df.to_csv(f'{wd}/negative_selection_redo_sketch_protein_using_faa/containment.protein.7.csv')

