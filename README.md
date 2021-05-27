# Tutorial

## Generating PSSM profiles

To generate PSSM profiles for protein sequences, the helper function *create_pssm_profile* can be used. However, before using the function, the steps mentioned below must be followed.

- Download a blast database: For eg. uniref50 database can be downloaded using this link 

    http://ftp.ebi.ac.uk/pub/databases/uniprot/uniref/uniref50/ 
    
    
- Download blast executables using this link preferably version 2.9.0. The psiblast program used to create the PSSM profiles will be downloaded as well along with makeblastdb program to be used in the next step

    https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/
    
    
- Make a local blast database using uniref50 fasta file and the blast executable (makeblastdb). The following command can be used for that purpose.

    \$makeblastdb -in uniref50.fasta -dbtype prot -out uniref50
    
For further information, refer here: 

https://quickgrid.blogspot.com/2018/10/Python-Sub-Process-Local-Psi-Blast-PSSM-Generation-from-FASTA-in-Directory-using-Uniref50-Database-in-Pycharm.html


Once the above steps have been followed and there is an indexed blast database on your local machine, the *create_pssm_profile* function can be used. It requires the following arguments:

1. A comma separated protein sequence file where each line contains the name of the protein followed by its sequence separated by a comma.

2. The output directory where the user would like to store the pssm profiles.

3. The path of the psiblast program executable downloaded in step 2 described above.

4. The path of the indexed blast database directory created in step 3 described above.
    


```python
# Usage example

from pssmpro.features import create_pssm_profile

# The comma separated protein sequence file
protein_sequence_file = "./pssmpro_test_data/test_seq.csv"
# Output directory where the pssm profiles will be stored
output_dir = "./pssmpro_test_data/pssm_profiles/"
# the path to the psiblast program executable downloaded as part of the blast program suite 
psiblast_executable_path = "/opt/aci/sw/ncbi-rmblastn/2.9/0_gcc-8.3.1-bxy/bin/psiblast"
# prefix of the indexed blast database files created using makeblastdb
blast_db_prefix = "./pssmpro_test_data/uniref50/uniref50db"
# number of cores to be used while creating the pssm profiles
number_of_cores = 8


create_pssm_profile(protein_sequence_file, output_dir, psiblast_executable_path,
                    blast_db_prefix, number_of_cores)
```

## Generating PSSM features

*pssmpro* contains 21 features which are capable of numerically encoding protein sequences using their pssm profiles. They are:

1. aac_pssm
2. aadp_pssm
3. aatp
4. ab_pssm
5. d_fpssm
6. dp_pssm
7. dpc_pssm
8. edp
9. eedp
10. k_separated_bigrams_pssm
11. medp
12. pse_pssm
13. pssm_ac
14. pssm_cc
15. pssm_composition
16. rpm_pssm
17. rpssm
18. s_fpssm
19. smoothed_pssm
20. tpc
21. tri_gram_pssm

**For a detailed description of the features, refer to the Supplementary Documents of the paper (link to be added)**


```python
# Usage example

# To create any one of the 21 features one can use the "get_feature" function

pssm_dir_path = "./pssmpro_test_data/pssm_profiles/"
feature_type = "aac_pssm"
output_dir_path = "./pssmpro_test_data/features/"

get_feature(pssm_dir_path, feature_type, output_dir_path)
```


```python
# To create all 21 features at once, one can use the "get_all_features" function

get_all_features(pssm_dir_path, output_dit_path)
```
