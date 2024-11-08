RBASE_TAG = '1.0'
PYBASE_TAG = '0.3.2'
MCMCGLMM_TAG = '1.0'
EMAIL = '' # Add your email address to download sequences via Entrez

rule all:
    input:
        'data/species_list_1st_select.tsv',
        'data/species_list_2nd_select.tsv',
        'data/genome/OG_summary_by_polar_inc_polar_1st_select_summary.tsv',
        'data/genome/OG_summary_by_polar_inc_polar_2nd_select_summary.tsv',
        'data/figure/n_gene_latitude/',
        'data/figure/n_gene_latitude_summary.pdf',
        # 'data/figure/n_gene_latitude_summary_renamed.pdf',
        'data/figure/tree_inc_polar_OG_2nd_select_all_species.pdf',
        'data/figure/tree_inc_polar_OG_1st_select_all_species.pdf',
        'data/figure/tree_inc_polar_OG_1st_select_downsampling.pdf',
        'data/species_list_1st_select_with_selected_OG_sort_by_family.tsv',
        'data/figure/tree_inc_polar_OG_2nd_select_pMCMC_all_species.pdf',
        # 'data/figure/tree_inc_polar_OG_2nd_select_all_species_renamed.pdf',
        # 'data/figure/tree_inc_polar_OG_2nd_select_pMCMC_all_species_renamed.pdf',
        'data/genome/OG_summary_by_polar_inc_polar_2nd_select_summary_GO_molecular_function.tsv',
        'data/genome/OG_summary_by_polar_inc_polar_2nd_select_summary_GO_biological_process.tsv',
        'data/genome/OG_summary_by_polar_inc_polar_2nd_select_summary_GO_cellular_component.tsv',
        'data/genome/OG_summary_by_polar_inc_polar_2nd_select_summary_pMCMC_filtered_GO_molecular_function.tsv',
        'data/genome/OG_summary_by_polar_inc_polar_2nd_select_summary_pMCMC_filtered_GO_biological_process.tsv',
        'data/genome/OG_summary_by_polar_inc_polar_2nd_select_summary_pMCMC_filtered_GO_cellular_component.tsv',
        'data/genome/191156at7898_Dm/gene_tree.pdf',
        'data/genome/428128at7898_Dm/gene_tree.pdf',
        'data/genome/409056at7898_Dm/gene_tree.pdf',
        'data/genome/489137at7898/gene_tree.pdf',
        'data/genome/428128at7898/gene_tree.pdf',
        'data/genome/299357at7898/gene_tree.pdf',
        'data/genome/465589at7898/gene_tree.pdf',
        'data/genome/409056at7898/gene_tree.pdf',
        'data/genome/172968at7898/gene_tree.pdf',
        'data/genome/29671at7898/gene_tree.pdf',
        'data/genome/191156at7898/gene_tree.pdf',
        'data/genome/498742at7898/gene_tree.pdf',
        'data/genome/454335at7898/gene_tree.pdf',
        'data/genome/488885at7898/gene_tree.pdf',
        'data/genome/428606at7898/gene_tree.pdf',
        'data/genome/469863at7898/gene_tree.pdf',
        'data/genome/485260at7898/gene_tree.pdf',
        'data/genome/482472at7898/gene_tree_labelright.pdf',
        'data/genome/438657at7898/gene_tree.pdf',
        'data/genome/500231at7898/gene_tree.pdf',
        'data/genome/307673at7898/gene_tree_labelright.pdf',
        'data/genome/493746at7898/gene_tree.pdf',
        'data/genome/482447at7898/gene_tree.pdf'

# Selection and classification of species into polar and non-polar fish
rule download_orthodb:
    output:
        'data/raw_data/orthodb/odb11v0_all_fasta.tab',
        # 'data/raw_data/orthodb/odb11v0_all_og_fasta.tab'
        'data/raw_data/orthodb/odb11v0_levels.tab',
        'data/raw_data/orthodb/odb11v0_species.tab',
        'data/raw_data/orthodb/odb11v0_level2species.tab',
        'data/raw_data/orthodb/odb11v0_genes.tab',
        'data/raw_data/orthodb/odb11v0_gene_xrefs.tab',
        'data/raw_data/orthodb/odb11v0_OGs.tab',
        'data/raw_data/orthodb/odb11v0_OG2genes.tab',
        'data/raw_data/orthodb/odb11v0_OG_xrefs.tab',
        'data/raw_data/orthodb/odb11v0_OG_pairs.tab',
        'data/raw_data/orthodb/README.txt'
    log:
        stdout = 'log/download_orthodb.stdout',
        stderr = 'log/download_orthodb.stderr'
    container:
        f'docker://aurelia01/dup_fish_rbase:{RBASE_TAG}'
    shell:
        'bash scripts/download_orthodb.sh > {log.stdout} 2> {log.stderr}'

rule download_fishtree:
    output:
        'data/raw_data/tree/fishtree.nwk'
    log:
        stdout = 'log/download_fishtree.stdout',
        stderr = 'log/download_fishtree.stderr'
    container:
        f'docker://aurelia01/dup_fish_rbase:{RBASE_TAG}'
    shell:
        'Rscript scripts/download_fishtree.R > {log.stdout} 2> {log.stderr}'

rule list_actinopterygii_sp_from_orthodb:
    input:
        levels = 'data/raw_data/orthodb/odb11v0_levels.tab',
        species = 'data/raw_data/orthodb/odb11v0_species.tab',
        level2species = 'data/raw_data/orthodb/odb11v0_level2species.tab',
        tree = 'data/raw_data/tree/fishtree.nwk'
    output:
        'data/species_list_tmp1.tsv'
    log:
        stdout = 'log/list_actinopterygii_sp_from_orthodb.stdout',
        stderr = 'log/list_actinopterygii_sp_from_orthodb.stderr'
    container:
        f'docker://aurelia01/dup_fish_rbase:{RBASE_TAG}'
    shell:
        'Rscript scripts/list_actinopterygii_sp_from_orthodb.R --levels {input.levels} --species {input.species} --level2species {input.level2species} --tree {input.tree} --output {output} > {log.stdout} 2> {log.stderr}'

rule select_proteomes_from_orthodb:
    input:
        sp_list = 'data/species_list_tmp1.tsv',
        fasta = 'data/raw_data/orthodb/odb11v0_all_fasta.tab'
    output:
        directory('data/genome/proteome/')
    log:
        stdout = 'log/select_proteomes_from_orthodb.stdout',
        stderr = 'log/select_proteomes_from_orthodb.stderr'
    container:
        f'docker://aurelia01/dup_fish_pybase:{PYBASE_TAG}'
    shell:
        'python scripts/select_proteomes_from_orthodb.py --species {input.sp_list} --fasta {input.fasta} --output {output} > {log.stdout} 2> {log.stderr}'

rule run_busco:
    input:
        'data/genome/proteome/'
    output:
        'data/genome/busco/batch_summary.txt'
    log:
        stdout = 'log/run_busco.stdout',
        stderr = 'log/run_busco.stderr'
    container:
        'docker://ezlabgva/busco:v5.6.1_cv1'
    threads:
        8
    shell:
        'busco -i {input} -o {output} -m proteins -l actinopterygii_odb10 -c {threads} > {log.stdout} 2> {log.stderr}'

rule join_busco_results:
    input:
        sp_list = 'data/species_list_tmp1.tsv',
        busco = 'data/genome/busco/batch_summary.txt'
    output:
        'data/species_list_tmp2.tsv'
    log:
        stdout = 'log/join_busco_results.stdout',
        stderr = 'log/join_busco_results.stderr'
    container:
        f'docker://aurelia01/dup_fish_rbase:{RBASE_TAG}'
    shell:
        'Rscript scripts/join_busco_results.R --sp_list {input.sp_list} --busco {input.busco} --output {output} > {log.stdout} 2> {log.stderr}'

rule join_polar_fish:
    input:
        sp_list = 'data/species_list_tmp2.tsv',
        polar_fish = 'polar_fish.csv'
    output:
        'data/species_list_tmp3.tsv'
    log:
        stdout = 'log/join_polar_fish.stdout',
        stderr = 'log/join_polar_fish.stderr'
    container:
        f'docker://aurelia01/dup_fish_rbase:{RBASE_TAG}'
    shell:
        'Rscript scripts/join_polar_fish.R --sp_list {input.sp_list} --polar_fish {input.polar_fish} --output {output} > {log.stdout} 2> {log.stderr}'

rule join_latitude:
    input:
        sp_list = 'data/species_list_tmp3.tsv',
        latitude = 'fb_latitude.tsv'
    output:
        'data/species_list_tmp4.tsv'
    log:
        stdout = 'log/join_latitude.stdout',
        stderr = 'log/join_latitude.stderr'
    container:
        f'docker://aurelia01/dup_fish_rbase:{RBASE_TAG}'
    shell:
        'Rscript scripts/join_latitude.R --sp_list {input.sp_list} --latitude {input.latitude} --output {output} > {log.stdout} 2> {log.stderr}'

rule mark_species_for_1st_filtering:
    input:
        'data/species_list_tmp4.tsv'
    output:
        'data/species_list_tmp5.tsv'
    log:
        stdout = 'log/mark_species_for_1st_filtering.stdout',
        stderr = 'log/mark_species_for_1st_filtering.stderr'
    container:
        f'docker://aurelia01/dup_fish_rbase:{RBASE_TAG}'
    shell:
        'Rscript scripts/mark_species_for_1st_filtering.R --input {input} --output {output} > {log.stdout} 2> {log.stderr}'

rule mark_species_for_2nd_filtering:
    input:
        'data/species_list_tmp5.tsv'
    output:
        'data/species_list_full.tsv'
    log:
        stdout = 'log/mark_species_for_2nd_filtering.stdout',
        stderr = 'log/mark_species_for_2nd_filtering.stderr'
    container:
        f'docker://aurelia01/dup_fish_rbase:{RBASE_TAG}'
    shell:
        'Rscript scripts/mark_species_for_2nd_filtering.R --input {input} --output {output} > {log.stdout} 2> {log.stderr}'

rule select_1st_selected_species:
    input:
        'data/species_list_full.tsv'
    output:
        'data/species_list_1st_select.tsv'
    log:
        stdout = 'log/select_1st_selected_species.stdout',
        stderr = 'log/select_1st_selected_species.stderr'
    container:
        f'docker://aurelia01/dup_fish_rbase:{RBASE_TAG}'
    shell:
        'Rscript scripts/select_1st_selected_species.R --input {input} --output {output} > {log.stdout} 2> {log.stderr}'

rule select_2nd_selected_species:
    input:
        'data/species_list_full.tsv'
    output:
        'data/species_list_2nd_select.tsv'
    log:
        stdout = 'log/select_2nd_selected_species.stdout',
        stderr = 'log/select_2nd_selected_species.stderr'
    container:
        f'docker://aurelia01/dup_fish_rbase:{RBASE_TAG}'
    shell:
        'Rscript scripts/select_2nd_selected_species.R --input {input} --output {output} > {log.stdout} 2> {log.stderr}'

# # Comprehensive detection of OGs with higher copy numbers in polar fish
rule extract_gene_from_orthodb:
    input:
        sp_list = 'data/species_list_full.tsv',
        genes = 'data/raw_data/orthodb/odb11v0_genes.tab',
        OGs = 'data/raw_data/orthodb/odb11v0_OGs.tab',
        OG2genes = 'data/raw_data/orthodb/odb11v0_OG2genes.tab'
    output:
        'data/genome/gene_list_full.tsv'
    log:
        stdout = 'log/extract_gene_from_orthodb.stdout',
        stderr = 'log/extract_gene_from_orthodb.stderr'
    container:
        f'docker://aurelia01/dup_fish_rbase:{RBASE_TAG}'
    shell:
        'Rscript scripts/extract_gene_from_orthodb.R --sp_list {input.sp_list} --genes {input.genes} --OGs {input.OGs} --OG2genes {input.OG2genes} --output {output} > {log.stdout} 2> {log.stderr}'

rule count_gene_by_OG:
    input:
        'data/genome/gene_list_full.tsv'
    output:
        'data/genome/n_gene_by_OG_full.tsv'
    log:
        stdout = 'log/count_gene_by_OG.stdout',
        stderr = 'log/count_gene_by_OG.stderr'
    container:
        f'docker://aurelia01/dup_fish_rbase:{RBASE_TAG}'
    shell:
        'Rscript scripts/count_gene_by_OG.R --input {input} --output {output} > {log.stdout} 2> {log.stderr}'

rule summarize_OG_by_polar:
    input:
        sp_list = 'data/species_list_{n_step}_select.tsv',
        n_gene_by_OG = 'data/genome/n_gene_by_OG_full.tsv'
    output:
        'data/genome/OG_summary_by_polar_{n_step}_select.tsv'
    log:
        stdout = 'log/summarize_OG_by_polar_{n_step}_select.stdout',
        stderr = 'log/summarize_OG_by_polar_{n_step}_select.stderr'
    container:
        f'docker://aurelia01/dup_fish_rbase:{RBASE_TAG}'
    shell:
        'Rscript scripts/summarize_OG_by_polar.R --sp_list {input.sp_list} --n_gene_by_OG {input.n_gene_by_OG} --output {output} > {log.stdout} 2> {log.stderr}'

rule summarize_OG_increased_in_polar_1st_select:
    input:
        select1 = 'data/genome/OG_summary_by_polar_1st_select.tsv',
        select2 = 'data/genome/OG_summary_by_polar_2nd_select.tsv'
    output:
        'data/genome/OG_summary_by_polar_inc_polar_1st_select_summary.tsv'
    log:
        stdout = 'log/summarize_OG_increased_in_polar_1st_select.stdout',
        stderr = 'log/summarize_OG_increased_in_polar_1st_select.stderr'
    container:
        f'docker://aurelia01/dup_fish_rbase:{RBASE_TAG}'
    shell:
        'Rscript scripts/summarize_OG_increased_in_polar_1st_select.R --select1 {input.select1} --select2 {input.select2} --output {output} > {log.stdout} 2> {log.stderr}'

rule summarize_OG_increased_in_polar_2nd_select:
    input:
        select1 = 'data/genome/OG_summary_by_polar_1st_select.tsv',
        select2 = 'data/genome/OG_summary_by_polar_2nd_select.tsv'
    output:
        'data/genome/OG_summary_by_polar_inc_polar_2nd_select_summary.tsv'
    log:
        stdout = 'log/summarize_OG_increased_in_polar_2nd_select.stdout',
        stderr = 'log/summarize_OG_increased_in_polar_2nd_select.stderr'
    container:
        f'docker://aurelia01/dup_fish_rbase:{RBASE_TAG}'
    shell:
        'Rscript scripts/summarize_OG_increased_in_polar_2nd_select.R --select1 {input.select1} --select2 {input.select2} --output {output} > {log.stdout} 2> {log.stderr}'

rule prune_tree:
    input:
        sp_list = 'data/species_list_{n_step}_select.tsv',
        tree = 'data/raw_data/tree/fishtree.nwk'
    output:
        'data/tree/fishtree_{n_step}_select.nwk'
    log:
        stdout = 'log/prune_tree_{n_step}_select.stdout',
        stderr = 'log/prune_tree_{n_step}_select.stderr'
    container:
        f'docker://aurelia01/dup_fish_rbase:{RBASE_TAG}'
    shell:
        'Rscript scripts/prune_tree.R --sp_list {input.sp_list} --tree {input.tree} --output {output} > {log.stdout} 2> {log.stderr}'

rule visualize_tree_OG_1st_select_all_species:
    input:
        sp_list = 'data/species_list_1st_select.tsv',
        gene = 'data/genome/gene_list_full.tsv',
        selected_OG = 'data/genome/OG_summary_by_polar_inc_polar_1st_select_summary.tsv',
        tree = 'data/tree/fishtree_1st_select.nwk'
    output:
        'data/figure/tree_inc_polar_OG_1st_select_all_species.pdf'
    log:
        stdout = 'log/visualize_tree_OG_1st_select_all_species.stdout',
        stderr = 'log/visualize_tree_OG_1st_select_all_species.stderr'
    container:
        f'docker://aurelia01/dup_fish_rbase:{RBASE_TAG}'
    shell:
        'Rscript scripts/visualize_tree_OG_1st_select_all_species.R --sp_list {input.sp_list} --gene {input.gene} --selected_OG {input.selected_OG} --tree {input.tree} --output {output} > {log.stdout} 2> {log.stderr}'

rule visualize_tree_OG_2nd_select_all_species:
    input:
        sp_list = 'data/species_list_1st_select.tsv',
        gene = 'data/genome/gene_list_full.tsv',
        selected_OG = 'data/genome/OG_summary_by_polar_inc_polar_2nd_select_summary.tsv',
        tree = 'data/tree/fishtree_1st_select.nwk'
    output:
        'data/figure/tree_inc_polar_OG_2nd_select_all_species.pdf'
    log:
        stdout = 'log/visualize_tree_OG_2nd_select_all_species.stdout',
        stderr = 'log/visualize_tree_OG_2nd_select_all_species.stderr'
    container:
        f'docker://aurelia01/dup_fish_rbase:{RBASE_TAG}'
    shell:
        'Rscript scripts/visualize_tree_OG_2nd_select_all_species.R --sp_list {input.sp_list} --gene {input.gene} --selected_OG {input.selected_OG} --tree {input.tree} --output {output} > {log.stdout} 2> {log.stderr}'

rule visualize_tree_OG_1st_select_downsampling:
    input:
        sp_list = 'data/species_list_2nd_select.tsv',
        gene = 'data/genome/gene_list_full.tsv',
        selected_OG = 'data/genome/OG_summary_by_polar_inc_polar_1st_select_summary.tsv',
        tree = 'data/tree/fishtree_2nd_select.nwk'
    output:
        'data/figure/tree_inc_polar_OG_1st_select_downsampling.pdf'
    log:
        stdout = 'log/visualize_tree_OG_1st_select_downsampling.stdout',
        stderr = 'log/visualize_tree_OG_1st_select_downsampling.stderr'
    container:
        f'docker://aurelia01/dup_fish_rbase:{RBASE_TAG}'
    shell:
        'Rscript scripts/visualize_tree_OG_1st_select_downsampling.R --sp_list {input.sp_list} --gene {input.gene} --selected_OG {input.selected_OG} --tree {input.tree} --output {output} > {log.stdout} 2> {log.stderr}'

rule visualize_n_gene_latitude:
    input:
        sp_list = 'data/species_list_1st_select.tsv',
        gene = 'data/genome/gene_list_full.tsv',
        selected_OG = 'data/genome/OG_summary_by_polar_inc_polar_2nd_select_summary_pMCMC_filtered.tsv'
    output:
        directory('data/figure/n_gene_latitude/')
    log:
        stdout = 'log/visualize_n_gene_latitude.stdout',
        stderr = 'log/visualize_n_gene_latitude.stderr'
    container:
        f'docker://aurelia01/dup_fish_rbase:{RBASE_TAG}'
    shell:
        'Rscript scripts/visualize_n_gene_latitude.R --sp_list {input.sp_list} --gene {input.gene} --selected_OG {input.selected_OG} --output {output} > {log.stdout} 2> {log.stderr}'

rule visualize_n_gene_latitude_summary:
    input:
        sp_list = 'data/species_list_1st_select.tsv',
        gene = 'data/genome/gene_list_full.tsv',
        selected_OG = 'data/genome/OG_summary_by_polar_inc_polar_2nd_select_summary_pMCMC_filtered.tsv'
    output:
        'data/figure/n_gene_latitude_summary.pdf'
    log:
        stdout = 'log/visualize_n_gene_latitude_summary.stdout',
        stderr = 'log/visualize_n_gene_latitude_summary.stderr'
    container:
        f'docker://aurelia01/dup_fish_rbase:{RBASE_TAG}'
    shell:
        'Rscript scripts/visualize_n_gene_latitude_summary.R --sp_list {input.sp_list} --gene {input.gene} --selected_OG {input.selected_OG} --output {output} > {log.stdout} 2> {log.stderr}'

# rule visualize_n_gene_latitude_summary_renamed:
#     input:
#         sp_list = 'data/species_list_1st_select.tsv',
#         gene = 'data/genome/gene_list_full.tsv',
#         selected_OG = 'data/genome/OG_summary_by_polar_inc_polar_2nd_select_summary_pMCMC_filtered_renamed.csv'
#     output:
#         'data/figure/n_gene_latitude_summary_renamed.pdf'
#     log:
#         stdout = 'log/visualize_n_gene_latitude_summary_renamed.stdout',
#         stderr = 'log/visualize_n_gene_latitude_summary_renamed.stderr'
#     container:
#         f'docker://aurelia01/dup_fish_rbase:{RBASE_TAG}'
#     shell:
#         'Rscript scripts/visualize_n_gene_latitude_summary_renamed.R --sp_list {input.sp_list} --gene {input.gene} --selected_OG {input.selected_OG} --output {output} > {log.stdout} 2> {log.stderr}'

rule mcmcglmm_n_gene_polar:
    input:
        sp_list = 'data/species_list_1st_select.tsv',
        gene_list = 'data/genome/gene_list_full.tsv',
        tree = 'data/tree/fishtree_1st_select.nwk',
        target_OG = 'data/genome/OG_summary_by_polar_inc_polar_2nd_select_summary.tsv'
    output:
        'data/genome/mcmcglmm/mcmcglmm_summary.tsv'
    log:
        stdout = 'log/mcmcglmm_n_gene_polar.stdout',
        stderr = 'log/mcmcglmm_n_gene_polar.stderr'
    container:
        f'docker://aurelia01/dup_fish_mcmcglmm:{MCMCGLMM_TAG}'
    threads:
        8
    shell:
        'Rscript scripts/mcmcglmm_n_gene_polar.R --sp_list {input.sp_list} --gene_list {input.gene_list} --tree {input.tree} --target_OG {input.target_OG} --output {output} --threads {threads} > {log.stdout} 2> {log.stderr}'

rule pmcmc_correction:
    input:
        'data/genome/mcmcglmm/mcmcglmm_summary.tsv'
    output:
        'data/genome/mcmcglmm/mcmcglmm_summary_BHcorrection.tsv'
    log:
        stdout = 'log/pmcmc_correction.stdout',
        stderr = 'log/pmcmc_correction.stderr'
    container:
        f'docker://aurelia01/dup_fish_mcmcglmm:{MCMCGLMM_TAG}'
    shell:
        'Rscript scripts/pmcmc_correction.R --input {input} --output {output} > {log.stdout} 2> {log.stderr}'

rule filter_OG_by_pMCMC:
    input:
        'data/genome/mcmcglmm/mcmcglmm_summary_BHcorrection.tsv'
    output:
        'data/genome/OG_summary_by_polar_inc_polar_2nd_select_summary_pMCMC_filtered.tsv'
    log:
        stdout = 'log/filter_OG_by_pMCMC.stdout',
        stderr = 'log/filter_OG_by_pMCMC.stderr'
    container:
        f'docker://aurelia01/dup_fish_mcmcglmm:{MCMCGLMM_TAG}'
    shell:
        'Rscript scripts/filter_OG_by_pMCMC.R --input {input} --output {output} > {log.stdout} 2> {log.stderr}'

rule visualize_tree_OG_2nd_select_pMCMC_all_species:
    input:
        sp_list = 'data/species_list_1st_select.tsv',
        gene = 'data/genome/gene_list_full.tsv',
        selected_OG = 'data/genome/OG_summary_by_polar_inc_polar_2nd_select_summary_pMCMC_filtered.tsv',
        tree = 'data/tree/fishtree_1st_select.nwk'
    output:
        'data/figure/tree_inc_polar_OG_2nd_select_pMCMC_all_species.pdf'
    log:
        stdout = 'log/visualize_tree_OG_2nd_select_pMCMC_all_species.stdout',
        stderr = 'log/visualize_tree_OG_2nd_select_pMCMC_all_species.stderr'
    container:
        f'docker://aurelia01/dup_fish_rbase:{RBASE_TAG}'
    shell:
        'Rscript scripts/visualize_tree_OG_2nd_select_pMCMC_all_species.R --sp_list {input.sp_list} --gene {input.gene} --selected_OG {input.selected_OG} --tree {input.tree} --output {output} > {log.stdout} 2> {log.stderr}'

# Formatting
rule merge_species_list_and_selected_OG:
    input:
        sp_list = 'data/species_list_1st_select.tsv',
        gene_list = 'data/genome/gene_list_full.tsv',
        selected_OG = 'data/genome/OG_summary_by_polar_inc_polar_2nd_select_summary.tsv'
    output:
        'data/species_list_1st_select_with_selected_OG.tsv'
    log:
        stdout = 'log/merge_species_list_and_selected_OG.stdout',
        stderr = 'log/merge_species_list_and_selected_OG.stderr'
    container:
        f'docker://aurelia01/dup_fish_rbase:{RBASE_TAG}'
    shell:
        'Rscript scripts/merge_species_list_and_selected_OG.R --sp_list {input.sp_list} --gene_list {input.gene_list} --selected_OG {input.selected_OG} --output {output} > {log.stdout} 2> {log.stderr}'

rule sort_by_family_name:
    input:
        'data/species_list_1st_select_with_selected_OG.tsv'
    output:
        'data/species_list_1st_select_with_selected_OG_sort_by_family.tsv'
    log:
        stdout = 'log/sort_by_family_name.stdout',
        stderr = 'log/sort_by_family_name.stderr'
    container:
        f'docker://aurelia01/dup_fish_rbase:{RBASE_TAG}'
    shell:
        'Rscript scripts/sort_by_family_name.R --input {input} --output {output} > {log.stdout} 2> {log.stderr}'

# rule visualize_tree_OG_2nd_select_all_species_renamed:
#     input:
#         sp_list = 'data/species_list_1st_select.tsv',
#         gene = 'data/genome/gene_list_full.tsv',
#         selected_OG = 'data/genome/OG_summary_by_polar_inc_polar_2nd_select_summary_renamed.csv',
#         tree = 'data/tree/fishtree_1st_select.nwk'
#     output:
#         'data/figure/tree_inc_polar_OG_2nd_select_all_species_renamed.pdf'
#     log:
#         stdout = 'log/visualize_tree_OG_2nd_select_all_species_renamed.stdout',
#         stderr = 'log/visualize_tree_OG_2nd_select_all_species_renamed.stderr'
#     container:
#         f'docker://aurelia01/dup_fish_rbase:{RBASE_TAG}'
#     shell:
#         'Rscript scripts/visualize_tree_OG_2nd_select_all_species_renamed.R --sp_list {input.sp_list} --gene {input.gene} --selected_OG {input.selected_OG} --tree {input.tree} --output {output} > {log.stdout} 2> {log.stderr}'

# rule visualize_tree_OG_2nd_select_pMCMC_all_species_renamed:
#     input:
#         sp_list = 'data/species_list_1st_select.tsv',
#         gene = 'data/genome/gene_list_full.tsv',
#         selected_OG = 'data/genome/OG_summary_by_polar_inc_polar_2nd_select_summary_pMCMC_filtered_renamed.csv',
#         tree = 'data/tree/fishtree_1st_select.nwk'
#     output:
#         'data/figure/tree_inc_polar_OG_2nd_select_pMCMC_all_species_renamed.pdf'
#     log:
#         stdout = 'log/visualize_tree_OG_2nd_select_pMCMC_all_species_renamed.stdout',
#         stderr = 'log/visualize_tree_OG_2nd_select_pMCMC_all_species_renamed.stderr'
#     container:
#         f'docker://aurelia01/dup_fish_rbase:{RBASE_TAG}'
#     shell:
#         'Rscript scripts/visualize_tree_OG_2nd_select_pMCMC_all_species_renamed.R --sp_list {input.sp_list} --gene {input.gene} --selected_OG {input.selected_OG} --tree {input.tree} --output {output} > {log.stdout} 2> {log.stderr}'

# Gene ontology (GO) enrichment analysis
rule GO_enrichment_analysis_full:
    input:
        GO = 'data/raw_data/orthodb/odb11v0_OG_xrefs.tab',
        OG_list = 'data/genome/OG_summary_by_polar_2nd_select.tsv',
        selected_OG = 'data/genome/OG_summary_by_polar_inc_polar_2nd_select_summary.tsv'
    output:
        'data/genome/OG_summary_by_polar_inc_polar_2nd_select_summary_GO_{GO_category}.tsv'
    log:
        stdout = 'log/GO_enrichment_analysis_full_{GO_category}.stdout',
        stderr = 'log/GO_enrichment_analysis_full_{GO_category}.stderr'
    container:
        f'docker://aurelia01/dup_fish_rbase:{RBASE_TAG}'
    shell:
        'Rscript scripts/GO_enrichment_analysis.R --GO {input.GO} --OG_list {input.OG_list} --selected_OG {input.selected_OG} --GO_category {wildcards.GO_category} --output {output} > {log.stdout} 2> {log.stderr}'

rule GO_enrichment_analysis_selected_by_pMCMC:
    input:
        GO = 'data/raw_data/orthodb/odb11v0_OG_xrefs.tab',
        OG_list = 'data/genome/OG_summary_by_polar_2nd_select.tsv',
        selected_OG = 'data/genome/OG_summary_by_polar_inc_polar_2nd_select_summary_pMCMC_filtered.tsv'
    output:
        'data/genome/OG_summary_by_polar_inc_polar_2nd_select_summary_pMCMC_filtered_GO_{GO_category}.tsv'
    log:
        stdout = 'log/GO_enrichment_analysis_selected_by_pMCMC_{GO_category}.stdout',
        stderr = 'log/GO_enrichment_analysis_selected_by_pMCMC_{GO_category}.stderr'
    container:
        f'docker://aurelia01/dup_fish_rbase:{RBASE_TAG}'
    shell:
        'Rscript scripts/GO_enrichment_analysis.R --GO {input.GO} --OG_list {input.OG_list} --selected_OG {input.selected_OG} --GO_category {wildcards.GO_category} --output {output} > {log.stdout} 2> {log.stderr}'

# Analysis at the sequence level
rule get_protein_and_nucleotide_sequences_of_OG:
    input:
        sp_list = 'data/species_list_2nd_select.tsv',
        gene_list = 'data/genome/gene_list_full.tsv',
        gene_protid = 'data/raw_data/orthodb/odb11v0_genes.tab'
    output:
        protein = 'data/raw_data/ncbi/{OG_id}_AA.fasta',
        nucleotide = 'data/raw_data/ncbi/{OG_id}_NT.fasta'
    log:
        stdout = 'log/get_protein_and_nucleotide_sequences_of_{OG_id}.stdout',
        stderr = 'log/get_protein_and_nucleotide_sequences_of_{OG_id}.stderr'
    container:
        f'docker://aurelia01/dup_fish_pybase:{PYBASE_TAG}'
    shell:
        'python scripts/get_protein_and_nucleotide_sequences_of_OG.py --og_id {wildcards.OG_id} --sp_list {input.sp_list} --gene_list {input.gene_list} --gene_protid {input.gene_protid} --output_protein {output.protein} --output_nucleotide {output.nucleotide} --email {EMAIL} > {log.stdout} 2> {log.stderr}'

rule prepare_protein_sequences_of_DmZPC5:
    output:
        'data/raw_data/ncbi/DmZPC5_AA.fasta'
    log:
        stdout = 'log/prepare_protein_sequences_of_DmZPC5_of_DmZPC5.stdout',
        stderr = 'log/prepare_protein_sequences_of_DmZPC5_of_DmZPC5.stderr'
    container:
        f'docker://aurelia01/dup_fish_pybase:{PYBASE_TAG}'
    shell:
        'python scripts/prepare_protein_sequences_of_DmZPC5.py --output {output} > {log.stdout} 2> {log.stderr}'

rule get_protein_and_nucleotide_sequences_of_DmZPAX1:
    output:
        protein = 'data/raw_data/ncbi/DmZPAX1_AA.fasta',
        nucleotide = 'data/raw_data/ncbi/DmZPAX1_NT.fasta'
    log:
        stdout = 'log/get_protein_and_nucleotide_sequences_of_DmZPAX1.stdout',
        stderr = 'log/get_protein_and_nucleotide_sequences_of_DmZPAX1.stderr'
    container:
        f'docker://aurelia01/dup_fish_pybase:{PYBASE_TAG}'
    shell:
        'python scripts/get_protein_and_nucleotide_sequences_of_DmZP.py --protein_id AIO03056.1 --description DmZPAX1 --output_protein {output.protein} --output_nucleotide {output.nucleotide} --email {EMAIL} > {log.stdout} 2> {log.stderr}'

rule get_protein_and_nucleotide_sequences_of_DmZPC1:
    output:
        protein = 'data/raw_data/ncbi/DmZPC1_AA.fasta',
        nucleotide = 'data/raw_data/ncbi/DmZPC1_NT.fasta'
    log:
        stdout = 'log/get_protein_and_nucleotide_sequences_of_DmZPC1.stdout',
        stderr = 'log/get_protein_and_nucleotide_sequences_of_DmZPC1.stderr'
    container:
        f'docker://aurelia01/dup_fish_pybase:{PYBASE_TAG}'
    shell:
        'python scripts/get_protein_and_nucleotide_sequences_of_DmZP.py --protein_id AJW66345.1 --description DmZPC1 --output_protein {output.protein} --output_nucleotide {output.nucleotide} --email {EMAIL} > {log.stdout} 2> {log.stderr}'

rule cp_fasta:
    input:
        'data/raw_data/ncbi/{OG_id}at7898_AA.fasta'
    output:
        'data/genome/{OG_id}at7898/AA.fasta'
    log:
        'log/cp_fasta_{OG_id}at7898.stderr'
    container:
        f'docker://aurelia01/dup_fish_pybase:{PYBASE_TAG}'
    shell:
        'cp {input} {output} 2> {log}'

rule concat_fasta_191156at7898:
    input:
        fasta1 = 'data/raw_data/ncbi/DmZPC5_AA.fasta',
        fasta2 = 'data/raw_data/ncbi/191156at7898_AA.fasta'
    output:
        'data/genome/191156at7898_Dm/AA.fasta'
    log:
        'log/concat_fasta_191156at7898.stderr'
    container:
        f'docker://aurelia01/dup_fish_pybase:{PYBASE_TAG}'
    shell:
        'cat {input.fasta1} {input.fasta2} > {output} 2> {log}'

rule concat_fasta_409056at7898:
    input:
        fasta1 = 'data/raw_data/ncbi/DmZPC1_AA.fasta',
        fasta2 = 'data/raw_data/ncbi/409056at7898_AA.fasta'
    output:
        'data/genome/409056at7898_Dm/AA.fasta'
    log:
        'log/concat_fasta_409056at7898.stderr'
    container:
        f'docker://aurelia01/dup_fish_pybase:{PYBASE_TAG}'
    shell:
        'cat {input.fasta1} {input.fasta2} > {output} 2> {log}'

rule concat_fasta_428128at7898:
    input:
        fasta1 = 'data/raw_data/ncbi/DmZPAX1_AA.fasta',
        fasta2 = 'data/raw_data/ncbi/428128at7898_AA.fasta'
    output:
        'data/genome/428128at7898_Dm/AA.fasta'
    log:
        'log/concat_fasta_428128at7898.stderr'
    container:
        f'docker://aurelia01/dup_fish_pybase:{PYBASE_TAG}'
    shell:
        'cat {input.fasta1} {input.fasta2} > {output} 2> {log}'

rule run_mafft:
    input:
        'data/genome/{OG_id}at7898/AA.fasta'
    output:
        'data/genome/{OG_id}at7898/aln_AA.fasta'
    log:
        'log/run_mafft_{OG_id}at7898.stderr'
    container:
        f'docker://aurelia01/dup_fish_pybase:{PYBASE_TAG}'
    shell:
        'mafft-linsi --amino --anysymbol {input} > {output} 2> {log}'

rule run_mafft_Dm:
    input:
        'data/genome/{OG_id}at7898_Dm/AA.fasta'
    output:
        'data/genome/{OG_id}at7898_Dm/aln_AA.fasta'
    log:
        'log/run_mafft_{OG_id}at7898_Dm.stderr'
    container:
        f'docker://aurelia01/dup_fish_pybase:{PYBASE_TAG}'
    shell:
        'mafft-linsi --amino --anysymbol {input} > {output} 2> {log}'

rule run_trimal:
    input:
        'data/genome/{OG_id}at7898/aln_AA.fasta'
    output:
        'data/genome/{OG_id}at7898/trimmed_aln_AA.fasta'
    log:
        stdout = 'log/run_trimal_{OG_id}at7898.stdout',
        stderr = 'log/run_trimal_{OG_id}at7898.stderr'
    container:
        f'docker://aurelia01/dup_fish_pybase:{PYBASE_TAG}'
    shell:
        'trimal -in {input} -out {output} -automated1 > {log.stdout} 2> {log.stderr}'

rule run_trimal_Dm:
    input:
        'data/genome/{OG_id}at7898_Dm/aln_AA.fasta'
    output:
        'data/genome/{OG_id}at7898_Dm/trimmed_aln_AA.fasta'
    log:
        stdout = 'log/run_trimal_{OG_id}at7898_Dm.stdout',
        stderr = 'log/run_trimal_{OG_id}at7898_Dm.stderr'
    container:
        f'docker://aurelia01/dup_fish_pybase:{PYBASE_TAG}'
    shell:
        'trimal -in {input} -out {output} -automated1 > {log.stdout} 2> {log.stderr}'

rule run_iqtree:
    input:
        'data/genome/{OG_id}at7898/trimmed_aln_AA.fasta'
    output:
        'data/genome/{OG_id}at7898/trimmed_aln_AA.fasta.treefile'
    log:
        stdout = 'log/run_iqtree_{OG_id}at7898.stdout',
        stderr = 'log/run_iqtree_{OG_id}at7898.stderr'
    container:
        f'docker://aurelia01/dup_fish_pybase:{PYBASE_TAG}'
    threads:
        8
    shell:
        'iqtree -s {input} -m LG+G -B 1000 -T {threads} > {log.stdout} 2> {log.stderr}'

rule run_iqtree_Dm:
    input:
        'data/genome/{OG_id}at7898_Dm/trimmed_aln_AA.fasta'
    output:
        'data/genome/{OG_id}at7898_Dm/trimmed_aln_AA.fasta.treefile'
    log:
        stdout = 'log/run_iqtree_{OG_id}at7898_Dm.stdout',
        stderr = 'log/run_iqtree_{OG_id}at7898_Dm.stderr'
    container:
        f'docker://aurelia01/dup_fish_pybase:{PYBASE_TAG}'
    threads:
        8
    shell:
        'iqtree -s {input} -m LG+G -B 1000 -T {threads} > {log.stdout} 2> {log.stderr}'

rule binarize_branch_length:
    input:
        'data/tree/fishtree_2nd_select.nwk'
    output:
        'data/tree/fishtree_2nd_select_tmp.nwk'
    log:
        stdout = 'log/binarize_branch_length.stdout',
        stderr = 'log/binarize_branch_length.stderr'
    container:
        f'docker://aurelia01/dup_fish_rbase:{RBASE_TAG}'
    shell:
        'Rscript scripts/binarize_branch_length.R --t {input}  --o {output} > {log.stdout} 2> {log.stderr}'

rule add_Dissostichus_mawsoni:
    input:
        'data/species_list_2nd_select.tsv',
    output:
        'data/species_list_2nd_select_Dm.tsv'
    log:
        stdout = 'log/add_Dissostichus_mawsoni.stdout',
        stderr = 'log/add_Dissostichus_mawsoni.stderr'
    container:
        f'docker://aurelia01/dup_fish_rbase:{RBASE_TAG}'
    shell:
        'Rscript scripts/add_Dissostichus_mawsoni.R --input {input} --output {output} > {log.stdout} 2> {log.stderr}'

rule prune_tree_Dm:
    input:
        sp_list = 'data/species_list_2nd_select_Dm.tsv',
        tree = 'data/raw_data/tree/fishtree.nwk'
    output:
        'data/tree/fishtree_2nd_select_Dm.nwk'
    log:
        stdout = 'log/prune_tree_2nd_select_Dm.stdout',
        stderr = 'log/prune_tree_2nd_select_Dm.stderr'
    container:
        f'docker://aurelia01/dup_fish_rbase:{RBASE_TAG}'
    shell:
        'Rscript scripts/prune_tree.R --sp_list {input.sp_list} --tree {input.tree} --output {output} > {log.stdout} 2> {log.stderr}'

rule binarize_branch_length_Dm:
    input:
        'data/tree/fishtree_2nd_select_Dm.nwk'
    output:
        'data/tree/fishtree_2nd_select_Dm_tmp.nwk'
    log:
        stdout = 'log/binarize_branch_length_Dm.stdout',
        stderr = 'log/binarize_branch_length_Dm.stderr'
    container:
        f'docker://aurelia01/dup_fish_rbase:{RBASE_TAG}'
    shell:
        'Rscript scripts/binarize_branch_length.R --t {input}  --o {output} > {log.stdout} 2> {log.stderr}'

rule run_generax:
    input:
        species_tree = 'data/tree/fishtree_2nd_select_tmp.nwk',
        gene_tree = 'data/genome/{OG_id}at7898/trimmed_aln_AA.fasta.treefile',
        aln = 'data/genome/{OG_id}at7898/trimmed_aln_AA.fasta'
    output:
        'data/genome/{OG_id}at7898/generax/results/{OG_id}at7898/geneTree.newick'
    log:
        stdout = 'log/run_generax_{OG_id}at7898.stdout',
        stderr = 'log/run_generax_{OG_id}at7898.stderr'
    container:
        f'docker://aurelia01/dup_fish_pybase:{PYBASE_TAG}'
    shell:
        'python scripts/run_generax.py --species_tree {input.species_tree} --gene_tree {input.gene_tree} --aln {input.aln} --output {output} > {log.stdout} 2> {log.stderr}'

rule run_generax_Dm:
    input:
        species_tree = 'data/tree/fishtree_2nd_select_Dm_tmp.nwk',
        gene_tree = 'data/genome/{OG_id}at7898_Dm/trimmed_aln_AA.fasta.treefile',
        aln = 'data/genome/{OG_id}at7898_Dm/trimmed_aln_AA.fasta'
    output:
        'data/genome/{OG_id}at7898_Dm/generax/results/{OG_id}at7898_Dm/geneTree.newick'
    log:
        stdout = 'log/run_generax_{OG_id}at7898_Dm.stdout',
        stderr = 'log/run_generax_{OG_id}at7898_Dm.stderr'
    container:
        f'docker://aurelia01/dup_fish_pybase:{PYBASE_TAG}'
    shell:
        'python scripts/run_generax.py --species_tree {input.species_tree} --gene_tree {input.gene_tree} --aln {input.aln} --output {output} > {log.stdout} 2> {log.stderr}'

rule calculate_global_DE_ratio:
    input:
        'data/genome/{OG_id}at7898/AA.fasta'
    output:
        'data/genome/{OG_id}at7898/global_DE_ratio.tsv'
    log:
        stdout = 'log/calculate_global_DE_ratio_{OG_id}at7898.stdout',
        stderr = 'log/calculate_global_DE_ratio_{OG_id}at7898.stderr'
    container:
        f'docker://aurelia01/dup_fish_pybase:{PYBASE_TAG}'
    shell:
        'python scripts/calculate_DE_ratio.py --input {input} --output {output} > {log.stdout} 2> {log.stderr}'

rule calculate_global_DE_ratio_Dm:
    input:
        'data/genome/{OG_id}at7898_Dm/AA.fasta'
    output:
        'data/genome/{OG_id}at7898_Dm/global_DE_ratio.tsv'
    log:
        stdout = 'log/calculate_global_DE_ratio_{OG_id}at7898_Dm.stdout',
        stderr = 'log/calculate_global_DE_ratio_{OG_id}at7898_Dm.stderr'
    container:
        f'docker://aurelia01/dup_fish_pybase:{PYBASE_TAG}'
    shell:
        'python scripts/calculate_DE_ratio.py --input {input} --output {output} > {log.stdout} 2> {log.stderr}'

rule summarize_sequence_information:
    input:
        sp_list = 'data/species_list_2nd_select.tsv',
        DE_ratio = 'data/genome/{OG_id}at7898/global_DE_ratio.tsv'
    output:
        'data/genome/{OG_id}at7898/seqinfo.tsv'
    log:
        stdout = 'log/summarize_sequence_information_{OG_id}at7898.stdout',
        stderr = 'log/summarize_sequence_information_{OG_id}at7898.stderr'
    container:
        f'docker://aurelia01/dup_fish_rbase:{RBASE_TAG}'
    shell:
        'Rscript scripts/summarize_sequence_information.R --sp_list {input.sp_list} --DE_ratio {input.DE_ratio} --output {output} > {log.stdout} 2> {log.stderr}'

rule summarize_sequence_information_Dm:
    input:
        sp_list = 'data/species_list_2nd_select_Dm.tsv',
        DE_ratio = 'data/genome/{OG_id}at7898_Dm/global_DE_ratio.tsv'
    output:
        'data/genome/{OG_id}at7898_Dm/seqinfo.tsv'
    log:
        stdout = 'log/summarize_sequence_information_{OG_id}at7898_Dm.stdout',
        stderr = 'log/summarize_sequence_information_{OG_id}at7898_Dm.stderr'
    container:
        f'docker://aurelia01/dup_fish_rbase:{RBASE_TAG}'
    shell:
        'Rscript scripts/summarize_sequence_information.R --sp_list {input.sp_list} --DE_ratio {input.DE_ratio} --output {output} > {log.stdout} 2> {log.stderr}'

rule visualize_gene_tree:
    input:
        tree = 'data/genome/{OG_id}at7898/generax/results/{OG_id}at7898/geneTree.newick',
        seqinfo = 'data/genome/{OG_id}at7898/seqinfo.tsv'
    output:
        'data/genome/{OG_id}at7898/gene_tree.pdf'
    log:
        stdout = 'log/visualize_gene_tree_{OG_id}at7898.stdout',
        stderr = 'log/visualize_gene_tree_{OG_id}at7898.stderr'
    container:
        f'docker://aurelia01/dup_fish_rbase:{RBASE_TAG}'
    shell:
        'Rscript scripts/visualize_gene_tree.R --tree {input.tree} --seqinfo {input.seqinfo} --output {output} > {log.stdout} 2> {log.stderr}'

rule visualize_gene_tree_labelright:
    input:
        tree = 'data/genome/{OG_id}at7898/generax/results/{OG_id}at7898/geneTree.newick',
        seqinfo = 'data/genome/{OG_id}at7898/seqinfo.tsv'
    output:
        'data/genome/{OG_id}at7898/gene_tree_labelright.pdf'
    log:
        stdout = 'log/visualize_gene_tree_{OG_id}at7898_labelright.stdout',
        stderr = 'log/visualize_gene_tree_{OG_id}at7898_labelright.stderr'
    container:
        f'docker://aurelia01/dup_fish_rbase:{RBASE_TAG}'
    shell:
        'Rscript scripts/visualize_gene_tree_labelright.R --tree {input.tree} --seqinfo {input.seqinfo} --output {output} > {log.stdout} 2> {log.stderr}'

rule visualize_gene_tree_Dm:
    input:
        tree = 'data/genome/{OG_id}at7898_Dm/generax/results/{OG_id}at7898_Dm/geneTree.newick',
        seqinfo = 'data/genome/{OG_id}at7898_Dm/seqinfo.tsv'
    output:
        'data/genome/{OG_id}at7898_Dm/gene_tree.pdf'
    log:
        stdout = 'log/visualize_gene_tree_{OG_id}at7898_Dm.stdout',
        stderr = 'log/visualize_gene_tree_{OG_id}at7898_Dm.stderr'
    container:
        f'docker://aurelia01/dup_fish_rbase:{RBASE_TAG}'
    shell:
        'Rscript scripts/visualize_gene_tree.R --tree {input.tree} --seqinfo {input.seqinfo} --output {output} > {log.stdout} 2> {log.stderr}'
