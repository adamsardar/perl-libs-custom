
=pod 

=head1 NAME

  Bio::EnsEMBL::Compara::PipeConfig::ProteinTrees_conf

=head1 SYNOPSIS

    #1. update ensembl-hive, ensembl and ensembl-compara CVS repositories before each new release

    #2. you may need to update 'schema_version' in meta table to the current release number in ensembl-hive/sql/tables.sql

    #3. make sure that all default_options are set correctly

    #4. Run init_pipeline.pl script:
        init_pipeline.pl Bio::EnsEMBL::Compara::PipeConfig::ProteinTrees_conf -password <your_password> -mlss_id <your_current_PT_mlss_id>

    #5. Sync and loop the beekeeper.pl as shown in init_pipeline.pl's output


=head1 DESCRIPTION  

    The PipeConfig file for ProteinTrees pipeline that should automate most of the pre-execution tasks.

=head2 rel.63 stats

    sequences to cluster:       1,198,678           [ SELECT count(*) from sequence; ]
    reused core dbs:            48                  [ SELECT count(*) FROM analysis JOIN job USING(analysis_id) WHERE logic_name='paf_table_reuse'; ]
    newly loaded core dbs:       5                  [ SELECT count(*) FROM analysis JOIN job USING(analysis_id) WHERE logic_name='load_fresh_members'; ]

    total running time:         8.7 days            [ SELECT (UNIX_TIMESTAMP(max(died))-UNIX_TIMESTAMP(min(born)))/3600/24 FROM worker;  ]  # NB: stable_id mapping phase not included
    blasting time:              1.9 days            [ SELECT (UNIX_TIMESTAMP(max(died))-UNIX_TIMESTAMP(min(born)))/3600/24 FROM worker JOIN analysis USING (analysis_id) WHERE logic_name='blastp_with_reuse'; ]

=head2 rel.62 stats

    sequences to cluster:       1,192,544           [ SELECT count(*) from sequence; ]
    reused core dbs:            46                  [ number of 'load_reuse_members' jobs ]
    newly loaded core dbs:       7                  [ number of 'load_fresh_members' jobs ]

    total running time:         6 days              [ SELECT (UNIX_TIMESTAMP(max(died))-UNIX_TIMESTAMP(min(born)))/3600/24 FROM hive;  ]
    blasting time:              2.7 days            [ SELECT (UNIX_TIMESTAMP(max(died))-UNIX_TIMESTAMP(min(born)))/3600/24 FROM hive JOIN analysis USING (analysis_id) WHERE logic_name='blastp_with_reuse'; ]

=head2 rel.61 stats

    sequences to cluster:       1,173,469           [ SELECT count(*) from sequence; ]
    reused core dbs:            46                  [ number of 'load_reuse_members' jobs ]
    newly loaded core dbs:       6                  [ number of 'load_fresh_members' jobs ]

    total running time:         6 days              [ SELECT (UNIX_TIMESTAMP(max(died))-UNIX_TIMESTAMP(min(born)))/3600/24 FROM hive;  ]
    blasting time:              1.4 days            [ SELECT (UNIX_TIMESTAMP(max(died))-UNIX_TIMESTAMP(min(born)))/3600/24 FROM hive JOIN analysis USING (analysis_id) WHERE logic_name like 'blast%' or logic_name like 'SubmitPep%'; ]

=head1 CONTACT

  Please contact ehive-users@ebi.ac.uk mailing list with questions/suggestions.

=cut

package Bio::EnsEMBL::Compara::PipeConfig::ProteinTrees_conf;

use strict;
use warnings;
use base ('Bio::EnsEMBL::Compara::PipeConfig::ComparaGeneric_conf');


sub default_options {
    my ($self) = @_;
    return {
        %{$self->SUPER::default_options},   # inherit the generic ones


    # parameters that are likely to change from execution to another:
#       'mlss_id'               => 40073,   # it is very important to check that this value is current (commented out to make it obligatory to specify)
        'release'               => '63',
        'rel_suffix'            => '',    # an empty string by default, a letter otherwise
        'ensembl_cvs_root_dir'  => $ENV{'ENSEMBL_CVS_ROOT_DIR'}, # make sure you have this variable defined & exported in your shell configs
        'email'                 => $ENV{'USER'}.'@ebi.ac.uk',    # NB: your EBI address may differ from the Sanger one!
        'work_dir'              => '/lustre/scratch101/ensembl/'.$ENV{'USER'}.'/protein_trees_'.$self->o('rel_with_suffix'),

    # dependent parameters:
        'rel_with_suffix'       => $self->o('release').$self->o('rel_suffix'),
        'pipeline_name'         => 'PT_'.$self->o('rel_with_suffix'),   # name the pipeline to differentiate the submitted processes
        'fasta_dir'             => $self->o('work_dir') . '/blast_db',  # affects 'dump_subset_create_blastdb' and 'blastp_with_reuse'
        'cluster_dir'           => $self->o('work_dir') . '/cluster',


    # blast parameters:
        'blast_options'             => '-filter none -span1 -postsw -V=20 -B=20 -sort_by_highscore -warnings -cpus 1',
        'blast_tmp_dir'             => '',  # if empty, will use Blast Analysis' default

    # clustering parameters:
        'outgroups'                     => [119],   # affects 'hcluster_dump_input_per_genome'
        'clustering_max_gene_halfcount' => 750,     # (half of the previously used 'clutering_max_gene_count=1500) affects 'hcluster_run'

    # tree building parameters:
        'tree_max_gene_count'       => 400,     # affects 'mcoffee' and 'mcoffee_himem'
        'use_exon_boundaries'       => 0,       # affects 'mcoffee' and 'mcoffee_himem'
        'use_genomedb_id'           => 0,       # affects 'njtree_phyml' and 'ortho_tree'
        'species_tree_input_file'   => '',      # you can define your own species_tree for 'njtree_phyml' and 'ortho_tree'

    # homology_dnds parameters:
        'codeml_parameters_file'    => $self->o('ensembl_cvs_root_dir').'/ensembl-compara/scripts/pipeline/protein_trees.codeml.ctl.hash',      # used by 'homology_dNdS'
        'taxlevels'                 => ['Theria', 'Sauria', 'Tetraodontiformes'],
        'filter_high_coverage'      => 1,   # affects 'group_genomes_under_taxa'

    # executable locations:
        'wublastp_exe'              => 'wublastp',
        'hcluster_exe'              => '/software/ensembl/compara/hcluster/hcluster_sg',
        'mcoffee_exe'               => '/software/ensembl/compara/tcoffee-7.86b/t_coffee',
        'mafft_exe'                 => '',  # if empty, will use module built-in default
        'mafft_binaries'            => '',  # if empty, will use module built-in default
        'sreformat_exe'             => '/usr/local/ensembl/bin/sreformat',
        'treebest_exe'              => '/software/ensembl/compara/treebest/treebest',
        'quicktree_exe'             => '/software/ensembl/compara/quicktree_1.1/bin/quicktree',
        'buildhmm_exe'              => '/software/ensembl/compara/hmmer3/hmmer-3.0/src/hmmbuild',
        'codeml_exe'                => '/usr/local/ensembl/bin/codeml',


    # connection parameters to various databases:

        'pipeline_db' => {                      # the production database itself (will be created)
            -host   => 'compara4',
            -port   => 3306,
            -user   => 'ensadmin',
            -pass   => $self->o('password'),                    
            -dbname => $ENV{'USER'}.'_compara_homology_'.$self->o('rel_with_suffix'),
        },

        'master_db' => {                        # the master database for synchronization of various ids
            -host   => 'compara1',
            -port   => 3306,
            -user   => 'ensro',
            -pass   => '',
            -dbname => 'sf5_ensembl_compara_master',
        },

        'staging_loc1' => {                     # general location of half of the current release core databases
            -host   => 'ens-staging',
            -port   => 3306,
            -user   => 'ensro',
            -pass   => '',
        },

        'staging_loc2' => {                     # general location of the other half of the current release core databases
            -host   => 'ens-staging2',
            -port   => 3306,
            -user   => 'ensro',
            -pass   => '',
        },

        'livemirror_loc' => {                   # general location of the previous release core databases (for checking their reusability)
            -host   => 'ens-livemirror',
            -port   => 3306,
            -user   => 'ensro',
            -pass   => '',
        },


        # "production mode"
        'reuse_core_sources_locs'   => [ $self->o('livemirror_loc') ],
        'curr_core_sources_locs'    => [ $self->o('staging_loc1'), $self->o('staging_loc2'), ],
        'prev_release'              => 0,   # 0 is the default and it means "take current release number and subtract 1"
        'reuse_db' => {   # usually previous release database on compara1
           -host   => 'compara1',
           -port   => 3306,
           -user   => 'ensro',
           -pass   => '',
           -dbname => 'sf5_ensembl_compara_62',
        },

        ## mode for testing the non-Blast part of the pipeline: reuse all Blasts
        #'reuse_core_sources_locs' => [ $self->o('staging_loc1'), $self->o('staging_loc2'), ],
        #'curr_core_sources_locs'  => [ $self->o('staging_loc1'), $self->o('staging_loc2'), ],
        #'prev_release'            => $self->o('release'),
        #'reuse_db' => {   # current release if we are testing after production
        #    -host   => 'compara1',
        #    -port   => 3306,
        #    -user   => 'ensro',
        #    -pass   => '',
        #    -dbname => 'sf5_ensembl_compara_61',
        #},

    };
}


sub pipeline_wide_parameters {  # these parameter values are visible to all analyses, can be overridden by parameters{} and input_id{}
    my ($self) = @_;
    return {
        %{$self->SUPER::pipeline_wide_parameters},          # here we inherit anything from the base class

        'email'             => $self->o('email'),           # for (future) automatic notifications (may be unsupported by your Meadows)
    };
}


sub pipeline_create_commands {
    my ($self) = @_;
    return [
        @{$self->SUPER::pipeline_create_commands},  # here we inherit creation of database, hive tables and compara tables
        
        'mkdir -p '.$self->o('fasta_dir'),
        'lfs setstripe '.$self->o('fasta_dir').' -c -1',    # stripe
        'mkdir -p '.$self->o('cluster_dir'),
    ];
}


sub resource_classes {
    my ($self) = @_;
    return {
         0 => { -desc => 'default',          'LSF' => '' },
         1 => { -desc => 'hcluster_run',     'LSF' => '-C0 -M25000000 -q hugemem -R"select[mycompara2<500 && mem>25000] rusage[mycompara2=10:duration=10:decay=1:mem=25000]"' },
         2 => { -desc => 'mcoffee_himem',    'LSF' => '-C0 -M7500000 -R"select[mem>7500] rusage[mem=7500]"' },
    };
}


sub pipeline_analyses {
    my ($self) = @_;
    return [

# ---------------------------------------------[rename PAF tables]-----------------------------------------------------------------------

        {   -logic_name    => 'rename_paf_tables',
            -module        => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
            -parameters    => {
                'sql' => [ 'RENAME TABLE peptide_align_feature TO peptide_align_feature_orig',
                           'RENAME TABLE peptide_align_feature_prod TO peptide_align_feature',
                ],
            },
            -input_ids  => [
                { },
            ],
            -flow_into => {
                1 => [ 'copy_table_factory' ],
            },
        },

# ---------------------------------------------[copy tables from master and fix the offsets]---------------------------------------------

        {   -logic_name => 'copy_table_factory',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
            -parameters => {
                'db_conn'      => $self->o('master_db'),
                'inputlist'    => [ 'method_link', 'species_set', 'method_link_species_set', 'ncbi_taxa_name', 'ncbi_taxa_node' ],
                'column_names' => [ 'table' ],
                'input_id'     => { 'src_db_conn' => '#db_conn#', 'table' => '#table#' },
                'fan_branch_code' => 2,
            },
            -flow_into => {
                2 => [ 'copy_table'  ],
                1 => [ 'innodbise_table_factory' ],     # backbone
            },
        },

        {   -logic_name    => 'copy_table',
            -module        => 'Bio::EnsEMBL::Hive::RunnableDB::MySQLTransfer',
            -parameters    => {
                'mode'          => 'overwrite',
            },
            -hive_capacity => 10,
        },

# ---------------------------------------------[turn all tables except 'genome_db' to InnoDB]---------------------------------------------

        {   -logic_name => 'innodbise_table_factory',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
            -parameters => {
                'inputquery'      => "SELECT table_name FROM information_schema.tables WHERE table_schema ='".$self->o('pipeline_db','-dbname')."' AND table_name!='genome_db' AND engine='MyISAM' ",
                'fan_branch_code' => 2,
            },
            -wait_for  => [ 'copy_table' ],
            -flow_into => {
                2 => [ 'innodbise_table'  ],
                1 => [ 'generate_reuse_ss' ],           # backbone
            },
        },

        {   -logic_name    => 'innodbise_table',
            -module        => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
            -parameters    => {
                'sql'         => "ALTER TABLE #table_name# ENGINE=InnoDB",
            },
            -hive_capacity => 10,
            -can_be_empty  => 1,
        },

# ---------------------------------------------[generate an empty species_set for reuse (to be filled in at a later stage) ]---------

        {   -logic_name => 'generate_reuse_ss',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
            -parameters => {
                'sql' => [  "INSERT INTO species_set VALUES ()",   # inserts a dummy pair (auto_increment++, 0) into the table
                            "DELETE FROM species_set WHERE species_set_id=#_insert_id_0#", # will delete the row previously inserted, but keep the auto_increment
                ],
            },
            -wait_for  => [ 'innodbise_table_factory', 'innodbise_table' ], # have to wait for both, because subfan can be empty
            -flow_into => {
                2 => { 'mysql:////meta' => { 'meta_key' => 'reuse_ss_id', 'meta_value' => '#_insert_id_0#' } },     # dynamically record it as a pipeline-wide parameter
                1 => [ 'load_genomedb_factory' ],       # backbone
            },
        },

# ---------------------------------------------[load GenomeDB entries from master+cores]---------------------------------------------

        {   -logic_name => 'load_genomedb_factory',
            -module     => 'Bio::EnsEMBL::Compara::RunnableDB::LoadGenomedbFactory',
            -parameters => {
                'compara_db'    => $self->o('master_db'),   # that's where genome_db_ids come from
                'mlss_id'       => $self->o('mlss_id'),
            },
            -flow_into => {
                2 => [ 'load_genomedb' ],
                1 => { 'make_species_tree' => undef,
                       'accumulate_reuse_ss' => undef,  # backbone
                },
            },
        },

        {   -logic_name => 'load_genomedb',
            -module     => 'Bio::EnsEMBL::Compara::RunnableDB::LoadOneGenomeDB',
            -parameters => {
                'registry_dbs'  => $self->o('curr_core_sources_locs'),
            },
            -hive_capacity => 1,    # they are all short jobs, no point doing them in parallel
            -flow_into => {
                1 => [ 'check_reusability' ],   # each will flow into another one
            },
        },

# ---------------------------------------------[load species tree]-------------------------------------------------------------------

        {   -logic_name    => 'make_species_tree',
            -module        => 'Bio::EnsEMBL::Compara::RunnableDB::MakeSpeciesTree',
            -parameters    => {
                'species_tree_input_file' => $self->o('species_tree_input_file'),   # empty by default, but if nonempty this file will be used instead of tree generation from genome_db
            },
            -wait_for => [ 'load_genomedb' ],  # a funnel
            -hive_capacity => -1,   # to allow for parallelization
            -flow_into  => {
                3 => { 'mysql:////protein_tree_tag' => { 'node_id' => 1, 'tag' => 'species_tree_string', 'value' => '#species_tree_string#' } },
            },
        },

# ---------------------------------------------[filter genome_db entries into reusable and non-reusable ones]------------------------

        {   -logic_name => 'check_reusability',
            -module     => 'Bio::EnsEMBL::Compara::RunnableDB::ProteinTrees::CheckGenomedbReusability',
            -parameters => {
                'reuse_db'      => $self->o('reuse_db'),
                'registry_dbs'  => $self->o('reuse_core_sources_locs'),
                'release'       => $self->o('release'),
                'prev_release'  => $self->o('prev_release'),
            },
            -hive_capacity => 10,    # allow for parallel execution
            -flow_into => {
                2 => { 'subset_table_reuse'         => undef,
                       'paf_table_reuse'            => undef,
                       'mysql:////species_set'      => { 'genome_db_id' => '#genome_db_id#', 'species_set_id' => '#reuse_ss_id#' },
                },
                3 => [ 'load_fresh_members', 'paf_create_empty_table' ],
            },
        },

        {   -logic_name    => 'accumulate_reuse_ss',
            -module        => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',     # a non-standard use of JobFactory for iterative insertion
            -parameters => {
                'inputquery'      => 'SELECT "reuse_ss_csv" meta_key, GROUP_CONCAT(genome_db_id) meta_value FROM species_set WHERE species_set_id=#reuse_ss_id#',
                'fan_branch_code' => 3,
            },
            -wait_for => [ 'load_genomedb', 'check_reusability' ],
            -hive_capacity => -1,   # to allow for parallelization
            -flow_into => {
                3 => [ 'mysql:////meta' ],
            },
        },

# ---------------------------------------------[reuse members and pafs]--------------------------------------------------------------

        {   -logic_name => 'subset_table_reuse',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::MySQLTransfer',
            -parameters => {
                'src_db_conn'   => $self->o('reuse_db'),
                'table'         => 'subset',
                'mode'          => 'insertignore',
                'where'         => 'description LIKE "gdb:#genome_db_id# %"',
            },
            -wait_for => [ 'accumulate_reuse_ss' ], # to make sure some fresh members won't start because they were dataflown first (as this analysis can_be_empty)
            -hive_capacity => 4,
            -can_be_empty  => 1,
            -flow_into => {
                 1 => [ 'subset_member_table_reuse' ],
            },
        },

        {   -logic_name => 'subset_member_table_reuse',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
            -parameters => {
                            'db_conn'    => $self->o('reuse_db'),
                            'inputquery' => "SELECT sm.* FROM subset_member sm JOIN subset USING (subset_id) WHERE member_id<=100000000 AND description LIKE 'gdb:#genome_db_id# %'",
                            'fan_branch_code' => 2,
            },
            -hive_capacity => 4,
            -can_be_empty  => 1,
            -flow_into => {
                1 => [ 'member_table_reuse' ],
                2 => [ 'mysql:////subset_member' ],
            },
        },

        {   -logic_name => 'member_table_reuse',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::MySQLTransfer',
            -parameters => {
                'src_db_conn'   => $self->o('reuse_db'),
                'table'         => 'member',
                'where'         => 'member_id<=100000000 AND genome_db_id = #genome_db_id#',
                'mode'          => 'insertignore',
            },
            -wait_for => [ 'accumulate_reuse_ss' ], # to make sure some fresh members won't start because they were dataflown first (as this analysis can_be_empty)
            -hive_capacity => 4,
            -can_be_empty  => 1,
            -flow_into => {
                1 => [ 'sequence_table_reuse' ],
            },
        },

        {   -logic_name => 'sequence_table_reuse',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
            -parameters => {
                            'db_conn'    => $self->o('reuse_db'),
                            'inputquery' => 'SELECT s.* FROM sequence s JOIN member USING (sequence_id) WHERE sequence_id<=100000000 AND genome_db_id = #genome_db_id#',
                            'fan_branch_code' => 2,
            },
            -hive_capacity => 4,
            -can_be_empty  => 1,
            -flow_into => {
                1 => [ 'store_sequences_factory', 'dump_subset_create_blastdb' ],
                2 => [ 'mysql:////sequence' ],
            },
        },

        {   -logic_name => 'paf_table_reuse',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::MySQLTransfer',
            -parameters => {
                'src_db_conn'   => $self->o('reuse_db'),
                'table'         => 'peptide_align_feature_#per_genome_suffix#',
                'where'         => 'hgenome_db_id IN (#reuse_ss_csv#)',
            },
            -wait_for   => [ 'accumulate_reuse_ss' ],     # have to wait until reuse_ss_csv is computed
            -hive_capacity => 4,
            -can_be_empty  => 1,
        },

# ---------------------------------------------[load the rest of members]------------------------------------------------------------

        {   -logic_name => 'load_fresh_members',
            -module     => 'Bio::EnsEMBL::Compara::RunnableDB::ProteinTrees::ReuseOrLoadMembers',
            -parameters => { },
            -wait_for => [ 'accumulate_reuse_ss', 'subset_table_reuse', 'subset_member_table_reuse', 'member_table_reuse', 'sequence_table_reuse' ],
            -hive_capacity => -1,
            -can_be_empty  => 1,
            -flow_into => {
                1 => [ 'dump_subset_create_blastdb', 'store_sequences_factory' ],
            },
        },

        {   -logic_name => 'paf_create_empty_table',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
            -parameters => {
                'sql' => [  'CREATE TABLE IF NOT EXISTS peptide_align_feature_#per_genome_suffix# like peptide_align_feature',
                            'ALTER TABLE peptide_align_feature_#per_genome_suffix# DISABLE KEYS',
                ],
            },
            -can_be_empty  => 1,
        },


# ---------------------------------------------[create and populate blast analyses]--------------------------------------------------

        {   -logic_name => 'dump_subset_create_blastdb',
            -module     => 'Bio::EnsEMBL::Compara::RunnableDB::ProteinTrees::DumpSubsetCreateBlastDB',
            -parameters => {
                'fasta_dir'                 => $self->o('fasta_dir'),
            },
            -batch_size    =>  20,  # they can be really, really short
            -flow_into => {
                1 => [ 'blast_factory' ],
            },
        },

        {   -logic_name => 'blast_factory',
            -module     => 'Bio::EnsEMBL::Compara::RunnableDB::ProteinTrees::SubsetMemberFactory',
            -parameters => {
                'input_id' => { 'member_id' => '#_range_start#' },
                'fan_branch_code' => 2,
            },
            -hive_capacity => 10,
            -flow_into => {
                2 => [ 'blastp_with_reuse' ],
                1 => [ 'hcluster_dump_input_per_genome' ],
            },
        },

        {   -logic_name         => 'blastp_with_reuse',
            -module             => 'Bio::EnsEMBL::Compara::RunnableDB::ProteinTrees::BlastpWithReuse',
            -program_file       => $self->o('wublastp_exe'),
            -parameters         => {
                'mlss_id'                   => $self->o('mlss_id'),
                'reuse_db'                  => $self->o('reuse_db'),
                'blast_options'             => $self->o('blast_options'),
                'blast_tmp_dir'             => $self->o('blast_tmp_dir'),
                'fasta_dir'                 => $self->o('fasta_dir'),
            },
            -wait_for => [ 'load_fresh_members', 'dump_subset_create_blastdb', 'paf_table_reuse', 'paf_create_empty_table' ],
            -batch_size    =>  40,
            -hive_capacity => 450,
        },

# ---------------------------------------------[clustering step]---------------------------------------------------------------------

        {   -logic_name => 'hcluster_dump_input_per_genome',
            -module     => 'Bio::EnsEMBL::Compara::RunnableDB::ProteinTrees::HclusterPrepare',
            -parameters => {
                'mlss_id'       => $self->o('mlss_id'),
                'outgroups'     => $self->o('outgroups'),
                'cluster_dir'   => $self->o('cluster_dir'),
            },
            -wait_for => [ 'blastp_with_reuse' ],
            -hive_capacity => 4,
            -flow_into => {
                1 => [ 'per_genome_clusterset_qc' ],
            },
        },

        {   -logic_name    => 'hcluster_merge_inputs',
            -module        => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters    => {
                'cluster_dir'               => $self->o('cluster_dir'),
            },
            -input_ids => [
                { 'cmd' => 'cat #cluster_dir#/*.hcluster.txt > #cluster_dir#/hcluster.txt' },
                { 'cmd' => 'cat #cluster_dir#/*.hcluster.cat > #cluster_dir#/hcluster.cat' },
            ],
            -wait_for => [ 'hcluster_dump_input_per_genome' ],
            -hive_capacity => -1,   # to allow for parallelization
        },

        {   -logic_name    => 'hcluster_run',
            -module        => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters    => {
                'clustering_max_gene_halfcount' => $self->o('clustering_max_gene_halfcount'),
                'cluster_dir'                   => $self->o('cluster_dir'),
                'hcluster_exe'                  => $self->o('hcluster_exe'),
                'cmd'                           => '#hcluster_exe# -m #clustering_max_gene_halfcount# -w 0 -s 0.34 -O -C #cluster_dir#/hcluster.cat -o #cluster_dir#/hcluster.out #cluster_dir#/hcluster.txt',
            },
            -input_ids => [
                { },    # backbone
            ],
            -wait_for => [ 'hcluster_merge_inputs' ],
            -hive_capacity => -1,   # to allow for parallelization
            -flow_into => {
                1 => [ 'hcluster_parse_output' ],
            },
            -rc_id => 1,
        },

        {   -logic_name => 'hcluster_parse_output',
            -module     => 'Bio::EnsEMBL::Compara::RunnableDB::ProteinTrees::HclusterParseOutput',
            -parameters => {
                'mlss_id'                   => $self->o('mlss_id'),
                'cluster_dir'               => $self->o('cluster_dir'),
            },
            -hive_capacity => -1,
            -flow_into => {
                1 => [ 'overall_clusterset_qc' ],   # backbone 
                2 => [ 'mcoffee' ],                 # per-cluster
            },
        },

# ---------------------------------------------[sequence caching step]---------------------------------------------------------------

        {   -logic_name => 'store_sequences_factory',
            -module     => 'Bio::EnsEMBL::Compara::RunnableDB::PeptideMemberGroupingFactory',
            -parameters => { },
            -hive_capacity => -1,
            -flow_into => {
                2 => [ 'store_sequences' ],
            },
        },

        {   -logic_name => 'store_sequences',
            -module     => 'Bio::EnsEMBL::Compara::RunnableDB::ProteinTrees::FlowMemberSeq',
            -parameters => { },
            -hive_capacity => 200,
            -flow_into => {
                2 => [ 'mysql:////sequence_cds' ],
                3 => [ 'mysql:////sequence_exon_bounded' ],
            },
        },

# ---------------------------------------------[a QC step before main loop]----------------------------------------------------------

        {   -logic_name => 'overall_clusterset_qc',
            -module     => 'Bio::EnsEMBL::Compara::RunnableDB::ProteinTrees::OverallGroupsetQC',
            -parameters => {
                'reuse_db'                  => $self->o('reuse_db'),
                'cluster_dir'               => $self->o('cluster_dir'),
                'groupset_tag'              => 'ClustersetQC',
            },
            -hive_capacity => 3,
            -flow_into => {
                1 => [ 'overall_genetreeset_qc' ],
            },
        },

        {   -logic_name => 'per_genome_clusterset_qc',
            -module     => 'Bio::EnsEMBL::Compara::RunnableDB::ProteinTrees::PerGenomeGroupsetQC',
            -parameters => {
                'reuse_db'                  => $self->o('reuse_db'),
                'groupset_tag'              => 'ClustersetQC',
            },
            -wait_for => [ 'hcluster_parse_output' ],
            -hive_capacity => 3,
            -flow_into => {
                1 => [ 'per_genome_genetreeset_qc' ],
            },
        },

# ---------------------------------------------[main tree creation loop]-------------------------------------------------------------

        {   -logic_name => 'mcoffee',
            -module     => 'Bio::EnsEMBL::Compara::RunnableDB::ProteinTrees::MCoffee',
            -program_file       => $self->o('mcoffee_exe'),
            -parameters => {
                'method'                    => 'cmcoffee',      # presumably, at the moment it refers to the 'initial' method
                'use_exon_boundaries'       => $self->o('use_exon_boundaries'),
                'max_gene_count'            => $self->o('tree_max_gene_count'),
                'mafft_exe'                 => $self->o('mafft_exe'),
                'mafft_binaries'            => $self->o('mafft_binaries'),
            },
            -wait_for => [ 'store_sequences', 'overall_clusterset_qc', 'per_genome_clusterset_qc' ],    # funnel
            -hive_capacity        => 600,
            -flow_into => {
               -1 => [ 'mcoffee_himem' ],
                1 => [ 'njtree_phyml' ],
                3 => [ 'quick_tree_break' ],
            },
        },

        {   -logic_name => 'mcoffee_himem',
            -module     => 'Bio::EnsEMBL::Compara::RunnableDB::ProteinTrees::MCoffee',
            -program_file       => $self->o('mcoffee_exe'),
            -parameters => {
                'method'                    => 'cmcoffee',      # presumably, at the moment it refers to the 'initial' method
                'use_exon_boundaries'       => $self->o('use_exon_boundaries'),
                'max_gene_count'            => $self->o('tree_max_gene_count'),
                'mafft_exe'                 => $self->o('mafft_exe'),
                'mafft_binaries'            => $self->o('mafft_binaries'),
            },
            -hive_capacity        => 600,
            -can_be_empty         => 1,
            -flow_into => {
                1 => [ 'njtree_phyml' ],
                3 => [ 'quick_tree_break' ],
            },
            -rc_id => 2,
        },

        {   -logic_name => 'njtree_phyml',
            -module     => 'Bio::EnsEMBL::Compara::RunnableDB::ProteinTrees::NJTREE_PHYML',
            -program_file       => $self->o('treebest_exe'),
            -parameters => {
                'cdna'                      => 1,
                'bootstrap'                 => 1,
                'use_genomedb_id'           => $self->o('use_genomedb_id'),
            },
            -hive_capacity        => 400,
            -failed_job_tolerance => 5,
            -flow_into => {
                1 => [ 'ortho_tree' ],
                2 => [ 'njtree_phyml' ],
                # 3 => [ 'quick_tree_break' ],      # in most cases we do not want quick_tree_break to happen automatically!
            },
        },

        {   -logic_name => 'ortho_tree',
            -module     => 'Bio::EnsEMBL::Compara::RunnableDB::ProteinTrees::OrthoTree',
            -parameters => {
                'use_genomedb_id'           => $self->o('use_genomedb_id'),
            },
            -hive_capacity        => 200,
            -failed_job_tolerance => 5,
            -flow_into => {
                1 => [ 'build_HMM_aa', 'build_HMM_cds' ],
            },
        },

        {   -logic_name => 'build_HMM_aa',
            -module     => 'Bio::EnsEMBL::Compara::RunnableDB::ProteinTrees::BuildHMM',
            -program_file       => $self->o('buildhmm_exe'),
            -parameters => {
                'sreformat_exe'     => $self->o('sreformat_exe'),
            },
            -hive_capacity        => 200,
            -failed_job_tolerance => 5,
        },

        {   -logic_name => 'build_HMM_cds',
            -module     => 'Bio::EnsEMBL::Compara::RunnableDB::ProteinTrees::BuildHMM',
            -program_file       => $self->o('buildhmm_exe'),
            -parameters => {
                'cdna'              => 1,
                'sreformat_exe'     => $self->o('sreformat_exe'),
            },
            -hive_capacity        => 200,
            -failed_job_tolerance => 5,
        },

        {   -logic_name => 'quick_tree_break',
            -module     => 'Bio::EnsEMBL::Compara::RunnableDB::ProteinTrees::QuickTreeBreak',
            -program_file       => $self->o('quicktree_exe'),
            -parameters => {
                'sreformat_exe'     => $self->o('sreformat_exe'),
                'mlss_id'           => $self->o('mlss_id'),
            },
            -hive_capacity        => 1, # this one seems to slow the whole loop down; why can't we have any more of these?
            -can_be_empty         => 1,
            -failed_job_tolerance => 5,
            -flow_into => {
                1 => [ 'other_paralogs' ],
                2 => [ 'mcoffee' ],
            },
        },

        {   -logic_name => 'other_paralogs',
            -module     => 'Bio::EnsEMBL::Compara::RunnableDB::ProteinTrees::OtherParalogs',
            -parameters => { },
            -wait_for => [ 'mcoffee', 'mcoffee_himem', 'njtree_phyml', 'ortho_tree', 'build_HMM_aa', 'build_HMM_cds', 'quick_tree_break' ],
            -hive_capacity        => 50,
            -failed_job_tolerance => 5,
        },

# ---------------------------------------------[a QC step after main loop]----------------------------------------------------------

        {   -logic_name => 'overall_genetreeset_qc',
            -module     => 'Bio::EnsEMBL::Compara::RunnableDB::ProteinTrees::OverallGroupsetQC',
            -parameters => {
                'reuse_db'                  => $self->o('reuse_db'),
                'cluster_dir'               => $self->o('cluster_dir'),
                'groupset_tag'              => 'GeneTreesetQC',
            },
            -wait_for => [ 'mcoffee', 'mcoffee_himem', 'njtree_phyml', 'ortho_tree', 'quick_tree_break' ],
            -hive_capacity => 3,
            -flow_into => {
                1 => [ 'group_genomes_under_taxa' ],    # backbone
            },
        },

        {   -logic_name => 'per_genome_genetreeset_qc',
            -module     => 'Bio::EnsEMBL::Compara::RunnableDB::ProteinTrees::PerGenomeGroupsetQC',
            -parameters => {
                'reuse_db'                  => $self->o('reuse_db'),
                'groupset_tag'              => 'GeneTreesetQC',
            },
            -wait_for => [ 'mcoffee', 'mcoffee_himem', 'njtree_phyml', 'ortho_tree', 'quick_tree_break' ],
            -hive_capacity => 3,
        },

# ---------------------------------------------[homology step]-----------------------------------------------------------------------

        {   -logic_name => 'group_genomes_under_taxa',
            -module     => 'Bio::EnsEMBL::Compara::RunnableDB::ProteinTrees::GroupGenomesUnderTaxa',
            -parameters => {
                'mlss_id'               => $self->o('mlss_id'),
                'taxlevels'             => $self->o('taxlevels'),
                'filter_high_coverage'  => $self->o('filter_high_coverage'),
            },
            -wait_for => [ 'per_genome_genetreeset_qc' ],   # funnel
            -hive_capacity => -1,
            -flow_into => {
                2 => [ 'homology_dNdS_factory' ],
            },
        },

        {   -logic_name => 'homology_dNdS_factory',
            -module     => 'Bio::EnsEMBL::Compara::RunnableDB::HomologyGroupingFactory',
            -parameters => {
		    'method_link_types'  => ['ENSEMBL_ORTHOLOGUES', 'ENSEMBL_PARALOGUES'],
		},
            -hive_capacity => -1,
            -flow_into => {
                1 => [ 'threshold_on_dS' ],
                2 => [ 'homology_dNdS' ],
            },
        },

        {   -logic_name => 'homology_dNdS',
            -module     => 'Bio::EnsEMBL::Compara::RunnableDB::ProteinTrees::Homology_dNdS',
            -program_file       => $self->o('codeml_exe'),
            -parameters => {
                'codeml_parameters_file'    => $self->o('codeml_parameters_file'),
            },
            -hive_capacity        => 200,
            -failed_job_tolerance => 2,
        },

        {   -logic_name => 'threshold_on_dS',
            -module     => 'Bio::EnsEMBL::Compara::RunnableDB::ProteinTrees::Threshold_on_dS',
            -parameters => { },
            -wait_for => [ 'homology_dNdS' ],
            -hive_capacity => -1,
        },

# ---------------------------------------------[homology duplications QC step]-------------------------------------------------------

        {   -logic_name => 'homology_duplications_factory',
            -module     => 'Bio::EnsEMBL::Compara::RunnableDB::MLSSfactory',
            -parameters => {
                'input_id' => { 'type' => '#short_type#', 'mlss_id' => '#_range_start#' },
                'fan_branch_code' => 2,
            },
            -wait_for => [ 'threshold_on_dS' ],
            -input_ids => [
                { 'method_link_type' => 'ENSEMBL_ORTHOLOGUES', 'short_type' => 'orthologues' },
                { 'method_link_type' => 'ENSEMBL_PARALOGUES',  'short_type' => 'paralogues'  },
            ],
            -hive_capacity => -1,
            -flow_into => {
                2 => [ 'homology_duplications' ],
            },
        },

        {   -logic_name => 'homology_duplications',
            -module     => 'Bio::EnsEMBL::Compara::RunnableDB::ProteinTrees::HDupsQC',
            -parameters => { },
            -hive_capacity => 10,
        },

    ];
}

1;

