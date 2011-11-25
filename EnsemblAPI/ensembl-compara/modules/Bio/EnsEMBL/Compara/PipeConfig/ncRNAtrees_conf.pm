
=pod 

=head1 NAME

  Bio::EnsEMBL::Compara::PipeConfig::ncRNAtrees_conf

=head1 SYNOPSIS

    init_pipeline.pl Bio::EnsEMBL::Compara::PipeConfig::ncRNAtrees_conf -password <your_password>

=head1 DESCRIPTION  

    This is an experimental PipeConfig file for ncRNAtrees pipeline (work in progress)

=head1 CONTACT

  Please contact ehive-users@ebi.ac.uk mailing list with questions/suggestions.

=cut

package Bio::EnsEMBL::Compara::PipeConfig::ncRNAtrees_confI;

use strict;
use warnings;
use base ('Bio::EnsEMBL::Compara::PipeConfig::ComparaGeneric_conf');

sub default_options {
    my ($self) = @_;
    return {
        %{$self->SUPER::default_options},

        'mlss_id'           => 40074,
        'max_gene_count'    => 1500,

        'release'           => '63',
        'rel_suffix'        => 'b',    # an empty string by default, a letter otherwise
        'rel_with_suffix'   => $self->o('release').$self->o('rel_suffix'),

        'ensembl_cvs_root_dir' => $ENV{'ENSEMBL_CVS_ROOT_DIR'},
        'work_dir'             => $ENV{'HOME'}.'/ncrna_trees_'.$self->o('rel_with_suffix'),

        'email'             => $ENV{'USER'}.'@ebi.ac.uk',    # NB: your EBI address may differ from the Sanger one!

        'species_tree_input_file'   => '',  # empty value means 'create using genome_db+ncbi_taxonomy information'; can be overriden by a file with a tree in it

        'pipeline_db' => {                                  # connection parameters
                          -driver => 'mysql',
                          -host   => 'compara3',
                          -port   => 3306,
                          -user   => 'ensadmin',
                          -pass   => $self->o('password'),
                          -dbname => $ENV{'USER'}.'_compara_nctrees_'.$self->o('rel_with_suffix'),
        },

        'reg1' => {
            -host   => 'ens-staging',
            -port   => 3306,
            -user   => 'ensro',
            -pass   => '',
        },

        'reg2' => {
            -host   => 'ens-staging2',
            -port   => 3306,
            -user   => 'ensro',
            -pass   => '',
        },

        'master_db' => {
            -host   => 'compara1',
            -port   => 3306,
            -user   => 'ensro',
            -pass   => '',
            -dbname => 'sf5_ensembl_compara_master', # 'sf5_ensembl_compara_master',
        },

        'epo_db' => {   # ideally, the current release database with epo pipeline results already loaded
            -host   => 'compara1',
            -port   => 3306,
            -user   => 'ensro',
            -pass   => '',
            -dbname => 'sf5_ensembl_compara_62',
        },
    };
}


sub pipeline_wide_parameters {  # these parameter values are visible to all analyses, can be overridden by parameters{} and input_id{}
    my ($self) = @_;
    return {
        'pipeline_name'     => 'NCT_'.$self->o('rel_with_suffix'),  # name the pipeline to differentiate the submitted processes
        'email'             => $self->o('email'),                   # for automatic notifications (may be unsupported by your Meadows)
        'work_dir'          => $self->o('work_dir'),                # data directories and filenames
    };
}

sub pipeline_create_commands {
    my ($self) = @_;
    return [
        @{$self->SUPER::pipeline_create_commands},  # here we inherit creation of database, hive tables and compara tables
        'mkdir -p '.$self->o('work_dir'),
    ];
}

sub resource_classes {
    my ($self) = @_;
    return {
            0 => { -desc => 'default', 'LSF' => '' },
            1 => { -desc => 'himem'  , 'LSF' => '-q hugemem -M15000000 -R"select[mem>15000] rusage[mem=15000]"' },
            2 => { -desc => 'long'   , 'LSF' => '-q long' },
           };
}

sub pipeline_analyses {
    my ($self) = @_;
    return [

# ---------------------------------------------[copy tables from master and fix the offsets]---------------------------------------------

        {   -logic_name => 'copy_table_factory',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
            -parameters => {
                'db_conn'   => $self->o('master_db'),
                'inputlist' => [ 'method_link', 'species_set', 'method_link_species_set', 'ncbi_taxa_name', 'ncbi_taxa_node' ],
                'column_names' => [ 'table' ],
                'input_id'  => { 'src_db_conn' => '#db_conn#', 'table' => '#table#' },
                'fan_branch_code' => 2,
            },
            -input_ids => [
                {},
            ],
            -flow_into => {
                2 => [ 'copy_table'  ],
                1 => [ 'offset_tables' ],  # backbone
            },
        },

        {   -logic_name    => 'copy_table',
            -module        => 'Bio::EnsEMBL::Hive::RunnableDB::MySQLTransfer',
            -parameters    => {
                'mode'          => 'overwrite',
            },
            -hive_capacity => 10,
        },

        {   -logic_name => 'offset_tables',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
            -parameters => {
                'sql'   => [
                    'ALTER TABLE member   AUTO_INCREMENT=200000001',
                    'ALTER TABLE sequence AUTO_INCREMENT=200000001',
                    'ALTER TABLE homology AUTO_INCREMENT=100000001',
                ],
            },
            -wait_for => [ 'copy_table_factory', 'copy_table' ],    # have to wait until the tables have been copied
            -flow_into => {
                1 => [ 'innodbise_table_factory' ],
            },
        },

# ---------------------------------------------[turn all tables except 'genome_db' to InnoDB]---------------------------------------------

        {   -logic_name => 'innodbise_table_factory',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
            -parameters => {
                'inputquery'      => "SELECT table_name FROM information_schema.tables WHERE table_schema ='".$self->o('pipeline_db','-dbname')."' AND table_name!='genome_db' AND engine='MyISAM' ",
                'fan_branch_code' => 2,
            },
            -flow_into => {
                2 => [ 'innodbise_table'  ],
                1 => [ 'load_genomedb_factory' ],
            },
        },

        {   -logic_name    => 'innodbise_table',
            -module        => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
            -parameters    => {
                'sql'         => "ALTER TABLE #table_name# ENGINE=InnoDB",
            },
            -hive_capacity => 10,
        },

# ---------------------------------------------[load GenomeDB entries from master+cores]---------------------------------------------

        {   -logic_name => 'load_genomedb_factory',
            -module     => 'Bio::EnsEMBL::Compara::RunnableDB::LoadGenomedbFactory',
            -parameters => {
                'compara_db'    => $self->o('master_db'),   # that's where genome_db_ids come from
                'mlss_id'       => $self->o('mlss_id'),
            },
            -wait_for  => [ 'innodbise_table_factory', 'innodbise_table' ],
            -flow_into => {
                2 => [ 'load_genomedb' ],
                1 => [ 'make_species_tree', 'create_lca_species_set', 'load_rfam_models' ],
            },
        },

        {   -logic_name => 'load_genomedb',
            -module     => 'Bio::EnsEMBL::Compara::RunnableDB::LoadOneGenomeDB',
            -parameters => {
                'registry_dbs'  => [ $self->o('reg1'), $self->o('reg2'), ],
            },
            -hive_capacity => 1,    # they are all short jobs, no point doing them in parallel
            -flow_into => {
                1 => [ 'load_members_factory' ],   # each will flow into another one
            },
        },

# ---------------------------------------------[load species tree]-------------------------------------------------------------------


        {   -logic_name    => 'make_species_tree',
            -module        => 'Bio::EnsEMBL::Compara::RunnableDB::MakeSpeciesTree',
            -parameters    => {
                'species_tree_input_file' => $self->o('species_tree_input_file'),   # empty by default, but if nonempty this file will be used instead of tree generation from genome_db
                'multifurcation_deletes_node'           => [ 33316, 129949, 314146 ],
                'multifurcation_deletes_all_subnodes'   => [  9347, 186625,  32561 ],
            },
            -wait_for => [ 'load_genomedb' ],  # have to wait for both to complete (so is a funnel)
            -hive_capacity => -1,   # to allow for parallelization
            -flow_into  => {
                3 => { 'mysql:////nc_tree_tag' => { 'node_id' => 1, 'tag' => 'species_tree_string', 'value' => '#species_tree_string#' } },
            },
        },


# ---------------------------------------------[create the low-coverage-assembly species set]-----------------------------------------

        {   -logic_name => 'create_lca_species_set',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
            -parameters => {
                'sql' => [  "INSERT INTO species_set VALUES ()",   # insert a dummy pair (auto_increment++, 0) into the table
                            "DELETE FROM species_set WHERE species_set_id IN (#_insert_id_0#)",     # delete the previously inserted row, but keep the auto_increment
                ],
            },
            -hive_capacity => -1,   # to allow for parallelization
            -flow_into => {
                2 => {
                    'store_lca_species_set'     => { 'lca_species_set_id' => '#_insert_id_0#' },     # pass it on to the query
                    'mysql:////species_set_tag' => { 'species_set_id' => '#_insert_id_0#', 'tag' => 'name', 'value' => 'low-coverage-assembly' },   # record the id in ss_tag table
                },
            },
        },

        {   -logic_name => 'store_lca_species_set',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',    # another non-stardard use of JobFactory for iterative insertion
            -parameters => {
                'db_conn'         => $self->o('epo_db'),
                'inputquery'      => "SELECT DISTINCT g.genome_db_id FROM genome_db g JOIN species_set ss USING(genome_db_id) JOIN method_link_species_set mlss USING(species_set_id) WHERE assembly_default AND mlss.name LIKE '%EPO_LOW_COVERAGE%' AND g.genome_db_id NOT IN (SELECT DISTINCT(g2.genome_db_id) FROM genome_db g2 JOIN species_set ss2 USING(genome_db_id) JOIN method_link_species_set mlss2 USING(species_set_id) WHERE assembly_default AND mlss2.name LIKE '%EPO')",
                'input_id'        => { 'species_set_id' => '#lca_species_set_id#', 'genome_db_id' => '#genome_db_id#' },
                'fan_branch_code' => 3,
            },
            -hive_capacity => -1,   # to allow for parallelization
            -flow_into => {
                3 => [ 'mysql:////species_set' ],
            },
        },

# ---------------------------------------------[load ncRNA and gene members and subsets]---------------------------------------------

        {   -logic_name    => 'load_members_factory',
            -module        => 'Bio::EnsEMBL::Compara::RunnableDB::ncRNAtrees::GenomePrepareNCMembers',
            -hive_capacity => 10,
            -flow_into => {
                2 => [ 'load_members' ],   # per-genome fan
            },
        },

        {   -logic_name    => 'load_members',
            -module        => 'Bio::EnsEMBL::Compara::RunnableDB::ncRNAtrees::GeneStoreNCMembers',
            -hive_capacity => 30,
            -flow_into => {
                3 => [ 'mysql:////subset_member' ],   # every ncrna member is added to the corresponding subset
                4 => [ 'mysql:////subset_member' ],   # every gene  member is added to the corresponding subset
            },
        },

# ---------------------------------------------[load RFAM models]---------------------------------------------------------------------

        {   -logic_name    => 'load_rfam_models',
            -module        => 'Bio::EnsEMBL::Compara::RunnableDB::ncRNAtrees::RFAMLoadModels',
            -hive_capacity => -1,   # to allow for parallelization
            -flow_into => {
                1 => [ 'rfam_classify' ],
            },
        },

# ---------------------------------------------[run RFAM classification]--------------------------------------------------------------

        {   -logic_name    => 'rfam_classify',
            -module        => 'Bio::EnsEMBL::Compara::RunnableDB::ncRNAtrees::RFAMClassify',
            -parameters    => {
                'mlss_id'        => $self->o('mlss_id'),
            },
            -wait_for => [ 'make_species_tree', 'create_lca_species_set', 'store_lca_species_set', 'load_members_factory', 'load_members' ], # mega-funnel
            -flow_into => {
                           2 => [ 'recover_epo' ],
                           1 => ['db_snapshot_after_Rfam_classify']
            },
        },

# ---------------------------------------------[by-cluster branches]----------------------------------------------------------------------

            {   -logic_name => 'db_snapshot_after_Rfam_classify',
                -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
                -parameters => {
                                'cmd'      => 'mysqldump '.$self->dbconn_2_mysql('pipeline_db', 0).' '.$self->o('pipeline_db','-dbname').' >#filename#',
                                'filename'  => $ENV{'HOME'}.'/db_snapshot_after_Rfam_classify',
                               },
            },

#--------------------------------------------------------------------------------

        {   -logic_name    => 'recover_epo',
            -module        => 'Bio::EnsEMBL::Compara::RunnableDB::ncRNAtrees::NCRecoverEPO',
            -parameters    => {
                'mlss_id'        => $self->o('mlss_id'),
                'epo_db'         => $self->o('epo_db'),
            },
            -hive_capacity => 100,
            -wait_for => ['db_snapshot_after_Rfam_classify'],
            -flow_into => {
                           1 => [ 'infernal','genomic_alignment' ],
            },
        },

#         {   -logic_name    => 'recover_search',
#             -module        => 'Bio::EnsEMBL::Compara::RunnableDB::ncRNAtrees::NCRecoverSearch',
#             -batch_size    => 5,
#             -hive_capacity => -1,
#             -flow_into => {
#                 1 => [ 'infernal' ],
#             },
#         },

        {   -logic_name    => 'infernal',
            -module        => 'Bio::EnsEMBL::Compara::RunnableDB::ncRNAtrees::Infernal',
            -hive_capacity => 200,
            -failed_job_tolerance => 10,    # that many per cent jobs are allowed to fail
            -flow_into => {
                           1 => ['pre_sec_struct_tree' ],
                          },
        },

            {   -logic_name    => 'pre_sec_struct_tree', ## pre_sec_struct_tree
            -module        => 'Bio::EnsEMBL::Compara::RunnableDB::ncRNAtrees::PrepareSecStructModels',  ## PrepareRAxMLSecModels -- rename
            -hive_capacity => 200,
            #            -failed_job_tolerance =>  5,   ### No failed_job_tolerance here
            -flow_into => {
                           1 => [ 'treebest_mmerge' ],
                           2 => [ 'sec_struct_model_tree'],
                           -1 => [ 'pre_sec_struct_tree_himem' ],
                           -2 => [ 'sec_struct_model_tree_himem' ],
                          },
        },

        {
         -logic_name => 'pre_sec_struct_tree_himem',
         -module => 'Bio::EnsEMBL::Compara::RunnableDB::ncRNAtrees::PrepareSecStructModels',
         -hive_capacity => 200,
         -flow_into => {
                        1 => [ 'treebest_mmerge' ],
                        2 => [ 'sec_struct_model_tree_himem' ],
                       },
         -can_be_empty => 1,
         -rc_id => 1,
        },

        {   -logic_name    => 'sec_struct_model_tree', ## sec_struct_model_tree
            -module        => 'Bio::EnsEMBL::Compara::RunnableDB::ncRNAtrees::SecStructModelTree', ## SecStrucModels
            -hive_capacity => 200,
            -failed_job_tolerance => 3,
            -flow_into => {
                           -1 => [ 'sec_struct_model_tree_himem' ],
                           -2 => [ 'sec_struct_model_tree_himem' ],
                          },
        },

        {
         -logic_name => 'sec_struct_model_tree_himem',
         -module => 'Bio::EnsEMBL::Compara::RunnableDB::ncRNAtrees::SecStructModelTree',
         -hive_capacity => 200,
         -failed_job_tolerance => 3,
         -can_be_empty => 1,
         -rc_id => 1,
        },

        {   -logic_name    => 'genomic_alignment',
            -module        => 'Bio::EnsEMBL::Compara::RunnableDB::ncRNAtrees::NCGenomicAlignment',
            -hive_capacity => 200,
            -failed_job_tolerance => 5,    # that many per cent jobs are allowed to fail
            -flow_into => {
                           -2 => ['genomic_alignment_long'],
                          },
        },

        {
         -logic_name => 'genomic_alignment_long',
         -module => 'Bio::EnsEMBL::Compara::RunnableDB::ncRNAtrees::NCGenomicAlignment',
         -hive_capacity => 200,
         -failed_job_tolerance => 5,
         -can_be_empty => 1,
         -rc_id => 2,
        },

        {   -logic_name    => 'treebest_mmerge',
            -module        => 'Bio::EnsEMBL::Compara::RunnableDB::ncRNAtrees::NCTreeBestMMerge',
            -hive_capacity => 400,
            -wait_for      => ['recover_epo', 'infernal','genomic_alignment', 'genomic_alignment_long', 'sec_struct_model_tree','sec_struct_model_tree_himem', ],
            -failed_job_tolerance => 5,
            -flow_into => {
                           1 => [ 'orthotree', 'ktreedist' ],
                           -1 => [ 'treebest_mmerge_himem' ],
            },
        },

        {
         -logic_name => 'treebest_mmerge_himem',
         -module => 'Bio::EnsEMBL::Compara::RunnableDB::ncRNAtrees::NCTreeBestMMerge',
         -hive_capacity => 400,
         -failed_job_tolerance => 5,
         -flow_into => {
                        1 => [ 'orthotree', 'ktreedist' ],
                       },
         -rc_id => 1,
        },

        {   -logic_name    => 'orthotree',
            -module        => 'Bio::EnsEMBL::Compara::RunnableDB::ncRNAtrees::NCOrthoTree',
            -hive_capacity => 200,
            -flow_into => {
                           -1 => ['orthotree_himem' ],
                          },
        },

        {
         -logic_name => 'orthotree_himem',
         -module => 'Bio::EnsEMBL::Compara::RunnableDB::ncRNAtrees::NCOrthotree',
         -hive_capacity => 200,
         -failed_job_tolerance => 5,
         -rc_id => 1,
        },

        {   -logic_name    => 'ktreedist',
            -module        => 'Bio::EnsEMBL::Compara::RunnableDB::ncRNAtrees::Ktreedist',
            -hive_capacity => -1,
            -failed_job_tolerance =>  5,    # that many per cent jobs are allowed to fail
            -flow_into => {
                           -1 => [ 'ktreedist_himem' ],
                          },
        },

        {
         -logic_name => 'ktreedist_himem',
         -module => 'Bio::EnsEMBL::Compara::RunnableDB::ncRNAtrees::Ktreedist',
         -hive_capacity => -1,
         -failed_job_tolerance => 5,
        },

    ];
}

1;

