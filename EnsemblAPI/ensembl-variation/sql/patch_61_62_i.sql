# change the consequence_type column in variation_feature to use SO accessions

ALTER TABLE variation_feature MODIFY consequence_type SET (
    'intergenic_variant',
    'splice_acceptor_variant',
    'splice_donor_variant',
    'complex_change_in_transcript', 
    'stop_lost',
    'coding_sequence_variant',
    'non_synonymous_codon',
    'stop_gained',
    'synonymous_codon',
    'frameshift_variant',
    'nc_transcript_variant',
    'mature_miRNA_variant',
    'NMD_transcript_variant',
    '5_prime_UTR_variant',
    '3_prime_UTR_variant',
    'incomplete_terminal_codon_variant',
    'intron_variant',
    'splice_region_variant',
    '5KB_downstream_variant',
    '500B_downstream_variant',
    '5KB_upstream_variant',
    '2KB_upstream_variant',
    'initiator_codon_change',
    'stop_retained_variant',
    'inframe_codon_gain',
    'inframe_codon_loss',
    'miRNA_target_site_variant',
    'pre_miRNA_variant',
    'regulatory_region_variant',
    'increased_binding_affinity',
    'decreased_binding_affinity',
    'binding_site_variant'
) DEFAULT 'intergenic_variant' NOT NULL;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_61_62_i.sql|change the consequence_type column in variation_feature to use SO terms');
