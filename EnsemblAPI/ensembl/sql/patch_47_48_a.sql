# patch_47_48_a.sql
#
# title: update schema version
#
# description:
# Update schema_version in meta table to 48

UPDATE meta SET meta_value='48' WHERE meta_key='schema_version';

# patch identifier
INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_47_48_a.sql|schema_version');


