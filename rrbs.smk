import os
from pathlib import Path
import time
import re
import shutil

# this is required to be set before importing pandas and numpy to avoid "OpenBLAS blas_thread_init: pthread_create failed for thread N of Ns: Resource temporarily unavailable"
# source: https://stackoverflow.com/questions/52026652/openblas-blas-thread-init-pthread-create-resource-temporarily-unavailable
os.environ['OPENBLAS_NUM_THREADS'] = '2'  # '1'

import pandas as pd
import numpy as np
import yagmail
from dotenv import load_dotenv


configfile: str(Path(workflow.basedir) / "config/config.yaml")
include: str(Path(workflow.basedir) / "scripts/common_functions.py")
# include: str(Path(workflow.basedir) / "scripts/spladder_functions.py")

snakemake_path = workflow.basedir  
load_dotenv(str(Path(workflow.basedir) / "config/.env"))

abort_execution = False
pre_run_messages = {'info': [], 'warning': []}

help = str_to_bool(replace_cfg_with_environment('help', get_config_value ('help')))
if help:
    print_help_info()
    sys.exit(1)

default_max_memory = get_config_value ('resources/max_memory')  # default max memory to be allocated for a rule
default_rule_memory = 10000  # default memory to be allocated for a rule; this could be overwritten by the config value of "rules/{rule_name}/memory" variable
default_retries = 0  # default retries to be applied for a rule; this could be overwritten by the config value of "rules/{rule_name}/retries" variable
# spladder_events_ase_edge_limit_decrease_step_default = 100
# spladder_events_ase_edge_limit_default = 300

db_cfg = get_config_value ('database/connection') # Extracts DB connection related setting from config to be passed to the steps requiring DB connectivity

debug = str_to_bool(replace_cfg_with_environment('debug', get_config_value ('debug')))
# run_alternative_splicing = str_to_bool(replace_cfg_with_environment('run_alternative_splicing', get_config_value ('run_alternative_splicing')))
run_rsem = str_to_bool(replace_cfg_with_environment('run_rsem', get_config_value ('run_rsem')))
run_rmats_turbo = str_to_bool(replace_cfg_with_environment('run_rmats_turbo', get_config_value ('run_rmats_turbo')))
run_rmats_novel = str_to_bool(replace_cfg_with_environment('run_rmats_novel', get_config_value ('run_rmats_novel')))
# run_leafcutter = str_to_bool(replace_cfg_with_environment('run_leafcutter', get_config_value ('run_leafcutter')))
run_spladder = str_to_bool(replace_cfg_with_environment('run_spladder', get_config_value ('run_spladder')))
run_spladder_long_execution_events = str_to_bool(replace_cfg_with_environment('run_spladder_long_execution_events', get_config_value ('run_spladder_long_execution_events')))
run_djexpress = str_to_bool(replace_cfg_with_environment('run_djexpress', get_config_value ('run_djexpress')))
run_ngscheckmate = str_to_bool(replace_cfg_with_environment('run_ngscheckmate', get_config_value ('run_ngscheckmate')))
run_feature_counts = str_to_bool(replace_cfg_with_environment('run_feature_counts', get_config_value ('run_feature_counts')))
run_telescope = str_to_bool(replace_cfg_with_environment('run_telescope', get_config_value ('run_telescope')))

create_rds_for_rscripts = str_to_bool(replace_cfg_with_environment('create_rds_for_rscripts', get_config_value ('create_rds_for_rscripts')))
print_prerun_info = str_to_bool(replace_cfg_with_environment('print_prerun_info', get_config_value ('print_prerun_info')))

data_path = replace_cfg_with_environment('data_path', get_config_value ('data_path'))
data_path_required_folders = get_config_value ('data_path_required_folders')
raw_data_file_ext = replace_cfg_with_environment('raw_data_file_ext', get_config_value ('raw_data_file_ext'))
# leafcutter_library_strandedness = replace_cfg_with_environment('leafcutter_library_strandedness', get_config_value ('leafcutter_library_strandedness'))
# leafcutter_contrast_db_mapping_code = replace_cfg_with_environment('leafcutter_contrast_db_mapping_code', get_config_value ('leafcutter_contrast_db_mapping_code'))
# leafcutter_contrast_group_count_min = replace_cfg_with_environment('leafcutter_contrast_group_count_min', get_config_value ('leafcutter_contrast_group_count_min'))

# get user defined study_id and center_id and validate the values
metadata_db_study_id = replace_cfg_with_environment('metadata_db_study_id', get_config_value ('metadata_db_study_id'))
if not isinstance(metadata_db_study_id, int):
    if not metadata_db_study_id is None:
        add_pre_run_warning('Warning: Not int value ({}) was provided for "metadata_db_study_id" and it will be ignored'.format(metadata_db_study_id))
    metadata_db_study_id = None
metadata_db_center_id = replace_cfg_with_environment('metadata_db_center_id', get_config_value ('metadata_db_center_id'))
if not isinstance(metadata_db_center_id, int):
    if not metadata_db_center_id is None:
        add_pre_run_warning('Warning: Not int value ({}) was provided for "metadata_db_study_id" and it will be ignored'.format(metadata_db_center_id))
    metadata_db_center_id = None

conda_envs_path=replace_cfg_with_environment ('conda_envs_path', get_config_value ('conda_envs_path'))

samplesToIgnore_str = replace_cfg_with_environment ('samplesToIgnore', get_config_value ('samplesToIgnore'))
samplesToIgnore = get_samples_to_ignore (samplesToIgnore_str)

# file names for the global messages (info and warnings)
pipeline_info_file_path = str(Path(data_path) / replace_cfg_with_environment('pipeline_info_file', get_config_value ('pipeline_info_file')))
pipeline_warning_file_path = str(Path(data_path) / replace_cfg_with_environment('pipeline_warning_file', get_config_value ('pipeline_warning_file')))

add_pre_run_info ('Pipeline pre-run started at {}'.format(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))

add_pre_run_info ('help = {} (set to True to display the pipeline help info)'.format(help))
add_pre_run_info ('debug = {} (set to True to preserve "temp" objects)'.format(debug))
# add_pre_run_info ('run_alternative_splicing = {} (set to False to skip Alternative Splicing related rules. Default: True)'.format(run_alternative_splicing))
# if not run_alternative_splicing:
#     add_pre_run_warning ('Warning: Alternative Splicing processing was turned off based on the value supplied for "run_alternative_splicing" parameter. All tools considered part of the Alternative Splicing will be turned off despite their dedicated on/off parameters.')
# add_pre_run_info ('run_rsem = {} (set to False to skip RSEM related rules. Default: True)'.format(run_rsem))
# if not run_rsem:
#     add_pre_run_warning ('Warning: Processing RSEM rule was turned off based on the value supplied for "run_rsem" parameter.')
# add_pre_run_info ('run_rmats_turbo = {} (set to True to use rmats-turbo tool. Otherwise original rmats tool will be used. Default: False)'.format(run_rmats_turbo))
# add_pre_run_info ('run_rmats_novel = {} (set to True to supply "--novelSS" argument to the rmats command line. Otherwise the argument will be skipped. Default: False)'.format(run_rmats_novel))

# add_pre_run_info ('run_leafcutter = {} (set to True to use leafcutter tool. Default: False)'.format(run_leafcutter))
# add_pre_run_info ('run_spladder = {} (set to False to skip Spladder tool related steps. Default: True)'.format(run_spladder))
# if not run_spladder:
#     add_pre_run_warning ('Warning: Processing Spladder related rules was turned off based on the value supplied for "run_spladder" parameter.')

# add_pre_run_info ('run_spladder_long_execution_events = {} (If set to True, additional events will be used to create the spladder output, but it can take a prolonged time. See the list of these events in the "spladder_events_long_execution" variable below. Default: False)'.format(run_spladder_long_execution_events))

# add_pre_run_info ('run_djexpress = {} (set to False to skip DJExpress tool. Default: True. This parameter will be ignored if run_alternative_splicing is set to False)'.format(run_djexpress))
# if not run_djexpress:
#    add_pre_run_warning ('Warning: Processing DJExpress related rules was turned off based on the value supplied for "run_djexpress" parameter.')

# add_pre_run_info ('run_ngscheckmate = {} (set to False to skip NGSCheckMate related rules. Default: True)'.format(run_ngscheckmate))
# if not run_ngscheckmate:
#     add_pre_run_warning ('Warning: Processing NGSCheckMate related rules was turned off based on the value supplied for "run_ngscheckmate" parameter.')

# add_pre_run_info ('run_feature_counts = {} (set to True to execute the featureCounts tool. Default: False)'.format(run_feature_counts))
# if not run_feature_counts:
#     add_pre_run_warning ('Warning: Processing FeatureCount related rules was turned off based on the value supplied for "run_feature_counts" parameter.')
    
# add_pre_run_info ('run_telescope = {} (set to True to execute the telescope related tools. Default: True)'.format(run_telescope))
# if not run_telescope:
#     add_pre_run_warning ('Warning: Processing telescope related rules was turned off based on the value supplied for "run_telescope" parameter.')

add_pre_run_info ('create_rds_for_rscripts = {} (if set to True, RDS files (being used for manual testing) will be created for all participating R scripts that are executed at the current run)'.format(create_rds_for_rscripts))
add_pre_run_info ('data_path = {}'.format(data_path))
add_pre_run_info ('data_path_required_folders = {}'.format(data_path_required_folders))
add_pre_run_info ('raw_data_file_ext = {}'.format(raw_data_file_ext))
# add_pre_run_info ('leafcutter_library_strandedness = {} (possible values "XS", "FR"; used for "intronMotifBAM2junc" and "leafcutter_cluster_regtools" rules)'.format(leafcutter_library_strandedness))

# add_pre_run_info ('leafcutter_contrast_db_mapping_code = {} (String value is expected. Mapping code assigned to pull data from the database; used for "leafcutter_contrast_file" rules)'.format(leafcutter_contrast_db_mapping_code))
# add_pre_run_info ('leafcutter_contrast_group_count_min = {} (Integer value is expected. Defines minimum count of rows in a contrast group; used for "leafcutter_contrast_file" rules)'.format(leafcutter_contrast_group_count_min))
add_pre_run_info ('metadata_db_study_id = {} (Integer value is expected. Represents study id to be used from the database. If it is provided, it will overwrite a default value associated with an aliquot.)'.format(metadata_db_study_id))
add_pre_run_info ('metadata_db_center_id = {} (Integer value is expected. Represents center id to be used from the database. If it is provided, it will overwrite a default value associated with an aliquot)'.format(metadata_db_center_id))

add_pre_run_info ('snakemake_path = {}'.format(snakemake_path))
add_pre_run_info ('conda_envs_path = {}'.format(conda_envs_path))

# load email related variables for sending status email on completion
smtp_server = replace_cfg_with_environment('pipeline_smtp_server', None)
smtp_server_port = replace_cfg_with_environment('pipeline_smtp_server_port', None)
email_from = replace_cfg_with_environment('pipeline_email_from', None)
emails_to = replace_cfg_with_environment('pipeline_emails_to', None)  # it is expected to be a comma delimited list of emails or a list object of email addresses
email_subject_template = replace_cfg_with_environment('pipeline_email_subject', get_config_value ('email/email_subject'))
email_attachment_max_size = replace_cfg_with_environment('pipeline_email_attachment_max_size', get_config_value ('email/email_attachment_max_size'))
if email_attachment_max_size:
    email_attachment_max_size_bytes = 1024 * 1024 * email_attachment_max_size
email_body_success = get_config_value ('email/email_body_success_template')
email_body_error = get_config_value ('email/email_body_error_template')

if emails_to:
    add_pre_run_info ('emails_to = {}'.format(emails_to))
else:
    add_pre_run_info ('emails_to = {}; this value can be set through "pipeline_emails_to" environment variable'.format(emails_to))

if not data_path is None:
    data_path = str(os.path.realpath(Path(data_path)))
    if os.path.isdir(data_path):
        if data_path_required_folders and isinstance(data_path_required_folders, list):
            for req_dir in data_path_required_folders:
                if not os.path.isdir(str(os.path.realpath(Path(data_path) / req_dir))):
                    abort_execution = True
                    add_pre_run_warning  ('ERROR: The "data_path" directory is expected to have "{}" sub-directory which does not exists or not accessible.'.format(req_dir))
    else:
        abort_execution = True
        add_pre_run_warning  ('ERROR: The "data_path" parameter is not a directory.')

if abort_execution:
    email_subject = email_subject_template.format('with ERRORS')
    cur_warnings = get_pre_run_warnings()
    add_pre_run_info('Workflow finished with errors. Attempting to send a status email to: {}.'.format(emails_to))
    
    email_body = email_body_error.format('', data_path, 'N/A', email_attachment_max_size) + '\n\nProduced Warnings:\n' + cur_warnings
    
    send_yagmail( \
        emails_to = emails_to, \
        subject = email_subject, \
        message = email_body, \
        email_from=email_from, \
        attachments=None, \
        smtp_server=smtp_server, \
        smtp_server_port=smtp_server_port, \
        attachment_max_bytes = email_attachment_max_size_bytes \
        )
        
    print_pre_run_info()
    print_pre_run_warnings()
    sys.exit(1)

# global ref data path
ref_data_path = replace_cfg_with_environment('ref_data_path', get_config_value ('ref_data_path'))
add_pre_run_info ('ref_data_path = {}'.format(ref_data_path))

motrpac_ref_data = replace_cfg_with_environment('motrpac_ref_data', get_config_value ('motrpac_ref_data'))
add_pre_run_info ('motrpac_ref_data = {}'.format(motrpac_ref_data))

# motrpac ref data as per config file
motrpac_ref_data_path = str(os.path.realpath(Path(ref_data_path) / motrpac_ref_data))
# if envir variable was provided, overwrite the config value with that
motrpac_ref_data_path = replace_cfg_with_environment('motrpac_ref_data_path', motrpac_ref_data_path)
add_pre_run_info ('motrpac_ref_data_path = {}'.format(motrpac_ref_data_path))
motrpac_ref_data_genom_dir = replace_cfg_with_environment('motrpac_ref_data_genom_dir', get_config_value ('motrpac_ref_data_genom_dir'))
add_pre_run_info ('motrpac_ref_data_genom_dir = {}'.format(motrpac_ref_data_genom_dir))

#alternative slicing ref data based on the config
alt_sl_ref_data_path = str(os.path.realpath(Path(ref_data_path) / get_config_value ('alt_sl_ref_data')))
# if envir variable was provided, overwrite the config value with that
alt_sl_ref_data_path = replace_cfg_with_environment('alt_sl_ref_data_path', alt_sl_ref_data_path)
add_pre_run_info ('alt_sl_ref_data_path = {}'.format(alt_sl_ref_data_path))

alt_sl_ref_data_genom_dir = replace_cfg_with_environment('alt_sl_ref_data_genom_dir', get_config_value ('alt_sl_ref_data_genom_dir')) 
if alt_sl_ref_data_genom_dir is None or len(alt_sl_ref_data_genom_dir) == 0:
    # if no config/envir value is provided, use the name given for the motrpac_ref_data_genom_dir
    alt_sl_ref_data_genom_dir = motrpac_ref_data_genom_dir
add_pre_run_info ('alt_sl_ref_data_genom_dir = {}'.format(alt_sl_ref_data_genom_dir))

# verify existens of main reference data paths
if not ref_data_path is None:
    ref_data_path = str(os.path.realpath(Path(ref_data_path)))
    if not os.path.isdir(ref_data_path):
        abort_execution = True
        add_pre_run_warning  ('ERROR: The "ref_data_path" parameter is not a directory.')
if not motrpac_ref_data_path is None:
    motrpac_ref_data_path = str(os.path.realpath(Path(motrpac_ref_data_path)))
    if not os.path.isdir(motrpac_ref_data_path):
        abort_execution = True
        add_pre_run_warning  ('ERROR: The "motrpac_ref_data_path" parameter is not a directory.')
if not alt_sl_ref_data_path is None:
    alt_sl_ref_data_path = str(os.path.realpath(Path(alt_sl_ref_data_path)))
    if not os.path.isdir(alt_sl_ref_data_path):
        abort_execution = True
        add_pre_run_warning  ('ERROR: The "alt_sl_ref_data_path" parameter is not a directory.')

if abort_execution:
    email_subject = email_subject_template.format('with ERRORS')
    cur_warnings = get_pre_run_warnings()
    add_pre_run_info('Workflow finished with errors. Attempting to send a status email to: {}.'.format(emails_to))
    
    email_body = email_body_error.format('', data_path, 'N/A', email_attachment_max_size) + '\n\nProduced Warnings:\n' + cur_warnings
    
    send_yagmail( \
        emails_to = emails_to, \
        subject = email_subject, \
        message = email_body, \
        email_from=email_from, \
        attachments=None, \
        smtp_server=smtp_server, \
        smtp_server_port=smtp_server_port, \
        attachment_max_bytes = email_attachment_max_size_bytes \
        )
        
    print_pre_run_info()
    print_pre_run_warnings()
    sys.exit(1)

# check if kallisto_genom_dir has a special assignment, if not use the alt_sl_ref_data_genom_dir
# kallisto_genom_dir = replace_cfg_with_environment('kallisto_genom_dir', get_config_value ('rules/kallisto/genom_dir'))
# if kallisto_genom_dir is None or len(kallisto_genom_dir) == 0:
#     kallisto_genom_dir = alt_sl_ref_data_genom_dir
# add_pre_run_info ('kallisto_genom_dir = {}'.format(kallisto_genom_dir))

# check if suppa_genom_dir has a special assignment, if not use the alt_sl_ref_data_genom_dir
# suppa_genom_dir = replace_cfg_with_environment('suppa_genom_dir', get_config_value ('rules/suppa_kallisto/genom_dir'))
# if suppa_genom_dir is None or len(suppa_genom_dir) == 0:
#     suppa_genom_dir = alt_sl_ref_data_genom_dir
# add_pre_run_info ('suppa_genom_dir = {}'.format(suppa_genom_dir))

# suppa_events = get_config_value ('suppa_events')
# add_pre_run_info ('suppa_events = {}'.format(suppa_events))
# rmats_events = get_config_value ('rmats_events')
# add_pre_run_info ('rmats_events = {}'.format(rmats_events))

# suppa_file_ext = get_config_value ('suppa_file_ext') 
# rmats_file_ext = get_config_value ('rmats_file_ext') 
# add_pre_run_info ('suppa_file_ext = {}'.format(suppa_file_ext))
# add_pre_run_info ('rmats_file_ext = {}'.format(rmats_file_ext))

## rmats_events_compiled_file_ext = get_config_value('rmats_events_compiled_file_ext')
## suppa_events_compiled_file_ext = get_config_value('suppa_events_compiled_file_ext')

# spladder_events = get_config_value ('spladder_events')
# add_pre_run_info ('spladder_events = {}'.format(spladder_events))

# spladder_events_long_execution = get_config_value ('spladder_events_long_execution')
# add_pre_run_info ('spladder_events_long_execution = {}'.format(spladder_events_long_execution))
# if run_spladder and not run_spladder_long_execution_events and spladder_events_long_execution:
#     add_pre_run_warning ('Warning: The following Spladder events will not be used during this run, based on the values of the "run_spladder_long_execution_events" parameter: {}'.format(spladder_events_long_execution))


file_types = get_config_value ('file_types')
file_type_groups = get_config_value ('file_type_groups')
add_pre_run_info ('file_types = {}'.format(file_types))
add_pre_run_info ('file_type_groups = {}'.format(file_type_groups))

samples_by_type = {}
samples_by_id = {}
all_samples = []
sample_groups = {}

# prefix for the alternative groups that contains only R1 samples 
# (i.e. for the R1R2 group it will contain R1 samlpes that have R1 files but not the R2 files)
alt_R1_grp_prefix = 'R1_not_in_'  

for ft in file_types:
    samples_by_type[ft['name']], = glob_wildcards(data_path+"/fastq_raw/{sample,[^/]+}_" + ft['name'] + raw_data_file_ext)  # ".fastq.gz"
    samples_by_type[ft['name']].sort()
    all_samples.extend(samples_by_type[ft['name']])

all_samples = list(set(all_samples))  # to remove duplicated aliquot ids (that come from different file types)

# create a dictionary with aliquot id as a key and list of present file types
# verify that all required samples are present
for smpl in all_samples:
    samples_by_id[smpl] = []
    for ft in file_types:
        if smpl in samples_by_type[ft['name']]:
            samples_by_id[smpl].append(ft['name'])
        else:
            if ft['required']:
                add_pre_run_warning ('Warning: Requied file type "{}" for aliquot id "{}" is not present.'.format(ft['name'], smpl))

# create a list of group samples
for smpl in all_samples:
    for ftg in file_type_groups.keys():
        if not ftg in sample_groups:
            sample_groups[ftg] = []
        qualify = False
        for t in file_type_groups[ftg]:
            if smpl in samples_by_type[t]:
                qualify = True
            else:
                qualify = False
                break
        if qualify:
            sample_groups[ftg].append(smpl)
        else:
            if smpl in samples_by_type['R1']:
                grp = alt_R1_grp_prefix + ftg
                if not grp in sample_groups:
                    sample_groups[grp] = []
                sample_groups[grp].append(smpl)

add_pre_run_info ('sample_groups: {}'.format(sample_groups))

# check if each R1 file has a corresponding R2 file. Based on that update samples_to_process and samples_R1R2_present variables. 
# these variables will define how some of the rules will be used (i.e. trim rule)
if not alt_R1_grp_prefix + 'R1R2' in sample_groups or not sample_groups[alt_R1_grp_prefix + 'R1R2']:
    # all R1 sampels have a corresponding R2 sample
    samples_to_process = sample_groups['R1R2']
    samples_R1R2_present = True
else:
    # not all R1 sampels have a corresponding R2 sample
    samples_to_process = samples_by_type['R1']
    samples_R1R2_present = False

# check if each R1 file has a corresponding I1 file. Based on that update samples_R1I1_present variable. 
# this variable will define how some of the rules will be used (i.e. UMI_attach rule)
if not alt_R1_grp_prefix + 'R1I1' in sample_groups or not sample_groups[alt_R1_grp_prefix + 'R1I1']:
    # all R1 sampels have a corresponding I1 sample
    samples_R1I1_present = True
else:
    # not all R1 sampels have a corresponding I1 sample
    samples_R1I1_present = False

samplesToIgnore_notMatched = []
# remove any samples that were requested (in samplesToIgnore) to be removed
samples_reduced = list(set(samples_to_process).intersection(samplesToIgnore))  # determines if any elements are common between samples_to_process and samplesToIgnore lists
if len(samples_reduced) < len(samplesToIgnore):
    # some to ignore samlpes were not matched to the samples to process
    samplesToIgnore_notMatched = [s for s in samplesToIgnore if not s in samples_reduced]

samples_reduced_bool = bool(samples_reduced)  # bool, 
samples_to_process = list(set(samples_to_process).difference(samplesToIgnore)) # removes matching items of samplesToIgnore list from samples_to_process

samples_to_process.sort()  # sort samples to make sure aliquots come in the same order in case of rerun

add_pre_run_info('samples_by_type = {}'.format(samples_by_type))
add_pre_run_info('all_samples (distinct) count = {}'.format(len(all_samples)))
add_pre_run_info('all_samples (distinct) = {}'.format(all_samples))
add_pre_run_info('sample_groups = {}'.format(sample_groups))
add_pre_run_info('samplesToIgnore count = {}'.format(len(samplesToIgnore)))
add_pre_run_info('samplesToIgnore = {}'.format(samplesToIgnore))
add_pre_run_info('samples reduced by samplesToIgnore count = {} {}'.format( \
                len(samples_reduced), \
                ('(check warning section for not matched ignore samples!)' if samplesToIgnore_notMatched else '') \
                ))
add_pre_run_info('samples_to_process count = {}'.format(len(samples_to_process)))
add_pre_run_info('samples_to_process {} = {}' \
    .format(( \
        ('({}reduced by samplesToIgnore)'.format('' if samples_reduced_bool else 'not ')) \
        if samplesToIgnore else ''), \
        samples_to_process))
# will record warning if some samplesToIgnore were provided but some of them were not matched to the samples_to_process
if samplesToIgnore_notMatched:
    add_pre_run_warning('Warning: The following {} sample(s) from the ignore list (samplesToIgnore) were no found in the list of samlpes (samples_to_process): {}'.format(len(samplesToIgnore_notMatched), samplesToIgnore_notMatched))
add_pre_run_info('samples_R1R2_present = {}'.format (samples_R1R2_present))
add_pre_run_info('samples_R1I1_present = {}'.format (samples_R1I1_present))

# validate selected samples_to_process
if not len(samples_to_process) > 0:
    abort_execution = True
    add_pre_run_warning  ('ERROR: No samples were selected for processing.')

# loop throug all samples and report any not accessable files    
for smp in samples_to_process:
    fileR1 = get_R1_file_path(data_path, 'fastq_raw').replace('{sample}', smp)
    if not os.path.isfile(os.path.realpath(fileR1)):
        abort_execution = True
        add_pre_run_warning  ('ERROR: Sample: {}. The following requested R1 file does not exist or is not accessable; real path of the file: {})'.format(smp, os.path.realpath(fileR1)))
    
    if samples_R1R2_present:
        fileR2 = get_R2_file_path(data_path, 'fastq_raw').replace('{sample}', smp)
        if not os.path.isfile(os.path.realpath(fileR2)):
            abort_execution = True
            add_pre_run_warning  ('ERROR: Sample: {}. The following requested R2 file does not exist or is not accessable; real path of the file: {})'.format(smp, os.path.realpath(fileR2)))
    
if abort_execution:
    email_subject = email_subject_template.format('with ERRORS')
    cur_warnings = get_pre_run_warnings()
    add_pre_run_info('Workflow finished with errors. Attempting to send a status email to: {}.'.format(emails_to))
    
    str_samples = ' for {} samples'.format(len(samples_to_process))
    email_body = email_body_error.format(str_samples, data_path, 'N/A', email_attachment_max_size) + '\n\nProduced Warnings:\n' + cur_warnings
    
    send_yagmail( \
        emails_to = emails_to, \
        subject = email_subject, \
        message = email_body, \
        email_from=email_from, \
        attachments=None, \
        smtp_server=smtp_server, \
        smtp_server_port=smtp_server_port, \
        attachment_max_bytes = email_attachment_max_size_bytes \
        )
    print_pre_run_info()
    print_pre_run_warnings()
    sys.exit(1)    

# check for the rmats contrast metadata file
# rmats_contrast_file_expected_path = str(os.path.realpath(Path(data_path) / get_config_value ('rmats_contrast_file')))  # expected location of the metadata contrast file
# rmats_contrast_file = replace_cfg_with_environment('rmats_contrast_file', get_config_value ('rmats_contrast_file')) 
# rmats_default_contrast = replace_cfg_with_environment('rmats_default_contrast', get_config_value ('rmats_default_contrast'))
# rmats_min_b1_qty = replace_cfg_with_environment('rmats_min_b1_qty', get_config_value ('rmats_min_b1_qty'))
# rmats_min_b2_qty = replace_cfg_with_environment('rmats_min_b2_qty', get_config_value ('rmats_min_b2_qty'))

# get rmats_contrast_file_path
# if os.path.isabs(rmats_contrast_file):
#     rmats_contrast_file_path = rmats_contrast_file
# else:
#     rmats_contrast_file_path = str(os.path.realpath(Path(data_path) / rmats_contrast_file))
# add_pre_run_info ('rmats_contrast_file_path provided from config or environment  = {}'.format (rmats_contrast_file_path))

"""
# if run_alternative_splicing:
    #validate the given path for the rmats contrast file. If the file is not present, create it using the default contrast assignment for all samples
    rmats_contrast_file_valid, meta_file_created = validate_rmats_contrast_file (rmats_contrast_file_path, rmats_contrast_file_expected_path, samples_to_process, 
                                                                rmats_default_contrast, rmats_min_b1_qty, rmats_min_b2_qty)
    rmats_contrast_file_path = rmats_contrast_file_expected_path

    add_pre_run_info ('rmats_contrast_file_path = {}'.format (rmats_contrast_file_path))
    add_pre_run_info ('rmats_min_b1_qty = {}'.format (rmats_min_b1_qty))
    add_pre_run_info ('rmats_min_b2_qty = {}'.format (rmats_min_b2_qty))
    add_pre_run_info ('rmats_default_contrast = {}'.format (rmats_default_contrast))
    add_pre_run_info ('rmats_contrast_file_valid = {}'.format (rmats_contrast_file_valid))
    if not rmats_contrast_file_valid:
        add_pre_run_warning ('Warning: Provided rmats contrast file is not valid. All rmats related rules will be skipped.')
    add_pre_run_info ('rmats_contrast_file auto created = {}'.format (meta_file_created))
else:
    rmats_contrast_file_valid = None
"""
"""
if run_spladder:
    # define chunk size related information for processing spladder_merge steps utilizing chunking approach
    spladder_chunk_size = replace_cfg_with_environment('spladder_chunk_size', get_config_value ('spladder_chunk_size'))
    spladder_chunk_size = convert_to_int(spladder_chunk_size)  # convert to int, if not convertable, None will be returned
    if spladder_chunk_size is None:
        spladder_chunk_size = 2
    spl_merge_chunks = get_spladder_merge_chanks(len(samples_to_process), spladder_chunk_size)
    add_pre_run_info ('spladder_chunk_size = {}'.format (spladder_chunk_size))
    add_pre_run_info ('spl_merge_chunks = {}'.format (spl_merge_chunks))
    
    # defines ase_edge_limit value to be used with spladder events rules
    spladder_events_ase_edge_limit = replace_cfg_with_environment('spladder_events_ase_edge_limit', get_config_value ('spladder_events_ase_edge_limit'))
    spladder_events_ase_edge_limit_min = get_config_value ('spladder_events_ase_edge_limit_min')
    spladder_events_ase_edge_limit_decrease_step = replace_cfg_with_environment('spladder_events_ase_edge_limit_decrease_step', get_config_value ('spladder_events_ase_edge_limit_decrease_step'))
    spladder_events_ase_edge_limit = convert_to_int(spladder_events_ase_edge_limit)  # convert to int, if not convertable, None will be returned
    if spladder_events_ase_edge_limit is None:
        spladder_events_ase_edge_limit = 500
    spladder_events_ase_edge_limit_decrease_step = convert_to_int(spladder_events_ase_edge_limit_decrease_step)  # convert to int, if not convertable, None will be returned
    if spladder_events_ase_edge_limit_decrease_step is None:
        spladder_events_ase_edge_limit_decrease_step = 100
    if spladder_events_ase_edge_limit < spladder_events_ase_edge_limit_min:
        spladder_events_ase_edge_limit = spladder_events_ase_edge_limit_min
        add_pre_run_warning('Warning: Requested spladder_events_ase_edge_limit value ({}) is less than set minimum and will be replaced with the minimum value: {}'.format(spladder_events_ase_edge_limit, spladder_events_ase_edge_limit_min))
    add_pre_run_info ('spladder_events_ase_edge_limit = {}'.format (spladder_events_ase_edge_limit))
    add_pre_run_info ('spladder_events_ase_edge_limit_min = {}'.format (spladder_events_ase_edge_limit_min))
    add_pre_run_info ('spladder_events_ase_edge_limit_decrease_step = {} (This is a decriment value that be used to lower spladder_events_ase_edge_limit value ({}) in case of retrying attempts of the associated rules'.format (spladder_events_ase_edge_limit_decrease_step, spladder_events_ase_edge_limit))
else:
    spladder_events_ase_edge_limit = 500
    spladder_chunk_size = 0
    spl_merge_chunks = {}
"""

# display collected global info before the actual execution of the pipeline
# NOTE: the following 2 rows have to be commented to create a DAG diagram
if print_prerun_info:
    print_pre_run_info()
    print_pre_run_warnings()

onstart:
    # this will reset the info and warning files, if those exist from the previous run
    add_info (message = '', mode = 'w')
    add_warning (message = '', mode = 'w')  
    # record pre-run messages to the pipeline message files
    for msg in pre_run_messages['info']:
        add_info (msg)
    for msg in pre_run_messages['warning']:
        add_warning (msg)
    add_info ('\nPipeline run started at {}'.format(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))

onsuccess:
    log_path = list({log})[0]  # {log} is a built-in snakemake variable which contains the path to a logfile
    attachments = [log_path, pipeline_info_file_path, pipeline_warning_file_path]
    email_subject = email_subject_template.format('SUCCESSFULLY')
    cur_warnings = get_global_warnings()
    add_info('Workflow finished successfully at {}. Attempting to send a status email to: {}.'.format(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), emails_to))
    
    email_body = email_body_success.format(len(samples_to_process), data_path, log_path, email_attachment_max_size) + '\n\nProduced Warnings:\n' + cur_warnings
    
    send_yagmail( \
        emails_to = emails_to, \
        subject = email_subject, \
        message = email_body, \
        email_from=email_from, \
        attachments=attachments, \
        smtp_server=smtp_server, \
        smtp_server_port=smtp_server_port, \
        attachment_max_bytes = email_attachment_max_size_bytes \
        )
    
    print_global_info()
    print_global_warnings()

onerror:
    log_path = list({log})[0]  # {log} is a built-in snakemake variable which contains the path to a logfile
    attachments = [log_path, pipeline_info_file_path, pipeline_warning_file_path]
    email_subject = email_subject_template.format('with ERRORS')
    cur_warnings = get_global_warnings()
    add_info('Workflow finished with errors at {}. Attempting to send a status email to: {}.'.format(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), emails_to))
    
    str_samples = ' for {} samples'.format(len(samples_to_process))
    email_body = email_body_error.format(str_samples, data_path, log_path, email_attachment_max_size) + '\n\nProduced Warnings:\n' + cur_warnings
    
    send_yagmail( \
        emails_to = emails_to, \
        subject = email_subject, \
        message = email_body, \
        email_from=email_from, \
        attachments=attachments, \
        smtp_server=smtp_server, \
        smtp_server_port=smtp_server_port, \
        attachment_max_bytes = email_attachment_max_size_bytes \
        )
    
    print_global_info()
    print_global_warnings()

rule all:
    input:
        fastq_R1R2_validate = expand(str(Path(data_path) / "fastq_raw_validate/valid_{sample}"), sample=samples_to_process) if samples_R1R2_present else [],
        fastq_I1_validate = expand(str(Path(data_path) / "fastq_raw_validate/valid_{sample}_I1"), sample=samples_to_process) if samples_R1I1_present else [],
        # umi_attach_R1 = expand(str(Path(data_path) / "fastq_attach/{sample}_R1.fastq.gz"), sample=samples_to_process) if samples_R1I1_present else [],
        # umi_attach_R2 = expand(str(Path(data_path) / "fastq_attach/{sample}_R2.fastq.gz"), sample=samples_to_process) if samples_R1R2_present and samples_R1I1_present else [],
        # the following is commented to avoid reprocessing pipeline just because the trim files were deleted at the end of the process
        # trim = expand(get_R1_file_path(data_path, 'fastq_trim'), sample=samples_to_process),
        # star_align_bam = expand(str(Path(data_path) / "star_align/{sample}/Aligned.sortedByCoord.out.bam"), sample=samples_to_process),
        # star_align_validate = expand(str(Path(data_path) / "star_align_validate/valid_{sample}"), sample=samples_to_process),
        # star_align_qc_log = expand(str(Path(data_path) / "star_align/{sample}/Log.final.out"), sample=samples_to_process),
        # star_align_intronMotif_bam = expand(str(Path(data_path) / "star_align_intronMotif/{sample}/Aligned.sortedByCoord.out.bam"), sample=samples_to_process) if run_leafcutter else [],
        # star_align_intronMotif_qc_log = expand(str(Path(data_path) / "star_align_intronMotif/{sample}/Log.final.out"), sample=samples_to_process) if run_leafcutter else [],
        # star_align_intronMotif_bai = expand(str(Path(data_path) / "star_align_intronMotif/{sample}/Aligned.sortedByCoord.out.bam.bai"), sample=samples_to_process) if run_leafcutter else [],
        # intronMotif_junc = expand(str(Path(data_path) / "leafcutter/regtools/{sample}_Aligned.sortedByCoord.out.bam.junc"), sample=samples_to_process) if run_leafcutter else [],
        # juncfiles_txt = str(Path(data_path) / "leafcutter/regtools/_juncfiles.txt") if run_leafcutter else [],
        # leafcutter_cluster_regtools = str(Path(data_path) / ("leafcutter/cluster/test_contrast_perind_numers.counts.gz")) if run_leafcutter else [],
        # leafcutter_contrast = str(Path(data_path) / "metadata/leafcutter_contrast.txt") if run_leafcutter else [],
        # ds_cluster_significance = str(Path(data_path) / ("leafcutter/ds/leafcutter_ds_cluster_significance.txt")) if run_leafcutter else [],

        # bismark_ok_file = expand(str(Path(data_path) / "bismark/{sample}_ok.txt"), sample=samples_to_process),
        bismark_bam_out = expand(str(Path(data_path) / "bismark/{sample}_R1_bismark_bt2_pe.bam"), sample=samples_to_process) if samples_R1R2_present else expand(str(Path(data_path) / "bismark/{sample}_R1_bismark_bt2.bam"), sample=samples_to_process),
        # bismark_sort_ok_file = expand(str(Path(data_path) / "bismark/{sample}_bismark_sort_ok.txt"), sample=samples_to_process),
        # this row will be in use only when I1 files are present
        bismark_pre_dup_bam_umi = expand(str(Path(data_path) / "bismark/{sample}_R1_bismark_bt2_pe_umi.bam") if samples_R1R2_present else str(Path(data_path) / "bismark/{sample}_R1_bismark_bt2_umi.bam"), sample=samples_to_process) if samples_R1I1_present else [],
        bismark_dedup_txt = expand(str(Path(data_path) / "bismark/{sample}_dedup.txt"), sample=samples_to_process),
        # bismark_methylation_extractor_ok_file = expand(str(Path(data_path) / "bismark/{sample}_bismark_methylation_extractor_ok.txt"), sample=samples_to_process),
        bismark_methylation_extractor_CHG_context = \
            expand( \
            str(Path(data_path) / "bismark/CHG_context_{sample}_R1_bismark_bt2_pe.deduplicated.txt") \
            if samples_R1R2_present else 
            str(Path(data_path) / "bismark/CHG_context_{sample}_R1_bismark_bt2.deduplicated.txt"), \
            sample=samples_to_process),
        # bismark_report_ok_file = expand(str(Path(data_path) / "bismark/{sample}_bismark_report_ok.txt"), sample=samples_to_process),
        bismark_report_html = expand(str(Path(data_path) / "bismark/{sample}_report.html"), sample=samples_to_process),
        
        # spladder_single_graphs_count_file = expand(str(Path(data_path) / "spladder/spladder/genes_graph_conf3.{sample}_Aligned.sortedByCoord.out.count.hdf5"), sample=samples_to_process) if run_spladder else [],
        # spladder_alignments = str(Path(data_path) / "spladder/alignments.txt") if run_spladder else [],
        # spladder_merged_graph_chunck_top_level = expand(str(Path(data_path) / "spladder/status_files/spladder_merged_graph_chunk_{chunk}.txt"), chunk=[spl_merge_chunks['level'][len(spl_merge_chunks['level'])][0]]) if run_spladder else [],
        # spladder_quantification1_count_file = expand(str(Path(data_path) / "spladder/spladder/genes_graph_conf3.merge_graphs.{sample}_Aligned.sortedByCoord.out.count.hdf5"), sample=samples_to_process) if run_spladder else [],
        # spladder_quantification2_status_file = str(Path(data_path) / "spladder/status_files/spladder_quantification2.txt") if run_spladder else [],
        # spladder_call_events_merge_graphs_pickle = expand(str(Path(data_path) / "spladder/merge_graphs_{event}_C3.pickle"), event=spladder_events) if run_spladder else [],
        # spladder_call_events_long_execution_merge_graphs_pickle = expand(str(Path(data_path) / "spladder/merge_graphs_{event}_C3.pickle"), event=spladder_events_long_execution) if run_spladder and run_spladder_long_execution_events else [],
        
        # star_align_djexpress_sj_file=expand(str(Path(data_path) / "star_align_djexpress/{sample}/SJ.out.tab"), sample=samples_to_process) if run_djexpress and run_alternative_splicing else [],
        # star_align_djexpress_named_junc_file=expand(str(Path(data_path) / "star_align_djexpress/junct_quant/{sample}_SJ.out.tab"), sample=samples_to_process) if run_djexpress and run_alternative_splicing else [],
        
        fastqc_raw = expand(get_R1_file_path(data_path, sub_dir_path = 'fastqc_raw', file_name_format = '{sample}_|FN|_fastqc.html'), sample=samples_to_process),
        # fastqc_trim = expand(get_R1_file_path(data_path, sub_dir_path = 'fastqc_trim', file_name_format = '{sample}_|FN|_fastqc.html'), sample=samples_to_process),
        # chr_info = expand(str(Path(data_path) / "star_align/{sample}/chr_info.txt"), sample=samples_to_process),
        # rRNA = expand(str(Path(data_path) / "rRNA/{sample}.txt"), sample=samples_to_process),
        # phix = expand(str(Path(data_path) / "phix/{sample}.txt"), sample=samples_to_process),
        # globin = expand(str(Path(data_path) / "globin/{sample}.txt"), sample=samples_to_process),
        
        # kallisto_abd_file1 = expand(str(Path(data_path) / "kallisto/{sample}/abundance.tsv"), sample=samples_to_process) if rmats_contrast_file_valid and run_alternative_splicing else [],
        # kallisto_abd_file2 = expand(str(Path(data_path) / "kallisto/{sample}/abundance.h5"), sample=samples_to_process) if rmats_contrast_file_valid and run_alternative_splicing else [],
        # kallisto_abd_tpm_file = expand(str(Path(data_path) / "kallisto/{sample}/abundance.tsv_tpm"), sample=samples_to_process) if rmats_contrast_file_valid and run_alternative_splicing else [],
        # qc53 = expand(str(Path(data_path) / "qc53/{sample}.RNA_Metrics"), sample=samples_to_process),
        # mark_dup_metrics = expand(str(Path(data_path) / "mark_dup/{sample}.dup_metrics"), sample=samples_to_process),
        # star_align_merge_all = str(Path(data_path) / "star_align/star_QC.txt"),
        # star_align_most_common_read_length = str(Path(data_path) / "star_align/most_common_read_length.txt"),
        
        multiqc_fastqc_raw_html = str(Path(data_path) / "multiqc/fastqc_raw.html"),
        
        # multiqc_fastqc_trim_html = str(Path(data_path) / "multiqc/fastqc_trim.html"),
        # multiqc_post_align = str(Path(data_path) / "multiqc/post_align.html"),
        # post_align_qc53 = str(Path(data_path) / "qc53.txt"),
        # suppa_kallisto = expand(str(Path(data_path) / "suppa/kallisto/{sample}/{sample}_{event}.psi"), sample=samples_to_process, event=suppa_events) if rmats_contrast_file_valid and run_alternative_splicing else [],
        
        # rmats_b1_file = str(Path(data_path) / ("rmats/b1.txt")),
        # rmats_b2_file = str(Path(data_path) / ("rmats/b2.txt")),
        # include rmats rule only if rmats contrast metadata file is valid
        # rmats_fromGTF=expand(str(Path(data_path) / "rmats/fromGTF.{event}.txt"), event=rmats_events) if rmats_contrast_file_valid and run_alternative_splicing else [] ,
        # rmats_jc_raw = expand(str(Path(data_path) / "rmats/JC.raw.input.{event}.txt"), event=rmats_events) if rmats_contrast_file_valid and run_alternative_splicing else [],
        # the following row implements original Frank's assignments of pkl and txt extension based on the event type, not in use as of now
        # rmats_compiled = expand(str(Path(data_path) / "compiled/jct_{event}.{ext}"), zip, event=rmats_events, ext=rmats_events_compiled_file_ext) if rmats_contrast_file_valid and run_alternative_splicing else [],
        # rmats_compiled= expand(str(Path(data_path) / "compiled/jct_{event}.{ext}"), event=rmats_events, ext=rmats_file_ext) if rmats_contrast_file_valid and run_alternative_splicing else [],
        # rmats_compiled_includes = expand(str(Path(data_path) / "compiled/includes/jct_{event}.Includes.{ext}"), event=rmats_events, ext=rmats_file_ext) if rmats_contrast_file_valid and run_alternative_splicing else [],
        # rmats_compiled_skips = expand(str(Path(data_path) / "compiled/skips/jct_{event}.Skips.{ext}"), event=rmats_events, ext=rmats_file_ext) if rmats_contrast_file_valid and run_alternative_splicing else [],
        
        # rmats_compiled_psi= expand(str(Path(data_path) / "compiled/psi/jct_{event}.PSItable.{ext}"), event=rmats_events, ext=rmats_file_ext) if rmats_contrast_file_valid and run_alternative_splicing else [],
        # rmats_compiled_psi_combined = str(Path(data_path) / "compiled/psi/jct_allEvents.txt") if rmats_contrast_file_valid and run_alternative_splicing else [],
        
        # the following row implements original Frank's assignments of pkl and txt extension based on the event type, not in use as of now
        # suppa_compiled = expand(str(Path(data_path) / "compiled/txr_{event}.{ext}"), zip, event=suppa_events, ext=suppa_events_compiled_file_ext),
        # suppa_compiled = expand(str(Path(data_path) / "compiled/txr_{event}.{ext}"), event=suppa_events, ext=suppa_file_ext) if rmats_contrast_file_valid and run_alternative_splicing else [],
        # suppa_compiled_includes = expand(str(Path(data_path) / "compiled/includes/txr_{event}.Includes.{ext}"), event=suppa_events, ext=suppa_file_ext) if rmats_contrast_file_valid and run_alternative_splicing else [],
        # suppa_compiled_skips = expand(str(Path(data_path) / "compiled/skips/txr_{event}.Skips.{ext}"), event=suppa_events, ext=suppa_file_ext) if rmats_contrast_file_valid and run_alternative_splicing else [],
        
        # suppa_compiled_psi = expand(str(Path(data_path) / "compiled/psi/txr_{event}.PSItable.{ext}"), event=suppa_events, ext=suppa_file_ext) if rmats_contrast_file_valid and run_alternative_splicing else [],
        # suppa_compiled_psi_combined = str(Path(data_path) / "compiled/psi/txr_allEvents.txt") if rmats_contrast_file_valid and run_alternative_splicing else [],
        
        # rsem_genes = expand(str(Path(data_path) / "rsem/{sample}.genes.results"), sample=samples_to_process) if run_rsem else [],
        # rsem_isoforms = expand(str(Path(data_path) / "rsem/{sample}.isoforms.results"), sample=samples_to_process) if run_rsem else [],
        # suppa_rsem = expand(str(Path(data_path) / "suppa/rsem/{sample}/{sample}_{event}.psi"), sample=samples_to_process, event=suppa_events) if run_rsem and run_alternative_splicing else [],
        # rsem_genes_count = str(Path(data_path) / "rsem_genes_count.txt") if run_rsem else [],
        # rsem_genes_tpm = str(Path(data_path) / "rsem_genes_tpm.txt") if run_rsem else [],
        # rsem_genes_fpkm = str(Path(data_path) / "rsem_genes_fpkm.txt") if run_rsem else [],
        # featureCounts = expand(str(Path(data_path) / "featureCounts/{sample}"), sample=samples_to_process) if run_feature_counts else [],
        # featureCounts_merge_all = str(Path(data_path) / "featureCounts.txt") if run_feature_counts else [],
        # umi_dup = expand(str(Path(data_path) / "star_align/{sample}/dup_log.txt"), sample=samples_to_process) if samples_R1I1_present else [],
        # qc_final = str(Path(data_path) / "qc_info.csv"),
        # NGSCheckMate_prep = expand(str(Path(data_path) / "NGSCheckMate/bam/{sample}.bam"), sample=samples_to_process) if run_ngscheckmate else [],
        # NGSCheckMate_vcf = expand(str(Path(data_path) / "NGSCheckMate/vcf/{sample}.vcf"), sample=samples_to_process) if run_ngscheckmate else [],
        # NGSCheckMate_matched = str(Path(data_path) / "NGSCheckMate/output_matched.txt") if run_ngscheckmate else [],
        # NGSCheckMate_matrix = str(Path(data_path) / "NGSCheckMate/output_corr_matrix.txt") if run_ngscheckmate else [],
        # NGSCheckMate_all = str(Path(data_path) / "NGSCheckMate/output_all.txt") if run_ngscheckmate else [],
        # incoherent_groups_out = str(Path(data_path) / "NGSCheckMate/incoherent_groups_out.json") if run_ngscheckmate else [],
        

rule fastqc_raw:
    input:
        fileR1 = get_R1_file_path(data_path, 'fastq_raw'),
        fileR2 = get_R2_file_path(data_path, 'fastq_raw') if samples_R1R2_present else []
    output:
        # outR1 = str(Path(data_path) / "fastqc_raw/{sample}_R1_fastqc.html"),
        outR1 = get_R1_file_path(data_path, sub_dir_path = 'fastqc_raw', file_name_format = '{sample}_|FN|_fastqc.html'),
        # outR2 = str(Path(data_path) / "fastqc_raw/{sample}_R2_fastqc.html") if samples_R1R2_present else [],
        outR2 = get_R2_file_path(data_path, sub_dir_path = 'fastqc_raw', file_name_format = '{sample}_|FN|_fastqc.html') if samples_R1R2_present else [],
    conda:
        get_conda_env(),
    log:
        str(Path(data_path) / "logs/fastqc_raw/fastqc_raw_{sample}.log")
    resources:
        # mem = next((el for el in [get_config_value ('rules/fastqc_raw/memory')] if el is not None), default_rule_memory)  # 40000
        cl_resources = lambda wildcards, attempt : get_cluster_resources_by_attempt (attempt, 'fastqc_raw'),
        cl_job_suffix = lambda wildcards : wildcards.sample,
    params:
        tool = "fastqc",
        out_dir = str(Path(data_path) / "fastqc_raw")
    shell:
        '''
        {params.tool} --outdir {params.out_dir} {input.fileR1} {input.fileR2} 2>&1 | tee {log}
        '''
        
rule bismark:
    input:
        fileR1 = get_R1_file_path(data_path, 'fastq_raw'),
        fileR2 = get_R2_file_path(data_path, 'fastq_raw') if samples_R1R2_present else [],
        valid_file_R1R2 = str(Path(data_path) / "fastq_raw_validate/valid_{sample}") if samples_R1R2_present else [],
    output:
        # ok_file = str(Path(data_path) / "bismark/{sample}_ok.txt")
        """        
        Temporarily commemnted to avoid re-runing "bismark" step
        tmp_dir = temp_debugcheck(directory(str(Path(data_path) / "bismark/tmp_{sample}"))),
        bam = temp_debugcheck(str(Path(data_path) / "bismark/work/{sample}_R1_bismark_bt2_pe.bam") if samples_R1R2_present else str(Path(data_path) / "bismark/{sample}_R1_bismark_bt2.bam")),
        align_report= str(Path(data_path) / "bismark/work/{sample}_R1_bismark_bt2_PE_report.txt") if samples_R1R2_present else str(Path(data_path) / "bismark/{sample}_R1_bismark_bt2_SE_report.txt"),
        """
    conda:
        get_conda_env('bismark'),
    log:
        str(Path(data_path) / "logs/bismark/bismark_{sample}.log")
    threads: get_rule_threads ('bismark')
    resources:
        cl_resources = lambda wildcards, attempt : get_cluster_resources_by_attempt (attempt, 'bismark'),  # required for loading memory 
        cl_job_suffix = lambda wildcards : wildcards.sample,
        walltime = get_rule_walltime('bismark'),
    params:
        genomeDir = str(Path(motrpac_ref_data_path) / (motrpac_ref_data_genom_dir)),
        # out_dir = str(Path(data_path) / "bismark/work"),
        # tmp_dir = str(Path(data_path) / "bismark/tmp_{sample}"),
        file1_arg_name = "-1" if samples_R1R2_present else "",
        file2_arg_name = "-2" if samples_R1R2_present else "",
    shell:
        '''
        cd "$(dirname "{output.bam}")"  # get to the parent dir of the output file
        
        bismark \
        {params.genomeDir} \
        --non_directional \
        --multicore {threads} \
        --temp_dir {output.tmp_dir} \
        {params.file1_arg_name} {input.fileR1} {params.file2_arg_name} {input.fileR2} \
        2>&1 | tee {log}
        '''
        # cd {params.out_dir}
        
rule bismark_sort:
    input:
        bam = str(Path(data_path) / "bismark/work/{sample}_R1_bismark_bt2_pe.bam") if samples_R1R2_present else str(Path(data_path) / "bismark/work/{sample}_R1_bismark_bt2.bam"),
    output:
        tmp_dir = temp_debugcheck(directory(str(Path(data_path) / "bismark/tmp_{sample}"))),
        sam = temp_debugcheck(str(Path(data_path) / "bismark/work/{sample}.sam")),
        bam_out = temp_debugcheck(str(Path(data_path) / "bismark/{sample}_R1_bismark_bt2_pe.bam") if samples_R1R2_present else str(Path(data_path) / "bismark/{sample}_R1_bismark_bt2.bam")),
        
        # ok_file = str(Path(data_path) / "bismark/{sample}_bismark_sort_ok.txt"),
    conda:
        get_conda_env('bismark'),
    log:
        str(Path(data_path) / "logs/bismark_sort/bismark_sort_{sample}.log")
    threads: get_rule_threads ('bismark_sort')
    resources:
        cl_job_suffix = lambda wildcards : wildcards.sample,
    params:
        # tmp_dir = temp_debugcheck(directory(str(Path(data_path) / "bismark/tmp_{sample}"))),
        # sam = temp_debugcheck(str(Path(data_path) / "bismark/work/{sample}.sam")),
        # bam_out = str(Path(data_path) / "bismark/{sample}_R1_bismark_bt2_pe.bam") if samples_R1R2_present else str(Path(data_path) / "bismark/{sample}_R1_bismark_bt2.bam"),
    shell:
        '''
        mkdir -p {output.tmp_dir}
        
        samtools view -@ {threads} -H {input.bam} > {output.sam}  2>&1 | tee {log}
        samtools view -@ {threads} {input.bam} |sort -k5,5nr -k1,1 -s -S30G -T {output.tmp_dir} >> {output.sam} 2>&1 | tee -a {log}
        samtools view -@ {threads} -b {output.sam} -o {output.bam_out} 2>&1 | tee -a {log}
        '''
        

# this rule performed only if I1 files are present
rule bismark_umi_format:
    input: 
        bam_bismark = str(Path(data_path) / "bismark/{sample}_R1_bismark_bt2_pe.bam") if samples_R1R2_present else str(Path(data_path) / "bismark/{sample}_R1_bismark_bt2.bam"),
    output:
        pre_dup_bam_umi = \
            temp_debugcheck( \
            str(Path(data_path) / "bismark/umi/{sample}_R1_bismark_bt2_pe.bam") \
            if samples_R1R2_present else \
            str(Path(data_path) / "bismark/umi/{sample}_R1_bismark_bt2.bam")),
    conda:
        get_conda_env('bismark'),
    log:
        str(Path(data_path) / "logs/bismark_umi_format/bismark_umi_format_{sample}.log"),
    resources:
        cl_resources = lambda wildcards, attempt : get_cluster_resources_by_attempt (attempt, 'bismark_umi_format'),
        cl_job_suffix = lambda wildcards : wildcards.sample,
    params:
        sam = str(Path(data_path) / "bismark/{sample}_R1_bismark_bt2_pe.sam") if samples_R1R2_present else str(Path(data_path) / "bismark/{sample}_R1_bismark_bt2.sam"),
        umi_script = str(Path(snakemake_path) / "scripts/bismark_bam_UMI_format.awk"),
    shell:
        '''
        samtools view -h {input.bam_bismark} | {params.umi_script} > {params.sam} 2>&1 | tee {log}
        samtools view -b -o {output.pre_dup_bam_umi} {params.sam} 2>&1 | tee -a {log}
        '''

rule bismark_dedup:
    input: 
        # input vary based on the presense of the I1 and R2 files
        bam_bismark = \
            (str(Path(data_path) / "bismark/umi/{sample}_R1_bismark_bt2_pe.bam") if samples_R1R2_present else str(Path(data_path) / "bismark/umi/{sample}_R1_bismark_bt2.bam")) \
            if samples_R1I1_present else \
            (str(Path(data_path) / "bismark/{sample}_R1_bismark_bt2_pe.bam") if samples_R1R2_present else str(Path(data_path) / "bismark/{sample}_R1_bismark_bt2.bam")),
    output:
        #(str(Path(data_path) / "bismark/{sample}_R1_bismark_bt2_pe_umi.deduplicated.bam") if samples_R1R2_present else str(Path(data_path) / "bismark/{sample}_R1_bismark_bt2_umi.deduplicated.bam")) \
        #    if samples_R1I1_present else \
        dedup_bam = \
            str(Path(data_path) / "bismark/{sample}_R1_bismark_bt2_pe.deduplicated.bam") if samples_R1R2_present else str(Path(data_path) / "bismark/{sample}_R1_bismark_bt2.deduplicated.bam"),
        #(str(Path(data_path) / "bismark/{sample}_R1_bismark_bt2_pe_umi.deduplication_report.txt") if samples_R1R2_present else str(Path(data_path) / "bismark/{sample}_R1_bismark_bt2_umi.deduplication_report.txt")) \
        #    if samples_R1I1_present else \
        dedup_report = \
            (str(Path(data_path) / "bismark/{sample}_R1_bismark_bt2_pe.deduplication_report.txt") if samples_R1R2_present else str(Path(data_path) / "bismark/{sample}_R1_bismark_bt2.deduplication_report.txt")),
        dedup_txt = str(Path(data_path) / "bismark/{sample}_dedup.txt"),
    conda:
        get_conda_env('bismark'),
    log:
        str(Path(data_path) / "logs/bismark_dedup/bismark_dedup_{sample}.log"),
    resources:
        cl_resources = lambda wildcards, attempt : get_cluster_resources_by_attempt (attempt, 'bismark_dedup'),
        walltime = get_rule_walltime('bismark_dedup'),
        cl_job_suffix = lambda wildcards : wildcards.sample,
    params:
        pair_option = '-p' if samples_R1R2_present else '-s',
        barcode = '--barcode' if samples_R1I1_present else '',
    shell:
        '''
        cd "$(dirname "{output.dedup_txt}")"
        
        deduplicate_bismark {params.pair_option} {params.barcode} --bam {input.bam_bismark} 2>&1 | tee {log}
        # copy the created log file to the output directory since it will be used for collecting stats
        cp {log[0]} {output.dedup_txt}
        '''

rule bismark_methylation_extractor:
    input: 
        dedup_bam = \
            str(Path(data_path) / "bismark/{sample}_R1_bismark_bt2_pe.deduplicated.bam") if samples_R1R2_present else str(Path(data_path) / "bismark/{sample}_R1_bismark_bt2.deduplicated.bam"),
    output:
        CHG_context = \
            str(Path(data_path) / "bismark/CHG_context_{sample}_R1_bismark_bt2_pe.deduplicated.txt") \
            if samples_R1R2_present else 
            str(Path(data_path) / "bismark/CHG_context_{sample}_R1_bismark_bt2.deduplicated.txt"),
        CHH_context = \
            str(Path(data_path) / "bismark/CHH_context_{sample}_R1_bismark_bt2_pe.deduplicated.txt") \
            if samples_R1R2_present else 
            str(Path(data_path) / "bismark/CHH_context_{sample}_R1_bismark_bt2.deduplicated.txt"),
        CpG_context = \
            str(Path(data_path) / "bismark/CpG_context_{sample}_R1_bismark_bt2_pe.deduplicated.txt") \
            if samples_R1R2_present else 
            str(Path(data_path) / "bismark/CpG_context_{sample}_R1_bismark_bt2.deduplicated.txt"),
        M_bias = \
            str(Path(data_path) / "bismark/{sample}_R1_bismark_bt2_pe.deduplicated.M-bias.txt") \
            if samples_R1R2_present else 
            str(Path(data_path) / "bismark/{sample}_R1_bismark_bt2.deduplicated.M-bias.txt"),
        splitting_report = \
            str(Path(data_path) / "bismark/{sample}_R1_bismark_bt2_pe.deduplicated_splitting_report.txt") \
            if samples_R1R2_present else 
            str(Path(data_path) / "bismark/{sample}_R1_bismark_bt2.deduplicated_splitting_report.txt"),
        bedGraph = \
            str(Path(data_path) / "bismark/{sample}_R1_bismark_bt2_pe.deduplicated.bedGraph.gz") \
            if samples_R1R2_present else 
            str(Path(data_path) / "bismark/{sample}_R1_bismark_bt2.deduplicated.bedGraph.gz"),
        bismark_cov = \
            str(Path(data_path) / "bismark/{sample}_R1_bismark_bt2_pe.deduplicated.bismark.cov.gz") \
            if samples_R1R2_present else 
            str(Path(data_path) / "bismark/{sample}_R1_bismark_bt2.deduplicated.bismark.cov.gz"),
        # ok_file = str(Path(data_path) / "bismark/{sample}_bismark_methylation_extractor_ok.txt"),
    conda:
        get_conda_env('bismark'),
    log:
        str(Path(data_path) / "logs/bismark_methylation_extractor/bismark_methylation_extractor_{sample}.log"),
    threads: get_rule_threads ('bismark_methylation_extractor')
    resources:
        cl_resources = lambda wildcards, attempt : get_cluster_resources_by_attempt (attempt, 'bismark_methylation_extractor'),
        walltime = get_rule_walltime('bismark_methylation_extractor'),
        cl_job_suffix = lambda wildcards : wildcards.sample,
    shell:
        '''
        cd "$(dirname "{output.CHG_context}")"  # get to the parent dir of the output file
        
        bismark_methylation_extractor --multicore {threads} --comprehensive --bedgraph {input.dedup_bam} 2>&1 | tee {log}
        '''
        # echo 'OK' > {output.ok_file}

rule bismark_report:
    input: 
        align_report = \
            str(Path(data_path) / "bismark/work/{sample}_R1_bismark_bt2_PE_report.txt") \
            if samples_R1R2_present else \
            str(Path(data_path) / "bismark/{sample}_R1_bismark_bt2_SE_report.txt"),
    output:
        # ok_file = str(Path(data_path) / "bismark/{sample}_bismark_report_ok.txt"),
        html = str(Path(data_path) / "bismark/{sample}_report.html"),
    conda:
        get_conda_env('bismark'),
    log:
        str(Path(data_path) / "logs/bismark_report/bismark_report_{sample}.log"),
    # threads: get_rule_threads ('bismark_report')
    resources:
        cl_resources = lambda wildcards, attempt : get_cluster_resources_by_attempt (attempt, 'bismark_report'),
        walltime = get_rule_walltime('bismark_report'),
        cl_job_suffix = lambda wildcards : wildcards.sample,
    params:
        # html = str(Path(data_path) / "bismark/{sample}_report.html"),
    shell:
        '''
        cd "$(dirname "{output.html}")"  # get to the parent dir of the output file
        
        bismark2report -o "{output.html}" -a "{input.align_report}" 2>&1 | tee {log}
        '''
        # echo 'OK' > {output.ok_file}

if samples_R1R2_present:
    rule fastq_raw_validate:
        # this rule asserts that number of rows in R1 and R2 files is the same
        input:
            fileR1 = get_R1_file_path(data_path, 'fastq_raw'),
            fileR2 = get_R2_file_path(data_path, 'fastq_raw'),
            # to make sure fastqc_raw rule is completed before starting this one, the following 2 rows are included here
            fastqc_raw_outR1 = get_R1_file_path(data_path, sub_dir_path = 'fastqc_raw', file_name_format = '{sample}_|FN|_fastqc.html'),
            fastqc_raw_outR2 = get_R2_file_path(data_path, sub_dir_path = 'fastqc_raw', file_name_format = '{sample}_|FN|_fastqc.html'),
        output:
            valid_file = str(Path(data_path) / "fastq_raw_validate/valid_{sample}"),
        conda:
            get_conda_env(),
        log:
            str(Path(data_path) / "logs/fastq_raw_validate/fastq_raw_validate_{sample}.log")
        resources:
            cl_job_suffix = lambda wildcards : wildcards.sample,
        params:
            outdir = str(Path(data_path) / "fastq_raw_validate"),
            pipeline_warning_file_path = pipeline_warning_file_path,
        shell:
            '''
            file1="{input.fileR1}"
            file2="{input.fileR2}"
            echo "R1 file: $file1" 2>&1 | tee {log} # do not append here, since this is a first entry to the log file
            echo "R2 file: $file2" 2>&1 | tee -a {log}
            
            r1=$(zcat $file1 | wc -l) # 2>&1 | tee -a {log}
            r2=$(zcat $file2 | wc -l) # 2>&1 | tee -a {log}
            echo "R1 count => $r1" 2>&1 | tee -a {log}
            echo "R2 count => $r2" 2>&1 | tee -a {log}
            
            content="R1=$r1 | R2=$r2"
            
            # Check if r1 and r2 are equal
            if [ "$r1" = "$r2" ] && [ -n "$r1" ] && [ -n "$r2" ]; then
                output_file="valid_{wildcards.sample}"
                warning_msg=""
            else
                output_file="invalid_{wildcards.sample}"
                warning_msg="Warning: fastq_raw_validate step. Sample {wildcards.sample} is invalid - $content"
            fi
            
            out_path={params.outdir}/$output_file
            echo "output file name: $out_path" 2>&1 | tee -a {log}

            # Create the output file with the appropriate name and content
            echo "$content" > "$out_path" # 2>&1 | tee -a {log}
            
            if [ -n "$warning_msg" ]; then
                # If it's not blank, add warning message to the pipeline warning file
                echo $warning_msg >> {params.pipeline_warning_file_path}
            fi
            '''

if samples_R1I1_present:
    # these rules will be in use only if I1 files are present.
    
    rule validate_fastqI1_umi_sequence_lenght:
        input: 
            fastq_file = get_I1_file_path(data_path, 'fastq_raw'),
        output:
            # status_file = str(Path(data_path) / "fastq_raw_validate/status/{sample}_I1"),
            valid_file = str(Path(data_path) / "fastq_raw_validate/valid_{sample}_I1"),
        conda:
            get_conda_env(),
        log:
            str(Path(data_path) / "logs/validate_fastqI1_umi_sequence_lenght/validate_fastqI1_umi_sequence_lenght_{sample}.log"),
        resources:
            # mem = next((el for el in [get_config_value ('rules/validate_fastqI1_umi_sequence_lenght/memory')] if el is not None), default_rule_memory)  # 10000
            cl_resources = lambda wildcards, attempt : get_cluster_resources_by_attempt (attempt, 'validate_fastqI1_umi_sequence_lenght'),
            cl_job_suffix = lambda wildcards : wildcards.sample,
        params:
            invalid_file = str(Path(data_path) / "fastq_raw_validate/invalid_{sample}_I1"), 
            default_seq_len = get_config_value ('rules/validate_fastqI1_umi_sequence_lenght/umi_length'),
            pipeline_warning_file_path = pipeline_warning_file_path,
        script:
            "scripts/validate_fastq_file.py"
   
    rule UMI_attach_R1:
        input:
            fileR1 = get_R1_file_path(data_path, 'fastq_raw'),
            valid_file_R12 = str(Path(data_path) / "fastq_raw_validate/valid_{sample}") if samples_R1R2_present else [],
            valid_file_I1 = str(Path(data_path) / "fastq_raw_validate/valid_{sample}_I1"),
        output:
            attach_fastq = temp_debugcheck(str(Path(data_path) / "fastq_attach/{sample}_R1.fastq.gz")),
            tmp_umi_attach= temp_debugcheck(str(Path(data_path) / "fastq_attach/{sample}_R1.fastq")),
        conda:
            get_conda_env(),
        log:
            str(Path(data_path) / "logs/UMI_attach_R1/UMI_attach_R1_{sample}.log"),
        resources:
            # mem = next((el for el in [get_config_value ('rules/UMI_attach_R1/memory')] if el is not None), default_rule_memory)  # 20000
            cl_resources = lambda wildcards, attempt : get_cluster_resources_by_attempt (attempt, 'UMI_attach_R1'),
            cl_job_suffix = lambda wildcards : wildcards.sample,
        params:
            fileI1 = get_I1_file_path(data_path, 'fastq_raw'),
            umi_script = str(Path(snakemake_path) / "scripts/UMI_attach.awk"),
            # original shell script:
            # zcat {input} | {params.umi_script} -v Ifq={params.fileI1}|gzip -c>{output.attach_fastq} # copy from Y's code; not suited for this implementation
            # zcat {input.fileR1}|{params.umi_script} -v Ifq={params.fileI1}|gzip -c>{output.attach_fastq}
            # modified script
            # zcat {input.fileR1} | {params.umi_script} -v Ifq={params.fileI1} > {output.tmp_umi_attach} 2>&1 | tee {log}
            # gzip -c {output.tmp_umi_attach} > {output.attach_fastq} 2>&1 | tee -a {log} 
        shell:
            '''
            zcat {input.fileR1} | {params.umi_script} -v Ifq={params.fileI1} > {output.tmp_umi_attach} 2>&1 | tee {log}
            gzip -c {output.tmp_umi_attach} > {output.attach_fastq} 2>&1 | tee -a {log}
            '''
    
    if samples_R1R2_present:
        rule UMI_attach_R2:
            input:
                fileR2 = get_R2_file_path(data_path, 'fastq_raw'),
                valid_file_R1R2 = str(Path(data_path) / "fastq_raw_validate/valid_{sample}") if samples_R1R2_present else [],
                valid_file_I1 = str(Path(data_path) / "fastq_raw_validate/valid_{sample}_I1"),
            output:
                attach_fastq = temp_debugcheck(str(Path(data_path) / "fastq_attach/{sample}_R2.fastq.gz")),
                tmp_umi_attach= temp_debugcheck(str(Path(data_path) / "fastq_attach/{sample}_R2.fastq")),
            conda:
                get_conda_env(),
            log:
                str(Path(data_path) / "logs/UMI_attach_R2/UMI_attach_R2_{sample}.log"),
            resources:
                # mem = next((el for el in [get_config_value ('rules/UMI_attach_R2/memory')] if el is not None), default_rule_memory)  # 20000
                cl_resources = lambda wildcards, attempt : get_cluster_resources_by_attempt (attempt, 'UMI_attach_R2'),
                cl_job_suffix = lambda wildcards : wildcards.sample,
            params:
                fileI1 = get_I1_file_path(data_path, 'fastq_raw'),
                umi_script = str(Path(snakemake_path) / "scripts/UMI_attach.awk"),
            shell:
                '''
                zcat {input.fileR2} | {params.umi_script} -v Ifq={params.fileI1} > {output.tmp_umi_attach} 2>&1 | tee {log}
                gzip -c {output.tmp_umi_attach} > {output.attach_fastq} 2>&1 | tee -a {log}
                '''
            

rule trim:
    input:
        fileR1 = get_R1_file_path(data_path, ('fastq_attach' if samples_R1I1_present else 'fastq_raw')),
        fileR2 = get_R2_file_path(data_path, ('fastq_attach' if samples_R1I1_present else 'fastq_raw')) if samples_R1R2_present else [],
        valid_file_R1R2 = str(Path(data_path) / "fastq_raw_validate/valid_{sample}") if samples_R1R2_present else [],
    output:
        trimR1 = temp_debugcheck(get_R1_file_path(data_path, 'fastq_trim')),
        trimR2 = temp_debugcheck(get_R2_file_path(data_path, 'fastq_trim')) if samples_R1R2_present else [],
        trimR1_tooshort = temp_debugcheck(get_R1_file_path(data_path, 'fastq_trim/tooshort')),
        trimR2_tooshort = temp_debugcheck(get_R2_file_path(data_path, 'fastq_trim/tooshort')) if samples_R1R2_present else [],
        trim_log_file = temp_debugcheck(str(Path(data_path) / "logs/trim/tool_logs/log.{sample}")), # required for multiqc_fastqc_trim rule
    conda:
        get_conda_env('cutadapt'),
    log:
        str(Path(data_path) / "logs/trim/log.{sample}")
    resources:
        cl_job_suffix = lambda wildcards : wildcards.sample,
    params:
        tool = "cutadapt",
        data_path = data_path,
        tooshort_path = str(Path(data_path) / "fastq_trim/tooshort"),
        index_adapter = get_config_value ('rules/trim/index_adapter'),
        # the following creates a partial bash script that will be inserted to the shell part
        univ_adapter = "-A " + get_config_value ('rules/trim/univ_adapter') if samples_R1R2_present else "",
        output = str(Path(data_path) / ("fastq_trim/" + os.path.basename(get_R1_file_path(data_path, 'fastq_raw')))),
        # the following creates a partial bash script that will be inserted to the shell part
        paired_output = "--paired-output " + str(Path(data_path) / ("fastq_trim/" + os.path.basename(get_R2_file_path(data_path, 'fastq_raw')))) if samples_R1R2_present else "",
        minimum_length = get_config_value ('rules/trim/minimum_length'),
        too_short_output = str(Path(data_path) / ("fastq_trim/tooshort/" + os.path.basename(get_R1_file_path(data_path, 'fastq_raw')))),
        # the following creates a partial bash script that will be inserted to the shell part
        too_short_paired_output = "--too-short-paired-output " + str(Path(data_path) / ("fastq_trim/tooshort/" + os.path.basename(get_R2_file_path(data_path, 'fastq_raw')))) if samples_R1R2_present else "",
    shell:
        '''
        # cd {params.data_path} # navigate to the data folder location
        mkdir -p {params.tooshort_path} # pre-create needed dir
        {params.tool} \
            --adapter {params.index_adapter} \
            {params.univ_adapter} \
            --output {params.output} \
            {params.paired_output} \
            --minimum-length {params.minimum_length} \
            --too-short-output {params.too_short_output} \
            {params.too_short_paired_output} \
            {input.fileR1} {input.fileR2} 2>&1 | tee {log}
        # copy created log file to the dedicated directory
        cp {log[0]} {output.trim_log_file}
        '''

rule fastqc_trim:
    input:
        fileR1 = get_R1_file_path(data_path, 'fastq_trim'),
        fileR2 = get_R2_file_path(data_path, 'fastq_trim') if samples_R1R2_present else []
    output:
        # outR1 = str(Path(data_path) / "fastqc_trim/{sample}_R1_fastqc.html"),
        outR1 = get_R1_file_path(data_path, sub_dir_path = 'fastqc_trim', file_name_format = '{sample}_|FN|_fastqc.html'),
        # outR2 = str(Path(data_path) / "fastqc_trim/{sample}_R2_fastqc.html") if samples_R1R2_present else [],
        outR2 = get_R2_file_path(data_path, sub_dir_path = 'fastqc_trim', file_name_format = '{sample}_|FN|_fastqc.html') if samples_R1R2_present else [],
    conda:
        get_conda_env(),
    log:
        str(Path(data_path) / "logs/fastqc_trim/fastqc_trim_{sample}.log")
    resources:
        # mem = next((el for el in [get_config_value ('rules/fastqc_trim/memory')] if el is not None), default_rule_memory)  # 40000
        cl_resources = lambda wildcards, attempt : get_cluster_resources_by_attempt (attempt, 'fastqc_trim'),
        cl_job_suffix = lambda wildcards : wildcards.sample,
    params:
        tool = "fastqc",
        out_dir = str(Path(data_path) / "fastqc_trim")
    shell:
        '''
        {params.tool} --outdir {params.out_dir} {input.fileR1} {input.fileR2} 2>&1 | tee {log}
        '''

rule star_align:
    # TODO: There is a possibility of running star_aliqn based on the original fastq files. Future implementation
    # should allow running star_align without trim rule. Current idea is to have a boolean configuration variable
    # that will define the path to be taken and based on this input files for this rule will differ.
    input:
        trimR1 = get_R1_file_path(data_path, 'fastq_trim'),
        trimR2 = get_R2_file_path(data_path, 'fastq_trim') if samples_R1R2_present else []
    output:
        bam=str(Path(data_path) / "star_align/{sample}/Aligned.sortedByCoord.out.bam"),
        rsem_bam=temp_debugcheck(str(Path(data_path) / "star_align/{sample}/Aligned.toTranscriptome.out.bam")),
        qc_log=str(Path(data_path) / "star_align/{sample}/Log.final.out"),
        sample_dir=directory(str(Path(data_path) / "star_align/{sample}")),
    conda:
        get_conda_env(),
    log:
        str(Path(data_path) / "logs/star_align/star_align_{sample}.log")
    threads: get_rule_threads ('star_align')
    resources:
        # mem = next((el for el in [get_config_value ('rules/star_align/memory')] if el is not None), default_rule_memory)  # 7000
        cl_resources = lambda wildcards, attempt : get_cluster_resources_by_attempt (attempt, 'star_align'),
        cl_job_suffix = lambda wildcards : wildcards.sample,
    params:
        tool = "STAR",
        genomeDir = str(Path(motrpac_ref_data_path) / (motrpac_ref_data_genom_dir + '/star_index')),
        sjdbOverhang = get_config_value ('rules/star_align/sjdbOverhang') if get_config_value ('rules/star_align/sjdbOverhang') else 100,
        outFileNamePrefix = str(Path(data_path) / "star_align/{sample}"),
        readFilesCommand = get_config_value ('rules/star_align/readFilesCommand') if get_config_value ('rules/star_align/readFilesCommand') else "zcat",
        outSAMattributes = get_config_value ('rules/star_align/outSAMattributes') if get_config_value ('rules/star_align/outSAMattributes') else "NH HI AS NM MD nM",
        outSAMtype = get_config_value ('rules/star_align/outSAMtype') if get_config_value ('rules/star_align/outSAMtype') else "BAM SortedByCoordinate",
        outFilterType = get_config_value ('rules/star_align/outFilterType') if get_config_value ('rules/star_align/outFilterType') else "BySJout",
        quantMode = get_config_value ('rules/star_align/quantMode') if get_config_value ('rules/star_align/quantMode') else "TranscriptomeSAM",
        outSAMmultNmax = get_config_value ('rules/star_align/outSAMmultNmax') if get_config_value ('rules/star_align/outSAMmultNmax') else -1, 
        outMultimapperOrder = get_config_value ('rules/star_align/outMultimapperOrder') if get_config_value ('rules/star_align/outMultimapperOrder') else "Old_2.4", 
        outFilterScoreMinOverLread = get_config_value ('rules/star_align/outFilterScoreMinOverLread') if get_config_value ('rules/star_align/outFilterScoreMinOverLread') else "0.66", 
        outFilterMatchNminOverLread = get_config_value ('rules/star_align/outFilterMatchNminOverLread') if get_config_value ('rules/star_align/outFilterMatchNminOverLread') else "0.66", 
        outTmpDir = str(Path(data_path) / ("star_align/tmpdir_{sample}")),
        # get the current memory allocation, deduct 1 GB (1000 MB) and convert to bytes
        limitBAMsortRAM = (int(next((el for el in [get_config_value ('rules/star_align/memory')] if el is not None), default_rule_memory)) - 1000) * 1000000,  
    shell:
        '''
        rm -rf {params.outTmpDir} || true # removes temp dir if it left from the previous run; it does nothing otherwise
        {params.tool} \
        --genomeDir {params.genomeDir} \
        --sjdbOverhang {params.sjdbOverhang} \
        --readFilesIn {input.trimR1} {input.trimR2} \
        --outFileNamePrefix {params.outFileNamePrefix}/ \
        --readFilesCommand {params.readFilesCommand} \
        --outSAMattributes {params.outSAMattributes} \
        --runThreadN {threads} \
        --outSAMtype {params.outSAMtype} \
        --outFilterType {params.outFilterType} \
        --quantMode {params.quantMode} \
        --outSAMmultNmax {params.outSAMmultNmax} \
        --outMultimapperOrder {params.outMultimapperOrder} \
        --outFilterScoreMinOverLread {params.outFilterScoreMinOverLread} \
        --outFilterMatchNminOverLread {params.outFilterMatchNminOverLread} \
        --outTmpDir {params.outTmpDir} \
        --limitBAMsortRAM {params.limitBAMsortRAM} \
        2>&1 | tee {log}
        '''

rule star_align_bai:
    input:
        bam=str(Path(data_path) / "star_align/{sample}/Aligned.sortedByCoord.out.bam"),
    output:
        bai = str(Path(data_path) / "star_align/{sample}/Aligned.sortedByCoord.out.bam.bai"),
    conda:
        get_conda_env(),
    log:
        str(Path(data_path) / "logs/star_align_bai/star_align_bai_{sample}.log")
    resources:
        cl_job_suffix = lambda wildcards : wildcards.sample,
    params:
        tool = "samtools",
        cmd = "index",
    shell:
        '''
        {params.tool} {params.cmd} {input.bam} 2>&1 | tee {log}
        '''

rule star_align_stats:
    input: 
        bam_files = expand(str(Path(data_path) / "star_align/{sample}/Aligned.sortedByCoord.out.bam"), sample=samples_to_process),
        bam_bai_files = expand(str(Path(data_path) / "star_align/{sample}/Aligned.sortedByCoord.out.bam.bai"), sample=samples_to_process),
    output:
        bam_read_counts = temp_debugcheck(expand(str(Path(data_path) / "star_align/{sample}/Aligned.sortedByCoord.out.bam.reads_count"), sample=samples_to_process)),
        star_align_stats = temp_debugcheck(str(Path(data_path) / "star_align/star_align_stats.txt")),
    conda:
        get_conda_env(),
    log:
        str(Path(data_path) / "logs/star_align_stats/star_align_stats.log")
    threads: get_rule_threads ('star_align_stats')
    resources:
        # mem = next((el for el in [get_config_value ('rules/star_align_stats/memory')] if el is not None), default_rule_memory)
        cl_resources = lambda wildcards, attempt : get_cluster_resources_by_attempt (attempt, 'star_align_stats'),
        walltime = get_rule_walltime('star_align_stats'), # "24:00",
    params:
        # pipeline_info_file_path = pipeline_info_file_path,
        pipeline_warning_file_path = pipeline_warning_file_path,
        parallel_count = next((el for el in [get_config_value ('rules/star_align_stats/parallel_count')] if el is not None), True),
        count_mapped_only = next((el for el in [get_config_value ('rules/star_align_stats/count_mapped_only')] if el is not None), True),
    script:
        "scripts/star_align_stats.py"

rule star_align_validate:
    input:
        bam_read_count = str(Path(data_path) / "star_align/{sample}/Aligned.sortedByCoord.out.bam.reads_count"),
        star_align_stats = str(Path(data_path) / "star_align/star_align_stats.txt"),
    output:
        valid_file = temp_debugcheck(str(Path(data_path) / "star_align_validate/valid_{sample}")),
    # conda:
    #     get_conda_env(),
    log:
        str(Path(data_path) / "logs/star_align_validate/star_align_validate_{sample}.log")
    resources:
        cl_job_suffix = lambda wildcards : wildcards.sample,
    params:
        out_dir = str(Path(data_path) / "star_align_validate"),
        outlier_threshold = next((el for el in [get_config_value ('rules/star_align_validate/outlier_threshold')] if el is not None), -2.0),
        min_reads_count = next((el for el in [get_config_value ('rules/star_align_validate/min_reads_count')] if el is not None), 50000),
        validate_z_score_max_read_count = next((el for el in [get_config_value ('rules/star_align_validate/validate_z_score_max_read_count')] if el is not None), 15000000),
        pipeline_warning_file_path = pipeline_warning_file_path,
    script:
        "scripts/star_align_validate.py"

rule star_align_intronMotif:
    input:
        trimR1 = get_R1_file_path(data_path, 'fastq_trim'),
        trimR2 = get_R2_file_path(data_path, 'fastq_trim') if samples_R1R2_present else []
    output:
        bam=str(Path(data_path) / "star_align_intronMotif/{sample}/Aligned.sortedByCoord.out.bam"), # Aligned.out.bam  #Aligned.sortedByCoord.out.bam
        qc_log=str(Path(data_path) / "star_align_intronMotif/{sample}/Log.final.out"),  # Log.progress.out
    conda:
        get_conda_env(),
    log:
        str(Path(data_path) / "logs/star_align_intronMotif/star_align_{sample}.log")
    threads: get_rule_threads ('star_align_intronMotif')
    resources:
        # mem = next((el for el in [get_config_value ('rules/star_align_intronMotif/memory')] if el is not None), default_rule_memory)  # 7000
        cl_resources = lambda wildcards, attempt : get_cluster_resources_by_attempt (attempt, 'star_align_intronMotif'),
        cl_job_suffix = lambda wildcards : wildcards.sample,
    params:
        tool = "STAR",
        genomeDir = str(Path(motrpac_ref_data_path) / (motrpac_ref_data_genom_dir + '/star_index')),
        outFileNamePrefix = str(Path(data_path) / "star_align_intronMotif/{sample}"),
        readFilesCommand = (get_config_value ('rules/star_align_intronMotif/readFilesCommand') or "zcat"),
        outSAMtype = (get_config_value ('rules/star_align_intronMotif/outSAMtype') or "BAM SortedByCoordinate"),
        twopassMode = (get_config_value ('rules/star_align_intronMotif/twopassMode') or "Basic"),
        outSAMstrandField = (get_config_value ('rules/star_align_intronMotif/outSAMstrandField') or "intronMotif"),
        outTmpDir = str(Path(data_path) / ("star_align_intronMotif/tmpdir_{sample}")),
    shell:
        '''
        rm -rf {params.outTmpDir} || true # removes temp dir if it left from the previous run; it does nothing otherwise
        {params.tool} \
        --genomeDir {params.genomeDir} \
        --readFilesIn {input.trimR1} {input.trimR2} \
        --outFileNamePrefix {params.outFileNamePrefix}/ \
        --readFilesCommand {params.readFilesCommand} \
        --runThreadN {threads} \
        --outSAMtype {params.outSAMtype} \
        --twopassMode Basic \
        --outSAMstrandField intronMotif \
        --outTmpDir {params.outTmpDir} \
        2>&1 | tee {log}
        '''

rule star_align_intronMotif_bai:
    input:
        bam=str(Path(data_path) / "star_align_intronMotif/{sample}/Aligned.sortedByCoord.out.bam"),
    output:
        bai = str(Path(data_path) / "star_align_intronMotif/{sample}/Aligned.sortedByCoord.out.bam.bai"),
    conda:
        get_conda_env(),
    log:
        str(Path(data_path) / "logs/star_align_intronMotif_bai/star_align_intronMotif_bai_{sample}.log")
    resources:
        cl_job_suffix = lambda wildcards : wildcards.sample,
    params:
        tool = "samtools",
        cmd = "index",
    shell:
        '''
        {params.tool} {params.cmd} {input.bam} 2>&1 | tee {log}
        '''

rule intronMotifBAM2junc:
    input:
        bam=str(Path(data_path) / "star_align_intronMotif/{sample}/Aligned.sortedByCoord.out.bam"),
        # index file (bai) is used implicitly by the tool, thus it is required to be available
        bai = str(Path(data_path) / "star_align_intronMotif/{sample}/Aligned.sortedByCoord.out.bam.bai"),
    output:
        intronMotif_junc = str(Path(data_path) / "leafcutter/regtools/{sample}_Aligned.sortedByCoord.out.bam.junc"),
    conda:
        get_conda_env(),
    log:
        str(Path(data_path) / "logs/intronMotifBAM2junc/intronMotifBAM2junc_{sample}.log")
    resources:
        cl_job_suffix = lambda wildcards : wildcards.sample,
    params:
        tool = "regtools",
        cmd = "junctions extract",
        anchor_length = (get_config_value ('rules/intronMotifBAM2junc/anchor_length') or "8"),
        min_intron_len = (get_config_value ('rules/intronMotifBAM2junc/min_intron_len') or "50"),
        max_intron_len = (get_config_value ('rules/intronMotifBAM2junc/max_intron_len') or "500000"),
#        library_strandedness = leafcutter_library_strandedness,
    shell:
        # grep -v "^GL" "{output.intronMotif_junc}" > {output.juncfile_txt}
        '''
        {params.tool} {params.cmd} \
        -a {params.anchor_length} \
        -m {params.min_intron_len} \
        -M {params.max_intron_len} \
        -s {params.library_strandedness} \
        {input.bam} \
        -o {output.intronMotif_junc} \
        2>&1 | tee {log}
        '''

rule intronMotifBAM2junc_list:
    input:
        intronMotif_junc_files = expand(str(Path(data_path) / "leafcutter/regtools/{sample}_Aligned.sortedByCoord.out.bam.junc"), sample=samples_to_process),
    output:
        juncfiles_txt = str(Path(data_path) / "leafcutter/regtools/_juncfiles.txt"),
    conda:
        get_conda_env(),
    log:
        str(Path(data_path) / "logs/intronMotifBAM2junc_list/intronMotifBAM2junc_list.log")
    # localrule: True
    params:
        # pipeline_info_file_path = pipeline_info_file_path,
        pipeline_warning_file_path = pipeline_warning_file_path,
    script:
        "scripts/intronMotifBAM2junc_list.py"

rule leafcutter_cluster_regtools:
    input:
        intronMotif_junc_files = expand(str(Path(data_path) / "leafcutter/regtools/{sample}_Aligned.sortedByCoord.out.bam.junc"), sample=samples_to_process),
        juncfiles_txt = str(Path(data_path) / "leafcutter/regtools/_juncfiles.txt"),
    output:
        # cluster_regtools = str(Path(data_path) / ("leafcutter/cluster/" + leafcutter_library_strandedness + "_perind_numers.counts.gz")),
        cluster_regtools = str(Path(data_path) / ("leafcutter/cluster/test_contrast_perind_numers.counts.gz")),
    conda:
        get_conda_env(),
    log:
        str(Path(data_path) / "logs/leafcutter_cluster_regtools/leafcutter_cluster_regtools.log")
    params:
        python = "python",
        tool = str(Path(snakemake_path) / "scripts/leafcutter/clustering/leafcutter_cluster_regtools.py"),
        split_reads_required = (get_config_value ('rules/leafcutter_cluster_regtools/split_reads_required') or "50"),
        intron_length_allowed = (get_config_value ('rules/leafcutter_cluster_regtools/intron_length_allowed') or "500000"),
        outputFileName_prefix = "test_contrast",  # leafcutter_library_strandedness,  # "leafcutter_out",  # leafcutter_library_strandedness,
        wrk_dir = str(Path(data_path) / "leafcutter/cluster"),
        # --nochromcheck=NOCHROMCHECK
    shell:
        '''
        cd {params.wrk_dir}
        {params.python} {params.tool} \
        -j {input.juncfiles_txt} \
        -m {params.split_reads_required} \
        -o {params.outputFileName_prefix} \
        -l {params.intron_length_allowed} \
        2>&1 | tee {log}
        '''

rule leafcutter_contrast_file:
    input:
        cluster_regtools = str(Path(data_path) / ("leafcutter/cluster/test_contrast_perind_numers.counts.gz")),
    output:
        leafcutter_contrast = str(Path(data_path) / "metadata/leafcutter_contrast.txt"),
        valid_file = str(Path(data_path) / "leafcutter/valid_leafcutter_contrast.txt"),
    # conda:
    #     get_conda_env(),
    log:
        str(Path(data_path) / "logs/leafcutter_contrast_file/leafcutter_contrast_file.log")
    localrule: True  # this step uses connection to the SQL Server that is currently not supported from the cluster
    params:
        db_cfg = db_cfg,
        samples = samples_to_process,
        db_study_id = metadata_db_study_id,
        db_center_id = metadata_db_center_id,
#        db_mapping_code = leafcutter_contrast_db_mapping_code,
#        contrast_group_count_min = leafcutter_contrast_group_count_min,
        pipeline_info_file_path = pipeline_info_file_path,
        pipeline_warning_file_path = pipeline_warning_file_path,
        rule_name = "leafcutter_contrast_file",
        invalid_file = str(Path(data_path) / "leafcutter/invalid_leafcutter_contrast.txt"),
    script:
        "scripts/create_leafcutter_contrast.py"

# TODO: this rule currently produces errors, fix is required
rule leafcutter_differential_splicing:
    input:
        cluster_regtools = str(Path(data_path) / ("leafcutter/cluster/test_contrast_perind_numers.counts.gz")),
        metadata = str(Path(data_path) / "metadata/leafcutter_contrast.txt"),
    output:
       ds_cluster_significance = str(Path(data_path) / ("leafcutter/ds/leafcutter_ds_cluster_significance.txt")),
    conda:
        get_conda_env(),
    log:
        str(Path(data_path) / "logs/leafcutter_differential_splicing/leafcutter_differential_splicing.log")
    conda:
        get_conda_env('leafcutter')
    threads: get_rule_threads ('leafcutter_differential_splicing') 
    params:
        r = "Rscript",
        tool = str(Path(snakemake_path) / "scripts/leafcutter/scripts/leafcutter_ds.R"),
        wrk_dir = str(Path(data_path) / "leafcutter/ds"),
    shell:
        '''
        cd {params.wrk_dir}
        {params.r} --vanilla \
        {params.tool} \
        --num_threads {threads} \
        {input.cluster_regtools} \
        {input.metadata} \
        2>&1 | tee {log}
        '''
        # -i 2 -g 2 \  # can be used for testing only, cannot be passed to production

rule star_align_djexpress:
    input:
        trimR1 = get_R1_file_path(data_path, 'fastq_trim'),
        trimR2 = get_R2_file_path(data_path, 'fastq_trim') if samples_R1R2_present else []
    output:
        bam=temp_debugcheck(str(Path(data_path) / "star_align_djexpress/{sample}/Aligned.sortedByCoord.out.bam")), 
        sj_file=str(Path(data_path) / "star_align_djexpress/{sample}/SJ.out.tab"),  # Log.progress.out
    conda:
        get_conda_env(),
    log:
        str(Path(data_path) / "logs/star_align_djexpress/star_align_djexpress_{sample}.log")
    threads: get_rule_threads ('star_align_djexpress')
    resources:
        # mem = next((el for el in [get_config_value ('rules/star_align_djexpress/memory')] if el is not None), default_rule_memory)  # 7000
        cl_resources = lambda wildcards, attempt : get_cluster_resources_by_attempt (attempt, 'star_align_djexpress'),
        cl_job_suffix = lambda wildcards : wildcards.sample,
    params:
        tool = "STAR",
        genomeDir = str(Path(motrpac_ref_data_path) / (motrpac_ref_data_genom_dir + '/star_index')),
        genomeGtf = str(Path(motrpac_ref_data_path) / (motrpac_ref_data_genom_dir + '/genome.gtf')), 
        outFileNamePrefix = str(Path(data_path) / "star_align_djexpress/{sample}"),
        readFilesCommand = (get_config_value ('rules/star_align_djexpress/readFilesCommand') or "zcat"),
        outSAMtype = (get_config_value ('rules/star_align_djexpress/outSAMtype') or "BAM SortedByCoordinate"),
        twopassMode = (get_config_value ('rules/star_align_djexpress/twopassMode') or "Basic"),
        outSAMstrandField = (get_config_value ('rules/star_align_djexpress/outSAMstrandField') or "intronMotif"),
        outTmpDir = str(Path(data_path) / ("star_align_djexpress/tmpdir_{sample}")),
    shell:
        '''
        rm -rf {params.outTmpDir} || true # removes temp dir if it left from the previous run; it does nothing otherwise
        {params.tool} \
        --genomeDir {params.genomeDir} \
        --readFilesIn {input.trimR1} {input.trimR2} \
        --outFileNamePrefix {params.outFileNamePrefix}/ \
        --readFilesCommand {params.readFilesCommand} \
        --runThreadN {threads} \
        --outSAMtype {params.outSAMtype} \
        --twopassMode Basic \
        --outSAMstrandField intronMotif \
        --outFilterMultimapScoreRange 1 \
        --outFilterMultimapNmax 20 \
        --outFilterMismatchNmax 10 \
        --alignIntronMax 500000 \
        --alignMatesGapMax 1000000 \
        --sjdbScore 2 \
        --alignSJDBoverhangMin 1 \
        --outFilterMatchNminOverLread 0.33 \
        --outFilterScoreMinOverLread 0.33 \
        --sjdbOverhang 100 \
        --outSAMattributes NH HI NM MD AS XS \
        --sjdbGTFfile {params.genomeGtf} \
        --limitSjdbInsertNsj 2000000 \
        --outSAMunmapped None \
        --outSAMheaderHD @HD VN:1.4 \
        --outSAMattrRGline ID:: \
        --outSAMmultNmax 1 \
        --outTmpDir {params.outTmpDir} \
        2>&1 | tee {log}
        '''

rule star_align_djexpress_named_junc:
    input:
        junc_file=str(Path(data_path) / "star_align_djexpress/{sample}/SJ.out.tab"),
    output:
        named_junc_file = str(Path(data_path) / "star_align_djexpress/junct_quant/{sample}_SJ.out.tab"),
    log:
        str(Path(data_path) / "logs/star_align_djexpress_named_junc/star_align_djexpress_named_junc_{sample}.log")
    # localrule: True
    resources:
        cl_resources = lambda wildcards, attempt : get_cluster_resources_by_attempt (attempt, 'star_align_djexpress_named_junc'),
        cl_job_suffix = lambda wildcards : wildcards.sample,
    shell:
        '''
        ln -srf "{input.junc_file}" "{output.named_junc_file}" 2>&1 | tee {log}
        '''

rule star_align_named_bams:
    input:
        bam=str(Path(data_path) / "star_align/{sample}/Aligned.sortedByCoord.out.bam"),
        star_align_valid_file = str(Path(data_path) / "star_align_validate/valid_{sample}"),
    output:
        bam = str(Path(data_path) / "star_align/named_files/{sample}_Aligned.sortedByCoord.out.bam"),
    log:
        str(Path(data_path) / "logs/star_align_named_bams/star_align_named_bams_{sample}.log")
    # localrule: True
    resources:
        cl_resources = lambda wildcards, attempt : get_cluster_resources_by_attempt (attempt, 'star_align_named_bams'),
        cl_job_suffix = lambda wildcards : wildcards.sample,
    shell:
        '''
        ln -srf "{input.bam}" "{output.bam}" 2>&1 | tee {log}
        '''

rule star_align_spladder_bai:
    input:
        # bam=str(Path(data_path) / "star_align/{sample}/Aligned.sortedByCoord.out.bam"),
        bam = str(Path(data_path) / "star_align/named_files/{sample}_Aligned.sortedByCoord.out.bam"),
        star_align_valid_file = str(Path(data_path) / "star_align_validate/valid_{sample}"),
    output:
        # bai = str(Path(data_path) / "star_align/{sample}/Aligned.sortedByCoord.out.bam.bai"),
        bai = str(Path(data_path) / "star_align/named_files/{sample}_Aligned.sortedByCoord.out.bam.bai"),
    conda:
        get_conda_env(),
    log:
        str(Path(data_path) / "logs/star_align_spladder_bai/star_align_spladder_bai_{sample}.log")
    resources:
        cl_job_suffix = lambda wildcards : wildcards.sample,
    params:
        tool = "samtools",
        cmd = "index",
    shell:
        '''
        {params.tool} {params.cmd} {input.bam} 2>&1 | tee {log}
        '''

rule spladder_single_graphs:
    input:
        bam = str(Path(data_path) / "star_align/named_files/{sample}_Aligned.sortedByCoord.out.bam"),
        bai = str(Path(data_path) / "star_align/named_files/{sample}_Aligned.sortedByCoord.out.bam.bai"),
    output:
        count_file = str(Path(data_path) / "spladder/spladder/genes_graph_conf3.{sample}_Aligned.sortedByCoord.out.count.hdf5"), 
        gene_exp = str(Path(data_path) / "spladder/spladder/genes_graph_conf3.{sample}_Aligned.sortedByCoord.out.gene_exp.hdf5"),
        pickle = str(Path(data_path) / "spladder/spladder/genes_graph_conf3.{sample}_Aligned.sortedByCoord.out.pickle"),
    conda:
        get_conda_env(),
    log:
        str(Path(data_path) / "logs/spladder_single_graphs/spladder_single_graphs_{sample}.log")
    threads: get_rule_threads ('spladder_single_graphs')
    resources:
        # mem = next((el for el in [get_config_value ('rules/spladder_single_graphs/memory')] if el is not None), default_rule_memory),
        cl_resources = lambda wildcards, attempt : get_cluster_resources_by_attempt (attempt, 'spladder_single_graphs'),
        cl_job_suffix = lambda wildcards : wildcards.sample,
        walltime = get_rule_walltime('spladder_single_graphs'), # "24:00",
    params:
        out_dir = str(Path(data_path) / "spladder"),
        annotation = str(Path(alt_sl_ref_data_path) / ('spladder/' + motrpac_ref_data_genom_dir + '/genome.gtf')),
        debug = "-d" if debug else "",
        tmp_dir = str(Path(data_path) / "spladder/tmp"),
    shell:
        '''
        mkdir -p {params.tmp_dir}
        
        spladder build \
            {params.debug} \
            -v \
            --parallel {threads} \
            -o {params.out_dir} \
            -a {params.annotation} \
            -b {input.bam} \
            --merge-strat single \
            --no-extract-ase \
            2>&1 | tee {log}
        '''
        
rule star_align_spladder_bam_files_list:
    input:
        input_files = expand(str(Path(data_path) / "star_align/named_files/{sample}_Aligned.sortedByCoord.out.bam"), sample=samples_to_process),
    output:
        combined_files = str(Path(data_path) / "spladder/alignments.txt"),
    conda:
        get_conda_env(),
    log:
        str(Path(data_path) / "logs/star_align_spladder_bam_files_list/star_align_spladder_bam_files_list.log")
    # localrule: True
    params:
        # pipeline_info_file_path = pipeline_info_file_path,
        pipeline_warning_file_path = pipeline_warning_file_path,
    script:
        "scripts/combine_files_to_list.py"
            
rule spladder_merged_graph_chuncks:
    input:
        bam_files = expand(str(Path(data_path) / "star_align/named_files/{sample}_Aligned.sortedByCoord.out.bam"), sample=samples_to_process),
        count_file = expand(str(Path(data_path) / "spladder/spladder/genes_graph_conf3.{sample}_Aligned.sortedByCoord.out.count.hdf5"), sample=samples_to_process),
        alignments = str(Path(data_path) / "spladder/alignments.txt"),
#        chunks = spladder_chunks_input  # function to dynamically assign values at run time
    output:
        spladder_merged_graph_chunck = str(Path(data_path) / "spladder/status_files/spladder_merged_graph_chunk_{chunk}.txt"),
    conda:
        get_conda_env(),
    log:
        str(Path(data_path) / "logs/spladder_merged_graph_chuncks/spladder_merged_graph_chunk_{chunk}.log"),
    # localrule: True
    threads: get_rule_threads ('spladder_merged_graph_chuncks')
    retries: get_rule_retries('spladder_merged_graph_chuncks'), 
    resources:
        # lambda is used to trigger callable function to retrieve the "attempt" value
        # mem = lambda wildcards, attempt : get_memory_by_attempt_and_threads( \
        #                                 attempt_num = attempt, \
        #                                 memory = next((el for el in [get_config_value ('rules/spladder_merged_graph_chuncks/memory')] if el is not None), default_rule_memory), \
        #                                 threads = next((el for el in [get_config_value ('rules/spladder_merged_graph_chuncks/threads')] if el is not None), 8), \
        #                                 max_memory = next((el for el in [get_config_value ('resources/max_memory')] if el is not None), default_max_memory) \
        #                                 ),
        cl_resources = lambda wildcards, attempt : get_cluster_resources_by_attempt (attempt, 'spladder_merged_graph_chuncks'),
        attempt = get_attempt,
        # define walltime values to be used by LSF
#        walltime = get_walltime_by_samples_count( \
#                        num_samples = spladder_chunk_size, \
#                        time_coef = get_config_value ('rules/spladder_merged_graph_chuncks/time_coef'), \
#                        min_walltime = get_config_value ('rules/spladder_merged_graph_chuncks/min_walltime'), \
#                        max_regular_walltime = get_config_value ('resources/max_regular_walltime'), \
#                        max_long_walltime = get_config_value ('resources/max_long_walltime')),
        # define name of the queue to be used with LSF
#        lsf_queue = get_queue_by_samples_count( \
#                        num_samples = spladder_chunk_size, \
#                        time_coef = get_config_value ('rules/spladder_merged_graph_chuncks/time_coef'), \
#                        min_walltime = get_config_value ('rules/spladder_merged_graph_chuncks/min_walltime'), \
#                        max_regular_walltime = get_config_value ('resources/max_regular_walltime'), \
#                        max_long_walltime = get_config_value ('resources/max_long_walltime'), \
#                        normal_queue = get_config_value ('resources/regular_queue'), \
#                        long_queue = get_config_value ('resources/long_queue')),
        cl_job_suffix = lambda wildcards : wildcards.chunk,
    params:
        out_dir = str(Path(data_path) / "spladder"),
        annotation = str(Path(alt_sl_ref_data_path) / ('spladder/' + motrpac_ref_data_genom_dir + '/genome.gtf')),
        debug = "-d" if debug else "",
        chunk_params = '{chunk}',  # .replace('_', ' '),
#        chunk_size = str(spladder_chunk_size),
    shell:
        '''
        chunk_args=$(echo "{params.chunk_params}" | tr '_' ' ')  # replace "_"s with " "s
        echo "chunk_args = "$chunk_args 2>&1 | tee {log}
        
        spladder build \
            {params.debug} \
            --parallel {threads} \
            -o {params.out_dir} \
            -a {params.annotation} \
            -b {input.alignments} \
            --merge-strat merge_graphs \
            --no-extract-ase \
            --no-quantify-graph \
            --chunksize {params.chunk_size} \
            --chunked-merge $chunk_args  \
            2>&1 | tee -a {log}
            
            echo 'OK' > {output.spladder_merged_graph_chunck}
        '''

rule spladder_quantification1:
    input:
        bam = str(Path(data_path) / "star_align/named_files/{sample}_Aligned.sortedByCoord.out.bam"),
        # spladder_merged_graph_count_file = str(Path(data_path) / "spladder/spladder/genes_graph_conf3.merge_graphs.count.hdf5"),  # this to utilize original spladder merge approach
        spladder_merged_graph_chunck_top_level = expand(str(Path(data_path) / "spladder/status_files/spladder_merged_graph_chunk_{chunk}.txt"), chunk=[spl_merge_chunks['level'][len(spl_merge_chunks['level'])][0]]) if run_spladder else [],  # this is to utilize spladder chunk merge approach
    output:
        count_file = str(Path(data_path) / "spladder/spladder/genes_graph_conf3.merge_graphs.{sample}_Aligned.sortedByCoord.out.count.hdf5"), 
        gene_exp = str(Path(data_path) / "spladder/spladder/genes_graph_conf3.merge_graphs.{sample}_Aligned.sortedByCoord.out.gene_exp.hdf5"),
    conda:
        get_conda_env(),
    log:
        str(Path(data_path) / "logs/spladder_quantification1/spladder_quantification1_{sample}.log")
    threads: get_rule_threads ('spladder_quantification1')
    retries: get_rule_retries('spladder_quantification1'), 
    resources:
        # mem = next((el for el in [get_config_value ('rules/spladder_quantification1/memory')] if el is not None), default_rule_memory),
        # mem = lambda wildcards, attempt : get_memory_by_attempt_and_threads( \
        #                                     attempt_num = attempt, \
        #                                     memory = next((el for el in [get_config_value ('rules/spladder_quantification1/memory')] if el is not None), default_rule_memory), \
        #                                     threads = next((el for el in [get_config_value ('rules/spladder_quantification1/threads')] if el is not None), 1), \
        #                                     max_memory = next((el for el in [get_config_value ('resources/max_memory')] if el is not None), default_max_memory) \
        #                                     ),
        cl_resources = lambda wildcards, attempt : get_cluster_resources_by_attempt (attempt, 'spladder_quantification1'),
        cl_job_suffix = lambda wildcards : wildcards.sample,
        attempt = get_attempt,
        walltime = get_rule_walltime('spladder_quantification1'),  # "24:00",
    params:
        out_dir = str(Path(data_path) / "spladder"),
        annotation = str(Path(alt_sl_ref_data_path) / ('spladder/' + motrpac_ref_data_genom_dir + '/genome.gtf')),
        debug = "-d" if debug else "",
    shell:
        '''
        spladder build \
            {params.debug} \
            -v \
            --parallel {threads} \
            -o {params.out_dir} \
            -a {params.annotation} \
            -b {input.bam} \
            --merge-strat merge_graphs \
            --no-extract-ase \
            --quantify-graph \
            --qmode single \
            2>&1 | tee {log}
        '''

rule spladder_quantification2:
    input:
        bam_files = expand(str(Path(data_path) / "star_align/named_files/{sample}_Aligned.sortedByCoord.out.bam"), sample=samples_to_process),
        spladder_quantification1_count_file = expand(str(Path(data_path) / "spladder/spladder/genes_graph_conf3.merge_graphs.{sample}_Aligned.sortedByCoord.out.count.hdf5"), sample=samples_to_process),
        alignments = str(Path(data_path) / "spladder/alignments.txt"),
    output:
        # this rule does not produce any visible output, so the formal status file is used to identify the step completion
        status_file = str(Path(data_path) / "spladder/status_files/spladder_quantification2.txt"),
    conda:
        get_conda_env(),
    log:
        str(Path(data_path) / "logs/spladder_quantification2/spladder_quantification2.log")
    threads: get_rule_threads ('spladder_quantification2')
    resources:
        # mem = next((el for el in [get_config_value ('rules/spladder_quantification2/memory')] if el is not None), default_rule_memory)
        cl_resources = lambda wildcards, attempt : get_cluster_resources_by_attempt (attempt, 'spladder_quantification2'),
    params:
        out_dir = str(Path(data_path) / "spladder"),
        annotation = str(Path(alt_sl_ref_data_path) / ('spladder/' + motrpac_ref_data_genom_dir + '/genome.gtf')),
        debug = "-d" if debug else "",
    shell:
        '''
        spladder build \
            {params.debug} \
            -v \
            --parallel {threads} \
            -o {params.out_dir} \
            -a {params.annotation} \
            -b {input.alignments} \
            --merge-strat merge_graphs \
            --no-extract-ase \
            --quantify-graph \
            --qmode collect \
            2>&1 | tee {log}
            
            echo 'OK' > {output.status_file}
        '''

rule spladder_call_events:
    input:
        bam_files = expand(str(Path(data_path) / "star_align/named_files/{sample}_Aligned.sortedByCoord.out.bam"), sample=samples_to_process),
        spladder_quantification2_status_file = str(Path(data_path) / "spladder/status_files/spladder_quantification2.txt"),
        alignments = str(Path(data_path) / "spladder/alignments.txt"),
    output:
        merge_graphs_confirmed_txt = str(Path(data_path) / "spladder/merge_graphs_{event}_C3.confirmed.txt.gz"),
        merge_graphs_confirmed_pickle = str(Path(data_path) / "spladder/merge_graphs_{event}_C3.confirmed.pickle"),
        merge_graphs_confirmed_gff3 = str(Path(data_path) / "spladder/merge_graphs_{event}_C3.confirmed.gff3"),
        merge_graphs_counts = str(Path(data_path) / "spladder/merge_graphs_{event}_C3.counts.hdf5"),
        merge_graphs_pickle = str(Path(data_path) / "spladder/merge_graphs_{event}_C3.pickle"),
    conda:
        get_conda_env(),
    log:
        str(Path(data_path) / "logs/spladder_call_events/spladder_events_{event}.log")
    # wildcard_constraints were added to help snakemake to distinguish between events used in spladder_call_events and spladder_call_events_long_execution rules
    wildcard_constraints:
        event="[a-zA-Z0-9]+_[a-zA-Z0-9]+"  # stands for 2 parts separated by an underscore
    threads: get_rule_threads ('spladder_call_events')
    retries: get_rule_retries('spladder_call_events'), 
    resources:
        # lambda is used to trigger callable function to retrieve the "attempt" value
        # mem = lambda wildcards, attempt : get_memory_by_attempt_and_threads( \
        #                                 attempt_num = attempt, \
        #                                 memory = next((el for el in [get_config_value ('rules/spladder_call_events/memory')] if el is not None), default_rule_memory), \
        #                                 threads = next((el for el in [get_config_value ('rules/spladder_call_events/threads')] if el is not None), 8), \
        #                                 max_memory = next((el for el in [get_config_value ('resources/max_memory')] if el is not None), default_max_memory) \
        #                                 ),
        cl_resources = lambda wildcards, attempt : get_cluster_resources_by_attempt (attempt, 'spladder_call_events'),
        cl_job_suffix = lambda wildcards : wildcards.event,
        attempt = get_attempt,
        walltime = get_rule_walltime('spladder_call_events'),  # "48:00",
        # ase-edge-limit value is calculated on a fly for each attempt
        ase_edge_limit = lambda wildcards, attempt : get_ase_edge_limit_by_attempt( \
                                        attempt_num = attempt, \
                                        ase_edge_limit = next((el for el in [spladder_events_ase_edge_limit] if el is not None), 500), \
                                        ase_edge_limit_decrease_step = next((el for el in [spladder_events_ase_edge_limit_decrease_step] if el is not None), spladder_events_ase_edge_limit_decrease_step_default), \
                                        ase_edge_limit_min = next((el for el in [spladder_events_ase_edge_limit_min] if el is not None), spladder_events_ase_edge_limit_default) \
                                        ),
    params:
        out_dir = str(Path(data_path) / "spladder"),
        annotation = str(Path(alt_sl_ref_data_path) / ('spladder/' + motrpac_ref_data_genom_dir + '/genome.gtf')),
        debug = "-d" if debug else "",
        cur_event = "{event}",
    shell:
        '''
        spladder build \
            {params.debug} \
            -v \
            --parallel {threads} \
            -o {params.out_dir} \
            -a {params.annotation} \
            -b {input.alignments} \
            --event-types {params.cur_event} \
            --ase-edge-limit {resources.ase_edge_limit} \
            2>&1 | tee {log}
        '''
        
rule spladder_call_events_long_execution:
    input:
        bam_files = expand(str(Path(data_path) / "star_align/named_files/{sample}_Aligned.sortedByCoord.out.bam"), sample=samples_to_process),
        spladder_quantification2_status_file = str(Path(data_path) / "spladder/status_files/spladder_quantification2.txt"),
        alignments = str(Path(data_path) / "spladder/alignments.txt"),
    output:
        merge_graphs_confirmed_txt = str(Path(data_path) / "spladder/merge_graphs_{event}_C3.confirmed.txt.gz"),
        merge_graphs_confirmed_pickle = str(Path(data_path) / "spladder/merge_graphs_{event}_C3.confirmed.pickle"),
        merge_graphs_confirmed_gff3 = str(Path(data_path) / "spladder/merge_graphs_{event}_C3.confirmed.gff3"),
        merge_graphs_counts = str(Path(data_path) / "spladder/merge_graphs_{event}_C3.counts.hdf5"),
        merge_graphs_pickle = str(Path(data_path) / "spladder/merge_graphs_{event}_C3.pickle"),
    conda:
        get_conda_env(),
    log:
        str(Path(data_path) / "logs/spladder_call_events_long_execution/spladder_events_long_execution_{event}.log")
        #str(Path(data_path) / "logs/spladder_call_events_long_execution/spladder_events_long_execution_mult_exon_skip.log")
    # wildcard_constraints were added to help snakemake to distinguish between events used in spladder_call_events and spladder_call_events_long_execution rules
    wildcard_constraints:
        event="[a-zA-Z0-9]+_[a-zA-Z0-9]+_[a-zA-Z0-9]+"  # stands for 3 parts separated by an underscore
    threads: get_rule_threads ('spladder_call_events_long_execution')
    retries: get_rule_retries('spladder_call_events_long_execution'), 
    resources:
        # mem = next((el for el in [get_config_value ('rules/spladder_call_events/memory')] if el is not None), default_rule_memory),
        # lambda is used to trigger callable function to retrieve the "attempt" value
        # mem = lambda wildcards, attempt : get_memory_by_attempt_and_threads( \
        #                                 attempt_num = attempt, \
        #                                 memory = next((el for el in [get_config_value ('rules/spladder_call_events/memory')] if el is not None), default_rule_memory), \
        #                                 threads = next((el for el in [get_config_value ('rules/spladder_call_events/threads')] if el is not None), 8), \
        #                                 max_memory = next((el for el in [get_config_value ('resources/max_memory')] if el is not None), default_max_memory) \
        #                                 ),
        cl_resources = lambda wildcards, attempt : get_cluster_resources_by_attempt (attempt, 'spladder_call_events_long_execution'),
        cl_job_suffix = lambda wildcards : wildcards.event,
        attempt = get_attempt,
        # define walltime values to be used by LSF
#        walltime = get_walltime_by_samples_count( \
#                        num_samples = len(samples_to_process), \
#                        time_coef = get_config_value ('rules/spladder_call_events_long_execution/time_coef'), \
#                        min_walltime = get_config_value ('rules/spladder_call_events_long_execution/min_walltime'), \
#                        max_regular_walltime = get_config_value ('resources/max_regular_walltime'), \
#                        max_long_walltime = get_config_value ('resources/max_long_walltime')),
        # define name of the queue to be used with LSF
#        lsf_queue = get_queue_by_samples_count( \
#                        num_samples = len(samples_to_process), \
#                        time_coef = get_config_value ('rules/spladder_call_events_long_execution/time_coef'), \
#                        min_walltime = get_config_value ('rules/spladder_call_events_long_execution/min_walltime'), \
#                        max_regular_walltime = get_config_value ('resources/max_regular_walltime'), \
#                        max_long_walltime = get_config_value ('resources/max_long_walltime'), \
#                        normal_queue = get_config_value ('resources/regular_queue'), \
#                        long_queue = get_config_value ('resources/long_queue')),
        # ase-edge-limit value is calculated on a fly for each attempt
        ase_edge_limit = lambda wildcards, attempt : get_ase_edge_limit_by_attempt( \
                                        attempt_num = attempt, \
                                        ase_edge_limit = next((el for el in [spladder_events_ase_edge_limit] if el is not None), 500), \
                                        ase_edge_limit_decrease_step = next((el for el in [spladder_events_ase_edge_limit_decrease_step] if el is not None), spladder_events_ase_edge_limit_decrease_step_default), \
                                        ase_edge_limit_min = next((el for el in [spladder_events_ase_edge_limit_min] if el is not None), spladder_events_ase_edge_limit_default) \
                                        ),
    params:
        out_dir = str(Path(data_path) / "spladder"),
        annotation = str(Path(alt_sl_ref_data_path) / ('spladder/' + motrpac_ref_data_genom_dir + '/genome.gtf')),
        debug = "-d" if debug else "",
        cur_event = "{event}",
    shell:
        '''
        spladder build \
            {params.debug} \
            -v \
            --parallel {threads} \
            -o {params.out_dir} \
            -a {params.annotation} \
            -b {input.alignments} \
            --event-types {params.cur_event} \
            --ase-edge-limit {resources.ase_edge_limit} \
            2>&1 | tee {log}
        '''

rule chr_info:
    input:
        bam = str(Path(data_path) / "star_align/{sample}/Aligned.sortedByCoord.out.bam"),
        # star_align_valid_file = str(Path(data_path) / "star_align_validate/valid_{sample}"),
    output:
        chr_info = str(Path(data_path) / "star_align/{sample}/chr_info.txt"),
        tmp_dir = temp_debugcheck(directory(str(Path(data_path) / "star_align/{sample}/chr_info"))),
    conda:
        get_conda_env(),
    log:
        str(Path(data_path) / "logs/chr_info/chr_info_{sample}.log")
    resources:
        cl_job_suffix = lambda wildcards : wildcards.sample,
    params:
        F_argument = get_config_value ('rules/chr_info/F_argument') if get_config_value ('rules/chr_info/F_argument') else "0x900",
        bam2 = str(Path(data_path) / "star_align/{sample}/chr_info/{sample}_primary.bam"),
    shell:
        '''
        mkdir -p {output.tmp_dir}
        samtools view -b -F {params.F_argument} {input.bam} -o {params.bam2} 2>&1 | tee {log}
        samtools index {params.bam2} 2>&1 | tee -a {log}
        samtools idxstats {params.bam2} > {output.chr_info} 2>&1 | tee -a {log}
        '''

rule star_align_merge_all:
    input: 
        star_align_qc_log=expand(str(Path(data_path) / "star_align/{sample}/Log.final.out"), sample=samples_to_process),
    output:
        out_file = str(Path(data_path) / "star_align/star_QC.txt")
    log:
        str(Path(data_path) / "logs/star_align_merge_all/star_align_merge_all.log")
    # localrule: True
    resources:
        cl_resources = lambda wildcards, attempt : get_cluster_resources_by_attempt (attempt, 'star_align_merge_all'),
        # cl_job_suffix = '',
    params:
        # defines parameters for the merge_files_by_sample function
        samples = samples_to_process,  # list of samples to be processed
        star_align_path = str(Path(data_path) / "star_align"),  # location of the "star_align" dir
        star_align_sub_dir_path_map = '|sample|',  # defines path under the "star_align" dir to find needed file
        input_file_name = 'Log.final.out',  # name of the files to be merged
        pipeline_warning_file_path = pipeline_warning_file_path,
    script:
        "scripts/star_align_merge_all.py"

rule rRNA:
    input:
        # it must be the same input as for the star_aliqn rule
        trimR1 = get_R1_file_path(data_path, 'fastq_trim'),
        trimR2 = get_R2_file_path(data_path, 'fastq_trim') if samples_R1R2_present else []
    output:
        rRNA = str(Path(data_path) / "rRNA/{sample}.txt"),
        sam_file = temp_debugcheck(str(Path(data_path) / "rRNA/{sample}.sam")),
    conda:
        get_conda_env(),
    log:
        str(Path(data_path) / "logs/rRNA/rRNA_{sample}.log")
    threads: get_rule_threads ('rRNA')
    resources:
        # mem = next((el for el in [get_config_value ('rules/rRNA/memory')] if el is not None), default_rule_memory)  # 6000
        cl_resources = lambda wildcards, attempt : get_cluster_resources_by_attempt (attempt, 'rRNA'),
        cl_job_suffix = lambda wildcards : wildcards.sample,
    params:
        tool = "bowtie2",
        # to prepare ref_dir - it will parse motrpac_ref_data_genom_dir config value, replace any numbers, split by "_", get the first item from the list and concatenate with "_rRNA"
        # if genom dir is "hg38_gencode_v30", it will return "hg_rRNA"
        ref_dir = get_ref_index_dir(motrpac_ref_data_genom_dir, 'rRNA'),  # re.sub('[\d]','', motrpac_ref_data_genom_dir).split('_')[0] + "_rRNA",
        ref_path = str(Path(motrpac_ref_data_path) / "misc_data"),
        sam_file = str(Path(data_path) / "rRNA/{sample}.sam"),
        # the following creates a partial bash script that will be inserted to the shell part
        input_files_bash_script = "-1 " + get_R1_file_path(data_path, 'fastq_trim') + \
                           " -2 " + get_R2_file_path(data_path, 'fastq_trim') \
                           if samples_R1R2_present else \
                           "-U " + get_R1_file_path(data_path, 'fastq_trim'),
    shell:
        '''
        {params.tool} \
        -p {threads} \
        {params.input_files_bash_script} \
        -x {params.ref_path}/{params.ref_dir}/{params.ref_dir} \
        --local \
        -S {params.sam_file} 2>&1 | tee {log} {output.rRNA}
        '''
        
rule phix:
    input:
        # it must be the same input as for the star_aliqn rule
        trimR1 = get_R1_file_path(data_path, 'fastq_trim'),
        trimR2 = get_R2_file_path(data_path, 'fastq_trim') if samples_R1R2_present else []
    output:
        phix = str(Path(data_path) / "phix/{sample}.txt"),
        sam_file = temp_debugcheck(str(Path(data_path) / "phix/{sample}.sam"))
    conda:
        get_conda_env(),
    log:
        str(Path(data_path) / "logs/phix/phix_{sample}.log")
    threads: get_rule_threads ('phix')
    resources:
        # mem = next((el for el in [get_config_value ('rules/phix/memory')] if el is not None), default_rule_memory)  # 6000,
        cl_resources = lambda wildcards, attempt : get_cluster_resources_by_attempt (attempt, 'phix'),
        cl_job_suffix = lambda wildcards : wildcards.sample,
    params:
        tool = "bowtie2",
        ref_dir = "phix",
        ref_path = str(Path(motrpac_ref_data_path) / "misc_data"),
        sam_file = str(Path(data_path) / "phix/{sample}.sam"),
        # the following creates a partial bash script that will be inserted to the shell part
        input_files_bash_script = "-1 " + get_R1_file_path(data_path, 'fastq_trim') + \
                           " -2 " + get_R2_file_path(data_path, 'fastq_trim') \
                           if samples_R1R2_present else \
                           "-U " + get_R1_file_path(data_path, 'fastq_trim'),
    shell:
        '''
        {params.tool} \
        -p {threads} \
        {params.input_files_bash_script} \
        -x {params.ref_path}/{params.ref_dir}/{params.ref_dir} \
        --local \
        -S {params.sam_file} 2>&1 | tee {log} {output.phix}
        '''

rule globin:
    input:
        # it must be the same input as for the star_aliqn rule
        trimR1 = get_R1_file_path(data_path, 'fastq_trim'),
        trimR2 = get_R2_file_path(data_path, 'fastq_trim') if samples_R1R2_present else []
    output:
        globin = str(Path(data_path) / "globin/{sample}.txt"),
        sam_file = temp_debugcheck(str(Path(data_path) / "globin/{sample}.sam"))
    conda:
        get_conda_env(),
    log:
        str(Path(data_path) / "logs/globin/globin_{sample}.log")
    threads: get_rule_threads ('globin')
    resources:
        # mem = next((el for el in [get_config_value ('rules/globin/memory')] if el is not None), default_rule_memory)  # 6000
        cl_resources = lambda wildcards, attempt : get_cluster_resources_by_attempt (attempt, 'globin'),
        cl_job_suffix = lambda wildcards : wildcards.sample,
    params:
        tool = "bowtie2",
        # to prepare ref_dir - it will parse motrpac_ref_data_genom_dir config value, replace any numbers, split by "_", get the first item from the list and concatenate with "_globin"
        # if genom dir is "hg38_gencode_v30", it will return "hg_globin"
        ref_dir = get_ref_index_dir(motrpac_ref_data_genom_dir, 'globin'), 
        ref_path = str(Path(motrpac_ref_data_path) / "misc_data"),
        sam_file = str(Path(data_path) / "globin/{sample}.sam"),
        # the following creates a partial bash script that will be inserted to the shell part
        input_files_bash_script = "-1 " + get_R1_file_path(data_path, 'fastq_trim') + \
                           " -2 " + get_R2_file_path(data_path, 'fastq_trim') \
                           if samples_R1R2_present else \
                           "-U " + get_R1_file_path(data_path, 'fastq_trim'),
    shell:
        '''
        {params.tool} \
        -p {threads} \
        {params.input_files_bash_script} \
        -x {params.ref_path}/{params.ref_dir}/{params.ref_dir} \
        --local \
        -S {params.sam_file} 2>&1 | tee {log} {output.globin}
        '''

rule telescope_sam:
    input:
        # it must be the same input as for the star_aliqn rule
        trimR1 = get_R1_file_path(data_path, 'fastq_trim'),
        trimR2 = get_R2_file_path(data_path, 'fastq_trim') if samples_R1R2_present else []
    output:
        telescop_bowtie_txt = temp_debugcheck(str(Path(data_path) / "telescope/bam/{sample}_bowtie.txt")),
        sam_file = temp_debugcheck(str(Path(data_path) / "telescope/bam/{sample}.sam")),
    conda:
        get_conda_env(),
    log:
        str(Path(data_path) / "logs/telescope_sam/telescope_sam_{sample}.log")
    threads: get_rule_threads ('telescope_sam')
    resources:
        cl_resources = lambda wildcards, attempt : get_cluster_resources_by_attempt (attempt, 'telescope_sam'),  # required for loading memory 
        cl_job_suffix = lambda wildcards : wildcards.sample,
        walltime = get_rule_walltime('telescope_sam'),
    params:
        tool = "bowtie2",
        ref_dir = "bowtie2_index/genome", 
        ref_path = str(Path(motrpac_ref_data_path) / "hg38_gencode_v38"),
        # the following creates a partial bash script that will be inserted to the shell part
        input_files_bash_script = "-1 " + get_R1_file_path(data_path, 'fastq_trim') + \
                           " -2 " + get_R2_file_path(data_path, 'fastq_trim') \
                           if samples_R1R2_present else \
                           "-U " + get_R1_file_path(data_path, 'fastq_trim'),
    shell:
        '''
        {params.tool} \
        -p {threads} \
        {params.input_files_bash_script} \
        -x {params.ref_path}/{params.ref_dir} \
        --local \
        -S {output.sam_file} 2>&1 | tee {log} {output.telescop_bowtie_txt}
        '''

rule telescope_to_bam_bai:
    input:
        sam_file = str(Path(data_path) / "telescope/bam/{sample}.sam"),
    output:
        bam_file = temp_debugcheck(str(Path(data_path) / "telescope/bam/{sample}.bam")),
        sorted_bam = temp_debugcheck(str(Path(data_path) / "telescope/bam/{sample}.sorted.bam")),
        bai_file = temp_debugcheck(str(Path(data_path) / "telescope/bam/{sample}.sorted.bam.bai")),
    conda:
        get_conda_env(),
    log:
        str(Path(data_path) / "logs/telescope_to_bam_bai/telescope_to_bam_bai_{sample}.log")
    resources:
        cl_job_suffix = lambda wildcards : wildcards.sample,
    # params:
    shell:
        '''
        samtools view -b {input.sam_file} -o {output.bam_file} 2>&1 | tee {log}
        samtools sort {output.bam_file} -o {output.sorted_bam} 2>&1 | tee -a {log}
        samtools index {output.sorted_bam} 2>&1 | tee -a {log}
        '''

rule telescope:
    input:
        sorted_bam = str(Path(data_path) / "telescope/bam/{sample}.sorted.bam"),
        bai_file = str(Path(data_path) / "telescope/bam/{sample}.sorted.bam.bai"),
    output:
        telescope_TE_counts = str(Path(data_path) / "telescope/{sample}/telescope-TE_counts.tsv"),
        telescope_run_stats = temp_debugcheck(str(Path(data_path) / "telescope/{sample}/telescope-run_stats.tsv")),
        telescope_checkpoint = temp_debugcheck(str(Path(data_path) / "telescope/{sample}/telescope-checkpoint.npz")),

    conda:
        get_conda_env('telescope'),  # contains installation of telescope ran from a fork of main project with a fix of the issue
    log:
        str(Path(data_path) / "logs/telescope/telescope_{sample}_stdout.log")
    resources:
        cl_resources = lambda wildcards, attempt : get_cluster_resources_by_attempt (attempt, 'telescope'),  # required for loading memory 
        cl_job_suffix = lambda wildcards : wildcards.sample,
        walltime = get_rule_walltime('telescope'),
    params:
        ref_path = str(Path(motrpac_ref_data_path) / "misc_data"),
        ref_dir = "TE/HERV_rmsk.hg38.v2/transcripts.gtf", 
        wrk_dir = str(Path(data_path) / "telescope"),
        out_dir = "{sample}",
    shell:
        '''
        telescope assign \
        {input.sorted_bam} \
        {params.ref_path}/{params.ref_dir} \
        --max_iter 200 --theta_prior 200000 \
        --outdir {params.wrk_dir}/{params.out_dir} \
        2>&1 | tee {log}
        '''

rule telescope_merge_all_count:
    input: 
        telescope_TE_counts = expand(str(Path(data_path) / "telescope/{sample}/telescope-TE_counts.tsv"), sample=samples_to_process),
    output:
        out_file = str(Path(data_path) / "telescope_TE_counts_all.txt")
    conda:
        get_conda_env(),
    log:
        str(Path(data_path) / "logs/telescope_merge_all_count/telescope_merge_all_count.log")
    params:
        # defines parameters for the telescope_merge_all script
        samples = samples_to_process,  # list of samples to be processed
        collect_column = "count",  # name of the columns to retrieved from the file
        pipeline_warning_file_path = pipeline_warning_file_path,
    script:
        "scripts/telescope_merge_all.py"


# include rmats related rules only if rmats contrast metadata file is valid
# if rmats_contrast_file_valid and run_alternative_splicing:
rule star_align_most_common_read_length:
    input: 
        bam_files=expand(str(Path(data_path) / "star_align/{sample}/Aligned.sortedByCoord.out.bam"), sample=samples_to_process),
        star_align_valid_file = expand(str(Path(data_path) / "star_align_validate/valid_{sample}"), sample=samples_to_process),
    output:
        out_file = str(Path(data_path) / "star_align/most_common_read_length.txt")
    conda:
        get_conda_env(),
    log:
        str(Path(data_path) / "logs/star_align_most_common_read_length/star_align_most_common_read_length.log")
    resources:
        # mem = next((el for el in [get_config_value ('rules/star_align_most_common_read_length/memory')] if el is not None), default_rule_memory)  # 12000
        cl_resources = lambda wildcards, attempt : get_cluster_resources_by_attempt (attempt, 'star_align_most_common_read_length'),
    params:
        max_sample_entries = get_config_value ('rules/star_align_most_common_read_length/max_sample_entries'),
        pipeline_info_file_path = pipeline_info_file_path,
        pipeline_warning_file_path = pipeline_warning_file_path,
    script:
        "scripts/star_align_most_common_read_length.py"
    
rule qc53:
    input:
        bam=str(Path(data_path) / "star_align/{sample}/Aligned.sortedByCoord.out.bam"),
        # star_align_valid_file = str(Path(data_path) / "star_align_validate/valid_{sample}"),
    output:
        qc53 = str(Path(data_path) / "qc53/{sample}.RNA_Metrics"),
    conda:
        get_conda_env(),
    log:
        str(Path(data_path) / "logs/qc53/qc53_{sample}.log")
    retries: get_rule_retries('qc53'), 
    resources:
        # mem = next((el for el in [get_config_value ('rules/qc53/memory')] if el is not None), default_rule_memory)  # 20000
        # lambda is used to trigger callable function to retrieve the "attempt" value
        # mem = lambda wildcards, attempt : get_memory_by_attempt_and_threads( \
        #                                 attempt_num = attempt, \
        #                                 memory = next((el for el in [get_config_value ('rules/qc53/memory')] if el is not None), default_rule_memory), \
        #                                 max_memory = next((el for el in [get_config_value ('resources/max_memory')] if el is not None), default_max_memory) \
        #                                 ),
        cl_resources = lambda wildcards, attempt : get_cluster_resources_by_attempt (attempt, 'qc53'),
        cl_job_suffix = lambda wildcards : wildcards.sample,
        attempt = get_attempt,
    params:
        tool = "java",
        # identify location of the picard.jar file based on the current installation of the picard tool
        # picard = str(Path(os.path.realpath(shutil.which ('picard'))).parent.absolute()/"picard.jar"),
        command = "CollectRnaSeqMetrics",
        MINIMUM_LENGTH = get_config_value ('rules/qc53/MINIMUM_LENGTH'),
        ref_flat = str(Path(motrpac_ref_data_path) / (motrpac_ref_data_genom_dir + '/qc53_ref/genome_flat.txt')),
        strand = "FIRST_READ_TRANSCRIPTION_STRAND", # The MOP specify FIRST_READ_TRANSCRIPTION_STRAND
        RRNA_FRAGMENT_PERCENTAGE = get_config_value ('rules/qc53/RRNA_FRAGMENT_PERCENTAGE'),
    shell:
        '''
        # get location of picard.jar file at the time of execution, so it will be retrieved in the proper environment
        picard="$(dirname "$(realpath "$(command -v picard)")")/picard.jar"
        
        {params.tool} \
        -jar $picard {params.command} \
        I={input.bam} \
        O={output}\
        MINIMUM_LENGTH={params.MINIMUM_LENGTH} \
        REF_FLAT={params.ref_flat} \
        STRAND_SPECIFICITY={params.strand} \
        RRNA_FRAGMENT_PERCENTAGE={params.RRNA_FRAGMENT_PERCENTAGE} \
        2>&1 | tee {log}
        '''
        
rule mark_dup:
    input:
        bam=str(Path(data_path) / "star_align/{sample}/Aligned.sortedByCoord.out.bam"),
        # star_align_valid_file = str(Path(data_path) / "star_align_validate/valid_{sample}"),
    output:
        mark_dup_metrics = str(Path(data_path) / "mark_dup/{sample}.dup_metrics"),
        mark_dup_bam = temp_debugcheck(str(Path(data_path) / "mark_dup/{sample}_markedDup.bam")),
        # mark_dup_bai = str(Path(data_path) / "mark_dup/{sample}_markedDup.bai"),  # this file is not being created if parameter params.CREATE_INDEX = false
    conda:
        get_conda_env(),
    log:
        str(Path(data_path) / "logs/mark_dup/mark_dup_{sample}.log"),
    resources:
        # mem = next((el for el in [get_config_value ('rules/mark_dup/memory')] if el is not None), default_rule_memory)  # 60000
        cl_resources = lambda wildcards, attempt : get_cluster_resources_by_attempt (attempt, 'mark_dup'),
        cl_job_suffix = lambda wildcards : wildcards.sample,
    params:
        tool = "java",
        # picard = str(Path(os.path.realpath(shutil.which ('picard'))).parent.absolute()/"picard.jar"),
        command = "MarkDuplicates",
        # if cluster config was supplied, it will use the memory allocation from the config file, otherwise it won't be set
        memory_allocation = get_config_value ('rules/mark_dup/memory_allocation'), 
        CREATE_INDEX = get_config_value ('rules/mark_dup/CREATE_INDEX'),
        VALIDATION_STRINGENCY = get_config_value ('rules/mark_dup/VALIDATION_STRINGENCY'),
        ASSUME_SORT_ORDER = get_config_value ('rules/mark_dup/ASSUME_SORT_ORDER'),
        REMOVE_DUPLICATES = get_config_value ('rules/mark_dup/REMOVE_DUPLICATES'),
    shell:
        '''
        # get location of picard.jar file at the time of execution, so it will be retrieved in the proper environment
        picard="$(dirname "$(realpath "$(command -v picard)")")/picard.jar"
        
        {params.tool} \
        {params.memory_allocation} \
        -jar $picard {params.command} \
        I={input.bam} \
        O={output.mark_dup_bam}\
        CREATE_INDEX={params.CREATE_INDEX} \
        VALIDATION_STRINGENCY={params.VALIDATION_STRINGENCY} \
        ASSUME_SORT_ORDER={params.ASSUME_SORT_ORDER} \
        M={output.mark_dup_metrics} \
        REMOVE_DUPLICATES={params.REMOVE_DUPLICATES} \
        2>&1 | tee {log}
        '''

rule multiqc_fastqc_raw:
    input:
        multiqc_fastqc_raw=expand(get_R1_file_path(data_path, sub_dir_path = 'fastqc_raw', file_name_format = '{sample}_|FN|_fastqc.html'), sample=samples_to_process),
    output:
        multiqc_fastqc_raw_html = str(Path(data_path) / "multiqc/fastqc_raw.html"),
        multiqc_fastqc_raw_txt = str(Path(data_path) / "multiqc/fastqc_raw_data/multiqc_fastqc.txt"),
    conda:
        get_conda_env('multiqc'), 
    log:
        str(Path(data_path) / "logs/multiqc/fastqc_raw.log"),
    retries: get_rule_retries('multiqc_fastqc_raw'), 
    resources:
        cl_resources = lambda wildcards, attempt : get_cluster_resources_by_attempt (attempt, 'multiqc_fastqc_raw'),
        attempt = get_attempt,
        # cl_job_suffix = '',
    params:
        tool = "multiqc",
        filename = "fastqc_raw",
        outdir = str(Path(data_path) / "multiqc"),
        inputdir = str(Path(data_path) / "fastqc_raw"),
    shell:
        '''
        {params.tool} \
        --dirs --force \
        --filename {params.filename} \
        --outdir {params.outdir} \
        {params.inputdir} \
        2>&1 | tee {log}
        '''

rule multiqc_fastqc_trim:
    input:
        multiqc_fastqc_trim=expand(get_R1_file_path(data_path, sub_dir_path = 'fastqc_trim', file_name_format = '{sample}_|FN|_fastqc.html'), sample=samples_to_process),
        trim_log_file = expand(str(Path(data_path) / "logs/trim/tool_logs/log.{sample}"), sample=samples_to_process),
    output:
        multiqc_fastqc_trim_html = str(Path(data_path) / "multiqc/fastqc_trim.html"),
        multiqc_fastqc_trim_txt = str(Path(data_path) / "multiqc/fastqc_trim_data/multiqc_fastqc.txt"),
        multiqc_fastqc_cutadapt_txt = str(Path(data_path) / "multiqc/fastqc_trim_data/multiqc_cutadapt.txt"),
    conda:
        get_conda_env('multiqc'),  # multiqc  multiqc_fastqc_trim
    log:
        str(Path(data_path) / "logs/multiqc/fastqc_trim.log"),
    retries: get_rule_retries('multiqc_fastqc_trim'), 
    resources:
        # lambda is used to trigger callable function to retrieve the "attempt" value
        # mem = lambda wildcards, attempt : get_memory_by_attempt_and_threads( \
        #                                 attempt_num = attempt, \
        #                                 memory = next((el for el in [get_config_value ('rules/multiqc_fastqc_trim/memory')] if el is not None), default_rule_memory), \
        #                                 max_memory = next((el for el in [get_config_value ('resources/max_memory')] if el is not None), default_max_memory) \
        #                                 ),
        cl_resources = lambda wildcards, attempt : get_cluster_resources_by_attempt (attempt, 'multiqc_fastqc_trim'),
        attempt = get_attempt,
        # cl_job_suffix = '',
    params:
        tool = "multiqc",
        filename = "fastqc_trim",
        outdir = str(Path(data_path) / "multiqc"),
        inputdir1 = str(Path(data_path) / "fastqc_trim"),
        inputdir2 = str(Path(data_path) / "logs/trim/tool_logs"),
        ignore_files = '--ignore "*_cluster.out" --ignore "*_cluster.err"'
    shell:
        '''
        {params.tool} \
        --dirs --force \
        --filename {params.filename} \
        --outdir {params.outdir} \
        {params.inputdir1} {params.inputdir2} \
        {params.ignore_files} \
        2>&1 | tee {log}
        '''

rule multiqc_post_align:
    input:
        star_align=expand(str(Path(data_path) / "star_align/{sample}/Log.final.out"), sample=samples_to_process),
        rsem_genes = expand(str(Path(data_path) / "rsem/{sample}.genes.results"), sample=samples_to_process) if run_rsem else [],
        rsem_isoforms = expand(str(Path(data_path) / "rsem/{sample}.isoforms.results"), sample=samples_to_process) if run_rsem else [],
        featureCounts=expand(str(Path(data_path) / "featureCounts/{sample}"), sample=samples_to_process) if run_feature_counts else [],
        featureCounts_summary=expand(str(Path(data_path) / "featureCounts/{sample}.summary"), sample=samples_to_process) if run_feature_counts else [],
    output:
        multiqc_post_align_html = str(Path(data_path) / "multiqc/post_align.html"),
    conda:
        get_conda_env('multiqc'),
    log:
        str(Path(data_path) / "logs/multiqc/post_align.log"),
    resources:
        # cl_job_suffix = '',
    params:
        tool = "multiqc",
        filename = "post_align",
        outdir = str(Path(data_path) / "multiqc"),
        inputdir1 = str(Path(data_path) / "star_align"),
        inputdir2 = str(Path(data_path) / "rsem") if run_rsem else "",
        inputdir3 = str(Path(data_path) / "featureCounts") if run_feature_counts else "",
    shell:
        '''
        {params.tool} \
        --dirs --force \
        --filename {params.filename} \
        --outdir {params.outdir} \
        {params.inputdir1} {params.inputdir2} {params.inputdir3}\
        2>&1 | tee {log}
        '''

rule post_align_qc53:
    input:
        qc53_RNA_Metrics=expand(str(Path(data_path) / "qc53/{sample}.RNA_Metrics"), sample=samples_to_process),
    output:
        post_align_qc53_txt = str(Path(data_path) / "qc53.txt"),
        post_align_qc53_plot = str(Path(data_path) / "qc53.pdf"),
    conda:
        get_conda_env(),
    log:
        str(Path(data_path) / "logs/post_align_qc53/post_align_qc53.log"),
    params:
        samples = samples_to_process,
        output_file_name = "qc53"
    script:
        "scripts/qc53_as.R"

rule rsem_validate:
    input:
        rsem_bam=str(Path(data_path) / "star_align/{sample}/Aligned.toTranscriptome.out.bam"),
    output:
        valid_file = temp_debugcheck(str(Path(data_path) / "rsem_validate/valid_{sample}")),
    conda:
        get_conda_env(),
    log:
        str(Path(data_path) / "logs/rsem_validate/rsem_validate_{sample}.log")
    retries: get_rule_retries('rsem_validate'),
    resources:
        # mem = next((el for el in [get_config_value ('rules/rsem_validate/memory')] if el is not None), default_rule_memory)  # 12000
        # lambda is used to trigger callable function to retrieve the "attempt" value
        # mem = lambda wildcards, attempt : get_memory_by_attempt_and_threads( \
        #                                 attempt_num = attempt, \
        #                                 memory = next((el for el in [get_config_value ('rules/rsem_validate/memory')] if el is not None), default_rule_memory), \
        #                                 max_memory = next((el for el in [get_config_value ('resources/max_memory')] if el is not None), default_max_memory) \
        #                                 ),
        cl_resources = lambda wildcards, attempt : get_cluster_resources_by_attempt (attempt, 'rsem_validate'),
        cl_job_suffix = lambda wildcards : wildcards.sample,
        attempt = get_attempt,
    params:
        outdir = directory(str(Path(data_path) / "rsem_validate")),
        pipeline_warning_file_path = pipeline_warning_file_path,
    shell:
        '''
        rsem_valid_output=$(rsem-sam-validator "{input.rsem_bam}" 2>&1)
        echo $rsem_valid_output 2>&1 | tee {log}
        
        if [[ $rsem_valid_output == *"The input file is valid"* ]]; then
            # File is valid
            output_file="valid_{wildcards.sample}"
            warning_msg=""
        else
            # File is not valid
            output_file="invalid_{wildcards.sample}"
            warning_msg="Warning: rsem_validate step. Sample {wildcards.sample} is invalid - $rsem_valid_output"
        fi
        out_path={params.outdir}/$output_file
        echo "output file name: $out_path" 2>&1 | tee -a {log}
        
        # Create the output file with the appropriate name and content
        echo "$rsem_valid_output" > "$out_path" # 2>&1 | tee -a {log}
        
        if [ -n "$warning_msg" ]; then
            # If it's not blank, add warning message to the pipeline warning file
            echo $warning_msg >> {params.pipeline_warning_file_path}
        fi
        '''

rule rsem:
    input:
        rsem_bam=str(Path(data_path) / "star_align/{sample}/Aligned.toTranscriptome.out.bam"),
        rsem_valid_file = str(Path(data_path) / "rsem_validate/valid_{sample}"),
    output:
        rsem_outdir = temp_debugcheck(directory(str(Path(data_path) / "rsem/{sample}.stat"))),
        rsem_genes = str(Path(data_path) / "rsem/{sample}.genes.results"),
        rsem_isoforms = str(Path(data_path) / "rsem/{sample}.isoforms.results"),
    conda:
        get_conda_env(),
    log:
        str(Path(data_path) / "logs/rsem/rsem_{sample}.log")
    threads: get_rule_threads ('rsem')
    retries: get_rule_retries('rsem'),
    resources:
        # mem = next((el for el in [get_config_value ('rules/rsem/memory')] if el is not None), default_rule_memory)  # 6000
        # lambda is used to trigger callable function to retrieve the "attempt" value
        # mem = lambda wildcards, attempt : get_memory_by_attempt_and_threads( \
        #                                 attempt_num = attempt, \
        #                                 memory = next((el for el in [get_config_value ('rules/rsem/memory')] if el is not None), default_rule_memory), \
        #                                 threads = next((el for el in [get_config_value ('rules/rsem/threads')] if el is not None), 1), \
        #                                 max_memory = next((el for el in [get_config_value ('resources/max_memory')] if el is not None), default_max_memory) \
        #                                 ),
        cl_resources = lambda wildcards, attempt : get_cluster_resources_by_attempt (attempt, 'rsem'),
        cl_job_suffix = lambda wildcards : wildcards.sample,
        attempt = get_attempt,
        walltime = get_rule_walltime('rsem'),  # "72:00",
    params:
        tool = "rsem-calculate-expression",
        pairopt = "--paired-end" if samples_R1R2_present else "",
        genomeDir = str(Path(motrpac_ref_data_path) / (motrpac_ref_data_genom_dir + '/rsem_index/genome')),
        sample = "{sample}",
        tempdir = directory(str(Path(data_path) / "rsem/{sample}_tmpdir")),
    shell:
        '''
        {params.tool} \
        {params.pairopt} \
        --num-threads {threads} \
        --no-bam-output \
        --forward-prob 0.5 \
        --seed 12345 \
        --bam {input.rsem_bam} \
        --temporary-folder {params.tempdir} \
        {params.genomeDir} \
        rsem/{params.sample} \
        2>&1 | tee {log}
        '''

rule rsem_tpm:
    input:
        isoforms_file = str(Path(data_path) / "rsem/{sample}.isoforms.results"),
    output:
        isoforms_tpm_file = str(Path(data_path) / "rsem/{sample}.isoforms.results_tpm"),
    log:
        str(Path(data_path) / "logs/rsem_tpm/rsem_tpm_{sample}.log")
    resources:
        cl_job_suffix = lambda wildcards : wildcards.sample,
    params:
        keep_columns = [0,5], # represents the first and last (6th) columns of the file 
        first_row_value = "{sample}",
        # pipeline_info_file_path = pipeline_info_file_path,
        pipeline_warning_file_path = pipeline_warning_file_path,
    script:
        "scripts/rsem_tpm.py"

rule suppa_rsem:
        input:
            rsem_tpm_file = str(Path(data_path) / "rsem/{sample}.isoforms.results_tpm"),
        output:
            event_file = str(Path(data_path) / "suppa/rsem/{sample}/{sample}_{event}.psi"),
        conda:
            get_conda_env(),
        log:
            stdout = str(Path(data_path) / "logs/suppa_rsem/suppa_{sample}_{event}.out"),
            stderr = str(Path(data_path) / "logs/suppa_rsem/suppa_{sample}_{event}.err"),
        resources:
            # mem = next((el for el in [get_config_value ('rules/suppa_rsem/memory')] if el is not None), default_rule_memory)  # 10000
            cl_resources = lambda wildcards, attempt : get_cluster_resources_by_attempt (attempt, 'suppa_rsem'),
            cl_job_suffix = lambda wildcards : wildcards.sample + '_' + wildcards.event,
        params:
            # identify location of the suppa.py in the bin directory of conda env (using python to get conda's location)
            # tool = str(os.path.realpath(Path(shutil.which ('python')).parent.absolute()/"suppa.py")), 
            command = get_config_value ('rules/suppa_rsem/command'), # "psiPerEvent"
            suppa_event_index = str(Path(alt_sl_ref_data_path) / \
                    (get_config_value ('rules/suppa_rsem/index_dir') \
                    + '/' + motrpac_ref_data_genom_dir \
                    + '/' + motrpac_ref_data_genom_dir + get_config_value ('rules/suppa_rsem/event_index_name_map'))),
            output = str(Path(data_path) / "suppa/rsem/{sample}/{sample}_{event}")
        shell:
            '''
            # get location of suppa file at the time of execution, so it will be retrieved in the proper environment
            tool=$(realpath "$(dirname "$(command -v python)")/suppa.py")
            
            python $tool {params.command} --ioe-file {params.suppa_event_index} --expression-file {input.rsem_tpm_file} -o {params.output} \
            > >(tee {log.stdout}) 2> {log.stderr} # this redirectds stdout to the "out" file and screen, and stderr to err file only.
            '''

rule rsem_merge_all_count:
    input: 
        rsem_gene_files = expand(str(Path(data_path) / "rsem/{sample}.genes.results"), sample=samples_to_process),
    output:
        out_file = str(Path(data_path) / "rsem_genes_count.txt")
    conda:
        get_conda_env(),
    log:
        str(Path(data_path) / "logs/rsem_merge_all_count/rsem_merge_all_count.log")
    params:
        # defines parameters for the merge_files_by_gene_id function
        samples = samples_to_process,  # list of samples to be processed
        collect_column = "expected_count",  # name of the columns to retrieved from the file
        pipeline_warning_file_path = pipeline_warning_file_path,
    script:
        "scripts/rsem_merge_all.py"

rule rsem_merge_all_tpm:
    input: 
        rsem_gene_files = expand(str(Path(data_path) / "rsem/{sample}.genes.results"), sample=samples_to_process),
    output:
        out_file = str(Path(data_path) / "rsem_genes_tpm.txt")
    conda:
        get_conda_env(),
    log:
        str(Path(data_path) / "logs/rsem_merge_all_tpm/rsem_merge_all_tpm.log")
    params:
        # defines parameters for the merge_files_by_gene_id function
        samples = samples_to_process,  # list of samples to be processed
        collect_column = "TPM",  # name of the columns to retrieved from the file
        pipeline_warning_file_path = pipeline_warning_file_path,
    script:
        "scripts/rsem_merge_all.py"
        
rule rsem_merge_all_fpkm:
    input: 
        rsem_gene_files = expand(str(Path(data_path) / "rsem/{sample}.genes.results"), sample=samples_to_process),
    output:
        out_file = str(Path(data_path) / "rsem_genes_fpkm.txt")
    conda:
        get_conda_env(),
    log:
        str(Path(data_path) / "logs/rsem_merge_all_fpkm/rsem_merge_all_fpkm.log")
    params:
        # defines parameters for the merge_files_by_gene_id function
        samples = samples_to_process,  # list of samples to be processed
        collect_column = "FPKM",  # name of the columns to retrieved from the file
        pipeline_warning_file_path = pipeline_warning_file_path,
    script:
        "scripts/rsem_merge_all.py"

rule featureCounts:
    input:
        bam=str(Path(data_path) / "star_align/{sample}/Aligned.sortedByCoord.out.bam"),
        # star_align_valid_file = str(Path(data_path) / "star_align_validate/valid_{sample}"),
    output:
        sam = temp_debugcheck(str(Path(data_path) / "featureCounts/{sample}_tmpdir/{sample}.sam")),
        sam_sorted= temp_debugcheck(str(Path(data_path) / "featureCounts/{sample}_tmpdir/{sample}_sorted.sam")),
        sam_full = temp_debugcheck(str(Path(data_path) / "featureCounts/{sample}_tmpdir/{sample}_full.sam")),
        tempdir = temp_debugcheck(directory(str(Path(data_path) / "featureCounts/{sample}_tmpdir"))),
        featureCounts=str(Path(data_path) / "featureCounts/{sample}"),
        featureCounts_summary=str(Path(data_path) / "featureCounts/{sample}.summary"),
    conda:
        get_conda_env(),
    log:
        str(Path(data_path) / "logs/featureCounts/featureCounts_{sample}.log")
    threads: get_rule_threads ('featureCounts')
    resources:
        # mem = next((el for el in [get_config_value ('rules/featureCounts/memory')] if el is not None), default_rule_memory)  # 6000
        cl_resources = lambda wildcards, attempt : get_cluster_resources_by_attempt (attempt, 'featureCounts'),
        cl_job_suffix = lambda wildcards : wildcards.sample,
    params:
        pairopt="-p" if samples_R1R2_present else "",
        genome_gtf_file=str(Path(motrpac_ref_data_path) / (motrpac_ref_data_genom_dir + '/genome.gtf')),
    shell:
        '''
        mkdir -p {output.tempdir}
        samtools view -@ {threads} {input.bam} >{output.sam} 2>&1 | tee {log}
        sleep 5
        export LC_ALL=C; export LC_LANG=C
        sort {output.sam} -k1,1 -o {output.sam_sorted} -T {output.tempdir} --parallel={threads} -S 5G 2>&1 | tee -a {log}
        sleep 5
        samtools view -H {input.bam} >{output.sam_full} 2>&1 | tee -a {log}
        sleep 5
        cat {output.sam_sorted} >>{output.sam_full} 2>&1 | tee -a {log}
        sleep 5
        featureCounts -T {threads} --tmpDir {output.tempdir}  -a {params.genome_gtf_file} -o {output.featureCounts} {params.pairopt} -M --fraction {output.sam_full} --donotsort 2>&1 | tee -a {log}
        '''

rule featureCounts_merge_all:
    input: 
        feature_counts_gene_files = expand(str(Path(data_path) / "featureCounts/{sample}"), sample=samples_to_process),
    output:
        out_file = str(Path(data_path) / "featureCounts.txt")
    conda:
        get_conda_env(),
    log:
        str(Path(data_path) / "logs/featureCounts_merge_all/featureCounts_merge_all.log")
    # localrule: True
    params:
        # defines parameters for the merge_files_by_gene_id function
        samples = samples_to_process,  # list of samples to be processed
        collect_column_num = 6,  # 0-indexed column number to be collected (7th column)
        pipeline_warning_file_path = pipeline_warning_file_path,
    script:
        "scripts/feature_counts_merge_all.py"

if samples_R1I1_present:
    rule UMI_dup:
        input:
            bam=str(Path(data_path) / "star_align/{sample}/Aligned.sortedByCoord.out.bam"),
            # star_align_valid_file = str(Path(data_path) / "star_align_validate/valid_{sample}"),
        output:
            # dup = str(Path(data_path) / "star_align/{sample}_dup/{sample}_dup_log.txt"),
            dup = str(Path(data_path) / "star_align/{sample}/dup_log.txt"),
            wrk_dir = temp_debugcheck(directory(str(Path(data_path) / "star_align/{sample}/umi_dup"))),
            tempdir = temp_debugcheck(directory(str(Path(data_path) / "star_align/{sample}/umi_dup/umi_dup_tmp"))),
        log:
            str(Path(data_path) / "logs/UMI_dup/{sample}_UMI_dup.log")
        conda:
            get_conda_env('python2'),  # python2_env_name
            # "conda_envs/ngscheckmate_python2.yml"
        threads: get_rule_threads ('UMI_dup')
        resources:
            # mem = next((el for el in [get_config_value ('rules/UMI_dup/memory')] if el is not None), default_rule_memory)  # 36000
            cl_resources = lambda wildcards, attempt : get_cluster_resources_by_attempt (attempt, 'UMI_dup'),
            cl_job_suffix = lambda wildcards : wildcards.sample,
        params:
            # to use with python 3
            # python = "python",
            # tool = str(Path(snakemake_path) / "scripts/nugen/nudup_py3.py"),
            # to use with python 2
            # python = python2_path,
            python = "python",
            tool = str(Path(snakemake_path) / "scripts/nugen/nudup_py2.py"),
            pairopt = "-2" if samples_R1R2_present else "",
            outopt = "{sample}",
            dup_log_tmp = str(Path(data_path) / "star_align/{sample}/umi_dup/{sample}_dup_log.txt"),
        shell:
            '''
            set +e # head might has a problem with the following command, avoiding interruptions
            len=$(samtools view {input.bam} |head -1 |awk '{{umi=gensub("^.*:","",1,$1); print length(umi)}}')
            set -e
            echo "len=$len"
            
            mkdir -p {output.tempdir}
            cd {output.wrk_dir}
            {params.python} {params.tool} {params.pairopt} -s $len -l $len --rmdup-only --out {params.outopt} -T {output.tempdir} {input.bam} 2>&1 | tee {log}
            mv {params.dup_log_tmp} {output.dup}
            '''

rule ngscheckmate_prep:
    input:
        bam=str(Path(data_path) / "star_align/{sample}/Aligned.sortedByCoord.out.bam"),
        star_align_valid_file = str(Path(data_path) / "star_align_validate/valid_{sample}"),
    output:
        out_file = str(Path(data_path) / "NGSCheckMate/bam/{sample}.bam")
    log:
        str(Path(data_path) / "logs/ngscheckmate_prep/ngscheckmate_prep_{sample}.log")
    # localrule: True
    resources:
        cl_resources = lambda wildcards, attempt : get_cluster_resources_by_attempt (attempt, 'ngscheckmate_prep'),
        cl_job_suffix = lambda wildcards : wildcards.sample,
    params:
        relative = get_config_value ('rules/ngscheckmate_prep/relative'),
        pipeline_warning_file_path = pipeline_warning_file_path,
    script:
        "scripts/ngscheckmate_prep.py"

rule ngscheckmate_bai:
    input:
        ngs_prep_bam = str(Path(data_path) / "NGSCheckMate/bam/{sample}.bam")
    output:
        ngs_bai = str(Path(data_path) / "NGSCheckMate/bam/{sample}.bam.bai")
    conda:
        get_conda_env(),
    log:
        str(Path(data_path) / "logs/ngscheckmate_bai/ngscheckmate_bai_{sample}.log")
    resources:
        cl_job_suffix = lambda wildcards : wildcards.sample,
    shell:
        '''
        samtools index {input.ngs_prep_bam}
        '''

rule ngscheckmate_vcf:
    input:
        ngs_prep_bam = str(Path(data_path) / "NGSCheckMate/bam/{sample}.bam"),
        ngs_bai = str(Path(data_path) / "NGSCheckMate/bam/{sample}.bam.bai"),
    output:
        ngs_vcf = str(Path(data_path) / "NGSCheckMate/vcf/{sample}.vcf")
    log:
        str(Path(data_path) / "logs/ngscheckmate_vcf/ngscheckmate_vcf_{sample}.log")
    resources:
        cl_job_suffix = lambda wildcards : wildcards.sample,
    conda: 
        get_conda_env('python2'),  # python2_env_name
    params:
        genome_ga = str(Path(motrpac_ref_data_path) / (motrpac_ref_data_genom_dir + '/genome.fa')),
        snp_bed = "SNP_GRCh38_hg38_wChr.bed",
    shell:
        '''
        snp_bed=$(dirname "$(which python)")/../NGSCheckMate/SNP/{params.snp_bed} # identify location of the snp.bed file in the conda's env installation
        samtools mpileup -I -uf {params.genome_ga} -l $snp_bed {input.ngs_prep_bam} | bcftools call -c - > {output.ngs_vcf} # version taking snp_bed from the conda's env installation
        '''

rule ngscheckmate:
    input:
        vcf_files = expand(str(Path(data_path) / "NGSCheckMate/vcf/{sample}.vcf"), sample=samples_to_process),
    output:
        matched = str(Path(data_path) / "NGSCheckMate/output_matched.txt"),
        matrix = str(Path(data_path) / "NGSCheckMate/output_corr_matrix.txt"),
        all = str(Path(data_path) / "NGSCheckMate/output_all.txt"),
        r_script_r = temp_debugcheck(str(Path(data_path) / "NGSCheckMate/r_script.r")),
        r_script_r_Rout = temp_debugcheck(str(Path(data_path) / "NGSCheckMate/r_script.r.Rout")),
        # RData = temp_debugcheck(str(Path(data_path) / "NGSCheckMate/.RData")),
    log:
        str(Path(data_path) / "logs/ngscheckmate/ngscheckmate.log")
    conda:
        get_conda_env('python2'),  # python2_env_name
    params:
        snp_bed = "SNP_GRCh38_hg38_wChr.bed",
        nsm = "ncm.py",
        out_dir = str(Path(data_path) / "NGSCheckMate"),
    shell:
        '''
        cd {params.out_dir}
        snp_bed=$(dirname "$(which python)")/../NGSCheckMate/SNP/{params.snp_bed}
        nsm=$(dirname "$(which python)")/../NGSCheckMate/{params.nsm}
        python $nsm -V -d vcf -bed $snp_bed -O {params.out_dir}
        '''

rule ngscheckmate_parse_output:
    input:
        output_matched = str(Path(data_path) / "NGSCheckMate/output_matched.txt"),
    output:
        incoherent_groups_out = str(Path(data_path) / "NGSCheckMate/incoherent_groups_out.json"),
    conda:
        get_conda_env(),
    log:
        str(Path(data_path) / "logs/ngscheckmate_parse_output/ngscheckmate_parse_output.log"),
    localrule: True  # this step uses connection to the SQL Server which is currently not supported from the cluster
    params:
        debug = True if debug or create_rds_for_rscripts else False,
        pipeline_warning_file_path = pipeline_warning_file_path,
        env_var_db_driver = get_config_value ('database/connection/env_db_driver'),  # "SNAKEMAKE_DB_DRIVER",
        env_var_db_server = get_config_value ('database/connection/env_db_server'),  # "SNAKEMAKE_DB_SERVER",
        env_var_db_database = get_config_value ('database/connection/env_db_name'),  # "SNAKEMAKE_DB_DATABASE",
        env_var_db_user = get_config_value ('database/connection/env_db_user_name'), # "SNAKEMAKE_DB_USER",
        env_var_db_pwd = get_config_value ('database/connection/env_db_user_pwd'),   # "SNAKEMAKE_DB_PWD",
        include_r_script = str(Path(workflow.basedir) / "scripts/DBfunctions.R"),
    script:
        "scripts/parseCheckmateOutput.R"

rule qc_final:
    input:
        star_align_merge_all = str(Path(data_path) / "star_align/star_QC.txt"),
        multiqc_fastqc_raw_txt = str(Path(data_path) / "multiqc/fastqc_raw_data/multiqc_fastqc.txt"),
        multiqc_fastqc_trim_txt = str(Path(data_path) / "multiqc/fastqc_trim_data/multiqc_fastqc.txt"),
        multiqc_fastqc_cutadapt_txt = str(Path(data_path) / "multiqc/fastqc_trim_data/multiqc_cutadapt.txt"),
        post_align_qc53_txt = str(Path(data_path) / "qc53.txt"),
        rRNA = expand(str(Path(data_path) / "rRNA/{sample}.txt"), sample=samples_to_process),
        phix = expand(str(Path(data_path) / "phix/{sample}.txt"), sample=samples_to_process),
        globin = expand(str(Path(data_path) / "globin/{sample}.txt"), sample=samples_to_process),
        chr_info = expand(str(Path(data_path) / "star_align/{sample}/chr_info.txt"), sample=samples_to_process),
        mark_dup_metrics = expand(str(Path(data_path) / "mark_dup/{sample}.dup_metrics"), sample=samples_to_process),
        umi_dup = expand(str(Path(data_path) / "star_align/{sample}/dup_log.txt"), sample=samples_to_process) if samples_R1I1_present else [],
        trim_log = expand(str(Path(data_path) / "logs/trim/tool_logs/log.{sample}"), sample=samples_to_process),
        
    output:
        qc_info_csv = str(Path(data_path) / "qc_info.csv"),
        # rds_file = str(Path(data_path) / "snakemake.RData") if debug or create_rds_for_rscripts else [],
    conda:
        get_conda_env(),
    log:
        str(Path(data_path) / "logs/qc_final/qc_final.log"),
    # localrule: True
    params:
        rds_file = str(Path(data_path) / "snakemake.RData"),
        samples = samples_to_process,
        debug = True if debug or create_rds_for_rscripts else False,
        pipeline_warning_file_path = pipeline_warning_file_path,
        samples_R1I1_present = samples_R1I1_present,
        rRNA = str(Path(data_path) / "rRNA"),
        phix = str(Path(data_path) / "phix"),
        globin = str(Path(data_path) / "globin"),
        chr_info = str(Path(data_path) / "star_align"),
        mark_dup_metrics = str(Path(data_path) / "mark_dup"),
        umi_dup = str(Path(data_path) / "star_align") if samples_R1I1_present else "",
        trim_log = str(Path(data_path) / "logs/trim/tool_logs"),
        samples_R1R2_present = samples_R1R2_present,
    script:
        "scripts/qc.R"
