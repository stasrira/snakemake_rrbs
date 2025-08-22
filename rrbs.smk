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

create_rds_for_rscripts = str_to_bool(replace_cfg_with_environment('create_rds_for_rscripts', get_config_value ('create_rds_for_rscripts')))
print_prerun_info = str_to_bool(replace_cfg_with_environment('print_prerun_info', get_config_value ('print_prerun_info')))

data_path = replace_cfg_with_environment('data_path', get_config_value ('data_path'))
data_path_required_folders = get_config_value ('data_path_required_folders')
raw_data_file_ext = replace_cfg_with_environment('raw_data_file_ext', get_config_value ('raw_data_file_ext'))

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

add_pre_run_info ('create_rds_for_rscripts = {} (if set to True, RDS files (being used for manual testing) will be created for all participating R scripts that are executed at the current run)'.format(create_rds_for_rscripts))
add_pre_run_info ('data_path = {}'.format(data_path))
add_pre_run_info ('data_path_required_folders = {}'.format(data_path_required_folders))
add_pre_run_info ('raw_data_file_ext = {}'.format(raw_data_file_ext))

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

# ruleorder: trim>trim_umi

rule all:
    input:
        # fastq_R1R2_validate = expand(str(Path(data_path) / "fastq_raw_validate/valid_{sample}"), sample=samples_to_process) if samples_R1R2_present else [],
        # fastq_I1_validate = expand(str(Path(data_path) / "fastq_raw_validate/valid_{sample}_I1"), sample=samples_to_process) if samples_R1I1_present else [],
        # umi_attach_R1 = expand(str(Path(data_path) / "fastq_attach/{sample}_R1.fastq.gz"), sample=samples_to_process) if samples_R1I1_present else [],
        # umi_attach_R2 = expand(str(Path(data_path) / "fastq_attach/{sample}_R2.fastq.gz"), sample=samples_to_process) if samples_R1R2_present and samples_R1I1_present else [],
        # trim_normalized_R1 = expand(get_R1_file_path(data_path, 'fastq_trim'), sample=samples_to_process),
        bismark_work_bam = expand(str(Path(data_path) / "bismark_work/{sample}_R1_bismark_bt2_pe.bam") if samples_R1R2_present else str(Path(data_path) / "bismark_work/{sample}_R1_bismark_bt2.bam"), sample=samples_to_process),
        bismark_sort_bam = expand(str(Path(data_path) / "bismark/{sample}_R1_bismark_bt2_pe.bam") if samples_R1R2_present else str(Path(data_path) / "bismark/{sample}_R1_bismark_bt2.bam"), sample=samples_to_process),
        bismark_lambda = expand(str(Path(data_path) / "bismark_lambda/{sample}_R1_bismark_bt2_pe.bam") if samples_R1R2_present else str(Path(data_path) / "bismark_lambda/{sample}_R1_bismark_bt2.bam"), sample=samples_to_process),
        
        # bismark_pre_dup_bam_umi = expand(str(Path(data_path) / "bismark_umi/{sample}_R1_bismark_bt2_pe.bam") if samples_R1R2_present else str(Path(data_path) / "bismark_umi/{sample}_R1_bismark_bt2.bam"), sample=samples_to_process),
        # bismark_dedup_txt = expand(str(Path(data_path) / "bismark/{sample}_dedup.txt"), sample=samples_to_process),
        # bismark_methylation_extractor_ok_file = expand(str(Path(data_path) / "bismark/{sample}_bismark_methylation_extractor_ok.txt"), sample=samples_to_process),
        bismark_methylation_extractor_CHG_context = \
            expand( \
            str(Path(data_path) / "bismark/CHG_context_{sample}_R1_bismark_bt2_pe.deduplicated.txt") \
            if samples_R1R2_present else 
            str(Path(data_path) / "bismark/CHG_context_{sample}_R1_bismark_bt2.deduplicated.txt"), \
            sample=samples_to_process),
        bismark_report_html = expand(str(Path(data_path) / "bismark/{sample}_report.html"), sample=samples_to_process),
        bismark_chr_info = expand(str(Path(data_path) / "chr_info/{sample}_chr_info.txt"), sample=samples_to_process),
        bismark_summary_html = str(Path(data_path) / "bismark/bismark_summary_report.html"),
        bismark_summary_txt = str(Path(data_path) / "bismark/bismark_summary_report.txt"),
        bismark_lambda_summary_html = str(Path(data_path) / "bismark_lambda/bismark_summary_report.html"),
        bismark_lambda_summary_txt = str(Path(data_path) / "bismark_lambda/bismark_summary_report.txt"),
        bismark_4strand_report_ = str(Path(data_path) / "bismark/bismark_4strand.tsv"),
        bismark_lambda_4strand_report_ = str(Path(data_path) / "bismark_lambda/bismark_4strand.tsv"),
        
        fastqc_raw = expand(get_R1_file_path(data_path, sub_dir_path = 'fastqc_raw', file_name_format = '{sample}_|FN|_fastqc.html'), sample=samples_to_process),
        fastqc_trim_html = expand(str(Path(data_path) / "fastqc_trim/{sample}_R1_fastqc.html"), sample=samples_to_process),
        phix = expand(str(Path(data_path) / "phix/{sample}.txt"), sample=samples_to_process),
        multiqc_fastqc_raw_html = str(Path(data_path) / "multiqc/fastqc_raw.html"),
        
        # multiqc_fastqc_trim_html = str(Path(data_path) / "multiqc/fastqc_trim.html"),
        # multiqc_post_align = str(Path(data_path) / "multiqc/post_align.html"),

        # qc_final = str(Path(data_path) / "qc_info.csv"),

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
        get_conda_env('trim_galore'),
    log:
        str(Path(data_path) / "logs/fastqc_raw/fastqc_raw_{sample}.log")
    resources:
        cl_resources = lambda wildcards, attempt : get_cluster_resources_by_attempt (attempt, 'fastqc_raw'),
        cl_job_suffix = lambda wildcards : wildcards.sample,
    params:
        tool = "fastqc",
        out_dir = str(Path(data_path) / "fastqc_raw")
    shell:
        '''
        {params.tool} --outdir {params.out_dir} {input.fileR1} {input.fileR2} 2>&1 | tee {log}
        '''
        
rule bismark_work:
    input:
        fileR1 = get_R1_file_path(data_path, 'fastq_trim'),
        fileR2 = get_R2_file_path(data_path, 'fastq_trim') if samples_R1R2_present else [],
    output:
        # ok_file = str(Path(data_path) / "bismark_work/{sample}_ok.txt")
        tmp_dir = temp_debugcheck(directory(str(Path(data_path) / "bismark_work/tmp_{sample}"))),
        bam = temp_debugcheck(str(Path(data_path) / "bismark_work/{sample}_R1_bismark_bt2_pe.bam") if samples_R1R2_present else str(Path(data_path) / "bismark_work/{sample}_R1_bismark_bt2.bam")),
        align_report= temp_debugcheck(str(Path(data_path) / "bismark_work/{sample}_R1_bismark_bt2_PE_report.txt") if samples_R1R2_present else str(Path(data_path) / "bismark_work/{sample}_R1_bismark_bt2_SE_report.txt")),
    conda:
        get_conda_env('bismark'),
    log:
        str(Path(data_path) / "logs/bismark_work/bismark_work_{sample}.log")
    threads: get_rule_threads ('bismark_work')
    resources:
        cl_resources = lambda wildcards, attempt : get_cluster_resources_by_attempt (attempt, 'bismark_work'),  # required for loading memory 
        cl_job_suffix = lambda wildcards : wildcards.sample,
        walltime = get_rule_walltime('bismark_work'),
    params:
        genomeDir = str(Path(motrpac_ref_data_path) / (motrpac_ref_data_genom_dir)),
        # tmp_dir = str(Path(data_path) / "bismark_work/tmp_{sample}"),
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

rule bismark_lambda:
    input:
        fileR1 = get_R1_file_path(data_path, 'fastq_trim'),
        fileR2 = get_R2_file_path(data_path, 'fastq_trim') if samples_R1R2_present else [],
    output:
        # ok_file = str(Path(data_path) / "bismark_lambda/{sample}_ok.txt"),
        tmp_dir = temp_debugcheck(directory(str(Path(data_path) / "bismark_lambda/tmp_{sample}"))),
        bam = str(Path(data_path) / "bismark_lambda/{sample}_R1_bismark_bt2_pe.bam") if samples_R1R2_present else str(Path(data_path) / "bismark_lambda/{sample}_R1_bismark_bt2.bam"),
        align_report= str(Path(data_path) / "bismark_lambda/{sample}_R1_bismark_bt2_PE_report.txt") if samples_R1R2_present else str(Path(data_path) / "bismark_lambda/{sample}_R1_bismark_bt2_SE_report.txt"),
    conda:
        get_conda_env('bismark'),
    log:
        str(Path(data_path) / "logs/bismark_lambda/bismark_lambda_{sample}.log")
    threads: get_rule_threads ('bismark_lambda')
    resources:
        cl_resources = lambda wildcards, attempt : get_cluster_resources_by_attempt (attempt, 'bismark_lambda'),  # required for loading memory 
        cl_job_suffix = lambda wildcards : wildcards.sample,
        walltime = get_rule_walltime('bismark_lambda'),
    params:
        genomeDir = str(Path(motrpac_ref_data_path) / 'misc_data/lambda'),
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
        # echo 'OK' > {output.ok_file}
        
rule bismark_sort:
    input:
        bam = str(Path(data_path) / "bismark_work/{sample}_R1_bismark_bt2_pe.bam") if samples_R1R2_present else str(Path(data_path) / "bismark_work/{sample}_R1_bismark_bt2.bam"),
        align_report= str(Path(data_path) / "bismark_work/{sample}_R1_bismark_bt2_PE_report.txt") if samples_R1R2_present else str(Path(data_path) / "bismark_work/{sample}_R1_bismark_bt2_SE_report.txt"),
    output:
        tmp_dir = temp_debugcheck(directory(str(Path(data_path) / "bismark_work/tmp_{sample}"))),
        sam = temp_debugcheck(str(Path(data_path) / "bismark_work/{sample}.sam")),
        bam_out = temp_debugcheck(str(Path(data_path) / "bismark/{sample}_R1_bismark_bt2_pe.bam") if samples_R1R2_present else str(Path(data_path) / "bismark/{sample}_R1_bismark_bt2.bam")),
        align_report= str(Path(data_path) / "bismark/{sample}_R1_bismark_bt2_PE_report.txt") if samples_R1R2_present else str(Path(data_path) / "bismark/{sample}_R1_bismark_bt2_SE_report.txt"),
        
        # ok_file = str(Path(data_path) / "bismark/{sample}_bismark_sort_ok.txt"),
    conda:
        get_conda_env('bismark'),
    log:
        str(Path(data_path) / "logs/bismark_sort/bismark_sort_{sample}.log")
    threads: get_rule_threads ('bismark_sort')
    resources:
        cl_resources = lambda wildcards, attempt : get_cluster_resources_by_attempt (attempt, 'bismark_sort'),
        walltime = get_rule_walltime('bismark_sort'),
        cl_job_suffix = lambda wildcards : wildcards.sample,
    params:
        # tmp_dir = temp_debugcheck(directory(str(Path(data_path) / "bismark_work/tmp_{sample}"))),
        # sam = temp_debugcheck(str(Path(data_path) / "bismark_work/{sample}.sam")),
        # bam_out = str(Path(data_path) / "bismark/{sample}_R1_bismark_bt2_pe.bam") if samples_R1R2_present else str(Path(data_path) / "bismark/{sample}_R1_bismark_bt2.bam"),
    shell:
        '''
        mkdir -p {output.tmp_dir}
        # sort bam file
        samtools view -@ {threads} -H {input.bam} > {output.sam}  2>&1 | tee {log}
        samtools view -@ {threads} {input.bam} |sort -k5,5nr -k1,1 -s -S30G -T {output.tmp_dir} >> {output.sam} 2>&1 | tee -a {log}
        samtools view -@ {threads} -b {output.sam} -o {output.bam_out} 2>&1 | tee -a {log}
        
        # move/copy align report to the main folder (bismark)
        cp {input.align_report} {output.align_report} 
        '''
 
# this rule performed only if I1 files are present
rule bismark_umi:
    input: 
        bam_bismark = str(Path(data_path) / "bismark/{sample}_R1_bismark_bt2_pe.bam") if samples_R1R2_present else str(Path(data_path) / "bismark/{sample}_R1_bismark_bt2.bam"),
    output:
        pre_dup_bam_umi = \
            temp_debugcheck( \
            str(Path(data_path) / "bismark_umi/{sample}_R1_bismark_bt2_pe.bam") \
            if samples_R1R2_present else \
            str(Path(data_path) / "bismark_umi/{sample}_R1_bismark_bt2.bam")),
        sam = temp_debugcheck(str(Path(data_path) / "bismark_umi/{sample}_R1_bismark_bt2_pe.sam") if samples_R1R2_present else str(Path(data_path) / "bismark_umi/{sample}_R1_bismark_bt2.sam")),
    conda:
        get_conda_env('bismark'),
    log:
        str(Path(data_path) / "logs/bismark_umi/bismark_umi_{sample}.log"),
    resources:
        cl_resources = lambda wildcards, attempt : get_cluster_resources_by_attempt (attempt, 'bismark_umi'),
        cl_job_suffix = lambda wildcards : wildcards.sample,
    params:
        sam = str(Path(data_path) / "bismark/{sample}_R1_bismark_bt2_pe.sam") if samples_R1R2_present else str(Path(data_path) / "bismark/{sample}_R1_bismark_bt2.sam"),
        umi_script = str(Path(snakemake_path) / "scripts/bismark_bam_UMI_format.awk"),
    shell:
        '''
        samtools view -h {input.bam_bismark} | {params.umi_script} > {output.sam} 2>&1 | tee {log}
        samtools view -b -o {output.pre_dup_bam_umi} {output.sam} 2>&1 | tee -a {log}
        '''

rule bismark_dedup:
    input: 
        # input vary based on the presense of the I1 and R2 files
        bam_bismark = \
            (str(Path(data_path) / "bismark_umi/{sample}_R1_bismark_bt2_pe.bam") if samples_R1R2_present else str(Path(data_path) / "bismark_umi/{sample}_R1_bismark_bt2.bam")) \
            if samples_R1I1_present else \
            (str(Path(data_path) / "bismark/{sample}_R1_bismark_bt2_pe.bam") if samples_R1R2_present else str(Path(data_path) / "bismark/{sample}_R1_bismark_bt2.bam")),
    output:
        dedup_bam = \
            str(Path(data_path) / "bismark/{sample}_R1_bismark_bt2_pe.deduplicated.bam") \
            if samples_R1R2_present else \
            str(Path(data_path) / "bismark/{sample}_R1_bismark_bt2.deduplicated.bam"),
        dedup_report = \
            (str(Path(data_path) / "bismark/{sample}_R1_bismark_bt2_pe.deduplication_report.txt") \
            if samples_R1R2_present else \
            str(Path(data_path) / "bismark/{sample}_R1_bismark_bt2.deduplication_report.txt")),
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
            str(Path(data_path) / "bismark/{sample}_R1_bismark_bt2_pe.deduplicated.bam") \
            if samples_R1R2_present else \
            str(Path(data_path) / "bismark/{sample}_R1_bismark_bt2.deduplicated.bam"),
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
            str(Path(data_path) / "bismark/{sample}_R1_bismark_bt2_PE_report.txt") \
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
    shell:
        '''
        cd "$(dirname "{output.html}")"  # get to the parent dir of the output file
        
        bismark2report -o "{output.html}" -a "{input.align_report}" 2>&1 | tee {log}
        '''
        # echo 'OK' > {output.ok_file}

rule bismark_summary:
    input: 
        align_report = \
            expand(str(Path(data_path) / "bismark/{sample}_R1_bismark_bt2_PE_report.txt") \
            if samples_R1R2_present else \
            str(Path(data_path) / "bismark/{sample}_R1_bismark_bt2_SE_report.txt"), sample=samples_to_process),
    output:
        # ok_file = str(Path(data_path) / "bismark/bismark_summary_ok.txt"),
        txt = str(Path(data_path) / "bismark/bismark_summary_report.txt"),
        html = str(Path(data_path) / "bismark/bismark_summary_report.html"),
    conda:
        get_conda_env('bismark'),
    log:
        str(Path(data_path) / "logs/bismark_summary/bismark_summary.log"),
    # threads: get_rule_threads ('bismark_summary')
    resources:
        cl_resources = lambda wildcards, attempt : get_cluster_resources_by_attempt (attempt, 'bismark_summary'),
        walltime = get_rule_walltime('bismark_summary'),
        # cl_job_suffix = lambda wildcards : wildcards.sample,
    shell:
        '''
        cd "$(dirname "{output.html}")"  # get to the parent dir of the input file
        
        bismark2summary 2>&1 | tee {log}
        '''

rule bismark_lambda_summary:
    input: 
        align_report = \
            expand(str(Path(data_path) / "bismark_lambda/{sample}_R1_bismark_bt2_PE_report.txt") \
            if samples_R1R2_present else \
            str(Path(data_path) / "bismark_lambda/{sample}_R1_bismark_bt2_SE_report.txt"), sample=samples_to_process),
    output:
        # ok_file = str(Path(data_path) / "bismark/bismark_summary_ok.txt"),
        txt = str(Path(data_path) / "bismark_lambda/bismark_summary_report.txt"),
        html = str(Path(data_path) / "bismark_lambda/bismark_summary_report.html"),
    conda:
        get_conda_env('bismark'),
    log:
        str(Path(data_path) / "logs/bismark_lambda_summary/bismark_lambda_summary.log"),
    # threads: get_rule_threads ('bismark_summary')
    resources:
        cl_resources = lambda wildcards, attempt : get_cluster_resources_by_attempt (attempt, 'bismark_summary'),
        walltime = get_rule_walltime('bismark_summary'),
        # cl_job_suffix = lambda wildcards : wildcards.sample,
    shell:
        '''
        cd "$(dirname "{output.html}")"  # get to the parent dir of the input file
        
        bismark2summary 2>&1 | tee {log}
        '''

rule bismark_4strand_report:
    input: 
        align_reports = \
            expand(str(Path(data_path) / "bismark/{sample}_R1_bismark_bt2_PE_report.txt") \
            if samples_R1R2_present else \
            str(Path(data_path) / "bismark/{sample}_R1_bismark_bt2_SE_report.txt"), sample=samples_to_process),
    output:
        out_file = str(Path(data_path) / "bismark/bismark_4strand.tsv"),
    conda:
        get_conda_env(),
    log:
        str(Path(data_path) / "logs/bismark_4strand_report/bismark_4strand_report.log"),
    resources:
        # mem = next((el for el in [get_config_value ('rules/bismark_4strand_report/memory')] if el is not None), default_rule_memory)  # 10000
        cl_resources = lambda wildcards, attempt : get_cluster_resources_by_attempt (attempt, 'bismark_4strand_report'),
        # cl_job_suffix = lambda wildcards : wildcards.sample,
    params:
        samples = samples_to_process,
        rule_name = 'bismark_4strand_report',  # required for reporting error from within the script
        pipeline_warning_file_path = pipeline_warning_file_path,
    script:
        "scripts/bismark_4strand.py"

rule bismark_lambda_4strand_report:
    input: 
        align_reports = \
            expand(str(Path(data_path) / "bismark_lambda/{sample}_R1_bismark_bt2_PE_report.txt") \
            if samples_R1R2_present else \
            str(Path(data_path) / "bismark_lambda/{sample}_R1_bismark_bt2_SE_report.txt"), sample=samples_to_process),
    output:
        out_file = str(Path(data_path) / "bismark_lambda/bismark_4strand.tsv"),
    conda:
        get_conda_env(),
    log:
        str(Path(data_path) / "logs/bismark_lambda_4strand_report/bismark_lambda_4strand_report.log"),
    resources:
        cl_resources = lambda wildcards, attempt : get_cluster_resources_by_attempt (attempt, 'bismark_lambda_4strand_report'),
        # cl_job_suffix = lambda wildcards : wildcards.sample,
    params:
        samples = samples_to_process,
        rule_name = 'bismark_lambda_4strand_report',  # required for reporting error from within the script
        pipeline_warning_file_path = pipeline_warning_file_path,
    script:
        "scripts/bismark_4strand.py"

rule bismark_chr_info:
    input: 
        dedup_bam = \
            str(Path(data_path) / "bismark/{sample}_R1_bismark_bt2_pe.deduplicated.bam") \
            if samples_R1R2_present else \
            str(Path(data_path) / "bismark/{sample}_R1_bismark_bt2.deduplicated.bam"),
    output:
        chr_info = str(Path(data_path) / "chr_info/{sample}_chr_info.txt"),
    conda:
        get_conda_env('bismark'),
    log:
        str(Path(data_path) / "logs/bismark_chr_info/bismark_chr_info_{sample}.log"),
    threads: get_rule_threads ('bismark_chr_info')
    resources:
        cl_resources = lambda wildcards, attempt : get_cluster_resources_by_attempt (attempt, 'bismark_chr_info'),
        walltime = get_rule_walltime('bismark_chr_info'),
        cl_job_suffix = lambda wildcards : wildcards.sample,
    shell:
        '''
        samtools view {input.dedup_bam}|cut -f 3|sort |uniq -c > {output.chr_info} 2>&1 | tee {log}
        '''
        # samtools view $bam|cut -f 3|sort |uniq -c > chr_info/${SID}.txt
        # echo 'OK' > {output.ok_file}

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
    shell:
        '''
        zcat {input.fileR1} | {params.umi_script} -v Ifq={params.fileI1} > {output.tmp_umi_attach} 2>&1 | tee {log}
        gzip -c {output.tmp_umi_attach} > {output.attach_fastq} 2>&1 | tee -a {log}
        '''

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
        # ok_file = str(Path(data_path) / "fastq_trim_temp/{sample}_ok.txt") if samples_R1I1_present else str(Path(data_path) / "fastq_trim/{sample}_ok.txt"),
        trim_R1 = temp_debugcheck(str(Path(data_path) / "fastq_trim_temp/{sample}_R1_val_1.fq.gz")) if samples_R1I1_present else str(Path(data_path) / "fastq_trim/{sample}_R1_val_1.fq.gz"), 
        report_R1 = str(Path(data_path) / "fastq_trim_temp/{sample}_R1.fastq.gz_trimming_report.txt") if samples_R1I1_present else str(Path(data_path) / "fastq_trim/{sample}_R1.fastq.gz_trimming_report.txt"),
        trim_R2 = (temp_debugcheck(str(Path(data_path) / "fastq_trim_temp/{sample}_R2_val_2.fq.gz")) if samples_R1I1_present else str(Path(data_path) / "fastq_trim/{sample}_R2_val_2.fq.gz")) if samples_R1R2_present else [],
        report_R2 = (str(Path(data_path) / "fastq_trim_temp/{sample}_R2.fastq.gz_trimming_report.txt") if samples_R1I1_present else str(Path(data_path) / "fastq_trim/{sample}_R2.fastq.gz_trimming_report.txt"))if samples_R1R2_present else [],
    conda:
        get_conda_env('trim_galore'),
    log:
        str(Path(data_path) / "logs/trim/{sample}.log")
    resources:
        cl_job_suffix = lambda wildcards : wildcards.sample,
    params:
        index_adapter = get_config_value ('rules/trim/index_adapter'),

        # the following creates a partial bash script that will be inserted to the shell part
        # note: it adds an addtional argument and adjusts the index2_adapter
        pair_options = \
            "--paired -a2 {}".format( \
            get_config_value ('rules/trim/index2_adapter') if samples_R1I1_present else get_config_value ('rules/trim/index_adapter')) \
            if samples_R1R2_present else "",
            
        # rename non-paired output of ".gz" to the format of the paired output 
        rename_single_trim = 'mv {} {}'.format( \
            str(Path(data_path) / "fastq_trim_temp/{sample}_R1_trimmed.fq.gz") if samples_R1I1_present else str(Path(data_path) / "fastq_trim/{sample}_R1_trimmed.fq.gz"), \
            str(Path(data_path) / "fastq_trim_temp/{sample}_R1_val_1.fq.gz") if samples_R1I1_present else str(Path(data_path) / "fastq_trim/{sample}_R1_val_1.fq.gz") \
            ) if not samples_R1R2_present else "",
        # move_report_R1 = 'mv {} {}'.format()
    shell:
        '''
        # this dir will define the detination of the output of the trim_galore
        trim_dir="$(dirname "{output.trim_R1}")"  # get trim dir from the path of the output file
        
        trim_galore -a {params.index_adapter} {params.pair_options} -o $trim_dir {input.fileR1} {input.fileR2} 2>&1 | tee {log}

        # rename non-paired output file in the next line (can be blank if R1 and R2 files are present)
        {params.rename_single_trim}
        
        '''
        # echo 'OK' > {output.ok_file}

rule trim_umi:
    input:
        trim_R1 = str(Path(data_path) / "fastq_trim_temp/{sample}_R1_val_1.fq.gz"),
        trim_R2 = str(Path(data_path) / "fastq_trim_temp/{sample}_R2_val_2.fq.gz") if samples_R1R2_present else [],
    output:
        # ok_file = str(Path(data_path) / "fastq_trim_temp/{sample}_trim_umi_ok.txt"),
        trim_R1 = temp_debugcheck(str(Path(data_path) / "fastq_trim/{sample}_R1_val_1.fq.gz")),
        trim_R2 = temp_debugcheck(str(Path(data_path) / "fastq_trim/{sample}_R2_val_2.fq.gz")) if samples_R1R2_present else [],
        trim_report_R1 = str(Path(data_path) / "fastq_trim/{sample}_R1.fastq.gz_trimming_report.txt"),
        trim_report_R2 = str(Path(data_path) / "fastq_trim/{sample}_R2.fastq.gz_trimming_report.txt") if samples_R1R2_present else [],
    log:
        str(Path(data_path) / "logs/trim_umi/trim_umi_{sample}.log")
    conda:
        get_conda_env('python2'),  # python2_env_name
    threads: get_rule_threads ('trim_umi')
    resources:
        cl_resources = lambda wildcards, attempt : get_cluster_resources_by_attempt (attempt, 'trim_umi'),
        cl_job_suffix = lambda wildcards : wildcards.sample,
    params:
        python = "python",
        tool = str(Path(snakemake_path) / "scripts/nugen/trimRRBSdiversityAdaptCustomers.py"),
        pairopt = "-2" if samples_R1R2_present else "",
        
        # trimmed fastq files that needs to be moved one level up to the trim folder
        trim_umi_R1 = str(Path(data_path) / "fastq_trim_temp/{sample}_R1_val_1.fq_trimmed.fq.gz"),
        trim_umi_R2 = str(Path(data_path) / "fastq_trim_temp/{sample}_R2_val_2.fq_trimmed.fq.gz") if samples_R1R2_present else [],
        
        # trimmed report files that needs to be moved one level up to the trim folder
        trim_report_R1 = str(Path(data_path) / "fastq_trim_temp/{sample}_R1.fastq.gz_trimming_report.txt"),
        trim_report_R2 = str(Path(data_path) / "fastq_trim_temp/{sample}_R2.fastq.gz_trimming_report.txt") if samples_R1R2_present else [],
    shell:
        '''
        cd "$(dirname "{input.trim_R1}")"  # get to the parent dir of the input file
        
        {params.python} {params.tool} -1 {input.trim_R1} {params.pairopt} {input.trim_R2} 2>&1 | tee {log}
        
        # get file names that needs to be moved to trim folder (one level up)
        trim_umi_R1="{params.trim_umi_R1}"
        trim_umi_R2="{params.trim_umi_R2}"
        trim_report_R1="{params.trim_report_R1}"
        trim_report_R2="{params.trim_report_R2}"

        # validate that file variables not blank and specified file is present and then move it one folder up
        [ -n "$trim_umi_R1" ] && [ -f "$trim_umi_R1" ] && mv "$trim_umi_R1" "{output.trim_R1}"
        [ -n "$trim_umi_R2" ] && [ -f "$trim_umi_R2" ] && mv "$trim_umi_R2" "{output.trim_R2}"
        [ -n "$trim_report_R1" ] && [ -f "$trim_report_R1" ] && mv "$trim_report_R1" "{output.trim_report_R1}"
        [ -n "$trim_report_R2" ] && [ -f "$trim_report_R2" ] && mv "$trim_report_R2" "{output.trim_report_R2}"
        '''
        # echo 'OK' > {output.ok_file}
        # mv {params.trim_umi_R1} {output.trim_umi_R1}
        # {params.mv_pairopt} {params.trim_umi_R2} {output.trim_umi_R2}
        # [ -n "$trim_report_R1" ] && [ -f "$trim_report_R1" ] && mv "$trim_report_R1" "$(dirname "$trim_report_R1")/.."
        # [ -n "$trim_report_R2" ] && [ -f "$trim_report_R2" ] && mv "$trim_report_R2" "$(dirname "$trim_report_R2")/.."

rule trim_name_normalized:
    input:
        trim_R1 = str(Path(data_path) / "fastq_trim/{sample}_R1_val_1.fq.gz"),
        trim_R2 = str(Path(data_path) / "fastq_trim/{sample}_R2_val_2.fq.gz") if samples_R1R2_present else [],
    output: 
        fileR1 = get_R1_file_path(data_path, 'fastq_trim'),
        fileR2 = get_R2_file_path(data_path, 'fastq_trim') if samples_R1R2_present else [],
    conda:
        get_conda_env(),
    log:
        str(Path(data_path) / "logs/trim_name_normalized/trim_name_normalized_{sample}.log")
    resources:
        cl_resources = lambda wildcards, attempt : get_cluster_resources_by_attempt (attempt, 'trim_name_normalized'),
        cl_job_suffix = lambda wildcards : wildcards.sample,
    params:
        rename_R2 = 'cp' if samples_R1R2_present else '',
        # rename_R2 = 'mv' if samples_R1R2_present else '',
    shell:
        '''
        # always rename R1 file
        cp {input.trim_R1} {output.fileR1}
        # rename R2 only if variable samples_R1R2_present is True
        {params.rename_R2} {input.trim_R2} {output.fileR2}
        '''
        # mv {input.trim_R1} {output.fileR1}
        
rule fastqc_trim:
    input:
        # trim_R1 = str(Path(data_path) / "fastq_trim/{sample}_R1.fq.gz"),
        # trim_R2 = str(Path(data_path) / "fastq_trim/{sample}_R2.fq.gz") if samples_R1R2_present else [],
        trim_R1 = get_R1_file_path(data_path, 'fastq_trim'),
        trim_R2 = get_R2_file_path(data_path, 'fastq_trim') if samples_R1R2_present else [],
    output:
        outR1 = str(Path(data_path) / "fastqc_trim/{sample}_R1_fastqc.zip"),
        outR1_html = str(Path(data_path) / "fastqc_trim/{sample}_R1_fastqc.html"),
        outR2 = str(Path(data_path) / "fastqc_trim/{sample}_R2_fastqc.zip") if samples_R1R2_present else [],
        outR2_html = str(Path(data_path) / "fastqc_trim/{sample}_R2_fastqc.html") if samples_R1R2_present else [],
    conda:
        get_conda_env('trim_galore'),
    log:
        str(Path(data_path) / "logs/fastqc_trim/fastqc_trim_{sample}.log")
    resources:
        cl_resources = lambda wildcards, attempt : get_cluster_resources_by_attempt (attempt, 'fastqc_trim'),
        cl_job_suffix = lambda wildcards : wildcards.sample,
    params:
        tool = "fastqc",
        out_dir = str(Path(data_path) / "fastqc_trim")
    shell:
        '''
        {params.tool} --outdir {params.out_dir} {input.trim_R1} {input.trim_R2} 2>&1 | tee {log}
        '''

rule phix:
    input:
        # trim_R1 = str(Path(data_path) / "fastq_trim/{sample}_R1.fq.gz"),
        # trim_R2 = str(Path(data_path) / "fastq_trim/{sample}_R2.fq.gz") if samples_R1R2_present else [],
        trim_R1 = get_R1_file_path(data_path, 'fastq_trim'),
        trim_R2 = get_R2_file_path(data_path, 'fastq_trim') if samples_R1R2_present else [],
    output:
        phix = str(Path(data_path) / "phix/{sample}.txt"),
        sam_file = temp_debugcheck(str(Path(data_path) / "phix/{sample}.sam"))
    conda:
        get_conda_env(),
    log:
        str(Path(data_path) / "logs/phix/phix_{sample}.log")
    threads: get_rule_threads ('phix')
    resources:
        cl_resources = lambda wildcards, attempt : get_cluster_resources_by_attempt (attempt, 'phix'),
        cl_job_suffix = lambda wildcards : wildcards.sample,
    params:
        tool = "bowtie2",
        ref_dir = "phix",
        ref_path = str(Path(motrpac_ref_data_path) / "misc_data"),
        # sam_file = str(Path(data_path) / "phix/{sample}.sam"),
        # the following creates a partial bash script that will be inserted to the shell part
        input_files_bash_script = "-1 " + str(Path(data_path) / "fastq_trim/{sample}_R1.fastq.gz") + \
                           " -2 " + str(Path(data_path) / "fastq_trim/{sample}_R2.fastq.gz") \
                           if samples_R1R2_present else \
                           "-U " + str(Path(data_path) / "fastq_trim/{sample}_R1.fastq.gz"),
    shell:
        '''
        {params.tool} \
        -p {threads} \
        {params.input_files_bash_script} \
        -x {params.ref_path}/{params.ref_dir}/{params.ref_dir} \
        --local \
        -S {output.sam_file} 2>&1 | tee {log} {output.phix}
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
"""
rule multiqc_post_align:
    input:
        star_align=expand(str(Path(data_path) / "star_align/{sample}/Log.final.out"), sample=samples_to_process),
        # rsem_genes = expand(str(Path(data_path) / "rsem/{sample}.genes.results"), sample=samples_to_process) if run_rsem else [],
        # rsem_isoforms = expand(str(Path(data_path) / "rsem/{sample}.isoforms.results"), sample=samples_to_process) if run_rsem else [],
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
        # inputdir2 = str(Path(data_path) / "rsem") if run_rsem else "",
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
"""

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
