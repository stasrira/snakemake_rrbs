import os
from pathlib import Path
import re
import pandas as pd
import yagmail


# this function is used to retrieve configuration values
def get_config_value(yaml_path, conf_file = None, delim='/'):
    try:
        config_local = config
    except NameError:
        config_local = None
    if conf_file is None:
        conf_file = config_local  # assign reference to the standard cofiguration object
    
    path_elems = yaml_path.split(delim)

    # loop through the path to get the required key
    val = conf_file
    for el in path_elems:
        # make sure "val" is not None and continue checking if "el" is part of "val"
        if val and el in val:
            try:
                val = val[el]
            except Exception:
                val = None
                break
        else:
            val = None

    return val

# returns the name of conda env yaml file or a path to the conda environment (based on the corresponding entry in the config)
def get_conda_env(env_alias = None):
    if env_alias is None:
        env_alias = 'main_tools'
    return get_config_value ('conda_envs/{}'.format(env_alias))
    
# this function parses the name of the ref_data_genom_dir and concatenates it with the provided postfix
# the outcome will correspond to the folder name where index data should be located
def get_ref_index_dir (ref_data_genom_dir, postfix = None):
    # to get speceis_id, it parses ref_data_genom_dir, replaces any numbers, splits by "_", gets the first item from the list
    species_id = re.sub('[\\d]','', ref_data_genom_dir).split('_')[0]
    if postfix:
        return species_id  + "_" + postfix
    else:
        return species_id

def get_raw_file_path(data_path, sub_dir_path = None, F_id = None, F_type = None, file_name_format = None):
    """
    F_id: expected values 1 or 2, which stands for R1 or R2 or I1
    F_type: expected values: R or I, which stands for R1 or I1
    """
    if F_type is None: 
        F_type = 'R'
    # define default_file_name_format based on the existense of the global variable raw_data_file_ext
    try: 
        raw_data_file_ext
    except NameError: 
        default_file_name_format = '{sample}_|FN|.fastq.gz'
    else:
        default_file_name_format = '{sample}_|FN|' + raw_data_file_ext
        
    # default_file_name_format = '{sample}_|FN|.fastq.gz'    
    if sub_dir_path is None:
        sub_dir_path = 'fastq_raw'
    if F_id is None: 
        F_id = '1'
    else:
        F_id = str(F_id)
    if file_name_format is None:
        file_name_format = default_file_name_format
    file_name = file_name_format.replace('|FN|', F_type + F_id)
    
    if len(sub_dir_path) > 0:
        out = str(Path(data_path) / (sub_dir_path + '/' + file_name))
    else:
        out = str(Path(data_path) / file_name)
    return out

def get_R1_file_path(data_path, sub_dir_path = None, file_name_format = None):
    return get_raw_file_path(data_path, sub_dir_path, '1', 'R', file_name_format)

def get_R2_file_path(data_path, sub_dir_path = None, file_name_format = None):
    return get_raw_file_path(data_path, sub_dir_path, '2', 'R', file_name_format)

def get_I1_file_path(data_path, sub_dir_path = None, file_name_format = None):
    return get_raw_file_path(data_path, sub_dir_path, '1', 'I', file_name_format)

def replace_cfg_with_environment(var_name, cur_value):
    env_val = os.environ.get(var_name)
    if env_val:
        return env_val
    else:
        return cur_value

def get_and_validate_env_var(var_name):
    try:
        env_val = os.environ[var_name]
        if env_val:
            return env_val
        else:
            print (f'Error: Required environment variable "{var_name}" was not set with any value!')
            return ''
    except KeyError:
        print (f'Error: Required environment variable "{var_name}" was not found!')
        return ''
    except Exception as ex:
        print (ex)
        return None

def validate_rmats_contrast_file (cf_path, cf_path_expected, samples_lst,
                                  contrast_default_val = None, rmats_min_b1_qty = None, rmats_min_b2_qty = None):
    if contrast_default_val is None:
        contrast_default_val = 'b1'
    rmats_min_b1_qty = convert_to_int(rmats_min_b1_qty)
    rmats_min_b2_qty = convert_to_int(rmats_min_b2_qty)
    
    # check rmats_min_b1_qty and rmats_min_b2_qty for None and assign a default value
    # if not rmats_min_b1_qty or not rmats_min_b1_qty >= 0:
    if rmats_min_b1_qty is None:
        rmats_min_b1_qty = 0
    # if not rmats_min_b2_qty or not rmats_min_b2_qty >= 0:
    if rmats_min_b2_qty is None:
        rmats_min_b2_qty = 0
    
    expected_contrast_values = ['b1','b2']  # list of expected contrast values in the contrast file
    out_val = None
    meta_file_created = False
    
    if os.path.isfile(cf_path):
        # contrast file exists, validate it
        ext = Path(cf_path).suffix
        if ext == '.csv':
            # read the comma delimited file
            cf = pd.read_csv (cf_path, sep=",", header=0)  # read the file with headers
        else:
            # read the tab delimited file
            cf = pd.read_csv (cf_path, sep="\t", header=0)  # read the file with headers
        if len(cf) >= (rmats_min_b1_qty + rmats_min_b2_qty):
            if len(cf.columns) >= 2:
                cf_samples = cf.iloc[:, 0].astype("string").tolist()  # get first column to the list
                cf_contrast_vals = cf.iloc[:, 1].astype("string").unique().tolist()  # get second column to the list (unique values only)
                # filter out contrast file samples that are not in samples_lst before checking counts for b1 and b2
                cf_qualified = cf[cf.iloc[:, 0].astype("string").isin(samples_lst)]  # leaves only rows where samlpe from the first row is in samples_lst
                if (cf_qualified.iloc[:, 1] == 'b1').sum() < rmats_min_b1_qty:
                    add_pre_run_warning(
                        'ERROR: The rmats contrast file does not contain the expected minimum ({}) of "b1" entries (only entries matching list of samlpes to be processed were counted)'
                        .format(rmats_min_b1_qty))
                    out_val = False

                if (cf_qualified.iloc[:, 1] == 'b2').sum() < rmats_min_b2_qty:
                    add_pre_run_warning(
                        'ERROR: The rmats contrast file does not contain the expected minimum ({}) of "b2" entries (only entries matching list of samlpes to be processed were counted)'
                        .format(rmats_min_b2_qty))
                    out_val = False

                no_match_samples = [sample for sample in samples_lst if not sample in cf_samples]  # get all samples from cf_samples that are not in the given samples_lst
                no_match_contrasts = [ctr for ctr in cf_contrast_vals if not ctr in expected_contrast_values]  # get all contrast values from cf_contrast_vals that are not in the expected list of values (expected_contrast_values)

                if no_match_samples:
                    add_pre_run_warning ('WARNING: the following samples have no corresponding entries in the given rmats contrast file: {}'.format(no_match_samples))
                if no_match_contrasts:
                    add_pre_run_warning ('WARNING: the following contrast values in the given rmats contrast file do not match the expected values for this file: {}. The list of expected contrast values: {}'.format(no_match_contrasts, expected_contrast_values))
                if out_val is None:
                    # if out_val was not updated yet, set it to True
                    out_val = True
                if cf_path != cf_path_expected:
                    # if the contras file is not in the expected location of the data folder, create it there
                    save_rmats_config_file(cf, cf_path_expected)
            else:
                add_pre_run_warning ('ERROR: the following contrast file is in wrong format - it has less than 2 columns : {}'.format(cf_path))
                out_val = False
        else:
            add_pre_run_warning(
                'ERROR: the following contrast file does not have a minimal number ({}) of records: {}'.format(
                    (rmats_min_b1_qty + rmats_min_b2_qty), cf_path))
            out_val = False
    else:
        # contrast file is not present, create one
        cf_content = {'sample_id':[], 'contrast':[]}
        b2_count = 0
        b_counts = {'b1':0, 'b2':0}
        for sample in samples_lst:
            cf_content['sample_id'].append(sample)
            # fill minimum b2 entries
            if b_counts['b2'] < rmats_min_b2_qty:
                cf_content['contrast'].append('b2')
                b_counts['b2'] += 1
                continue
            # fill minimum b1 entries
            if b_counts['b1'] < rmats_min_b1_qty:
                cf_content['contrast'].append('b1')
                b_counts['b1'] += 1
                continue
            # fill all other entries using the provided default value
            cf_content['contrast'].append(contrast_default_val)

        cf_df = pd.DataFrame(cf_content)
        save_rmats_config_file(cf_df, cf_path_expected)
        # validate if the created file has a minimal number of records
        if len(cf_df) >= (rmats_min_b1_qty + rmats_min_b2_qty):
            out_val = True
        else:
            add_pre_run_warning(
                'ERROR: the following contrast file does not have a minimal number ({}) of records: {}'.format(
                    (rmats_min_b1_qty + rmats_min_b2_qty), cf_path))
            out_val = False
        meta_file_created = True
    return out_val, meta_file_created

def save_rmats_config_file(cf_dataframe, cf_path_expected):
    # check if required directory exists, create one if needed
    dir_path = os.path.dirname(cf_path_expected)
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)
    cf_dataframe.to_csv(cf_path_expected, sep ='\t', index=False)  # save contrast file 

def get_samples_to_ignore(samplesToIgnore_str):
    # print ('samplesToIgnore_str = '.format(samplesToIgnore_str))
    
    if samplesToIgnore_str is None:
        samplesToIgnore = []
    elif os.path.isfile(samplesToIgnore_str):
        # if a valid file path is provided, get files content
        # the file is expected to have 1 column with the list of sample ids and no headers
        igf = pd.read_csv (samplesToIgnore_str, sep=",", header=None)  # read the file without headers
        if len(igf.columns) >= 1:
            samplesToIgnore = igf.iloc[:, 0].astype("string").tolist()  # get first column to the list
        else:
            samplesToIgnore = []
        
    elif len(samplesToIgnore_str.strip()) > 0:
        # if string is provided, it is assumed to be comma delimited string
        samplesToIgnore = [smp.strip() for smp in samplesToIgnore_str.split(',')]
    else:
        samplesToIgnore = []
        
    return samplesToIgnore

def send_yagmail(emails_to, subject, message, email_from=None, attachments=None, smtp_server=None, smtp_server_port=None, attachment_max_bytes=None):
    emails_to_list = None
    attachment_lst = []
    if not email_from:
        email_from = os.environ.get('SNAKEMAKE_FROM_EMAIL')  # from the local config/.env file
    if not smtp_server:
        smtp_server = os.environ.get('SNAKEMAKE_SMTP_SERVER')  # from the local config/.env file
    if not smtp_server_port:
        smtp_server_port = os.environ.get('SNAKEMAKE_SMTP_SERVER_PORT')  # from the local config/.env file
    
    # if a comma delimeted list of emails was provided, convert it to the list
    if isinstance(emails_to, str):
        emails_to_list = emails_to.split(',')
    elif isinstance(emails_to, list):
        emails_to_list = emails_to
    
    if isinstance(attachments, str):
        # a single path was sent as an attachemnt
        if os.path.isfile(attachments):
            attachments = [attachemnts]
    elif isinstance(attachments, list):
        attachments = [atch for atch in attachments if os.path.isfile(atch)]
            
    # verify that the attachment file exists, otherwise set to None
    # if attachment_path:
        # if not os.path.isfile(attachment_path):
            # attachment_path = None
    
    if attachments:
        # the first attachment is the snakemake log file
        log_file_size = round(os.stat(attachments[0]).st_size/1024/1024, 2)  # get size in MB
        message = message.replace('|log_size|', '{} MB'.format(log_file_size))
        if not attachment_max_bytes is None:
            # validate that the attachemnt file is not over max size, otherwise remove the attachemnt from the list
            for atch in attachments:
                statinfo = os.stat(atch)
                if statinfo.st_size <= attachment_max_bytes:
                    attachment_lst.append(atch)
        else:
            attachment_lst = attachments
    else:
        message = message.replace('|log_size|', 'N/A')
    
    body = message
    
    if emails_to_list is None or email_from is None or smtp_server is None or smtp_server_port is None:
        add_warning('WARNING: Cannot send an email since one or more of the key values are not defined '
            '(current settings: smtp_server: {}, smtp_server_port: {}, email_from: {}, emails_to: {})'
            .format(smtp_server, smtp_server_port, email_from, emails_to_list)
            )
    else:
        yag = yagmail.SMTP(email_from,
                           host=smtp_server,
                           smtp_skip_login=True,
                           smtp_ssl=False,
                           soft_email_validation=False,
                           port=smtp_server_port)
        if attachment_lst:
            yag.send(
                to=emails_to_list,
                subject=subject,
                contents=body,
                attachments=attachment_lst,
            )
        else:
            yag.send(
                to=emails_to_list,
                subject=subject,
                contents=body,
            )

def add_pre_run_info(message, lst = None):
    if lst is None:
        try:
            lst = pre_run_messages['info']
        except NameError:
            lst = None
    add_pre_run_message (message, lst)

def add_pre_run_warning(message, lst = None):
    if lst is None:
        try:
            lst = pre_run_messages['warning']
        except NameError:
            lst = None
    add_pre_run_message (message, lst)

def add_pre_run_message(message, lst):
    if isinstance(lst, list):
        lst.append(message)
    else: 
        print ('Cannot save the following pre-run message since the object to save it in is not of list type: {}'.format(message))

def print_pre_run_info(lst = None):
    if lst is None:
        try:
            lst = pre_run_messages['info']
        except NameError:
            lst = None
    header = 'The following are the pre-run pipeline settings and sample info that will be used for the current run:'
    print_pre_run_messages (lst, header)
    
def print_pre_run_warnings(lst = None):
    if lst is None:
        try:
            lst = pre_run_messages['warning']
        except NameError:
            lst = None
    header = 'The following are the warnings generated during pre-run of the pipeline:'
    print_pre_run_messages (lst, header, sep = '*', add_blank_row_after_section = True)

def get_pre_run_warnings(lst = None, output_format = None):
    if lst is None:
        try:
            lst = pre_run_messages['warning']
        except NameError:
            lst = None
    if output_format is None:
        output_format = 'string'
    
    out = None
    
    if lst:
        if output_format == 'string':
            out = '\n'.join(lst)
        elif output_format == 'list':
            out = lst
    else:
        out = 'No messages to report'
    return out
    
def print_pre_run_messages(lst, header = None, sep = None, add_blank_row_after_section = None):
    sep_len = 0
    if add_blank_row_after_section is None:
        add_blank_row_after_section = False
    if isinstance(lst, list):
        if header: 
            sep_len = print_header (header, sep)
        
            
        if len(lst) > 0:
            for msg in lst:
                print(msg)
        else:
            print ('No messages to report')
        
        if header:
            print_info_separator(sep_len, sep, add_blank_row_after_section)
    else:
        print('WARNING: Provided list of messages is not the "list" type. Retrieving messages has failed.')

# adds a message to the pipeline info file
def add_info(message, save_to_path = None, mode = None):
    # print('add_info ========>{}'.format(message))
    file_type = 'info'
    add_message_to_file (message, save_to_path, mode, file_type)

# adds a message to the pipeline warning file
def add_warning(message, save_to_path = None, mode = None):
    file_type = 'warning'
    add_message_to_file (message, save_to_path, mode, file_type)

def add_message_to_file(message, save_to_path = None, mode = None, file_type = None):
    if file_type is None:
        file_type = 'info'
    if save_to_path is None:
        try:
            if file_type == 'info':
                save_to_path = pipeline_info_file_path
            elif file_type == 'warning':
                save_to_path = pipeline_warning_file_path
            else:
                save_to_path = None
                print('WARNING: Unexpected message type "{}" was provided! Expected types are: "info", "warning".'.format(file_type))
        except NameError:
            save_to_path = None
    # print('add_message_to_file --> save_to_path = {}'.format(save_to_path))
    if mode is None:
        mode = 'a'  # append to the file 
    if save_to_path:
        with open(save_to_path, mode) as fw:
            if isinstance(message, str):
                if len(message) > 0:
                    # add return character for the not blank messages
                    fw.write(message + '\n')
                else:
                    # blank message in combination with the mode = 'w' can be used to reset the file
                    fw.write(message)
            elif isinstance(message, list):
                msg_lst = [msg + '\n' for msg in message]
                fw.writelines(msg_lst)
    else:
        print('WARNING: No pipeline message file path was provided! The following message was not recorded: {}'.format(message))

def print_global_info(info_file_path = None):
    if info_file_path is None:
        try:
            info_file_path = pipeline_info_file_path
        except NameError:
            info_file_path = None
    if info_file_path:
        header = 'The following are the pipeline settings and sample info used for the current run:'
        print_global_messages (info_file_path, header)
    else:
        print('WARNING: No pipeline file_info path was provided! Retrieving info messages has failed.')
        
def print_global_warnings(warning_file_path = None):
    if warning_file_path is None:
        try:
            warning_file_path = pipeline_warning_file_path
        except NameError:
            warning_file_path = None
    if os.path.isfile(warning_file_path):
        header = 'The following are the warnings generated during run of the pipeline:'
        print_global_messages (warning_file_path, header, sep = '*', add_blank_row_after_section = True)
    else:
        print('No warnings are available since the pipeline file_warning file does not exist (which is normal).')
 
def get_global_warnings(msg_file_path = None):
    if msg_file_path is None:
        try:
            msg_file_path = pipeline_warning_file_path
        except NameError:
            msg_file_path = None
    
    out = None
    
    if os.path.isfile(msg_file_path):
        with open(msg_file_path, 'r') as f:
            file_content = f.read()
            if len(file_content.strip()) > 0:
                out = file_content.strip()
            else:
                out = 'No messages to report'
    else:
        out = 'WARNING: No pipeline message file path was provided or the file does not exist!Retrieving messages has failed.'
    return out
 
def print_global_messages(msg_file_path, header = None, sep = None, add_blank_row_after_section = None):
    sep_len = 0
    if add_blank_row_after_section is None:
        add_blank_row_after_section = False
    if os.path.isfile(msg_file_path):
        with open(msg_file_path, 'r') as f:
            if header: 
                sep_len = print_header (header, sep)
            file_content = f.read()
            if len(file_content.strip()) > 0:
                print(file_content[:-1])  # print string retrieved from the file without the last character (which is a return character)
            else:
                print ('No messages to report')
            if header:
                print_info_separator(sep_len, sep, add_blank_row_after_section)
    else:
        print('WARNING: No pipeline message file path was provided or the file does not exist! Retrieving messages has failed.')

def print_header(message, sep = None):
    m_len = len(message)
    print  ('\n')
    print_info_separator(m_len, sep)
    print  (message)
    print_info_separator(m_len, sep)
    return m_len

def print_info_separator(m_len, sep = None, add_blank_row_after = None):
    if sep is None:
        sep = '-'
    print  ('{}'.format(sep * m_len))
    if add_blank_row_after:
        print  ('\n')
    
def temp_debugcheck(str):
    # checks if global variable "debug' is set to True, do not mark provided arguments as snakemake's "temp" object
    # otherwise return the provided argument "str" as snakemake's temp
    if debug:
        return str
    else:
        return temp(str)

def str_to_bool(value):
    """Convert string value to boolean."""
    valid = {'true': True, 't': True, '1': True,
         'false': False, 'f': False, '0': False,
         }

    if isinstance(value, bool):
        return value

    if not isinstance(value, str):
        return None

    lower_value = value.lower()
    if lower_value in valid:
        return valid[lower_value]
    else:
        raise None

def convert_to_int(s):
    try:
        int_value = int(s)
        return int_value
    except ValueError:
        return None

# this prepares suffix to be added to the cluster's job name (where {rule_name} is always included as the first part of the job name)
def get_job_suffix(wildcards, custom_val = None):
    if not custom_val is None:
        return custom_val
    elif hasattr(wildcards, 'sample'):
        # if sample attribute exists
        return wildcards.sample
    elif hasattr(wildcards, 'event'):
        # if event attribute exists, but no sample proccessed
        return {wildcards.event}
    elif hasattr(wildcards, 'chunk'):
        # if chunk attribute exists
        return wildcards.chunk
    else:
        # provide no suffix, so just {rule_name} will be used for naming the cluster job 
        return ''

# function to get a dynamic attempt number to be passed as an additional resource parameter 
# It is used for simple visual representation of the current attempt number in the log file
def get_attempt(wildcards, attempt):
    return attempt

# this is a wrapper for prepare_cluster_resources function 
# it is used if the allocated memory has to be identified dynamically based on the attempt number
def get_cluster_resources_by_attempt (attempt_num, rule_name, config_file = None):
    memory = get_memory_by_attempt_and_threads(
                                attempt_num = attempt_num,
                                memory = get_rule_memory(rule_name = rule_name, conf_file = config_file),
                                threads = get_rule_threads(rule_name = rule_name, conf_file = config_file),
                                config_file = config_file
                                )
    return prepare_cluster_resources(
                                memory = memory,
                                single_host = get_rule_single_host(rule_name = rule_name, conf_file = config_file),
                                himem = get_rule_himem(rule_name = rule_name, conf_file = config_file)
                                )

# common function to prepare resource string for submitting to cluster
def prepare_cluster_resources(memory, single_host = None, himem = None):
    if single_host is None:
        single_host = False
    if himem is None:
        himem = False
        
    single_host_val = ' span[hosts=1]' if single_host else ''
    himem_val = ' select[himem]' if himem else ''
    
    resources = f'"rusage[mem={memory}]{single_host_val}{himem_val}"'
    return resources

# retrievs memory assigned for a rule        
def get_rule_memory(rule_name, conf_file = None, default_rule_memory = None):
    default_mem = 1024
    if default_rule_memory is None:
        # attempt to retrive default_rule_memory from config, but use default_mem, if None is received
        default_rule_memory = get_config_value ('default_rule_memory', conf_file)
        default_rule_memory = default_rule_memory if default_rule_memory is not None else default_mem
    
    # get rule's memory from config, but use default_rule_memory if None is received
    rule_memory = get_config_value ('rules/{}/memory'.format(rule_name), conf_file)
    return rule_memory if rule_memory is not None else default_rule_memory

# retrievs single_host directive assigned for a rule    
def get_rule_single_host(rule_name, conf_file = None, default_single_host = None):
    if default_single_host is None:
        default_single_host = False
    # get rule's single_host value from config, but use default_single_host if None is received
    rule_single_host = get_config_value ('rules/{}/single_host'.format(rule_name), conf_file)
    return rule_single_host if rule_single_host is not None else default_single_host

# retrievs walltime directive assigned for a rule    
def get_rule_walltime(rule_name, conf_file = None, default_walltime = None):
    if default_walltime is None:
        default_walltime = "12:00"
    # get rule's walltime value from config, but use default_walltime if None is received
    rule_walltime = get_config_value ('rules/{}/walltime'.format(rule_name), conf_file)
    return rule_walltime if rule_walltime is not None else default_walltime

# retrievs himem directive assigned for a rule    
def get_rule_himem(rule_name, conf_file = None, default_himem = None):
    if default_himem is None:
        default_himem = False
    # get rule's single_host value from config, but use default_single_host if None is received
    rule_himem = get_config_value ('rules/{}/himem'.format(rule_name), conf_file)
    return rule_himem if rule_himem is not None else default_himem

# retrievs threads assigned for a rule    
def get_rule_threads (rule_name, conf_file = None, default_threads = None):
    if default_threads is None:
        default_threads = 1
    # get rule's threads from config, but use default_threads if None is received
    rule_threads = get_config_value ('rules/{}/threads'.format(rule_name), conf_file)
    return rule_threads if rule_threads is not None else default_threads
    
# retrievs retries assigned for a rule    
def get_rule_retries (rule_name, conf_file = None, default_retries = None):
    if default_retries is None:
        default_retries = 1
    # get rule's threads from config, but use default_threads if None is received
    rule_retries = get_config_value ('rules/{}/retries'.format(rule_name), conf_file)
    return rule_retries if rule_retries is not None else default_retries

# def get_memory_by_attempt(attempt_num, memory, max_memory = None):
    # if max_memory is None:
        # max_memory = 300000  # 300GB is the default maximum memory allocation
        
    # for i in range(attempt_num):
        # if i > 0:
            # cur_mem = cur_mem * 2  # double memory allocation for each next attempt
        # else:
            # cur_mem = memory  # for first time use the assigned memory allocation
    
    # if cur_mem > max_memory:
        # cur_mem = max_memory  # if cur_mem is over max_memory allocation, reset it to the set maximum
        
    # return cur_mem
    
def get_memory_by_attempt_and_threads(attempt_num, memory, threads = None, max_memory = None, config_file = None):
    default_max_memory = 1400000  # 1.4 TB
    if max_memory is None:
        max_memory = get_config_value ('resources/max_memory', config_file)
        max_memory = max_memory if max_memory is not None else default_max_memory
    
    if threads is None:
        threads = 1
    
    for i in range(attempt_num):
        if i > 0:
            cur_mem = cur_mem * 2  # double memory allocation for each next attempt
        else:
            cur_mem = memory  # for first time use the assigned memory allocation
    
    if cur_mem * threads > max_memory:
        cur_mem = max_memory / threads # if cur_mem is over max_memory allocation, use the maximum divided by number of threads
        
    return cur_mem

def get_memory_by_samples_count(num_samples, memory, samples_high_count = None, high_memory_multiplier = None):
    if samples_high_count is None:
        samples_high_count = 500
    if high_memory_multiplier is None:
        high_memory_multiplier = 2
        
    if num_samples < samples_high_count:
        return memory
    else:
        return memory * high_memory_multiplier
        
def get_retries_by_samples_count(num_samples, retries, samples_high_count = None, high_memory_divider = None):
    import math
    if samples_high_count is None:
        samples_high_count = 500
    if high_memory_divider is None:
        high_memory_divider = 2
        
    if num_samples < samples_high_count:
        return retries
    else:
        return int(math.floor(retries / high_memory_divider))  # round to an up integer number

def get_queue_by_samples_count (num_samples, time_coef, \
                                min_walltime = None, max_regular_walltime = None, \
                                max_long_walltime = None, \
                                normal_queue = None, long_queue = None):
    if normal_queue is None:
        normal_queue = 'premium'
    if long_queue is None:
        long_queue = 'long'
    
    wall_time = get_walltime_by_samples_count(num_samples, time_coef, \
                                            min_walltime, max_regular_walltime, \
                                            max_long_walltime, True)
    if wall_time <= max_regular_walltime:
        return normal_queue
    else:
        return long_queue
        
def get_walltime_by_samples_count(num_samples, time_coef, 
                                min_walltime = None, max_regular_walltime = None, \
                                max_long_walltime = None, int_value = None):
    import math
    if min_walltime is None:
        min_walltime = 12
    if max_regular_walltime is None:
        max_regular_walltime = 144
    if max_long_walltime is None:
        max_long_walltime = 336
    if int_value is None:
        int_value = False
    
    num_samples = convert_to_int(num_samples)
    if num_samples is None:
        num_samples = 0
    out_wall_time = 12
    
    walltime = int(math.floor(num_samples * time_coef))
    if walltime < min_walltime:
        out_wall_time =  min_walltime
    elif walltime > max_long_walltime: 
        out_wall_time = max_long_walltime
    else:
        out_wall_time = walltime
    
    if int_value:
        return out_wall_time
    else:
        return "{}:00".format(out_wall_time)  # convert to expected format with hours/minutes

def print_help_info():
    help_msg = []
    header = '---------------------- Pipeline Help info ----------------------'
    help_msg.append ('This pipeline help info is displayed, because "help" environment variable was set to True. Set it to False or unset the variable to avoid this behaviour.')
    help_msg.append ('Listed below are the configuration variables (with the default values provided) that can be overwritten by providing environmental variables with the same names.')
    help_msg.append ('Among these variables the following are configuration variables that require additional information:')
    help_msg.append ('    - help: displays this message and aborts pipeline execution.')
    help_msg.append ('    - debug: if set to True will instruct snakemake not to delete any files/folder marked as "temp".')
    help_msg.append ('    - data_path: (required) It must contain a path to the data directory that is being processed.')
    help_msg.append ('    - ref_data_path: (optional) It has to contain a path to the main location of the reference data (see defult value from the config below).')
    help_msg.append ('')
    help_msg.append ('Here is the full list of overwritable variables:')
    print_all_overwritable_vars (config, help_msg)
    print_pre_run_messages (help_msg, header, sep = '*', add_blank_row_after_section = True)

def print_all_overwritable_vars(cnf, msg_list, prefix = '    - '):
    if cnf:
        for item in cnf:
            # print ('isinstance(cnf[item], dict) = {}'.format(isinstance(cnf[item], dict)))
            # print ('isinstance(cnf[item], list) = {}'.format(isinstance(cnf[item], list)))
            # print ('isinstance(cnf[item], str) = {}'.format(isinstance(cnf[item], str)))
            # print ('isinstance(cnf[item], bool) = {}'.format(isinstance(cnf[item], bool)))
            if not isinstance(cnf[item], dict) and not isinstance(cnf[item], list):
                if isinstance(msg_list, list):
                    msg_list.append('{}{} = {}'.format(prefix, item, cnf[item]))
                else:
                    print('{}{} = {}'.format(prefix, item, cnf[item]))
 
