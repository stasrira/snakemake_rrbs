from Bio import SeqIO
import gzip
from common_functions import *
import os


# it goes through the given fastq file (I1 file is expected), identifies sequence for each entry and compares its length to the given default_len
def validate_fastq_file(fastq_path, output_valid_path, output_invalid_path, default_len = None, pipeline_warning_file_path = None): # status_file_path, 
    if default_len is None:
        default_len = 8
    
    umi_length_dict = {}  # Dictionary to store UMIs of different lengths
    
    # read fastq file
    # with open(fastq_path, 'r') as fastq_file:
    with gzip.open(fastq_path, 'rt') as fastq_file:
        for record in SeqIO.parse(fastq_file, 'fastq'):
            umi = str(record.seq)
            umi_length = len(umi)

            # Check if the UMI length is not equal to deault (8)
            if umi_length != default_len:
                # Collect UMIs of non-standard lengths in the dictionary
                if umi_length not in umi_length_dict:
                    umi_length_dict[umi_length] = []
                umi_length_dict[umi_length].append(umi)
    
    # delete all existing validate output files before creating new ones
    if os.path.isfile(output_valid_path):
        os.remove(output_valid_path)
    if os.path.isfile(output_invalid_path):
        os.remove(output_invalid_path)
    
    # Save results to files
    if not umi_length_dict:
        # If umi_length_dict is empty, save a valid file
        # msg_status = 'valid'
        with open(output_valid_path, 'w') as valid_file:
            msg = 'All sequences are of the expected length {}'.format(default_len)
            valid_file.write(msg)
            print(msg)
    else:
        # If umi_length_dict is not empty, save an invalid file
        # msg_status = 'invalid'
        with open(output_invalid_path, 'w') as invalid_file:
            msg = 'The following are invalid UMI sequence lengths (different from {}) with the associated counts of found entries:'.format(default_len)
            invalid_file.write(msg + '\n')
            print(msg)
            msg_prep = 'Outcome of validating of fastq file: {}'.format(fastq_path)
            add_warning (msg_prep, pipeline_warning_file_path) # fastq_path
            add_warning (msg, pipeline_warning_file_path)
            
            # for length, umi_list in umi_length_dict.items():
            for length in sorted(umi_length_dict):
                umi_list = umi_length_dict[length]
                msg = 'Length {}: {} entries'.format(length, len(umi_list))
                invalid_file.write(msg + '\n')
                print(msg)
                add_warning (msg, pipeline_warning_file_path)
                

# open log file and redirect there stdout and stderr outputs    
with open(snakemake.log[0], 'w') as f:
    sys.stderr = sys.stdout = f

    fastq_path = snakemake.input.fastq_file
    # status_file_path = snakemake.output.status_file
    output_valid_path = snakemake.output.valid_file
    output_invalid_path = snakemake.params.invalid_file
    default_len = snakemake.params.default_seq_len
    # pipeline_info_file_path = snakemake.params.pipeline_info_file_path
    pipeline_warning_file_path = snakemake.params.pipeline_warning_file_path

    # Replace 'your_fastq_file.fastq' with the actual Fastq file name
    validate_fastq_file(
        fastq_path = fastq_path, 
        output_valid_path = output_valid_path, 
        output_invalid_path = output_invalid_path, 
        # status_file_path = status_file_path,
        default_len = default_len,
        pipeline_warning_file_path = pipeline_warning_file_path
        )

