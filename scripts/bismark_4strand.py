import sys
import re
import traceback

# Mapping from the 4 Bismark triplet codes to the desired columns
CODE_TO_COL = {
    "CT/GA/CT": "OT",    # ((converted) top strand)
    "CT/GA/GA": "OB",    # ((converted) bottom strand)
    "GA/CT/CT": "CTOT",  # (complementary to (converted) top strand)
    "GA/CT/GA": "CTOB",  # (complementary to (converted) bottom strand)
}

# Regex to capture lines like:
#   CT/GA/CT:    4534    ((converted) top strand)
LINE_RE = re.compile(r'^(CT/GA/CT|CT/GA/GA|GA/CT/CT|GA/CT/GA):\s+(\d+)\b')

def parse_bismark_report(report_path, rule_name):
    """
    Parse a Bismark *_report.txt file and return a dict with keys OT, OB, CTOT, CTOB (ints).
    We scan the whole file and pick the first occurrence of each of the 4 codes.
    """
    counts = {"OT": None, "OB": None, "CTOT": None, "CTOB": None}
    try:
        with open(report_path, 'r', encoding='utf-8', errors='replace') as fh:
            for line in fh:
                m = LINE_RE.match(line)
                if not m:
                    continue
                code, num = m.group(1), int(m.group(2))
                col = CODE_TO_COL[code]
                if counts[col] is None:
                    counts[col] = num
    except FileNotFoundError:
        msg = f'Error (({rule_name} rule)): Expected to exists report not found: {report_path}'
        print(msg)
        add_warning (msg, pipeline_warning_file_path)
        raise
    
    # check 'counts' dict for not populated values
    missing = [k for k, v in counts.items() if v is None]
    if missing:
        # If any strands were missing, report an warning
        msg = f'Warning (({rule_name} rule)): Missing strand counts {missing} in report: {report_path}'
        print(msg)
        add_warning (msg, pipeline_warning_file_path)

    return counts

def create_strand_report(samples, align_reports, out_file, rule_name, pipeline_warning_file_path):
    # Parse and write output
    with open(out_file, 'w', encoding='utf-8') as out:
        out.write("samples\tOT\tOB\tCTOT\tCTOB\n")
        for sample, rpt in zip(samples, align_reports):
            counts = parse_bismark_report(rpt, rule_name)
            out.write(f"{sample}\t{counts['OT']}\t{counts['OB']}\t{counts['CTOT']}\t{counts['CTOB']}\n")

def report_error(ex):
    msg = 'Unexpected Error "{}" occurred while executing bismark_4strand.py scirpt. \n{}'.format(ex, traceback.format_exc())
    print (msg)
    if pipeline_warning_file_path:
        add_warning (msg, pipeline_warning_file_path)

try: 
    pipeline_warning_file_path = None
    
    log_file = snakemake.log[0]
    with open(log_file, "w") as f:
        sys.stderr = sys.stdout = f  # allows creating log file for the step
        
        samples = snakemake.params.samples
        align_reports = snakemake.input.align_reports
        out_file = snakemake.output.out_file
        rule_name = snakemake.params.rule_name
        pipeline_warning_file_path = snakemake.params.pipeline_warning_file_path
        
        print (f'samples = {samples}')
        print (f'align_reports = {align_reports}')
        print (f'out_file = {out_file}')
        print (f'rule_name = {rule_name}')
        print (f'pipeline_warning_file_path = {pipeline_warning_file_path}')
        
        create_strand_report ( \
            samples = samples, \
            align_reports = align_reports, \
            out_file = out_file, \
            rule_name = rule_name, \
            pipeline_warning_file_path = pipeline_warning_file_path)
                
except Exception as ex:
    # report unexpected error 
    report_error(ex)
    raise
