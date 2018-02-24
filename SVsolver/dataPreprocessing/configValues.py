import json

class ConfigValues():

    def __init__(self, config):
        #sample
        self.bam_file = config['sample']['bam_file']
        self.coverage_directory = config['sample']['coverage_directory']
        self.coverage_file_prefix = config['sample']['coverage_file_prefix']

        #reference genom
        self.cyto_band = config['reference_genom']['cyto_band']

        #filtering
        self.contigs_names = config['filtering']['contigs_names'].split(',')
        self.contigs = config['filtering']['contigs'].split(',')
        self.min_aligment_length = int(config['filtering']['min_aligment_length'])
        self.min_read_aligment = float(config['filtering']['min_read_aligment'])

        #micro_sv
        self.max_alignment_distance = int(config['micro_sv']['max_alignment_distance'])
        self.max_micro_deletion_length = int(config['micro_sv']['max_micro_deletion_length'])
        self.min_micro_insertion_length = int(config['micro_sv']['min_micro_insertion_length'])
        self.min_micro_sv_support = int(config['micro_sv']['min_micro_sv_support'])
        self.max_sv_distance_to_merge = int(config['micro_sv']['max_sv_distance_to_merge'])

        #graph
        self.include_not_mapped = json.loads(config['graph']['include_not_mapped'])
        self.min_breakpoint_support = int(config['graph']['min_breakpoint_support'])
        self.max_breakpoint_distance_to_merge = int(config['graph']['max_breakpoint_distance_to_merge'])

        #optimisation
        self.mean_coverage = int(config['optimisation']['mean_coverage'])
        self.copy_numbers = json.loads(config['optimisation']['copy_numbers'])
        self.prior_breakpoint = json.loads(config['optimisation']['prior_breakpoint'])
        self.prior_segment = json.loads(config['optimisation']['prior_segment'])
        self.penalty_coefficient = float(config['optimisation']['penalty_coefficient'])
        self.penalty_coefficient_multiplier = float(config['optimisation']['penalty_coefficient_multiplier'])
        self.max_iteration = int(config['optimisation']['max_iteration'])
        self.min_x_change = float(config['optimisation']['min_x_change'])
        self.min_x_change_iterations = int(config['optimisation']['min_x_change_iterations'])

        #output
        self.output_path = config['output']['output_path']