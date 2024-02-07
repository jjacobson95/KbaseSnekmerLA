# -*- coding: utf-8 -*-
#BEGIN_HEADER
# The header block is where all import statments should live
import logging
import os
from pprint import pformat
import subprocess 
from Bio import SeqIO
import sys


from installed_clients.AssemblyUtilClient import AssemblyUtil
from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.WorkspaceClient import Workspace as workspaceService
#END_HEADER


class SnekmerLearnApply:
    '''
    Module Name:
    SnekmerLearnApply

    Module Description:
    A KBase module: SnekmerLearnApply
This sample module contains one small method that filters contigs.
This will have to be changed soon.
    '''

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "0.0.1"
    GIT_URL = "https://github.com/jjacobson95/KbaseSnekmerLA.git"
    GIT_COMMIT_HASH = "1421a376a1e6a3f23c863257f976721e20b6acce"

    #BEGIN_CLASS_HEADER
    # Class variables and functions can be defined in this block
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        
        # Any configuration parameters that are important should be parsed and
        # saved in the constructor.
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        self.shared_folder = config['scratch']
        self.workspaceURL = config['workspace-url']
        self.wsClient = workspaceService(self.workspaceURL)
        logging.basicConfig(format='%(created)s %(levelname)s: %(message)s',
                            level=logging.INFO)
        #END_CONSTRUCTOR
        pass


    def run_SnekmerLearnApply(self, ctx, params):
        """
        run_Snekmer_model accepts some of the model params for now, and returns results in a KBaseReport
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN run_SnekmerLearnApply


        logging.info('Starting run_Snekmer_search function. Params=' + pformat(params))

        logging.info('Validating parameters.')

        # check inputs
        workspace_name = params['workspace_name']
        if 'input_seqs' not in params:
            raise ValueError('Parameter kmer is not set in input arguments')
        input_seqs = params['input_seqs']


        text_message = '\n'.join(params['input_seqs'])
        
        protein_seq_set= self.wsClient.get_objects2({'objects': [{'ref': text_message}]})

        # report = KBaseReport(self.callback_url)
        # text_message = '\n'.join(params['input_seqs'])
        # report_info = report.create({'report': {'objects_created': [],
        #                                 'text_message': text_message},
        #                                 'workspace_name': params['workspace_name']})
        
        logging.info(text_message)
        logging.info(sys.version)
        # logging.info(protein_seq_set)
    
        output_file_path = "output1.fasta"

        with open(output_file_path, 'w') as fasta_file:
            # Assuming we want the first item in the list which contains 'sequences'
            sequences_data = protein_seq_set['data'][0]['data']['sequences']  # Adjusted access here
            # Iterate over each sequence
            for seq in sequences_data:
                # Construct the header with sequence ID and description using str.format()
                header = ">{} {}".format(seq['id'], seq['description'])
                # Write the header and sequence to the FASTA file
                fasta_file.write("{}\n{}\n".format(header, seq['sequence']))

        sequences = list(SeqIO.parse(output_file_path, 'fasta'))

        # Example: Print each sequence's ID and its first 60 nucleotides
        for seq_record in sequences[0:40]:
            logging.info("ID: {}".format(seq_record.id))
            logging.info("Description: {}".format(seq_record.description))
            logging.info("Sequence: {}\n".format(str(seq_record.seq)[:60]))
    
        report_params = {
            'message': text_message,
            'workspace_name': workspace_name,
            'objects_created': []}

        report_client = KBaseReport(self.callback_url)
        report_info = report_client.create_extended_report(report_params)
        
        
        logging.info(os.getcwd())
        logging.info(os.listdir())
        
        cmd_string = "snekmer apply"
        cmd_process = subprocess.Popen(cmd_string, stdout=subprocess.PIPE,
                                       stderr=subprocess.STDOUT, cwd=self.shared_folder,
                                       shell=True)
        cmd_process.wait()
        logging.info('return code: ' + str(cmd_process.returncode))
        logging.info("="*80)
        output, errors = cmd_process.communicate()
        logging.info("output: " + str(output) + '\n')
        logging.info("errors: " + str(errors) + '\n')
    
        output = {
            'report_name': report_info['name'],
            'report_ref': report_info['ref'],
        }

        #END run_SnekmerLearnApply

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method run_Snekmer_model return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]
            

    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK",
                     'message': "",
                     'version': self.VERSION,
                     'git_url': self.GIT_URL,
                     'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
