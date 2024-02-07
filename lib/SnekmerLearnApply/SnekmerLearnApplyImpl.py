# -*- coding: utf-8 -*-
#BEGIN_HEADER
# The header block is where all import statments should live
import logging
import os
from pprint import pformat
import subprocess 
from Bio import SeqIO

from installed_clients.AssemblyUtilClient import AssemblyUtil
from installed_clients.KBaseReportClient import KBaseReport
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
        logging.basicConfig(format='%(created)s %(levelname)s: %(message)s',
                            level=logging.INFO)
        #END_CONSTRUCTOR
        pass



    def run_SnekmerLearnApply(self, ctx, params):
        # ctx is the context object
        # return variables are: returnVal
        # BEGIN run_SnekmerLearnApply

        # Initialize logging
        logging.basicConfig(level=logging.INFO)
        log_messages = []

        # Step 1 - Validate and extract the parameters
        log_messages.append('Starting run_SnekmerLearnApply function. Params=' + pformat(params))
        if 'input_seqs' not in params:
            error_message = 'Parameter input_seqs is not set in input arguments'
            log_messages.append(error_message)
            logging.error(error_message)
            raise ValueError(error_message)
        input_seqs_ref = params['input_seqs']

        # Step 2 - Retrieve the ProteinSequenceSet object
        log_messages.append(f'Retrieving ProteinSequenceSet with ref: {input_seqs_ref}')
        try:
            # Placeholder for actual code to retrieve the ProteinSequenceSet object using input_seqs_ref
            # This should be replaced with actual data retrieval logic
            protein_seq_set = {
                'id': 'example_id',
                'md5': 'example_md5',
                'description': 'example_description',
                'sequences': [
                    {'id': 'seq1', 'sequence': 'MEEPQSDPSV', 'md5': 'example_md5_1'},
                    {'id': 'seq2', 'sequence': 'EEPQSDPSVE', 'md5': 'example_md5_2'}
                ]
            }
            log_messages.append(f'ProteinSequenceSet details: {pformat(protein_seq_set)}')
        except Exception as e:
            error_message = 'Error retrieving ProteinSequenceSet: ' + str(e)
            log_messages.append(error_message)
            logging.error(error_message)
            raise

        # Log steps (for illustration purposes, replace with actual steps in your function)
        # ...

        # Step X - Build a Report and return (assuming KBaseReport is properly initialized elsewhere)
        reportObj = {
            'text_message': '\n'.join(log_messages)
        }
        # Assuming KBaseReport is properly initialized and available
        report = KBaseReport(self.callback_url)
        report_info = report.create({'report': reportObj, 'workspace_name': params['workspace_name']})

        # Construct the output to send back
        output = {
            'report_name': report_info['name'],
            'report_ref': report_info['ref']
        }

        # END run_SnekmerLearnApply

        if not isinstance(output, dict):
            raise ValueError('Method run_SnekmerLearnApply return value output is not type dict as required.')
        return [output]
            
        # At some point might do deeper type checking...
    
        #For me. I think this is where i make my script
        #Steps:
        # 1 Gather inputs
        # 2 Build cofig file
        # 3 Pull this  - annotation option from dropdown. 
        # 4 Run snekmer learn / apply.
        # For learn - this might be stored as KBaseExperiments.CorrelationMatrix
        # For apply - the output might be KBaseSequences.ProteinSequenceSet
        # Annotations should be built into the docker image. Remove all but the one on dropdown

        # At some point might do deeper type checking...
    
    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK",
                     'message': "",
                     'version': self.VERSION,
                     'git_url': self.GIT_URL,
                     'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
