# -*- coding: utf-8 -*-
#BEGIN_HEADER
# The header block is where all import statments should live
import logging
import os
from pprint import pformat
import subprocess 
from Bio import SeqIO
import sys
import zipfile
import uuid
from datetime import datetime
import pandas as pd

from installed_clients.AssemblyUtilClient import AssemblyUtil
from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.WorkspaceClient import Workspace as workspaceService
from installed_clients.cb_annotation_ontology_apiClient import CBAnnotationOntologyAPI

# from installed_clients.DataFileUtilClient import DataFileUtil
# from installed_clients.GenomeFileUtilClient import GenomeFileUtil
# from installed_clients.KBaseDataObjectToFileUtilsClient import KBaseDataObjectToFileUtils
# from installed_clients.GenomeAnnotationAPIClient import GenomeAnnotationAPI
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


        object_refs = [{'ref': ref} for ref in input_seqs]

        # Fetch the objects using get_objects2
        protein_seq_set = self.wsClient.get_objects2({'objects': object_refs})

        # report = KBaseReport(self.callback_url)
        text_message = '\n'.join(params['input_seqs'])
        # report_info = report.create({'report': {'objects_created': [],
        #                                 'text_message': text_message},
        #                                 'workspace_name': params['workspace_name']})
        
        logging.info(object_refs)
        logging.info(sys.version)
        # logging.info(protein_seq_set)
    
    
        output_dir = "input"
        output_file_name = "output1.fasta"
        output_file_path = os.path.join(output_dir, output_file_name)

        # Ensure the output directory exists
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
    
        with open(output_file_path, 'w') as fasta_file:
            # Iterate over each fetched object in the response
            for obj in protein_seq_set['data']:
                # Extract sequences data for each object
                sequences_data = obj['data']['sequences']
                # Iterate over each sequence
                for seq in sequences_data:
                    # Construct the header with sequence ID and description
                    header = ">{} {}".format(seq['id'], seq['description'])
                    # Write the header and sequence to the FASTA file
                    fasta_file.write("{}\n{}\n".format(header, seq['sequence']))

        # Read the FASTA file and log the first 40 sequences
        sequences = list(SeqIO.parse(output_file_path, 'fasta'))
        for seq_record in sequences[0:40]:
            logging.info("ID: {}".format(seq_record.id))
            logging.info("Description: {}".format(seq_record.description))
            logging.info("Sequence: {}\n".format(str(seq_record.seq)[:60]))
            
        
        cwd = os.getcwd()
        cwd_contents = os.listdir(cwd)
        
        logging.info(os.getcwd())
        logging.info(cwd_contents)
        
        
        cmd_string = "snekmer apply --cores=4"
        cmd_process = subprocess.Popen(cmd_string, stdout=subprocess.PIPE,
                                    stderr=subprocess.STDOUT, cwd=os.getcwd(),
                                    shell=True)
        
        output, errors = cmd_process.communicate()
        logging.info('return code: ' + str(cmd_process.returncode))
        logging.info("="*80)
        logging.info("output: " + str(output) + '\n')
        logging.info("errors: " + str(errors) + '\n')
                
        cwd_contents = os.listdir(os.getcwd())
        logging.info(str(cwd_contents))
        
        apply_dir_path = os.path.join(cwd, "output", "apply")
        logging.info(os.listdir(apply_dir_path))
    
        # New Stuff start
        specific_file_path = os.path.join(cwd, "output", "apply", "kmer-summary-output1.csv")
        logging.info(specific_file_path)

        # Setup output directory for the ZIP file
        output_directory = os.path.join(self.shared_folder, str(uuid.uuid4()))
        os.makedirs(output_directory)
        run_date = datetime.now().strftime("%Y.%m.%d-%I:%M:%S%p")
        result_name = "KmerSummaryOutput" + run_date + ".zip"
        result_file = os.path.join(output_directory, result_name)


        ### new stuff
        # Read the results into a DataFrame
        annotations_df = pd.read_csv(result_file)

        # Prepare the ontology events based on your result_file structure
        ontology_events = []
        for index, row in annotations_df.iterrows():
            ontology_event = {
                "event_id": "your_event_id",  # This should be a unique ID for the event
                "description": "Protein annotation event",
                "ontology_id": row['Prediction'],  # Assuming this is the ontology ID
                "method": "SnekmerLearnApply_annotation",
                "method_version": "1.0",
                "timestamp": datetime.now().isoformat(),
                "feature_types": {},  # Fill in if you have feature types
                "ontology_terms": {
                    row['index']: [  # The protein ID
                        {
                            "term": row['Prediction'],
                            "modelseed_ids": [],  # Fill in if you have ModelSEED IDs
                            "evidence_only": 0,
                            "evidence": {
                                "reference": ("protein_annotation", "SnekmerLearnApply"),
                                "scores": {"confidence": row['Confidence']}  # Or any other scores
                            }
                        }
                    ]
                }
            }
            ontology_events.append(ontology_event)

        # Use the cb_annotation_ontology_api to add the ontology events
        ontology_api = CBAnnotationOntologyAPI(url=self.callback_url, token=os.environ.get('KB_AUTH_TOKEN'))
        params = {
            "input_ref": object_refs[0],
            "events": ontology_events,
            "output_name": "updated_protein_set_with_annotations",
            "output_workspace": params['workspace_name'],
        }
        # Call the method to add ontology events
        update_result = ontology_api.add_annotation_ontology_events(params)
        updated_protein_set_ref = update_result.get('output_ref')
        # New Stuff end


        # Zip the specific output file
        with zipfile.ZipFile(result_file, 'w', zipfile.ZIP_DEFLATED, allowZip64=True) as zip_file:
            # Add the specific file to the zip archive, adjust the arcname to change its name within the archive if needed
            zip_file.write(specific_file_path, os.path.basename(specific_file_path))




        # Prepare report parameters with the zipped file
        report_params = {
            'message': text_message,
            'workspace_name': workspace_name,
            'objects_created': [{
                'ref': updated_protein_set_ref,
                'description': 'Updated protein set with new ontologies'}],
            'file_links': [
                {
                    'path': result_file,
                    'name': os.path.basename(result_file),  
                    'label': 'Kmer Summary Output Archive',
                    'description': 'Zipped archive of Kmer analysis summary output'
                }
            ]
        }

        # Create the report
        report_client = KBaseReport(self.callback_url)
        report_info = report_client.create_extended_report(report_params)

        # Prepare output
        output = {
            'report_name': report_info['name'],
            'report_ref': report_info['ref'],
            'protein_set_ref': updated_protein_set_ref
        }

        # Log the details for debugging
        logging.info("Zipped file directory: " + output_directory)
        logging.info("=" * 80)
        logging.info("Zipped result file: " + result_file)


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
