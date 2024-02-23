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
import csv

from installed_clients.AssemblyUtilClient import AssemblyUtil
from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.WorkspaceClient import Workspace as workspaceService
from installed_clients.cb_annotation_ontology_apiClient import cb_annotation_ontology_api

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
        # logging.info(protein_seq_set)
        logging.info(sys.version)
        logging.info(protein_seq_set)
    
    
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
    
        specific_file_path = os.path.join(cwd, "output", "apply", "kmer-summary-output1.csv")
        logging.info(specific_file_path)

    

        ### This might work for a genome object but not for a proteinsequenceset object
        # Read the results into a DataFrame
        # with open(specific_file_path, 'r') as csvfile:
        #     csvreader = csv.DictReader(csvfile)
        #     ontology_events = []
        #     for row in csvreader:
        #         ontology_event = {
        #             "event_id": "your_event_id",
        #             "description": "Protein annotation event",
        #             "ontology_id": row['Prediction'],
        #             "method": "SnekmerLearnApply_annotation",
        #             "method_version": "1.0",
        #             "timestamp": datetime.now().isoformat(),
        #             "feature_types": {},
        #             "ontology_terms": {
        #                 row['index']: [
        #                     {
        #                         "term": row['Prediction'],
        #                         "modelseed_ids": [],
        #                         "evidence_only": 0,
        #                         "evidence": {
        #                             "reference": ("protein_annotation", "SnekmerLearnApply"),
        #                             "scores": {"confidence": row['Confidence']}
        #                         }
        #                     }
        #                 ]
        #             }
        #         }
        #         ontology_events.append(ontology_event)
        # logging.info(ontology_events)




        # Use the cb_annotation_ontology_api to add the ontology events
        # ontology_api = cb_annotation_ontology_api(url=self.callback_url, token=os.environ.get('KB_AUTH_TOKEN'))
        # params = {
        #     "input_ref": object_refs[0]['ref'],
        #     "events": ontology_events,
        #     "output_name": "updated_protein_set_with_annotations",
        #     "output_workspace": params['workspace_name'],
        #     "save": 1,
        # }
        # # Call the method to add ontology events
        # update_result = ontology_api.add_annotation_ontology_events(params)
        # updated_protein_set_ref = update_result.get('output_ref')
        
        
        # New Stuff End


        #Attempt 2:
            
        # def generate_event_id():
        #     return "event_" + datetime.now().strftime("%Y%m%d%H%M%S")

        # # Read predictions and update ProteinSequence objects
        # with open(specific_file_path, 'r') as csvfile:
        #     csvreader = csv.DictReader(csvfile)
        #     for row in csvreader:
        #         # Find the matching ProteinSequence by 'id'
        #         for protein_sequence in protein_seq_set['sequences']:
        #             if protein_sequence['id'] == row['id']:  # Assuming 'id' column matches
        #                 # Update the ontology_terms for this ProteinSequence
        #                 ontology_term = row['Prediction']  # Assuming 'Prediction' column contains the ontology term
        #                 if ontology_term not in protein_sequence['ontology_terms']:
        #                     protein_sequence['ontology_terms'][ontology_term] = []
        #                 # Append an event index to the list for this term
        #                 # Here, you might need to adapt this based on how you handle ontology_events
        #                 protein_sequence['ontology_terms'][ontology_term].append(generate_event_id())

        # # Optionally, add a generic ontology event to the ontology_events list if needed
        # ontology_event = {
        #     "id": generate_event_id(),
        #     "ontology_ref": "your_ontology_ref",  # Replace with actual reference
        #     "method": "SnekmerLearnApply_annotation",
        #     "method_version": "1.0",
        #     "timestamp": datetime.now().isoformat(),
        #     "description": "Protein annotation based on SnekmerLearnApply predictions"
        # }
        # protein_seq_set['ontology_events'].append(ontology_event)

        # THIS WORKS - Do Not REMOVE
        
        #         
        # with open(specific_file_path, 'r') as csvfile:
        #     csvreader = csv.DictReader(csvfile)
        #     all_predictions = {}
        #     for row in csvreader:
        #         all_predictions[row['index']] = {
        #             "prediction": row['Prediction'],
        #             "score":row['Score'],
        #             "delta":row['delta'],
        #             "confidence":row['Confidence']
        #         }
        # sequences = protein_seq_set['data'][0]['data']['sequences']


        # for i,item in enumerate(sequences):
        #     if item["id"] in all_predictions:
        #     # Assuming all_predictions[item["id"]]["prediction"] gives a string like "Ribulokinase (EC 2.7.1.16)"
        #         prediction = all_predictions[item["id"]]["prediction"]
        #         ref_id = str(params['workspace_id']) + "." + str(i)
        #         # item["ontology_terms"] = {ref_id: {"term": [prediction]}}
        #         item["ontology_terms"] = {prediction: {"term": []}}
                
                
        # modified_data = protein_seq_set['data'][0]['data']


        # logging.info("New Protein Set \n\n\n")
        # logging.info(protein_seq_set)


        # object_name = 'my_protein_fasta2_with_Snekmer_annotations'
        # object_type = 'KBaseSequences.ProteinSequenceSet-1.0'


        # # Save the modified object
        # save_params_mod = {
        #     'workspace': params['workspace_name'],
        #     'objects': [{
        #         'type': object_type,
        #         'data': modified_data,
        #         'name': object_name
        #     }]
        # }



        # result = self.wsClient.save_objects(save_params_mod)
        # logging.info("Object saved successfully:")
        # logging.info("result")
        
        # saved_object_info = result[0]
        # object_id = str(saved_object_info[0])  # Object ID
        # workspace_id = params['workspace_id']  # Assuming this is the workspace ID
        # workspace_ref = "{}/{}".format(workspace_id, object_id)  # Workspace reference
        # ABOVE CODE WORKS - DO NOT REMOVE
    
        
        with open(specific_file_path, 'r') as csvfile:
            csvreader = csv.DictReader(csvfile)
            all_predictions = {}
            for row in csvreader:
                all_predictions[row['index']] = {
                    "prediction": row['Prediction'],
                    "score":row['Score'],
                    "delta":row['delta'],
                    "confidence":row['Confidence']
                }
        sequences = protein_seq_set['data'][0]['data']['sequences']
        
        logging.info(sequences)
        
        events = []
        for i,item in enumerate(sequences):
            if item["id"] in all_predictions:
                prediction = all_predictions[item["id"]]["prediction"]
                confidence = all_predictions[item["id"]]["confidence"]
                index = item["id"]
                # ref_id = str(params['workspace_id']) + "." + str(i)
                # item["ontology_terms"] = {ref_id: {"term": [prediction]}}
                item["ontology_terms"] = {prediction: {"term": []}}
                events.append({
                    "ontology_id" : "TIGR",
                    "description" : "TIGR annotations with Snekmer Apply",
                    "method_version" : "1.0",
                    "method" : "Snekmer Apply",
                    "timestamp" : datetime.now().strftime("%Y.%m.%d-%I:%M:%S%p"),
                    "ontology_terms":{ index : [
                        {
                            "term" : prediction,
                            "evidence" : {"scores":{"probability":confidence}}
                        }
                    ]
                    }
                })
        
        
        output = cb_annotation_ontology_api.add_annotation_ontology_events(self,params={
            "input_ref": object_refs['ref'], #Name of your input object
            "input_workspace":params['workspace_id'],#Workspace with your input object
            "output_name":"New_ProteinSetObj_CH_method",#Name to which the modified object should be saved
            "output_workspace":params['workspace_id'],#Workspace where output should be saved
            "clear_existing":0,#Set to 1 to clear existing annotations (don’t do this)
            "overwrite_matching":1,#Overwrites annotations for matching event IDs
            "save":1,#Set to one to save the output object
            "events":events
            })
                
             
        logging.info(output)   
        logging.info(output[0])
        saved_object_info = output[0]
        object_id = str(saved_object_info[0])
        workspace_ref = "{}/{}".format(params['workspace_id'], object_id)  # Workspace reference

    
        # Setup output directory for the ZIP file
        output_directory = os.path.join(self.shared_folder, str(uuid.uuid4()))
        os.makedirs(output_directory)
        run_date = datetime.now().strftime("%Y.%m.%d-%I:%M:%S%p")
        result_name = "KmerSummaryOutput" + run_date + ".zip"
        result_file = os.path.join(output_directory, result_name)

        # Zip the specific output file
        with zipfile.ZipFile(result_file, 'w', zipfile.ZIP_DEFLATED, allowZip64=True) as zip_file:
            # Add the specific file to the zip archive, adjust the arcname to change its name within the archive if needed
            zip_file.write(specific_file_path, os.path.basename(specific_file_path))


        # Prepare report parameters with the zipped file
        report_params = {
            'message': text_message,
            'workspace_name': workspace_name,
            'objects_created': [
                {
                    'ref': workspace_ref,
                    'description': 'Updated protein set with new ontologies and annotations'
                }
            ],
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
