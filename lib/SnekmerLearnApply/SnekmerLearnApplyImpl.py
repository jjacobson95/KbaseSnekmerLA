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
        run_type = []
        if 'protein_input' not in params and 'genome_input' not in params:
            raise ValueError('Neither protein_input nor genome_input found')
        if 'protein_input' in params:
            run_type.append("protein")
            protein_input = params['protein_input']
        if 'genome_input' in params:
            run_type.append("genome")
            genome_input  = params['genome_input']
            
        # protein_input = params['protein_input']
        # genome_input  = params['genome_input']

        ###
        ###  Protein Path - Parse into fasta file(s)
        ### 
        
        fasta_index = 0
        if "protein" in run_type:
            object_refs = [{'ref': ref} for ref in protein_input]
            # Fetch the objects using get_objects2
            protein_seq_set = self.wsClient.get_objects2({'objects': object_refs})
            text_message = '\n'.join(params['protein_input'])
            
            logging.info(object_refs)
            # logging.info(protein_seq_set)
            logging.info(sys.version)
            logging.info(protein_seq_set)
        
            output_dir = "input"
            # output_file_name = "output1.fasta"
            # output_file_path = os.path.join(output_dir, output_file_name)

            # Ensure the output directory exists
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
        
            #old working version for concat all into 1
            # with open(output_file_path, 'w') as fasta_file:
            #     # Iterate over each fetched object in the response
            #     for obj in protein_seq_set['data']:
            #         # Extract sequences data for each object
            #         sequences_data = obj['data']['sequences']
            #         # Iterate over each sequence
            #         for seq in sequences_data:
            #             # Construct the header with sequence ID and description
            #             header = ">{} {}".format(seq['id'], seq['description'])
            #             # Write the header and sequence to the FASTA file
            #             fasta_file.write("{}\n{}\n".format(header, seq['sequence']))
            
            for index, ref in enumerate(protein_input):
                # Fetch the object for the current reference
                protein_seq_set = self.wsClient.get_objects2({'objects': [ref]})
                
                # Generate a unique output file name for the current input
                fasta_index += 1
                output_file_name = "output_" + str(fasta_index) + ".fasta"
                output_file_path = os.path.join(output_dir, output_file_name)
                
                with open(output_file_path, 'w') as fasta_file:
                    # Iterate over each fetched object in the response (should be only one in this case)
                    for obj in protein_seq_set['data']:
                        # Extract sequences data for the object
                        sequences_data = obj['data']['sequences']
                        # Iterate over each sequence
                        for seq in sequences_data:
                            # Construct the header with sequence ID and description
                            header = ">{} {}".format(seq['id'], seq['description'])
                            # Write the header and sequence to the FASTA file
                            fasta_file.write("{}\n{}\n".format(header, seq['sequence']))

            # Read the FASTA file and log the first 40 sequences
            logging.info("Genome Object(s) Received and parsed into Fasta. ")
            sequences = list(SeqIO.parse(output_file_path, 'fasta'))
            for seq_record in sequences[0:40]:
                logging.info("ID: {}".format(seq_record.id))
                logging.info("Description: {}".format(seq_record.description))
                logging.info("Sequence: {}\n".format(str(seq_record.seq)[:40]))
            
            
        ###
        ###  Genome Path - Parse into fasta file(s)
        ### 
        
        if "genome" in run_type:
            # object_refs = [{'ref': ref} for ref in genome_input]
            
            # logging.info(object_refs)
            # # logging.info(protein_seq_set)
            # logging.info(sys.version)
            # logging.info(protein_seq_set)
            text_message = '\n'.join(params['genome_input'])
            output_dir = "input"
            # output_file_name = "output1.fasta"
            # output_file_path = os.path.join(output_dir, output_file_name)

            # Ensure the output directory exists
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
        
        
            for index, ref in enumerate(genome_input):
                    # Fetch the object for the current reference
                    genome_seq_set = self.wsClient.get_objects2({'objects': [ref]})
                    
                    # Generate a unique output file name for the current input
                    fasta_index += 1
                    output_file_name = "output_" + str(fasta_index) + ".fasta"
                    output_file_path = os.path.join(output_dir, output_file_name)
            
                    with open(output_file_path, 'w') as fasta_file:
                        # Iterate over each fetched object in the response (should be only one in this case)
                        for obj in genome_seq_set['data']:
                            # Extract sequences data for the object
                            sequences_data = obj['data']["cdss"]
                            # Iterate over each sequence
                            for seq in sequences_data:
                                # Construct the header with sequence ID and description
                                header = ">{} {}".format(seq['id'], seq['function'])
                                # Write the header and sequence to the FASTA file
                                fasta_file.write("{}\n{}\n".format(header, seq['protein_translation']))

            # Read the FASTA file and log the first 40 sequences
            logging.info("Genome Object(s) Received and parsed into Fasta. ")
            sequences = list(SeqIO.parse(output_file_path, 'fasta'))
            for seq_record in sequences[0:40]:
                logging.info("ID: {}".format(seq_record.id))
                logging.info("Description: {}".format(seq_record.description))
                logging.info("Sequence: {}\n".format(str(seq_record.seq)[:40]))
            
        
        # Run Snekmer 
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
    
        # specific_file_path = os.path.join(cwd, "output", "apply", "kmer-summary-output1.csv")
        target_dir = os.path.join(cwd, "output", "apply")
        file_paths = [f for f in os.listdir(target_dir) if os.path.isfile(os.path.join(target_dir, f))]
        
        logging.info(file_paths)

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
        
                
                
                
                
        for file in file_paths:
            with open(file, 'r') as csvfile:
                csvreader = csv.DictReader(csvfile)
                all_predictions = {}
                for row in csvreader:
                    all_predictions[row['index']] = {
                        "prediction": row['Prediction'],
                        "score":row['Score'],
                        "delta":row['delta'],
                        "confidence":row['Confidence']
                    }
                    
        # sequences = protein_seq_set['data'][0]['data']['sequences']
        saved_object_info_list = []
        object_id_list = []
        workspace_ref_list = []
        if "protein" in run_type:
            for seq_obj_num,sequences in enumerate(protein_seq_set['data']):

                for i,item in enumerate(sequences['data']['sequences'] ):
                    if item["id"] in all_predictions:
                        prediction = all_predictions[item["id"]]["prediction"]
                        confidence = all_predictions[item["id"]]["confidence"]
                        index = item["id"]
                    # Assuming all_predictions[item["id"]]["prediction"] gives a string like "Ribulokinase (EC 2.7.1.16)"
                        prediction = all_predictions[item["id"]]["prediction"]
                        ref_id = str(params['workspace_id']) + "." + str(i)
                        item["ontology_terms"] = {prediction: {"evidence" : confidence}}
                        
                        
                        # item["ontology_terms"] = {index:                       {
                        #             "term" : prediction,
                        #             "evidence" : {"scores":{"probability":confidence}}
                        #         }}
                        
                        
                        # "ontology_terms":{ index : [
                        #         {
                        #             "term" : prediction,
                        #             "evidence" : {"scores":{"probability":confidence}}
                        #         }
                        
                        
                        
                modified_data = protein_seq_set['data'][seq_obj_num]['data']


            logging.info("New Protein Set \n\n\n")
            logging.info(modified_data)
            logging.info(protein_seq_set['data'][seq_obj_num]['info'][1])

            object_name = protein_seq_set['data'][seq_obj_num]['info'][1]
            object_type = 'KBaseSequences.ProteinSequenceSet-1.0'

            # Save the modified object
            save_params_mod = {
                'workspace': params['workspace_name'],
                'objects': [{
                    'type': object_type,
                    'data': modified_data,
                    'name': object_name
                }]
            }


            result = self.wsClient.save_objects(save_params_mod)
            logging.info("Object saved successfully:")
            logging.info("result")
            
            saved_object_info_list.append(result[0])
            object_id_list.append(str(result[0][0]))
            # object_id = str(result[0][0])  # Object ID
            workspace_id = params['workspace_id']  # Assuming this is the workspace ID
            # workspace_ref = "{}/{}".format(workspace_id, result[0][0])  # Workspace reference
            workspace_ref_list.append("{}/{}".format(workspace_id, result[0][0]))  # Workspace reference

        # ABOVE CODE WORKS - DO NOT REMOVE
    
        
    
        #### Below code works on my side - but chris needs to fix something. "Location" Keyerror - I think this works for genome but not protein.
        # ## Do not delete
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
        
        # logging.info(sequences)
        
        # events = []
        # for i,item in enumerate(sequences):
        #     if item["id"] in all_predictions:
        #         prediction = all_predictions[item["id"]]["prediction"]
        #         confidence = all_predictions[item["id"]]["confidence"]
        #         index = item["id"]
        #         # ref_id = str(params['workspace_id']) + "." + str(i)
        #         # item["ontology_terms"] = {ref_id: {"term": [prediction]}}
        #         item["ontology_terms"] = {prediction: {"term": []}}
        #         events.append({
        #             "ontology_id" : "TIGR",
        #             "description" : "TIGR annotations with Snekmer Apply",
        #             "method_version" : "1.0",
        #             "method" : "Snekmer Apply",
        #             "timestamp" : datetime.now().strftime("%Y.%m.%d-%I:%M:%S%p"),
        #             "ontology_terms":{ index : [
        #                 {
        #                     "term" : prediction,
        #                     "evidence" : {"scores":{"probability":confidence}}
        #                 }
        #             ]
        #             }
        #         })
        
        # ontology_api = cb_annotation_ontology_api(url=self.callback_url, token=os.environ.get('KB_AUTH_TOKEN'))
        # output = ontology_api.add_annotation_ontology_events(params={
        #     "input_ref": object_refs[0]['ref'], #Name of your input object
        #     "input_workspace":params['workspace_id'],#Workspace with your input object
        #     "output_name":"New_ProteinSetObj_CH_method",#Name to which the modified object should be saved
        #     "output_workspace":params['workspace_id'],#Workspace where output should be saved
        #     "clear_existing":0,#Set to 1 to clear existing annotations (don’t do this)
        #     "overwrite_matching":1,#Overwrites annotations for matching event IDs
        #     "save":1,#Set to one to save the output object
        #     "events":events
        #     })
                
             
        # logging.info(output)   
        # logging.info(output[0])
        # saved_object_info = output[0]
        # object_id = str(saved_object_info[0])
        # workspace_ref = "{}/{}".format(params['workspace_id'], object_id)  # Workspace reference
        
        #DO NOT DELETE above Code - needs fix from chris
    
    
        # Setup output directory for the ZIP file
        output_directory = os.path.join(self.shared_folder, str(uuid.uuid4()))
        os.makedirs(output_directory)
        run_date = datetime.now().strftime("%Y.%m.%d-%I:%M:%S%p")
        result_name = "KmerSummaryOutput" + run_date + ".zip"
        result_file = os.path.join(output_directory, result_name)
        
        def zipdir(path, ziph):
            # ziph is zipfile handle
            for root, dirs, files in os.walk(path):
                for file in files:
                    # Create a relative path for files to keep the directory structure
                    relative_path = os.path.relpath(os.path.join(root, file), os.path.join(path, '..'))
                    ziph.write(os.path.join(root, file), relative_path)

        # Zip the specific output file
        with zipfile.ZipFile(result_file, 'w', zipfile.ZIP_DEFLATED, allowZip64=True) as zip_file:
            # Add the specific file to the zip archive, adjust the arcname to change its name within the archive if needed
            zipdir(target_dir, zip_file)

        object_list = []
        for workspace_ref in workspace_ref_list:
            object_list.append({
                    'ref': workspace_ref,
                    'description': 'Updated protein set with new ontologies and annotations'
                })

        # Prepare report parameters with the zipped file
        report_params = {
            'message': text_message,
            'workspace_name': workspace_name,
            'objects_created': object_list,
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
