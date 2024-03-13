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
import gzip
import shutil

from installed_clients.AssemblyUtilClient import AssemblyUtil
from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.WorkspaceClient import Workspace as workspaceService
from installed_clients.cb_annotation_ontology_apiClient import cb_annotation_ontology_api

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

        logging.info('Starting Snekmer Apply function. Params=' + pformat(params))
        logging.info('Validating parameters.')

        # check inputs
        workspace_name = params['workspace_name']
        run_type = []
        
        if 'protein_input' not in params and 'genome_input' not in params:
            raise ValueError('Neither protein_input nor genome_input found')
        if len(params['protein_input']) > 0:
            run_type.append("protein")
            protein_input = params['protein_input']
        if len(params['genome_input']) > 0:
            run_type.append("genome")
            genome_input  = params['genome_input']
        object_refs = []
            
        
        fasta_index = 0
        if "protein" in run_type and "protein" == params["input_type"]:
            object_refs.extend([{'ref': ref} for ref in protein_input])
            text_message = '\n'.join(params['protein_input'])

            logging.info(sys.version)
        
            output_dir = "input"

            # Ensure the output directory exists
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
            
            for index, ref in enumerate(protein_input):
                # Fetch the object for the current reference
                protein_seq_set = self.wsClient.get_objects2({'objects': [{"ref": ref}]})
                
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
        
        if "genome" in run_type and "genome" == params["input_type"]:
            object_refs.extend([{'ref': ref} for ref in genome_input])
            
            text_message = '\n'.join(params['genome_input'])
            output_dir = "input"

            # Ensure the output directory exists
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
    
            for index, ref in enumerate(genome_input):
                    # Fetch the object for the current reference
                    genome_seq_set = self.wsClient.get_objects2({'objects': [{"ref": ref}]})

                    # Generate a unique output file name for the current input
                    fasta_index += 1
                    output_file_name = "output_" + str(fasta_index) + ".fasta"
                    output_file_path = os.path.join(output_dir, output_file_name)

                    with open(output_file_path, 'w') as fasta_file:
                        for obj in genome_seq_set['data']:
                            sequences_data = obj['data']["cdss"]
                            for seq in sequences_data:
                                if 'functions' in seq:
                                    function_or_functions = seq['functions']
                                elif 'function' in seq:
                                    function_or_functions = seq['function']
                                else:
                                    function_or_functions = 'Unknown Function' 
                                header = ">{} {}".format(seq['id'], function_or_functions)
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
        
        

        def decompress_and_move(source_path, destination_path):
            # Adjust the destination path to have the correct filename (remove .gz)
            destination_path_unzipped = destination_path.rstrip('.gz')
            
            # Decompress the .gz file and write the contents to the new location
            with gzip.open(source_path, 'rb') as f_in:
                with open(destination_path_unzipped, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)

        if params["family"] in ["TIGRFams", "Pfam", "PANTHER"]:
            source_confidence = "data/{0}-global-confidence-scores.csv.gz".format(params["family"])
            destination_confidence = "confidence/{0}-global-confidence-scores.csv".format(params["family"])
            source_counts = "data/{0}-kmer-counts-total.csv.gz".format(params["family"])
            destination_counts = "counts/{0}-kmer-counts-total.csv".format(params["family"])

        # Make dirs if not present
        os.makedirs(os.path.dirname(destination_confidence))
        os.makedirs(os.path.dirname(destination_counts))
            
        if os.path.exists(source_confidence):
            decompress_and_move(source_confidence, destination_confidence)
        
        # Decompress and move kmer counts total file
        if os.path.exists(source_counts):
            decompress_and_move(source_counts, destination_counts)
                     
                
                
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


        for file in file_paths:
            with open(os.path.join(cwd, "output", "apply", file), 'r') as csvfile:
                csvreader = csv.DictReader(csvfile)
                all_predictions = {}
                for row in csvreader:
                    all_predictions[row['index']] = {
                        "prediction": row['Prediction'],
                        "score":row['Score'],
                        "delta":row['delta'],
                        "confidence":row['Confidence']
                    }
                
        
        logging.info("Check object_refs here")
        logging.info(object_refs)
        object_id_list = []
        workspace_ref_list = []    
        ontology_api = cb_annotation_ontology_api(url=self.callback_url, token=os.environ.get('KB_AUTH_TOKEN'))
        
        
        
        fam_map = {"TIGRFams": "TIGR", "Pfam": "PF", "PANTHER":"PTHR"}
        family_type = params["family"]
        ontology_id = fam_map.get(family_type)
        description_prefix = family_type + " annotations with Snekmer Apply"

        if "protein" in run_type and "protein" == params["input_type"]:
            for seq_obj_num, ref in enumerate(protein_input):
                    # Fetch the object for the current reference
                    protein_seq_set = self.wsClient.get_objects2({'objects': [{"ref": ref}]})
                    object_name = protein_seq_set['data'][0]['info'][1] + "_Annotated_with_Snekmer_Apply"
                    sequences = protein_seq_set['data'][0]['data']['sequences']
                    
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
                                "ontology_id" : ontology_id,
                                "description" : description_prefix,
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
                
                    result = ontology_api.add_annotation_ontology_events(params={
                        "input_ref": object_refs[seq_obj_num]['ref'], #Name of your input object
                        "input_workspace":params['workspace_id'],#Workspace with your input object
                        "output_name":object_name,#Name to which the modified object should be saved
                        "output_workspace":params['workspace_id'],#Workspace where output should be saved
                        "clear_existing":0,#Set to 1 to clear existing annotations (don’t do this)
                        "overwrite_matching":1,#Overwrites annotations for matching event IDs
                        "save":1,#Set to one to save the output object
                        "events":events
                        })

                    
                    logging.info("Object saved successfully:")
                    logging.info(result)
                    
                    object_id_list.append(str(result['output_ref']))
                    workspace_ref_list.append(str(result['output_ref']))  # Workspace reference


        
        fam_map = {"TIGRFams": "TIGR", "Pfam": "PF", "PANTHER":"PTHR"}
        family_type = params["family"]
        ontology_id = fam_map.get(family_type)
        description_prefix = family_type + " annotations with Snekmer Apply"

        if "genome" in run_type and "genome" == params["input_type"]:
            for seq_obj_num, ref in enumerate(genome_input):
                    # Fetch the object for the current reference
                    genome_seq_set = self.wsClient.get_objects2({'objects': [{"ref": ref}]})
                    object_name = genome_seq_set['data'][0]['info'][1] + "_Annotated_with_Snekmer_Apply"
                    sequences = genome_seq_set['data'][0]['data']["cdss"]
                    logging.info(object_name)
                    logging.info(sequences)
                    
                    events = []
                    for i,item in enumerate(sequences):
                        if item["id"] in all_predictions:
                            prediction = all_predictions[item["id"]]["prediction"]
                            confidence = all_predictions[item["id"]]["confidence"]
                            index = item["id"]
                            item["ontology_terms"] = {prediction: {"term": []}}
                            events.append({
                                "ontology_id" : ontology_id,
                                "description" : description_prefix,
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
                
                    result = ontology_api.add_annotation_ontology_events(params={
                        "input_ref": object_refs[seq_obj_num]['ref'], #Name of your input object
                        "input_workspace":params['workspace_id'],#Workspace with your input object
                        "output_name":object_name,#Name to which the modified object should be saved
                        "output_workspace":params['workspace_id'],#Workspace where output should be saved
                        "clear_existing":0,#Set to 1 to clear existing annotations (don’t do this)
                        "overwrite_matching":1,#Overwrites annotations for matching event IDs
                        "save":1,#Set to one to save the output object
                        "events":events
                        })
                    logging.info("Object saved successfully:")
                    logging.info(result)
                    
                    object_id_list.append(str(result['output_ref']))
                    workspace_ref_list.append(str(result['output_ref']))  # Workspace reference

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
                    'description': 'Updated object with new ontologies and annotations'
                })
    
        text_message = "<b>All object(s) above have been successfully annotated with {0} using Snekmer Apply.</b><br><br>The Zipfile below contains output files with the following columns:<br><b>Index</b>: Coding Sequence ID<br><b>Prediction</b>: Predicted {0} Annotation<br><b>Score</b>: Cosine Similarity Score between coding sequence and nearest annotation in the Kmer-Association Matrix.<br><b>Delta</b>: Difference between the top two cosine similarity scores. A greater difference indicates a higher resolution and confidence.<br><b>Confidence</b>: Approximate probability of the prediction correctness.".format(params["family"])
        # text_message = "<b>All object(s) above have been successfully annotated with TIGRFAMs using Snekmer Apply.</b><br><br>The Zipfile below contains output files with the following columns:<br><b>Index</b>: Coding Sequence ID<br><b>Prediction</b>: Predicted TIGRFAMS Annotation<br><b>Score</b>: Cosine Similarity Score between coding sequence and nearest annotation in the Kmer-Association Matrix.<br><b>Delta</b>: Difference between the top two cosine similarity scores. A greater difference indicates a higher resolution and confidence.<br><b>Confidence</b>: Approximate probability of the prediction correctness."
        
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
        logging.info(report_params)

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
