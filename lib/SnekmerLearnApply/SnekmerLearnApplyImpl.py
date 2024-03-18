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
        Executes the Snekmer Apply function, which is designed to annotate biological sequences with ontology terms (TIGRFams, Pfam, PANTHER)
        based on their k-mer profiles. The method supports processing both protein and genome inputs. It performs several key operations:
        
        1. Validates input parameters to ensure the presence of either protein or genome inputs.
        2. Dynamically constructs object references based on the input type and fetches corresponding sequence data from a workspace.
        3. Writes the fetched sequences to FASTA files, which are then used as input for the Snekmer application.
        4. Decompresses necessary reference files for Snekmer analysis based on the specified family (TIGRFams, Pfam, PANTHER).
        5. Executes the Snekmer application to generate k-mer based annotations for the input sequences.
        6. Parses the Snekmer output to annotate the original sequence objects with the new ontology terms.
        7. Packages the Snekmer output files into a ZIP archive for easy download.
        8. Generates a detailed report through the KBaseReport service, summarizing the annotation results and providing links to the output files.
        
        Parameters:
            ctx (dict): A context object containing information about the runtime environment, including user authentication.
            params (dict): A dictionary containing method input parameters. Expected keys include:
                - workspace_name: Name of the workspace where the output will be saved.
                - protein_input or genome_input: Lists of object references to protein or genome data.
                - input_type: A string indicating the type of input ('protein' or 'genome').
                - family: The ontology family to use for annotation (e.g., 'TIGRFams', 'Pfam', 'PANTHER').
        
        Returns:
            dict: A dictionary with keys 'report_name' and 'report_ref', referring to the name and reference of the generated report.
        
        Raises:
            ValueError: If neither protein_input nor genome_input are provided in the params, or if other expected parameters are missing.
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN run_SnekmerLearnApply

        logging.info('Starting Snekmer Apply function. Params=' + pformat(params))
        logging.info('Validating parameters.')

        # check inputs
        workspace_name = params['workspace_name']
        run_type = []
        
        # Check for the presence of input keys
        if 'protein_input' not in params and 'genome_input' not in params:
            raise ValueError('Neither protein_input nor genome_input found in parameters.')

        # Validate that inputs match the specified input_type
        if 'input_type' not in params:
            raise ValueError("Input type not specified. Please include 'input_type' with value 'protein' or 'genome'.")

        input_type = params['input_type']
        if input_type not in ['protein', 'genome']:
            raise ValueError("Invalid 'input_type' specified. It must be either 'protein' or 'genome'.")

        if input_type == 'protein':
            if 'protein_input' not in params or not params['protein_input']:
                # protein_input must not be empty if input_type is protein
                raise ValueError("'protein_input' is required and cannot be empty when 'input_type' is 'protein'.")
        elif input_type == 'genome':
            if 'genome_input' not in params or not params['genome_input']:
                # genome_input must not be empty if input_type is genome
                raise ValueError("'genome_input' is required and cannot be empty when 'input_type' is 'genome'.")

        if len(params['protein_input']) > 0:
            run_type.append("protein")
            protein_input = params['protein_input']
        if len(params['genome_input']) > 0:
            run_type.append("genome")
            genome_input  = params['genome_input']
        object_refs = []
            
        
        # Here we constuct a fasta file from a Protein Sequence Set Object
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
        
        
        # Here we constuct a fasta file from a Genome Object
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
            
        
        # Here we Prepare Snekmer Apply based on the selected Family Type.
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
                     
                
        # Here we run Snekmer Apply on the newly Generated Fasta File
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
    
        target_dir = os.path.join(cwd, "output", "apply")
        file_paths = [f for f in os.listdir(target_dir) if os.path.isfile(os.path.join(target_dir, f))]
        
        logging.info(file_paths)

        all_predictions = {}
        for file in file_paths:
            with open(os.path.join(cwd, "output", "apply", file), 'r') as csvfile:
                csvreader = csv.DictReader(csvfile)
                for row in csvreader:
                    all_predictions[row['index']] = {
                        "prediction": row['Prediction'],
                        "score":row['Score'],
                        "delta":row['delta'],
                        "confidence":row['Confidence']
                    }
                
        
        # Prepare for returning Object
        logging.info("Check object_refs here")
        logging.info(object_refs)
        object_id_list = []
        workspace_ref_list = []            
        fam_map = {"TIGRFams": "TIGR", "Pfam": "PF", "PANTHER":"PTHR"}
        family_type = params["family"]
        ontology_id = fam_map.get(family_type)
        description_prefix = family_type + " annotations with Snekmer Apply"

        # Here we return the results to a Protein Sequence Set Object
        if "protein" in run_type and "protein" == params["input_type"]:
            for seq_obj_num, ref in enumerate(protein_input):
                    ontology_api = None
                    ontology_api = cb_annotation_ontology_api(url=self.callback_url, token=os.environ.get('KB_AUTH_TOKEN'))
                    # Fetch the object for the current reference
                    protein_seq_set = self.wsClient.get_objects2({'objects': [{"ref": ref}]})
                    object_name = protein_seq_set['data'][0]['info'][1] + "_Snekmer_LA_" + family_type
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
                        "input_ref": ref, #Name of your input object
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

        # Here we return the results to a Genome Object
        if "genome" in run_type and "genome" == params["input_type"]:
            for seq_obj_num, ref in enumerate(genome_input):
                    logging.info("Genome ref (from iterating through line 323 for loop): " + ref)
                    ontology_api = None
                    ontology_api = cb_annotation_ontology_api(url=self.callback_url, token=os.environ.get('KB_AUTH_TOKEN'))
                    # Fetch the object for the current reference
                    genome_seq_set = self.wsClient.get_objects2({'objects': [{"ref": ref}]})
                    object_name = genome_seq_set['data'][0]['info'][1] + "_Snekmer_LA_" + family_type
                    sequences = genome_seq_set['data'][0]['data']["cdss"]
                    logging.info("object_name: ")
                    logging.info(object_name)
                    # logging.info(sequences)
                    
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
                        "input_ref": ref, #Name of your input object
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
        
        # Prepare the Report Parameters
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
