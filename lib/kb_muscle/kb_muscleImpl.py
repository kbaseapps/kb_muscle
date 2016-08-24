#BEGIN_HEADER
import os
import sys
import shutil
import hashlib
import subprocess
import requests
import re
import traceback
import uuid
from datetime import datetime
from pprint import pprint, pformat
import numpy as np
import gzip

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein
from biokbase.workspace.client import Workspace as workspaceService
from requests_toolbelt import MultipartEncoder  # added
from biokbase.AbstractHandle.Client import AbstractHandle as HandleService  # added

# KBase Data API
#from doekbase.data_api.annotation.genome_annotation.api import GenomeAnnotationAPI as GenomeAnnotationAPI
    
# Standard setup for accessing Data API
#services = {"workspace_service_url": "https://ci.kbase.us/services/ws/",
#            "shock_service_url": "https://ci.kbase.us/services/shock-api/"}
#token = os.environ["KB_AUTH_TOKEN"]


# silence whining
import requests
requests.packages.urllib3.disable_warnings()

#END_HEADER


class kb_muscle:
    '''
    Module Name:
    kb_muscle

    Module Description:
    ** A KBase module: kb_muscle
**
** This module runs MUSCLE to make MSAs of either DNA or PROTEIN sequences.  "MUSCLE nuc" will build nucleotide alignments, even for protein coding genes.  "MUSCLE prot" will build protein sequence alignments, and will ignore any features that do not code for proteins.
**
    '''

    ######## WARNING FOR GEVENT USERS #######
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    #########################################
    #BEGIN_CLASS_HEADER
    workspaceURL = None
    shockURL = None
    handleURL = None

    MUSCLE_bin = '/kb/module/muscle/bin/muscle'

    # target is a list for collecting log messages
    def log(self, target, message):
        # we should do something better here...
        if target is not None:
            target.append(message)
        print(message)
        sys.stdout.flush()

    def get_single_end_read_library(self, ws_data, ws_info, forward):
        pass

    def get_feature_set_seqs(self, ws_data, ws_info):
        pass

    def KBase_data2file_GenomeAnnotation2Fasta(self, ws_data, ws_info):
        pass

    def get_genome_set_feature_seqs(self, ws_data, ws_info):
        pass

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.workspaceURL = config['workspace-url']
        self.shockURL = config['shock-url']
        self.handleURL = config['handle-service-url']
        self.scratch = os.path.abspath(config['scratch'])
        # HACK!! temporary hack for issue where megahit fails on mac because of silent named pipe error
        #self.host_scratch = self.scratch
        self.scratch = os.path.join('/kb','module','local_scratch')
        # end hack
        if not os.path.exists(self.scratch):
            os.makedirs(self.scratch)

        #END_CONSTRUCTOR
        pass


    # Helper script borrowed from the transform service, logger removed
    #
    def upload_file_to_shock(self,
                             console,  # DEBUG
                             shock_service_url = None,
                             filePath = None,
                             ssl_verify = True,
                             token = None):
        """
        Use HTTP multi-part POST to save a file to a SHOCK instance.
        """
        self.log(console,"UPLOADING FILE "+filePath+" TO SHOCK")

        if token is None:
            raise Exception("Authentication token required!")

        #build the header
        header = dict()
        header["Authorization"] = "Oauth {0}".format(token)
        if filePath is None:
            raise Exception("No file given for upload to SHOCK!")

        dataFile = open(os.path.abspath(filePath), 'rb')
        m = MultipartEncoder(fields={'upload': (os.path.split(filePath)[-1], dataFile)})
        header['Content-Type'] = m.content_type

        #logger.info("Sending {0} to {1}".format(filePath,shock_service_url))
        try:
            response = requests.post(shock_service_url + "/node", headers=header, data=m, allow_redirects=True, verify=ssl_verify)
            dataFile.close()
        except:
            dataFile.close()
            raise
        if not response.ok:
            response.raise_for_status()
        result = response.json()
        if result['error']:
            raise Exception(result['error'][0])
        else:
            return result["data"]


    def upload_SingleEndLibrary_to_shock_and_ws (self,
                                                 ctx,
                                                 console,  # DEBUG
                                                 workspace_name,
                                                 obj_name,
                                                 file_path,
                                                 provenance,
                                                 sequencing_tech):

        self.log(console,'UPLOADING FILE '+file_path+' TO '+workspace_name+'/'+obj_name)

        # 1) upload files to shock
        token = ctx['token']
        forward_shock_file = self.upload_file_to_shock(
            console,  # DEBUG
            shock_service_url = self.shockURL,
            filePath = file_path,
            token = token
            )
        #pprint(forward_shock_file)
        self.log(console,'SHOCK UPLOAD DONE')

        # 2) create handle
        self.log(console,'GETTING HANDLE')
        hs = HandleService(url=self.handleURL, token=token)
        forward_handle = hs.persist_handle({
                                        'id' : forward_shock_file['id'], 
                                        'type' : 'shock',
                                        'url' : self.shockURL,
                                        'file_name': forward_shock_file['file']['name'],
                                        'remote_md5': forward_shock_file['file']['checksum']['md5']})

        
        # 3) save to WS
        self.log(console,'SAVING TO WORKSPACE')
        single_end_library = {
            'lib': {
                'file': {
                    'hid':forward_handle,
                    'file_name': forward_shock_file['file']['name'],
                    'id': forward_shock_file['id'],
                    'url': self.shockURL,
                    'type':'shock',
                    'remote_md5':forward_shock_file['file']['checksum']['md5']
                },
                'encoding':'UTF8',
                'type':'fasta',
                'size':forward_shock_file['file']['size']
            },
            'sequencing_tech':sequencing_tech
        }
        self.log(console,'GETTING WORKSPACE SERVICE OBJECT')
        ws = workspaceService(self.workspaceURL, token=ctx['token'])
        self.log(console,'SAVE OPERATION...')
        new_obj_info = ws.save_objects({
                        'workspace':workspace_name,
                        'objects':[
                            {
                                'type':'KBaseFile.SingleEndLibrary',
                                'data':single_end_library,
                                'name':obj_name,
                                'meta':{},
                                'provenance':provenance
                            }]
                        })
        self.log(console,'SAVED TO WORKSPACE')

        return new_obj_info[0]

    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.workspaceURL = config['workspace-url']
        self.shockURL = config['shock-url']
        self.handleURL = config['handle-service-url']
        self.scratch = os.path.abspath(config['scratch'])
        # HACK!! temporary hack for issue where megahit fails on mac because of silent named pipe error
        #self.host_scratch = self.scratch
        self.scratch = os.path.join('/kb','module','local_scratch')
        # end hack
        if not os.path.exists(self.scratch):
            os.makedirs(self.scratch)

        #END_CONSTRUCTOR
        pass

    def MUSCLE_nuc(self, ctx, params):
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN MUSCLE_nuc
        console = []
        invalid_msgs = []
        self.log(console,'Running MUSCLE_nuc with params=')
        self.log(console, "\n"+pformat(params))
        report = ''
#        report = 'Running MUSCLE_nuc with params='
#        report += "\n"+pformat(params)


        #### do some basic checks
        #
        if 'workspace_name' not in params:
            raise ValueError('workspace_name parameter is required')
        if 'input_name' not in params:
            raise ValueError('input_name parameter is required')
        if 'output_name' not in params:
            raise ValueError('output_name parameter is required')


        #### Get the input_name object
        ##
        input_forward_reads_file_compression = None
        sequencing_tech = 'N/A'
        try:
            ws = workspaceService(self.workspaceURL, token=ctx['token'])
            objects = ws.get_objects([{'ref': params['workspace_name']+'/'+params['input_name']}])
            data = objects[0]['data']
            info = objects[0]['info']
            input_type_name = info[2].split('.')[1].split('-')[0]

            if input_type_name == 'SingleEndLibrary':
                input_type_namespace = info[2].split('.')[0]
                if input_type_namespace == 'KBaseAssembly':
                    file_name = data['handle']['file_name']
                elif input_type_namespace == 'KBaseFile':
                    file_name = data['lib']['file']['file_name']
                else:
                    raise ValueError('bad data type namespace: '+input_type_namespace)
                #self.log(console, 'INPUT_FILENAME: '+file_name)  # DEBUG
                if file_name[-3:] == ".gz":
                    input_forward_reads_file_compression = 'gz'
                if 'sequencing_tech' in data:
                    sequencing_tech = data['sequencing_tech']

        except Exception as e:
            traceback.format_exc()
            raise ValueError('Unable to fetch input_name object from workspace: ' + str(e))


        # Handle overloading (input_name can be SingleEndLibrary or FeatureSet)
        #
        if input_type_name == 'SingleEndLibrary':

            # DEBUG
            #for k in data:
            #    self.log(console,"SingleEndLibrary ["+k+"]: "+str(data[k]))

            try:
                if 'lib' in data:
                    input_forward_reads = data['lib']['file']
                elif 'handle' in data:
                    input_forward_reads = data['handle']
                else:
                    self.log(console,"bad structure for 'input_forward_reads'")
                    raise ValueError("bad structure for 'input_forward_reads'")

                ### NOTE: this section is what could be replaced by the transform services
                input_forward_reads_file_path = os.path.join(self.scratch,input_forward_reads['file_name'])
                input_forward_reads_file_handle = open(input_forward_reads_file_path, 'w', 0)
                self.log(console, 'downloading reads file: '+str(input_forward_reads_file_path))
                headers = {'Authorization': 'OAuth '+ctx['token']}
                r = requests.get(input_forward_reads['url']+'/node/'+input_forward_reads['id']+'?download', stream=True, headers=headers)
                for chunk in r.iter_content(1024):
                    input_forward_reads_file_handle.write(chunk)
                input_forward_reads_file_handle.close();
                self.log(console, 'done')
                ### END NOTE


                # remove carriage returns
                new_file_path = input_forward_reads_file_path+"-CRfree"
                new_file_handle = open(new_file_path, 'w', 0)
                input_forward_reads_file_handle = open(input_forward_reads_file_path, 'r', 0)
                for line in input_forward_reads_file_handle:
                    line = re.sub("\r","",line)
                    new_file_handle.write(line)
                input_forward_reads_file_handle.close();
                new_file_handle.close()
                input_forward_reads_file_path = new_file_path

                # convert FASTQ to FASTA (if necessary)
                new_file_path = input_forward_reads_file_path+".fna"
                new_file_handle = open(new_file_path, 'w', 0)
                if input_forward_reads_file_compression == 'gz':
                    input_forward_reads_file_handle = gzip.open(input_forward_reads_file_path, 'r', 0)
                else:
                    input_forward_reads_file_handle = open(input_forward_reads_file_path, 'r', 0)
                header = None
                last_header = None
                last_seq_buf = None
                last_line_was_header = False
                was_fastq = False
                for line in input_forward_reads_file_handle:
                    if line.startswith('>'):
                        break
                    elif line.startswith('@'):
                        was_fastq = True
                        header = line[1:]
                        if last_header != None:
                            new_file_handle.write('>'+last_header)
                            new_file_handle.write(last_seq_buf)
                        last_seq_buf = None
                        last_header = header
                        last_line_was_header = True
                    elif last_line_was_header:
                        last_seq_buf = line
                        last_line_was_header = False
                    else:
                        continue
                if last_header != None:
                    new_file_handle.write('>'+last_header)
                    new_file_handle.write(last_seq_buf)

                new_file_handle.close()
                input_forward_reads_file_handle.close()
                if was_fastq:
                    input_forward_reads_file_path = new_file_path

            except Exception as e:
                print(traceback.format_exc())
                raise ValueError('Unable to download single-end read library files: ' + str(e))

        # FeatureSet
        #
        elif input_type_name == 'FeatureSet':
            # retrieve sequences for features
            input_featureSet = data

            genome2Features = {}
            features = input_featureSet['elements']
            for fId in features.keys():
                genomeRef = features[fId][0]
                if genomeRef not in genome2Features:
                    genome2Features[genomeRef] = []
                genome2Features[genomeRef].append(fId)

            # export features to FASTA file
            input_forward_reads_file_path = os.path.join(self.scratch, params['input_name']+".fasta")
            self.log(console, 'writing fasta file: '+input_forward_reads_file_path)
            records = []
            for genomeRef in genome2Features:
                genome = ws.get_objects([{'ref':genomeRef}])[0]['data']
                these_genomeFeatureIds = genome2Features[genomeRef]
                for feature in genome['features']:
                    if feature['id'] in these_genomeFeatureIds:
                        #self.log(console,"kbase_id: '"+feature['id']+"'")  # DEBUG
                        record = SeqRecord(Seq(feature['dna_sequence']), id=feature['id'], description=genome['id'])
                        records.append(record)
            SeqIO.write(records, input_forward_reads_file_path, "fasta")


        # Missing proper input_input_type
        #
        else:
            raise ValueError('Cannot yet handle input_name type of: '+type_name)


        ### Construct the command
        #
        #  e.g. muscle -in <fasta_in> -out <fasta_out> -maxiters <n> -haxours <h>
        #
        muscle_cmd = [self.MUSCLE_bin]

        # check for necessary files
        if not os.path.isfile(self.MUSCLE_bin):
            raise ValueError("no such file '"+self.MUSCLE_bin+"'")
        if not os.path.isfile(input_forward_reads_file_path):
            raise ValueError("no such file '"+input_forward_reads_file_path+"'")
        elif not os.path.getsize(input_forward_reads_file_path) > 0:
            raise ValueError("empty file '"+input_forward_reads_file_path+"'")

        # set the output path
        timestamp = int((datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()*1000)
        output_dir = os.path.join(self.scratch,'output.'+str(timestamp))
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        output_aln_file_path = os.path.join(output_dir, params['input_name']+'-MSA.fasta');
        file_extension = ''

        muscle_cmd.append('-in')
        muscle_cmd.append(input_forward_reads_file_path)
        muscle_cmd.append('-out')
        muscle_cmd.append(output_aln_file_path)

        # options
        if 'maxiters' in params and params['maxiters'] != None:
            muscle_cmd.append('-maxiters')
            muscle_cmd.append(str(params['maxiters']))
        if 'maxhours' in params and params['maxhours'] != None:
            muscle_cmd.append('-maxhours')
            muscle_cmd.append(str(params['maxhours']))


        # Run MUSCLE, capture output as it happens
        #
        self.log(console, 'RUNNING MUSCLE:')
        self.log(console, '    '+' '.join(muscle_cmd))
#        report += "\n"+'running MUSCLE:'+"\n"
#        report += '    '+' '.join(muscle_cmd)+"\n"

        p = subprocess.Popen(muscle_cmd, \
                             cwd = self.scratch, \
                             stdout = subprocess.PIPE, \
                             stderr = subprocess.STDOUT, \
                             shell = False)

        while True:
            line = p.stdout.readline()
            if not line: break
            self.log(console, line.replace('\n', ''))

        p.stdout.close()
        p.wait()
        self.log(console, 'return code: ' + str(p.returncode))
        if p.returncode != 0:
            raise ValueError('Error running MUSCLE, return code: '+str(p.returncode) + 
                '\n\n'+ '\n'.join(console))


        # Parse the FASTA MSA output
        #
        self.log(console, 'PARSING MUSCLE MSA FASTA OUTPUT')
        if not os.path.isfile(output_aln_file_path):
            raise ValueError("failed to create MUSCLE output: "+output_aln_file_path)
        elif not os.path.getsize(output_aln_file_path) > 0:
            raise ValueError("created empty file for MUSCLE output: "+output_aln_file_path)
        output_aln_file_handle = open (output_aln_file_path, "r", 0)

        row_order = []
        alignment = {}
        alignment_length = None

        last_header = None
        header = None
        last_seq = ''
        leading_chars_pattern = re.compile("^\S+")
        for line in output_aln_file_handle:
            line = line.rstrip('\n')
            if line.startswith('>'):
                header = line[1:]

                if last_header != None:
                    last_id = leading_chars_pattern.findall(last_header)[0]
                    row_order.append(last_id)
                    #self.log(console,"ID: '"+last_id+"'\nALN: '"+last_seq+"'")  # DEBUG
                    report += last_id+"\t"+last_seq+"\n"
                    alignment[last_id] = last_seq
                    if alignment_length == None:
                        alignment_length = len(last_seq)
                    elif alignment_length != len(last_seq):
                        raise ValueError ("unequal alignment row for "+last_header+": '"+last_seq+"'")
                last_header = header
                last_seq = ''
            else:
                last_seq += line
        if last_header != None:
            last_id = leading_chars_pattern.findall(last_header)[0]
            row_order.append(last_id)
            #self.log(console,"ID: '"+last_id+"'\nALN: '"+last_seq+"'")  # DEBUG
            report += last_id+"\t"+last_seq+"\n"
            alignment[last_id] = last_seq
            if alignment_length == None:
                alignment_length = len(last_seq)
            elif alignment_length != len(last_seq):
                raise ValueError ("unequal alignment row for "+last_header+": '"+last_seq+"'")
        
        output_aln_file_handle.close()


        # load the method provenance from the context object
        #
        self.log(console,"SETTING PROVENANCE")  # DEBUG
        provenance = [{}]
        if 'provenance' in ctx:
            provenance = ctx['provenance']
        # add additional info to provenance here, in this case the input data object reference
        provenance[0]['input_ws_objects'] = []
        provenance[0]['input_ws_objects'].append(params['workspace_name']+'/'+params['input_name'])
        provenance[0]['service'] = 'kb_muscle'
        provenance[0]['method'] = 'MUSCLE_nuc'


        # Upload results
        #
        if len(invalid_msgs) == 0:
            self.log(console,"UPLOADING RESULTS")  # DEBUG

            MSA_name = params['output_name']
            MSA_description = params['desc']
            sequence_type = 'dna'
            ws_refs = None  # may add these later from FeatureSet
            kb_refs = None
            #alignment_length  # already have
            #row_order  # already have
            #alignment  # already have
            # NO trim_info
            # NO alignment_attributes
            # NO default_row_labels
            # NO parent_msa_ref

#            if input_type_name == 'FeatureSet':
#                features = featureSet['elements']
#                genome2Features = {}
#                for fId in row_order:
#                    genomeRef = features[fId][0]
#                    if genomeRef not in genome2Features:
#                        genome2Features[genomeRef] = []
#                    genome2Features[genomeRef].append(fId)
#
#                for genomeRef in genome2Features:
#                    genome = ws.get_objects([{'ref':genomeRef}])[0]['data']
#                    these_genomeFeatureIds = genome2Features[genomeRef]
#                    for feature in genome['features']:
#                        if feature['id'] in these_genomeFeatureIds:

            output_MSA = {
                      'name': MSA_name,
                      'description': MSA_description,
                      'sequence_type': sequence_type,
                      'alignment_length': alignment_length,
                      'row_order': row_order,
                      'alignment': alignment
                     }

            new_obj_info = ws.save_objects({
                            'workspace': params['workspace_name'],
                            'objects':[{
                                    'type': 'KBaseTrees.MSA',
                                    'data': output_MSA,
                                    'name': params['output_name'],
                                    'meta': {},
                                    'provenance': provenance
                                }]
                        })


        # build output report object
        #
        self.log(console,"BUILDING REPORT")  # DEBUG
        if len(invalid_msgs) == 0:
#            self.log(console,"sequences in many set: "+str(seq_total))
#            self.log(console,"sequences in hit set:  "+str(hit_total))
#            report += 'sequences in many set: '+str(seq_total)+"\n"
#            report += 'sequences in hit set:  '+str(hit_total)+"\n"

            reportObj = {
                'objects_created':[{'ref':params['workspace_name']+'/'+params['output_name'], 'description':'MUSCLE_nuc MSA'}],
                'text_message':report
                }
        else:
            report += "FAILURE:\n\n"+"\n".join(invalid_msgs)+"\n"
            reportObj = {
                'objects_created':[],
                'text_message':report
                }

        reportName = 'muscle_report_'+str(hex(uuid.getnode()))
        ws = workspaceService(self.workspaceURL, token=ctx['token'])
        report_obj_info = ws.save_objects({
#                'id':info[6],
                'workspace':params['workspace_name'],
                'objects':[
                    {
                        'type':'KBaseReport.Report',
                        'data':reportObj,
                        'name':reportName,
                        'meta':{},
                        'hidden':1,
                        'provenance':provenance
                    }
                ]
            })[0]


        self.log(console,"BUILDING RETURN OBJECT")
        returnVal = { 'report_name': reportName,
                      'report_ref': str(report_obj_info[6]) + '/' + str(report_obj_info[0]) + '/' + str(report_obj_info[4]),
                      }
        self.log(console,"MUSCLE_nuc DONE")

        #END MUSCLE_nuc

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method MUSCLE_nuc return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]

    def MUSCLE_prot(self, ctx, params):
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN MUSCLE_prot
        console = []
        invalid_msgs = []
        self.log(console,'Running MUSCLE_prot with params=')
        self.log(console, "\n"+pformat(params))
        report = ''
#        report = 'Running MUSCLE_prot with params='
#        report += "\n"+pformat(params)


        #### do some basic checks
        #
        if 'workspace_name' not in params:
            raise ValueError('workspace_name parameter is required')
        if 'input_name' not in params:
            raise ValueError('input_name parameter is required')
        if 'output_name' not in params:
            raise ValueError('output_name parameter is required')


        #### Get the input_name object
        ##
#        input_forward_reads_file_compression = None
#        sequencing_tech = 'N/A'
        try:
            ws = workspaceService(self.workspaceURL, token=ctx['token'])
            objects = ws.get_objects([{'ref': params['workspace_name']+'/'+params['input_name']}])
            data = objects[0]['data']
            info = objects[0]['info']
            input_type_name = info[2].split('.')[1].split('-')[0]

#            if input_type_name == 'SingleEndLibrary':
#                input_type_namespace = info[2].split('.')[0]
#                if input_type_namespace == 'KBaseAssembly':
#                    file_name = data['handle']['file_name']
#                elif input_type_namespace == 'KBaseFile':
#                    file_name = data['lib']['file']['file_name']
#                else:
#                    raise ValueError('bad data type namespace: '+input_type_namespace)
#                #self.log(console, 'INPUT_FILENAME: '+file_name)  # DEBUG
#                if file_name[-3:] == ".gz":
#                    input_forward_reads_file_compression = 'gz'
#                if 'sequencing_tech' in data:
#                    sequencing_tech = data['sequencing_tech']

        except Exception as e:
            traceback.format_exc()
            raise ValueError('Unable to fetch input_name object from workspace: ' + str(e))


        # Handle overloading (input_name can be SingleEndLibrary or FeatureSet)
        #
        """
        if input_type_name == 'SingleEndLibrary':

            # DEBUG
            #for k in data:
            #    self.log(console,"SingleEndLibrary ["+k+"]: "+str(data[k]))

            try:
                if 'lib' in data:
                    input_forward_reads = data['lib']['file']
                elif 'handle' in data:
                    input_forward_reads = data['handle']
                else:
                    self.log(console,"bad structure for 'input_forward_reads'")
                    raise ValueError("bad structure for 'input_forward_reads'")

                ### NOTE: this section is what could be replaced by the transform services
                input_forward_reads_file_path = os.path.join(self.scratch,input_forward_reads['file_name'])
                input_forward_reads_file_handle = open(input_forward_reads_file_path, 'w', 0)
                self.log(console, 'downloading reads file: '+str(input_forward_reads_file_path))
                headers = {'Authorization': 'OAuth '+ctx['token']}
                r = requests.get(input_forward_reads['url']+'/node/'+input_forward_reads['id']+'?download', stream=True, headers=headers)
                for chunk in r.iter_content(1024):
                    input_forward_reads_file_handle.write(chunk)
                input_forward_reads_file_handle.close();
                self.log(console, 'done')
                ### END NOTE


                # remove carriage returns
                new_file_path = input_forward_reads_file_path+"-CRfree"
                new_file_handle = open(new_file_path, 'w', 0)
                input_forward_reads_file_handle = open(input_forward_reads_file_path, 'r', 0)
                for line in input_forward_reads_file_handle:
                    line = re.sub("\r","",line)
                    new_file_handle.write(line)
                input_forward_reads_file_handle.close();
                new_file_handle.close()
                input_forward_reads_file_path = new_file_path

                # convert FASTQ to FASTA (if necessary)
                new_file_path = input_forward_reads_file_path+".fna"
                new_file_handle = open(new_file_path, 'w', 0)
                if input_forward_reads_file_compression == 'gz':
                    input_forward_reads_file_handle = gzip.open(input_forward_reads_file_path, 'r', 0)
                else:
                    input_forward_reads_file_handle = open(input_forward_reads_file_path, 'r', 0)
                header = None
                last_header = None
                last_seq_buf = None
                last_line_was_header = False
                was_fastq = False
                for line in input_forward_reads_file_handle:
                    if line.startswith('>'):
                        break
                    elif line.startswith('@'):
                        was_fastq = True
                        header = line[1:]
                        if last_header != None:
                            new_file_handle.write('>'+last_header)
                            new_file_handle.write(last_seq_buf)
                        last_seq_buf = None
                        last_header = header
                        last_line_was_header = True
                    elif last_line_was_header:
                        last_seq_buf = line
                        last_line_was_header = False
                    else:
                        continue
                if last_header != None:
                    new_file_handle.write('>'+last_header)
                    new_file_handle.write(last_seq_buf)

                new_file_handle.close()
                input_forward_reads_file_handle.close()
                if was_fastq:
                    input_forward_reads_file_path = new_file_path

            except Exception as e:
                print(traceback.format_exc())
                raise ValueError('Unable to download single-end read library files: ' + str(e))
        """

        # FeatureSet
        #
#        elif input_type_name == 'FeatureSet':
        if input_type_name == 'FeatureSet':
            # retrieve sequences for features
            input_featureSet = data

            genome2Features = {}
            features = input_featureSet['elements']
            for fId in features.keys():
                genomeRef = features[fId][0]
                if genomeRef not in genome2Features:
                    genome2Features[genomeRef] = []
                genome2Features[genomeRef].append(fId)

            # export features to FASTA file
            input_forward_reads_file_path = os.path.join(self.scratch, params['input_name']+".fasta")
            self.log(console, 'writing fasta file: '+input_forward_reads_file_path)
            records = []
            proteins_found = 0
            for genomeRef in genome2Features:
                genome = ws.get_objects([{'ref':genomeRef}])[0]['data']
                these_genomeFeatureIds = genome2Features[genomeRef]
                for feature in genome['features']:
                    if feature['id'] in these_genomeFeatureIds:
                        #self.log(console,"kbase_id: '"+feature['id']+"'")  # DEBUG
                        #record = SeqRecord(Seq(feature['dna_sequence']), id=feature['id'], description=genome['id'])
                        if feature['type'] != 'CDS':
                            self.log(console,"attempt to include non-CDS Feature "+feature['id'])
                            self.log(invalid_msgs,"attempt to include non-CDS Feature "+feature['id'])
                            continue
                        elif 'protein_translation' not in feature or feature['protein_translation'] == None:
                            self.log(console,"bad CDS Feature "+feature['id']+": no protein_translation found")
                            self.log(invalid_msgs,"bad CDS Feature "+feature['id']+": no protein_translation found")
                            continue
                        else:
                            record = SeqRecord(Seq(feature['protein_translation']), id=feature['id'], description=genome['id'])
                            proteins_found += 1
                            records.append(record)

            if proteins_found < 2:
                self.log(invalid_msgs,"Less than 2 protein Features (CDS) found.  exiting...")
            else:
                SeqIO.write(records, input_forward_reads_file_path, "fasta")

        # Missing proper input_input_type
        #
        else:
            raise ValueError('Cannot yet handle input_name type of: '+type_name)


        ### Construct the command
        #
        #  e.g. muscle -in <fasta_in> -out <fasta_out> -maxiters <n> -haxours <h>
        #
        if len(invalid_msgs) == 0:
            muscle_cmd = [self.MUSCLE_bin]

            # check for necessary files
            if not os.path.isfile(self.MUSCLE_bin):
                raise ValueError("no such file '"+self.MUSCLE_bin+"'")
            if not os.path.isfile(input_forward_reads_file_path):
                raise ValueError("no such file '"+input_forward_reads_file_path+"'")
            elif not os.path.getsize(input_forward_reads_file_path) > 0:
                raise ValueError("empty file '"+input_forward_reads_file_path+"'.  May have not provided any protein coding genes?")

            # set the output path
            timestamp = int((datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()*1000)
            output_dir = os.path.join(self.scratch,'output.'+str(timestamp))
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
            output_aln_file_path = os.path.join(output_dir, params['input_name']+'-MSA.fasta');
            file_extension = ''

            muscle_cmd.append('-in')
            muscle_cmd.append(input_forward_reads_file_path)
            muscle_cmd.append('-out')
            muscle_cmd.append(output_aln_file_path)

            # options
            if 'maxiters' in params and params['maxiters'] != None:
                muscle_cmd.append('-maxiters')
                muscle_cmd.append(str(params['maxiters']))
            if 'maxhours' in params and params['maxhours'] != None:
                muscle_cmd.append('-maxhours')
                muscle_cmd.append(str(params['maxhours']))


            # Run MUSCLE, capture output as it happens
            #
            self.log(console, 'RUNNING MUSCLE:')
            self.log(console, '    '+' '.join(muscle_cmd))
#            report += "\n"+'running MUSCLE:'+"\n"
#            report += '    '+' '.join(muscle_cmd)+"\n"

            p = subprocess.Popen(muscle_cmd, \
                             cwd = self.scratch, \
                             stdout = subprocess.PIPE, \
                             stderr = subprocess.STDOUT, \
                             shell = False)

            while True:
                line = p.stdout.readline()
                if not line: break
                self.log(console, line.replace('\n', ''))

            p.stdout.close()
            p.wait()
            self.log(console, 'return code: ' + str(p.returncode))
            if p.returncode != 0:
                raise ValueError('Error running MUSCLE, return code: '+str(p.returncode) + 
                                 '\n\n'+ '\n'.join(console))


            # Parse the FASTA MSA output
            #
            self.log(console, 'PARSING MUSCLE MSA FASTA OUTPUT')
            if not os.path.isfile(output_aln_file_path):
                raise ValueError("failed to create MUSCLE output: "+output_aln_file_path)
            elif not os.path.getsize(output_aln_file_path) > 0:
                raise ValueError("created empty file for MUSCLE output: "+output_aln_file_path)
            output_aln_file_handle = open (output_aln_file_path, "r", 0)

            row_order = []
            alignment = {}
            alignment_length = None

            last_header = None
            header = None
            last_seq = ''
            leading_chars_pattern = re.compile("^\S+")
            for line in output_aln_file_handle:
                line = line.rstrip('\n')
                if line.startswith('>'):
                    header = line[1:]

                    if last_header != None:
                        last_id = leading_chars_pattern.findall(last_header)[0]
                        row_order.append(last_id)
                        #self.log(console,"ID: '"+last_id+"'\nALN: '"+last_seq+"'")  # DEBUG
                        report += last_id+"\t"+last_seq+"\n"
                        alignment[last_id] = last_seq
                        if alignment_length == None:
                            alignment_length = len(last_seq)
                        elif alignment_length != len(last_seq):
                            raise ValueError ("unequal alignment row for "+last_header+": '"+last_seq+"'")
                    last_header = header
                    last_seq = ''
                else:
                    last_seq += line
            if last_header != None:
                last_id = leading_chars_pattern.findall(last_header)[0]
                row_order.append(last_id)
                #self.log(console,"ID: '"+last_id+"'\nALN: '"+last_seq+"'")  # DEBUG
                report += last_id+"\t"+last_seq+"\n"
                alignment[last_id] = last_seq
                if alignment_length == None:
                    alignment_length = len(last_seq)
                elif alignment_length != len(last_seq):
                    raise ValueError ("unequal alignment row for "+last_header+": '"+last_seq+"'")
        
            output_aln_file_handle.close()


        # load the method provenance from the context object
        #
        self.log(console,"SETTING PROVENANCE")  # DEBUG
        provenance = [{}]
        if 'provenance' in ctx:
            provenance = ctx['provenance']
        # add additional info to provenance here, in this case the input data object reference
        provenance[0]['input_ws_objects'] = []
        provenance[0]['input_ws_objects'].append(params['workspace_name']+'/'+params['input_name'])
        provenance[0]['service'] = 'kb_muscle'
        provenance[0]['method'] = 'MUSCLE_prot'


        # Upload results
        #
        if len(invalid_msgs) == 0:
            self.log(console,"UPLOADING RESULTS")  # DEBUG

            MSA_name = params['output_name']
            MSA_description = params['desc']
            sequence_type = 'protein'
            ws_refs = None  # may add these later from FeatureSet
            kb_refs = None
#            alignment_length  # already have
#            row_order  # already have
#            alignment  # already have
#            NO trim_info
#            NO alignment_attributes
#            NO default_row_labels
#            NO parent_msa_ref

#            if input_type_name == 'FeatureSet':
#                features = featureSet['elements']
#                genome2Features = {}
#                for fId in row_order:
#                    genomeRef = features[fId][0]
#                    if genomeRef not in genome2Features:
#                        genome2Features[genomeRef] = []
#                    genome2Features[genomeRef].append(fId)
#
#                for genomeRef in genome2Features:
#                    genome = ws.get_objects([{'ref':genomeRef}])[0]['data']
#                    these_genomeFeatureIds = genome2Features[genomeRef]
#                    for feature in genome['features']:
#                        if feature['id'] in these_genomeFeatureIds:

            output_MSA = {
                      'name': MSA_name,
                      'description': MSA_description,
                      'sequence_type': sequence_type,
                      'alignment_length': alignment_length,
                      'row_order': row_order,
                      'alignment': alignment
                     }

            new_obj_info = ws.save_objects({
                            'workspace': params['workspace_name'],
                            'objects':[{
                                    'type': 'KBaseTrees.MSA',
                                    'data': output_MSA,
                                    'name': params['output_name'],
                                    'meta': {},
                                    'provenance': provenance
                                }]
                        })


        # build output report object
        #
        self.log(console,"BUILDING REPORT")  # DEBUG
        if len(invalid_msgs) == 0:
#            self.log(console,"sequences in many set: "+str(seq_total))
#            self.log(console,"sequences in hit set:  "+str(hit_total))
#            report += 'sequences in many set: '+str(seq_total)+"\n"
#            report += 'sequences in hit set:  '+str(hit_total)+"\n"
            reportObj = {
                'objects_created':[{'ref':params['workspace_name']+'/'+params['output_name'], 'description':'MUSCLE_prot MSA'}],
                'text_message':report
                }
        else:
            report += "FAILURE:\n\n"+"\n".join(invalid_msgs)+"\n"
            reportObj = {
                'objects_created':[],
                'text_message':report
                }

        reportName = 'muscle_report_'+str(hex(uuid.getnode()))
        ws = workspaceService(self.workspaceURL, token=ctx['token'])
        report_obj_info = ws.save_objects({
#                'id':info[6],
                'workspace':params['workspace_name'],
                'objects':[
                    {
                        'type':'KBaseReport.Report',
                        'data':reportObj,
                        'name':reportName,
                        'meta':{},
                        'hidden':1,
                        'provenance':provenance
                    }
                ]
            })[0]


        self.log(console,"BUILDING RETURN OBJECT")
        returnVal = { 'report_name': reportName,
                      'report_ref': str(report_obj_info[6]) + '/' + str(report_obj_info[0]) + '/' + str(report_obj_info[4]),
                      }
        self.log(console,"MUSCLE_prot DONE")

        #END MUSCLE_prot

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method MUSCLE_prot return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]
