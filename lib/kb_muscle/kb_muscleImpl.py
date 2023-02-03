 -*- coding: utf-8 -*-
#BEGIN_HEADER
import gzip
import os
import re
import subprocess
import sys
import traceback
import uuid
import json
from datetime import datetime
from pprint import pformat

import requests
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from requests_toolbelt import MultipartEncoder

from installed_clients.AbstractHandleClient import AbstractHandle as HandleService
from installed_clients.DataFileUtilClient import DataFileUtil as DFUClient
from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.WorkspaceClient import Workspace as workspaceService
from installed_clients.SetAPIServiceClient import SetAPI
from installed_clients.AssemblyUtilClient import AssemblyUtil
from installed_clients.kb_ObjectUtilitiesClient import kb_ObjectUtilities
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

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "1.1.1"
    GIT_URL = "https://github.com/kbaseapps/kb_muscle"
    GIT_COMMIT_HASH = "d25d4d112be6c3fce5d879734dabdf5cc524ea2f"

    #BEGIN_CLASS_HEADER
    workspaceURL     = None
    shockURL         = None
    handleURL        = None
    serviceWizardURL = None
    callbackURL      = None
    scratch          = None

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


    # Translation
    def TranslateNucToProtSeq(self, ctx, params):
        if 'nuc_seq' not in params or params['nuc_seq'] == None or params['nuc_seq'] == '':
            raise ValueError('Method TranslateNucToProtSeq() requires nuc_seq parameter')
        if 'genetic_code' not in params or params['genetic_code'] == None or params['genetic_code'] == '':
            params['genetic_code'] = '11'

        if params['genetic_code'] != '11':
            raise ValueError('Method TranslateNucToProtSeq() only knows genetic code 11')
        
        nuc_seq = params['nuc_seq'].upper()
        prot_seq = ''

        genetic_code = params['genetic_code']
        genetic_code_table = dict()
        genetic_code_table['11'] = {
            'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
            'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
            'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
            'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
            'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
            'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
            'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
            'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
            'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
            'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
            'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
            'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
            'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
            'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
            'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
            'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'
            }
        if genetic_code not in genetic_code_table:
            raise ValueError ("genetic code '"+str(genetic_code)+"' not configured in genetic_code_table")

        prot_seq = ''.join([genetic_code_table[genetic_code].get(nuc_seq[3*i:3*i+3],'X') for i in range(len(nuc_seq)//3)])
        if prot_seq.endswith('_'):
            prot_seq = prot_seq.rstrip('_')

        return prot_seq

    # AMA_METHODS
    #def _get_ama_features_as_json (self, features_handle_ref, gff_handle_ref, protein_handle_ref):
    def _get_ama_features_as_json (self, features_handle_ref):
        this_id = str(uuid.uuid4())
        this_scratch_dir = os.path.join (self.scratch, this_id)
        json_features_file_path = os.path.join (this_scratch_dir, 'features.json')
        #gff_file_path = os.path.join (this_scratch_dir, 'genes.gff')
        #protein_file_path = os.path.join (this_scratch_dir, 'protein.fasta')

        try:
            dfu = DFUClient (self.callbackURL)
        except Exception as e:
            raise ValueError('Unable to connect to DFU: ' + str(e))

        try:
            dfu.shock_to_file({'handle_id': features_handle_ref,
                               'file_path': json_features_file_path+'.gz',
                               'unpack': 'uncompress'
                           })
        except Exception as e:
            raise ValueError('Unable to fetch AnnotatedMetagenomeAssembly features from SHOCK: ' + str(e))
        """
        try:
            dfu.shock_to_file({'handle_id': gff_handle_ref,
                               'file_path': gff_file_path+'.gz',
                               'unpack': 'uncompress'
                           })
        except Exception as e:
            raise ValueError('Unable to fetch AnnotatedMetagenomeAssembly gffs from SHOCK: ' + str(e))

        try:
            dfu.shock_to_file({'handle_id': protein_handle_ref,
                               'file_path': protein_file_path+'.gz',
                               'unpack': 'uncompress'
                           })
        except Exception as e:
            raise ValueError('Unable to fetch AnnotatedMetagenomeAssembly protein FASTA from SHOCK: ' + str(e))
        """

        # DEBUG
        """
        print ("SCRATCH CONTENTS")
        sys.stdout.flush()
        for this_file in os.listdir (this_scratch_dir):
            print ("\t"+this_file)
            sys.stdout.flush()

        buf = []
        #with open(json_features_file_path, 'r') as f:
        with open(protein_file_path, 'r') as f:
            for line in f.readlines():
                buf.append (line)
            #features_json = json.load(f)

        print ("FEATURES_JSON:\n"+"\n".join(buf))
        sys.stdout.flush()
        """
        
        with open(json_features_file_path, 'r') as f:
            features_json = json.load(f)


        os.remove(json_features_file_path+'.gz')
        os.remove(json_features_file_path)
        #os.remove(gff_file_path+'.gz')
        #os.remove(gff_file_path)
        #os.remove(protein_file_path+'.gz')
        #os.remove(protein_file_path)

        return features_json


    def _get_features_from_AnnotatedMetagenomeAssembly(self, ctx, ama_ref):

        # get ama object
        try:
            ws = workspaceService(self.workspaceURL, token=ctx['token'])
            ama_object = ws.get_objects2({'objects':[{'ref':ama_ref}]})['data'][0]
            ama_object_data = ama_object['data']
            ama_object_info = ama_object['info']
        except Exception as e:
            raise ValueError('Unable to fetch AnnotatedMetagenomeAssembly object from workspace: ' + str(e))
        #to get the full stack trace: traceback.format_exc()

        # get features from json
        features_handle_ref = ama_object_data['features_handle_ref']
        #gff_handle_ref = ama_object_data['gff_handle_ref']
        #protein_handle_ref = ama_object_data['protein_handle_ref']
        #features_json = self._get_ama_features_as_json (features_handle_ref, gff_handle_ref, protein_handle_ref)
        features_json = self._get_ama_features_as_json (features_handle_ref)

        return features_json

    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.workspaceURL = config['workspace-url']
        self.shockURL = config['shock-url']
        self.handleURL = config['handle-service-url']
        self.serviceWizardURL = config['srv-wiz-url']

        self.callbackURL = os.environ.get('SDK_CALLBACK_URL')
        if self.callbackURL == None:
            raise ValueError ("SDK_CALLBACK_URL not set in environment")

        self.scratch = os.path.abspath(config['scratch'])
        # HACK!! temporary hack for issue where megahit fails on mac because of silent named pipe error
        #self.host_scratch = self.scratch
        #self.scratch = os.path.join('/kb','module','local_scratch')
        # end hack
        if not os.path.exists(self.scratch):
            os.makedirs(self.scratch)

        #END_CONSTRUCTOR
        pass


    def MUSCLE_nuc(self, ctx, params):
        """
        Methods for MSA building of either DNA or PROTEIN sequences
        **
        **    overloading as follows:
        **        input_ref: SingleEndLibrary (just MUSCLE_nuc), FeatureSet (both)
        **        output_name: MSA
        :param params: instance of type "MUSCLE_Params" (MUSCLE Input Params
           ** ** MUSCLE_prot(): input_ref must be FeatureSet ** MUSCLE_nuc():
           input_ref must be FeatureSet, SingleEndLibrary, or AssemblySet) ->
           structure: parameter "workspace_name" of type "workspace_name" (**
           The workspace object refs are of form: ** **    objects =
           ws.get_objects([{'ref':
           params['workspace_id']+'/'+params['obj_name']}]) ** ** "ref" means
           the entire name combining the workspace id and the object name **
           "id" is a numerical identifier of the workspace or object, and
           should just be used for workspace ** "name" is a string identifier
           of a workspace or object.  This is received from Narrative.),
           parameter "desc" of String, parameter "input_ref" of type
           "data_obj_ref", parameter "output_name" of type "data_obj_name",
           parameter "genome_disp_name_config" of String, parameter
           "maxiters" of Long, parameter "maxhours" of Double
        :returns: instance of type "MUSCLE_Output" (MUSCLE Output) ->
           structure: parameter "report_name" of type "data_obj_name",
           parameter "report_ref" of type "data_obj_ref"
        """
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
        [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I,
         WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = range(11)  # object_info tuple
        
        row_labels = {}


        #### do some basic checks
        #
        if 'workspace_name' not in params:
            raise ValueError('workspace_name parameter is required')
        if 'input_ref' not in params:
            raise ValueError('input_ref parameter is required')
        if 'output_name' not in params:
            raise ValueError('output_name parameter is required')


        #### Get the input_ref object
        ##
        input_forward_reads_file_compression = None
        sequencing_tech = 'N/A'
        try:
            ws = workspaceService(self.workspaceURL, token=ctx['token'])
            objects = ws.get_objects([{'ref': params['input_ref']}])
            data = objects[0]['data']
            info = objects[0]['info']
            input_name = info[1]
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
            raise ValueError('Unable to fetch input_ref object from workspace: ' + str(e))


        # Handle overloading (input_ref can be SingleEndLibrary or FeatureSet)
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
                input_forward_reads_file_handle = open(input_forward_reads_file_path, 'w')
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
                new_file_handle = open(new_file_path, 'w')
                input_forward_reads_file_handle = open(input_forward_reads_file_path, 'r')
                for line in input_forward_reads_file_handle:
                    line = re.sub("\r","",line)
                    new_file_handle.write(line)
                input_forward_reads_file_handle.close();
                new_file_handle.close()
                input_forward_reads_file_path = new_file_path

                # convert FASTQ to FASTA (if necessary)
                new_file_path = input_forward_reads_file_path+".fna"
                new_file_handle = open(new_file_path, 'w')
                if input_forward_reads_file_compression == 'gz':
                    input_forward_reads_file_handle = gzip.open(input_forward_reads_file_path, 'r')
                else:
                    input_forward_reads_file_handle = open(input_forward_reads_file_path, 'r')
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
            
            genome_id_feature_id_delim = '.f:'

            # retrieve sequences for features
            input_featureSet = data

            genomeObjName = {}
            genomeObjVer  = {}
            genomeSciName = {}
            genome2Features = {}
            new_id = {}
            featureSet_elements = input_featureSet['elements']
            if 'element_ordering' in input_featureSet and input_featureSet['element_ordering']:
                feature_order = input_featureSet['element_ordering']
            else:
                feature_order = sorted(featureSet_elements.keys())
            for fId in feature_order:
                genomeRef = featureSet_elements[fId][0]
                if genomeRef not in genome2Features:
                    genome2Features[genomeRef] = []
                    new_id[genomeRef] = {}
                if genome_id_feature_id_delim in fId:
                    [genome_id, feature_id] = fId.split(genome_id_feature_id_delim)
                else:
                    feature_id = fId
                genome2Features[genomeRef].append(feature_id)
                this_id = genomeRef + genome_id_feature_id_delim + feature_id
                new_id[genomeRef][fId] = this_id

            # export features to FASTA file
            input_forward_reads_file_path = os.path.join(self.scratch, input_name+".fasta")
            self.log(console, 'writing fasta file: '+input_forward_reads_file_path)
            records_by_fid = dict()
            for genomeRef in genome2Features:
                genome_obj = ws.get_objects([{'ref':genomeRef}])[0]
                genome_type = re.sub('-[0-9]+\.[0-9]+$', "", genome_obj['info'][TYPE_I])
                genomeObjName[genomeRef] = genome_obj['info'][NAME_I]
                genomeObjVer[genomeRef]  = genome_obj['info'][VERSION_I]
                these_genomeFeatureIds = genome2Features[genomeRef]

                # Genome
                if genome_type == 'KBaseGenomes.Genome':
                    genome = genome_obj['data']
                    genomeSciName[genomeRef] = genome['scientific_name']
                    for feature in genome['features']:
                        if feature['id'] in these_genomeFeatureIds:
                            #self.log(console,"kbase_id: '"+feature['id']+"'")  # DEBUG
                            this_id = genomeRef + genome_id_feature_id_delim + feature['id']
                            short_feature_id = re.sub("^.*\.([^\.]+)\.([^\.]+)$", r"\1.\2", feature['id'])
                            genome_disp_name = ''
                            if 'obj_name' in params.get('genome_disp_name_config'):
                                genome_disp_name += genomeObjName[genomeRef]
                                if 'ver' in params.get('genome_disp_name_config'):
                                    genome_disp_name += '.v'+str(genomeObjVer[genomeRef])

                                if genome_type == "KBaseGenomes.Genome" and \
                                   'sci_name' in params.get('genome_disp_name_config'):
                                    genome_disp_name += ': '+genomeSciName[genomeRef]
                            else:
                                genome_disp_name = genomeObjName[genomeRef]
                
                            row_labels[this_id] = genome_disp_name+' - '+short_feature_id

                            #record = SeqRecord(Seq(feature['dna_sequence']), id=feature['id'], description=genome['id'])
                            record = SeqRecord(Seq(feature['dna_sequence']), id=this_id, description=genome['id'])
                            records_by_fid[this_id] = record

                            # DEBUG
                            print("THIS_ID: '"+this_id+"'")
                            print("RECORD: '"+record+"'")
                            
                # AnnotatedMetagenomeAssembly
                elif genome_type == 'KBaseMetagenomes.AnnotatedMetagenomeAssembly':
                    ama_features = self._get_features_from_AnnotatedMetagenomeAssembly (ctx, genomeRef)
                    for feature in ama_features:
                        if feature['id'] in these_genomeFeatureIds:
                            if not feature.get('dna_sequence'):
                                raise ValueError("bad feature "+feature['id']+": No dna_sequence field.")
                            this_id = genomeRef + genome_id_feature_id_delim + feature['id']
                            short_feature_id = re.sub("^.*\.([^\.]+)\.([^\.]+)$", r"\1.\2", feature['id'])
                            genome_disp_name = genomeObjName[genomeRef]
                            row_labels[this_id] = genome_disp_name+' - '+short_feature_id

                            record = SeqRecord(Seq(feature['dna_sequence']), id=this_id, description=genomeObjName[genomeRef])
                            records_by_fid[this_id] = record

                else:
                    raise ValueError ("unable to handle feature from object type: "+genome_type)


            records = []
            for fId in feature_order:
                genomeRef = featureSet_elements[fId][0]

                # DEBUG
                print("GENOME_REF: '"+genomeRef+"'")
                print("FID: '"+fId+"'")
                print("NEW_ID: '"+new_id[genomeRef][fId]+"'")
                            
                records.append(records_by_fid[new_id[genomeRef][fId]])
            SeqIO.write(records, input_forward_reads_file_path, "fasta")

        # Missing proper input_input_type
        #
        else:
            raise ValueError('Cannot yet handle input_ref type of: '+input_type_name)

        """
        # AssemblySet
        #
        elif input_type_name == 'AssemblySet':
            try:
                SetAPI_Client = SetAPI(self.serviceWizardURL, token=ctx['token'])
            except Exception as e:
                raise ValueError ("unable to instantiate SetAPI Client")
            try:
                auClient = AssemblyUtil(self.callbackURL, token=ctx['token'])
            except Exception as e:
                raise ValueError ("unable to instantiate AssemblyUtil Client")

            # HERE
        """
        

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
        output_aln_file_path = os.path.join(output_dir, params['output_name']+'-MSA.fasta');
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
            line = p.stdout.readline().decode()
            if not line: break
            self.log(console, line.replace('\n', ''))

        p.stdout.close()
        p.wait()
        self.log(console, 'return code: ' + str(p.returncode))
        if p.returncode != 0:
            raise ValueError('Error running MUSCLE, return code: '+str(p.returncode) + 
                '\n\n'+ '\n'.join(console))


        # Parse the FASTA MSA output and replace id for txt upload
        #
        self.log(console, 'PARSING MUSCLE MSA FASTA OUTPUT')
        if not os.path.isfile(output_aln_file_path):
            raise ValueError("failed to create MUSCLE output: "+output_aln_file_path)
        elif not os.path.getsize(output_aln_file_path) > 0:
            raise ValueError("created empty file for MUSCLE output: "+output_aln_file_path)
        output_aln_file_handle = open (output_aln_file_path, 'r')

        output_fasta_buf = []
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
                if row_labels:
                    this_id = leading_chars_pattern.findall(header)[0]
                    this_row_label = re.sub ('\s', '_', row_labels[this_id])
                    output_fasta_buf.append('>'+this_row_label)
                else:
                    output_fasta_buf.append(line)

                if last_header != None:
                    last_id = leading_chars_pattern.findall(last_header)[0]
                    row_order.append(last_id)
                    #self.log(console,"ID: '"+last_id+"'\nALN: '"+last_seq+"'")  # DEBUG
                    #report += last_id+"\t"+last_seq+"\n"
                    alignment[last_id] = last_seq
                    if alignment_length == None:
                        alignment_length = len(last_seq)
                    elif alignment_length != len(last_seq):
                        raise ValueError ("unequal alignment row for "+last_header+": '"+last_seq+"'")
                last_header = header
                last_seq = ''
            else:
                last_seq += line
                output_fasta_buf.append(line)
        if last_header != None:
            last_id = leading_chars_pattern.findall(last_header)[0]
            row_order.append(last_id)
            #self.log(console,"ID: '"+last_id+"'\nALN: '"+last_seq+"'")  # DEBUG
            #report += last_id+"\t"+last_seq+"\n"
            alignment[last_id] = last_seq
            if alignment_length == None:
                alignment_length = len(last_seq)
            elif alignment_length != len(last_seq):
                raise ValueError ("unequal alignment row for "+last_header+": '"+last_seq+"'")
        
        output_aln_file_handle.close()

        # write remapped ids
        with open(output_aln_file_path, 'w') as output_aln_file_handle:
            output_aln_file_handle.write("\n".join(output_fasta_buf)+"\n")


        # load the method provenance from the context object
        #
        self.log(console,"SETTING PROVENANCE")  # DEBUG
        provenance = [{}]
        if 'provenance' in ctx:
            provenance = ctx['provenance']
        # add additional info to provenance here, in this case the input data object reference
        provenance[0]['input_ws_objects'] = []
        provenance[0]['input_ws_objects'].append(params['input_ref'])
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
            if row_labels:
                output_MSA['default_row_labels'] = row_labels

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


            # create CLW formatted output file
            max_row_width = 60
            id_aln_gap_width = 1
            gap_chars = ''
            for sp_i in range(id_aln_gap_width):
                gap_chars += ' '
# DNA
            strong_groups = { 'AG': True,
                              'CTU': True
                              }
            weak_groups = None
# PROTEINS
#        strong_groups = { 'AST':  True,
#                          'EKNQ': True,
#                          'HKNQ': True,
#                          'DENQ': True,
#                          'HKQR': True,
#                          'ILMV': True,
#                          'FILM': True,
#                          'HY':   True,
#                          'FWY':  True
#                          }
#        weak_groups = { 'ACS':    True,
#                        'ATV':    True,
#                        'AGS':    True,
#                        'KNST':   True,
#                        'APST':   True,
#                        'DGNS':   True,
#                        'DEKNQS': True,
#                        'DEHKNQ': True,
#                        'EHKNQR': True,
#                        'FILMV':  True,
#                        'FHY':    True
#                        }

            clw_buf = []
            clw_buf.append ('CLUSTALW format of MUSCLE alignment '+MSA_name+': '+MSA_description)
            clw_buf.append ('')

            long_id_len = 0
            aln_pos_by_id = dict()
            for row_id in row_order:
                aln_pos_by_id[row_id] = 0
                if row_labels:
                    row_id_disp = row_labels[row_id]
                else:
                    row_id_disp = row_id
                if long_id_len < len(row_id_disp):
                    long_id_len = len(row_id_disp)

            full_row_cnt = alignment_length // max_row_width
            if alignment_length % max_row_width == 0:
                full_row_cnt -= 1
            for chunk_i in range (full_row_cnt + 1):
                for row_id in row_order:
                    if row_labels:
                        row_id_disp = re.sub('\s', '_', row_labels[row_id])
                    else:
                        row_id_disp = row_id
                    for sp_i in range (long_id_len-len(row_id_disp)):
                        row_id_disp += ' '

                    aln_chunk_upper_bound = (chunk_i+1)*max_row_width
                    if aln_chunk_upper_bound > alignment_length:
                        aln_chunk_upper_bound = alignment_length
                    aln_chunk = alignment[row_id][chunk_i*max_row_width:aln_chunk_upper_bound]
                    for c in aln_chunk:
                        if c != '-':
                            aln_pos_by_id[row_id] += 1

                    clw_buf.append (row_id_disp+gap_chars+aln_chunk+' '+str(aln_pos_by_id[row_id]))

                # conservation line
                cons_line = ''
                for pos_i in range(chunk_i*max_row_width, aln_chunk_upper_bound):
                    col_chars = dict()
                    seq_cnt = 0
                    for row_id in row_order:
                        char = alignment[row_id][pos_i]
                        if char != '-':
                            seq_cnt += 1
                            col_chars[char] = True
                    if seq_cnt <= 1:
                        cons_char = ' '
                    elif len(col_chars.keys()) == 1:
                        cons_char = '*'
                    else:
                        strong = False
                        for strong_group in strong_groups.keys():
                            this_strong_group = True
                            for seen_char in col_chars.keys():
                                if seen_char not in strong_group:
                                    this_strong_group = False
                                    break
                            if this_strong_group:
                                strong = True
                                break
                        if not strong:
                            weak = False
                            if weak_groups is not None:
                                for weak_group in weak_groups.keys():
                                    this_weak_group = True
                                    for seen_char in col_chars.keys():
                                        if seen_char not in weak_group:
                                            this_strong_group = False
                                            break
                                    if this_weak_group:
                                        weak = True
                        if strong:
                            cons_char = ':'
                        elif weak:
                            cons_char = '.'
                        else:
                            cons_char = ' '
                    cons_line += cons_char

                lead_space = ''
                for sp_i in range(long_id_len):
                    lead_space += ' '
                lead_space += gap_chars

                clw_buf.append(lead_space+cons_line)
                clw_buf.append('')

            # write clw to file
            clw_buf_str = "\n".join(clw_buf)+"\n"
            output_clw_file_path = os.path.join(output_dir, input_name+'-MSA.clw')
            with open (output_clw_file_path, 'w') as output_clw_file_handle:
                output_clw_file_handle.write(clw_buf_str)


            # upload MUSCLE FASTA output to SHOCK for file_links
            dfu = DFUClient(self.callbackURL)
            try:
                output_upload_ret = dfu.file_to_shock({'file_path': output_aln_file_path,
# DEBUG
#                                                      'make_handle': 0,
#                                                      'pack': 'zip'})
                                                       'make_handle': 0})
            except:
                raise ValueError ('error loading aln_out file to shock')

            # upload MUSCLE CLW output to SHOCK for file_links
            try:
                output_clw_upload_ret = dfu.file_to_shock({'file_path': output_clw_file_path,
# DEBUG
#                                                      'make_handle': 0,
#                                                      'pack': 'zip'})
                                                           'make_handle': 0})
            except:
                raise ValueError ('error loading clw_out file to shock')


            # make HTML reports
            #
            # HERE



            # build output report object
            #
            self.log(console,"BUILDING REPORT")  # DEBUG

            reportName = 'muscle_report_'+str(uuid.uuid4())
            reportObj = {
                'objects_created':[{'ref':params['workspace_name']+'/'+params['output_name'],
                                    'description':'MUSCLE_nuc MSA'}],
                #'message': '',
                'message': clw_buf_str,
                'file_links': [],
                'workspace_name': params['workspace_name'],
                'report_object_name': reportName
                }
            reportObj['file_links'] = [{'shock_id': output_upload_ret['shock_id'],
                                        'name': params['output_name']+'-MUSCLE_nuc.FASTA',
                                        'label': 'MUSCLE_nuc FASTA'
                                        },
                                       {'shock_id': output_clw_upload_ret['shock_id'],
                                        'name': params['output_name']+'-MUSCLE_nuc.CLW',
                                        'label': 'MUSCLE_nuc CLUSTALW'
                                        }]

            # save report object
            #
            SERVICE_VER = 'release'
            reportClient = KBaseReport(self.callbackURL, token=ctx['token'], service_ver=SERVICE_VER)
            #report_info = report.create({'report':reportObj, 'workspace_name':params['workspace_name']})
            report_info = reportClient.create_extended_report(reportObj)                                       

        else:  # len(invalid_msgs) > 0
            reportName = 'muscle_report_'+str(uuid.uuid4())
            report += "FAILURE:\n\n"+"\n".join(invalid_msgs)+"\n"
            reportObj = {
                'objects_created':[],
                'text_message':report
                }

            ws = workspaceService(self.workspaceURL, token=ctx['token'])
            report_obj_info = ws.save_objects({
                    #'id':info[6],
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

            report_info = dict()
            report_info['name'] = report_obj_info[1]
            report_info['ref'] = str(report_obj_info[6])+'/'+str(report_obj_info[0])+'/'+str(report_obj_info[4])


        # done
        returnVal = { 'report_name': report_info['name'],
                      'report_ref': report_info['ref']
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
        """
        :param params: instance of type "MUSCLE_Params" (MUSCLE Input Params
           ** ** MUSCLE_prot(): input_ref must be FeatureSet ** MUSCLE_nuc():
           input_ref must be FeatureSet, SingleEndLibrary, or AssemblySet) ->
           structure: parameter "workspace_name" of type "workspace_name" (**
           The workspace object refs are of form: ** **    objects =
           ws.get_objects([{'ref':
           params['workspace_id']+'/'+params['obj_name']}]) ** ** "ref" means
           the entire name combining the workspace id and the object name **
           "id" is a numerical identifier of the workspace or object, and
           should just be used for workspace ** "name" is a string identifier
           of a workspace or object.  This is received from Narrative.),
           parameter "desc" of String, parameter "input_ref" of type
           "data_obj_ref", parameter "output_name" of type "data_obj_name",
           parameter "genome_disp_name_config" of String, parameter
           "maxiters" of Long, parameter "maxhours" of Double
        :returns: instance of type "MUSCLE_Output" (MUSCLE Output) ->
           structure: parameter "report_name" of type "data_obj_name",
           parameter "report_ref" of type "data_obj_ref"
        """
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
        [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I,
         WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = range(11)  # object_info tuple
        
        row_labels = {}


        #### do some basic checks
        #
        if 'workspace_name' not in params:
            raise ValueError('workspace_name parameter is required')
        if 'input_ref' not in params:
            raise ValueError('input_ref parameter is required')
        if 'output_name' not in params:
            raise ValueError('output_name parameter is required')


        #### Get the input_ref object
        ##
#        input_forward_reads_file_compression = None
#        sequencing_tech = 'N/A'
        try:
            ws = workspaceService(self.workspaceURL, token=ctx['token'])
            objects = ws.get_objects([{'ref': params['input_ref']}])
            data = objects[0]['data']
            info = objects[0]['info']
            input_name = info[1]
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
            raise ValueError('Unable to fetch input_ref object from workspace: ' + str(e))


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
                input_forward_reads_file_handle = open(input_forward_reads_file_path, 'w')
                self.log(console, 'downloading reads file: '+str(input_forward_reads_file_path))
                headers = {'Authorization': 'OAuth '+ctx['token']}
                r = requests.get(input_forward_reads['url']+'/node/'+input_forward_reads['id']+'?download', stream=True, headers=headers)
                for chunk in r.iter_content(1024):
                    input_forward_reads_file_handle.write(chunk)
                input_forward_reads_file_handle.close()
                self.log(console, 'done')
                ### END NOTE


                # remove carriage returns
                new_file_path = input_forward_reads_file_path+"-CRfree"
                new_file_handle = open(new_file_path, 'w')
                input_forward_reads_file_handle = open(input_forward_reads_file_path, 'r')
                for line in input_forward_reads_file_handle:
                    line = re.sub("\r","",line)
                    new_file_handle.write(line)
                input_forward_reads_file_handle.close()
                new_file_handle.close()
                input_forward_reads_file_path = new_file_path

                # convert FASTQ to FASTA (if necessary)
                new_file_path = input_forward_reads_file_path+".fna"
                new_file_handle = open(new_file_path, 'w')
                if input_forward_reads_file_compression == 'gz':
                    input_forward_reads_file_handle = gzip.open(input_forward_reads_file_path, 'r')
                else:
                    input_forward_reads_file_handle = open(input_forward_reads_file_path, 'r')
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
            
            genome_id_feature_id_delim = '.f:'

            # retrieve sequences for features
            input_featureSet = data

            genomeObjName = {}
            genomeObjVer  = {}
            genomeSciName = {}
            genome2Features = {}
            new_id = {}
            featureSet_elements = input_featureSet['elements']
            if 'element_ordering' in input_featureSet and input_featureSet['element_ordering']:
                feature_order = input_featureSet['element_ordering']
            else:
                feature_order = sorted(featureSet_elements.keys())
            for fId in feature_order:
                genomeRef = featureSet_elements[fId][0]
                if genomeRef not in genome2Features:
                    genome2Features[genomeRef] = []
                    new_id[genomeRef] = {}
                if genome_id_feature_id_delim in fId:
                    [genome_id, feature_id] = fId.split(genome_id_feature_id_delim)
                else:
                    feature_id = fId
                genome2Features[genomeRef].append(feature_id)
                this_id = genomeRef + genome_id_feature_id_delim + feature_id
                new_id[genomeRef][fId] = this_id

            # export features to FASTA file
            input_forward_reads_file_path = os.path.join(self.scratch, input_name+".fasta")
            self.log(console, 'writing fasta file: '+input_forward_reads_file_path)
            records_by_fid = dict()
            proteins_found = 0
            for genomeRef in genome2Features:
                genome_obj = ws.get_objects([{'ref':genomeRef}])[0]
                genome_type = re.sub('-[0-9]+\.[0-9]+$', "", genome_obj['info'][TYPE_I])
                genomeObjName[genomeRef] = genome_obj['info'][NAME_I]
                genomeObjVer[genomeRef]  = genome_obj['info'][VERSION_I]
                these_genomeFeatureIds = genome2Features[genomeRef]

                # Genome
                if genome_type == 'KBaseGenomes.Genome':
                    genome = genome_obj['data']
                    genomeSciName[genomeRef] = genome['scientific_name']

                    for feature in genome['features']:
                        if feature['id'] in these_genomeFeatureIds:
                            if 'protein_translation' not in feature or feature['protein_translation'] == None:
                                self.log(invalid_msgs,"bad CDS Feature "+feature['id']+": no protein_translation found")
                                continue
                            else:
                                #self.log(console,"kbase_id: '"+feature['id']+"'")  # DEBUG
                                this_id = genomeRef + genome_id_feature_id_delim + feature['id']
                                this_id = re.sub ('\s', '_', this_id)
                                short_feature_id = re.sub("^.*\.([^\.]+)\.([^\.]+)$", r"\1.\2", feature['id'])

                                genome_disp_name = ''
                                if 'obj_name' in params.get('genome_disp_name_config'):
                                    genome_disp_name += genomeObjName[genomeRef]
                                    if 'ver' in params.get('genome_disp_name_config'):
                                        genome_disp_name += '.v'+str(genomeObjVer[genomeRef])

                                    if genome_type == "KBaseGenomes.Genome" and \
                                       'sci_name' in params.get('genome_disp_name_config'):
                                        genome_disp_name += ': '+genomeSciName[genomeRef]
                                else:
                                    genome_disp_name = genomeObjName[genomeRef]
                
                                row_labels[this_id] = genome_disp_name+' - '+short_feature_id

                                #record = SeqRecord(Seq(feature['protein_translation']), id=feature['id'], description=genome['id'])
                                record = SeqRecord(Seq(feature['protein_translation']), id=this_id, description=genome['id'])
                                proteins_found += 1
                                records_by_fid[this_id] = record

                # AnnotatedMetagenomeAssembly
                elif genome_type == 'KBaseMetagenomes.AnnotatedMetagenomeAssembly':
                    ama_features = self._get_features_from_AnnotatedMetagenomeAssembly (ctx, genomeRef)
                    for feature in ama_features:
                        if feature['id'] in these_genomeFeatureIds:
                            if not feature.get('type'):
                                raise ValueError ("No type for AMA feature "+feature['id'])
                            if feature['type'] != 'CDS':
                                self.log ("skipping non-CDS AMA feature "+feature['id'])
                                continue
                            if not feature.get('protein_translatkon'):
                                self.log(console,"AMA CDS Feature "+feature['id']+": no protein_translation found.  Auto-translatiing from dna_sequence")
                                prot_translation = self.TranslateNucToProtSeq(ctx,
                                                                              {'nuc_seq': feature['dna_sequence'],
                                                                               'genetic_code': '11'})
                            else:
                                prot_translation = feature['protein_translation']
                            this_id = genomeRef + genome_id_feature_id_delim + feature['id']
                            short_feature_id = re.sub("^.*\.([^\.]+)\.([^\.]+)$", r"\1.\2", feature['id'])
                            genome_disp_name = genomeObjName[genomeRef]
                            row_labels[this_id] = genome_disp_name+' - '+short_feature_id

                            record = SeqRecord(Seq(prot_translation), id=this_id, description=genomeObjName[genomeRef])
                            proteins_found += 1
                            records_by_fid[this_id] = record

                                
            if proteins_found < 2:
                self.log(invalid_msgs,"Less than 2 protein Features (CDS) found.  exiting...")
            else:
                records = []
                for fId in feature_order:
                    genomeRef = featureSet_elements[fId][0]
                    records.append(records_by_fid[new_id[genomeRef][fId]])
                SeqIO.write(records, input_forward_reads_file_path, "fasta")

        # Missing proper input_input_type
        #
        else:
            raise ValueError('Cannot yet handle input_ref type of: '+input_type_name)


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
            output_aln_file_path = os.path.join(output_dir, params['output_name']+'-MSA.fasta')
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
                line = p.stdout.readline().decode()
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
            output_aln_file_handle = open (output_aln_file_path, 'r')

            output_fasta_buf = []
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
                    if row_labels:
                        this_id = leading_chars_pattern.findall(header)[0]
                        this_row_label = re.sub ('\s', '_', row_labels[this_id])
                        output_fasta_buf.append('>'+this_row_label)
                    else:
                        output_fasta_buf.append(line)

                    if last_header != None:
                        last_id = leading_chars_pattern.findall(last_header)[0]
                        row_order.append(last_id)
                        #self.log(console,"ID: '"+last_id+"'\nALN: '"+last_seq+"'")  # DEBUG
                        #report += last_id+"\t"+last_seq+"\n"
                        alignment[last_id] = last_seq
                        if alignment_length == None:
                            alignment_length = len(last_seq)
                        elif alignment_length != len(last_seq):
                            raise ValueError ("unequal alignment row for "+last_header+": '"+last_seq+"'")
                    last_header = header
                    last_seq = ''
                else:
                    last_seq += line
                    output_fasta_buf.append(line)
            if last_header != None:
                last_id = leading_chars_pattern.findall(last_header)[0]
                row_order.append(last_id)
                #self.log(console,"ID: '"+last_id+"'\nALN: '"+last_seq+"'")  # DEBUG
                #report += last_id+"\t"+last_seq+"\n"
                alignment[last_id] = last_seq
                if alignment_length == None:
                    alignment_length = len(last_seq)
                elif alignment_length != len(last_seq):
                    raise ValueError ("unequal alignment row for "+last_header+": '"+last_seq+"'")
        
            output_aln_file_handle.close()

            # write remapped ids
            with open(output_aln_file_path, 'w') as output_aln_file_handle:
                output_aln_file_handle.write("\n".join(output_fasta_buf)+"\n")


        # load the method provenance from the context object
        #
        self.log(console,"SETTING PROVENANCE")  # DEBUG
        provenance = [{}]
        if 'provenance' in ctx:
            provenance = ctx['provenance']
        # add additional info to provenance here, in this case the input data object reference
        provenance[0]['input_ws_objects'] = []
        provenance[0]['input_ws_objects'].append(params['input_ref'])
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
            if row_labels:
                output_MSA['default_row_labels'] = row_labels

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


            # create CLW formatted output file
            max_row_width = 60
            id_aln_gap_width = 1
            gap_chars = ''
            for sp_i in range(id_aln_gap_width):
                gap_chars += ' '
# DNA
#            strong_groups = { 'AG': True,
#                              'CTU': True
#                              }
#            weak_groups = None
# PROTEINS
            strong_groups = { 'AST':  True,
                              'EKNQ': True,
                              'HKNQ': True,
                              'DENQ': True,
                              'HKQR': True,
                              'ILMV': True,
                              'FILM': True,
                              'HY':   True,
                              'FWY':  True
                              }
            weak_groups = { 'ACS':    True,
                            'ATV':    True,
                            'AGS':    True,
                            'KNST':   True,
                            'APST':   True,
                            'DGNS':   True,
                            'DEKNQS': True,
                            'DEHKNQ': True,
                            'EHKNQR': True,
                            'FILMV':  True,
                            'FHY':    True
                            }

            clw_buf = []
            clw_buf.append ('CLUSTALW format of MUSCLE alignment '+MSA_name+': '+MSA_description)
            clw_buf.append ('')

            long_id_len = 0
            aln_pos_by_id = dict()
            for row_id in row_order:
                aln_pos_by_id[row_id] = 0
                if row_labels:
                    row_id_disp = row_labels[row_id]
                else:
                    row_id_disp = row_id
                if long_id_len < len(row_id_disp):
                    long_id_len = len(row_id_disp)

            full_row_cnt = alignment_length // max_row_width
            if alignment_length % max_row_width == 0:
                full_row_cnt -= 1
            for chunk_i in range (full_row_cnt + 1):
                for row_id in row_order:
                    if row_labels:
                        row_id_disp = re.sub('\s', '_', row_labels[row_id])
                    else:
                        row_id_disp = row_id
                    for sp_i in range (long_id_len-len(row_id_disp)):
                        row_id_disp += ' '

                    aln_chunk_upper_bound = (chunk_i+1)*max_row_width
                    if aln_chunk_upper_bound > alignment_length:
                        aln_chunk_upper_bound = alignment_length
                    aln_chunk = alignment[row_id][chunk_i*max_row_width:aln_chunk_upper_bound]
                    for c in aln_chunk:
                        if c != '-':
                            aln_pos_by_id[row_id] += 1

                    clw_buf.append (row_id_disp+gap_chars+aln_chunk+' '+str(aln_pos_by_id[row_id]))

                # conservation line
                cons_line = ''
                for pos_i in range(chunk_i*max_row_width, aln_chunk_upper_bound):
                    col_chars = dict()
                    seq_cnt = 0
                    for row_id in row_order:
                        char = alignment[row_id][pos_i]
                        if char != '-':
                            seq_cnt += 1
                            col_chars[char] = True
                    if seq_cnt <= 1:
                        cons_char = ' '
                    elif len(col_chars.keys()) == 1:
                        cons_char = '*'
                    else:
                        strong = False
                        for strong_group in strong_groups.keys():
                            this_strong_group = True
                            for seen_char in col_chars.keys():
                                if seen_char not in strong_group:
                                    this_strong_group = False
                                    break
                            if this_strong_group:
                                strong = True
                                break
                        if not strong:
                            weak = False
                            if weak_groups is not None:
                                for weak_group in weak_groups.keys():
                                    this_weak_group = True
                                    for seen_char in col_chars.keys():
                                        if seen_char not in weak_group:
                                            this_strong_group = False
                                            break
                                    if this_weak_group:
                                        weak = True
                        if strong:
                            cons_char = ':'
                        elif weak:
                            cons_char = '.'
                        else:
                            cons_char = ' '
                    cons_line += cons_char

                lead_space = ''
                for sp_i in range(long_id_len):
                    lead_space += ' '
                lead_space += gap_chars

                clw_buf.append(lead_space+cons_line)
                clw_buf.append('')

            # write clw to file
            clw_buf_str = "\n".join(clw_buf)+"\n"
            output_clw_file_path = os.path.join(output_dir, input_name+'-MSA.clw')
            with open (output_clw_file_path, 'w') as output_clw_file_handle:
                output_clw_file_handle.write(clw_buf_str)


            # upload MUSCLE FASTA output to SHOCK for file_links
            dfu = DFUClient(self.callbackURL)
            try:
                output_upload_ret = dfu.file_to_shock({'file_path': output_aln_file_path,
# DEBUG
#                                                      'make_handle': 0,
#                                                      'pack': 'zip'})
                                                       'make_handle': 0})
            except:
                raise ValueError ('error loading aln_out file to shock')

            # upload MUSCLE CLW output to SHOCK for file_links
            try:
                output_clw_upload_ret = dfu.file_to_shock({'file_path': output_clw_file_path,
# DEBUG
#                                                      'make_handle': 0,
#                                                      'pack': 'zip'})
                                                           'make_handle': 0})
            except:
                raise ValueError ('error loading clw_out file to shock')


            # make HTML reports
            #
            # HERE



            # build output report object
            #
            self.log(console,"BUILDING REPORT")  # DEBUG

            reportName = 'muscle_report_'+str(uuid.uuid4())
            reportObj = {
                'objects_created':[{'ref':params['workspace_name']+'/'+params['output_name'],
                                    'description':'MUSCLE_prot MSA'}],
                #'message': '',
                'message': clw_buf_str,
                'file_links': [],
                'workspace_name': params['workspace_name'],
                'report_object_name': reportName
                }
            reportObj['file_links'] = [{'shock_id': output_upload_ret['shock_id'],
                                        'name': params['output_name']+'-MUSCLE_prot.FASTA',
                                        'label': 'MUSCLE_prot FASTA'
                                        },
                                       {'shock_id': output_clw_upload_ret['shock_id'],
                                        'name': params['output_name']+'-MUSCLE_prot.CLW',
                                        'label': 'MUSCLE_prot CLUSTALW'
                                        }]

            # save report object
            #
            SERVICE_VER = 'release'
            reportClient = KBaseReport(self.callbackURL, token=ctx['token'], service_ver=SERVICE_VER)
            #report_info = report.create({'report':reportObj, 'workspace_name':params['workspace_name']})
            report_info = reportClient.create_extended_report(reportObj)                                       

        else:  # len(invalid_msgs) > 0
            reportName = 'muscle_report_'+str(uuid.uuid4())
            report += "FAILURE:\n\n"+"\n".join(invalid_msgs)+"\n"
            reportObj = {
                'objects_created':[],
                'text_message':report
                }

            ws = workspaceService(self.workspaceURL, token=ctx['token'])
            report_obj_info = ws.save_objects({
                    #'id':info[6],
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

            report_info = dict()
            report_info['name'] = report_obj_info[1]
            report_info['ref'] = str(report_obj_info[6])+'/'+str(report_obj_info[0])+'/'+str(report_obj_info[4])


        # done
        returnVal = { 'report_name': report_info['name'],
                      'report_ref': report_info['ref']
                      }
        self.log(console,"MUSCLE_prot DONE")

        #END MUSCLE_prot

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method MUSCLE_prot return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]
    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK",
                     'message': "",
                     'version': self.VERSION,
                     'git_url': self.GIT_URL,
                     'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
