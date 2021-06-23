# -*- coding: utf-8 -*-
import os  # noqa: F401
import time
import unittest
import uuid
import json
import shutil
from os import environ

try:
    from ConfigParser import ConfigParser  # py2
except:
    from configparser import ConfigParser  # py3

from pprint import pprint  # noqa: F401

from installed_clients.WorkspaceClient import Workspace as workspaceService
from installed_clients.GenomeFileUtilClient import GenomeFileUtil
from installed_clients.authclient import KBaseAuth as _KBaseAuth
from kb_muscle.kb_muscleImpl import kb_muscle
from kb_muscle.kb_muscleServer import MethodContext


class kb_muscleTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        token = environ.get('KB_AUTH_TOKEN', None)
        config_file = environ.get('KB_DEPLOYMENT_CONFIG', None)
        cls.cfg = {}
        config = ConfigParser()
        config.read(config_file)
        for nameval in config.items('kb_muscle'):
            cls.cfg[nameval[0]] = nameval[1]
        # Getting username from Auth profile for token
        authServiceUrl = cls.cfg['auth-service-url']
        auth_client = _KBaseAuth(authServiceUrl)
        user_id = auth_client.get_user(token)
        # WARNING: don't call any logging methods on the context object,
        # it'll result in a NoneType error
        cls.ctx = MethodContext(None)
        cls.ctx.update({'token': token,
                        'user_id': user_id,
                        'provenance': [
                            {'service': 'kb_muscle',
                             'method': 'please_never_use_it_in_production',
                             'method_params': []
                             }],
                        'authenticated': 1})
        cls.wsURL = cls.cfg['workspace-url']
        cls.wsClient = workspaceService(cls.wsURL)
        cls.serviceImpl = kb_muscle(cls.cfg)
        cls.scratch = cls.cfg['scratch']
        cls.callback_url = os.environ['SDK_CALLBACK_URL']

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')

    def getWsClient(self):
        return self.__class__.wsClient

    def getWsName(self):
        if hasattr(self.__class__, 'wsName'):
            return self.__class__.wsName
        suffix = int(time.time() * 1000)
        wsName = "test_kb_muscle_" + str(suffix)
        ret = self.getWsClient().create_workspace({'workspace': wsName})  # noqa
        self.__class__.wsName = wsName
        return wsName

    def getImpl(self):
        return self.__class__.serviceImpl

    def getContext(self):
        return self.__class__.ctx


    # call this method to get the WS object info of a Genome
    #   (will upload the example data if this is the first time the method is called during tests)
    def getGenomeInfo(self, genome_basename, item_i=0):
        if hasattr(self.__class__, 'genomeInfo_list'):
            try:
                info = self.__class__.genomeInfo_list[item_i]
                name = self.__class__.genomeName_list[item_i]
                if info != None:
                    if name != genome_basename:
                        self.__class__.genomeInfo_list[item_i] = None
                        self.__class__.genomeName_list[item_i] = None
                    else:
                        return info
            except:
                pass

        # 1) transform genbank to kbase genome object and upload to ws
        shared_dir = "/kb/module/work/tmp"
        genome_data_file = 'data/genomes/'+genome_basename+'.gbff.gz'
        genome_file = os.path.join(shared_dir, os.path.basename(genome_data_file))
        shutil.copy(genome_data_file, genome_file)

        SERVICE_VER = 'release'
        #SERVICE_VER = 'dev'
        GFU = GenomeFileUtil(os.environ['SDK_CALLBACK_URL'],
                             token=self.getContext()['token'],
                             service_ver=SERVICE_VER
                         )
        print ("UPLOADING genome: "+genome_basename+" to WORKSPACE "+self.getWsName()+" ...")
        genome_upload_result = GFU.genbank_to_genome({'file': {'path': genome_file },
                                                      'workspace_name': self.getWsName(),
                                                      'genome_name': genome_basename
                                                  })
#                                                  })[0]
        pprint(genome_upload_result)
        genome_ref = genome_upload_result['genome_ref']
        new_obj_info = self.getWsClient().get_object_info_new({'objects': [{'ref': genome_ref}]})[0]

        # 2) store it
        if not hasattr(self.__class__, 'genomeInfo_list'):
            self.__class__.genomeInfo_list = []
            self.__class__.genomeName_list = []
        for i in range(item_i+1):
            try:
                assigned = self.__class__.genomeInfo_list[i]
            except:
                self.__class__.genomeInfo_list.append(None)
                self.__class__.genomeName_list.append(None)

        self.__class__.genomeInfo_list[item_i] = new_obj_info
        self.__class__.genomeName_list[item_i] = genome_basename
        return new_obj_info

    # call this method to get the WS object info of an AnnotatedMetagenomeAssembly
    #   (will upload the example data if this is the first time the method is called during tests)
    def getAMAInfo(self, ama_basename, item_i=0):
        if hasattr(self.__class__, 'amaInfo_list'):
            try:
                info = self.__class__.amaInfo_list[item_i]
                name = self.__class__.amaName_list[item_i]
                if info != None:
                    if name != ama_basename:
                        self.__class__.amaInfo_list[item_i] = None
                        self.__class__.amaName_list[item_i] = None
                    else:
                        return info
            except:
                pass

        # 1) transform GFF+FNA to kbase AMA object and upload to ws
        shared_dir = "/kb/module/work/tmp"
        ama_gff_srcfile = 'data/amas/'+ama_basename+'.gff'
        ama_fna_srcfile = 'data/amas/'+ama_basename+'.fa'
        ama_gff_dstfile = os.path.join(shared_dir, os.path.basename(ama_gff_srcfile))
        ama_fna_dstfile = os.path.join(shared_dir, os.path.basename(ama_fna_srcfile))
        shutil.copy(ama_gff_srcfile, ama_gff_dstfile)
        shutil.copy(ama_fna_srcfile, ama_fna_dstfile)

        try:
            SERVICE_VER = 'release'
            #SERVICE_VER = 'dev'
            GFU = GenomeFileUtil(os.environ['SDK_CALLBACK_URL'],
                                 token=self.getContext()['token'],
                                 service_ver=SERVICE_VER
            )
        except:
            raise ValueError ("unable to obtain GenomeFileUtil client")
        print ("UPLOADING AMA: "+ama_basename+" to WORKSPACE "+self.getWsName()+" ...")
        ama_upload_params = {
            "workspace_name": self.getWsName(),
            "genome_name": ama_basename,
            "fasta_file": {"path": ama_fna_dstfile},
            "gff_file": {"path": ama_gff_dstfile},
            "source": "GFF",
            "scientific_name": "TEST AMA",
            "generate_missing_genes": "True"
        }        
        try:
            ama_upload_result = GFU.fasta_gff_to_metagenome(ama_upload_params)
        except:
            raise ValueError("unable to upload test AMA data object")
        print ("AMA UPLOADED")
        pprint(ama_upload_result)

        ama_ref = ama_upload_result['metagenome_ref']
        new_obj_info = self.getWsClient().get_object_info_new({'objects': [{'ref': ama_ref}]})[0]

        # 2) store it
        if not hasattr(self.__class__, 'amaInfo_list'):
            self.__class__.amaInfo_list = []
            self.__class__.amaName_list = []
        for i in range(item_i+1):
            try:
                assigned = self.__class__.amaInfo_list[i]
            except:
                self.__class__.amaInfo_list.append(None)
                self.__class__.amaName_list.append(None)

        self.__class__.amaInfo_list[item_i] = new_obj_info
        self.__class__.amaName_list[item_i] = ama_basename
        return new_obj_info


    ##############
    # UNIT TESTS #
    ##############

    #### test_MUSCLE_nuc_01()
    ##
    # HIDE @unittest.skip("skipped test_MUSCLE_nuc_01()")  # uncomment to skip
    ##
    def test_MUSCLE_nuc_01(self):
        obj_out_name = 'test_MUSCLE_nuc'
        obj_out_type = 'KBaseTrees.MSA'

        [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = range(11)  # object_info tuple
        
        genomeInfo_0 = self.getGenomeInfo('GCF_001566335.1_ASM156633v1_genomic', 0)  # E. coli K-12 MG1655
        genomeInfo_1 = self.getGenomeInfo('GCF_000021385.1_ASM2138v1_genomic', 1)    # D. vulgaris str. 'Miyazaki F
        genomeInfo_2 = self.getGenomeInfo('GCF_001721825.1_ASM172182v1_genomic', 2)  # Pseudomonas aeruginosa
        genomeInfo_3 = self.getGenomeInfo('GCF_000020845.1_ASM2084v1_genomic', 3)    # Vibrio fischeri MJ11
        genome_ref_0 = '/'.join([str(genomeInfo_0[WSID_I]),
                                 str(genomeInfo_0[OBJID_I]),
                                 str(genomeInfo_0[VERSION_I])])
        genome_ref_1 = '/'.join([str(genomeInfo_1[WSID_I]),
                                 str(genomeInfo_1[OBJID_I]),
                                 str(genomeInfo_1[VERSION_I])])
        genome_ref_2 = '/'.join([str(genomeInfo_2[WSID_I]),
                                 str(genomeInfo_2[OBJID_I]),
                                 str(genomeInfo_2[VERSION_I])])
        genome_ref_3 = '/'.join([str(genomeInfo_3[WSID_I]),
                                 str(genomeInfo_3[OBJID_I]),
                                 str(genomeInfo_3[VERSION_I])])

        genome_refs = [ genome_ref_0, genome_ref_1, genome_ref_2, genome_ref_3 ]
        
        # RplF
        feature_ids = [ 
            'AWN69_RS09265',
            'DVMF_RS00495',
            'A6701_RS03520',
            'VFMJ11_RS08285'
        ]

        featureSet_obj = { 'description': 'test genomeSet',
                           'element_ordering': [],
                           'elements': {}
                        }
        for i,genome_ref in enumerate(genome_refs):
            feature_id = feature_ids[i]
            featureSet_obj['element_ordering'].append(feature_id)
            featureSet_obj['elements'][feature_id] = [genome_ref]

        provenance = [{}]
        featureSet_info = self.getWsClient().save_objects({
            'workspace': self.getWsName(), 
            'objects': [
                {
                    'type': 'KBaseCollections.FeatureSet',
                    'data': featureSet_obj,
                    'name': 'test_featureSet',
                    'meta': {},
                    'provenance': provenance
                }
            ]})[0]

        [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = list(range(11))  # object_info tuple
        featureSet_ref = str(featureSet_info[WSID_I])+'/'+str(featureSet_info[OBJID_I])+'/'+str(featureSet_info[VERSION_I])

        output_dir = os.path.join(self.scratch,'fasta_out.'+str(uuid.uuid4()))
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        parameters = {
                'desc':           'test MUSCLE nuc',
                'input_ref':      featureSet_ref,
                'output_name':    obj_out_name,
                'genome_disp_name_config':   'obj_name_ver_sci_name',
                'maxiters':       '16',
                'maxhours':       '0.5',
                'workspace_name': self.getWsName()
                }
        ret = self.getImpl().MUSCLE_nuc(self.getContext(), parameters)[0]
        self.assertIsNotNone(ret['report_ref'])

        # check created obj
        #report_obj = self.getWsClient().get_objects2({'objects':[{'ref':ret['report_ref']}]})[0]['data']
        report_obj = self.getWsClient().get_objects([{'ref':ret['report_ref']}])[0]['data']
        self.assertIsNotNone(report_obj['objects_created'][0]['ref'])

        created_obj_0_info = self.getWsClient().get_object_info_new({'objects':[{'ref':report_obj['objects_created'][0]['ref']}]})[0]
        self.assertEqual(created_obj_0_info[NAME_I], obj_out_name)
        self.assertEqual(created_obj_0_info[TYPE_I].split('-')[0], obj_out_type)
        pass
        

    #### test_MUSCLE_prot_01()
    ##
    # HIDE @unittest.skip("skipped test_MUSCLE_prot_01()")  # uncomment to skip
    ##
    def test_MUSCLE_prot_01(self):
        obj_out_name = 'test_MUSCLE_prot'
        obj_out_type = 'KBaseTrees.MSA'

        [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = range(11)  # object_info tuple

        genomeInfo_0 = self.getGenomeInfo('GCF_001566335.1_ASM156633v1_genomic', 0)  # E. coli K-12 MG1655
        genomeInfo_1 = self.getGenomeInfo('GCF_000021385.1_ASM2138v1_genomic', 1)    # D. vulgaris str. 'Miyazaki F
        genomeInfo_2 = self.getGenomeInfo('GCF_001721825.1_ASM172182v1_genomic', 2)  # Pseudomonas aeruginosa
        genomeInfo_3 = self.getGenomeInfo('GCF_000020845.1_ASM2084v1_genomic', 3)    # Vibrio fischeri MJ11
        genome_ref_0 = '/'.join([str(genomeInfo_0[WSID_I]),
                                 str(genomeInfo_0[OBJID_I]),
                                 str(genomeInfo_0[VERSION_I])])
        genome_ref_1 = '/'.join([str(genomeInfo_1[WSID_I]),
                                 str(genomeInfo_1[OBJID_I]),
                                 str(genomeInfo_1[VERSION_I])])
        genome_ref_2 = '/'.join([str(genomeInfo_2[WSID_I]),
                                 str(genomeInfo_2[OBJID_I]),
                                 str(genomeInfo_2[VERSION_I])])
        genome_ref_3 = '/'.join([str(genomeInfo_3[WSID_I]),
                                 str(genomeInfo_3[OBJID_I]),
                                 str(genomeInfo_3[VERSION_I])])

        genome_refs = [ genome_ref_0, genome_ref_1, genome_ref_2, genome_ref_3 ]
        
        # RplF
        feature_ids = [ 
            'AWN69_RS09265',
            'DVMF_RS00495',
            'A6701_RS03520',
            'VFMJ11_RS08285'
        ]

        featureSet_obj = { 'description': 'test genomeSet',
                           'element_ordering': [],
                           'elements': {}
                        }
        for i,genome_ref in enumerate(genome_refs):
            feature_id = feature_ids[i]
            featureSet_obj['element_ordering'].append(feature_id)
            featureSet_obj['elements'][feature_id] = [genome_ref]

        provenance = [{}]
        featureSet_info = self.getWsClient().save_objects({
            'workspace': self.getWsName(), 
            'objects': [
                {
                    'type': 'KBaseCollections.FeatureSet',
                    'data': featureSet_obj,
                    'name': 'test_featureSet',
                    'meta': {},
                    'provenance': provenance
                }
            ]})[0]

        [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = list(range(11))  # object_info tuple
        featureSet_ref = str(featureSet_info[WSID_I])+'/'+str(featureSet_info[OBJID_I])+'/'+str(featureSet_info[VERSION_I])

        output_dir = os.path.join(self.scratch,'fasta_out.'+str(uuid.uuid4()))
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        parameters = {
                'desc':           'test MUSCLE nuc',
                'input_ref':      featureSet_ref,
                'output_name':    obj_out_name,
                'genome_disp_name_config':   'obj_name_ver_sci_name',
                'maxiters':       '16',
                'maxhours':       '0.5',
                'workspace_name': self.getWsName()
                }
        ret = self.getImpl().MUSCLE_prot(self.getContext(), parameters)[0]
        self.assertIsNotNone(ret['report_ref'])

        # check created obj
        #report_obj = self.getWsClient().get_objects2({'objects':[{'ref':ret['report_ref']}]})[0]['data']
        report_obj = self.getWsClient().get_objects([{'ref':ret['report_ref']}])[0]['data']
        self.assertIsNotNone(report_obj['objects_created'][0]['ref'])

        created_obj_0_info = self.getWsClient().get_object_info_new({'objects':[{'ref':report_obj['objects_created'][0]['ref']}]})[0]
        self.assertEqual(created_obj_0_info[NAME_I], obj_out_name)
        self.assertEqual(created_obj_0_info[TYPE_I].split('-')[0], obj_out_type)
        pass


    #### test_MUSCLE_prot_02():  test FeatureSets with AMA features
    ##
    # HIDE @unittest.skip("skipped test_MUSCLE_prot_02()")  # uncomment to skip
    ##
    def test_MUSCLE_prot_02(self):
        obj_out_name = 'test_MUSCLE_prot_02'
        obj_out_type = 'KBaseTrees.MSA'

        [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = range(11)  # object_info tuple

        genomeInfo_0 = self.getGenomeInfo('GCF_001566335.1_ASM156633v1_genomic', 0)  # E. coli K-12 MG1655
        genomeInfo_1 = self.getGenomeInfo('GCF_000021385.1_ASM2138v1_genomic', 1)    # D. vulgaris str. 'Miyazaki F
        genomeInfo_2 = self.getGenomeInfo('GCF_001721825.1_ASM172182v1_genomic', 2)  # Pseudomonas aeruginosa
        genomeInfo_3 = self.getGenomeInfo('GCF_000020845.1_ASM2084v1_genomic', 3)    # Vibrio fischeri MJ11
        genome_ref_0 = '/'.join([str(genomeInfo_0[WSID_I]),
                                 str(genomeInfo_0[OBJID_I]),
                                 str(genomeInfo_0[VERSION_I])])
        genome_ref_1 = '/'.join([str(genomeInfo_1[WSID_I]),
                                 str(genomeInfo_1[OBJID_I]),
                                 str(genomeInfo_1[VERSION_I])])
        genome_ref_2 = '/'.join([str(genomeInfo_2[WSID_I]),
                                 str(genomeInfo_2[OBJID_I]),
                                 str(genomeInfo_2[VERSION_I])])
        genome_ref_3 = '/'.join([str(genomeInfo_3[WSID_I]),
                                 str(genomeInfo_3[OBJID_I]),
                                 str(genomeInfo_3[VERSION_I])])

        amaInfo_0 = self.getAMAInfo("ama_with_2_rplF_genes", 0)
        ama_ref_0 = '/'.join([str(amaInfo_0[WSID_I]),
                              str(amaInfo_0[OBJID_I]),
                              str(amaInfo_0[VERSION_I])])

        genome_refs = [ genome_ref_0, genome_ref_1, genome_ref_2, genome_ref_3, ama_ref_0, ama_ref_0 ]
        
        # RplF
        feature_ids = [ 
            'AWN69_RS09265',
            'DVMF_RS00495',
            'A6701_RS03520',
            'VFMJ11_RS08285',
            'NFFIPACM_29235',
            'NFFIPACM_75173'
        ]

        featureSet_obj = { 'description': 'test genomeSet',
                           'element_ordering': [],
                           'elements': {}
                        }
        for i,genome_ref in enumerate(genome_refs):
            feature_id = feature_ids[i]
            featureSet_obj['element_ordering'].append(feature_id)
            featureSet_obj['elements'][feature_id] = [genome_ref]

        provenance = [{}]
        featureSet_info = self.getWsClient().save_objects({
            'workspace': self.getWsName(), 
            'objects': [
                {
                    'type': 'KBaseCollections.FeatureSet',
                    'data': featureSet_obj,
                    'name': 'test_featureSet_02',
                    'meta': {},
                    'provenance': provenance
                }
            ]})[0]

        [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = list(range(11))  # object_info tuple
        featureSet_ref = str(featureSet_info[WSID_I])+'/'+str(featureSet_info[OBJID_I])+'/'+str(featureSet_info[VERSION_I])

        output_dir = os.path.join(self.scratch,'fasta_out_02.'+str(uuid.uuid4()))
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        parameters = {
                'desc':           'test MUSCLE nuc',
                'input_ref':      featureSet_ref,
                'output_name':    obj_out_name,
                'genome_disp_name_config':   'obj_name_ver_sci_name',
                'maxiters':       '16',
                'maxhours':       '0.5',
                'workspace_name': self.getWsName()
                }
        ret = self.getImpl().MUSCLE_prot(self.getContext(), parameters)[0]
        self.assertIsNotNone(ret['report_ref'])

        # check created obj
        #report_obj = self.getWsClient().get_objects2({'objects':[{'ref':ret['report_ref']}]})[0]['data']
        report_obj = self.getWsClient().get_objects([{'ref':ret['report_ref']}])[0]['data']
        self.assertIsNotNone(report_obj['objects_created'][0]['ref'])

        created_obj_0_info = self.getWsClient().get_object_info_new({'objects':[{'ref':report_obj['objects_created'][0]['ref']}]})[0]
        self.assertEqual(created_obj_0_info[NAME_I], obj_out_name)
        self.assertEqual(created_obj_0_info[TYPE_I].split('-')[0], obj_out_type)
        pass
