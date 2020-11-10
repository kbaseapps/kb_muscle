# -*- coding: utf-8 -*-
import os  # noqa: F401
import time
import unittest
import uuid
from os import environ

try:
    from ConfigParser import ConfigParser  # py2
except:
    from configparser import ConfigParser  # py3

from installed_clients.WorkspaceClient import Workspace as workspaceService
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


    #### MUSCLE_nuc_01()
    ##
    def test_MUSCLE_nuc_01(self):
        obj_out_name = 'test_MUSCLE_nuc'
        obj_out_type = 'KBaseTrees.MSA'

        reference_prok_genomes_WS = 'ReferenceDataManager'  # PROD and CI
        genome_refs = [
            reference_prok_genomes_WS+'/GCF_001566335.1/1',  # E. coli K-12 MG1655
            reference_prok_genomes_WS+'/GCF_000021385.1/1',  # D. vulgaris str. 'Miyazaki F'
            reference_prok_genomes_WS+'/GCF_001721825.1/1',  # Pseudomonas aeruginosa
            reference_prok_genomes_WS+'/GCF_000020845.1/1'  # Vibrio fischeri MJ11
        ]
        # DnaA
        feature_ids = [ 
            'AWN69_RS07145',
            'DVMF_RS00005',
            'A6701_RS00005',
            'VFMJ11_RS07145'
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
        

    #### MUSCLE_prot_01()
    ##
    def test_MUSCLE_prot_01(self):
        obj_out_name = 'test_MUSCLE_prot'
        obj_out_type = 'KBaseTrees.MSA'

        reference_prok_genomes_WS = 'ReferenceDataManager'  # PROD and CI
        genome_refs = [
            reference_prok_genomes_WS+'/GCF_001566335.1/1',  # E. coli K-12 MG1655
            reference_prok_genomes_WS+'/GCF_000021385.1/1',  # D. vulgaris str. 'Miyazaki F'
            reference_prok_genomes_WS+'/GCF_001721825.1/1',  # Pseudomonas aeruginosa
            reference_prok_genomes_WS+'/GCF_000020845.1/1'  # Vibrio fischeri MJ11
        ]
        # DnaA
        feature_ids = [ 
            'AWN69_RS07145',
            'DVMF_RS00005',
            'A6701_RS00005',
            'VFMJ11_RS07145'
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
