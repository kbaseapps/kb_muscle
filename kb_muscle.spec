/*
** A KBase module: kb_muscle
**
** This module runs MUSCLE to make MSAs of either DNA or PROTEIN sequences.  "MUSCLE nuc" will build nucleotide alignments, even for protein coding genes.  "MUSCLE prot" will build protein sequence alignments, and will ignore any features that do not code for proteins.
** 
*/

module kb_muscle {

    /* 
    ** The workspace object refs are of form:
    **
    **    objects = ws.get_objects([{'ref': params['workspace_id']+'/'+params['obj_name']}])
    **
    ** "ref" means the entire name combining the workspace id and the object name
    ** "id" is a numerical identifier of the workspace or object, and should just be used for workspace
    ** "name" is a string identifier of a workspace or object.  This is received from Narrative.
    */
    typedef string workspace_name;
    typedef string data_obj_name;
    typedef string data_obj_ref;


    /* MUSCLE Input Params
    */
    typedef structure {
        workspace_name workspace_name;
	string         desc;
	data_obj_name  input_ref;
        data_obj_name  output_name;
	int            maxiters;
	float          maxhours;
    } MUSCLE_Params;


    /* MUSCLE Output
    */
    typedef structure {
	data_obj_name report_name;
	data_obj_ref  report_ref;
/*       data_obj_ref  output_filtered_ref;
*
*        int n_initial_seqs;
*        int n_seqs_matched;
*        int n_seqs_notmatched;
*/
    } MUSCLE_Output;
	

    /*  Methods for MSA building of either DNA or PROTEIN sequences
    **
    **    overloading as follows:
    **        input_name: SingleEndLibrary, FeatureSet
    **        output_name: MSA
    */
    funcdef MUSCLE_nuc (MUSCLE_Params params)  returns (MUSCLE_Output) authentication required;
    funcdef MUSCLE_prot (MUSCLE_Params params)  returns (MUSCLE_Output) authentication required;
};
