#BEGIN_HEADER
#END_HEADER


class kb_muscle:
    '''
    Module Name:
    kb_muscle

    Module Description:
    ** A KBase module: kb_muscle
**
** This module runs MUSCLE to make MSAs of either DNA or PROTEIN sequences
**
    '''

    ######## WARNING FOR GEVENT USERS #######
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    #########################################
    #BEGIN_CLASS_HEADER
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        #END_CONSTRUCTOR
        pass

    def MUSCLE_nuc(self, ctx, params):
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN MUSCLE_nuc
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
        #END MUSCLE_prot

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method MUSCLE_prot return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]
