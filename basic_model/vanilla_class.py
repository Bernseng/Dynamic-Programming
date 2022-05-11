from types import SimpleNamespace # model infrastructure
from model_functions import settings, setup, allocate, solve, simulate

class VanillaClass:
    def __init__(self,name='',**kwargs):
        """ defines default attributes """

        # a. name
        self.name = name

        # b. unpack namespaces from pre-defined functions

        # i. settings
        assert hasattr(self,'settings'), 'the model must have defined a .settings() method'
        self.settings()

        for ns in self.namespaces:
            setattr(self,ns,SimpleNamespace())

        # ii setup (parameters)
        assert hasattr(self,'setup'), 'the model must have defined a .setup() method'
        self.setup() # call setup method to pass parameters to the model

        # iii. update
        self.__update(kwargs)
            
        # vi. allocate
        assert hasattr(self,'allocate'), 'the model must have defined an .allocate() method'
        self.allocate()
    

    def __update(self,upd_dict):
        """ update """

        for nskey,values in upd_dict.items():
            assert nskey in self.namespaces, f'{nskey} is not a namespace'
            assert type(values) is dict, f'{nskey} must be dict'
            for key,value in values.items():
                assert hasattr(getattr(self,nskey),key), f'{key} is not in {nskey}' 
                setattr(getattr(self,nskey),key,value) 
    
    settings = settings
    setup = setup
    allocate = allocate
    solve = solve
    simulate = simulate