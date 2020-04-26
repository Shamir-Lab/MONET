from module import Module
class Globals:
    """
    The class contains the state of a MONET run - all patients, modules, omic graphs, and configuration.
    """

    def __init__(self, len):
        self.patients = {}
        self.modules = {}
        self.converged_modules = {}
        self.active_modules = {}
        self.omics = {}
        self.actions = [0 for i in range(len)]
        self.index = 0
        self.min_mod_size = 1
        self.max_pats_per_action = 10
        self.gmm_params = {}
        return

    def kill_module(self, mod):
        pats = mod.get_patients().copy().items()
        for name, pat in pats:
            mod.remove_patient(pat)
        if isinstance(mod, int):
            del self.modules[mod]
        elif isinstance(mod, Module):
            for name, module in self.modules.items():
                if module == mod:
                    break
            del self.modules[name]
        else:
            raise Exception('unknown module id')
        return self

    def add_module(self, mod):
        self.modules.update({self.index : mod})
        self.index += 1
        return
