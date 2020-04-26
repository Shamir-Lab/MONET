import time
class Module:
    """
    Represents a single module (set of patients and omics).
    Contains functions to manipulate modules which are used by MONET actions
    (e.g. split module).
    """

    def __init__(self, glob_var, patients={}, omics={}, weight=0):
        self.patients = {}
        self.omics = {}
        self.weight = 0
        for name, pat in patients.items():
            self.add_patient(pat)
        for omic in omics:
            self.add_omic(omic, glob_var)
        if weight:
            if abs(self.weight - weight) > 0.0001:
		# This error should never occur!
                print("self.weight != weight")
                exit(1)
        glob_var.add_module(self)
        return

    def is_patient_in_module(self, pat):
        return pat in self.patients

    def get_patients(self):
        return self.patients

    def get_size(self):
        return len(self.patients)

    def add_patient(self, pat):
        self.patients.update({pat.get_name():pat})
        pat_lst = self.get_patients_names_as_list()
        for omic in self.omics.values():
            edges = omic.graph.subgraph(pat_lst).edges(pat.get_name(), data=True)
            for e in edges:
                self.weight += e[2]['weight']
        pat.set_module(self)
        return self.weight

    def remove_patient(self, pat):
        if isinstance(pat, str):
            pat = self.get_patients()[pat]
        pat.remove_module(self)
        pat_lst = self.get_patients_names_as_list()
        for omic in self.omics.values():
            edges = omic.graph.subgraph(pat_lst).edges(pat.get_name(), data=True)
            for e in edges:
                self.weight -= e[2]['weight']
        try:
            self.patients.pop(pat.get_name())
        except:
	    # This error should never occur!
            import pdb;pdb.set_trace()
            print("error in remove patient {} from module {}.".format(pat, self))
            exit(1)
        return self.weight

    def get_omics(self):
        return self.omics

    def add_weight(self, weight):
        self.weight += weight
        return self.weight

    def dec_weight(self, weight):
        self.weight -= weight
        return self.weight

    def get_weight(self):
        return self.weight

    # returns but does not change the object
    def calc_module_weight(self):
        weight = 0
        for omic in self.omics.values():
            g1 = omic.graph.subgraph(list(self.patients.keys()))
            weight += g1.size('weight')
        return weight

    def get_patients_names_as_list(self):
        return [name for name, node in self.get_patients().items()]

    def add_omic(self, omic, glob_var=None):
        if isinstance(omic, str):
            omic_name = omic
            omic = glob_var.omics[omic]
        else:  # omic is an omic object
            for name, obj in glob_var.omics.items():
                if obj == omic:
                    omic_name = name
                    break
        if omic_name not in self.omics:
            self.weight += omic.graph.subgraph(self.get_patients_names_as_list()).size("weight")
            self.omics[omic_name] = omic
        return self.weight

    def remove_omic(self, omic, glob_var=None):
        if isinstance(omic, str):
            omic_name = omic
            omic = self.omics[omic_name]
        else:
            for name, obj in glob_var.omics.items():
                if obj == omic:
                    omic_name = name
                    break
        if omic_name in self.omics:
            self.weight -= omic.graph.subgraph(self.get_patients_names_as_list()).size("weight")
            del self.omics[omic_name]
        return self.weight

    def merge_with_module_union(self, mod1, glob_var):
        for omic in mod1.get_omics().values():
            if omic not in self.omics:
                self.add_omic(omic, glob_var)
        for name in mod1.get_patients_names_as_list():
            mod1.remove_patient(glob_var.patients[name])
            self.add_patient(glob_var.patients[name])
        glob_var.kill_module(mod1)
        return glob_var

    def eat_module(self, mod, glob_var):
        lst = [(name, pat) for name, pat in mod.get_patients().items()]
        for (name, pat) in lst:
            mod.remove_patient(pat)
            self.add_patient(pat)
        glob_var.kill_module(mod)
        return glob_var

    def merge_me_into(self, mod, glob_var):
        mods_pats = [(name, pat) for name, pat in mod.get_patients().items()]
        selfs_omics_list = list(self.get_omics().items())
        # selfs_omics_list.sort(key=lambda x: x[0])
        for name, pat in mods_pats:
            mod.remove_patient(pat)
            self.add_patient(pat)
        for name, omic in selfs_omics_list:
            self.remove_omic(omic, glob_var)
        for name, omic in mod.get_omics().items():
            self.add_omic(omic, glob_var)
        glob_var.kill_module(mod)
        return glob_var

    def merge_to_intersection_omics(self, mod, glob_var):
        lst = [(name, pat) for name, pat in mod.get_patients().items()]
        for name, pat in lst:
            mod.remove_patient(pat)
            self.add_patient(pat)
        omics_lst = list(self.get_omics().items())
        omics_lst.sort(key=lambda x: x[0])
        for name, omic in omics_lst:
            if omic not in mod.get_omics().values():
                self.remove_omic(omic, glob_var)
        glob_var.kill_module(mod)
        return glob_var

    def split_module(self, sub_nodes, glob_var):
        pat_list = self.get_patients_names_as_list()
        left_out = set(sub_nodes).symmetric_difference(set(pat_list))
        new_mod = Module(glob_var)
        for omic in self.get_omics():
            new_mod.add_omic(omic, glob_var)
        for pat in left_out:
            self.remove_patient(glob_var.patients[pat])
            new_mod.add_patient(glob_var.patients[pat])
        return glob_var

    def split_and_add_omic(self, omic, sub_nodes, glob_var):
        pats = self.get_patients_names_as_list()
        left_out = set(pats).symmetric_difference(set(sub_nodes))
        if left_out:  # not empty
            new_mod = Module(glob_var=glob_var, omics=self.get_omics())
            for pat in left_out:
                self.remove_patient(pat)
                new_mod.add_patient(glob_var.patients[pat])
        self.add_omic(omic, glob_var)
        return glob_var
    
    def split_and_remove_omic(self, omic, sub_nodes, glob_var):
        pats = self.get_patients_names_as_list()
        pats_in_new_module = set(sub_nodes)
        pats_remain_in_this_module = set(pats).symmetric_difference(set(pats_in_new_module))
        if pats_in_new_module and (pats_remain_in_this_module or len(self.omics) > 1):  # not empty
            new_mod = Module(glob_var=glob_var, omics=[omic])
            for pat in pats_in_new_module:
                self.remove_patient(pat)
                new_mod.add_patient(glob_var.patients[pat])
        return glob_var

    def merge_with_module(self, params, glob_var):
        if params[1][0] == 0:
            return self.merge_with_module_union(params[1][1][1], glob_var)
        if params[1][0] == 1:
            return self.eat_module(params[1][1][1], glob_var)
        if params[1][0] == 2:
            return self.merge_me_into(params[1][1][1], glob_var)
        if params[1][0] == 3:
            return self.merge_to_intersection_omics(params[1][1][1], glob_var)


    def spread_module(self, params, glob_var):
        for pat, new_mod in params.items():
            self.remove_patient(pat)
            if new_mod is not None:
                new_mod.add_patient(pat)
        glob_var.kill_module(self)

