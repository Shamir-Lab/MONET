class Patient:
    """
    Represents a single patient (or any kind of sample in the data).
    """


    def __init__(self, name):
        self.module = None
        self.name = name

    def set_module(self, module):
        if self.is_in_module():
            print(self.name)
            raise NameError('trying to add patient to 2 modules', self)
        self.module = module
        return

    def remove_module(self, module):
        self.module = None
        return module

    def get_module(self):
        return self.module

    def is_in_module(self):
        return self.module is not None

    def __str__(self):
        return self.name

    def get_name(self):
        return self.name
