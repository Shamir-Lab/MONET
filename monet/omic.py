class Omic:
    """
    Represents a single omic, and contains the omic's graph.
    """

    def __init__(self, graph=None, name=""):
        self.graph = graph
        self.name = name
        return

    def set_graph(self, g):
        self.graph = g
        return

    def get_name(self):
        return self.name

    def __repr__(self):
        return self.name

