import pandas as pd

class QuantModel:
    def __init__(self):
        self.ratio = None

    # dname is the root directory of salmon output
    def from_file(self, dname):
        import os
        quant_name = os.path.sep.join([dname, 'quant.sf'])
        # Convert the tab separted data from quant.sf in to a pandas dataframe
        df = pd.read_csv(quant_name, sep='\t', header=(0))
        # Calculate the the list of all Effective to Actual length ratios
        self.ratio = (df['EffectiveLength']/df['Length']).tolist()
        return True