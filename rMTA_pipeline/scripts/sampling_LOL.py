import cobra
import pandas as pd
import numpy as np
from cobra.sampling import sample
import time

import os

os.chdir(os.path.dirname(__file__))
model = cobra.io.load_matlab_model('tissueModel_imat_late.mat')



# SAMPLING
print("Starting sampling at " + str(time.strftime('%d-%m-%Y__%H_%M')) )
print()
sampler = sample(model, 2000, method='achr')

sampler.to_csv(f'Sampler_ACHR_Recon3D_LOL.csv', index=False, )

print("Sampling done at " + str(time.strftime('%d-%m-%Y__%H_%M')) )
print()





