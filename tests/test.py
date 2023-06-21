from src.orient_index import orient_index
import numpy as np

H = np.arange(-5.1,5.1, 0.01)
K = np.arange(-5.1,5.1, 0.01)
L = np.arange(-9.1,9.1, 0.02)

orient_index(H, K, L, data_dir='K:/pokharel-3470-a/PrCd3P3/BRO7/', orm_temp='300', index_temps=None,
             stacks_list=[1, 2, 3])