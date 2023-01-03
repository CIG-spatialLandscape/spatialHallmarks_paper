import multiprocessing as mp
import pandas as pd
import random
import numpy as np
import time

coord = pd.read_csv("../../matrix_coordinates_4A.csv", index_col=0)
coord = coord.transpose()
prediction = pd.read_csv("../../matrix_4A_prediction.csv", index_col=0)
#sub = random.sample(range(0,2772), 200)
#prediction = prediction.iloc[:,sub]
#prediction = prediction.iloc[[2,4,6],:]


def compute(spot1):
    score = np.array([0.0,0.0,0.0,0.0,0.0])
    score[0] = float(prediction[[spot1]].iloc[4]) #Epithelial - CAF
    score[1] = float(prediction[[spot1]].iloc[6]) #Malignant - CAF
    score[2] = float(prediction[[spot1]].iloc[4]) #Epithelial - T cells
    score[3] = float(prediction[[spot1]].iloc[6]) #Malignant - T cells
    score[4] = float(prediction[[spot1]].iloc[2]) #CAFs - T cells


    value = np.array([0.0,0.0,0.0,0.0,0.0])
    for spot2 in prediction:
        CAF = float(prediction[[spot2]].iloc[2])
        Tcell = float(prediction[[spot2]].iloc[7])
        d = np.sqrt((float(coord[[spot2]].iloc[3])-float(coord[[spot1]].iloc[3]))**2+(float(coord[[spot2]].iloc[4])-float(coord[[spot1]].iloc[4]))**2)
        if d == 0:
            d = 1
        value[0] += CAF/d  # Epithelial - CAF  
        value[1] += CAF/d    # Malignant - CAF
        value[2] += Tcell/d   # Epithelial - T cells
        value[3] += Tcell/d    # Malignant -  T cells
        value[4] += Tcell/d    # CAFs -  T cells
    for i in range(0,5):
        score[i] = score[i]*value[i]
    return(score)


import time
start_time = time.time()
pool= mp.Pool(mp.cpu_count())

from operator import add
print(sum(pool.map(compute, list(prediction.columns))))

pool.close()


print("--- %s seconds ---" % (time.time() - start_time))
#print(float(prediction[['TAAGTAAATGTGCCGC-1']].iloc[2]))
#print(compute('TAAGTAAATGTGCCGC-1'))