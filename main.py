import numpy as np
from scipy.optimize import minimize
import subprocess
import csv
import pandas as pd
#from interruptingcow import timeout
#function to be optimized
#order: T,R,a,b,r
def f(y):
    with open('constants.csv', mode='w') as file:
        writer = csv.writer(file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        writer.writerow(y)
    process = subprocess.Popen("java geodesic", shell=True, stdout=subprocess.PIPE)
    process.wait()
    data = pd.read_csv('data.csv')
    data = data.iloc[1:]['r']
    print(data.var())
    return data.var()
with open("optimal_parameters.csv","w+") as ans:
    writer = csv.writer(ans,delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    for search in range(10):
        try:
            #with timeout(60, exception = RuntimeError):
            print("Step {}".format(search))
            sol = minimize(f,np.random.random_sample(5)*100.0,method='Nelder-Mead')
            if(sol.success):
                writer.writerow(sol.x)
                print('Step was successful')
        except:
            print("Step wasn't successful")
            pass
    ans.close()
