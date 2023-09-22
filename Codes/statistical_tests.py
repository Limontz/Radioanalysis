# A class made of multiple statistical tests
# commonly use to compare two samples 
# LIMITED TO 2 SAMPLES ONLY

import numpy as np
from scipy import stats as sts
import math


class StatisticalTests:
    def __init__(self, data1, data2):
        self.data1 = data1
        self.data2 = data2

    def mann_whitney_u_test(self, method='asymptotic'): #choose method="exact" for smaller samples
        
        """ Perform the mann whitney U test to check 
            whether two sample come from the 
            same distribution 

            Params
            method: string
                  Specify which method you want to use. 
                  Common choices are "asymptotic" and "exact"
            output : 3 scalar
                  Return the U statistic, the z statistic and the p-value

            References
            https://datatab.net/tutorial/mann-whitney-u-test
        """

        n1 = len(self.data1)
        n2 = len(self.data2)
        U, p_value = sts.mannwhitneyu(self.data1, self.data2, method=method)
        expected_U = n1 * n2 * 0.5
        simga_U = np.sqrt( n1*n2*(n1+n2+1)/12 )
        z = (U-expected_U)/simga_U
        return U, z, p_value
    
    def two_sample_ks_test(self):

         """Perform the KS test to check 
            whether two sample come from the 
            same distribution 

            Params
            output : 2 scalar
                  Return the D statistic and the p-value

            References
            https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.ks_2samp.html
        """

         D, p_value = sts.ks_2samp(self.data1, self.data2)
         return D, p_value
    
    def two_sample_t_test(self):
        t_statistic, p_value = sts.ttest_ind(self.data1, self.data2)
        return t_statistic, p_value
    

class CircularStatistics:
    def __init__(self, data):
        self.data = data

    def Rayleigh_test(self):
        
        """Perform the Rayleigh test to check 
           the uniformity of circular data

            Params
            input : np.ndarray
                    Array of angles in degrees 
            output : 2 scalar
                  Return the statistic and the p-value

            References
            http://palaeo.spb.ru/pmlibrary/pmbooks/mardia&jupp_2000.pdf (pag.94)
        """


        n = len(self.data)
        r = np.sqrt(((1/n) * np.sum(np.cos(2 * self.data * math.pi / 180)))**2 + ((1/n) * np.sum(np.sin(2 * self.data * math.pi / 180)))**2)
        statistic = (1 - (2/n)) * 2 * n * r**2 + (n * r**4) / 2
        p_value = 1 - sts.chi2.cdf(x=statistic, df=2, loc=0, scale=1)
        return statistic, p_value