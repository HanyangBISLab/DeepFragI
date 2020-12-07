
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import math

import numpy as np
import pandas as pd

from scipy.stats import pearsonr
from scipy.stats import spearmanr
from scipy.spatial.distance import cosine

import os
import sys
sys.path.append(os.path.dirname(os.path.abspath(os.path.dirname(__file__))))

from utils import peptide2
from utils import mass_TMT

from keras.models import load_model

import time
import datetime

FEATURE_SIZE = 96
AA_FEATURE_SIZE = 96
target_max_peplen = 30
FRAGMENT_SIZE = 29
ION_TYPE_SIZE = 8 # b and y ion, netural loss
iTRAQLabelMass = 229.163 #144.102
carbamidomethylMass = 57.021
phosphoMass = 79.966

def getSequence(seq, modiLocation):
    output = '+' + str(iTRAQLabelMass)
    i=1
    for c in seq:
        output += c
        if c=='C':
            output += '+' + str(carbamidomethylMass)
        elif c=='K':
            output += '+' + str(iTRAQLabelMass)
        elif str(i) in modiLocation:
            output += '+' + str(phosphoMass)
        i += 1
    return output



def zero_padding(features, max_frag_site_length):
    ''' zero padding the feature vector '''
    # get max feature vector length
    max_feature_vector_len = max_frag_site_length * AA_FEATURE_SIZE  # max_frag_site_length = 20 -1 = 19 이므로 max_feature_vector_len = 19 * 92 가 된다.

    padding_size = max_feature_vector_len - len(features)  # 펩타이드의 길이가 20보다 작을경우, 예를들면 길이가 15 일경우 features 벡터는 길이가 14*92 이므로 5*92만큼 제로패딩해준다.
    features.extend([0] * padding_size)




    return features



if __name__ == "__main__":

    modelList = os.listdir('../model')
    deep_frag = load_model('../model/'+modelList[0])

    i = open('uniprot_anno_psppept_list.txt')

    o = open('output.txt','w')
    t = 0

    start_time = time.time()
    while(True):
        line = i.readline()

        if not line: break
        t+=1

        split = line.split('\t')

        seq = split[0]
        charge = split[1]
        modiLocation = split[2][:-1].split(';')  ## \n은 빼기위해서 [:-1]로 한다.

        precursorMass = mass_TMT.getPrecursorMass(seq, modiLocation,int(charge))
        modiSeq = getSequence(seq,modiLocation)  ## C+57.~~ 처럼 modification의 mass가 같이 적혀있는 sequence

        # match = split[3].split(';')
        # inten = match.split(';')
        #
        #
        #
        # realValues = []
        # for i in len(match):
        #     realValues.append([match[i] , inten[i]])
        #
        # realValues = sorted(realValues , key = lambda x : x[0])



        modiLocationFirst = int(modiLocation[0])

        pept = peptide2.Peptide(seq, int(charge), modiLocation)

        feature = pept.get_features()

        feature = zero_padding(feature, target_max_peplen - 1)  # target_max_peplen = 20

        X = np.array(feature)
        # reshape to (batch, time, feature)
        X = X.reshape(1, FRAGMENT_SIZE, FEATURE_SIZE)

        predicted = [[0 for x in range(ION_TYPE_SIZE)] for j in range(len(seq) - 1)] ##['b(2+)P', 'b(2+)P-loss', 'bP', 'bP-loss', 'y(2+)P', 'y(2+)P-loss', 'yP', 'yP-loss']
        predicted = np.array(predicted)  #print(predicted)  # [[b1,b1-H2O,b1-NH3,y3,y3-H2O,y3-NH3],[b2,b2-H2O,b2-NH3,y2,y2-H2O,y2-NH3],[b3,b3-H2O,b3-NH3,y1,y1-H2O,y1-NH3],]

        predicted = deep_frag.predict(X)[0, :len(seq)-1, :]
        predicted[predicted < 0] = 0

        massPeakPair = []


        b_modi = mass_TMT.bion_phospho_mass_list(seq, modiLocation)
        #b_modi_iso = mass_TMT.bion_phospho_mass_isotope_list(seq, modiLocation)

        y_modi = mass_TMT.yion_phospho_mass_list(seq, modiLocation)
        y_modi.reverse()
        #y_modi_iso = mass_TMT.yion_phospho_mass_isotope_list(seq, modiLocation)
        #y_modi_iso.reverse()

        b_modi_loss = mass_TMT.bion_phospho_loss_mass_list(seq, modiLocation)
        #b_modi_loss_iso = mass_TMT.bion_phospho_loss_mass_isotope_list(seq, modiLocation)

        y_modi_loss = mass_TMT.yion_phospho_loss_mass_list(seq, modiLocation)
        y_modi_loss.reverse()
        #y_modi_loss_iso = mass_TMT.yion_phospho_loss_mass_isotope_list(seq, modiLocation)
        #y_modi_loss_iso.reverse()

        b2_modi = mass_TMT.bion2_phospho_mass_list(seq, modiLocation)
        #b2_modi_iso = mass_TMT.bion2_phospho_mass_isotope_list(seq, modiLocation)

        y2_modi = mass_TMT.yion2_phospho_mass_list(seq, modiLocation)
        y2_modi.reverse()
        #y2_modi_iso = mass_TMT.yion2_phospho_mass_isotope_list(seq, modiLocation)
        #y2_modi_iso.reverse()

        b2_modi_loss = mass_TMT.bion2_phospho_loss_mass_list(seq, modiLocation)
        #b2_modi_loss_iso = mass_TMT.bion2_phospho_loss_mass_isotope_list(seq, modiLocation)

        y2_modi_loss = mass_TMT.yion2_phospho_loss_mass_list(seq, modiLocation)
        y2_modi_loss.reverse()
        #y2_modi_loss_iso = mass_TMT.yion2_phospho_loss_mass_isotope_list(seq, modiLocation)
        #y2_modi_loss_iso.reverse()

        modiLocationReverse = len(seq) - modiLocationFirst + 1

        for j in range(len(seq) - modiLocationReverse):
            massPeakPair.append([round(y_modi_loss[j],7),round(predicted[j][7],7)])  ##y(n-1)부터 채운다.
            #massPeakPair.append([round(y_modi_loss_iso[j], 7), round((5.42e-4 * y_modi_loss_iso[j] + 1e-30)*predicted[j][7], 7)])

            massPeakPair.append([round(y2_modi_loss[j],7),round(predicted[j][5],7)])
            #massPeakPair.append([round(y2_modi_loss_iso[j],7),round((5.42e-4 * y2_modi_loss_iso[j] + 1e-30)*predicted[j][5],7)])




        k = 0
        for j in range(modiLocationFirst -1, len(seq)-1):
            massPeakPair.append([round(b_modi_loss[k],7) , round(predicted[j][3],7)])
            #massPeakPair.append([round(b_modi_loss_iso[k],7) , round((5.42e-4 * b_modi_loss_iso[k] + 1e-30)*predicted[j][3],7)])

            massPeakPair.append([round(b2_modi_loss[k],7) , round(predicted[j][1],7)])
            #massPeakPair.append([round(b2_modi_loss_iso[k],7) , round((5.42e-4 * b2_modi_loss_iso[k] + 1e-30)*predicted[j][1],7)])
            k+=1

        for j in range(len(seq) - 1):
            massPeakPair.append([round(b_modi[j],7), round(predicted[j][2],7)])
            #massPeakPair.append([round(b_modi_iso[j],7), round((5.42e-4 * b_modi_iso[j] + 1e-30)*predicted[j][2],7)])

            massPeakPair.append([round(b2_modi[j],7), round(predicted[j][0],7)])
            #massPeakPair.append([round(b2_modi_iso[j],7), round((5.42e-4 * b2_modi_iso[j] + 1e-30)*predicted[j][0],7)])

            massPeakPair.append([round(y_modi[j],7), round(predicted[j][6],7)])
            #massPeakPair.append([round(y_modi_iso[j],7), round((5.42e-4 * y_modi_iso[j] + 1e-30)*predicted[j][6],7)])

            massPeakPair.append([round(y2_modi[j],7), round(predicted[j][4],7)])
            #massPeakPair.append([round(y2_modi[j],7), round((5.42e-4 * y2_modi_iso[j] + 1e-30)*predicted[j][4],7)])

        massPeakPair = sorted(massPeakPair , key = lambda x : x[0])

        o.write('BEGIN IONS\n')
        o.write('TITLE=' + str(t) + '\n')
        o.write('SEQ=' + modiSeq + '\n')
        o.write('PEPMASS=' + str(precursorMass) + '\n')
        o.write('CHARGE=' + str(charge) + '+\n')


        for massPeak in massPeakPair:
            if massPeak[1] > 0:
                o.write(str(massPeak[0]) + ' ' + str(massPeak[1]) + '\n')

        o.write('END IONS\n')

        if t%1000 == 0:
            print(t)

    i.close()
    o.close()

    end_time = time.time()

    print(str(datetime.datetime.fromtimestamp(start_time).strftime('%Y-%m-%d_%H:%M:%S')))
    print(str(datetime.datetime.fromtimestamp(end_time).strftime('%Y-%m-%d_%H:%M:%S')))






































