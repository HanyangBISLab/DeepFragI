"""
Feature extraction from peptides
"""

import re
import numpy as np

from . import aa_table as aa
from . import hydph_table as hydph
from . import mass_table as mass

FEATURE_SIZE = 96
MAX_PEPLEN = 30
num_ion_type = 8



def load_dataset(filepath):
    data_list = []

    data_count = 0
    with open(filepath) as f:
        while (True):
            line = f.readline()
            data_count += 1
            if not line: break

            split = line.split(" ")
            # fragment length is 19
            pep_seq = ""
            frag_len = int(split[-1])

            for i in range(frag_len):
                start = FEATURE_SIZE * i
                start_from = start + 46  # 46 is for other features
                end = start_from + 20

                # 40 x 1 is (20: left AA + 20: right AA)
                featurized_left_aa = split[start_from:end]
                #        print ("feature aa len {}".format(len(featurized_left_aa)))
                aminoacid = aa.get_aa_from_index(featurized_left_aa.index('1'))
                #        print (aminoacid)
                pep_seq += aminoacid

                if i == (frag_len - 1):
                    featurized_right_aa = split[start_from + 20:end + 20]
                    aminoacid = aa.get_aa_from_index(featurized_right_aa.index('1'))
                    #          print (aminoacid)
                    pep_seq += aminoacid

            feature = list(map(float, split[:-((MAX_PEPLEN-1) * num_ion_type + 1)]))  # 데이터에서 각 라인은 feature vector(19*92개),intensity값들(ion_type수*19개),펩타이드 서열길이 로 되어있으므로 19*num_ion_type에 +1을 해야한다.
            #print(len(feature))
            label = list(map(float, split[FEATURE_SIZE * (MAX_PEPLEN-1):-1]))
            #print(len(label))
            charge = split[FEATURE_SIZE-6:FEATURE_SIZE].index('1') + 1

            data_list.append(Data(pep_seq, charge, feature, label))  #class인 Data로 된 리스트

    return data_list


class Data:
    def __init__(self, pep_seq, charge, feature, label):
        self.pep_seq = pep_seq
        self.charge = charge
        self.feature = feature
        self.label = label


####################################################
## Featurize sequence to vector                   ##
####################################################
class Peptide:
    def __init__(self, peptide, charge, pick_location):
        self.peptide = self.get_strip_sequence(peptide)
        self.charge = charge
        self.pick_location = pick_location
        self.length = len(self.peptide)
        self.min = self.min_mass()
        self.max = self.max_mass()

    # get strip sequence
    def get_strip_sequence(self, peptide):
        p = re.compile("[^a-zA-Z]") #a-z는 제외, A-Z만 포함
        return p.sub('', peptide)  # peptide가 그대로 나옴


    # min(b1 ion, y1 ion)
    def min_mass(self):
        b = mass.get_aa_mass(self.peptide[0])
        if self.peptide[0] == 'C':
            b += 57.02146

        y = mass.get_aa_mass(self.peptide[self.length - 1])
        if self.peptide[self.length - 1] == 'C':
            y += 57.02146
        return min(b, y)

    # peptide의 length가 11일경우 max(b10 ion, y10 ion)을 리턴한다.
    # def max_mass(self):
    #     b = 0
    #     y = 0
    #     for i in range(self.length - 1):
    #         b += mass.get_aa_mass(self.peptide[i])
    #         if self.peptide[i] == 'C':
    #             b += 57.02146
    #
    #         y += mass.get_aa_mass(self.peptide[i + 1])
    #         if self.peptide[i + 1] == 'C':
    #             y += 57.02146
    #     return max(b, y)

    def max_mass(self):
        precursor = 0
        for i in range(self.length):
            precursor += mass.get_aa_mass(self.peptide[i])
            if self.peptide[i] == 'C':
                precursor += 57.02146
        return min(precursor, 2000)

    def get_features(self):
        ''' get featurized vector from peptide sequence

    Returns:
      92 dimensional vector

    Description:
      Peak Location Feature (1)
      Peak Composition Feature (40)
      Hydrophobicity (5)
        Sum of residue hydrophobicities (1)
        HYDC_1 (1)
        HYDN_1 (1)
        HYDN (1)
        HYDC (1)
      Amino acid sequence (40)
        left_spanning (20)
        right_sapnning (20)
     Charge state (6) one-hot

      Total 92 (1 + 40 + 5 + 40 + 6) dimensions vector

     Fragment site definition:

      peptide: ACDEDGGD
      peptide length: 8
      e.g.) A C D E D G G D
             | | | | | | |
      site   1 2 3 4 5 6 7

      To get the left side AA of fragmentation site //AA는 아미노산을 말한다.
      peptide[site - 1]
      To get the right side AA of fragmentation site
      peptide[site]
    '''
        features_vector = []

        for frag_site in range(1, self.length):
            features_vector.extend(self.get_peak_location_features(frag_site))

            features_vector.extend(self.get_composition_features(frag_site))

            features_vector.extend(self.get_hyd_features(frag_site, distance=1))

            features_vector.extend(self.get_spanning_aa_features(frag_site))

            left = 0
            right = 0
            if str(frag_site) in self.pick_location:
                left += 1
            if str(frag_site+1) in self.pick_location:
                right += 1
            features_vector.extend([left,right])

            total_num_modi = len(self.pick_location)
            left = 0
            for m in self.pick_location:
                m = int(m)
                if m <= frag_site:
                    left += 1
            right = total_num_modi - left
            features_vector.extend([left,right])

            # if frag_site == self.pick_location - 1:
            #     features_vector.extend([0,1])
            # elif frag_site == self.pick_location:
            #     features_vector.extend([1,0])
            # else:
            #     features_vector.extend([0,0])
            #
            # if frag_site < self.pick_location:
            #     features_vector.extend([0,1])
            # else:
            #     features_vector.extend([1,0])



            # 새로운 feature는 get_spanning_aa_features(frag_site) 와 get_charge_features(self.charge) 사이에 넣기
            features_vector.extend(self.get_charge_features(self.charge))
            # features_vector.extend(self.get_peptide_common_features())
        return features_vector

    """
  Peak location features
  """

    # (1 * 1)
    def get_peak_location_features(self, fragmentation_site):
        # [float, float, float]

        dist_from_min = self.min_peak_mass(fragmentation_site)
        dist_from_max = self.max_peak_mass(fragmentation_site)
        scaling_factor = self.max_mass() - self.min_mass()

        # Write the smaller one (the neariest to termini site)
        #return [dist_from_min / scaling_factor]
        if (dist_from_min < dist_from_max):
            return [dist_from_min / scaling_factor]
        else:
            return [dist_from_max / scaling_factor]

    # (1 * 1)
    def min_peak_mass(self, fragmentation_site):
        if fragmentation_site <= 0 or fragmentation_site >= self.length:
            return 0

        peak_mass = 0
        for i in range(fragmentation_site):
            # Fixed modification (Carbamidometyl) // proteometools 논문 참고
            if self.peptide[i] == 'C':
                peak_mass += 57.02146
            peak_mass += mass.get_aa_mass(self.peptide[i])  # mass.get_aa_mass(self.peptide[i]) => peptide[i]에 해당하는 아미노산의 mass를 리턴한다.

        # peak_mass = peptide의 N-term부터 fragmention site까지의 theoretical mass의 합.

        if peak_mass > 2000:
            peak_mass = 2000

        return peak_mass - self.min

    # (1 * 1)
    def max_peak_mass(self, fragmentation_site):
        if fragmentation_site <= 0 or fragmentation_site >= self.length:
            return 0

        peak_mass = 0
        for i in range(fragmentation_site):
            if self.peptide[i] == 'C':
                peak_mass += 57.02146
            peak_mass += mass.get_aa_mass(self.peptide[i])
        if peak_mass > 2000:
            peak_mass = 2000

        return self.max - peak_mass
    #

    # LEGACY FEATURE
    # (1 * 1)

    """
  Peptide composition features
  """

    # (40 * 1)
    def get_composition_features(self, fragmentation_site):
        # [] + []
        return self.n_x(fragmentation_site) + self.c_x(fragmentation_site)

    # number of x from n-term to fragmentation_site
    # (20 * 1)
    def n_x(self, fragmentation_site):
        if fragmentation_site <= 0 or fragmentation_site >= self.length:
            return 0

        vector = [0] * 20
        # amino acids before fragmentation site
        for i in range(fragmentation_site):
            vector[aa.get_index(self.peptide[i])] += 1
        return vector

    # number of x from i + 3 to fragmentation_site
    # (20 * 1)
    def c_x(self, fragmentation_site):
        if fragmentation_site <= 0 or fragmentation_site >= self.length:
            return 0

        vector = [0] * 20
        # amino acids after fragmentation site
        for i in range(fragmentation_site, self.length):
            vector[aa.get_index(self.peptide[i])] += 1
        return vector

    """
  Peptide Hydrophobicity Features
  """

    # (5 * 1)
    def get_hyd_features(self, fragmentation_site, distance):
        # float, float, float, float, float
        return [self.hydf(fragmentation_site),
                self.hydn_x(fragmentation_site, distance),
                self.hydc_x(fragmentation_site, distance),
                self.hydn_fragmentation_site(fragmentation_site),
                self.hydc_fragmentation_site(fragmentation_site)]

    # Hydrophobicity of fragment = sum of residue hydrophobicities
    # (1 * 1)
    def hydf(self, fragmentation_site):
        if fragmentation_site <= 0 or fragmentation_site >= self.length:
            return 0

        h_left = hydph.get_aa_hydph(self.peptide[fragmentation_site - 1])
        h_right = hydph.get_aa_hydph(self.peptide[fragmentation_site])
        return h_left + h_right

    # LEAGACY FEATURE
    # The average hydrophobicity of the amino acids
    # on both sides of fragmentation site
    # (1 * 1)
    def hydpra(self, fragmentation_site):
        if fragmentation_site <= 0 or fragmentation_site >= self.length:
            return 0

        return (hydph.get_aa_hydph(self.peptide[fragmentation_site])
                + hydph.get_aa_hydph(self.peptide[fragmentation_site - 1])) / 2

    # LEAGACY FEATURE
    # Differences in the hydrophobicity of amino acids
    # on both sides of fragmentation site
    # (1 * 1)
    def hydprd(self, fragmentation_site):
        if fragmentation_site <= 0 or fragmentation_site >= self.length:
            return 0

        return hydph.get_aa_hydph(self.peptide[fragmentation_site]) \
               - hydph.get_aa_hydph(self.peptide[fragmentation_site - 1])

    # Hydrophobicity of the amino acid at the x-distance
    # from the fragmentation site to the n-term
    # (1 * 1)
    def hydn_x(self, fragmentation_site, distance):
        if fragmentation_site - distance <= 0 \
                or fragmentation_site - distance >= self.length:
            return 0

        return hydph.get_aa_hydph(self.peptide[fragmentation_site - distance])

    # Hydrophobicity of the amino acid at the x-distance
    # from the fragmentation site to the c-term
    # (1 * 1)
    def hydc_x(self, fragmentation_site, distance):
        if fragmentation_site + distance <= 0 \
                or fragmentation_site + distance >= self.length:
            return 0

        return hydph.get_aa_hydph(self.peptide[fragmentation_site + distance])  # fragmentation_site + distance -1 로 해야하지 않나?

    # Hydrophobicity of the amino acid at the fragmentation site to n-term
    # (1 * 1)
    def hydn_fragmentation_site(self, fragmentation_site):
        if fragmentation_site <= 0 or fragmentation_site >= self.length:
            return 0

        return hydph.get_aa_hydph(self.peptide[fragmentation_site - 1])

    # Hydrophobicity of the amino acid at the fragmentation site to c-term
    # (1 * 1)
    def hydc_fragmentation_site(self, fragmentation_site):
        if fragmentation_site <= 0 or fragmentation_site >= self.length:
            return 0

        return hydph.get_aa_hydph(self.peptide[fragmentation_site])

    """
    Aminoa Acids feature
    (40 x 1) Get spanning AA vector
  """

    def get_spanning_aa_features(self, frag_site):
        left_vector = [0] * 20
        right_vector = [0] * 20

        if (frag_site < 1) or (frag_site > len(self.peptide) - 1):
            print("error" + str(frag_site))
            return 0

        left_vector[aa.get_index(self.peptide[frag_site - 1])] = 1
        right_vector[aa.get_index(self.peptide[frag_site])] = 1

        return left_vector + right_vector

    """
    Charge state feature
    (6 x 1)
  """

    def get_charge_features(self, charge):
        charge_vector = [0] * 6  # +1 to +6
        charge_vector[charge - 1] = 1

        return charge_vector

    """
  Peptide common features
  """

    # (41 + 21 * (peptide length)
    def get_peptide_common_features(self):
        # [] + [] + [float] + [] + []
        return self.nterm_is_x() + self.cterm_is_x() + [self.hydp()] \
               + self.sequence_info() + self.sequence_hydph_info()

    # n term is x
    # (20 * 1)
    def nterm_is_x(self):
        vector = [0] * 20
        vector[aa.get_index(self.peptide[0])] = 1
        return vector

    # c term is x
    # (20 * 1)
    def cterm_is_x(self):
        vector = [0] * 20
        vector[aa.get_index(self.peptide[self.length - 1])] = 1
        return vector

    # Hydph sum of peptide
    # (1 * 1)
    def hydp(self):
        s = 0
        for aa in self.peptide:
            s += hydph.get_aa_hydph(aa)
        return s

    # Intuitive sequence information
    # (20 * peptide length)
    def sequence_info(self):
        # peptide length limit is 11
        vector = []
        for i in range(self.length):
            v = [0] * 20
            v[aa.get_index(self.peptide[i])] = 1
            vector.extend(v)
        return vector

    # Intuitive peptide hydrophobicity information
    # (peptide length * 1)
    def sequence_hydph_info(self):
        vector = []
        for i in range(self.length):
            vector.append(hydph.get_aa_hydph(self.peptide[i]))
        return vector


if __name__ == "__main__":
    # Test
    print("TEST")
    pep = Peptide("ACCDDDEF", "2", "y")
    vec = pep.get_spanning_aa_features(1)
    print(vec)

    feature_vec = pep.get_features()

    print(feature_vec)

    print("peptide sequence: " + str("ACCDDDEF"))
    print("peptide length: " + str(len("ACCDDDEF")))
    print("feature vec length:" + str((len(feature_vec))))
    print("feature vec length / 66 : " + str(len(feature_vec) / 86))
