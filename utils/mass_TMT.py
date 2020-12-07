LabelMass = 229.162932
proton = 1.007276
isotope = 1.00235
h2o = 18.01528
phospho = 79.966
h2oProton = h2o + proton

def mass(ch):
    return aatable[ch]

def bion_mass(ch):
    return aatable[ch] + 1.0078

def yion_mass(ch):
    return aatable[ch] + 19.0184

def aion_mass(ch):
    return aatable[ch] - 26.9871

def xion_mass(ch):
    return aatable[ch] + 44.9977

def zion_mass(ch):
    return aatable[ch] + 1.9918

def getPrecursorMass(seq, modiLocation, charge):
    mass = LabelMass
    numModi = len(modiLocation)
    for c in seq:
        mass += aatable[c]
        if c=='C':
            mass += 57.02146
    mass = (mass + 79.966*numModi + 18.01528 + 1.0078*charge)/charge

    return mass

def bion_phospho_mass_list(seqeunce ,modi_location):
    masses = list()
    pre_mass = LabelMass

    for i in range(len(seqeunce)-1):
        ch = seqeunce[i]
        if ch == 'C':
            pre_mass += 57.02146
        elif str(i+1) in modi_location:
            pre_mass += 79.966
        pre_mass += aatable[ch]
        masses.append(pre_mass + 1.007276)


    return masses

def bion_phospho_mass_isotope_list(seqeunce ,modi_location):
    masses = list()
    pre_mass = LabelMass

    for i in range(len(seqeunce)-1):
        ch = seqeunce[i]
        if ch == 'C':
            pre_mass += 57.02146
        elif str(i+1) in modi_location:
            pre_mass += 79.966
        pre_mass += aatable[ch]
        masses.append(pre_mass + 1.007276 + isotope )


    return masses


def bion2_phospho_mass_list(seqeunce ,modi_location):
    masses = list()
    pre_mass = LabelMass

    for i in range(len(seqeunce)-1):
        ch = seqeunce[i]
        if ch == 'C':
            pre_mass += 57.02146
        elif str(i+1) in modi_location:
            pre_mass += 79.966
        pre_mass += aatable[ch]
        masses.append((pre_mass + 2*proton)/2)



    return masses

def bion2_phospho_mass_isotope_list(seqeunce ,modi_location):
    masses = list()
    pre_mass = LabelMass

    for i in range(len(seqeunce)-1):
        ch = seqeunce[i]
        if ch == 'C':
            pre_mass += 57.02146
        elif str(i+1) in modi_location:
            pre_mass += 79.966
        pre_mass += aatable[ch]
        masses.append((pre_mass + 2*proton + isotope )/2)

    return masses

def bion_phospho_loss_mass_list(seqeunce ,modi_location):
    masses = list()
    pre_mass = LabelMass
    modi_location_first = int(modi_location[0])

    for i in range(modi_location_first - 1):
        ch = seqeunce[i]
        if ch == 'C':
            pre_mass += 57.02146
        elif str(i+1) in modi_location:
            pre_mass += 79.966
        pre_mass += aatable[ch]
        # masses.append(pre_mass + 1.0078)###


    pre_mass = pre_mass - 18.01528 - phospho

    for i in range(modi_location_first - 1,len(seqeunce)-1):
        ch = seqeunce[i]
        if ch == 'C':
            pre_mass += 57.02146
        elif str(i+1) in modi_location:
            pre_mass += 79.966
        pre_mass += aatable[ch]
        masses.append(pre_mass + proton)  #############################################

    return masses

def bion_phospho_loss_mass_isotope_list(seqeunce ,modi_location):
    masses = list()
    pre_mass = LabelMass
    modi_location_first = int(modi_location[0])

    for i in range(modi_location_first - 1):
        ch = seqeunce[i]
        if ch == 'C':
            pre_mass += 57.02146
        elif str(i+1) in modi_location:
            pre_mass += 79.966
        pre_mass += aatable[ch]
        # masses.append(pre_mass + 1.0078)###


    pre_mass = pre_mass - 18.01528 - phospho

    for i in range(modi_location_first - 1,len(seqeunce)-1):
        ch = seqeunce[i]
        if ch == 'C':
            pre_mass += 57.02146
        elif str(i+1) in modi_location:
            pre_mass += 79.966
        pre_mass += aatable[ch]
        masses.append(pre_mass + proton + isotope )  #############################################

    return masses

def bion2_phospho_loss_mass_list(seqeunce ,modi_location):
    masses = list()
    pre_mass = LabelMass
    modi_location_first = int(modi_location[0])

    for i in range(modi_location_first- 1):
        ch = seqeunce[i]
        if ch == 'C':
            pre_mass += 57.02146
        pre_mass += aatable[ch]


    pre_mass = pre_mass - 18.01528 - phospho

    for i in range(modi_location_first - 1,len(seqeunce)-1):
        ch = seqeunce[i]
        if ch == 'C':
            pre_mass += 57.02146
        elif str(i+1) in modi_location:
            pre_mass += 79.966
        pre_mass += aatable[ch]
        masses.append((pre_mass + 2*proton)/2)

    return masses

def bion2_phospho_loss_mass_isotope_list(seqeunce ,modi_location):
    masses = list()
    pre_mass = LabelMass
    modi_location_first = int(modi_location[0])

    for i in range(modi_location_first- 1):
        ch = seqeunce[i]
        if ch == 'C':
            pre_mass += 57.02146
        pre_mass += aatable[ch]


    pre_mass = pre_mass - 18.01528 - phospho

    for i in range(modi_location_first - 1,len(seqeunce)-1):
        ch = seqeunce[i]
        if ch == 'C':
            pre_mass += 57.02146
        elif str(i+1) in modi_location:
            pre_mass += 79.966
        pre_mass += aatable[ch]
        masses.append((pre_mass + 2*proton + isotope )/2)

    return masses

def bion_mass_list(seqeunce):
    masses = list()
    pre_mass = 0
    for i in range(len(seqeunce) - 1):
        ch = seqeunce[i]
        if ch == 'C':
            pre_mass += 57.02146
        pre_mass += aatable[ch]
        masses.append(pre_mass + proton)
    return masses

def bion_H2O_mass_list(seqeunce):
    masses = list()
    pre_mass = -18.01528
    for i in range(len(seqeunce) - 1):
        ch = seqeunce[i]
        if ch == 'C':
            pre_mass += 57.02146
        pre_mass += aatable[ch]
        masses.append(pre_mass + proton)
    return masses

def bion_NH3_mass_list(seqeunce):
    masses = list()
    pre_mass = -17.031
    for i in range(len(seqeunce) - 1):
        ch = seqeunce[i]
        if ch == 'C':
            pre_mass += 57.02146
        pre_mass += aatable[ch]
        masses.append(pre_mass + proton)
    return masses

def bion2_mass_list(seqeunce): ## charge 2+ b fragmentation ion mass list
    masses = list()
    pre_mass = 0
    for i in range(len(seqeunce) - 1):
        ch = seqeunce[i]
        if ch == 'C':
            pre_mass += 57.02146
        pre_mass += aatable[ch]
        masses.append((pre_mass + 2*proton)/2)
    return masses

def yion_phospho_mass_list(seqeunce,modi_location):
    masses = list()
    pre_mass = 0
    modi_location_first = int(modi_location[0])
    seqLen = len(seqeunce)
    modi_location2 = []
    for m in modi_location:
        modi_location2.append(str(seqLen-int(m)+1))

    # modi_location_reverse = len(seqeunce) - modi_location + 1

    for i in range(len(seqeunce)-1):
        ch = seqeunce[::-1][i] # Reverse order
        if ch == 'C':
            pre_mass += 57.02146
        elif str(i+1) in modi_location2:
            pre_mass += 79.966
        pre_mass += aatable[ch]
        masses.append(pre_mass + h2oProton)


    return masses

def yion_phospho_mass_isotope_list(seqeunce,modi_location):
    masses = list()
    pre_mass = 0
    modi_location_first = int(modi_location[0])
    seqLen = len(seqeunce)
    modi_location2 = []
    for m in modi_location:
        modi_location2.append(str(seqLen-int(m)+1))

    # modi_location_reverse = len(seqeunce) - modi_location + 1

    for i in range(len(seqeunce)-1):
        ch = seqeunce[::-1][i] # Reverse order
        if ch == 'C':
            pre_mass += 57.02146
        elif str(i+1) in modi_location2:
            pre_mass += 79.966
        pre_mass += aatable[ch]
        masses.append(pre_mass + h2oProton + isotope )


    return masses

def yion2_phospho_mass_list(seqeunce,modi_location):
    masses = list()
    pre_mass = 0
    # modi_location_reverse = len(seqeunce) - modi_location + 1

    seqLen = len(seqeunce)
    modi_location2 = []
    for m in modi_location:
        modi_location2.append(str(seqLen-int(m)+1))

    for i in range(seqLen - 1):
        ch = seqeunce[::-1][i] # Reverse order
        if ch == 'C':
            pre_mass += 57.02146
        elif str(i+1) in modi_location2:
            pre_mass += 79.966
        pre_mass += aatable[ch]
        masses.append((pre_mass + h2oProton + proton)/2)



    return masses

def yion2_phospho_mass_isotope_list(seqeunce,modi_location):
    masses = list()
    pre_mass = 0
    # modi_location_reverse = len(seqeunce) - modi_location + 1

    seqLen = len(seqeunce)
    modi_location2 = []
    for m in modi_location:
        modi_location2.append(str(seqLen-int(m)+1))

    for i in range(seqLen - 1):
        ch = seqeunce[::-1][i] # Reverse order
        if ch == 'C':
            pre_mass += 57.02146
        elif str(i+1) in modi_location2:
            pre_mass += 79.966
        pre_mass += aatable[ch]
        masses.append((pre_mass + h2oProton + proton + isotope )/2)



    return masses


def yion_phospho_loss_mass_list(seqeunce, modi_location):
    masses = list()
    pre_mass = 0
    modi_location_reverse = len(seqeunce) - int(modi_location[0]) + 1

    seqLen = len(seqeunce)
    modi_location2 = []
    for m in modi_location:
        modi_location2.append(str(seqLen-int(m)+1))

    for i in range(modi_location_reverse - 1):
        ch = seqeunce[::-1][i]  # Reverse order
        if ch == 'C':
            pre_mass += 57.02146
        elif str(i+1) in modi_location2:
            pre_mass += 79.966
        pre_mass += aatable[ch]
        # masses.append(pre_mass + 19.0184) ##


    pre_mass -= h2o

    for i in range(modi_location_reverse - 1, len(seqeunce) - 1):
        ch = seqeunce[::-1][i]  # Reverse order
        if ch == 'C':
            pre_mass += 57.02146
        pre_mass += aatable[ch]
        masses.append(pre_mass + h2oProton)

    return masses

def yion_phospho_loss_mass_isotope_list(seqeunce, modi_location):
    masses = list()
    pre_mass = 0
    modi_location_reverse = len(seqeunce) - int(modi_location[0]) + 1

    seqLen = len(seqeunce)
    modi_location2 = []
    for m in modi_location:
        modi_location2.append(str(seqLen-int(m)+1))

    for i in range(modi_location_reverse - 1):
        ch = seqeunce[::-1][i]  # Reverse order
        if ch == 'C':
            pre_mass += 57.02146
        elif str(i+1) in modi_location2:
            pre_mass += 79.966
        pre_mass += aatable[ch]
        # masses.append(pre_mass + 19.0184) ##


    pre_mass -= h2o

    for i in range(modi_location_reverse - 1, len(seqeunce) - 1):
        ch = seqeunce[::-1][i]  # Reverse order
        if ch == 'C':
            pre_mass += 57.02146
        pre_mass += aatable[ch]
        masses.append(pre_mass + h2oProton + isotope )

    return masses

def yion2_phospho_loss_mass_list(seqeunce, modi_location):
    masses = list()
    pre_mass = 0

    modi_location_reverse = len(seqeunce) - int(modi_location[0]) + 1

    seqLen = len(seqeunce)
    modi_location2 = []
    for m in modi_location:
        modi_location2.append(str(seqLen-int(m)+1))

    for i in range(modi_location_reverse - 1):
        ch = seqeunce[::-1][i]  # Reverse order
        if ch == 'C':
            pre_mass += 57.02146
        elif str(i+1) in modi_location2:
            pre_mass += 79.966
        pre_mass += aatable[ch]


    pre_mass -= h2o

    for i in range(modi_location_reverse - 1, len(seqeunce) - 1):
        ch = seqeunce[::-1][i]  # Reverse order
        if ch == 'C':
            pre_mass += 57.02146
        pre_mass += aatable[ch]
        masses.append((pre_mass + h2oProton + proton)/2)

    return masses

def yion2_phospho_loss_mass_isotope_list(seqeunce, modi_location):
    masses = list()
    pre_mass = 0

    modi_location_reverse = len(seqeunce) - int(modi_location[0]) + 1

    seqLen = len(seqeunce)
    modi_location2 = []
    for m in modi_location:
        modi_location2.append(str(seqLen-int(m)+1))

    for i in range(modi_location_reverse - 1):
        ch = seqeunce[::-1][i]  # Reverse order
        if ch == 'C':
            pre_mass += 57.02146
        elif str(i+1) in modi_location2:
            pre_mass += 79.966
        pre_mass += aatable[ch]


    pre_mass -= h2o

    for i in range(modi_location_reverse - 1, len(seqeunce) - 1):
        ch = seqeunce[::-1][i]  # Reverse order
        if ch == 'C':
            pre_mass += 57.02146
        pre_mass += aatable[ch]
        masses.append((pre_mass + h2oProton + proton + isotope )/2)

    return masses

def yion_mass_list(seqeunce):
    masses = list()
    pre_mass = 0
    for i in range(len(seqeunce) - 1):
        ch = seqeunce[::-1][i] # Reverse order
        if ch == 'C':
            pre_mass += 57.02146
        pre_mass += aatable[ch]
        masses.append(pre_mass + h2oProton)
    return masses

def yion_H2O_mass_list(seqeunce):
    masses = list()
    pre_mass = -18.01528
    for i in range(len(seqeunce) - 1):
        ch = seqeunce[::-1][i] # Reverse order
        if ch == 'C':
            pre_mass += 57.02146
        pre_mass += aatable[ch]
        masses.append(pre_mass + h2oProton)
    return masses

def yion_NH3_mass_list(seqeunce):
    masses = list()
    pre_mass = -17.031
    for i in range(len(seqeunce) - 1):
        ch = seqeunce[::-1][i] # Reverse order
        if ch == 'C':
            pre_mass += 57.02146
        pre_mass += aatable[ch]
        masses.append(pre_mass + h2oProton)
    return masses

def yion2_mass_list(seqeunce):
    masses = list()
    pre_mass = 0
    for i in range(len(seqeunce) - 1):
        ch = seqeunce[::-1][i] # Reverse order
        if ch == 'C':
            pre_mass += 57.02146
        pre_mass += aatable[ch]
        masses.append((pre_mass + h2oProton + proton)/2)
    return masses

aatable = {
    'A': 71.037114,
    # 'B'
    'C': 103.00919,
    'D': 115.02694,
    'E': 129.04259,
    'F': 147.06841,
    'G': 57.021464,
    'H': 137.05891,
    'I': 113.08406,
    # 'J'
    'K': 357.257892, ## K의 mass에 TMT 레이블의 mass를 더해서 357.257892가 된다.
    'L': 113.08406,
    'M': 131.04048,
    'N': 114.04293,
    # 'O'
    'P': 97.052764,
    'Q': 128.05858,
    'R': 156.10111,
    'S': 87.032029,
    'T': 101.04768,
    'U': 151.00919,
    'V': 99.068414,
    'W': 186.07931,
    # 'X'
    'Y': 163.06333,
}

# a = 'AAEAAGGKYRSTVSKSKD'
# b = ['9','11']
# print(bion_phospho_mass_list(a,b))
# print(yion_phospho_mass_list(a,b))
# print(bion_phospho_loss_mass_list(a,b))
# print(yion_phospho_loss_mass_list(a,b))
# print(yion2_phospho_loss_mass_list(a,b))
# print(bion2_phospho_loss_mass_list(a,b))
# print(yion2_phospho_mass_list(a,b))
# print(bion2_phospho_mass_list(a,b))
#
# a = aatable['A'] + aatable['A'] + aatable['E'] + aatable['A']+ aatable['A'] +aatable['G'] + aatable['G'] + aatable['K'] + aatable['Y'] -18.02 + 229.163
# b = a + aatable['R'] + aatable['S'] + 79.966
#
# print(a)
# print(b)
