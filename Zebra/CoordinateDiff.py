import numpy as np
#MyCal
Bead11 = np.array([743.2362,  45.6734, 153.5677])
Bead21= np.array([741.0815,  48.6720, 158.6730])
#GPT
Bead12 = np.array([743.2362, 45.6732, 153.5676])
Bead22 = np.array([741.0814, 48.6720, 158.6729])

#Gemini
Bead13 = np.array([743.2364, 45.6732, 153.5674])
Bead23 =  np.array([741.0813, 48.6720, 158.6730])
#Claude
Bead14 = np.array([743.4934, 45.6126, 153.6911])
Bead24 = np.array([741.1606, 48.7168, 158.6621])

print('GPT')
print(Bead11-Bead12)
print(Bead21-Bead22)
print('Gemini')
print(Bead11-Bead13)
print(Bead21-Bead23)
print('Claude')
print(Bead11-Bead14)
print(Bead11-Bead14)