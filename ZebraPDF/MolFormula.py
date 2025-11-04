
CStart = 8
HStart = 15
NStart = 1
OStart = 2

nmonomer = 50

Ctotal = CStart*nmonomer
Htotal = HStart*nmonomer
Ntotal = NStart*nmonomer
Ototal = OStart*nmonomer

Ctotal+=2
Htotal+=6

print(f'C{Ctotal}H{Htotal}N{Ntotal}O{Ototal}')


file = open('./ZebraPDF/CoilPDB.pdb','r')

Ctot=0
Ntot=0
Otot=0
Htot=0
for line in file:
    if not line.startswith('ATOM'):
        continue
    else:
        lsp = line.split()
        atom = lsp[-1]
        if atom =='C':
            Ctot+=1
        if atom =='N':
            Ntot+=1
        if atom =='O':
            Otot+=1
        if atom =='H':
            Htot+=1

print('##########################################')
print('Ctotal = ', Ctot)
print('Htotal = ', Htot)
print('Ntotal = ', Ntot)
print('Ototal = ', Otot)
file.close()
print('##########################################')

file = open('./ZebraPDF/Monomer.pdb','r')

Ctot=0
Ntot=0
Otot=0
Htot=0
Formula = ''
for line in file:
    if not line.startswith('ATOM'):
        continue
    else:
        lsp = line.split()
        atom = lsp[-1]
        if atom =='C':
            Ctot+=1
            Formula+=atom
        if atom =='N':
            Ntot+=1
            Formula+=atom
        if atom =='O':
            Otot+=1
            Formula+=atom
        if atom =='H':
            Htot+=1
            Formula+=atom
            
print('Ctotal = ', Ctot)
print('Htotal = ', Htot)
print('Ntotal = ', Ntot)
print('Ototal = ', Otot)
print(Formula)