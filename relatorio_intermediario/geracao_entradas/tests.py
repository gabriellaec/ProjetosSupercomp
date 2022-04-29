# -*- coding: utf-8 -*-
import random

i=10
index=0
while i < 3000:
     n = i # tamanho da primeira sequência
     m = i # tamanho da segunda sequência
     file = "../in_range_maior/dna{}.seq".format(index) # nome do arquivo a ser gerado
     f = open(file, 'w')
     seq=[str(n)+'\n',
          str(m)+'\n',
          ''.join(random.choices(['A','T','C','G','-'],k=n))+'\n',
          ''.join(random.choices(['A','T','C','G','-'],k=m))]
     f.writelines(seq)
     f.close()
     print(''.join(seq))
     i+=100
     index+=1
