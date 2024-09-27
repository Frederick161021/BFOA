from fastaReader import fastaReader
import random
import numpy
import copy
from evaluadorBlosum import evaluadorBlosum

class bacteria():


    def __init__(self, path):
        self.matrix = fastaReader(path)
        self.blosumScore = 0
        self.fitness = 0
        self.interaction =0
        self.NFE = 0

    def showGenome(self):
     for seq in self.matrix.seqs:
        print(seq)

    def clonar(self, path):
        newBacteria = bacteria(path)
        newBacteria.matrix.seqs = numpy.array(copy.deepcopy(self.matrix.seqs))
        return newBacteria

    def tumboNado(self, numGaps):

        self.cuadra()
        matrixCopy = copy.deepcopy(self.matrix.seqs)
        """convierto a lista para poder modificar"""
        matrixCopy = matrixCopy.tolist()
        gapRandomNumber = random.randint(0,numGaps)  #numero de gaps a insertar
        for i in range(gapRandomNumber):                    #cilco de gaps
            seqnum = random.randint(0, len(matrixCopy)-1)   #selecciono secuencia
            pos = random.randint(0, len(matrixCopy[0]))     #determina de forma alatoria la posicion del gap entre un numero de 0 a la longitud de la fila
            part1 = matrixCopy[seqnum][:pos]    #divide la fila sin incluir el indice pos
            part2 = matrixCopy[seqnum][pos:]    #divide la fila incluyendo el indice pos
            temp = "-".join([part1, part2])     #inserto gap
            matrixCopy[seqnum] = temp
        matrixCopy = numpy.array(matrixCopy)   #convierto a numpy array de regreso para fijar tama�os
        self.matrix.seqs = matrixCopy

        self.cuadra2()
        self.limpiaColumnas()

    def cuadra(self):
        """rellena con gaps las secuencias mas cortas"""
        import numpy
        seq = self.matrix.seqs
        maxLen = len(max(seq, key=len))
        for i in range(len(seq)):
            if len(seq[i]) < maxLen:
                seq[i] = seq[i] + "-"*(maxLen-len(seq[i]))
        self.matrix.seqs = numpy.array(seq)

    #Este meotodo sera modificado con el objetivo de no tener una enorme cantidad de gaps al final
    def cuadra2(self):
        """rellena con gaps las secuencias mas cortas"""
        import numpy
        seq = self.matrix.seqs
        maxLen = len(max(seq, key=len))
        # print (maxLen)
        for i in range(len(seq)):
            # if len(seq[i]) < maxLen: #original
            #     seq[i] = seq[i] + "-"*(maxLen-len(seq[i])) #original
            diferencia = maxLen-len(seq[i])
            if diferencia > 0:
                for j in range(diferencia):
                    posicion = random.randint(0, len(seq[i]))
                    parte1 = seq[i][:posicion]
                    parte2 = seq[i][posicion:]
                    temp = "-".join([parte1, parte2])
                    seq[i] = temp

        self.matrix.seqs = numpy.array(seq)


    """metodo para saber si alguna columna de self.matrix tiene  gap en todos los elementos"""
    def gapColumn(self, col):
        for i in range(len(self.matrix.seqs)):
            if self.matrix.seqs[i][col] != "-":
                return False
        return True



    """metodo que recorre la matriz y elimina las columnas con gaps en todos los elementos"""
    def limpiaColumnas(self):
        i = 0
        while i < len(self.matrix.seqs[0]):
            if self.gapColumn(i):
                self.deleteCulmn(i)
            else:
                i += 1


        """metodo para eliminar un elemento especifico en cada secuencia"""
    def deleteCulmn(self, pos):
        for i in range(len(self.matrix.seqs)):
            self.matrix.seqs[i] = self.matrix.seqs[i][:pos] + self.matrix.seqs[i][pos+1:]





        """metodo para obtener una lista con los elementos de cada columna"""
    def getColumn(self, col):
        column = []
        for i in range(len(self.matrix.seqs)):
            column.append(self.matrix.seqs[i][col])
        return column



        """metodo para evaluar columnas"""
    def autoEvalua(self):
        evaluador = evaluadorBlosum()
        score = 0
        for i in range(len(self.matrix.seqs[0])):
            column = self.getColumn(i)
            """cuenta gaps de columna"""
            gapCount = column.count("-")
            """eliminar gaps de columna"""
            column = [x for x in column if x != "-"]
            """metodo para recorrer todos los pares unicos y enviarlos a evaluador"""
            pares = self.obtener_pares_unicos(column)
            for par in pares:
                score += evaluador.getScore(par[0], par[1])
            """si hay gaps en la columna, se penaliza"""
            score -= gapCount*2
        self.blosumScore = score
        self.NFE += 1


    def obtener_pares_unicos(self, columna):
        pares_unicos = set()
        for i in range(len(columna)):
            for j in range(i+1, len(columna)):
                par = tuple(sorted([columna[i], columna[j]]))
                pares_unicos.add(par)
        return list(pares_unicos)

