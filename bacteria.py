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
        """Rellena las secuencias más cortas con gaps, asegurando que todas tengan la misma longitud."""
        seq = self.matrix.seqs
        maxLen = len(max(seq, key=len))  # Longitud máxima de las secuencias

        for i in range(len(seq)):
            # Calcula la diferencia entre la longitud máxima y la longitud actual
            diferencia = maxLen - len(seq[i])

            # Inserta gaps solo si la secuencia es más corta que la longitud máxima
            if diferencia > 0:
                # Crea una lista de posiciones para insertar gaps
                for _ in range(diferencia):
                    # Determina una posición aleatoria para insertar el gap
                    posicion = random.randint(0, len(seq[i]))  # Puede ser igual a la longitud

                    # Inserta el gap en la posición seleccionada
                    seq[i] = seq[i][:posicion] + "-" + seq[i][posicion:]

        self.matrix.seqs = numpy.array(seq)  # Asegúrate de que todas las secuencias tengan la misma longitud




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
        conserved_columns = []

        # Evaluar cada columna
        for i in range(len(self.matrix.seqs[0])):
            column = self.getColumn(i)
            gapCount = column.count("-")
            column_no_gaps = [x for x in column if x != "-"]

            # Contar la frecuencia de cada aminoácido/nucleótido
            frequency = {}
            for char in column_no_gaps:
                if char in frequency:
                    frequency[char] += 1
                else:
                    frequency[char] = 1

            # Evaluar la conservación
            most_common = max(frequency.values(), default=0)
            conservation_score = (most_common / len(column_no_gaps)) if column_no_gaps else 0

            # Penalizar si hay gaps
            score -= gapCount * 4
            score += conservation_score * 3  # Aumentar el score por conservación

            # Evaluar pares únicos
            pares = self.obtener_pares_unicos(column_no_gaps)
            for par in pares:
                score += evaluador.getScore(par[0], par[1])

        self.blosumScore = score
        self.NFE += 1


    def obtener_pares_unicos(self, columna):
        pares_unicos = set()
        for i in range(len(columna)):
            for j in range(i+1, len(columna)):
                par = tuple(sorted([columna[i], columna[j]]))
                pares_unicos.add(par)
        return list(pares_unicos)

