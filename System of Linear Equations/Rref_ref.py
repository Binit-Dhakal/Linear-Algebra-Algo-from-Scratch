#This is the algorithm of row reduced echelon form(Gauss Jordan) and row echelon form(Gauss elimination)

import numpy as np

class GaussElimination:
    def __init__(self,matrix):
        self.matrix = matrix
        self.row, self.col = matrix.shape[0], matrix.shape[1]
        
    def solve(self):
        R, C = 0, 0  # pivot row number, pivot column number
        for i in range(self.row - 1):
            #handling the case where pivot has zero element
            k = R     #for iterating rows to find the non zero pivot
            while ( self.matrix[k][C] == 0 ):
                # for checking if the search has reached at the final row and final column entry 
                if (k + 1 == self.row and C == self.col - 1 ):
                    return self.matrix
                
                # for checking if the final row number has been reached: if it is reached 
                #we start our search to next column of initial row
                if(k + 1 == self.row ):
                    C += 1
                    k = R
                    continue
                
                k += 1     
                if (self.matrix[k][C] != 0):
                    self.matrix[[k-1,k],:] = self.matrix[[k,k-1],:]  #swap rows
                    break
                    
            for j in range(R, self.row - 1):     
                if (matrix[j + 1][C] != 0):
                    mul = self.matrix[j + 1][C] / self.matrix[R][C]
                    self.matrix[j + 1] = self.matrix[j + 1] -  self.matrix[R] * mul
            R += 1
            C += 1
            
        return self.matrix




class GaussJordan:
    def __init__(self,matrix):
        self.matrix = matrix
        self.row, self.col = matrix.shape[0], matrix.shape[1]
        
    def solve(self):
        R, C = 0, 0  # pivot row number, pivot column number
        for i in range(self.row-1):
            #handling the case where pivot has zero element
            k = R     #for iterating rows to find the non zero pivot
            while ( self.matrix[k][C] == 0 ):
               
                # for checking if the search has reached at the final row and final column entry       
                if (k == self.row - 1 and C == self.col - 1 ):
                    return self.matrix
                
                # for checking if the final row number has been reached: if it is reached we start our search to next column of initial row
                if(k == self.row - 1):
                    C += 1
                    k = R
                    continue
                
                k += 1     
                if (self.matrix[k][C] != 0):
                    self.matrix[[k-1,k],:] = self.matrix[[k,k-1],:]  #swap rows
                    break
            
            # making the pivot entry 1 for gauss jordan
            if( self.matrix[R][C] != 1):
                self.matrix[R] = self.matrix[R] / self.matrix[R][C]
                
            for j in range(0, self.row):     
                if (self.matrix[j][C] != 0 and j != R):
                    mul = self.matrix[j][C] / self.matrix[R][C]
                    self.matrix[j] = self.matrix[j] -  self.matrix[R] * mul

            R += 1
            C += 1
            
        return self.matrix