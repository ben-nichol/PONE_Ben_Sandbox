#!/usr/bin/env python

import numpy as np

class FunctionFromTable:

    def __init__(self, inputs, outputs):
        self.inputs = inputs
        self.outputs = outputs

    def getValue(self, inputVal):
        # spacing might not be even just search through

        if inputVal > self.inputs[len(self.inputs)-1] or inputVal < self.inputs[0]:
            raise RuntimeError("Error: input value not in domain")
        
        index = 0
        for i in range(len(self.outputs)): 
            if( inputVal <= self.inputs[i+1]):
                index = i
                break
        fraction = (inputVal - self.inputs[index])/(self.inputs[index+1] - self.inputs[index])

        return self.outputs[index] + (self.outputs[index+1]-self.outputs[index])*fraction

class Polynomial:
    
    def __init__(self, coefficients, xmin, xmax):
        if len(coefficients) == 0:
            raise RuntimeError("no proper coefficient model passed")
    
        self.coefficients = coefficients
        self.xmin = xmin
        self.xmax = xmax

        self.normalize()

    def normalize(self):
        x = np.linspace(self.xmin, self.xmax, 1000)
        dx = x[1] - x[0]
        integral = 0.0

        for i in range(len(x)):
            y = 0
            for j in range(len(self.coefficients)):
                y += x[i] ** j * self.coefficients[j]
        
            integral += y*dx
    
        self.coefficients = [coefficient/integral for coefficient in self.coefficients]

    def getValue(self, inputVal):
        if inputVal > self.xmax or inputVal < self.xmin:
            raise RuntimeError("Error: input value not in domain")
        
        sumVal = 0
        for i in range(len(self.coefficients)):
            sumVal += inputVal**i * self.coefficients[i]
    
        return sumVal