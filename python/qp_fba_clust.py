# -*- coding: utf-8 -*-
"""
Created on Tue Jul  9 10:28:04 2013

Does the optimization with gurobi

@author: vitoz
"""

#!/usr/bin/python

import numpy as np
import re
from gurobipy import *
import scipy
from scipy import linalg, matrix

def null(A, eps=1e-15):
    u, s, vh = scipy.linalg.svd(A)
    null_mask = (s <= eps)
    null_space = scipy.compress(null_mask, vh, axis=0)
    return scipy.transpose(null_space)

def dense_optimize(rows, cols, c, Q, A, sense, rhs, lb, ub, vtype,
                   rownames,colnames):
    model = Model()
    model.setParam('BarConvTol',1e-16)
    model.setParam('OptimalityTol',1e-9)    
    model.setParam('FeasibilityTol',1e-9) 
    model.setParam('Method',1) 
    # Add variables to model
    for j in range(cols):
        model.addVar(lb=lb[j], ub=ub[j], vtype=vtype[j],name=colnames[j])
        #model.addVar(lb=lb[j], ub=ub[j], vtype=vtype[j])
    model.update()
    vars = model.getVars()
    
    # Populate A matrix
    for i in range(rows):
        expr = LinExpr()
        for j in range(cols):
            if A[i][j] != 0:
                expr += A[i][j]*vars[j]
        model.addConstr(expr, sense[i], rhs[i],name=rownames[i])
        #model.addConstr(expr, sense[i], rhs[i])
        
    # Populate objective
    obj = QuadExpr()
    for i in range(cols):
        for j in range(cols):
            if Q[i][j] != 0:
                #obj += Q[i][j]*vars[i]*vars[j]
                #the /2 is arbitrary, but the R version does it and we want consistency                
                obj += Q[i][j]*vars[i]*vars[j]
    for j in range(cols):
        if c[j] != 0:
            obj += c[j]*vars[j]
    model.setObjective(obj)
    model.setParam('LogFile',"")
    # Write model to a file
    model.update()
    #model.write('dense.lp')
    
    # Solve
    model.optimize()
    return model


def import_Rtab_numpy(filepath):
    # imports a table exported from R with rownames and colnames
    with open(filepath,'r') as tsv:
        tabinput = [line.strip().replace('\"','').split('\t') for line in tsv]
    
    outTab = dict()
    outTab['colNames'] = tabinput[0]
    tabinput = tabinput[1:]
   
   # Check wether the content is numeric
    try:
        float(tabinput[1][1])
        isNum =1
    except ValueError:
        isNum = 0
        
    if isNum == 1:
        outTab['data'] = np.zeros((len(tabinput), len(tabinput[1])-1))
    else:
        outTab['data'] = np.chararray((len(tabinput), len(tabinput[1])-1))
    
    outTab['rowNames'] = list()
    for i in range(len(tabinput)):
        outTab['rowNames'].append(tabinput[i][0])
        for j in range(len(tabinput[0])-1):
            outTab['data'][i,j] = tabinput[i][j+1]
    return outTab
    

def import_Rtab(filepath):
    # imports a table exported from R with rownames and colnames
    with open(filepath,'r') as tsv:
        tabinput = [line.strip().replace('\"','').split('\t') for line in tsv]
        
    outTab = dict()
    outTab['v'] = tabinput[0]
    tabinput = tabinput[1:]
    outTab['data'] = list()
    outTab['rowNames'] = list()
     # Check wether the content is numeric
    try:
        float(tabinput[1][1])
        isNum =1
    except ValueError:
        isNum = 0
    for i in range(len(tabinput)):
        outTab['rowNames'].append(tabinput[i][0])
        if isNum:
            outTab['data'].append(map(float,tabinput[i][1:]))
        else:
            outTab['data'].append(tabinput[i][1:])
    # if the outTab data is just a vector, flatten it
    if len(tabinput[1][1:]) == 1:
        outTab['data'] = [y for x in outTab['data'] for y in x]
    return outTab

def flatten(inputList):
    return([y for x in inputList for y in x])

def readRpythonDat(pythonDat):
    inDat = dict()
    for line in pythonDat:
        line = line.replace('\n','')
        line = line.replace('"','')
        line = line.split('\t')
        if not(line[0] in inDat.keys()):   
            inDat[line[0]] = list()    
        try:
            inDat[line[0]].append(map(float,line[1:]))
        except ValueError:
            inDat[line[0]].append(line[1:])
            
    return inDat

def fva_qp(model,objvals):
    # does a variability analysis, where the objective function is constraint to
    # the optimal value and all the variables are max resp minimized given this
    # constraint.
    # works only fast with the model without the integer constraints
    obj = model.getObjective()
    
    if isinstance(objvals,float):
        objvals = [objvals]
    
    # set the obj function as constraint
    model.addQConstr(obj,GRB.LESS_EQUAL,objvals[0],name='formerObj')
    model.update()
    Qconstr = model.getQConstrs()
    idx = [idx for idx in range(Qconstr.__len__()) if Qconstr[idx].getAttr('QCname') == 'formerObj']
    objConstr = Qconstr[idx[0]]
    
    model.setParam('OutputFlag',0)
    model.setParam('Presolve',0) # for continuous models only
    model.setParam('BarHomogeneous',1) # for continuous models only
    model.setParam('TimeLimit',120)
    model.setParam('preqlinearize',1) # linearize the quadratic constraint
    model.setParam('BarQCPConvTol',1e-10)
    
    # add the variables in varList as objective
    # use the objective max/min: r_000_F - r_000_R
    modvars = model.getVars()
    allnames = model.getAttr('VarName',modvars)
    
    # only reactions
    fvaList = [] #name, min, max
    for val in objvals:
        names = [rxn for rxn in allnames if rxn[-2:] in ['_F','_R']]
        objConstr.setAttr('qcrhs',float(val))
        
        while names.__len__() >0:
            curVar = names[0]
            # check if it is a forward reaction
            if curVar[-2:] == '_F':
                if (curVar[:-2]+'_R') in names:
                    tobj = model.getVarByName(curVar) - model.getVarByName(curVar[:-2]+'_R')
                    names.remove(curVar[:-2]+'_R')
                else: tobj = model.getVarByName(curVar)
                rxnName =curVar[:-2]
            elif curVar[-2:] == '_R':
                tobj = -model.getVarByName(curVar)
                rxnName =curVar[:-2]
            else:
                tobj = model.getVarByName(curVar)
                rxnName =curVar
                
            names.remove(curVar)
            
            # determine the minimum flux through reaction which is consistent with the QP solution
            
            model.setObjective(tobj,GRB.MINIMIZE )
            model.optimize()
            minCode = model.status
            
            # if an optimal solution (2) is found continue
            # if a suboptimal solution (13) is found, save and rerun
            # otherwise rerun up to 10 times
            
            if minCode == 2:
                minV= model.getObjective().getValue()
            else:
                count = 0
                bestSol = float("-inf")
                while (count < 9 and minCode != 2):
                    if minCode == 13:
                        bestSol = model.getObjective().getValue()
                    model.reset()
                    model.optimize()
                    minCode = model.status
                    count += 1
                minV = bestSol
            
            fvaList.append([rxnName+'_FVA_'+str(val) +'_min',minV, 'status_'+str(minCode), 'viol_'+str(model.BoundSVioSum + model.ConstrVioSum)])
            
            # determine the maximum flux through reaction which is consistent with the QP solution
            
            model.setObjective(tobj,GRB.MAXIMIZE)
            model.optimize()
            maxCode = model.status
            
            if maxCode == 2:
                maxV= model.getObjective().getValue()
            else:
                count = 0
                bestSol = float("inf")
                while (count < 9 and maxCode != 2):
                    if maxCode == 13:
                        bestSol = model.getObjective().getValue()
                    model.reset()
                    model.optimize()
                    maxCode = model.status
                    count += 1
                maxV = bestSol
            
            fvaList.append([rxnName+'_FVA_'+str(val) +'_max',maxV, 'status_'+str(maxCode), 'viol_'+str(model.BoundSVioSum + model.ConstrVioSum)])
            
            if len(fvaList) > 5:
                print fvaList
                quit()
    return fvaList
    
    

## Start script

#pythonDat = sys.stdin.readlines()
pythonDat = open('./test_files/pythonDat_13.txt','r').readlines()
inDat = readRpythonDat(pythonDat)

## calc some input
nrow = inDat['A'].__len__()
ncol = inDat['A'][1].__len__()
vtype = [GRB.CONTINUOUS] * ncol
for i in range(len(inDat['ub'][0])):
    inDat['ub'][0][i] = min(inDat['ub'][0][i],GRB.INFINITY)

# do the optimization:
if inDat['mode'][0][0] == 'simple': 
    model = dense_optimize(nrow,ncol,inDat['obj'][0], inDat['Q'], 
                           inDat['A'], inDat['sense'][0], inDat['rhs'][0], 
    inDat['lb'][0],inDat['ub'][0], vtype, inDat['rownamesA'][0], inDat['colnamesA'][0])
    model.setParam('Presolve',0) # for continuous models only -> for FVA
    model.setParam('BarHomogeneous',1) # for continuous models only -> for FVA

# return the optimized values
vars = model.getVars()
values = model.getAttr('x', vars)
names = model.getAttr('VarName',vars)
optstatus = [''] * len(names)
optviol = [''] * len(names)

if inDat['FVA'][0][0] == 'T':
    obj = model.getObjective()
    objvals = obj.getValue()
    fvaRes = fva_qp(model,objvals)
    for el in fvaRes:
        names.append(el[0])
        values.append(el[1])
        optstatus.append(el[2])
        optviol.append(el[3])

sys.stdout.write('output_start'+'\n')
for i in range(len(values)):
    sys.stdout.write(names[i] + '\t' + str(values[i]) + '\t' + str(optstatus[i]) + '\t' + str(optviol[i]) + '\n')
sys.stdout.write('output_end'+'\n')
