#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  8 07:23:08 2022

@author: kirill
"""

for i in range(90):
    summ=0
    counter=0
    sum_2=0
    file = open("energy."+str(i+1)+".log", "r")
    for line in file:
        try:
            c = list(map(float, line.split()))
            summ+=c[1]
            counter+=1
        except:
            pass
    av=summ/counter
    file.close()
    file = open("energy."+str(i+1)+".log", "r")
    for line in file:
        try:
            c = list(map(float, line.split()))
            sum_2+=c[1]**2
        except:
            pass
    av_2=sum_2/counter
    print(*[1+0.1*i, av, av**2, av_2])
    file.close()