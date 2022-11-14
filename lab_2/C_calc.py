#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 11 15:09:15 2022

@author: kirill
"""
diff = 0
file_input = open("dump.txt", "r")
for line in file_input:
    if (diff == 0):
        c_prev = list(map(float, line.split()))
        diff = 1
    else:
        c_cur = list(map(float, line.split()))
        diff = (c_cur[1] - c_prev[1])/(c_cur[0] - c_prev[0])
        diff_1 = (c_cur[3] - c_cur[2])/(c_cur[0]**2)
        print(*[c_cur[0], diff, diff_1])
        c_prev=[c_cur[i] for i in range(len(c_cur))]
        