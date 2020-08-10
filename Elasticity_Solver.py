#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 15 15:15:03 2020

@author: merlenye
"""
from __future__ import division
import numpy as np
from technical_tools_PC_3_alt import Tech_total
from economic_tools_elasticity_copy import Econ_total
import xlsxwriter
import math
import pandas as pd
def Verify_Solver(coeff, loadkWh):
    LoadKW_MAK = pd.read_excel('LoadKW_MAK.xlsx',index_col=None, header=None)
    loadkWh_orig = sum(LoadKW_MAK[0])
    workbook = xlsxwriter.Workbook('Verify_Solver.xlsx')
    worksheet = workbook.add_worksheet()
    worksheet.write(0,1, "Coefficient")
    worksheet.write(1,1, coeff)
    worksheet.write(0,2, "Original_8760_sum")
    worksheet.write(1,2, loadkWh_orig)
    worksheet.write(0,3, "New_8760_sum")
    worksheet.write(1,3, (loadkWh*12))
    worksheet.write(0,4, "Original_8760_by_coeff")
    worksheet.write(1,4, (loadkWh_orig*coeff))
    worksheet.write(0,5, "degree_difference")
    worksheet.write(1,5, ((loadkWh_orig*coeff)/(loadkWh*12)))
    workbook.close()
    
    
def Elasticity(implicit_tariff, price_elasticity, BattGuess, PVGuess):
    print(implicit_tariff)
    print(price_elasticity)
    guesses= np.array([BattGuess, PVGuess])
    Propane_ec, Batt_SOC, LoadkW, P_gen, P_PV, P_batt, P_dump,Limit_charge, Limit_discharge, BattkW, Batt_kWh_tot_ec,loadkWh,peakload, lifecycle = Tech_total(guesses[0],guesses[1],0)
    baseline_run = np.array(Econ_total(Propane_ec,guesses[1]*peakload,guesses[0]*peakload,Batt_kWh_tot_ec,peakload,loadkWh, lifecycle))
    tariff = baseline_run[8]
    target_tariff = implicit_tariff
    print(target_tariff)
    coeff=1
    while (target_tariff != tariff):
        print("tariff" + str(tariff))
        target_tariff = (target_tariff+tariff)/2
        coeff = 1.69388687849908 + 0.233593301304739*(np.log(implicit_tariff)) + 1.96372615122544*implicit_tariff*price_elasticity + 0.330369124083302*target_tariff**2 - price_elasticity - price_elasticity*target_tariff - np.log(implicit_tariff)*math.atanh(math.atanh(0.384738252918903*price_elasticity*target_tariff)) - 0.932075208105459*target_tariff
        #Equation derived using eurequa and extensive tariff and price elasticity data set. R^2 = 0.972
        print("coeff"+str(coeff))
        if coeff <= 1:
            Propane_ec, Batt_SOC, LoadkW, P_gen, P_PV, P_batt, P_dump,Limit_charge, Limit_discharge, BattkW, Batt_kWh_tot_ec,loadkWh,price_elasticityakload, lifecycle = Tech_total(guesses[0],guesses[1],coeff)
            tech_values = np.array([Propane_ec, Batt_SOC, LoadkW, P_gen, P_PV, P_batt, P_dump,Limit_charge, Limit_discharge, BattkW, Batt_kWh_tot_ec,loadkWh,peakload, lifecycle])
            generation_solution = np.array(Econ_total(Propane_ec,guesses[1]*peakload,guesses[0]*peakload,Batt_kWh_tot_ec,peakload,loadkWh, lifecycle))
            if(tariff-target_tariff) <0.05: #.05 of a cent accuracy. Increase to imrpove accuracy or decrease to improve runtime
                target_tariff =tariff
    Verify_Solver(coeff, loadkWh)
    return np.array([tech_values[0], tech_values[1], tech_values[2], tech_values[3], tech_values[4], tech_values[5], tech_values[6],tech_values[7], tech_values[8], tech_values[9], tech_values[10],tech_values[11],tech_values[12], tech_values[13], generation_solution[0], generation_solution[1], generation_solution[2], generation_solution[3], generation_solution[4], generation_solution[5], generation_solution[6], generation_solution[7], generation_solution[8], generation_solution[9]])
    
    
