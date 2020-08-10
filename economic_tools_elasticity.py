#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 10 08:16:48 2020

@author: merlenye
"""

from __future__ import division
import numpy as np
import math
import pandas as pd

## MCASHFLOW FUNCTION =====================================================================================================
def mcashflow(tariff_hillclimb_multiplier,lifetime,f_pv,a_pv,f,a,Batt_life_yrs, equity_pct, term, loadkwh_lec, interest_rate, construction_period, grace_period, loanfactor, PVkW, BattKWh, LEC, C1_pv, C1_LPG, Cost_bank, Cost_Propane_yr, min_irr, dsc_ratio, price_elasticity, T_sat, construction_rate):
#Removed all thermal system variables and calculations

    #construction_period: period of time starting from month t=1 to t=construction_period when there is no O&M cost
    #grace_period: period of time starting from month t=1 to t=grace_period when there is no finance cost
    monthly_interest_rate=(1+interest_rate)**(1/12)-1#interest_rate is annual
    Batt_life_mths = Batt_life_yrs*12
    grace_period=int(grace_period)
    Cost_Propane_mth = Cost_Propane_yr/12
    construction_period = int(construction_period)
    loanfactor = 1

    #Initialize output variables
    lifetime = int(lifetime)
    LoanPrincipal = np.zeros(lifetime)
    Interest = np.zeros(lifetime)
    year = np.zeros(lifetime)
    month = np.zeros(lifetime)
    Cost = np.zeros(lifetime)
    Revenue = np.zeros(lifetime)
    BalanceCost=np.zeros(lifetime)
    BattCost=np.zeros(lifetime)
    CashonHand = np.zeros(lifetime)
    Balance = np.zeros(lifetime)
    M = np.zeros(lifetime)
    O = np.zeros(lifetime)
    debt_service=np.zeros(lifetime)
    income_debt_ratio=np.zeros(lifetime)#net operating income/debt service for a given month
    tariff = np.copy(LEC) # "USD/kWh tariff for electricity"
    Batt_penalty = 0

    debt_pct=1-equity_pct
    #ILD: "Interannual Load Dynamic"
    ild=np.zeros(lifetime)
    for t in range(1,lifetime):
        #ILD[j] represents demand at t=j months as a fraction of loadkwh_lec
        #the ILD function f is a piecewise. It is a sigmoid while t<=24, such that f(24)=1. On t>=24, f grows at 4% per year.
        #However, before the constructin period, demand is 0
        a=.2#adjusts steepness of curve
        b=(T_sat/2)#inflection point at t=3 months
        #Inflection point should be T_SAT time to saturation of demand (month #) b should be half that input
        c=.04#growth rate after sigmoid

        m=float(1+c)**float(1/12)
        if t<construction_period:
            ild[t]=0
        elif t<=2*b+construction_period:
            ild[t]=1/(1+math.exp(-a*(t-construction_period-b)))
        else:
            ild[t]=m*ild[t-1]
    adj_factor = True
    while adj_factor == True: #While loop that accounts for loan factor
        loanfactor+=.01
        capex = ((C1_pv+C1_LPG)*loanfactor)
        LoanPrincipal[0] = ((C1_pv+C1_LPG)*loanfactor)*(debt_pct)   #"The amount of project finance needed to cover capex and initial opex"
        LoanPrincipal[1] = np.copy(LoanPrincipal[0])
        equity=equity_pct*((C1_pv+C1_LPG)*loanfactor)
        CashonHand[0] = equity + LoanPrincipal[0]-(C1_pv+C1_LPG)
    
    
        if Batt_life_mths == 0:  #"If the battery needs to be replaced more often than once a year, inflict a heavy penalty"
            Batt_life_mths = 12
            Batt_penalty = 1
    
        for i in range(1,lifetime):
            month[i]=i
            j=int(i/12)+1
            year[i]=j
    
            x=i-construction_period#months after construction
    
            if(x>0):
                y=int(x/12)+1#k is years after the constructin period
            #"Operating cost is backup fuel - some ITC charges could also be added depending on the metering situation"
    
    
            #IF NOT CONSTRUCTION PERIOD, then there are O&M costs:
            #"Maintenance cost is a function of the type of equipment"
            if i<=construction_period:
                M[i]=0
                O[i]=0
            else:
                M[i] = ((f_pv*a_pv*C1_pv/lifetime)+(f_pv*(1-a_pv)*C1_pv/lifetime**2)*(2*y-1)+(f*a*C1_LPG/lifetime)+(f*(1-a)*C1_LPG/lifetime**2)*(2*y-1))/12
                O[i] =  (np.copy(Cost_Propane_mth))
    
            #cost of battery replacement every n=batt_life_yrs"
            Modulo = i % Batt_life_mths #find the remainder
            if Modulo == 0:
                BattCost[i] = Cost_bank + 2 * Cost_bank * Batt_penalty
    
    
            #IF IN GRACE PERIOD:
            #Cost does not include finance
            #Do not debit finance from the principal
            if i > 1:
                LoanPrincipal[i] = LoanPrincipal[i-1] + Interest[i-1] - debt_service[i-1]
    
            #INTEREST
            if i > construction_period + grace_period:
                Interest[i]=LoanPrincipal[i]*monthly_interest_rate#Interest to be applied to the next loanprincipal calculation
            elif i<= construction_period+ grace_period:
                Interest[i]=LoanPrincipal[i]*monthly_interest_rate*construction_rate
    
            amortized=LoanPrincipal[int(grace_period)]*(monthly_interest_rate+monthly_interest_rate/((1+monthly_interest_rate)**term-1))#Amortization formula
            if LoanPrincipal[i]>0:
                if i<grace_period:
                    debt_service[i]=0
                elif amortized>=LoanPrincipal[i]+Interest[i]:
                    debt_service[i]=LoanPrincipal[i]+Interest[i]
                else:#If not in grace period
                    debt_service[i]=amortized
            else:
                debt_service[i]=0
    
    
            #"Financial outlays in each year of operation include loan repayment and O&M"
            Cost[i] = Interest[i] + O[i] + M[i] + BattCost[i]
            BalanceCost[i] = O[i] + M[i]
    	    #"if the loan is paid off THEN there no finance charge"
    
        acceptable_irr = False
        if (equity_pct>0):
            #Calculate IRR
            present_value_terms = [0 for i in range(lifetime+1)]
            construction_period = int(construction_period)
            baseline = True
            while (not acceptable_irr  or not (income_debt_ratio[construction_period+grace_period:lifetime:1]>=dsc_ratio).all() or LoanPrincipal[lifetime-1]>0): #continue loop until all values in CashonHand[1:] are greater than 0.
                #removed the check for cashonhand being above 0 in final step and instead put check in for loop
                if adj_factor == True and baseline == False:
                    break
                if adj_factor == False:
                    tariff = tariff*tariff_hillclimb_multiplier #" Increase the tariff until the cash flows are positive "
                for j in range(1,lifetime):
                    Revenue[j]= loadkwh_lec * ild[j] * tariff
                    Balance[j] = Revenue[j] - BalanceCost[j]
                    CashonHand[j] = CashonHand[j-1] + Revenue[j] - Cost[j] - BattCost[j]
                    #if cash on hand goes negative
                    if CashonHand[j]<0 & j < grace_period + construction_period: #If cash on hand goes negative during the construction or grace period, break the loop
                        baseline = False
                        break
                    if j > grace_period + construction_period:
                        adj_factor= False
                    if CashonHand[j]<0 & adj_factor == False: #If cash on hand goes negative after 
                        break
                    if(debt_service[j]>0):
                        income_debt_ratio[j]=Balance[j]/debt_service[j]
                    else:
                        income_debt_ratio[j]=10
                    present_value_terms[i]=  Balance[i]/(min_irr**i)
                npv = sum(present_value_terms)
                acceptable_irr=(npv>capex)#npv-capex>0
    
        else:#if there is no equity, just make sure that the project doesn't lose money
            print("else")
            while (CashonHand[lifetime-1]<0) or not (income_debt_ratio>=dsc_ratio).all() or LoanPrincipal[lifetime-1]>0:
                  #removed the check for cashonhand being above 0 in final step and instead put check in for loop
                if adj_factor == True and baseline == False:
                    break
                if adj_factor == False:
                    tariff = tariff*tariff_hillclimb_multiplier #" Increase the tariff until the cash flows are positive "
                for j in range(1,lifetime):
                    Revenue[j]= loadkwh_lec * ild[j] * tariff
                    Balance[j] = Revenue[j] - BalanceCost[j]
                    CashonHand[j] = CashonHand[j-1] + Revenue[j] - Cost[j] - BattCost[j]
                    #if cash on hand goes negative
                    if CashonHand[j]<0 & j < grace_period + construction_period: #If cash on hand goes negative during the construction or grace period, break the loop
                        baseline = False
                        break
                    if j > grace_period + construction_period:
                        adj_factor= False
                    if CashonHand[j]<0 & adj_factor == False: #If cash on hand goes negative after 
                        break
                    if(debt_service[j]>0):
                        income_debt_ratio[j]=Balance[j]/debt_service[j]
                    else:
                        income_debt_ratio[j]=10
                    present_value_terms[i]=  Balance[i]/(min_irr**i)
                npv = sum(present_value_terms)
                acceptable_irr=(npv>capex)#npv-capex>0
    sum_payments=LoanPrincipal[0]+sum(Interest)
    sum_paid = sum(debt_service)
    print(loadkwh_lec)
    print("Sum Owed: " + str(sum_payments) + ", Sum Paid: " + str(sum_paid), " diff: " + str(sum_paid-sum_payments))
    print("tariff: " + str(tariff))
    writer = pd.ExcelWriter("economic_tools_output.xlsx", engine='xlsxwriter')
    d={'Month': month, 'Year': year, 'LoanPrincipal': LoanPrincipal, 'Interest': Interest, 'Cost': Cost, 'Revenue': Revenue, 'CashOnHand': CashonHand, 'Balance': Balance, 'Debt_Service': debt_service, 'Maintenance': M, 'Operation': O, 'ILD': ild, 'Income_Debt_Ratio': income_debt_ratio,'Tariff' : tariff, 'Loan_Factor' : loanfactor}
    cashflow_df=pd.DataFrame(d, columns =['Month', 'Year', 'LoanPrincipal', 'Interest', 'Cost', 'Revenue', 'CashOnHand', 'Balance', 'Debt_Service', 'Maintenance', 'Operation', 'ILD', 'Income_Debt_Ratio', 'Tariff', 'Loan_Factor'])
    cashflow_df.to_excel(writer, sheet_name='Cashflow')

    workbook  = writer.book
    worksheet = writer.sheets['Cashflow']
    worksheet.set_landscape()
    LoanChart = workbook.add_chart({'type': 'line'})
    LoanChart.add_series({
    #start_row, start_column, end_row, end_column
    'name':       ['Cashflow', 0, 3],#cost name cell
    'categories': ['Cashflow', 1, 1, lifetime-1, 1],#should be from 1 to lifetime instead of 3 to lifetime-1. changed for debugging purposes
    'values':     ['Cashflow', 1, 3, lifetime-1, 3],#lifetime-1 should be last row
    })
    LoanChart.set_x_axis({
    'name': 'Month'
    })
    LoanChart.set_y_axis({
    'name': '$'
    })
    LoanChart.set_legend({'position': 'none'})
    worksheet.insert_chart('P1', LoanChart)

    FinanceChart = workbook.add_chart({'type': 'line'})
    FinanceChart.add_series({
        'name':       ['Cashflow', 0, 5],#cost name cell
        'categories': ['Cashflow', 1, 1, lifetime-1, 1],
        'values':     ['Cashflow', 1, 5, lifetime-1, 5],#cost value cells
    })
    FinanceChart.add_series({
        'name':       ['Cashflow', 0, 6],#revenue name cell
        'categories': ['Cashflow', 1, 1, lifetime-1, 1],
        'values':     ['Cashflow', 1, 6, lifetime-1, 6],#revenue value cells
    })
    FinanceChart.set_x_axis({
    'name': 'Month'
    })
    FinanceChart.set_y_axis({
    'name': '$'
    })
    worksheet.insert_chart('P16', FinanceChart)

    CostBreakdownChart = workbook.add_chart({'type': 'line'})
    CostBreakdownChart.add_series({#Debt Service
        'name':       ['Cashflow', 0, 9],#cost name cell
        'categories': ['Cashflow', 1, 1, lifetime-1, 1],
        'values':     ['Cashflow', 1, 9, lifetime-1, 9],#cost value cells
    })
    CostBreakdownChart.add_series({#Maintenance
        'name':       ['Cashflow', 0, 10],#revenue name cell
        'categories': ['Cashflow', 1, 1, lifetime-1, 1],
        'values':     ['Cashflow', 1, 10, lifetime-1, 10],#revenue value cells
    })
    CostBreakdownChart.add_series({#Operations
        'name':       ['Cashflow', 0, 11],#revenue name cell
        'categories': ['Cashflow', 1, 1, lifetime-1, 1],
        'values':     ['Cashflow', 1, 11, lifetime-1, 11],#revenue value cells
    })
    CostBreakdownChart.set_x_axis({
    'name': 'Month'
    })
    CostBreakdownChart.set_y_axis({
    'name': '$'
    })

    worksheet.insert_chart('X1', CostBreakdownChart)


    COH_Chart = workbook.add_chart({'type': 'line'})
    COH_Chart.add_series({#Debt Service
        'name':       ['Cashflow', 0, 7],#cost name cell
        'categories': ['Cashflow', 1, 1, lifetime-1, 1],
        'values':     ['Cashflow', 1, 7, lifetime-1, 7],#cost value cells
    })
    COH_Chart.set_x_axis({
    'name': 'Month'
    })
    COH_Chart.set_y_axis({
    'name': '$'
    })
    worksheet.insert_chart('X16', COH_Chart)
    writer.save()
    return LoanPrincipal, year, Cost, Revenue, CashonHand, Balance, M, O, tariff
###================================================================================================


    

## Econ_total function ===========================================================================
def Econ_total(propane, PVkW,BattKWh,Batt_kWh_tot,peakload,loadkwh_lec, batt_lifecycle):

    #Load all Econ input
    Econ_Parameters = pd.read_excel('uGrid_Input.xlsx', sheet_name = 'Econ')
    #"factors for distributing maintenance costs in time as a function of capex see Orosz IMechE"
    f_pv= Econ_Parameters['f_pv'][0]
    a_pv=Econ_Parameters['a_pv'][0]
    f=Econ_Parameters['f'][0]
    a=Econ_Parameters['a'][0]

    #"Set the financial return and period"
    interest_rate=Econ_Parameters['interest_rate'][0]
    construction_period=Econ_Parameters['construction_period'][0]
    grace_period=Econ_Parameters['grace_period'][0]
    dsc_ratio = Econ_Parameters['dsc_ratio'][0]
    loanfactor=Econ_Parameters['loanfactor'][0]
    T_sat=Econ_Parameters['T_sat'][0] #Time to saturation
    construction_rate=Econ_Parameters['construction_rate']

    min_irr=Econ_Parameters['min_irr'][0]
    price_elasticity=Econ_Parameters['price_elasticity'][0]
    equity_pct=Econ_Parameters['equity_pct'][0]
    lifetime = Econ_Parameters['lifetime'][0]
    term = Econ_Parameters['term'][0]
    tariff_hillclimb_multiplier = Econ_Parameters['tariff_hillclimb_multiplier'][0]

    #"Convert battery throughput into lifetime"
    Batt_life_yrs = np.floor((BattKWh*batt_lifecycle)/(Batt_kWh_tot+0.01))  #"Years of battery life before replacement is necessary, rounded down to an integer"

    #"Cost functions"
    Pole_num=Econ_Parameters['Dist_km'][0] /0.050   #"1 pole for every 50m distribution we"
    Cost_panels=PVkW*Econ_Parameters['Cost_panel_per_kW'][0]  #"PV Price via Alibaba 2016"
    Cost_charge_controllers=Econ_Parameters['Cost_charge_controllers_per_kW'][0]*PVkW
    Cost_Smartmeter=65*Econ_Parameters['node_num'][0]  #"Iometer"
    Cost_MPesa = Econ_Parameters['Cost_Mpesa_per_kWLoad'][0]*peakload  #"Estimate for merchant services with vodacom"
    Cost_inv=peakload*Econ_Parameters['Cost_inv_per_kWLoad'][0] #"[$/kW peak]"
    Cost_EPC_tracker=Econ_Parameters['Cost_EPC_tracker_per_kW'][0]*PVkW

    #"Cost aggregators"
    C1_LPG = (-10354.1143  + 6192.606 * math.log(peakload))   #"Propane Genset costs"  "Based on generac lineup"

    Cost_bank = BattKWh * Econ_Parameters['Cost_batt'][0]   #"[NREL, USAID Tetratech, health mgmt., PIH] "

    Cost_Propane_yr = propane*1.3  #"USD/kg"  "RSA prices 2016"

    Cost_Dist = Econ_Parameters['Cost_Dist_wire'][0] * Econ_Parameters['Dist_km'][0] + Econ_Parameters['Cost_Step_up_Trans'][0] * Econ_Parameters['Step_up_Trans_num'][0] + Econ_Parameters['Cost_Pole_Trans'][0] * Econ_Parameters['Pole_Trans_num'][0] + Econ_Parameters['Cost_Pole'][0] * Pole_num

    Cost_BOS = Cost_bank + Cost_inv + Econ_Parameters['Cost_control'][0] + Cost_Dist + Cost_Smartmeter + Cost_MPesa + Cost_charge_controllers   #"Balance of System"

    Cost_EPC = Cost_EPC_tracker + Econ_Parameters['Cost_EPC_LPG_tank'][0] + Econ_Parameters['Cost_EPC_Power_house'][0] + Econ_Parameters['Cost_EPC_Labor_Plant'][0] + Econ_Parameters['Cost_EPC_Labor_Dist'][0]

    Cost_Dev = Econ_Parameters['Cost_Dev_land'][0] + Econ_Parameters['Cost_Dev_EIA'][0] + Econ_Parameters['Cost_Dev_connection'][0] + Econ_Parameters['Cost_Dev_ICT'][0] + Econ_Parameters['Cost_Dev_contingency'][0] + Econ_Parameters['Cost_Dev_overhead'][0] + Econ_Parameters['Cost_taxes'][0]

    C1_pv = Cost_panels + Cost_BOS + Cost_EPC + Cost_Dev

    LEC = 0.1 #this is a starting point for LEC. This could potentially be done without a hill-climb and be dectly solved

    return mcashflow(tariff_hillclimb_multiplier,lifetime,f_pv,a_pv,f,a,Batt_life_yrs, equity_pct, term, loadkwh_lec, interest_rate, construction_period, grace_period, loanfactor, PVkW, BattKWh, LEC, C1_pv, C1_LPG, Cost_bank, Cost_Propane_yr, min_irr, dsc_ratio, price_elasticity, T_sat, construction_rate)
    

##===============================================================================================

#Run economic code standalone
if __name__ == "__main__":
    #These are for running as standalone
    propane=2254.555
    PVkW=196.166
    BattKWh=248.19
    Batt_kWh_tot=73215.76
    LoadKW_MAK = pd.read_excel('LoadKW_MAK.xlsx',index_col=None, header=None)
    peakload_buffer = 1.2
    peakload=44.87
    #max(LoadKW_MAK[0])*peakload_buffer
    loadkwh_lec = (121221.56/12)
    batt_lifecycle = 1200
    #sum(LoadKW_MAK[0])/12 #Estimated monthly load consumption

    LoanPrincipal, year, Cost, Revenue, CashonHand, Balance, M, O, tariff, Batt_life_yrs = Econ_total(propane,PVkW,BattKWh,Batt_kWh_tot,peakload,loadkwh_lec, batt_lifecycle)