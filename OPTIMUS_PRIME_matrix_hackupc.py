# -*- coding: utf-8 -*-
"""
Created on Sat Apr 30 00:03:32 2022

@author: VSU
"""

import pandas as pd
import numpy as np
import datetime


#read CSV
test = pd.read_csv("test.csv")
geo_params = pd.read_csv("geo_params.csv")
test = pd.read_csv("test.csv")
sales = pd.read_csv("sales.csv")
sku = pd.read_csv("sku.csv")


sales.count()


n_sku = len(sku)
n_geo = len(geo_params)
vdates = sales['date'].unique()
vdates.sort() 
n_dates = len(vdates)
n_sales = len(sales)

# functions to obtain partial data
def get_day(sku_id, date_str):
    vector_input = np.zeros(n_geo)
    slice_day = sales[sales['date'] == date_str]
    slice_sales = slice_day[slice_day['SKU'] == sku_id]
    for i_geo in range(n_geo):
        geo_id = geo_params['geoCluster'][i_geo]
        point = slice_sales[slice_sales['geoCluster'] == geo_id]['sales']
        if len(point) == 1: 
            if not np.isnan(list(point)[0]):
                vector_input[i_geo] = list(point)[0]
        
    return vector_input

def get_day_all(date_str):
    vector_day = np.zeros((n_sku,n_geo))
    slice_day = sales[sales['date'] == date_str]
    for i_sku in range(n_sku):
        sku_id = sku['SKU'][i_sku]
        vector_input = np.zeros(n_geo)
        slice_sales = slice_day[slice_day['SKU'] == sku_id]
        for i_geo in range(n_geo):
            geo_id = geo_params['geoCluster'][i_geo]
            point = slice_sales[slice_sales['geoCluster'] == geo_id]['sales']
            if len(point) == 1: 
                if not np.isnan(list(point)[0]):
                    vector_input[i_geo] = list(point)[0]
        vector_day[i_sku,:] = vector_input
        #vector_day[i_sku,:] = get_day(sku_id, date_str)
    
    return vector_day

# We will use index for SKU and geoCluster
dict_sku = {}
for i_sku in range(n_sku):
    sku_id = sku['SKU'][i_sku]
    dict_sku[sku_id] = i_sku
    
dict_geo = {}
for i_geo in range(n_geo):
    geo_id = geo_params['geoCluster'][i_geo]
    dict_geo[geo_id] = i_geo


def get_all():   # Data into a matrix
    vector_all = np.zeros((n_dates, n_sku, n_geo, 2))
    vector_all[:,:,:,:] = np.nan    # Useful later
    for index, row in sales.iterrows():
        if index%10000==0:
            print(index)
        geo_id = row['geoCluster']
        sku_id = row['SKU']
        date_str = row['date']
        price_value = row['price']
        sales_value = row['sales']
        
        i_geo = dict_geo[geo_id]
        i_sku = dict_sku[sku_id]
        l_date = list(map(int,date_str.split('-')))
        
        i_date = (datetime.datetime(l_date[0],l_date[1],l_date[2])-datetime.datetime(2020,1,1)).days
        #NANs in sales become zeroes
 
        if np.isnan(sales_value) or sales_value == 0:
            price_value = 700 #lol
            sales_value = 0
               
        vector_all[i_date, i_sku, i_geo, 0] = price_value
        vector_all[i_date, i_sku, i_geo, 1] = sales_value
    for i_date in range(n_dates):
        for i_sku in range(n_sku):
            for i_geo in range(n_geo):
                if vector_all[i_date,i_sku,i_geo,0] == 700:
                    vector_all[i_date,i_sku,i_geo,0] = 0    # This solves some weird issue
    return vector_all

def get_ampli(vector_all):     # Extends data with test values
    vector_ampli = np.zeros((n_dates+14, n_sku, n_geo, 2))
    vector_ampli[0:n_dates,:,:,:] = vector_all[:,:,:,:]
    for index, row in test.iterrows():
        if index%10000==0:
            print(index)
        geo_id = row['geoCluster']
        sku_id = row['SKU']
        date_str = row['date']
        price_value = row['price_filled']
        #sales_value = row['sales']
        
        i_geo = dict_geo[geo_id]
        i_sku = dict_sku[sku_id]
        l_date = list(map(int,date_str.split('-')))
        
        i_date = (datetime.datetime(l_date[0],l_date[1],l_date[2])-datetime.datetime(2020,1,1)).days
        
        vector_ampli[i_date, i_sku, i_geo, 0] = price_value
        
    return vector_ampli

def populate(vector_all):
    # Now we will populate prices with neighbour values
    n = len(vector_all)
    for i_sku in range(n_sku):
        print('i_sku: ', i_sku)
        for i_geo in range(n_geo):
            f_store_product = vector_all[:, i_sku, i_geo, 0]
            for i_date in range(1,n):
                if f_store_product[i_date] == 0 and f_store_product[i_date-1]>0:
                    f_store_product[i_date] = f_store_product[i_date - 1]
            for i_date in range(1,n):
                if f_store_product[n - i_date-1] == 0 and f_store_product[n - i_date]>0:
                    f_store_product[n - i_date-1] = f_store_product[n - i_date]
                    
            vector_all[:, i_sku, i_geo, 0] = f_store_product 
                        
            
    return vector_all

# We obtain matrixes and save them into .npy files, that will be used
# by our models
vector_all = get_all()
vector_all = populate(vector_all)
np.save('sales.npy', vector_all)

vector_ampli = get_ampli(vector_all)
vector_ampli = populate(vector_ampli)
np.save('sales_expanded.npy', vector_all)


