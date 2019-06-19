"""
The purpose of this script is to download ICESat-2 data in a defined bounding box and temporal range. 
This script only needs to be edited in the user input sections. 

The script gathers IceSat-2 processed data with the user defined temporal range and spatial bounding box. 
It will place data in the AWS S3 cloud (aws s3 ls s3://pangeo-data-upload-oregon/icesat2/Snowblower/).
You will be prompted for your NASA earth data login

To run the scipt, open a terminal and execute the following command:
$ ipython get_IC2_Data.py

"""

# Import Packages
import requests
import getpass
import socket
import json
import zipfile
import io
import math
import os
import shutil
import pprint
import time
import geopandas as gpd
import matplotlib.pyplot as plt
import fiona
import h5py
import re
# To read KML files with geopandas, we will need to enable KML support in fiona (disabled by default)
fiona.drvsupport.supported_drivers['LIBKML'] = 'rw'
from shapely.geometry import Polygon, mapping
from shapely.geometry.polygon import orient
from statistics import mean
from requests.auth import HTTPBasicAuth

############################################################ Begin User Input

# Download note
download_note = "Novo_bbox"

# Input data set ID (e.g. ATL06) of interest here, also known as "short name".
short_name = 'ATL06'

# Input bounding box
# Input lower left longitude in decimal degrees
LL_lon = '9'
# Input lower left latitude in decimal degrees
LL_lat = '-72'
# Input upper right longitude in decimal degrees
UR_lon = '15'
# Input upper right latitude in decimal degrees
UR_lat = '-70'

# Input temporal range 
# Input start date in yyyy-MM-dd format
start_date = '2018-10-15'
# Input start time in HH:mm:ss format
start_time = '00:00:00'
# Input end date in yyyy-MM-dd format
end_date = '2018-10-31'
# Input end time in HH:mm:ss format
end_time = '23:59:59'


# Earthdata Login credentials
# Enter your Earthdata Login user name
uid = 'eric_keenan'
# Enter your email address associated with your Earthdata Login account
email = 'eric.keenan@colorado.edu'
pswd = getpass.getpass('Earthdata Login:')

############################################################ End User Input

# Update download note
download_note = download_note + "_" + LL_lon + "_" + LL_lat + "_" + UR_lon + "_" + UR_lat + "_" + start_date + \
    "_" + start_time + "_" + end_date + "_" + end_time
print("Download Note: " + download_note)

# Request token from Common Metadata Repository using Earthdata credentials
token_api_url = 'https://cmr.earthdata.nasa.gov/legacy-services/rest/tokens'
hostname = socket.gethostname()
ip = socket.gethostbyname(hostname)

data = {
    'token': {
        'username': uid,
        'password': pswd,
        'client_id': 'NSIDC_client_id',
        'user_ip_address': ip
    }
}
headers={'Accept': 'application/json'}
response = requests.post(token_api_url, json=data, headers=headers)
token = json.loads(response.content)['token']['id']
# print(token)

# Get json response from CMR collection metadata and print results. This provides high-level metadata on a data set or "collection", provide in json format.

params = {
    'short_name': short_name
}

cmr_collections_url = 'https://cmr.earthdata.nasa.gov/search/collections.json'
response = requests.get(cmr_collections_url, params=params)
results = json.loads(response.content)
# pprint.pprint(results)

# Find all instances of 'version_id' in metadata and print most recent version number
versions = [i['version_id'] for i in results['feed']['entry']]
latest_version = max(versions)
# print(latest_version)

temporal = start_date + 'T' + start_time + 'Z' + ',' + end_date + 'T' + end_time + 'Z'
print("Temporal Range: " + temporal)

# # Commenting for tutorial since we will be walking through option 3 (spatial file input) together
# # Bounding Box spatial parameter in 'W,S,E,N' format

bounding_box = LL_lon + ',' + LL_lat + ',' + UR_lon + ',' + UR_lat
# aoi value used for CMR params below
aoi = '1'
print("Bounding Box: " + bounding_box)

#Create CMR parameters used for granule search. Modify params depending on bounding_box or polygon input.

if aoi == '1':
# bounding box input:
    params = {
    'short_name': short_name,
    'version': latest_version,
    'temporal': temporal,
    'page_size': 100,
    'page_num': 1,
    'bounding_box': bounding_box
    }
else:
    
# If polygon input (either via coordinate pairs or shapefile/KML/KMZ):
    params = {
    'short_name': short_name,
    'version': latest_version,
    'temporal': temporal,
    'page_size': 100,
    'page_num': 1,
    'polygon': polygon,
    }

# print('CMR search parameters: ', params)

# Query number of granules using our (paging over results)

granule_search_url = 'https://cmr.earthdata.nasa.gov/search/granules'

granules = []
while True:
    response = requests.get(granule_search_url, params=params, headers=headers)
    results = json.loads(response.content)

    if len(results['feed']['entry']) == 0:
        # Out of results, so break out of loop
        break

    # Collect results and increment page_num
    granules.extend(results['feed']['entry'])
    params['page_num'] += 1

    
# Get number of granules over my area and time of interest
print("Number of granules to download: " + str(len(granules)))

granule_sizes = [float(granule['granule_size']) for granule in granules]

# Average size of granules in MB
print("Mean Granule file size: " + str(mean(granule_sizes)) + " MB")

# Total volume in MB
print("Total file size: " + str(sum(granule_sizes)) + " MB")



# Query service capability URL 

from xml.etree import ElementTree as ET

capability_url = f'https://n5eil02u.ecs.nsidc.org/egi/capabilities/{short_name}.{latest_version}.xml'

# print(capability_url)

# Create session to store cookie and pass credentials to capabilities url

session = requests.session()
s = session.get(capability_url)
response = session.get(s.url,auth=(uid,pswd))

root = ET.fromstring(response.content)

bbox = bounding_box
timevar = start_date + 'T' + start_time + ',' + end_date + 'T' + end_time
# print(timevar)

coverage = '/ancillary_data/atlas_sdp_gps_epoch,\
/gt1l/land_ice_segments/atl06_quality_summary,\
/gt1l/land_ice_segments/delta_time,\
/gt1l/land_ice_segments/h_li,\
/gt1l/land_ice_segments/geophysical,\
/gt1r/land_ice_segments/geophysical,\
/gt2l/land_ice_segments/geophysical,\
/gt2r/land_ice_segments/geophysical,\
/gt3l/land_ice_segments/geophysical,\
/gt3r/land_ice_segments/geophysical,\
/gt1l/land_ice_segments/h_li_sigma,\
/gt1l/land_ice_segments/latitude,\
/gt1l/land_ice_segments/longitude,\
/gt1l/land_ice_segments/segment_id,\
/gt1l/land_ice_segments/sigma_geo_h,\
/gt1r/land_ice_segments/atl06_quality_summary,\
/gt1r/land_ice_segments/delta_time,\
/gt1r/land_ice_segments/h_li,\
/gt1r/land_ice_segments/h_li_sigma,\
/gt1r/land_ice_segments/latitude,\
/gt1r/land_ice_segments/longitude,\
/gt1r/land_ice_segments/segment_id,\
/gt1r/land_ice_segments/sigma_geo_h,\
/gt2l/land_ice_segments/atl06_quality_summary,\
/gt2l/land_ice_segments/delta_time,\
/gt2l/land_ice_segments/h_li,\
/gt2l/land_ice_segments/h_li_sigma,\
/gt2l/land_ice_segments/latitude,\
/gt2l/land_ice_segments/longitude,\
/gt2l/land_ice_segments/segment_id,\
/gt2l/land_ice_segments/sigma_geo_h,\
/gt2r/land_ice_segments/atl06_quality_summary,\
/gt2r/land_ice_segments/delta_time,\
/gt2r/land_ice_segments/h_li,\
/gt2r/land_ice_segments/h_li_sigma,\
/gt2r/land_ice_segments/latitude,\
/gt2r/land_ice_segments/longitude,\
/gt2r/land_ice_segments/segment_id,\
/gt2r/land_ice_segments/sigma_geo_h,\
/gt3l/land_ice_segments/atl06_quality_summary,\
/gt3l/land_ice_segments/delta_time,\
/gt3l/land_ice_segments/h_li,\
/gt3l/land_ice_segments/h_li_sigma,\
/gt3l/land_ice_segments/latitude,\
/gt3l/land_ice_segments/longitude,\
/gt3l/land_ice_segments/segment_id,\
/gt3l/land_ice_segments/sigma_geo_h,\
/gt3r/land_ice_segments/atl06_quality_summary,\
/gt3r/land_ice_segments/delta_time,\
/gt3r/land_ice_segments/h_li,\
/gt3r/land_ice_segments/h_li_sigma,\
/gt3r/land_ice_segments/latitude,\
/gt3r/land_ice_segments/longitude,\
/gt3r/land_ice_segments/segment_id,\
/gt3r/land_ice_segments/sigma_geo_h,\
/orbit_info/cycle_number,\
/orbit_info/rgt,\
/orbit_info/orbit_number' 

base_url = 'https://n5eil02u.ecs.nsidc.org/egi/request'

# Set number of granules requested per order, which we will initially set to 10.
page_size = 10

#Determine number of pages basd on page_size and total granules. Loop requests by this value
page_num = math.ceil(len(granules)/page_size)

#Set request mode. 
request_mode = 'async'

# Determine how many individual orders we will request based on the number of granules requested

# print(page_num)

subset_params = {
    'short_name': short_name, 
    'version': latest_version, 
    'temporal': temporal, 
    'time': timevar, 
    'bounding_box': bounding_box, 
    'bbox': bbox, 
    'Coverage': coverage, 
    'request_mode': request_mode, 
    'page_size': page_size,  
    'token': token, 
    'email': email, 
    }
# # print(subset_params)

# Create and clear target directory
get_ipython().system('mkdir -p Outputs/$download_note')
get_ipython().system('rm -r Outputs/$download_note')

# Create output directory
path = str(os.getcwd() + '/Outputs')
if not os.path.exists(path):
    os.mkdir(path)
    
# Request data service for each page number, and unzip outputs

for i in range(page_num):
    page_val = i + 1
#     print('Order: ', page_val)
    subset_params.update( {'page_num': page_val} )
    
# For all requests other than spatial file upload, use get function
    request = session.get(base_url, params=subset_params)
    
#     print('Request HTTP response: ', request.status_code)

# Raise bad request: Loop will stop for bad response code.
    request.raise_for_status()
#     print('Order request URL: ', request.url)
    esir_root = ET.fromstring(request.content)
#     print('Order request response XML content: ', request.content)

#Look up order ID
    orderlist = []   
    for order in esir_root.findall("./order/"):
        orderlist.append(order.text)
    orderID = orderlist[0]
#     print('order ID: ', orderID)

#Create status URL
    statusURL = base_url + '/' + orderID
#     print('status URL: ', statusURL)

#Find order status
    request_response = session.get(statusURL)    
#     print('HTTP response from order response URL: ', request_response.status_code)
    
# Raise bad request: Loop will stop for bad response code.
    request_response.raise_for_status()
    request_root = ET.fromstring(request_response.content)
    statuslist = []
    for status in request_root.findall("./requestStatus/"):
        statuslist.append(status.text)
    status = statuslist[0]
    print('Data request ', page_val, ' is submitting...')
    print('Initial request status is ', status)

#Continue loop while request is still processing
    while status == 'pending' or status == 'processing': 
        print('Status is not complete. Trying again.')
        time.sleep(10)
        loop_response = session.get(statusURL)

# Raise bad request: Loop will stop for bad response code.
        loop_response.raise_for_status()
        loop_root = ET.fromstring(loop_response.content)

#find status
        statuslist = []
        for status in loop_root.findall("./requestStatus/"):
            statuslist.append(status.text)
        status = statuslist[0]
#         print('Retry request status is: ', status)
        if status == 'pending' or status == 'processing':
            continue

#Order can either complete, complete_with_errors, or fail:
# Provide complete_with_errors error message:
    if status == 'complete_with_errors' or status == 'failed':
        messagelist = []
        for message in loop_root.findall("./processInfo/"):
            messagelist.append(message.text)
#         print('error messages:')
#         pprint.pprint(messagelist)

# Download zipped order if status is complete or complete_with_errors
    if status == 'complete' or status == 'complete_with_errors':
        downloadURL = 'https://n5eil02u.ecs.nsidc.org/esir/' + orderID + '.zip'
        print('Zip download URL: ', downloadURL)
        print('Beginning download of zipped output...')
        zip_response = session.get(downloadURL)
        # Raise bad request: Loop will stop for bad response code.
        zip_response.raise_for_status()
        with zipfile.ZipFile(io.BytesIO(zip_response.content)) as z:
            z.extractall(path)
        print('Data request', page_val, 'is complete.')
    else: print('Request failed.')

# Move data into target directory         
get_ipython().system('mkdir Outputs/$download_note')
get_ipython().system('mv Outputs/*/* Outputs/$download_note/')
get_ipython().system('aws s3 sync Outputs/$download_note/ s3://pangeo-data-upload-oregon/icesat2/Snowblower/')
