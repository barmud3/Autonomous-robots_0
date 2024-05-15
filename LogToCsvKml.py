import sys, os, csv
from datetime import datetime, timedelta
import pandas as pd
import numpy as np
import navpy
from gnssutils import EphemerisManager
import simplekml

WEEKSEC = 604800
LIGHTSPEED = 2.99792458e8

parent_directory = os.getcwd()
ephemeris_data_directory = os.path.join(parent_directory, 'data')
sys.path.insert(0, parent_directory)
# Get path to sample file in data directory, which is located in the parent directory of this notebook
input_filepath = os.path.join(parent_directory, 'data', 'sample', 'Walking.txt')

with open(input_filepath) as csvfile:
    reader = csv.reader(csvfile)
    for row in reader:
        if row[0][0] == '#':
            if 'Fix' in row[0]:
                android_fixes = [row[1:]]
            elif 'Raw' in row[0]:
                measurements = [row[1:]]
        else:
            if row[0] == 'Fix':
                android_fixes.append(row[1:])
            elif row[0] == 'Raw':
                measurements.append(row[1:])

android_fixes = pd.DataFrame(android_fixes[1:], columns = android_fixes[0])
measurements = pd.DataFrame(measurements[1:], columns = measurements[0])

# Format satellite IDs
measurements.loc[measurements['Svid'].str.len() == 1, 'Svid'] = '0' + measurements['Svid']
measurements.loc[measurements['ConstellationType'] == '1', 'Constellation'] = 'G'
measurements.loc[measurements['ConstellationType'] == '3', 'Constellation'] = 'R'
measurements['SatPRN (ID)'] = measurements['Constellation'] + measurements['Svid']

# Remove all non-GPS measurements
measurements = measurements.loc[measurements['Constellation'] == 'G']

# Convert columns to numeric representation
measurements['Cn0DbHz'] = pd.to_numeric(measurements['Cn0DbHz'])
measurements['TimeNanos'] = pd.to_numeric(measurements['TimeNanos'])
measurements['FullBiasNanos'] = pd.to_numeric(measurements['FullBiasNanos'])
measurements['ReceivedSvTimeNanos']  = pd.to_numeric(measurements['ReceivedSvTimeNanos'])
measurements['PseudorangeRateMetersPerSecond'] = pd.to_numeric(measurements['PseudorangeRateMetersPerSecond'])
measurements['ReceivedSvTimeUncertaintyNanos'] = pd.to_numeric(measurements['ReceivedSvTimeUncertaintyNanos'])

# A few measurement values are not provided by all phones
# We'll check for them and initialize them with zeros if missing
if 'BiasNanos' in measurements.columns:
    measurements['BiasNanos'] = pd.to_numeric(measurements['BiasNanos'])
else:
    measurements['BiasNanos'] = 0
if 'TimeOffsetNanos' in measurements.columns:
    measurements['TimeOffsetNanos'] = pd.to_numeric(measurements['TimeOffsetNanos'])
else:
    measurements['TimeOffsetNanos'] = 0

measurements['GpsTimeNanos'] = measurements['TimeNanos'] - (measurements['FullBiasNanos'] - measurements['BiasNanos'])
gpsepoch = datetime(1980, 1, 6, 0, 0, 0)
measurements['UnixTime'] = pd.to_datetime(measurements['GpsTimeNanos'], utc = True, origin=gpsepoch)
measurements['UnixTime'] = measurements['UnixTime']

# Split data into measurement epochs
measurements['Epoch'] = 0
measurements.loc[measurements['UnixTime'] - measurements['UnixTime'].shift() > timedelta(milliseconds=200), 'Epoch'] = 1
measurements['Epoch'] = measurements['Epoch'].cumsum()

# This should account for rollovers since it uses a week number specific to each measurement
measurements['tRxGnssNanos'] = measurements['TimeNanos'] + measurements['TimeOffsetNanos'] - (measurements['FullBiasNanos'].iloc[0] + measurements['BiasNanos'].iloc[0])
measurements['GpsWeekNumber'] = np.floor(1e-9 * measurements['tRxGnssNanos'] / WEEKSEC)
measurements['tRxSeconds'] = 1e-9*measurements['tRxGnssNanos'] - WEEKSEC * measurements['GpsWeekNumber']
measurements['tTxSeconds'] = 1e-9*(measurements['ReceivedSvTimeNanos'] + measurements['TimeOffsetNanos'])
# Calculate pseudorange in seconds
measurements['prSeconds'] = measurements['tRxSeconds'] - measurements['tTxSeconds']
# Conver to meters
measurements['PrM'] = LIGHTSPEED * measurements['prSeconds']
measurements['PrSigmaM'] = LIGHTSPEED * 1e-9 * measurements['ReceivedSvTimeUncertaintyNanos']

manager = EphemerisManager(ephemeris_data_directory)

epoch = 0
num_sats = 0
while num_sats < 5 :
    one_epoch = measurements.loc[(measurements['Epoch'] == epoch) & (measurements['prSeconds'] < 0.1)].drop_duplicates(subset='SatPRN (ID)')
    timestamp = one_epoch.iloc[0]['UnixTime'].to_pydatetime(warn=False)
    one_epoch.set_index('SatPRN (ID)', inplace=True)
    num_sats = len(one_epoch.index)
    epoch += 1

sats = one_epoch.index.unique().tolist()
ephemeris = manager.get_ephemeris(timestamp, sats)

# Reorder the columns to include 'SvName' as the first column
one_epoch_selected = one_epoch[['Cn0DbHz', 'UnixTime', 'tTxSeconds']]
one_epoch_selected['tTxSeconds'] = pd.to_timedelta(one_epoch_selected['tTxSeconds'], unit='s')
one_epoch_selected['GPS time'] = one_epoch_selected.apply(lambda row: row['UnixTime'] + row['tTxSeconds'], axis=1)
# Drop the original 'UnixTime' and 'tTxSeconds' columns
one_epoch_selected.drop(['UnixTime', 'tTxSeconds'], axis=1, inplace=True)
# Specify the output file path for the selected dataframe
output_filepath = os.path.join(parent_directory, 'data', 'sample', 'output4_measurements.csv')


def least_squares(xs, measured_pseudorange, x0, b0):
    dx = 100*np.ones(3)
    b = b0
    # set up the G matrix with the right dimensions. We will later replace the first 3 columns
    # note that b here is the clock bias in meters equivalent, so the actual clock bias is b/LIGHTSPEED
    G = np.ones((measured_pseudorange.size, 4))
    iterations = 0
    while np.linalg.norm(dx) > 1e-3:
        r = np.linalg.norm(xs - x0, axis=1)
        phat = r + b0
        deltaP = measured_pseudorange - phat
        G[:, 0:3] = -(xs - x0) / r[:, None]
        sol = np.linalg.inv(np.transpose(G) @ G) @ np.transpose(G) @ deltaP
        dx = sol[0:3]
        db = sol[3]
        x0 = x0 + dx
        b0 = b0 + db
    norm_dp = np.linalg.norm(deltaP)
    return x0, b0, norm_dp

def ecef_to_lla(location):
    lla = navpy.ecef2lla(location)  # Specify latitude and longitude units
    latitude = lla[0]
    longitude = lla[1]
    altitude = lla[2]
    return latitude, longitude, altitude

def calculate_satellite_position(ephemeris, transmit_time):
    mu = 3.986005e14
    OmegaDot_e = 7.2921151467e-5
    F = -4.442807633e-10
    sv_position = pd.DataFrame()
    sv_position['sv']= ephemeris.index
    sv_position.set_index('sv', inplace=True)
    sv_position['t_k'] = transmit_time - ephemeris['t_oe']
    A = ephemeris['sqrtA'].pow(2)
    n_0 = np.sqrt(mu / A.pow(3))
    n = n_0 + ephemeris['deltaN']
    M_k = ephemeris['M_0'] + n * sv_position['t_k']
    E_k = M_k
    err = pd.Series(data=[1]*len(sv_position.index))
    i = 0
    while err.abs().min() > 1e-8 and i < 10:
        new_vals = M_k + ephemeris['e']*np.sin(E_k)
        err = new_vals - E_k
        E_k = new_vals
        i += 1
        
    sinE_k = np.sin(E_k)
    cosE_k = np.cos(E_k)
    delT_r = F * ephemeris['e'].pow(ephemeris['sqrtA']) * sinE_k
    delT_oc = transmit_time - ephemeris['t_oc']
    sv_position['delT_sv'] = ephemeris['SVclockBias'] + ephemeris['SVclockDrift'] * delT_oc + ephemeris['SVclockDriftRate'] * delT_oc.pow(2)

    v_k = np.arctan2(np.sqrt(1-ephemeris['e'].pow(2))*sinE_k,(cosE_k - ephemeris['e']))

    Phi_k = v_k + ephemeris['omega']

    sin2Phi_k = np.sin(2*Phi_k)
    cos2Phi_k = np.cos(2*Phi_k)

    du_k = ephemeris['C_us']*sin2Phi_k + ephemeris['C_uc']*cos2Phi_k
    dr_k = ephemeris['C_rs']*sin2Phi_k + ephemeris['C_rc']*cos2Phi_k
    di_k = ephemeris['C_is']*sin2Phi_k + ephemeris['C_ic']*cos2Phi_k

    u_k = Phi_k + du_k

    r_k = A*(1 - ephemeris['e']*np.cos(E_k)) + dr_k

    i_k = ephemeris['i_0'] + di_k + ephemeris['IDOT']*sv_position['t_k']

    x_k_prime = r_k*np.cos(u_k)
    y_k_prime = r_k*np.sin(u_k)

    Omega_k = ephemeris['Omega_0'] + (ephemeris['OmegaDot'] - OmegaDot_e)*sv_position['t_k'] - OmegaDot_e*ephemeris['t_oe']
    pr = one_epoch['PrM'] + LIGHTSPEED * sv_position['delT_sv']
    pr = pr.to_numpy()
    sv_position['Psuedo-range'] = pr

    sv_position['x_k'] = x_k_prime*np.cos(Omega_k) - y_k_prime*np.cos(i_k)*np.sin(Omega_k)
    sv_position['y_k'] = x_k_prime*np.sin(Omega_k) + y_k_prime*np.cos(i_k)*np.cos(Omega_k)
    sv_position['z_k'] = y_k_prime*np.sin(i_k)
    
    return sv_position


# Run the function and check out the results:

sv_position = calculate_satellite_position(ephemeris, one_epoch['tTxSeconds'])

xs = sv_position[['x_k', 'y_k', 'z_k']].to_numpy()
pr = sv_position['Psuedo-range']
b0 = 0
x0 = np.array([0, 0, 0])

x, b, dp = least_squares(xs, pr, x0, b0)
    
sv_position['Pos.X'] = x[0]
sv_position['Pos.Y'] = x[1]
sv_position['Pos.Z'] = x[2]
pos_lla = ecef_to_lla(x)
sv_position['Lat'] ,sv_position['Lon'] , sv_position['Alt'] = pos_lla[0] , pos_lla[1] , pos_lla[2]
sv_position.drop(columns=['delT_sv'], inplace=True)
sv_position.drop(columns=['t_k'], inplace=True)
# Assuming you've already run the first and second code snippets

# Merge the dataframes based on the index
merged_df = pd.merge(one_epoch_selected, sv_position, left_index=True, right_index=True)

# Specify the output file path for the merged dataframe
merged_output_filepath = os.path.join(parent_directory, 'CsvFile.csv')

merged_df = merged_df.rename(columns={'x_k': 'Sat.X'})
merged_df = merged_df.rename(columns={'y_k': 'Sat.Y'})
merged_df = merged_df.rename(columns={'z_k': 'Sat.Z'})
merged_df = merged_df.rename(columns={'Cn0DbHz': 'CN0'})

# Save the merged dataframe to a CSV file

merged_df.to_csv(merged_output_filepath, index=True)


def create_kml_from_lla(lla_df, output_filepath):
    """
    Create a KML file containing points from the given DataFrame of LLA coordinates.

    Args:
        lla_df (pd.DataFrame): DataFrame containing columns for Latitude, Longitude, and Altitude.
        output_filepath (str): Path to save the output KML file.
    """
    kml = simplekml.Kml()

    for index, row in lla_df.iterrows():
        point = kml.newpoint(
            name=f"Point {index}",
            coords=[(row['Longitude'], row['Latitude'], row['Altitude'])]  # (lon, lat, alt)
        )
        point.style.iconstyle.icon.href = 'http://maps.google.com/mapfiles/kml/shapes/placemark_circle.png'

    kml.save(output_filepath)
    print(f"KML file saved at: {output_filepath}")


b0 = 0
x0 = np.array([0, 0, 0])
xs = sv_position[['x_k', 'y_k', 'z_k']].to_numpy()

# Apply satellite clock bias to correct the measured pseudorange values
pr = sv_position["Psuedo-range"]
x, b, dp = least_squares(xs, pr, x0, b0)

ecef_list = []
for epoch in measurements['Epoch'].unique():
    one_epoch = measurements.loc[(measurements['Epoch'] == epoch) & (measurements['prSeconds'] < 0.1)] 
    one_epoch = one_epoch.drop_duplicates(subset='SatPRN (ID)').set_index('SatPRN (ID)')
    if len(one_epoch.index) > 4:
        timestamp = one_epoch.iloc[0]['UnixTime'].to_pydatetime(warn=False)
        sats = one_epoch.index.unique().tolist()
        ephemeris = manager.get_ephemeris(timestamp, sats)
        sv_position = calculate_satellite_position(ephemeris, one_epoch['tTxSeconds'])

        xs = sv_position[['x_k', 'y_k', 'z_k']].to_numpy()
        pr = sv_position['Psuedo-range']
        x, b, dp = least_squares(xs, pr, x, b)
        ecef_list.append(x)

# Perform coordinate transformations using the Navpy library

# convert x y z position to Latitude, Longitude, Altitude
ecef_array = np.stack(ecef_list, axis=0)
lla_array = np.stack(ecef_to_lla(ecef_array), axis=1)
ref_lla = lla_array[0, :]

# Create DataFrames
lla_df = pd.DataFrame(lla_array, columns=['Latitude', 'Longitude', 'Altitude'])

create_kml_from_lla(lla_df, 'KmlFile.kml')

