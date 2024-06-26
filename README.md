# Project Overview:

This code is designed to process GPS measurements and generate key information such as satellite positions and pseudorange corrections. 
It primarily focuses on extracting and processing raw GNSS measurements from Android devices, performing necessary conversions and calculations, and then exporting the results for further analysis or visualization.

## Dependencies:

* Python 3.x
* Required Libraries: sys, os, csv, datetime, pandas, numpy, navpy, gnssutils, simplekml

#### - Ensure that all required dependencies are installed. You can install missing dependencies using pip:

 ```` pip install pandas numpy navpy gnssutils simplekml ````

## Useage
* Place your GNSS data file (e.g., driving.txt) in data -> sample -> drop your raw file here.
*  Modify the script's main function to specify the input file path.
*  Ensure that the file contains the necessary information, including Fix and Raw data. Open the Python script containing the code in your preferred Python environment.
* The script will process the GNSS data, calculate satellite and receiver positions, and generate output files including CSV files containing satellite and receiver positions, as well as a KML file.
  
  **How to run:**
```` python LogToCsvKml.py ````

## Data Preprocessing
The code begins by reading the input file and extracting Fix and Raw data into Pandas DataFrames for further processing. It then cleans and formats the data, ensuring consistency and correctness.
Satellite Position Calculation: The code utilizes ephemeris data to calculate satellite positions for each measurement epoch. It applies necessary corrections and least squares optimization to improve accuracy.
Coordinate Transformation: Using the Navpy library, the code transforms satellite positions from ECEF (Earth-Centered, Earth-Fixed) coordinates to Latitude, Longitude, and Altitude (LLA) coordinates.
Output Generation: The code exports the processed data and results into CSV files for further analysis. Additionally, it generates a KML file containing plotted points based on the calculated LLA coordinates.

## Explanation of the code

First of all, the code reads GNSS (Global Navigation Satellite System) data from a CSV file, separates it into Android fixes and raw measurements, processes and filters the data, performs time-related calculations,
selects satellite measurements based on criteria such as signal strength and pseudorange, retrieves ephemeris data for the satellites selected, and saves the processed data to a CSV file.
Essentially, it prepares GNSS data for further analysis by formatting, filtering and organizing it.

### Functionality
#### least_squares
This function implements the least squares method to estimate the receiver's position (x0) and clock bias (b0) based on measured pseudorange data.
It iteratively updates these estimates until the difference between predicted and measured pseudorange values converges.
The function returns the estimated position, clock bias, and the norm of the difference between predicted and measured pseudorange values.
#### ecef_to_lla
This function converts Earth-Centered, Earth-Fixed (ECEF) coordinates (x, y, z) to Latitude, Longitude, and Altitude (LLA) coordinates.
It uses the NavPy library to perform the conversion and returns the result as a Pandas DataFrame with columns for Latitude, Longitude, and Altitude.
If the input ECEF coordinates are one-dimensional, it reshapes the output LLA array accordingly before converting it to a DataFrame.
#### calculate_satellite_position
This function calculates the position of satellites based on ephemeris data and transmit time.
It involves a series of calculations, including solving for eccentric anomaly, computing satellite clock corrections, determining satellite position in ECEF coordinates, performing least squares estimation to refine the position,
and finally converting the ECEF coordinates to Latitude, Longitude, and Altitude (LLA) coordinates. The resulting satellite positions are returned as a DataFrame containing various satellite information, including pseudorange, position in ECEF and LLA coordinates, and more.
#### create_kml_from_lla
This function defines a function `create_kml_from_lla` that generates a KML file from a DataFrame containing Latitude, Longitude, and Altitude (LLA) coordinates.
It uses the SimpleKML library to create KML objects for each point in the DataFrame and saves the resulting KML file to the specified output filepath.
Finally, it prints a message confirming the file's creation. When called with appropriate arguments, it creates a KML file named 'points.kml' containing the plotted points.

### Tests
This repository includes a suite of tests to ensure the accuracy and reliability of the implemented algorithms. The tests cover key functionalities such as satellite position estimation and coordinate conversion.

#### Running Tests
To run the tests, execute the TestAlgorithm class located in the LogToCsvKml module using the unittest framework. Ensure that all dependencies are installed before running the tests.

``` python -m unittest LogToCsvKml.TestAlgorithm ```
#### Test Results
Successful execution of all tests indicates that the algorithms are likely functioning correctly.
Test failures suggest potential issues with the implementation that require further investigation and debugging.

## Output Files 
CsvFile.csv
KmlFile.kml
### Examples
#### Kml File : 
![image](https://github.com/barmud3/Autonomous-robots_0/assets/130641348/50a6bfad-1111-4c4d-911f-82d2be6e6a59)

#### Csv File : 
![image](https://github.com/barmud3/Autonomous-robots_0/assets/130641348/f6a53788-fe8a-47c5-ae88-49b7a06b4e09)

## Crew members 
* Bar Alayof - 206840621
* Avihay Finish - 208907113
* Amit Rovshitz - 207701426
* Shoval Zohar - 318284668
