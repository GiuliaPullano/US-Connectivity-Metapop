📂 DATA
counties_fips_index.csv
This file provides a mapping of U.S. counties' FIPS codes to the numerical indices used in the metapopulation model and connectivity networks.

📁 US_connectivity_network_2020_county/
This folder contains 12 data files representing the daily average connectivity network between U.S. counties, aggregated at a monthly level. These files capture intercounty mobility patterns and are used to model the spatial spread of diseases in the metapopulation framework.

File Structure
Each file corresponds to a specific month and follows the naming convention:

📄 monthly_matrix_sd_county_monthX.csv

Description: Represents the connectivity network for month X.

Columns:
Column 1: Index of the origin county
Column 2: Index of the destination county
Column 3: Month
Column 4: Connectivity value
Additional File
📄 network_in-out_degree.txt

This file provides in-degree and out-degree metrics for each county in the network.

Columns:
Column 1: Index of the county
Column 2: Month
Column 3: In-degree (number of incoming connections)
Column 4: Out-degree (number of outgoing connections)

