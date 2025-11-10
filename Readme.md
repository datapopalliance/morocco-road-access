# 🛣️ Road Access and Isolation Analysis in Morocco

This project estimates the **percentage of isolated territories in Morocco** by analyzing the relationship between **road networks**, **population distribution**, and **accessibility** to key services or population centers.

---

## 🌍 Overview

Transportation accessibility is a crucial factor in equitable development and service delivery. Using **OpenStreetMap**, **WorldPop**, and **OpenRouteService**, this project identifies and quantifies territories in Morocco that experience poor road connectivity or limited access to essential services.

---

## 🎯 Objectives
- Map and analyze the **road network** across Morocco.
- Compute **isochrones** (travel time areas) using OpenRouteService.
- Overlay **population data** to estimate the share of isolated inhabitants.
- Quantify the **percentage of isolated territories** (e.g., beyond 60-minute access to the nearest major center).

---

## 🗂️ Project Structure
P144-morocco-road-access/  
├── Data/ # Input geoparquet, or road network data  
├── output/ # output folder  
├── Road_Access.ipynb/ # Jupyter notebooks for analysis and map export  
├── src/ # Python scripts and utilities  
├── requirements.txt # Environment dependencies   
└── README.md # Project documentation  
└── ors-config.yml # Configuration to set-up ORS Local Server  

## 🧠 Methodology
1. **Data Collection:** Import road and territorial data (e.g., from OpenStreetMap or Moroccan government sources).
2. **Processing:** Use tools like `GeoPandas`, `Shapely`, and `NetworkX` to analyze accessibility.
3. **Analysis:** Compute isolation metrics such as travel time to nearest all weather roads.
4. **Visualization:** Plot maps and charts showing isolated areas.

## Installation

Code implemented in Python 3.12

### Setting up environment

Clone and go to repository

```
$ git clone https://github.com/datapopalliance/P144-morocco-road-access.git
$ cd morocco-road-access
```

Create, activate environment and Install dependancies

```
$ (new_env) pip install -r requirements.txt
```

## Getting Started

Install local OpenRouteService on computer following instructions from [here](https://giscience.github.io/openrouteservice/run-instance/running-with-docker).

For this specific case, after installing [Docker](https://docs.docker.com/get-docker/) on your configuration, users need to download an extract of OpenStreetMap data on their study area. For Morocco, the extract can be found [here](https://download.geofabrik.de/africa/morocco.html). The OSM extract needs to be placed within an ORS folder along the "ors-config.yml" file provided in this reposirory. 

The server can be configured like this (after the replacing the path)

```
docker run -it --name ors-app -p 8080:8082   -v /home/mikmeh01/ors/ors-config.yml:/home/ors/config/ors-config.yml   -v /home/mikmeh01/ors/ors-docker/elevation_cache:/home/ors/elevation_cache   -v /home/mikmeh01/ors/ors-docker/graphs:/home/ors/graphs   -v /home/mikmeh01/ors/ors-docker/logs:/home/ors/logs   -v /home/mikmeh01/ors/morocco-latest.osm.pbf:/home/ors/morocco-latest.osm.pbf   -e "JAVA_OPTS=-Xms4g -Xmx8g"   openrouteservice/openrouteservice:latest
```

and then run from Jupyter file "road_access.ipynb"

Here is an example of **[Interactive Map](https://rawcdn.githack.com/datapopalliance/morocco-road-access/8b0fc0cc6cecfe5b6c22903c33e3cfe3c85f8c95/enclavement_commune_50.html)** for Morocco that can be generated