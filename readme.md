## Explore the datasets and preprocessing code for the 3D-PV-Locator paper

This repository provides you with three demo notebooks.

### Explore_Dataset.ipynb

This notebook allows you to visualize the locations of all the images in our classification dataset on an interactive map. 

For example, the locations of non-PV images in our training dataset are illustrated below:

![Exemplary Pipeline Output](https://github.com/kdmayer/3D-PV-Locator/blob/master/exemplary_pipeline_output.png)

### Download_Images_from_OpenNRW.ipynb

This notebook allows you to create your own dataset by downloading images from the openNRW server directly.

### Extract_Rooftop_Information_from_CityGML.ipynb

This notebook demonstrates the code to extract 3D rooftop information from the CityGML files provided by the state of North Rhine-Westphalia ([here](https://www.bezreg-koeln.nrw.de/brk_internet/geobasis/3d_gebaeudemodelle/index.html))

For example, after processing the exemplary .gml files in data/GML/, we can load the extracted rooftop polygons with their respective attributes in QGIS:



## Dependencies

The repository requires the packages GeoPandas, Scipy, Geopy, Pandas, Numpy, and related typical packages used.

## BibTex Citation:

Please cite our work as

    @article{MAYER2022,
    title = {3D-PV-Locator: Large-scale detection of rooftop-mounted photovoltaic systems in 3D},
    journal = {Applied Energy},
    volume = {310},
    pages = {118469},
    year = {2022},
    issn = {0306-2619},
    doi = {https://doi.org/10.1016/j.apenergy.2021.118469},
    url = {https://www.sciencedirect.com/science/article/pii/S0306261921016937},
    author = {Kevin Mayer and Benjamin Rausch and Marie-Louise Arlt and Gunther Gust and Zhecheng Wang and Dirk Neumann and Ram Rajagopal},
    keywords = {Solar panels, Renewable energy, Image recognition, Deep learning, Computer vision, 3D building data, Remote sensing, Aerial imagery},
    }

## License
[MIT](https://choosealicense.com/licenses/mit/)